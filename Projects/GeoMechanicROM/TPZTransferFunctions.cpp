//
//  TPZTransfer.cpp
//  PZ
//
//  Created by Omar on 3/1/17.
//
//

#include "TPZTransferFunctions.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


/** @brief Default constructor */
TPZTransferFunctions::TPZTransferFunctions(){
    
    fSimulationData = NULL;
    
}

/** @brief Default desconstructor */
TPZTransferFunctions::~TPZTransferFunctions(){
    
}

/** @brief Copy constructor $ */
TPZTransferFunctions::TPZTransferFunctions(const TPZTransferFunctions &copy)
{
    fSimulationData = copy.fSimulationData;
}

/** @brief Copy assignemnt operator $ */
TPZTransferFunctions & TPZTransferFunctions::operator=(const TPZTransferFunctions &other)
{
    if (this != & other) // prevent self-assignment
    {
        fSimulationData = other.fSimulationData;
    }
    return *this;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Matrices Initialization Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TPZTransferFunctions::Initialize_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || !fCmeshRB_projections || fgeomechanic_galerkinp_cindexes.size() == 0) {
        std::cout << "There is no computational meshes, cmesh_multiphysics||fCmeshRB_projections = Null." << std::endl;
        DebugStop();
    }
#endif
    
    
    fCmeshRB_projections->LoadReferences();
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    // Compute destination index scatter by element (Omega and Gamma)
    fphi_u_dof_scatter.Resize(nel);
    
    // Block size structue including (Omega and Gamma)
    TPZVec< std::pair<long, long> > blocks_dimensions_phi(nel);
    TPZVec< std::pair<long, long> > blocks_dimensions_grad_phi(nel);
    
    long gel_index, geomech_index, gp_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;

        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(geomec_cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = gp_cel->Index();
        
        TPZInterpolationSpace * gp_intel = dynamic_cast<TPZInterpolationSpace * >(gp_cel);
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        int gel_dim = gel->Dimension();
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        this->ElementDofIndexes(gp_intel, dof_indexes);
        fphi_u_dof_scatter[element_index] = dof_indexes;
        blocks_dimensions_phi[element_index].first = int_point_indexes.size()*n_var_dim;
        blocks_dimensions_phi[element_index].second = dof_indexes.size();
        
        blocks_dimensions_grad_phi[element_index].first = int_point_indexes.size()*n_var_dim*gel_dim;
        blocks_dimensions_grad_phi[element_index].second = dof_indexes.size();
    }
    
    // Initialize the matrix
    fphi_u_To_Geomechanic.Initialize(blocks_dimensions_phi);
    fgrad_phi_u_To_Geomechanic.Initialize(blocks_dimensions_grad_phi);
    
}


void TPZTransferFunctions::Fill_RB_basis_To_Geomechanic(TPZCompMesh * cmesh_multiphysics){
    
    // It verify the consistency of dynamic_cast operations and mesh structure, and  finally it initialize diagonal matrix blocks
    Initialize_RB_basis_To_Geomechanic(cmesh_multiphysics);
    
    
    fCmeshRB_projections->LoadReferences();
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int n_var_dim = fCmeshRB_projections->Reference()->Dimension(); // vectorial
    long element_index = 0;
    
    std::pair<long, long> block_dim;
    
    long gel_index, geomech_index, gp_index;
    for (long iel = 0; iel < nel; iel++) {
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        

        
        TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(geomec_cel);
#ifdef PZDEBUG
        if(!mf_cel)
        {
            DebugStop();
        }
#endif
        element_index = gp_cel->Index();
        TPZInterpolationSpace * gp_intel = dynamic_cast<TPZInterpolationSpace * >(gp_cel);
        
        // Getting local integration index
        TPZManVector<long> int_point_indexes(0,0);
        TPZManVector<long> dof_indexes(0,0);
        
        mf_cel->GetMemoryIndices(int_point_indexes);
        dof_indexes = fphi_u_dof_scatter[element_index];
        
        block_dim.first = int_point_indexes.size();
        block_dim.second = dof_indexes.size();
        
        // Computing the local integration points indexes
        const TPZIntPoints & int_points_geomechanic = mf_cel->GetIntegrationRule();
        int np_cel = int_points_geomechanic.NPoints();
        
#ifdef PZDEBUG
        if (int_point_indexes.size() != np_cel) {
            DebugStop();
        }
#endif
        
        // Computing over all integration points of the compuational element cel
        int gel_dim = gel->Dimension();
        TPZFMatrix<double> block_phi, block_grad_phi;
        
        block_phi.Resize(block_dim.first*n_var_dim,block_dim.second);
        block_grad_phi.Resize(block_dim.first*n_var_dim*gel_dim,block_dim.second);

        // for derivatives in real space
        int nshape = gp_intel->NShapeF();
        TPZFNMatrix<220> phi(nshape,1);
        TPZFNMatrix<660> dphi(gel_dim,nshape),dphix_axes(gel_dim,nshape);
        TPZFMatrix<double> dphidx;
        TPZFNMatrix<9,STATE> jacobian(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> jacinv(gel_dim,gel_dim);
        TPZFNMatrix<9,STATE> axes;
        REAL detjac;
        
        for (int ip = 0; ip < block_dim.first ; ip++)
        {
            TPZManVector<REAL,3> qsi(gel_dim,0.0);
            STATE w;
            int_points_geomechanic.Point(ip, qsi, w);
            
            // Get the phi and dphix for H1 elasticity
            gp_intel->Shape(qsi, phi, dphi);
            gel->Jacobian( qsi, jacobian, axes, detjac , jacinv);
            
            switch(gel_dim) {
                case 0:
                    break;
                case 1:
                    dphix_axes = dphi;
                    dphix_axes *= (1./detjac);
                    break;
                case 2:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq);
                    }
                    break;
                case 3:
                    for(int ieq = 0; ieq < nshape; ieq++) {
                        dphix_axes(0,ieq) = jacinv(0,0)*dphi(0,ieq) + jacinv(1,0)*dphi(1,ieq) + jacinv(2,0)*dphi(2,ieq);
                        dphix_axes(1,ieq) = jacinv(0,1)*dphi(0,ieq) + jacinv(1,1)*dphi(1,ieq) + jacinv(2,1)*dphi(2,ieq);
                        dphix_axes(2,ieq) = jacinv(0,2)*dphi(0,ieq) + jacinv(1,2)*dphi(1,ieq) + jacinv(2,2)*dphi(2,ieq);
                    }
                    break;
                default:
                    std::stringstream sout;
                    sout << "pzintel.c please implement the " << gel_dim << "d Jacobian and inverse\n";
                    LOGPZ_ERROR(logger,sout.str());
            }
            
            TPZAxesTools<STATE>::Axes2XYZ(dphix_axes, dphidx, axes);
            
#ifdef PZDEBUG
            if(block_dim.second != phi.Rows() * n_var_dim){
                DebugStop();
            }
#endif
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    block_phi(ip*n_var_dim+id,jp*n_var_dim+id) = phi(jp,0);
                }
            }
            
            for (int jp = 0; jp < phi.Rows(); jp++) {
                for (int id = 0; id < n_var_dim; id++) {
                    for (int jd = 0; jd < gel_dim; jd++) {
                        if(gel_dim == n_var_dim){
                            block_grad_phi(ip*n_var_dim*gel_dim + id*gel_dim + jd,jp*n_var_dim+id) = dphidx(jd,jp);
                        }
                        else{
                            block_grad_phi(ip*n_var_dim*gel_dim+id*gel_dim + jd,jp*n_var_dim+id) = 0.0;
                        }
                    }
                }
            }
            
        }
        
        fphi_u_To_Geomechanic.SetBlock(element_index, block_phi);
        fgrad_phi_u_To_Geomechanic.SetBlock(element_index, block_grad_phi);
    }
    
}

void TPZTransferFunctions::RB_basis_To_Geomechanic_Memory(TPZCompMesh * cmesh_multiphysics){
    
    
    this->Fill_RB_basis_To_Geomechanic(cmesh_multiphysics);
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || fgeomechanic_galerkinp_cindexes.size() == 0) {
        DebugStop();
    }
#endif
    

    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int dim = cmesh_multiphysics->Dimension();
    int n_rb = fCmeshRB_projections->Solution().Cols();
    

    // Step zero scatter
    TPZFMatrix<STATE> ScatterDisplacements(fphi_u_To_Geomechanic.Cols(),n_rb,0.0);

    //    for all RB functions
    for(int i_rb = 0; i_rb < n_rb; i_rb++){
        long pos = 0;
        for (int el = 0; el < nel; el++) {
            for(int ip = 0; ip < fphi_u_dof_scatter[el].size(); ip++) {
                ScatterDisplacements(pos,i_rb) = fCmeshRB_projections->Solution()(fphi_u_dof_scatter[el][ip],i_rb);
                pos++;
            }
        }
        
    }
    
    
    // Step two
    TPZFMatrix<STATE> u_at_intpoints,grad_u_at_intpoints;
    fphi_u_To_Geomechanic.Multiply(ScatterDisplacements,u_at_intpoints);
    fgrad_phi_u_To_Geomechanic.Multiply(ScatterDisplacements, grad_u_at_intpoints);
    
    
    long iblock = 0;
    long first_point_phi = 0;
    long first_point_dphi = 0;
    std::pair<long, long> b_size_phi, b_size_dphi;
    b_size_phi.first = 0;
    b_size_phi.second = 0;
    b_size_dphi.first = 0;
    b_size_dphi.second = 0;
    
    long gel_index, geomech_index, gp_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        first_point_phi += b_size_phi.first;
        first_point_dphi += b_size_dphi.first;
        b_size_phi = fphi_u_To_Geomechanic.GetSizeofBlock(iblock);
        b_size_dphi = fgrad_phi_u_To_Geomechanic.GetSizeofBlock(iblock);
        iblock++;
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> *>(material);
            
            long np_cmesh = associated_material->GetMemory().NElements();

            for(long i = 0; i <  np_cmesh; i++) {
                associated_material->GetMemory()[i].phi_u().Resize(n_rb,dim);
                associated_material->GetMemory()[i].grad_phi_u().Resize(n_rb,dim*dim);
            }
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);

            // Transfering values
            int n_points = int_point_indexes.size();
            TPZFMatrix<STATE> phi_u(n_rb,dim,0.0);
            TPZFMatrix<STATE> grad_phi_u(n_rb,dim*dim);
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                
                //    for all RB functions
                for(int i_rb = 0; i_rb < n_rb; i_rb++){
                    for (int id = 0; id < dim ; id++) {
                        phi_u(i_rb,id) = u_at_intpoints(first_point_phi + ip*dim + id,i_rb);
                    }
                    
                    int c_pos = 0;
                    for (int id = 0; id < dim ; id++) {
                        for (int jd = 0; jd < dim ; jd++) {
                            grad_phi_u(i_rb,c_pos)= grad_u_at_intpoints(first_point_dphi + ip*dim*dim + id*dim + jd,i_rb);
                            c_pos++;
                        }
                    }
                }
                associated_material->GetMemory()[int_point_indexes[ip]].Set_phi_u(phi_u);
                associated_material->GetMemory()[int_point_indexes[ip]].Set_grad_phi_u(grad_phi_u);
            }
            
        }
        else{
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> *>(material);
            
            long np_cmesh = associated_material->GetMemory().NElements();
            
            for(long i = 0; i <  np_cmesh; i++) {
                associated_material->GetMemory()[i].phi_u().Resize(n_rb,dim);
            }
            
            TPZManVector<long, 30> int_point_indexes;
            cmesh_multiphysics->LoadReferences();
            TPZCompEl * mf_cel = gel->Reference();
            GlobalPointIndexes(mf_cel, int_point_indexes);
            
            // Transfering values
            int n_points = int_point_indexes.size();
            TPZFMatrix<STATE> phi_u(n_rb,dim,0.0);
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                
                //    for all RB functions
                for(int i_rb = 0; i_rb < n_rb; i_rb++){
                    for (int id = 0; id < dim ; id++) {
                        phi_u(i_rb,id) = u_at_intpoints(first_point_phi + ip*dim + id,i_rb);
                    }
                }
                associated_material->GetMemory()[i_pos].Set_phi_u(phi_u);
            }
            
        }

    }
    
}

/** @brief Transfer the RB Solution to multiphysics mesh  */
void TPZTransferFunctions::RB_Solution_To_Geomechanic(TPZCompMesh * cmesh_multiphysics, TPZFMatrix<STATE> & rb_solution){
    
#ifdef PZDEBUG
    if (!cmesh_multiphysics || fgeomechanic_galerkinp_cindexes.size() == 0) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    long nel = fgeomechanic_galerkinp_cindexes.size();
    int dim = cmesh_multiphysics->Dimension();
    int n_gp = fCmeshRB_projections->Solution().Cols();
    
#ifdef PZDEBUG
    if (n_gp != rb_solution.Rows()) {
        DebugStop();
    }
#endif
    
    
    std::pair<long, long> block_size;
    block_size.first = 0;
    block_size.second = 0;
    
    
    long gel_index, geomech_index, gp_index;
    for (int iel = 0; iel < nel; iel++) {
        
        
        gel_index = fgeomechanic_galerkinp_cindexes[iel].first;
        TPZGeoEl * gel = geometry->Element(gel_index);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        
        if (gel->HasSubElement()) {
            continue;
        }
        
        geomech_index   = fgeomechanic_galerkinp_cindexes[iel].second.first;
        gp_index        = fgeomechanic_galerkinp_cindexes[iel].second.second;
        
        TPZCompEl * geomec_cel = cmesh_multiphysics->Element(geomech_index);
        TPZCompEl * gp_cel = fCmeshRB_projections->Element(gp_index);
        
#ifdef PZDEBUG
        if (!geomec_cel || !gp_cel) {
            DebugStop();
        }
#endif
        
        
        //  Getting the total integration point of the destination cmesh
        int matd_id = gel->MaterialId();
        if(matd_id == 1){ // The volumetric ones!
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZDiscontinuousGalerkin> *>(material);
            
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);
            
            // Transfering values
            long n_points = int_point_indexes.size();
            TPZFNMatrix<3,REAL> u(dim,1);
            TPZFNMatrix<9,REAL> grad_u(dim,dim);
            
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                u.Zero();
                grad_u.Zero();

                 for (int igp = 0; igp <n_gp; igp++) {
                     int c_pos = 0;
                     for (int id = 0; id <dim; id++) {
                        u(id,0) += associated_material->GetMemory()[i_pos].phi_u()(igp,id) * rb_solution(igp,0);
                        for (int jd = 0; jd <dim; jd++) {
                            grad_u(id,jd) += associated_material->GetMemory()[i_pos].grad_phi_u()(igp,c_pos) * rb_solution(igp,0);
                            c_pos++;
                        }
                    }
                }
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[i_pos].Set_u_n(u);
                    associated_material->GetMemory()[i_pos].Set_grad_u_n(grad_u);
                }
                else{
                    associated_material->GetMemory()[i_pos].Set_u(u);
                    associated_material->GetMemory()[i_pos].Set_grad_u(grad_u);                    
                }
            }
            
        }
        else{
            
            TPZMaterial * material = cmesh_multiphysics->FindMaterial(matd_id);
            TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> * associated_material = dynamic_cast<TPZMatWithMem<TPZPoroPermMemory,TPZBndCond> *>(material);
            
            TPZManVector<long, 30> int_point_indexes;
            GlobalPointIndexes(geomec_cel, int_point_indexes);
            
            // Transfering values
            long n_points = int_point_indexes.size();
            TPZFNMatrix<3,REAL> u(dim,1);
            
            long i_pos;
            for(long ip = 0; ip <  n_points; ip++){
                i_pos = int_point_indexes[ip];
                u.Zero();
                for (int igp = 0; igp <n_gp; igp++) {
                    int c_pos = 0;
                    for (int id = 0; id <dim; id++) {
                        u(id,0) += associated_material->GetMemory()[i_pos].phi_u()(igp,id) * rb_solution(igp,0);
                    }
                }
                if (fSimulationData->IsCurrentStateQ()) {
                    associated_material->GetMemory()[i_pos].Set_u_n(u);
                }
                else{
                    associated_material->GetMemory()[i_pos].Set_u(u);
                }
            }
            
        }
        
    }
    
    
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Utility Methods
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Get Global integration point indexes associaded  */
void TPZTransferFunctions::GlobalPointIndexes(TPZCompEl * cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsElement * mf_cel = dynamic_cast<TPZMultiphysicsElement * >(cel);
    
#ifdef PZDEBUG
    if(!mf_cel)
    {
        DebugStop();
    }
#endif
    
    mf_cel->GetMemoryIndices(int_point_indexes);
    
}

/** @brief Get Global integration point indexes associaded  */
void TPZTransferFunctions::GlobalPointIndexesInterface(TPZCompEl * int_cel, TPZManVector<long,30> &int_point_indexes){
    
    TPZMultiphysicsInterfaceElement * mf_int_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(int_cel);
    
#ifdef PZDEBUG
    if(!mf_int_cel)
    {
        DebugStop();
    }
#endif
    
    mf_int_cel->GetMemoryIndices(int_point_indexes);
    
}

bool TPZTransferFunctions::IdentifyFace(int &side, TPZGeoEl * vol, TPZGeoEl * face){
    
    int volu_nsides = vol->NSides();
    int face_nsides = face->NSides();
    side = -1;
    TPZGeoElSide face_itself =  face->Neighbour(face_nsides-1);
    bool IsMembershipQ = false;
    
    for (int iside = 0; iside < volu_nsides; iside++) {
        IsMembershipQ = bool(vol->NeighbourExists(iside, face_itself));
        if (IsMembershipQ) {
            side = iside;
            break;
        }
    }
    
    TPZGeoElSide vol_itself =  vol->Neighbour(volu_nsides-1);
    
    if(!IsMembershipQ){
        TPZGeoElSide neigh = face->Neighbour(face_nsides-1);
        
        if(!neigh.Element()){
            DebugStop();
        }
        
        TPZGeoElSide neigh_father = neigh.Father2();
        
        if(!neigh_father.Element()){
            DebugStop();
        }
        bool IsNeighQ = false;
        for (int iside = vol_itself.Element()->NNodes(); iside < volu_nsides; iside++) {
            
            IsNeighQ = bool(vol_itself.Element()->NeighbourExists(iside, neigh_father));
            
            if (IsNeighQ) {
                side = iside;
                IsMembershipQ = IsNeighQ;
                break;
            }
        }
    }
    
    return IsMembershipQ;
}

/** @brief Compute parametric transformation form origin to tarfet (xinverse based) */
void TPZTransferFunctions::ComputeTransformation(TPZGeoEl * face_gel_origin, TPZGeoEl * gel_origin , TPZGeoEl * gel_target, TPZVec<REAL> & origin, TPZVec<REAL> & target){
    
#ifdef PZDEBUG
    if (!face_gel_origin || !gel_origin || !gel_target) {
        DebugStop();
    }
    
    if (gel_origin->Dimension() != gel_target->Dimension()) {
        DebugStop();
    }
    
#endif
    
    bool IsmemberQ;
    int face_side;
    TPZManVector<REAL,3> origin_vol(gel_origin->Dimension());
    target.Resize(gel_origin->Dimension(),0.0);
    
    // Transformation step one
    TPZTransform<REAL> tr_face_to_vol(gel_origin->Dimension());
    IsmemberQ =  IdentifyFace(face_side, gel_origin, face_gel_origin);
    
    
#ifdef PZDEBUG
    if (!IsmemberQ) {
        DebugStop();
    }
#endif
    
    tr_face_to_vol = gel_origin->SideToSideTransform(face_side, gel_origin->NSides()-1);
    tr_face_to_vol.Apply(origin, origin_vol);
    gel_origin->TransformSonToFather(gel_target, origin_vol, target);
    
    //    std::cout << "origin = "        << origin << std::endl;
    //    std::cout << "origin_vol = "    << origin_vol << std::endl;
    //    std::cout << "target = "        << target << std::endl;
    
}

/** @brief Compute indices associated to faces on 3D topologies */
void TPZTransferFunctions::ComputeFaceIndex(TPZGeoEl * gel , TPZVec<int> &sides){
    
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            
        }
            break;
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute sides associated to faces on 3D topologies */
void TPZTransferFunctions::ComputeFaceNormals(TPZGeoEl * gel , TPZVec<int> &sides, TPZFMatrix<STATE> &normals){
    
    //  @omar:: Just for linear mapping
    
    TPZFMatrix<REAL> mat_normals;
    TPZVec<int> v_sides;
    gel->ComputeNormals(mat_normals, v_sides);
    
    switch (gel->Type()) {
        case ECube:
        {
            int nfaces = 6;
            sides.Resize(nfaces);
            sides[0] = 20;
            sides[1] = 21;
            sides[2] = 22;
            sides[3] = 23;
            sides[4] = 24;
            sides[5] = 25;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETetraedro:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 10;
            sides[1] = 11;
            sides[2] = 12;
            sides[3] = 13;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case EQuadrilateral:
        {
            int nfaces = 4;
            sides.Resize(nfaces);
            sides[0] = 4;
            sides[1] = 5;
            sides[2] = 6;
            sides[3] = 7;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
        case ETriangle:
        {
            int nfaces = 3;
            sides.Resize(nfaces);
            sides[0] = 3;
            sides[1] = 4;
            sides[2] = 5;
            int iside = 0;
            normals.Resize(3, nfaces);
            
            for (int i = 0 ; i < v_sides.size(); i++) {
                if (nfaces <= iside) {
                    break;
                }
                if(v_sides[i] ==  sides[iside]){
                    normals(0,iside) = mat_normals(0,i);
                    normals(1,iside) = mat_normals(1,i);
                    normals(2,iside) = mat_normals(2,i);
                    iside++;
                }
            }
            
            
        }
            break;
            
        default:
        {
            std::cout << "Element not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
}

/** @brief Compute left and right geometric element indexes */
void TPZTransferFunctions::ComputeLeftRight(TPZCompMesh * transport_mesh){
    
    DebugStop();
    
}


/** @brief Compute left and right geometric element indexes */
void TPZTransferFunctions::ComputeLeftRightII(TPZCompMesh * transport_mesh){
    
    DebugStop();
    
}


/** @brief Dimensionla Measure of the elemnt */
REAL TPZTransferFunctions::DimensionalMeasure(TPZGeoEl * gel){
    
#ifdef PZDEBUG
    if (!gel) {
        DebugStop();
    }
#endif
    REAL measure = 0.0;
    int order = 10;
    int element_itself  = gel->NSides() - 1;
    TPZIntPoints * int_points = gel->CreateSideIntegrationRule(element_itself, order);
    REAL detjac, w;
    TPZVec<REAL> par(gel->Dimension(),0.0);
    TPZFMatrix<REAL> jac;
    TPZFMatrix<REAL> axes;
    TPZFMatrix<REAL> jacinv;
    for (int i = 0; i < int_points->NPoints(); i++) {
        int_points->Point(i, par, w);
        gel->Jacobian(par, jac, axes, detjac, jacinv);
        measure += w * detjac;
    }
    
    return measure;
    
}


void TPZTransferFunctions::ElementDofIndexes(TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
//    int dfvar = block.Size(dfseq);
//    long pos = block.Position(dfseq);
//    for(int jn=0; jn<dfvar; jn++) {
//        for (long is=0; is<numbersol; is++) {
//            sol[is][iv%nstate] += (STATE)phi(iv/nstate,0)*MeshSol(pos+jn,is);
//            for(d=0; d<dphix.Rows(); d++){
//                dsol[is](d,iv%nstate) += (STATE)dphix(d,iv/nstate)*MeshSol(pos+jn,is);
//            }
//        }
//        iv++;
//    }
    
    TPZStack<long> index(0,0);
    int nconnect = intel->NConnects();
    
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = intel->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = intel->Mesh()->Block().Position(seqnumber);
        long b_size = intel->Mesh()->Block().Size(seqnumber);
            for(int jb=0; jb<b_size; jb++) {
                index.Push(position + jb);
            }
    }
    
    dof_indexes = index;
    return;
}

void TPZTransferFunctions::ElementDofIndexes(TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el) {
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intel_vol = dynamic_cast<TPZInterpolationSpace * >(m_el->Element(0));
    
#ifdef PZDEBUG
    if (!intel_vol) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    int nconnect = intel_vol->NConnects();
    for (int icon = 0; icon < nconnect; icon++) {
        TPZConnect  & con = m_el->Connect(icon);
        long seqnumber = con.SequenceNumber();
        long position = m_el->Mesh()->Block().Position(seqnumber);
        int nshape = con.NShape();
        for (int ish=0; ish < nshape; ish++) {
            index.Push(position+ ish);
        }
    }
    
    dof_indexes = index;
    return;
    
}

void TPZTransferFunctions::ElementDofFaceIndexes(int connect_index,TPZInterpolationSpace * &intel, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!intel) {
        DebugStop();
    }
#endif
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = intel->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = intel->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}

void TPZTransferFunctions::ElementDofFaceIndexes(int connect_index, TPZMultiphysicsElement * &m_el, TPZVec<long> &dof_indexes){
    
    
#ifdef PZDEBUG
    if (!m_el && connect_index > 4) {
        DebugStop();
    }
#endif
    
    
    TPZStack<long> index(0,0);
    TPZConnect  & con = m_el->Connect(connect_index);
    long seqnumber = con.SequenceNumber();
    long position = m_el->Mesh()->Block().Position(seqnumber);
    int nshape = con.NShape();
    for (int ish=0; ish < nshape; ish++) {
        index.Push(position+ ish);
    }
    
    dof_indexes = index;
    return;
}


/** @brief Compute compuational mesh pair (geomechanic, gp mesh) indexed by geometric element index */
void TPZTransferFunctions::FillGeomechanicElPairs(TPZCompMesh * cmesh_mphysics){
    

    fgeomechanic_galerkinp_cindexes.Resize(0);
    
#ifdef PZDEBUG
    if (!cmesh_mphysics || !fCmeshRB_projections) {
        DebugStop();
    }
    
#endif

    fCmeshRB_projections->LoadReferences();
    TPZGeoMesh * geometry = fCmeshRB_projections->Reference();
    int dimension = geometry->Dimension();
    
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif

    std::pair<long, std::pair <long,long> > gel_indexes;
    
    for (long i = 0; i < geometry->NElements(); i++) {
        TPZGeoEl * gel = geometry->Element(i);
        
#ifdef PZDEBUG
        if (!gel) {
            DebugStop();
        }
#endif
        if (gel->HasSubElement()) {
            continue;
        }
        gel_indexes.first = gel->Index();
        gel_indexes.second.first = -1;
        gel_indexes.second.second = -1;
        fgeomechanic_galerkinp_cindexes.Push(gel_indexes);
        
    }
    
    // counting volumetric elements
    long n_elements = fgeomechanic_galerkinp_cindexes.size();
    fgeomechanic_galerkinp_cindexes.Resize(n_elements);
    
    // inserting geomechanic indexes
    cmesh_mphysics->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * cel = geometry->Element(fgeomechanic_galerkinp_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!cel) {
            //continue;
            DebugStop();
        }
#endif
        fgeomechanic_galerkinp_cindexes[iel].second.first = cel->Index();
        
    }
    
    // inserting gp indexes
    fCmeshRB_projections->LoadReferences();
    for(long iel = 0; iel < n_elements; iel++){
        
        TPZCompEl * gp_cel = geometry->Element(fgeomechanic_galerkinp_cindexes[iel].first)->Reference();
        
#ifdef PZDEBUG
        if (!gp_cel) {
            DebugStop();
        }
#endif
        fgeomechanic_galerkinp_cindexes[iel].second.second = gp_cel->Index();
        
    }
    
//    for (int k = 0; k < fgeomechanic_galerkinp_cindexes.size(); k++) {
//        std::cout << " volume k : " << k <<std::endl;
//        std::cout << " volume gel : " << fgeomechanic_galerkinp_cindexes[k].first <<std::endl;
//        std::cout << " volume cgeo : " << fgeomechanic_galerkinp_cindexes[k].second.first <<std::endl;
//        std::cout << " volume cgp : " << fgeomechanic_galerkinp_cindexes[k].second.second <<std::endl;
//    }
    
}

