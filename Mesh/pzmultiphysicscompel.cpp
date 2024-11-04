/**
 * @file
 * @brief Contains the implementation of the Multiphysics computational element methods.
 */

#include "pzmultiphysicscompel.h"

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pztrnsform.h"
#include "TPZMaterial.h"
#include "TPZMatLoadCases.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "tpzautopointer.h"
#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"
#include "TPZMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZElementMatrixT.h"
#include "pzconnect.h"
#include "pzinterpolationspace.h"
#include "pzlog.h"
#include "pzcompelwithmem.h"

#include "TPZBndCond.h"

#include <set>

using namespace pzgeom;


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzmultiphysiccompEl");
static TPZLogger logerror("pz.material.darcy");
#endif

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl() : TPZRegisterClassId(&TPZMultiphysicsCompEl::ClassId),
TPZMultiphysicsElement(), fElementVec(0){
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh, const TPZMultiphysicsCompEl<TGeometry> &copy) : TPZRegisterClassId(&TPZMultiphysicsCompEl::ClassId),
TPZMultiphysicsElement(mesh,copy),
fConnectIndexes(copy.fConnectIndexes),  fIntRule(copy.fIntRule)
{
#ifdef PZDEBUG
    TPZMaterial * material = Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        return;
    }
#endif
    
}

template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh,
                                                        const TPZMultiphysicsCompEl<TGeometry> &copy,
                                                        std::map<int64_t,int64_t> & gl2lcConMap,
                                                        std::map<int64_t,int64_t> & gl2lcElMap) : TPZRegisterClassId(&TPZMultiphysicsCompEl::ClassId){
    
    DebugStop(); // if this is called, withmem should not work
}


template <class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::TPZMultiphysicsCompEl(TPZCompMesh &mesh, TPZGeoEl *ref) :TPZRegisterClassId(&TPZMultiphysicsCompEl::ClassId),
TPZMultiphysicsElement(mesh, ref), fElementVec(0) {
}

template<class TGeometry>
TPZMultiphysicsCompEl<TGeometry>::~TPZMultiphysicsCompEl(){
    int nc = NConnects();
    for (int ic=0; ic<nc; ic++) {
        Connect(ic).RemoveDepend();
    }
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::AffineTransform(TPZVec<TPZTransform<> > &trVec) const
{
    int64_t nel;
    int side, dim, dimmf;
    nel=fElementVec.size();
    trVec.Resize(nel);
    TPZGeoEl *gelmf = Reference();
    dimmf = gelmf->Dimension();
    side = gelmf->NSides()-1;
    TPZGeoEl  *geoel;
    for (int64_t i = 0; i<nel; i++) {
        if (!fElementVec[i]) {
            continue;
        }
        geoel = fElementVec[i].Element()->Reference();
        dim =  geoel->Dimension();
        if (dim == dimmf) {
            TPZTransform<> tr(dim);
            TPZGeoElSide gelside(geoel,geoel->NSides()-1);
            TPZGeoElSide gelmfside(gelmf,side);
            if (gelside.NeighbourExists(gelmfside))
            {
                trVec[i] = tr;
                gelmfside.SideTransform3(gelside, trVec[i]);
            }
            else
            {
                trVec[i] = gelmf->BuildTransform2(side, geoel, tr);
            }
        }
        else
        {
            TPZTransform<> LocalTransf(dimmf);
            TPZGeoElSide thisgeoside(gelmf,gelmf->NSides()-1);
            TPZGeoElSide neighgeoside = fElementVec[i].Reference();
            thisgeoside.SideTransform3(neighgeoside, LocalTransf);
            TPZGeoElSide highdim(neighgeoside.Element(), neighgeoside.Element()->NSides()-1);
            trVec[i] = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
            
        }
    }
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::GetReferenceIndexVec(TPZManVector<TPZCompMesh *> cmeshVec, std::set<int64_t> &refIndexVec){
    
    if(cmeshVec.NElements() == 0) return;
    TPZCompMesh *cmesh = cmeshVec[0];
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    int64_t isub;
    int64_t ncm = cmeshVec.NElements();
    for (isub=0; isub<ncm; isub++) {
        cmeshVec[isub]->LoadReferences();
    }
    int64_t ncel;
    TPZStack<TPZCompElSide> sidevec;
    for(int64_t i = 0; i< ncm; i++){
        ncel = cmeshVec[i]->NElements();
        for (int64_t j=0; j<ncel; j++) {
            TPZCompEl * cel = cmeshVec[i]->ElementVec()[j];
            if(cel){
                TPZGeoEl *geoel = cel->Reference();
                
#ifdef PZDEBUG
                if (!geoel){
                    PZError << "Error at " << __PRETTY_FUNCTION__ << " Geometry element null!\n";
                    DebugStop();
                }
#endif
                
                int ns = geoel->NSides();
                TPZGeoElSide *geoside = new TPZGeoElSide(geoel,ns-1);
                sidevec.Resize(0);
                geoside->HigherLevelCompElementList2(sidevec, 1,1);
                int64_t nel = sidevec.NElements();
                if (nel==0){
                    //std::cout << "Incluindo elemento " << geoel->Index() << std::endl;
                    refIndexVec.insert(geoel->Index());
                }
            }
        }
    }
    
    
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout << "Number of elements : " << refIndexVec.size() << std::endl;
        sout <<"Reference index of elements : "<< std::endl;
        std::set<int64_t>::iterator it;
        for (it=refIndexVec.begin() ; it != refIndexVec.end(); it++ )
            sout << " " << *it;
        sout << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
}

/*{
 //ComputeNodElCon();
 out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
 out << "TITLE-> " << fName << "\n\n";
 
 out << "number of connects            = " << NConnects() << std::endl;
 out << "number of elements            = " << NElements() << std::endl;
 out << "number of materials           = " << NMaterials() << std::endl;
 //  out << "number of nodal bound cond    = " << NBCConnects() << endl;
 
 out << "\n\t Connect Information:\n\n";
 int i, nelem = NConnects();
 for(i=0; i<nelem; i++) {
 if(fConnectVec[i].SequenceNumber() == -1) {
 if(fConnectVec[i].HasDependency()) {
 cout << "TPZCompMesh::Print inconsistency of connect\n";
 cout << "Index " << i << ' ';
 fConnectVec[i].Print(*this,std::cout);
 }
 continue;
 }
 out << "Index " << i << ' ';
 fConnectVec[i].Print(*this,out);
 }
 out << "\n\t Computable Element Information:\n\n";
 nelem = NElements();
 for(i=0; i<nelem; i++) {
 if(!fElementVec[i]) continue;
 TPZCompEl *el = fElementVec[i];
 out << "Index " << i << ' ';
 el->Print(out);
 if(!el->Reference()) continue;
 out << "\tReference Index = " << el->Reference()->Index() << std::endl << std::endl;
 }
 out << "\n\tMaterial Information:\n\n";
 std::map<int, TPZMaterial * >::const_iterator mit;
 nelem = NMaterials();
 for(mit=fMaterialVec.begin(); mit!= fMaterialVec.end(); mit++) {
 TPZMaterial *mat = mit->second;
 if (!mat) {
 DebugStop();
 }
 mat->Print(out);
 }
 
 }*/

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Print(std::ostream & out) const {
    
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompEl::Print(out);
    out << "Integration rule ";
    GetIntegrationRule().Print(out);
    TPZGeoEl *gref = Reference();
    if(gref)
    {
        out << "\nCenter coordinate: ";
        TPZVec< REAL > centerMaster( this->Reference()->Dimension(),0. );
        TPZVec< REAL > centerEuclid( 3,0.);
        this->Reference()->CenterPoint(this->Reference()->NSides()-1,centerMaster);
        this->Reference()->X(centerMaster,centerEuclid);
        out << centerEuclid << std::endl;
    }
    int matid = 0;
    
    if(gref)
    {
        out << "Material id " << gref->MaterialId() << "\n";
    }
    else {
        out << "No material\n";
    }
    
    out << "Number of connects = " << NConnects();
    out << "Connect indexes : ";
    int nod;
    for(nod=0; nod< NConnects(); nod++)
    {
        out << ConnectIndex(nod) <<  ' ' ;
    }
    if(this->Reference())
    {
        out << "Reference Index = " << this->Reference()->Index() << std::endl;
    }
    std::list<TPZOneShapeRestraint> restraints = this->GetShapeRestraints();
    if(restraints.size())
    {
        out << "One shape restraints\n";
        for (std::list<TPZOneShapeRestraint>::const_iterator it = restraints.begin(); it != restraints.end(); it++) {
            it->Print(out);
        }
    }
    
    out << "\nComputational elements within multi-physics element: \n";
    int64_t nmesh  = fElementVec.size();
    for(int64_t iel = 0; iel< nmesh; iel++ ){
        
        out << "\nComputational element belonging to the mesh " << iel <<":\n";
        TPZCompEl *cel = fElementVec[iel].Element();
        if(!cel){
            out << "\n There is not element to computational mesh " << iel <<"\n";
            continue;
        }
        cel->Print(out);
        if(!cel->Reference()) continue;
        out << "\tReference Index = " << cel->Reference()->Index();
        
        TPZManVector<TPZTransform<> > tr;
        AffineTransform(tr);
        out << "\n\tAffine transformation of the multiphysics element for this computational element:"<<"\n";
        out << std::endl;
        out <<"\t" << tr[iel];
    }
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::ShortPrint(std::ostream & out) const {

    out << __PRETTY_FUNCTION__ << std::endl;
    out << "\n\tMulti-physics element:: \n";
    TPZCompEl::Print(out);

    std::list<TPZOneShapeRestraint> restraints = this->GetShapeRestraints();
    if(restraints.size())
    {
        out << "One shape restraints\n";
        for (std::list<TPZOneShapeRestraint>::const_iterator it = restraints.begin(); it != restraints.end(); it++) {
            it->Print(out);
        }
    }

    out << "\n\tComputational elements within this multi-physics element: \n\n";
    int64_t nmesh  = fElementVec.size();
    for(int64_t iel = 0; iel< nmesh; iel++ ){

        out << "Computational element belonging to the mesh " << iel <<":\n";
        TPZCompEl *cel = fElementVec[iel].Element();
        if(!cel){
            out << "There is not element to computational mesh " << iel <<"\n";
            continue;
        }
        cel->ShortPrint(out);
        if(!cel->Reference()) continue;
        out << "\tReference Index = " << cel->Reference()->Index();

        TPZManVector<TPZTransform<> > tr;
        AffineTransform(tr);
        out << "\n\tAffine transformation of the multiphysics element for this computational element:"<<"\n";
        out << std::endl;
        out <<"\t" << tr[iel];
    }
    out << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}


template <class TGeometry>
TPZCompEl * TPZMultiphysicsCompEl<TGeometry>::Clone(TPZCompMesh &mesh) const {
    
    PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
    DebugStop();
    return 0;
}

template <class TGeometry>
TPZCompEl* TPZMultiphysicsCompEl<TGeometry>::ClonePatchEl(TPZCompMesh &mesh,
                                                          std::map<int64_t,int64_t> & gl2lcConMap,
                                                          std::map<int64_t,int64_t> & gl2lcElMap) const {
    PZError << "Error at " << __PRETTY_FUNCTION__ << " method is not implementedl!\n";
    DebugStop();
    return 0;
}

template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::NConnects() const {
    return fConnectIndexes.NElements();
}

template <class TGeometry>
int64_t TPZMultiphysicsCompEl<TGeometry>::ConnectIndex(int i) const {
    return fConnectIndexes[i];
}

template <class TGeometry>
int64_t TPZMultiphysicsCompEl<TGeometry>::ConnectIndex(int elem, int connect) const {
    
    int first = 0;
    for(int64_t el=0; el<elem; el++){
        TPZCompEl *cel = fElementVec[el].Element();
        if (cel && fActiveApproxSpace[el]) {
            first+=cel->NConnects();
        }
    }
    return fConnectIndexes[first+connect];
}


template <class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::Dimension() const {
    for(int el = 0; el < fElementVec.NElements(); el++)
    {
        if(fElementVec[el])
        {
            return fElementVec[el].Element()->Dimension();
        }
    }
    return -1;
}



template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::Integrate(int variable, TPZVec<STATE> & value){
    TPZMaterial * material = this->Material();
    if(!material){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " : no material for this element\n";
        return;
    }
    if (!this->Reference()){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " : no reference element\n";
        return;
    }
    
    int64_t nref = fElementVec.size();
    TPZVec<TPZMaterialDataT<STATE>> datavec;
    datavec.resize(nref);
    
#ifdef PZDEBUG
    if (nref != datavec.size()) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
        DebugStop();
    }
#endif
    
    TPZVec<int> nshape(nref);
    for (int64_t iref = 0; iref<nref; iref++)
    {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        msp->InitMaterialData(datavec[iref]);
    }
    
    
    //	TPZManVector<REAL, 3> intpoint(dim,0.);
    //	const int varsize = material->NSolutionVariables(variable);
    //	value.Resize(varsize);
    value.Fill(0.);
    //
    //	const TPZIntPoints &intrule = this->GetIntegrationRule();
    //	int npoints = intrule.NPoints(), ip, iv;
    //	TPZManVector<REAL> sol(varsize);
    //	for(ip=0;ip<npoints;ip++){
    //		intrule.Point(ip,intpoint,weight);
    //		sol.Fill(0.);
    //		this->Solution(intpoint, variable, sol);
    //		//Tiago: Next call is performed only for computing detcaj. The previous method (Solution) has already computed jacobian.
    //		//       It means that the next call would not be necessary if I wrote the whole code here.
    //		this->Reference()->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
    //		weight *= fabs(data.detjac);
    //		for(iv = 0; iv < varsize; iv++) {
    //#if !BUILD_COMPLEX_PROJECTS
    //			DebugStop();
    //#else
    //			value[iv] += sol[iv]*weight;
    //#endif
    //		}//for iv
    //	}//for ip
}//method


template<class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::SolutionT(TPZVec<REAL> &qsi, int var,TPZVec<TVar> &sol)
{
    
    if (var >= 99) {
        TPZCompEl::Solution(qsi, var, sol);
        return;
    }
    
    auto *material =
        dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(this->Material());
    if(!material){
        sol.Resize(0);
        return;
    }
    
    //    TPZManVector<REAL,3> xi(qsi.size());
    //    Reference()->CenterPoint(Reference()->NSides()-1,xi);
    //    for (int i=0; i<xi.size(); i++) {
    //        qsi[i] += 0.001*(xi[i]-qsi[i]);
    //    }
    
    
    TPZManVector<TPZTransform<> > trvec;
    AffineTransform(trvec);
    
    TPZManVector<REAL,3> myqsi(qsi);
    myqsi.resize(qsi.size());
    
    int64_t nref = fElementVec.size();
    TPZManVector<TPZMaterialDataT<TVar>,4> datavec;
    datavec.resize(nref);
    this->InitMaterialData(datavec);
    material->FillDataRequirements(datavec);

    TPZGeoEl * ref = this->Reference();

    
    if (!ref){
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
        return;
    }
    ref->Jacobian(qsi, datavec[0].jacobian, datavec[0].axes,
                  datavec[0].detjac , datavec[0].jacinv);
    for (int64_t iref = 0; iref<nref; iref++)
    {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        trvec[iref].Apply(qsi, myqsi);
        datavec[iref].p = msp->MaxOrder();
        
        TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
        if(shapetype==datavec[iref].EVecandShape){
            msp->ComputeRequiredData(datavec[iref], myqsi);
        }
        else
        {
            datavec[iref].jacobian = datavec[0].jacobian;
            datavec[iref].axes = datavec[0].axes;
            datavec[iref].detjac = datavec[0].detjac;
            datavec[iref].jacinv = datavec[0].jacinv;
            msp->ComputeShape(myqsi,datavec[iref]);
        }
        constexpr bool hasPhi{true};
        msp->ComputeSolution(myqsi, datavec[iref],hasPhi);
        
        datavec[iref].x.Resize(3);
        msp->Reference()->X(myqsi, datavec[iref].x);
    }
    
    material->Solution(datavec, var, sol);
    for (int64_t iref = 0; iref<nref; iref++)
    {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        msp->CleanupMaterialData(datavec[iref]);
    }
}



template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::SetConnectIndex(int inode, int64_t index){
    fConnectIndexes[inode] = index;
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef)
{
    ek.fMesh = Mesh();
    ef.fMesh = Mesh();
    ek.fType = TPZElementMatrix::EK;
    ef.fType = TPZElementMatrix::EF;
    const int ncon = this->NConnects();
    int numeq = 0;
    int ic;
    
    for(ic=0; ic<ncon; ic++)
    {
        TPZConnect &c = Connect(ic);
        int64_t neqThisConn = c.NDof(*Mesh());
        numeq += neqThisConn;
    }
    
    int64_t nref = this->fElementVec.size();
    int nstate = 0;
    //nstate=1;
    
    const int numloadcases = [this](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(this->Material()); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
    for (int64_t iref=0; iref<nref; iref++) {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if (! msp) {
            continue;
        }
        TPZMaterial *mat = msp->Material();
        nstate += mat->NStateVariables();
    }
    
    const int numstate = nstate;
    ek.Matrix().Redim(numeq,numeq);
    ef.Matrix().Redim(numeq,numloadcases);
    ek.Block().SetNBlocks(ncon);
    ef.Block().SetNBlocks(ncon);
    
    int i;
    for(i=0; i<ncon; i++){
        unsigned int ndof = Connect(i).NDof(*Mesh());
#ifdef PZDEBUG
        TPZConnect &c = Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
        ek.Block().Set(i,ndof);
        ef.Block().Set(i,ndof);
    }
    ek.fConnect.Resize(ncon);
    ef.fConnect.Resize(ncon);
    for(i=0; i<ncon; i++){
        (ek.fConnect)[i] = ConnectIndex(i);
        (ef.fConnect)[i] = ConnectIndex(i);
    }
}//void

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::InitializeElementMatrix(TPZElementMatrix &ef)
{
    const int ncon = this->NConnects();
    ef.fMesh = Mesh();
    ef.fType = TPZElementMatrix::EF;
    int numeq = 0;
    int ic;
    
    for(ic=0; ic<ncon; ic++)
    {
        int64_t neqThisConn = Connect(ic).NDof(*Mesh());
        numeq += neqThisConn;
    }
    
    int64_t nref = this->fElementVec.size();
    int nstate = 0;
    //nstate=1;
    const int numloadcases = [this](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(this->Material()); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
    for (int64_t iref=0; iref<nref; iref++) {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if (! msp) {
            continue;
        }
        TPZMaterial *mat = msp->Material();
        nstate += mat->NStateVariables();
    }
    
    const int numstate = nstate;
    ef.Matrix().Redim(numeq,numloadcases);
    ef.Block().SetNBlocks(ncon);
    
    int i;
    for(i=0; i<ncon; i++){
        unsigned int ndof = Connect(i).NDof(*Mesh());
#ifdef PZDEBUG
        TPZConnect &c = Connect(i);
        if (c.NShape()*c.NState() != ndof) {
            DebugStop();
        }
#endif
        ef.Block().Set(i,ndof);
    }
    ef.fConnect.Resize(ncon);
    for(i=0; i<ncon; i++){
        (ef.fConnect)[i] = ConnectIndex(i);
    }
}//void

template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::InitMaterialDataT(TPZVec<TPZMaterialDataT<TVar> > &dataVec, TPZVec<int64_t> *indices)
{
    int64_t nref = this->fElementVec.size();
    
#ifdef PZDEBUG
    if (nref != dataVec.size()) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
        DebugStop();
    }
#endif
    if(indices){
        int64_t nindices = indices->size();
        TPZVec<int> nshape(nindices);
        for (int64_t iref = 0; iref <nindices; iref++) {
            int64_t indiciref = indices->operator[](iref);
            if(fElementVec[indiciref])
            {
                dataVec[indiciref].gelElId = fElementVec[indiciref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[indiciref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[indiciref]);
        }
    }else{
        TPZVec<int> nshape(nref);
        for (int64_t iref = 0; iref < nref; iref++)
        {
            if(fElementVec[iref])
            {
                dataVec[iref].gelElId = fElementVec[iref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[iref]);
        }
    }
    
    int n_active_approx_spaces = fActiveApproxSpace.size();
    if (n_active_approx_spaces == 0) { /// it preserves the integrity for old version of multiphycis codes.
        fActiveApproxSpace.Resize(nref, 1);
    }
    
    for (int64_t iref = 0; iref < nref; iref++) {
        dataVec[iref].fActiveApproxSpace = fActiveApproxSpace[iref];
		if(!fElementVec[iref].Element()){
			dataVec[iref].fActiveApproxSpace = 0;
		}
    }
    auto * mat =
        dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(this->Material());
    mat->FillDataRequirements(dataVec);
    
}

template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::InitMaterialDataT(std::map<int, TPZMaterialDataT<TVar> > &dataVec, TPZVec<int64_t> *indices)
{
    int64_t nref = this->fElementVec.size();
    
//#ifdef PZDEBUG
//    if (nref != dataVec.size()) {
//        PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
//        DebugStop();
//    }
//#endif
    if(indices){
        int64_t nindices = indices->size();
        TPZVec<int> nshape(nindices);
        for (int64_t iref = 0; iref <nindices; iref++) {
            int64_t indiciref = indices->operator[](iref);
            if(fElementVec[indiciref])
            {
                dataVec[indiciref].gelElId = fElementVec[indiciref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[indiciref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[indiciref]);
        }
    }else{
        TPZVec<int> nshape(nref);
        for (int64_t iref = 0; iref < nref; iref++)
        {
            if(fElementVec[iref])
            {
                dataVec[iref].gelElId = fElementVec[iref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[iref]);
        }
    }
    
    int n_active_approx_spaces = fActiveApproxSpace.size();
    if (n_active_approx_spaces == 0) { /// it preserves the integrity for old version of multiphycis codes.
        fActiveApproxSpace.Resize(nref, 1);
    }
    
    for (auto &it : dataVec) {
        it.second.fActiveApproxSpace = fActiveApproxSpace[it.first];
    }
    auto * mat =
        dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(this->Material());
    mat->FillDataRequirements(dataVec);
    
}

/**
 * @brief Initialize a material data vector and its attributes based on element dimension, number
 * of state variables and material definitions
 */
template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::CleanupMaterialDataT(TPZVec<TPZMaterialDataT<TVar>> &dataVec)
{

    int64_t nref = this->fElementVec.size();

    TPZVec<int> nshape(nref);
    for (int64_t iref = 0; iref < nref; iref++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if (!msp) {
            continue;
        }
        // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
        msp->CleanupMaterialData(dataVec[iref]);
    }
}

template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::CalcStiffT(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef)
{
    TPZMaterial * material = Material();
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
    if(!material || !matCombined){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ek.Reset();
        ef.Reset();
        return;
    }
    
    auto *nullmat = dynamic_cast<TPZNullMaterialCS<TVar> *>(material);
    if(nullmat)
    {
        ek.Reset();
        ef.Reset();
        ek.fType = TPZElementMatrix::EK;
        ef.fType = TPZElementMatrix::EF;
        return;
    }
    InitializeElementMatrix(ek,ef);
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    
    TPZManVector<TPZMaterialDataT<TVar>,6> datavec;
    const int64_t nref = fElementVec.size();
    datavec.resize(nref);
    InitMaterialData(datavec);
    
    TPZManVector<TPZTransform<> > trvec;
    AffineTransform(trvec);
    
    int dim = Dimension();
    TPZAutoPointer<TPZIntPoints> intrule;
    TPZGeoEl *ref = this->Reference();

    if(fIntRule.NPoints() == 1) {
        
        TPZManVector<int,4> ordervec;
        //ordervec.resize(nref);
        for (int64_t iref=0;  iref<nref; iref++)
        {
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
            int svec;
            if(msp)
            {
                ordervec.Resize(ordervec.size()+1);
                svec = ordervec.size();
            }
            else
            {
                continue;
            }
            datavec[iref].p = msp->MaxOrder();
            ordervec[svec-1] = datavec[iref].p;
        }
        
        int order = matCombined->IntegrationRuleOrder(ordervec);
        
        intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
        
        TPZManVector<int,4> intorder(dim,order);
        intrule->SetOrder(intorder);
        if(intrule->NPoints() > 1000) {
            DebugStop();
        }
    } else {
        intrule = fIntRule.Clone();
    }
    
    int intrulepoints = intrule->NPoints();
    TPZManVector<REAL,4> intpointtemp(TGeometry::Dimension,0.);
    REAL weight = 0.;

    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule->Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        for (int i = 0; i < fElementVec.size(); i++) {
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[i].Element());
            if (!msp) {
                continue;
            }
            datavec[i].intLocPtIndex = int_ind;
        }
        
        this->ComputeRequiredData(intpointtemp,trvec,datavec);
        
        matCombined->Contribute(datavec,weight,ek.fMat,ef.fMat);
    }//loop over integration points
    
    CleanupMaterialData(datavec);
}//CalcStiff

template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::CalcResidualT(TPZElementMatrixT<TVar> &ef)
{
    TPZMaterial * material = Material();
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
    if(!material || !matCombined){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ef.Reset();
        return;
    }
    
    InitializeElementMatrix(ef);
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    
    TPZManVector<TPZMaterialDataT<TVar>,3> datavec;
    const int64_t nref = fElementVec.size();
    datavec.resize(nref);
    InitMaterialData(datavec);
    
    TPZManVector<TPZTransform<> > trvec;
    AffineTransform(trvec);
    
    int dim = Dimension();
    TPZAutoPointer<TPZIntPoints> intrule;
    
    TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
    REAL weight = 0.;
    
    TPZManVector<int,3> ordervec;
    //ordervec.resize(nref);
    for (int64_t iref=0;  iref<nref; iref++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        int svec;
        if(msp)
        {
            ordervec.Resize(ordervec.size()+1);
            svec = ordervec.size();
        }
        else
        {
            continue;
        }
        datavec[iref].p = msp->MaxOrder();
        ordervec[svec-1] = datavec[iref].p;
    }
    int order = matCombined->IntegrationRuleOrder(ordervec);
    
    TPZGeoEl *ref = this->Reference();
    intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
    
    TPZManVector<int,3> intorder(dim,order);
    intrule->SetOrder(intorder);
    int intrulepoints = intrule->NPoints();
    if(intrulepoints > 1000) {
        DebugStop();
    }
    
    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    int nmeshes = datavec.size();
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule->Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        for (int imesh = 0; imesh < nmeshes; imesh++) {
            datavec[imesh].intLocPtIndex = int_ind;
        }
        
        
        this->ComputeRequiredData(intpointtemp,trvec,datavec);
        
        matCombined->Contribute(datavec,weight,ef.fMat);
    }//loop over integratin points
    
    CleanupMaterialData(datavec);
    
}//CalcStiff

/**
 * @brief Compute the integral of a variable
 */
template <class TGeometry>
TPZVec<STATE> TPZMultiphysicsCompEl<TGeometry>::IntegrateSolution(int var) const
{
    TPZManVector<STATE> result;
    TPZMaterial * material = Material();
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(material);
    auto *matBC = dynamic_cast<TPZBndCondT<STATE>*>(material);
    if(!material || !matCombined){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        return result;
    }
    if (matBC) {
        material = matBC->Material();
        matCombined =dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(matBC->Material());
    }
    
    if (this->NConnects() == 0) return result;//boundary discontinuous elements have this characteristic
    
    TPZMultiphysicsCompEl<TGeometry> *thisnonconst = (TPZMultiphysicsCompEl<TGeometry> *) this;
    
    TPZManVector<TPZMaterialDataT<STATE>,4> datavec;
    const int64_t nref = fElementVec.size();
    datavec.resize(nref);
    thisnonconst->InitMaterialData(datavec);
    
    TPZManVector<TPZTransform<> > trvec;
    AffineTransform(trvec);
    
    int dim = Dimension();
    
    const TPZIntPoints &intrule = GetIntegrationRule();
    TPZManVector<REAL,3> intpoint(dim,0.), intpointtemp(dim,0.);
    REAL weight = 0.;
    const TPZCompElWithMem<TPZMultiphysicsCompEl<TGeometry>> *CompElWithMem = dynamic_cast<const TPZCompElWithMem<TPZMultiphysicsCompEl<TGeometry>> *>(this);
    TPZManVector<int64_t> intpoint_indices;
    if(CompElWithMem)
    {
        CompElWithMem->GetMemoryIndices(intpoint_indices);
#ifdef PZDEBUG
        if(intrule.NPoints() != intpoint_indices.size()) DebugStop();
#endif
    }
    //    int intrulepoints = intrule->NPoints();

//    TPZManVector<int,3> ordervec;
//    //ordervec.resize(nref);
   for (int64_t iref=0;  iref<nref; iref++)
   {
       datavec[iref].fNeedsSol = true;
//        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
//        int svec;
//        if(msp)
//        {
//            ordervec.Resize(ordervec.size()+1);
//            svec = ordervec.size();
//        }
//        else
//        {
//            continue;
//        }
//        datavec[iref].p = msp->MaxOrder();
//        ordervec[svec-1] = datavec[iref].p;
   }
//    int order = matCombined->IntegrationRuleOrder(ordervec);
//
//    TPZGeoEl *ref = this->Reference();
//    intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
//
//    TPZManVector<int,3> intorder(dim,order);
//    intrule->SetOrder(intorder);
//    int intrulepoints = intrule->NPoints();
//    if(intrulepoints > 1000) {
//        DebugStop();
//    }
    TPZGeoEl *ref = this->Reference();
    int intrulepoints = intrule.NPoints();
    int nvar = material->NSolutionVariables(var);
    TPZManVector<STATE> solout(var);
    result.Resize(nvar, 0.);
    
    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule.Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        datavec[0].intLocPtIndex = int_ind;
        datavec[1].intLocPtIndex = int_ind;
        
        
        thisnonconst->ComputeRequiredData(intpointtemp,trvec,datavec);
        if(CompElWithMem) {
            datavec[0].intGlobPtIndex = intpoint_indices[int_ind];
        }
        matCombined->Solution(datavec, var, solout);
        
        for (int iv=0; iv<nvar; iv++) {
            result[iv] += weight*solout[iv];
        }
        
    }//loop over integratin points
    
    return result;
}

template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::SetIntegrationRule(int int_order) {
    
#ifdef PZDEBUG
    if (!this->Reference()) {
        DebugStop();
    }
#endif
    int dim = this->Reference()->Dimension();
    TPZManVector<int,3> ord(dim,int_order);
    if (TPZCompEl::fIntegrationRule) {
        TPZCompEl::fIntegrationRule->SetOrder(ord);
    }else{
        fIntRule.SetOrder(ord);
    }
    
}

/// After adding the elements initialize the integration rule
template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::InitializeIntegrationRule()
{
    if (fElementVec.NElements() == 0) {
        DebugStop(); // This case should be treated before
    }
    
    auto basemat = this->Material();
    auto *material =
        dynamic_cast<TPZMatCombinedSpaces*>(basemat);
    if(!material){
        std::cout << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        if(basemat) std::cout << "basemat \n";
        basemat->Print();
        DebugStop();
    }
    
    const int64_t nref = fElementVec.size();
    TPZGeoEl *gel = this->Reference();
    int dim = gel->Dimension();
    TPZStack<int> ordervec;
    for (int64_t iref=0;  iref<nref; iref++)
    {
        if(fActiveApproxSpace[iref] == 0){
            continue;
        }
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(msp)
        {
            ordervec.Push(msp->MaxOrder());
        }
    }
    int order = material->IntegrationRuleOrder(ordervec);
    TPZManVector<int,3> orderdim(dim,order);
    fIntRule.SetOrder(orderdim);
    
}



template <class TGeometry>
TPZIntPoints & TPZMultiphysicsCompEl<TGeometry>::GetIntegrationRule()
{
    if (TPZCompEl::fIntegrationRule) {
        return *TPZCompEl::fIntegrationRule;
    }
    return fIntRule;
}

template <class TGeometry>
const TPZIntPoints & TPZMultiphysicsCompEl<TGeometry>::GetIntegrationRule() const
{
    if (TPZCompEl::fIntegrationRule) {
        return *TPZCompEl::fIntegrationRule;
    }
    return fIntRule;
}


template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::EvaluateError(TPZVec<REAL> &errors, bool store_errors) {
    if (this->Mesh()->GetSolType() ==  EReal){
        return EvaluateErrorT<STATE>(errors,store_errors);
    }else{
        return EvaluateErrorT<CSTATE>(errors,store_errors);
    }
}

template<class TGeometry>
template<class TVar>
void TPZMultiphysicsCompEl<TGeometry>::EvaluateErrorT(TPZVec<REAL> &errors, bool store_errors) {
  errors.Fill(0.);
  auto *nullmat = dynamic_cast<TPZNullMaterialCS<TVar> *>(this->Material());
    if(nullmat) return;
  auto *mat =
      dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(this->Material());
  auto *matError =
      dynamic_cast<TPZMatErrorCombinedSpaces<TVar>*>(mat);
  
  //TPZMaterial * matptr = material.operator->();
  if (!mat) {
      PZError << __PRETTY_FUNCTION__;
      PZError<<" no material for this element\n";
      return;
  }
  if(!matError){
//      PZError << __PRETTY_FUNCTION__;
//      PZError<<" no error interface for this element.\n";
//      PZError<<"See TPZMatErrorCombinedSpaces\n";
      return;
  }
#ifdef PZ_LOG
    if(logerror.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Element index " << Index();
        LOGPZ_DEBUG(logerror,sout.str())
    }
#endif
    

  const int NErrors = matError->NEvalErrors();
  errors.Resize(NErrors);
  errors.Fill(0.);

  // Return if BC
  if (this->NConnects() == 0) return;

  // Adjust the order of the integration rule
  int dim = Dimension();
  TPZAutoPointer<TPZIntPoints> intrule = this->GetIntegrationRule().Clone();
  int maxIntOrder = intrule->GetMaxOrder();
  TPZManVector<int, 3> prevorder(dim), maxorder(dim, maxIntOrder);
  intrule->GetOrder(prevorder);
  const int order_limit = 15;
  if (maxIntOrder > order_limit) {
    if (prevorder[0] > order_limit) {
      maxIntOrder = prevorder[0];
    } else {
      maxIntOrder = order_limit;
    }
  }

  intrule->SetOrder(maxorder);

  const int ndof = this->Material()->NStateVariables();
  TPZManVector<STATE, 10> u_exact(ndof);
  TPZFNMatrix<3, STATE> du_exact(dim, ndof);
  TPZManVector<REAL, 10> intpoint(Reference()->Dimension()), values(NErrors);
  values.Fill(0.0);
  REAL weight;

  TPZManVector<TPZMaterialDataT<TVar>, 6> datavec;
  const int64_t nref = fElementVec.size();
  datavec.resize(nref);
  InitMaterialData(datavec);
  int iactive = -1;
  for (unsigned int i = 0; i < nref; ++i) {
    if (datavec[i].fShapeType != TPZMaterialData::EEmpty) {
      datavec[i].fNeedsSol = true;
      iactive = i;
    }
  }

  if (iactive == -1) DebugStop();

  TPZManVector<TPZTransform<> > trvec;
  AffineTransform(trvec);

  int nintpoints = intrule->NPoints();
  TPZFNMatrix<9, REAL> jac, axe, jacInv;
  REAL detJac;
  TPZGeoEl *ref = this->Reference();

  for (int nint = 0; nint < nintpoints; nint++) {

    intrule->Point(nint, intpoint, weight);

    ref->Jacobian(intpoint, jac, axe, detJac, jacInv);
    this->ComputeRequiredData(intpoint, trvec, datavec);

    weight *= fabs(detJac);
    matError->Errors(datavec, values);

    for (int ier = 0; ier < NErrors; ier++) {
      errors[ier] += values[ier] * weight;
    }

  }//fim for : integration rule
  //Norma sobre o elemento
  for (int ier = 0; ier < NErrors; ier++) {
      errors[ier] = sqrt(errors[ier]);
  }//for ier

  if (store_errors) {
    int64_t index = Index();
    TPZFMatrix<STATE> &elvals = Mesh()->ElementSolution();
    if (elvals.Cols() < NErrors) {
      std::cout << "The element solution of the mesh should be resized before EvaluateError\n";
      DebugStop();
    }
    for (int ier = 0; ier < NErrors; ier++) {
      elvals(index, ier) = errors[ier];
    }
  }
  intrule->SetOrder(prevorder);
}


/** Returns the maximum interpolation order of all connected elements */
template <class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::PolynomialOrder(TPZVec<int> &order) const {
    const int64_t nref = fElementVec.size();
    order.resize(nref);
    for (int64_t iref = 0; iref < nref; iref++) {
        
        TPZInterpolationSpace *msp = dynamic_cast<TPZInterpolationSpace *> (fElementVec[iref].Element());
        order[iref] = msp ? msp->MaxOrder() : 0;
    }
    
    //    TPZMaterial * material = Material();
    //    if(!material)
    //    {
    //        int matid = Reference()->MaterialId();
    //        std::cout << "Matid " << matid << " for multiphysics element index " << Index() << " does not exist\n";
    //        DebugStop();
    //    }
    //    int order = material->IntegrationRuleOrder(ordervec);
    //    return order;
}

template<class TGeometry>
int TPZMultiphysicsCompEl<TGeometry>::ClassId() const{
    return Hash("TPZMultiphysicsCompEl") ^ TPZMultiphysicsElement::ClassId() << 1 ^ TGeometry().ClassId() << 2;
}

#include "pzgraphel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "pzgraphel1d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pztrigraph.h"
#include "tpzgraphelt2dmapped.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel.h"
#include "TPZBndCond.h"
#include "pzgraphmesh.h"

template<class TGeometry>
void TPZMultiphysicsCompEl<TGeometry>::CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension)
{
    
    
    TPZGeoEl *ref = Reference();
    if (ref->Dimension() != dimension) {
        return;
    }
    TPZMaterial * material = Material();
    
    TPZBndCond * BDC = dynamic_cast<TPZBndCond * > (material);
    
    if (BDC) {
        return;
    }
    int matid = material->Id();
    int nsides = ref->NSides();
    bool to_postpro = grmesh.Material_Is_PostProcessed(matid);
    if(dimension == 2 && to_postpro){
        if(nsides == 9){
            new TPZGraphElQ2dd(this,&grmesh);
            return;
        }
        if(nsides == 7){
            new TPZGraphElT2dMapped(this,&grmesh);
            return;
        }
    }//2d
    
    if(dimension == 3 && to_postpro){
        if(nsides == 27){
            new TPZGraphElQ3dd(this,&grmesh);
            return;
        }//cube
        if(nsides == 21){
            new TPZGraphElPrismMapped(this,&grmesh);
            return;
        }//prism
        if(nsides == 15){
            new TPZGraphElT3d(this,&grmesh);
            return;
        }//tetra
        if(nsides == 19){
            new TPZGraphElPyramidMapped(this,&grmesh);
            return;
        }//pyram
    }//3d
    
    if(dimension == 1  && to_postpro){
        new TPZGraphEl1dd(this,&grmesh);
    }//1d
}

//---------------------------------------------------------------	
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoCube>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>;
template class TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid>;


/*
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoCube> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >;
 template class TPZCompElWithMem<TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid> >;
 */


TPZCompEl * CreateMultiphysicsPointEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint>(mesh, gel);
}

TPZCompEl * CreateMultiphysicsLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear>(mesh,gel);
}

TPZCompEl * CreateMultiphysicsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle >(mesh,gel);
}

TPZCompEl * CreateMultiphysicsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad>(mesh,gel);
}

TPZCompEl * CreateMultiphysicsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoCube >(mesh,gel);
}

TPZCompEl * CreateMultiphysicsPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism>(mesh,gel);
}

TPZCompEl * CreateMultiphysicsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra>(mesh,gel);
}

TPZCompEl * CreateMultiphysicsPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh) {
    return new TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid >(mesh,gel);
}

//--------------------- WITH MEMORY ----------------------

TPZCompEl *CreateMultiphysicsPointElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPoint> >(mesh,gel) ;
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsLinearElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoLinear> >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsQuadElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoQuad> >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsTriangleElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoTriangle > >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsCubeElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoCube > >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsPrismElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPrism> >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsPyramElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoPyramid > >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
TPZCompEl *CreateMultiphysicsTetraElWithMem(TPZGeoEl *gel,TPZCompMesh &mesh) {
    //	if(!gel->Reference() && gel->NumInterfaces() == 0)
    return new TPZCompElWithMem < TPZMultiphysicsCompEl<pzgeom::TPZGeoTetrahedra> >(mesh,gel);
    //	index = -1;
    //	return NULL;
}
