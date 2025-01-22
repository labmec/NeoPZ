//
//  TPZMultiphysicsCompMesh.cpp
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicscompel.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mulptiphysicscompmesh");
#endif


TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh() : TPZCompMesh(){
    
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(TPZGeoMesh * gmesh, bool isComplex) : TPZCompMesh(gmesh, isComplex){
    
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(TPZAutoPointer<TPZGeoMesh>  gmesh, bool isComplex) : TPZCompMesh(gmesh, isComplex),
    m_active_approx_spaces(), m_mesh_vector(){
    
}

TPZMultiphysicsCompMesh::TPZMultiphysicsCompMesh(const TPZMultiphysicsCompMesh &other) : TPZCompMesh(other) {
    m_active_approx_spaces  = other.m_active_approx_spaces;
    m_mesh_vector           = other.m_mesh_vector;
}


TPZMultiphysicsCompMesh::~TPZMultiphysicsCompMesh()
{
    m_active_approx_spaces.Resize(0);
    m_mesh_vector.Resize(0);
    int64_t nc = fConnectVec.NElements();
    for (int ic=0; ic<nc; ic++) {
        fConnectVec[ic].RemoveDepend();
    }
}

TPZMultiphysicsCompMesh & TPZMultiphysicsCompMesh::operator=(const TPZMultiphysicsCompMesh &other){
    
    if (this != & other) // prevent self-assignment
    {
        TPZCompMesh::operator=(other);
        m_active_approx_spaces  = other.m_active_approx_spaces;
        m_mesh_vector           = other.m_mesh_vector;
    }
    return *this;
}

TPZVec<TPZCompMesh *> & TPZMultiphysicsCompMesh::MeshVector() {
    return  m_mesh_vector;
}

TPZVec<MSpaceConfig> &  TPZMultiphysicsCompMesh::GetActiveApproximationSpaces(){
    return m_active_approx_spaces;
}

/// Set active approximation spaces
void TPZMultiphysicsCompMesh::BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector)
{
    TPZManVector<MSpaceConfig> active(mesh_vector.size(),ETestTrial);
    BuildMultiphysicsSpace(active,mesh_vector);
}

void TPZMultiphysicsCompMesh::BuildMultiphysicsSpace(TPZVec<MSpaceConfig> & active_approx_spaces, const TPZVec<TPZCompMesh * > & mesh_vector){

    m_active_approx_spaces = active_approx_spaces;
    m_mesh_vector          = mesh_vector;
    if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
        std::cout<< "TPZMultiphysicsCompMesh:: The vector provided should have the same size." << std::endl;
        DebugStop();
    }
    
    int n_approx_spaces = m_mesh_vector.size();
  
    SetNMeshes(n_approx_spaces);
    Reference()->ResetReference();
    if (ApproxSpace().Style() == TPZCreateApproximationSpace::EMultiphysics)
    {
        SetAllCreateFunctionsMultiphysicElem();
    }
    else if (ApproxSpace().Style() == TPZCreateApproximationSpace::EMultiphysicsSBFem)
    {
        SetAllCreateFunctionsSBFemMultiphysics();
    }
    else {
        DebugStop(); // Please set one of the above available spaces
    }
    // delete all elements and connects in the mesh
    CleanElementsConnects();
    TPZCompMesh::AutoBuild();
    AddElements();
    AddConnects();
    LoadSolutionFromMeshes();
    // reorder the connect sequence numbers in case there are hanging nodes
    CleanUpUnconnectedNodes();
}

void TPZMultiphysicsCompMesh::BuildMultiphysicsSpaceWithMemory(TPZVec<MSpaceConfig> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector){
    m_active_approx_spaces = active_approx_spaces;
    m_mesh_vector          = mesh_vector;
    if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
        std::cout<< "TPZMultiphysicsCompMesh:: The vector provided should have the same size." << std::endl;
        DebugStop();
    }
    
    int n_approx_spaces = m_mesh_vector.size();
    
    SetNMeshes(n_approx_spaces);
    Reference()->ResetReference();
    SetAllCreateFunctionsMultiphysicElemWithMem();
    ApproxSpace().CreateWithMemory(true);
    // delete all elements and connects in the mesh
    CleanElementsConnects();
    TPZCompMesh::AutoBuild();
    AddElements();
    AddConnects();
    LoadSolutionFromMeshes();
    
    int nel_res = NElements();
    for (long el = 0; el < nel_res; el++) {
        TPZCompEl *cel = Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->PrepareIntPtIndices();
    }
    
}
void TPZMultiphysicsCompMesh::BuildMultiphysicsSpaceWithMemory(TPZVec<MSpaceConfig> & active_approx_spaces, TPZVec<TPZCompMesh * > & mesh_vector, std::set<int> matsIdWithMem, std::set<int> matsIdNoMem){
    m_active_approx_spaces = active_approx_spaces;
    m_mesh_vector          = mesh_vector;
    if (m_mesh_vector.size() != m_active_approx_spaces.size()) {
        std::cout<< "TPZMultiphysicsCompMesh:: The vector provided should have the same size." << std::endl;
        DebugStop();
    }
    
    int n_approx_spaces = m_mesh_vector.size();
    
    SetNMeshes(n_approx_spaces);
    Reference()->ResetReference();
    SetAllCreateFunctionsMultiphysicElemWithMem();
    ApproxSpace().CreateWithMemory(true);
    // delete all elements and connects in the mesh
    CleanElementsConnects();
    TPZCompMesh::AutoBuild(matsIdWithMem);
    SetAllCreateFunctionsMultiphysicElem();
    TPZCompMesh::AutoBuild(matsIdNoMem);
    SetAllCreateFunctionsMultiphysicElemWithMem();
    AddElements();
    AddConnects();
    LoadSolutionFromMeshes();
    
    int nel_res = NElements();
    for (long el = 0; el < nel_res; el++) {
        TPZCompEl *cel = Element(el);
        TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfcel) {
            continue;
        }
        mfcel->PrepareIntPtIndices();
    }
    
}

void TPZMultiphysicsCompMesh::BuildMultiphysicsSpace(TPZVec<TPZCompMesh * > & mesh_vector, const TPZVec<int64_t> &gelindexes){
    
    m_mesh_vector          = mesh_vector;
    m_active_approx_spaces.Resize(m_mesh_vector.size());
    for(int64_t i = 0; i< m_active_approx_spaces.size(); i++) m_active_approx_spaces[i] = ETestTrial;
    int n_approx_spaces = m_mesh_vector.size();
    SetNMeshes(n_approx_spaces);
    Reference()->ResetReference();
    // if(ApproxSpace().Style() != TPZCreateApproximationSpace::EMultiphysics)
    // {
    //     std::cout << __PRETTY_FUNCTION__ << " Modifying the style of approximation space "
    //     " to multiphysics\n";
    //     SetAllCreateFunctionsMultiphysicElem();
    // }
    // delete all elements and connects in the mesh
    CleanElementsConnects();
    TPZCompMesh::AutoBuild(gelindexes);
    AddElements();
    AddConnects();
    LoadSolutionFromMeshes();
}


void TPZMultiphysicsCompMesh::AutoBuild(){
    
    std::cout << __PRETTY_FUNCTION__ << " has not been implemented. Use BuildMultiphysicsSpace instead\n";
    DebugStop();
}

void TPZMultiphysicsCompMesh::AddElements(){

    TPZGeoMesh * geometry = Reference();
//    geometry->ResetReference();
    int64_t n_cels = NElements();
    int64_t n_gels = geometry->NElements();
    TPZVec<TPZCompEl *> Referred(n_gels,0);
    int n_approx_spaces = m_mesh_vector.size();
    for(int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        TPZCompMesh *atom = m_mesh_vector[i_as];
        if(!atom) continue;
        // load computational element references in the Referred vector
        {
            auto &elemvec = atom->ElementVec();
            auto &gelvec = geometry->ElementVec();
            int64_t nel = elemvec.NElements();
            for(int64_t el=0; el<nel; el++)
            {
                TPZCompEl *cel = elemvec[el];
                if(!cel) continue;
                int64_t gelindex = cel->ReferenceIndex();
                if(gelindex < 0) continue;
                Referred[gelindex] = cel;
            }
        }
//        m_mesh_vector[i_as]->LoadReferences();
        int64_t icel;
        // loop over the multiphysics elements
        for(icel=0; icel < n_cels; icel++)
        {
            TPZCompEl * cel = ElementVec()[icel];
            TPZMultiphysicsElement * mfcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
            if(mfcel)
            {
                int64_t found = 0;
                int64_t gelindex = mfcel->ReferenceIndex();
                TPZCompEl *celatom = Referred[gelindex];
                
                if (celatom) {
                    mfcel->AddElement(celatom, i_as);
                    continue;
                }
                else
                {
                    TPZGeoEl *gel = geometry->Element(gelindex);
                    TPZGeoEl *gelF = gel;
                    while(gelF->Father())
                    {
                        gelF = gelF->Father();
                        int gelFindex = gelF->Index();
                        if (Referred[gelFindex]) {
#ifdef PZDEBUG
                            if (gelF->MaterialId() != gel->MaterialId()) {
                                DebugStop();
                            }
#endif
                            mfcel->AddElement(Referred[gelFindex], i_as);
                            found = true;
                            break;
                        }
                    }
                }
                if (!found) {
                    mfcel->AddElement(0, i_as);
                }
            }
            else {
                DebugStop();
            }
        }
        Referred.Fill(0);
    }
    
    for (int64_t icel = 0; icel < n_cels; icel++) {
        TPZCompEl *cel = Element(icel);
        TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mfel) {
            continue;
        }
#ifdef PZDEBUG
        {
            int ncontained = 0;
            for(int i=0; i<n_approx_spaces; i++)
            {
                TPZCompEl *atcel = mfel->Element(i);
                if(atcel) ncontained++;
            }
            if(ncontained == 0)
            {
                TPZGeoEl *gel = cel->Reference();
                std::cout << "Multiphysics element " << icel << " with matid " <<
                    gel->MaterialId() << " does not refer to any elements "
                << " geometric index " << gel->Index() << std::endl;
                TPZCompMesh *flux = this->m_mesh_vector[0];
                fReference->ResetReference();
                flux->LoadReferences();
                std::cout << "Flux reference " << (void *) gel->Reference() << std::endl;
                DebugStop();
            }
        }
#endif
        mfel->SetActiveApproxSpaces(m_active_approx_spaces);
        mfel->InitializeIntegrationRule();
    }
    
}

void TPZMultiphysicsCompMesh::AddConnects(){
    
    int n_approx_spaces = m_mesh_vector.size();

    TPZVec<int64_t> FirstConnect(n_approx_spaces,0);
    int64_t nconnects = 0;
    for (int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        FirstConnect[i_as] = nconnects;
        nconnects += m_mesh_vector[i_as]->ConnectVec().NElements();
    }
    ConnectVec().Resize(nconnects);
    Block().SetNBlocks(nconnects);
    
    int64_t counter = 0;
    int64_t seqnum = 0;
    for (int i_as = 0; i_as < n_approx_spaces; i_as++)
    {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ic;
        int64_t nc = m_mesh_vector[i_as]->ConnectVec().NElements();
        for (ic=0; ic<nc; ic++)
        {
            TPZConnect &refcon =  m_mesh_vector[i_as]->ConnectVec()[ic];
            ConnectVec()[counter] = refcon;
            if (refcon.SequenceNumber() >= 0) {
                ConnectVec()[counter].SetSequenceNumber(seqnum);
                ConnectVec()[counter].SetNState(refcon.NState());
                ConnectVec()[counter].SetNShape(refcon.NShape());
                ConnectVec()[counter].SetLagrangeMultiplier(refcon.LagrangeMultiplier());
                int ndof = refcon.NDof(*m_mesh_vector[i_as]);
                Block().Set(seqnum,ndof);
                seqnum++;
            }
            counter++;
        }

        for (ic=0; ic<nc; ic++)
        {
            TPZConnect &cn = ConnectVec()[FirstConnect[i_as]+ic];
            if (cn.HasDependency())
            {
                TPZConnect::TPZDependBase *dep = cn.FirstDepend();
                while (dep) {
                    dep->fDepConnectIndex = dep->fDepConnectIndex+FirstConnect[i_as];
                    dep = dep->fNext;
                }
            }
        }
    }
    
    Block().SetNBlocks(seqnum);
    ExpandSolution();
    int64_t iel;
    int64_t nelem = NElements();
    for (iel = 0; iel < nelem; iel++)
    {
        
        
        TPZCompEl *celorig = ElementVec()[iel];
        TPZMultiphysicsElement *cel = dynamic_cast<TPZMultiphysicsElement *> (celorig);
        TPZMultiphysicsInterfaceElement *interfacee = dynamic_cast<TPZMultiphysicsInterfaceElement *> (celorig);
        if (interfacee) {
            continue;
        }
        if (!cel) {
            continue;
        }
        TPZStack<int64_t> connectindexes;
        int64_t imesh;
//        std::list<TPZOneShapeRestraint> oneshape;
        for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
            
            if (m_active_approx_spaces[i_as] == 0) {
                continue;
            }
            
            TPZCompEl *celref = cel->ReferredElement(i_as);
            if (!celref) {
                continue;
            }
//            std::list<TPZOneShapeRestraint> celrest;
//            celrest = celref->GetShapeRestraints();
//            for (std::list<TPZOneShapeRestraint>::iterator it = celrest.begin(); it != celrest.end(); it++) {
//                TPZOneShapeRestraint rest = *it;
//                TPZOneShapeRestraint convertedrest(rest);
//                for(int face = 0; face < rest.fFaces.size(); face++)
//                {
//                    int ic = rest.fFaces[face].first;
//                    convertedrest.fFaces[face].first = ic+FirstConnect[i_as];
//                }
//                oneshape.push_back(convertedrest);
//            }
            int64_t ncon = celref->NConnects();
            int64_t ic;
            for (ic=0; ic<ncon; ic++) {
                connectindexes.Push(celref->ConnectIndex(ic)+FirstConnect[i_as]);
            }
        }
        cel->SetConnectIndexes(connectindexes);
//        for (std::list<TPZOneShapeRestraint>::iterator it = oneshape.begin(); it != oneshape.end(); it++) {
//            cel->AddShapeRestraint(*it);
//        }
    }
    
}

void TPZMultiphysicsCompMesh::LoadSolutionFromMeshes()
{
    if(fSolType == EReal){
        return LoadSolutionFromMeshesInternal<STATE>();
    }
    if(fSolType == EComplex){
        return LoadSolutionFromMeshesInternal<CSTATE>();
    }
    PZError<<__PRETTY_FUNCTION__<<'\n';
    PZError<<"Invalid type! Aborting...\n";
    DebugStop();
}

template<class TVar>
void TPZMultiphysicsCompMesh::LoadSolutionFromMeshesInternal()
{
    int n_approx_spaces = m_active_approx_spaces.size();
    TPZManVector<int64_t> FirstConnectIndex(n_approx_spaces+1,0);
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        FirstConnectIndex[i_as+1] = FirstConnectIndex[i_as]+m_mesh_vector[i_as]->NConnects();
    }
    TPZBlock &blockMF = Block();
    TPZFMatrix<TVar> &solMF = Solution();
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ncon = m_mesh_vector[i_as]->NConnects();
        TPZBlock &block = m_mesh_vector[i_as]->Block();
        TPZFMatrix<TVar> &sol = m_mesh_vector[i_as]->Solution();
        int64_t ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = m_mesh_vector[i_as]->ConnectVec()[ic];
            int64_t seqnum = con.SequenceNumber();
            if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = ConnectVec()[FirstConnectIndex[i_as]+ic];
            int64_t seqnumMF = conMF.SequenceNumber();
            if(seqnumMF < 0) continue;
            for (int idf=0; idf<blsize; idf++) {
                solMF(blockMF.Index(seqnumMF, idf)) = sol(block.Index(seqnum, idf));
            }
        }
    }
}

void TPZMultiphysicsCompMesh::LoadSolutionFromMultiPhysics()
{
    if(fSolType == EReal){
        return LoadSolutionFromMultiPhysicsInternal<STATE>();
    }
    if(fSolType == EComplex){
        return LoadSolutionFromMultiPhysicsInternal<CSTATE>();
    }
    PZError<<__PRETTY_FUNCTION__<<'\n';
    PZError<<"Invalid type! Aborting...\n";
    DebugStop();
}

template<class TVar>
void TPZMultiphysicsCompMesh::LoadSolutionFromMultiPhysicsInternal()
{
    int n_approx_spaces = m_active_approx_spaces.size();
    TPZManVector<int64_t> FirstConnectIndex(n_approx_spaces+1,0);
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        FirstConnectIndex[i_as+1] = FirstConnectIndex[i_as]+m_mesh_vector[i_as]->NConnects();
    }
    TPZBlock &blockMF = Block();
    TPZFMatrix<TVar> &solMF = Solution();
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        
        int64_t ncon = m_mesh_vector[i_as]->NConnects();
        TPZBlock &block = m_mesh_vector[i_as]->Block();
        TPZFMatrix<TVar> &sol = m_mesh_vector[i_as]->Solution();
        int64_t ic;
        for (ic=0; ic<ncon; ic++) {
            TPZConnect &con = m_mesh_vector[i_as]->ConnectVec()[ic];
            int64_t seqnum = con.SequenceNumber();
            if(seqnum<0) continue;       /// Whether connect was deleted by previous refined process
            int blsize = block.Size(seqnum);
            TPZConnect &conMF = ConnectVec()[FirstConnectIndex[i_as]+ic];
            int nelconnected = conMF.NElConnected();
            if (nelconnected == 0) {
                continue;
            }
            int64_t seqnumMF = conMF.SequenceNumber();
            int idf;
            for (idf=0; idf<blsize; idf++) {
                sol(block.Index(seqnum, idf)) = solMF(blockMF.Index(seqnumMF, idf));
            }
        }
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Solutions of the referred meshes";
    }
#endif
    for (int i_as = 0; i_as < n_approx_spaces; i_as++) {
        if (m_active_approx_spaces[i_as] == 0) {
            continue;
        }
        m_mesh_vector[i_as]->LoadSolution((m_mesh_vector[i_as]->Solution()));
    }
}

/// delete the elements and connects
void TPZMultiphysicsCompMesh::CleanElementsConnects()
{
    int64_t nel = NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = Element(el);
        if(cel)
        {
            delete cel;
            fElementVec[el] = 0;
        }
    }
    fElementVec.Resize(0);
    nel = fConnectVec.NElements();
    for (int64_t el=0; el<nel; el++) {
        fConnectVec[el].RemoveDepend();
    }
    fConnectVec.Resize(0);
}
