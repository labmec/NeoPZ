/**
 * @file pzcondensedcompel.cpp
 * @brief Contains the implementations of the TPZCondensedCompEl methods.
 */

#include "pzcondensedcompel.h"
#include "pzlog.h"
#include "pzstepsolver.h"
#include "pzelementgroup.h"
#include "pzcmesh.h"
#include "TPZElementMatrixT.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzcondensedcompel");
#endif

#include "TPZLapack.h"


TPZCondensedCompEl::~TPZCondensedCompEl()
{
    if (fMesh->ElementVec()[fIndex] == this) {
        fMesh->ElementVec()[fIndex] = fReferenceCompEl;
        delete fReferenceCompEl;
    }
}


/** @brief unwrap the condensed element from the computational element and delete the condensed element */
void TPZCondensedCompEl::Unwrap()
{
    int64_t myindex = fIndex;
    fMesh->ElementVec()[myindex] = 0;
    TPZCompEl *ReferenceEl = fReferenceCompEl;
    int ncon = NConnects();
    for (int ic=0; ic<ncon ; ic++) {
        Connect(ic).SetCondensed(false);
    }
    TPZGeoEl *gel = Reference();
    if (gel) {
        gel->ResetReference();
    }
    delete this;
    ReferenceEl->Mesh()->ElementVec()[myindex] = ReferenceEl;
    if (gel) {
        gel->SetReference(ReferenceEl);
    }
}

/**
 * @brief Set the index i to node inode
 * @param inode node to set index
 * @param index index to be seted
 */
void TPZCondensedCompEl::SetConnectIndex(int inode, int64_t index)
{
    LOGPZ_ERROR(logger,"SetConnectIndex should never be called")
    DebugStop();
}


/** @brief Verifies if the material associated with the element is contained in the set */
bool TPZCondensedCompEl::HasMaterial(const std::set<int> &materialids) const {
    if(!fReferenceCompEl){
        return false;
    }
    bool has_material_Q = fReferenceCompEl->HasMaterial(materialids);
    return has_material_Q;
}

/**
 * @brief Prints element data
 * @param out Indicates the device where the data will be printed
 */
void TPZCondensedCompEl::Print(std::ostream &out) const
{
    out << "Output for a condensed element\n";
    TPZCompEl::Print(out);
    out << "Connect indices (Condensed or not)\n";
    for(int i=0; i<NConnects(); i++){
        TPZConnect &c = Connect(i);
        out << ConnectIndex(i) << "/" << c.IsCondensed() << ' ';
    }
    out << std::endl;
    TPZElementGroup *eg = dynamic_cast<TPZElementGroup *>(fReferenceCompEl);
    if(eg)
    {
        out << "Index of grouped elements: ";
        int nel = eg->GetElGroup().size();
        for(int i=0; i<nel-1; i++){
            out << eg->GetElGroup()[i]->Index() <<", ";
        }
        out << eg->GetElGroup()[nel-1]->Index() <<std::endl;
        out << "Connect indexes of the contained elements \n";
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = eg->GetElGroup()[i];
            TPZGeoEl *gel = cel->Reference();
            out << "cel index " << cel->Index() << " cindex ";
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                out << cel->ConnectIndex(ic) << " ";
            }
            if (gel) {
                out << "\ngelindex " << gel->Index() << " ";
                out << "matid " << gel->MaterialId();
                TPZManVector<REAL,3> xi(gel->Dimension()), xco(3);
                gel->CenterPoint(gel->NSides()-1, xi);
                gel->X(xi,xco);
                out << " center x " << xco;
            }
            out << std::endl;
        }
    }
    else
    {
        if (fReferenceCompEl) {
            fReferenceCompEl->Print(out);
        }
        else
        {
            DebugStop();
        }
    }
    out << "Internal index resequencing: " << fIndexes << std::endl;
//    fCondensed.Print("Condensed matrix",out,EMathematicaInput);
}


/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZCondensedCompEl::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    int nc = NConnects();
    std::set<int64_t> refconn;
    fReferenceCompEl->BuildCornerConnectList(refconn);
    for (int ic=0; ic<nc ; ic++) {
        TPZConnect &c = Connect(ic);
        if (!c.IsCondensed() || !c.HasDependency()) {
            int64_t index = ConnectIndex(ic);
            if (refconn.find(index) != refconn.end()) {
                connectindexes.insert(index);
            }
        }
    }
}

int TPZCondensedCompEl::ClassId() const{
    return Hash("TPZCondensedCompEl") ^ TPZCompEl::ClassId() << 1;
}


template<class TVar>
TPZCondensedCompElT<TVar>::TPZCondensedCompElT(TPZCompEl *ref, bool keepmatrix) :
TPZRegisterClassId(&TPZCondensedCompElT<TVar>::ClassId)
{
    fKeepMatrix = keepmatrix;
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    TPZGeoEl *gel = ref->Reference();
    if (gel) {
        SetReference(gel->Index());
    }
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(ref);
    if(elgr)
    {
        elgr->ReorderConnects();
    }

    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int ic=0; ic< ncon ; ic++) {
        fIndexes[ic] = ic;
    }
    Resequence();
#ifdef PZDEBUG
    if(elgr)
    {
        for (int ic=0; ic<ncon; ic++) {
            if(fIndexes[ic] != ic) DebugStop();
        }
        for (int ic=0; ic<fCondensedConnectIndexes.size(); ic++) {
            if(fCondensedConnectIndexes[ic] != elgr->ConnectIndex(ic))
            {
                DebugStop();
            }
        }
        int ncondense = fCondensedConnectIndexes.size();
        for (int ic=0; ic<fActiveConnectIndexes.size(); ic++) {
            if(fActiveConnectIndexes[ic] != elgr->ConnectIndex(ncondense+ic))
            {
                DebugStop();
            }
        }
    }
#endif
}

/** @brief create a copy of the condensed computational element in the other mesh */
template<class TVar>
TPZCondensedCompElT<TVar>::TPZCondensedCompElT(const TPZCondensedCompElT<TVar> &copy,
                                               TPZCompMesh &mesh)
    : TPZRegisterClassId(&TPZCondensedCompElT<TVar>::ClassId)
{
    TPZCompEl *ref = fReferenceCompEl->Clone(mesh);
    if(!ref)
    {
        DebugStop();
    }
    fReferenceCompEl = ref;
    fMesh = ref->Mesh();
    fKeepMatrix = copy.fKeepMatrix;
    SetReference(ref->Reference()->Index());
    SetIndex(ref->Index());
    fMesh->ElementVec()[fIndex] = this;
    int ncon = ref->NConnects();
    fIndexes.resize(ncon);
    for (int i=0; i<ncon ; ++i) {
        fIndexes[i] = i;
    }
}

/**
 * @brief Method for creating a copy of the element in a patch mesh
 * @param mesh Patch clone mesh
 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
 * @param gl2lcElMap map the computational elements
 */
/**
 * Otherwise of the previous clone function, this method don't
 * copy entire mesh. Therefore it needs to map the connect index
 * from the both meshes - original and patch
 */
template<class TVar>
TPZCompEl *TPZCondensedCompElT<TVar>::ClonePatchEl(TPZCompMesh &mesh,
                                std::map<int64_t,int64_t> & gl2lcConMap,
                                std::map<int64_t,int64_t> & gl2lcElMap) const
{
    TPZCompEl *cel = fReferenceCompEl->ClonePatchEl(mesh,gl2lcConMap,gl2lcElMap);
    auto *result = new TPZCondensedCompElT<TVar>(cel);
    result->fIndexes = fIndexes;
    result->fCondensed = fCondensed;
    return result;
}

template<class TVar>
void TPZCondensedCompElT<TVar>::Resequence()
{
    TPZStack<int> condensed;
    TPZStack<int> notcondensed;
    int nint=0,next=0;
    std::set<int64_t> depreceive;
    int ncon = fIndexes.size();
    // if a connect of the condensed element receives a contribution through connect dependency
    // this will be tracked by depreceive
    for (int ic=0; ic<ncon; ic++) {
        TPZConnect &c = fReferenceCompEl->Connect(ic);
        TPZConnect::TPZDependBase * dep = c.FirstDepend();
        while (dep) {
            depreceive.insert(dep->fDepConnectIndex);
            dep = dep->fNext;
        }
    }
    for (int i=0; i<ncon ; ++i) {
        TPZConnect &c = fReferenceCompEl->Connect(i);
        int64_t cindex = fReferenceCompEl->ConnectIndex(i);
        if(c.NElConnected() == 1 && c.HasDependency() == 0 && depreceive.find(cindex) == depreceive.end())
        {
            c.SetCondensed(true);
        }
        if (c.IsCondensed()) {
            if(c.HasDependency())
            {
                std::cout << "Not Implemented yet\n";
                DebugStop();
            }
            condensed.Push(i);
            nint += c.NDof();
        }
        else
        {
            notcondensed.Push(i);
            next += c.NDof();
        }
    }
    fNumInternalEqs = nint;
    fNumTotalEqs = nint+next;
    
    int ncond = condensed.size();
    fCondensedConnectIndexes.Resize(ncond);
    for (int i=0; i<ncond; ++i) {
        fIndexes[i] = condensed[i];
        fCondensedConnectIndexes[i] = fReferenceCompEl->ConnectIndex(condensed[i]);
    }
    fActiveConnectIndexes.Resize(notcondensed.size());
    for (int i=0; i<notcondensed.size(); ++i) {
        fIndexes[i+ncond] = notcondensed[i];
        fActiveConnectIndexes[i] = fReferenceCompEl->ConnectIndex(notcondensed[i]);
    }
    //TPZAutoPointer<TPZMatrix<TVar> > k00 = new TPZFMatrix<TVar>(nint,nint,0.);
    if (fKeepMatrix == false) {
        nint = 0;
    }
	TPZAutoPointer<TPZMatrix<TVar> > k00 = new TPZFMatrix<TVar>(nint, nint, 0.);
    //TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(k00);
    TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(k00);
    if(0)
    {
        TPZAutoPointer<TPZMatrix<TVar> > mat2 = k00->Clone();
        TPZStepSolver<TVar> *gmrs = new TPZStepSolver<TVar>(mat2);
        step->SetReferenceMatrix(mat2);
        gmrs->SetGMRES(20, 20, *step, 1.e-20, 0);
        TPZAutoPointer<TPZMatrixSolver<TVar> > autogmres = gmrs;
    }
    step->SetDirect(ELU);
    TPZAutoPointer<TPZMatrixSolver<TVar> > autostep = step;
    

    
    if(fKeepMatrix == true)
    {
//        fCondensed.Redim(nint+next,nint);
        fCondensed.Redim(fNumTotalEqs,fNumInternalEqs);
    }
    fCondensed.SetSolver(autostep);
}

/// Assemble the stiffness matrix in locally kept datastructure
template<class TVar>
void TPZCondensedCompElT<TVar>::Assemble()
{
    fCondensed.K00()->Redim(fNumInternalEqs, fNumInternalEqs);
    fCondensed.Redim(fNumTotalEqs, fNumInternalEqs);

    fCondensed.Zero();
    
    TPZElementMatrixT<TVar> ek,ef;
    
    fReferenceCompEl->CalcStiff(ek,ef);
    ek.PermuteGather(fIndexes);
    ef.PermuteGather(fIndexes);

    int64_t dim = ek.fMat.Rows();
    for (int64_t i=0; i<dim ; ++i) {
        for (int64_t j=0; j<dim ; ++j) {
            fCondensed(i,j) = ek.fMat(i,j);
        }
    }
    
    fCondensed.SetF(ef.fMat);
    fCondensed.SetReduced();
}


/**
 * @brief Computes the element stifness matrix and right hand side
 * @param ek element stiffness matrix
 * @param ef element load vector
 */
template<class TVar>
void TPZCondensedCompElT<TVar>::CalcStiff(TPZElementMatrixT<TVar> &ekglob,
                                          TPZElementMatrixT<TVar> &efglob)
{

    auto* material_ss =
        dynamic_cast<TPZMatSingleSpaceT<TVar> *>(this->Material());
    auto* material_cs =
        dynamic_cast<TPZMatCombinedSpacesT<TVar> *>(this->Material());
    /*
      apparently you can have a condensed group that will not have an
      associated material but still needs to compute a stiff matrix
     */
    const bool has_mat = material_ss || material_cs;
    
    auto *nullmat = dynamic_cast<TPZNullMaterial<TVar> *>(material_ss);
    auto *nullmatcs = dynamic_cast<TPZNullMaterialCS<TVar> *>(material_cs);
    if(has_mat && (nullmat || nullmatcs))
    {
        ekglob.Reset();
        efglob.Reset();
        ekglob.fType = TPZElementMatrix::EK;
        efglob.fType = TPZElementMatrix::EF;
        return;
    }
    
    if(fKeepMatrix == false)
    {
        fKeepMatrix = true;
        fCondensed.K00()->Redim(fNumInternalEqs, fNumInternalEqs);
        fCondensed.Redim(fNumTotalEqs, fNumInternalEqs);
        fKeepMatrix = false;
    }
    InitializeElementMatrix(ekglob, efglob);
    TPZElementMatrixT<TVar> ek, ef;
    
    fReferenceCompEl->CalcStiff(ek,ef);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        Print(sout);
        sout << "Connect indices of element stiffness" << ek.fConnect << std::endl;
        ek.fMat.Print("EKOrig = ",sout,EMathematicaInput);
        ef.fMat.Print("EFOrig = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(fReferenceCompEl);
    if(!elgr)
    {
        ek.PermuteGather(fIndexes);
        ef.PermuteGather(fIndexes);
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Element connects\n";
        int nc = NConnects();
        for (int ic=0; ic<nc; ic++) {
            sout << "ic = " << ic << ' ' << " index " << ConnectIndex(ic) << ' ';
            Connect(ic).Print(*Mesh(),sout);
        }
        sout << "Permutations " << fIndexes << std::endl;
        sout << "Connect indices " << ek.fConnect << std::endl;

        ek.fMat.Print("EKPermute = ",sout,EMathematicaInput);
        ef.fMat.Print("EFPermute = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
//    ek.fMat.Print("ek = ",std::cout,EMathematicaInput);
//    ef.fMat.Print("ef = ",std::cout,EMathematicaInput);
    
#ifdef USING_DGER
    
    // Initialization TPZMatRed structure
    fCondensed.Zero();
    {
        TPZFMatrix<TVar> * K00_temp = dynamic_cast<TPZFMatrix<TVar> * >(fCondensed.K00().operator->());
        K00_temp->InitializePivot();
    }
    int64_t dim0 = fCondensed.Dim0();
    int64_t dim1 = fCondensed.Dim1();
    int64_t rows = ek.fMat.Rows();
    int64_t cols = ek.fMat.Cols()+ef.fMat.Cols();
    
    TPZFMatrix<TVar> KF(rows,cols); //  Local object
    for(int64_t i=0; i<rows;i++)
        for (int64_t j=0; j<cols; j++)
        {
            if (j<rows)
                KF(i,j) = ek.fMat(i,j);
            else
                KF(i,j) = ef.fMat(i,j-rows);
        }
    
    for (int64_t i=0; i<dim1; i++)
        for (int64_t j=0; j<dim0; j++)
            fCondensed.K10().operator()(i,j)=KF(i+dim0,j);
    
    for (int64_t i=0; i<dim1; i++)
        for (int64_t j=0; j<dim1; j++)
            fCondensed(i+dim0,j+dim0)=KF(i+dim0,j+dim0);
    
    fCondensed.SetF(ef.fMat);
    
    // LDLt Decomposition
    for (int64_t i=0; i<rows-dim1; i++)
    {
        for(int64_t j=i+1;j<cols;j++)
        {
            if (j<rows)
            {
                KF(j,i)/=KF(i,i);
                KF(i,j)/=KF(i,i);
            }
            else
                KF(i,j)/=KF(i,i);
        }
        
#ifdef TVardouble
        cblas_dger (CblasColMajor, rows-i-1, cols-i-1,
                    -KF(i,i), &KF(i+1,i), 1,
                    &KF(i,i+1), rows, &KF(i+1,i+1), rows);
#else
        DebugStop();
#endif
    }
    
    for (int64_t i=dim0; i< rows; i++)
    {
        ef.fMat(i,0) = KF.GetVal(i,fCondensed.Rows());
        for (int64_t j=dim0; j< fCondensed.Rows(); j++)
        {
            fCondensed(i,j) = KF.GetVal(i,j);
        }
    }
    
    TPZAutoPointer<TPZMatrix<TVar> > K00 = fCondensed.K00();
    for (int64_t i=0; i<dim0; i++)
        for (int64_t j=0; j<dim0; j++)
            K00->operator()(i, j) = KF(i,j);
    
    fCondensed.K00()->SetIsDecomposed(ELDLt);
    
    //    if(0){// Implementation not complete
    //
    //        for (int64_t i=0; i<dim0; i++){ // Substituindo valores obtidos para K01 usando o BLAS
    //            for (int64_t j=0; j<dim1; j++){
    //                fCondensed.K01().operator()(i,j) = KF(i,j+dim0)/KF(i,i);
    //            }
    //
    //        }
    //
    //        double scale = 1.0;
    //        cblas_dtrsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasUnit,dim0,dim1,scale,&fCondensed.K00().operator->()->operator()(0,0),dim0,&fCondensed.K01()(0,0),dim0);
    //    }else{
    fCondensed.K00()->SolveDirect(fCondensed.K01(), ELDLt);
    fCondensed.K00()->SolveDirect(fCondensed.F0(), ELDLt);
    //    }
    fCondensed.SetF(ef.fMat);
    fCondensed.SetK01IsComputed(true);
    fCondensed.SetF0IsComputed(true);
    fCondensed.SetReduced();// Directive that instructs the object to behave as reduced matrix.

    int64_t dim = ek.fMat.Rows();
    for (int64_t i = dim0; i<dim; i++) {
        ef.fMat(i,0) = fCondensed.F1()(i-dim0,0);
        for (int64_t j = dim0; j<dim; j++) {
            ek.fMat(i,j) = fCondensed.K11()(i-dim0,j-dim0);
        }
    }
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        fCondensed.Print("Condensed = ",sout);
        ek.fMat.Print("Rigidez",sout);
        ef.fMat.Print("Carga",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#else
    
    fCondensed.Zero();
    
    int64_t dim = ek.fMat.Rows();
    for (int64_t i=0; i<dim ; ++i) {
        for (int64_t j=0; j<dim ; ++j) {
            fCondensed(i,j) = ek.fMat.g(i,j);
        }
    }
    
    fCondensed.SetF(ef.fMat);
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "BEFORE CONDENSING THE EQUATIONS\n";
        sout << "Index = " << Index() << std::endl;
        fCondensed.Print("Reduced = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    int64_t dim1 = fCondensed.Dim1();
    TPZFNMatrix<200,TVar> K11(dim1,dim1),F1(dim1,ef.fMat.Cols());
    //const TPZFMatrix<REAL> &k11 = fCondensed.K11Red();
    
	fCondensed.K11Reduced(K11, F1);

#ifdef PZDEBUG
    {
        REAL normk01 = Norm(fCondensed.K01());
        REAL normf0 = Norm(fCondensed.F0());
        REAL normK11 = Norm(K11);
        if(std::isnan(normk01) || std::isnan(normf0) || std::isnan(normK11))
        {
            Print();
            DebugStop();
        }
    }
#endif
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Index = " << Index() << std::endl;
        fCondensed.Print("Reduced = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fCondensed.SetReduced();
    
    //const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    // this will change. Only fill in the active connects
    dim = K11.Rows();
#ifdef PZDEBUG
    if(dim != ekglob.fMat.Rows())
    {
        DebugStop();
    }
#endif
    for (int64_t i=0; i<dim; i++) {
        efglob.fMat(i,0) = F1.g(i,0);
        for (int64_t j=0; j<dim; j++) {
            ekglob.fMat.g(i,j) = K11.g(i,j);
        }
    }
    
    
#endif
    
#ifdef USING_DGER
#ifdef USING_LAPACK
    TPZFMatrix<TVar> * K00_temp = dynamic_cast<TPZFMatrix<TVar> * >(fCondensed.K00().operator->());
    K00_temp->InitializePivot();
#endif
#endif
    
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << std::endl;
        int nc = NConnects();
        for (int ic=0; ic<nc; ic++) {
            sout << "ic = " << ic << ' ';
            Connect(ic).Print(*Mesh(),sout);
        }
        ekglob.fMat.Print("EK11Reduced",sout,EMathematicaInput);
        efglob.fMat.Print("EF11Reduced",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

//    fCondensed.Print("Condensed = ",std::cout);
//    fCondensed.K11().Print("kbar = ",std::cout,EMathematicaInput);
//    fCondensed.F1().Print("fbar = ",std::cout,EMathematicaInput);
//    fCondensed.K01().Print("k01 = ",std::cout,EMathematicaInput);
//    fCondensed.F0().Print("f0 = ",std::cout,EMathematicaInput);
    
    if (fKeepMatrix == false) {
        fCondensed.K00()->Redim(0, 0);
        fCondensed.K10().Redim(0,0);
        fCondensed.K11().Redim(0,0);
    }
}


/**
 * @brief Computes the element right hand side
 * @param ef element load vector(s)
 */
template<class TVar>
void TPZCondensedCompElT<TVar>::CalcResidual(TPZElementMatrixT<TVar> &ef)
{
    // we need the stiffness matrix computed to compute the residual
    if (fKeepMatrix == false) {
        DebugStop();
        fKeepMatrix = true;
        fKeepMatrix = false;
    }
    fReferenceCompEl->CalcResidual(ef);
    ef.PermuteGather(fIndexes);
    fCondensed.SetF(ef.fMat);
    //const TPZFMatrix<REAL> &f1 = fCondensed.F1Red();
    TPZFNMatrix<100,TVar> f1(fCondensed.Dim1(),ef.fMat.Cols());
	fCondensed.F1Red(f1);
    int64_t dim1 = f1.Rows();
    int64_t dim = ef.fMat.Rows();
    int64_t dim0 = dim-dim1;
    for (int64_t i= dim0; i<dim; i++) {
        ef.fMat(i,0) = f1.GetVal(i-dim0,0);
    }
    if (fKeepMatrix == false) {
        fCondensed.Redim(0,0);
    }
}

/** @brief Loads the solution within the internal data structure of the element */
/**
 * Is used to initialize the solution of connect objects with dependency \n
 * Is also used to load the solution within SuperElements
 */
template<class TVar>
void TPZCondensedCompElT<TVar>::LoadSolution()
{
//    if (fKeepMatrix == false) {
//        fKeepMatrix = true;
//        fCondensed.K00()->Redim(fNumInternalEqs, fNumInternalEqs);
//        fCondensed.Redim(fNumTotalEqs, fNumInternalEqs);
//        TPZElementMatrix ek,ef;
//        CalcStiff(ek, ef);
//        fKeepMatrix = false;
//    }
    // initialize the solution of the constrained connects
    TPZCompEl::LoadSolution();
    
    // if the matrix has not been condensed then nothing to do
    if (fCondensed.Rows() != fCondensed.Dim1())
    {
        return;
    }
    // compute the solution of the internal equations
    int dim0=0, dim1=0;
    int nc = fIndexes.size(),nc0 = fCondensedConnectIndexes.size(), nc1 = fActiveConnectIndexes.size();
//    int ic;
//    for (ic=0; ic<nc ; ic++) {
//        int64_t connectindex = fIndexes[ic];
//        TPZConnect &con = fReferenceCompEl->Connect(connectindex);
//        int sz = con.NShape()*con.NState();
//        if (con.IsCondensed()) {
//#ifdef PZDEBUG
//            if (dim1) {
//                DebugStop();
//            }
//#endif
//            dim0 += sz;
//        }
//        else
//        {
//            dim1 += sz;
//        }
//    }
    dim0 = fNumInternalEqs;
    dim1 = fNumTotalEqs-dim0;
    //TPZBlock &bl = Mesh()->Block();
	TPZBlock &bl = Mesh()->Block();
    TPZFMatrix<TVar> &sol = Mesh()->Solution();
    int64_t count = 0;
    //TPZFMatrix<REAL> u1(dim1,1,0.);
	TPZFMatrix<TVar> u1(dim1,1,0.);
    //TPZFMatrix<REAL> elsol(dim0+dim1,1,0.);
	TPZFNMatrix<60,TVar> elsol(dim0+dim1,1,0.);
    for (int ic=0; ic<nc1 ; ic++) {
        TPZConnect &c = Connect(ic);
        int64_t seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        if(blsize)
        {
            int64_t firsteq = bl.Position(seqnum);
            for (int ibl=0; ibl<blsize; ibl++) {
                u1(count++,0) = sol(firsteq+ibl,0);
            }
        }
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        TPZManVector<int64_t,60> u1eq(dim1);
        int64_t count = 0;
        for (int ic=0; ic<nc1 ; ic++) {
            TPZConnect &c = Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            int blsize = bl.Size(seqnum);
            if(blsize)
            {
                int64_t firsteq = bl.Position(seqnum);
                for (int ibl=0; ibl<blsize; ibl++) {
                    u1eq[count++] = firsteq+ibl;
                }
            }
        }
        std::stringstream sout;
        sout << "Computing UGlobal Index " << Index();
        sout << " Norm fK01 " << Norm(fCondensed.K01()) << std::endl;
        TPZVec<TVar> u1vec(dim1);
        for(int i=0; i<u1vec.size(); i++) u1vec[i] = u1(i,0);
        sout << "u1 " << u1vec << std::endl;
        sout << "u1 eq " << u1eq << std::endl;
        for(int i=0; i<dim1; i++) sout << sol(u1eq[i],0) << ' ';
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fCondensed.UGlobal(u1, elsol);
    count = 0;
    for (int ic=0; ic<nc0 ; ic++) {
        int64_t cindex = fCondensedConnectIndexes[ic];
        TPZConnect &c = Mesh()->ConnectVec()[cindex];
        int64_t seqnum = c.SequenceNumber();
        int blsize = bl.Size(seqnum);
        for (int ibl=0; ibl<blsize; ibl++) {
            sol(bl.Index(seqnum,ibl),0) = elsol(count++,0);
        }
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "After Computing UGlobal Index " << Index() ;
        sout << " Norm fK01 " << Norm(fCondensed.K01()) << std::endl;
        TPZVec<TVar> u1vec(dim1+dim0);
        for(int i=0; i<u1vec.size(); i++) u1vec[i] = elsol(i,0);
        sout << "elsol " << u1vec;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fReferenceCompEl->LoadSolution();
//    if (fKeepMatrix == false) {
//        fCondensed.Redim(0,0);
//        fCondensed.K00()->Redim(0, 0);
//    }

}

template<class TVar>
void TPZCondensedCompElT<TVar>::PermuteActiveConnects(TPZManVector<int64_t> &perm)
{
    auto ncon = NConnects();
    auto nint = fNumInternalEqs;
    for (int64_t ic = 0; ic < ncon; ic++)
    {
        fActiveConnectIndexes[ic] = perm[ic];
    }
    
    if (fKeepMatrix == false)
    {
        nint = 0;
    }
	TPZAutoPointer<TPZMatrix<TVar> > k00 = new TPZFMatrix<TVar>(nint, nint, 0.);
    TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(k00);
    step->SetDirect(ELU);
    TPZAutoPointer<TPZMatrixSolver<TVar> > autostep = step;
    
    fCondensed.SetSolver(autostep);
    if(fKeepMatrix == true)
    {
        fCondensed.Redim(fNumTotalEqs,fNumInternalEqs);
    }
}

template<class TVar>
int TPZCondensedCompElT<TVar>::ClassId() const{
    return Hash("TPZCondensedCompElT")
        ^ ClassIdOrHash<TVar>() << 1
        ^ TPZCondensedCompEl::ClassId() << 2;
}

template class TPZCondensedCompElT<STATE>;
template class TPZCondensedCompElT<CSTATE>;
