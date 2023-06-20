/**
 * @file
 * @brief Contains implementations of the TPZSubMeshAnalysis methods: implementation of the TPZSubMeshAnalysis class.
 */

#include "pzsmanal.h"
#include "pzsubcmesh.h"
#include "pzfmatrix.h"
#include "TPZMatrixSolver.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.analysis.pzsmanalysis");

#endif

using namespace std;

// Construction/Destruction

TPZSubMeshAnalysis::TPZSubMeshAnalysis(TPZSubCompMesh *mesh) : TPZRegisterClassId(&TPZSubMeshAnalysis::ClassId),
TPZLinearAnalysis(mesh), fReducableStiff(0){
	fMesh = mesh;
    if (fMesh)
    {
        fReferenceSolution.Redim(fCompMesh->NEquations(),1);
    }
}

TPZSubMeshAnalysis::~TPZSubMeshAnalysis()
{
	
}

/** @brief Set the computational mesh of the analysis. */
void TPZSubMeshAnalysis::SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth)
{
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(mesh);
    if (submesh) {
        fMesh = submesh;
    }
    else
    {
        DebugStop();
    }
    TPZLinearAnalysis::SetCompMesh(mesh, mustOptimizeBandwidth);
    if (fCompMesh) {
        fReferenceSolution.Redim(fCompMesh->NEquations(), 1);
    }
}


template<class TVar>
void TPZSubMeshAnalysis::AssembleInternal()
{
    TPZCompMesh *mesh = Mesh();
    auto &mySolver = MatrixSolver<TVar>();
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *>(mesh);
    if (!submesh) {
        DebugStop();
    }
    TPZCompMesh *fathermesh = submesh->Mesh();

	int numeq = fCompMesh->NEquations();
	int numinternal = fMesh->NumInternalEquations();
	fReferenceSolution.Redim(numeq,1);
	fRhs.Redim(numeq,1);
    if(!fReducableStiff) 
    {
        fReducableStiff = new TPZMatRed<TVar, TPZFMatrix<TVar> > ();
    }
    if(!fStructMatrix->HasRange()){
        fReducableStiff->Redim(numeq,numinternal);
    }
	TPZMatRed<TVar, TPZFMatrix<TVar> > *matred = dynamic_cast<TPZMatRed<TVar, TPZFMatrix<TVar> > *> (fReducableStiff.operator->());
    if(!mySolver.Matrix())
    {
        if (fStructMatrix->HasRange()) {
            const int64_t nintactiveeq = fStructMatrix->EquationFilter().NActiveEquations();
            const int64_t nexteq = numeq - numinternal;
            const int64_t neqfiltered = nintactiveeq + nexteq;
            fReducableStiff->Redim(neqfiltered,nintactiveeq);
            
            // Creates the sparse structure with only the active equations of the internal dofs
            mySolver.SetMatrix(fStructMatrix->Create());
            
            // Adding external equations to active equations of filter
            TPZVec<int64_t> active(neqfiltered);
            int count = 0;
            for (auto eq : fStructMatrix->EquationFilter().GetActiveEquations()) {
                active[count++] = eq;
            }
            for (int ieq = numinternal; ieq < numeq ; ieq++) {
                active[count++] = ieq;
            }
            fStructMatrix->EquationFilter().Reset();
            fStructMatrix->EquationFilter().SetActiveEquations(active);
            fStructMatrix->EquationFilter().SetNumEq(numeq); // it does not hurt to set it again in case it was forgotten after condensing...
        }
        else{
            // fReducableStiff->Redim(numeq,numinternal);
            fStructMatrix->SetEquationRange(0, numinternal);
            mySolver.SetMatrix(fStructMatrix->Create());
            fStructMatrix->EquationFilter().Reset();
        }
    }
    matred->SetMaxNumberRigidBodyModes(fMesh->NumberRigidBodyModes());
	//	fReducableStiff.SetK00(fSolver->Matrix());
	// this will initialize fK00 too
	matred->SetSolver(dynamic_cast<TPZMatrixSolver<STATE> *>(mySolver.Clone()));
	//	TPZStructMatrix::Assemble(fReducableStiff,fRhs, *fMesh);
//	time_t before = time (NULL);
	fStructMatrix->Assemble(fReducableStiff,fRhs,fGuiInterface);
    
    if (fStructMatrix->HasRange()) {
        TPZFMatrix<TVar> rhsloc(fStructMatrix->EquationFilter().NActiveEquations(),fRhs.Cols(),0.);
        fStructMatrix->EquationFilter().Gather(fRhs, rhsloc);
        matred->SetF(rhsloc);
    }
    else{
        matred->SetF(fRhs);
    }
//	time_t after = time(NULL);
//	double diff = difftime(after, before);
//	std::cout << __PRETTY_FUNCTION__ << " tempo " << diff << std::endl;
}
void TPZSubMeshAnalysis::Assemble(){
	//TODOCOMPLEX
    return AssembleInternal<STATE>();
}

void TPZSubMeshAnalysis::Run(std::ostream &out){
	
	//fReducableStiff.Print("Reducable stiff before assembled");
	fReferenceSolution = fSolution;
//	time_t tempo = time(NULL);
	Assemble();
//	time_t tempodepois = time(NULL);
//	double elapsedtime = difftime(tempodepois, tempo);
	
//	std::cout << "Tempo para assemblagem " << elapsedtime << std::endl;
	if (!fReducableStiff) {
        DebugStop();
    }
	TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
    if(!matred)
    {
        DebugStop();
    }
}
void TPZSubMeshAnalysis::CondensedSolution(TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
//	time_t tempo = time(NULL);
	if (!fReducableStiff) {
        DebugStop();
    }
	TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        matred->Print("Before = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    matred->K11Reduced(ek, ef);
	
//	time_t tempodepois = time(NULL);
//	double elapsedtime = difftime(tempodepois, tempo);
	
//	std::cout << "Tempo para inversao " << elapsedtime << std::endl;
	
}

/** @brief compute the reduced right hand side using the current stiffness. Abort if there is no stiffness computed */
void TPZSubMeshAnalysis::ReducedRightHandSide(TPZFMatrix<STATE> &rhs)
{
    if (!fReducableStiff) {
        DebugStop();
    }
    TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        matred->Print("Before = ",sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    matred->SetF(fRhs);
    matred->F1Red(rhs);
    
}


void TPZSubMeshAnalysis::LoadSolution(const TPZFMatrix<STATE> &sol)
{
	
	//	sol.Print("sol");
	int numinter = fMesh->NumInternalEquations();
	int numeq = fMesh->TPZCompMesh::NEquations();
	TPZFMatrix<STATE> soltemp(numeq-numinter,1,0.);
	int i;
	for(i=0; i<numeq-numinter; i++) {
		soltemp(i,0) = sol.GetVal(numinter+i,0)-fReferenceSolution(numinter+i,0);
	}
	TPZFMatrix<STATE> uglobal(numeq,1,0.);
    if(fReducableStiff)
    {
        TPZMatRed<STATE, TPZFMatrix<STATE> > *matred = dynamic_cast<TPZMatRed<STATE, TPZFMatrix<STATE> > *> (fReducableStiff.operator->());
        if(this->StructMatrix()->EquationFilter().IsActive()){
            // If there is equations filter we need to scatter the solution to the "complete" solution vector
            TPZFMatrix<STATE> uglobaltemp(numeq,1,0.);
            matred->UGlobal(soltemp,uglobaltemp);
            this->StructMatrix()->EquationFilter().Scatter(uglobaltemp, uglobal);
        }
        else{
            matred->UGlobal(soltemp,uglobal);
        }
                
        fSolution = fReferenceSolution + uglobal;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        soltemp.Print("External DOF Solution",sout);
        uglobal.Print("Expanded solution",sout);
        fSolution.Print("fSolution",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	TPZLinearAnalysis::LoadSolution();
}

int TPZSubMeshAnalysis::ClassId() const{
    return Hash("TPZSubMeshAnalysis") ^ TPZLinearAnalysis::ClassId() << 1;
}
