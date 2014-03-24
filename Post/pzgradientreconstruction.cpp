//
//  pzgradientreconstruction.cpp
//  PZ
//
//  Created by Agnaldo Farias on 4/10/13.
//
//

#include "pzgradientreconstruction.h"
#include "pzgradient.h"
#include "tpzintpoints.h"
#include "pzmultiphysicselement.h"
#include "pzmaterial.h"
#include "pzskylstrmatrix.h"
#include "pzintel.h"
#include "pzgnode.h"
#include "pzstepsolver.h"
#include <cmath>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.pzgradientreconstruction"));
#endif


using namespace std;

TPZGradientReconstruction::TPZGradientReconstruction(bool distmesh, REAL param){
    
    if (param != 1. && param !=2.) {
        DebugStop();
    }
    fGradData = new TPZGradientData;
    fDistortedMesh = distmesh;
    fparam = param;
}

TPZGradientReconstruction::~TPZGradientReconstruction(){
    
}

TPZGradientReconstruction::TPZGradientReconstruction(const TPZGradientReconstruction &cp){
    
    fGradData = cp.fGradData;
    fDistortedMesh = cp.fDistortedMesh;
    fparam = cp.fparam;
}

TPZGradientReconstruction &TPZGradientReconstruction::operator=(const TPZGradientReconstruction &copy){
    
    fGradData = copy.fGradData;
    fDistortedMesh = copy.fDistortedMesh;
    fparam = copy.fparam;
    return *this;
}

void TPZGradientReconstruction::ProjectionL2GradientReconstructed(TPZCompMesh *cmesh, int matidl2proj)
{
    if (cmesh->Reference()->Reference() != cmesh) {
        DebugStop();
    }
    
    // Redimensionando a matriz dos dados da reconstruca de gradientes
    int dim  = cmesh->Dimension();
    int nelem = cmesh->NElements();
    
    bool useweight;
    REAL paramK;
    GetDataDistortedMesh(useweight, paramK);
    
    //criar ponteiro para TPZFunction
    TPZGradient *pGrad = new TPZGradient;
    TPZAutoPointer<TPZFunction<STATE> > fp(pGrad);
    
    //Criar matrix de rigidez e vetor de carga
    int numloadcases=0;
    std::map<int, TPZMaterial * >::const_iterator mit;
	for(mit=cmesh->MaterialVec().begin(); mit!= cmesh->MaterialVec().end(); mit++) {
        TPZMaterial *mat = mit->second;
        if (!mat) {
            DebugStop();
        }
        numloadcases = mat->NumLoadCases();
        break;
    }
    
    int neq = cmesh->NEquations();
    TPZFMatrix<STATE> rhs;
    rhs.Redim(neq,numloadcases);
    
    TPZSkylineStructMatrix stmatrix(cmesh);
    TPZMatrix<STATE> *stiffmatrix = stmatrix.Create();
    
    int matid;
    for(int i=0; i<nelem; i++)
    {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel || cel->Dimension()!=dim) continue;
        
        TPZElementMatrix ek(cel->Mesh(), TPZElementMatrix::EK);
        TPZElementMatrix ef(cel->Mesh(), TPZElementMatrix::EF);
        
        fGradData->SetCel(cel, useweight, paramK);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            fGradData->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
        //set data of the gradient reconstructed
        TPZManVector<REAL,3> centerPoint, gradient;
        STATE cellAverage, slopeLimiter;
        fGradData->GetData(centerPoint, gradient, cellAverage, slopeLimiter);
        pGrad->SetData(centerPoint,gradient,cellAverage,slopeLimiter);
        
        //change material id current to material id of the L2Projection
        matid = cel->Material()->Id();
        ChangeMaterialIdIntoCompElement(cel, matid, matidl2proj);
        
        //set forcing function of l2 projection material
        TPZMaterial *mat = cel->Material();
        mat->SetForcingFunction(fp);
        
        //load the matrix ek and vector ef of the element
        cel->CalcStiff(ek,ef);
        
        //        ek.fMat.Print("ek = ");
        //        ef.fMat.Print("ef = ");
        
        //assemble pos l2 projection
        AssembleGlobalMatrix(cel, ek, ef, *stiffmatrix, rhs);
        
        //        stiffmatrix->Print("Matriz de Rigidez: ");
        //        rhs.Print("Right Handside: ");
        
        //Return for original material and current solution of the mesh
        ChangeMaterialIdIntoCompElement(cel, matidl2proj, matid);
    }
    
    //Solve linear system and transfer the solution to computational mesh
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    step.SetMatrix(stiffmatrix);
    TPZFMatrix<STATE> result;
    step.Solve(rhs, result);
    cmesh->Solution().Zero();
    cmesh->LoadSolution(result);
    
    //    stiffmatrix->Print("MatKRG = ");
    //    rhs.Print("FComRG = ");
    //    result.Print("SolComRG = ");
}

void TPZGradientReconstruction::ChangeMaterialIdIntoCompElement(TPZCompEl *cel, int oldmatid, int newmatid) {
    
    // Changes material Id only elements with required id (matid)
    if(cel->Material()->Id() != oldmatid) return;
    
    //mudar o material id
    TPZGeoEl *gel;
    gel = cel->Reference();
    gel->SetMaterialId(newmatid);
}


void TPZGradientReconstruction::AssembleGlobalMatrix(TPZCompEl *el, TPZElementMatrix &ek, TPZElementMatrix &ef,TPZMatrix<STATE> & stiffmatrix, TPZFMatrix<STATE> &rhs){
    
    if(!el->HasDependency()) {
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
        
    } else {
        // the element has dependent nodes
        ek.ApplyConstraints();
        ef.ApplyConstraints();
        ek.ComputeDestinationIndices();
        
        stiffmatrix.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
        rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
    }
}

TPZGradientReconstruction::TPZGradientData::TPZGradientData()
{
    fSolCellAndNeighbors.resize(0);
    fCenterPointCellAndNeighbors.resize(0);
    fCelAndNeighbors.resize(0);
    fCenterPointInterface.resize(0);
    
    this->fForcingFunctionExact = NULL;
    this->fUseForcinfFuncion = false;
    
    fWeightsGrad.resize(0);
    fGradient.resize(0);
    fdim=0;
    fSlopeLimiter=0.;
    fUseWeight=false;
    fparamK = 0;
    
}

TPZGradientReconstruction::TPZGradientData::~TPZGradientData()
{
    
}


TPZGradientReconstruction::TPZGradientData::TPZGradientData(const TPZGradientData &cp)
{
    fdim = cp.fdim;
    fCelAndNeighbors = cp.fCelAndNeighbors;
    fWeightsGrad = cp.fWeightsGrad;
    fUseWeight = cp.fUseWeight;
    fparamK = cp.fparamK;
    
    fSolCellAndNeighbors = cp.fSolCellAndNeighbors;
    fCenterPointCellAndNeighbors = cp.fCenterPointCellAndNeighbors;
    fCenterPointInterface = cp.fCenterPointInterface;
    fGradient = cp.fGradient;
    fSlopeLimiter = cp.fSlopeLimiter;
    
    fForcingFunctionExact = cp.fForcingFunctionExact;
    fUseForcinfFuncion = cp.fUseForcinfFuncion;
}

TPZGradientReconstruction::TPZGradientData & TPZGradientReconstruction::TPZGradientData::operator=(const TPZGradientData &copy)
{
    fdim = copy.fdim;
    fCelAndNeighbors = copy.fCelAndNeighbors;
    fWeightsGrad = copy.fWeightsGrad;
    fUseWeight = copy.fUseWeight;
    fparamK = copy.fparamK;
    
    fSolCellAndNeighbors = copy.fSolCellAndNeighbors;
    fCenterPointCellAndNeighbors = copy.fCenterPointCellAndNeighbors;
    fCenterPointInterface = copy.fCenterPointInterface;
    fGradient = copy.fGradient;
    fSlopeLimiter = copy.fSlopeLimiter;
    
    fForcingFunctionExact = copy.fForcingFunctionExact;
    fUseForcinfFuncion = copy.fUseForcinfFuncion;
    return *this;
}


void TPZGradientReconstruction::TPZGradientData::SetCel(TPZCompEl * cel, bool useweight, REAL paramK)
{
    if(!cel || cel->Dimension()!=cel->Mesh()->Dimension())
    {
        std::cout << "\nError: Element does not exist or has dimension different to the geometric mesh\n";
		DebugStop();
    }
    
    fSolCellAndNeighbors.resize(0);
    fCenterPointCellAndNeighbors.resize(0);
    fCelAndNeighbors.resize(0);
    fCenterPointInterface.resize(0);
    
    fdim = cel->Dimension();
    fUseWeight = useweight;
    fparamK = paramK;
    
    InitializeGradData(cel);
    
    ComputeGradient();
    
    ComputeSlopeLimiter2();
}

void TPZGradientReconstruction::TPZGradientData::Print(std::ostream &out) const
{
    const char *name1 = "Cell average";
    const char *name2 = "Center Point";
    const char *name3 = "Gradient";
    const char *name4 = "Slope limiter";
    
    const char *name5 = "Cel";
    const char *name6 = "Neigh";
    const char *name7 = "Center Point Interface";
    
    char string[256];
    
    out<<"\n\n";
    sprintf(string, "\t%s", name1);
    out << string ;
    sprintf(string, "\t%s", name2);
    out << string ;
    sprintf(string, "\t\t%s", name3);
    out << string ;
    sprintf(string, "\t\t%s", name4);
    out << string ;
    out<<"\n";
    
    sprintf(string, "%s", name5);
    out << string ;
    
    sprintf(string, "\t%f", fSolCellAndNeighbors[0]);
    out << string ;
    
    int i, j;
    out<<"\t(";
    for(j=0; j<3; j++ ){
        sprintf(string, "%f", fCenterPointCellAndNeighbors[0][j]);
        out << string;
        if (j<2) out<<",";
    }
    out <<")";
    
    out<<"\t(";
    for(i=0; i<3; i++ ){
        sprintf(string, "%f", fGradient[i]);
        out << string;
        if (i<2) out<<",";
    }
    out <<")";
    
    sprintf(string, "\t%f", fSlopeLimiter);
    out << string<<"\n";
    
    
    for (i=1; i<fSolCellAndNeighbors.size(); i++)
    {
        sprintf(string, "%s", name6);
        out << string <<i;
        
        sprintf(string, "\t%f", fSolCellAndNeighbors[i]);
        out << string ;
        
        out<<"\t(";
        for(j=0; j<3; j++ ){
            sprintf(string, "%f", fCenterPointCellAndNeighbors[i][j]);
            out << string;
            if (j<2) out<<",";
        }
        out <<")\n";
    }
    
    //Interface
    out <<"\n\n";
    sprintf(string, "%s", name7);
    out << string <<"\n";
    for (i=0; i<fCenterPointInterface.size(); i++)
    {
        out<<"(";
        for(j=0; j<3; j++ ){
            sprintf(string, "%f", fCenterPointInterface[i][j]);
            out << string;
            if (j<2) out<<",";
        }
        out <<")\n";
    }
    
}


void TPZGradientReconstruction::TPZGradientData::GetCenterPointAndCellAveraged(TPZCompEl *cel, TPZManVector<REAL,3> &xcenter, STATE &solcel)
{
    // ---------- calculating center point -----------
    TPZGeoEl* gel = cel->Reference();
    TPZManVector<REAL> centerpsi(3,0.0);
    gel->CenterPoint(gel->NSides()-1,centerpsi);
    xcenter.Fill(0.);
    gel->X(centerpsi,xcenter);
    
    //-------- calculating cell averaged ------------
    int intOrder = cel->GetgOrder();
    TPZIntPoints *pointIntRule = ((TPZInterpolatedElement*)cel)->Reference()->CreateSideIntegrationRule((cel->Reference()->NSides())-1,intOrder);
    int it, npoints = pointIntRule->NPoints();
    REAL integral = 0.0;
    TPZManVector<REAL> point(3,0.);
    TPZManVector<REAL> xpoint(3,0.);
	REAL weight;
    REAL Area = gel->RefElVolume();
    
    for(it=0;it<npoints;it++)
    {
		pointIntRule->Point(it,point,weight);
		weight /= Area;
		cel->Reference()->X(point,xpoint);
        
        TPZVec<STATE> sol;
        if (this->HasForcingFunctionExact()){
            sol.Resize(1, 0.);
            this->fForcingFunctionExact->Execute(xpoint, sol);
        }
        else{
            cel->Solution(xpoint, 1, sol);
        }
        if(sol.size()!=1) {
            PZError << "TPZGradientReconstruction::TPZDataGradient: The number of solutions variable can not be other than 1.\n";
            DebugStop();
        }
#ifdef STATE_COMPLEX
        integral += weight*sol[0].real();
#else
        integral += weight*sol[0];
#endif
    }
    
    solcel = integral;
}

void TPZGradientReconstruction::TPZGradientData::InitializeGradData(TPZCompEl *cel)
{
    if(!cel){
        std::cout << "\nError: Element does not exist\n";
		DebugStop();
    }
    
    fCelAndNeighbors.Push(cel);
    
    TPZManVector<REAL,3> xcenter(3);
    STATE cellaveraged;
    
    //------ Solution and center point of the cel -----------------
    GetCenterPointAndCellAveraged(cel,xcenter,cellaveraged);
    fCenterPointCellAndNeighbors.Push(xcenter);
    fSolCellAndNeighbors.Push(cellaveraged);
    
    
    //-------- Solution and center point of the neighbors, and center point of the interfaces ---------
    TPZStack<TPZCompElSide> neighequal,neighsmaller;
    TPZCompElSide neighbigger;
    int nneighs=0;
    TPZManVector<REAL,3> point(3,0.);
    TPZManVector<REAL,3> xpoint(3,0.);
    
    int oldneighs=0;
    int newneighs=0;
    
    std::set<TPZCompEl *> interfaces;
    for(int side = cel->Reference()->NSides()-2; side > cel->Reference()->NCornerNodes()-1; side--)
    {
        neighequal.Resize(0);
        neighsmaller.Resize(0);
        
        TPZCompElSide celside(cel,side);
        celside.EqualLevelElementList(neighequal, 0, 0);//(neighs,0,0);
        neighbigger = celside.LowerLevelElementList(0);
        
        if(neighbigger){
            neighequal.Push(neighbigger);
        }
        
        celside.HigherLevelElementList(neighsmaller, 1, 1);
        
        nneighs = neighequal.NElements();
        int nneighsmaller = neighsmaller.size();
        if(nneighs && nneighsmaller)
        {
            DebugStop();
        }
        
        if(nneighs != 0)
        {
            //Loop on neighboring elements greater or equal level
            for(int i =0; i<nneighs; i++)
            {
                TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighequal[i].Element());
                if(!InterpEl) continue;
                
                //Do not assume neighbors by node
                if(neighequal[i].Side() <neighequal[i].Reference().Element()->NCornerNodes()) continue;
                
                //verificando se eh elemento de contorno
                if(neighequal[i].Element()->Dimension() == fdim-1)
                {
                    TPZGeoElSide gelside = celside.Reference();
                    gelside.CenterPoint(point);
                    gelside.X(point,xpoint);
                    fCenterPointInterface.Push(xpoint);
                    continue;
                }
                
                oldneighs = interfaces.size();
                interfaces.insert(neighequal[i].Element());
                newneighs = interfaces.size();
                
                if(newneighs > oldneighs)
                {
                    fCelAndNeighbors.Push(neighequal[i].Element());
                    
                    GetCenterPointAndCellAveraged(neighequal[i].Element(),xcenter,cellaveraged);
                    fCenterPointCellAndNeighbors.Push(xcenter);
                    fSolCellAndNeighbors.Push(cellaveraged);
                    
                    TPZGeoElSide gelside = celside.Reference();
                    gelside.CenterPoint(point);
                    gelside.X(point,xpoint);
                    fCenterPointInterface.Push(xpoint);
                }
            }
            
        }
        
        //Loop on neighboring lower level elements
        if(nneighsmaller != 0)
        {
            
            for(int i =0; i<nneighsmaller; i++)
            {
                
                TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighsmaller[i].Element());
                if(!InterpEl) continue;
                
                //Do not assume neighbors by node
                if(neighsmaller[i].Side() <neighsmaller[i].Reference().Element()->NCornerNodes()) continue;
                
                //verificando se eh elemento de contorno
                if(neighequal[i].Element()->Dimension() == fdim-1)
                {
                    TPZGeoElSide gelside = neighsmaller[i].Reference();
                    gelside.CenterPoint(point);
                    gelside.X(point,xpoint);
                    fCenterPointInterface.Push(xpoint);
                    continue;
                }
                
                oldneighs = interfaces.size();
                interfaces.insert(neighsmaller[i].Element());
                newneighs = interfaces.size();
                
                if(newneighs > oldneighs)
                {
                    fCelAndNeighbors.Push(neighsmaller[i].Element());
                    
                    GetCenterPointAndCellAveraged(neighsmaller[i].Element(),xcenter,cellaveraged);
                    fCenterPointCellAndNeighbors.Push(xcenter);
                    fSolCellAndNeighbors.Push(cellaveraged);
                    
                    TPZGeoElSide gelside = neighsmaller[i].Reference();
                    gelside.CenterPoint(point);
                    gelside.X(point,xpoint);
                    fCenterPointInterface.Push(xpoint);
                }
            }
        }
    }
}


#include <stdio.h>
#ifdef USING_LAPACK
#include <clapack.h>
#endif
#ifdef USING_BLAS
#include <cblas.h>
#endif

void TPZGradientReconstruction::TPZGradientData::ComputeGradient()
{
    if (fCenterPointInterface.size()<1 || fSolCellAndNeighbors.size()<1)
    {
        DebugStop();
    }
    
    int i, j;
	fGradient.Resize(3, 0.);
	
	//matrices to apply the method of least squares
	TPZFMatrix<REAL> DeltaXcenter;//xcenter cell menos xcenter neigh
	TPZFMatrix<REAL> DeltaXcenterTranspose;
	TPZFMatrix<REAL> DifSol;//sol cell menos sol neigh
    
    int nneighs = fCenterPointCellAndNeighbors.size()-1;
    
    DeltaXcenter.Redim(nneighs,fdim);
    DeltaXcenterTranspose.Redim(fdim,nneighs);
    DifSol.Redim(nneighs,1);
    
    for(i = 0; i < nneighs; i++)
    {
        for(j=0; j<fdim; j++)
        {
            DeltaXcenter(i,j) = fCenterPointCellAndNeighbors[0][j] - fCenterPointCellAndNeighbors[i+1][j];
        }
        DifSol(i,0) = fSolCellAndNeighbors[0] - fSolCellAndNeighbors[i+1];
    }
    
    //insert weight
    if(fUseWeight==true)
    {
        ComputeWeights(fparamK);
        InsertWeights(DeltaXcenter, DifSol);
    }
    
    TPZFMatrix<REAL> grad;
    grad.Redim(fdim, 1);
    
    if(nneighs==fdim)
    {
        DeltaXcenter.SolveDirect(DifSol,ELU);
        grad= DifSol;
    }
    else
    {
#ifdef USING_LAPACK
        //QR factorization
        QRFactorization(DeltaXcenter,DifSol);
        DeltaXcenter.SolveDirect(DifSol,ELU);
        grad=DifSol;
#else
        
        // Solving the system by least squares: DeltaXcenter_T * DifSol = DeltaXcenter_T * DeltaXcenter*Grad(u)
        TPZFMatrix<REAL> A(fdim,fdim);
        DeltaXcenter.Transpose(&DeltaXcenterTranspose);
        A = DeltaXcenterTranspose*DeltaXcenter;
        grad = DeltaXcenterTranspose*DifSol;
        A.SolveDirect(grad,ELU);
#endif
    }
    
    for (i = 0; i<grad.Rows(); i++)
    {
        fGradient[i] = grad(i,0);
    }
}

void TPZGradientReconstruction::TPZGradientData::QRFactorization(TPZFMatrix<REAL> &matA,TPZFMatrix<REAL> &vecb)
{
    //-------- Compute a QR factorization of a real	M-by-N matrix A = QR ---------
    //On exit, the elements on and above the diagonal of the array
    //contain the min(M,N)-by-N upper trapezoidal matrix R (R is
    //upper triangular if m >= n); the elements below the diagonal,
    //with the array TAU, represent the orthogonal matrix Q as a
    //product of min(m,n) elementary reflectors

#ifdef USING_LAPACK
    
    int m = matA.Rows();
    int n = matA.Cols();
    
    int lda = m; //the leading dimension of the matA
    double *tau = new double[n];//The scalar factors of the elementary reflectors
    int lwork = n;
    double *work = new double[n];
    int info;
    
    double *A = new double[m*n];
    for(int j = 0; j<n; j++){
        for(int i=0; i<m; i++){
            A[i+m*j] = matA(i,j);
        }
    }
    
    dgeqrf_(&m, &n, A, &lda,tau,work,&lwork,&info);
    
    //matrix R: upper triangular
    matA.Redim(n, n);
    for(int j = 0; j<n; j++)
    {
        for(int i=0; i<n; i++)
        {
            if(i<=j) matA(i,j)=A[i+m*j];
        }
    }
    //matA.Print("\n\n \tmtR = ");
    
    //metodo que retorna a matrix Q
    int kk=n;
    double *Q = new double[m*n];
    Q=A;

    dorgqr_(&m, &n, &kk, Q, &lda,tau,work,&lwork,&info);
    
    TPZFMatrix<REAL> matQ;
    matQ.Redim(m, n);
    for(int j = 0; j<n; j++){
        for(int i=0; i<m; i++){
            matQ(i,j) = Q[i+m*j];
        }
    }
    //matQ.Print("\n\n \tmQ = ");
    TPZFMatrix<REAL> matQt;
    matQ.Transpose(&matQt);
    TPZFMatrix<REAL> res;
    
    matQt.Multiply(vecb, res);
    vecb.Redim(res.Rows(), res.Cols());
    vecb = res;
    
#else
    DebugStop();
#endif
    
    //------- Product of the transpose of matrix Q by b: Qˆt*vecb -------
    /*
     char side = 'L';
     char Trans = 'T';
     int nb = vecb.Cols();
     int k = n;
     
     int ldc = m;
     double work2[nb];
     int lwork2 =nb;
     
     double *b = new double[m];
     for (int i = 0; i<m; i++) {
     b[i]=vecb(i,0);
     }
     
     //o metodo nao explicita a matriz Q. Retorna Qˆt*b diretamente;
     dormqr_(&side, &Trans, &m, &nb, &k, A, &lda, tau, b, &ldc, work2, &lwork2, &info);
     
     vecb.Redim(n, 1);
     for(int j = 0; j<n; j++){
     vecb(j,0)=b[j];
     }
     vecb.Print("\n\n \tmQˆt*b = ");
     */
}

//Finite volume methods:foundation and analysis (chapter 3.3),Timothy Barth and Mario Ohlberge-2004
void TPZGradientReconstruction::TPZGradientData::ComputeSlopeLimiter()
{
    if(fGradient.size()==0 || fSolCellAndNeighbors.size()==0 || fCenterPointInterface.size()==0){
        PZError << "TPZGradientReconstruction::TPZDataGradient: gradient has size equal zero.\n";
        DebugStop();
    }
    
    int nneigh = fSolCellAndNeighbors.size()-1;
    if(nneigh==0) DebugStop();
    int i, j;
    
    //-------- getting min and max solution of neighbors ----
    STATE solKmax, solKmin;
    solKmin = fSolCellAndNeighbors[1];
    solKmax = fSolCellAndNeighbors[1];
    
    for(i = 2; i<= nneigh; i++)
    {
        if(solKmin > fSolCellAndNeighbors[i])
        {
            solKmin = fSolCellAndNeighbors[i];
            continue;
        }
        
        if(solKmax < fSolCellAndNeighbors[i])
        {
            solKmax = fSolCellAndNeighbors[i];
        }
    }
    
    if(solKmax < solKmin) DebugStop();
    
    // ----------- Calculating slope limiter --------------
    int ninterf = fCenterPointInterface.size();
    STATE solKside;
    STATE temp;
    TPZStack<STATE> alphavec;
    STATE solcel = fSolCellAndNeighbors[0];
    
    for (i = 0; i<ninterf; i++)
    {
        solKside = solcel;
        for(j=0; j<fdim; j++)
        {
            solKside += (STATE)(fCenterPointInterface[i][j] - fCenterPointCellAndNeighbors[0][j])*fGradient[j];
        }
        
        
        if(IsZero(solKside - solKmax) || IsZero(solKside - solKmin)) {
            temp = 1.;
        }
        else if((solKside - solKmax) > 1.e-12)
        {
            temp = (solKmax - solcel)/(solKside-solcel);
            if(temp>1.) temp = 1.;
            if(temp<0.) temp = fabs(temp);
        }
        else if((solKside - solKmin) < 1.e-12)
        {
            temp = (solKmin - solcel)/(solKside-solcel);
            if(temp>1.) temp = 1.;
            if(temp<0.) temp = fabs(temp);
        }
        else temp = 1.;
        
        if(temp < 0. || temp > 1.) DebugStop();
        alphavec.Push(temp);
    }
    
    //getting min slope limiter
    STATE alphaK = alphavec[0];
    for (long j=1; j<alphavec.size(); j++)
    {
        if(alphaK > alphavec[j])
        {
            alphaK = alphavec[j];
        }
    }
    fSlopeLimiter = alphaK;
    
    if(alphaK<0. || alphaK>1.) {
        DebugStop();
    }
}

//Paper: Development of a cell centred upwind finite volume algorithm for a new conservation law formulation in structural dynamics
void TPZGradientReconstruction::TPZGradientData::ComputeSlopeLimiter2()
{
    if(fGradient.size()==0 || fSolCellAndNeighbors.size()==0 || fCenterPointInterface.size()==0){
        PZError << "TPZGradientReconstruction::TPZDataGradient: gradient has size equal zero.\n";
        DebugStop();
    }
    
    int nneigh = fSolCellAndNeighbors.size()-1;
    if(nneigh==0) DebugStop();
    int i, j;
    
    //-------- getting min and max solution of neighbors ----
    STATE solKmax, solKmin;
    solKmin = fSolCellAndNeighbors[0];
    solKmax = fSolCellAndNeighbors[0];
    
    for(i = 1; i<= nneigh; i++)
    {
        if(solKmin > fSolCellAndNeighbors[i])
        {
            solKmin = fSolCellAndNeighbors[i];
            continue;
        }
        
        if(solKmax < fSolCellAndNeighbors[i])
        {
            solKmax = fSolCellAndNeighbors[i];
        }
    }
    
    if(solKmax < solKmin) DebugStop();
    
    // ----------- Calculating slope limiter --------------
    int ninterf = fCenterPointInterface.size();
    STATE solKside;
    STATE temp;
    TPZStack<STATE> alphavec;
    STATE solcel = fSolCellAndNeighbors[0];
    
    for (i = 0; i<ninterf; i++)
    {
        solKside = solcel;
        for(j=0; j<fdim; j++)
        {
            solKside += (STATE)(fCenterPointInterface[i][j] - fCenterPointCellAndNeighbors[0][j])*fGradient[j];
        }
        
        
        if(IsZero(solKside - solcel)) {
            temp = 1.;
        }
        else if((solKside - solcel) > 1.e-12)
        {
            temp = (solKmax - solcel)/(solKside-solcel);
            if(temp>1.) temp = 1.;
        }
        else if(solKside - solcel < 1.e-12)
        {
            temp = (solKmin-solcel)/(solKside-solcel);
            if(temp>1.) temp = 1.;
            if(temp<0.) temp = fabs(temp);
        }
        else temp =1.;
        
        if(temp < 0. || temp > 1.) DebugStop();
        alphavec.Push(temp);
    }
    
    //getting min slope limiter
    STATE alphaK = alphavec[0];
    for (long j=1; j<alphavec.size(); j++)
    {
        if(alphaK > alphavec[j])
        {
            alphaK = alphavec[j];
        }
    }
    fSlopeLimiter = alphaK;
    
    if(alphaK<0. || alphaK>1.) {
        DebugStop();
    }
}

//Venkatakrishnan V (1993). On the accuracy of limiters and convergence to steady state solutions
//Christopher Michalak, Carl Ollivier-Gooch (2009)-Accuracy preserving limiter for the high-order accurate solution of the Euler equations
void TPZGradientReconstruction::TPZGradientData::ComputeSlopeLimiter3()
{
    if(fGradient.size()==0 || fSolCellAndNeighbors.size()==0 || fCenterPointInterface.size()==0){
        PZError << "TPZGradientReconstruction::TPZDataGradient: gradient has size equal zero.\n";
        DebugStop();
    }
    
    int nneigh = fSolCellAndNeighbors.size()-1;
    if(nneigh==0) DebugStop();
    int i, j;
    
    //-------- getting min and max solution of neighbors ----
    STATE solKmax, solKmin;
    solKmin = fSolCellAndNeighbors[0];
    solKmax = fSolCellAndNeighbors[0];
    
    for(i = 1; i<= nneigh; i++)
    {
        if(solKmin > fSolCellAndNeighbors[i])
        {
            solKmin = fSolCellAndNeighbors[i];
            continue;
        }
        
        if(solKmax < fSolCellAndNeighbors[i])
        {
            solKmax = fSolCellAndNeighbors[i];
        }
    }
    
    if(solKmax < solKmin) DebugStop();
    
    // ----------- Calculating slope limiter --------------
    int ninterf = fCenterPointInterface.size();
    STATE solKside;
    STATE temp;
    TPZStack<STATE> alphavec;
    STATE solcel = fSolCellAndNeighbors[0];
    
    for (i = 0; i<ninterf; i++)
    {
        solKside = solcel;
        for(j=0; j<fdim; j++)
        {
            solKside += (STATE)(fCenterPointInterface[i][j] - fCenterPointCellAndNeighbors[0][j])*fGradient[j];
        }
        
        
        if(IsZero(solKside - solcel)) {
            temp = 1.;
        }
        else if((solKside - solcel) > 1.e-12)
        {
            temp = (solKmax - solcel)/(solKside-solcel);
            temp = (temp*temp + 2.*temp)/(temp*temp + temp + 2.);
            
        }
        else if(solKside - solcel < 1.e-12)
        {
            temp = (solKmin-solcel)/(solKside-solcel);
            temp = (temp*temp + 2.*temp)/(temp*temp + temp + 2.);
        }
        else temp =1.;
        
        if(temp < 0.) DebugStop();
        alphavec.Push(temp);
    }
    
    //getting min slope limiter
    STATE alphaK = alphavec[0];
    for (long j=1; j<alphavec.size(); j++)
    {
        if(alphaK > alphavec[j])
        {
            alphaK = alphavec[j];
        }
    }
    fSlopeLimiter = alphaK;
}

void TPZGradientReconstruction::TPZGradientData::ComputeWeights(REAL paramk)
{
    
    if(paramk<1 || paramk>2) DebugStop();
    
    TPZManVector<REAL,3> nodecelX;
    
    //Node more closer of the cel barycenter
    NodeCloserCenterX(nodecelX);
    
    long nneighs = fCelAndNeighbors.size()-1;
    TPZManVector<REAL,30> dist(nneighs,0.);
    TPZManVector<REAL,3> centerneigh(3,0.0);
    
    REAL sum=0.;
    for(int in = 0; in < nneighs; in++)
    {
        TPZGeoEl * gel = fCelAndNeighbors[in+1]->Reference();
        if(!gel) DebugStop();
        
        for(int k=0; k<fdim; k++) centerneigh[k]=fCenterPointCellAndNeighbors[in+1][k];
        
        dist[in]=gel->Distance(centerneigh, nodecelX);
        
        sum += 1./pow((REAL)dist[in],paramk);
    }
    
    fWeightsGrad.Resize(dist.size(), 0.);
    REAL temp;
    
    for (int id = 0; id<dist.size(); id++)
    {
        temp = 1./pow((REAL)dist[id],paramk);
        fWeightsGrad[id] = temp/sum;
        if(fWeightsGrad[id]<0 || fWeightsGrad[id]>1) DebugStop();
    }
}

void TPZGradientReconstruction::TPZGradientData::NodeCloserCenterX(TPZManVector<REAL,3> &nodecelX)
{
    TPZGeoEl* gel = fCelAndNeighbors[0]->Reference();
    if(!gel) DebugStop();
    
    TPZManVector<REAL,3> centercelX(3,0.);
    TPZManVector<REAL,3> coordX(3,0.);
    nodecelX.Resize(3, 0.);
    
    for(int k=0; k<fdim; k++) centercelX[k]=fCenterPointCellAndNeighbors[0][k];
    
    int nnodes = gel->NCornerNodes();
    TPZManVector<REAL,8> dist(nnodes,0.);
    
    
    for (int in = 0; in<nnodes; in++)
    {
        TPZGeoNode geonode = gel->Node(in);
        geonode.GetCoordinates(coordX);
        
        dist[in]=gel->Distance(centercelX, coordX);
    }
    
    int id = 0;
    for(int j=1; j<dist.size(); j++){
        
        if(dist[id] > dist[j]){
            id = j;
        }
    }
    
    gel->Node(id).GetCoordinates(nodecelX);
}

void TPZGradientReconstruction::TPZGradientData::InsertWeights(TPZFMatrix<REAL> &DeltaH, TPZFMatrix<REAL> &DifSol){
    
    if(fWeightsGrad[0]==0) DebugStop();
    if(DeltaH.Rows()!=fWeightsGrad.size()) DebugStop();
    if(DifSol.Rows()!=fWeightsGrad.size()) DebugStop();
    
    long ncH = DeltaH.Cols();
    long ncD = DifSol.Cols();
    for (long i = 0; i<fWeightsGrad.size(); i++) {
        
        for (long j = 0; j<ncH; j++){
            DeltaH(i,j) = fWeightsGrad[i]*DeltaH(i,j);
        }
        
        for (long k = 0; k<ncD; k++){
            DifSol(i,k) = fWeightsGrad[i]*DifSol(i,k);
        }
    }
}


