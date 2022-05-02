/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixLCC methods.
 */

#include "pzstrmatrixLCC.h"

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pzsubcmesh.h"
#include "pzanalysis.h"
#include "pzsfulmat.h"

#include "pzgnode.h"
#include "TPZTimer.h"


#include "pzcheckconsistency.h"
#include "TPZMaterial.h"

using namespace std;

#include "pzlog.h"

#include "run_stats_table.h"
#include <thread>

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixLCC");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif

RunStatsTable stat_ass_graph_LCC("-ass_graph_ot", "Run statistics table for the graph creation and coloring TPZStructMatrixLCC.");


TPZStructMatrixLCC::TPZStructMatrixLCC(TPZCompMesh *mesh) : TPZStructMatrixBase(mesh) {
    stat_ass_graph_LCC.start();
    TPZManVector<int64_t> ElementOrder;
    
    
    TPZStructMatrixLCC::OrderElement(this->Mesh(), ElementOrder);
    TPZVec<int64_t> elcolors;
    TPZStructMatrixLCC::ElementColoring(this->Mesh(), ElementOrder, fElSequenceColor, fElBlocked, elcolors);
    stat_ass_graph_LCC.stop();
    

}

TPZStructMatrixLCC::TPZStructMatrixLCC(TPZAutoPointer<TPZCompMesh> cmesh) : TPZStructMatrixBase(cmesh) {
    stat_ass_graph_LCC.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixLCC::OrderElement(this->Mesh(), ElementOrder);
    TPZVec<int64_t> elcolors;
    TPZStructMatrixLCC::ElementColoring(this->Mesh(), ElementOrder, fElSequenceColor, fElBlocked, elcolors);
    stat_ass_graph_LCC.stop();
    
    
}

TPZStructMatrixLCC::TPZStructMatrixLCC(const TPZStructMatrixLCC &copy) : TPZStructMatrixBase(copy) {
    fElSequenceColor = copy.fElSequenceColor;
    fElBlocked = copy.fElBlocked;
    fEquationFilter = copy.fEquationFilter;
    fElementsComputed = copy.fElementsComputed;
    fElementCompleted = copy.fElementCompleted;
    fSomeoneIsSleeping = copy.fSomeoneIsSleeping;
    
}


TPZMatrix<STATE> *TPZStructMatrixLCC::Create() {
    cout << "TPZStructMatrixLCC::Create should never be called\n";
    return 0;
}

TPZStructMatrixLCC *TPZStructMatrixLCC::Clone() {
    cout << "TPZStructMatrixLCC::Clone should never be called\n";
    return 0;
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

void TPZStructMatrixLCC::Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_stiff.start();
    if (fEquationFilter.IsActive()) {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<STATE> rhsloc(neqcondense,rhs.Cols(),0.);
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhsloc,guiInterface);
        }
        else{
            this->Serial_Assemble(stiffness,rhsloc,guiInterface);
        }
        
        fEquationFilter.Scatter(rhsloc, rhs);
    }
    else
    {
        if(this->fNumThreads){
            this->MultiThread_Assemble(stiffness,rhs,guiInterface);
        }
        else{
            this->Serial_Assemble(stiffness,rhs,guiInterface);
        }
    }
    ass_stiff.stop();
}

void TPZStructMatrixLCC::Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface){
    ass_rhs.start();
    if(fEquationFilter.IsActive())
    {
        int64_t neqcondense = fEquationFilter.NActiveEquations();
        int64_t neqexpand = fEquationFilter.NEqExpand();
        if(rhs.Rows() != neqexpand || Norm(rhs) != 0.)
        {
            DebugStop();
        }
        TPZFMatrix<STATE> rhsloc(neqcondense,1,0.);
        if(this->fNumThreads)
        {
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhsloc);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhsloc,guiInterface);
#ifdef HUGEDEBUG
            rhsserial -= rhsloc;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else
        {
            this->Serial_Assemble(rhsloc,guiInterface);
        }
        fEquationFilter.Scatter(rhsloc,rhs);
    }
    else
    {
        if(this->fNumThreads){
#ifdef HUGEDEBUG
            TPZFMatrix<STATE> rhsserial(rhs);
            this->Serial_Assemble(rhsserial, guiInterface);
#endif
            this->MultiThread_Assemble(rhs,guiInterface);
#ifdef HUGEDEBUG
            REAL normrhs = Norm(rhs);
            REAL normrhsserial = Norm(rhsserial);
            std::cout << "normrhs = " << normrhs << " normrhsserial " << normrhsserial << std::endl;
            rhsserial -= rhs;
            REAL norm = Norm(rhsserial);
            std::cout << "difference between serial and parallel " << norm << std::endl;
#endif
        }
        else{
            this->Serial_Assemble(rhs,guiInterface);
        }
    }
    ass_rhs.stop();
}



void TPZStructMatrixLCC::Serial_Assemble(TPZMatrix<STATE> & stiffness, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface ){
    
    if(!fMesh){
        LOGPZ_ERROR(logger,"Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef PZ_LOG
    if(dynamic_cast<TPZSubCompMesh * >(fMesh))
    {
        std::stringstream sout;
        sout << "AllEig = {};";
        LOGPZ_DEBUG(loggerelmat,sout.str())
        
    }
#endif
#ifdef PZDEBUG
    if (rhs.Rows() != fEquationFilter.NActiveEquations()) {
        DebugStop();
    }
#endif
    
    int64_t iel;
    int64_t nelem = fMesh->NElements();
    TPZElementMatrix ek(fMesh, TPZElementMatrix::EK),ef(fMesh, TPZElementMatrix::EF);
#ifdef PZ_LOG
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    int64_t count = 0;
    for(iel=0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        int matidsize = fMaterialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(fMaterialIds) == false) continue;
            }
            else
            {
                int matid = mat->Id();
                if (this->ShouldCompute(matid) == false) continue;
            }
        }
        
        count++;
        if(!(count%1000))
        {
            std::cout << '*';
            std::cout.flush();
        }
        if(!(count%20000))
        {
            std::cout << "\n";
        }
        calcstiff.start();
        
        el->CalcStiff(ek,ef);
        
        if(guiInterface) if(guiInterface->AmIKilled()){
            return;
        }
        
#ifdef PZ_LOG
        if(dynamic_cast<TPZSubCompMesh * >(fMesh))
        {
            std::stringstream objname;
            objname << "Element" << iel;
            std::string name = objname.str();
            objname << " = ";
            std::stringstream sout;
            ek.fMat.Print(objname.str().c_str(),sout,EMathematicaInput);
            sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";
            
            LOGPZ_DEBUG(loggerelmat,sout.str())
            /*		  if(iel == 133)
             {
             std::stringstream sout2;
             el->Reference()->Print(sout2);
             el->Print(sout2);
             LOGPZ_DEBUG(logger,sout2.str())
             }
             */
        }
        
#endif
        
#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if(!result)
        {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck,sout.str())
        }
        
#endif
        
        calcstiff.stop();
        assemble.start();
        
        if(!el->HasDependency()) {
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<STATE> test(stiffness);
            //			TPZFMatrix<STATE> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);
            
#ifdef PZ_LOG
            if(loggerel.isDebugEnabled())
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                if(gel)
                {
                    TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                    gel->CenterPoint(gel->NSides()-1, center);
                    gel->X(center, xcenter);
                    sout << "Stiffness for computational element index " << el->Index() << std::endl;
                    sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                }
                else {
                    sout << "Stiffness for computational element without associated geometric element\n";
                }
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            fEquationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
#ifdef PZ_LOG
            if(loggerel.isDebugEnabled() && ! dynamic_cast<TPZSubCompMesh *>(fMesh))
            {
                std::stringstream sout;
                TPZGeoEl *gel = el->Reference();
                TPZManVector<REAL> center(gel->Dimension()),xcenter(3,0.);
                gel->CenterPoint(gel->NSides()-1, center);
                gel->X(center, xcenter);
                sout << "Stiffness for geometric element " << gel->Index() << " center " << xcenter << std::endl;
                ek.Print(sout);
                ef.Print(sout);
                LOGPZ_DEBUG(loggerel,sout.str())
            }
#endif
        }
        
        assemble.stop();
    }//fim for iel
    if(count > 20) std::cout << std::endl;
    
#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        stiffness.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
    
#endif
    
}

void TPZStructMatrixLCC::Serial_Assemble(TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface){
    
    int64_t iel;
//    int64_t nelem = fMesh->NElements();
    
    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");
    
    TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();
    
    stat_ass_graph_LCC.start();
    TPZManVector<int64_t> ElementOrder;
    TPZStructMatrixLCC::OrderElement(this->Mesh(), ElementOrder);
    TPZVec<int64_t> elcolors;
    TPZStructMatrixLCC::ElementColoring(this->Mesh(), ElementOrder, fElSequenceColor, fElBlocked, elcolors);
    stat_ass_graph_LCC.stop();

    int64_t elseqsize = fElSequenceColor.size();
    for (int64_t index = 0; index < elseqsize; index++) {
        iel = fElSequenceColor[index];
        TPZCompEl *el = elementvec[iel];
        if(!el) continue;
        
        TPZMaterial * mat = el->Material();
        if (!mat) continue;
        int matid = mat->Id();
        if (this->ShouldCompute(matid) == false) continue;
        
        TPZElementMatrix ef(fMesh, TPZElementMatrix::EF);
        
        calcresidual.start();
        
        el->CalcResidual(ef);
        
        calcresidual.stop();
        
        assemble.start();
        
        if(!el->HasDependency()) {
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            fEquationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
        }
        
        assemble.stop();
        
    }//fim for iel
#ifdef PZ_LOG
    {
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
    }
#endif
    //std::cout << std::endl;
}

TPZMatrix<STATE> * TPZStructMatrixLCC::CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    TPZMatrix<STATE> *stiff = Create();
    
    //int64_t neq = stiff->Rows();
    int64_t cols = MAX(1, rhs.Cols());
    rhs.Redim(fEquationFilter.NEqExpand(),cols);
    Assemble(*stiff,rhs,guiInterface);
    
#ifdef PZ_LOG2
    if(loggerel.isDebugEnabled())
    {
        std::stringstream sout;
        stiff->Print("Stiffness matrix",sout);
        rhs.Print("Right hand side", sout);
        LOGPZ_DEBUG(loggerel,sout.str())
    }
#endif
    return stiff;
    
}


void TPZStructMatrixLCC::MultiThread_Assemble(TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs, TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    
    const int numthreads = this->fNumThreads;
    std::cout << "Assemble numthreads = " << numthreads << std::endl;
    std::vector<std::thread> allthreads(numthreads);
    int itr;
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return;
        }
    }

    this->fCurrentIndex = 0;
    
    fElementCompleted = -1;
    fElementsComputed.Resize(fMesh->NElements());
    fElementsComputed.Fill(0);
    fSomeoneIsSleeping = 0;
    TPZManVector<ThreadData*> allthreaddata(numthreads);
#ifdef PZDEBUG
    {
        for (int64_t i=1; i<fElBlocked.size(); i++) {
            if (fElBlocked[i] < fElBlocked[i-1]) {
                std::cout << "i = " << i << " fElBlocked[i-1] " << fElBlocked[i-1] << " fElBlocked[i] " << fElBlocked[i] << std::endl;
            }
        }
    }
#endif

    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] = new ThreadData(this, itr, mat, rhs, fMaterialIds, guiInterface);
        ThreadData &threaddata = *allthreaddata[itr];
        
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;
        threaddata.fSomeoneIsSleeping = &fSomeoneIsSleeping;
        allthreads.emplace_back(ThreadData::ThreadWork, &threaddata);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreads[itr].join();
    }

    for(itr=0; itr<numthreads; itr++)
    {
        delete allthreaddata[itr];
    }

#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matriz de Rigidez: ",sout,EMathematicaInput);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
}


void TPZStructMatrixLCC::MultiThread_Assemble(TPZFMatrix<STATE> & rhs,TPZAutoPointer<TPZGuiInterface> guiInterface)
{
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads(numthreads);
    int itr;
    if(guiInterface){
        if(guiInterface->AmIKilled()){
            return;
        }
    }
    

    this->fCurrentIndex = 0;
    fElementCompleted = -1;
    fElementsComputed.Resize(fMesh->NElements());
    fElementsComputed.Fill(0);
    fSomeoneIsSleeping = 0;
    TPZManVector<ThreadData*> allthreaddata(numthreads);

    for(itr=0; itr<numthreads; itr++)
    {
        allthreaddata[itr] = new ThreadData(this, itr, rhs, fMaterialIds, guiInterface);
        ThreadData &threaddata = *allthreaddata[itr];
        threaddata.fElBlocked=&fElBlocked;
        threaddata.fElSequenceColor=&fElSequenceColor;
        threaddata.fElementCompleted = &fElementCompleted;
        threaddata.fComputedElements = &fElementsComputed;
        threaddata.fSomeoneIsSleeping = &fSomeoneIsSleeping;

        allthreads.emplace_back(ThreadData::ThreadWork, &threaddata);
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        allthreads[itr].join();
    }
    
    for(itr=0; itr<numthreads; itr++)
    {
        delete allthreaddata[itr];
    }


#ifdef PZ_LOG
    if(loggerCheck.isDebugEnabled())
    {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        rhs.Print("Right Handside", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck,sout.str())
    }
#endif
    
}



TPZStructMatrixLCC::ThreadData::ThreadData(TPZStructMatrixLCC *strmat, int seqnum, TPZMatrix<STATE> &mat,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat), fGuiInterface(guiInterface), fGlobMatrix(&mat), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
    fCurrentIndex = &strmat->fCurrentIndex;
    
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

TPZStructMatrixLCC::ThreadData::ThreadData(TPZStructMatrixLCC *strmat, int seqnum,
                                          TPZFMatrix<STATE> &rhs,
                                          std::set<int> &MaterialIds,
                                          TPZAutoPointer<TPZGuiInterface> guiInterface)
: fStruct(strmat),fGuiInterface(guiInterface), fGlobMatrix(0), fGlobRhs(&rhs), fThreadSeqNum(seqnum)
{
    this->fCurrentIndex = &strmat->fCurrentIndex;
    /*	sem_t *sem_open( ... );
     int sem_close(sem_t *sem);
     int sem_unlink(const char *name);
     */
    /*
     #ifdef MACOSX
     std::stringstream sout;
     static int counter = 0;
     sout << "AssemblySem" << counter++;
     fAssembly = sem_open(sout.str().c_str(), O_CREAT,777,1);
     if(fAssembly == SEM_FAILED)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     DebugStop();
     }
     #else
     int sem_result = sem_init(&fAssembly,0,0);
     if(sem_result != 0)
     {
     std::cout << __PRETTY_FUNCTION__ << " could not open the semaphore\n";
     }
     #endif
     */
}

//#define DRY_RUN

void *TPZStructMatrixLCC::ThreadData::ThreadWork(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
//    TPZStructMatrixLCC *strmat = data->fStruct;
    TPZVec<int64_t> &ComputedElements = *(data->fComputedElements);
    TPZVec<int64_t> &ElBlocked = *(data->fElBlocked);

    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ek(cmesh,TPZElementMatrix::EK);
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    int64_t numelements = data->fElSequenceColor->size();
    int64_t index = data->fCurrentIndex->fetch_add(1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef HUGEDEBUG
    data->fMutexAccessElement.lock();
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    data->fMutexAccessElement.unlock();
#endif
    
    while (index < numelements)
    {
        
        int64_t iel = data->fElSequenceColor->operator[](index);
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Computing element " << index;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
#ifdef PZ_LOG
        std::stringstream sout;
        sout << "Element " << index << " elapsed time ";
        TPZTimer timeforel(sout.str());
        timeforel.start();
#endif

        if (iel >= 0){
            TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
            el->CalcStiff(ek,ef);
#else
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
            intel->InitializeElementMatrix(ek,ef);
#endif

            if(guiInterface) if(guiInterface->AmIKilled()){
                break;
            }
            
#ifndef DRY_RUN
            if(!el->HasDependency()) {
                
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
                
#ifdef PZ_LOG
                if(loggerel.isDebugEnabled())
                {
                    std::stringstream sout;
                    ek.fMat.Print("Element stiffness matrix",sout);
                    ef.fMat.Print("Element right hand side", sout);
                    LOGPZ_DEBUG(loggerel,sout.str())
                }
#endif
            } else {
                // the element has dependent nodes
                ek.ApplyConstraints();
                ef.ApplyConstraints();
                ek.ComputeDestinationIndices();
                data->fStruct->FilterEquations(ek.fSourceIndex,ek.fDestinationIndex);
            }
#endif
            
#ifdef PZ_LOG
            timeforel.stop();
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << timeforel.processName() <<  timeforel;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif

            int64_t localcompleted = *(data->fElementCompleted);
            bool localupdated = false;
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            int64_t needscomputed = ElBlocked[index];
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                if (localupdated) {
                    sout << "Localcompleted updated without thread lock\n";
                }
                
                sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            
#ifdef HUGEDEBUG
            data->fMutexAccessElement.lock();
            std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
            data->fMutexAccessElement.unlock();
#endif
            
            bool hadtowait = false;
            while (needscomputed > localcompleted) {
                // block the thread till the element needed has been assembled
              std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
                SomeoneIsSleeping = 1;
                hadtowait = true;
#ifdef HUGEDEBUG
                std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                std::cout.flush();
#endif
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " cannot be assembled - going to sleep";
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                data->fConditionVar.wait(lock);
                
                localcompleted = *data->fElementCompleted;
                localupdated = false;
                while (ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "thread wakeup for element index " << index << std::endl;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
            
#ifdef HUGEDEBUG
            if (hadtowait) {
              std::scoped_lock lock(data->fMutexAccessElement);
                std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
            }
#endif
            
#ifndef DRY_RUN
            auto globMatrix =
                    dynamic_cast<TPZMatrix<STATE> *>(data->fGlobMatrix);
                auto globRhs =
                    dynamic_cast<TPZFMatrix<STATE> *>(data->fGlobRhs);
            if(globMatrix){
                // Assemble the matrix
                if(!ek.HasDependency())
                {
                    globMatrix->AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                    globRhs->AddFel(ef.fMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                else
                {
                    globMatrix->AddKel(ek.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                    globRhs->AddFel(ef.fConstrMat,ek.fSourceIndex,ek.fDestinationIndex);
                }
                
            }
#endif
            
            localupdated = false;
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            if (localcompleted == index-1) {
                localcompleted++;
            }
            while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                localcompleted++;
                localupdated = true;
            }
            bool elementcompletedupdate = false;
            if (*data->fElementCompleted < localcompleted) {
//                std::cout << "Updating element completed " << localcompleted << std::endl;
                *data->fElementCompleted = localcompleted;
                elementcompletedupdate = true;
            }
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Element " << index << " has been assembled ";
                if (localupdated) {
                    sout << "\nLocalcompleted updated without thread lock\n";
                }
                if (elementcompletedupdate) {
                    sout << "\nfElementCompleted updated to localcompleted\n";
                }
                sout << "local completed " << localcompleted;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            ComputedElements[index] = 1;
            if (SomeoneIsSleeping) {
              std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
#ifdef HUGEDEBUG
                std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                std::cout.flush();
#endif
                
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    LOGPZ_DEBUG(logger, "condition broadcast")
                }
#endif
                SomeoneIsSleeping = 0;
                data->fConditionVar.notify_all();
            }

        }
        else
        {
            std::cout << "the element in ElColorSequence is negative???\n";
            DebugStop();
        }
        index = data->fCurrentIndex->fetch_add(1);
        
    }
    // just make sure threads that were accidentally blocked get woken up
    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
    bool localupdated = false;
    int64_t localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    data->fConditionVar.notify_all();
    SomeoneIsSleeping = 0;
    return 0;
}

void *TPZStructMatrixLCC::ThreadData::ThreadWorkResidual(void *datavoid)
{
    ThreadData *data = (ThreadData *) datavoid;
    //    TPZStructMatrixLCC *strmat = data->fStruct;
    TPZVec<int64_t> &ComputedElements = *(data->fComputedElements);
    TPZVec<int64_t> &ElBlocked = *(data->fElBlocked);
    
    int &SomeoneIsSleeping = *(data->fSomeoneIsSleeping);
    
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    TPZAutoPointer<TPZGuiInterface> guiInterface = data->fGuiInterface;
    TPZElementMatrix ef(cmesh,TPZElementMatrix::EF);
    int64_t numelements = data->fElSequenceColor->size();
    int64_t index = data->fCurrentIndex->fetch_add(1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
        sout << "index = " << index << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef HUGEDEBUG
    data->fMutexAccessElement.lock();
    std::cout << "ThreadData starting with " << data->fThreadSeqNum << " total elements " << numelements << std::endl;
    std::cout << "index = " << index << std::endl;
    std::cout.flush();
    data->fMutexAccessElement.unlock();
#endif
    
        while (index < numelements)
        {
            
            int64_t iel = data->fElSequenceColor->operator[](index);
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Computing element " << index;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
#ifdef PZ_LOG
            std::stringstream sout;
            sout << "Element " << index << " elapsed time ";
            TPZTimer timeforel(sout.str());
            timeforel.start();
#endif
            
            if (iel >= 0){
                TPZCompEl *el = cmesh->ElementVec()[iel];
#ifndef DRY_RUN
                el->CalcResidual(ef);
#else
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(el);
                intel->InitializeElementMatrix(ef);
#endif
                if(guiInterface) if(guiInterface->AmIKilled()){
                    break;
                }
                
#ifndef DRY_RUN
                if(!el->HasDependency()) {
                    
                    ef.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                    
#ifdef PZ_LOG
                    if(loggerel.isDebugEnabled())
                    {
                        std::stringstream sout;
                        ef.fMat.Print("Element right hand side", sout);
                        LOGPZ_DEBUG(loggerel,sout.str())
                    }
#endif
                } else {
                    // the element has dependent nodes
                    ef.ApplyConstraints();
                    ef.ComputeDestinationIndices();
                    data->fStruct->FilterEquations(ef.fSourceIndex,ef.fDestinationIndex);
                }
#endif
                
#ifdef PZ_LOG
                timeforel.stop();
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << timeforel.processName() <<  timeforel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
                int64_t localcompleted = *(data->fElementCompleted);
                bool localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                int64_t needscomputed = ElBlocked[index];
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    if (localupdated) {
                        sout << "Localcompleted updated without thread lock\n";
                    }
                    
                    sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                
#ifdef HUGEDEBUG
                data->fMutexAccessElement.lock();
                std::cout << "threadEK " << data->fThreadSeqNum << " index " << index << " localcompleted " << localcompleted << " needscomputed " << needscomputed << std::endl;
                data->fMutexAccessElement.unlock();
#endif
                
                bool hadtowait = false;
#ifdef PZ_LOG
                if (logger.isInfoEnabled() && needscomputed > localcompleted)
                {
                    std::stringstream sout;
                    sout << "Element " << index << " cannot be assembled - going to sleep";
                    LOGPZ_INFO(logger, sout.str())
                }
#endif
                while (needscomputed > localcompleted) {
                    // block the thread till the element needed has been assembled
                    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
                    SomeoneIsSleeping = 1;
                    hadtowait = true;
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " going to sleep waiting for " << needscomputed << std::endl;
                    std::cout.flush();
#endif
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Element " << index << " cannot be assembled - going to sleep";
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    data->fConditionVar.wait(lock);
                    
                    localcompleted = *data->fElementCompleted;
                    localupdated = false;
                    while (ComputedElements[localcompleted+1] == 1) {
                        localcompleted++;
                        localupdated = true;
                    }
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "thread wakeup for element index " << index << std::endl;
                        if (localupdated) {
                            sout << "Localcompleted updated without thread lock\n";
                        }
                        
                        sout << "Element " << index << " is computed and can assemble if required " << needscomputed << " is smaller than localcompleted " << localcompleted;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                

#ifdef PZ_LOG
                if (logger.isInfoEnabled() && hadtowait)
                {
                    std::stringstream sout;
                    sout << "thread wakeup for element index " << index << std::endl;
                    LOGPZ_INFO(logger, sout.str())
                }
#endif


#ifdef HUGEDEBUG
                if (hadtowait) {
                  std::scoped_lock lock (data->fMutexAccessElement);
                    std::cout << "threadEK " <<data->fThreadSeqNum << " Index " << index << " continuing\n";
                }
#endif
                
#ifndef DRY_RUN
                auto globRhs =
                    dynamic_cast<TPZFMatrix<STATE> *>(data->fGlobRhs);
                if(globRhs){
                    // Assemble the matrix
                    if(!ef.HasDependency())
                    {
                        globRhs->AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);
                    }
                    else
                    {
                        globRhs->AddFel(ef.fConstrMat,ef.fSourceIndex,ef.fDestinationIndex);
                    }
                    
                }
#endif
                
                localupdated = false;
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                if (localcompleted == index-1) {
                    localcompleted++;
                }
                while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
                    localcompleted++;
                    localupdated = true;
                }
                bool elementcompletedupdate = false;
                if (*data->fElementCompleted < localcompleted) {
                    //                std::cout << "Updating element completed " << localcompleted << std::endl;
                    *data->fElementCompleted = localcompleted;
                    elementcompletedupdate = true;
                }
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Element " << index << " has been assembled ";
                    if (localupdated) {
                        sout << "\nLocalcompleted updated without thread lock\n";
                    }
                    if (elementcompletedupdate) {
                        sout << "\nfElementCompleted updated to localcompleted\n";
                    }
                    sout << "local completed " << localcompleted;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                ComputedElements[index] = 1;
                if (SomeoneIsSleeping) {
                  std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
#ifdef HUGEDEBUG
                    std::cout << "threadEK " <<data->fThreadSeqNum <<  " Computed index " << index << " Waking up ElementsCompleted " << *data->fElementCompleted << std::endl;
                    std::cout.flush();
#endif
                    
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        LOGPZ_DEBUG(logger, "condition broadcast")
                    }
#endif
                    SomeoneIsSleeping = 0;
                    data->fConditionVar.notify_all();
                }
                
            }
            else
            {
                std::cout << "the element in ElColorSequence is negative???\n";
                DebugStop();
            }
            index = data->fCurrentIndex->fetch_add(1);
            
        }
    // just make sure threads that were accidentally blocked get woken up
    std::unique_lock<std::mutex> lock(data->fMutexAccessElement);
    bool localupdated = false;
    int64_t localcompleted = *(data->fElementCompleted);
    while (localcompleted < numelements-1 && ComputedElements[localcompleted+1] == 1) {
        localcompleted++;
        localupdated = true;
    }
    if (localcompleted == index-1) {
        localcompleted++;
    }
    bool elementcompletedupdate = false;
    if (*data->fElementCompleted < localcompleted) {
        //                std::cout << "Updating element completed " << localcompleted << std::endl;
        *data->fElementCompleted = localcompleted;
        elementcompletedupdate = true;
    }
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        LOGPZ_DEBUG(logger, "Finishing up")
        if (localupdated) {
            LOGPZ_DEBUG(logger, "updated localcompleted")
        }
        if (elementcompletedupdate) {
            LOGPZ_DEBUG(logger, "updated fElementCompleted")
        }
        if (localupdated || elementcompletedupdate) {
            std::stringstream sout;
            sout << "localcompleted " << localcompleted;
            LOGPZ_DEBUG(logger, sout.str())
        }
        LOGPZ_DEBUG(logger, "finishing and condition broadcast")
    }
#endif
    data->fConditionVar.notify_all();
    SomeoneIsSleeping = 0;
    return 0;
}

//static bool CanAssemble(TPZStack<int64_t> &connectlist, TPZVec<int> &elContribute)
//{
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        if (elContribute[connectlist[i]] >= 0){
//            return false;
//        }
//    }
//    return true;
//}

static void AssembleColor(int64_t el, TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute) {
    for (auto connect : connectlist) {
        elContribute[connect] = el;
    }
}

static int64_t WhoBlockedMe(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &elSequenceColorInv) {
    int64_t el = -1;
    for (auto connect : connectlist) {
        auto elBlocked = elContribute[connect];
        if (elBlocked == -1) continue;
        auto elBlockedIndex = elSequenceColorInv[elBlocked];
        if (el == -1) el = elBlockedIndex;
        if (elBlockedIndex > el) el = elBlockedIndex;
    }
    return el;
}

//static void RemoveEl(int el,TPZCompMesh *cmesh,TPZVec<int> &elContribute,int elSequence)
//{
//    TPZCompEl *cel = cmesh->ElementVec()[el];
//    if(!cel) DebugStop();
//    TPZStack<int64_t> connectlist;
//    cel->BuildConnectList(connectlist);
//    for (int i = 0 ; i < connectlist.NElements() ; i++)
//    {
//        int conindex = connectlist[i];
//        if (elContribute[conindex] != elSequence){
//            DebugStop();
//        }
//        elContribute[conindex] = -1;
//    }
//}

static int64_t MinPassIndex(TPZStack<int64_t> &connectlist, TPZVec<int64_t> &elContribute, TPZVec<int64_t> &passIndex) {
    int64_t minPassIndex = -1;
    for (int64_t i = 0; i < connectlist.NElements(); i++) {
        int64_t elcont = elContribute[connectlist[i]];
        int64_t passindex = -1;
        if (elcont != -1) {
            passindex = passIndex[elcont];
            if (minPassIndex == -1) minPassIndex = passindex;
        }
        if (minPassIndex < passindex) minPassIndex = passindex;
    }
    return minPassIndex;
}

//Coloring computational elements 2022
void TPZStructMatrixLCC::orderElements(){
    int64_t nconnects = fMesh->NConnects();
    int64_t nelement = fMesh->ElementVec().NElements();
    fElVecColor.Resize(nelement);
    fElVecColor.Fill(-1);
    TPZVec<int> used(nconnects,0);
    bool haswork = true;
    int color = 0;
    
    while (haswork){
        haswork = false;
        for (int64_t iel = 0; iel < nelement; iel++){
            TPZCompEl *compel = fMesh->Element(iel);
            if (!compel) continue;
            if (fElVecColor[iel] != -1) continue;
            
            TPZStack<int64_t> connectlist;
            compel->BuildConnectList(connectlist);
            
            bool canColor = true;
            for (auto ic : connectlist){
                if (used[ic]){
                    canColor = false;
                    haswork = true;
                }
            }
            if (canColor){
                for (auto ic : connectlist) used[ic] = 1;
                fElVecColor[iel] = color;
            }
        }
        if (haswork) {
            color++;
            used.Fill(0);
        }
    }
}

/**
 * elSequence (input) element sequence according to the connect sequence numbers
 * elSequenceColor (output) the colored element sequence
 * elBlocked the element index which needs to have been computed before assembling the element
 * elColors (output) number of elements in each color
 */ 
void TPZStructMatrixLCC::ElementColoring(TPZCompMesh *cmesh, TPZVec<int64_t> &elSequence, TPZVec<int64_t> &elSequenceColor,
        TPZVec<int64_t> &elBlocked, TPZVec<int64_t> &NumelColors) {

    const int64_t nnodes = cmesh->NConnects();
    const int64_t nel = cmesh->NElements();

    if (nel == 0) return;

    NumelColors.Resize(nel, -1);

    // elContribute contains the element index which last contributed to the node
    // passIndex contains the color of the element
    TPZManVector<int64_t> elContribute(nnodes, -1), passIndex(nel, -1), elSequenceColorInv(nel, -1);
    int64_t elsequencesize = elSequence.size();
    elSequenceColor.Resize(elsequencesize);
    elSequenceColor.Fill(-1);
    elBlocked.Resize(elsequencesize);
    elBlocked.Fill(-1);
    int64_t nelProcessed = 0;
    int64_t currentPassIndex = 0;
    while (nelProcessed < elSequence.NElements()) {
        for (auto elindex : elSequence) {
            // if this element hasn't been computed in a previous pass
            if (elSequenceColorInv[elindex] == -1) {
                TPZCompEl *cel = cmesh->Element(elindex);
                if (!cel) continue;
                TPZStack<int64_t> connectlist;
                cel->BuildConnectList(connectlist);
                // compute the lowest color (pass index) of the elements that have contributed to this set of nodes
                int64_t minPass = MinPassIndex(connectlist, elContribute, passIndex);
                // no element has ever seen any of these nodes
                if (minPass == -1) {
                    passIndex[elindex] = currentPassIndex;
                    // the element index is put into elContribute (from there it is possible to know the color as well)
                    AssembleColor(elindex, connectlist, elContribute);
                    // initialize the data structures
                    elSequenceColor[nelProcessed] = elindex;
                    elSequenceColorInv[elindex] = nelProcessed;
                    nelProcessed++;
                } // this element cannot be computed as it is connected to another element of the current color
                else if (minPass == currentPassIndex) {
                } // the elements connected to this node are from a previous color
                else if (minPass < currentPassIndex) {
                    // the element with largest index which contributes to the set of nodes
                    // el is given in the new sequence order
                    const int64_t el = WhoBlockedMe(connectlist, elContribute, elSequenceColorInv);
                    // elblocked means the future element the element el will block
                    if (elBlocked[nelProcessed] != -1) DebugStop();
                    elBlocked[nelProcessed] = el;
                    passIndex[elindex] = currentPassIndex;
                    AssembleColor(elindex, connectlist, elContribute);
                    elSequenceColor[nelProcessed] = elindex;
                    elSequenceColorInv[elindex] = nelProcessed;
                    nelProcessed++;
                } else {
                    DebugStop();
                }
            }
        }
        NumelColors[currentPassIndex] = nelProcessed;
        currentPassIndex++;
    }

    NumelColors[currentPassIndex] = NumelColors[currentPassIndex - 1] + 1;
    NumelColors.Resize(currentPassIndex + 1);

#ifdef PZDEBUG
    std::ofstream toto("../ColorMeshDebug.txt");
    toto << "elSequence\n" << elSequence << std::endl;
    toto << "elSequenceColor\n" << elSequenceColor << std::endl;
    toto << "elSequenceColorInv\n" << elSequenceColorInv << std::endl;
    toto << "elBlocked\n" << elBlocked << std::endl;
    toto << "elContribute\n" << elContribute << std::endl;
    toto << "passIndex\n" << passIndex << std::endl;
    toto << "NumelColors\n" << NumelColors << std::endl;
    toto.close();
#endif
}

void TPZStructMatrixLCC::OrderElement(TPZCompMesh *cmesh, TPZVec<int64_t> &ElementOrder) {
    int64_t numelconnected = 0;
    int64_t nconnect = cmesh->NConnects();
    //firstelconnect contains the first element index in the elconnect vector
    TPZVec<int64_t> firstelconnect(nconnect + 1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        numelconnected += cmesh->ConnectVec()[ic].NElConnected();
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "numelconnected " << numelconnected << endl;
    //cout << "firstelconnect ";
    //  for(ic=0; ic<nconnect; ic++) cout << firstelconnect[ic] << ' ';
    TPZVec<int64_t> elconnect(numelconnected, -1);
    int64_t el;
    TPZCompEl *cel;

#ifdef NOORDER
    int64_t count = 0;
    std::cout << __PRETTY_FUNCTION__ << " no element order\n";
    ElementOrder.Resize(cmesh->NElements(), -1);
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        ElementOrder[count] = el;
        count++;
    }
    ElementOrder.Resize(count);
    return;
#endif
    for (el = 0; el < cmesh->ElementVec().NElements(); el++) {
        cel = cmesh->ElementVec()[el];
        if (!cel) continue;
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t nc = connectlist.NElements();
        for (int64_t ic = 0; ic < nc; ic++) {
            int64_t cindex = connectlist[ic];
            elconnect[firstelconnect[cindex]] = el;
            firstelconnect[cindex]++;
        }
    }
    //  for(ic=0; ic<numelconnected; ic++) cout << elconnect[ic] << endl;
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        firstelconnect[ic + 1] = firstelconnect[ic] + cmesh->ConnectVec()[ic].NElConnected();
    }
    //cout << "elconnect\n";
    //  int no;
    //  for(no=0; no< fMesh->ConnectVec().NElements(); no++) {
    //cout << "no numero " << no << ' ' << " seq num " << fMesh->ConnectVec()[no].SequenceNumber() << ' ';
    //       for(ic=firstelconnect[no]; ic<firstelconnect[no+1];ic++) cout << elconnect[ic] << ' ';
    //cout << endl;
    //  }

    ElementOrder.Resize(cmesh->ElementVec().NElements(), -1);
    ElementOrder.Fill(-1);
    TPZVec<int64_t> nodeorder(cmesh->ConnectVec().NElements(), -1);
    firstelconnect[0] = 0;
    for (int64_t ic = 0; ic < nconnect; ic++) {
        int64_t seqnum = cmesh->ConnectVec()[ic].SequenceNumber();
        if (seqnum >= 0) nodeorder[seqnum] = ic;
    }
    //  cout << "nodeorder ";
    /*  for(ic=0; ic<fMesh->ConnectVec().NElements(); ic++) cout << nodeorder[ic] << ' ';
     cout << endl;
     cout.flush();*/
    int64_t elsequence = 0;
    TPZVec<int> elorderinv(cmesh->ElementVec().NElements(), -1);
    for (int64_t seq = 0; seq < nconnect; seq++) {
        int64_t ic = nodeorder[seq];
        if (ic == -1) continue;
        int64_t firstind = firstelconnect[ic];
        int64_t lastind = firstelconnect[ic + 1];
        for (int64_t ind = firstind; ind < lastind; ind++) {
            el = elconnect[ind];
            if (el == -1) {
                continue;
            }
            if (elorderinv[el] == -1) elorderinv[el] = elsequence++;
        }
    }
    //  cout << "elorderinv ";
    //  for(seq=0;seq<fMesh->ElementVec().NElements();seq++) cout << elorderinv[seq] << ' ';
    //  cout << endl;
    elsequence = 0;
    for (int64_t seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (elorderinv[seq] == -1) continue;
        ElementOrder[elorderinv[seq]] = seq;
    }
    int64_t seq;
    for (seq = 0; seq < cmesh->ElementVec().NElements(); seq++) {
        if (ElementOrder[seq] == -1) break;
    }

    ElementOrder.Resize(seq);
}

int TPZStructMatrixLCC::ClassId() const{
    return Hash("TPZStructMatrixLCC") ^ TPZStructMatrixBase::ClassId() << 1;
}


void TPZStructMatrixLCC::Read(TPZStream& buf, void* context) {
    TPZStructMatrixBase::Read(buf, context);
    
    buf.Read(fElBlocked);
    buf.Read(fElSequenceColor);
    buf.Read(&fElementCompleted);
    buf.Read(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    buf.Read(&fCurrentIndexLong);
    fCurrentIndex = fCurrentIndexLong;
}

void TPZStructMatrixLCC::Write(TPZStream& buf, int withclassid) const {
    TPZStructMatrixBase::Write(buf, withclassid);
    
    buf.Write(fElBlocked);
    buf.Write(fElSequenceColor);
    buf.Write(&fElementCompleted);
    buf.Write(&fSomeoneIsSleeping);
    int64_t fCurrentIndexLong;
    fCurrentIndexLong = fCurrentIndex;
    buf.Write(&fCurrentIndexLong);
}

template class TPZRestoreClass<TPZStructMatrixLCC>;
