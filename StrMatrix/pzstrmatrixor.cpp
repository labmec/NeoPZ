/**
 * @file
 * @brief Contains the implementation of the TPZStructMatrixOR methods.
 */
#include "pzstrmatrixor.h"
#include "TPZStructMatrix.h"
#include "pzsubcmesh.h"
#include "TPZElementMatrixT.h"
#include "TPZTimer.h"
#include "run_stats_table.h"
#include "fpo_exceptions.h"
#include "TPZThreadPool.h"
#include <thread>
#include <vector>
using namespace std;

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix.TPZStructMatrixOR");
static TPZLogger loggerel("pz.strmatrix.element");
static TPZLogger loggerel2("pz.strmatrix.elementinterface");
static TPZLogger loggerelmat("pz.strmatrix.elementmat");
static TPZLogger loggerCheck("pz.strmatrix.checkconsistency");
static TPZLogger loggerGlobStiff("pz.strmatrix.globalstiffness");
#endif

#ifdef CHECKCONSISTENCY
static TPZCheckConsistency stiffconsist("ElementStiff");
#endif


template<class TVar>
TPZStructMatrixOR<TVar>::TPZStructMatrixOR() : TPZStrMatParInterface(){
    this->SetNumThreads(TPZThreadPool::globalInstance().threadCount());
}

static RunStatsTable ass_stiff("-ass_stiff", "Assemble Stiffness");
static RunStatsTable ass_rhs("-ass_rhs", "Assemble Stiffness");

template<class TVar>
void 
TPZStructMatrixOR<TVar>::Assemble(TPZBaseMatrix & stiffness, TPZBaseMatrix & rhs) {
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    ass_stiff.start();
    if (equationFilter.IsActive()) {
        int64_t neqcondense = equationFilter.NActiveEquations();
#ifdef PZDEBUG
        if (stiffness.Rows() != neqcondense) {
            DebugStop();
        }
#endif
        TPZFMatrix<TVar> rhsloc;
        if(ComputeRhs()) rhsloc.Redim(neqcondense, rhs.Cols());
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhsloc);
        } else {
            this->Serial_Assemble(stiffness, rhsloc);
        }

        if(ComputeRhs()) equationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(stiffness, rhs);
        } else {
            this->Serial_Assemble(stiffness, rhs);
        }
    }
    ass_stiff.stop();
}

template<class TVar>
void
TPZStructMatrixOR<TVar>::Assemble(TPZBaseMatrix & rhs) {
    const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
    ass_rhs.start();
    if (equationFilter.IsActive()) {
        auto rhsState = dynamic_cast<const TPZFMatrix<TVar> *>(&rhs);
            
        int64_t neqcondense = equationFilter.NActiveEquations();
        int64_t neqexpand = equationFilter.NEqExpand();
        if (rhs.Rows() != neqexpand || Norm(*rhsState) != 0.) {
            DebugStop();
        }
        TPZFMatrix<TVar> rhsloc(neqcondense, 1, 0.);
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhsloc);
        } else {
            this->Serial_Assemble(rhsloc);
        }
        equationFilter.Scatter(rhsloc, rhs);
    } else {
        if (this->fNumThreads) {
            this->MultiThread_Assemble(rhs);
        } else {
            this->Serial_Assemble(rhs);
        }
    }
    ass_rhs.stop();
}

template<class TVar>
void 
TPZStructMatrixOR<TVar>::Serial_Assemble(TPZBaseMatrix & stiff_base, TPZBaseMatrix & rhs_base) {
    if(
        !dynamic_cast<TPZMatrix<TVar>*>(&stiff_base)||
        !dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)
       ){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" Incompatible type. Aborting...\n";
        DebugStop();
    }
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    auto *cmesh = myself->Mesh();
    auto &materialIds = myself->MaterialIds();
    TPZMatrix<TVar> &stiffness = dynamic_cast<TPZMatrix<TVar>&>(stiff_base);
    TPZFMatrix<TVar> &rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
    if (!cmesh) {
        LOGPZ_ERROR(logger, "Serial_Assemble called without mesh")
        DebugStop();
    }
#ifdef PZ_LOG
    if (loggerelmat.isDebugEnabled()) {
        if (dynamic_cast<TPZSubCompMesh *> (cmesh)) {
            std::stringstream sout;
            sout << "AllEig = {};";
            LOGPZ_DEBUG(loggerelmat, sout.str())
        }
    }
#endif

#ifdef PZDEBUG
    if (ComputeRhs() && myself->EquationFilter().IsActive() &&
        rhs.Rows() != myself->EquationFilter().NActiveEquations()) {
        DebugStop();
    }
#endif

    int64_t iel;
    int64_t nelem = cmesh->NElements();
    TPZElementMatrixT<TVar> ek(cmesh, TPZElementMatrix::EK), ef(cmesh, TPZElementMatrix::EF);
#ifdef PZ_LOG
    bool globalresult = true;
    bool writereadresult = true;
#endif
    TPZTimer calcstiff("Computing the stiffness matrices");
    TPZTimer assemble("Assembling the stiffness matrices");
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();

    int64_t count = 0;
    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = materialIds.size();
        if(matidsize){
            if(!el->NeedsComputing(materialIds)) continue;
        }

        count++;
        if (!(count % 1000)) {
            std::cout << '*';
            std::cout.flush();
        }
        if (!(count % 20000)) {
            std::cout << "\n";
        }
        calcstiff.start();
        ek.Reset();
        ef.Reset();
        el->CalcStiff(ek, ef);

#ifdef PZ_LOG
        if (loggerelmat.isDebugEnabled()) {
            if (dynamic_cast<TPZSubCompMesh *> (cmesh)) {
                std::stringstream objname;
                objname << "Element" << iel;
                std::string name = objname.str();
                objname << " = ";
                std::stringstream sout;
                ek.fMat.Print(objname.str().c_str(), sout, EMathematicaInput);
                sout << "AppendTo[AllEig,Eigenvalues[" << name << "]];";

                LOGPZ_DEBUG(loggerelmat, sout.str())
                        /*		  if(iel == 133)
                         {
                         std::stringstream sout2;
                         el->Reference()->Print(sout2);
                         el->Print(sout2);
                         LOGPZ_DEBUG(logger,sout2.str())
                         }
                         */
            }
        }
#endif

#ifdef CHECKCONSISTENCY
        //extern TPZCheckConsistency stiffconsist("ElementStiff");
        stiffconsist.SetOverWrite(true);
        bool result;
        result = stiffconsist.CheckObject(ek.fMat);
        if (!result) {
            globalresult = false;
            std::stringstream sout;
            sout << "element " << iel << " computed differently";
            LOGPZ_ERROR(loggerCheck, sout.str())
        }

#endif

        calcstiff.stop();
        const auto &equationFilter =
        (dynamic_cast<TPZStructMatrix*>(this))->EquationFilter();
        assemble.start();

        if (!ek.HasDependency()) {
            ek.ComputeDestinationIndices();
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            //			TPZSFMatrix<TVar> test(stiffness);
            //			TPZFMatrix<TVar> test2(stiffness.Rows(),stiffness.Cols(),0.);
            //			stiffness.Print("before assembly",std::cout,EMathematicaInput);
            stiffness.AddKel(ek.fMat, ek.fSourceIndex, ek.fDestinationIndex);
#ifdef PZDEBUG
            REAL rhsnorm = Norm(ef.fMat);
            REAL eknorm = Norm(ek.fMat);
            if (std::isnan(rhsnorm) || std::isnan(eknorm)) {
                std::cout<<__PRETTY_FUNCTION__
                         <<"\nERROR:\n"
                         << "\telement " << iel << " has norm " << rhsnorm << std::endl;
                el->Print();
                ek.fMat.Print("ek",std::cout);
                ef.fMat.Print("ef",std::cout);
            }
#endif
            if(ComputeRhs())rhs.AddFel(ef.fMat, ek.fSourceIndex, ek.fDestinationIndex);
            //			stiffness.Print("stiffness after assembly STK = ",std::cout,EMathematicaInput);
            //			rhs.Print("rhs after assembly Rhs = ",std::cout,EMathematicaInput);
            //			test2.AddKel(ek.fMat,ek.fSourceIndex,ek.fDestinationIndex);
            //			test -= stiffness;
            //			test.Print("matriz de rigidez diference",std::cout);
            //			test2.Print("matriz de rigidez interface",std::cout);

        } else {
            // the element has dependent nodes
            ek.ApplyConstraints();
            ef.ApplyConstraints();
            ek.ComputeDestinationIndices();
            equationFilter.Filter(ek.fSourceIndex, ek.fDestinationIndex);
            stiffness.AddKel(ek.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
            if(ComputeRhs())rhs.AddFel(ef.fConstrMat, ek.fSourceIndex, ek.fDestinationIndex);
        }

#ifdef PZ_LOG
        if (loggerel.isDebugEnabled()) {
            std::stringstream sout;
            TPZGeoEl *gel = el->Reference();
            if (gel) {
                TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                gel->CenterPoint(gel->NSides() - 1, center);
                gel->X(center, xcenter);
                sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                sout << "Geometric element index " << gel->Index() << " center " << xcenter << std::endl;
            } else {
                sout << "Stiffness for computational element without associated geometric element index " << el->Index() << "\n";
            }
            if(ek.HasDependency()){
                sout << "source index " << ek.fSourceIndex << " destination index " << ek.fDestinationIndex << std::endl;
                ek.fConstrMat.Print(sout);
                ef.fConstrMat.Print(sout);
            }else{
                ek.fMat.Print("ek", sout, EMathematicaInput);
                ef.fMat.Print("ef", sout, EMathematicaInput);
            }
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif
        // tototototo
        //        GK.Multiply(Sol, GKSol);
        //        GKSol -= GF;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GKSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        stiffness.Multiply(Sol, GKSol);
        //        GKSol -= rhs;
        //        GKSol.Transpose();
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "StiffSol" << iel << " = ";
        //            GKSol.Print(str.str().c_str(),sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        {
        //            std::stringstream sout;
        //            sout << "Element " << iel << std::endl;
        //            std::stringstream str;
        //            str << "GK" << iel << " = ";
        //            GK.Print(str.str().c_str(),sout,EMathematicaInput);
        //            std::stringstream str2;
        //            str2 << "ST" << iel << " = ";
        //            stiffness.Print(str2.str().c_str(),sout,EMathematicaInput);
        //            sout << "GK-ST\n";
        //            LOGPZ_DEBUG(logger, sout.str())
        //        }
        //        
        //        stiffness.Zero();
        //        rhs.Zero();
        assemble.stop();
    }//fim for iel
    if (count > 1000) std::cout << std::endl;

#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        sout << "The comparaison results are : consistency check " << globalresult << " write read check " << writereadresult;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
    if (loggerGlobStiff.isDebugEnabled())
    {
        std::stringstream sout;
        stiffness.Print("GK = ",sout,EMathematicaInput);
        if(ComputeRhs())rhs.Print("GR = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(loggerGlobStiff,sout.str())
    }

#endif

}

template<class TVar>
void 
TPZStructMatrixOR<TVar>::Serial_Assemble(TPZBaseMatrix & rhs_base) {
    if(!dynamic_cast<TPZFMatrix<TVar>*>(&rhs_base)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<": incompatible types. Aborting...\n";
        DebugStop();
    }
    auto &rhs = dynamic_cast<TPZFMatrix<TVar>&>(rhs_base);
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    const auto &equationFilter = myself->EquationFilter();
    auto *cmesh = myself->Mesh();
    auto &materialIds = myself->MaterialIds();
    int64_t iel;
    int64_t nelem = cmesh->NElements();

    TPZTimer calcresidual("Computing the residual vector");
    TPZTimer assemble("Assembling the residual vector");

    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();

    for (iel = 0; iel < nelem; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el || el->IsInterface()) continue;

        int matid = 0;
        TPZGeoEl *gel = el->Reference();
        if (gel) {
            matid = gel->MaterialId();
        }
        int matidsize = materialIds.size();
        if(matidsize){
            TPZMaterial * mat = el->Material();
            TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (el);
            if (!mat)
            {
                if (!submesh) {
                    continue;
                }
                else if(submesh->NeedsComputing(materialIds) == false) continue;
            }
            else
            {
                if (myself->ShouldCompute(matid) == false) continue;
            }
        }
                
        TPZElementMatrixT<TVar> ef(cmesh, TPZElementMatrix::EF);

        calcresidual.start();

        el->CalcResidual(ef);

        calcresidual.stop();
#ifdef PZ_LOG
        if(loggerel.isDebugEnabled())
        {
            std::stringstream sout;
            TPZGeoEl *gel = el->Reference();
            if (gel) {
                TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                gel->CenterPoint(gel->NSides() - 1, center);
                gel->X(center, xcenter);
                sout << "Residual for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                sout << "Residual for geometric element " << gel->Index() << " center " << xcenter << std::endl;
            } else {
                sout << "Residual for computational element without associated geometric element index " << el->Index() << "\n";
            }
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif
        assemble.start();

        if (!ef.HasDependency()) {
            ef.ComputeDestinationIndices();
            equationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
        } else {
            // the element has dependent nodes
            ef.ApplyConstraints();
            ef.ComputeDestinationIndices();
            equationFilter.Filter(ef.fSourceIndex, ef.fDestinationIndex);
            rhs.AddFel(ef.fConstrMat, ef.fSourceIndex, ef.fDestinationIndex);
        }
#ifdef PZ_LOG
        if(loggerel.isDebugEnabled())
        {
            std::stringstream sout;
            if(ef.HasDependency()){
                ef.fConstrMat.Print(sout);
            }else{
                ef.fMat.Print(sout);
            }
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif

        assemble.stop();

    }//fim for iel
#ifdef PZ_LOG
    {
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << calcresidual.processName() << " " << calcresidual << std::endl;
            sout << assemble.processName() << " " << assemble;
            LOGPZ_DEBUG(logger, sout.str().c_str());
        }
    }
#endif
    //std::cout << std::endl;
}

template<class TVar>
void
TPZStructMatrixOR<TVar>::MultiThread_Assemble(TPZBaseMatrix & mat, TPZBaseMatrix & rhs) {
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    if(!myself){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Run-time error.Aborting...\n";
        DebugStop();
    }
    ThreadData threaddata(myself,mat,rhs,myself->MaterialIds(),ComputeRhs());
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads;
    int itr;
    for (itr = 0; itr < numthreads; itr++) {
      allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
    }

    ThreadData::ThreadAssembly(&threaddata);

    for (itr = 0; itr < numthreads; itr++) {
      allthreads[itr].join();
    }

#ifdef PZ_LOG2
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        //stiffness.Print("Matriz de Rigidez: ",sout);
        mat.Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        rhs.Print("Right Handside", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
}

template<class TVar>
void
TPZStructMatrixOR<TVar>::MultiThread_Assemble(TPZBaseMatrix & rhs) {
    auto *myself = dynamic_cast<TPZStructMatrix*>(this);
    ThreadData threaddata(myself, rhs, myself->MaterialIds(),ComputeRhs());
    const int numthreads = this->fNumThreads;
    std::vector<std::thread> allthreads;
    int itr;

    for (itr = 0; itr < numthreads; itr++) {
      allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
    }

    ThreadData::ThreadAssembly(&threaddata);

    for (itr = 0; itr < numthreads; itr++) {
      allthreads[itr].join();
    }
}

template<class TVar>
int TPZStructMatrixOR<TVar>::ClassId() const{
    return Hash("TPZStructMatrixOR") ^
        TPZStrMatParInterface::ClassId() << 1 ^
        ClassIdOrHash<TVar>() << 2;
}

template<class TVar>
void TPZStructMatrixOR<TVar>::Read(TPZStream& buf, void* context) {
    TPZStrMatParInterface::Read(buf,context);
}

template<class TVar>
void TPZStructMatrixOR<TVar>::Write(TPZStream& buf, int withclassid) const {
    TPZStrMatParInterface::Write(buf,withclassid);
}

template<class TVar>
TPZStructMatrixOR<TVar>::ThreadData::~ThreadData() {
}

template<class TVar>
TPZStructMatrixOR<TVar>::ThreadData::ThreadData(
  TPZStructMatrix *strmat,
  TPZBaseMatrix &mat, TPZBaseMatrix &rhs,
  const std::set<int> &MaterialIds,
  bool computeRhs)
: fStruct(strmat), fGlobMatrix(&mat),
  fGlobRhs(&rhs), fNextElement(0), fComputeRhs(computeRhs) {

}

template<class TVar>
TPZStructMatrixOR<TVar>::ThreadData::ThreadData(
  TPZStructMatrix *strmat,
  TPZBaseMatrix &rhs,
  const std::set<int> &MaterialIds,
  bool computeRhs)
: fStruct(strmat), fGlobMatrix(0),
  fGlobRhs(&rhs), fNextElement(0), fComputeRhs(computeRhs) {
}

#ifdef USING_MKL
#include <mkl.h>
#endif
//#define DRY_RUN
template<class TVar>
void *
TPZStructMatrixOR<TVar>::ThreadData::ThreadWork(void *datavoid) {
#ifdef PZDEBUG
    TExceptionManager activateExceptions;
#endif
#ifdef USING_MKL
    //we keep the previous value
    auto mkl_threads = mkl_set_num_threads_local(1);
#endif
    ThreadData *data = (ThreadData *) datavoid;
    const bool computeRhs = data->fComputeRhs;
    // compute the next element (this method is threadsafe)
    int64_t iel = data->NextElement();

    TPZCompMesh *cmesh = data->fStruct->Mesh();
    int64_t nel = cmesh->NElements();
    while (iel < nel) {

        TPZAutoPointer<TPZElementMatrixT<TVar>> ek;
        TPZAutoPointer<TPZElementMatrixT<TVar>> ef =
            new TPZElementMatrixT<TVar>(cmesh, TPZElementMatrix::EF);
        if (data->fGlobMatrix) {
            ek = new TPZElementMatrixT<TVar>(cmesh, TPZElementMatrix::EK);
        } else {
            ek = ef;
        }

        TPZCompEl *el = cmesh->ElementVec()[iel];

#ifndef DRY_RUN
        if (data->fGlobMatrix) {
            el->CalcStiff(ek, ef);
        } else {
            el->CalcResidual(ef);
        }
#else
        {
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (el);
            if (!intel) {
                DebugStop();
            }
            if (data->fGlobMatrix) {
                intel->InitializeElementMatrix(ekr, efr);
            } else {
                intel->InitializeElementMatrix(efr);
            }
        }
#endif


        if (!ek->HasDependency()) {
            ek->ComputeDestinationIndices();

            if (data->fStruct->HasRange()) {
                data->fStruct->FilterEquations(ek->fSourceIndex, ek->fDestinationIndex);
            }
        } else {
            // the element has dependent nodes
            if (data->fGlobMatrix) {
                ek->ApplyConstraints();
            }
            ef->ApplyConstraints();
            ek->ComputeDestinationIndices();
            if (data->fStruct->HasRange()) {
                data->fStruct->FilterEquations(ek->fSourceIndex, ek->fDestinationIndex);
            }
        }

#ifdef PZ_LOG
        if (loggerel2.isDebugEnabled() && el->Reference() && el->Reference()->MaterialId() == 1 && el->IsInterface()) {
            std::stringstream sout;
            el->Reference()->Print(sout);
            el->Print(sout);
            ek->Print(sout);
            //			ef->Print(sout);
            LOGPZ_DEBUG(loggerel2, sout.str())
        }
        if (loggerel.isDebugEnabled()) {
            std::stringstream sout;
            TPZGeoEl *gel = el->Reference();
            if (gel) {
                TPZManVector<REAL> center(gel->Dimension()), xcenter(3, 0.);
                gel->CenterPoint(gel->NSides() - 1, center);
                gel->X(center, xcenter);
                sout << "Stiffness for computational element index " << el->Index() << " material id " << gel->MaterialId() << std::endl;
                sout << "Geometric element index " << gel->Index() << " center " << xcenter << std::endl;
            } else {
                sout << "Stiffness for computational element without associated geometric element index " << el->Index() << "\n";
            }
            if(ek->HasDependency()){
                ek->fConstrMat.Print(sout);
                if(computeRhs){ef->fConstrMat.Print(sout);}
            }else{
                ek->fMat.Print(sout);
                if(computeRhs){ef->fMat.Print(sout);}
            }
            
            LOGPZ_DEBUG(loggerel, sout.str())
        }
#endif
        // put the elementmatrices on the stack to be assembled (threadsafe)
        data->ComputedElementMatrix(iel, ek, ef);
        // compute the next element (this method is threadsafe)
        iel = data->NextElement();
    }

    {
      std::scoped_lock lock(data->fMutexAccessElement);
      data->fAssembly.Post();
    }
#ifdef USING_MKL
    //we restore the previous value
    mkl_set_num_threads_local(mkl_threads);
#endif
    return 0;
}

// The function which will compute the assembly
template<class TVar>
void *
TPZStructMatrixOR<TVar>::ThreadData::ThreadAssembly(void *threaddata) {
    ThreadData *data = (ThreadData *) threaddata;
    const bool computeRhs = data->fComputeRhs;
    TPZCompMesh *cmesh = data->fStruct->Mesh();
    int64_t nel = cmesh->NElements();
    data->fMutexAccessElement.lock();
    int64_t nextel = data->fNextElement;
    int numprocessed = data->fProcessed.size();
    while (nextel < nel || numprocessed) {
        bool keeplooking = false;
        if (data->fSubmitted.size() && data->fProcessed.size()) {
            auto itavail = data->fSubmitted.begin();
            auto itprocess = data->fProcessed.begin();
            if (itavail->first == *itprocess) {
                // make sure we come back to look for one more element
                keeplooking = true;
                // Get a hold of the data
#ifdef PZ_LOG
                int iel = *itprocess;
#endif
                data->fProcessed.erase(itprocess);
                TPZAutoPointer<TPZElementMatrixT<TVar>> ek = itavail->second.first;
                TPZAutoPointer<TPZElementMatrixT<TVar>> ef = itavail->second.second;
                data->fSubmitted.erase(itavail);
#ifdef PZ_LOG
                if (logger.isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Assembling element " << iel;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
#ifdef CHECKCONSISTENCY
                //static TPZCheckConsistency stiffconsist("ElementStiff");
                stiffconsist.SetOverWrite(true);
                bool result;
                result = stiffconsist.CheckObject(ek->fMat);
                if (!result) {
                    globalresult = false;
                    std::stringstream sout;
                    sout << "element " << iel << " computed differently";
                    LOGPZ_ERROR(loggerCheck, sout.str())
                }
#endif

                // Release the mutex
                data->fMutexAccessElement.unlock();
#ifndef DRY_RUN
                auto globMatrix =
                    dynamic_cast<TPZMatrix<TVar> *>(data->fGlobMatrix);
                auto globRhs =
                    dynamic_cast<TPZFMatrix<TVar> *>(data->fGlobRhs);
                // Assemble the matrix
                if (!ek->HasDependency()) {
                    if (globMatrix) {
                        globMatrix->AddKel(ek->fMat, ek->fSourceIndex, ek->fDestinationIndex);
                    }
                    if(computeRhs)
                        globRhs->AddFel(ef->fMat, ek->fSourceIndex, ek->fDestinationIndex);
                } else {
                    if (globMatrix) {
                        globMatrix->AddKel(ek->fConstrMat, ek->fSourceIndex, ek->fDestinationIndex);
                    }
                    if(computeRhs)
                        globRhs->AddFel(ef->fConstrMat, ek->fSourceIndex, ek->fDestinationIndex);
                }
#endif
                // acquire the mutex
                data->fMutexAccessElement.lock();
            }
        }
        if (!keeplooking) {
          data->fMutexAccessElement.unlock();
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger, "Going to sleep within assembly")
            }
#endif
            // wait for a signal
            data->fAssembly.Wait();
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger, "Waking up for assembly")
            }
#endif
            data->fMutexAccessElement.lock();
        }
        nextel = data->fNextElement;
        numprocessed = data->fProcessed.size();

    }
    //	std::cout << std::endl;
#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        sout << "nextel = " << nextel << " numprocessed = " << numprocessed << " submitted " << data->fSubmitted.size() << std::endl;
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
#ifdef PZ_LOG
    if (loggerCheck.isDebugEnabled()) {
        std::stringstream sout;
        if (data->fGlobMatrix) {
            data->fGlobMatrix->Print("Matriz de Rigidez: ", sout, EMathematicaInput);
        }
        data->fGlobRhs->Print("Right hand side", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerCheck, sout.str())
    }
#endif
    data->fMutexAccessElement.unlock();
    return 0;
}

template<class TVar>
int64_t 
TPZStructMatrixOR<TVar>::ThreadData::NextElement() {
    fMutexAccessElement.lock();
    int64_t iel;
    int64_t nextel = fNextElement;
    TPZCompMesh *cmesh = fStruct->Mesh();
    TPZAdmChunkVector<TPZCompEl *> &elementvec = cmesh->ElementVec();
    int64_t nel = elementvec.NElements();
    for (iel = fNextElement; iel < nel; iel++) {
        TPZCompEl *el = elementvec[iel];
        if (!el) continue;
        if (fStruct->MaterialIds().size() == 0) break;
        if (el->NeedsComputing(fStruct->MaterialIds())) break;
    }
    fNextElement = iel + 1;
    nextel = iel;
    if (iel < nel) fProcessed.insert(iel); //AQUIBORIN pelo que percebi, aqui tem que acontecer antes do unlock no caso sem Critical Section
    fMutexAccessElement.unlock();
#ifdef PZ_LOG
    {
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " returning " << nextel << " fNextElement " << fNextElement;
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
#endif
    return nextel;
}


// put the computed element matrices in the map

template<class TVar>
void 
TPZStructMatrixOR<TVar>::ThreadData::ComputedElementMatrix(int64_t iel, TPZAutoPointer<TPZElementMatrixT<TVar>> &ek, TPZAutoPointer<TPZElementMatrixT<TVar>> &ef) {
    std::scoped_lock lock(fMutexAccessElement);
    std::pair< TPZAutoPointer<TPZElementMatrixT<TVar>>, TPZAutoPointer<TPZElementMatrixT<TVar>> > el(ek, ef);
    fSubmitted[iel] = el;
    fAssembly.Post();
}

template<class TVar>
bool 
TPZStructMatrixOR<TVar>::ThreadData::ShouldCompute(int matid) const{
    return fStruct->ShouldCompute(matid);
}

template class TPZStructMatrixOR<STATE>;
template class TPZRestoreClass<TPZStructMatrixOR<STATE>>;
template class TPZStructMatrixOR<CSTATE>;
template class TPZRestoreClass<TPZStructMatrixOR<CSTATE>>;
