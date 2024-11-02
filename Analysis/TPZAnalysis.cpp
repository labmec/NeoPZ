/**
 * \file
 * @brief Contains implementations of the TPZAnalysis methods.
 */

#include "TPZAnalysis.h"
#include <math.h>                          // for sqrt, fabs
#include <thread>
#include <vector>
#include <stdio.h>                         // for NULL
#include <string.h>                        // for strcpy, strlen
#ifdef MACOSX
//#include <__functional_base>               // for less
#include <__tree>                          // for __tree_const_iterator, ope...
#endif
#include <list>                            // for list, __list_iterator, lis...
#include <map>                             // for __map_iterator, map, map<>...
#include <string>                          // for allocator, basic_string
#include <utility>                         // for pair
#include "TPZLagrangeMultiplier.h"         // for TPZLagrangeMultiplier
#include "TPZSloanRenumbering.h"           // for TPZSloanRenumbering
#include "TPZCutHillMcKee.h"
#ifdef PZ_USING_METIS
#include "pzmetis.h"
#endif
#include "pzadmchunk.h"                    // for TPZAdmChunkVector
#include "pzbdstrmatrix.h"                 // for TPZBlockDiagonalStructMatrix
#include "pzblock.h"                       // for TPZBlock
#include "pzblockdiag.h"                   // for TPZBlockDiagonal
#include "TPZBndCond.h"                     // for TPZBndCond
#include "TPZChunkVector.h"                       // for TPZChunkVector
#include "pzcmesh.h"                       // for TPZCompMesh
#include "pzcompel.h"                      // for TPZCompEl
#include "pzconnect.h"                     // for TPZConnect
#include "pzdxmesh.h"                      // for TPZDXGraphMesh
#include "TPZEquationFilter.h"              // for TPZEquationFilter
#include "pzgeoel.h"                       // for TPZGeoEl
#include "pzgmesh.h"                       // for TPZGeoMesh
#include "pzgraphmesh.h"                   // for TPZGraphMesh
#include "pzlog.h"                         // for glogmutex, LOGPZ_DEBUG
#include "pzmanvector.h"                   // for TPZManVector
#include "TPZMaterial.h"                    // for TPZMaterial
#include "TPZMatLoadCases.h"               // for TPZMatLoadCases
#include "pzmvmesh.h"                      // for TPZMVGraphMesh
#include "pzseqsolver.h"                   // for TPZSequenceSolver
#include "TPZMatrixSolver.h"                       // for TPZMatrixSolver, TPZSolver
#include "pzstack.h"                       // for TPZStack
#include "pzstepsolver.h"                  // for TPZStepSolver
#include "TPZStructMatrix.h"               // for TPZStructMatrix
#include "pzv3dmesh.h"                     // for TPZV3DGraphMesh
#include "pzvec.h"                         // for TPZVec, operator<<
#include "pzvtkmesh.h"                     // for TPZVTKGraphMesh
#include "tpznodesetcompute.h"             // for TPZNodesetCompute
#include "tpzsparseblockdiagonal.h"        // for TPZSparseBlockDiagonal
#include "TPZMatError.h"
#include "TPZSimpleTimer.h"
#include "pzelementgroup.h"
#include "TPZJacobiPrecond.h"
#ifdef WIN32
#include "pzsloan.h"                       // for TPZSloan
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.analysis");
static TPZLogger loggerError("pz.analysis.error");
static TPZLogger loggerPrecond("pz.analysis.precondgraph");
#endif

//@orlandini: does anyone know if boost renumbering still works?
//@shauer: I don't think so. This #define structure is also now deprecated. I vote for deleting
#undef USE_BOOST_RENUMBERING
#if defined(USING_BOOST) && defined(USE_BOOST_RENUMBERING)
#include "TPZBoostGraph.h"
/**
 * @brief Renumbering will use boost library.
 * @ingroup analysis
 */
#define RENUMBER TPZBoostGraph(TPZBoostGraph::KMCExpensive))
#else
/**
 * @brief Renumbering will use sloan library.
 * @ingroup analysis
 */

//#ifdef PZ_USING_METIS
//#define RENUMBER TPZMetis()
//#else
//#define RENUMBER TPZSloanRenumbering()
//#endif
//#define RENUMBER TPZCutHillMcKee() //TPZCutHillMcKee usually performs worse than metis and sloan. Use carefully

#endif

using namespace std;

void TPZAnalysis::SetStructuralMatrix(TPZStructMatrix &strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix.Clone());
    /*let us leave this check even in release for a while, until we
      can be sure that our implementation of complex vars is ok*/
// #ifdef PZDEBUG
    auto* tmp =
        dynamic_cast<TPZStructMatrixT<STATE> *>(&strmatrix);
    if(tmp){//real struct matrix
        if(fCompMesh && fCompMesh->GetSolType() == EReal){
            return;
        }
    }else{//complex struct matrix
        if(fCompMesh && fCompMesh->GetSolType() == EComplex){
            return;
        }
    }

    PZError<<__PRETTY_FUNCTION__;
    PZError<<" Incompatible types. Aborting\n...";
    DebugStop();
// #endif
}

void TPZAnalysis::SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix){
	fStructMatrix = TPZAutoPointer<TPZStructMatrix>(strmatrix->Clone());
// #ifdef PZDEBUG
    auto* tmp =
        dynamic_cast<TPZStructMatrixT<STATE> *>((strmatrix.operator->()));
    if(tmp){//real struct matrix
        if(fCompMesh && fCompMesh->GetSolType() == EReal){
            return;
        }
    }else{//complex struct matrix
        if(fCompMesh && fCompMesh->GetSolType() == EComplex){
            return;
        }
    }
// #endif
}
TPZAnalysis::TPZAnalysis() :
    TPZRegisterClassId(&TPZAnalysis::ClassId){
}

TPZAnalysis::TPZAnalysis(TPZCompMesh *mesh, const RenumType& renumtype, std::ostream &out) :
TPZRegisterClassId(&TPZAnalysis::ClassId),
fSolType(mesh->GetSolType()),
// fRhs(fSolType == EComplex ? true : false),
fSolution(fSolType == EComplex ? true : false)
{
  //we must not call virtual methods in constructor
  bool mustOptimizeBandwidth = false;
  if (renumtype != RenumType::ENone){
    CreateRenumberObject(renumtype);
    mustOptimizeBandwidth = true;
  }
  this->SetCompMeshInit(mesh, mustOptimizeBandwidth);
}

TPZAnalysis::TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh, const RenumType& renumtype, std::ostream &out) :
TPZRegisterClassId(&TPZAnalysis::ClassId),
fSolType(mesh->GetSolType()),
// fRhs(fSolType == EComplex ? true : false),
fSolution(fSolType == EComplex ? true : false)
{
  //we must not call virtual methods in constructor
  bool mustOptimizeBandwidth = false;
  if (renumtype != RenumType::ENone){
    CreateRenumberObject(renumtype);
    mustOptimizeBandwidth = true;
  }
	this->SetCompMeshInit(mesh.operator ->(), mustOptimizeBandwidth);
}

void TPZAnalysis::CreateRenumberObject(const RenumType& renumtype) {
  switch (renumtype) {
    case RenumType::ENone:
      fRenumber = nullptr;
      return;
    case RenumType::ECutHillMcKee:
      fRenumber = new TPZCutHillMcKee(true);
      return;
    case RenumType::ECutHillMcKeeFast:
      fRenumber = new TPZCutHillMcKee(false);
      return;
    case RenumType::EDefault:
    case RenumType::EMetis:
#ifdef PZ_USING_METIS
      fRenumber = new TPZMetis;
      return;
#else
      if(renumtype==RenumType::EMetis){
        //maybe we ended up here because of EDefault,
        //no need to print the message then
        std::cout<<__PRETTY_FUNCTION__
                 <<":\nMetis not available, setting Sloan for renumbering"
                 <<std::endl;
      }
#endif
    case RenumType::ESloan:
      fRenumber = new TPZSloanRenumbering;
      return;
  }
  DebugStop();
}

void TPZAnalysis::SetCompMeshInit(TPZCompMesh *mesh, bool mustOptimizeBandwidth)
{
  if(mesh)
    {
        fCompMesh = mesh;
        fSolType = fCompMesh->GetSolType();
        fGeoMesh = mesh->Reference();
        if(fGraphMesh[0]){
            delete fGraphMesh[0];
            fGraphMesh[0] = nullptr;
        }
        if(fGraphMesh[1]){
            delete fGraphMesh[1];
            fGraphMesh[1] = nullptr;
        }
        if(fGraphMesh[2]){
            delete fGraphMesh[2];
            fGraphMesh[2] = nullptr;
        }
        if(fSolver) fSolver->ResetMatrix();
//        fCompMesh->InitializeBlock();
        int64_t neq = fCompMesh->NEquations();
        if(neq > 20000 && mustOptimizeBandwidth)
        {
            std::cout << __PRETTY_FUNCTION__ << " optimizing bandwidth neq = " << neq << "\n";
            std::cout.flush();
        }
        if(mustOptimizeBandwidth)
        {
            TPZSimpleTimer bd("Optimize bandwidth",false);
            OptimizeBandwidth();
        }
        if(neq > 20000 && mustOptimizeBandwidth)
        {
            std::cout << __PRETTY_FUNCTION__ << " optimizing bandwidth finished\n";
            std::cout.flush();
        }
        fSolution = fCompMesh->Solution();
        fSolution.Resize(neq,1);
    }
    else
    {
        CleanUp();
    }
    fStep = 0;
}

void TPZAnalysis::SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth) {
  SetCompMeshInit(mesh,mustOptimizeBandwidth);
}

TPZAnalysis::~TPZAnalysis(void){
    CleanUp();
}

/// deletes all data structures
void TPZAnalysis::CleanUp()
{
    if(fSolver) {
        delete fSolver;
        fSolver = nullptr;
    }
    int dim;
    for(dim=0; dim<3; dim++) {
        if(fGraphMesh[dim]) delete fGraphMesh[dim];
        fGraphMesh[dim] = nullptr;
        fScalarNames[dim].resize(0);
        fVectorNames[dim].resize(0);
    }
    fCompMesh = nullptr;
    fGeoMesh = nullptr;
    fSolution.Redim(0,0);
    fStructMatrix = nullptr;
    fRenumber = nullptr;
    
}


void TPZAnalysis::OptimizeBandwidth() {
	//enquanto nao compilamos o BOOST no windows, vai o sloan antigo
#ifdef WIN32
	if(!fCompMesh) return;
//    fCompMesh->InitializeBlock();
	TPZVec<int64_t> perm,iperm;
	
	TPZStack<int64_t> elgraph;
	TPZStack<int64_t> elgraphindex;
	int64_t nindep = fCompMesh->NIndependentConnects();
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int64_t nel = elgraphindex.NElements()-1;
	TPZSloan sloan(nel,nindep);
	sloan.SetElementGraph(elgraph,elgraphindex);
	sloan.Resequence(perm,iperm);
	fCompMesh->Permute(perm);
#else
	if(!fCompMesh) return;
//    fCompMesh->InitializeBlock();
	
	TPZVec<int64_t> perm,iperm;
	
	TPZStack<int64_t> elgraph,elgraphindex;
	int64_t nindep = fCompMesh->NIndependentConnects();
    
    /// if there are no connects, there is no bandwidth to be optimized
    if(nindep == 0) return;
    
	fCompMesh->ComputeElGraph(elgraph,elgraphindex);
	int64_t nel = elgraphindex.NElements()-1;
	int64_t el,ncel = fCompMesh->NElements();
	int maxelcon = 0;
	for(el = 0; el<ncel; el++)
	{
		TPZCompEl *cel = fCompMesh->ElementVec()[el];
		if(!cel) continue;
		std::set<int64_t> indepconlist,depconlist;
		cel->BuildConnectList(indepconlist,depconlist);
		int64_t locnindep = indepconlist.size();
		maxelcon = maxelcon < locnindep ? locnindep : maxelcon;
	}
	fRenumber->SetElementsNodes(nel,nindep);
	fRenumber->SetElementGraph(elgraph,elgraphindex);
#ifdef PZ_LOG2
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        fRenumber->Print(elgraph, elgraphindex, "Elgraph of submesh", sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    bool shouldstop = false;
    try {
        fRenumber->Resequence(perm,iperm);
    } catch(...)
    {
        fRenumber->PlotElementGroups("ElementGroups.vtk", fCompMesh);
        std::ofstream out("fElementGroups.txt");
        fCompMesh->Print(out);
        shouldstop = true;
    }
    if(shouldstop) DebugStop();
    
	fCompMesh->Permute(perm);
    if (nel > 100000) {
        std::cout << "Applying Saddle Permute\n";
    }
    fCompMesh->SaddlePermute();
//    fCompMesh->SaddlePermute2();
	
#endif
	
}

/** @brief Determine the number of load cases from the material objects and return its value */
/**
 * this method will modify the material objects so that they have all the same number of load cases
 * the number of load cases is the maximum value of load cases of all material objects
 */
int TPZAnalysis::ComputeNumberofLoadCases()
{
    int res = 1;
    if(!fCompMesh) 
    {
        return res;
    }
    bool everyMatHasLoadCase{false};
    for(auto &it : fCompMesh->MaterialVec()){
        if(auto *matLoad = dynamic_cast<TPZMatLoadCasesBase*>(it.second);
           !matLoad) {
            everyMatHasLoadCase = false;
            break;
        }
    }
    if(everyMatHasLoadCase){
        for(auto &it : fCompMesh->MaterialVec()){
            auto *matLoad = dynamic_cast<TPZMatLoadCasesBase*>(it.second);
            const int matNCases = matLoad->MinimumNumberofLoadCases();
            res = res < matNCases ? matNCases : res;
        }
        for(auto &it : fCompMesh->MaterialVec()){
            auto *matLoad = dynamic_cast<TPZMatLoadCasesBase*>(it.second);
            matLoad->SetNumLoadCases(res);
        }
    }

    return res;
}


void TPZAnalysis::LoadSolution() {	
	if(fCompMesh) {
		fCompMesh->LoadSolution(fSolution);
	}
}

void TPZAnalysis::Print( const std::string &name, std::ostream &out) {
	out<<endl<<name<<endl;
    if (!fCompMesh) {
        out << "No computational mesh\n";
    }
    else
    {
        int64_t i,nelements = fCompMesh->ConnectVec().NElements();
        for(i=0;i<nelements;i++) {
            TPZConnect &gnod = fCompMesh->ConnectVec()[i];
            if(gnod.SequenceNumber()!=-1) {
                out << "Connect index " << i << endl;
                gnod.Print(*fCompMesh,out);
            }
        }
    }
    TPZBaseMatrix &sol = this->Solution();
	sol.Print("fSolution",out);
    if (fCompMesh) {
        fCompMesh->ConnectSolution(out);
    }
}


void TPZAnalysis::PostProcess(TPZVec<REAL> &ervec, std::ostream &out) {
	int64_t i;
	//TPZVec<REAL> ux((int) neq);
	//TPZVec<REAL> sigx((int) neq);
	TPZManVector<REAL,10> values(10,0.);
	TPZManVector<REAL,10> values2(10,0.);
	fCompMesh->LoadSolution(fSolution);
	TPZAdmChunkVector<TPZCompEl *> &elvec = fCompMesh->ElementVec();
	TPZManVector<REAL,10> errors(10);
	errors.Fill(0.0);
	int64_t nel = elvec.NElements();
	int matId0 = 0;
	for(i=0;i<nel;i++) {
        TPZCompEl *cel = elvec[i];
        if(!cel) continue;
		matId0=cel->Material()->Id();
		if(matId0 > 0)
			break;
	}
	bool lastEl=false;
	for(i=0;i<nel;i++) {
		TPZCompEl *el = (TPZCompEl *) elvec[i];
		if(el) {
			errors.Fill(0.0);
      TPZGeoEl *gel = el->Reference();
      if(gel->Dimension() != fCompMesh->Dimension()) continue;
      auto *matError =
          dynamic_cast<TPZMatError<STATE>*>(el->Material());
      if (matError || matError->HasExactSol()) {
        el->EvaluateError(errors, 0);
      } else {
          PZError<<__PRETTY_FUNCTION__;
          PZError<<" the material has no associated exact solution\n";
          PZError<<"Aborting...";
          DebugStop();
      }
      if(matId0==el->Material()->Id()){
				for(int ier = 0; ier < errors.NElements(); ier++) 	values[ier] += errors[ier] * errors[ier];
				lastEl=false;
			}
			else {
				for(int ier = 0; ier < errors.NElements(); ier++)	values2[ier] += errors[ier] * errors[ier];
				lastEl=true;
			}
		}
	}
	int nerrors = errors.NElements();
	ervec.Resize(2*nerrors);
	ervec.Fill(-10.0);
	
	if (nerrors==4) {
		if(lastEl){
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values2[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values2[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values2[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values2[3])  <<endl;

			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values[2])  <<endl;
		}
		else{
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values2[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values2[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values2[2])  <<endl;
		}
	}
	else {
		if(lastEl){
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values2[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values2[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values2[2])  <<endl;
		}
		else{
			out << endl << "############" << endl;
			out << endl << "L2 Norm for pressure  = "  << sqrt(values2[0]) << endl;
			out << endl << "L2 Norm for flux = "    << sqrt(values2[1]) << endl;
			out << endl << "L2 Norm for divergence = "    << sqrt(values2[2])  <<endl;
			out << endl << "Hdiv Norm for flux = "    << sqrt(values2[3])  <<endl;
			
			out << endl << "############" << endl;
			out << endl << "true_error (Norma H1) = "  << sqrt(values[0]) << endl;
			out << endl << "L2_error (Norma L2) = "    << sqrt(values[1]) << endl;
			out << endl << "estimate (Semi-norma H1) = "    << sqrt(values[2])  <<endl;
		}
		
		for(i=0;i<nerrors;i++) {
			ervec[i] = sqrt(values[i]);
		}
		for(i=nerrors;i<2*nerrors;i++){
			ervec[i] = sqrt(values2[i-nerrors]);
		}
	}
	return;
}

void TPZAnalysis::PostProcessError(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){
  if(!fNthreadsError){
    PostProcessErrorSerial(ervec, store_error, out);
  }
  else{
    PostProcessErrorParallel(ervec, store_error, out);
  }
}

void *TPZAnalysis::ThreadData::ThreadWork(void *datavoid)
{
  ThreadData *data = (ThreadData *) datavoid;
  const int64_t nelem = data->fElvec.NElements();
  TPZManVector<REAL,10> errors(10);
 
  // Getting unique id for each thread
  const int64_t myid = [&]()
  {
    std::scoped_lock lock(data->fMutexThreadId);
    const int64_t myid_loc = data->ftid;
    data->ftid++;
    return myid_loc;
  }();

  
  
  do{

    const int64_t iel = [&]()
    {
      std::scoped_lock lock(data->fMutexAccessEl);
      const int64_t iel_loc = data->fNextElement;
      data->fNextElement++;
      return iel_loc;
    }();
    
    // For all the elements it tries to get after the last one
    if ( iel >= nelem ) continue;
    
    TPZCompEl *cel = data->fElvec[iel];

    if(!cel) continue;
    if(cel->Dimension() != cel->Mesh()->Dimension()) continue;
    
    data->fIsUsed[myid] = true;
      
    cel->EvaluateError(errors, data->fStoreError);
    
    const int nerrors = errors.NElements();
    data->fvalues[myid].Resize(nerrors, 0.);


    for(int ier = 0; ier < nerrors; ier++)
    {
      (data->fvalues[myid])[ier] += errors[ier] * errors[ier];
    }
    
    
  } while (data->fNextElement < nelem);
  return data;
}

void TPZAnalysis::PostProcessErrorParallel(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){ //totto

  fCompMesh->LoadSolution(fSolution);
  const int numthreads = this->fNthreadsError;
  std::vector<std::thread> allthreads;

  TPZAdmChunkVector<TPZCompEl *> & elvec = fCompMesh->ElementVec();
  
  ThreadData threaddata(elvec,store_error);
  threaddata.fvalues.Resize(numthreads);
  threaddata.fIsUsed.Resize(numthreads,false);
  for(int iv = 0 ; iv < numthreads ; iv++){
      threaddata.fvalues[iv].Resize(10);
      threaddata.fvalues[iv].Fill(0.0);
  }


  {
      TPZSimpleTimer t("ThreadWork",true);
  
      for(int itr=0; itr<numthreads; itr++)
          {
              allthreads.push_back(std::thread(ThreadData::ThreadWork, &threaddata));
          }
  
      for(int itr=0; itr<numthreads; itr++)
          {
              allthreads[itr].join();
          }
  }
  
  
  // Sanity check. There should be number of ids equal to number of threads
  if(threaddata.ftid != numthreads){
    DebugStop();
  }

  
  TPZManVector<REAL,10> values;
  // Assuming the first is equal to the others
  int nerrors = threaddata.fvalues[0].NElements();
    for(int it = 0 ; it < numthreads ; it++)
    {
        if(!threaddata.fIsUsed[it]) continue;
        int locmin = threaddata.fvalues[it].NElements();
        nerrors = nerrors > locmin ? nerrors : locmin;
    }
  values.Resize(nerrors,0);
  // Summing up all the values of all threads
  for(int it = 0 ; it < numthreads ; it++){
      if(!threaddata.fIsUsed[it]) continue;
      for(int ir = 0 ; ir < nerrors ; ir++){
          values[ir] += (threaddata.fvalues[it])[ir];
      }
  }
  
  //const int nerrors = values.NElements();
  ervec.Resize(nerrors);
  ervec.Fill(-10.0);
    
  if (nerrors < 3) {
    PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
    out<<endl<<"############"<<endl;

#ifdef PZDEBUG
    for(int ier = 0; ier < nerrors; ier++)
      out << endl << "error " << ier << "  = " << sqrt(values[ier]);
#endif
  }
  else{
#ifdef PZDEBUG
    out << "############" << endl;
    out <<"Norma H1 or L2 -> p = "  << sqrt(values[0]) << endl;
    out <<"Norma L2 or L2 -> u = "    << sqrt(values[1]) << endl;
    out << "Semi-norma H1 or L2 -> div = "    << sqrt(values[2])  <<endl;
    for(int ier = 3; ier < nerrors; ier++)
      out << "other norms = " << sqrt(values[ier]) << endl;
#endif
  }
  
  // Returns the calculated errors.
  for(int i=0;i<nerrors;i++)
    ervec[i] = sqrt(values[i]);
  return;
  
}

#include "pzsubcmesh.h"

void TPZAnalysis::PostProcessErrorSerial(TPZVec<REAL> &ervec, bool store_error, std::ostream &out ){

    fCompMesh->EvaluateError(store_error, ervec);


    const int nerrors = ervec.NElements();

    if (nerrors < 3) {
        PZError << endl << "TPZAnalysis::PostProcess - At least 3 norms are expected." << endl;
        out<<endl<<"############"<<endl;
        for(int ier = 0; ier < nerrors; ier++)
            out << endl << "error " << ier << "  = " << ervec[ier];
    }
    else {
#ifdef PZDEBUG
        out << "############" << endl;
        out << "Norma H1 or L2 -> p = " << ervec[0] << std::endl;
        out << "Norma L2 or L2 -> u = " << ervec[1] << std::endl;
        out << "Semi-norma H1 or L2 -> div = " << ervec[2] << std::endl;
        for (int ier = 3; ier < ervec.size(); ier++)
            out << "other norms = " << ervec[ier] << std::endl;
#endif
    }
}

void TPZAnalysis::PostProcessTable( TPZFMatrix<REAL> &,std::ostream & )//pos,out
{
	TPZAdmChunkVector<TPZCompEl *> elvec = fCompMesh->ElementVec();
	int64_t nel = elvec.NElements();
	for(int64_t i=0;i<nel;i++) {
		//TPZCompEl *el = (TPZCompEl *) elvec[i];
		//if(el)	el->PrintTable(fExact,pos,out);
	}
	return;
}


void TPZAnalysis::SetExact(std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)>f,int p){
    
    if(fCompMesh->GetSolType() == ESolType::EReal){
        TPZAnalysis::SetExactInternal<STATE>(f,p);
    }else{
        PZError<<__PRETTY_FUNCTION__
               <<"\nError: incompatible type. Trying to set real solution in complex mesh.\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
}


void TPZAnalysis::SetExact(std::function<void (const TPZVec<REAL> &loc, TPZVec<CSTATE> &result, TPZFMatrix<CSTATE> &deriv)>f,int p){
    
    if(fCompMesh->GetSolType() == ESolType::EComplex){
        TPZAnalysis::SetExactInternal<CSTATE>(f,p);
    }else{
        PZError<<__PRETTY_FUNCTION__
               <<"\nError: incompatible type. Trying to set complex solution in real mesh.\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
}

template<class TVar>
void TPZAnalysis::SetExactInternal(std::function<void (const TPZVec<REAL> &loc, TPZVec<TVar> &result, TPZFMatrix<TVar> &deriv)>f,int p){
    if(!fCompMesh){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" there is no associated computational mesh.\n";
        PZError<<"Aborting...\n";
        DebugStop();
    }
    for(auto imat : fCompMesh->MaterialVec()){
        auto *matError =
            dynamic_cast<TPZMatError<TVar>*>(imat.second);
        if(matError) matError->SetExactSol(f,p);
    }
}
template<class TVar>
void TPZAnalysis::ShowShapeInternal(
    const TPZStack<std::string> &scalnames,
    const TPZStack<std::string> &vecnames,
    const std::string &plotfile,
    TPZVec<int64_t> &equationindices){
  DefineGraphMesh(fCompMesh->Dimension(), scalnames, vecnames, plotfile);
  int porder = fCompMesh->GetDefaultOrder();

  int neq = equationindices.size();
  TPZFMatrix<TVar> solkeep(fSolution);
  TPZFMatrix<TVar> &mysol = fSolution;
  mysol.Zero();
  for (int ieq = 0; ieq < neq; ieq++) {
    mysol(equationindices[ieq], 0) = 1.;
    LoadSolution();
    Mesh()->TransferMultiphysicsSolution();
    PostProcess(porder + 1);
    fSolution.Zero();
  }
  fSolution = solkeep;
  LoadSolution();
}
void TPZAnalysis::ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices)
{
	
    SetStep(1);
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = fCompMesh->MaterialVec().begin()->second;
    int nstate = mat->NStateVariables();
    if (nstate == 1) {
        scalnames.Push("Solution");
    }
    else
    {
        vecnames.Push("Solution");
    }
    //TODOCOMPLEX
    ShowShapeInternal<STATE>(scalnames,vecnames,plotfile,equationindices);
}

void TPZAnalysis::ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices, int matid, const TPZVec<std::string> &varname)
{
    TPZMaterial *mat = fCompMesh->FindMaterial(matid);
    if(!mat) DebugStop();
    
    int n_varnames = varname.size();
    TPZStack<std::string> scalnames,vecnames;
    for(int i_var = 0 ; i_var<n_varnames ; i_var++){
        int varindex = mat->VariableIndex(varname[i_var]);
        if(varindex == -1) DebugStop();
        SetStep(1);
        
        int nstate = mat->NSolutionVariables(varindex);
        if (nstate == 1) {
            scalnames.Push(varname[i_var]);
        }
        else
        {
            vecnames.Push(varname[i_var]);
        }
    }
    //TODOCOMPLEX
    ShowShapeInternal<STATE>(scalnames,vecnames,plotfile,equationindices);
}

void TPZAnalysis::Run(std::ostream &out)
{
    int64_t neq = fCompMesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }

    {
        TPZSimpleTimer t("Time for assembly",true);
        Assemble();
    }
    
    if(neq > 20000)
    {
        std::cout << "Entering Solve\n";
        std::cout.flush();
    }
    
    {
        TPZSimpleTimer t("Time for solving",true);
        Solve();
    }
}

/** @brief Define GrapMesh as V3D, DX, MV or VTK depending on extension of the file */
void TPZAnalysis::DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile){
    std::set<int> matids;
    IdentifyPostProcessingMatIds(dimension, matids);
    this->DefineGraphMesh(dimension, matids, scalnames, vecnames, plotfile);
    
}
/** @brief Define GrapMesh as VTK with tensorial names depending on extension of the file */
void TPZAnalysis::DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile){
    std::set<int> matids;
    IdentifyPostProcessingMatIds(dimension, matids);
    this->DefineGraphMesh(dimension, matids, scalnames, vecnames, tensnames, plotfile);
}

void TPZAnalysis::DefineGraphMesh(int dim, const std::set<int> & matids , const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile) {

    int dim1 = dim - 1;
    if (!fCompMesh || matids.size() == 0) {
        cout << "TPZAnalysis::DefineGraphMesh nothing to post-process\n";
        return;
    }

    if (fGraphMesh[dim1]) delete fGraphMesh[dim1];
    fScalarNames[dim1] = scalnames;
    fVectorNames[dim1] = vecnames;
    int posplot = plotfile.rfind(".plt");
    int64_t filelength = plotfile.size();
    if (filelength - posplot == 3) {
        if(tensnames.size() != 0)
        {
            std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
        }
        fGraphMesh[dim1] = new TPZV3DGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
    } else {
        int posdx = plotfile.rfind(".dx");
        if (filelength - posdx == 3) {
            if(tensnames.size() != 0) {
                std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
                
            }
            fGraphMesh[dim1] = new TPZDXGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
        } else {
            int pospos = plotfile.rfind(".pos");
            if (filelength - pospos == 3) {
                if(tensnames.size() != 0) {
                    std::cout << "Tensors are not beint post-process for file = " << plotfile << std::endl;
                }
                fGraphMesh[dim1] = new TPZMVGraphMesh(fCompMesh, dim, matids, scalnames, vecnames);
            } else {
                int posvtk = plotfile.rfind(".vtk");
                if (filelength - posvtk == 4) {
                    fGraphMesh[dim1] = new TPZVTKGraphMesh(fCompMesh, dim, matids, scalnames, vecnames, tensnames);
                } else {
                    cout << "grafgrid was not created\n";
                    fGraphMesh[dim1] = 0;
                }
            }
        }
    }
    if (fGraphMesh[dim1]) {
        fGraphMesh[dim1]->SetFileName(plotfile);
    }
}

void TPZAnalysis::DefineGraphMesh(int dim, const std::set<int> & matids, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile){
    TPZVec<std::string> tensnames;
    this->DefineGraphMesh(dim,matids, scalnames, vecnames, tensnames, plotfile);
}



void TPZAnalysis::CloseGraphMesh(){
	for(int i = 0; i < 3; i++){
		if ( this->fGraphMesh[i] ){
			delete this->fGraphMesh[i];
			this->fGraphMesh[i] = NULL;
		}
	}
}

void TPZAnalysis::PostProcess(int resolution) {
	int dim;
	for(dim=1; dim<=3; dim++) {
		PostProcess(resolution, dim);
	}
}

/** @brief Fill mat ids with materials with dimension dim wich are not boundary conditions or interface  */
void TPZAnalysis::IdentifyPostProcessingMatIds(int dimension, std::set<int> & matids){
    std::map<int, TPZMaterial * >::iterator matit;
    for (matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++) {
        
        TPZMaterial *mat = matit->second;
        auto *lag = dynamic_cast<TPZLagrangeMultiplierBase *> (mat);
        TPZBndCond *bc = dynamic_cast<TPZBndCond *> (mat);
        if (mat && !lag && mat->Dimension() == dimension){
            matids.insert(mat->Id());
        }
    }
    if (matids.size() == 0) {
        std::cout << __PRETTY_FUNCTION__ << " No materials were identified to perform post-processing. \n";
        DebugStop();
        return;
    }
}

void TPZAnalysis::PostProcess(int resolution, int dimension){
	int dim1 = dimension-1;
	if(!fGraphMesh[dim1]) return;

	fGraphMesh[dim1]->SetCompMesh(fCompMesh,fGraphMesh[dim1]->MaterialIds());
	fGraphMesh[dim1]->SetResolution(resolution);
	fGraphMesh[dim1]->DrawMesh(1);
	fGraphMesh[dim1]->DrawSolution(fStep,0);
	fStep++;
}

int TPZAnalysis::HighestDimension(){
	int dim = 0;
	std::map<int, TPZMaterial * >::iterator matit;
	for(matit = fCompMesh->MaterialVec().begin(); matit != fCompMesh->MaterialVec().end(); matit++)
	{
		if(!matit->second) continue;
		dim = matit->second->Dimension() > dim ? matit->second->Dimension() : dim;
	}
	return dim;
}

TPZAnalysis::TTablePostProcess::TTablePostProcess() :
fGeoElId(), fCompElPtr(), fDimension(-1), fLocations(), fVariableNames() {
}

int TPZAnalysis::TTablePostProcess::ClassId() const{
    return Hash("TPZAnalysis::TTablePostProcess");
}

void TPZAnalysis::TTablePostProcess::Write(TPZStream &buf, int withclassid) const {
    buf.Write(fGeoElId);
    buf.WritePointers(fCompElPtr);
    buf.Write(&fDimension);
    buf.Write(fLocations);
    buf.Write(fVariableNames);
}

void TPZAnalysis::TTablePostProcess::Read(TPZStream &buf, void *context){
    buf.Read(fGeoElId);
    buf.ReadPointers(fCompElPtr);
    buf.Read(&fDimension);
    buf.Read(fLocations);
    buf.Read(fVariableNames);
}


TPZAnalysis::TTablePostProcess::~TTablePostProcess() {
	fDimension = -1;
}

void TPZAnalysis::DefineElementTable(int dimension, TPZVec<int64_t> &GeoElIds, TPZVec<REAL> &points) {
	fTable.fDimension = dimension;
	fTable.fGeoElId = GeoElIds;
	fTable.fLocations = points;
	int64_t numel = GeoElIds.NElements();
	fTable.fCompElPtr.Resize(numel);
	int64_t iel;
	for(iel=0; iel<numel; iel++) {
		TPZGeoEl *gel = (TPZGeoEl *) fGeoMesh->FindElement(GeoElIds[iel]);
		fTable.fCompElPtr[iel] = (gel) ? gel->Reference() : 0;
	}
}

void TPZAnalysis::SetTableVariableNames(TPZVec<std::string> varnames) {
    fTable.fVariableNames = varnames;
}

void TPZAnalysis::PrePostProcessTable(std::ostream &out_file){
	TPZCompEl *cel;
	int numvar = fTable.fVariableNames.NElements();
	for(int iv=0; iv<numvar; iv++) {
		int64_t numel = fTable.fCompElPtr.NElements();
		for(int64_t iel=0; iel<numel; iel++) {
			cel = (TPZCompEl *) fTable.fCompElPtr[iel];
			if(cel) cel->PrintTitle(fTable.fVariableNames[iv].c_str(),out_file);
		}
	}
	out_file << endl;
	int dim;
	TPZVec<REAL> point(fTable.fDimension);
	for(dim=1; dim<fTable.fDimension+1; dim++) {
		for(int iv=0; iv<numvar; iv++) {
			int64_t numel = fTable.fCompElPtr.NElements();
			for(int64_t iel=0; iel<numel; iel++) {
				int d;
				for(d=0; d<fTable.fDimension; d++) {
					point[d] = fTable.fLocations[iel*fTable.fDimension+d];
				}
				cel = (TPZCompEl *) fTable.fCompElPtr[iel];
				if(cel) cel->PrintCoordinate(point,dim,out_file);
			}
		}
		out_file << endl;
	}
}

void TPZAnalysis::PostProcessTable(std::ostream &out_file) {
	TPZVec<REAL> point(fTable.fDimension);
	int numvar = fTable.fVariableNames.NElements();
	TPZCompEl *cel;
	for(int iv=0; iv<numvar; iv++) {
		int64_t numel = fTable.fCompElPtr.NElements();
		for(int64_t iel=0; iel<numel; iel++) {
			int d;
			for(d=0; d<fTable.fDimension; d++) {
				point[d] = fTable.fLocations[iel*fTable.fDimension+d];
			}
			cel = (TPZCompEl *) fTable.fCompElPtr[iel];
			if(cel) cel->PrintSolution(point,fTable.fVariableNames[iv].c_str(),out_file);
		}
	}
	out_file << endl;
}

template<class TVar>
TPZMatrixSolver<TVar> *TPZAnalysis::BuildPreconditioner(Precond::Type preconditioner,
                                                        bool overlap)
{
  TPZSimpleTimer build("BuildPreconditioner");
  auto mySolver = dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
	if(!mySolver || !mySolver->Matrix())
	{
#ifndef BORLAND
		cout << __FUNCTION__ << " called with uninitialized stiffness matrix\n";
#else
		cout << "TPZMatrixSolver *TPZAnalysis::BuildPreconditioner" << " called with uninitialized stiffness matrix\n";
#endif
		
	}
  /*
    First we treat two special cases in which the precond
    is simpler to create
   */
	if(preconditioner == Precond::Jacobi){
    return new TPZJacobiPrecond<TVar>(mySolver->Matrix());
	}
  else if(preconditioner == Precond::BlockJacobi && overlap){
    TPZSimpleTimer build("BuildPreconditioner::Create");
    TPZBlockDiagonalStructMatrix<TVar> blstr(fCompMesh);
    //now we check if equation filter is active
    if(this->StructMatrix()->EquationFilter().IsActive()){
      const int64_t nTotalEq = fCompMesh->NEquations();
      auto active_eqs = fStructMatrix->EquationFilter().GetActiveEquations();
      blstr.EquationFilter().SetActiveEquations(active_eqs);
    }
    TPZVec<int> blocksizes;
    blstr.BlockSizes(blocksizes);
    TPZBlockDiagonal<TVar> *sp = new TPZBlockDiagonal<TVar>();
    sp->Initialize(blocksizes);
    TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(sp);
    step->SetDirect(ELU);
    step->SetReferenceMatrix(mySolver->Matrix());
    return step;
  }
  else{
		TPZNodesetCompute nodeset;
		TPZStack<int64_t> elementgraph,elementgraphindex;
		int64_t nindep = fCompMesh->NIndependentConnects();
		int64_t neq = fCompMesh->NEquations();
    /*
      note: this graph will have even elements with no active connects.
      if such an element is at position i, then elemengraphindex[i] = elementgraphindex[i+1]
     */
		fCompMesh->ComputeElGraph(elementgraph,elementgraphindex);
		int64_t nel = elementgraphindex.NElements()-1;
		TPZRenumbering renum(nel,nindep);
        //this call will generate, based on the element graph, a graph that illustrates
        //connectivity between nodes
		renum.ConvertGraph(elementgraph,elementgraphindex,nodeset.Nodegraph(),nodeset.Nodegraphindex());
		nodeset.AnalyseGraph();

		TPZStack<int64_t> blockgraph,blockgraphindex;
		switch(preconditioner)
		{
    case Precond::Jacobi:
        DebugStop();
    case Precond::BlockJacobi:
				nodeset.BuildNodeGraph(blockgraph,blockgraphindex);
				break;
    case Precond::Element:
                blockgraph = elementgraph;
                blockgraphindex = elementgraphindex;
                //@orlandini: BuildElementGraph is broken. elgraph should temporarily suffice
				//nodeset.BuildElementGraph(blockgraph,blockgraphindex);
				break;
    case Precond::NodeCentered:
				nodeset.BuildVertexGraph(blockgraph,blockgraphindex);
				break;
		}

    
		TPZStack<int64_t> expblockgraph,expblockgraphindex;
		/*
      blockgraph and blockgraph indexes refer to CONNECTS (or equation blocks)
      now we expand them to equation numbers
     */
        {
            TPZManVector<int64_t> removed_blocks;
            nodeset.ExpandGraph(blockgraph,blockgraphindex,fCompMesh->Block(),
                                expblockgraph,expblockgraphindex,removed_blocks);
            TPZNodesetCompute::UpdateGraph(blockgraph,blockgraphindex,removed_blocks);
        }
    
#ifdef PZ_LOG
#ifdef PZDEBUG2
        if (logger.isDebugEnabled())
        {
            std::map<int64_t,int64_t> blocksizes;
            int64_t i;
            int64_t totalsize = 0;
            for(i=0; i< expblockgraphindex.NElements()-1;i++)
            {
                int64_t bls = expblockgraphindex[i+1]-expblockgraphindex[i];
                blocksizes[bls]++;
                totalsize += bls*bls;
            }
            std::map<int64_t,int64_t>::iterator it;
            std::stringstream sout;
            sout << __PRETTY_FUNCTION__ << " total size of allocation " << totalsize << std::endl;
            for(it=blocksizes.begin(); it != blocksizes.end(); it++)
            {
                sout << "block size " << (*it).first << " number of blocks " << (*it).second << std::endl;
            }
            LOGPZ_DEBUG(logger,sout.str().c_str());
        }
#endif
#endif
    if(this->StructMatrix()->EquationFilter().IsActive()){
      //now we need to convert the graphs as we have a filtered system,
      //therefore the eq nums have changed
      TPZManVector<int64_t> removed_blocks;
      TPZNodesetCompute::FilterGraph(this->StructMatrix()->EquationFilter(),
                                     expblockgraph,expblockgraphindex, removed_blocks);
      TPZNodesetCompute::UpdateGraph(blockgraph, blockgraphindex, removed_blocks);
      neq = fStructMatrix->NReducedEquations();
    }
      
    std::cout<<"Building "<<expblockgraphindex.size()-1<<" blocks"<<std::endl;
    if(overlap)
		{
#ifdef PZ_LOG
      if(loggerPrecond.isDebugEnabled()){
        std::stringstream sout;
        sout << "Precond type :" <<Precond::Name(preconditioner)<<'\n';
        sout << "Overlap : true\n";
        TPZNodesetCompute::Print(sout,blockgraphindex,blockgraph);
        LOGPZ_DEBUG(loggerPrecond,sout.str());
      }
#endif
      TPZSimpleTimer build("BuildPreconditioner::Create");
      TPZSparseBlockDiagonal<TVar> *sp = new TPZSparseBlockDiagonal<TVar>(expblockgraph,expblockgraphindex,neq);
			TPZStepSolver<TVar> *step = new TPZStepSolver<TVar>(sp);
			step->SetDirect(ELU);
      //this will allow the TPZAnalysis::Solve to call UpdateFrom and insert values
			step->SetReferenceMatrix(mySolver->Matrix());
			return step;
		}
    else
		{
      TPZSimpleTimer build("BuildPreconditioner::Create");
      TPZVec<int> blockcolor;
      //we must use the number of INDEPENDENT connects!
      int numcolors = nodeset.ColorGraph(blockgraph,blockgraphindex,nindep,blockcolor);
#ifdef PZ_LOG
      if(loggerPrecond.isDebugEnabled()){
        std::stringstream sout;
        sout << "Precond type :" <<Precond::Name(preconditioner)<<'\n';
        sout << "Overlap : false\n";
        sout << "Number of precond colors :"<<numcolors<<'\n';
        TPZNodesetCompute::Print(sout,blockgraphindex,blockgraph,blockcolor);
        LOGPZ_DEBUG(loggerPrecond,sout.str());
      }else if(loggerPrecond.isInfoEnabled()){
        std::stringstream sout;
        sout << "Precond type :" <<Precond::Name(preconditioner)<<'\n';
        sout << "Number of precond colors :"<<numcolors<<'\n';
        LOGPZ_INFO(loggerPrecond,sout.str());
      }
#endif
      std::cout<<"Blocks divided into "<<numcolors<<" colors"<<std::endl;
			return BuildSequenceSolver<TVar>(expblockgraph,expblockgraphindex,neq,numcolors,blockcolor);
		}
	}
	return 0;
}

/** @brief Build a sequence solver based on the block graph and its colors */
template<class TVar>
TPZMatrixSolver<TVar> *TPZAnalysis::BuildSequenceSolver(const TPZVec<int64_t> &graph,
                                                        const TPZVec<int64_t> &graphindex,
                                                        const int64_t neq, const int numcolors,
                                                        const TPZVec<int> &colors)
{
	TPZVec<TPZMatrix<TVar> *> blmat(numcolors);
	TPZVec<TPZStepSolver<TVar> *> steps(numcolors);
	int c;
    auto mySolver =
        dynamic_cast<TPZMatrixSolver<TVar>*>(fSolver);
	for(c=0; c<numcolors; c++)
	{
		blmat[c] = new TPZSparseBlockDiagonal<TVar>(graph,graphindex, neq, c, colors);
		steps[c] = new TPZStepSolver<TVar>(blmat[c]);
		steps[c]->SetDirect(ELU);
		steps[c]->SetReferenceMatrix(mySolver->Matrix());
	}
	if(numcolors == 1) return steps[0];
	TPZSequenceSolver<TVar> *result = new TPZSequenceSolver<TVar>;
  result->SetMatrix(mySolver->Matrix());
  for(c=0; c<numcolors; c++)
	{
		result->AppendSolver(*steps[c]);
	}
  for(c=numcolors-2; c>=0; c--)
	{
    /*
      append solver will copy the solver
      steps[c] has already been added to the
      sequence solver
      there is no need to update the matrix
      twice

      setting SetReferenceMatrix to nullptr
      will prevent this duplicate update
     */
		steps[c]->SetReferenceMatrix(nullptr);
		result->AppendSolver(*steps[c]);
	}
	for(c=0; c<numcolors; c++)
	{
		delete steps[c];
	}
	return result;
}

/// Integrate the postprocessed variable name over the elements included in the set matids
TPZVec<STATE> TPZAnalysis::Integrate(const std::string &varname, const std::set<int> &matids)
{
    // the postprocessed index of the varname for each material id
    std::map<int,int> variableids;
    int nvars = 0;
    
    std::map<int,TPZMaterial *>::iterator itmap;
    for (itmap = fCompMesh->MaterialVec().begin(); itmap != fCompMesh->MaterialVec().end(); itmap++) {
        if (matids.find(itmap->first) != matids.end()) {
            variableids[itmap->first] = itmap->second->VariableIndex(varname);
            nvars = itmap->second->NSolutionVariables(variableids[itmap->first]);
        }
    }
    TPZManVector<STATE,3> result(nvars,0.);
    int64_t nelem = fCompMesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = fCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            continue;
        }
        int matid = gel->MaterialId();
        if (matids.find(matid) == matids.end()) {
            continue;
        }
        TPZManVector<STATE,3> locres(nvars,0.);
        locres = cel->IntegrateSolution(variableids[matid]);
        for (int iv=0; iv<nvars; iv++)
        {
            result[iv] += locres[iv];
        }
    }
    return result;
}

/// extract the values corresponding to the connect from the vector
static void ConnectSolution(int64_t cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob, TPZVec<STATE> &sol)
{
    int64_t seqnum = cmesh->ConnectVec()[cindex].SequenceNumber();
    int blsize = cmesh->Block().Size(seqnum);
    int position = cmesh->Block().Position(seqnum);
    sol.resize(blsize);
    for (int64_t i=position; i< position+blsize; i++) {
        sol[i-position] = glob(i,0);
    }
}

static STATE ConnectNorm(int64_t cindex, TPZCompMesh *cmesh, TPZFMatrix<STATE> &glob)
{
    TPZManVector<STATE,20> cvec;
    ConnectSolution(cindex, cmesh, glob, cvec);
    STATE norm = 0.;
    for (int64_t i=0; i<cvec.size(); i++) {
        norm += cvec[i]*cvec[i];
    }
    norm = sqrt(norm);
    return norm;
}

/// Print the residual vector for those elements with entry above a given tolerance
void TPZAnalysis::PrintVectorByElement(std::ostream &out, TPZFMatrix<STATE> &vec, REAL tol)
{
    int64_t nel = fCompMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCompMesh->Element(el);
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        int ic;
        for (ic=0; ic<nc; ic++) {
            int64_t cindex = cel->ConnectIndex(ic);
            TPZConnect &c = cel->Connect(ic);
            if(c.HasDependency() || c.IsCondensed())
            {
                continue;
            }
#ifndef STATE_COMPLEX
            if (ConnectNorm(cindex, fCompMesh, vec) >= tol) {
                break;
            }
#endif
        }
        // if all connects have norm below tol, do not print the element
        if (ic == nc) {
            continue;
        }
        bool hasgeometry = false;
        TPZManVector<REAL> xcenter(3,0.);
        TPZGeoEl *gel = cel->Reference();
        // if a geometric element exists, extract its center coordinate
        if (gel) {
            hasgeometry = true;
            TPZManVector<REAL,3> xicenter(gel->Dimension(),0.);
            gel->CenterPoint(gel->NSides()-1, xicenter);
            gel->X(xicenter,xcenter);
        }
        out << "CompEl " << el;
        if (hasgeometry) {
            out << " Gel " << gel->Index() << " matid " << gel->MaterialId() << " Center " << xcenter << std::endl;
        }
        for (ic = 0; ic<nc; ic++) {
            TPZManVector<STATE> connectsol;
            int64_t cindex = cel->ConnectIndex(ic);
            TPZConnect &c = fCompMesh->ConnectVec()[cindex];
            if(c.HasDependency() || c.IsCondensed()) {
                out << "connect " << ic << " is restrained\n";
                continue;
            }
            ConnectSolution(cindex, fCompMesh, vec, connectsol);
            for (int i=0; i<connectsol.size(); i++) {
                if (fabs(connectsol[i]) < tol) {
                    connectsol[i] = 0.;
                }
            }
            out << ic << " index " << cindex << " values " << connectsol << std::endl;
        }
    }

}

int TPZAnalysis::ClassId() const{
    return Hash("TPZAnalysis");
}

void TPZAnalysis::Write(TPZStream &buf, int withclassid) const{
    TPZPersistenceManager::WritePointer(fGeoMesh, &buf);
    TPZPersistenceManager::WritePointer(fCompMesh, &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[0], &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[1], &buf);
    TPZPersistenceManager::WritePointer(fGraphMesh[2], &buf);
    buf.Write((int)fSolType);
    fSolution.Write(buf,withclassid);
    TPZPersistenceManager::WritePointer(fSolver, &buf);
    buf.Write(fScalarNames[0]);
    buf.Write(fScalarNames[1]);
    buf.Write(fScalarNames[2]);
    buf.Write(fVectorNames[0]);
    buf.Write(fVectorNames[1]);
    buf.Write(fVectorNames[2]);
    buf.Write(&fStep);
    buf.Write(&fNthreadsError);
    TPZPersistenceManager::WritePointer(fStructMatrix.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fRenumber.operator ->(), &buf);
    fTable.Write(buf,withclassid);
    buf.Write(fTensorNames[0]);
    buf.Write(fTensorNames[1]);
    buf.Write(fTensorNames[2]);
}

void TPZAnalysis::Read(TPZStream &buf, void *context){
    fGeoMesh = dynamic_cast<TPZGeoMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fCompMesh = dynamic_cast<TPZCompMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[0] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[1] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fGraphMesh[2] = dynamic_cast<TPZGraphMesh*>(TPZPersistenceManager::GetInstance(&buf));
    fSolType = [&buf](){
        int tmp;
        buf.Read(&tmp);
        return (ESolType) tmp;
    }();
    fSolution.Read(buf,context);
    fSolver = dynamic_cast<TPZMatrixSolver<STATE>*>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(fScalarNames[0]);
    buf.Read(fScalarNames[1]);
    buf.Read(fScalarNames[2]);
    buf.Read(fVectorNames[0]);
    buf.Read(fVectorNames[1]);
    buf.Read(fVectorNames[2]);
    buf.Read(&fStep);
    buf.Read(&fNthreadsError);
    fStructMatrix = TPZAutoPointerDynamicCast<TPZStructMatrix>(TPZPersistenceManager::GetAutoPointer(&buf));
    fRenumber = TPZAutoPointerDynamicCast<TPZRenumbering>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTable.Read(buf,context);
    buf.Read(fTensorNames[0]);
    buf.Read(fTensorNames[1]);
    buf.Read(fTensorNames[2]);
}

TPZAnalysis::ThreadData::ThreadData(TPZAdmChunkVector<TPZCompEl *> &elvec, bool store_error) : fNextElement(0), fvalues(0), fStoreError(store_error), ftid(0), fElvec(elvec){
}


TPZAnalysis::ThreadData::~ThreadData()
{
}

template
TPZMatrixSolver<STATE> *TPZAnalysis::BuildPreconditioner<STATE>(
    Precond::Type preconditioner,bool overlap);
template
TPZMatrixSolver<CSTATE> *TPZAnalysis::BuildPreconditioner<CSTATE>(
    Precond::Type preconditioner,bool overlap);
