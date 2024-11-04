/**
 * @file
 * @brief Contains the implementation of the TPZCompMesh methods.
 */

#include "pzcmesh.h"
#ifdef MACOSX
//#include <__functional_base>               // for less
#include <__tree>                          // for __tree_const_iterator, ope...
#endif
#include <cmath>                           // for fabs, sqrt, abs
#include <iterator>                        // for operator!=, reverse_iterator
#include <map>                             // for map, __map_iterator, opera...
#include <set>                             // for set, set<>::reverse_iterator
#include <string>                          // for char_traits, allocator
#include <utility>                         // for pair
#include "TPZCompElDisc.h"                 // for TPZCompElDisc
#include "TPZInterfaceEl.h"                // for TPZInterfaceElement
#include "pzadmchunk.h"                    // for TPZAdmChunkVector
#include "pzblock.h"                       // for TPZBlock
#include "TPZBndCond.h"                     // for TPZBndCond
#include "pzcompel.h"                      // for TPZCompEl, TPZCompElSide
#include "pzcondensedcompel.h"             // for TPZCondensedCompEl
#include "pzconnect.h"                     // for TPZConnect
#include "pzelementgroup.h"                // for TPZElementGroup
#include "pzeltype.h"                      // for MElementType::EAgglomerate
#include "pzerror.h"                       // for PZError, DebugStop
#include "TPZStream.h"                     // for TPZStream
#include "pzgeoel.h"                       // for TPZGeoEl
#include "pzgeoelside.h"                   // for TPZGeoElSide
#include "pzgmesh.h"                       // for TPZGeoMesh
#include "pzgnode.h"                       // for TPZGeoNode
#include "pzintel.h"                       // for TPZInterpolatedElement
#include "pzinterpolationspace.h"          // for TPZInterpolationSpace
#include "pzlog.h"                         // for glogmutex, LOGPZ_DEBUG
#include "pzmanvector.h"                   // for TPZManVector
#include "TPZMaterial.h"                    // for TPZMaterial
#include "TPZMaterialDataT.h"                // for TPZSolVec
#include "pzmatrix.h"                      // for TPZFMatrix, TPZMatrix
#include "pzmultiphysicselement.h"         // for TPZMultiphysicsElement
#include "pzsubcmesh.h"                    // for TPZSubCompMesh
#include "pztransfer.h"                    // for TPZTransfer
#include "pztrnsform.h"                    // for TPZTransform
#include "pzvec.h"                         // for TPZVec, operator<<
#include "TPZMatError.h"

#ifndef STATE_COMPLEX
	#include "TPZAgglomerateEl.h" // for TPZAgglomerateElement
#endif

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzcompmesh");
static TPZLogger aloclogger("pz.allocate");
#endif
using namespace std;

TPZCompMesh::TPZCompMesh () :
    TPZRegisterClassId(&TPZCompMesh::ClassId),
    fElementVec(0),fConnectVec(0),fMaterialVec(),
    fSolType(EUndefined),
    fSolution(),fElementSolution(),fSolN(){
    
#ifdef PZ_LOG
    if (aloclogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Allocate TPZCompMesh this = " << (void *)this;
        LOGPZ_DEBUG(aloclogger, sout.str())
    }
#endif
	fDefaultOrder = TPZCompEl::GetgOrder();
	
	//Initializing class members
	fDimModel = 0;
	fReference = nullptr;
	//  fChecked = 0;
	//fName[0] = '\0';
	//fName[126] = '\0';
  SetName( "Computational mesh");
}


TPZCompMesh::TPZCompMesh (TPZGeoMesh* gr, bool isComplex) :
    TPZRegisterClassId(&TPZCompMesh::ClassId),
    fElementVec(0),fConnectVec(0),fMaterialVec(),
    fSolType(isComplex ? EComplex : EReal),
    fSolution(0,1,isComplex),
    fElementSolution(0,1,isComplex),
    fSolN(0,1,isComplex){
    
#ifdef PZ_LOG
    if (aloclogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Allocate TPZCompMesh this = " << (void *)this;
        LOGPZ_DEBUG(aloclogger, sout.str())
    }
#endif
	fDefaultOrder = TPZCompEl::GetgOrder();
	
	//Initializing class members
	fDimModel = 0;
	fReference = gr;
	//  fChecked = 0;
	//fName[0] = '\0';
	//fName[126] = '\0';
	if(gr) {
		SetName( gr->Name() );
		gr->ResetReference();
		gr->SetReference(this);
        SetDimModel(gr->Dimension());
	}
    else {
        SetName( "Computational mesh");
    }
    TPZBaseMatrix &sol = fSolution;
	fBlock.SetMatrix(&sol);
	fSolutionBlock.SetMatrix(&sol);
    
    fNmeshes = 0;
}


TPZCompMesh::TPZCompMesh(TPZAutoPointer<TPZGeoMesh> &gmesh,
                         bool isComplex) :
    TPZRegisterClassId(&TPZCompMesh::ClassId),
    fGMesh(gmesh),fElementVec(0),fConnectVec(0),fMaterialVec(),
    fSolType(isComplex ? EComplex : EReal),
    fSolution(0,1,isComplex),
    fElementSolution(0,1,isComplex),
    fSolN(0,1,isComplex)
{
#ifdef PZ_LOG
    if (aloclogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Allocate TPZCompMesh this = " << (void *)this;
        LOGPZ_DEBUG(aloclogger, sout.str())
    }
#endif

    fDefaultOrder = TPZCompEl::GetgOrder();
    
    //Initializing class members
    fDimModel = gmesh->Dimension();
    fReference = gmesh.operator->();
    //  fChecked = 0;
    //fName[0] = '\0';
    //fName[126] = '\0';
    if(fReference) {
        SetName( fReference->Name() );
        fReference->ResetReference();
        fReference->SetReference(this);
    }
    else {
        SetName( "Computational mesh");
    }
    TPZBaseMatrix &sol = fSolution;
    fBlock.SetMatrix(&sol);
    fSolutionBlock.SetMatrix(&sol);
    
    fNmeshes = 0;
}

TPZCompMesh::~TPZCompMesh() {
	
#ifdef PZ_LOG
    if (aloclogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Delete TPZCompMesh this = " << (void *)this;
        LOGPZ_DEBUG(aloclogger, sout.str())
    }
#endif

#ifdef PZ_LOG2
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	// THIS NEEDS TO INCLUDE THE DELETION ROUTINES OF ALL ITEMS
	this->CleanUp();
	TPZGeoMesh * ref = this->Reference();
	if (ref){
		if(ref->Reference() == this){
			ref->ResetReference();
		}//if == this
	}//if (ref)
}//~

void TPZCompMesh::CleanUp() {
	
	// THIS ROUTINE NEEDS TO INCLUDE THE DELETION OF THE LIST POINTERS
	TPZGeoMesh *ref = Reference();
//	if (ref){
//		ref->ResetReference();
//		this->LoadReferences();
//	}

	int64_t i, nelem = this->NElements();

	//deleting subcompmesh
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZSubCompMesh * subc = dynamic_cast<TPZSubCompMesh*>(el);
		if(subc){
			delete subc;
		}
	}
	
    
	//unwrapping condensed compel
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(el);
		if(cond){
			cond->Unwrap();
		}
	}
	
	//unwrapping element groups
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZElementGroup * group = dynamic_cast<TPZElementGroup*>(el);
		if(group){
            group->Unwrap();
		}
	}
	
	//deleting interface elements
	for(i=0; i<nelem; i++){
		TPZCompEl *el = fElementVec[i];
		TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(el);
		if(face){
            delete el;
		}
	}
	
	//deleting other elements
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		delete el;
	}
	
	fElementVec.Resize(0);
	fElementVec.CompactDataStructure(fElementVec.NOW);
	fConnectVec.Resize(0);
	fConnectVec.CompactDataStructure(fConnectVec.NOW);
	nelem = NMaterials();
    std::map<int, TPZMaterial *>::iterator it;
    for (it = fMaterialVec.begin(); it != fMaterialVec.end(); it++) {
        delete it->second;
    }
	fMaterialVec.clear();
	
	fBlock.SetNBlocks(0);
	fSolutionBlock.SetNBlocks(0);
}

void TPZCompMesh::SetName (const string &nm) {
	fName = nm;
}

void TPZCompMesh::Print (std::ostream & out) const {
	
	//ComputeNodElCon();
	out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
	out << "TITLE-> " << fName << "\n\n";
	
	out << "number of connects            = " << NConnects() << std::endl;
	out << "number of elements            = " << NElements() << std::endl;
	out << "number of materials           = " << NMaterials() << std::endl;
	out << "dimension of the mesh         = " << this->Dimension() << std::endl;
	
	out << "\n\t Connect Information:\n\n";
	int64_t i, nelem = NConnects();
	for(i=0; i<nelem; i++) {
		if(fConnectVec[i].SequenceNumber() == -1) {
			if(fConnectVec[i].HasDependency()) {
				cout << " TPZCompMesh::Print inconsistency of connect\n";
				cout << " Index " << i << ' ';
				fConnectVec[i].Print(*this,std::cout);
			}
			continue;
		}
		out << " Index " << i << ' ';
		fConnectVec[i].Print(*this,out);
	}
	out << "\n\t Computable Element Information:\n\n";
	nelem = NElements();
	for(i=0; i<nelem; i++) {
		if(!fElementVec[i]) continue;
		TPZCompEl *el = fElementVec[i];
		out << "\n Index " << i << ' ';
		el->Print(out);
        TPZMultiphysicsElement *mpel = dynamic_cast<TPZMultiphysicsElement *>(el);
        if(!mpel){
            if(!el->Reference()) continue;
            out << "\tReference Index = " << el->Reference()->Index() << std::endl << std::endl;
        }
	}
	out << "\n\t Material Information:\n\n";
	std::map<int, TPZMaterial * >::const_iterator mit;
	nelem = NMaterials();
	for(mit=fMaterialVec.begin(); mit!= fMaterialVec.end(); mit++) {
        TPZMaterial *mat = mit->second;
        if (!mat) {
          DebugStop();
        }
        int matindex = mit->first;
        out << "material index " << matindex << std::endl;
		mat->Print(out);
        out << std::endl;
	}
}

void TPZCompMesh::ShortPrint(std::ostream &out) const {
    //ComputeNodElCon();
    out << "\n\t\tCOMPUTABLE GRID INFORMATIONS:\n\n";
    out << "TITLE-> " << fName << "\n\n";

    out << "number of connects            = " << NConnects() << std::endl;
    out << "number of elements            = " << NElements() << std::endl;
    out << "number of materials           = " << NMaterials() << std::endl;
    out << "dimension of the mesh         = " << this->Dimension() << std::endl;

    out << "\n\t Connect Information:\n\n";
    int64_t i, nelem = NConnects();
    for(i=0; i<nelem; i++) {
        if(fConnectVec[i].SequenceNumber() == -1) {
            if(fConnectVec[i].HasDependency()) {
                cout << " TPZCompMesh::Print inconsistency of connect\n";
                cout << " Index " << i << ' ';
                fConnectVec[i].Print(*this,std::cout);
            }
            continue;
        }
        out << " Index " << i << ' ';
        fConnectVec[i].Print(*this,out);
        out << "\n";
    }
    out << "\n\t Computable Element Information:\n\n";
    nelem = NElements();
    for(i=0; i<nelem; i++) {
        if(!fElementVec[i]) continue;
        TPZCompEl *el = fElementVec[i];
        el->ShortPrint(out);
        TPZMultiphysicsElement *mpel = dynamic_cast<TPZMultiphysicsElement *>(el);
        if(!mpel){
            if(!el->Reference()) continue;
            out << "\tReference Index = " << el->Reference()->Index() << std::endl << std::endl;
        }
    }
    out << "\n\t Material Information:\n\n";
    std::map<int, TPZMaterial * >::const_iterator mit;
    nelem = NMaterials();
    for(mit=fMaterialVec.begin(); mit!= fMaterialVec.end(); mit++) {
        TPZMaterial *mat = mit->second;
        if (!mat) {
            DebugStop();
        }
        mat->Print(out);
    }
}

/**Insert a material object in the datastructure*/
int TPZCompMesh::InsertMaterialObject(TPZMaterial * mat) {
	if(!mat) return -1;
	int matid = mat->Id();
    if (fMaterialVec.find(matid) != fMaterialVec.end()) {
        DebugStop();
    }
	fMaterialVec[matid] = mat;
	return fMaterialVec.size();
}
/*Due to the TPZMaterialRefactor, the type TPZBndCond
  no longer inherits from TPZMaterial. But every possible instance
  will be a TPZMaterial, so this dynamic_cast is not expected
  to fail
*/
int TPZCompMesh::InsertMaterialObject(TPZBndCond * mat) {
    return InsertMaterialObject(dynamic_cast<TPZMaterial*>(mat));
}

TPZMaterial * TPZCompMesh::FindMaterial(int matid){	// find the material object with id matid
	TPZMaterial * result = 0;
	std::map<int, TPZMaterial * >::iterator mit;
	mit = fMaterialVec.find(matid);
	if(mit != fMaterialVec.end())
	{
		result = mit->second;
	}
	return result;
}

void TPZCompMesh::AutoBuild(const std::set<int> *MaterialIDs) {
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = Reference();
        TPZCheckGeom check(gmesh);
        check.CheckIds();
    }
#endif
    if(MaterialIDs)
    {
        fCreate.BuildMesh(*this, *MaterialIDs);        
    }
    else {
        fCreate.BuildMesh(*this);
    }
	
}

void TPZCompMesh::AutoBuildContDisc(const TPZVec<TPZGeoEl*> &continuous, const TPZVec<TPZGeoEl*> &discontinuous) {
#ifdef PZDEBUG
    {
        TPZGeoMesh *gmesh = Reference();
        TPZCheckGeom check(gmesh);
        check.CheckIds();
    }
#endif
	
	TPZAdmChunkVector<TPZGeoEl *> &elvec = Reference()->ElementVec();
	int64_t nelem = elvec.NElements();
	
	int64_t neltocreate = 0;
	for(int64_t i=0; i<nelem; i++) {
		TPZGeoEl *gel = elvec[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			neltocreate++;
		}
	}
	
	int64_t nbl = fBlock.NBlocks();
	if(neltocreate > nbl) fBlock.SetNBlocks(neltocreate);
	fBlock.SetNBlocks(nbl);
	
	//Creating continuous elements
	fCreate.SetAllCreateFunctionsContinuous();
	int64_t ncont = continuous.NElements();
	for(int64_t i = 0; i < ncont; i++){
		TPZGeoEl *gel = continuous[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int printing = 0;
			if (printing) {
				gel->Print(cout);
			}
			
			if(gel->NumInterfaces() == 0){
				CreateCompEl(gel);
			}
		}
	}
	
	//Creating discontinuous elements
	fCreate.SetAllCreateFunctionsDiscontinuous();
	int64_t ndisc = discontinuous.NElements();
	for(int64_t i = 0; i < ndisc; i++){
		TPZGeoEl *gel = discontinuous[i];
		if(!gel) continue;
		if(!gel->HasSubElement()) {
			int printing = 0;
			if (printing) {
				gel->Print(cout);
			}
			
			if(gel->NumInterfaces() == 0){
				CreateCompEl(gel);
			}
		}
	}
	
	this->InitializeBlock();
}

void TPZCompMesh::InitializeBlock() {
	ExpandSolution();
	CleanUpUnconnectedNodes();
}


void TPZCompMesh::ExpandSolution(){
    try{
        if(fSolType==EReal)
            ExpandSolutionInternal<STATE>(fSolution);
        else
            ExpandSolutionInternal<CSTATE>(fSolution);
    }
    catch(...){
        PZError << "Incompatible matrix type in ";
        PZError << __PRETTY_FUNCTION__;
        PZError << std::endl;
        DebugStop();
    }
}

template<class TVar>
void TPZCompMesh::ExpandSolutionInternal(TPZFMatrix<TVar> &sol) {
	fBlock.Resequence();
	int64_t ibl,nblocks = fBlock.NBlocks();
	
	TPZFMatrix<TVar> OldSolution(sol);
	
	int64_t cols = sol.Cols();
	sol.Redim(fBlock.Dim(),cols);
	int64_t minblocks = nblocks < fSolutionBlock.NBlocks() ? nblocks : fSolutionBlock.NBlocks();
	
	int64_t ic;
	for(ic=0; ic<cols; ic++) {
		for(ibl = 0;ibl<minblocks;ibl++) {
			int64_t oldsize = fSolutionBlock.Size(ibl);
			int64_t oldposition = fSolutionBlock.Position(ibl);
			int64_t newsize = fBlock.Size(ibl);
			int64_t newposition = fBlock.Position(ibl);
			int64_t minsize = (oldsize < newsize) ? oldsize : newsize;
			int64_t ieq;
			for(ieq=0; ieq<minsize; ieq++) {
				sol.PutVal(newposition+ieq,ic,OldSolution(oldposition+ieq,ic));
			}
		}
	}
	fSolutionBlock = fBlock;
}

void TPZCompMesh::LoadSolution(const TPZSolutionMatrix &mat){
    /*
      The TPZLinearAnalysis class will store the solution associated with the independent
      equations, i.e., the solution returned from the solver.
      Meanwhile, the TPZCompMesh class stores the *full* solution, therfore it needs
      extra room for the dependent dofs
    */
    fSolution.ExpandAndSetSol(mat, fSolution.Rows());
    
    const auto nelem = NElements();
	for(auto i=0; i<nelem; i++) {
		TPZCompEl *cel = fElementVec[i];
		if(!cel) continue;
		cel->LoadSolution();
	}
}

void TPZCompMesh::TransferMultiphysicsSolution()
{
    int64_t nel = this->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZCompEl *cel = this->Element(iel);
        if (!cel) {
            continue;
        }
        cel->TransferMultiphysicsElementSolution();
    }
}

void TPZCompMesh::LoadReferences() {
	
	//	Reference()->ResetReference();
	Reference()->SetReference(this);
	int64_t i, nelem = NElements();
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		el->LoadElementReference();
		/*  TPZGeoEl *gel = el->Reference();
		 if(!gel) continue;
		 gel->SetReference(el);
		 */
	}
}

void TPZCompMesh::CleanUpUnconnectedNodes() {
    // we assume the sequence numbers of the connects are distinct and that each
    // connect with a sequence number >= 0 has a corresponding block
    // the values of the solution vector will be permuted toghether with the sequence
    // numbers
	ComputeNodElCon();
	int64_t nconnects = NConnects();
    int64_t ndepblocks = 0, nvalidblocks = 0, nremoved = 0, ncondensed = 0, maxseq = -1, numnowithseq = 0;
    int64_t nblocks = fBlock.NBlocks();

	for (int64_t i=0;i<nconnects;i++)
    {
		TPZConnect &no = fConnectVec[i];
		int64_t seq = no.SequenceNumber();
        // a connect with a sequence number needs a representation in the fblock
        // data structure
        if(seq >= 0 && seq >= nblocks) DebugStop();
        if(seq > maxseq) maxseq = seq;
        if(seq >= 0) numnowithseq++;
        if(seq < 0 && (no.NElConnected() || no.IsCondensed()))
        {
            DebugStop();
        }
		if(!no.NElConnected() && !no.IsCondensed() && seq != -1)
		{
			nremoved++;
		}
		else if(!no.HasDependency() && no.NElConnected() && !no.IsCondensed()) nvalidblocks++;
        else if(!no.HasDependency() && no.IsCondensed()) ncondensed++;
		else if(no.HasDependency() && no.NElConnected()) ndepblocks++;
#ifdef PZDEBUG
        if(no.HasDependency() && no.IsCondensed())
        {
            std::cout << __PRETTY_FUNCTION__ << " I dont understand\n";
            DebugStop();
        }
#endif
    }
	int need = 0;
    if(maxseq != numnowithseq-1)
    {
        need = 1;
    }
	for (int64_t i=0;i<nconnects;i++) {
		TPZConnect &no = fConnectVec[i];
		if (no.SequenceNumber() == -1) continue;
		if (no.HasDependency() && no.NElConnected() == 0) {
			PZError << "TPZCompMesh::CleanUpUnconnectedNodes node has dependency\n";
			continue;
		}
		if (!no.NElConnected() && !no.IsCondensed() && no.SequenceNumber() != -1)
		{
			need = 1;
			break;
		}
		if (no.HasDependency() && no.SequenceNumber() < nvalidblocks)
		{
			need = 1;
			break;
		} else if(!no.HasDependency() && !no.IsCondensed() && no.SequenceNumber() >= nvalidblocks)
		{
			need = 1;
			break;
		} else if(!no.HasDependency() && no.IsCondensed() && no.SequenceNumber() < nvalidblocks)
		{
			need = 1;
			break;
		}
	}
    TPZManVector<int64_t> permute(nblocks,-1);
    // the down datastructure will be equal 1 for the sequence numbers that
    // have been assigned and zero for the others
    // the down data structure allows for sequence numbers with "holes"
    TPZManVector<int64_t> down(nblocks,0);
	int64_t idepblocks = 0, iremovedblocks= 0, icondensed = 0;
	
	if (need) {
		for(int64_t i=0; i<nconnects; i++) {
			TPZConnect &no = fConnectVec[i];
			if(no.SequenceNumber() == -1) continue;
			int seq = no.SequenceNumber();
            // if the value of permute is not -1, then there are connects with
            // duplicate sequence number
            if(permute[seq] != -1) DebugStop();
            // if the sequence number is larger than the number of blocks
            // the datastructure is inconsistent
            if(seq >= nblocks) DebugStop();
            // a condensed connect cannot be removed
			if(no.NElConnected() == 0 && !no.IsCondensed())
			{
				permute[seq] = nvalidblocks+ndepblocks+iremovedblocks+ncondensed;
				down[seq] = 1;
				fBlock.Set(seq,0);
                no.Reset();
				fConnectVec.SetFree(i);
				iremovedblocks++;
			}
			else if(no.HasDependency()) {
				permute[seq] = nvalidblocks+ncondensed+idepblocks;
				down[seq] = 1;
				idepblocks++;
			}
            else if(no.IsCondensed())
            {
				permute[seq] = nvalidblocks+icondensed;
				down[seq] = 1;
				icondensed++;
            }
		}
		for(int64_t i=1; i<nblocks; i++) down[i] += down[i-1];
		for(int64_t i=0; i<nblocks; i++)
		{
			if(permute[i] == -1)
			{
				permute[i] = i-down[i];
			}
		}
	}
#ifdef PZ_LOG
	if(need)
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "permute to put the free connects to the back\n";
            if(nblocks < 50)
            {
                sout << "original sequence numbers|nelconected\n";
                int64_t nel = fConnectVec.NElements();
                for (int64_t el=0; el<nel; el++) {
                    TPZConnect &c = fConnectVec[el];
                    int64_t seqnum = c.SequenceNumber();
                    sout << seqnum << '|' << c.NElConnected() << " ";
                }
                sout << std::endl;
            }
            if(nblocks < 50) {
                for (int64_t i=0;i<nblocks;i++) sout << permute[i] << ' ';
                sout << std::endl;
            }
            sout << "need = " << need << endl;
            LOGPZ_DEBUG(logger,sout.str());
        }
#endif
	
	if (need) {
#ifdef PZDEBUG
		std::set<int64_t> check;
		nconnects = permute.NElements();
		for(int64_t i=0; i<nconnects; i++)
        {
            if(permute[i] < 0 || permute[i] >= nconnects)
            {
                std::cout << "value of permute " << permute[i] << " is larger than " << nconnects << std::endl;
            }
            else
            {
                check.insert(permute[i]);
            }
        }
		if(static_cast<int>(check.size()) != nconnects)
		{
			cout << __PRETTY_FUNCTION__ << " The permutation vector is not a permutation!\n" << permute << endl;
			DebugStop();
		}
#endif
		Permute(permute);
        
#ifdef PZ_LOG
        if (logger.isDebugEnabled() && nblocks < 50)
        {
            if(nblocks < 50)
            {
                std::stringstream sout;
                sout << "after permute sequence numbers|nelconected\n";
                int64_t nel = fConnectVec.NElements();
                for (int64_t el=0; el<nel; el++) {
                    TPZConnect &c = fConnectVec[el];
                    int64_t seqnum = c.SequenceNumber();
                    sout << seqnum << '|' << c.NElConnected() << " ";
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
            
        }
#endif
        int64_t nel = fConnectVec.NElements();
        for (int64_t i=0;i<nel;i++) {
            TPZConnect &no = fConnectVec[i];
            if (no.NElConnected() == 0 && !no.IsCondensed() && no.SequenceNumber() >= nblocks-nremoved) {
                no.Reset();
                fConnectVec.SetFree(i);
            }
            else if(no.NElConnected() == 0 && !no.IsCondensed() && no.SequenceNumber() != -1)
            {
                DebugStop();
            }
        }
#ifdef PZDEBUG
        {
            int64_t nel = fConnectVec.NElements();
            for (int64_t el=0; el<nel; el++) {
                TPZConnect &c = fConnectVec[el];
                int64_t seqnum = c.SequenceNumber();
                if (seqnum > nblocks-nremoved) {
                    DebugStop();
                }
            }
        }
#endif

		fBlock.SetNBlocks(nblocks-nremoved);
	}
}

void TPZCompMesh::ComputeNodElCon() {
	
	int64_t i, nelem = NConnects();
	for(i=0; i<nelem; i++) {
		TPZConnect &no = fConnectVec[i];
		if(no.SequenceNumber() == -1) continue;
		no.ResetElConnected();
	}
	
	
	TPZStack<int64_t> nodelist;
	int64_t numnod;
	// modified Philippe 22/7/97
	// in order to account for constrained nodes
	nelem = NElements();
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		nodelist.Resize(0);
		el->BuildConnectList(nodelist);
		numnod = nodelist.NElements();
		for (int64_t in=0; in<numnod; ++in) {
			int64_t dfnindex = nodelist[in];
#ifdef PZDEBUG
            if(dfnindex < 0)
            {
                std::cout << "element index i " << i << " node in " << in <<
                " has negative node index " << dfnindex;
                DebugStop();
            }
#endif
			TPZConnect *dfn = &fConnectVec[dfnindex];
			dfn->IncrementElConnected();
		}
	}
}

void TPZCompMesh::ComputeNodElCon(TPZVec<int> &nelconnected ) const {
	
	int64_t i, nelem = NConnects();
	nelconnected.Resize(nelem);
	nelconnected.Fill(0);
	TPZStack<int64_t> nodelist;
	int64_t numnod;
	// modified Philippe 22/7/97
	// in order to account for constrained nodes
	nelem = NElements();
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		nodelist.Resize(0);
		el->BuildConnectList(nodelist);
		numnod = nodelist.NElements();
		for (int64_t in=0; in<numnod; ++in) {
			int64_t dfnindex = nodelist[in];
			nelconnected[dfnindex]++;
		}
	}
}


int64_t TPZCompMesh::NEquations() {
	int64_t neq = 0;
	int64_t i, ncon = NConnects();
	for(i=0; i<ncon; i++) {
		TPZConnect &df = fConnectVec[i];
        if(df.HasDependency() || df.IsCondensed() || !df.NElConnected() || df.SequenceNumber() == -1){
            continue;
        }
        
        int dofsize = df.NShape()*df.NState();
#ifdef PZDEBUG
        // check the consistency between the block size and the data structure of the connect
        {
            int64_t seqnum = df.SequenceNumber();
            int64_t blsize = fBlock.Size(seqnum);
            if (blsize != dofsize) {
                DebugStop();
            }
        }
#endif
        neq += dofsize;
	}
	return neq;
}

/** Este metodo nao e apropriado para aproximantes descontinuas */
int TPZCompMesh::BandWidth() {
	
	int bw = 0;
	TPZStack<int64_t> connectlist;
	// modified Philippe 24/7/97
	// in order to take dependent nodes into account
	
	int64_t i, nelem = NElements();
	for(i=0; i<nelem; i++) {
		connectlist.Resize(0);
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		el->BuildConnectList(connectlist);
		int64_t nnod = connectlist.NElements();
		if(!nnod) continue;
        // look for a node which has equations associated with it
		int64_t ifirstnode = 0;
		TPZConnect *np = &fConnectVec[connectlist[ifirstnode++]];
		while(ifirstnode < nnod && (np->HasDependency() || np->IsCondensed() || !fBlock.Size(np->SequenceNumber()))) {
			np = &fConnectVec[connectlist[ifirstnode++]];
		}
		int64_t ibl = np->SequenceNumber();
		int64_t loweq = fBlock.Position(ibl);
		int64_t higheq = loweq+fBlock.Size(ibl)-1;
		for(int64_t n=ifirstnode;n<nnod;n++) {
			np = &fConnectVec[connectlist[n]];
			if(np->HasDependency() || np->IsCondensed() ) continue;
			int64_t ibl = np->SequenceNumber();
			if(!fBlock.Size(ibl)) continue;
			int64_t leq = fBlock.Position(ibl);
			int64_t heq = leq+fBlock.Size(ibl)-1;
			loweq = (loweq > leq) ? leq : loweq;
			higheq = (higheq < heq) ? heq : higheq;
		}
		int64_t elbw = higheq - loweq;
		bw = (bw < elbw) ? elbw : bw;
	}
	return bw;
}

void TPZCompMesh::Skyline(TPZVec<int64_t> &skyline) {
	
	TPZStack<int64_t> connectlist;
	// modified Philippe 24/7/97
	// in order to take dependent nodes into account
	
#ifdef PZDEBUG
  TPZAdmChunkVector <TPZConnect > & connectVec = this->ConnectVec();
  int maxSequenceNumberIndependentConnect = 0;
  TPZVec<int> depConInd(0,0);//list of connects that have dependencies
  for (int i = 0; i < connectVec.NElements(); i++) {
    if (connectVec[i].HasDependency() ){
      int oldSize = depConInd.size();
      depConInd.Resize( oldSize + 1 );
      depConInd[oldSize] = i;
			continue;
    }
    if (connectVec[i].SequenceNumber() > maxSequenceNumberIndependentConnect && !connectVec[i].IsCondensed())
    {
      maxSequenceNumberIndependentConnect  = connectVec[i].SequenceNumber();
    }
  }
  for (int i = 0; i < depConInd.size(); i++) {
    if (connectVec[ depConInd[i] ].SequenceNumber() < maxSequenceNumberIndependentConnect ) {
      
      std::cout<<"A connect that has dependency has "
      <<"a sequency number smaller than an independent connect. "
      <<"have you tried this->CleanUpUnconnectedNodes() after "
      <<"building the computational mesh?"<<std::endl;
      DebugStop();
    }
  }
#endif
  
	int64_t neq = NEquations();
	skyline.Resize(neq);
    if (neq == 0) {
        return;
    }
//	cout << "Element skyline\n";
	//int eleq=0;
	int64_t i, n, l, nelem = NElements();
	for(i=0; i<neq; i++) skyline[i] = i;
	for(i=0; i<nelem; i++) {
		TPZCompEl *el = fElementVec[i];
		if(!el) continue;
		//      if(!el) continue;
		connectlist.Resize(0);
		el->BuildConnectList(connectlist);
		int64_t nnod = connectlist.NElements();
		if(!nnod) continue;
        // look for a connect with global equations associated to it
		int64_t ifirstnode = 0;
		TPZConnect *np = &fConnectVec[connectlist[0]];
		while(ifirstnode < nnod && (np->HasDependency() || np->IsCondensed()) ) {
			ifirstnode++;
            if (ifirstnode == nnod) {
                break;
            }
			np = &fConnectVec[connectlist[ifirstnode]];
		}
		int64_t ibl = np->SequenceNumber();
		int64_t loweq = fBlock.Position(ibl);
		int64_t higheq = loweq+fBlock.Size(ibl)-1;
		for(n=ifirstnode;n<nnod;n++) {
			np = &fConnectVec[connectlist[n]];
			if(np->HasDependency() || np->IsCondensed()) continue;
			int64_t ibl = np->SequenceNumber();
			int64_t leq = fBlock.Position(ibl);
			int64_t heq = leq+fBlock.Size(ibl)-1;
			//for(int _eq=leq; _eq<= heq; _eq++) {
			//   if((eleq%20==0)) cout << endl;
			//   cout << _eq << ' ';
			//   eleq++;
			//}
			loweq = (loweq > leq) ? leq : loweq;
			higheq = (higheq < heq) ? heq : higheq;
		}
//        std::cout << "Element " << i << " loweq " << loweq << " higheq " << higheq << std::endl;
//		cout << "Equations ";
		for(n=ifirstnode;n<nnod;n++) {
			np = &fConnectVec[connectlist[n]];
			if(np->HasDependency() || np->IsCondensed()) continue;
			int64_t ibl = np->SequenceNumber();
			int64_t leq = fBlock.Position(ibl);
			int64_t heq = leq+fBlock.Size(ibl);
			for(l=leq;l<heq;l++) {
				skyline[l] = skyline[l] < loweq ? skyline[l] : loweq;
//                cout << l << "/" << skyline[l] << " ";
			}
		}
//        cout << endl;
	}
//    std::cout << "Skyline " << skyline << std::endl;
}

void TPZCompMesh::BuildTransferMatrix(TPZCompMesh &coarsemesh, TPZTransfer<STATE> &transfer) {
	
	//TPZBlock &localblock = Block();
	TPZBlock &localblock = Block();
	//TPZBlock &coarseblock = coarsemesh.Block();
	TPZBlock &coarseblock = coarsemesh.Block();
	// adapt the block size of the blocks, dividing by the number of variables
	//  of the material
	int nmat = NMaterials();
	if(!nmat) {
		PZError << "TPZCompMesh::BuildTransferMatrix, no material object found\n";
		return;
	}
	TPZMaterial * mat;
	mat = fMaterialVec.begin()->second;
	int nvar = mat->NStateVariables();
	int dim = mat->Dimension();
	
	transfer.SetBlocks(localblock,coarseblock,nvar,NIndependentConnects(),coarsemesh.NIndependentConnects());
	Reference()->ResetReference();
	coarsemesh.LoadReferences();
	int64_t nelem = NElements();
	int64_t i;
	for(i=0; i<nelem; i++) {
		if(!fElementVec[i]) continue;
		TPZInterpolationSpace * finecel = dynamic_cast<TPZInterpolationSpace *> (fElementVec[i]);
		if(!finecel) continue;
		if(finecel->Dimension() != dim) continue;
		TPZGeoEl *finegel = finecel->Reference();
		TPZGeoEl *coarsegel = finegel;
		if(!finegel) {
			cout << "TPZCompMesh::BuildTransferMatrix is not implemented for super elements\n";
			continue;
		}
		
		while(coarsegel && !coarsegel->Reference()) {
			coarsegel = coarsegel->Father();
		}
		if(!coarsegel) {
			cout << "TPZCompMesh::BuildTransferMatrix corresponding coarse element not found\n";
			finecel->Print(cout);
			continue;
		}
		
		TPZInterpolationSpace * coarsel = dynamic_cast<TPZInterpolationSpace *> ( coarsegel->Reference() );
		if(!coarsel) continue;
		
		if(coarsel->Mesh() != &coarsemesh) {
			cout << "TPZCompMesh::BuildTransferMatrix is not implemented for transfers"
			" between superelements\n";
			continue;
		}
		TPZTransform<STATE> t(coarsel->Dimension());
		t=finegel->BuildTransform2(finegel->NSides()-1,coarsegel,t);
		finecel->BuildTransferMatrix(*coarsel,t,transfer);
	}
}

void TPZCompMesh::BuildTransferMatrixDesc(TPZCompMesh &transfermesh,
										  TPZTransfer<STATE> &transfer) {
#ifdef STATE_COMPLEX
	DebugStop();
#else
	//TPZBlock &localblock = Block();
	TPZBlock &localblock = Block();
	//TPZBlock &transferblock = transfermesh.Block();
	TPZBlock &transferblock = transfermesh.Block();
	// adapt the block size of the blocks, dividing by the number of variables
	//  of the material
	int nmat = NMaterials();
	if(!nmat) {
		PZError << "TPZCompMesh::BuildTransferMatrixDesc no material object found\n";
		return;
	}
	TPZMaterial * mat;
	mat = fMaterialVec.begin()->second;
	int nvar = mat->NStateVariables();
	int dim = mat->Dimension();
	//o seguinte �igual ao nmero de conects da malha
	int64_t ncon = NIndependentConnects(),coarncon = transfermesh.NIndependentConnects();
	transfer.SetBlocks(localblock,transferblock,nvar,ncon,coarncon);
	Reference()->ResetReference();//geom�ricos apontam para nulo
	transfermesh.LoadReferences();
	//geom�ricos apontam para computacionais da malha atual
	TPZAgglomerateElement *aggel = 0;
	TPZAdmChunkVector<TPZCompEl *> &elvec = transfermesh.ElementVec();
	int64_t nelem = elvec.NElements();
	int64_t i;
	for(i=0; i<nelem; i++) {
		TPZCompEl *comp = elvec[i];
		if(!comp) continue;
		if(comp->Dimension() != dim) continue;
		if(comp->Type() != EAgglomerate){
			PZError << "TPZCompMesh::BuildTransferMatrixDesc mesh agglomerated"
			<< " with element of volume not agglomerated\n";
			continue;
		}
		aggel = dynamic_cast<TPZAgglomerateElement *>(comp);
		//TPZStack<int> elvec;
		//retorna todos os descont�uos aglomerados por aggel
		//aggel->IndexesDiscSubEls(elvec);
		//int size = elvec.NElements(),i;
		int64_t size = aggel->NIndexes(),i;
		for(i=0;i<size;i++){
			//TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(fElementVec[elvec[i]]);
			//      TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(aggel->FineElement(i));
			TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(aggel->SubElement(i));
			if(!disc){
				PZError << "TPZCompMesh::BuildTransferMatrixDesc index with null"
				<< " elemento\n";
				continue;
			}
			if(disc->Type() != EDiscontinuous) {
				PZError << "TPZCompMesh::BuildTransferMatrixDesc index of not"
				<< " discontinous element\n";
				continue;
			}
			disc->BuildTransferMatrix(*aggel,transfer);
		}
	}
#endif
}

/*
 void TPZCompMesh::CreateConnectBC() {
 TPZGeoMesh *geo = Reference();
 TPZAdmChunkVector<TPZGeoNodeBC> *geobndcondvec = &geo->BCNodeVec();
 int ibc, nnodebc = geobndcondvec->NElements();
 for(ibc=0; ibc<nnodebc; ibc++) {
 TPZGeoNodeBC &gbc = (*geobndcondvec)[ibc];
 int bcnumber = gbc.fBCId;
 TPZMaterial *bc = FindMaterial(bcnumber);
 if(!bc) continue;
 // find a geometric element which has a
 //      corresponding computational element
 TPZGeoElSide gel(gbc.fGeoEl,gbc.fGeoElSide);
 TPZStack<TPZGeoElSide> neighbourset;
 gel.AllNeighbours(neighbourset);
 int in=0,nneigh;
 nneigh = neighbourset.NElements();
 while(in<nneigh && !neighbourset[in].Reference().Exists()) in++;
 //     neighbour = gel.Neighbour();
 //     while(neighbour.Exists() && neighbour != gel && !neighbour.Reference().Exists()) {
 //       neighbour = neighbour.Neighbour();
 //     }
 //     if(!neighbour.Exists() || ! neighbour.Reference().Exists()) continue;
 if(in == nneigh) continue;
 TPZCompElSide cel = neighbourset[in].Reference();
 TPZConnect *df = &cel.Element()->Connect(neighbourset[in].Side());
 if(!df) continue;
 TPZConnectBC dfbc(df,(TPZBndCond *) bc);
 int key = BCConnectVec().AllocateNewElement();
 BCConnectVec()[key] = dfbc;
 }
 }
 */

/*
 void TPZCompMesh::ComputeConnecttoElGraph(TPZVec<int> &firstel, TPZVec<int> &connectelgraph){
 int connectstackstore[50];
 TPZStack<int> connectstack(connectstackstore,50);
 int i, ncon = NConnects();
 firstel.Resize(ncon+1);
 firstel[0] = 0;
 for(i=0; i<ncon; i++) {
 TPZConnect &c = fConnectVec[i];
 int seqnum = c.SequenceNumber();
 if(seqnum == -1) {
 firstel[i+1] = firstel[i];
 } else {
 firstel[i+1] = firstel[i] + c.NElConnected();
 }
 }
 connectelgraph.Resize(firstel[ncon]);
 connectelgraph.Fill(-1);
 int nelem = NElements();
 for(i=0; i<nelem; i++) {
 TPZCompEl *el = fElementVec[i];
 if(!el) continue;
 connectstack.Resize(0);
 el->BuildConnectList(connectstack);
 int in;
 ncon = connectstack.NElements();
 for(in=0; in<ncon; in++) {
 int ic = connectstack[in];
 int first = firstel[ic];
 int last = firstel[ic+1];
 while(connectelgraph[first] != -1 && first < last) first++;
 if(first == last) {
 PZError << "TPZCompMesh::ComputeConnecttoElGraph wrong data structure\n";
 continue;
 }
 connectelgraph[first] = i;
 }
 }
 }
 */
int64_t TPZCompMesh::NIndependentConnects() {
	int64_t i, ncon = NConnects();
	int64_t NIndependentConnects = 0;
	for(i=0; i<ncon; i++) {
		TPZConnect &c = fConnectVec[i];
		if(c.HasDependency() || c.IsCondensed() || c.SequenceNumber() == -1) continue;
		NIndependentConnects++;
	}
	return NIndependentConnects;
}

void TPZCompMesh::ComputeElGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex){
	
    std::set<int> mat_ids;
    for (auto item : fMaterialVec) {
        mat_ids.insert(item.first);
    }
    ComputeElGraph(elgraph, elgraphindex, mat_ids);
}

void TPZCompMesh::ComputeElGraph(TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, std::set<int> & mat_ids){
    
    int64_t i, ncon;
    TPZStack<int64_t> connectstack;
    int64_t nelem = NElements();
    elgraphindex.Resize(nelem+1);
    elgraphindex[0] = 0;
    elgraph.Resize(0);
    int64_t curel=0;
    int64_t nindep = this->NIndependentConnects();
    for(i=0; i<nelem; i++) {
        TPZCompEl *el = fElementVec[i];
        /// Apply the filter by material identifiers
        bool has_material_Q = true;
        if (el) {
            has_material_Q = el->HasMaterial(mat_ids);
        }
        if(!el || !has_material_Q){
            elgraphindex[curel+1]=elgraph.NElements();
            curel++;
            continue;
        }
        

        
        connectstack.Resize(0);
        el->BuildConnectList(connectstack);
        int64_t in;
        ncon = connectstack.NElements();
        for(in=0; in<ncon; in++) {
            int ic = connectstack[in];
            TPZConnect &c = fConnectVec[ic];
            if(c.HasDependency() || c.IsCondensed()) continue;
#ifdef  PZDEBUG
            if (c.SequenceNumber() >= nindep) {
                DebugStop();
            }
#endif
            elgraph.Push(c.SequenceNumber());
        }
        elgraphindex[curel+1]=elgraph.NElements();
        curel++;
    }
    elgraphindex.Resize(curel+1);
}

void TPZCompMesh::Divide(int64_t index,TPZVec<int64_t> &subindex,int interpolate) {
	
	TPZCompEl * el = fElementVec[index];
	if (!el) {
		PZError << "TPZCompMesh::Divide element not found index = " << index << endl;
		subindex.Resize(0);
		return;
	}
	if(index != el->Index()){
		PZError << "TPZCompMesh::Divide - element->Index() != index " << endl;
		subindex.Resize(0);
		return;
	}
	
	
	el->Divide(index,subindex,interpolate);
}

void TPZCompMesh::Coarsen(TPZVec<int64_t> &elements, int64_t &index, bool CreateDiscontinuous) {
	int64_t i;
	const int64_t nelem = elements.NElements();
	
	if(!nelem) {
		index = -1;
		return;
	}//if
	
	TPZGeoEl *father = 0;
	TPZCompEl *cel = fElementVec[elements[0]];
	if (!cel) {
		index = -1;
		return;
	}//if
	
	if(cel) father = cel->Reference()->Father();
	if(!father) {
		index = -1;
		return;
	}//if
	
	for(i=1;i<nelem;i++) {
		if(!fElementVec[elements[i]]) {
			index = -1;
			return;
		}//if
		
		TPZGeoEl *father2 = fElementVec[elements[i]]->Reference()->Father();
		if(!father2 || father != father2) {
			index = -1;
			return;
		}//if
		
	}//for i
	
	if(nelem != father->NSubElements()){
		cout << "TPZCompEl::Coarsen : incomplete list of elements sons\n";
	}//if
	
	for(i=0; i<nelem; i++) {
		TPZInterpolationSpace * cel = dynamic_cast<TPZInterpolationSpace*>(this->ElementVec()[elements[i]]);
		if (!cel) continue;
		cel->RemoveInterfaces();
		TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		if (intel){
			intel->RemoveSideRestraintsII(TPZInterpolatedElement::EDelete);
		}//if (intel)
		
		cel->Reference()->ResetReference();
	}//for
	
	
	for(i=0; i<nelem; i++) {
		TPZCompEl * cel = this->ElementVec()[elements[i]];
		if (cel) delete cel;
	}
	
	if (CreateDiscontinuous) fCreate.SetAllCreateFunctionsDiscontinuous();
	else fCreate.SetAllCreateFunctionsContinuous();
	
	TPZCompEl * newcel = CreateCompEl(father);
    index = newcel->Index();
	
	TPZCompElDisc * newdisc = dynamic_cast<TPZCompElDisc*>(newcel);
	if (newdisc){
		newdisc->SetDegree( this->GetDefaultOrder() );
	}
	
}//method


void TPZCompMesh::RemakeAllInterfaceElements(){
	
	int64_t n = this->ElementVec().NElements();
	
	for(int64_t i = 0; i < n; i++){
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>( this->ElementVec()[i] );
		if (!disc) continue;
		disc->RemoveInterfaces();
	}//for
	
#ifdef PZDEBUG
	{
		n = this->ElementVec().NElements();
		for(int64_t i = 0; i < n; i++){
			TPZCompEl * cel = this->ElementVec()[i];
			if (!cel) continue;
			//    MElementType type = cel->Type();
			TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
			if(face){
				PZError << __PRETTY_FUNCTION__ << " - At this point no TPZInterfaceElement may exist.\n";
			}
		}
	}
#endif
	
	n = this->ElementVec().NElements();
	for(int64_t i = 0; i < n; i++){
		TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>( this->ElementVec()[i] );
		if (!disc) continue;
		disc->CreateInterfaces();
	}//for
	
}//method

#ifdef PZ_LOG
#include "pzsubcmesh.h"
#endif

/**ExpandSolution must be called before calling this*/
// it is a gather permutation
void TPZCompMesh::Permute(TPZVec<int64_t> &permute) {

  try {
      if(fSolType==EReal)
          PermuteInternal<STATE>(fSolution, permute);
      else if (fSolType==EComplex)
          PermuteInternal<CSTATE>(fSolution, permute);
    
  } catch (...) {
    PZError << "Incompatible matrix type in ";
    PZError << __PRETTY_FUNCTION__ << '\n';
    PZError << std::endl;
    DebugStop();
  }
}

template<class TVar>
void TPZCompMesh::PermuteInternal(TPZFMatrix<TVar> &sol,TPZVec<int64_t> &permute) {
	
	ExpandSolution();
	//   if (permute.NElements() != fBlock.NBlocks()) {
	//     PZError << "TPZCompMesh::Permute : permute vector size not equal to fBlock size\n";
	//   }
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (this);
		if (submesh) {
			sout << "Index = " << submesh->Index() << " ";
		}
		sout << "Permutation " << permute;
        std::set<int64_t> permset;
        if (permute.size() != 0) {
            permset.insert(&permute[0],(&permute[permute.size()-1]+1));
        }
        sout << " permute size " << permute.size() << " no dist elements " << permset.size();
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	int64_t i,j;
	int64_t permutenel = permute.NElements();
	for (i = 0; i < permutenel; i++) fBlock.Set(permute[i],fSolutionBlock.Size(i));
	fBlock.Resequence();
	if (sol.Rows() != 0) {
		//TPZFMatrix<REAL>	newsol(fSolution);
		TPZFMatrix<TVar> newsol(sol);
		for (i=0;i<fBlock.NBlocks();i++) {
			int64_t oldpos = fSolutionBlock.Position(i);
			int64_t newpos;
			if(i < permutenel) {
				newpos = fBlock.Position(permute[i]);
			} else {
				newpos = fBlock.Position(i);
			}
			for (j=0;j<fSolutionBlock.Size(i);j++) sol.PutVal(newpos+j,0,newsol(oldpos+j,0));
		}    //a sol. inicial esta em newsol
	}
	
	fSolutionBlock = fBlock;
	int64_t ncon = NConnects();
	for(i=0; i<ncon; i++) {
		TPZConnect &df = fConnectVec[i];
		int64_t seqnum = df.SequenceNumber();
		if(seqnum == -1) continue;
		if(seqnum < permutenel) df.SetSequenceNumber(permute[seqnum]);
	}
}


template<class TVar>
void TPZCompMesh::ConnectSolutionInternal(std::ostream &out, const TPZFMatrix<TVar>&sol) const{
    
	out << "\n\t\tCONNECT INDEX SOLUTION:\n\n";
	const int64_t ncon = NConnects();
	for(int i=0; i<ncon; i++) {
		out << i << ") ";
		TPZConnect &df = ConnectVec()[i];
		int64_t seqnum = df.SequenceNumber();
		if(df.NElConnected()==0) {
			out << "free node" << endl;
		} else if (seqnum < 0 || Block().Size(seqnum)==0) {
			out << "non solution connect" << endl;
		} else {
			int64_t pos = Block().Position(seqnum);
			for(int64_t j=0;j<Block().Size(seqnum);j++)
				out << sol.Get(pos+j,0) << "  ";
			out << std::endl;
		}
	}
}
void TPZCompMesh::ConnectSolution(std::ostream & out) {
	if(fSolType==EReal)
            ConnectSolutionInternal<STATE>(out,fSolution);
        else
            ConnectSolutionInternal<CSTATE>(out,fSolution);
    
}


void TPZCompMesh::EvaluateError(bool store_error, TPZVec<REAL> &errorSum) {
    std::set<int> matset;
    for(auto matpair : this->MaterialVec()){
        TPZMaterial* mat = matpair.second;
        if(!mat) DebugStop();
        if(mat->Dimension() == this->Dimension()){
            TPZBndCond* bnd = dynamic_cast<TPZBndCond*>(mat);
            if(!bnd){
                matset.insert(mat->Id());
            }
        }
    }
    this->EvaluateError(store_error, errorSum, matset);
}

void TPZCompMesh::AccountForElementError(TPZCompEl* cel, bool store_error, TPZVec<REAL>& true_error,
                                         TPZVec<REAL>& errorSum, std::set<int> &matset) {
    
    // Skipping cels that are not included in the set of materials to compute error
    const int celmatid = cel->Reference()->MaterialId();
    if(matset.find(celmatid) == matset.end()) return;
    
    cel->EvaluateError(true_error, store_error);

    int64_t nerrors = true_error.NElements();
    errorSum.Resize(nerrors, 0.);
    for (int64_t ii = 0; ii < nerrors; ii++)
        errorSum[ii] += true_error[ii] * true_error[ii];

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "true_errors: ";
        for (int ierr = 0; ierr < nerrors; ierr++) {
            sout << true_error[ierr] << " ";
        }
        sout << "\n";
        sout << "acc_errors^2: ";
        for (int ierr = 0; ierr < nerrors; ierr++) {
            sout << errorSum[ierr] << " ";
        }
        sout << "\n";
        sout << "GelID: ";
        if (cel->Reference()) {
            sout << cel->Reference()->Index();
            TPZGeoElSide side(cel->Reference());
            TPZManVector<REAL> coord(3);
            side.CenterX(coord);
            sout << " CenterCoord: " << coord;
            sout << " MatID: " << cel->Material()->Id();
        }
        sout << "\n";
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZCompMesh::EvaluateError(bool store_error, TPZVec<REAL> &errorSum, std::set<int> &matset) {
	
    if(!matset.size()) DebugStop();

	errorSum.Fill(0.);
	
	TPZManVector<REAL,7> true_error(errorSum.size());
	true_error.Fill(0.);
	
	TPZCompEl *cel;
	int gridDim = Dimension();
	
	//soma de erros sobre os elementos
	for(int64_t el=0;el< fElementVec.NElements();el++) {
        cel = fElementVec[el];
        if (!cel) continue;
        
        TPZCondensedCompEl* condcompel = dynamic_cast<TPZCondensedCompEl*>(cel);
        if(condcompel) {
            TPZElementGroup* elgr = dynamic_cast<TPZElementGroup*>(condcompel->ReferenceCompEl());
            if(elgr){ // loop over elements in the condensed element and account for the error of each one
                for(auto celgr : elgr->GetElGroup()){
                    if(!celgr) continue;
                    AccountForElementError(celgr,store_error,true_error,errorSum,matset);
                }
            }
            else{ // condcompel->ReferenceCompEl() is a single element
                AccountForElementError(condcompel->ReferenceCompEl(),store_error,true_error,errorSum,matset);
            }

        }
        else{ // It is a normal compel. Just account for the element error
          TPZElementGroup* elgr = dynamic_cast<TPZElementGroup*>(cel);
          if(elgr){ // loop over elements in the condensed element and account for the error of each one
            for(auto celgr : elgr->GetElGroup()){
              if(!celgr) continue;
              AccountForElementError(celgr,store_error,true_error,errorSum,matset);
            }
          }
          else {
              TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
              if(subcmesh) {
                  TPZCompMesh *cmesh = subcmesh;
                  cmesh->EvaluateError(store_error, errorSum, matset);
              }
              else{
                  AccountForElementError(cel,store_error,true_error,errorSum,matset);
              }
          }
        }
	}
	
	int64_t nerrors = errorSum.NElements();
	for(int64_t ii = 0; ii < nerrors; ii++)
		errorSum[ii] = sqrt(errorSum[ii]);
}

void TPZCompMesh::AdjustBoundaryElements() {
	int changed = 1;
	while(changed) {
		changed = 0;
		int64_t nel = fElementVec.NElements();
		int64_t el;
		TPZVec<int64_t> subelindex;
		for(el=0; el<nel; el++) {
			TPZStack<TPZCompElSide> elvec;
			TPZCompEl *elp = fElementVec[el];
			
			if(!elp || !dynamic_cast<TPZInterpolatedElement*>(elp) ) continue;
			
			TPZMaterial * mat = elp->Material();
			// this statement determines thata the element is associated with a boundary condition
			TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
			if(!bnd) continue;
			//      if(mat && mat->Id() >= 0) continue;
			int nsides = elp->Reference()->NSides();
			int is;
			int nc = elp->Reference()->NCornerNodes();
			for(is=nc; is<nsides; is++) {
				TPZCompElSide elpside(elp,is);
				// this should never be called
				if(elp->Reference()->SideDimension(is) == 0) continue;
				elvec.Resize(0);
				// verify if there are smaller elements connected to this boundary element
				elpside.HigherLevelElementList(elvec,0,0);
				if(elvec.NElements()) {
					// the first small element
					TPZGeoElSide fatherside = elvec[0].Reference();//el BC
					// elp is the element we will eventually refine
					int face = elp->Reference()->NSides() - 1;
					// this is the volume side of the element
					while(fatherside.Exists()) {
						if(elp->Reference()->NeighbourExists(face,fatherside)) break;
						fatherside = fatherside.Father2();
					}
					// fatherside is a neighbour of the current element
					// I wouldnt know when this test could fail
					if(fatherside.Exists()) {
#ifdef PZ_LOG
                        if(logger.isDebugEnabled())
						{
							std::stringstream sout;
							sout << "Dividing element " << el << " of type " << elp->Reference()->TypeName();
							LOGPZ_DEBUG(logger,sout.str().c_str());
						}
#endif
						Divide(el,subelindex,0);
						changed = 1;
						break;
					}
				}
				// we are working on the last side and nothing was divided
				if(is == nsides-1) {
					TPZInterpolatedElement *elpint = dynamic_cast <TPZInterpolatedElement *> (elp);
					if(!elpint) continue;
					int porder = elpint->PreferredSideOrder(is);
					int maxorder = 0;
					elvec.Resize(0);
					elpside.EqualLevelElementList(elvec,0,0);
					int64_t eq;
					for(eq=0; eq<elvec.NElements(); eq++) {
						TPZInterpolatedElement *eqel = dynamic_cast<TPZInterpolatedElement *> (elvec[eq].Element());
						int eqside = elvec[eq].Side();
						if(!eqel) continue;
						if(maxorder < eqel->PreferredSideOrder(eqside)) maxorder =  eqel->PreferredSideOrder(eqside);
					}
					// set the order to the largest order of all connecting elements
					if(porder < maxorder) {
#ifdef PZ_LOG
                        if(logger.isDebugEnabled())
						{
							std::stringstream sout;
							sout << "Refining element " << el << " to order " << maxorder;
							LOGPZ_DEBUG(logger,sout.str().c_str());
						}
#endif
						elpint->PRefine(maxorder);
						changed = 1;
					}
				}
			}
		}
	}
    InitializeBlock();
}

int64_t TPZCompMesh::PutinSuperMesh (int64_t local, TPZCompMesh *super){
	if (super != this) return -1;
	else return local;
}

int64_t TPZCompMesh::GetFromSuperMesh (int64_t superind, TPZCompMesh *super){
	if (super != this) return -1;
	else return superind;
}

REAL TPZCompMesh::CompareMesh(int var, char *matname){
	
	REAL error = 0.;
	int64_t i=0;
	for (i=0;i<fElementVec.NElements();i++){
		TPZCompEl *el = fElementVec[i];
		if(el) error+= el->CompareElement(var,matname);
	}
	return (error);
}


template<class TVar>
void TPZCompMesh::SetElementSolution(int64_t i, TPZVec<TVar> &sol){
    try{
        SetElementSolutionInternal<TVar>(fElementSolution,i,sol);
    } catch(...){
        PZError << "Incompatible matrix type in ";
        PZError << __PRETTY_FUNCTION__ << '\n';
        PZError << std::endl;
        DebugStop();
    }
}

template<class TVar>
void TPZCompMesh::SetElementSolutionInternal(TPZFMatrix<TVar> &mysol, int64_t i, TPZVec<TVar> &sol) {
	if(sol.NElements() != NElements()) {
		cout << __PRETTY_FUNCTION__<<" size of the vector doesn't match\n";
	}
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        TVar norm=0.;
        for (int64_t ii=0; ii<sol.size(); ii++) {
            norm += sol[ii];
        }
        norm = sqrt(norm);
        sout << "Norma da solucao " << i << " norma " << norm;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	if(mysol.Cols() <= i) mysol.Resize(NElements(),i+1);
	int64_t el,nel= NElements();
	for(el=0; el<nel; el++) {
		mysol.PutVal(el,i,sol[el]);
	}
}

void TPZCompMesh::GetRefPatches(std::set<TPZGeoEl *> &grpatch){
	int64_t i;
	int64_t nel = NElements();
	//	cout << "GetRefPatches\n" << nel << endl;
	for (i=0; i<nel; i++){
		if (fElementVec[i]){
			TPZGeoEl *gel = fElementVec[i]->Reference();
			if (gel) gel = fElementVec[i]->GetRefElPatch();
			if (gel)
			{
				grpatch.insert(gel);
			}
		}
	}
	return;
}


void  TPZCompMesh::GetNodeToElGraph(TPZVec<int64_t> &nodtoelgraph, TPZVec<int64_t> &nodtoelgraphindex, TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex){
	
	ComputeElGraph(elgraph,elgraphindex);
	
	//  int i,j;
	//  for (i=0; i<elgraphindex.NElements()-1; i++){
	//    cout << "Block: " << i << "\t";
	//    for (j = elgraphindex[i]; j<elgraphindex[i+1]; j++){
	//      cout << elgraph[j] << "\t";
	//    }
	//   cout << endl;
	//  }
	
	TPZRenumbering renum(elgraphindex.NElements() -1 ,NIndependentConnects());
	renum.SetElementGraph(elgraph,elgraphindex);
	renum.NodeToElGraph(elgraph,elgraphindex,nodtoelgraph, nodtoelgraphindex);
	/*   TPZRenumbering *re = &renum; */
	/*   re->Print(nodtoelgraph, nodtoelgraphindex , "Grapho de n� para elementos "); */
	
	
}


void TPZCompMesh::GetElementPatch(TPZVec<int64_t> nodtoelgraph, TPZVec<int64_t> nodtoelgraphindex, TPZStack<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex,int64_t elind ,TPZStack<int64_t> &patch){
	
	//  int aux =0;
	//TPZAVLMap<int,int> elconmap(aux);
	std::set<int64_t > elconmap;
	int64_t i,j;
	for (i= elgraphindex[elind]; i<elgraphindex[elind+1];i++){
		int64_t node = elgraph[i];
		for (j = nodtoelgraphindex[node];  j<nodtoelgraphindex[node+1]; j++){
			elconmap.insert(nodtoelgraph[j]);
			
		}
	}
	patch.Resize(0);
	
	//TPZPix iter = elconmap.First();
	set<int64_t >::iterator iter = elconmap.begin();
	while(iter!=elconmap.end()){
		//patch.Push(elconmap.Key(iter));
		patch.Push(*iter);//elconmap.Key(iter));
		iter++;//elconmap.Next(iter);
	}
}

TPZCompMesh::TPZCompMesh(const TPZCompMesh &copy) :
TPZRegisterClassId(&TPZCompMesh::ClassId),
fReference(copy.fReference),fConnectVec(copy.fConnectVec),
fMaterialVec(), fSolutionBlock(copy.fSolutionBlock),
fCreate(copy.fCreate), fBlock(copy.fBlock),
fSolution(copy.fSolution),
fElementSolution(copy.fElementSolution),fDimModel(copy.fDimModel),
fSolN(copy.fSolN),
fSolType(copy.fSolType)
{
#ifdef PZ_LOG
    if (aloclogger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Allocate TPZCompMesh this = " << (void *)this;
        LOGPZ_DEBUG(aloclogger, sout.str())
    }
#endif

	fDefaultOrder = copy.fDefaultOrder;
	fReference->ResetReference();
    TPZBaseMatrix &sol = fSolution;
	fBlock.SetMatrix(&sol);
    fSolutionBlock.SetMatrix(&sol);
	copy.CopyMaterials(*this);
	int64_t nel = copy.fElementVec.NElements();
	fElementVec.Resize(nel);
	int64_t iel;
	for(iel = 0; iel<nel; iel++) fElementVec[iel] = 0;
	for(iel = 0; iel<nel; iel++) {
		TPZCompEl *cel = copy.fElementVec[iel];
		if(cel && !dynamic_cast<TPZInterfaceElement* >(cel) )
		{
			TPZCompEl *clone  =  cel->Clone(*this);
			/*#ifdef PZ_LOG
			 {
			 std::stringstream sout;
			 sout << "original\n";
			 cel->Print(sout);
			 sout << "cloned\n";
			 clone->Print(sout);
			 LOGPZ_DEBUG(logger,sout.str())
			 }
			 #endif*/
		}
	}
	/** Update data into all the connects */
	ComputeNodElCon();
//	int nconn = fConnectVec.NElements();
//	for(iel=0;
	/*#ifdef PZ_LOG
	 {
	 std::stringstream sout;
	 Print(sout);
	 LOGPZ_DEBUG(logger,sout.str())
	 }
	 #endif*/
//	for(iel = 0; iel<nel; iel++) {
//		TPZCompEl *cel = copy.fElementVec[iel];
//		if(cel && dynamic_cast<TPZInterfaceElement* >(cel) ) cel->Clone(*this);
//	}
	fDimModel = copy.fDimModel;
	fName = copy.fName;
}

TPZCompMesh &TPZCompMesh::operator=(const TPZCompMesh &copy)
{
    CleanUp();
    fReference = copy.fReference;
    fReference->ResetReference();
    fConnectVec = copy.fConnectVec;
    copy.CopyMaterials(*this);    
    fElementSolution = copy.fElementSolution;
    fSolution = copy.fSolution;
    fSolType = copy.fSolType;
    TPZBaseMatrix &sol = fSolution;
    
    fSolutionBlock = copy.fSolutionBlock;
    fSolutionBlock.SetMatrix(&sol);
    fBlock = copy.fBlock;
    fBlock.SetMatrix(&sol);
    fSolN = copy.fSolN;
    fDefaultOrder = copy.fDefaultOrder;
    int64_t nel = copy.fElementVec.NElements();
    fElementVec.Resize(nel);
    int64_t iel;
    for(iel = 0; iel<nel; iel++) fElementVec[iel] = nullptr;
    for(iel = 0; iel<nel; iel++) {
        TPZCompEl *cel = copy.fElementVec[iel];
        if(cel && !dynamic_cast<TPZInterfaceElement* >(cel) )
        {
                cel->Clone(*this);
        }
    }
    for(iel = 0; iel<nel; iel++) 
    {
        TPZCompEl *cel = copy.fElementVec[iel];
        if(cel && dynamic_cast<TPZInterfaceElement* >(cel) ) cel->Clone(*this);
    }
    fDimModel = copy.fDimModel;
    fName = copy.fName;
    fCreate = copy.fCreate;

    return *this;
}

TPZCompMesh* TPZCompMesh::Clone() const {
	return new TPZCompMesh(*this);
}
/*
 TPZGeoMesh *gmesh = Reference();
 gmesh->ResetReference();
 TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
 CopyMaterials(cmesh);
 
 int iel, nel = fElementVec.NElements();
 for (iel=0;iel<nel;iel++){
 if (!fElementVec[iel]) continue;
 TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (fElementVec[iel]);
 
 TPZGeoEl *gel = fElementVec[iel]->Reference();
 if (!gel) continue;
 int indexcomp;
 gel->CreateCompEl(*cmesh,indexcomp);
 
 TPZInterpolatedElement *clintel = dynamic_cast<TPZInterpolatedElement *> (cmesh->ElementVec()[indexcomp]);
 int side, nsides = intel->NConnects();
 for (side =0;side<nsides;side++){
 int porder = intel->SideOrder(side);
 clintel->PRefine(side,porder);
 int clorder = clintel->SideOrder(side);
 if(porder != clorder) {
 cout << "PZCompMesh::Clone order was not honoured side " <<
 side << " porder " << porder << " clorder " << clorder <<
 endl;
 }
 int elndof = intel->Connect(side).NDof(*this);
 int clelndof = clintel->Connect(side).NDof(*cmesh);
 if(elndof != clelndof) {
 cout << "PZCompMesh::Clone ndof different side " <<
 side << " porder " << porder << " clorder " << clorder <<
 " elndof = "<< elndof << " clelndof " << clelndof << endl;
 }
 }
 }
 
 cmesh->InitializeBlock();
 
 gmesh->ResetReference();
 LoadReferences();
 
 
 int clnel = cmesh->ElementVec().NElements();
 //  if (clnel != nel){
 //    cout << "Clone mesh inconsistancy: clone elements failed!\n";
 //    return 0;
 //  }
 
 for (iel=0;iel<clnel;iel++){
 TPZCompEl *cel = cmesh->ElementVec()[iel];
 TPZGeoEl *gel = cel->Reference();
 TPZCompEl *el = gel->Reference();
 int ic,nc = el->NConnects();
 for (ic=0;ic<nc;ic++){
 int elseqnum = el->Connect(ic).SequenceNumber();
 int elndof = el->Connect(ic).NDof(*this);
 int clelseqnum = cel->Connect(ic).SequenceNumber();
 int clelndof = cel->Connect(ic).NDof(*cmesh);
 if (elndof != clelndof){
 cout << "Number of degree of freedom incompatible between clone and original mesh!\n";
 cout << "Clone connect id: " << ic << "  Clone dof: " << clelndof << "  Original dof: " << clelndof << endl;
 continue;
 }
 int idof;
 for (idof=0; idof<elndof; idof++){
 cmesh->fSolutionBlock(clelseqnum,0,idof,0) = Block()(elseqnum,0,idof,0);
 }
 }
 }
 return cmesh;
 }
 */

void TPZCompMesh::CopyMaterials(TPZCompMesh &meshto) const {
    // Clone volumetric mats
    for (auto it : fMaterialVec) {
        if (!dynamic_cast<TPZBndCond *> (it.second)) {
            it.second->Clone(meshto.fMaterialVec);
        }
    }
    // Clone BC mats
    for (auto it : fMaterialVec) {
        auto *bc = dynamic_cast<TPZBndCond *> (it.second);
        if (bc) {
            it.second->Clone(meshto.fMaterialVec);
            auto *cloned_mat = meshto.FindMaterial(bc->Material()->Id());
            auto *new_bc = dynamic_cast<TPZBndCond*>(meshto.FindMaterial(bc->Id()));
            if (!new_bc) DebugStop();
            new_bc->SetMaterial(cloned_mat);
        }
    }
}

void TPZCompMesh::DeleteMaterial(const int matId) {
    delete this->MaterialVec()[matId];
    this->MaterialVec().erase(matId);
}

REAL TPZCompMesh::DeltaX(){
	
	int64_t nel = ElementVec().NElements(),i,j;
	if(nel == 0) cout << "\nTPZCompMesh::DeltaX no elements\n";
	REAL maxdist = 0.0,dist=0.0;
	TPZVec<REAL> point0(3),point1(3);
	TPZGeoNode *node0,*node1;
	for(i=0;i<nel;i++){
		TPZCompEl *com = ElementVec()[i];
		if(!com) continue;
		if(com->Type() == EInterface || com->Material()->Id() < 0) continue;
		node0 = com->Reference()->NodePtr(0);
		node1 = com->Reference()->NodePtr(1);
		for(j=0;j<3;j++){
			point0[j] = node0->Coord(j);
			point1[j] = node1->Coord(j);
		}
		dist = TPZGeoEl::Distance(point0,point1);
		if(dist > maxdist) maxdist = dist;
	}
	return maxdist;
}

REAL TPZCompMesh::MaximumRadiusOfMesh(){
	
	int64_t nel = ElementVec().NElements(),i;
	if(nel == 0) cout << "\nTPZCompMesh::MaximumRadiusOfMesh no elements\n";
	REAL maxdist = 0.0,dist=0.0;
	TPZVec<REAL> point0(3),point1(3);
	for(i=0;i<nel;i++){
		TPZCompEl *com = ElementVec()[i];
		if(!com) continue;
		if(com->Type() == EInterface || com->Material()->Id() < 0) continue;
		dist = com->MaximumRadiusOfEl();
		if(dist > maxdist) maxdist = dist;
	}
	return maxdist;
}

REAL TPZCompMesh::LesserEdgeOfMesh(){
	
	int64_t nel = ElementVec().NElements(),i;
	if(nel == 0) cout << "\nTPZCompMesh::MaximumRadiusOfMesh no elements\n";
	REAL mindist =10000.0,dist=0.0;
	for(i=0;i<nel;i++){
		TPZCompEl *com = ElementVec()[i];
		if(!com) continue;
		int type = com->Type();
		if( type == EInterface || com->Material()->Id() < 0 ) continue;
		if(type == EAgglomerate) dist = dynamic_cast<TPZCompElDisc *>(com)->LesserEdgeOfEl();
		else dist = com->LesserEdgeOfEl();
		if(dist < mindist) mindist = dist;
	}
	return mindist;
}

/** This method will fill the matrix passed as parameter with a representation of the fillin of the global stiffness matrix, based on the sequence number of the connects
 @param resolution Number of rows and columns of the matrix
 @param fillin Matrix which is mapped onto the global system of equations and represents the fillin be assigning a value between 0. and 1. in each element */
void TPZCompMesh::ComputeFillIn(int64_t resolution, TPZFMatrix<REAL> &fillin){
	ComputeNodElCon();
	int64_t nequations = NEquations();
	int64_t divider = nequations/resolution;
	if(divider*resolution != nequations) divider++;
	REAL factor = 1./(divider*divider);
	fillin.Redim(resolution,resolution);
	
	TPZStack<int64_t> graphelindex, graphel, graphnodeindex, graphnode;
	this->ComputeElGraph(graphel,graphelindex);
	TPZRenumbering renum(fElementVec.NElements(),fConnectVec.NElements());
	renum.ConvertGraph(graphel,graphelindex,graphnode,graphnodeindex);
	std::map<int64_t,TPZConnect *> seqtoconnect;
	int ic,ncon = fConnectVec.NElements();
	for(ic=0; ic<ncon; ic++) {
		TPZConnect &c = fConnectVec[ic];
		if(c.HasDependency() || c.IsCondensed() || c.SequenceNumber() < 0) continue;
		seqtoconnect[c.SequenceNumber()] = &c;
	}
	int64_t iseqnum;
	for(iseqnum = 0; iseqnum < graphnodeindex.NElements()-1; iseqnum++) {
		if(!seqtoconnect.count(iseqnum)) continue;
		int64_t firstieq = Block().Position(iseqnum);
		int64_t lastieq = Block().Size(iseqnum)+firstieq;
		int64_t firstnode = graphnodeindex[iseqnum];
		int64_t lastnode = graphnodeindex[iseqnum+1];
		{
			int64_t ieq;
			for(ieq=firstieq; ieq<lastieq; ieq++) {
				int64_t rowp = ieq/divider;
				int64_t ieq2;
				for(ieq2=firstieq; ieq2<lastieq; ieq2++) {
					int64_t rowp2 = ieq2/divider;
					fillin(rowp,rowp2) += factor;
				}
			}
		}
		int64_t in;
		for(in=firstnode; in<lastnode; in++) {
			int64_t jseqnum = graphnode[in];
			int64_t firstjeq = Block().Position(jseqnum);
			int64_t lastjeq = Block().Size(jseqnum)+firstjeq;
			int64_t ieq;
			for(ieq=firstieq; ieq<lastieq; ieq++) {
				int64_t rowp = ieq/divider;
				int64_t jeq;
				for(jeq=firstjeq; jeq<lastjeq; jeq++) {
					int64_t colp = jeq/divider;
					fillin(rowp,colp) += factor;
				}
			}
		}
	}
}

template<class TVar>
void TPZCompMesh::ProjectSolution(TPZFMatrix<TVar> &projectsol) {
	
	//  * * A MALHA ATUAL DEVE SER AGLOMERADA * * *
	
	//   TPZBlock &localblock = Block();
	//   TPZBlock &transferblock = finemesh.Block();
	// adapt the block size of the blocks, dividing by the number of variables
	//  of the material
	int64_t neq = NEquations();
	projectsol.Redim(neq,1);
	projectsol.Zero();
	int nmat = NMaterials();
	if(!nmat) {
		PZError << "TPZCompMesh::BuildTransferMatrixDesc2 no material object found\n";
		return;
	}
	Reference()->ResetReference();//geom�ricos apontam para nulo
	LoadReferences();
	
#ifdef STATE_COMPLEX
	DebugStop();
#else
	TPZMaterial * mat = fMaterialVec.begin()->second;
	//geom�ricos apontam para computacionais da malha atual
    int dim = mat->Dimension();
	TPZAgglomerateElement *aggel = 0;
	TPZAdmChunkVector<TPZCompEl *> &elvec = ElementVec();
	int64_t nelem = elvec.NElements();

	for(int64_t i=0; i<nelem; i++) {
		TPZCompEl *comp = elvec[i];
		if(!comp) continue;
		if(comp->Dimension() != dim) continue;
		if(comp->Type() != EAgglomerate){
			PZError << "TPZCompMesh::BuildTransferMatrixDesc2 mesh agglomerated"
			<< " with element of volume not agglomerated\n";
			continue;
		}
		aggel = dynamic_cast<TPZAgglomerateElement *>(comp);
		aggel->ProjectSolution<TVar>(projectsol);
	}
#endif
}
/**
 * returns the unique identifier for reading/writing objects to streams
 */
int TPZCompMesh::ClassId() const{
    return Hash("TPZCompMesh");
}
/**
 Save the element data to a stream
 */
void TPZCompMesh::Write(TPZStream &buf, int withclassid) const { //ok
    TPZPersistenceManager::WritePointer(fReference,&buf);
    TPZPersistenceManager::WritePointer(fGMesh.operator->(), &buf);
    buf.Write((int)fSolType);
    buf.Write(&fName);
    buf.WritePointers(fElementVec);
    fConnectVec.Write(buf, withclassid);
    std::map<int, TPZMaterial*> internal_materials;
    std::map<int, TPZMaterial*> boundary_materials;
    for (auto mat_pair : fMaterialVec) {
        if (dynamic_cast<TPZBndCond*>(mat_pair.second)){
            boundary_materials.insert(mat_pair);
        } else {
            internal_materials.insert(mat_pair);
        }
    }
    buf.WritePointers(internal_materials);
    buf.WritePointers(boundary_materials);
    fSolutionBlock.Write(buf,0);
    fSolution.Write(buf,0);
    fSolN.Write(buf,0);
    fBlock.Write(buf,0);
    fElementSolution.Write(buf,0);
    buf.Write(&fDimModel);
    buf.Write(&fDefaultOrder);
    fCreate.Write(buf, withclassid);
    buf.Write(&fNmeshes);
	
}

/**
 Read the element data from a stream
 */
void TPZCompMesh::Read(TPZStream &buf, void *context) { //ok
    fReference = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fGMesh = TPZAutoPointerDynamicCast<TPZGeoMesh >(TPZPersistenceManager::GetAutoPointer(&buf));
    fSolType = [&buf](){
        int tmp;
        buf.Read(&tmp);
        return (ESolType) tmp;
    }();
    buf.Read(&fName);
    buf.ReadPointers(fElementVec);
    fConnectVec.Read(buf, context);
    buf.ReadPointers(fMaterialVec); //internal materials
    buf.ReadPointers(fMaterialVec); //boundary materials
    fSolutionBlock.Read(buf, NULL);
    fSolution.Read(buf,NULL);
    fSolN.Read(buf,NULL);
    fBlock.Read(buf, NULL);
    fElementSolution.Read(buf, NULL);
    buf.Read(&fDimModel);
    buf.Read(&fDefaultOrder);
    fCreate.Read(buf, context);
    buf.Read(&fNmeshes);
}

/// Integrate the postprocessed variable name over the elements included in the set matids
TPZVec<STATE> TPZCompMesh::Integrate(const std::string &varname, const std::set<int> &matids)
{
    // the postprocessed index of the varname for each material id
    std::map<int,int> variableids;
    int nvars = 0;
    
    std::map<int,TPZMaterial *>::iterator itmap;
    for (itmap = MaterialVec().begin(); itmap != MaterialVec().end(); itmap++) {
        if (matids.find(itmap->first) != matids.end()) {
            TPZMaterial *mat = itmap->second;
            TPZBndCond *bndcond = dynamic_cast<TPZBndCond *>(mat);
            int varindex = mat->VariableIndex(varname);
            if (varindex == -1 && bndcond) {
                mat = bndcond->Material();
                varindex= mat->VariableIndex(varname);
            }
            if (varindex == -1) {
                DebugStop();
            }
            variableids[itmap->first] = varindex;
            int nvarnew = mat->NSolutionVariables(variableids[itmap->first]);
            // the number of variables has to be the same for all materials
            if(nvars && nvars != nvarnew)
            {
                DebugStop();
            }
            nvars = nvarnew;
        }
    }
    TPZManVector<STATE,3> result(nvars,0.);
    int64_t nelem = NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZCompEl *cel = Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            TPZManVector<STATE,3> locres;
            locres = cel->IntegrateSolution(varname, matids);
            if (locres.size() == nvars) {
                for (int iv = 0; iv<nvars; iv++) {
                    result[iv] += locres[iv];
                }
            }
            else if(locres.size())
            {
                DebugStop();
            }
        }
        else
        {
            int matid = gel->MaterialId();
            if (matids.find(matid) == matids.end()) {
                continue;
            }
            TPZManVector<STATE,3> locres(nvars,0.);
            locres = cel->IntegrateSolution(variableids[matid]);
            if (locres.size() != nvars) {
                DebugStop();
            }
            for (int iv=0; iv<nvars; iv++)
            {
                result[iv] += locres[iv];
            }
            //        std::cout << "el = " << el << " integral " << locres << " result " << result << std::endl;
        }
    }
    return result;
}


/*
void TPZCompMesh::SaddlePermute()
{
    
#ifdef PZ_LOG
    {
        std::stringstream sout;
        sout<< "Implementando permutacao para problemas de ponto de sela"<< std::endl;
        LOGPZ_DEBUG(logger, sout.str().c_str());
    }
#endif
    TPZVec<int64_t> permute;
    int numinternalconnects = NIndependentConnects();
    permute.Resize(numinternalconnects,0);
    
    TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (this);
    if(submesh)
    {
        int nexternal = submesh->NConnects();
        numinternalconnects -= nexternal;
    }
    //	else {
    //		DebugStop();
    //	}
    
    int jperm=0;
    int nel=ElementVec().NElements();
    for (int jel=0; jel<nel; jel++) {
        
        for (int64_t ip=0; ip<permute.NElements(); ip++) {
            permute[ip]=ip;
        }
        
        TPZCompEl *cel= ElementVec()[jel];
        //	int idtroca=0;
        int eqmax=0;
        if(!cel)continue;
//        int ncon=elvec->NConnects();
        std::set<int> connects;
        cel->BuildConnectList(connects );
        //	if(ncon==1) continue;
        int pressureconectindex = cel->PressureConnectIndex();
        if(pressureconectindex == -1) continue;
        int64_t eqpress=cel->Connect(pressureconectindex).SequenceNumber();

        for (std::set<int>::const_iterator it= connects.begin(); it != connects.end(); it++) {
//        for (int icon=0; icon< ncon-1; icon++) {
            if(*it == pressureconectindex) continue;
            TPZConnect &coel= fConnectVec[*it];
            if(coel.HasDependency()) continue;
            int eqflux=coel.SequenceNumber();
            eqmax = max(eqmax,eqflux);
        }
        
        
        if(eqpress<eqmax) {
            
            permute[eqpress]=eqmax;
            
        }
        
        
        for ( jperm = eqpress+1; jperm<=eqmax; jperm++) {
            permute[jperm]=jperm-1;
            
        }
        
//         #ifdef PZ_LOG
//         {
//         std::stringstream sout;
//         sout << "vetor SaddlePermute  do elemento - "<<jel<< " - " <<permute;
//         LOGPZ_DEBUG(logger, sout.str().c_str());
//         }
//         #endif
         
        Permute(permute);
        
    }		
}
*/

static void switchEq(int64_t eqsmall, int64_t eqlarge, TPZVec<int64_t> &permutegather, TPZVec<int64_t> &permutescatter)
{
    int64_t eqkeep = permutegather[eqsmall];
    for (int64_t eq = eqsmall; eq< eqlarge; eq++) {
        permutegather[eq] = permutegather[eq+1];
    }
    permutegather[eqlarge] = eqkeep;
    for (int64_t eq = eqsmall; eq<= eqlarge; eq++) {
        permutescatter[permutegather[eq]] = eq;
    }
}


void TPZCompMesh::SaddlePermute()
{
    TPZVec<int64_t> permutegather,permutescatter;
    int64_t numinternalconnects = NIndependentConnects();
    permutegather.Resize(numinternalconnects,0);
    permutescatter.Resize(numinternalconnects,0);
    for (int64_t i=0; i<numinternalconnects; i++) {
        permutegather[i] = i;
        permutescatter[i] = i;
    }
    int64_t numconnects = ConnectVec().NElements();
    int64_t numindepconnects = NIndependentConnects();
    if (numconnects==0) {
        return;
    }
    int minlagrange = 0;
    int maxlagrange = 0;
    for (int64_t ic=0; ic<numconnects; ic++) {
        TPZConnect &c = ConnectVec()[ic];
        if(c.HasDependency() || c.IsCondensed()) continue;
        if (c.SequenceNumber() < 0) {
            continue;
        }
        minlagrange = c.LagrangeMultiplier();
        maxlagrange = c.LagrangeMultiplier();
        break;
    }
    for (int ic=0; ic<numconnects; ic++) {
        TPZConnect &c = ConnectVec()[ic];
        if(c.HasDependency() || c.IsCondensed()) continue;
        int lagrange = c.LagrangeMultiplier();
        minlagrange = Min(lagrange, minlagrange);
        maxlagrange = Max(lagrange,maxlagrange);
    }

    int64_t nel = NElements();
    for (int lagr = minlagrange+1; lagr <= maxlagrange; lagr++)
    {
        for (int64_t el = 0; el<nel ; el++) {
            TPZCompEl *cel = ElementVec()[el];
            if (!cel) {
                continue;
            }
            int nc = cel->NConnects();
            if (nc <= 1) {
                continue;
            }
        
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "el " << el << " Before renumbering : ";
                for (int ic=0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if (c.HasDependency() || c.IsCondensed()) {
                        continue;
                    }
                    sout << permutescatter[c.SequenceNumber()] << "/" << (int)c.LagrangeMultiplier() << " ";
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            // put all connects after the connect largest seqnum and lower lagrange number
            int64_t maxseq = -1;
            for (int ic=0; ic<nc ; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                int64_t eq = permutescatter[c.SequenceNumber()];
                if (!c.HasDependency() && c.LagrangeMultiplier() < lagr && eq > maxseq) {
                    maxseq = eq;
                }
            }
            std::set<int64_t> seteq;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                int eq = permutescatter[c.SequenceNumber()];
                if (c.LagrangeMultiplier() == lagr && eq < maxseq) {
#ifdef PZDEBUG
                    if (c.HasDependency()) {
                        DebugStop();
                    }
#endif
                    seteq.insert(eq);
                }
            }
            std::set<int64_t>::reverse_iterator it;
            int64_t count = 0;
            for (it = seteq.rbegin(); it != seteq.rend(); it++) {
                int64_t eq = *it;
#ifdef PZ_LOG
                if (logger.isDebugEnabled()) {
                    std::stringstream sout;
                    sout << "Switch ceq = " << eq << " with maxeq = " << maxseq-count;
                    //                    sout << permute;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                // we have to switch
                switchEq(eq, maxseq-count, permutegather, permutescatter);
                count++;
            }

#ifdef PZDEBUG
            for (int64_t i=0; i<numinternalconnects; i++) {
                if (permutescatter[permutegather[i]] != i) {
                    std::cout << "permutegather " << permutegather << std::endl;
                    std::cout << "permutescatter " << permutescatter << std::endl;
                    DebugStop();
                }
            }
            // put all connects after the connect largest seqnum and lower lagrange number
            maxseq = -1;
            for (int ic=0; ic<nc ; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.IsCondensed() || c.HasDependency()) {
                    continue;
                }
                int64_t eq = permutescatter[c.SequenceNumber()];
                if (c.LagrangeMultiplier() < lagr && eq < maxseq) {
                    maxseq = eq;
                }
            }
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                int eq = permutescatter[c.SequenceNumber()];
                if (c.LagrangeMultiplier() == lagr && eq < maxseq) {
                    // we have to switch
                    DebugStop();
                }
            }
#endif
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "el " << el << " After renumbering  : ";
                for (int ic=0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if (c.HasDependency() || c.IsCondensed()) {
                        continue;
                    }
                    sout << permutescatter[c.SequenceNumber()] << "/" << (int)c.LagrangeMultiplier() << " ";
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
        }
        
    }
#ifdef PZ_LOG2
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Saddle permute new permutation ";
        for (int i = 0; i< permutescatter.size(); i++) {
            int64_t jmax = i+10;
            if (jmax > permutescatter.size()) {
                jmax = permutescatter.size();
            }
            for (int64_t j=i; j< jmax; j++) {
                sout << permutescatter[j] << ' ';
            }
            sout << std::endl;
            i+= 10;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    Permute(permutescatter);
#ifdef PZDEBUG
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        LOGPZ_DEBUG(logger, "******************* AFTER PERMUTATION **************************")
    }
#endif
    
    for (int64_t i=0L; i<numinternalconnects; i++) {
        permutegather[i] = i;
        permutescatter[i] = i;
    }
    for (int64_t el = 0L; el<nel ; el++) {
        TPZCompEl *cel = ElementVec()[el];
        if (!cel) {
            continue;
        }
        int nc = cel->NConnects();
        if (nc == 0) {
            continue;
        }
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "el " << el << " Final numbering : ";
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                sout << permutescatter[c.SequenceNumber()] << "/" << (int)c.LagrangeMultiplier() << " ";
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZConnect &c0 = cel->Connect(0);
        int minlagrange = c0.LagrangeMultiplier();
        int maxlagrange = c0.LagrangeMultiplier();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if(c.HasDependency() || c.IsCondensed()) continue;
            int lagrange = c.LagrangeMultiplier();
            minlagrange = Min(lagrange, minlagrange);
            maxlagrange = Max(lagrange,maxlagrange);
        }
        for (int lagr = minlagrange+1; lagr <= maxlagrange; lagr++) {
            // put all connects after the connect largest seqnum and lower lagrange number
            int64_t maxseq = -1;
            for (int ic=0; ic<nc ; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                int64_t eq = permutescatter[c.SequenceNumber()];
                if (!c.HasDependency() && c.LagrangeMultiplier() < lagr && eq > maxseq) {
                    maxseq = eq;
                }
            }
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.HasDependency() || c.IsCondensed()) {
                    continue;
                }
                int eq = permutescatter[c.SequenceNumber()];
                if (c.LagrangeMultiplier() == lagr && eq < maxseq) {
                    // we have to switch
                    cel->Print(std::cout);
                    for (int iic=0; iic<nc; iic++) {
                        cel->Connect(iic).Print(*this, std::cout);
                        std::cout << "Destination seqnum " << permutegather[cel->Connect(iic).SequenceNumber()] << std::endl;
                    }
                    DebugStop();
                }
            }
        }
    }

#endif
}

void TPZCompMesh::SaddlePermute2()
{
    TPZVec<int64_t> permute;
    int64_t numinternalconnects = NIndependentConnects();
    permute.Resize(numinternalconnects,0);
    for (int64_t i=0L; i<numinternalconnects; i++) {
        permute[i] = i;
    }
    int64_t nel = NElements();
    for (int64_t el = 0L; el<nel ; el++) {
        TPZCompEl *cel = ElementVec()[el];
        if (!cel) {
            continue;
        }
        unsigned char minlagrange = 255, maxlagrange = 0;
        int nc = cel->NConnects();
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            int64_t seqnum = c.SequenceNumber();
            // if seqnum is larger than internal connects the equation is restrained, no permutation is necessary
            if (seqnum >= numinternalconnects) {
                continue;
            }
            unsigned char lagrange = c.LagrangeMultiplier();
            if (lagrange < minlagrange) {
                minlagrange = lagrange;
            }
            if (lagrange > maxlagrange) {
                maxlagrange = lagrange;
            }
        }
        for (unsigned char lagrange = minlagrange+1; lagrange <= maxlagrange; lagrange++) {
            int64_t maxeq = -1;
            for (int ic=0; ic<nc ; ic++) {
                TPZConnect &c = cel->Connect(ic);
                if (c.SequenceNumber() >= numinternalconnects) {
                    continue;
                }
                if (c.LagrangeMultiplier() < lagrange) {
                    int64_t origeq = c.SequenceNumber();
                    if(maxeq < permute[origeq])
                    {
                        maxeq = permute[origeq];
                    }
                }
            }
            if (maxeq < 0) {
                continue;
            }
            std::set<int64_t> lagrangeseqnum;
            for (int ic=nc-1; ic>=0 ; ic--) {
                TPZConnect &c = cel->Connect(ic);
                int clagrange = c.LagrangeMultiplier();
                int64_t ceqnum = c.SequenceNumber();
                if (ceqnum >= numinternalconnects) {
                    continue;
                }
                int64_t ceq = permute[ceqnum];
                if (clagrange == lagrange && ceq < maxeq) {
                    lagrangeseqnum.insert(ceqnum);
                }
            }
            std::set<int64_t>::reverse_iterator it;
            int64_t count = 0;
            for (it = lagrangeseqnum.rbegin(); it != lagrangeseqnum.rend(); it++) {
                int64_t ceq = permute[*it];
                ModifyPermute(permute, ceq, maxeq-count);
#ifdef PZ_LOG
                if(logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Switch ceq = " << ceq << " with maxeq = " << maxeq-count;
//                    sout << permute;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                count++;
            }
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Resequence for element " << el << std::endl;
                for (int ic=0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    c.Print(*this,sout);
                    if (c.SequenceNumber() < numinternalconnects) {
                        sout << "New seqnum = " << permute[c.SequenceNumber()] << std::endl;
                    }
                    else
                    {
                        sout << "Connect with restraint, seqnum " << c.SequenceNumber() << std::endl;
                    }
                }
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
        }
    }
#ifdef PZ_LOG2
    if(logger.isDebugEnabled())
    {
        for (int64_t el=0L; el<nel; el++) {
            TPZCompEl *cel = ElementVec()[el];
            if (!cel) {
                continue;
            }
            int nc = cel->NConnects();
            std::stringstream sout;
            sout << "Resequence for element " << el << std::endl;
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                c.Print(*this,sout);
                if (c.SequenceNumber() < numinternalconnects) {
                    sout << "New seqnum = " << permute[c.SequenceNumber()] << std::endl;
                }
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
    }
#endif
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Saddle permute old permutation ";
        for (int i = 0; i< permute.size(); i++) {
            int64_t jmax = i+10;
            if (jmax > permute.size()) {
                jmax = permute.size();
            }
            for (int64_t j=i; j< jmax; j++) {
                sout << permute[j] << ' ';
            }
            sout << std::endl;
            i+= 10;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Permute(permute );

}

void TPZCompMesh::GetEquationSetByMat(std::set<int64_t>& matidset, std::set<int64_t>& eqset) {
    
    if (!matidset.size()) {
        std::cout << "\nNo materials provided to TPZCompMesh::GetEquationSetByMat(). Returning..." << std::endl;
        return;
    }
    
    for (TPZCompEl* cel : this->ElementVec()) {
        if(!cel) continue;
        
        TPZGeoEl* gel = cel->Reference();
        if(!gel) continue; // It skips SubCompMeshes, ElementGroups, and CondensedElements with this check
                
        const int64_t gelmatid = gel->MaterialId();
        if(matidset.find(gelmatid) == matidset.end()) continue;
        
        const int ncon = cel->NConnects();
        for(int ic = 0 ; ic < ncon ; ic++){
            TPZConnect& con = cel->Connect(ic);
//            if(con.HasDependency() || con.IsCondensed()) continue; //CHECAR COM PHIL
            AddConnectEquationsToSet(con, eqset);
        }
    }
    
#ifdef PZDEBUG
    if(!eqset.size()){
        std::cout << "\n\n\t===> Warning! TPZCompMesh::GetEquationSetByMat() did not find any equations for the material set provided | ";
        std::cout << "matidset = ";
        for(auto const& id : matidset) std::cout << id << " ";
        std::cout << std::endl;
    }
#endif
}

void TPZCompMesh::AddConnectEquationsToSet(TPZConnect& con, std::set<int64_t>& eqset){
    const int64_t seqnum = con.SequenceNumber();
    if(seqnum < 0) DebugStop();
    const int64_t firsteq = this->Block().Position(seqnum);
    const int64_t blocksize = this->Block().Size(seqnum);
    for (int i = 0; i < blocksize; i++) {
        eqset.insert(firsteq+i);
    }
}

/// Modify the permute vector swapping the lagrangeq with maxeq and shifting the intermediate equations
void TPZCompMesh::ModifyPermute(TPZVec<int64_t> &permute, int64_t lagrangeq, int64_t maxeq)
{
    int64_t neq = permute.size();
#ifdef PZDEBUG
    if (lagrangeq < 0 || lagrangeq >= neq || maxeq < 0 || maxeq >= neq) {
        DebugStop();
    }
#endif
    // find the equation which maps to lagrangeq
    //int lagrangeqindex = permuteinv[lagrangeq];
    TPZVec<int64_t> accpermute(neq,0),input(permute);
    for (int64_t i=0; i<neq; i++) {
        accpermute[i] = i;
    }
    
    int64_t lagrangeqindex = lagrangeq;

    // this equation should never be sent forwards
    if (accpermute[lagrangeqindex] > lagrangeq) {
        DebugStop();
    }
    
    accpermute[lagrangeqindex] = maxeq;
    int64_t index = lagrangeqindex+1;
    while (index < neq && (accpermute[index] <= maxeq || accpermute[index] < index)) {
        accpermute[index] = accpermute[index]-1;
        index++;
    }
    for (int64_t i=0; i<neq; i++) {
        permute[i] = accpermute[input[i]];
    }
    
#ifdef PZDEBUG22
    {
        std::set<int64_t> acc;
        for (int64_t i=0; i<neq; i++) {
            acc.insert(permute[i]);
        }
        if (acc.size() != neq) {
            std::cout << "input " << input << std::endl;
            std::cout << "accpermute " << accpermute << std::endl;
            std::cout << "permute " << permute << std::endl;
            DebugStop();
        }
    }
#endif
}

/** @brief adds the connect indexes associated with base shape functions to the set */
void TPZCompMesh::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    int64_t nel = NElements();
    for (int64_t el=0; el<nel ; el++) {
        TPZCompEl *cel = ElementVec()[el];
        if (!cel) {
            continue;
        }
        cel->BuildCornerConnectList(connectindexes);
    }
}

TPZCompMesh * TPZCompMesh::CommonMesh(TPZCompMesh *mesh){
	
	TPZStack<TPZCompMesh *> s1, s2;
	int64_t pos1=0, pos2, comind;
	TPZCompMesh *father = FatherMesh();
	s1.Push(this);
	while (father){
		s1.Push((father));
		pos1++;
		father = s1[pos1]->FatherMesh();
	}
	pos2 = 0;
	s2.Push(mesh);
	father = mesh->FatherMesh();
	while (father){
		s2.Push(father);
		pos2++;
		father = s2[pos2]->FatherMesh();
	}
	if (s1[pos1] != s2[pos2]) return 0;
	comind=0; //The first mesh is common for all submeshes
	for (; pos1>=0 && pos2>=0; pos1--, pos2--) {
		if((s1[pos1])!=(s2[pos2])) {
			comind=pos1+1;
			return (s1[comind]);
		}
	}
	return (pos1 >=0 ) ? (s1[pos1+1]) : s2[pos2+1];
}

/** update the solution at the previous state with fSolution and
    set fSolution to the previous state */
template<class TVar>
void TPZCompMesh::UpdatePreviousState(TVar mult)
{
    if(fSolN.Rows() != fSolution.Rows() || fSolN.Cols() != fSolution.Cols())
    {
        fSolN.Redim(fSolution.Rows(),fSolution.Cols());
    }
    TPZFMatrix<TVar>& solN = fSolN;
    TPZFMatrix<TVar>& sol = fSolution;
    solN += mult*sol;
    fSolution = fSolN;
    int64_t nel = NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = Element(el);
        TPZSubCompMesh *subc = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subc) subc->UpdatePreviousState(mult);
    }
}


template class TPZRestoreClass<TPZCompMesh>;

/// extract the values corresponding to the connect from the vector
template<class TVar>
void TPZCompMesh::ConnectSolution(int64_t cindex, TPZCompMesh *cmesh, TPZFMatrix<TVar> &glob, TPZVec<TVar> &sol)
{
    int64_t seqnum = cmesh->ConnectVec()[cindex].SequenceNumber();
    int blsize = cmesh->Block().Size(seqnum);
    int position = cmesh->Block().Position(seqnum);
    sol.resize(blsize);
    for (int64_t i=position; i< position+blsize; i++) {
        sol[i-position] = glob(i,0);
    }
}

#define INSTANTIATE_METHODS(TVar) \
template \
void TPZCompMesh::UpdatePreviousState<TVar>(TVar); \
template \
void TPZCompMesh::SetElementSolution<TVar>(int64_t , TPZVec<TVar>&); \
template \
void TPZCompMesh::ConnectSolution<TVar>(int64_t , TPZCompMesh *, TPZFMatrix<TVar> &, TPZVec<TVar> &); \
template \
void TPZCompMesh::ProjectSolution<TVar>(TPZFMatrix<TVar> &);

INSTANTIATE_METHODS(STATE)
INSTANTIATE_METHODS(CSTATE)
#undef INSTANTIATE_METHODS
