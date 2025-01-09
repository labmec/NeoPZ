/**
 * @file
 * @brief Contains the implementation of the TPZSpStructMatrix methods. 
 */

#include "TPZSpStructMatrix.h"
#include "pzcmesh.h"
#include "TPZYSMPMatrix.h"
#include "TPZRenumbering.h"

#include "TPZTimer.h"

#ifdef USING_EIGEN
#include "TPZEigenSparseMatrix.h"
#endif

#ifdef USING_MKL
#include "TPZYSMPPardiso.h"
#endif

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar, class TPar>
TPZStructMatrix * TPZSpStructMatrix<TVar,TPar>::Clone(){
    return new TPZSpStructMatrix(*this);
}

template<class TVar, class TPar>
void TPZSpStructMatrix<TVar,TPar>::EndCreateAssemble(TPZBaseMatrix * mat){
    auto spMat = dynamic_cast<TPZFYsmpMatrix<TVar> *>(mat);
    spMat->ComputeDiagonal();
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSpStructMatrix<TVar,TPar>::Create(){
    int64_t neq = this->fMesh->NEquations();
    if(this->fMesh->FatherMesh()) {
//		PZError << "TPZSpStructMatrix should not be called with CreateAssemble for a substructure mesh\n";
    }
	
    /**
     *Longhin implementation
	 */
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> elgraphindex;
    //    int nnodes = 0;
    this->fMesh->ComputeElGraph(elgraph,elgraphindex);
    TPZMatrix<TVar> * mat = SetupMatrixData(elgraph,elgraphindex);
    return mat;
}

template<class TVar, class TPar>
TPZMatrix<TVar> *
TPZSpStructMatrix<TVar,TPar>::SetupMatrixData(TPZStack<int64_t> & elgraph,
                                              TPZVec<int64_t> &elgraphindex){
    
    const int64_t neq = this->fEquationFilter.NActiveEquations();
#ifdef USING_MKL
    TPZFYsmpMatrixPardiso<TVar> * mat = new TPZFYsmpMatrixPardiso<TVar>(neq,neq);
#elif USING_EIGEN
    TPZEigenSparseMatrix<TVar> * mat = new TPZEigenSparseMatrix<TVar>(neq,neq);
#else
    TPZFYsmpMatrix<TVar> *mat = new TPZFYsmpMatrix<TVar>(neq,neq);
    DebugStop();
#endif
    
    /**Creates a element graph*/
    TPZRenumbering graph;
    graph.SetElementsNodes(elgraphindex.NElements() -1 ,
                           this->fMesh->NIndependentConnects());
	
    TPZManVector<int64_t> nodegraph;
    TPZManVector<int64_t> nodegraphindex;
    /**
     *converts an element graph structure into a node graph structure
     *those vectors have size ZERO !!!
     */
    graph.ConvertGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
    /**vector sizes*/
    const int64_t nblock = nodegraphindex.NElements()-1;
    // number of values in the sparse matrix
    int64_t totalvar = 0;
    // number of equations
    int64_t totaleq = 0;
    for(int64_t i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive = this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
		totaleq += iblsize;
		const int64_t icfirst = nodegraphindex[i];
		const int64_t iclast = nodegraphindex[i+1];
		totalvar+=iblsize*iblsize;
		for(auto j=icfirst;j<iclast;j++) {
			const int64_t col = nodegraph[j];
			const int64_t colsize = this->fMesh->Block().Size(col);
			const int64_t colpos = this->fMesh->Block().Position(col);
            const int64_t numactive = this->fEquationFilter.NumActive(colpos, colpos+colsize);
            if (!numactive) {
                continue;
            }
			totalvar += iblsize*colsize;
		}
    }
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Number of equations " << totaleq << " number of nonzero s " << totalvar;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int64_t ieq = 0;
    // pos is the position where we will put the column value
    int64_t pos = 0;
	
    TPZManVector<int64_t,400> Eq(totaleq+1);
    TPZVec<int64_t> EqCol(totalvar);
    TPZVec<TVar> EqValue(totalvar);
    //lambda for avoid repeating code
    auto AddColEqs =
        [this,&EqCol,&EqValue,&pos](const int colsize, const int colpos)
        {
            TPZManVector<int64_t> destindices(colsize);
            for (int64_t i=0; i<colsize; i++) {
                destindices[i] = colpos+i;
            }
            this->fEquationFilter.Filter(destindices);
            for(auto jbleq=0; jbleq<destindices.size(); jbleq++) {
                EqCol[pos] = destindices[jbleq];
                EqValue[pos] = 0.;
                pos++;
            }
        };

    for(auto i=0;i<nblock;i++){
		const int64_t iblsize = this->fMesh->Block().Size(i);
		const int64_t iblpos = this->fMesh->Block().Position(i);
        const int64_t numactive =
            this->fEquationFilter.NumActive(iblpos, iblpos+iblsize);
        if (!numactive) {
            continue;
        }
        TPZManVector<int64_t> rowdestindices(iblsize);
        for (int64_t ij=0; ij<iblsize; ij++) {
            rowdestindices[ij] = iblpos+ij;
        }
        this->fEquationFilter.Filter(rowdestindices);
        // working equation by equation
        // rowdestindices contains the equation number of each element in the block number "i"
		for(auto ibleq=0; ibleq<rowdestindices.size(); ibleq++) {
            if (rowdestindices[ibleq] != ieq) {
                DebugStop();
            }
			Eq[ieq] = pos;
            //we need to add the diagonal block
    //since nodegraph does not contain the pointer to itself
            AddColEqs(this->fMesh->Block().Size(i),
                      this->fMesh->Block().Position(i));
			const int64_t icfirst = nodegraphindex[i];
			const int64_t iclast = nodegraphindex[i+1];
			for(auto j=icfirst;j<iclast;j++)
            {
                // col is the block linked to block "i"
				const int64_t col = nodegraph[j];
				const auto colsize = this->fMesh->Block().Size(col);
				const auto colpos = this->fMesh->Block().Position(col);
                const int64_t numactive =
                    this->fEquationFilter.NumActive(colpos, colpos+colsize);
                if (!numactive) {
                    continue;
                }
                AddColEqs(colsize,colpos);
			}
            //we need to sort the entries of ja otherwise pardiso will complain
            const auto firstentry = Eq[ieq];
            const auto nentries = pos-Eq[ieq];

            std::stable_sort(&EqCol[firstentry],
                             &EqCol[firstentry]+nentries);
            
			ieq++;
		}
    }
    Eq[ieq] = pos;
    if(pos != totalvar)
    {
        PZError<<__PRETTY_FUNCTION__<<'\n';
        PZError<<"nvars: "<<totalvar<<'\n';
        PZError<<"l pos: "<<pos<<'\n';
        DebugStop();
    }
    mat->SetData(std::move(Eq),std::move(EqCol),std::move(EqValue));
    return mat;
}

template<class TVar, class TPar>
int TPZSpStructMatrix<TVar,TPar>::ClassId() const{
    return Hash("TPZSpStructMatrix") ^
        TPZStructMatrixT<TVar>::ClassId() << 1 ^
        TPar::ClassId() << 2;
}

template<class TVar, class TPar>
void TPZSpStructMatrix<TVar,TPar>::Read(TPZStream& buf, void* context){
    TPZStructMatrix::Read(buf,context);
    TPar::Read(buf,context);
}

template<class TVar, class TPar>
void TPZSpStructMatrix<TVar,TPar>::Write(TPZStream& buf, int withclassid) const{
    TPZStructMatrix::Write(buf,withclassid);
    TPar::Write(buf,withclassid);
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"
#include "TPZStructMatrixOMPorTBB.h"

template class TPZSpStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSpStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZSpStructMatrix<STATE,TPZStructMatrixOMPorTBB<STATE>>;

// #ifndef USING_EIGEN
template class TPZSpStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSpStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSpStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;
template class TPZSpStructMatrix<CSTATE,TPZStructMatrixOMPorTBB<CSTATE>>;
// #endif
