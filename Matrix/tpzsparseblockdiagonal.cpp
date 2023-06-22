
/**
 * @file
 * @brief Contains the implementation of the TPZSparseBlockDiagonal methods.
 */

#include "tpzsparseblockdiagonal.h"
#include "pzfmatrix.h"

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.StrMatrix");
#endif

using namespace std;

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal() : TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId)
{
}
template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(const TPZVec<int64_t> &blockgraph,
													 const TPZVec<int64_t> &blockgraphindex,
													 const int64_t rows) :
	TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId), fBlock(blockgraph), fBlockIndex(blockgraphindex)
{
	int64_t numbl = blockgraphindex.NElements()-1;
	this->fBlockSize.Resize(numbl);
	int64_t ibl;
	for(ibl=0; ibl<numbl; ibl++)
	{
		this->fBlockSize[ibl] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
		fGlobalBlockIndex[ibl] = ibl;
	}
	this->Initialize(this->fBlockSize);
	//initialize from parent class will set fRow and fCol as the number of equations
	this->fRow = rows;
	this->fCol = rows;
}

template<class TVar>
TPZSparseBlockDiagonal<TVar>::TPZSparseBlockDiagonal(const TPZVec<int64_t> &blockgraph,
													 const TPZVec<int64_t> &blockgraphindex,
													 const int64_t rows, const int color,
													 const TPZVec<int> &colors) :
	TPZRegisterClassId(&TPZSparseBlockDiagonal::ClassId)
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled()){
		LOGPZ_DEBUG(logger, "Constructor of TPZSparseBlockDiagonal");
	}
#endif
	const int64_t numbl = blockgraphindex.NElements()-1;
#ifdef PZ_LOG
	if(numbl != colors.size()){
		PZError<<__PRETTY_FUNCTION__
					 <<"\nInvalid input! number of blocks "<<numbl
					 <<"\nSize of color vec: "<<colors.size()
					 <<std::endl;
		DebugStop();
	}
#endif
	this->fBlockSize.Resize(numbl);
	int64_t ibl,iblcount,graphsize = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			this->fBlockSize[iblcount] = blockgraphindex[ibl+1]-blockgraphindex[ibl];
			graphsize += this->fBlockSize[iblcount++];
		}
	}
	this->fBlockSize.Resize(iblcount);
	fBlock.Resize(graphsize);
	fBlockIndex.Resize(iblcount+1);
	fBlockIndex[0] = 0;
	for(ibl=0, iblcount=0; ibl<numbl; ibl++)
	{
		if(colors[ibl]==color) 
		{
			fGlobalBlockIndex[ibl] = iblcount;
			int64_t first = blockgraphindex[ibl];
			int64_t last = blockgraphindex[ibl+1];
			int64_t firstcp = fBlockIndex[iblcount];
			this->fBlockIndex[iblcount+1] = firstcp+this->fBlockSize[iblcount];
			//      int lastcp = fBlockIndex[iblcount+1];
			int64_t ieq,ieqcp;
			for(ieq=first,ieqcp=firstcp; ieq<last; ieq++,ieqcp++)
			{
				fBlock[ieqcp] = blockgraph[ieq];
			}
			iblcount++;
		}
	}
#ifdef PZDEBUG
	std::set<int64_t> eqset;
	for(auto eq : fBlock){
		eqset.insert(eq);
	}

	if(eqset.size() != fBlock.NElements()){
		std::cout<<__PRETTY_FUNCTION__
						 <<"\ncoloring is not correct!\n"
						 <<"set of equations has  "<<eqset.size()<<" eqs\n"
						 <<"fBlock has size "<<fBlock.NElements()
						 <<std::endl;
	}
#endif
	
	this->Initialize(this->fBlockSize);
	this->fRow = rows;
	this->fCol = rows;
}

template<class TVar>
const TVar TPZSparseBlockDiagonal<TVar>::Get(const int64_t row, const int64_t col) const
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return (TVar)0;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return (TVar)0;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
const TVar TPZSparseBlockDiagonal<TVar>::GetVal(const int64_t row, const int64_t col) const
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return (TVar) 0;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return (TVar) 0;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template <class TVar>
int TPZSparseBlockDiagonal<TVar>::Put(const int64_t row, const int64_t col, const TVar& value)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::PutVal(const int64_t row, const int64_t col, const TVar& value)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return -1;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return -1;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	this->fStorage[this->fBlockPos[rblock]+pos] = value;
	return 0;
}
template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::operator ( )(const int64_t row, const int64_t col)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
int TPZSparseBlockDiagonal<TVar>::Substitution(TPZFMatrix<TVar>* B) const
{
	TPZFNMatrix<1000,TVar > BG(fBlock.NElements(),B->Cols());
	Gather(*B,BG);
	int result = TPZBlockDiagonal<TVar>::Substitution(&BG);
	B->Zero();
	ScatterAdd(BG,*B);
	return result;
	
}

template<class TVar>
TVar& TPZSparseBlockDiagonal<TVar>::s(const int64_t row, const int64_t col)
{
	int64_t rblock,rblockindex,cblock,cblockindex;
	FindBlockIndex(row,rblock,rblockindex);
	if(rblock == -1) return this->gZero;
	FindBlockIndex(col,cblock,cblockindex);
	if(cblock != rblock) return this->gZero;
	int64_t pos = rblockindex + cblockindex*this->fBlockSize[rblock];
	return this->fStorage[this->fBlockPos[rblock]+pos];
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Print(const char* message, std::ostream& out, const MatrixOutputFormat format) const
{
	TPZBlockDiagonal<TVar>::Print(message, out, format);
	if(format == EFormatted)
	{
		out << "Equations for each block " << endl;
		int64_t nbl = fBlockIndex.NElements()-1;
		int64_t ibl;
		for(ibl = 0; ibl<nbl ; ibl++)
		{
			int64_t first = fBlockIndex[ibl];
			int64_t last = fBlockIndex[ibl+1];
			out << "Block " << ibl << " : ";
			int64_t i;
			for(i=first; i<last; i++) out << fBlock[i] << " ";
			out << endl;
		}
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::BuildFromMatrix(TPZMatrix<TVar>& matrix)
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled()){
		LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::BuildFromMatrix");
	}
#endif
	TPZManVector<int64_t> indices;
	TPZFNMatrix<10000,TVar> submat(0,0);
	int64_t ibl,nbl = fBlockIndex.NElements()-1;
	for(ibl=0; ibl<nbl; ibl++)
	{
		int64_t nel = this->fBlockSize[ibl];
		indices.Resize(nel);
		submat.Resize(nel,nel);
		int64_t iel,first = fBlockIndex[ibl];
		for(iel=0; iel<nel; iel++) indices[iel] = fBlock[first+iel];
		matrix.GetSub(indices,submat);
		this->SetBlock(ibl,submat);
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::GetBlock(int64_t i, TPZFMatrix<TVar>& block)
{
    TPZBlockDiagonal<TVar>::GetBlock(i, block);
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::MultAdd(const TPZFMatrix<TVar>& x, const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z, const TVar alpha, const TVar beta, const int opt) const
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled()){
		LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::MultAdd");
	}
#endif
	TPZFNMatrix<1000,TVar> xsc(0,0),ysc(0,0,0.),zsc(0,0);
	xsc.Resize(this->fBlock.NElements(),x.Cols());
	if(abs(beta) != 0.) ysc.Resize(fBlock.NElements(),y.Cols());
	zsc.Resize(fBlock.NElements(),z.Cols());
	Gather(x,xsc);
	if(abs(beta) != 0.) Gather(y,ysc);
	const auto rows = this->fRow;
	/*sorry for non-constness*/
	TPZSparseBlockDiagonal<TVar> *other = (TPZSparseBlockDiagonal<TVar> *)this;
	other->fRow = other->fCol = this->fBlock.NElements();
	TPZBlockDiagonal<TVar>::MultAdd(xsc, ysc, zsc, alpha, beta, opt);
	z.Zero();
	ScatterAdd(zsc,z);
	other->fRow = other->fCol = rows;
}

/*!
 \fn TPZSparseBlockDiagonal::FindBlockIndex(int glob, int &block, int &blockind)
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::FindBlockIndex(int64_t glob, int64_t &block, int64_t &blockind) const
{
    int64_t numbl = fBlockIndex.NElements()-1;
    int64_t ieq,ibl;
    for(ibl = 0; ibl<numbl; ibl++)
    {
		for(ieq = fBlockIndex[ibl];ieq<fBlockIndex[ibl+1];ieq++)
		{
			if(fBlock[ieq] == glob)
			{
				block = ibl;
				blockind = ieq-fBlockIndex[ibl];
				return;
			}
		}
    }
    block = -1;
    blockind = -1;
}

/*!
 \fn TPZSparseBlockDiagonal::Scatter(TPZFMatrix<>&in, TPZFMatrix<>&out) const
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::ScatterAdd(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const
{
    int64_t neq = fBlock.NElements();
    int64_t nc = in.Cols();
    int64_t ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(fBlock[ieq],ic) += in.GetVal(ieq,ic);
    }
}

/*!
 \fn TPZSparseBlockDiagonal::Gather(TPZFMatrix<>&in, TPZFMatrix<>&out) const
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::Gather(const TPZFMatrix<TVar> &in, TPZFMatrix<TVar> &out) const
{
    int64_t neq = fBlock.NElements();
    int64_t nc = in.Cols();
    int64_t ieq,ic;
    for(ic=0; ic<nc; ic++)
    {
		for(ieq=0; ieq<neq; ieq++) out(ieq,ic) = in.GetVal(fBlock[ieq],ic);
    }
}

/**
 * Updates the values of the matrix based on the values of the matrix
 */
template<class TVar>
void TPZSparseBlockDiagonal<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled()){
		LOGPZ_DEBUG(logger, "TPZSparseBlockDiagonal::UpdateFrom");
	}
#endif
	if(!mat) 
	{
		cout << __PRETTY_FUNCTION__ << " called with zero argument\n";
		return;
	}
	this->fDecomposed = ENoDecompose;
	int64_t nblock = this->fBlockSize.NElements();
	int64_t b,bsize,pos;
	TPZManVector<int64_t,1000> indices;
	for(b=0; b<nblock; b++) {
		bsize = this->fBlockSize[b];
		if(bsize){
			indices.Resize(bsize);
			int64_t r;
			pos = this->fBlockPos[b];
			for(r=0; r<bsize; r++) indices[r] = fBlock[fBlockIndex[b]+r]; 
			auto &block = *(this->fBlockMats[b]);
			mat->GetSub(indices,block);
		}
	}
	
}

template<class TVar>
int64_t TPZSparseBlockDiagonal<TVar>::HasBlock(const int64_t global) const
{
	if(fGlobalBlockIndex.count(global)){
		return fGlobalBlockIndex.at(global);
	}else{
		return -1;
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::GetBlockList(TPZVec<int64_t> &loc_ind, TPZVec<int64_t> &glob_ind) const
{
	const int nbl = fGlobalBlockIndex.size();
	loc_ind.Resize(nbl);
	glob_ind.Resize(nbl);
	int i = 0;
	for(auto it : fGlobalBlockIndex){
		glob_ind[i] = it.first;
		loc_ind[i++] = it.second;
	}
}

template<class TVar>
void TPZSparseBlockDiagonal<TVar>::GetBlockEqs(const int64_t global_ibl, TPZVec<int64_t> &eqs) const{
	const auto loc_ibl = HasBlock(global_ibl);
	if (loc_ibl < 0){eqs.Resize(0);}
	else{
		const auto first = fBlockIndex[loc_ibl];
		const auto last = fBlockIndex[loc_ibl+1];
		const auto sz = last-first;
		eqs.Resize(sz);
		for(int i = first; i < last; i++){
			eqs[i-first] = fBlock[i];
		}
	}
}

template <class TVar>
int TPZSparseBlockDiagonal<TVar>::ClassId() const{
    return Hash("TPZSparseBlockDiagonal") ^ TPZBlockDiagonal<TVar>::ClassId() << 1;
}
template class TPZSparseBlockDiagonal<float>;
template class TPZSparseBlockDiagonal<double>;
template class TPZSparseBlockDiagonal<long double>;

template class TPZSparseBlockDiagonal<std::complex<float> >;
template class TPZSparseBlockDiagonal<std::complex<double> >;
template class TPZSparseBlockDiagonal<std::complex<long double> >;
