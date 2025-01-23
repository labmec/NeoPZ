/**
 * @file pzadmchunk.h
 * @brief Declarates the TPZBlock<REAL>class which implements block matrices.
 */

#ifndef _TBLOCKHH_
#define _TBLOCKHH_

#include "pzmanvector.h"
#include "pzreal.h"
#include "TPZSavable.h"
#include "Hash/TPZHash.h"

template<class TVar>
class TPZFMatrix;

class TPZBaseMatrix;

/**
 * @brief Implements block matrices. \ref matrixutility "Matrix utility"
 * @author Misael Luis Santana Mandujano
 * @since 12/1994
 * @ingroup matrixutility
 */
class TPZBlock : public TPZSavable
{
public:
    TPZBlock() : TPZRegisterClassId(&TPZBlock::ClassId), fBlock(), fpMatrix(nullptr)
    {
        
    }
	/**
	 * @brief For each elements on matrix a size 1 block is created
	 * @param matrix_to_represent Indicates which matrix is to be represented
	 * @param num_of_blocks Indicates number of blocks
	 * @param initial_blocks_size Indicates initial block size, default value is 1
	 */
	TPZBlock(TPZBaseMatrix *const matrix_to_represent,const int num_of_blocks = 0,
			 const int initial_blocks_size = 1 );
	
	/**
	 * @brief Copy constructor
	 * @param bl New object is created based on bl
	 */
	TPZBlock(const TPZBlock &bl);
	
	/** @brief Simple Destrutor */
	virtual ~TPZBlock();
	
    /// Reset the data structure
    void Reset() 
    {
        fBlock.Resize(0);
        fpMatrix = nullptr;
    }

	/**
     * @brief Changes pointer to other
     * @param other New matrix to be pointed to
	 */
	void SetMatrix(TPZBaseMatrix *const other)
    {
        fpMatrix = other;
    }
	
	/** @brief Returns a pointer to current matrix */
    template <class TVar>
	TPZFMatrix<TVar> *Matrix(){
		if(auto tmp = dynamic_cast<TPZFMatrix<TVar>*>(fpMatrix); tmp){
            return tmp;
        }
        else{
            PZError<<"Incompatible matrix type in ";
            PZError<<__PRETTY_FUNCTION__<<std::endl;
            DebugStop();
            return nullptr;
        }
	}

    template <class TVar>
	const TPZFMatrix<TVar> *Matrix()const{
		if(auto tmp = dynamic_cast<const TPZFMatrix<TVar>*>(fpMatrix); tmp){
            return tmp;
        }
        else{
            PZError<<"Incompatible matrix type in ";
            PZError<<__PRETTY_FUNCTION__<<std::endl;
            DebugStop();
            return nullptr;
        }
	}
	
	/**
     * @brief Sets number of blocks on diagonal matrix
     * @param num_of_blocks Number of blocks
	 */
	int SetNBlocks(const int num_of_blocks );
	
	/**
     * @brief Modifies existing block dimensions or creates a new block with given index
     * @param index Given index to be redimensioned or created
     * @param dim New dimension
     * @param pos New position
	 */
	int Set(const int index,const int dim,const int pos = -1 );
	
	/**
     * @brief Computes blocks sequence
     * @param dimensions Contains blocks sequence
	 */
	int SetAll( TPZVec<int> & dimensions );

	/**
     * @brief Resequences blocks positioning
     * @param start Starting position
	 */
	int Resequence(const int start=0);
	
	/**
     * @brief Removes a block
     * @param index Index of the block to be removed
	 */
	int Remove(const int index );
	
	/**
     * @brief Verifies if blocks are sequential and does not overcome matrix size
	 */
	int Verify() const;
	
    /// Return the index in the blocked matrix
    int64_t Index(const int64_t block_row, const int r) const;
    
    /// Return the row-column index for the row and column block
    std::pair<int64_t, int64_t> at(const int block_row,const int block_col,const int r,const int c) const
    {
        return std::pair<int64_t, int64_t>(Index(block_row, r),Index(block_col,c));
    }
	
public:
    
	TPZBlock&operator=(const TPZBlock& ); 		
	
    void PrintStructure(std::ostream &out = std::cout);
	
	/** @brief Returns the max number of blocks on diagonal */
	int64_t MaxBlockSize() const {return fBlock.NElements();}
	/** @brief Returns number of blocks on diagonal */
	int64_t NBlocks() const {return fBlock.NElements();}
	
	/**
     * @brief Returns block dimension
     * @param block_diagonal Inquired block_diagonal
	 */
	int Size(const int64_t block_diagonal) const { return fBlock[block_diagonal].dim; }
	
	/**
     * @brief Returns the position of first element block dependent on matrix diagonal
     * @param block_diagonal Inquired block_diagonal
	 */
	int Position(const int64_t block_diagonal) const { return fBlock[block_diagonal].pos;}
	
	/** @brief Returns matrix dimension pointed by block */
	int Dim() const {return fBlock.NElements() ? fBlock[fBlock.NElements()-1].pos+fBlock[fBlock.NElements()-1].dim : 0; }

    template <class TVar>
    int PutBlock(const int bRow, const int bCol,
                 const TPZFMatrix<TVar> &block);

    template <class TVar>
    int GetBlock(const int bRow, const int bCol,
                 TPZFMatrix<TVar> &block) const;
        /** @brief returns the unique identifier for reading/writing objects to
         * streams
         */

        int ClassId() const override;

        /** @brief Save the element data to a stream */
        void Write(TPZStream &buf, int withclassid) const override;

        /** @brief Read the element data from a stream */
        void Read(TPZStream &buf, void *context) override;

      private:
	
	/**
     * @struct TNode
     * @brief Defines a node
	 */
	class TNode : public TPZSavable {
            public :
                
		int pos; /**< Position of node */
		int dim; /**< Dimension of node */
		
		TNode() {
			pos=0;
			dim=0;
		}
                
        int ClassId() const override{
            return Hash("TNode") ^ ClassIdOrHash<TPZBlock>() << 1;
        }
                
		void Read(TPZStream &buf, void *context) override;
		void Write(TPZStream &buf, int withclassid) const override;
	};
	
	/** @brief Nodes vector */
	TPZManVector<TNode>    fBlock;
	/**  @brief Pointer to TPZMatrix */
	TPZBaseMatrix *fpMatrix;
	static REAL gZero;//zero
	
};

extern template TPZFMatrix<float> * TPZBlock::Matrix();
extern template TPZFMatrix<double> * TPZBlock::Matrix();
extern template TPZFMatrix<long double> * TPZBlock::Matrix();

extern template TPZFMatrix<std::complex<float> > * TPZBlock::Matrix();
extern template TPZFMatrix<std::complex<double> > * TPZBlock::Matrix();
extern template TPZFMatrix<std::complex<long double> > * TPZBlock::Matrix();

extern template const TPZFMatrix<float> * TPZBlock::Matrix() const;
extern template const TPZFMatrix<double> * TPZBlock::Matrix() const;
extern template const TPZFMatrix<long double> * TPZBlock::Matrix() const;

extern template const TPZFMatrix<std::complex<float> > * TPZBlock::Matrix() const;
extern template const TPZFMatrix<std::complex<double> > * TPZBlock::Matrix() const;
extern template const TPZFMatrix<std::complex<long double> > * TPZBlock::Matrix() const;



extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<float> &);
extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<double> &);
extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<long double> &);
extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<float>> &);
extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<double>> &);
extern template
int TPZBlock::PutBlock(const int, const int,
                       const TPZFMatrix<std::complex<long double>> &);


extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<float> &) const;
extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<double> &) const;
extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<long double> &) const;
extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<float>> &) const;
extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<double>> &) const;
extern template
int TPZBlock::GetBlock(const int, const int,
                       TPZFMatrix<std::complex<long double>> &) const;


#endif
