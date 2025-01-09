/**
 * @file
 * @brief Contains the TPZFYsmpMatrix class which implements a non symmetric sparse matrix. \n
 * Purpose: Defines operations on non-symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */

#ifndef EIGENYSMPMATH
#define EIGENYSMPMATH

#include "pz_config.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSYSMPMatrix.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define USEACCELERATESUPPORT
#if defined(MACOSX) && defined(USEACCELERATESUPPORT)
#include <Eigen/AccelerateSupport>
#endif

/**
 * @brief Implements a non symmetric sparse matrix (Yale Sparse Matrix Storage). \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Defines operations on general sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
template<class TVar>
class TPZEigenSparseMatrix : public TPZSYsmpMatrix<TVar> {
    
    public :
    typedef Eigen::SparseMatrix<TVar,0,int64_t> SpMat; // declares a column-major sparse matrix type of double
    typedef Eigen::SparseMatrix<TVar,0> SpMatInt; // declares a column-major sparse matrix type of double with int indexer. Needed for mac accelerate
    
    typedef Eigen::Map<SpMat> EigenSparse;
    typedef Eigen::Map<SpMatInt> EigenSparseInt;
    typedef Eigen::SimplicialLLT<SpMat> EigenCholesky;
    
    typedef Eigen::SparseLU<SpMat> EigenLU;
    
#if defined(MACOSX) && defined(USEACCELERATESUPPORT)
    typedef Eigen::AccelerateLDLT<SpMatInt,Eigen::Lower | Eigen::Symmetric> EigenLDLT;
#else
    typedef Eigen::SimplicialLDLT<SpMat> EigenLDLT;
#endif
    

private:
    EigenSparse *fEigenMatrix = 0;
    EigenCholesky *fCholesky = 0;
    EigenLDLT *fLDLT = 0;
    EigenLU *fLU = 0;
    
public:
    
    TPZEigenSparseMatrix();
    TPZEigenSparseMatrix(const int64_t rows,const int64_t cols );
    TPZEigenSparseMatrix(const TPZEigenSparseMatrix<TVar>&) = default;
    TPZEigenSparseMatrix(const TPZSYsmpMatrix<TVar>&copy) : TPZSYsmpMatrix<TVar>(copy) {

    }
    TPZEigenSparseMatrix(TPZEigenSparseMatrix<TVar>&&) = default;
    
    
    TPZEigenSparseMatrix &operator=(const TPZEigenSparseMatrix<TVar> &copy) = default;
    TPZEigenSparseMatrix &operator=(TPZEigenSparseMatrix<TVar> &&copy) = default;
    
    inline TPZEigenSparseMatrix<TVar>*NewMatrix() const override {
        return new TPZEigenSparseMatrix<TVar>();
        }
    CLONEDEF(TPZEigenSparseMatrix)
    
    virtual ~TPZEigenSparseMatrix();
    
    /** @brief Creates a copy from another TPZFYsmpMatrix*/
    void CopyFrom(const TPZMatrix<TVar> *  mat) override
    {
        auto *from = dynamic_cast<const TPZEigenSparseMatrix<TVar> *>(mat);
        if (from) {
            *this = *from;
        }
        else
        {
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"\nERROR: Called with incompatible type\n.";
            PZError<<"Aborting...\n";
            DebugStop();
        }
    }
    
    
    /**@name Arithmetics */
    /** @{ */
    TPZEigenSparseMatrix<TVar> operator+(const TPZEigenSparseMatrix<TVar> & A) const;
    TPZEigenSparseMatrix<TVar> operator-(const TPZEigenSparseMatrix<TVar> & A) const;
    TPZEigenSparseMatrix<TVar> operator*(const TVar alpha) const;
    TPZEigenSparseMatrix<TVar> &operator+=(const TPZEigenSparseMatrix<TVar> &A );
    TPZEigenSparseMatrix<TVar> &operator-=(const TPZEigenSparseMatrix<TVar> &A );
    TPZMatrix<TVar> &operator*=(const TVar val) override;
    
    
    
    /** @brief Print the matrix along with a identification title */
    virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
    
    /**
     * @name Solvers
     * @brief Linear system solvers. \n
     */
    /** For symmetric decompositions lower triangular matrix is used. \n
     * Solves a system A*X = B returning X in B
     */
    
    /**
     * @brief Add a contribution of a stiffness matrix
     * putting it on destination indexes position
     */
    
    
    /// this is a class that doesn't implement direct decompostion
    /** @brief decompose the system of equations acording to the decomposition
     * scheme */
    virtual int Decompose(const DecomposeType dt) override;
    /**
     * @brief Solves the linear system using Direct methods
     * @param F The right hand side of the system and where the solution is stored.
     * @param dt Indicates type of decomposition
     */
    virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) override;
    virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const override;
    
    
    /// Reimplements the multadd function using eigen optimized operations
    /// Multadd is:
    /// z = beta * y + alpha * opt(this)*x
    /// z and x cannot share storage
    /// - Parameters:
    ///   - x: vector that the sparse matrix opt(this) will multiply
    ///   - y: a vector that multiplied by alpha will sum to opt(this)*x*
    ///   - z: the result of the operation
    virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                         const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const override;
    
    int ClassId() const override;
    
    
    
    
protected:
    void CheckTypeCompatibility(const TPZMatrix<TVar>*aPtr,
                                const TPZMatrix<TVar>*bPtr) const override;
    
protected:
    
    /**
     * @brief Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
     */
    /**
     * -fIA will contain the initial positions for all the equations
     * -fJA will contain (-1) on all its positions
     * -fA will contain 0 on all its value
     */
    void InitializeData();
};


#endif
