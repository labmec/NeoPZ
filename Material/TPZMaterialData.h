/**
 * @file
 * @brief Contains the TPZMaterialData class which implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods.
 */

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "TPZShapeData.h"

/// @brief Determine the type of contribution to the stiffness matrix
enum MSpaceConfig {ENotActive, ETestTrial, ETest, ETrial};

/**
 * @ingroup material
 * @brief This class implements a type-agnostic interface between TPZCompEl::CalcStiff and Contribute methods of the materials. \n
 * It requests to the material which attributes must be computed by the computational element and trigger their computation.
 */
class TPZMaterialData : public TPZShapeData, virtual public TPZSavable {

public:
    //! Default constructor.
    TPZMaterialData();
    //! Copy constructor.
    TPZMaterialData(const TPZMaterialData&) = default;
    //! Move constructor.
    TPZMaterialData(TPZMaterialData&&) = default;
    //! Destructor.
    virtual ~TPZMaterialData() = default;
    //! Copy assignment operator
    TPZMaterialData& operator=(const TPZMaterialData&) = default;
    //! Move assignment operator
    TPZMaterialData& operator=(TPZMaterialData&&) = default;

    /** @brief Computes the flux divergence values for HDiv approximation spaces*/
    [[deprecated("This call is useless")]] virtual void ComputeFunctionDivergence() = 0;
    //@{
    //! Read and Write methods
    int ClassId() const override;
    
    void Write(TPZStream &buf, int withclassid) const override;
    
    void Read(TPZStream &buf, void *context) override;
    //@}
    //! Shape function type as a string.
    std::string ShapeFunctionType() const;
    /** @brief Set all flags at once (except FAD directions)*/
    void SetAllRequirements(bool set);
    /** @brief Prints the data */
    void Print(std::ostream &out) const;
    /** @brief Prints the data in a format suitable for Mathematica */
    void PrintMathematica(std::ostream &out) const;    

    /** @brief Resize the solution vectors to their appropriate size
        @param[in] nSol number of solutions
        @param[in] uLen size of each solution
        @param[in] duRow number of rows of the solutions derivative
        @param[in] duCol number of columns of the solutions derivative*/
    virtual void SetSolSizes(const int nSol, const int uLen,
                             const int duRow, const int duCol) = 0;
    static constexpr int MatDataNumPhi{20};
    static constexpr int MatDataNumDPhi{60};
    static constexpr int MatDataNumDir{110};
    static constexpr int MatDataDimSol{10};
    static constexpr int MatDataNumSol{20};
    //! Auxiliary attribute for collapsed HDiv elements
    void* fUserData;
    /** @name Flags
     * @brief Flags indicating which attributes should be available when computing the FEM matrix.*/
    /** @{ */
    /** Whether the solution is needed for computing the FEM matrix */
    bool fNeedsSol{false};
    /** Whether the neighbour's solution is needed for computing the FEM matrix */
    bool fNeedsNeighborSol{false};
    /** Whether the element size is needed for computing the FEM matrix */
    bool fNeedsHSize{false};
    /** Whether the neighbour's center coordinate is needed for computing the FEM matrix */
    bool fNeedsNeighborCenter{false};
    /** Whether FAD directions are needed for computing the FEM matrix. Can be used for HDiv and HCurl approximation spaces */
    bool fNeedsDeformedDirectionsFad{false};
    /** Whether the normal vector is needed for computing the FEM matrix */
    bool fNeedsNormal{false};
    /** Whether the material has an active approximation space. */
    MSpaceConfig fActiveApproxSpace{ETestTrial};
    /** @} */
    
    /** @name Data
     * @brief Attributes to be computed in CalcStiff */
    /** @{ */
    
    /// Vector of shapefunctions (format is dependent on the value of shapetype)
    TPZFNMatrix<MatDataNumPhi, REAL> phi;
    /** @brief Values of the derivative of the shape functions in the deformed element.
        @note This derivative is calculated in the element axes. For obtaining the
        derivative in the xyz-coordinate system, one should call
        TPZAxesTools<REAL>::Axes2XYZ(dphix, dphix_xyz, data.axes)
        in the contribute method.
     */
    TPZFNMatrix<MatDataNumDPhi, REAL> dphix;
    /// Values of the divergence of the shape functions in the mapped element (only applicable to H(div) spaces)    
    TPZFNMatrix<MatDataNumPhi, REAL> divphi;
    /// Values of the curl of the shape functions in the mapped element (only applicable to H(curl) spaces)
    TPZFNMatrix<MatDataNumDPhi, REAL> curlphi;
    /// Axes between the R3 space and the element's coordinates. Used for the derivatives of the shape functions
    TPZFNMatrix<9,REAL> axes;
    /// Value of the jacobian at the integration point
    TPZFNMatrix<9,REAL> jacobian;
    /// Value of the inverse of the jacobian at the integration point
    TPZFNMatrix<9,REAL> jacinv;
    /// Normal to the element at the integration point
    TPZManVector<REAL,3> normal;
    /// Value of the coordinate at the integration point
    TPZManVector<REAL,3> x;
    /// Value of the coordinate at the integration point corresponding to the x-parametric coordinate (master element)
    TPZManVector<REAL,3> xParametric;
    /// Maximum polinomial order of the shape functions
    int p{-1};
    /// Measure of the size of the element
    REAL HSize = 0.;
    /// Determinant of the jacobian
    REAL detjac = 0.;
    /// Value of the coordinate at the center of the element
    TPZManVector<REAL,3> XCenter;
    // Id of associated geometric element
    int gelElId{-1};    
    /// Correspondence between direction vector index and index of the shape functions. Used for H(div) and H(curl) approximation spaces.
    TPZManVector<std::pair<int,int64_t>, MatDataNumPhi > fVecShapeIndex;
    /// Directions on the deformed element. Used for H(div) and H(curl) approximation spaces.
    TPZFNMatrix<MatDataNumDir> fDeformedDirections;
    /** @} */
    //Used for HDiv pressure
    int numberdualfunctions{0};
    /// Directions on the deformed element using Fad
    TPZFMatrix<Fad<REAL>> fDeformedDirectionsFad;
    

    /** @brief Index of the current integration point being evaluated **/
    /** Needed for materials with memory **/
    int intLocPtIndex = -1;
    
    /** @brief global point index */
    int intGlobPtIndex = -1;
    
    /** @brief Number of points in the integration rule */
    int NintPts = 0;
    
protected:
    //! Dummy function to force this class to be abstract.
    virtual bool HasSol()=0;
};
#endif

