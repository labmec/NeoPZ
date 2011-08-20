/** @file pzmaterial.h
 *
 * @brief Header file for abstract class TPZMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef PZMATERIALHPP
#define PZMATERIALHPP

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"
#include "pzsave.h"
#include "pzmaterialdata.h"

#include <iostream>
#include <string>
//#ifdef _AUTODIFF
//#include "fadType.h"
//#endif


class TPZBndCond;
class TPZMaterial;
class TPZMaterialData;
class TPZIntPoints;

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the TPZMaterial class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * TPZMaterial objects also need to implement the interface for post processing the results
 */
class  TPZMaterial : public TPZSaveable
{
private:
    int fId;
    
protected:
    void (*fForcingFunction)(TPZVec<REAL> &loc,TPZVec<REAL> &result);
    
    /// @brief Defines whether the equation context is linear solver or non linear
    /**
     * True means linear (default)
     * @since 08 oct 2010
     */
    bool fLinearContext;
    
public:
    
    static REAL gBigNumber;
    
    /** @brief Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
	 */
	/** 
	 * Upon return vectorindex
     *  contains the index of the material object within the
     *  vector
     */
    TPZMaterial(int id);
    
    /**
     * @brief Default constructor
     */
    TPZMaterial();
    
    /** @brief Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
	 */
	/**  Upon return vectorindex contains the index of the material
     *  object within the vector
     */
    TPZMaterial(const TPZMaterial &mat);
    
    virtual ~TPZMaterial();
    
    /** @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data);
    
    /** @brief returns the name of the material*/
    virtual std::string Name() { return "no_name"; }
    
    /** @brief returns the integrable dimension of the material*/
    virtual int Dimension() = 0;
    
    int Id() const { return fId; }
    void SetId(int id) {
        if(id == 0) {
            std::cout << "\n*** Material Id can't be ZERO! ***\n";
            std::cout << "*** This Will Be a Disaster!!! ***\n";
        }
        fId = id; }
    
    /** @brief returns the number of state variables associated with the material*/
    virtual int NStateVariables() = 0;
    
    /** @brief return the number of components which form the flux function*/
    virtual int NFluxes() {return 0;}
    
    /** @brief print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);
    
    /** @brief returns the variable index associated with the name*/
    virtual int VariableIndex(const std::string &name);
    
    /** @brief returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);
    
    /** @brief returns the solution associated with the var index based on
     * the finite element approximation*/
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
    
protected:
    /** @deprecated Deprecated interface for Solution method which must use material data. */
    virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
    
public:
    
    /** @brief compute the value of the flux function to be used by ZZ error
     * estimator*/
    virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol,
                      TPZFMatrix &DSol, TPZFMatrix &axes,
                      TPZVec<REAL> &flux) {}
    
    /** @brief Create an object TPZBndCond derived of TPZMaterial*/
    virtual TPZBndCond *CreateBC(TPZAutoPointer<TPZMaterial> &reference, int id, int typ, TPZFMatrix &val1,
                                 TPZFMatrix &val2);
    
    /** @name Metodos Contribute
	 * @{
	 */
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef) = 0;
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc) = 0;
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc);
    /** @} */
	
    /**Compute contribution to the energy at an integration point*/
    //      virtual void ContributeEnergy(TPZVec<REAL> &x,
    //			      TPZVec<FADFADREAL> &sol,
    //			      TPZVec<FADFADREAL> &dsol,
    //			      FADFADREAL &U,
    //			      REAL weight);
    
    //#endif
    
    
    //#ifdef _AUTODIFF
    
    /** Compute contribution of BC to the Energy*/
    //      virtual void ContributeBCEnergy(TPZVec<REAL> & x,
    //	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
    //	REAL weight, TPZBndCond &bc);
    
    //#endif
    /** 
	 * @brief Set a procedure as source function for the material.
	 * @param loc corresponds to the coordinate of the point where the source function is applied
	 */
    void SetForcingFunction(void (*fp)(TPZVec<REAL> &loc,
                                       TPZVec<REAL> &result))
    {
        fForcingFunction = fp;
    }
    
    virtual int HasForcingFunction() {return (fForcingFunction != 0);}
    
    /** 
     * @brief Get the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;
    
    /** Set the integration rule order based on the element
     *  @param p order of interpolation, its dimension and the characteristics
     *   of the material
     *
	 virtual void SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
	 int elPMaxOrder,
	 int elDimension);
	 */
    /**
	 * @brief Compute the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol,
                        TPZFMatrix &axes, TPZVec<REAL> &flux,
                        TPZVec<REAL> &uexact, TPZFMatrix &duexact,
                        TPZVec<REAL> &val){
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
    
    /**
     * @brief Returns the number of norm errors. Default is 3: energy, L2 and H1.
     */
    virtual int NEvalErrors() {return 3;}
    
    /** @brief To create another material of the same type*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial();
    
    /** @brief Read data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data);
    
    /**
     * @brief Create a copy of the material object and put it in the vector
     * which is passed on
     */
    virtual void Clone(std::map<int, TPZAutoPointer<TPZMaterial> > &matvec);
    
    /** @brief To return a numerical flux type to apply over the interfaces
     * of the elements*/
    virtual int FluxType() { return 2; }
    
    /** Factor to diffussive term*/
    //      virtual int IdBC(REAL *x) { return 5; }
    
    virtual void ContributeErrors(TPZMaterialData &data,
                                  REAL weight,
                                  TPZVec<REAL> &nk,
                                  int &errorid){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
    }
    
    /**
     * @brief Compute square of residual of the differential equation at one integration point.
     * @param X is the point coordinate (x,y,z)
     * @param sol is the solution vector
     * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
     */    
    virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<REAL> &sol, TPZFMatrix &dsol){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
        return -1.;
    }
    
    /**
     * @brief Unique identifier for serialization purposes
     */
    virtual int ClassId() const;
    
    /**
     * @brief Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /**
     * @brief Read the element data from a stream
     */
    virtual void Read(TPZStream &buf, void *context);
    
    /**
     * @brief Pushes a new entry in the context of materials with memory,
     * returning its index at the internal storage stack.
	 */
	/** To be implemented only in the proper materials.
     */
    virtual int PushMemItem(int sourceIndex = -1){ return -1; }
    
    /**
     * @brief Frees an entry in the material with memory internal history storage
     */
    virtual void FreeMemItem(int index){ return; }
    
    /** @brief Set fLinearContext attribute */
    void SetLinearContext(bool IsLinear);
    
    /** @brief Returns fLinearContext attribute */
    bool GetLinearContext() const;
    
};

inline bool TPZMaterial::GetLinearContext() const{
    return fLinearContext;
}

/** @brief Extern variable - Vector of force values */
extern TPZVec< void(*) ( TPZVec<REAL> &, TPZVec<REAL>& ) > GFORCINGVEC;

#endif

