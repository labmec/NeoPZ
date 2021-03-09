/**
 * @file pzmaterial.h
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
#include "TPZSavable.h"
#include "pzmaterialdata.h"
#include "pzfunction.h"
#include "pzcompel.h"

#include <iostream>
#include <string>

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
class  TPZMaterial : public virtual TPZSavable
{
private:
    int fId;
    
protected:
   
	/** @brief Pointer to forcing function, it is the right member at differential equation */
    TPZAutoPointer<TPZFunction<STATE> > fForcingFunction;
	
	/** @brief Pointer to exact solution function, needed to calculate exact error */
    TPZAutoPointer<TPZFunction<STATE> > fForcingFunctionExact;
	
public:
	/** @brief Pointer to time dependent forcing function, it is the right member at differential equation */
    TPZAutoPointer<TPZFunction<STATE> > fTimeDependentForcingFunction;
protected:
	/** @brief Pointer to time dependent exact solution function, needed to calculate exact error */
    TPZAutoPointer<TPZFunction<STATE> > fTimedependentFunctionExact;
    
    /** @brief Pointer to bc forcing function, it is a variable boundary condition at differential equation */
    TPZAutoPointer<TPZFunction<STATE> > fBCForcingFunction;
    
    /** @brief Pointer to time dependent bc forcing function, it is a variable boundary condition at differential equation */
    TPZAutoPointer<TPZFunction<STATE> > fTimedependentBCForcingFunction;    
    

    /**
	 * @brief Defines whether the equation context is linear solver or non linear
     * @return True means linear (default)
     * @since 08 oct 2010
     */
    bool fLinearContext;
    /** @brief Defines the number of load cases generated by this material */
    /**
     * The number of load cases will determine the number of columns of the right hand side vector
     * this variable defaults to one
     */
    int fNumLoadCases;
    
    /** @brief indicates which solution should be used for post processing */
    int fPostProcIndex;
    
public:
    /** @brief Big number to penalization method, used for Dirichlet conditions */
    static REAL gBigNumber;
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMaterial(int id);
    
    /** @brief Default constructor */
    TPZMaterial();
    
    /** @brief Creates a material object based on the referred object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */
    TPZMaterial(const TPZMaterial &mat);
    
    /// operator =
    TPZMaterial &operator=(const TPZMaterial &copy);
    
    /** @brief Default destructor */
    virtual ~TPZMaterial();
  
    /** @brief This method is responsible for setting
        the dimensions of the data structures used for
        computing the exact solution at an integration point.
        @param u_len length of the TPZVector<STATE> of the state variable
        @param du_row number of rows of the TPZFNMatrix<STATE> containing the appropriate state variable derivative
        @param du_col number of cols of the TPZFNMatrix<STATE> containing the appropriate state variable derivative*/
    virtual void GetExactSolDimensions(uint64_t &u_len,
                                            uint64_t &du_row,
                                            uint64_t &du_col);
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** 
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data);
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        // default is no specific data requirements
        if(type == 50)
        {
            data.fNeedsSol = true;
        }
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
    {
        // default is no specific data requirements
        int nref = datavec.size();
        if(type == 50)
        {
            for(int iref = 0; iref<nref; iref++){
                datavec[iref].fNeedsSol = true;
            }
        }
    }

    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data)
    {
        data.fNeedsNormal = false;
    }
    
    /** @brief This method defines which parameters need to be initialized in order to compute the contribution of interface elements */
    virtual void FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
    {
        data.SetAllRequirements(false);
        int nref_left = datavec_left.size();
        for(int iref = 0; iref<nref_left; iref++){
            datavec_left[iref].SetAllRequirements(false);
        }
        int nref_right = datavec_right.size();
        for(int iref = 0; iref<nref_right; iref++){
            datavec_right[iref].SetAllRequirements(false);
        }

    }
    
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "no_name"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const = 0;
    
    int Id() const { return fId; }
    void SetId(int id) {
/*        if(id == 0) {
            std::cout << "\n*** Material Id can't be ZERO! ***\n";
            std::cout << "*** This Will Be a Disaster!!! ***\n";
            DebugStop();
        }*/
        fId = id; }
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const = 0;
    
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 0;}
    
    /** @brief returns the number of load cases for this material object */
    int NumLoadCases()
    {
        return fNumLoadCases;
    }
    
    /** @brief returns the minimum number of load cases for this material */
    virtual int MinimumNumberofLoadCases()
    {
        return 1;
    }
    
    /** @brief changes the number of load cases for this material */
    void SetNumLoadCases(int numloadcases)
    {
        if(numloadcases <= 0)
        {
            std::cout << __PRETTY_FUNCTION__ << " numloadcases " << numloadcases << " cannot be less or equal to zero\n";
            DebugStop();
        }
        fNumLoadCases = numloadcases;
    }
    
    /** @brief indicates which variable should be post processed */
    void SetPostProcessIndex(int index)
    {
#ifdef PZDEBUG
        if (index < 0 || index >= fNumLoadCases)
        {
            DebugStop();
        }
#endif
        fPostProcIndex = index;
    }
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout);
    
    /** @name Post processing methods
     * @{
     */

    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * left, TPZCompEl * ritgh);	
    
protected:
    /** @deprecated Deprecated interface for Solution method which must use material data. */
    virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
    
public:
    
    /** @brief Computes the value of the flux function to be used by ZZ error estimator */
    virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol,
                      TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes,
                      TPZVec<STATE> &flux) {}
    
    /** @} */

    /** @brief Creates an object TPZBndCond derived of TPZMaterial*/
    virtual TPZBndCond *CreateBC(TPZMaterial *reference, int id, int typ, const TPZFMatrix<STATE> &val1,
                                 const TPZFMatrix<STATE> &val2);
    
    /** @name Contribute methods
	 * @{
	 */
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) = 0;
    
    /**
     * @brief It computes a contribution to the residual vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the residual vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef){
      TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
      this->Contribute(data, weight, fakeek, ef);
    }
    

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
    {
        TPZFMatrix<STATE> fakeek(ef.Rows(),ef.Rows(),0.);
        this->Contribute(datavec, weight, fakeek, ef);
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) = 0;
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
	 * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
    {
        TPZFMatrix<STATE> fakeek(ef.Rows(),ef.Rows(),0.);
        this->ContributeBC(datavec, weight, fakeek, ef, bc);
    }
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
      TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
      this->ContributeBC(data, weight, fakeek, ef, bc);
    }
	
    /** @} */
	
    /** 
	 * @brief Sets a procedure as source function for the material.
	 * @param fp pointer of the forces function
	 * @note Parameter loc corresponds to the coordinate of the point where the source function is applied
	 * @note Parameter result contains the forces resulting
	 */
    void SetForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
			fForcingFunction = fp;
    }
    void SetForcingFunction(void (*fp)(const TPZVec<REAL> &loc, TPZVec<STATE> &result), int porder )
		{
				if(fp)
                {
                    TPZDummyFunction<STATE> *loc = new TPZDummyFunction<STATE>(fp, porder);
                    loc->SetPolynomialOrder(porder);
                    fForcingFunction = loc;
                }
            
				else fForcingFunction = NULL;
		}
    
    void SetForcingFunction(void (*fp)(const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &gradu), int porder )
    {
        if(fp)
        {
            TPZDummyFunction<STATE> *loc = new TPZDummyFunction<STATE>(fp, porder);
            loc->SetPolynomialOrder(porder);
            fForcingFunction = loc;
        }
        
        else fForcingFunction = NULL;
    }

	/** @brief Returns a procedure as source function for the material */
	TPZAutoPointer<TPZFunction<STATE> > &ForcingFunction() {
		return fForcingFunction;
	}
	
    /** 
	 * @brief Sets a procedure as exact solution for the problem
	 * @param fp pointer of exact solution function
	 */
	void SetForcingFunctionExact(TPZAutoPointer<TPZFunction<STATE> > fp)
	{
		fForcingFunctionExact = fp;
	}
	
    /** @brief Returns a procedure as exact solution for the problem */
    TPZAutoPointer<TPZFunction<STATE> > &ForcingFunctionExact() {
        return fForcingFunctionExact;
    }
    
    /**
     * @brief Sets a procedure as time dependent source function for the material
     * @param fp pointer of the function
     */
    void SetTimeDependentForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
		fTimeDependentForcingFunction = fp;
    }
	
    /** @brief Returns a procedure as time dependent source function for the material */
    TPZAutoPointer<TPZFunction<STATE> > &TimeDependentForcingFunction() {
        return fTimeDependentForcingFunction;
    }
    
    /** 
	 * @brief Sets a procedure as time dependent exact solution for the problem
	 * @param fp pointer of the function
	 */
	void SetTimeDependentFunctionExact(TPZAutoPointer<TPZFunction<STATE> > fp)
	{
		fTimedependentFunctionExact = fp;
	}
    
    /** @brief Returns a procedure as time dependent exact solution for the problem */
    TPZAutoPointer<TPZFunction<STATE> > &TimedependentFunctionExact() {
        return fTimedependentFunctionExact;
    }
	
    /** 
     * @brief Sets a procedure as variable boundary condition
     * @param fp pointer of exact solution function
     */
    void SetBCForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fBCForcingFunction = fp;
    }
    
    /** @brief Returns a procedure as variable boundary condition */
    TPZAutoPointer<TPZFunction<STATE> > &BCForcingFunction() {
        return fBCForcingFunction;
    }
    
    /** 
     * @brief Sets a procedure as time variable boundary condition
     * @param fp pointer of exact solution function
     */
    void SetTimedependentBCForcingFunction(TPZAutoPointer<TPZFunction<STATE> > fp)
    {
        fTimedependentBCForcingFunction = fp;
    }
    
    /** @brief Returns a procedure as time variable boundary condition */
    TPZAutoPointer<TPZFunction<STATE> > &TimedependentBCForcingFunction() {
        return fTimedependentBCForcingFunction;
    }
    
    /** @brief Directive that gives true if the material has a forcing function   */
    virtual int HasForcingFunction() {return (fForcingFunction != 0);}
    
    /** @brief Directive that gives true if the material has a function exact  */
	virtual int HasForcingFunctionExact() {return (fForcingFunctionExact != 0);}
    
    /** @brief Directive that gives true if the material has a bc forcing function exact  */
    virtual int HasBCForcingFunction() {return (fBCForcingFunction != 0);}
    
    /** @brief Directive that gives true if the material has a time dependent function exact  */
    virtual int HasTimedependentFunctionExact() {return (fTimedependentFunctionExact != 0);}
    
    /** @brief Directive that gives true if the material has a time dependent forcing function   */
    virtual int HasTimedependentForcingFunction() {return (fTimeDependentForcingFunction != 0);}
    
    /** @brief Directive that gives true if the material has a time dependent bc forcing function   */
    virtual int HasTimedependentBCForcingFunction() {return (fTimedependentBCForcingFunction != 0);}
    
    
    /** @brief Gets the order of the integration rule necessary to integrate an element with polinomial order p */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;
	
	/** @brief Gets the order of the integration rule necessary to integrate an element multiphysic */
    virtual int IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const;
	
    void Errors(TPZMaterialData &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
    {
        TPZManVector<STATE,3> flux;
        Flux(data.x, data.sol[0], data.dsol[0], data.axes, flux);
        Errors(data.x, data.sol[0], data.dsol[0], data.axes, flux, u_exact, du_exact, errors );
    }
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors){
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
    
    virtual void ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc){
    
            PZError << __PRETTY_FUNCTION__ << std::endl;
            PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
        DebugStop();
        
    }

    /**
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
                        TPZVec<STATE> &uexact, TPZFMatrix<STATE> &duexact,
                        TPZVec<REAL> &val) {
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
	virtual	void ErrorsHdiv(TPZMaterialData &data, TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
		PZError << __PRETTY_FUNCTION__ << std::endl;
		PZError << "Nao sei o q fazer." << std::endl;
		
	}
    

    
    /** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
    virtual int NEvalErrors() {return 3;}
    
    /** @brief To create another material of the same type*/
    virtual TPZMaterial * NewMaterial();
    
    /** @brief Reads data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data);
    
    /** @brief Creates a copy of the material object and put it in the vector which is passed on */
    virtual void Clone(std::map<int, TPZMaterial * > &matvec);
    
    /** @brief To return a numerical flux type to apply over the interfaces of the elements */
    virtual int FluxType() { return 2; }
    
    /**
     * @brief Computes square of residual of the differential equation at one integration point.
     * @param X is the point coordinate (x,y,z)
     * @param sol is the solution vector
     * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
     */    
    virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
        return -1.;
    }
    
    /**
     * @brief Pushes a new entry in the context of materials with memory,
     * returning its index at the internal storage stack.
	 */
	/** To be implemented only in the proper materials. */
    virtual int PushMemItem(int sourceIndex = -1){ return -1; }
    
    /** @brief Frees an entry in the material with memory internal history storage */
    virtual void FreeMemItem(int index){ return; }
    
    /** @brief Sets fLinearContext attribute */
    void SetLinearContext(bool IsLinear);
    
    /** @brief Returns fLinearContext attribute */
    bool GetLinearContext() const {
        return fLinearContext;
    }
    
    /** @{
     * @name Save and Load methods
     */
    
    /** @brief Unique identifier for serialization purposes */
    public:
int ClassId() const override;

    
    /** @brief Saves the element data to a stream */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Reads the element data from a stream */
    void Read(TPZStream &buf, void *context) override;
    
    /** @} */
	
};

/** @brief Extern variable - Vector of force values */
extern TPZVec< void(*) (const TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

#endif

