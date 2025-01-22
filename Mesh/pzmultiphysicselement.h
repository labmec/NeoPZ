/**
 * @file
 * @brief Contains the declaration of the TPZMultiphysicsElement class. This class is abstract.
 */

#ifndef PZMULTIPHYSICSELEMENTH
#define PZMULTIPHYSICSELEMENTH 

#include <iostream>

#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "pzfunction.h"
#include "TPZMaterialDataT.h"

#include <map>

class TPZMultiphysicsInterfaceElement;

class TPZMultiphysicsElement : public TPZCompEl {
	
protected:
    
    /** @brief List of active approximation spaces */
    TPZManVector<MSpaceConfig,5> fActiveApproxSpace;

    template<class TVar>
    void TransferMultiphysicsElementSolutionT();
    template<class TVar>
    void ComputeRequiredDataT(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialDataT<TVar>> &datavec);
    template<class TVar>
    void ComputeRequiredDataT(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, std::map<int, TPZMaterialDataT<TVar>> &datavec);
public:
	/** @brief Default constructor */
	TPZMultiphysicsElement() : TPZCompEl()
	{
	}
	/**
	 * @brief Constructor
	 * @param mesh Multiphysics mesh where will be created the element
	 * @param ref geometric element reference
	 * @param index Index of the element created
	 */
	TPZMultiphysicsElement(TPZCompMesh &mesh, TPZGeoEl *ref) : TPZCompEl(mesh, ref)
	{
	}
  
  /** @brief Put a copy of the element in the referred mesh */
  TPZMultiphysicsElement(TPZCompMesh &mesh, const TPZMultiphysicsElement &copy);
  
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsElement()
	{
        TPZGeoEl *reference = Reference();
        if (reference) {
            reference->ResetReference();
        }
	}
	
	virtual void AddElement(TPZCompEl *cel, int64_t mesh) = 0;
    
    virtual void AddElement(const TPZCompElSide &cel, int64_t mesh) = 0;
    
    /// return the element associated with a mesh
    virtual TPZCompEl *Element(int64_t mesh) = 0;
	
    /** @brief Returns a reference to the element pointers vector */
    virtual TPZManVector<TPZCompElSide,5>   &ElementVec() = 0;
    
	virtual TPZCompEl *ReferredElement(int64_t mesh) = 0;
    
    virtual int64_t NMeshes() = 0;
	
	virtual void SetConnectIndexes(TPZVec<int64_t> &indexes) = 0;
    
    void GetConnectMeshPairs(TPZVec<std::pair<int64_t,int>> &connectpairs);
	
	virtual void AffineTransform(TPZVec<TPZTransform<> > &trVec) const = 0;
	
	virtual void InitMaterialData(TPZVec<TPZMaterialDataT<STATE>> &dataVec, TPZVec<int64_t> *indices = 0) = 0;
    virtual void InitMaterialData(TPZVec<TPZMaterialDataT<CSTATE>> &dataVec, TPZVec<int64_t> *indices = 0) = 0;
    
    virtual void InitMaterialData(std::map<int,TPZMaterialDataT<STATE>> &dataVec, TPZVec<int64_t> *indices = 0) = 0;
    virtual void InitMaterialData(std::map<int,TPZMaterialDataT<CSTATE>> &dataVec, TPZVec<int64_t> *indices = 0) = 0;
    
//    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec, TPZVec<int64_t> &indices);
    
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialDataT<STATE>> &datavec){
        if(datavec.size()) ComputeRequiredData(datavec[0],point);
        ComputeRequiredDataT(point,trvec,datavec);
    }
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialDataT<CSTATE>> &datavec){
        if(datavec.size()) ComputeRequiredData(datavec[0],point);
        ComputeRequiredDataT(point,trvec,datavec);
    }
    
    /** @brief Compute and fill data with requested attributes for each of the compels in fElementVec*/
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, std::map<int, TPZMaterialDataT<STATE>> &datavec){
        ComputeRequiredDataT(point,trvec,datavec);
    }
    virtual void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, std::map<int, TPZMaterialDataT<CSTATE>> &datavec){
        ComputeRequiredDataT(point,trvec,datavec);
    }
    //@{
    
    /** @brief TPZCompElWithMem can override this method. */
    virtual void ComputeRequiredData(TPZMaterialDataT<STATE> &data,
									 TPZVec<REAL> &qsi){

	}

	virtual void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data,
									 TPZVec<REAL> &qsi){

	}
    //@}
	/**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors (output) each norm or true error of the error of the solution at each physics
	 * @param flux (input) value of the interpolated flux values
	 */
    void EvaluateError(TPZVec<STATE> &errors, bool store_error) override;  

	virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override = 0 ;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override =0 ;
	
	//virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, std::set<int> dimension, std::set<int> dimension)=0;	
    
    void CreateInterfaces();
    
    bool ExistsInterface(int side);
    
    TPZMultiphysicsInterfaceElement * CreateInterface(int side);
    
    void RemoveInterfaces();
    
    void RemoveInterface(int side);
    
   	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const override = 0 ;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() = 0;
    
    /// After adding the elements initialize the integration rule
    virtual void InitializeIntegrationRule() override = 0 ;

    virtual void PolynomialOrder(TPZVec<int> &order) const = 0;
    
    virtual int ComputeIntegrationOrder() const override;
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        std::cout << "To Be Implemented\n";
        DebugStop();
    }

    virtual void TransferMultiphysicsElementSolution() override;
    
    /**
     * @brief Set the active approximation spaces
     * @param indexes List of the active approximation spaces
     */
    virtual void SetActiveApproxSpaces(TPZVec<MSpaceConfig> & active_approx_space)
    {
#ifdef PZDEBUG
        if(fActiveApproxSpace.size()!= ElementVec().size()){
            DebugStop();
        }
#endif
        fActiveApproxSpace = active_approx_space;
        
    }
    
    /**
     * @brief Set the active approximation spaces
     * @param indexes List of the active approximation spaces
     */
    virtual TPZVec<MSpaceConfig> & GetActiveApproxSpaces()
    {
        return fActiveApproxSpace;
        
    }
    
    virtual MSpaceConfig IsActiveApproxSpaces(int space_index)
    {
        return fActiveApproxSpace[space_index];
    }
    
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override  = 0 ;
    virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<CSTATE> &sol) override  = 0 ;
    

};

#endif
