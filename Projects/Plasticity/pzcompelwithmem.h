//$Id: pzcompelwithmem.h,v 1.12 2009-04-20 03:33:44 erick Exp $

#ifndef PZCOMPELWITHMEM_H
#define PZCOMPELWITHMEM_H

class TPZMaterialData;

#include "pzinterpolationspace.h"
#include "pzstack.h"
#include "pzcmesh.h"
#include "pzquad.h" //TPZIntPoints
#include "pzmaterial.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr CompElWMemlogger(Logger::getLogger("mesh.TPZCompElWithMem"));
#endif

/**
This class implements the TPZCompEl structure to enable material memory feature.
It should be instantiated using one of TPZCompEl bottom classes as parent in the template parameter.
 The computational element identifies for each of its points wich is the context of plastic history.
 This class is responsable to ask for storage space or instruct the instance of the class TPZMatElastoplastic to free a set of positions to storage or to destruct or create it. 
@since July, 23 2008
*/



template <class TBASE>

class TPZCompElWithMem : public TBASE
{
public:

  TPZCompElWithMem();

  virtual ~TPZCompElWithMem();
		
  TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy);

  /**
   * used to generate patch mesh... generates a map of connect index from
   * global mesh to clone mesh
   */
  TPZCompElWithMem(TPZCompMesh &mesh,
              const TPZCompElWithMem<TBASE> &copy,
              std::map<int,int> & gl2lcConMap,
              std::map<int,int> & gl2lcElMap);


  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZCompElWithMem<TBASE> (mesh, *this);
  }

  /**
   * Create a copy of the given element. The clone copy have the connect indexes
   * mapped to the local clone connects by the given map
   * @param mesh Patch clone mesh
   * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
   * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const
  {
    return new TPZCompElWithMem<TBASE> (mesh, *this, gl2lcConMap, gl2lcElMap);
  }
	
  void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &qsi);

protected:
		
  /**
   * PrepareIntPtIndices initializes the material damage varibles memory in the proper material class.
   */
  void PrepareIntPtIndices();
	
  /**
   * Frees the material damage varibles memory in the proper material class.
   */
  void SetFreeIntPtIndices();
	
  void CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy);
	
public:

  /**
   * Save the element data to a stream
   */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  virtual void Read(TPZStream &buf, void *context);

	
private:
		
  TPZStack<int,128> fIntPtIndices;
	
};

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem() : TBASE() {
	PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::~TPZCompElWithMem() {
	SetFreeIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
  TBASE(mesh, gel, index){
	PrepareIntPtIndices();
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh, const TPZCompElWithMem<TBASE> &copy) :
  TBASE(mesh, copy) {
	CopyIntPtIndicesFrom(copy);
}

template<class TBASE>
TPZCompElWithMem<TBASE>::TPZCompElWithMem(TPZCompMesh &mesh,
                                      const TPZCompElWithMem<TBASE> &copy,
                                      std::map<int,int> & gl2lcConMap,
                                      std::map<int,int> & gl2lcElMap) :
    TBASE(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	CopyIntPtIndicesFrom(copy);
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::PrepareIntPtIndices(){

  // This code was copied from TPZInterpolationSpace::CalcStiff with minor changes
  // regarding integration point index evaluation.
  // Inclusions were commented properly.
	
  if(fIntPtIndices.NElements())	
	{
#ifdef LOG4CXX
{
  std::stringstream sout;
  sout << __PRETTY_FUNCTION__ << " Attempting to add memory indices to an already configured TPZCompElWithMem";
  LOGPZ_ERROR(CompElWMemlogger,sout.str().c_str());
}
#endif
	}
	
  TPZMaterial * material = TBASE::Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    return;
  }

  int p = TBASE::MaxOrder();

  int dim = TBASE::Dimension();

	/**
	 * Só para teste! Não sabemos que tipo é (TPZIntCube3D ou outro tipo). PHIL HELP ME!!!
	 */									    //*/
	//TPZIntCube3D intrule = *(dynamic_cast<TPZIntCube3D *>(TBASE::GetIntegrationRule().Clone()));
	TPZIntPoints &intrule =*(TBASE::GetIntegrationRule().Clone());
	TPZManVector<int,3> p2(dim,p*2);
	intrule.SetOrder(p2);
  if(material->HasForcingFunction()) {
    TPZManVector<int,3> order(dim,intrule.GetMaxOrder());
    intrule.SetOrder(order);
  }

  int intrulepoints = intrule.NPoints();
  fIntPtIndices.Resize(intrulepoints);
	
  for(int int_ind = 0; int_ind < intrulepoints; ++int_ind){
	fIntPtIndices[int_ind] = this->Material()->PushMemItem();
	 // Pushing a new entry in the material memory
  } //Loop over integratin points generating a reference vector of memory
	//entries in the related pzmatwithmem for further use.

}


template <class TBASE>
inline void TPZCompElWithMem<TBASE>::SetFreeIntPtIndices(){

  TPZMaterial * material = TBASE::Material();
  if(!material){
    PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    return;
  }

  int n = fIntPtIndices.NElements();
	
  for(int i = 0; i < n; i++){
	this->Material()->FreeMemItem(fIntPtIndices[i]);
  }
	
  fIntPtIndices.Resize(0);

}

template <class TBASE>
void TPZCompElWithMem<TBASE>::CopyIntPtIndicesFrom(const TPZCompElWithMem<TBASE> & copy)
{
	
	TPZMaterial * material = TBASE::Material();
  	if(!material){
    	PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
    	return;
  	}
	
	int i, n = copy.fIntPtIndices.NElements();
	fIntPtIndices.Resize(n);
	
	for(i = 0; i < n; i++)fIntPtIndices[i] = material->PushMemItem(copy.fIntPtIndices[i]);
}

/**
 * Save the element data to a stream
 */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::ComputeRequiredData(TPZMaterialData &data,
                                                TPZVec<REAL> &qsi){
  TBASE::ComputeRequiredData(data, qsi); 
													
  data.intPtIndex = fIntPtIndices[ data.intPtIndex ]; // returning the 
  //material index for the n-th CompEl integration point									
 
}

template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Write(TPZStream &buf, int withclassid)
{
  TBASE::Write(buf,withclassid);
  TPZSaveable::WriteObjects(buf, fIntPtIndices);
}

/**
 * Read the element data from a stream
 */
template <class TBASE>
inline void TPZCompElWithMem<TBASE>::Read(TPZStream &buf, void *context)
{
  TBASE::Read(buf,context);
  TPZSaveable::ReadObjects(buf, fIntPtIndices);
}

#endif
