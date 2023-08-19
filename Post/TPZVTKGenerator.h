/*********************************************************************/
/* File:   TPZVTKGenerator.h                                         */
/* Author: Francisco Orlandini                                       */
/* Date:   5. July 2022                                              */
/* Adapted from: NGSolve's vtkoutput.hpp                             */
/*********************************************************************/

#ifndef _TPZVTKGENERATOR_H_
#define _TPZVTKGENERATOR_H_

#include "pzmanvector.h"
#include "tpzautopointer.h"
#include "pzeltype.h"
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <array>

class TPZCompMesh;



namespace TPZVTK{
  //! max number of nodes for a given cell (hexahedron)
  static constexpr int MAX_PTS{8};
  //! max number of sub elements after a geometric refinement (pyramid)
  static constexpr int MAX_SUBEL{10};

  //! cell types for .VTK format
  static int CellType(const MElementType el){
    switch (el){
    case EPoint:
      return 1;
    case EOned:
      return 3;
    case ETriangle:
      return 5;
    case EQuadrilateral:
      return 9;
    case ETetraedro:
      return 10;
    case EPiramide:
      return 14;
    case EPrisma:
      return 13;
    case ECube:
      return 12;
    default:
      std::cout << "TPZVTKGenerator Element Type "
                << MElementType_Name(el) << " not supported!"
                << std::endl;
      DebugStop();
      unreachable();
    }
  }
};
/**
   @brief Stores pointwise information on post-processed variables.
   This class is used internally by the TPZVTKGenerator
*/
class TPZVTKField : public TPZVec<STATE>
{
public:
  enum class Type{scal,vec,tens};
  TPZVTKField() { ; };
  TPZVTKField(Type t, int id, std::string name);
  void SetType(Type t) {fType = t;}
  Type GetType() const {return fType;}
  std::string TypeName() const{
    switch (fType){
    case Type::scal: return "SCALARS";
    case Type::vec: return "VECTORS";
    case Type::tens: return "TENSORS";
    }
    unreachable();
  }
  
  int Dimension() const {
    switch (fType){
    case Type::scal: return 1;
    case Type::vec: return 3;
    case Type::tens: return 9;
    }
    unreachable();
  }
  void SetName(std::string name) { fName = name; }
  std::string Name() const{ return fName; }
  void SetId(int id) {fId = id;}
  int Id() {return fId;}
private:
  Type fType = Type::scal;
  int fId = -1;
  std::string fName = "none";
};


/**
   @brief Manages exporting of post-processed quantities to .VTK files.
*/
class TPZVTKGenerator{
protected:
  //! Computational mesh
  TPZCompMesh* fCMesh = nullptr;
  //! Set of materials in which post-processing will take place
  std::set<int> fPostProcMats;
  //! Post-processing dimension
  int fPostProcDim{-1};
  //! File name (no extension)
  std::string fFilename = "";
  //! Number of subdivisions of each geometric element
  int fSubdivision{0};
  //! Post-processed quantities
  TPZManVector<TPZAutoPointer<TPZVTKField>> fFields;
  //! Geometric coordinates of points in which quantities will be processed
  TPZVec<TPZManVector<REAL,3>> fPoints;
  /** @brief Domain triangulation after subdivisions.
      Each cell is organised as type nnodes node0 node1 node2 ... nodeN*/
  TPZVec<std::array<int,TPZVTK::MAX_PTS+2>> fCells;//max 
  //! Used for exporting .VTK series (time-steps, different modes, etc)
  int fOutputCount = 0;
  //! Used for manually setting the step value (fOutputCount is ignored if fStep>0)
  int fStep = -1;
  //! Current file in which .VTK data will be written
  TPZAutoPointer<std::ofstream> fFileout{nullptr};
  //! Name (no extension) of last output file
  std::string fLastOutputName = "";
  //! All nodes in the reference element in which quantities are evaluated
  std::map<MElementType,TPZVec<TPZManVector<REAL,3>>> fRefVertices;
  //! All sub-elements resulting of dividing the original reference element
  std::map<MElementType,TPZVec<std::array<int,TPZVTK::MAX_PTS+2>>> fRefEls;
  //! (el,first pos) list of elements and the initial position of their solution in the field vec
  TPZVec<std::pair<TPZCompEl*,int>> fElementVec;
  //! number of threads to be used for computing the fields
  int fNThreads{4};

  //! Computes points (at reference element) for evaluating fields
  template<class TOPOL>
  void FillReferenceEl(TPZVec<TPZManVector<REAL,3>> &ref_coords,
                       TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems);
  //! Calls FillReferenceEl to fill fRefVertices and fRefEls for all relevant topologies
  void FillRefEls();
  //! Print all points in VTK Legacy format
  void PrintPointsLegacy();
  //! Print all cells in VTK Legacy format
  void PrintCellsLegacy();
  //! Print all cell types in VTK Legacy format
  void PrintCellTypesLegacy();
  //! Print all fields in VTK Legacy format
  void PrintFieldDataLegacy();
  //! Check if element should be processed
  bool IsValidEl(TPZCompEl *el);
  //! Initializes field data
  void InitFields(const TPZVec<std::string> &fields);
  /**@brief Compute all post-processing points and vtk cells. 
     It also creates list of valid computational elements for post-processing.*/
  void ComputePointsAndCells();
  //!Compute fields for a range of elements
  void ComputeFields(int first, int last);
public:
  /**
     @brief Creates instance for generating .vtk results for a given mesh
     @param[in] cmesh Computational mesh
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename (without extension)
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
     @param[in] dim Post-processing dimension (defaults to dimension of mesh)
  */
  TPZVTKGenerator(TPZCompMesh* cmesh,
                  const TPZVec<std::string> &fields,
                  std::string filename,
                  int vtkres,
                  int dim = -1);
  /**
     @brief Creates instance for generating .vtk results for a given mesh
     @param[in] cmesh Computational mesh
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename without extension
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
     @param[in] dim Post-processing dimension (defaults to dimension of mesh)
  */
  TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                  const TPZVec<std::string> &fields,
                  std::string filename,
                  int vtkres,
                  int dim = -1);

  /**
     @brief Creates instance for generating .vtk results for given materials in a given mesh
     @param[in] cmesh Computational mesh
     @param[in] mats identifiers of materials to be post-processed
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename (without extension)
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
     @note All materials should have the same dimension
  */
  TPZVTKGenerator(TPZCompMesh* cmesh,
                  std::set<int> mats,
                  const TPZVec<std::string> &fields,
                  std::string filename,
                  int vtkres);
  /**
     @brief Creates instance for generating .vtk results for given materials in a given mesh
     @param[in] cmesh Computational mesh
     @param[in] mats identifiers of materials to be post-processed
     @param[in] fields names of fields to be post-processed
     @param[in] filename filename without extension
     @param[in] vtkres resolution of vtk post-processing (number of el subdivision)
     @note All materials should have the same dimension
  */
  TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                  std::set<int> mats,
                  const TPZVec<std::string> &fields,
                  std::string filename,
                  int vtkres);
  //!Generates .vtk file for a given current solution of the mesh
  void Do(REAL time = -1);

  /** @brief Resets post-processing data structure. Call this function if
   the mesh has had changes between Do() calls (refinement, etc)*/
  void ResetArrays();

  //! Set number of threads to compute fields (0 for serial)
  void SetNThreads(int nt) {
    //let us avoid negative number of threads for the sake of it. 0 for serial
    if(nt > -1 ){fNThreads = nt;}
  }
  //! Get number of threads to compute fields (0 for serial)
  int NThreads() const {return fNThreads;}
  //! Sets manually the step of next post-processing round (must be >0)
  void SetStep(int st){fStep = st;}
  //! Gets the current step of next post-rpocessing round (if <0, fOutputCount is used)
  int Step() const {return fStep;}
    
private:
  // No copy constructor allowed until someone implements it properly
  TPZVTKGenerator(const TPZVTKGenerator& vtkgen);
    
  // No operator= allowed until someone implements it properly
  TPZVTKGenerator& operator=(const TPZVTKGenerator& vtkgen);
};

#endif /* _TPZVTKGENERATOR_H_ */
