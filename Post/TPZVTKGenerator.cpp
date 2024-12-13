/*********************************************************************/
/* File:   TPZVTKGenerator.cpp                                       */
/* Author: Francisco Orlandini                                       */
/* Date:   5. July 2022                                              */
/* Adapted from: NGSolve's vtkoutput.hpp                             */
/*********************************************************************/
#include "TPZVTKGenerator.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "TPZSimpleTimer.h"
#include "pzinterpolationspace.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMaterialDataT.h"
#include "pzmultiphysicselement.h"
#include "TPZGeoMeshTools.h"
#include "TPZRefPatternDataBase.h"


#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"



/**
@brief Auxiliary class for evaluating fields at computational elements.
Provides an abstraction of compel type (multiphysics or not). 
Template argument referes to solution type (real/complex).
This class allows for efficient computation of the fields by avoiding
repeated calls to InitMaterialData and ComputeRequiredData. Also,
it stores both TPZMaterialData<TVar> and a statically allocated 
vector of TPZMaterialDataT<TVar> in order to minimize dinamic allocation.

Whenever possible, the logic in this class prefers stack allocation.
*/
template<class TVar>
class TPZPostProcEl{
public:
  TPZPostProcEl(TPZCompEl *cel);
  //! Iniatilizes data for the computational element
  void InitData();
  //! Initializes data need at each integration point
  void ComputeRequiredData(TPZVec<REAL> &qsi);
  //! Evalutes the field id at the point qsi and stores result in sol
  void Solution(const TPZVec<REAL> &qsi, const int id, TPZVec<TVar> &sol);
private:
  bool fIsMultiphysics{false};
  TPZCompEl *fCel{nullptr};
  TPZMaterialDataT<TVar> fMatdata;
  TPZManVector<TPZMaterialDataT<TVar>,10> fDatavec; 
};

/**
   @brief Computes all relevant fields for a given computational element
   @param[in] cel Computational element
   @param[in] ref_vertices points (at reference el) to evaluate the fields
   @param[in] fields quantities to be evaluated
   @param[out]  initial position to store the results in the fields array
*/
template<class TVar>
void ComputeFieldAtEl(TPZCompEl *cel,
                      const TPZVec<TPZManVector<REAL,3>> &ref_vertices,
                      TPZVec<TPZAutoPointer<TPZVTKField>>& fields,
                      const TPZVec<int> &init_pos){
  // maximum size is 9 (3x3 tensor variable)
  TPZManVector<TVar,9> sol;

  const auto celdim = cel->Dimension();
  const TPZGeoEl *gel = cel->Reference();
  if(!gel) DebugStop();
  bool collapsed = false;
  int64_t node1 = -1;
  int64_t node2 = -1;
  int dir = -1;
  int leftright = 0;
  TPZManVector<REAL,2> scale(2,1.);
  if(gel->Type() == EQuadrilateral) {
    for(int in=0; in<4; in++) {
      node1 = gel->NodeIndex(in);
      node2 = gel->NodeIndex((in+1)%4);
      if(node1 == node2) {
        collapsed = true;
        if(in%2==0) {
          dir = 1;
          leftright = 2*(in/2)-1;
        } else {
          dir = 0;
          leftright = 2*(1-in/2)-1;
        }
        break;
      }
    }
  }

  TPZPostProcEl<TVar> graphel(cel);
  graphel.InitData();

  //vertex counter
  int iv = 0;
  //iterate through points in the reference element
  REAL shift_distance = 1.0e-2;
  for (auto ip : ref_vertices){
    ip.Resize(celdim);
#ifdef PZDEBUG
    if(ip.size() != celdim) DebugStop();
#endif
    if(collapsed && fabs(ip[dir]-leftright) < shift_distance) ip[dir] *= 1.-shift_distance;
    //computes all relevant data for a given integration point
    graphel.ComputeRequiredData(ip);
    const int nfields = fields.size();
    for (int i = 0; i < nfields; i++){
      auto &field = *(fields[i]);
      auto fdim = field.Dimension();
      //position will change depending on field dimension
      auto pos = init_pos[i] + fdim*iv;
      sol.Resize(fdim);
        sol.Fill(0.);
      graphel.Solution(ip, field.Id(), sol);

      /**TPZPostProcEl<TVar>::Solution might resize the sol array
         (i.e., 2d vectors). We thus make sure that we don't read
         out of bounds quantities and fill with zero if needed*/
      const auto sz = sol.size();
#ifdef PZDEBUG
      if(sz > fdim){
        //something wrong here
        PZError<<__PRETTY_FUNCTION__
               <<"\nsz: "<<sz<<" fdim: "<<fdim<<std::endl;
        DebugStop();
      }
#endif
      if constexpr (std::is_same_v<TVar,CSTATE>){
        for (int d = 0; d < sz; ++d){
          field[pos++] =  std::real(sol[d]);
        }
        for (int d = sz; d < fdim; ++d){
          field[pos++] =  0.0;
        }
      }else{
        for (int d = 0; d < sz; ++d){
          field[pos++] = sol[d];
        }
        for (int d = sz; d < fdim; ++d){
          field[pos++] = 0.0;
        }
      }
    }
    iv++;
  }
}


TPZVTKField::TPZVTKField(TPZVTKField::Type type, int id, std::string aname) :
  fType( type), fId(id), fName(aname) { ; }

TPZVTKGenerator::TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                                 const TPZVec<std::string> &fields,
                                 std::string filename,
                                 int vtkres, int dim)
  : TPZVTKGenerator(cmesh.operator->(),fields,filename,vtkres,dim)//delegates to other ctor
{
}



TPZVTKGenerator::TPZVTKGenerator(TPZCompMesh* cmesh,
                                 const TPZVec<std::string> &fields,
                                 std::string filename,
                                 int vtkres, int dim)
  : fCMesh(cmesh), fFilename(filename), fSubdivision(vtkres),
    fPostProcDim(cmesh->Dimension())
{
  if(dim > 0) {
    fPostProcDim = dim;
  }
  
    const auto nvars = fields.size();
  
  //let us check for valid post-processing matials
  for(auto [id,matp] : cmesh->MaterialVec()){
    auto bnd =
      dynamic_cast<TPZBndCond *>(matp);
      auto bnd2 = dynamic_cast<TPZMatCombinedSpacesBC<STATE> *>(matp);
    //we skip boundary materials
    if(matp && !bnd && !bnd2){
      if(matp->Dimension() != fPostProcDim) {continue;}
      bool foundAllVars{true};
      for(int i = 0; (i < nvars) && foundAllVars; i++){
        const auto &name = fields[i];        
        const auto index = matp->VariableIndex(name);
        foundAllVars = foundAllVars && (index > -1) ;
      }
      if(foundAllVars){
        fPostProcMats.insert(id);
      }
    }
  }

  if(fPostProcMats.size() == 0){
    std::cout<<"No post processing materials could be found!"<<std::endl;
    return;
  }
#ifdef PZDEBUG
  std::cout<<"The following materials will be post-processed:";
  for(auto id : fPostProcMats){std::cout<<" "<<id;}
  std::cout<<std::endl;
#endif
  InitFields(fields);
  
  //computes all points in the relevant reference element
  FillRefEls();
  /**computes all mapped points in the domain and resize all arrays
     to appropriate size*/
  ComputePointsAndCells();
}

TPZVTKGenerator::TPZVTKGenerator(TPZAutoPointer<TPZCompMesh> cmesh,
                                 std::set<int> mats,
                                 const TPZVec<std::string> &fields,
                                 std::string filename,
                                 int vtkres)
  : TPZVTKGenerator(cmesh.operator->(),mats, fields,filename,vtkres)//delegates to other ctor
{
}

TPZVTKGenerator::TPZVTKGenerator(TPZCompMesh* cmesh,
                                 std::set<int> mats,
                                 const TPZVec<std::string> &fields,
                                 std::string filename,
                                 int vtkres)
  : fCMesh(cmesh), fFilename(filename), fSubdivision(vtkres),
    fPostProcMats(std::move(mats))
{

  //checks post-processing dimension from given materials
  {
    const auto id = *(fPostProcMats.begin());
    auto *matp = cmesh->FindMaterial(id);
    if(!matp){
      PZError<<__PRETTY_FUNCTION__
             <<"\nError: material id "<<id<<" is not valid."
             <<std::endl;
      DebugStop();
    }
    fPostProcDim = matp->Dimension();
  }
  
  InitFields(fields);
  
  //computes all points in the relevant reference element
  FillRefEls();
  /**computes all mapped points in the domain and resize all arrays
     to appropriate size*/
  ComputePointsAndCells();
}


void TPZVTKGenerator::InitFields(const TPZVec<std::string> &fields)
{
  const int nvars = fields.size();
  fFields.resize(nvars);

  auto * matp = fCMesh->FindMaterial(*(fPostProcMats.begin()));
  for(int i = 0; i < nvars; i++){
    const auto &name = fields[i];
    const auto index = matp->VariableIndex(name);
    const auto dim = matp->NSolutionVariables(index);
    const TPZVTKField::Type type = [dim](){
      switch(dim){
      case 1:
        return TPZVTKField::Type::scal;
      case 2:
      case 3:
        return TPZVTKField::Type::vec;
      case 9:
        return TPZVTKField::Type::tens;
      default:
        PZError<<"TPZVTKGenerator:\n"
               <<"field dim "<<dim
               <<" not supported. It should be either "
               <<"\t1 (scalar)\n"
               <<"\t3 (vector)\n"
               <<"\t9 (tensor)\n";
        DebugStop();
        return TPZVTKField::Type::scal;
      }
    }();
    fFields[i] = new TPZVTKField(type,index,name);
  }
}

void TPZVTKGenerator::FillRefEls()
{
  const auto &dim = fPostProcDim;

  /**
     Uniform refinement patterns are used in order to compute
     the resulting points after the reference element has been divided
     fSubdivision times.

     The points are computed only for elements of correct dimension.
   */
  TPZSimpleTimer timer("FillRefEls");
  
  TPZManVector<TPZManVector<MElementType,3>,4> eltypes{
    {{}},
    {EOned},
    {ETriangle, EQuadrilateral},
    {ETetraedro, EPrisma, ECube, EPiramide}
  };

  for( auto type : eltypes[dim]){
    auto refp = gRefDBase.GetUniformRefPattern(type);
    if(!refp){gRefDBase.InitializeUniformRefPattern(type);}
    TPZVec<TPZManVector<REAL,3>> ref_coords;
    TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> ref_elems;
    switch(type){
    case EOned:
      FillReferenceEl<pztopology::TPZLine>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ETriangle:
      FillReferenceEl<pztopology::TPZTriangle>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EQuadrilateral:
      FillReferenceEl<pztopology::TPZQuadrilateral>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ETetraedro:
      FillReferenceEl<pztopology::TPZTetrahedron>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EPrisma:
      FillReferenceEl<pztopology::TPZPrism>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case ECube:
      FillReferenceEl<pztopology::TPZCube>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    case EPiramide:
      FillReferenceEl<pztopology::TPZPyramid>(ref_coords, ref_elems);
      fRefEls[type] = std::move(ref_elems);
      fRefVertices[type] = std::move(ref_coords);
      break;
    default:
      unreachable();
    }
  }
}

void TPZVTKGenerator::ComputePointsAndCells()
{

  /*
    Computes a list of mapped points in which the quantities will be evaluated.
    It also computes the resulting triangulation of the domain after the
    subdivisions were performed (vtk CELLS).
    After the points are computed, the fFields arrays are resized accordingly.
  **/
  fPoints.resize(0);

  /*
    vector of pairs of (compel,pos), where 
    - compel is a computational element apt for post-processing
    - pos is the position of its first post-processing node in the fPoints vector
  */
  fElementVec.resize(0);
  
  for (auto orig_cel : fCMesh->ElementVec()) {
    if(!orig_cel){continue;}
    /**
       orig_cel might be a TPZSubCompMesh, TPZElementGroup, etc.
       So we need to check exactly how many computational elements are "inside" it.
       note: the same logic does not apply for multiphysics elements, 
       where post-processing each element individually wouldn't make sense.
    */
    TPZStack<TPZCompEl*> cellist;
    orig_cel->GetCompElList(cellist);
    for(auto cel : cellist){
      if (! IsValidEl(cel)){continue;}
      //we need to add to valid elements for post-processing

      //we take note of the position of the first node
      const int offset = fPoints.size();
      AppendToVec(fElementVec,std::make_pair(cel,offset));
      
      const auto type = cel->Reference()->Type();
      TPZManVector<REAL, 3> pt(3, 0.);
      //append every post-processing node to fPoints vec
      for (auto &ip : fRefVertices[type]) {
        cel->Reference()->X(ip, pt);
        AppendToVec(fPoints, pt);
      }
      /*
        append every VTK cell corresponding to cel
        0-th position = element type
        1-th position = nnodes
        2--(nnodes+2) positions = node indexes
      */
      for (const auto &elem : fRefEls[type]) {
        std::array<int, TPZVTK::MAX_PTS + 2> new_elem = elem;
        
        const int npts = new_elem[1];
        for (int i = 0; i <= npts; ++i)
          new_elem[i+2] += offset;
        AppendToVec(fCells, new_elem);
      }
    }
    /**at this point we know how many post processing nodes are there, 
       so we can allocate memory just once.*/
    const int npts = fPoints.size();
    for(auto &f : fFields){
      //we first resize to zero to ensure that we allocate the exact size
      f->Resize(0);
      const int fdim = f->Dimension();
      f->Resize(npts*fdim);
    }
  }
}

/// Empty all fields, points and cells
void TPZVTKGenerator::ResetArrays()
{
  fElementVec.Resize(0);
  fPoints.Resize(0);
  fCells.Resize(0);
  for(auto field : fFields){
    field->Resize(0);
  }
}

template <class TOPOL>
void TPZVTKGenerator::FillReferenceEl(TPZVec<TPZManVector<REAL,3>> &ref_coords,
                                      TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &ref_elems){
  /*
    Computes every post-processing node at the REFERENCE element along with
    the VTK cells resulting from its subdivision
   */

  TPZAutoPointer<TPZGeoMesh> refpmesh = TPZGeoMeshTools::CreateGeoMeshSingleElT<TOPOL>(1,false,0);
  if(fSubdivision){
    const int nrefs = fSubdivision;

    TPZManVector<TPZGeoEl*,TPZVTK::MAX_SUBEL> sons(TPZVTK::MAX_SUBEL);
    for(int i = 0; i < nrefs; i++){
      for(auto gel : refpmesh->ElementVec()){
        if(gel->HasSubElement()==false){
          gel->Divide(sons);
        }
      }
    }
  }

  // fill nodes in the divided reference element
  const auto nnodes = refpmesh->NNodes();
  ref_coords.Resize(nnodes, {0, 0, 0});
  for (auto in = 0; in < nnodes; in++)
  {
    auto &node = refpmesh->NodeVec()[in];
    node.GetCoordinates(ref_coords[in]);
  }

  //fill cells in the divided reference element
  int ic = 0;//i-th cell
  for(TPZGeoEl* gel : refpmesh->ElementVec()){
    if(gel->HasSubElement()==false){
      const int gnnodes = gel->NNodes();
      AppendToVec(ref_elems,std::array<int,TPZVTK::MAX_PTS+2>{});
      //cell type
      ref_elems[ic][0] = TPZVTK::CellType(gel->Type());
      //nnodes
      ref_elems[ic][1] = gnnodes;
      //node indexes (reference element)
      for(int in = 0; in < gnnodes; in++){
        ref_elems[ic][in+2] = gel->NodeIndex(in);
      }
      ic++;
    }
  }
}

void TPZVTKGenerator::PrintPointsLegacy()
{
  /*print all post-processing nodes in legacy .VTK format*/
  TPZSimpleTimer timer("PrintPts");
  
  (*fFileout) << "POINTS " << fPoints.size() << " float" << std::endl;
  for(const auto &p : fPoints){
    for(const auto &x : p){
      *fFileout << x<<'\t';
    }
    *fFileout << std::endl;
  }
}


void TPZVTKGenerator::PrintCellsLegacy()
{
  /*print all resulting post-processing cells in legacy .VTK format*/
  TPZSimpleTimer timer("PrintCells");
  
  // count number of data for cells, one + number of vertices
  int ndata = 0;
  for (auto &c : fCells){
    ndata++;
    ndata += c[1];
  }
  *fFileout << "CELLS " << fCells.size() << " " << ndata << std::endl;
  for (const auto &c : fCells){
    const int nv = c[1];//n nodes
    *fFileout << nv << '\t';
    for (int i = 0; i < nv; i++)
      *fFileout << c[i + 2] << '\t';
    *fFileout << std::endl;
  }
}


void TPZVTKGenerator::PrintCellTypesLegacy()
{
  /*print all cell type info in legacy .VTK format*/
  TPZSimpleTimer timer("PrintCellTypes");
  
  *fFileout << "CELL_TYPES " << fCells.size() << std::endl;

  for (const auto &c : fCells){
    *fFileout << c[0] << '\n';
  }
  *fFileout << "CELL_DATA " << fCells.size() << std::endl;
  *fFileout << "POINT_DATA " << fPoints.size() << std::endl;
}


void TPZVTKGenerator::PrintFieldDataLegacy()
{
  /*print all post-processed quantities in legacy .VTK format*/
  TPZSimpleTimer timer("PrintField");
  
  for (auto field : fFields){
    *fFileout << field->TypeName() <<' ' << field->Name()
              << " float " << std::endl;
    if(field->GetType() == TPZVTKField::Type::scal){
      *fFileout << "LOOKUP_TABLE default" << std::endl;
      for (const auto &v : *field)
      {*fFileout << v << ' ';}
      *fFileout << '\n';
    }else{
      //both tensor and vector should be written as
      //(1,2,3)
      //(1,2,3)
      //etc
      const auto nfs = field->size();
      for (int i = 0; i < nfs; i++){
        *fFileout << (*field)[i] << ' ';
        if((i+1)%3 == 0){
          *fFileout << '\n';
        }
      }
    }
    
      if(field->TypeName()=="TENSORS"){
          *fFileout << "METADATA\n"
                    << "COMPONENT_NAMES\n"
                    << "XX\n"
                    << "XY\n"
                    << "XZ\n"
                    << "YX\n"
                    << "YY\n"
                    << "YZ\n"
                    << "ZX\n"
                    << "ZY\n"
                    << "ZZ\n\n";
      }
  }
}

bool TPZVTKGenerator::IsValidEl(TPZCompEl *cel)
{
  /*checks whether a given computational element should be post-processed.*/
  if(!cel || ! cel->Reference()){return false;}

  //check if element has an appropriate material
  if(fPostProcMats.count(cel->Material()->Id()) == 0) {return false;}
  //checks dimension and whether the element has been refined
  const auto gel = cel->Reference();
  const auto geldim = gel->Dimension();
  if(geldim != fPostProcDim){return false;}

  
  /* the condition below is not always applicable.
     for instance, if one uses the same geometric mesh, 
     with several refinements, but with one comp mesh for each
     refinement level.
     anyway, i will leave this here in case it provides a hint for debugging*/
  
  // if(gel->HasSubElement()){return false;}

  switch (gel->Type()){
    case EPoint:
    case EOned:
    case ETriangle:
    case EQuadrilateral:
    case ETetraedro:
    case EPiramide:
    case EPrisma:
    case ECube:
      return true;
    //it should really not reach any other type, but let us be on the safe side
    case EPolygonal: /*8*/	
    case EInterface: /*9*/	
    case EInterfacePoint: /*10*/	
    case EInterfaceLinear: /*11*/	
    case EInterfaceSurface: /*12*/	
    case ESubstructure: /*13*/	
    case EGlobLoc: /*14*/	
    case EDiscontinuous: /*15*/	
    case EAgglomerate: /*16*/	
    case ENoType: /*17*/	
      return false;
    }
  return false;
}

void TPZVTKGenerator::Do()
{
  if(fPostProcMats.size() == 0){
    std::cout<<"No post-processing materials found..."<<std::endl;
    return;
  }
  /*
    Main function for exporting one .VTK file.
    May be called repeatedly if the solution associated with 
    the computational mesh has changed. 
    NOTE: if the mesh has changed (i.e., geometric refinement), 
    ResetArrays() must be called.
*/
  TPZSimpleTimer timer("Do");

  std::ostringstream filenamefinal;
  std::stringstream appended;
  std::vector<int> datalength;
  int offs = 0;
  // create name of the current .vtk file
  filenamefinal << fFilename << '.';
  if(fStep > -1){
    filenamefinal << fStep;
  }else{
    filenamefinal << fOutputCount;
  }

  fLastOutputName = filenamefinal.str();

  filenamefinal << ".vtk";

  fFileout = new std::ofstream(filenamefinal.str());

  if(fPoints.size() == 0){
    /* A call to ResetArrays has been made. 
       Therefore all post-processing nodes must be recomputed and fFields
       arrays must be resized accordingly*/
    ComputePointsAndCells();
  }

  // header:
  *fFileout << "# vtk DataFile Version 3.0" << std::endl;
  *fFileout << "vtk output" << std::endl;
  *fFileout << "ASCII" << std::endl;
  *fFileout << "DATASET UNSTRUCTURED_GRID" << std::endl;

  
  const int nel = fElementVec.size();

  if(fNThreads == 0){
    
    ComputeFields(0,nel);
  }else{
    std::vector<std::thread> allthreads;
    const int nt = NThreads();
    //number of elements per thread (last one may have more elements)
    const int threadnel = nel / nt;
    for(int i = 0; i < nt; i++){
      int firstel = threadnel*i;
      int lastel = i == nt -1 ? nel : firstel + threadnel;
      allthreads.push_back(
        std::thread(
          [this] (int f, int l){this->ComputeFields(f,l);},firstel,lastel));
    }

    for(int i = 0; i < nt; i++){
      allthreads[i].join();
    }
  }

  
  
  
  //write everything to file
  PrintPointsLegacy();
  PrintCellsLegacy();
  PrintCellTypesLegacy();
  PrintFieldDataLegacy();

  fOutputCount++;
  if(fStep > -1){fStep++;}
  std::cout << " Done." << std::endl;

}

void TPZVTKGenerator::ComputeFields(int first, int last)
{

  //whether the mesh has real or cplx dofs
  const bool isCplxMesh = fCMesh->GetSolType() == ESolType::EComplex;

  //number of fields to be computed
  const int nfields = fFields.size();
  TPZManVector<int> posvec(nfields);

  //iterate through valid range
  for (int i = first; i < last; i++){
    auto [cel,pos] =  fElementVec[i];
    const auto eltype = cel->Reference()->Type();

    /*
      pos is the index of the first post-processing node of the element.
      To find the position in the fFields[i] vec, one must take into account
      the dimension of the each field.
     */
    for(auto f = 0; f < nfields; f++){
      posvec[f] = pos * fFields[f]->Dimension();
    }
    if (isCplxMesh) {
      ComputeFieldAtEl<CSTATE>(cel, fRefVertices[eltype], fFields,posvec);
    } else {
      ComputeFieldAtEl<STATE>(cel, fRefVertices[eltype], fFields,posvec);
    }
  }
}


template<class TVar>
TPZPostProcEl<TVar>::TPZPostProcEl(TPZCompEl *cel) : fCel(cel){
  auto test = dynamic_cast<TPZMultiphysicsElement*>(cel);
  if(test){fIsMultiphysics = true;}
}

template<class TVar>
void TPZPostProcEl<TVar>::InitData(){
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    const int64_t nref = mfcel->NMeshes();
    fDatavec.resize(nref);
    mfcel->InitMaterialData(fDatavec);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    intel->InitMaterialData(fMatdata);
  }
}

template<class TVar>
void TPZPostProcEl<TVar>::ComputeRequiredData(TPZVec<REAL> &qsi){
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    const int64_t nref = mfcel->NMeshes();
    for(int ir = 0; ir < nref; ir++){
      fDatavec[ir].fNeedsSol= true;
    }
    TPZManVector<TPZTransform<> > trvec;
    mfcel->AffineTransform(trvec);
    mfcel->ComputeRequiredData(qsi, trvec, fDatavec);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    fMatdata.fNeedsSol = true;
		intel->ComputeRequiredData(fMatdata, qsi);
  }
}


template<class TVar>
void TPZPostProcEl<TVar>::Solution(const TPZVec<REAL> &qsi, const int id, TPZVec<TVar> &sol)
{
    if(id >= 100) {
        TPZVec<REAL> &qsicopy = (TPZVec<REAL> &) qsi;
        fCel->Solution(qsicopy, id, sol);
        return;
    }
  TPZMaterial * material = fCel->Material();
  if(fIsMultiphysics){
    auto mfcel = dynamic_cast<TPZMultiphysicsElement*>(fCel);
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
    matCombined->Solution(fDatavec,id,sol);
  }else{
    auto intel = dynamic_cast<TPZInterpolationSpace*>(fCel);
    auto *matSingle =
       dynamic_cast<TPZMatSingleSpaceT<TVar>*>(material);
    if(matSingle) {
        matSingle->Solution(fMatdata,id,sol);
    } else {
      auto *matCombined = dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
      if(!matCombined) {
        DebugStop();
      }
      fDatavec.Resize(2);
      fDatavec[0] = fMatdata;
      matCombined->Solution(fDatavec,id,sol);
    }
  }
}


#define FILLREF(T)\
  template \
  void TPZVTKGenerator::FillReferenceEl<T>(TPZVec<TPZManVector<REAL,3>> &, \
                                           TPZVec<std::array<int,TPZVTK::MAX_PTS + 2>> &);

FILLREF(pztopology::TPZTriangle)
FILLREF(pztopology::TPZQuadrilateral)
FILLREF(pztopology::TPZTetrahedron)
FILLREF(pztopology::TPZCube)
FILLREF(pztopology::TPZPrism)
#undef FILLREF
