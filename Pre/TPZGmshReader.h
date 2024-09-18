//
//  TPZGmshReader.h
//  PZ
//
//  Created by Omar on 2/7/16.
//
//

#ifndef TPZGmshReader_h
#define TPZGmshReader_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "pzgmesh.h"
#include "SPZPeriodicData.h"

class TPZGeoMesh;


struct MaterialDataS {
    
    TPZStack<int> fMatID;
    TPZStack<std::pair<int ,std::string> >  fMaterial;
    
    MaterialDataS() : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(int num) : fMatID(), fMaterial(){
        
    }
    
    MaterialDataS(const MaterialDataS &copy) : fMatID(copy.fMatID),
    fMaterial(copy.fMaterial) {
    }
    
    MaterialDataS &operator=(const MaterialDataS &copy){
        fMatID = copy.fMatID;
        fMaterial = copy.fMaterial;
        return *this;
    }
    
};


/**
 * @brief Implement the interface between TPZGeoMesh and the files produced by Gmsh (version 3.0 or 4.0 ) in msh format.
 * @since January 16, 2017
 */

/** What is Gmsh ? Take a look on http://gmsh.info/
 * Gmsh is a free 3D finite element grid generator with a build-in CAD engine and post-processor. Its design goal is to provide a fast, light and user-friendly
 * meshing tool with parametric input and advanced visualization capabilities. Gmsh is built around four modules: geometry, mesh, solver and post-processing.
 * The specification of any input to these modules is done either interactively using the graphical user interface or in ASCII text files using Gmsh's own
 * scripting language.
 */

/** Note about the implementation for file format 4
 * The mandatory sections are considered MeshFormat, Entities, Nodes and Elements.
 * The optional section PhysicalName is considered the others (PartitionedEntities,Periodic,GhostElements,NodeData,ElementData,ElementNodeData) are just ignored.
 * To conclude a successful read of your *.msh file, you should have physical tags to be able to insert elements into a TPZGeoMesh object.
 */
class TPZGmshReader{
    
    /// gmsh file format version (supported versions = {3,4})
    std::string m_format_version;
    
    /// pointer to the geometric mesh being read
    TPZGeoMesh *m_gmesh = 0;
    
    /// Number of volumes
    int m_n_volumes = 0;
    
    /// Number of surfaces
    int m_n_surfaces = 0;
    
    /// Number of curves
    int m_n_curves = 0;
    
    /// Number of points
    int m_n_points = 0;
    
    /// Number of volumes with physical tag
    int m_n_physical_volumes = 0;
    
    /// Number of surfaces with physical tag
    int m_n_physical_surfaces = 0;
    
    /// Number of curves with physical tag
    int m_n_physical_curves = 0;
    
    /// Number of points with physical tag
    int m_n_physical_points = 0;
    
    /// Geometry dimension
    int m_dimension = -1;
    
    /// Characteristic length to apply a Scale affine transformation
    REAL m_characteristic_lentgh = 1;

    //////////// Members related to file format with version 4 ////////////
    /// Data structure of both: physical entities and names indexed by dimension
    TPZManVector<std::map<int,std::vector<int>>,4> m_dim_entity_tag_and_physical_tag;
    //////////// Members related to file format with version 3 ////////////
    
    /// Structure of both: physical entities and names indexed by dimension
    TPZManVector<std::map<int,std::string>,4> m_dim_physical_tag_and_name;
    
    /// Structure of both: names and physical id indexed by dimension
    TPZManVector<std::map<std::string,int>,4> m_dim_name_and_physical_tag;
    
    /// Structure of both: physical id and user defined physical tag indexed by dimension
    TPZManVector<std::map<int,int>,4> m_dim_physical_tag_and_physical_tag;
    
    /// Flag indicating whether elements without defined physical tag should be read
    bool m_read_undefined_physical_tag_elements = true;
    
    /// Entity index to which the element belongs
    TPZVec<int64_t> m_entity_index;

    /// Structure for periodic entities
    TPZAutoPointer<SPZPeriodicData> m_periodic_data;
    
    /// Number of hexahedra
    int m_n_hexahedron_els = 0;
    
    /// Number of tetrahedra
    int m_n_tetrahedron_els = 0;
    
    /// Number of prisms
    int m_n_prism_els = 0;
    
    /// Number of pyramids
    int m_n_pyramid_els = 0;
    
    /// Number of quadrilaterals
    int m_n_quadrilateral_els = 0;
    
    /// Number of triangles
    int m_n_triangle_els = 0;
    
    /// Number of triangles
    int m_n_line_els = 0;
    
    /// Number of points
    int m_n_point_els = 0;
    
    /// verbosity of output
    int m_verbose = 0;
    
    /// Convert a Gmsh *.msh file with format 4 to a TPZGeoMesh object
    TPZGeoMesh * GeometricGmshMesh4(const std::string &file_name, TPZGeoMesh *gmesh = NULL, bool addNonAssignedEls = true);

    /// Convert a Gmsh *.msh file with format 3 to a TPZGeoMesh object
    TPZGeoMesh * GeometricGmshMesh3(const std::string &file_name, TPZGeoMesh *gmesh = NULL);

    // read the version line at the top of the gmsh file
    void ReadVersion(const std::string &file_name);
    
    // read the version line at the top of the gmsh file
    void ReadVersion(std::istream &file_name);
    
public:
    /// read the physical property section of the gmsh file
    void ReadPhysicalProperties4(std::istream &input);

private:
    /// read the elements and nodes
    void ReadElements4(std::istream &input);
    
    /// read the data with respect to periodic nodes
    void ReadPeriodic4(std::istream &input);
    
    //! Fills m_periodic_els structure after creating the mesh
    void FillPeriodicData(
      TPZGeoMesh *gmesh,
      const TPZVec<std::map<int64_t, std::map<int64_t, int64_t>>>
      &entity_periodic_nodes,
      const TPZVec<std::map<int64_t, int64_t>> &periodic_entities);

  public:
    /// Default constructor
    TPZGmshReader();
    
    /// Default destructor
    ~TPZGmshReader() = default;
    
    /// setting copy constructor to default
    TPZGmshReader(const TPZGmshReader &copy) = default;
    
    /// setting = operator to default
    TPZGmshReader &operator=(const TPZGmshReader &copy) = default;
    
    /// Convert Gmsh msh files in a TPZGeoMesh object, detecting .msh version
    TPZGeoMesh * GeometricGmshMesh(const std::string &file_name, TPZGeoMesh *gmesh = NULL, bool addNonAssignedEls = true);

    
    /// Set the Characteristic length (before reading the mesh)
    void SetCharacteristiclength(REAL length);
    
    /// Print the partition summary after the reading process
    void PrintPartitionSummary(std::ostream & out);
    
    /// set the verbose output level (0 -> minimum output)
    void SetVerbose(int verbose) {
        m_verbose = verbose;
    }

    void InsertElement(TPZGeoMesh *gmesh, const int physical_identifier,
                       const int el_type, const int el_identifier,
                       const std::vector<int> &node_identifiers);

    int GetNumberofNodes(int & el_type);

    //! Get periodic data
    auto GetPeriodicData() const { return m_periodic_data; }
    /// Insert elements following msh file format */
    bool InsertElement(TPZGeoMesh * gmesh, std::ifstream & line);
    
    /// Get the structure dim -  physical tag - name
    TPZManVector<std::map<int,std::string>,4> & GetDimPhysicalTagName(){
        return m_dim_physical_tag_and_name;
    }
    
    /// Set the structure dim -  physical tag - name
    void SetDimPhysicalTagName(TPZManVector<std::map<int,std::string>,4> & dim_physical_tag_and_name){
        m_dim_physical_tag_and_name = dim_physical_tag_and_name;
    }
    
    /// Get the structure dim - name - physical tag
    void SetDimNamePhysical(TPZManVector<std::map<std::string,int>,4> & dim_name_and_physical_tag){
        m_dim_name_and_physical_tag = dim_name_and_physical_tag;
    }
    
    /// Get the structure dim - name - physical tag
    
    TPZManVector<std::map<std::string,int>,4> & GetDimNamePhysical(){
        return m_dim_name_and_physical_tag;
    }
    
    /// Get the structure dim - name - physical tag
    void GetDimNamePhysical(TPZVec<std::map<std::string,int>> & dim_name_and_physical_tag){
        m_dim_name_and_physical_tag = dim_name_and_physical_tag;
    }
    
    TPZVec<int64_t> &EntityIndex()
    {
        return m_entity_index;
    }
    
    /// Return the number of hexahedra created
    int NHexahedra(){
        return m_n_hexahedron_els;
    }
    
    /// Return the number of tetrahedra created
    int NTetrahera(){
        return m_n_tetrahedron_els;
    }
    
    /// Return the number of prisms created
    int NPrisms(){
        return m_n_prism_els;
    }
    
    /// Return the number of pyramids created
    int NPyramids(){
        return m_n_pyramid_els;
    }
    
    /// Return the number of quadrilaterals created
    int NQuadrilaterals(){
        return m_n_quadrilateral_els;
    }
    
    /// Return the number of triangles created
    int NTriangles(){
        return m_n_triangle_els;
    }
    
    /// Return the number of line created
    int NLines(){
        return m_n_line_els;
    }
    
    /// Return the number of line created
    int NPoints(){
        return m_n_point_els;
    }

    /// Return the dimension of the mesh
    int Dimension() {
        return m_dimension;
    }
};

#endif /* TPZGmshReader_h */
