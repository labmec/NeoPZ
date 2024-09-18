#ifndef _SPZPERIODICDATA_H_
#define _SPZPERIODICDATA_H_

#include "pzvec.h"
#include "tpzautopointer.h"

#include <map>


/** @brief Structure for storing periodic information
    All three structures have the same indexing
    @note For gmsh generated meshes, they are sorted by dimension
    of periodic entities (points, lines, surfaces)
*/
struct SPZPeriodicData{
  //! list of dependent periodic ids
  TPZVec<int> dep_mat_ids;
  //! list of independent periodic ids
  TPZVec<int> indep_mat_ids;
  //! list of maps of periodic nodes INDICES for the pairs (dep_id,indep_id)
  TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> nodes_map;
  //note that for periodic meshes node->Index() != node->Id() for dependent nodes
};
#endif /* _SPZPERIODICDATA_H_ */
