//$Id: pzgmesh.cpp,v 1.18 2004-06-23 15:57:29 phil Exp $

// -*- c++ -*-
/**File : pzgmesh.c

Method definition for class TPZGeoMesh.*/

#include "pzgmesh.h"
#include "pzvec.h"
//template class TPZVec<REAL>;
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
//#include "pzavlmap.h"

#include "pzelgt2d.h"
#include "pzelgq2d.h"
#include "pzelgt3d.h"
#include "pzelgpi3d.h"
#include <TPZRefPattern.h>
#include <tpzgeoelrefpattern.h>


TPZGeoMesh::TPZGeoMesh() : fElementVec(0), fNodeVec(0), fCosysVec(0),
  fBCElementVec(0) {

//  fName[0] = '\0';
  fReference = 0;
  fNodeMaxId = -1;
  fElementMaxId = -1;
}

TPZGeoMesh::~TPZGeoMesh() {
  CleanUp();
}

/**Delete element, nodes, Cosys, boundary elements and boundary nodes in list*/
void TPZGeoMesh::CleanUp() {
  int i, nel = fElementVec.NElements();
  for(i=0; i<nel; i++) {
    TPZGeoEl *el = fElementVec[i];
    if(el) delete el;
    fElementVec[i] = 0;
  }
  fElementVec.Resize(0);
  fElementVec.CompactDataStructure(1);
  fNodeVec.Resize(0);
  fNodeVec.CompactDataStructure(1);
  fCosysVec.Resize(0);
  fCosysVec.CompactDataStructure(1);
  fBCElementVec.Resize(0);
  fBCElementVec.CompactDataStructure(1);
//  fBCNodeVec.Resize(0);
//  fBCNodeVec.CompactDataStructure(1);

  map< int,map<string,TPZRefPattern *> >::iterator first = fRefPatterns.begin();
  map< int,map<string,TPZRefPattern *> >::iterator last = fRefPatterns.end();
  map< int,map<string,TPZRefPattern *> >::iterator  iter;
  for (iter = first; iter != last; iter++){
    map<string,TPZRefPattern *> &map_el = (*iter).second;

    map<string,TPZRefPattern *>::iterator name_first = map_el.begin();
    map<string,TPZRefPattern *>::iterator name_last = map_el.end();
    map<string,TPZRefPattern *>::iterator name_iter;    
    for(name_iter = name_first; name_iter != name_last; name_iter++){
      TPZRefPattern *refpat = (*name_iter).second;
      delete refpat;
    }
  }
}

void TPZGeoMesh::SetName (char *nm) {
  fName = nm;
//  if(nm != NULL) {
//    strncpy(fName,nm,62);
//    fName[62] = '\0';
//  }
}


void TPZGeoMesh::Print (ostream & out) {
  out << "\n\t\t GEOMETRIC TPZGeoMesh INFORMATIONS:\n\n";
  out << "TITLE-> " << fName << "\n\n";
  out << "number of nodes               = " << fNodeVec.NElements() << "\n";
  out << "number of elements            = " << fElementVec.NElements() << "\n";

  out << "\n\tGeometric Node Information:\n\n";
  int i;
  int nnodes = fNodeVec.NElements();
  for(i=0; i<nnodes; i++) {
    fNodeVec[i].Print(out);
    out << "\n";
  }
  out << "\n\tGeometric Element Information:\n\n";
  int nelem = fElementVec.NElements();
  for(i=0; i<nelem; i++) {
    if(fElementVec[i]) fElementVec[i]->Print(out);
    out << "\n";
  }

  out << "\nBoundary Element Information : \n\n";
  for(i=0; i<fBCElementVec.NElements () ; i++) {
	  fBCElementVec[i].Print (out);
  }
}

void TPZGeoMesh::GetNodePtr(TPZVec<int> &nos,TPZVec<TPZGeoNode *> &nodep) {

  int i,nnodes=nos.NElements();
  for(i=0;i<nnodes;i++) nodep[i]=&fNodeVec[nos[i]];
}

void  TPZGeoMesh::ResetReference() {

  TPZGeoEl *elp;
  int i,nelements=fElementVec.NElements();
  for(i=0;i<nelements;i++) {
    elp = fElementVec[i];
    if(elp) elp->ResetReference();
  }
  fReference = 0;
}

void TPZGeoMesh::RestoreReference(TPZCompMesh *cmesh) {

  ResetReference();
  fReference = cmesh;
  TPZGeoEl *gel;
  TPZCompEl *cel;
  int i,nelem = cmesh->ElementVec().NElements();
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(cel) {
      gel = cel->Reference();
      if(!gel) {
	PZError << "RestoreReference incomplete. Exist computational element with geometrical\n";
	PZError << "element not belongs to the current geometrical mesh.\n";
	return;
      }
      gel->SetReference(cel);
    }
  }
}

// GetBoundaryElements returns all elements beweeen NodFrom and NodTo counterclock wise
//		this method uses the connectivity of the elements
//		BuildConnectivity should be called to initialize the connectivity information
// 	this method will only work for grid with 2-D topology
//		the current version will only work for a grid with only one level
void TPZGeoMesh::GetBoundaryElements(int NodFrom, int NodTo,TPZStack<TPZGeoEl *> &ElementVec,TPZStack<int> &Sides) {
  // Find a first element whose first node on the side is NodFrom
//  TPZGeoEl *def = 0;
  //TPZAVLMap<int,TPZGeoEl *> elmap(def);
  map<int,TPZGeoEl *> elmap;
  int i,nelements=NElements();
  for(i=0;i<nelements;i++) {
    TPZGeoEl *el = fElementVec[i];
    if(el) elmap[el->Id()]=fElementVec[i];
  }

  int currentnode = NodFrom;
  TPZGeoEl *candidate = 0;
  int candidateside = 0;
  while(currentnode != NodTo) {
    // put all elements connected to currentnode in elmap, eliminate the elements
    //		from elmap which do not contain the node
    BuildElementsAroundNode(currentnode,elmap);
    // find, within elmap the element which has currentnode as its first boundary side
    //  	node
    FindElement(elmap, currentnode, candidate, candidateside);
    //	if the element found is already contained in the list, we have a circular list
    // if no element was found, the topology may not be two dimensional
    if(!candidate) break;
    int index = 0;
    int nelvec = ElementVec.NElements();
    while(index<nelvec && ElementVec[index] != candidate) index++;
    if(index <nelvec && Sides[index]==candidateside) break;
    ElementVec.Push(candidate);
    Sides.Push(candidateside);
    elmap.erase(elmap.begin(), elmap.end());//CleanUp();
    elmap[candidate->Id()] = candidate;
    // initialize the list in which to look for connected elements
    currentnode = candidate->SideNodeIndex(candidateside,1);
  }
}

// Find all elements in elmap or neighbour of elements in elmap which contain a node
//void TPZGeoMesh::BuildElementsAroundNode(int currentnode,TPZAVLMap<int,TPZGeoEl*> &elmap){
void TPZGeoMesh::BuildElementsAroundNode(int currentnode,map<int,TPZGeoEl*> &elmap){
  // first eliminate all elements which do not contain currentnode
  //TPZPix iel = elmap.First();
	map<int, TPZGeoEl *>::iterator ielprev,iel=elmap.begin();
  TPZGeoEl *el;
  int i;
  while(iel!=elmap.end()) {
    el = iel->second;
	ielprev=iel;
    iel++;//elmap.Next(iel);
    int numnode = el->NNodes();
    for(i=0; i< numnode; i++) {
      if(el->NodeIndex(i) == currentnode) break;
    }
    if(i == numnode) {
		elmap.erase(ielprev);
	}
  }
  iel = elmap.begin();
  while(iel!=elmap.end()) {
    el = iel->second;//elmap.Contents(iel);
    iel++;//elmap.Next(iel);
    int nside = el->NSides();
    for(int is=0; is<nside; is++) {
      TPZGeoElSide neigh = el->Neighbour(is);
      if(!neigh.Exists()) continue;
      int numnode = neigh.Element()->NNodes();
      for(i=0; i< numnode; i++) {
	if(neigh.Element()->NodeIndex(i) == currentnode){
	  if(elmap.find(neigh.Element()->Id())!=elmap.end()) {
	    // this should be implemented as a stack, so that we dont have to
	    // 	go through the list again each time
	    elmap[neigh.Element()->Id()] = neigh.Element();
	    iel = elmap.begin();
	  }
	  break;	// get out of the loop over the nodes
	}
      }
    }
  }
}

// find, within elmap the element which has currentnode as its first boundary side
//  	node
//void TPZGeoMesh::FindElement(TPZAVLMap<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {
void TPZGeoMesh::FindElement(map<int,TPZGeoEl *> &elmap,int currentnode,TPZGeoEl* &candidate,int &candidateside) {

  candidate = 0;
  //TPZPix iel = elmap.First();
  map<int , TPZGeoEl *>::iterator ielprev, iel = elmap.begin();
  while(iel!=elmap.end()) {
    TPZGeoEl *el = iel->second;//elmap.Contents(iel);
	ielprev=iel;
    iel++;//elmap.Next(iel);
    int ns = el->NSides();
    int is = el->NCornerNodes();
    for(; is < ns; is++) {
      TPZGeoElSide neigh = el->Neighbour(is);
      TPZGeoElSide father = el->Father2(is);
      if(!neigh.Exists() && !father.Exists() && el->SideNodeIndex(is,0) == currentnode) {
	candidate = el;
	candidateside = is;
	return;
      }
    }
  }
}

TPZGeoNode *TPZGeoMesh::FindNode(TPZVec<REAL> &co) {

  int i=0, in, nnodes = fNodeVec.NElements();
  while(i<nnodes && fNodeVec[i].Id() == -1) i++;
  if(i == nnodes) return 0;
  TPZGeoNode *gnkeep = &fNodeVec[i];
  REAL distkeep = 0.;
  for(in=0;in<3;in++)
    distkeep += (co[in]-(gnkeep->Coord(in)))*(co[in]-(gnkeep->Coord(in)));
  while(i< nnodes) {
    TPZGeoNode *gn = &fNodeVec[i];
    REAL dist = 0.;
    for(in=0;in<3;in++)
      dist += (co[in]-gn->Coord(in))*(co[in]-gn->Coord(in));
    if(dist < distkeep) {
      gnkeep = gn;
      distkeep = dist;
    }
    i++;
    while(i<nnodes && fNodeVec[i].Id() == -1) i++;
  }
  return gnkeep;
}

void TPZGeoMesh::BuildConnectivity()
{
  TPZVec<int> SideNum(NNodes(),-1);
  TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
  int nelem = NElements();
  int iel = 0;
  for(iel=0; iel<nelem; iel++)
    {
      TPZGeoEl *gel = fElementVec[iel];
      if(!gel) continue;
      int ncor = gel->NCornerNodes();
      int in;
      for(in=0; in<ncor; in++) {
	int nod = gel->NodeIndex(in);
	if(SideNum[nod] == -1)
	  {
	    NeighNode[nod] = gel;
	    SideNum[nod] = in;
	    if(gel->SideIsUndefined(in)) gel->SetSideDefined(in);
	  } else
	    {
	      TPZGeoElSide neigh(NeighNode[nod],SideNum[nod]);
	      TPZGeoElSide gelside(gel,in);
              if(!neigh.NeighbourExists(gelside))
              {
                neigh.SetConnectivity(gelside);
              }
	    }
      }
    }
  for(iel=0; iel<nelem; iel++)
    {
      TPZGeoEl *gel = fElementVec[iel];
      if(!gel) continue;
      int ncor = gel->NCornerNodes();
      int nsides = gel->NSides();
      int is;
      for(is=ncor; is<nsides; is++)
	{
	  if( gel->SideIsUndefined(is))
	    {
	      gel->SetSideDefined(is);
	      TPZGeoElSide gelside(gel,is);
	      TPZStack<TPZGeoElSide> neighbours;
	      gelside.ComputeNeighbours(neighbours);
	      int nneigh = neighbours.NElements();
	      int in;
	      for(in=0; in<nneigh; in++) {
		if(neighbours[in].Side() == -1) 
		  {
		    cout << "TPZGeoMesh::BuildConnectivity : Inconsistent mesh detected!\n";
		    continue;
		  }
		gelside.SetConnectivity(neighbours[in]);
	      }
	    }
	}
    }
}

void TPZGeoMesh::BuildConnectivity2() {

  TPZVec<int> SideNum(NNodes(),-1);
  TPZVec<TPZGeoEl *> NeighNode(NNodes(),0);
  int nelem = NElements();
  int iel = 0;
  while(iel<nelem && fElementVec[iel] == 0) iel++;

  long numsearch =1;
  // if there are no elements, do nothing
  while(iel < nelem) {
    TPZGeoEl *el = fElementVec[iel];
    int numsides = el->NSides();
    int side;
    for(side = 0;side<numsides;side++) {

      // check whether all entries in NeighNode are equal

	int equalnode = 1;
    int numsidenodes = el->NSideNodes(side);
    int sidenode = el->SideNodeIndex(side,0);
    TPZGeoEl *neigh = NeighNode[sidenode];
    int sidenumber = SideNum[sidenode];
    for(int sn = 0;sn < numsidenodes; sn++) {
		sidenode = el->SideNodeIndex(side,sn);
		if (neigh != NeighNode[sidenode]){
			equalnode=0;
			break;
		}
	}

    if(equalnode && neigh == 0) {
		if(el->SideIsUndefined(side)) {
			int elloaded = 0;
			for(int in=0; in<el->NNodes(); in++) {
				if(NeighNode[el->NodeIndex(in)] == el) elloaded = 1;
			}
			// this element is not loaded and its side is undefined

			// load the element side in the NeighNode vector
			for(int sn=0;!elloaded && sn < numsidenodes; sn++) {
				sidenode = el->SideNodeIndex(side,sn);
				NeighNode[sidenode] = el;
				SideNum[sidenode] = side;
			}
			numsearch++;
		}
	} else if(equalnode && side == sidenumber && neigh == el) {
	// unload the element side
		for(int sn=0;sn < numsidenodes; sn++) {
			sidenode = el->SideNodeIndex(side,sn);
			NeighNode[sidenode] = 0;
			SideNum[sidenode] = -1;
		}
	// if no neighbouring element was detected during the loop
	//    define the element side as undefined
	TPZGeoElSide neighbour = el->Neighbour(side);
	if(!neighbour.Exists()) el->SetSideDefined(side);
	numsearch++;
      } else if(equalnode && neigh != el) {
	// we found a neigbour
	TPZManVector<int> SideNodes(numsidenodes);
	// detect which side of the neigbour is loaded witin NeighNode
	for(int sn=0;sn < numsidenodes; sn++) {
		sidenode = el->SideNodeIndex(side,sn);
		SideNodes[sn] = sidenode;
	}
	int neighside = neigh->WhichSide(SideNodes);
	TPZGeoElSide neighbour(neigh,neighside);
	// WhichSide will tell the side number which contains the vector
	//    of node ids SideNodes
	if(neighbour.Side() != -1 && !el->NeighbourExists(side,neighbour)){
	  TPZGeoElSide(el,side).SetConnectivity(neighbour);
	  numsearch++;
	}
      }
    } // loop over the sides
    iel++;
    while(iel<nelem && fElementVec[iel] == 0) iel++;
    if(iel==nelem && numsearch) {
      numsearch = 0;
      iel = 0;
      while(iel<nelem && fElementVec[iel] == 0) iel++;
    }
  }
}

//Cedric : 03/03/99
TPZGeoEl *TPZGeoMesh::FindElement(int elid) {

	int nel = fElementVec.NElements();
   TPZGeoEl *gel = 0;
	for(int i=0;i<nel;i++) {
    	gel = fElementVec[i];
      if(gel && gel->Id()==elid) break;
   }
   return gel;
}

int TPZGeoMesh::ElementIndex(TPZGeoEl *gel){
	int i=0;
	int numel = ElementVec().NElements();
	while ( i < numel ) {
		if (ElementVec()[i] == gel) break;
		i++;
	}
	if(i<numel) return i;
	else return -1;
}

int TPZGeoMesh::NodeIndex(TPZGeoNode *nod){
	int i=0;
	int numel = NodeVec().NElements();
	while ( i < numel ) {
		if (&NodeVec()[i] == nod) break;
		i++;
	}
	return i;
}

#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzgeopoint.h"
#include "pzrefpoint.h"
#include "pzshapepoint.h"

TPZGeoEl *TPZGeoMesh::CreateGeoElement(MElementType type,
                                       TPZVec<int>& nodeindexes,
                                       int matid,
                                       int& index,
                                       int reftype){
  if (!reftype)  switch( type ){
    case 0://point
      return new TPZGeoElement< TPZShapePoint, TPZGeoPoint, TPZRefPoint>(
                              nodeindexes, matid, *this, index );
    case 1://line
      return new TPZGeoElement< TPZShapeLinear, TPZGeoLinear, TPZRefLinear>(
                              nodeindexes, matid, *this, index );
    case 2://triangle
      return new TPZGeoElement< TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >(
                              nodeindexes, matid, *this, index );
      //return new TPZGeoElT2d(nodeindexes,matid,*this);
    case 3://quadrilatera
      return  new TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >(
                              nodeindexes, matid, *this, index );
      //return new TPZGeoElQ2d(nodeindexes,matid,*this);
    case 4://tetraedra
      //return new TPZGeoElT3d(nodeindexes,matid,*this);
      return new TPZGeoElement< TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >(
                              nodeindexes, matid, *this, index );
    case 5:
      //return new TPZGeoElPi3d(nodeindexes,matid,*this);
      return new TPZGeoElement< TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >(
                              nodeindexes, matid, *this, index );
    case 6:
      return new TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >(
                              nodeindexes, matid, *this, index );
    case 7:
      return new TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >(
                              nodeindexes, matid, *this, index );
    default:
      PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
              << " type = " << type << endl;
    return NULL;
  } else {
    switch( type ){
      case 0://point
        return new TPZGeoElRefPattern<TPZShapePoint, TPZGeoPoint>(
                                nodeindexes, matid, *this, index, 0 );
      case 1://line
        return new TPZGeoElRefPattern<TPZShapeLinear, TPZGeoLinear>(
                                nodeindexes, matid, *this, index, 0 );
      case 2://triangle
        return new TPZGeoElRefPattern<TPZShapeTriang, TPZGeoTriangle>(
                                nodeindexes, matid, *this, index, 0 );
      case 3://quadrilatera
        return  new TPZGeoElRefPattern<TPZShapeQuad, TPZGeoQuad>(
                                nodeindexes, matid, *this, index, 0 );
      case 4://tetraedra
        return new TPZGeoElRefPattern<TPZShapeTetra, TPZGeoTetrahedra>(
                                nodeindexes, matid, *this, index, 0 );
      case 5://pyramid
        return new TPZGeoElRefPattern<TPZShapePiram, TPZGeoPyramid>(
                                nodeindexes, matid, *this, index, 0 );
      case 6://prism
        return new TPZGeoElRefPattern<TPZShapePrism, TPZGeoPrism>(
                                nodeindexes, matid, *this, index, 0 );
      case 7://cube
        return new TPZGeoElRefPattern<TPZShapeCube, TPZGeoCube>(
                                nodeindexes, matid, *this, index, 0 );
      default:
        PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
                << " type = " << type << endl;
      return NULL;
    }
  }
  return NULL;
}

void TPZGeoMesh::InsertRefPattern(TPZRefPattern *refpat){
  if (!refpat) {
    PZError << "TPZGeoMesh::InsertRefPattern ERROR : NULL refinement pattern! " << endl;
    return;
  }
  int eltype = refpat->Element(0)->Type();
  string name = refpat->GetName();
  map<int,map<string,TPZRefPattern *> >::iterator eltype_iter = fRefPatterns.find(eltype);

  if (eltype_iter == fRefPatterns.end()) {
    fRefPatterns [eltype][name] = refpat;
    return;
  }
  
  map <string,TPZRefPattern *>::iterator name_iter = fRefPatterns[eltype].find(name);
  if (name_iter != fRefPatterns[eltype].end()) return; // ja existe...
  fRefPatterns [eltype][name] = refpat;
}

TPZRefPattern * TPZGeoMesh::GetRefPattern(int eltype, string name){
  TPZRefPattern *refpat = 0;
  map<int,map<string,TPZRefPattern *> >::iterator eltype_iter = fRefPatterns.find(eltype);
  if (eltype_iter == fRefPatterns.end()) return refpat;
  map <string,TPZRefPattern *>::iterator name_iter = fRefPatterns[eltype].find(name);
  if (name_iter != fRefPatterns[eltype].end()) refpat = fRefPatterns[eltype][name];
  return refpat;

}

/*
TPZGeoEl* TPZGeoMesh::CreateGeoElement( MElementType type, int* nodeindexes,
					int matid, int& index )
{
   switch( type )
   {
      case 0://point
	 return new TPZGeoElement< TPZShapeLinear, TPZGeoPoint, TPZRefPoint >(
	    nodeindexes, matid, *this, index );

      case 1://line
	 return new TPZGeoElement< TPZShapeLinear, TPZGeoLinear, TPZRefLinear >(
	    nodeindexes, matid, *this, index );

      case 2://triangle
	 return new TPZGeoElement<
	    TPZShapeTriang, TPZGeoTriangle, TPZRefTriangle >(
	       nodeindexes, matid, *this, index );

      case 3://quadrilatera
	 return  new TPZGeoElement< TPZShapeQuad, TPZGeoQuad, TPZRefQuad >(
	    nodeindexes, matid, *this, index );

      case 4://tetraedra
	 return new TPZGeoElement<
	    TPZShapeTetra, TPZGeoTetrahedra, TPZRefTetrahedra >(
	       nodeindexes, matid, *this, index );

      case 5:
	 return new TPZGeoElement<
	    TPZShapePiram, TPZGeoPyramid, TPZRefPyramid >(
	       nodeindexes, matid, *this, index );

      case 6:
	 return new TPZGeoElement< TPZShapePrism, TPZGeoPrism, TPZRefPrism >(
	    nodeindexes, matid, *this, index );

      case 7:
	 return new TPZGeoElement< TPZShapeCube, TPZGeoCube, TPZRefCube >(
	    nodeindexes, matid, *this, index );

      default:
	 PZError << "TPZGeoMesh::CreateGeoElement type element not exists:"
		 << " type = " << type << endl;
	 return NULL;
   }

   return NULL;
}

void TPZGeoMesh::DeleteElement(TPZGeoEl *gel,int index){ 
  if(index < 0 || gel != fElementVec[index]){ 
    index = ElementIndex(gel);
    if(index < 0) {
      PZError << "TPZGeoMesh::DeleteElement index error\n"; 
      return; 
    }
  } 
  if(gel) delete gel; 
  fElementVec[index] = NULL; 
  fElementVec.SetFree(index); 
} 
*/
/** Verifies if the side based refinement pattern exists. If the refinement pattern doesn't exists return a Null refinement Pattern. */
TPZRefPattern * TPZGeoMesh::GetRefPattern (TPZGeoEl *gel, int side){
  int type = gel->Type();
  //construct the named for the refinement pattern.
  string name = "SIDE_";
  TPZGeoElSide gelside (gel,side);
  int dimension = gelside.Dimension();
  switch (dimension){
    case (0) :{
      name += "NODE_";
      break;
    }
    case (1) :{
      name += "RIB_";
      break;
    }
    case (2) :{
      name += "FACE_";
      break;
    }
    case (3) :{
      name += "VOLUME_";
      break;
    }
    default :{
      name += "RIB_";
      break;
    }
  }
  int size = 2;
  if (side/10 == 0){
    name += "0" ;
    size = 1;
  }
  char aux[256];
  sprintf(aux,"%d",side);
  name += aux;

  switch (type){
    case (0) : {
      name += "_POINT";
      break;
    }
    case (1) : {
      name += "_LINE";
      break;
    }
    case (2) : {
      name += "_TRIANGLE";
      break;
    }
    case (3) : {
      name += "_QUAD";
      break;
    }
    case (4) : {
      name += "_TETRA";
      break;
    }
    case (5) : {
      name += "_PYRAMID";
      break;
    }
    case (6) : {
      name += "_PRISM";
      break;
    }
    case (7) : {
      name += "_HEXA";
      break;
    }
    default:{
      PZError << "TPZGeoMesh::GetRefPattern ERROR : Undefined element type " << type << endl;
      return 0;
    }
  }
  return GetRefPattern(type,name);
}

int TPZGeoMesh::ClassId() const {
  return TPZGEOMESHID;
}

void TPZGeoMesh::Read(TPZStream &buf, void *context)
{
  TPZSaveable::Read(buf,context);
  buf.Read(&fName,1);
  ReadObjects(buf,fNodeVec,this);
  ReadObjectPointers(buf,fElementVec,this);
  ReadObjects(buf,this->fBCElementVec,this);
  buf.Read(&fNodeMaxId,1);
  buf.Read(&fElementMaxId,1);
  int ninterfacemaps;
  buf.Read(&ninterfacemaps,1);
  int c;
  for(c=0; c< ninterfacemaps; c++) 
  {
    int vals[3];
    buf.Read(vals,3);
    fInterfaceMaterials[pair<int,int>(vals[0],vals[1])]=vals[2];
    }
    BuildConnectivity();
}

void TPZGeoMesh::Write(TPZStream &buf, int withclassid)
{
  TPZSaveable::Write(buf,withclassid);
  buf.Write(&fName,1);
  WriteObjects(buf,fNodeVec);
  WriteObjectPointers(buf,fElementVec);
  WriteObjects(buf,fBCElementVec);
  buf.Write(&fNodeMaxId,1);
  buf.Write(&fElementMaxId,1);
  int ninterfacemaps = fInterfaceMaterials.size();
  buf.Write(&ninterfacemaps,1);
  InterfaceMaterialsMap::iterator it = fInterfaceMaterials.begin();
  while(it != fInterfaceMaterials.end())
  {
    int vals[3];
    vals[0] = (it->first).first;
    vals[1] = (it->first).second;
    vals[2] = it->second;
    buf.Write(vals,3);
  }
}

