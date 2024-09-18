#include "TPZGeoMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZGenGrid2D.h"
#include "TPZGenGrid3D.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"
#include "tpzpyramid.h"
#include "pzgeoelbc.h"
#include "pzlog.h"
#include "pzvec_extras.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzgeomeshtools");
#endif

void TPZGeoMeshTools::DividePyramidsIntoTetra(TPZGeoMesh *gmesh) {
    const char buf[] =
            "3     5  "
            "37     PyramidsIntoTetra  "
            "-1.    -1.    0.  "
            " 1.    -1.    0.  "
            " 1.     1.    0.  "
            "-1.     1.    0.  "
            " 0.     0.    1.  "
            "7     5     0     1     2     3     4 "
            "4     4     0     1     2     4 "
            "4     4     0     2     3     4 "
    ;
    std::istringstream str(buf);
    TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(str);
    TPZAutoPointer<TPZRefPattern> refpatFound = gRefDBase.FindRefPattern(refpat);
    if(!refpatFound){
        gRefDBase.InsertRefPattern(refpat);
    }
    else{
        refpatFound->SetName(refpat->Name());
    }
    refpat->InsertPermuted();

    TPZManVector<TPZGeoEl *, 2> sons;
    for(auto &gel : gmesh->ElementVec()){
        if(gel->Type() == EPiramide || gel->NSubElements() == 0){
            gel->SetRefPattern(refpat);
            gel->Divide(sons);
        }
    }
    gmesh->BuildConnectivity();
}

TPZGeoMesh *
TPZGeoMeshTools::CreateGeoMesh1D(const REAL minX, const REAL maxX, const int nEls,
                const TPZVec<int> &matids, bool createBoundEls){
#ifdef PZDEBUG
    if(minX > maxX){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Invalid values of minX and maxX\n";
        PZError<<"minX: "<<minX<<" maxX: "<<maxX;
        DebugStop();
    }
    if(nEls < 1){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Invalid value of nEls\n";
        PZError<<"nEls: "<<nEls;
        DebugStop();
    }
#endif
    if((createBoundEls && matids.size()!=3)||
       (!createBoundEls && matids.size()!=1)){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Wrong size of matids vec!\n";
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    const int matId{matids[0]};
    const auto nnodes = nEls + 1;
    const auto elSize = (maxX-minX)/nEls;
    gmesh->NodeVec().Resize(nnodes);
    for (auto i = 0 ; i < nnodes; i++)
      {
        const REAL pos = minX + i*elSize;
        TPZManVector <REAL,3> coord= {pos,0.,0.};
        coord[0] = pos;
        gmesh->NodeVec()[i].SetCoord(coord);
        gmesh->NodeVec()[i].SetNodeId(i);
      }
    TPZManVector <int64_t,2> nodeVec(2);
    int64_t id;
    for (auto iel = 0; iel < nEls; iel++)
      {
        nodeVec={iel,iel+1};
        gmesh->CreateGeoElement(EOned, nodeVec, matId, id);
        gmesh->ElementVec()[id];
      }
    if(createBoundEls){
        const int bc0{matids[1]},bc1{matids[2]};
        TPZManVector<int64_t,1> nodeVec = {0};
        gmesh->CreateGeoElement(EPoint, nodeVec, bc0, id);
        nodeVec = {nEls};
        gmesh->CreateGeoElement(EPoint, nodeVec, bc1, id);
    }
    gmesh->BuildConnectivity();
    return gmesh;
}

TPZGeoMesh *
TPZGeoMeshTools::CreateGeoMeshOnGrid(int dim, const TPZVec<REAL> &minX, const TPZVec<REAL> &maxX, const TPZVec<int> &matids,
                                     const TPZVec<int> nDivs, MMeshType meshType, bool createBoundEls) {
#ifdef PZDEBUG
    const REAL tol{ZeroTolerance()};
    if(dim !=3 && dim!= 2){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Dimension = "<<dim<<" is not supported. Aborting...\n";
        DebugStop();
    }
    if(MMeshType_Dimension(meshType) != dim){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"Element type "<<meshType<<" is not supported. Aborting...\n";
        DebugStop();
    }
    if(minX.NElements() != 3 || maxX.NElements() != 3){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"minX and maxX must have size = 3!\n";
        PZError<<"size(minX) = "<<minX.NElements()<<"\n"
               <<"size(maxX) = "<<maxX.NElements()<<"\n";
        DebugStop();
    }
    if(nDivs.NElements() != dim){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"size(nDivs) != dim !\n";
        PZError<<"dim = "<<minX.NElements()<<"\n"
               <<"size(nDivs) = "<<nDivs.NElements()<<"\n";
        DebugStop();
    }
    if(matids.size() != dim*2 + 1 && createBoundEls){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"The number of material ids must be equal to the number of boundaries + 1\n"
               <<"# of matids: "<<matids.size()<<"\n"
               <<"# of boundaries + 1: "<<dim*2+1<<"\n";
        DebugStop();
    } else if(matids.size() != 1 && !createBoundEls){
        PZError<<__PRETTY_FUNCTION__<<" error\n";
        PZError<<"When no boundary is created, the number of material ids must be equal to one\n"
               <<"# of matids: "<<matids.size()<<"\n";
        DebugStop();
    }
    if(dim == 2){
        if(nDivs[0] < 1 || nDivs[1] < 1){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
            PZError<<"nel x: "<<nDivs[0]<<"\n"<<"nel y: "<<nDivs[1]<<"\n";
            DebugStop();
        }else if(maxX[0] < tol || maxX[1] < tol){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
        else if(maxX[0] < minX[0] || maxX[1] < minX[1]){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"min x: "<<minX[0]<<"\n"<<"min y: "<<minX[1]<<"\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n";
            DebugStop();
        }
    }
    else{//dim == 3
        if(nDivs[0] < 1 || nDivs[1] < 1 || nDivs[2] < 1){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"The read number of grid divisions is not allowed. The parameters are:\n";
            PZError<<"nel x: "<<nDivs[0]<<"\n"<<"nel y: "<<nDivs[1]<<"\n"<<"nel z: "<<nDivs[2]<<"\n";
            DebugStop();
        }else if(maxX[0] < tol || maxX[1] < tol || maxX[2] < tol){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n"<<"max z: "<<maxX[2]<<"\n";
            DebugStop();
        }
        else if(maxX[0] < minX[0] || maxX[1] < minX[1] || maxX[2] < minX[2]){
            PZError<<__PRETTY_FUNCTION__<<" error\n";
            PZError<<"Dimensions of grid not allowed. The parameters are:\n";
            PZError<<"min x: "<<minX[0]<<"\n"<<"min y: "<<minX[1]<<"\n"<<"min z: "<<minX[2]<<"\n";
            PZError<<"max x: "<<maxX[0]<<"\n"<<"max y: "<<maxX[1]<<"\n"<<"max z: "<<maxX[2]<<"\n";
            DebugStop();
        }
    }
#endif
    return [&](){
        switch(dim){
            case 2:{
                TPZGeoMesh *gmesh = new TPZGeoMesh();
                TPZGenGrid2D gengrid(nDivs, minX, maxX);
                gengrid.SetElementType(meshType);
                gengrid.Read(gmesh, matids[0]);
                if(createBoundEls){
                    gengrid.SetBC(gmesh, 4, matids[1]);
                    gengrid.SetBC(gmesh, 5, matids[2]);
                    gengrid.SetBC(gmesh, 6, matids[3]);
                    gengrid.SetBC(gmesh, 7, matids[4]);
                }
                gmesh->BuildConnectivity();
                return gmesh;
            }
            case 3:{
                TPZGeoMesh *gmesh = nullptr;
                TPZGenGrid3D genGrid3D(minX,maxX,nDivs,meshType);
                gmesh = genGrid3D.BuildVolumetricElements(matids[0]);
                if(createBoundEls){
                    gmesh = genGrid3D.BuildBoundaryElements(matids[1],matids[2],matids[3],matids[4],matids[5],matids[6]);
                }

                gmesh->BuildConnectivity();
                return gmesh;
            }
            default:{
                PZError<<__PRETTY_FUNCTION__;
                PZError<<"Invalid dimension, returning nullptr\n";
                return (TPZGeoMesh*)nullptr;
            }
        }
    }();

}

TPZGeoMesh *TPZGeoMeshTools::CreateGeoMeshSingleEl(const MMeshType meshType,
                                                   const int matid,
                                                   const bool createBoundEls,
                                                   const int matidbc)
{
  auto gmesh = [=]() ->TPZGeoMesh* {
    switch (meshType) {
    case MMeshType::ETriangular:
      return CreateGeoMeshSingleElT<pztopology::TPZTriangle>(
          matid, createBoundEls, matidbc);
    case MMeshType::EQuadrilateral:
      return CreateGeoMeshSingleElT<pztopology::TPZQuadrilateral>(
          matid, createBoundEls, matidbc);
    case MMeshType::ETetrahedral:
      return CreateGeoMeshSingleElT<pztopology::TPZTetrahedron>(
          matid, createBoundEls, matidbc);
    case MMeshType::EPyramidal:
      return CreateGeoMeshSingleElT<pztopology::TPZPyramid>(
          matid, createBoundEls, matidbc);
    case MMeshType::EPrismatic:
      return CreateGeoMeshSingleElT<pztopology::TPZPrism>(
          matid, createBoundEls, matidbc);
    case MMeshType::EHexahedral:
      return CreateGeoMeshSingleElT<pztopology::TPZCube>(
          matid, createBoundEls, matidbc);
    case MMeshType::EHexaPyrMixed:
      return nullptr;
    case MMeshType::ENoType:
      return nullptr;
    }
    unreachable();
  }();
  return gmesh;
}

template <class TGEOM>
TPZGeoMesh *TPZGeoMeshTools::CreateGeoMeshSingleElT(const int matid,
                                                    const bool createBoundEls,
                                                    const int matidbc)
{
  constexpr auto nNodes{TGEOM::NCornerNodes};
  constexpr auto dim{TGEOM::Dimension};
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  gmesh->SetDimension(dim);
  // Auxiliary vector to store coordinates:
  TPZManVector<REAL, 3> coords(3, 0.);
  for (int iNode = 0; iNode < nNodes; iNode++) {
    TGEOM::ParametricDomainNodeCoord(iNode, coords);
    coords.Resize(3);
    auto newindex = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[newindex].Initialize(coords, *gmesh);
  }
  TPZManVector<int64_t, TGEOM::NCornerNodes> nodesIdVec(TGEOM::NCornerNodes,
                                                        -1);
  for (int i = 0; i < TGEOM::NCornerNodes; i++)
    nodesIdVec[i] = i;
  int64_t index{-1};
  TPZGeoEl *gel =
      gmesh->CreateGeoElement(TGEOM::Type(), nodesIdVec, matid, index, 0);
  gmesh->BuildConnectivity();

  if (createBoundEls) {
    const auto fside = gel->FirstSide(dim - 1);
    const auto lastside = fside + gel->NSides(dim - 1);
    for (int iside = fside; iside < lastside; iside++)
      TPZGeoElBC(gel, iside, matidbc);
  }
  return gmesh;
}

TPZGeoEl* TPZGeoMeshTools::FindElementByMatId(TPZGeoMesh* gmesh, TPZVec<REAL> &x, TPZVec<REAL> & qsi, int64_t & InitialElIndex, const std::set<int>& matids) {
    if (InitialElIndex < 0 || InitialElIndex >= gmesh->NElements()) {
        DebugStop();
    }
    
    // FindApproxElement
    FindCloseElement(gmesh, x, InitialElIndex, matids);
    TPZGeoEl * gel = gmesh->ElementVec()[InitialElIndex]->LowestFather();
 
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "x coordinate " << x << std::endl;
        sout << "element index " << gel->Index() << std::endl;
        sout << "coordinates of the corner nodes of the element\n";
        for (int i=0; i<gel->NNodes(); i++) {
            TPZGeoNode *ptr = gel->NodePtr(i);
            for (int ic=0; ic<3; ic++) {
                sout << ptr->Coord(ic) << " ";
            }
            sout << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    qsi.Resize(gel->Dimension(), 0.);
    qsi.Fill(0.);
    
    if(matids.find(gel->MaterialId()) == matids.end())
    {
        DebugStop();
    }
    REAL zero;
    ZeroTolerance(zero);
    
    std::set<TPZGeoEl *> tested;
    // this method will call ComputeXInverse if the element dimension != 3
    bool memberQ = gel->ComputeXInverse(x, qsi,zero);
    if(memberQ)
    {
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            LOGPZ_DEBUG(logger, "Going into the FindSubElement alternative")
        }
#endif
        gel = gmesh->FindSubElement(gel, x, qsi, InitialElIndex);
        return gel;
    }
    tested.insert(gel);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Looking for x = " << x << std::endl;
        sout << "Tried out gel index " << gel->Index() << std::endl;
        sout << "found qsi " << qsi << std::endl;
        TPZManVector<REAL> locx(3);
        gel->X(qsi, locx);
        sout << "found x co " << locx << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZManVector<REAL,3> projection(gel->Dimension());
    int side = -1;
    side = gel->ProjectInParametricDomain(qsi, projection);

    
    TPZGeoElSide mySide(gel,side);
    TPZManVector<REAL,3> xclose(3);
    gel->X(projection, xclose);
//    mySide.X(projection, xclose);
    REAL mindist = 0.;
    for (int i=0; i<3; i++) {
        mindist += (xclose[i]-x[i])*(xclose[i]-x[i]);
    }
    TPZGeoElSide bestside(mySide);
    TPZManVector<REAL,3> bestproj(projection);
    
    
    TPZStack<TPZGeoElSide> allneigh;
    TPZGeoElSide neighSide(mySide.Neighbour());
    while (neighSide != mySide) {
        if (matids.find(neighSide.Element()->MaterialId()) != matids.end() && tested.find(neighSide.Element()) == tested.end()) {
            allneigh.Push(neighSide);
        }
        neighSide = neighSide.Neighbour();
    }
    
    while (allneigh.size())
    {
        TPZGeoElSide  gelside = allneigh.Pop();
        TPZGeoEl *locgel = gelside.Element();
        if (tested.find(locgel) != tested.end()) {
            continue;
        }
        qsi.Fill(0.);
        if(locgel->ComputeXInverse(x, qsi, zero*100.) == true)
        {
#ifdef PZ_LOG
            if (logger.isDebugEnabled())LOGPZ_DEBUG(logger, "FOUND ! Going into the FindSubElement alternative")
#endif
            gel = gmesh->FindSubElement(locgel, x, qsi, InitialElIndex);
            return gel;
        }
        
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "Looking for x = " << x << std::endl;
            sout << "Tried out gel index " << locgel->Index() << std::endl;
            sout << "found qsi " << qsi << std::endl;
            TPZManVector<REAL> locx(3);
            gel->X(qsi, locx);
            sout << "found x co " << locx << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

        
        tested.insert(locgel);
        int side = -1;
        side = locgel->ProjectInParametricDomain(qsi, projection);
        TPZManVector<REAL,3> xclose(3);
        TPZGeoElSide locgelside(locgel,side);
        locgel->X(projection, xclose);
        REAL dist = 0.;
        for (int i=0; i<3; i++) {
            dist += (xclose[i]-x[i])*(xclose[i]-x[i]);
        }
        if (dist < mindist) {
            mindist = dist;
            bestside = locgelside;
            bestproj = projection;
        }
        mySide = TPZGeoElSide(locgel,side);
        TPZGeoElSide neighSide(mySide.Neighbour());
        while (neighSide != mySide) {
            if (matids.find(neighSide.Element()->MaterialId()) != matids.end() && tested.find(neighSide.Element()) == tested.end()) {
                allneigh.Push(neighSide);
            }
            neighSide = neighSide.Neighbour();
        }
    }
    TPZGeoEl *bestgel = bestside.Element();
    qsi = bestproj;
    gel = gmesh->FindSubElement(bestgel, x, qsi, InitialElIndex);
    return gel;
}

TPZGeoEl * TPZGeoMeshTools::FindCloseElement(TPZGeoMesh* gmesh, TPZVec<REAL> &x, int64_t & InitialElIndex, const std::set<int>& matids)
{
    TPZManVector<REAL,3> xcenter(3);
    TPZManVector<TPZManVector<REAL, 3>, 8> cornercenter;
    // what happens if gelnext == 0 or if InitialIndex is out of scope??
    if (InitialElIndex >= gmesh->NElements()) {
        DebugStop();
    }
    TPZGeoEl *gelnext = gmesh->ElementVec()[InitialElIndex];
    if (matids.find(gelnext->MaterialId()) == matids.end()) {
        DebugStop();
    }
    
    if (gelnext == 0) {
        DebugStop();
    }
    
    TPZGeoEl *gel = 0;
    REAL geldist;
    std::map<REAL,int> cornerdist;
    while (gel != gelnext) {
        gel = gelnext;
        gelnext = 0;
        cornerdist.clear();
        // compute the corner coordinates
        for (int ic=0; ic<gel->NCornerNodes(); ic++) {
            TPZManVector<REAL,3> xcorner(3);
            gel->NodePtr(ic)->GetCoordinates(xcorner);
            REAL dist = 0.;
            for (int i=0; i<3; i++) {
                dist += (x[i]-xcorner[i])*(x[i]-xcorner[i]);
            }
            dist = sqrt(dist);
            cornerdist[dist]=ic;
        }
        // compute the distance of the center of the element
        {
            geldist = 0.;
            TPZManVector<REAL,3> xcenter(3);
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            gelside.CenterX(xcenter);
            geldist = dist(x,xcenter);
        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "FindClosestElement Tried element index " << gel->Index() << std::endl;
            sout << "Distance from the center " << geldist << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // find the closest corner node
        REAL closestcorner = cornerdist.begin()->first;
        // if the center node is closer than the cornernode, return the element
        if (geldist < closestcorner || closestcorner < 1.e-15) {
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                std::stringstream sout;
                sout << "Distance from the closest corner " << closestcorner << "bailing out " << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            InitialElIndex = gel->Index();
            return gel;
        }
        // look for all neighbours of the corner side
        TPZGeoElSide gelside(gel,cornerdist.begin()->second);
        TPZGeoElSide neighbour = gelside.Neighbour();
        std::map<REAL,int> distneigh;
        while (neighbour != gelside)
        {
            int neighmatid = neighbour.Element()->MaterialId();
            if (matids.find(neighmatid) != matids.end())
            {
                TPZManVector<REAL,3> center(3);
                TPZGeoElSide centerneigh(neighbour.Element(),neighbour.Element()->NSides()-1);
                centerneigh.CenterX(center);
                REAL distcenter = dist(center,x);
                distneigh[distcenter] = neighbour.Element()->Index();
            }
            neighbour = neighbour.Neighbour();
        }
        // choose the element whose center is closest to the coordinate
        REAL gelnextdist = 0.;
        if (distneigh.size() == 0) {
            gelnext = gel;
            gelnextdist = geldist;
        }
        else {
            gelnext = gmesh->ElementVec()[distneigh.begin()->second];
            gelnextdist = distneigh.begin()->first;
        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Closest element index " << gelnext->Index() << std::endl;
            sout << "Distance from the center " << gelnextdist << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        // return if its center distance is larger than the distance of the current element
        if (geldist <= gelnextdist) {
            InitialElIndex = gel->Index();
            return gel;
        }
    }
    InitialElIndex = gel->Index();
    return gel;
}

void
TPZGeoMeshTools::FindPeriodicElements(
  TPZGeoMesh *gmesh,
  const TPZVec<int> &dep_ids, const TPZVec<int> &indep_ids,
  const TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &periodic_nodes,
  TPZVec<TPZAutoPointer<std::map<int64_t,int64_t>>> &el_map)
{

    if((dep_ids.size() != indep_ids.size()) || dep_ids.size() != periodic_nodes.size()){
        DebugStop();
    }
    
    const int n_periodic_reg = dep_ids.size();
    el_map.Resize(n_periodic_reg,nullptr);

    /**
       now we want to map periodic ELEMENTS
       IMPORTANT: since we have changed the node ids,
       now we no longer have index == id
    */
    for(auto ireg = 0; ireg < n_periodic_reg; ireg++){
        el_map[ireg] = new std::map<int64_t,int64_t>{};
        const auto depmatid = dep_ids[ireg];
        const auto indepmatid = indep_ids[ireg];
        const auto &node_map = *periodic_nodes[ireg];
        for (auto depel : gmesh->ElementVec()) {
            if(!depel){continue;}
            if(depel->MaterialId() != depmatid){continue;}
            
            //we have found a dependent el
            const auto nnodes = depel->NNodes();
            const auto nsides = depel->NSides();
            const auto dep_type = depel->Type();
            constexpr int max_nnodes{8};
            TPZManVector<int64_t,max_nnodes> mapped_nodes(nnodes);
            for (auto in = 0; in < nnodes; in++) {
                const auto depnode = depel->NodeIndex(in);
                mapped_nodes[in] = node_map.at(depnode);
            }
            TPZGeoEl* indepel{nullptr};
            const int nelem = gmesh->ElementVec().NElements();
            for(auto gel : gmesh->ElementVec()){
                if(!gel){continue;}
                if(indepel){break;}
                if (gel->MaterialId() != indepmatid ||
                    gel->Type() != dep_type) {continue;}
            
                TPZManVector<int64_t,max_nnodes> indepnodes(nnodes);
                gel->GetNodeIndices(indepnodes);
                bool samenodes = true;
                for (auto in = 0; in < nnodes && samenodes; in++) {
                    const auto indepnode = indepnodes[in];
                    const bool hasnode =
                        std::find(mapped_nodes.begin(), mapped_nodes.end(),
                                  indepnode) != mapped_nodes.end();
                    samenodes = samenodes && hasnode;
                }
                if (!samenodes) {continue;}
                indepel=gel;
            }
            if(!indepel){
                PZError<<__PRETTY_FUNCTION__
                       <<"\nCould not find periodic element for:\n"
                       <<"el index "<<depel->Index()<<" matid:  "<<depmatid<<std::endl;
                DebugStop();
            }
#ifdef PZDEBUG
            if(el_map[ireg]->count(depel->Index())!=0){
                DebugStop();
            }
#endif
            el_map[ireg]->insert({depel->Index(),indepel->Index()});
        }
    }
}



#define CREATE_TEMPL(TTOPOL) \
    template \
    TPZGeoMesh * \
    TPZGeoMeshTools::CreateGeoMeshSingleElT<TTOPOL>(const int,  \
                                                    const bool, \
                                                    const int);

CREATE_TEMPL(pztopology::TPZLine);
CREATE_TEMPL(pztopology::TPZTriangle);
CREATE_TEMPL(pztopology::TPZQuadrilateral);
CREATE_TEMPL(pztopology::TPZTetrahedron);
CREATE_TEMPL(pztopology::TPZCube);
CREATE_TEMPL(pztopology::TPZPrism);
CREATE_TEMPL(pztopology::TPZPyramid);
