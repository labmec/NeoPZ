/**
 * @file BlendUnitTest.cpp
 * @brief Define a Unit Test using Boost for blend elements over all kind of geometries available in the library
 *
 */
#include <iostream>

#include "pzlog.h"
#include "TPZRefPatternDataBase.h"

#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeoprism.h"
#include "pzgeopyramid.h"

#include "pzgeoelrefless.h"

#include "tpzquadraticline.h"

#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"

#include "tpzquadratictetra.h"
#include "tpzquadraticcube.h"
#include "tpzquadraticprism.h"
#include "tpzquadraticpyramid.h"

#include "tpzgeoblend.h"

#include "TPZVTKGeoMesh.h"

#include "tpzgeoelrefpattern.h"

#include "tpzarc3d.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.testgeom"));
#endif

#ifdef _AUTODIFF
#include "fad.h"
#endif

// Using Unit Test of the Boost Library
#ifdef USING_BOOST

#ifndef WIN32
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MAIN pz blend_tests tests

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "boost/test/output_test_stream.hpp"

#endif

//#define NOISY_BLEND //outputs x and grad comparisons
#define NOISYVTK _BLEND//prints all elements in .vtk format

std::string dirname = PZSOURCEDIR;

#ifdef USING_BOOST


BOOST_AUTO_TEST_SUITE(blend_tests)

    namespace blendtest{
        const int pOrder = 10;
        const REAL tol = 1e-8;
        template <class TGeo>
        void CompareQuadraticAndBlendEls();
        void CreateGeoMesh2D(int nDiv);
        void CreateGeoMesh3D(int nDiv);
    }

    BOOST_AUTO_TEST_CASE(geoblend_tests) {
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoQuad>();
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoTriangle>();
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoTetrahedra>();
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoCube>();
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoPrism>();
        blendtest::CompareQuadraticAndBlendEls<pzgeom::TPZGeoPyramid>();

        gRefDBase.InitializeUniformRefPattern(EOned);
        gRefDBase.InitializeUniformRefPattern(ETriangle);
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        {
            const int nDiv = 4;
            blendtest::CreateGeoMesh2D(nDiv);
        }

        gRefDBase.InitializeUniformRefPattern(ETetraedro);
        gRefDBase.InitializeUniformRefPattern(EPiramide);
        gRefDBase.InitializeUniformRefPattern(EPrisma);
        gRefDBase.InitializeUniformRefPattern(ECube);
        {
            const int nDiv = 3;
            blendtest::CreateGeoMesh3D(nDiv);
        }
    }

    namespace blendtest{
        template <class TGeo>
        void CompareQuadraticAndBlendEls() {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            auto elType = TGeo::Type();

            const int64_t nCornerNodes = TGeo::NCornerNodes;
            const int64_t nEdges = TGeo::NumSides(1);
            const int64_t dim = TGeo::Dimension;
            const REAL sphereRadius = 1;

            TPZVec<REAL> phiPts(nCornerNodes,-1),thetaPts(nCornerNodes,-1); //r is always equal to sphereRadius
            switch(elType){
                case ETriangle:
                    for(int i = 0; i < nCornerNodes; i++){
                        thetaPts[i] = M_PI/2;
                        phiPts[i] = i * M_PI/3;
                    }
                    break;
                case EQuadrilateral:
                    for(int i = 0; i < nCornerNodes; i++){
                        thetaPts[i] = M_PI/2;
                        phiPts[i] = i * M_PI/2;
                    }
                    break;
                case ETetraedro:
                    for(int i = 0; i < 3; i++){
                        thetaPts[i] = 3 * M_PI/4;
                        phiPts[i] = i * M_PI/3;
                    }
                    thetaPts[3] = 0;
                    phiPts[3] = 0;
                    break;
                case EPrisma:
                    for(int i = 0; i < 3; i++){
                        thetaPts[i] = 3 * M_PI/4;
                        phiPts[i] = i * M_PI/3;
                    }
                    for(int i = 0; i < 3; i++){
                        thetaPts[i+3] = 1 * M_PI/4;
                        phiPts[i+3] = i * M_PI/3;
                    }
                    break;
                case ECube:
                    for(int i = 0; i < 4; i++){
                        thetaPts[i] = 3 * M_PI/4;
                        phiPts[i] = i * M_PI/2;
                    }
                    for(int i = 0; i < 4; i++){
                        thetaPts[4+i] = 1 * M_PI/4;
                        phiPts[4+i] = i * M_PI/2;
                    }
                    break;
                case EPiramide:
                    for(int i = 0; i < 4; i++){
                        thetaPts[i] = 3 * M_PI/4;
                        phiPts[i] = i * M_PI/2;
                    }
                    thetaPts[4] = 0;
                    phiPts[4] = 0;
                    break;
                default:
                    DebugStop();
            }


            TPZGeoMesh *gmesh = new TPZGeoMesh();

            TPZManVector<REAL,3> coordsOffset(3,0);
            coordsOffset[0] = 2;
            TPZManVector<int, 12> midSideNodesIndexVec(nEdges, -1);
            const int nNodesOriginalMesh = nCornerNodes + nEdges;
            {
                TPZManVector<REAL, 3> coord(3, 0.);

                ///CREATE NODES FOR QUADRATIC ELEMENT
                for (int64_t i = 0; i < nCornerNodes; i++) {
                    const REAL theta = ((REAL) i / nCornerNodes) * 2 * M_PI;
                    coord[0] = sphereRadius * sin(thetaPts[i]) * cos(phiPts[i]);
                    coord[1] = sphereRadius * sin(thetaPts[i]) * sin(phiPts[i]);
                    coord[2] = sphereRadius * cos(thetaPts[i]);
//                std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                    const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
                }//quadraticElCornerNodes


                //CREATE MIDSIDE NODES
                TPZStack<int> edgeStack;
                TGeo::LowerDimensionSides(TGeo::NSides - 1, edgeStack, 1);
                for (int64_t edgeIndex = 0;
                     edgeIndex < edgeStack.NElements(); edgeIndex++) {
                    const int64_t edge = edgeStack[edgeIndex];
                    const int64_t nNodesSide = TGeo::NSideNodes(edge);
                    coord[0] = 0;
                    coord[1] = 0;
                    coord[2] = 0;
                    for (int64_t node = 0; node < nNodesSide; node++) {
                        const int64_t nodeIndex = TGeo::SideNodeLocId(edge, node);
                        coord[0] += gmesh->NodeVec()[nodeIndex].Coord(0);
                        coord[1] += gmesh->NodeVec()[nodeIndex].Coord(1);
                        coord[2] += gmesh->NodeVec()[nodeIndex].Coord(2);
                    }
                    const REAL norm = sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
                    coord[0] *= sphereRadius / norm;
                    coord[1] *= sphereRadius / norm;
                    coord[2] *= sphereRadius / norm;
                    const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
                    midSideNodesIndexVec[edgeIndex] = newindex;
                }

                ///CREATE NODES FOR BLEND ELEMENT
                for (int i = 0; i < nNodesOriginalMesh; i++) {
                    coord[0] = gmesh->NodeVec()[i].Coord(0) + coordsOffset[0];
                    coord[1] = gmesh->NodeVec()[i].Coord(1) + coordsOffset[1];
                    coord[2] = gmesh->NodeVec()[i].Coord(2) + coordsOffset[2];
                    const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                    gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
                }
            }

            const int matIdVol = 1, matIdSphere = 2;
            TPZManVector<int64_t,8> nodesIdVec(1);
            int64_t elId = 0;

            ///////CREATE MAX DIM ELEMENTS
            nodesIdVec.Resize(nCornerNodes + nEdges);
            for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i;

            TPZGeoEl *quadraticEl = nullptr;
            switch(elType){
                case ETriangle:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                case EQuadrilateral:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                case ETetraedro:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticTetra>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                case EPrisma:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticPrism>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                case ECube:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticCube>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                case EPiramide:
                    quadraticEl = new TPZGeoElRefLess<pzgeom::TPZQuadraticPyramid>(elId,nodesIdVec,matIdVol,*gmesh);
                    break;
                default:
                    DebugStop();
                    break;
            }
            nodesIdVec.Resize(nCornerNodes);
            for(int i = 0; i < nodesIdVec.size(); i++ ) nodesIdVec[i] = i + nNodesOriginalMesh;
            TPZGeoEl *blendEl = new TPZGeoElRefLess<pzgeom::TPZGeoBlend<TGeo>>(elId,nodesIdVec,matIdVol,*gmesh);

#ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "gmesh_" + TGeo::TypeName();
                const size_t strlen = meshFileName.length();
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());
                gmesh->Print(outTXT);
                outTXT.close();
            }
#endif
            ///////CREATE MAX DIM - 1 ELEMENTS FOR BLENDING
            {
                TPZStack<int> containedSides;
                TGeo::LowerDimensionSides(TGeo::NSides - 1, containedSides, TGeo::Dimension - 1);
                for (int64_t subSideIndex = 0;
                     subSideIndex < containedSides.NElements(); subSideIndex++){
                    const int64_t subSide = containedSides[subSideIndex];
                    const int64_t nNodesSide = TGeo::NSideNodes(subSide);
                    TPZStack<int> containedEdges;
                    TGeo::LowerDimensionSides(subSide, containedEdges, 1);
                    if(TGeo::SideDimension(subSide) == 1)   containedEdges.Push(subSide);
                    const int64_t nEdgesSide = containedEdges.NElements();
                    nodesIdVec.Resize(nNodesSide+nEdgesSide);
                    TPZManVector<int64_t,8> nodesIdVecQuad(nNodesSide+nEdgesSide);
                    for (int64_t node = 0; node < nNodesSide; node++){
                        const int64_t nodeId = TGeo::SideNodeLocId(subSide,node);
                        nodesIdVec[node] = blendEl->Node(nodeId).Id();
                        nodesIdVecQuad[node] = quadraticEl->Node(nodeId).Id();
//                    std::cout<<nodesIdVec[node]<<std::endl;
                    }
                    for (int64_t edge = 0; edge < nEdgesSide; edge++){
                        const int64_t edgeId = containedEdges[edge];
                        nodesIdVec[nNodesSide + edge] = nNodesOriginalMesh + midSideNodesIndexVec[edgeId- nCornerNodes];
                        nodesIdVecQuad[nNodesSide + edge] = midSideNodesIndexVec[edgeId- nCornerNodes];
//                    std::cout<<nodesIdVec[nNodesSide + edge]<<std::endl;
                    }
                    switch (TGeo::Type(subSide)){
                        case EOned:
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticLine>(nodesIdVec,matIdSphere,*gmesh);
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticLine>(nodesIdVecQuad,matIdSphere,*gmesh);
                            break;
                        case ETriangle:
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(nodesIdVec,matIdSphere,*gmesh);
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticTrig>(nodesIdVecQuad,matIdSphere,*gmesh);
                            break;
                        case EQuadrilateral:
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(nodesIdVec,matIdSphere,*gmesh);
                            new TPZGeoElRefLess<pzgeom::TPZQuadraticQuad>(nodesIdVecQuad,matIdSphere,*gmesh);
                            break;
                        default:
                            DebugStop();
                    }
                }
            }
            gmesh->BuildConnectivity();


#ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "gmesh_" + TGeo::TypeName();
                const size_t strlen = meshFileName.length();
                meshFileName.append(".vtk");
                std::ofstream outVTK(meshFileName.c_str());
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());

                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
                gmesh->Print(outTXT);
                outTXT.close();
                outVTK.close();
            }
#endif
            //X COMPARE
            {
                TPZManVector<REAL,3> xi;
                REAL notUsedHere = -1;
                uint64_t errorsInterior = 0;
                TPZManVector<uint64_t,18> errorsSide(TGeo::NSides - TGeo::NNodes - 1,0);
                std::cout<<"Element: "<<MElementType_Name(TGeo::Type())<<std::endl;
                uint64_t nPoints = 0;
                for(int iSide = TGeo::NNodes; iSide < TGeo::NSides; iSide ++){
                    auto intRule = blendEl->CreateSideIntegrationRule(iSide, pOrder);
                    xi.Resize(dim,0);
                    for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                        bool hasAnErrorOccurred = false;
                        nPoints++;
                        TPZManVector<REAL,3> xiSide(TGeo::SideDimension(iSide),0);
                        intRule->Point(iPt,xiSide,notUsedHere);
                        auto transf = TGeo::TransformSideToElement(iSide);
                        transf.Apply(xiSide,xi);
                        TPZManVector<REAL,3> xBlend(3);
                        blendEl->X(xi, xBlend);
                        TPZManVector<REAL,3> xQuad(3);
                        quadraticEl->X(xi, xQuad);
                        std::ostringstream xBlendM, xQuadM;
                        xBlendM<<"x_blend:"<<std::endl;
                        xQuadM<<"x_quad:"<<std::endl;
                        const auto VAL_WIDTH = 15;
                        REAL diff = 0;
                        for(int i = 0; i < xBlend.size(); i++){
                            xBlendM<<std::setw(VAL_WIDTH) << std::right<<xBlend[i] - coordsOffset[i]<<"\t";
                            xQuadM <<std::setw(VAL_WIDTH) << std::right<<xQuad[i] <<"\t";
                            diff += (xBlend[i] - coordsOffset[i] - xQuad[i])*(xBlend[i] - coordsOffset[i] - xQuad[i]);
                        }
                        if(diff > tol){
                            if(!hasAnErrorOccurred && iSide < TGeo::NSides - 1) errorsSide[iSide - TGeo::NNodes]+=1;
                            if(!hasAnErrorOccurred && iSide == TGeo::NSides - 1) errorsInterior++;
                            hasAnErrorOccurred = true;
                        }
                        BOOST_CHECK(!hasAnErrorOccurred);
                        if(hasAnErrorOccurred){
                            std::cout<<std::flush;
                            std::cout<<xBlendM.str()<<std::endl;
                            if(iSide < TGeo::NSides - 1){
                                auto neigh = blendEl->Neighbour(iSide);
                                if(neigh.Id() != blendEl->Id()) {
                                    TPZGeoElSide thisside(blendEl, iSide);
                                    auto neighTransf = thisside.NeighbourSideTransform(neigh);
                                    TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                                    neighTransf.Apply(xiSide, neighXi);
                                    TPZManVector<REAL, 3> xNeigh(3);
                                    neigh.X(neighXi, xNeigh);
                                    std::cout<<"x_neigh_blend:"<<std::endl;
                                    for (int i = 0; i < xNeigh.size(); i++) {
                                        std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i] - coordsOffset[i]
                                                  << "\t";
                                    }
                                    std::cout<<std::endl;
                                }
                            }
                            std::cout<<xQuadM.str()<<std::endl;
                            if(iSide < TGeo::NSides - 1){
                                auto neigh = quadraticEl->Neighbour(iSide);
                                if(neigh.Id() != quadraticEl->Id()) {
                                    TPZGeoElSide thisside(quadraticEl, iSide);
                                    auto neighTransf = thisside.NeighbourSideTransform(neigh);
                                    TPZManVector<REAL, 3> neighXi(neigh.Dimension(), 0);
                                    neighTransf.Apply(xiSide, neighXi);
                                    TPZManVector<REAL, 3> xNeigh(3);
                                    neigh.X(neighXi, xNeigh);
                                    std::cout<<"x_neigh_quad:"<<std::endl;
                                    for (int i = 0; i < xNeigh.size(); i++) {
                                        std::cout << std::setw(VAL_WIDTH) << std::right << xNeigh[i]<< "\t";
                                    }
                                    std::cout<<std::endl;
                                }
                            }
                            std::cout<<"diff :"<<diff<<std::endl;
                        }
                    }
#ifdef NOISY_BLEND
                    if(iSide == TGeo::NSides - 1){
                        std::cout<<"\tNumber of interior points: "<<intRule->NPoints()<<"\tErrors: "<<errorsInterior<<std::endl;
                    }else{
                        std::cout<<"\tNumber of side points for side "<<iSide<<": "<<intRule->NPoints()<<"\tErrors: "<<errorsSide[iSide - TGeo::NNodes]<<std::endl;
                    }
#endif
                }
                uint64_t errorsTotal = errorsInterior;
                for(int i = 0; i< errorsSide.size(); i++){
                    errorsTotal +=errorsSide[i];
                }
                std::cout<<"\tNumber of points: "<<nPoints<<"\tErrors: "<<errorsTotal<<std::endl;
            }
            delete gmesh;
        }

        void CreateGeoMesh2D(int nDiv) {
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            TPZGeoMesh *gmesh = new TPZGeoMesh();

            TPZVec<REAL> coord(3, 0.);
            for (int64_t i = 0; i < 4; i++) {
                coord[0] = i / 2 == i % 2 ? -1 : 1;
                coord[1] = -1 + 2 * (i / 2);
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
            }//quad nodes

            TPZManVector<int64_t, 4> nodesIdArcsVec(4);
            for (int64_t i = 0; i < 4; i++) {
                const int xOff = i % 2;
                const int yOff = 1 - xOff;
                const int xSign = 1 - 2 * (i / 2);
                const int ySign = -1 * xSign;
                coord[0] = xSign * xOff * M_SQRT2;
                coord[1] = ySign * yOff * M_SQRT2;
//            std::cout<<"xOff:\t"<<xOff<<"\txSign:"<<xSign<<std::endl;
//            std::cout<<"yOff:\t"<<yOff<<"\tySign:"<<ySign<<std::endl;
//            std::cout<<coord[0]<<"\t"<<coord[1]<<std::endl;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord, *gmesh);
                nodesIdArcsVec[i] = newindex;
            }//quad nodes

            int matIdQuad = 1, matIdArc = 2;
            TPZManVector<int64_t, 4> nodesIdVec(1);
            int64_t elId = 0;

            nodesIdVec.Resize(4);
            TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > *quadEl = nullptr;
            {
                nodesIdVec[0] = 0;
                nodesIdVec[1] = 1;
                nodesIdVec[2] = 2;
                nodesIdVec[3] = 3;
                quadEl = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> >(nodesIdVec, matIdQuad,
                                                                                          *gmesh);
            }


            TPZAutoPointer<TPZRefPattern> refp;
            char buf[] =
                    "6 4 "
                    "-98 QuadIntoAll "
                    "-1 -1   0 "//0
                    " 1 -1   0 "//1
                    " 1  1   0 "//2
                    "-1  1   0 "//3
                    " 0 -1   0 "//4
                    " 0  1   0 "//5
                    "3 4  0 1 2 3 " //master quad
                    "3 4  0 4 5 3 " //left half of the cube
                    "2 3  4 1 2 "//triang1
                    "2 3  4 2 5 "; //triang2

            std::istringstream str(buf);
            refp = new TPZRefPattern(str);
            refp->GenerateSideRefPatterns();
            gRefDBase.InsertRefPattern(refp);
            const int firstEdge = pztopology::TPZQuadrilateral::NumSides(0);
            const int nEdges = pztopology::TPZQuadrilateral::NumSides(1);
            quadEl->SetRefPattern(refp);
            TPZGeoElRefPattern<pzgeom::TPZArc3D> *arc = nullptr;
            for (int edge = firstEdge; edge < firstEdge + nEdges; edge++) {
                const int nNodes = pztopology::TPZQuadrilateral::NSideNodes(edge);
                nodesIdVec.Resize(nNodes + 1);
                for (int node = 0; node < nNodes; node++) {
                    nodesIdVec[node] = pztopology::TPZQuadrilateral::SideNodeLocId(edge, node);
//                std::cout<<nodesIdVec[node]<<"\t";
                }
                nodesIdVec[nNodes] = nodesIdArcsVec[edge - firstEdge];
//            std::cout<<nodesIdVec[nNodes]<<std::endl;
                arc = new TPZGeoElRefPattern<pzgeom::TPZArc3D>(nodesIdVec, matIdArc, *gmesh);
                auto arcRefp = refp->SideRefPattern(edge);
                arc->SetRefPattern(arcRefp);
            }
            gmesh->BuildConnectivity();

            {
                TPZVec<TPZGeoEl *> sons;
                const int nel = gmesh->NElements();
                quadEl->Divide(sons);
                for (int iel = 1; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->Element(iel);
                    auto geoElSide = geo->Neighbour(geo->NSides() - 1);
                    if (geoElSide.NSubElements() > 1) {
                        geo->Divide(sons);
                    }
                }
            }

            //Create GeoBlend elements
            {
                TPZGeoMesh *newgmesh = new TPZGeoMesh();
                for (int64_t i = 0; i < gmesh->NNodes(); i++) {
                    coord[0] = gmesh->NodeVec()[i].Coord(0);
                    coord[1] = gmesh->NodeVec()[i].Coord(1);
                    coord[2] = gmesh->NodeVec()[i].Coord(2);
                    const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
                    newgmesh->NodeVec()[newindex].Initialize(coord, *newgmesh);
                }
                //TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >
                const int nel = gmesh->NElements();
                newgmesh->ElementVec().Resize(nel);
                for (int iel = 0; iel < nel; iel++) {
                    newgmesh->ElementVec()[iel] = nullptr;
                    TPZGeoEl *geo = gmesh->Element(iel);
                    nodesIdVec.resize(geo->NNodes());
                    for (int i = 0; i < nodesIdVec.size(); i++) {
                        nodesIdVec[i] = geo->NodeIndex(i);
                    }
                    const int matId = geo->MaterialId();
                    const int elId = geo->Index();
//                geo->Ind
                    switch (geo->Type()) {
                        case ETriangle:
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle>>(elId, nodesIdVec, matId,
                                                                                                *newgmesh);
                            break;
                        case EQuadrilateral:
                            if (geo->HasSubElement()) continue;//the original cube should not be copied
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad>>(elId, nodesIdVec, matId,
                                                                                            *newgmesh);
                            break;
                        default:
                            geo->Clone(*newgmesh);
                            break; //the element should not be deleted
                    }
                    //TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
                }
                delete gmesh;
                gmesh = newgmesh;
            }


            gmesh->ResetConnectivities();
            gmesh->BuildConnectivity();
#ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "blendmesh2D_partial";
                const size_t strlen = meshFileName.length();
                meshFileName.append(".vtk");
                std::ofstream outVTK(meshFileName.c_str());
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());

                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
                gmesh->Print(outTXT);
                outTXT.close();
                outVTK.close();
            }
#endif
            #ifdef _AUTODIFF
            {
                TPZManVector<REAL,3> xiReal;
                TPZManVector<Fad<REAL>,3> xiFad;
                REAL weight = -1;//useless
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
                        uint64_t errors = 0;
                        auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                        xiReal.Resize(geo->Dimension(),0);
                        xiFad.Resize(geo->Dimension(),0);
                        for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                            bool hasAnErrorOccurred = false;
                            intRule->Point(iPt,xiReal,weight);
                            for(int x = 0; x < geo->Dimension(); x++){
                                xiFad[x] = Fad<REAL>(geo->Dimension(),x,xiReal[x]);
                            }
//                        Fad<REAL> func;
//                        func = 9*xiFad[0];
//                        for(int i =0; i < func.size(); i++){
//                            std::cout<<"dx["<<i<<"]:\t"<<func.dx(i)<<std::endl;
//                        }
                            TPZManVector<Fad<REAL>,3> xFad(3);
                            geo->X(xiFad, xFad);
                            TPZFNMatrix<9,REAL> gradXreal(3,3,0);
                            geo->GradX(xiReal, gradXreal);
                            std::ostringstream xFadM, gradXM;
                            xFadM<<"xFad:"<<std::endl;
                            gradXM<<"grad x:"<<std::endl;
                            const auto VAL_WIDTH = 10;
                            for(int i = 0; i < gradXreal.Rows(); i++){
                                for(int j = 0; j < gradXreal.Cols(); j++){
                                    xFadM<<std::setw(VAL_WIDTH) << std::right<<xFad[i].dx(j)<<"\t";
                                    gradXM<<std::setw(VAL_WIDTH) << std::right<<gradXreal(i,j)<<"\t";
                                    const REAL diff = (xFad[i].dx(j)-gradXreal(i,j))*(xFad[i].dx(j)-gradXreal(i,j));
                                    if(diff > tol){
                                        if(!hasAnErrorOccurred) errors++;
                                        hasAnErrorOccurred = true;
                                    }
                                }
                                xFadM<<std::endl;
                                gradXM<<std::endl;
                            }
                            BOOST_CHECK(!hasAnErrorOccurred);
                            if(hasAnErrorOccurred){
                                std::cout<<std::flush;
                                std::cout<<xFadM.str()<<std::endl;
                                std::cout<<gradXM.str()<<std::endl;
                            }
                        }
                        if(errors > 0){
                            std::cout<<"============================"<<std::endl;
                            std::cout<<"Element: "<<geo->Id()<<"\tType: "<<MElementType_Name(geo->Type());
                            std::cout<<"\tIs blend? : "<<geo->IsGeoBlendEl()<<std::endl;
                            std::cout<<"\tNumber of points: "<<intRule->NPoints()<<"\tErrors: "<<errors<<std::endl;
                        }
                    }
                }
            }
            #endif


            {
                TPZVec<TPZGeoEl *> sons;
                std::vector<std::string> loading = {"-","/","|","\\"};
                for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                    #ifdef NOISY_BLEND
                    std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
                    #endif
                    const int nel = gmesh->NElements();
                    for (int iel = 0; iel < nel; iel++) {
                        #ifdef NOISY_BLEND
                        std::cout<<"\b"<<loading[iel%4]<<std::flush;
                        #endif
                        TPZGeoEl *geo = gmesh->ElementVec()[iel];
                        if (geo && !geo->HasSubElement()) {
                            geo->Divide(sons);
                        }
                    }
                    #ifdef NOISY_BLEND
                    std::cout<<"\b"<<std::endl;
                    #endif
                }
            }

            #ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "blendmesh2D";
                const size_t strlen = meshFileName.length();
                meshFileName.append(".vtk");
                std::ofstream outVTK(meshFileName.c_str());
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());

                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
                gmesh->Print(outTXT);
                outTXT.close();
                outVTK.close();
            }
            #endif
            delete gmesh;
        }

        void CreateGeoMesh3D(int nDiv){
            std::cout << __PRETTY_FUNCTION__ << std::endl;
            TPZGeoMesh *gmesh = new TPZGeoMesh();

            const int64_t nNodesHexa = 8;
            const int64_t nNodes= nNodesHexa;//+nNodesArc;

            const REAL radius = std::sqrt(3);
            const REAL cubeSide= 2;
            TPZVec<REAL> coord(3, 0.);
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///                                                                                                         ////
            ////the face containing the nodes 0,1,2 and 3 will be joined with a half-sphere                             ////
            ///                                                                                                         ////
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            for (int64_t i = 0; i < nNodesHexa; i++) {
                int iAux = i%4;
                bool xFactor = (bool) !( (iAux % 2) == (iAux / 2) );
                bool yFactor = (bool) ( i / 4 );
                bool zFactor = (bool) ( ( iAux /2 ) % 2 );
                coord[0] = 0.5 * cubeSide * (2 * (int)xFactor - 1) ;
                coord[1] = 0.5 * cubeSide  * (1 - 2 * (int)yFactor);
                coord[2] = 0.5 * cubeSide * (2 * (int)zFactor - 1) ;
                const int64_t newindex = gmesh->NodeVec().AllocateNewElement();
                gmesh->NodeVec()[newindex].Initialize(coord,*gmesh);
//            std::cout<<"x : "<<xFactor<<"\ty : "<<yFactor<<"\tz : "<<zFactor<<std::endl;
//            std::cout<<"x : "<<coord[0]<<"\ty : "<<coord[1]<<"\tz : "<<coord[2]<<std::endl;
            }


            int matIdVol = 1, matIdSphere = 2;
            TPZVec<int64_t> nodesIdVec(1);
            int64_t elId = 0;

            nodesIdVec.Resize(8);

            nodesIdVec[0] = 0;
            nodesIdVec[1] = 1;
            nodesIdVec[2] = 2;
            nodesIdVec[3] = 3;
            nodesIdVec[4] = 4;
            nodesIdVec[5] = 5;
            nodesIdVec[6] = 6;
            nodesIdVec[7] = 7;
            TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > *hexaEl =
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >(nodesIdVec, matIdVol, *gmesh);
            TPZAutoPointer<TPZRefPattern> refp;
            char buf[] =
                    "12 5 "
                    "-99 CubeIntoAll "
                    "-1 -1  -1 "
                    " 1 -1  -1 "
                    " 1  1  -1 "
                    "-1  1  -1 "
                    "-1 -1   0 "
                    " 1 -1   0 "
                    " 1  1   0 "
                    "-1  1   0 "
                    "-1 -1   1 "
                    " 1 -1   1 "
                    " 1  1   1 "
                    "-1  1   1 "
                    "7 8  0 1 2 3 8 9 10 11 " //master cube
                    "7 8  4 5 6 7 8 9 10 11 " //upper half of the cube
                    "6 6  0 5 4 3 6 7 "//prism
                    "5 5  0 3 6 5 2 " //pyramid
                    "4 4  1 2 5 0 "; //tetra

            std::istringstream str(buf);
            refp = new TPZRefPattern(str);
            refp->GenerateSideRefPatterns();
            gRefDBase.InsertRefPattern(refp);
            const int firstFace = pztopology::TPZCube::NumSides(0) + pztopology::TPZCube::NumSides(1);
            const int nFaces = pztopology::TPZCube::NumSides(2);
            hexaEl->SetRefPattern(refp);
            TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>> *sphere = nullptr;
            for(int face = firstFace; face < firstFace + nFaces; face++){
                const int nNodes = pztopology::TPZCube::NSideNodes(face);
                nodesIdVec.Resize(nNodes);
                for (int node = 0; node < nNodes; node++){
                    nodesIdVec[node] = pztopology::TPZCube::SideNodeLocId(face,node);
                }
                sphere = new TPZGeoElRefPattern<pzgeom::TPZQuadSphere<>>(nodesIdVec, matIdSphere, *gmesh);
                coord[0] = coord[1] = coord[2] = 0;
                sphere->Geom().SetData(radius, coord);
                auto sphereRefp = refp->SideRefPattern(face);
                sphere->SetRefPattern(sphereRefp);
            }
            gmesh->BuildConnectivity();

            {
                TPZVec<TPZGeoEl *> sons;
                const int nel = gmesh->NElements();
                hexaEl->Divide(sons);
                for (int iel = 1; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->Element(iel);
                    auto geoElSide = geo->Neighbour(geo->NSides()-1);
                    if(geoElSide.NSubElements() > 1)
                    {
                        geo->Divide(sons);
                    }
                }
            }

            //Create GeoBlend elements
            {
                TPZGeoMesh *newgmesh = new TPZGeoMesh();
                for (int64_t i = 0; i < gmesh->NNodes(); i++) {
                    coord[0] = gmesh->NodeVec()[i].Coord(0);
                    coord[1] = gmesh->NodeVec()[i].Coord(1);
                    coord[2] = gmesh->NodeVec()[i].Coord(2);
                    const int64_t newindex = newgmesh->NodeVec().AllocateNewElement();
                    newgmesh->NodeVec()[newindex].Initialize(coord,*newgmesh);
                }
                //TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> >
                const int nel = gmesh->NElements();
                newgmesh->ElementVec().Resize(nel);
                for (int iel = 0; iel < nel; iel++) {
                    newgmesh->ElementVec()[iel] = nullptr;
                    TPZGeoEl *geo = gmesh->Element(iel);
                    nodesIdVec.resize(geo->NNodes());
                    for(int i = 0; i < nodesIdVec.size(); i++){
                        nodesIdVec[i] = geo->NodeIndex(i);
                    }
                    const int matId = geo->MaterialId();
                    const int elId = geo->Index();
//                geo->Ind
                    switch(geo->Type()){
                        case ETetraedro:
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTetrahedra>>(elId,nodesIdVec,matId,*newgmesh);
                            break;
                        case EPiramide:
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPyramid>>(elId,nodesIdVec,matId,*newgmesh);
                            break;
                        case EPrisma:
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism>>(elId,nodesIdVec,matId,*newgmesh);
                            break;
                        case ECube:
                            if(geo->HasSubElement()) continue;//the original cube should not be copied
                            new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube>>(elId,nodesIdVec,matId,*newgmesh);
                            break;
                        default:
                            geo->Clone(*newgmesh);
                            break; //the element should not be deleted
                    }
                    //TPZGeoElRefPattern(int64_t id,TPZVec<int64_t> &nodeindexes,int matind,TPZGeoMesh &mesh);
                }
                delete gmesh;
                gmesh = newgmesh;
            }
            gmesh->ResetConnectivities();
            gmesh->BuildConnectivity();
#ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "blendmesh3D_partial";
                const size_t strlen = meshFileName.length();
                meshFileName.append(".vtk");
                std::ofstream outVTK(meshFileName.c_str());
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());

                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
                gmesh->Print(outTXT);
                outTXT.close();
                outVTK.close();
            }
#endif
#ifdef _AUTODIFF
            {
                TPZManVector<REAL,3> xiReal;
                TPZManVector<Fad<REAL>,3> xiFad;
                REAL weight = -1;//useless
                const int nel = gmesh->NElements();
                for (int iel = 0; iel < nel; iel++) {
                    TPZGeoEl *geo = gmesh->ElementVec()[iel];
                    if (geo && !geo->HasSubElement()) {
                        uint64_t errors = 0;
                        auto intRule = geo->CreateSideIntegrationRule(geo->NSides()-1, pOrder);
                        xiReal.Resize(geo->Dimension(),0);
                        xiFad.Resize(geo->Dimension(),0);
                        for(int iPt = 0; iPt < intRule->NPoints(); iPt++){
                            bool hasAnErrorOccurred = false;
                            intRule->Point(iPt,xiReal,weight);
                            for(int x = 0; x < geo->Dimension(); x++){
                                xiFad[x] = Fad<REAL>(geo->Dimension(),x,xiReal[x]);
                            }

                            TPZManVector<Fad<REAL>,3> xFad(3);
                            geo->X(xiFad, xFad);
                            TPZFNMatrix<9,REAL> gradXreal(3,3,0);
                            geo->GradX(xiReal, gradXreal);
                            std::ostringstream xFadM, gradXM;
                            xFadM<<"xFad:"<<std::endl;
                            gradXM<<"grad x:"<<std::endl;
                            const auto VAL_WIDTH = 10;
                            for(int i = 0; i < gradXreal.Rows(); i++){
                                for(int j = 0; j < gradXreal.Cols(); j++){
                                    xFadM<<std::setw(VAL_WIDTH) << std::right<<xFad[i].dx(j)<<"\t";
                                    gradXM<<std::setw(VAL_WIDTH) << std::right<<gradXreal(i,j)<<"\t";
                                    const REAL diff = (xFad[i].dx(j)-gradXreal(i,j))*(xFad[i].dx(j)-gradXreal(i,j));
                                    if(diff > tol){
                                        if(!hasAnErrorOccurred) errors++;
                                        hasAnErrorOccurred = true;
                                    }
                                }
                                xFadM<<std::endl;
                                gradXM<<std::endl;
                            }
                            BOOST_CHECK(!hasAnErrorOccurred);
                            if(hasAnErrorOccurred){
                                std::cout<<std::flush;
                                std::cout<<xFadM.str()<<std::endl;
                                std::cout<<gradXM.str()<<std::endl;
                            }
                        }
                        if(errors > 0){
                            std::cout<<"============================"<<std::endl;
                            std::cout<<"Element: "<<geo->Id()<<"\tType: "<<MElementType_Name(geo->Type());
                            std::cout<<"\tIs blend? : "<<geo->IsGeoBlendEl()<<std::endl;
                            std::cout<<"\tNumber of points: "<<intRule->NPoints()<<"\tErrors: "<<errors<<std::endl;
                        }
                    }
                }
            }
#endif
            {
                TPZVec<TPZGeoEl *> sons;
                std::vector<std::string> loading = {"-","/","|","\\"};
                for (int iDiv = 0; iDiv < nDiv; iDiv++) {
                    #ifdef NOISY_BLEND
                    std::cout<<"Performing "<<iDiv+1<<" ref step out of " << nDiv<<std::endl;
                    #endif
                    const int nel = gmesh->NElements();
                    for (int iel = 0; iel < nel; iel++) {
                        #ifdef NOISY_BLEND
                        std::cout<<"\b"<<loading[iel%4]<<std::flush;
                        #endif
                        TPZGeoEl *geo = gmesh->ElementVec()[iel];
                        if (geo && !geo->HasSubElement()) {
                            geo->Divide(sons);
                        }
                    }
                    #ifdef NOISY_BLEND
                    std::cout<<"\b"<<std::endl;
                    #endif
                }
            }

            #ifdef NOISYVTK_BLEND
            {
                std::string meshFileName = "blendmesh3D";
                const size_t strlen = meshFileName.length();
                meshFileName.append(".vtk");
                std::ofstream outVTK(meshFileName.c_str());
                meshFileName.replace(strlen, 4, ".txt");
                std::ofstream outTXT(meshFileName.c_str());

                TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outVTK, true);
                gmesh->Print(outTXT);
                outTXT.close();
                outVTK.close();
            }
            #endif
            delete gmesh;
        }



    }


BOOST_AUTO_TEST_SUITE_END()


#endif

