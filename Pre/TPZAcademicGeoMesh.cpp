//
//  TPZAcademicGeoMesh.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/29/16.
//
//

#include "TPZAcademicGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZGenGrid2D.h"
#include "TPZExtendGridDimension.h"
#include "TPZVTKGeoMesh.h"

#include "pzcheckmesh.h"

#include "TPZMaterial.h"
#include "TPZBndCond.h"

#include "TPZLinearAnalysis.h"
#include "pzstepsolver.h"

#include "TPZRefPatternTools.h"


#include <stdio.h>

#include "pzlog.h"

#include "pzgeoelbc.h"

/** Initialiazing file for Log4CXX for this project */
#ifdef PZ_LOG
static TPZLogger logger("pz.academicmesh");
#endif



TPZAcademicGeoMesh::TPZAcademicGeoMesh() : fMeshType(EHexa), fBCNumbers(6,-1), fShouldDeform(false), fNumberElements(1)
{
    REAL coord[8][3] = {
        //        {0,0,0},
        //        {1,0,0},
        //        {1,1,0},
        //        {0,1,0},
        //        {0,0,1},
        //        {1,0,1},
        //        {1,1,1},
        //        {0,1,1}
        {-0.5,-0.5,-0.5},
        {2,-1,-1},
        {1.1,1.1,-0.1},
        {-1,2,-1},
        {-1,-1,2},
        {1.2,-0.2,1.2},
        {2,2,2},
        {-0.5,1.5,1.5}
    };
    fDeformed.NodeVec().Resize(8);
    TPZManVector<int64_t,8> indices(8);
    for (int i=0; i<8; i++) {
        indices[i] = i;
        for (int c=0; c<3; c++) {
            fDeformed.NodeVec()[i].SetCoord(c, coord[i][c]);
        }
    }
    int64_t index;
    fDeformed.CreateGeoElement(ECube, indices, 1, index);
}

/// Deform the geometric mesh according to the coordinates of fDeformed
void TPZAcademicGeoMesh::DeformGMesh(TPZGeoMesh &gmesh)
{
    int64_t nnodes = gmesh.NodeVec().NElements();
    TPZManVector<REAL,3> xbefore(3),xafter(3);
    for (int64_t nod=0; nod<nnodes; nod++) {
        gmesh.NodeVec()[nod].GetCoordinates(xbefore);
        for (int i=0; i<3; i++) {
            xbefore[i] = 2.*xbefore[i]-1.;
        }
        fDeformed.ElementVec()[0]->X(xbefore, xafter);
        gmesh.NodeVec()[nod].SetCoord(xafter);
    }
}

static int pyramid[2][5]=
{
    {0,1,2,3,4},
    {4,5,6,7,2}
};
static int tetraedra[2][4]=
{
    {1,2,5,4},
    {4,7,3,2}
};
static int tetraedra_2[6][4]=
{
    {1,2,5,4},
    {4,7,3,2},
    {0,1,2,4},
    {0,2,3,4},
    {4,5,6,2},
    {4,6,7,2}
};

static int prism[2][6] =
{
    {0,1,2,4,5,6},
    {0,2,3,4,6,7}
};


void TPZAcademicGeoMesh::GenerateNodes(TPZGeoMesh *gmesh)
{
    int64_t nelem = fNumberElements;
    gmesh->NodeVec().Resize((nelem+1)*(nelem+1)*(nelem+1));
    for (int64_t i=0; i<=nelem; i++) {
        for (int64_t j=0; j<=nelem; j++) {
            for (int64_t k=0; k<=nelem; k++) {
                TPZManVector<REAL,3> x(3);
                x[0] = k*1./nelem;
                x[1] = j*1./nelem;
                x[2] = i*1./nelem;
                gmesh->NodeVec()[i*(nelem+1)*(nelem+1)+j*(nelem+1)+k].Initialize(x, *gmesh);
            }
        }
    }
}

TPZGeoMesh *TPZAcademicGeoMesh::PyramidalAndTetrahedralMesh()
{
    int64_t nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    const int dim = 3;
    gmesh->SetDimension(dim);
    GenerateNodes(gmesh);
  
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                if(logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Pyramid and tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<2; el++)
                {
                    TPZManVector<int64_t,5> elnodes(5);
                    for (int il=0; il<5; il++) {
                        elnodes[il] = nodes[pyramid[el][il]];
                    }
                    int64_t index;
                    gmesh->CreateGeoElement(EPiramide, elnodes, MaterialId, index);
                    elnodes.resize(4);
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
//    AddBoundaryElements(gmesh);
    AddBoundaryElementsByCoord(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}

TPZGeoMesh *TPZAcademicGeoMesh::RedBlackPyramidalAndHexagonalMesh(){
    
    int64_t n = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh);
    gmesh->SetDimension(3);
    
    int64_t s = n + 1;
    REAL dir = 1.0;
    REAL dx = 0.5/n;
    int64_t n_nodes = gmesh->NNodes();
    TPZStack<int> perm;
    perm.push_back(0);
    perm.push_back(1);
    perm.push_back(3);
    perm.push_back(2);
    perm.push_back(4);
    perm.push_back(5);
    perm.push_back(7);
    perm.push_back(6);
    
    TPZManVector<int,8> indexes(8);
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                
                int stride =  j*s + k*s*s;
                indexes[0] = i + stride;
                indexes[1] = i + 1 + stride;
                indexes[2] = i + s + stride;
                indexes[3] = i + 1 + s + stride;
                for (int l = 0; l < 4; l++) { // Vertical stride
                    indexes[l+4] = indexes[l]+s*s;
                }
                
                bool is_even_Q = (i+j+k)%2==0;
                if (is_even_Q) {
                    
                    TPZManVector<int64_t,8> nodes(8);
                    for (int inode=0; inode<8; inode++) {
                        nodes[inode] = gmesh->NodeVec()[indexes[perm[inode]]].Id();
                    }
                    
                    // Adding center point
                    TPZManVector<REAL,3> x(3);
                    gmesh->NodeVec()[nodes[0]].GetCoordinates(x);
                    for (int i = 0; i < 3; i++) {
                        x[i]+=dx*dir;
                    }
                    
                    gmesh->NodeVec().Resize(n_nodes+1);
                    gmesh->NodeVec()[n_nodes].Initialize(x, *gmesh);
                    int center_id = gmesh->NodeVec()[n_nodes].Id();
                    n_nodes++;
                    
                    // inserting Pyramids
                    int64_t index;
                    TPZManVector<int64_t,5> el_nodes(5);
                    
                    el_nodes[0] = nodes[0];
                    el_nodes[1] = nodes[1];
                    el_nodes[2] = nodes[2];
                    el_nodes[3] = nodes[3];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    el_nodes[0] = nodes[0];
                    el_nodes[1] = nodes[3];
                    el_nodes[2] = nodes[7];
                    el_nodes[3] = nodes[4];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    el_nodes[0] = nodes[0];
                    el_nodes[1] = nodes[1];
                    el_nodes[2] = nodes[5];
                    el_nodes[3] = nodes[4];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    el_nodes[0] = nodes[1];
                    el_nodes[1] = nodes[2];
                    el_nodes[2] = nodes[6];
                    el_nodes[3] = nodes[5];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    el_nodes[0] = nodes[3];
                    el_nodes[1] = nodes[2];
                    el_nodes[2] = nodes[6];
                    el_nodes[3] = nodes[7];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    el_nodes[0] = nodes[4];
                    el_nodes[1] = nodes[5];
                    el_nodes[2] = nodes[6];
                    el_nodes[3] = nodes[7];
                    el_nodes[4] = center_id;
                    gmesh->CreateGeoElement(EPiramide, el_nodes, MaterialId, index);
                    
                    
                }else{
                    TPZManVector<int64_t,8> el_nodes(8);
                    int64_t index;
                    for (int inode=0; inode<8; inode++) {
                        el_nodes[inode] = gmesh->NodeVec()[indexes[perm[inode]]].Id();
                    }
                    gmesh->CreateGeoElement(ECube, el_nodes, MaterialId, index);
                }
                
            }
        }
    }
    
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    
    return gmesh;
}

TPZGeoMesh *TPZAcademicGeoMesh::TetrahedralMesh()
{
    int64_t nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GenerateNodes(gmesh);
    gmesh->SetDimension(3);
  
    for (int64_t i=0; i<nelem; i++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t k=0; k<nelem; k++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                if(logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Tetrahedral nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                for (int el=0; el<6; el++)
                {
                    TPZManVector<int64_t,4> elnodes(4);
                    int64_t index;
                    for (int il=0; il<4; il++) {
                        elnodes[il] = nodes[tetraedra_2[el][il]];
                    }
                    gmesh->CreateGeoElement(ETetraedro, elnodes, MaterialId, index);
                }
            }
        }
    }
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}


TPZGeoMesh *TPZAcademicGeoMesh::HexahedralMesh()
{
    int64_t nelem = fNumberElements;
    int MaterialId = fMaterialId;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(3);
    GenerateNodes(gmesh);
    
    for (int64_t k=0; k<nelem; k++) {
        for (int64_t j=0; j<nelem; j++) {
            for (int64_t i=0; i<nelem; i++) {
                TPZManVector<int64_t,8> nodes(8,0);
                nodes[0] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[1] = k*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[2] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[3] = k*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
                nodes[4] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i;
                nodes[5] = (k+1)*(nelem+1)*(nelem+1)+j*(nelem+1)+i+1;
                nodes[6] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i+1;
                nodes[7] = (k+1)*(nelem+1)*(nelem+1)+(j+1)*(nelem+1)+i;
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Cube nodes " << nodes;
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                int64_t index;
                gmesh->CreateGeoElement(ECube, nodes, MaterialId, index);
            }
        }
    }
    gmesh->BuildConnectivity();
    AddBoundaryElements(gmesh);
    if (fShouldDeform) {
        DeformGMesh(*gmesh);
    }
    return gmesh;
}

/// verify if the faces without neighbour should be orthogonal to the main planes
void TPZAcademicGeoMesh::CheckConsistency(TPZGeoMesh *mesh)
{
    int64_t nel = mesh->NElements();
    for(int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = mesh->ElementVec()[el];
        int nsides = gel->NSides();
        for(int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.Dimension() != 2) {
                continue;
            }
            if(gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            gelside.CenterPoint(xi);
            TPZFNMatrix<6,REAL> axes(2,3);
            TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
            REAL detjac;
            gelside.Jacobian(xi, jac, axes, detjac, jacinv);
            TPZManVector<REAL,3> x(3,0.);
            gelside.X(xi, x);
            TPZManVector<REAL,3> normal(3);
            normal[0] = fabs(axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1));
            normal[1] = fabs(-axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0));
            normal[2] = fabs(axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0));
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numtol = 0;
            for(int i=0; i<3; i++) {
                if(xmin > x[i]) xmin = x[i];
                if(xmax < x[i]) {
                    xmax = x[i];
                }
                if(normal[i] > tol) {
                    numtol++;
                }
            }
            if(numtol != 1) {
                DebugStop();
            }
            if(xmin > tol && xmax < 1.-tol) {
                DebugStop();
            }
        }
    }
}


int TPZAcademicGeoMesh::AddBoundaryElements(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for(int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        int nsides = gel->NSides();
        for(int is=0; is<nsides; is++) {
            TPZGeoElSide gelside(gel,is);
            if(gelside.Dimension() != 2) {
                continue;
            }
            if(gelside.Neighbour() != gelside) {
                continue;
            }
            TPZManVector<REAL,2> xi(2,0.);
            TPZManVector<REAL,3> normal(3);
            gelside.Normal(xi, normal);
            REAL tol = 1.e-6;
            REAL xmin = 1., xmax = 0.;
            int numfound = 0;
            int dir = -5;
            for (int i=0; i<3; i++) {
                if (fabs(normal[i] - 1.) < tol) {
                    dir = 1+i;
                    numfound++;
                }
                if (fabs(normal[i] + 1.) < tol) {
                    dir = -1-i;
                    numfound++;
                }
            }
            if (dir == -5) {
                DebugStop();
            }
            if(numfound != 1) {
                DebugStop();
            }
            switch (dir) {
                case -3:
                    TPZGeoElBC(gelside,fBCNumbers[0]);
                    break;
                case 1:
                    TPZGeoElBC(gelside,fBCNumbers[1]);
                    break;
                case 2:
                    TPZGeoElBC(gelside,fBCNumbers[2]);
                    break;
                case -1:
                    TPZGeoElBC(gelside,fBCNumbers[3]);
                    break;
                case -2:
                    TPZGeoElBC(gelside,fBCNumbers[4]);
                    break;
                case 3:
                    TPZGeoElBC(gelside,fBCNumbers[5]);
                    break;
                default:
                    DebugStop();
                    break;
            }
        }
    }
	return 0;
}

int TPZAcademicGeoMesh::AddBoundaryElementsByCoord(TPZGeoMesh *gmesh)
{
  int64_t nel = gmesh->NElements();
  for(int64_t el=0; el<nel; el++) {
    TPZGeoEl *gel = gmesh->ElementVec()[el];
    int nsides = gel->NSides();
    for(int is=0; is<nsides; is++) {
      TPZGeoElSide gelside(gel,is);
      if(gelside.Dimension() != 2) {
        continue;
      }
      if(gelside.Neighbour() != gelside) {
        continue;
      }
      TPZManVector<REAL,2> xi(2,0.);
      TPZFNMatrix<6,REAL> axes(2,3);
      TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
      REAL detjac;
      gelside.Jacobian(xi, jac, axes, detjac, jacinv);
      TPZManVector<REAL,3> x(3,0.);

      REAL tol = 1.e-6;
      REAL coordmin = 0., coordmax = 0.;
      const int nnodes = gelside.NSideNodes();
      int nxbot = 0, nybot = 0, nzbot = 0, nxtop = 0, nytop = 0, nztop = 0;
      for (int i = 0; i < nnodes; i++) {
        const int nodeindex = gelside.SideNodeIndex(i);
        TPZManVector<REAL,3> coord(3,0);
        gmesh->NodeVec()[nodeindex].GetCoordinates(coord);
        const REAL diffxbot = fabs(coord[0]);
        const REAL diffybot = fabs(coord[1]);
        const REAL diffzbot = fabs(coord[2]);
        const REAL diffxtop = fabs(coord[0]-1.);
        const REAL diffytop = fabs(coord[1]-1.);
        const REAL diffztop = fabs(coord[2]-1.);
        if(diffxbot == 0){
          nxbot++;
        }
        if(diffybot == 0){
          nybot++;
        }
        if(diffzbot == 0){
          nzbot++;
        }
        if(diffxtop == 0){
          nxtop++;
        }
        if(diffytop == 0){
          nytop++;
        }
        if(diffztop == 0){
          nztop++;
        }
      }
      
      int numbcsfound = 0;
      if (nxbot >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[0]);
        numbcsfound++;
      }
      if (nybot >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[1]);
        numbcsfound++;
      }
      if (nzbot >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[2]);
        numbcsfound++;
      }
      if (nxtop >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[3]);
        numbcsfound++;
      }
      if (nytop >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[4]);
        numbcsfound++;
      }
      if (nztop >= 3) {
        TPZGeoElBC(gelside,fBCNumbers[5]);
        numbcsfound++;
      }

      if (numbcsfound != 1) {
        DebugStop();
      }
    }
  }
  return 0;
}




