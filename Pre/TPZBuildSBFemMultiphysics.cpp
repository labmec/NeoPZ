//
//  TPZBuildSBFemMultiphysics.cpp
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#include "TPZBuildSBFemMultiphysics.h"

#include "TPZCompElHDivSBFem.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZBndCondBase.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzgeoelbc.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "TPZVTKGeoMesh.h"

#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"

#include "TPZSBFemMultiphysicsElGroup.h"
#include "TPZSBFemVolumeL2.h"
#include "TPZSBFemVolumeHdiv.h"
#include "TPZSBFemVolumeMultiphysics.h"

#include "pzcondensedcompel.h"
#include "pzgeoelbc.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzbuildsbfem"));
static LoggerPtr loggersbfemhdivgeom(Logger::getLogger("buildsbfemhdivgeom"));
static LoggerPtr loggercmesh(Logger::getLogger("cmeshsbfemhdiv"));
#endif

// This function will create both TPZSBFemVolumeMultiphysics (and its atomic elements) and TPZSBFemMultiphysicsElGroup.
// Sequence of the method:
// 1. Define geometry with the external dim-1 elements and collapsed element;
// 2. Update atomic meshes (cmeshflux and cmeshpressure);
// 3. Update Multiphysics mesh;
// 4. Create interface elements;
// 5. Create TPZSBFemVolumeMultiphysics objects based on the multiphysics elements;
// 6. Update the multiphysics mesh to contemplate only the SBFem elements (ignore FE els);
// 7. Define SBFemMultiphysicsElGroups as a group of TPZSBFemVolumeMultiphysics;
// 8. Condense the DOFs of the SBFemMultiphysicsElGroups.

/// create geometric volumetric elements
void TPZBuildSBFemMultiphysics::CreateVolumetricElementsFromSkeleton(TPZCompMesh &cmesh, std::set<int> & matidstarget, std::set<int> & matids1d)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    cmesh.LoadReferences();
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmesh.FindMaterial(mat)) {
            matidstarget.insert(it->second);
        }
    }
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement() ) {
            continue;
        }
        if (!gel->Reference()) {
            continue;
        }
        if (fElementPartition[el] == -1) {
            continue;
        }
        if (gel->Dimension() > dim - 1) {
            DebugStop();
        }
        if (matids1d.find(gel->MaterialId()) == matids1d.end()) {
            continue;
        }
        int nsides = gel->NSides();
        int is = nsides - 1;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gelside(gel,is);
        int nnodes = gelside.NSideNodes();
        TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
        int matid = fMatIdTranslation[gel->MaterialId()];
        int64_t index;
        for (int in=0; in<nnodes; in++) {
            Nodes[in] = gelside.SideNodeIndex(in);
        }
        int elpartition = fElementPartition[el];
        for (int in=nnodes; in < 2*nnodes; in++) {
            Nodes[in] = fPartitionCenterNode[elpartition];
        }
        if (gelside.IsLinearMapping())
        {
            switch(nnodes)
            {
                case 1:
                    gmesh->CreateGeoElement(EOned, Nodes, matid, index);
                    break;
                case 2:
                    gmesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                    break;
                case 4:
                    gmesh->CreateGeoElement(ECube, Nodes, matid, index);
                    break;
                case 3:
                    gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                    break;
                default:
                    std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                    DebugStop();
            }
            
        }
        else
        {
            int64_t elementid = gmesh->NElements()+1;
            switch(nnodes)
            {
                case 2:
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);
                    break;
                case 4:
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmesh,index);
                    break;
                case 3:
                    gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                    break;
                default:
                    std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                    DebugStop();
            }
        }
        if (index >= fElementPartition.size()) {
            fElementPartition.resize(index+1);
        }
        fElementPartition[index] = elpartition;
        
        
    }
    gmesh->BuildConnectivity();
}

// /// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFemMultiphysics::BuildComputationalMeshFromSkeleton(TPZMultiphysicsCompMesh &cmeshm)
{
    TPZManVector<TPZCompMesh*, 2> cmeshvec = cmeshm.MeshVector();

    // The flux cmesh must have:
    // 1. Material ID defined in the map EGroup -> EMatVol
    auto cmeshflux = cmeshvec[0];    

    // The pressure cmesh must have:
    // 1. Material ID defining the map EGroup -> EMatVol.
    // 2. The Material ID related to the Skeleton - a Neumann BC - external average pressure.
    // 3. The boundary conditions.
    // All these materials must be defined in the cmeshm too.
    auto cmeshpressure = cmeshvec[1];

    // ********** CREATING THE GEOMETRY
    // Skeleton Elements + Collapsed Elements + External elements:

    // Before calling this function the Skeleton elements has been already created
    // So I just need to create the collapsed and external elements

    // Creating dim-1 elements CompEls - Skeleton
    int dim = cmeshm.Dimension();
    fGMesh->SetName("gmesh with collapsed els");
    cmeshpressure->SetName("cmesh pressure");

    set<int> matids1d;
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->Dimension() == fGMesh->Dimension()-1)
        {
            matids1d.insert(gel->MaterialId());
        }
    }
    cmeshpressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure->AutoBuild(matids1d);

    // Creating volumetric collapsed elements
    set<int> matidstarget; // matid of the collapsed element (output parameter)
    // The updated mesh will be fGMesh
    // But the comp mesh used must be cmeshpressure because the connectivity of the multiphysics mesh
    // hasn't been updated yet.
    // CreateCollapsedGeoEls(*cmeshpressure, matidstarget, matids1d);
    CreateVolumetricElementsFromSkeleton(*cmeshpressure, matidstarget, matids1d);

    // Creating dim-1 elements: External flux and external pressure.
    CreateExternalElements(fGMesh, matidstarget, matids1d);

    {
        std::ofstream gout("gmeshbuild.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(fGMesh, gout, true);
    }

    // ********** CREATING THE COMPUTATIONAL MESH
    // At this point the code created all geometric elements needed for ther SBFEM simulation and it's stored in fGMesh
    // Here it's created the SBFemVolumeHdiv elements -> CompElHdivElement -> Adjust connectivities
    // But firstly the atomic meshes must be properly created
    CreateCompElPressure(*cmeshpressure, matids1d);
    // Adding the volumetric elements
    CreateSBFemVolumePressure(*cmeshpressure, matids1d, matidstarget);
    cmeshpressure->LoadReferences();

    CreateCompElFlux(*cmeshflux, matidstarget, matids1d);
    // Adding the volumetric elements
    CreateSBFemVolumeFlux(*cmeshflux, matids1d, matidstarget);
    cmeshflux->LoadReferences();

    // Then, creating the multiphysics mesh
    // Next method will setup the multiphysics mesh
    // Deleting the Hdiv elements existing in the mesh.
    // The elements must be specified as TPZCompElHDivSBFem elements
    CreateSBFemMultiphysicsMesh(cmeshm,matidstarget);
    cmeshm.ApproxSpace().SetAllCreateFunctionsSBFemMultiphysics(dim);
    cmeshm.SetName("multiphysicssbfem");
    TPZManVector<MSpaceConfig> active(2,ETestTrial);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();
    {
        ofstream sout("cmeshmultiphysics.txt");
        cmeshm.Print(sout);
    }
    
    // Adding the interface elements
    AddInterfaceElements(cmeshm, matids1d);

    // Creating the SBFEMVolumeHdiv elements
    CreateSBFemMultiphysicsVol(cmeshm, matids1d, matidstarget);

    // Based on the multiphysics mesh, create SBFemElementGroups
    // The SBFemMultiphysicsElGroups = Group of collapsed flux elements + pressure compels
    // Then, condense the flux into the pressures.
    CreateSBFemMultiphysicsElGroups(cmeshm, matidstarget);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();
    
    GroupandCondense(cmeshm);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        ofstream mout("cmeshmultiphysicscondensed.txt");
        cmeshm.Print(mout);
    }
#endif

    auto nel = cmeshm.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshm.Element(el);
        if (!cel)
        {
            continue;
        }
        auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup * >(cel);
        if (sbfemgr)
        {
            continue;
        }
        if (cel->Reference())
        {
            auto matid = cel->Reference()->MaterialId();
            auto it = matids1d.find(matid);
            if (it != matids1d.end())
            {
                continue;
            }
        }
        cmeshm.ElementVec()[el] = 0;
    }
}

void TPZBuildSBFemMultiphysics::BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh & cmeshm)
{
    TPZManVector<TPZCompMesh*, 2> cmeshvec = cmeshm.MeshVector();

    // The flux cmesh must have:
    // 1. Material ID defined in the map EGroup -> EMatVol
    auto cmeshflux = cmeshvec[0];    

    // The pressure cmesh must have:
    // 1. Material ID defining the map EGroup -> EMatVol.
    // 2. The Material ID related to the Skeleton - a Neumann BC - external average pressure.
    // 3. The boundary conditions.
    // All these materials must be defined in the cmeshm too.
    auto cmeshpressure = cmeshvec[1];

    // ********** CREATING THE GEOMETRY
    // Skeleton Elements + Collapsed Elements + External elements:

    // Before calling this function the Skeleton elements has been already created
    // So I just need to create the collapsed and external elements

    // Creating dim-1 elements CompEls - Skeleton
    int dim = cmeshm.Dimension();
    fGMesh->SetName("gmesh with collapsed els");
    cmeshpressure->SetName("cmesh pressure");

    set<int> matids1d;
    for (auto gel : fGMesh->ElementVec())
    {
        if (!gel) continue;
        if (gel->Dimension() == fGMesh->Dimension()-1)
        {
            matids1d.insert(gel->MaterialId());
        }
    }
    cmeshpressure->ApproxSpace().CreateDisconnectedElements(true);
    cmeshpressure->AutoBuild(matids1d);

    // Creating volumetric collapsed elements
    set<int> matidstarget; // matid of the collapsed element (output parameter)
    // The updated mesh will be fGMesh
    // But the comp mesh used must be cmeshpressure because the connectivity of the multiphysics mesh
    // hasn't been updated yet.
    CreateCollapsedGeoEls(*cmeshpressure, matidstarget, matids1d);

    // Creating dim-1 elements: External flux and external pressure.
    CreateExternalElements(fGMesh, matidstarget, matids1d);

    // ********** CREATING THE COMPUTATIONAL MESH
    // At this point the code created all geometric elements needed for ther SBFEM simulation and it's stored in fGMesh
    // Here it's created the SBFemVolumeHdiv elements -> CompElHdivElement -> Adjust connectivities
    // But firstly the atomic meshes must be properly created
    CreateCompElPressure(*cmeshpressure, matids1d);
    // Adding the volumetric elements
    CreateSBFemVolumePressure(*cmeshpressure, matids1d, matidstarget);
    cmeshpressure->LoadReferences();

    CreateCompElFlux(*cmeshflux, matidstarget, matids1d);
    // Adding the volumetric elements
    CreateSBFemVolumeFlux(*cmeshflux, matids1d, matidstarget);
    cmeshflux->LoadReferences();

    // Then, creating the multiphysics mesh
    // Next method will setup the multiphysics mesh
    // Deleting the Hdiv elements existing in the mesh.
    // The elements must be specified as TPZCompElHDivSBFem elements
    CreateSBFemMultiphysicsMesh(cmeshm,matidstarget);
    cmeshm.ApproxSpace().SetAllCreateFunctionsSBFemMultiphysics(dim);
    cmeshm.SetName("multiphysicssbfem");
    TPZManVector<MSpaceConfig> active(2,ETestTrial);
    cmeshm.BuildMultiphysicsSpace(active, cmeshvec);
    cmeshm.LoadReferences();
    cmeshm.CleanUpUnconnectedNodes();
    {
        ofstream sout("cmeshmultiphysics.txt");
        cmeshm.Print(sout);
    }
    
    // Adding the interface elements
    AddInterfaceElements(cmeshm, matids1d);
    {
        ofstream sout("cmeshmultiphysics.txt");
        cmeshm.Print(sout);
    }

    // Creating the SBFEMVolumeHdiv elements
    CreateSBFemMultiphysicsVol(cmeshm, matids1d, matidstarget);

    // Based on the multiphysics mesh, create SBFemElementGroups
    // The SBFemMultiphysicsElGroups = Group of collapsed flux elements + pressure compels
    // Then, condense the flux into the pressures.
    CreateSBFemMultiphysicsElGroups(cmeshm, matidstarget);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();
    
    GroupandCondense(cmeshm);
    cmeshm.ExpandSolution();
    cmeshm.ComputeNodElCon();
    cmeshm.CleanUpUnconnectedNodes();

#ifdef PZDEBUG
    {
        ofstream mout("cmeshmultiphysicscondensed.txt");
        cmeshm.Print(mout);
    }
#endif

    auto nel = cmeshm.NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmeshm.Element(el);
        if (!cel)
        {
            continue;
        }
        auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup * >(cel);
        if (sbfemgr)
        {
            continue;
        }
        if (cel->Reference())
        {
            auto matid = cel->Reference()->MaterialId();
            auto it = matids1d.find(matid);
            if (it != matids1d.end())
            {
                continue;
            }
        }
        cmeshm.ElementVec()[el] = 0;
    }
#ifdef LOG4CXX
    if (loggercmesh->isDebugEnabled())
    {
        stringstream sout;
        cmeshm.Print(sout);
        sout << "flux\n";
        cmeshflux->Print(sout);
        sout << "pressure\n";
        ofstream pout("cmeshpressure.txt");
        cmeshpressure->Print(sout);
    }
#endif
}

// Creates geometric collapsed elements
void TPZBuildSBFemMultiphysics::CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget, set<int> & matids1d)
{
    // all computational elements for internal DOFs have been loaded for flux and pressure meshes
    set<int> matids;
    for (auto it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++)
    {
        int64_t mat = it->second;
        if (cmeshpressure.FindMaterial(mat))
        {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }

    // Creating geometric collapsed elements for the multiphysics mesh based on the pressure mesh
    auto gmeshpressure = cmeshpressure.Reference();
    auto dim = gmeshpressure->Dimension();

    for (auto gelpressure : fGMesh->ElementVec())
    {
        if (!gelpressure || gelpressure->HasSubElement())
        {
            continue;
        }
        // we create SBFemVolume elements by partitioning the volume elements
        auto el = gelpressure->Index();
        if (gelpressure->Dimension() != dim || fElementPartition[el] == -1)
        {
            continue;
        }
        // I am searching elements with dim dimension and with FE matids
        // to transform it in collapsed elements with matidtargets
        if (matids.find(gelpressure->MaterialId()) == matids.end())
        {
            continue;
        }
        int nsides = gelpressure->NSides();
        for (int is = 0; is<nsides; is++)
        {
            if (gelpressure->SideDimension(is) != dim-1)
            {
                continue;
            }            
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gelpressure,is);
            int onlyinterpolated = false;
            int removeduplicates = true;
            // we identify all computational elements connected to this element side
            gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);
            // we create a volume element based on all smaller elements linked to this side
            int ncelstack = celstack.NElements();
            for (int icel=0; icel<ncelstack; icel++) {
                TPZGeoElSide subgelside = celstack[icel].Reference();
                // we are only interested in faces
                if (subgelside.Dimension() != dim-1) {
                    continue;
                }
                int nnodes = subgelside.NSideNodes();
                TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
                int matid = fMatIdTranslation[gelpressure->MaterialId()];
                int64_t index;
                for (int in=0; in<nnodes; in++) {
                    Nodes[in] = subgelside.SideNodeIndex(in);
                }
                int elpartition = fElementPartition[el];
                for (int in=nnodes; in < 2*nnodes; in++) {
                    Nodes[in] = fPartitionCenterNode[elpartition];
                }
                if (subgelside.IsLinearMapping())
                {
                    switch(nnodes)
                    {
                        case 2:
                            gmeshpressure->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                            break;
                        case 4:
                            gmeshpressure->CreateGeoElement(ECube, Nodes, matid, index);
                            break;
                        case 3:
                            gmeshpressure->CreateGeoElement(EPrisma, Nodes, matid, index);
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                    
                }
                else
                {
                    int64_t elementid = gmeshpressure->NElements()+1;
                    switch(nnodes)
                    {
                        case 2:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmeshpressure,index);
                            break;
                        case 4:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmeshpressure,index);
                            break;
                        case 3:
                            gmeshpressure->CreateGeoElement(EPrisma, Nodes, matid, index);
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                }
                auto gel = fGMesh->ElementVec()[index];
                TPZVec<REAL> qsi(3,0);
                TPZFMatrix<REAL> jacobian, axes, jacinv;
                REAL detjac = 0;
                gel->JacobianXYZ(qsi, jacobian, axes, detjac, jacinv);
                if (detjac < 0)
                {
                    if(0 )cout << "The element " << index << " is not countclockwise\n";
                    TPZManVector<int64_t> nodescopy(Nodes);
                    gel->SetNodeIndex(0,Nodes[1]);
                    gel->SetNodeIndex(1,Nodes[0]);
                }
                if (index >= fElementPartition.size()) {
                    fElementPartition.resize(index+1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    fGMesh->ResetConnectivities();
    fGMesh->BuildConnectivity();

#ifdef PZDEBUG
    {
        ofstream gout("gmeshwithvol.txt");
        fGMesh->Print(gout);
        std::ofstream outvtk("GeometrySBFEM.vtk");
        TPZVTKGeoMesh vtk;
        vtk.PrintGMeshVTK(fGMesh, outvtk, true);
    }
#endif
}

// The order is: fDifPressure - fFluxRight - fFluxLeft - fAverPressure - fSkeleton
// fSkeleton is the internal flux and pressures.
void TPZBuildSBFemMultiphysics::CreateExternalElements(TPZAutoPointer<TPZGeoMesh> & gmesh, set<int> & matidtarget, set<int> &matids1d)
{
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel) continue;
        if (fElementPartition[gel->Index()] == -1) continue;

        // matidtarget contains the material id of the collapsed elements
        auto it = matidtarget.find(gel->MaterialId());
        if (it == matidtarget.end()) continue;
        if (gel->Dimension() < fGMesh->Dimension()) continue;

        fElementPartition.resize(fElementPartition.size()+6);
        auto idcollapsed = fElementPartition[gel->Index()];

        // getting the side in which the 1d elements will be constructed as neighbours
        auto iside = GetSideCollapsedEl(gel);

        TPZGeoElBC(gel, iside, fRightFlux);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;
        
        TPZGeoElBC(gel, iside, fInternal);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;

        TPZGeoElBC(gel, iside, fLeftFlux);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;

        TPZGeoElBC(gel, iside, fDifPressure);
        fElementPartition[gmesh->NElements()-1] = idcollapsed;
        
        TPZGeoElSide gelside(gel,iside);
        auto gelsideskeleton = gelside.HasNeighbour(matids1d);

        for (int i = 0; i < 2; i++)
        {
            gmesh->ElementVec()[gmesh->NElements()-1]->SetNodeIndex(0,gelsideskeleton.Element()->NodeIndex(0));
            gmesh->ElementVec()[gmesh->NElements()-1]->SetNodeIndex(1,gelsideskeleton.Element()->NodeIndex(1));
        }
    }
}

void TPZBuildSBFemMultiphysics::CreateCompElPressure(TPZCompMesh &cmeshpressure, set<int> & matids1d)
{
    // getting the materials
    TPZCompMesh pressurecopy(cmeshpressure);
    cmeshpressure.SetReference(fGMesh);
    auto matvec = cmeshpressure.MaterialVec();
    // map<int, TPZMaterial*> matveccopy(matvec);
    // matveccopy.insert(matvec.begin(), matvec.end());
    int d1 = fGMesh->Dimension()-1;
    auto nstate = cmeshpressure.MaterialVec().begin()->second->NStateVariables();

    // for (auto const& [key, val] : matvec)
    // {
    //     auto valbndcnd = dynamic_cast<TPZBndCondT<STATE> *>(val);
    //     if(!valbndcnd)
    //     {
    //         matveccopy[key] = val->NewMaterial();
    //     } else 
    //     {

    //         val->Clone(matveccopy);
    //     }
    // }

    auto dim = cmeshpressure.Dimension()-1; // materials with dim-1 dimensional elements
    auto matext = cmeshpressure.MaterialVec().begin()->second;
    auto POrder = cmeshpressure.GetDefaultOrder();

    auto matint = matext->NewMaterial();
    matint->SetId(fInternal);

    cmeshpressure.CleanUp();

    cmeshpressure.SetReference(fGMesh);
    cmeshpressure.SetDefaultOrder(POrder);
    cmeshpressure.SetDimModel(dim);
    cmeshpressure.SetAllCreateFunctionsContinuous();
    cmeshpressure.ApproxSpace().CreateDisconnectedElements(true);
    pressurecopy.CopyMaterials(cmeshpressure);
    // for (auto const& [key, val] : matveccopy)
    // {
    //     cmeshpressure.InsertMaterialObject(val);
    // }
    cmeshpressure.InsertMaterialObject(matint);
    set<int> matids = {fInternal};
    cmeshpressure.AutoBuild(matids);
    
    auto matleft = new TPZNullMaterial(fDifPressure, dim, nstate);
    cmeshpressure.InsertMaterialObject(matleft);

    // Pressure mesh will be composed of:
    // Internal DOFs - Materials already created!

    // 2rd. Dif. pressure
    matids = {fDifPressure};
    cmeshpressure.AutoBuild(matids);

    // 3rd. Aver. pressure
    cmeshpressure.AutoBuild(matids1d); 
    
    for(auto newnod : cmeshpressure.ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    cmeshpressure.LoadReferences();

}

void TPZBuildSBFemMultiphysics::CreateSBFemVolumePressure(TPZCompMesh & cmeshpressure, set<int> &matids1d, set<int> & matidtarget)
{
    cmeshpressure.Reference()->ResetReference();
    cmeshpressure.LoadReferences();
    // This for creates the volumetric multiphysics element, working as an AutoBuild()
    for (auto gelcollapsed : fGMesh->ElementVec())
    {
        // Now I'm searching for the collapsed element, in which the matid is in matidtarget.
        if (!gelcollapsed)
        {
            continue;
        }
        auto it = matidtarget.find(gelcollapsed->MaterialId());
        if (it == matidtarget.end())
        {
            continue;
        }

        auto idvol = fElementPartition[gelcollapsed->Index()];
        if (idvol == -1)
        {
            DebugStop();
        }
        
        auto cel = CreateSBFemPressureCompEl(cmeshpressure, gelcollapsed);
        int64_t index = cel->Index();

        // Getting the TPZSBFemVolumeHdiv compel
        auto celsbfem = cmeshpressure.Element(index);
        auto sbfem = dynamic_cast<TPZSBFemVolumeL2 * >(celsbfem);
        if(!sbfem)
        {
            DebugStop();
        }

        // Adding multiphysics compels:
        // The neighbour will be in the surface identified as side. So I get the gelside and search for its neighbour.
        // I'll include all computational elements related to geo collapsed element
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelside(gelcollapsed, side);

        // 1st element will be fDifPressure
        auto neigh = gelside.Neighbour();
        auto gel1d = neigh.Element();
        while (gel1d->MaterialId() != fDifPressure || fElementPartition[gel1d->Index()] != idvol)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
        }
        if(!(gel1d->Reference()))
        {
            DebugStop();
        }
        auto cel1d = gel1d->Reference();
        sbfem->AddElement1D(cel1d, 0);
        int count = 1;

        while (gelside != neigh)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
            int gel1dmat = gel1d->MaterialId();
            auto it = matids1d.find(gel1dmat);
            if(gel1d->Reference() && (gel1d->MaterialId() == fInternal) )
            {
                if(fElementPartition[gel1d->Index()] == idvol)
                {
                    cel1d = gel1d->Reference();
                    sbfem->AddElement1D(cel1d, 1);
                }
            }
            if(gel1d->Reference() && (it != matids1d.end() ) )
            {
                cel1d = gel1d->Reference();
                sbfem->AddElement1D(cel1d, 2);
                break;
            }
        }
    }
}

void TPZBuildSBFemMultiphysics::CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget, set<int> & matid1d)
{
    // Now creating the material for the external elements:
    auto dim = cmeshflux.Dimension()-1; // materials with dim-1 dimensional elements
    auto nstate = cmeshflux.MaterialVec().begin()->second->NStateVariables();
    auto matext = cmeshflux.MaterialVec().begin()->second;
    
    auto matinternal = new TPZNullMaterial(fInternal, dim, nstate);
    cmeshflux.InsertMaterialObject(matinternal);
    
    auto matboundleft = new TPZNullMaterial(fLeftFlux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundleft);
    
    auto matboundright = new TPZNullMaterial(fRightFlux, dim, nstate);
    cmeshflux.InsertMaterialObject(matboundright);

    map<int64_t,TPZCompEl *> geltocel;

    fGMesh->ResetReference();
    int previouspartition = -1;
    for (auto gelcollapsed : fGMesh->ElementVec())
    {
        if (!gelcollapsed) continue;
        auto it = matidtarget.find(gelcollapsed->MaterialId());
        if (it == matidtarget.end()) continue;

        int partition = fElementPartition[gelcollapsed->Index()];
        if (previouspartition != partition)
        {
            cmeshflux.Reference()->ResetReference();
        }
        
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelsidecollapsed(gelcollapsed, side);

        auto gelsideint = gelsidecollapsed.Neighbour();
        while (gelsideint.Element()->MaterialId() != fInternal)
        {
            gelsideint = gelsideint.Neighbour();
            if (gelsideint == gelsidecollapsed)
            {
                DebugStop();
            }
        }
        auto gel1d = gelsideint.Element();

        if(fElementPartition[gel1d->Index()] != fElementPartition[gelcollapsed->Index()])
        {
            DebugStop();
        }

        auto celhdivc = new TPZCompElHDivSBFem<pzshape::TPZShapeLinear>(cmeshflux, gel1d, gelsidecollapsed);
        geltocel[gel1d->Index()] = celhdivc;

        previouspartition = partition;
    }
    cmeshflux.ExpandSolution();
    cmeshflux.Reference()->ResetReference();

    // Now setting the geoelement for the HdivBound
    for (auto gelleft : fGMesh->ElementVec())
    {
        if (!gelleft || gelleft->MaterialId() != fLeftFlux) continue;

        TPZGeoElSide gelsideleft(gelleft, gelleft->NSides()-1);

        auto gelsideinternal = gelsideleft.Neighbour(); // This element if fInternal
        auto gelsideright = gelsideinternal.Neighbour(); // This element is fRightFlux
#ifdef PZDEBUG
        if(gelsideright.Element()->MaterialId() != fRightFlux || gelsideinternal.Element()->MaterialId() != fInternal)
        {
            DebugStop();
        }
#endif

        auto cel  = geltocel[gelsideinternal.Element()->Index()];
        TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * celhdivc = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(cel);
#ifdef PZDEBUG
        if (!celhdivc)
        {
            DebugStop();
        }
#endif
        if(fElementPartition[gelleft->Index()] != fElementPartition[gelsideinternal.Element()->Index()])
        {
            DebugStop();
        }
        if(fElementPartition[gelsideright.Element()->Index()] != fElementPartition[gelsideinternal.Element()->Index()])
        {
            DebugStop();
        }
        

        auto hdivboundleft = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gelsideleft.Element());
        celhdivc->SetConnectIndex(3,hdivboundleft->ConnectIndex(0));
        celhdivc->SetCompElFlux(hdivboundleft);
        gelsideleft.Element()->ResetReference();

        auto hdivboundright = new TPZCompElHDivBound2<pzshape::TPZShapeLinear>(cmeshflux,gelsideright.Element());
        celhdivc->SetConnectIndex(4,hdivboundright->ConnectIndex(0));
        gelsideright.Element()->ResetReference();
    }

    cmeshflux.LoadReferences();
    cmeshflux.CleanUpUnconnectedNodes();
    cmeshflux.ExpandSolution();

#ifdef PZDEBUG
    ofstream sout("cmeshflux0.txt");
    cmeshflux.Print(sout);
    ofstream gout("gmeshflux0.txt");
    cmeshflux.Reference()->Print(gout);
#endif
}

void TPZBuildSBFemMultiphysics::CreateSBFemVolumeFlux(TPZCompMesh & cmesh, set<int> &matids1d, set<int> & matidtarget)
{
    // This for creates the volumetric multiphysics element, working as an AutoBuild()
    for (auto gelcollapsed : fGMesh->ElementVec())
    {
        // Now I'm searching for the collapsed element, in which the matid is in matidtarget.
        if (!gelcollapsed)
        {
            continue;
        }
        auto it = matidtarget.find(gelcollapsed->MaterialId());
        if (it == matidtarget.end())
        {
            continue;
        }

        auto idvol = fElementPartition[gelcollapsed->Index()];
        if (idvol == -1)
        {
            DebugStop();
        }
        
        auto cel = CreateSBFemFluxCompEl(cmesh, gelcollapsed);
        const int64_t index = cel->Index();

        // Getting the TPZSBFemVolumeHdiv compel
        auto celsbfem = cmesh.Element(index);
        auto sbfem = dynamic_cast<TPZSBFemVolumeHdiv * >(celsbfem);
        if(!sbfem)
        {
            DebugStop();
        }

        // Adding multiphysics compels:
        // The neighbour will be in the surface identified as side. So I get the gelside and search for its neighbour.
        // I'll include all computational elements related to geo collapsed element
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelside(gelcollapsed, side);

        // 1st element will be fDifPressure
        auto neigh = gelside.Neighbour();
        auto gel1d = neigh.Element();
        while (gel1d->MaterialId() != fLeftFlux || fElementPartition[gel1d->Index()] != idvol)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
        }
        if(!(gel1d->Reference()))
        {
            DebugStop();
        }
        auto cel1d = gel1d->Reference();
        sbfem->AddElement1D(cel1d, 0);
        int count = 1;

        while (gelside != neigh)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
            if(gel1d->Reference() && (gel1d->MaterialId() == fInternal || gel1d->MaterialId() == fRightFlux))
            {
                if(fElementPartition[gel1d->Index()] == idvol)
                {
                    cel1d = gel1d->Reference();
                    sbfem->AddElement1D(cel1d, count);
                    count++;
                    if(count == 3) break;
                }
            }
        }
    }
}

void TPZBuildSBFemMultiphysics::CreateSBFemMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidstarget)
{
    // The materials for ESkeleton and Emat1 were already created before.
    auto dim = cmeshm.Dimension();
    auto id = cmeshm.MaterialVec().begin()->second->Id();
    auto nstate = cmeshm.MaterialVec().begin()->second->NStateVariables();
    auto matint = cmeshm.MaterialVec().begin()->second;
    
    {
        auto mat = new TPZNullMaterialCS<STATE>(fDifPressure, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = matint->NewMaterial();
        mat->SetId(fInternal);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZNullMaterialCS<STATE>(fLeftFlux, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZNullMaterialCS<STATE>(fRightFlux, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }

    {
        auto mat = new TPZLagrangeMultiplierCS<STATE>(fInterface, dim, nstate);
        cmeshm.InsertMaterialObject(mat);
    }
}

// Order of elements to add the interface:
// 1. fDifPressure
// 2. fLeftFlux
// 3. fInternal
// 4. fRightFlux
// 5. Skeleton/BC (External average pressure)
void TPZBuildSBFemMultiphysics::AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> &matids1d)
{
    for (auto geldifpr : fGMesh->ElementVec())
    {
        if (!geldifpr || geldifpr->MaterialId() != fDifPressure) continue;

        auto side = GetSideSkeletonEl(geldifpr);
        TPZGeoElSide gelsidedifpr(geldifpr, side);

        auto gelsideleft = gelsidedifpr.Neighbour(); // fLeftFlux
#ifdef PZDEBUG
        if(gelsideleft.Element()->MaterialId() != fLeftFlux) DebugStop();
#endif
        

        TPZCompElSide celsidedifpr = gelsidedifpr.Reference();
        TPZCompElSide celsideleft = gelsideleft.Reference();
#ifdef PZDEBUG
        if (!celsidedifpr || !celsideleft) DebugStop();
        if (fElementPartition[gelsidedifpr.Element()->Index()] != fElementPartition[gelsideleft.Element()->Index()])
        {
            DebugStop();
        }
#endif
        TPZGeoElBC gelbc(gelsideleft, fInterface);
        TPZMultiphysicsInterfaceElement *intl = new TPZMultiphysicsInterfaceElement(cmeshm, gelbc.CreatedElement(),celsidedifpr,celsideleft);
        fElementPartition[gelbc.CreatedElement()->Index()] = fElementPartition[gelsidedifpr.Element()->Index()];
    }

        for (auto gelright : fGMesh->ElementVec())
    {
        if (!gelright || gelright->MaterialId() != fRightFlux) continue;

        auto side = GetSideSkeletonEl(gelright);
        TPZGeoElSide gelsideright(gelright, side);

        auto gelsideaverpr = gelsideright.Neighbour();
        auto it = matids1d.find(gelsideaverpr.Element()->MaterialId());
        while (it == matids1d.end())
        {
            gelsideaverpr = gelsideaverpr.Neighbour();
            it = matids1d.find(gelsideaverpr.Element()->MaterialId());
#ifdef PZDEBUG
            if (gelsideaverpr.Element() == gelsideright.Element()) DebugStop();
#endif
        }

        TPZCompElSide celsideaverpr = gelsideaverpr.Reference();
        TPZCompElSide celsideright = gelsideright.Reference();
#ifdef PZDEBUG
        if (!celsideaverpr || !celsideright)
        {
            DebugStop();
        }
#endif
        TPZGeoElBC gelbc(gelsideright, fInterface);
        auto intr = new TPZMultiphysicsInterfaceElement(cmeshm, gelbc.CreatedElement(),celsideright,celsideaverpr);
        fElementPartition[gelbc.CreatedElement()->Index()] = fElementPartition[gelsideright.Element()->Index()];
    }
}

// Here SBFEM Multiphysics Volume elements are created.
// The Multiphysics element must be fSkeleton
// But I still need to pass the Volumetric Collapsed element so the TPZSBFemVolumeHdiv class will be able to perform integrations.
void TPZBuildSBFemMultiphysics::CreateSBFemMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matids1d, set<int> & matidtarget)
{
    // This for creates the volumetric multiphysics element, working as an AutoBuild()
    for (auto cel : cmeshm.ElementVec())
    {
        // Now I'm searching for the collapsed element, in which the matid is in matidtarget.
        if (!cel) continue;
        auto gelcollapsed = cel->Reference();
        
        auto it = matidtarget.find(gelcollapsed->MaterialId());
        if (it == matidtarget.end()) continue;

        auto idvol = fElementPartition[gelcollapsed->Index()];
        if (idvol == -1) DebugStop();

        // int64_t index;

        // Getting the TPZSBFemVolumeMultiphysics compel
        // auto celsbfem = cmeshm.Element(index);

        auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> * >(cel);
        if(!sbfem) continue;


        // Adding multiphysics compels:
        // The neighbour will be in the surface identified as side. So I get the gelside and search for its neighbour.
        // I'll include all computational elements related to geo collapsed element
        auto side = GetSideCollapsedEl(gelcollapsed);
        TPZGeoElSide gelside(gelcollapsed, side);

        // 1st element will be fDifPressure
        auto neigh = gelside.Neighbour();
        auto gel1d = neigh.Element();
        if(!gel1d || !(gel1d->Reference())) DebugStop();
        
        while (gel1d->MaterialId() != fDifPressure || fElementPartition[gel1d->Index()] != idvol)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
        }
        auto cel1d = gel1d->Reference();
        TPZManVector<int64_t> indexes(0);
        sbfem->SetConnectIndexes(indexes);
        sbfem->AddElement1D(cel1d, 0);
        int count = 1;

        while (gelside != neigh)
        {
            neigh = neigh.Neighbour();
            gel1d = neigh.Element();
            if(gel1d->Reference() && gel1d->Dimension() == fGMesh->Dimension()-1)
            {
                if(fElementPartition[gel1d->Index()] == idvol && count != 6)
                {
                    cel1d = gel1d->Reference();
                    sbfem->AddElement1D(cel1d, count);
                    count++;
                    if(count == 7) break;
                }
                else
                {
                    auto it = matids1d.find(gel1d->MaterialId());
                    if((fElementPartition[gel1d->Index()] == idvol || fElementPartition[gel1d->Index()] == -1) && it != matids1d.end() && count == 6)
                    {
                        cel1d = gel1d->Reference();
                        sbfem->AddElement1D(cel1d, count);
                        count++;
                        if(count == 7) break;
                    }
                    if(gel1d->MaterialId() == fSkeletonMatId && count == 6)
                    {
                        cel1d = gel1d->Reference();
                        sbfem->AddElement1D(cel1d, count);
                        count++;
                        if(count == 7) break;
                    }
                }
            }
        }
    }
    
}

void TPZBuildSBFemMultiphysics::GroupandCondense(TPZMultiphysicsCompMesh & cmeshm)
{
    for (auto cel : cmeshm.ElementVec())
    {
        if (!cel)
        {
            continue;
        }
        auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(cel);
        if (!sbfemgr)
        {
            continue;
        }
        sbfemgr->GroupandCondense(fCondensedMatids);
        cmeshm.ComputeNodElCon();
        cmeshm.CleanUpUnconnectedNodes();
    }
}

// ******** UTILITY FUNCTIONS

int TPZBuildSBFemMultiphysics::GetSideSkeletonEl(TPZGeoEl * gel)
{
    auto side = -1;
    switch (gel->Type())
    {
    case EOned:
        side = 2;
        break;
    case ETriangle:
        side = 6;
        break;
    case EQuadrilateral:
        side = 8;
    default:
        break;
    }
    return side;
}

int TPZBuildSBFemMultiphysics::GetSideCollapsedEl(TPZGeoEl * gel)
{
    auto side = -1;
    switch (gel->Type())
    {
    case EQuadrilateral:
        side = 4;
        break;
    case EPrisma:
        side = 15;
        break;
    case ECube:
        side = 20;
        break;
    default:
        break;
    }
    return side;
}

void TPZBuildSBFemMultiphysics::CreateSBFemMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmesh, set<int> & matidtarget)
{
    int64_t numgroups = fPartitionCenterNode.size();
    int64_t groupelementindices(numgroups);
    
    TPZManVector<int64_t> elementgroupindices(numgroups); 
    
    for (int64_t el=0; el<numgroups; el++)
    {
        
        TPZCompEl* cel = new TPZSBFemMultiphysicsElGroup(cmesh);
        elementgroupindices[el] = cel->Index();
    }
    int dim = cmesh.Dimension();
    for (auto cel : cmesh.ElementVec())
    {
        if (!cel) continue;
        
        auto sbfem = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cel);
        if (sbfem)
        {
            TPZGeoEl *gel = sbfem->Reference();
            if(gel->Dimension() != fGMesh->Reference()->Dimension()) continue;
            int64_t gelindex = gel->Index();

            int side = GetSideCollapsedEl(gel);
            TPZGeoElSide gelside(gel,side);
            int geldim = gel->Dimension();
            int nsidenodes = gel->NSideNodes(side);
            TPZManVector<int64_t,8> cornernodes(nsidenodes);
            for (int node = 0; node<nsidenodes; node++)
            {
                cornernodes[node] = gel->SideNodeIndex(side, node);
            }
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension() == geldim-1 && neighbour.Element()->Reference())
                {
                    int nsidenodesneighbour = neighbour.Element()->NCornerNodes();
                    if (nsidenodesneighbour == nsidenodes)
                    {
                        TPZManVector<int64_t,8> neighnodes(nsidenodesneighbour);
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            neighnodes[i] = neighbour.Element()->NodeIndex(i);
                        }
                        int numequal = 0;
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            if (neighnodes[i] == cornernodes[i]) {
                                numequal++;
                            }
                        }
                        if (numequal == nsidenodesneighbour) {
                            break;
                        }
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                // we are not handling open sides (yet)
                DebugStop();
            }
            int64_t skelindex = neighbour.Element()->Reference()->Index();
            sbfem->SetSkeleton(skelindex);
            
            int64_t gelgroup = fElementPartition[gelindex];
            if (gelgroup == -1) {
                DebugStop();
            }
            int64_t celgroupindex = elementgroupindices[gelgroup];
            auto celgr = cmesh.Element(celgroupindex);
            auto sbfemgr = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(celgr);
            if (!sbfemgr) {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
        }
    }
    
    for (auto index : elementgroupindices)
    {
        auto cel = cmesh.Element(index);
        auto sbfemgroup = dynamic_cast<TPZSBFemMultiphysicsElGroup *>(cel);
        if (!sbfemgroup)
        {
            DebugStop();
        }
        const TPZVec<TPZCompEl *> &subgr = sbfemgroup->GetElGroup();
        int64_t nsub = subgr.NElements();
        for (auto cels : subgr)
        {
            auto femvol = dynamic_cast<TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad> *>(cels);
            if (!femvol) {
                DebugStop();
            }
            femvol->SetElementGroupIndex(index);
        }
    }
}
