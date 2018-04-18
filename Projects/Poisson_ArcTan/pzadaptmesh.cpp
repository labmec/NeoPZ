/* Generated by Together */

#include "pzadaptmesh.h"
#include "pzgclonemesh.h"
#include "pzcclonemesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzfstrmatrix.h"
#include "pzstepsolver.h"
#include "pzintel.h"
#include "pzquad.h"
#include "TPZMaterial.h"
#include "pzonedref.h"

#include "TPZVTKGeoMesh.h"

using namespace std;

// Save information of the current mesh to compare with cloned mesh (geometric mesh plus computational mesh)
void SaveCompMesh(TPZCompMesh *cmesh, int timessave,TPZCompMesh *cmeshmodified=NULL,bool check=false);


TPZAdaptMesh::TPZAdaptMesh(int maxorder) {
    fReferenceCompMesh = 0;
    fGeoRef.Resize(0);
    fPatch.Resize(0);
    fPatchIndex.Resize(0);
    fElementError.Resize(0);
    fCloneMeshes .Resize(0);
    fFineCloneMeshes .Resize(0);
    fMaxP = maxorder;
}

TPZAdaptMesh::~TPZAdaptMesh() {
    CleanUp();
}

void TPZAdaptMesh::SetCompMesh(TPZCompMesh * mesh) {
    if(!mesh) {
        cout <<"TPZAdaptMesh::Error:\n computational reference mesh must not be NULL!\n";
        return;
    }
    CleanUp();
    fReferenceCompMesh = mesh;
    fReferenceCompMesh->SetDimModel(mesh->Dimension());
    int nel = fReferenceCompMesh->ElementVec().NElements();
    fElementError.Resize(nel);
    fElementError.Fill(0.);
}

void TPZAdaptMesh::SetMaxP(int maxp) {
    if(maxp < 1) {
        cout << "TPZAdaptMesh::Error : SetMaxP - maximum p order must be greater than 0... You given " << maxp << ". Trying to set maximum p to new value " << 1 << endl;
        maxp = 1;
    }
    if(maxp > 12) maxp = 12;
    fMaxP = maxp;
}
void TPZAdaptMesh::SetOneDMaxP(int maxp) {
    if(maxp < 1) {
        cout << "TPZAdaptMesh::Error : SetMaxP - maximum p order must be greater than 0... You given " << maxp << ". Trying to set maximum p to new value " << 1 << endl;
        maxp = 1;
    }
    if(maxp > 12) maxp = 12;
    TPZOneDRef::gMaxP = maxp;
}

void TPZAdaptMesh::CleanUp() {
    int i;
    for(i=0;i<fCloneMeshes .NElements();i++) {
        TPZGeoCloneMesh *gmesh = dynamic_cast<TPZGeoCloneMesh *> (fCloneMeshes [i]->Reference());
        gmesh->ResetReference();
        //Cesar July 2003 ->
        //If some reference element is not used to analyse error its fine clone mesh is not created!
        if(fFineCloneMeshes[i]) {
            fFineCloneMeshes[i]->LoadReferences();
            RemoveCloneBC(fFineCloneMeshes[i]);
            DeleteElements(fFineCloneMeshes[i]);
            delete fFineCloneMeshes[i];
        }
        
        fCloneMeshes[i]->LoadReferences();
        RemoveCloneBC(fCloneMeshes[i]);
        DeleteElements(fCloneMeshes[i]);
        delete fCloneMeshes[i];
        delete gmesh;
    }
    fCloneMeshes.Resize(0);
    fFineCloneMeshes.Resize(0);
}

void PrintGeoMeshAsCompMeshInVTKWithElementData(TPZGeoMesh *gmesh,char *filename,TPZVec<TPZVec<REAL> > &elData) {
	if(!gmesh || !elData.NElements())
		return;
	int i, j, size = gmesh->NElements();

	TPZVec<TPZVec<REAL> > DataElement;
	DataElement.Resize(size);
	int ndata = elData[0].NElements();
	for(i=0;i<size;i++) {
		DataElement[i].Resize(ndata);
		DataElement[i].Fill(0.0);
	}
	// Making dimension of the elements as data element
	for(i=0;i<size;i++) {
		TPZGeoEl *gel = gmesh->ElementVec()[i];
		if(gel && gel->Reference())
			for(j=0;j<ndata;j++)
				DataElement[i][j] = elData[gel->Reference()->Index()][j];
		else
			for(j=0;j<ndata;j++)
				DataElement[i][j] = -1.0;
	}
	// Printing geometric mesh to visualization in Paraview
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filename, DataElement);
}

TPZCompMesh * TPZAdaptMesh::GetAdaptedMesh(REAL &error, REAL &truerror, TPZVec<REAL> &ervec, 
                                           void(*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),
                                           TPZVec<REAL> &truervec,TPZVec<REAL> &effect,ofstream &out,int use_trueerror,MElementType eltype,int print) {
    int i;
    //clone analysis
    int cliter;
    int nelmesh = fReferenceCompMesh->ElementVec().NElements();
    fElementError.Resize(nelmesh);
    effect.Resize(nelmesh);
    truervec.Resize(nelmesh);
    ervec.Resize(nelmesh);
    ervec.Fill(0.);
    truervec.Fill(0.);
    effect.Fill(0.);
    fElementError.Fill(0.);
    
    //Used to evaluate the error when the true solution is provided======
    TPZManVector<int64_t> perm(nelmesh,0);
    TPZVec<REAL> auxerrorvec(nelmesh,0.);
    REAL minerror = 0.;
    //===================================================================
    
    //gets the geometric reference elements that will generate the patch
    GetReferenceElements();
    
    if (use_trueerror && f){
        for (i=0;i<fReferenceCompMesh->ElementVec().NElements();i++){
            TPZInterpolatedElement *el = dynamic_cast<TPZInterpolatedElement *> (fReferenceCompMesh->ElementVec()[i]);
            if (el) 
                auxerrorvec[i] = UseTrueError (el , f);
            else
                auxerrorvec[i] = 0.;
        }
        minerror = SortMinError(auxerrorvec,perm,0.65);
    } 
    
    //Generates the patch - put the data into fPatch and fPatchIndexes
    BuildReferencePatch();
    out << "NElements " << nelmesh << ". Number of Patchs " << fPatchIndex.NElements() - 1 << std::endl;
    std::cout << "NElements " << nelmesh << ". Number of Patchs " << fPatchIndex.NElements() - 1 << std::endl;
    //Creates the patch clones; fReferenceCompMesh is the original computational mesh
	
	// First, asserting connects in computational mesh
    fReferenceCompMesh->ComputeNodElCon();
    
    int printing = 0;
    if(printing) {
        ofstream test("testAdaptMesh.txt",ios::app);
        fReferenceCompMesh->Print(test);
        fReferenceCompMesh->Reference()->Print(test);
        test.close();
    }
    
    // Create the clone meshes and put the data in fCloneMeshes
    // Each element of the vector contains a computational mesh
    CreateClones();
    int ncl = fCloneMeshes.NElements();
    out << "Number of cloned meshes " << ncl << std::endl;
    std::cout << "Number of cloned meshes " << ncl << std::endl;
    fFineCloneMeshes.Resize(ncl);
	static int count = 0;

	bool gDebug = false;
	char saida[512];
	TPZVec<TPZVec<REAL> > ErrorVecData;
	ErrorVecData.Resize(nelmesh);
	for(i=0;i<nelmesh;i++)
		ErrorVecData[i].Resize(2);

    //Creates an uniformly refined mesh and evaluates the error
    for (cliter = 0; cliter<ncl; cliter++) {
		ErrorVecData.clear();
        //Análise dos Clones
        TPZGeoCloneMesh *gcmesh = dynamic_cast<TPZGeoCloneMesh *> (fCloneMeshes[cliter]->Reference());
        //       if (gcmesh && gcmesh->GetMeshRootElement()->MaterialId() < 0) continue;
        if(gcmesh->ReferenceElement(0)->MaterialId() <  0 || (use_trueerror && f && !HasTrueError(cliter,minerror,auxerrorvec))) {
            fFineCloneMeshes [cliter] = 0;
            continue;
        }
        fFineCloneMeshes[cliter] = fCloneMeshes[cliter]->UniformlyRefineMesh(fMaxP,print);
		printing = 0;
        if(printing) {
            sprintf(saida,"OutputE%d_Cliter%d.txt",((int)eltype),cliter);
            ofstream outtemp(saida);
            fCloneMeshes[cliter]->Print(outtemp);
			fFineCloneMeshes[cliter]->Print(outtemp);
            outtemp.close();
        }
        std::cout << cliter << " " << fCloneMeshes[cliter]->NElements() << " " << fCloneMeshes[cliter]->NEquations() << " ";
        fCloneMeshes[cliter]->MeshError(fFineCloneMeshes [cliter],fElementError,f,truervec,out);
		// Printing errors on geometric mesh to validate
		memset(saida,0,512);
		if(gDebug) {
			for(i=0;i<nelmesh;i++) {
				ErrorVecData[i][0] = fElementError[i];
				ErrorVecData[i][1] = truervec[i];
			}
			
			sprintf(saida,"ErroOnCMesh%02dDFrom_Patch%d_Iter%d.vtk",fReferenceCompMesh->Dimension(),cliter,count);
			PrintGeoMeshAsCompMeshInVTKWithElementData((TPZGeoMesh *)gcmesh,saida,ErrorVecData);
		}
    }
	// To validate we print all the calculated errors by geo mesh
	ErrorVecData.clear();
	ErrorVecData.Resize(nelmesh);
	for(i=0;i<nelmesh;i++)
		ErrorVecData[i].Resize(3);
	for(i=0;i<nelmesh;i++) {
		ErrorVecData[i][0] = fElementError[i];
		ErrorVecData[i][1] = truervec[i];
		if(IsZero(truervec[i]))
			ErrorVecData[i][2] = fElementError[i];
		else
			ErrorVecData[i][2] = fElementError[i]/truervec[i];
	}
    
	// Printing calculated errors for each element in computational mesh
	sprintf(saida,"ErroOnCMesh%02dDFrom_Iter%d.vtk",fReferenceCompMesh->Dimension(),count++);
	PrintGeoMeshAsCompMeshInVTKWithElementData(fReferenceCompMesh->Reference(),saida,ErrorVecData);

    //Ordena o vetor de erros
    for(i=0;i<nelmesh;i++) {
        perm[i] = i;
        ervec[i]=fElementError[i];
    }
    Sort(fElementError,perm);
    //somatorio dos componentes do vetor de erro
    for(i=0;i<nelmesh;i++) error += fElementError[i];
    
    // Determining minimum error for applying p-refinement or h-refinement     // It's WRONG
    REAL minimumerror = 1.0, auxerror = 0.0;
    int counter = 0;
    for(i=0;i<nelmesh;i++) {
        auxerror += fElementError[perm[i]];
        if(!IsZero(fElementError[perm[i]])) {
            counter++;
//            minimumerror = (minimumerror < fElementError[perm[i]] ? minimumerror : fElementError[perm[i]]);
        }
    }
    if(counter)
        minimumerror = fabs(auxerror/counter) - 100*ZeroTolerance();
    
    if(f) {
        for(i=0; i<nelmesh; i++) {
            truerror += truervec[i];
        }
    }
    
    //inicializa effect com o tamanho de trueeerror
    effect.Resize(truervec.NElements());
    effect.Fill(0.);
    if(f) {
        for(i=0; i<nelmesh; i++) {
            if(truervec[i] >= 1.e-4*truerror && truervec[i] >= 5e-20) {
                effect[i] = ervec[i]/truervec[i];
            }
            else {
                effect[i]=ervec[i];
            }
        }
    }
    
    TPZStack <TPZGeoEl*> gelstack;
    TPZStack <int> porder;
    
    //Analyse clone element error and, if necessary, analyse element and changes its refinement pattern
    //std::cout << "Aplying ref pattern after analisis of the error on the clone meshes.\n";
    for (i=0;i<ncl;i++) {
        if(!fFineCloneMeshes[i]) continue;
        fCloneMeshes[i]->ApplyRefPattern(minimumerror,fElementError,fFineCloneMeshes[i],gelstack,porder);
    }
    
    TPZCompMesh *adapted = CreateCompMesh(fReferenceCompMesh,gelstack,porder);
    // recording in disk the computational mesh data
    static int countermesh = 0;
    if(print > 0) {
        SaveCompMesh(fReferenceCompMesh,countermesh++);
    }
    
	ErrorVecData.clear();
    std::cout << endl;
    return adapted;
}

void TPZAdaptMesh::GetReferenceElements() {
    if(!fReferenceCompMesh) {
        cout << "TPZAdaptMesh::Error:\n computational mesh must be initialized to call GetReferenceElements!\n";
        return;
    }
    std::set<TPZGeoEl*> georef;
    //  fReferenceCompMesh->GetRefPatches(fGeoRef);
    fReferenceCompMesh->GetRefPatches(georef);
    
#ifndef CLONEBCTOO
    //This will exclude geometric elements associated to bc from clone creation
    std::set<TPZGeoEl *>::iterator it;
    for(it=georef.begin();it!=georef.end();it++) {
        int id = (*it)->MaterialId();
        if (id > 0) {
            fGeoRef.Push(*it);
        }
    }
#else
    std::set<TPZGeoEl *>::iterator it;
    for(it=georef.begin();it!=georef.end();it++) {
        fGeoRef.Push(*it);
    }
#endif
}

void TPZAdaptMesh::BuildReferencePatch() {
    
    // the fGeoRef elements are a partition of the computational domain (should be)
    // create a computational element based on each reference element
    TPZGeoMesh *gmesh = fReferenceCompMesh->Reference();
    gmesh->ResetReference();
    TPZCompMesh *tmpcmesh = new TPZCompMesh (gmesh);
    int64_t i,j;
    for (i=0;i<fGeoRef.NElements();i++){
        int64_t index;
        tmpcmesh->CreateCompEl(fGeoRef[i],index);
    } 
    tmpcmesh->CleanUpUnconnectedNodes();
	tmpcmesh->ExpandSolution();
    TPZStack <int64_t> patchelindex;
    TPZStack <TPZGeoEl *> toclonegel;
    TPZStack<int64_t> elgraph;
    TPZVec<int64_t> n2elgraph;
    TPZVec<int64_t> n2elgraphid;
    TPZVec<int64_t> elgraphindex;

    tmpcmesh->GetNodeToElGraph(n2elgraph,n2elgraphid,elgraph,elgraphindex);
    // we use the  node to elgraph structure to decide which elements will be included
    int64_t clnel = tmpcmesh->NElements();
    // clnel corresponds to the number of patches
    // fPatch and fPatchIndex form a compacted list which form the patches.
    // Boundary elements will be added to each patch.
    fPatchIndex.Push(0);
    for (int64_t ipatch=0; ipatch<clnel; ipatch++){
        tmpcmesh->GetElementPatch(n2elgraph,n2elgraphid,elgraph,elgraphindex,ipatch,patchelindex);
        for (j=0; j<patchelindex.NElements(); j++){
            TPZGeoEl *gel = tmpcmesh->ElementVec()[patchelindex[j]]->Reference();
            if(gel) fPatch.Push(gel);
        }
        int64_t sum = fPatch.NElements();
        fPatchIndex.Push(sum);
    }
	
#ifdef PZDEBUG2 
	// CAJU TOOL
	{
		std::string filename("cMeshVtk.");
		{
			std::stringstream finalname;
			finalname << filename << 0 << ".vtk";
			ofstream file(finalname.str().c_str());
			/** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
			//static void PrintCMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = false);
			TPZVTKGeoMesh::PrintCMeshVTK(gmesh,file,true);
		}
		for (int64_t ip=0; ip<clnel; ip++) {
			int64_t firstindex = fPatchIndex[ip];
			int64_t lastindex = fPatchIndex[ip+1];
			gmesh->ResetReference();
			tmpcmesh->LoadReferences();
			std::set<TPZGeoEl *> loaded;
			for (int64_t ind=firstindex; ind<lastindex; ind++) {
				TPZGeoEl *gel = fPatch[ind];
				loaded.insert(gel);
			}
			int64_t ngel = gmesh->NElements();
			for (int64_t el=0; el<ngel; el++) {
				TPZGeoEl *gel = gmesh->ElementVec()[el];
				if (!gel) {
					continue;
				}
				if (gel->Reference() && loaded.find(gel) == loaded.end()) {
					gel->ResetReference();
				}
			}
			std::stringstream finalname;
			finalname << filename << ip+1 << ".vtk";
			ofstream file(finalname.str().c_str());
			/** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
			//static void PrintCMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = false);
			TPZVTKGeoMesh::PrintCMeshVTK(gmesh,file,true);
		}
	}
#endif
	// cleaning reference to computational elements into temporary cmesh
    gmesh->ResetReference();
    delete tmpcmesh;
	// loading references between geometric and computational meshes (originals)
    fReferenceCompMesh->LoadReferences();
}

void TPZAdaptMesh::CreateClones(){
	// asserting references of the original meshes
    fReferenceCompMesh->Reference()->ResetReference();
    fReferenceCompMesh->LoadReferences();
    TPZGeoMesh *geomesh = fReferenceCompMesh->Reference();
    
    int64_t clid,elid;
    for (clid=0; clid<fPatchIndex.NElements()-1;clid++) {
		// making clone of the original geometric mesh, only to construct computational clone
        TPZGeoCloneMesh *geoclone = new TPZGeoCloneMesh(geomesh);
        TPZStack<TPZGeoEl*> patch;
        for (elid=fPatchIndex[clid];elid<fPatchIndex[clid+1];elid++){
            patch.Push(fPatch[elid]);
        }
        geoclone->SetElements(patch,fGeoRef[clid]);

		int printing = 0;
        if(printing) {
            ofstream out("testAdaptMesh.txt",ios::app);
            geoclone->Print(out);
        }
        
        TPZCompCloneMesh *clonecompmesh = new TPZCompCloneMesh(geoclone,fReferenceCompMesh);
        clonecompmesh->AutoBuild();
/*#ifdef LOG4CXX
        {
            TPZFileStream fstr;
            std::stringstream sout;
            sout << "LOG/" << (void*) clonecompmesh << ".txt";
            std::string filename(sout.str());
            fstr.OpenWrite(filename);
            clonecompmesh->Reference()->Write(fstr,1);
            clonecompmesh->Write(fstr,1);
        }
#endif */
		// Computational mesh clone is stored
        fCloneMeshes.Push(clonecompmesh);    
    }
}

void TPZAdaptMesh::Sort(TPZVec<REAL> &vec, TPZVec<int64_t> &perm) {
    int64_t i,j;
    int imin = 0;
    int64_t imax = vec.NElements();
    for(i=imin; i<imax; i++) {
        for(j=i+1; j<imax; j++) {
            if(vec[perm[i]] < vec[perm[j]]) {
                int64_t kp = perm[i];
                perm[i] = perm[j];
                perm[j] = kp;
            }
        }
    }
}

void TPZAdaptMesh::HeapSort(TPZVec<REAL> &sol, TPZVec<int64_t> &perm){
    
    int64_t nelem = perm.NElements();
    int64_t i,j;
    for(i=0; i<nelem; i++) perm[i] = i;
    
    if(nelem == 1) return;
    int64_t l, ir,ind;
    REAL q;
    l= nelem/2;
    ir = nelem-1;
    while(l>0 && ir>0) {
        if(l> 0) {
            l--;
            ind = perm[l];
            q=sol[ind];
        } else {
            ind = perm[ir];
            q = sol[ind];
            perm[ir] = perm[0];
            ir--;
        }
        i=l;
        j=l+l+1;
        while(j<=ir) {
            if(j<ir && sol[perm[j]] < sol[perm[j+1]]) j++;
            if(q < sol[perm[j]]) {
                perm[i] = perm[j];
                i=j;
                j= i+i+1;
            } else {
                break;
            }
        }
        perm[i] = ind;
    }
}

TPZCompMesh *TPZAdaptMesh::CreateCompMesh (TPZCompMesh *mesh,                                          //malha a refinar
                                           TPZVec<TPZGeoEl *> &gelstack,   //
                                           TPZVec<int> &porders) {
    
    //Cria um ponteiro para a malha geométrica de mesh
    TPZGeoMesh *gmesh = mesh->Reference();
    if(!gmesh) {
        cout << "TPZAdaptMesh::CreateCompMesh encountered no geometric mesh\n";
        return 0;
    }
    
    //Reseta as referências do ponteiro para a malha geométrica criada
    //e cria uma nova malha computacional baseada nesta malha geométrica
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(mesh->Dimension());
    
    //Cria um clone do vetor de materiais da malha mesh
    mesh->CopyMaterials(*cmesh);

    //Idenifica o vetor de elementos computacionais de mesh
    int64_t el,nelem = gelstack.NElements();
    for(el=0L; el<nelem; el++) {
        
        //identifica os elementos geométricos passados em gelstack
        TPZGeoEl *gel = gelstack[el];
        if(!gel) {
            cout << "TPZAdaptMesh::CreateCompMesh encountered an null element\n";
            continue;
        }
        int64_t celindex;
        
        //Cria um TPZIntel baseado no gel identificado
        TPZInterpolatedElement *csint;
        csint = dynamic_cast<TPZInterpolatedElement *> (cmesh->CreateCompEl(gel,celindex));
        if(!csint) continue;
        
        csint->PRefine(porders[el]);
    }
#ifndef CLONEBCTOO
    nelem = gmesh->NElements();
    for (el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if (!gel || gel->Reference()) {
            continue;
        }
        int matid = gel->MaterialId();
        if (matid < 0) {
            TPZStack<TPZCompElSide> celstack;
            int ns = gel->NSides();
            TPZGeoElSide gelside(gel,ns-1);
            gelside.HigherLevelCompElementList2(celstack, 1, 1);
            if (celstack.size()) {
                TPZStack<TPZGeoEl *> subels;
                gel->Divide(subels);
            }
        }
    }
    nelem = gmesh->NElements();
    for (el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->ElementVec()[el];
        if (!gel || gel->Reference()) {
            continue;
        }
        int matid = gel->MaterialId();
        if (matid < 0) {
            TPZStack<TPZCompElSide> celstack;
            int ns = gel->NSides();
            TPZGeoElSide gelside(gel,ns-1);
            gelside.EqualLevelCompElementList(celstack, 1, 0);
            if (celstack.size()) {
                int64_t index;
                cmesh->CreateCompEl(gel, index);
            }
        }
    }
#endif

    cmesh->AdjustBoundaryElements();
    return cmesh;
}

void TPZAdaptMesh::RemoveCloneBC(TPZCompMesh *mesh) {
    int64_t nelem = mesh->NElements();
    for(int64_t iel=0L; iel<nelem; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if(!cel) continue;
		if(!cel->Material()) {
			delete cel;
			cel = 0;
			continue;
		}
        int matid = cel->Material()->Id();
        if(matid == -1000) delete cel;
    }
}

void TPZAdaptMesh::DeleteElements(TPZCompMesh *mesh)
{
    int64_t nelem = mesh->NElements();
    int64_t iel;
    for(iel=0; iel<nelem; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if(!cel) continue;
        TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
        if(!cint) continue;
        while(cint->HasDependency()) {
            TPZInterpolatedElement *large = LargeElement(cint);
            TPZInterpolatedElement *nextlarge = LargeElement(large);
            while(nextlarge != large) {
                large = nextlarge;
                nextlarge = LargeElement(large);
            }
            large->RemoveSideRestraintsII(TPZInterpolatedElement::EDelete);
            delete large;
        }
        cint->RemoveSideRestraintsII(TPZInterpolatedElement::EDelete);
        delete cint;
    }
}

TPZInterpolatedElement * TPZAdaptMesh::LargeElement(TPZInterpolatedElement *cint)
{
    int64_t nc = cint->NConnects();
    int side;
    TPZInterpolatedElement *result = cint;
    for(side=0; side<nc; side++) {
        if(cint->Connect(side).HasDependency()) {
            TPZCompElSide cintside(cint,side);
            TPZCompElSide large = cintside.LowerLevelElementList(1);
            if(!large.Exists()) {
                cout << "TPZAdaptMesh::DeleteElements I dont understand\n";
                large = cintside.LowerLevelElementList(1);
				return cint;
            }
            result = dynamic_cast<TPZInterpolatedElement *> (large.Element());
            break;
        }
    }
    return result;
}


REAL TPZAdaptMesh::UseTrueError(TPZInterpolatedElement *coarse, 
                                void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)){
    if (coarse->Material()->Id() < 0) return 0.0;
    
    REAL error = 0.;
    
    //  REAL loclocmatstore[500] = {0.},loccormatstore[500] = {0.};
    // TPZFMatrix loccormat(locmatsize,cormatsize,loccormatstore,500);
    
    //Cesar 25/06/03 - Uso a ordem máxima???
    TPZAutoPointer<TPZIntPoints> intrule = coarse->GetIntegrationRule().Clone();
    
    int dimension = coarse->Dimension();
    int numdof = coarse->Material()->NStateVariables();
        
    //TPZSolVec corsol;
    //TPZGradSolVec cordsol;
    TPZGradSolVec cordsolxy;
//    TPZVec<REAL> corsol(numdof);
//    TPZFNMatrix<9,REAL> cordsol(dimension,numdof),cordsolxy(dimension,numdof);
    
    TPZManVector<int> order(dimension,20);
    intrule->SetOrder(order);
    
    // derivative of the shape function
    // in the master domain
    TPZManVector<REAL,3> coarse_int_point(dimension);
//    TPZFNMatrix<9,REAL> jaccoarse(dimension,dimension),jacinvcoarse(dimension,dimension);
//    TPZFNMatrix<9,REAL> axescoarse(3,3);
//    TPZManVector<REAL,3> xcoarse(3);
    TPZFNMatrix<9,REAL> axesinner(3,3);
    
    
    TPZMaterialData datacoarse;
    coarse->InitMaterialData(datacoarse);
    
//    REAL jacdetcoarse;
    int numintpoints = intrule->NPoints();
    REAL weight;
    
    TPZVec<STATE> truesol(numdof);
    TPZFMatrix<STATE> truedsol(dimension,numdof);
    for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
        intrule->Point(int_ind,coarse_int_point,weight);
        //coarse->Reference()->X(coarse_int_point, xcoarse);
        coarse->Reference()->X(coarse_int_point, datacoarse.x);
        //if(f) f(xcoarse,truesol,truedsol);
        if(f) f(datacoarse.x,truesol,truedsol);

//        coarse->Reference()->Jacobian(coarse_int_point, jaccoarse, axescoarse, jacdetcoarse, jacinvcoarse);
        coarse->Reference()->Jacobian(coarse_int_point, datacoarse.jacobian, datacoarse.axes, datacoarse.detjac, datacoarse.jacinv);
        //weight *= fabs(jacdetcoarse);
        weight *= fabs(datacoarse.detjac);
//Er        int iv=0;
//        corsol[0].Fill(0.);
//        cordsol[0].Zero();

        //coarse->ComputeSolution(coarse_int_point, corsol, cordsol, axescoarse);
        coarse->ComputeShape(coarse_int_point,datacoarse);
        coarse->ComputeSolution(coarse_int_point,datacoarse);
        
        //int nc = cordsol[0].Cols();
        int64_t nc = datacoarse.dsol[0].Cols();
        for (int64_t col=0; col<nc; col++)
        {
            for (int d=0; d<dimension; d++) {
                REAL deriv = 0.;
                for (int d2=0; d2<dimension; d2++) {
                    deriv += datacoarse.dsol[0](d2,col)*datacoarse.axes(d2,d);
                }
               // cordsolxy[0](d,col) = deriv;
            }
        }
        int64_t jn;
        for(jn=0; jn<numdof; jn++) {
            for(int d=0; d<dimension; d++) {
                error += (datacoarse.dsol[0](d,jn)-truedsol(d,jn))*(datacoarse.dsol[0](d,jn)-truedsol(d,jn))*weight;
            }
        }
    }
    return error;
}


REAL TPZAdaptMesh::SortMinError (TPZVec<REAL> errvec, TPZVec<int64_t> perm, REAL percenterror){
    //Ordena o vetor de erros
    int64_t nelem = errvec.NElements();
    int64_t i;
    REAL error = 0.0;
    for(i=0; i<nelem; i++) {
        perm[i] = i;
        //ervec[i]=fElementError[i];
    }
    Sort(fElementError,perm);
    //somatório dos componentes do vetor de erro
    for(i=0; i<nelem; i++) error += fElementError[i];
    REAL ninetyfivepercent = 0.,auxerror = 0.;
    for(i=0;i<nelem;i++){
        auxerror += fElementError[perm[i]];
        if (auxerror >=  percenterror*error){
            ninetyfivepercent = fElementError[perm[i]];
            break;
        }
    }
    return ninetyfivepercent;
}

int TPZAdaptMesh::HasTrueError(int64_t clindex, REAL &minerror, TPZVec<REAL> &ervec){
    
    TPZGeoCloneMesh *gmesh = dynamic_cast<TPZGeoCloneMesh *> (fCloneMeshes [clindex]->Reference());
    int nref = gmesh->NReference();
    int i;
    for (i=0;i<nref;i++){
        TPZGeoEl *gel = gmesh->ReferenceElement(i);
        TPZCompEl *cel = gel->Reference();
        if (!cel) continue;
        int celindex = cel->Index();
        if (ervec[celindex] >= minerror) return 1;
    }
    return 0;
}
