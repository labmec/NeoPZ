// -*- c++ -*-
/* Generated by Together */

#include "pzstrmatrix.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzelmat.h"
#include "pzcompel.h"

#include "pzelgq2d.h"
#include "pzgnode.h"
#include "pzmat2dlin.h"

TPZStructMatrix::~TPZStructMatrix() {}

TPZMatrix *TPZStructMatrix::CreateAssemble(TPZFMatrix &rhs) {
  cout << "TPZStructMatrix::CreateAssemble should never be called\n";
  return 0;
}

TPZMatrix *TPZStructMatrix::Create() {
  cout << "TPZStructMatrix::Create should never be called\n";
  return 0;
}

TPZStructMatrix *TPZStructMatrix::Clone() {
  cout << "TPZStructMatrix::Clone should never be called\n";
  return 0;
}

void TPZStructMatrix::Assemble(TPZMatrix & stiffness, TPZFMatrix & rhs){

  int iel;
  int numel = 0, nelem = fMesh->NElements();
  TPZElementMatrix ek,ef;
  TPZManVector<int> destinationindex(0);
  TPZManVector<int> sourceindex(0);
  REAL stor1[1000],stor2[1000],stor3[100],stor4[100];
  ek.fMat = new TPZFMatrix(0,0,stor1,1000);
  ek.fConstrMat = new TPZFMatrix(0,0,stor2,1000);
  ef.fMat = new TPZFMatrix(0,0,stor3,100);
  ef.fConstrMat = new TPZFMatrix(0,0,stor4,100);

  TPZAdmChunkVector<TPZCompEl *> &elementvec = fMesh->ElementVec();

  for(iel=0; iel < nelem; iel++) {
    TPZCompEl *el = elementvec[iel];
    if(!el) continue;
    //		int dim = el->NumNodes();
    el->CalcStiff(ek,ef);
      //ek.fMat->Print(out);
    //ef.fMat->Print();
/*    if(!(numel%20)) cout << endl << numel;   //Jorge 8/7/2001
//    cout << '*';
//	cout.flush(); */
    numel++;

    if(!el->HasDependency()) {
      //ek.fMat->Print("stiff has no constraint",test);
      //ef.fMat->Print("rhs has no constraint",test);
      //test.flush();
      destinationindex.Resize(ek.fMat->Rows());
      int destindex = 0;
      int numnod = ek.NConnects();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.ConnectIndex(in);
         TPZConnect &np = fMesh->ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fMesh->Block().Position(blocknumber);
         int ndf = fMesh->Block().Size(blocknumber);
	 //	 if (numnod == 27){
	 //   cout << "First equation " << firsteq <<"\t ndf " << ndf << endl;
	 //	 }
         for(int idf=0; idf<ndf; idf++) {
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      stiffness.AddKel(*ek.fMat,destinationindex);
      rhs.AddFel(*ef.fMat,destinationindex);                 //  ??????????? Erro
    } else {
      // the element has dependent nodes
      el->ApplyConstraints(ek,ef);
      //ek.fMat->Print("stif no constraint",test);
      //ek.fConstrMat->Print("stif constrained",test);
      //ef.fMat->Print("rhs no constraint",test);
      //ef.fConstrMat->Print("rhs constrained",test);
      //test.flush();
      //test << "sum of columns\n";
      int destindex = 0;
      int fullmatindex = 0;
      destinationindex.Resize(ek.fConstrMat->Rows());
      sourceindex.Resize(ek.fConstrMat->Rows());
      int numnod = ek.fConstrConnect.NElements();
      for(int in=0; in<numnod; in++) {
         int npindex = ek.fConstrConnect[in];
         TPZConnect &np = fMesh->ConnectVec()[npindex];
         int blocknumber = np.SequenceNumber();
         int firsteq = fMesh->Block().Position(blocknumber);
         int ndf = fMesh->Block().Size(blocknumber);
         if(np.HasDependency()) {
           fullmatindex += ndf;
           continue;
         }
         for(int idf=0; idf<ndf; idf++) {
           sourceindex[destindex] = fullmatindex++;
           destinationindex[destindex++] = firsteq+idf;
         }
      }
      sourceindex.Resize(destindex);
      destinationindex.Resize(destindex);
      stiffness.AddKel(*ek.fConstrMat,sourceindex,destinationindex);
      rhs.AddFel(*ef.fConstrMat,sourceindex,destinationindex);
/*
if(ek.fConstrMat->Decompose_LU() != -1) {
    el->ApplyConstraints(ek,ef);
    ek.Print(*this,check);
    check.flush();
}
*/
    }
  }//fim for iel
//  cout << endl;
}

TPZStructMatrix::TPZStructMatrix(TPZCompMesh *mesh){
    fMesh = mesh;
}

