#ifndef GRAFNODEH
#define GRAFNODEH

#include <iostream>

using namespace std;

#include <string.h>
#include "pzconnect.h"
#include "pzgraphmesh.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphEl;
class TPZBlock;

class TPZGraphNode {

public:

   TPZGraphNode();
	TPZGraphNode(TPZConnect *cn, TPZGraphMesh *gm);
	~TPZGraphNode(void);

   //int ElIndex();
   int SequenceNumber() {return fSequenceNumber;}
   void SetSequenceNumber(int seqnum) {fSequenceNumber = seqnum;}
	void SetElement(TPZGraphEl *gel);
   void SetConnect(TPZConnect *connect);
   void SetGraphMesh(TPZGraphMesh *mesh);
	int NPoints();
	void SetPointNumber(long num);
	void DrawCo(TPZDrawStyle st = EDXStyle);
	void DrawSolution(int solutionid, TPZDrawStyle st = EDXStyle);
   void DrawSolution(TPZVec<int> &solutionid, TPZDrawStyle st= EDXStyle);
	void DrawSolution(TPZBlock &sol, TPZDrawStyle st = EDXStyle);

	long FirstPoint();

	void Print(ostream &out);

protected:
	TPZConnect *fConnect;
	TPZGraphMesh *fGraphMesh;
	TPZGraphEl *fGraphEl;
	long fPointNum;
   //int fElIndex

private:
   int fSequenceNumber;
};

#endif
