#ifndef PZGRAPHELQ2D
#define PZGRAPHELQ2D

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 */
/// To export a graphical two dimensional element
class TPZGraphElQ2d : public TPZGraphEl {
public:

  TPZGraphElQ2d(TPZCompEl *cel, TPZGraphMesh *gmesh);

  virtual ~TPZGraphElQ2d(void);

  virtual int NConnects(){ return 9;}

  virtual MElementType Type() { return EQuadrilateral;}
	
  virtual int ExportType(TPZDrawStyle st);
	
  virtual int NNodes();
	
  virtual TPZGraphNode *Connect(int i);

  virtual int NPoints(TPZGraphNode *n);

  virtual int NElements();

  virtual void SetNode(int i,TPZGraphNode *n);

  virtual long EqNum(TPZVec<int> &co);

  virtual void Connectivity(TPZDrawStyle st = EDXStyle);


protected:

virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr);

virtual void NextIJ(int connect,TPZVec<int> &co, int incr);

//	void ComputeSequence(TPZGraphNode *n,int *ibound,int *incr);

protected :   

	TPZGraphNode *fConnects[9];


};

#endif
