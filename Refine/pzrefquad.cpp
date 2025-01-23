/**
 * @file
 * @brief Contains the implementation of the TPZRefQuad methods. 
 */

#include "pzrefquad.h"
#include "pzgeoquad.h"
#include "pzshapequad.h"
#include "TPZGeoElement.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

using namespace pzshape;
using namespace std;

namespace pzrefine {
	
	static int nsubeldata[9] = {1,1,1,1,3,3,3,3,9};
	
	static int subeldata[9][9][2] =
	{
		{{0,0}},/*side=0 {isub0{0,1},isub1{0,1},isub2{0,1},...}*/
		{{1,1}},/*side=1*/
		{{2,2}},/*side=2*/
		{{3,3}},/*side=3*/
		{{0,1},{0,4},{1,4}},/*side=4*/  // {{0,4},{0,1},{1,4}},/*side=4*/
		{{1,2},{1,5},{2,5}},/*side=5*/
		{{2,3},{2,6},{3,6}},/*side=6*/
		{{3,0},{3,7},{0,7}},/*side=7*/         //{0,2},{0,8},{1,8},{2,8},{3,8},{0,5},{0,6},{2,4},{2,7}}/*side=8*/
		{{0,2},{0,5},{0,6},{2,4},{2,7},{0,8},{1,8},{2,8},{3,8}}/*side=8*/  
	};
	
	static int MidSideNodes[5][2] = {
		{0,1},//side 4
		{1,2},//side 5
		{2,3},//side 6
		{3,0},//side 7
		{0,2}	//side 8
	};
	
	
	static REAL MidCoord[5][2] = { 
		{0.,-1.},//side 4
		{1.,0.}, //side 5
		{0.,1.}, //side 6
		{-1.,0.},//side 7
		{0.,0.}  //side 8
	};
	
	/**
	 * define as conectividades entre sub-elementos
	 * linha i é filho i, {a,b,c} = {lado do filho atual,
	 * irmão vizinho,lado do vizinho}
	 */
	
	const int NumInNeigh = 5;
	static int InNeigh[4][NumInNeigh][3] = {
		{{1,1,0},{2,1,3},{3,3,0},{5,1,7},{6,3,4}},
		{{0,0,1},{2,2,1},{3,2,0},{6,2,4},{7,0,5}},
		{{0,3,1},{1,1,2},{3,3,2},{4,1,6},{7,3,5}},
		{{0,0,3},{1,0,2},{2,2,3},{4,0,6},{5,2,7}}
	};
	
	/**
	 * define os cantos locais dos fihos
	 */
	static int CornerSons[4][4] = {
		{0,4,8,7},
		{4,1,5,8},
		{8,5,2,6},
		{7,8,6,3}
	};
	
	static REAL buildt[4][9][3][2] = {
		/*S0*/
		{/*0*/{{0.,0.},{0.,0.},{-1.,-1.}},//por colunas
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{0.,0.},{0.,0.},{0.,0.}},
			/*4*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*5*/{{0.,.5},{0.,0.},{0.,-.5}},
			/*6*/{{-.5,0.},{0.,0.},{-.5,0.}},
			/*7*/{{.5,0.},{0.,0.},{.5,0.}},
			/*8*/{{.5,0.},{0.,.5},{-.5,-.5}}},
		/*S1*/
		{/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{1.,-1.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{0.,0.},{0.,0.},{0.,0.}},
			/*4*/{{.5,0.},{0.,0.},{.5,0.}},
			/*5*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*6*/{{-.5,0.},{0.,0.},{.5,0.}},
			/*7*/{{0.,-.5},{0.,0.},{0.,-.5}},
			/*8*/{{.5,0},{0.,.5},{.5,-.5}}},
		/*S2*/
		{/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{1.,1.}},
			/*3*/{{0.,0.},{0.,0.},{0.,0.}},
			/*4*/{{.5,0.},{0.,0.},{.5,0.}},
			/*5*/{{.5,0.},{0.,0.},{.5,0.}},
			/*6*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*7*/{{0.,-.5},{0.,0.},{0.,.5}},
			/*8*/{{.5,0.},{0.,.5},{.5,.5}}},
		/*S3*/
		{/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{0.,0.},{0.,0.},{-1.,1.}},
			/*4*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*5*/{{0.,.5},{0.,0.},{0.,.5}},
			/*6*/{{.5,0.},{0.,0.},{.5,0.}},
			/*7*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*8*/{{.5,0.},{0.,.5},{-.5,.5}}}
	};
	
	static int fatherside[TPZRefQuad::NSubEl][TPZShapeQuad::NSides] = {
		/*0*/ {0,4,8,7,4,8,8,7,8},
		/*1*/ {4,1,5,8,4,5,8,8,8},
		/*2*/ {8,5,2,6,8,5,6,8,8},
		/*3*/ {7,8,6,3,8,8,6,7,8}
	};
	
	
	//Into Divides is necesary to consider the connectivity with the all neighboards
	void TPZRefQuad::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
		int i;
		if(geo->HasSubElement()) {
			SubElVec.Resize(NSubEl);
			for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
			return;//If exist fSubEl return this sons
		}
		int j,sub,matid=geo->MaterialId();
		int64_t index;
		int np[TPZShapeQuad::NSides];//guarda conectividades dos 4 subelementos
		
		for(j=0;j<TPZShapeQuad::NCornerNodes;j++) np[j] = geo->NodeIndex(j);
		for(sub=TPZShapeQuad::NCornerNodes;sub<TPZShapeQuad::NSides;sub++) {
			NewMidSideNode(geo,sub,index);
			np[sub] = index;
		}
		// creating new subelements
		for(i=0;i<TPZShapeQuad::NCornerNodes;i++) {
			TPZManVector<int64_t>  cornerindexes(TPZShapeQuad::NCornerNodes);
			for(int j=0;j<TPZShapeQuad::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
			int64_t index;
			TPZGeoEl *subel = geo->Mesh()->CreateGeoElement(EQuadrilateral,cornerindexes,matid,index,0);
			geo->SetSubElement(i , subel);
		}
		
		SubElVec.Resize(NSubEl);
		for(sub=0;sub<NSubEl;sub++) {
			SubElVec[sub] = geo->SubElement(sub);
			SubElVec[sub]->SetFather(geo);
			SubElVec[sub]->SetFatherIndex(geo->Index());
		}
		for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
			for(j=0;j<NumInNeigh;j++) {        //lado do subel              numero do filho viz.             lado do viz.
				geo->SubElement(i)->SetNeighbour(InNeigh[i][j][0],TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
			}
		}
		
		geo->SetSubElementConnectivities();
	}
	
	void TPZRefQuad::NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index) {
		
		MidSideNodeIndex(gel,side,index);
		if(index < 0) {
			TPZGeoElSide gelside = gel->Neighbour(side);
			if(gelside.Element()) {
				while(gelside.Element() != gel) {
					gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
					if(index!=-1) return;
					gelside = gelside.Neighbour();
				}
			}
			TPZVec<REAL> par(pztopology::TPZQuadrilateral::Dimension,0.);
			TPZVec<REAL> coord(3,0.);
			if(side < TPZShapeQuad::NCornerNodes) {
				index = gel->NodeIndex(side); 
				return;
			}
			//aqui side = 8 a 26
			side-=TPZShapeQuad::NCornerNodes;//0,1,..,18
			par[0] = MidCoord[side][0];
			par[1] = MidCoord[side][1];
			gel->X(par,coord);
			index = gel->Mesh()->NodeVec().AllocateNewElement();
			gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
		}
	}
	
	void TPZRefQuad::MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index) {
		index = -1;
		if(side<0 || side>TPZShapeQuad::NSides-1) {
			PZError << "TPZRefQuad::MidSideNodeIndex. Bad parameter side = " << side << endl;
			return;
		}
		//sides 0 a 3
		if(side<TPZShapeQuad::NCornerNodes) {//o nó medio do lado 0 é o 0 etc.
			index = (gel)->NodeIndex(side);
			return; 
		}
		//o no medio da face eh o centro e o no medio do centro e o centro
		//como no de algum filho se este existir
		//caso tenha filhos eh o canto de algum filho, se nao tiver filhos retorna -1
		if(gel->HasSubElement()) {
			side-=TPZShapeQuad::NCornerNodes;
			index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
		}
	}
	
	void TPZRefQuad::GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel){
		
//		subel.Resize(0);
		if(side<0 || side>TPZShapeQuad::NSides || !father->HasSubElement()){
			PZError << "TPZRefQuad::GetSubelements called with error arguments\n";
			return;
		}
		int nsub = NSideSubElements(side);
		for(int i=0;i<nsub;i++)
			subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),subeldata[side][i][1]));
	}
	
	int TPZRefQuad::NSideSubElements(int side) {  
		if(side<0 || side>TPZShapeQuad::NSides-1){
			PZError << "TPZRefQuad::NSideSubelements called with error arguments\n";
			return -1;
		}
		return nsubeldata[side];
	}
	
	TPZTransform<> TPZRefQuad::GetTransform(int side,int whichsubel){
		
		if(side<0 || side>(TPZShapeQuad::NSides-1)){
			PZError << "TPZRefQuad::GetTransform side out of range or father null\n";
			return TPZTransform<>(0,0);
		}
		int smalldim = TPZShapeQuad::SideDimension(side);
		int fatherside = FatherSide(side,whichsubel);
		int largedim = TPZShapeQuad::SideDimension(fatherside);
		TPZTransform<> trans(largedim,smalldim);
		int i,j;
		for(i=0; i<largedim; i++) {
			for(j=0; j<smalldim; j++) {
				trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
			}
			trans.Sum() (i,0) = buildt[whichsubel][side][2][i];
		}
		return trans;
	}
	
	int TPZRefQuad::FatherSide(int side,int whichsubel){
		
		if(side<0 || side>TPZShapeQuad::NSides-1){
			PZError << "TPZRefQuad::Father2 called error" << endl;
			return -1;
		}
		return fatherside[whichsubel][side];
	}
        
        int TPZRefQuad::ClassId() const{
            return Hash("TPZRefQuad");
        }
	
};
