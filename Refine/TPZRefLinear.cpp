/**
 * @file
 * @brief Contains the implementation of the TPZRefLinear methods. 
 */
#include "TPZRefLinear.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZGeoElement.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

using namespace pzshape;
using namespace std;

namespace pzrefine {
	static int InNeigh[2][1][3] = {
		{{1,1,0}},{{0,0,1}}
	};
	
	static int CornerSons[2][2] = {
		{0,2},{2,1}
	};
	
	static int MidSideNodes[1][2]  = {
		{0,1}
	};
	
	static REAL MidCoord[1][1] = { {0.} };
	
	static int subeldata[3][3][2] =
	{
		{{0,0}},/*side=0 { isub0{0,1},isub1{0,1} }*/
		{{1,1}},/*side=1*/
       	{{0,2},{0,1},{1,2}}/*side=2*/
		//	{{0,1},{0,2},{1,2}}/*side=2*/
	};
	
	static int nsubeldata[3] = {1,1,3};
	
	static REAL buildt[2][3][2][1] = {//por colunas
		/*S0*/
		{/*0*/{{0.},{-1.}},
			/*1*/{{0.},{0.}},
			/*2*/{{.5},{-.5}}
		},
		/*S1*/
		{/*0*/{{0.},{0.}},
			/*1*/{{0.},{1.}},
			/*2*/{{.5},{.5}}
		}
	};
	
	static int fatherside[2][3] = {
		{0,2,2},{2,1,2}
	};
	
	
	void TPZRefLinear::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
		
		int i;
		if(geo->HasSubElement()) {
			SubElVec.Resize(NSubEl);
			for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
			return;//If exist fSubEl return this sons
		}
		int j,sub,matid=geo->MaterialId();
		int64_t index;
		int np[TPZShapeLinear::NSides];//guarda conectividades dos 8 subelementos
		
		for(j=0;j<TPZShapeLinear::NCornerNodes;j++) np[j] = geo->NodeIndex(j);
		for(sub=TPZShapeLinear::NCornerNodes;sub<TPZShapeLinear::NSides;sub++) {
			NewMidSideNode(geo,sub,index);
			np[sub] = index;
		}
		// creating new subelements
		for(i=0;i<TPZShapeLinear::NCornerNodes;i++) {
			TPZManVector<int64_t> cornerindexes(TPZShapeLinear::NCornerNodes);
			for(int j=0;j<TPZShapeLinear::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
			int64_t index;
			TPZGeoEl *subel = geo->Mesh()->CreateGeoElement(EOned,cornerindexes,matid,index);
			geo->SetSubElement(i , subel);
		}
		
		SubElVec.Resize(NSubEl);
		for(sub=0;sub<NSubEl;sub++) {
			SubElVec[sub] = geo->SubElement(sub);
			SubElVec[sub]->SetFather(geo);
			SubElVec[sub]->SetFatherIndex(geo->Index());
		}
		for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
			for(j=0;j<1;j++) {        //lado do subel                                          numero do filho viz.             lado do viz.
				geo->SubElement(i)->SetNeighbour(InNeigh[i][j][0],TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
			}
		}
		
		geo->SetSubElementConnectivities();
		
	}
	void TPZRefLinear::MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index){
		index = -1;
		if(side<0 || side>TPZShapeLinear::NSides-1) {
			PZError << "TPZRefCube::MidSideNodeIndex. Bad parameter side = " << side << endl;
			return;
		}
		//sides 0 a 7
		if(side<TPZShapeLinear::NCornerNodes) {//o nó medio do lado 0 é o 0 etc.
			index = (gel)->NodeIndex(side);
			return; 
		}
		//o nó medio da face é o centro e o nó medio do centro é o centro
		//como nó de algum filho se este existir
		//caso tenha filhos é o canto de algum filho, se não tiver filhos retorna -1
		if(gel->HasSubElement()) {
			side-=TPZShapeLinear::NCornerNodes;
			index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
		}
	}
	void TPZRefLinear::NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index) {
		
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
			TPZVec<REAL> par(pztopology::TPZLine::Dimension,0.);
			TPZVec<REAL> coord(3,0.);
			if(side < TPZShapeLinear::NCornerNodes) {
				index = gel->NodeIndex(side); 
				return;
			}
			//aqui side = 8 a 26
			side-=TPZShapeLinear::NCornerNodes;//0,1,..,18
			par[0] = MidCoord[side][0];
			gel->X(par,coord);
			index = gel->Mesh()->NodeVec().AllocateNewElement();
			gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
		}
		
	}
	
	void TPZRefLinear::GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel) {
//		subel.Resize(0);
		if(side<0 || side>TPZShapeLinear::NSides || !father->HasSubElement()){
			PZError << "TPZRefCube::GetSubelements2 called with error arguments\n";
			return;
		}
		int nsub = NSideSubElements(side);//nsubeldata[side];
		for(int i=0;i<nsub;i++)
			subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),subeldata[side][i][1]));
		
	}
	int TPZRefLinear::NSideSubElements(int side) {
		if(side<0 || side>TPZShapeLinear::NSides-1){
			PZError << "TPZRefCube::NSideSubelements2 called with error arguments\n";
			return -1;
		}
		return nsubeldata[side];
		
	}
	
	TPZTransform<> TPZRefLinear::GetTransform(int side,int whichsubel) {
		if(side<0 || side>TPZShapeLinear::NSides-1){
			PZError << "TPZRefLinear::GetTransform side out of range or father null\n";
			return TPZTransform<>(0,0);
		}
		int smalldim = TPZShapeLinear::SideDimension(side);
		int fatherside = FatherSide(side,whichsubel);
		int largedim = TPZShapeLinear::SideDimension(fatherside);
		TPZTransform<> trans(largedim,smalldim);
		int i,j;
		for(i=0; i<largedim; i++) {
			for(j=0; j<smalldim; j++) {
				trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
			}
			trans.Sum() (i,0) = buildt[whichsubel][side][1][i];
		}
		return trans;
		
	}
	
	int TPZRefLinear::FatherSide(int side,int whichsubel) {
		if(side<0 || side>TPZShapeLinear::NSides-1){
			PZError << "TPZRefCube::Father2 called error" << endl;
			return -1;
		}
		return fatherside[whichsubel][side];
		
	}
        
        int TPZRefLinear::ClassId() const{
            return Hash("TPZRefLinear");
        }
};
