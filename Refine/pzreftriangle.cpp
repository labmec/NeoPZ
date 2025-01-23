/**
 * @file
 * @brief Contains the implementation of the TPZRefTriangle methods. 
 */

#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapetriang.h"
#include "TPZGeoElement.h"
#include "pzgeoelside.h"
#include "pzgeoel.h"
#include "pzgmesh.h"

using namespace pzshape;
using namespace std;

namespace pzrefine {
	
	static int nsubeldata[7] = {1,1,1,3,3,3,7};
	
	static int subeldata[7][7][2] = {
		{{0,0}},/*side=0 {isub0{0,1},isub1{0,1},isub2{0,1},...}*/
		{{1,1}},/*side=1*/
		{{2,2}},/*side=2*/
		//  {{0,3},{0,1},{1,3}},/*side=3*/
		{{0,1},{0,3},{1,3}},/*side=3*/
		//  {{1,4},{1,2},{2,4}},/*side=4*/
		{{1,2},{1,4},{2,4}},/*side=4*/
		//  {{2,5},{2,0},{0,5}},/*side=5*/
		{{2,0},{2,5},{0,5}},/*side=5*/
		//{{0,6},{1,6},{2,6},{3,6},{2,3},{0,4},{1,5}}/*side=6*/
		{{2,3},{0,4},{1,5},{0,6},{1,6},{2,6},{3,6}}/*side=6*/
	};
	
	static int MidSideNodes[3][2] = {
		{0,1},
		{1,2},
		{2,0}
	};
	
	static REAL MidCoord[3][2] = { 
		{0.5,0.},  //side 3
		{0.5,0.5}, //side 4
		{0.,0.5}   //side 5
	};
	
	/**
	 * define as conectividades entre sub-elementos
	 * linha i eh filho i, {a,b,c} = {lado do filho atual,
	 * irmao vizinho,lado do vizinho}
	 */
	const int NumInNeigh = 6;
	static int InNeigh[4][NumInNeigh][3] = {
		{{1,1,0},{2,3,1},{4,3,4},{-1},{-1},{-1}},
		{{0,3,2},{2,2,1},{5,3,5},{-1},{-1},{-1}},
		{{0,0,2},{1,3,0},{3,3,3},{-1},{-1},{-1}},
		{{0,1,2},{1,2,0},{2,0,1},{3,2,3},{4,0,4},{5,1,5}}
	};
	
	/**
	 * define os cantos locais dos fihos
	 */
	static int CornerSons[4][3] = {
		{0,3,5},
		{3,1,4},
		{5,4,2},
		{4,5,3}
	};
	
	static REAL buildt[4][7][3][2] = {//por colunas
		/*S0*/ {
			/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*4*/{{-.25,.25},{0.,0.},{.25,.25}},
			/*5*/{{.5,0.},{0.,0.},{.5,0.}},
			/*6*/{{.5,0.},{0.,.5},{0.,0.}}},
		/*S1*/ {
			/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,1.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{.5,0.},{0.,0.},{.5,0.}},
			/*4*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*5*/{{0.,-.25},{0.,0.},{.5,.25}},
			/*6*/{{.5,0.},{0.,.5},{.5,0.}}},
		/*S2*/ {
			/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{0.,1.}},
			/*3*/{{.25,0.},{0.,0.},{.25,.5}},
			/*4*/{{.5,0.},{0.,0.},{.5,0.}},
			/*5*/{{.5,0.},{0.,0.},{-.5,0.}},
			/*6*/{{.5,0.},{0.,.5},{0.,.5}}},
		/*S3*/ {
			/*0*/{{0.,0.},{0.,0.},{0.,0.}},
			/*1*/{{0.,0.},{0.,0.},{0.,0.}},
			/*2*/{{0.,0.},{0.,0.},{0.,0.}},
			/*3*/{{-.25,0.},{0.,0.},{.25,.5}},
			/*4*/{{.25,-.25},{0.,0.},{.25,.25}},
			/*5*/{{0.,.25},{0.,0.},{.5,.25}},
			/*6*/{{-.5,0.},{0.,-.5},{.5,.5}}}
	};
	
	static int fatherside[TPZRefTriangle::NSubEl][TPZShapeTriang::NSides] = {
		/*0*/ {0,3,5,3,6,5,6},
		/*1*/ {3,1,4,3,4,6,6},
		/*2*/ {5,4,2,6,4,5,6},
		/*3*/ {4,5,3,6,6,6,6}
	};
	
	//Into Divides is necesary to consider the connectivity with the all neighboards
	void TPZRefTriangle::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
		int i;
		if(geo->HasSubElement()) {
			SubElVec.Resize(NSubEl);
			for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
			return;//If exist fSubEl return this sons
		}
		int j,sub,matid=geo->MaterialId();
		int64_t index;
		int np[TPZShapeTriang::NSides];//guarda conectividades dos 4 subelementos
		
		for(j=0;j<TPZShapeTriang::NCornerNodes;j++) np[j] = geo->NodeIndex(j);
		for(sub=TPZShapeTriang::NCornerNodes;sub<TPZShapeTriang::NSides-1;sub++) {
			NewMidSideNode(geo,sub,index);
			np[sub] = index;
		}
		// creating new subelements
		for(i=0;i<NSubEl;i++) {
			TPZManVector<int64_t> cornerindexes(TPZShapeTriang::NCornerNodes);
			for(int j=0;j<TPZShapeTriang::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
			int64_t index;
			TPZGeoEl *subel = geo->CreateGeoElement(ETriangle,cornerindexes,matid,index);
			geo->SetSubElement(i ,subel);
		}
		
		SubElVec.Resize(NSubEl);
		for(sub=0;sub<NSubEl;sub++) {
			SubElVec[sub] = geo->SubElement(sub);
			SubElVec[sub]->SetFather(geo);
			SubElVec[sub]->SetFatherIndex(geo->Index());
		}
		for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
			for(j=0;j<NumInNeigh;j++) {        //lado do subel numero do filho viz.             lado do viz.
				if(InNeigh[i][j][0]<0) continue;
				geo->SubElement(i)->SetNeighbour(InNeigh[i][j][0],TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
			}
		}
		geo->SetSubElementConnectivities();
	}
	
	void TPZRefTriangle::NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index) {
		
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
			TPZVec<REAL> par(pztopology::TPZTriangle::Dimension,0.);
			TPZVec<REAL> coord(3,0.);
			if(side < TPZShapeTriang::NCornerNodes) {
				index = gel->NodeIndex(side); 
				return;
			}
			//aqui side = 8 a 26
			side-=TPZShapeTriang::NCornerNodes;//0,1,..,18
			par[0] = MidCoord[side][0];
			par[1] = MidCoord[side][1];
			gel->X(par,coord);
			index = gel->Mesh()->NodeVec().AllocateNewElement();
			gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
		}
	}
	
	void TPZRefTriangle::MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index) {
		index = -1;
		if(side<0 || side>TPZShapeTriang::NSides-1) {
			PZError << "TPZRefTriangle::MidSideNodeIndex. Bad parameter side = " << side << endl;
			return;
		}
		//sides 0 a 3
		if(side<TPZShapeTriang::NCornerNodes) {//o n� medio do lado 0 � o 0 etc.
			index = (gel)->NodeIndex(side);
			return; 
		}
		//o n� medio da face � o centro e o n� medio do centro � o centro
		//como n� de algum filho se este existir
		//caso tenha filhos � o canto de algum filho, se n�o tiver filhos retorna -1
		if(gel->HasSubElement()) {
			side-=TPZShapeTriang::NCornerNodes;
			index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
		}
	}
	
	void TPZRefTriangle::GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel){
//		subel.Resize(0);
		if(side<0 || side>TPZShapeTriang::NSides || !father->HasSubElement()){
			PZError << "TPZRefTriangle::GetSubelements2 called with error arguments\n";
			return;
		}
		int nsub = NSideSubElements(side);//nsubeldata[side];
		for(int i=0;i<nsub;i++)
			subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),subeldata[side][i][1]));
	}
	
	int TPZRefTriangle::NSideSubElements(int side) {  
		if(side<0 || side>TPZShapeTriang::NSides-1){
			PZError << "TPZRefTriangle::NSideSubelements2 called with error arguments\n";
			return -1;
		}
		return nsubeldata[side];
	}
	
	TPZTransform<> TPZRefTriangle::GetTransform(int side,int whichsubel){
		if(side<0 || side>(TPZShapeTriang::NSides-1)){
			PZError << "TPZRefTriangle::GetTransform side out of range or father null\n";
			return TPZTransform<>(0,0);
		}
		int smalldim = TPZShapeTriang::SideDimension(side);
		int fatherside = FatherSide(side,whichsubel);
		int largedim = TPZShapeTriang::SideDimension(fatherside);
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
	
	int TPZRefTriangle::FatherSide(int side,int whichsubel){
		if(side<0 || side>TPZShapeTriang::NSides-1){
			PZError << "TPZRefTriangle::Father2 called error" << endl;
			return -1;
		}
		return fatherside[whichsubel][side];
	}
        
    int TPZRefTriangle::ClassId() const{
        return Hash("TPZRefTriangle");
    }	
};
