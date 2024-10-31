/**
 * @file
 * @brief Contains the implementation of the TPZGeoElSide and TPZGeoElSideIndex methods.
 */

#include "pzgeoelside.h"
#include "pzgeoel.h"
#include "pztrnsform.h"
#include "pzstack.h"
#include "pzvec_extras.h"
#include "pzquad.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzcompel.h"
#include "pzintel.h"
#include "pznumeric.h"

#include "pzmultiphysicscompel.h"

using namespace pzshape;
using namespace std;

#include "pzlog.h"

#include <sstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzgeoelside");
#endif

// Implementation of the TPZGeoElSideIndex methods

TPZGeoElSideIndex::TPZGeoElSideIndex(TPZGeoEl *gel,int side){  
    if (gel) this->fGeoElIndex = gel->Index();
    else this->fGeoElIndex = -1;
    this->fSide = side;
}

TPZGeoElSideIndex::TPZGeoElSideIndex(const TPZGeoElSide &side) {
    TPZGeoEl * gel = side.Element();
    if (gel) this->fGeoElIndex = gel->Index();
    else this->fGeoElIndex = -1;
    this->fSide = side.Side();
}

void TPZGeoElSideIndex::SetElement(TPZGeoEl* geoel){
    if (geoel) this->fGeoElIndex = geoel->Index();
    else this->fGeoElIndex = -1;
}

int TPZGeoElSideIndex::ClassId() const {
    return Hash("TPZGeoElSideIndex");
}

void TPZGeoElSideIndex::Read(TPZStream& buf, void* context) { //ok
    buf.Read(&fGeoElIndex);
    buf.Read(&fSide);
}

void TPZGeoElSideIndex::Write(TPZStream& buf, int withclassid) const { //ok
    buf.Write(&fGeoElIndex);
    buf.Write(&fSide);
}


// Implementation of the TPZGeoElSide methods

TPZGeoElSide::TPZGeoElSide(TPZGeoEl *gel, std::set<int64_t> &sideCornerNodes)
{
	fGeoEl = 0; fSide = -1;
	
	std::set<int64_t> nodes;
	for(int s = 0; s < gel->NSides(); s++)
	{
		nodes.clear();
		TPZGeoElSide actualSide(gel, s);
		if(actualSide.NSideNodes() != (int)sideCornerNodes.size())
		{
			continue;
		}
		for(int as = 0; as < actualSide.NSideNodes(); as++)
		{
			nodes.insert(actualSide.SideNodeIndex(as));
		}
		if(nodes == sideCornerNodes)
		{
			fGeoEl = gel;
			fSide = s;
			
			return;
		}
	}
}

TPZGeoElSide::TPZGeoElSide(TPZGeoEl *gel) : fGeoEl(gel), fSide(-1)
{
    if(fGeoEl) fSide = fGeoEl->NSides()-1;
}


bool TPZGeoElSide::IsAncestor(const TPZGeoElSide &other) const {
	if(*this == other) return true;
	TPZGeoElSide father = this->Father2();
	if(father.Element()){
		if(father.Element() == other.Element()){
			return true;
		}
		else{
			if(father.IsAncestor(other)){
				return true;
			}
		}
	}
	return false;
}

bool TPZGeoElSide::IsRelative(const TPZGeoElSide &other) const {
	if( this->IsAncestor(other) ) return true;
	if( other.IsAncestor(*this) ) return true;
	return false;
}



void TPZGeoElSide::X(TPZVec< REAL > &loc, TPZVec< REAL > &result) const {
	
	TPZManVector< REAL,3 > locElement(fGeoEl->Dimension(), 0.);
    result.Resize(3);
    QsiElement(loc, locElement);
	fGeoEl->X(locElement, result);
}

/** @brief X coordinate of a point loc of the side */
void TPZGeoElSide::QsiElement(TPZVec< REAL > &qsi_side, TPZVec< REAL > &qsi_element) const{
    
    qsi_element.Resize(fGeoEl->Dimension());
    TPZTransform<> ElementDim = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
    ElementDim.Apply(qsi_side, qsi_element);
}

/** @brief X coordinate of a point loc of the side */
void TPZGeoElSide::GradX(TPZVec<REAL> &loc, TPZFMatrix<REAL> &gradx) const{
    
#ifdef PZDEBUG
    if(!fGeoEl) return;
#endif
    
    int dim = fGeoEl->Dimension();
    TPZFNMatrix<9,REAL> gradx_vol(3,dim);
    
    TPZManVector< REAL,3 > locElement(dim, 0.);
    
    TPZTransform<> Transformation = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);

    Transformation.Apply(loc, locElement);
    TPZFNMatrix<9,REAL> trans_mult(Transformation.Mult().Cols(),Transformation.Mult().Rows());
    Transformation.Mult().Transpose(&trans_mult);
    fGeoEl->GradX(locElement, gradx_vol);
    gradx_vol.Multiply(Transformation.Mult(), gradx);
    
}

/** @brief GradX loc of the side */
void TPZGeoElSide::GradX(TPZVec< Fad<REAL> > &loc, TPZFMatrix< Fad<REAL> > &gradx) const{
    
    TPZManVector< Fad<REAL> ,3 > locElement(fGeoEl->Dimension(), 0.);
    gradx.Resize(3,Dimension());
    int dim = fGeoEl->Dimension();
    TPZFNMatrix<9,Fad<REAL>> gradx_vol(3,dim);

    TPZTransform<> ElementDimR = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
    TPZTransform<Fad<REAL> > ElementDim;
    ElementDim.CopyFrom(ElementDimR);
    ElementDim.Apply(loc, locElement);
    fGeoEl->GradX(locElement, gradx_vol);
    gradx_vol.Multiply(ElementDim.Mult(), gradx);
}


/** @brief X coordinate of a point loc of the side */
void TPZGeoElSide::X(TPZVec< Fad<REAL> > &loc, TPZVec< Fad<REAL> > &result) const
{
    TPZManVector<Fad<REAL>,3 > locElement(fGeoEl->Dimension(), 0.);
    result.Resize(3);
    
    TPZTransform<> ElementDimR = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
    TPZTransform<Fad<REAL> > ElementDim;
    ElementDim.CopyFrom(ElementDimR);
    
    ElementDim.Apply(loc, locElement);
    
    fGeoEl->X(locElement, result);

}





void TPZGeoElSide::Jacobian(TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const {
	
	if(!fGeoEl) return;
	int DIM = fGeoEl->Dimension();
	
	TPZTransform<> ThisTransf = fGeoEl->SideToSideTransform(fSide, fGeoEl->NSides()-1);
	TPZManVector< REAL,3 > paramElement(DIM,0.);
	
	ThisTransf.Apply(param,paramElement);
	REAL DetElement;
	
	TPZFNMatrix<9> JCn(3,DIM,0.), JacElement(DIM,DIM,0.), AxesElement(DIM,3,0.), Temp(3,DIM,0.), InvElement(DIM,DIM,0.), axest(3,DIM);
	
	fGeoEl->Jacobian(paramElement,JacElement,AxesElement,DetElement,InvElement);
	
	AxesElement.Transpose();
	AxesElement.Multiply(JacElement,JCn);
	JCn.Multiply(ThisTransf.Mult(),Temp);
	
	Temp.GramSchmidt(axest,jacobian); 
	axest.Transpose(&axes);
	jacinv.Resize(jacobian.Rows(),jacobian.Cols());
	if(axes.Rows() == 1) 
	{ 
		detjac = jacobian(0,0); 
		if(detjac)
		{
			jacinv(0,0) = 1./jacobian(0,0); 
		}
		else
		{
			jacinv(0,0) = 0.;
		}
	}
	if(axes.Rows() == 2) {
		detjac = jacobian(0,0)*jacobian(1,1) - jacobian(1,0)*jacobian(0,1);
		if(detjac)
		{
			jacinv(0,0) =  jacobian(1,1)/detjac; jacinv(0,1) = -jacobian(0,1)/detjac;
			jacinv(1,0) = -jacobian(1,0)/detjac; jacinv(1,1) =  jacobian(0,0)/detjac;
		}
		else
		{
			jacinv.Zero();
		}
	}
	if(axes.Rows() == 3) {
		detjac = -(jacobian(0,2)*jacobian(1,1)*jacobian(2,0)) + jacobian(0,1)*jacobian(1,2)*jacobian(2,0) + jacobian(0,2)*jacobian(1,0)*jacobian(2,1) - jacobian(0,0)*jacobian(1,2)*jacobian(2,1) - jacobian(0,1)*jacobian(1,0)*jacobian(2,2) + jacobian(0,0)*jacobian(1,1)*jacobian(2,2);
		jacinv(0,0) = (-(jacobian(1,2)*jacobian(2,1)) + jacobian(1,1)*jacobian(2,2))/detjac;
		jacinv(0,1) = (jacobian(0,2)*jacobian(2,1) - jacobian(0,1)*jacobian(2,2))/detjac;
		jacinv(0,2) = (-(jacobian(0,2)*jacobian(1,1)) + jacobian(0,1)*jacobian(1,2))/detjac;
		jacinv(1,0) = (jacobian(1,2)*jacobian(2,0) - jacobian(1,0)*jacobian(2,2))/detjac;
		jacinv(1,1) = (-(jacobian(0,2)*jacobian(2,0)) + jacobian(0,0)*jacobian(2,2))/detjac;
		jacinv(1,2) = (jacobian(0,2)*jacobian(1,0) - jacobian(0,0)*jacobian(1,2))/detjac;
		jacinv(2,0) = (-(jacobian(1,1)*jacobian(2,0)) + jacobian(1,0)*jacobian(2,1))/detjac;
		jacinv(2,1) = (jacobian(0,1)*jacobian(2,0) - jacobian(0,0)*jacobian(2,1))/detjac;
		jacinv(2,2) = (-(jacobian(0,1)*jacobian(1,0)) + jacobian(0,0)*jacobian(1,1))/detjac;
	}
}

/// Returns the number of sides in which the current side can be decomposed
int TPZGeoElSide::NSides() const
{
    TPZStack<int> lower;
    fGeoEl->LowerDimensionSides(fSide,lower);
    return lower.NElements()+1;
}


/// Area associated with the side
REAL TPZGeoElSide::Area() const
{
	TPZManVector<REAL,3> sideparam(Dimension(),0);
	REAL detjac;
	TPZFNMatrix<9> jacinv(3,3),jacobian(3,3),axes(3,3);
    //supondo jacobiano constante: X linear
	CenterPoint(sideparam);
	Jacobian(sideparam,jacobian,axes,detjac,jacinv);
    TPZIntPoints *intrule = fGeoEl->CreateSideIntegrationRule(fSide, 0);
    REAL RefElVolume = 0.;
    int np = intrule->NPoints();
    TPZManVector<REAL,3> points(Dimension());
    REAL weight;
    for (int ip=0; ip<np; ip++) {
        intrule->Point(ip, points, weight);
        RefElVolume += weight;
    }
	return (RefElVolume*detjac);//RefElVolume(): volume do elemento de refer�ncia
	
}


int TPZGeoElSide::NNeighbours() const
{
	int nneighbours = 0;
	TPZGeoElSide neigh = this->Neighbour();
	
	while(neigh.Element() != this->Element())
	{
		nneighbours++;
		neigh = neigh.Neighbour();
	}
	
	return nneighbours;
}

/** @brief Get number of neighbours of a given dimension */
int TPZGeoElSide::NNeighbours(const int dimfilter) const{
	if(dimfilter < 0) return NNeighbours();
	int nneighbours = 0;

	TPZGeoElSide neig = this->Neighbour();
	while(neig != *this){
		nneighbours += neig.Element()->Dimension() == dimfilter;
		neig = neig.Neighbour();
	}

	return nneighbours;
}

int TPZGeoElSide::NNeighbours(const int dimfilter, const std::set<int>& matids) const{
	if(matids.size() == 0) return NNeighbours(dimfilter);
	int nneighbours = 0;

	TPZGeoElSide neig = this->Neighbour();
	auto end = matids.end();
	while(neig != *this){
		nneighbours += (dimfilter<0 || neig.Element()->Dimension() == dimfilter) && matids.find(neig.Element()->MaterialId()) != end;
		neig = neig.Neighbour();
	}

	return nneighbours;
}

int TPZGeoElSide::NNeighboursButThisElem(TPZGeoEl *thisElem) const
{
	int nneighbours = 0;
	TPZGeoElSide neigh = this->Neighbour();
	
	while(neigh.Element() != this->Element())
	{
		if(neigh.Element() != thisElem)
		{
			nneighbours++;	
		}
		neigh = neigh.Neighbour();
	}
	
	return nneighbours;
}

void TPZGeoElSide::RemoveConnectivity(){
	
	if(!Exists()) return;
	if(fSide < 0 || fSide >= fGeoEl->NSides()) {
		PZError << "TPZGeoElSide::RemoveConnectivity Index out of bound\n";
	}
	//it removes the connectivity of the cycle where this inserted one: 
	//neighpre->this->neighpos => neighpre->neighpos
	TPZGeoElSide neighpre,neigh = Neighbour();
	if(neigh.Element() == NULL || neigh.Side() == -1){
		PZError << "TPZGeoElSide::RemoveConnectivity trying to remove null or inexistent connection";
	}
	TPZGeoElSide neighpos = neigh;
	while(neigh.Exists() && neigh != *this){
		neighpre = neigh;
		neigh = neigh.Neighbour();
	}
	if(neigh == *this){
		this->SetNeighbour(TPZGeoElSide());
		if (neighpre.Exists()) neighpre.SetNeighbour(neighpos);
	} else {
		PZError << "TPZGeoElSide::RemoveConnectivity neighbourhood cycle error";
	}
}

bool TPZGeoElSide::ResetBlendConnectivity(const int64_t &index){
    return fGeoEl->ResetBlendConnectivity(fSide, index);
}
using namespace std;

/**
 * This method inserts the element/side and all lowerdimension sides into the connectivity loop
 */
void TPZGeoElSide::InsertConnectivity(TPZGeoElSide &neighbour)
{
	if(!fGeoEl || !neighbour.fGeoEl) return;
	
	TPZStack<int> mylowerdimension, neighbourlowerdimension;
	fGeoEl->LowerDimensionSides(fSide,mylowerdimension);
	neighbour.fGeoEl->LowerDimensionSides(neighbour.fSide,neighbourlowerdimension);
	SetConnectivity(neighbour);
	int ns = mylowerdimension.NElements();
	int is;
	for(is=0; is<ns; is++)
	{
		TPZGeoElSide myl(fGeoEl,mylowerdimension[is]);
		TPZGeoElSide neighl(neighbour.fGeoEl,neighbourlowerdimension[is]);
#ifdef PZDEBUG
        std::set<int64_t> mynodes,neighnodes;
        int nsidenodes = myl.NSideNodes();
        for (int is = 0; is<nsidenodes; is++) {
            mynodes.insert(myl.SideNodeIndex(is));
            neighnodes.insert(neighl.SideNodeIndex(is));
        }
        if(mynodes != neighnodes) DebugStop();
#endif
		myl.SetConnectivity(neighl);
	}
}

/**
 * This method inserts the element/side and all lowerdimension sides into the connectivity loop
 */
void TPZGeoElSide::InsertConnectivity(TPZGeoElSide &neighbour, const TPZVec<int> &mapsides)
{
    if(!fGeoEl || !neighbour.fGeoEl) return;
    
    TPZStack<int> mylowerdimension, neighbourlowerdimension;
    fGeoEl->LowerDimensionSides(fSide,mylowerdimension);
    neighbour.fGeoEl->LowerDimensionSides(neighbour.fSide,neighbourlowerdimension);
    SetConnectivity(neighbour);
    int ns = mylowerdimension.NElements();

    for(int is=0; is<ns; is++)
    {
        TPZGeoElSide myl(fGeoEl,mylowerdimension[mapsides[is]]);
        TPZGeoElSide neighl(neighbour.fGeoEl,neighbourlowerdimension[is]);
#ifdef PZDEBUG
        std::set<int64_t> mynodes,neighnodes;
        int nsidenodes = myl.NSideNodes();
        for (int ss = 0; ss<nsidenodes; ss++) {
            mynodes.insert(myl.SideNodeIndex(ss));
            neighnodes.insert(neighl.SideNodeIndex(ss));
        }
        if(mynodes != neighnodes) DebugStop();
#endif
        myl.SetConnectivity(neighl);
    }
}

void TPZGeoElSide::SetConnectivity(const TPZGeoElSide &neighbour) const{
	
	if(!Exists()) return;
	if(fSide >= fGeoEl->NSides()) {
		PZError << "ERROR(TPZGeoElSide::SetConnectivity)-> Index greater than number of sides.\n";
		PZError << " fNumSides = " << fGeoEl->NSides() << " side = " << fSide << "\n";
	}
	if(!neighbour.Exists()) {
		fGeoEl->SetSideDefined(fSide);
	}
	
	TPZGeoElSide neighneigh, currentneigh;
	neighneigh = neighbour.Neighbour();
	currentneigh = Neighbour();
	if (!neighneigh.Exists() && !currentneigh.Exists()) {
		SetNeighbour(neighbour);
		neighbour.SetNeighbour(*this);
	} else if (neighneigh.Exists() && currentneigh.Exists()) {
		// It would be convenient to check the consistency of both loops
		// insert the connectivity between two independent loops
		int a = NeighbourExists(neighbour);
		int b = neighbour.NeighbourExists(*this);
		if(
		   (a && !b) ||
		   (!a && b)
		   )
		{
			std::cout << "This element side : " << fSide << std::endl;
			this->Element()->Print(std::cout);
			cout << "\nNeighbour side :"   << neighbour.Side() << std::endl;
			neighbour.Element()->Print(std::cout);
			PZError << "TPZGeoElSide::SetConnectivity Fourth untreated case, wrong data structure\n";
		} else if(!a) {
			SetNeighbour(neighneigh);
			neighbour.SetNeighbour(currentneigh);
		}
		/**The neighbouring element has already a connectivity loop, insert in his loop*/
	} else if (neighneigh.Exists()) {
		// the neighbour already has a loop, insert this into the loop
		SetNeighbour(neighneigh);
		neighbour.SetNeighbour(*this);
		/**The current element has already a connectivity loop and the neighbour doesnt*/
	} else if (currentneigh.Exists()) {
		// this is already inserted in a loop, insert neighbour in the loop
		SetNeighbour(neighbour);
		neighbour.SetNeighbour(currentneigh);
	}
}

void TPZGeoElSide::CenterPoint(TPZVec<REAL> &center) const
{
	if(!fGeoEl) return;
    TPZManVector<REAL,3> gelcenter(fGeoEl->Dimension());
	fGeoEl->CenterPoint(fSide,gelcenter);
    TPZTransform<> tr(Dimension(),fGeoEl->Dimension());
    tr = fGeoEl->SideToSideTransform(fGeoEl->NSides()-1,fSide);
    tr.Apply(gelcenter, center);
}

/** @brief return the coordinates of the center of the side in real space */
void TPZGeoElSide::CenterX(TPZVec<REAL> &Xcenter) const
{
	if(!fGeoEl) return;
    TPZManVector<REAL,3> gelcenter(fGeoEl->Dimension());
	fGeoEl->CenterPoint(fSide,gelcenter);
    fGeoEl->X(gelcenter, Xcenter);
}

void TPZGeoElSide::ComputeNeighbours(TPZStack<TPZGeoElSide> &compneigh) {
	if(fSide < fGeoEl->NCornerNodes())
    {
		AllNeighbours(compneigh);
		return;
    }
	int nsnodes = NSideNodes();
	TPZStack<TPZGeoElSide> GeoElSideSet;
	TPZStack<int> GeoElSet[27];
	int in;
	TPZManVector<int64_t> nodeindexes(nsnodes);
	for(in=0; in<nsnodes; in++)
	{
		nodeindexes[in] = SideNodeIndex(in);
		int locnod = fGeoEl->SideNodeLocIndex(fSide,in);
		GeoElSideSet.Resize(0);
		TPZGeoElSide locside(fGeoEl,locnod);
		locside.AllNeighbours(GeoElSideSet);
		int nel = GeoElSideSet.NElements();
		int el;
		for(el=0; el<nel; el++) {
			GeoElSet[in].Push(GeoElSideSet[el].Element()->Index());
		}
		Sort<int>(GeoElSet[in]);
	}
	TPZStack<int> result;
	switch(nsnodes) {
		case 1:
		{
			result = GeoElSet[0];
		}
			break;
		case 2:
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[1],result);
			break;
		case 3:
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[1],GeoElSet[2],result);
			break;
		case 4:
		{
			TPZStack<int> inter1, inter2;
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[0],GeoElSet[2],inter1);
			if(inter1.NElements()==0) break;
			Intersect<int,DEFAULTVEC_ALLOC>(GeoElSet[1],GeoElSet[3],inter2);
			if(inter2.NElements()==0) break;
			Intersect<int,DEFAULTVEC_ALLOC>(inter1,inter2,result);
		}
			break;
		default:
		{
			TPZStack<int> inter1, inter2;
			inter1 = GeoElSet[0];
			for(in=0; in<nsnodes-1; in++) {
				inter2.Resize(0);
				Intersect<int,DEFAULTVEC_ALLOC>(inter1,GeoElSet[in+1],inter2);
				if(inter2.NElements() == 0) break;
				inter1 = inter2;
			}
			result = inter2;
		}
	}
	int el,nel = result.NElements();
	TPZGeoMesh * geoMesh = fGeoEl->Mesh();
	for(el=0; el<nel; el++) {
		TPZGeoEl * gelResult = geoMesh->ElementVec()[result[el]];
		int whichSd = gelResult->WhichSide(nodeindexes);
		if(whichSd > 0)
		{
			compneigh.Push(TPZGeoElSide( gelResult, whichSd));
		}
	}
}


TPZTransform<> TPZGeoElSide::NeighbourSideTransform(const TPZGeoElSide &neighbour) {
	
#ifdef PZDEBUG
	if(!NeighbourExists(neighbour))
	{
		stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Neighbour does not exist : expect trouble";
		LOGPZ_ERROR(logger,sout.str());
		std::cout << "fSide " << fSide << " neighbour side" << neighbour.Side() << 
		" neighbour index " << neighbour.Element()->Index() << std::endl;
		std::cout << "fGeoEl\n";
		this->Element()->Print(std::cout);
		DebugStop();
		TPZTransform<> toto;
		return toto;
	}
#endif
	int sidedimension = Dimension();
	TPZTransform<> tside(sidedimension);//transforma�o local
	switch (sidedimension) {
		case 0://canto para canto viz
			
			break;
			
		case 1://aresta para aresta viz
			if (SideNodeIndex(0) == neighbour.SideNodeIndex(0))
			{
				tside.Mult()(0,0) =  1.;
			}
			else
			{
				tside.Mult()(0,0) = -1.;
			}
			break;
			
		case 2://transformacoes entre faces viz
			int i;
			//TPZCompEl *cel = Element()->Reference();
			TPZVec<int64_t> idto(0),idfrom(0);
			if(Element()->NSideNodes(Side()) == 4) {//faces quadrilaterais
				idto.Resize(4);
				idfrom.Resize(4);
				tside.Mult()(0,0) = 0;
				tside.Mult()(1,1) = 0;
				for(i=0;i<4;i++) idto[i]=neighbour.Element()->SideNodeIndex(neighbour.Side(),i);
				for(i=0;i<4;i++) idfrom[i]=Element()->SideNodeIndex(Side(),i);
				
				int transid = Element()->GetTransformId2dQ(idfrom,idto);
				tside.Mult()(0,0) = TPZShapeQuad::gTrans2dQ[transid][0][0];//cel->gTrans2dQ[transid][0][0];
				tside.Mult()(0,1) = TPZShapeQuad::gTrans2dQ[transid][0][1];
				tside.Mult()(1,0) = TPZShapeQuad::gTrans2dQ[transid][1][0];
				tside.Mult()(1,1) = TPZShapeQuad::gTrans2dQ[transid][1][1];
			} else if(Element()->NSideNodes(Side()) == 3) {//faces triangulares
				idto.Resize(3);
				idfrom.Resize(3);
				tside.Mult()(0,0) = 0.;
				tside.Mult()(1,1) = 0.;
				for(i=0;i<3;i++) idfrom[i] = Element()->SideNodeIndex(Side(),i);
				for(i=0;i<3;i++) idto[i] = neighbour.Element()->SideNodeIndex(neighbour.Side(),i);
				int transid = Element()->GetTransformId2dT(idfrom,idto);
				tside.Mult()(0,0) = TPZShapeTriang::gTrans2dT[transid][0][0];
				tside.Mult()(0,1) = TPZShapeTriang::gTrans2dT[transid][0][1];
				tside.Mult()(1,0) = TPZShapeTriang::gTrans2dT[transid][1][0];
				tside.Mult()(1,1) = TPZShapeTriang::gTrans2dT[transid][1][1];
				tside.Sum()(0,0) = TPZShapeTriang::gVet2dT[transid][0];
				tside.Sum()(1,0) = TPZShapeTriang::gVet2dT[transid][1];
			} else {
				PZError << "TPZGeoElSide::NeighbourSideTransform : elemento desconhecido" << std::endl;
			}
			break;
	}
	return tside;
}

int TPZGeoElSide::NeighbourExists(const TPZGeoElSide &gel) const {
	if(gel == *this) return 1;
	TPZGeoElSide neighbour = Neighbour();
	if(!neighbour.Exists()) return 0;
	while(neighbour != *this) {
		if(gel == neighbour) return 1;
		neighbour = neighbour.Neighbour();
	}
	return 0;
}

TPZCompElSide TPZGeoElSide::Reference() const {
	if (!fGeoEl) return TPZCompElSide();
	TPZCompElSide compside(fGeoEl->Reference(),fSide);
	return compside;
}

int TPZGeoElSide::Dimension() const {
    if (!fGeoEl) {
        PZError << "TPZGeoElSide::Dimension : null element\n";
        DebugStop();
    }
    return fGeoEl->SideDimension(fSide);
}

void TPZGeoElSide::SideTransform3(TPZGeoElSide neighbour,TPZTransform<> &t)	{
	//t : atual -> neighbour
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << __FUNCTION__ << " this = \n";
        Print(sout);
        sout << "neighbour\n";
        neighbour.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	TPZGeoElSide father(*this);
	if(!father.Exists()) {
		PZError << "TPZGeoElSide::SideTransform3 I dont understand\n";
		return;
	}
	while(father.Exists())
	{
		if(father.NeighbourExists(neighbour)) {
			TPZTransform<> Temp =  father.NeighbourSideTransform(neighbour);
			//       t =  NeighbourSideTransform(neighbour).Multiply(t);
			t =  Temp.Multiply(t);
			return;
		} else {
			TPZGeoElSide nextfather = father.StrictFather();
			if(nextfather.Exists())
			{
				t = father.Element()->BuildTransform2(father.Side(),nextfather.Element(),t);
			}
			father = nextfather;
		}
	}  
	TPZGeoElSide start,neighbourwithfather,neighfather;
	start = *this;
	neighbourwithfather = *this;
	do {
		neighfather = neighbourwithfather.StrictFather();
		if(!neighfather.Exists()) neighbourwithfather = neighbourwithfather.Neighbour();
	} while(!neighfather.Exists() && neighbourwithfather.Exists() && neighbourwithfather != start);
	int secondcase = 0;
	
	if(neighfather.Exists()) {
		if(neighbourwithfather != start) 
		{
			secondcase++;
			t = start.NeighbourSideTransform(neighbourwithfather).Multiply(t);
		}
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "neighbourwithfather\n";
            neighbourwithfather.Print(sout);
            sout << "neighbour\n";
            neighbour.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
		neighbourwithfather.SideTransform3(neighbour,t);
		return;
	}
    Print(std::cout);
    neighbour.Print(std::cout);
    fGeoEl->Print();
    neighbour.Element()->Print();
	PZError << "TPZGeoElSide:SideTranform3 did not find the neighbour\n";
	return;
}

void TPZGeoElSide::ConnectedCompElementList(TPZStack<TPZCompElSide> &ellist,
											int onlyinterpolated, int removeduplicates) {
	if(!fGeoEl) return;
	TPZCompElSide father = LowerLevelCompElementList2(onlyinterpolated);
	if(father.Exists()) ellist.Push(father);
	EqualLevelCompElementList(ellist,onlyinterpolated,removeduplicates);
	HigherLevelCompElementList2(ellist,onlyinterpolated,removeduplicates);
}

void TPZGeoElSide::EqualLevelCompElementList(TPZStack<TPZCompElSide> &elsidevec,
											 int onlyinterpolated, int removeduplicates) {
	
	TPZGeoElSide neighbour;
	TPZCompElSide ref;
	neighbour = Neighbour();
	if(!neighbour.Exists()) return;
	
	while(neighbour != *this) {
		ref = neighbour.Reference();
		if(ref.Element() && ref != Reference() && (!onlyinterpolated || dynamic_cast<TPZInterpolatedElement*>(ref.Element()) )) {
			elsidevec.Push(ref);
			if(removeduplicates) return;
		}
		neighbour = neighbour.Neighbour();
	}
}

void TPZGeoElSide::EqualLevelCompElementList3(TPZStack<TPZCompElSide> &elsidevec,
											 int onlymultiphysicelement, int removeduplicates) {
	
	TPZGeoElSide neighbour;
	TPZCompElSide ref;
	neighbour = Neighbour();
	if(!neighbour.Exists()) return;
	
	while(neighbour.Element() != this->Element()) {
		ref = neighbour.Reference();
		if(ref.Element() && ref.Element() != Reference().Element() && (!onlymultiphysicelement || dynamic_cast<TPZMultiphysicsElement *>(ref.Element()) )) {
			elsidevec.Push(ref);
			if(removeduplicates) return;
		}
		neighbour = neighbour.Neighbour();
	}
    if(neighbour.Element() == this->Element())
    {
        ref = neighbour.Reference();
        if(ref.Exists()) elsidevec.Push(ref);
    }
}


void TPZGeoElSide::HigherDimensionElementList(TPZStack<TPZCompElSide> &elsidevec, int onlyinterpolated) {
	
	TPZStack<TPZGeoElSide> gelsides;
	fGeoEl->AllHigherDimensionSides(fSide,2,gelsides);
	int il,cap = gelsides.NElements();
	for(il=0; il<cap; il++) {
		TPZCompElSide cels = gelsides[il].Reference();
		if(onlyinterpolated) {
			TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *> (cels.Element());
			if(!cel) continue;
			//int64_t conind = cel->ConnectIndex(cels.Side());
            int locconind = cel->MidSideConnectLocId(cels.Side());
            int64_t conind = cel->ConnectIndex(locconind);
			if(conind < 0) continue;
		}
		elsidevec.Push(cels);
	}
	
}

int64_t TPZGeoElSide::Id() {
	return fGeoEl->Id();
}

/**fill in the data structure for the neighbouring information*/
void TPZGeoElSide::SetNeighbour(const TPZGeoElSide &neighbour) const {
    if(!fGeoEl) {
		PZError << "TPZGeoElSide::SetNeighbour. Don't exist geometrical element.\n";
		return;
    }
    fGeoEl->SetNeighbour(fSide,neighbour);
}


TPZTransform<> TPZGeoElSide::SideToSideTransform(TPZGeoElSide &higherdimensionside) {
    if(fGeoEl != higherdimensionside.fGeoEl) {
		PZError << "TPZGeoElSide::SideToSideTransform inconsistent id1 = " << fGeoEl->Id() << 
		" id2 = " << higherdimensionside.fGeoEl->Id() << std::endl;
    }
    return fGeoEl->SideToSideTransform(fSide,higherdimensionside.fSide);
}

TPZGeoElSide TPZGeoElSide::LowestFatherSide()
{
    TPZGeoEl * actGel = this->fGeoEl;
    TPZGeoElSide side(*this);
    while(actGel->Father())
    {
        side = actGel->Father2(side.Side());
        actGel = actGel->Father();
    }
    return side;
}

/// return the TPZGeoElSide element that contains the current element/side
TPZGeoElSide TPZGeoElSide::LowerLevelSide() const
{
    TPZGeoElSide lower = StrictFather();
    if(lower) return lower;
    TPZGeoElSide neighbour = Neighbour();
    while(neighbour != *this)
    {
        lower = neighbour.StrictFather();
        if(lower) return lower;
        neighbour = neighbour.Neighbour();
    }
    return lower;
}



void TPZGeoElSide::GetAllSiblings(TPZStack<TPZGeoElSide> &sonSides)
{
#ifdef PZDEBUG
	PZError << __PRETTY_FUNCTION__ << "is deprecated. Use TPZGeoEl::YoungestChildren instead \n";
	DebugStop();
#endif // PZDEBUG
    if(this->Element()->HasSubElement() == false)
    {
        sonSides.Push(*this);
    }
    int dim = Dimension();
    TPZStack<TPZGeoElSide> lowerSubelements;
    fGeoEl->GetSubElements2(fSide,lowerSubelements,dim);
    int nsub = lowerSubelements.size();
    for (int s=0; s<nsub; s++) {
        lowerSubelements[s].GetAllSiblings(sonSides);
    }
}

void TPZGeoElSide::YoungestChildren(TPZStack<TPZGeoElSide> &sonSides)
{
    if(this->Element()->HasSubElement() == false)
    {
        sonSides.Push(*this);
    }
    int dim = Dimension();
    TPZStack<TPZGeoElSide> lowerSubelements;
    fGeoEl->GetSubElements2(fSide,lowerSubelements,dim);
    int nsub = lowerSubelements.size();
    for (int s=0; s<nsub; s++) {
        lowerSubelements[s].YoungestChildren(sonSides);
    }
}

/*return 1 if the element has subelements along side*/
int TPZGeoElSide::HasSubElement() {
    if(!fGeoEl) return 0;
    return fGeoEl->HasSubElement();
}
/*return the number of nodes for a particular side*/
int TPZGeoElSide::NSideNodes() const {
    if(!fGeoEl) return 0;
    return fGeoEl->NSideNodes(fSide);
}

/**returns the index of the nodenum node of side*/
int64_t TPZGeoElSide::SideNodeIndex(int nodenum) const {
    if(!fGeoEl) return -1;
    return ( fGeoEl->SideNodeIndex(fSide,nodenum) );
}

/**returns the index of the local nodenum  node of side*/
int64_t TPZGeoElSide::SideNodeLocIndex(int nodenum) const {
    if(!fGeoEl) return -1;
    return ( fGeoEl->SideNodeLocIndex(fSide,nodenum) );
}


TPZCompElSide TPZGeoElSide::LowerLevelCompElementList2(int onlyinterpolated)
{
	// This method was modified to look for the father of any neighbouring element
	// It is not sufficient to look for the father of the current element only, because a neighbour
	// might have a father where the current element doesn t. This happens in the clone meshes. It probably
	// will happen when working with interface elements or any situation where an element is inserted in an already refined mesh
	TPZGeoElSide father,neighbour,start;
	start = *this;
	neighbour =  Neighbour();
	father = StrictFather();
	int secondcase = 0;
	while(!father.Exists() && neighbour.Exists() && neighbour != start) {
		father = neighbour.StrictFather();
		neighbour = neighbour.Neighbour();
		secondcase++;
	}
	if(!father.Exists()) return TPZCompElSide();
	TPZStack<TPZCompElSide> equal;
	//2003-12-03
	if (father.Reference().Exists())  equal.Push(father.Reference());
	
	father.EqualLevelCompElementList(equal,onlyinterpolated,1);
	
	
	while(father.Exists() && equal.NElements() == 0) {
		neighbour = father.Neighbour();
		start = father;
		father = father.StrictFather();
		while(!father.Exists() && neighbour.Exists() && neighbour != start) {
			secondcase++;
			father = neighbour.StrictFather();
			neighbour = neighbour.Neighbour();
		}
		father.EqualLevelCompElementList(equal,onlyinterpolated,1);
		//2003-12-03
		if (father.Reference().Exists())  equal.Push(father.Reference());
		
	}
	if(equal.NElements()) return equal[0];
	return TPZCompElSide();
}

TPZGeoElSide TPZGeoElSide::Father2() const
{
	if(!fGeoEl) return TPZGeoElSide();
	return fGeoEl->Father2(fSide);
}

TPZGeoElSide TPZGeoElSide::StrictFather() const
{
	TPZGeoElSide father = Father2();
	int nfathsub = 0;
	if(father.Exists()) nfathsub = father.fGeoEl->NSideSubElements(father.fSide);
	while(father.Exists() && nfathsub == 1) {
		father = father.Father2();
		if(father.Exists()) nfathsub = father.fGeoEl->NSideSubElements(father.fSide);
	}
	return father;
}

void TPZGeoElSide::GetSubElements2(TPZStack<TPZGeoElSide> &subelements)
{
	if(!fGeoEl || !fGeoEl->HasSubElement()) {   // Jorge 10/01/2000
// comentei para poder acumular subelementos
//		subelements.Resize(0);
		return;
	}
	fGeoEl->GetSubElements2(fSide,subelements);
}

void TPZGeoElSide::HigherLevelCompElementList2(TPZStack<TPZCompElSide> &elvec, int onlyinterpolated, int removeduplicates) {
	
	if(!Dimension()) return;
	TPZGeoElSide neighbour(*this);
	TPZStack<TPZGeoElSide> subel;
	do {
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
			neighbour.GetSubElements2(subel);
			int nsub = subel.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				subel[is].EqualorHigherCompElementList2(elvec,onlyinterpolated,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);
}

void TPZGeoElSide::HigherLevelCompElementList3(TPZStack<TPZCompElSide> &elvec, int onlymultiphysicelement, int removeduplicates) {
	
	if(!Dimension()) return;
	TPZGeoElSide neighbour(*this);
	TPZStack<TPZGeoElSide> subel;
	do {
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
			neighbour.GetSubElements2(subel);
			int nsub = subel.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				subel[is].EqualorHigherCompElementList3(elvec,onlymultiphysicelement,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);
}


void TPZGeoElSide::EqualorHigherCompElementList2(TPZStack<TPZCompElSide> &celside, int onlyinterpolated, int removeduplicates){
	
	
    int ncelsides = celside.NElements();
	if(Reference().Exists()) {
		celside.Push(Reference());
		if(removeduplicates) {
			return;
		}
	}
	this->EqualLevelCompElementList(celside,onlyinterpolated,removeduplicates);
	if(ncelsides != celside.NElements()) return;
	TPZStack<TPZGeoElSide> gelsides;
	TPZGeoElSide neighbour(*this);
	do {
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
			neighbour.GetSubElements2(gelsides);
			int nsub = gelsides.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				gelsides[is].EqualorHigherCompElementList2(celside,onlyinterpolated,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);
	
}

void TPZGeoElSide::EqualorHigherCompElementList3(TPZStack<TPZCompElSide> &celside, int onlymultiphysicelement, int removeduplicates){
	
	
	int ncelsides = celside.NElements();
	if(Reference().Exists()) {
		celside.Push(Reference());
		if(removeduplicates) {
			return;
		}
	}
	this->EqualLevelCompElementList3(celside,onlymultiphysicelement,removeduplicates);
	if(ncelsides != celside.NElements()) return;
	TPZStack<TPZGeoElSide> gelsides;
	TPZGeoElSide neighbour(*this);
	do {
		if(neighbour.HasSubElement() && neighbour.NSubElements() > 1) {
			neighbour.GetSubElements2(gelsides);
			int nsub = gelsides.NElements();
			int is;
			for(is=0; is<nsub; is++) {
				gelsides[is].EqualorHigherCompElementList3(celside,onlymultiphysicelement,removeduplicates);
			}
			neighbour = *this;
		} else {
			neighbour = neighbour.Neighbour();
		}
		if(!neighbour.Exists()) break;
	} while(neighbour != *this);
	
}



int TPZGeoElSide::NSubElements()
{
	
	if(!Exists()) return -1;
	return fGeoEl->NSideSubElements(fSide);
}

void TPZGeoElSide::BuildConnectivities(TPZVec<TPZGeoElSide> &sidevec,TPZVec<TPZGeoElSide> &neighvec){
	/**
	 os vetores trazem a partic� do lado comum a  
	 dois vizinhos segundo os seus proprios padr�s de
	 refinamento, a divis� �identica para este lado comum*/ //cout << "Sao iguais: acertar as vizinhancas!!!\n";
	static int count = 0;

	int64_t size = sidevec.NElements();
	int64_t neighsize = neighvec.NElements();
	if(size!=neighsize || !size){
		PZError << "TPZGeoElSide::BuildConnectivities wrong vectors: abort!!!\n";
		DebugStop();
	}
	int64_t iv,ivn,side,neighside,sidedim,neighsidedim;
	TPZGeoElSide subside,neighsubside;
	for(iv=0;iv<size;iv++){
		subside = sidevec[iv];
		side = subside.Side();
		sidedim = subside.Dimension();
		TPZGeoEl *elside = subside.Element();
		for(ivn=0;ivn<neighsize;ivn++){
			neighsubside = neighvec[ivn];
			neighside = neighsubside.Side();
			neighsidedim = neighsubside.Dimension();
			TPZGeoEl *elneigh = neighsubside.Element();
			if(neighsidedim != sidedim) continue;
			//if(temp.Neighbour().Element()) continue;//?????? MELHORAR ESTA LINHA: pode ser viz. do irm�
			int in[4],face,i,j,num;
			int im[4],neighface;      
			switch(sidedim){
				case 0://canto	    
					if(elside->SideNodeIndex(side,0) == elneigh->SideNodeIndex(neighside,0)){
						if(subside.NeighbourExists(neighsubside)) {
							if(count < 5) cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?\n";
							count++;
						} else {
							subside.SetConnectivity(neighsubside);
						}
					}
					break;
				case 1://aresta
					in[0] = elside->SideNodeIndex(side,0);
					in[1] = elside->SideNodeIndex(side,1);
					im[0] = elneigh->SideNodeIndex(neighside,0);
					im[1] = elneigh->SideNodeIndex(neighside,1);
					if( (in[0] == im[0] && in[1] == im[1]) || (in[0] == im[1] && in[1] == im[0]) ){
						if(subside.NeighbourExists(neighsubside)) {
							if(count < 5) cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?\n";
							count++;
						} else {
							subside.SetConnectivity(neighsubside);
						}
					}
					break;
				case 2://face
					//face = NSideNodes(side);//original substituida pela seguinte
					face = elside->NSideNodes(side);
					neighface = elneigh->NSideNodes(neighside);
					if(face!=neighface) break;
					for(i=0;i<face;i++){
						in[i] = elside->SideNodeIndex(side,i);
						im[i] = elneigh->SideNodeIndex(neighside,i);
					}
					if(face==3) in[3] = im[3] = -1;
					num = 0;
					for(i=0;i<4;i++) for(j=0;j<4;j++) if(in[i]==im[j]) num++;
					if(num==4){
						if(subside.NeighbourExists(neighsubside)) {
							if(count < 5) cout << "TPZGeoElSide::BuildConnectivities the neighbour already exists?\n";
							count++;
						} else {
							subside.SetConnectivity(neighsubside);
						}
					}
					break;
				default:
					PZError << "TPZGeoElSide::BuildConnectivities error !!!\n";
					
			}
		}
	}
}

std::ostream &operator << (std::ostream & out,const TPZGeoElSide &geoside){
	out << "TPZGeoElSide : side = " << geoside.Side() << std::endl ;
	if (geoside.Element()) geoside.Element()->Print(out);
	return out;
}

bool TPZGeoElSide::IsLinearMapping() const
{
	if(!fGeoEl) return false;
	return fGeoEl->IsLinearMapping();
}


/** @brief compute the normal to the point from left to right neighbour */
void TPZGeoElSide::Normal(TPZVec<REAL> &point, TPZGeoEl *LeftEl, TPZGeoEl *RightEl, TPZVec<REAL> &normal) const
{
    normal.Resize(3,0.);
	normal.Fill(0.);
	int faceleft,faceright;
    int Leftdim = LeftEl->Dimension();
    int Rightdim = RightEl->Dimension();
    
	TPZManVector<REAL, 3> centleft(Leftdim),centright(Rightdim),result(3,0.),xint(3),xvolleft(3),xvolright(3),vec(3),rib(3);
	REAL normalize;
	int i;
	
	int InterfaceDimension = Dimension();
    TPZFNMatrix<9,REAL> axes(InterfaceDimension,3), jacobian(InterfaceDimension,InterfaceDimension),invjacobian(InterfaceDimension,InterfaceDimension);
    REAL detjac;
    this->Jacobian(point,jacobian,axes,detjac,invjacobian);
	faceleft = LeftEl->NSides()-1;//lado interior do elemento esquerdo
	faceright = RightEl->NSides()-1; // lado interior do element direito
	LeftEl->CenterPoint(faceleft,centleft);//ponto centro do elemento de volume
	RightEl->CenterPoint(faceright,centright);
	LeftEl->X(centleft,xvolleft);
	RightEl->X(centright,xvolright);
	for(i=0;i<3;i++) vec[i] = xvolright[i]-xvolleft[i];//nao deve ser nulo
	
	
	REAL vecnorm = sdot(vec, vec);
	if(vecnorm < 1.e-10)
	{
		LOGPZ_ERROR(logger,"Left and Right element centers coincide")
        vec[0]=1.;
	}
	
	
	switch(InterfaceDimension){
		case 0:
			normal[0] = vec[0];// a normal sempre aponta direcao positiva do eixo
			normal[1] = vec[1];
			normal[2] = vec[2];
            
            normalize = 0.;
			for(i=0;i<3;i++) normalize += normal[i]*normal[i];
            normalize = sqrt(normalize);
            if(!IsZero(normalize))
                for(i=0;i<3;i++) normal[i] = normal[i]/normalize;
            
			break;
		case 1:
			for(i=0;i<3;i++) rib[i] = axes(0,i);//direcao da aresta
            TPZNumeric::ProdVetorial(rib, vec, result);
            TPZNumeric::ProdVetorial(result, rib, normal);
			//normalizando a normal
			normalize = 0.;
			for(i=0;i<3;i++) normalize += normal[i]*normal[i];
			if(normalize == 0.0)
			{
				PZError << __PRETTY_FUNCTION__ << " null normal vetor\n";
#ifdef PZ_LOG
				{
					std::stringstream sout;
					Print(sout);
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
				
				DebugStop();
			}
			normalize = sqrt(normalize);
			for(i=0;i<3;i++) normal[i] = normal[i]/normalize;
			break;
		case 2:{
			TPZManVector<REAL,3> axes1(3), axes2(3);
			for(int iax = 0; iax < 3; iax++){
				axes1[iax] = axes(0,iax);
				axes2[iax] = axes(1,iax);
			}
            TPZNumeric::ProdVetorial(axes1,axes2,normal);
		}
			break;
		default:
			PZError << "TPZInterfaceElement::NormalToFace in case that not treated\n";
			normal.Resize(0);
			DebugStop();
			return;
	}
	
	//to guarantee the normal points from left to right neighbours:
	REAL dot = 0.;
	for(i=0; i<3; i++) dot += normal[i]*vec[i];
	if(dot < 0.) {
		for(i=0; i<3; i++) normal[i] = -normal[i];
	}

}

/** @brief compute the normal to the point */
void  TPZGeoElSide::Normal(TPZVec<REAL> &qsi_side, TPZVec<REAL> &normal) const{
    
    if(Dimension() == 2 && fGeoEl->Dimension() == 2) {
        TPZFNMatrix<9,REAL> gradx(3,2);
        fGeoEl->GradX(qsi_side, gradx);
        TPZFNMatrix<6,REAL> jac(2,2),jacinv(2,2),axes(2,3);
        REAL detjac;
        fGeoEl->Jacobian(gradx, jac, axes, detjac, jacinv);
        for(int i=0; i<3; i++) normal[i] = axes(0,(i+1)%3)*axes(1,(i+2)%3)-axes(0,(i+2)%3)*axes(1,(i+1)%3);
        return;
    }
    if (Dimension() != fGeoEl->Dimension()-1) {
        DebugStop();
    }
    
    normal.Resize(3,0.);
    normal.Fill(0.);

    TPZManVector<REAL, 3> center_dir(3,0.0), vol_center_x(3,0.0), x(3,0.0);
    
    TPZGeoElSide vol_side(fGeoEl,fGeoEl->NSides()-1);
    vol_side.CenterX(vol_center_x);
    this->X(qsi_side, x);
    center_dir = x - vol_center_x;
    
    int vol_dim   = vol_side.Dimension();
    int side_dim = Dimension();
    
    TPZManVector<REAL, 3> qsi_vol(vol_dim,0.0);
    this->QsiElement(qsi_side, qsi_vol);
    
    TPZFMatrix<REAL> side_jac(side_dim,side_dim), side_inv_jac(side_dim,side_dim), side_axes(side_dim,3);
    TPZFMatrix<REAL> vol_jac(vol_dim,vol_dim), vol_inv_jac(vol_dim,vol_dim), vol_axes(vol_dim,3);
    REAL detjac;
    
    this->Jacobian(qsi_side, side_jac, side_axes, detjac, side_inv_jac);
    vol_side.Jacobian(qsi_vol, vol_jac, vol_axes, detjac, vol_inv_jac);
    
    switch (side_dim) {
        case 0:
        {
            for (unsigned int i = 0; i < 3; i++) {
                normal[i] = vol_axes(i,0);
            }
        }
            break;
        case 1:
        {
            TPZManVector<REAL, 3> v1(3,0.0),v2(3,0.0),v3(3,0.0),v4(3,0.0);
            for (unsigned int i = 0; i < 3 ; i++) {
                v1[i] = vol_axes(0,i);
                v2[i] = vol_axes(1,i);
                v3[i] = side_axes(0,i);
            }
            
            TPZNumeric::ProdVetorial(v1, v2, v4);
            TPZNumeric::ProdVetorial(v4, v3, normal);
            
        }
            break;
        case 2:
        {
            TPZManVector<REAL, 3> v1(3,0.0),v2(3,0.0);
            for (unsigned int i = 0; i < 3 ; i++) {
                v1[i] = side_axes(0,i);
                v2[i] = side_axes(1,i);
            }
            TPZNumeric::ProdVetorial(v1, v2, normal);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
    TPZNumeric::NormalizeVetor(normal);
    
    //to guarantee the normal points from left to right neighbours:
    REAL dot = 0.;
    for(unsigned int i=0; i<3; i++) dot += normal[i]*center_dir[i];
    if(dot < 0.) {
        for(unsigned int i=0; i<3; i++) normal[i] = -normal[i];
    }
    
}

/** @brief print geometric characteristics of the element/side */
void TPZGeoElSide::Print(std::ostream &out) const
{
    if(! fGeoEl)
    {
        out << "Null TPZGeoElSide\n";
        return;
    }
    out << "Element index " << fGeoEl->Index() << " Side " << fSide << " SideNode indexes " ;
    TPZManVector<REAL,3> center(Dimension(),0.), centerX(3,0.);
    CenterPoint(center);
    X(center, centerX);
    for (int i=0; i<NSideNodes(); i++) {
        out << SideNodeIndex(i) << " ";
    }
    out << "Center coordinate " << centerX << std::endl;
}

TPZIntPoints *TPZGeoElSide::CreateIntegrationRule(int order) {
    return this->Element()->CreateSideIntegrationRule(this->Side(), order);
}

int TPZGeoElSide::GelLocIndex(int index) const
{
    if (!fGeoEl) {
        DebugStop();
    }
    return fGeoEl->SideNodeLocIndex(fSide,index);
}

int TPZGeoElSide::ClassId() const {
    return Hash("TPZGeoElSide");
}

void TPZGeoElSide::Read(TPZStream& buf, void* context) { //ok
    fGeoEl = dynamic_cast<TPZGeoEl*>(TPZPersistenceManager::GetInstance(&buf));
    buf.Read(&fSide);
}

void TPZGeoElSide::Write(TPZStream& buf, int withclassid) const { //ok
    TPZPersistenceManager::WritePointer(fGeoEl, &buf);
    buf.Write(&fSide);
}

TPZGeoEl *TPZGeoElSideIndex::Element(const TPZGeoMesh *mesh) const{
    if (this->fSide == -1 || this->fGeoElIndex == -1){
		return NULL;
    }
    return mesh->ElementVec()[this->fGeoElIndex];
}

TPZGeoElSide TPZGeoElSide::HasNeighbour(int materialid) const
{
    if(!fGeoEl) return TPZGeoElSide();
    if(fGeoEl->MaterialId() == materialid)
    {
        return (*this);
    }
    TPZGeoElSide neighbour = Neighbour();
    while(neighbour != *this)
    {
		int matid = neighbour.Element()->MaterialId();
        if(matid == materialid) return neighbour;
        neighbour = neighbour.Neighbour();
    }
    return TPZGeoElSide();
}

TPZGeoElSide TPZGeoElSide::HasNeighbour(const std::set<int> &matIDs) const
{
    for (const auto &it : matIDs) {
        TPZGeoElSide neigh = HasNeighbour(it);
        if (neigh.Exists())
        {
            return neigh;
        }
    }
    return TPZGeoElSide();
}

/** verifiy if a larger (lower level) neighbour exists with the given material id
 */
TPZGeoElSide TPZGeoElSide::HasLowerLevelNeighbour(int materialid) const
{
    if(!fGeoEl) return TPZGeoElSide();
    TPZGeoElSide lower = LowerLevelSide();
    while(lower)
    {
        TPZGeoElSide neighbour = lower.HasNeighbour(materialid);
        if(neighbour) return neighbour;
        lower = lower.LowerLevelSide();
    }
    return lower;
}

/** verifiy if a larger (lower level) neighbour exists with the given material id
 */
TPZGeoElSide TPZGeoElSide::HasLowerLevelNeighbour(const std::set<int> &matids) const
{
    if(!fGeoEl) return TPZGeoElSide();
    TPZGeoElSide lower = LowerLevelSide();
    while(lower)
    {
        TPZGeoElSide neighbour = lower.HasNeighbour(matids);
        if(neighbour) return neighbour;
        lower = lower.LowerLevelSide();
    }
    return lower;
}
/// Compute if the neighbour along a face has normal pointing outward
bool TPZGeoElSide::IsNeighbourCounterClockWise(TPZGeoElSide &neighbour)
{
#ifdef PZDEBUG
    // the side must have dimension dim-1
    if(!fGeoEl) DebugStop();
    if(Dimension() != fGeoEl->Dimension()-1) DebugStop();
    if(!IsNeighbour(neighbour)) DebugStop();
    if(neighbour.Side() != neighbour.Element()->NSides()-1) DebugStop();
#endif
    int nnodes = fGeoEl->NSideNodes(fSide);
    TPZManVector<int64_t> thisids(nnodes),neighids(nnodes);
    for(int in = 0; in<nnodes; in++) {
        thisids[in] = fGeoEl->SideNodeIndex(fSide, in);
        neighids[in] = neighbour.Element()->SideNodeIndex(neighbour.Side(), in);
    }
    int trid = 0;
    if(nnodes == 2)
    {
        trid = thisids[0] == neighids[0] ? 0 : 1;
    } else if (nnodes == 3) {
        trid = Element()->GetTransformId2dT(thisids,neighids) %2;
    } else if (nnodes == 4) {
        trid = Element()->GetTransformId2dQ(thisids,neighids) %2;
    } else {
        DebugStop();
    }
    int face = fSide-fGeoEl->FirstSide(2);
    int faceorient = fGeoEl->GetFaceOrientation(face);
    if(faceorient == -1) trid = (trid+1)%2;
    return trid == 0;
}
