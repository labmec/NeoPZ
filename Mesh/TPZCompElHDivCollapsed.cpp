/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDivCollapsed methods.
 */


#include "TPZCompElHDivCollapsed.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "pzlog.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzelchdiv.h"
#include "TPZShapeHDivCollapsed.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDivCollapsed");
#else
static int logger;
#endif

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,gel,index), fBottom(mesh,gel,index), fTop(mesh,gel,index)
{
    index = this->fIndex;
    mesh.ElementVec().SetFree(fTop.Index());
    mesh.ElementVec().SetFree(fBottom.Index());
    int64_t bottom_c_index = fBottom.ConnectIndex(0);
    int64_t top_c_index = fTop.ConnectIndex(0);
    fBottom.SetIndex(-1);
    fTop.SetIndex(-1);
    this->Reference()->SetReference(this);
    

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	 {
         std::stringstream sout;
         sout << "Finalizando criacao do elemento ";
         this->Print(sout);
         LOGPZ_DEBUG(logger,sout.str())
	 }
#endif
	 
}

template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh, const TPZCompElHDivCollapsed<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy),fBottom(copy.fBottom), fTop(copy.fTop)
{
}

// NAO TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed(TPZCompMesh &mesh,
												 const TPZCompElHDivCollapsed<TSHAPE> &copy,
												 std::map<int64_t,int64_t> & gl2lcConMap,
												 std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap),fBottom(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap),
fTop(mesh,copy.fBottom,gl2lcConMap,gl2lcElMap)
{
	
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<=TSHAPE::NFacets;i++)
	{
		int64_t lcIdx = -1;
		int64_t glIdx = copy.fConnectIndexes[i];
		if(glIdx == -1)
		{
			// nothing to clone
			this->fConnectIndexes[i] = -1;
			continue;
		}
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end())
		{
			lcIdx = gl2lcConMap[glIdx];
		}
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
    // write the code when this constructor is called
    DebugStop();
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::TPZCompElHDivCollapsed() :
TPZRegisterClassId(&TPZCompElHDivCollapsed::ClassId),
TPZCompElHDiv<TSHAPE>(),fBottom(),fTop()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NFacets+1;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

// TESTADO
template<class TSHAPE>
TPZCompElHDivCollapsed<TSHAPE>::~TPZCompElHDivCollapsed(){
    TPZGeoEl *gel = this->Reference();
    if(!gel) return;
    if (gel && gel->Reference() != this) {
        return;
    }
    int side = TSHAPE::NSides-1;
    TPZGeoElSide gelside(this->Reference(),side);
    TPZStack<TPZCompElSide> celstack;
    TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
    if (largecel) {
        int cindex = this->SideConnectLocId(0, side);
        TPZConnect &c = this->Connect(cindex);
        c.RemoveDepend();
    }
    gelside.HigherLevelCompElementList3(celstack, 0, 1);
    int64_t ncel = celstack.size();
    for (int64_t el=0; el<ncel; el++) {
        TPZCompElSide celsidesmall = celstack[el];
        TPZGeoElSide gelsidesmall = celsidesmall.Reference();
        if (gelsidesmall.Dimension() != gel->Dimension()) {
            continue;
        }
        TPZCompEl *cel = celsidesmall.Element();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if (!intel) {
            DebugStop();
        }
        int cindex = intel->SideConnectLocId(0, celsidesmall.Side());
        TPZConnect &c = intel->Connect(cindex);
        c.RemoveDepend();
    }
    if (gel){
        gel->ResetReference();
    }

}

// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides + 1) {
        DebugStop();
    }


    if(side < TSHAPE::NSides-1) TPZCompElHDiv<TSHAPE>::SetSideOrient(side, sideorient);
    else if(side == TSHAPE::NSides-1 ) fBottom.SetSideOrient(side, sideorient);
    else fTop.SetSideOrient(side-1, sideorient);
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::GetSideOrient(int side)
{
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides + 1) {
        DebugStop();
    }

    if (side < TSHAPE::NSides - 1) {
        return TPZCompElHDiv<TSHAPE>::GetSideOrient(side);
    }
    else if(side == TSHAPE::NSides-1) return fBottom.GetSideOrient(side);
    else return fTop.GetSideOrient(side-1);
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnects() const {
	
	return TPZCompElHDiv<TSHAPE>::NConnects()+2;
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	if(i <= TSHAPE::NFacets)
	{
        SetConnectIndex(i, connectindex);
	}
    else if(i == TSHAPE::NFacets+1)
    {
        fBottom.SetConnectIndex(0, connectindex);
    }
    else if(i == TSHAPE::NFacets+2)
    {
        fTop.SetConnectIndex(0, connectindex);
    }
    else
    {
        DebugStop();
    }
	
}

template<class TSHAPE>
int TPZCompElHDivCollapsed<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
    if(connect <= TSHAPE::NFacets)
    {
        return TPZCompElHDiv<TSHAPE>::NConnectShapeF(connect, connectorder);
    }
    else if(connect == TSHAPE::NFacets+1)
    {
        return fBottom.NConnectShapeF(0, connectorder);
    }
    else if(connect == TSHAPE::NFacets+2)
    {
        return fTop.NConnectShapeF(0, connectorder);
    }
    else DebugStop();
    return -1;
}


/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialData(TPZMaterialData &data){
    auto *tmp =
        dynamic_cast<TPZMaterialDataT<STATE>*>(&data);
    if(tmp){
        InitMaterialDataT(*tmp);
    }else{
        auto *tmp =
        dynamic_cast<TPZMaterialDataT<CSTATE>*>(&data);
        if(tmp){
            InitMaterialDataT(*tmp);
        }
    }
}
template<class TSHAPE>
template<class TVar>
void TPZCompElHDivCollapsed<TSHAPE>::InitMaterialDataT(TPZMaterialDataT<TVar> &data)
{
	TPZCompElHDiv<TSHAPE>::InitMaterialData(data);
    if(data.fUserData) DebugStop();
    auto datapair = new std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>>;
    data.fUserData = datapair;
    TPZMaterialData &datatop = datapair->second, &databottom = datapair->first;
    fTop.InitMaterialData(datatop);
    fBottom.InitMaterialData(databottom);
    // expand the shape vector and normal vector
//    int nvecshape = data.fVecShapeIndex.size();
//    int nscalar = data.phi.Rows();
//    int nscalartop = datatop.phi.Rows();
//    int nscalarbottom = databottom.phi.Rows();
//    int nvec = data.fVecShapeIndex.size(); // same as nvecshape
//    const int dim = TSHAPE::Dimension;
//    data.fMasterDirections.Resize(dim+1, nvec+2);
//    data.fMasterDirections(dim,nvec) = 1.;
//    data.fMasterDirections(dim,nvec+1) = -1.;
    // nvec is the number of vector shapes in this element. Then we need to sum the vec shapes of top and bottom
//    data.fDeformedDirections.Resize(3, nvec+nscalartop+nscalarbottom);
//    data.fVecShapeIndex.Resize(nvecshape+nscalartop+nscalarbottom, {0,0});
//    for(int i=0; i<nscalartop; i++) data.fVecShapeIndex[nvecshape+i] = std::pair<int,int64_t>(nvec,nscalar+i);
//    for(int i=0; i<nscalarbottom; i++) data.fVecShapeIndex[nvecshape+nscalartop+i] = std::pair<int,int64_t>(nvec+1,nscalar+nscalartop+i);
//    data.phi.Resize(nscalar+nscalartop+nscalarbottom, 1);
//    data.fDPhi.Resize(dim+1,nscalar+nscalartop+nscalarbottom);
//    data.divphi.Resize(nvecshape+nscalartop+nscalarbottom,1);
    if(data.fNeedsDeformedDirectionsFad) DebugStop();

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
        std::stringstream sout;
        sout << "After InitMaterialData\n";
        data.Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
    
}

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension || point.size() != TSHAPE::Dimension ){
		DebugStop() ;
	}
    fBottom.SideShapeFunction(side, point, phi, dphi);
    
    return;
}

/** Compute the shape function at the integration point */
//template<class TSHAPE>
//void TPZCompElHDivCollapsed<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
//{
//    fBottom.Shape(pt, phi, dphi);
//}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZIntelGen<TSHAPE>::Read(buf,context);
    fBottom.Read(buf,context);
    fTop.Read(buf, context);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZIntelGen<TSHAPE>::Write(buf,withclassid);
    fBottom.Write(buf, false);
    fTop.Write(buf, false);
}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZCompElHDiv<TSHAPE>::Print(out);
    fBottom.Print(out);
    fTop.Print(out);
    
}

#include "pzvec_extras.h"

static void ExpandAxes(TPZFMatrix<REAL> &axinput, TPZMatrix<REAL> &axout)
{
    int dim = axinput.Rows();
    switch (dim) {
        case 1:
        {
            REAL norms[3];
            TPZManVector<REAL,3> v1(3),try1[3];
            for(int i=0; i<3; i++) v1[i] = axinput(0,i);
            for(int i=0;i<3;i++)
            {
                TPZManVector<REAL,3> uni(3,0.);
                uni[i] = 1.;
                try1[i].resize(3);
                Cross(v1,uni,try1[i]);
                norms[i] = Norm<REAL>(try1[i]);
            }
            int maxi = 0;
            if(norms[1]>norms[0]) maxi = 1;
            if(norms[2]>norms[maxi]) maxi = 2;
            for(int i=0; i<3; i++)
            {
                axout(0,i) = axinput(0,i);
                axout(1,i) = try1[maxi][i]/norms[maxi];
            }
        }
            break;
        case 2:
        {
            TPZManVector<REAL,3> v1(3),v2(3),v3(3);
            for(int i=0; i<3; i++)
            {
                v1[i] = axinput(0,i);
                v2[i] = axinput(1,i);
            }
            Cross<REAL>(v1,v2,v3);
            for(int i=0; i<3; i++)
            {
                axout(0,i) = v1[i];
                axout(1,i) = v2[i];
                axout(2,i) = v3[i];
            }

        }
            break;
        default:
            DebugStop();
            break;
    }
}

/*
template<class TSHAPE>
template<class TVar>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){
    
//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZCompElHDiv<TSHAPE>::ComputeRequiredData(data, qsi);
    data.fNeedsSol = needsol;
    TPZFNMatrix<9,REAL> axeslocal(TSHAPE::Dimension+1,3);
    ExpandAxes(data.axes, axeslocal);
    data.axes = axeslocal;
    const int dim = TSHAPE::Dimension+1;
    const int nvecshapestd = data.fDeformedDirections.Cols();
    TPZManVector<REAL> topdir(dim,0.), botdir(dim,0.); // top and bot directions in the deformed element
    TPZManVector<REAL,3> vecup={0,0,0}, vecdown={0,0,0};
    vecup[dim-1] = 1.;
    vecdown[dim-1] = -1.;
    // compute the deformed directions for the two additional vectors
    {
        for(int i=0; i<3; i++){
                topdir[i] = data.axes(dim-1,i);
                botdir[i] = -data.axes(dim-1,i);
        }
    }
    
    std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>> *datapair = (std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>> *) data.fUserData;
    TPZMaterialDataT<TVar> &datatop = datapair->second, &databottom = datapair->first;
    
    // compute the divergence of the top and bottom elements
    // the value is the value of the shape function times the sign of the vector in master direction
    {
        fTop.ComputeRequiredData(datatop, qsi);
        fBottom.ComputeRequiredData(databottom, qsi);
        int64_t numvec = data.divphi.Rows();
        int64_t numphi = data.phi.Rows();
        int64_t nvec_top = datatop.phi.Rows();
        int64_t nvec_bottom = databottom.phi.Rows();
        int64_t nvec_hdiv = numvec;
      
        // fDeformedDirections (for now) represents the H1 shape functions
        // times the element vectors. So, it is already the hdiv shape function itself.
        // Its size is, therefore, the size for a standard 2d hdiv element, plus
        // the shape functions related to the top and bottom connect that communicate
        // with the adjacent 3D elements
        const int64_t nvecshapecollpased = nvecshapestd+nvec_top+nvec_bottom;
        data.fDeformedDirections.Resize(dim,nvecshapecollpased);
        data.fVecShapeIndex.Resize(nvecshapecollpased);
        data.divphi.Resize(nvecshapecollpased,1);
        data.phi.Resize(nvecshapecollpased,1);
        for(int i=numvec; i<nvecshapecollpased; i++)
        {
            data.fVecShapeIndex[i] = std::pair<int,int64_t>(i,i);
            data.phi(i) = 1.;
        }

        // First we append the bottom shapes and then the top shapes
        for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
            for (int d = 0; d < 3; d++) {
                data.fDeformedDirections(d,i) = databottom.phi(i-nvec_hdiv)*botdir[d];
            }
        }
        for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
            for (int d = 0; d < 3; d++) {
                data.fDeformedDirections(d,i) = datatop.phi(i-nvec_hdiv-nvec_bottom)*topdir[d];
            }
        }
        // Same for divphi
        for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
            data.divphi(i,0) = -databottom.phi(i-nvec_hdiv);
        }
        for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
            data.divphi(i,0) = datatop.phi(i-nvec_hdiv-nvec_bottom);
        }
    }
    
    if (data.fNeedsSol) {
        TPZCompElHDiv<TSHAPE>::ReallyComputeSolution(data);
    }


#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        data.fDeformedDirections.Print("Normal Vectors " , sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    

}//void
*/

template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) {

    TPZShapeData &shapedata = data;

    TPZFMatrix<REAL> auxPhi;
    TPZShapeHDiv<TSHAPE>::Shape(qsi, shapedata, auxPhi, data.divphi);
    
    TPZFMatrix<REAL> gradx(3,TSHAPE::Dimension,0.);
    this->Reference()->GradX(qsi, gradx);
//    TPZFMatrix<REAL> phiSHdiv;
    gradx.Multiply(auxPhi,data.fDeformedDirections);

    
    TPZFNMatrix<9,REAL> axeslocal(TSHAPE::Dimension+1,3);
    ExpandAxes(data.axes, axeslocal);
    data.axes = axeslocal;
    const int dim = TSHAPE::Dimension+1;
    const int nvecshapestd = data.fDeformedDirections.Cols();
    TPZManVector<REAL> topdir(dim,0.), botdir(dim,0.); // top and bot directions in the deformed element
    TPZManVector<REAL,3> vecup={0,0,0}, vecdown={0,0,0};
    vecup[dim-1] = 1.;
    vecdown[dim-1] = -1.;
    // compute the deformed directions for the two additional vectors
    {
        for(int i=0; i<3; i++){
                topdir[i] = data.axes(dim-1,i);
                botdir[i] = data.axes(dim-1,i);
        }
    }
    
    std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *datapair = (std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *) data.fUserData;
    TPZMaterialDataT<STATE> &datatop = datapair->second, &databottom = datapair->first;
    
    // compute the divergence of the top and bottom elements
    // the value is the value of the shape function times the sign of the vector in master direction
    fTop.ComputeRequiredData(datatop, qsi);
    fBottom.ComputeRequiredData(databottom, qsi);
    int toporder = fTop.Connect(0).Order();
    int bottomorder = fBottom.Connect(0).Order();
    int64_t numvec = data.divphi.Rows();
    int64_t numphi = data.phi.Rows();
    int64_t nvec_top = datatop.phi.Rows();
    int64_t nvec_bottom = databottom.phi.Rows();
    int64_t nvec_hdiv = numvec;
  
    // fDeformedDirections (for now) represents the H1 shape functions
    // times the element vectors. So, it is already the hdiv shape function itself.
    // Its size is, therefore, the size for a standard 2d hdiv element, plus
    // the shape functions related to the top and bottom connect that communicate
    // with the adjacent 3D elements
    const int64_t nvecshapecollpased = nvecshapestd+nvec_top+nvec_bottom;
    data.fDeformedDirections.Resize(dim,nvecshapecollpased);
    data.divphi.Resize(nvecshapecollpased,1);

    data.fDeformedDirections *= 1./data.detjac;
    data.divphi *= 1/data.detjac;
    // First we append the bottom shapes and then the top shapes
    for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
        for (int d = 0; d < 3; d++) {
            data.fDeformedDirections(d,i) = databottom.phi(i-nvec_hdiv)*botdir[d];
        }
    }
    for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
        for (int d = 0; d < 3; d++) {
            data.fDeformedDirections(d,i) = datatop.phi(i-nvec_hdiv-nvec_bottom)*topdir[d];
        }
    }
    // Same for divphi
    for (int64_t i= nvec_hdiv; i<nvecshapecollpased-nvec_top; i++) {
        data.divphi(i,0) = databottom.phi(i-nvec_hdiv);
    }
    for (int64_t i= nvec_hdiv+nvec_bottom; i<nvecshapecollpased; i++) {
        data.divphi(i,0) = datatop.phi(i-nvec_hdiv-nvec_bottom);
    }

    
    // alternative formulation
    {
        // adding a column do gradx pointing to the normal direction
        TPZFNMatrix<9> gradxlocal(3,dim);
        for (int i=0; i<3; i++) {
            for (int d=0; d<dim-1; d++) {
                gradxlocal(i,d) = gradx(i,d);
            }
            gradxlocal(i,dim-1) = axeslocal(dim-1,i);
        }
        TPZShapeData localshapedata;
        TPZShapeHDivCollapsed<TSHAPE> shape;
        TPZFNMatrix<27> locphi(dim,nvecshapecollpased),locdivphi(nvecshapecollpased,1);
        TPZFNMatrix<27> locdefphi(dim,nvecshapecollpased);
        TPZManVector<int> connectorders(data.fHDivConnectOrders);
        int ncon = NConnects();
        connectorders.Resize(ncon, 0);
        connectorders[ncon-2] = bottomorder;
        connectorders[ncon-1] = toporder;
        TPZManVector<int> sideorient(data.fSideOrient);
        sideorient.Resize(TSHAPE::NFacets+2,0);
        sideorient[TSHAPE::NFacets] = fBottom.GetSideOrient(TSHAPE::NSides-1);
        sideorient[TSHAPE::NFacets+1] = fTop.GetSideOrient(TSHAPE::NSides-1);
        shape.Initialize(data.fCornerNodeIds, connectorders, sideorient, localshapedata);
        shape.Shape(qsi,localshapedata,locphi,locdivphi);
        locdivphi *= 1./data.detjac;
        gradxlocal.Multiply(locphi,locdefphi);
        locdefphi *= 1./data.detjac;
        
        TPZFNMatrix<50> differencephi,differencedivphi;
        differencephi = data.fDeformedDirections - locdefphi;
        differencedivphi = data.divphi - locdivphi;
        auto normdiffphi = Norm(differencephi);
        auto normdiffdivphi = Norm(differencedivphi);
        std::cout << "normdiffphi " << normdiffphi << " normdiffdivphi " << normdiffdivphi << std::endl;
        if(normdiffphi > 1.e-10 || normdiffdivphi > 1.e-10)
        {
            differencephi.Print("diffphi",std::cout,EMathematicaInput);
            differencedivphi.Print("diffdiv",std::cout,EMathematicaInput);
            DebugStop();
        }
    }

    data.fVecShapeIndex.Resize(nvecshapecollpased);
    data.phi.Resize(nvecshapecollpased,1);
    for(int i=numvec; i<nvecshapecollpased; i++)
    {
        data.fVecShapeIndex[i] = std::pair<int,int64_t>(i,i);
        data.phi(i) = 1.;
    }
}



template<class TSHAPE>
int64_t TPZCompElHDivCollapsed<TSHAPE>::ConnectIndex(int con) const
{
    if(con <= TSHAPE::NFacets) return TPZCompElHDiv<TSHAPE>::ConnectIndex(con);
    if(con > TSHAPE::NFacets + 2) DebugStop();
    if(con == TSHAPE::NFacets+1) return fBottom.ConnectIndex(0);
    if(con == TSHAPE::NFacets+2) return fTop.ConnectIndex(0);
    DebugStop();
    return -1;
}

/**
 * @brief Destroy internally allocated data structures
 */
//auto datapair = new std::pair<TPZMaterialDataT<TVar>,TPZMaterialDataT<TVar>>;
template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::CleanupMaterialData(TPZMaterialData &data)
{
    TPZMaterialDataT<STATE> *dataS = dynamic_cast<TPZMaterialDataT<STATE> *>(&data);
    if(dataS)
    {
        std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *userdataS = (std::pair<TPZMaterialDataT<STATE>,TPZMaterialDataT<STATE>> *) data.fUserData;
        if(userdataS) delete userdataS;
    }
    TPZMaterialDataT<CSTATE> *dataC = dynamic_cast<TPZMaterialDataT<CSTATE> *>(&data);
    if(dataC)
    {
        std::pair<TPZMaterialDataT<CSTATE>,TPZMaterialDataT<CSTATE>> *userdataC = (std::pair<TPZMaterialDataT<CSTATE>,TPZMaterialDataT<CSTATE>> *) data.fUserData;
        if(userdataC) delete userdataC;
    }
    data.fUserData = nullptr;
}




template<class TSHAPE>
void TPZCompElHDivCollapsed<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsHDiv();
}

#include "pzshapetriang.h"
#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"

using namespace pzshape;

template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElHDivCollapsed<TPZShapeQuad>>;


template class TPZCompElHDivCollapsed<TPZShapeTriang>;
template class TPZCompElHDivCollapsed<TPZShapeLinear>;
template class TPZCompElHDivCollapsed<TPZShapeQuad>;

