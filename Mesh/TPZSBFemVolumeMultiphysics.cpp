
//
//  TPZSBFemVolumeMultiphysics.cpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#include "TPZSBFemVolumeL2.h"
#include "TPZSBFemVolumeHdiv.h"
#include "TPZSBFemVolumeMultiphysics.h"
#include "pzgeoelside.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "TPZNullMaterialCS.h"
#include "pzcmesh.h"
#include "pzmultiphysicscompel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
#endif

template<class TGeometry>
TPZSBFemVolumeMultiphysics<TGeometry>::TPZSBFemVolumeMultiphysics(TPZCompMesh & mesh, TPZGeoEl * gel) : TPZMultiphysicsCompEl<TGeometry>(mesh, gel)
{
    fElementVec1D.Resize(7);
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::AddElement1D(TPZCompEl * cel, int localindex)
{
    fElementVec1D[localindex] = cel;
    auto ncon = fConnectIndexes.size();
    auto nconcel =cel->NConnects();
    fConnectIndexes.Resize(ncon+nconcel);
    for (auto i = 0; i < nconcel; i++)
    {
        fConnectIndexes[i+ncon] = cel->ConnectIndex(i);
    }
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::AddElement(TPZCompEl *cel, int64_t meshindex)
{
    if (fElementVec.size() <= meshindex) 
    {
        fElementVec.resize(meshindex+1);
        TPZMultiphysicsElement::fActiveApproxSpace.Resize(meshindex+1, ETestTrial);
    }
    if (cel)
    {
        TPZGeoEl *gel = cel->Reference();
        TPZCompElSide celside(cel,gel->NSides()-1);
        fElementVec[meshindex] = celside;
    }
    else
    {
        fElementVec[meshindex] = TPZCompElSide();
    }
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::SetElementGroupIndex(int64_t index)
{
    fElementGroupIndex = index;
    TPZCompEl *celgr = this->Mesh()->Element(index);
    fElementGroup = celgr;
}

template<class TGeometry>
TPZCompEl * TPZSBFemVolumeMultiphysics<TGeometry>::Element(int64_t elindex)
{
    return fElementVec[elindex].Element();
}

template<class TGeometry>
TPZCompEl * TPZSBFemVolumeMultiphysics<TGeometry>::ReferredElement(int64_t mesh)
{
#ifdef PZDEBUG
    if (fElementVec.size() <= mesh) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " index does not exist!\n";
        DebugStop();
    };
#endif
    return fElementVec[mesh].Element();
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::SetSkeleton(int64_t skeleton)
{
    fSkeleton = skeleton;
}

template<class TGeometry>
int64_t TPZSBFemVolumeMultiphysics<TGeometry>::SkeletonIndex()
{
    return fSkeleton;
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd)
{
    auto sbfemflux = dynamic_cast<TPZSBFemVolumeHdiv * >(fElementVec[0].Element());
#ifdef PZDEBUG
    if(!sbfemflux) DebugStop();
#endif
    sbfemflux->LoadCoef(coef, coefd);

    auto sbfeml2 = dynamic_cast<TPZSBFemVolumeL2 * >(fElementVec[1].Element());
#ifdef PZDEBUG
    if(!sbfeml2) DebugStop();
#endif
    sbfeml2->LoadCoef(coef, coefd);
    
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::SetLocalIndices(TPZManVector<int64_t> &localindices, TPZManVector<int64_t> &localindicesint, TPZManVector<int64_t> &localindicesflux)
{
    auto sbfemflux = dynamic_cast<TPZSBFemVolumeHdiv * >(fElementVec[0].Element());
#ifdef PZDEBUG
    if(!sbfemflux) DebugStop();
#endif
    sbfemflux->SetLocalIndicesFlux(localindicesint,localindicesflux);

    auto sbfeml2 = dynamic_cast<TPZSBFemVolumeL2 * >(fElementVec[1].Element());
#ifdef PZDEBUG
    if(!sbfeml2) DebugStop();
#endif
    sbfeml2->SetLocalIndices(localindices);
}

template<class TGeometry>
int TPZSBFemVolumeMultiphysics<TGeometry>::NConnects() const
{
    return fConnectIndexes.size();   
}

template<class TGeometry>
int64_t TPZSBFemVolumeMultiphysics<TGeometry>::ConnectIndex(int i) const
{
    if (i > fConnectIndexes.size())
    {
        DebugStop();
    }        
    return fConnectIndexes[i];
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::CreateGraphicalElement(TPZGraphMesh &graphmesh, int dimension)
{
    TPZGeoEl *ref = this->Reference();
    if (ref->Dimension() != dimension) {
        return;
    }
    MElementType ty = ref->Type();
    if (ty == EQuadrilateral) {
        new TPZGraphElQ2dd(this, &graphmesh);
    } else if (ty == ECube) {
        new TPZGraphElQ3dd(this, &graphmesh);
    } else if (ty == EPrisma) {
        new TPZGraphElPrismMapped(this, &graphmesh);
    } else {
        DebugStop();
    }
}

template<class TGeometry>
TPZIntPoints & TPZSBFemVolumeMultiphysics<TGeometry>::GetIntegrationRule() const
{
    if(!fIntRule) DebugStop();
    return *fIntRule;
}

template<class TGeometry>
TPZIntPoints & TPZSBFemVolumeMultiphysics<TGeometry>::GetIntegrationRule()
{
    if(!fIntRule) InitializeIntegrationRule();
    return *fIntRule;
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    auto celskeleton = fElementVec1D[6];
    celskeleton->BuildCornerConnectList(connectindexes);
}

template<class TGeometry>
TPZCompEl * TPZSBFemVolumeMultiphysics<TGeometry>::Element(int elindex)
{
    return fElementVec1D[elindex];
}

template<class TGeometry>
TPZManVector<TPZCompEl * ,7> & TPZSBFemVolumeMultiphysics<TGeometry>::SBFemElementVec()
{
    return fElementVec1D;
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::SetPhiEigVal(TPZFMatrix<std::complex<double>> &phi, TPZManVector<std::complex<double>> &eigval)
{
    fEigenvalues = eigval;
    
    // int nrow = fLocalIndices.size();
    // fPhi.Resize(nrow, phi.Cols());
    // for (int i = 0; i < nrow; i++) {
    //     for (int j = 0; j < phi.Cols(); j++) {
    //         fPhi(i, j) = phi(fLocalIndices[i], j);
    //     }
    // }
    

    auto sbfeml2 = dynamic_cast<TPZSBFemVolumeL2 * >(fElementVec[1].Element());
#ifdef PZDEBUG
    if(!sbfeml2) DebugStop();
#endif
    sbfeml2->SetPhiEigVal(phi, eigval);
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::SetPhiFlux(TPZFMatrix<std::complex<double>> &phiflux, TPZManVector<std::complex<double>> &eigval)
{
    auto sbfemhdiv = dynamic_cast<TPZSBFemVolumeHdiv * >(fElementVec[0].Element());
#ifdef PZDEBUG
    if(!sbfemhdiv) DebugStop();
#endif
    sbfemhdiv->SetPhiFlux(phiflux, eigval);
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::AffineTransform(TPZVec<TPZTransform<> > &trVec) const
{
    int64_t nel;
    int side, dim, dimmf;
    nel=fElementVec.size();
    trVec.Resize(nel);
    TPZGeoEl *gelmf = this->Reference();
    dimmf = gelmf->Dimension();
    side = gelmf->NSides()-1;
    TPZGeoEl  *geoel;
    for (int64_t i = 0; i<nel; i++) {
        if (!fElementVec[i]) {
            continue;
        }
        geoel = fElementVec[i].Element()->Reference();
        dim =  geoel->Dimension();
        if (dim == dimmf) {
            TPZTransform<> tr(dim);
            TPZGeoElSide gelside(geoel,geoel->NSides()-1);
            TPZGeoElSide gelmfside(gelmf,side);
            if (gelside.NeighbourExists(gelmfside))
            {
                trVec[i] = tr;
                gelmfside.SideTransform3(gelside, trVec[i]);
            }
            else
            {
                trVec[i] = gelmf->BuildTransform2(side, geoel, tr);
            }
        }
        else
        {
            TPZTransform<> LocalTransf(dimmf);
            TPZGeoElSide thisgeoside(gelmf,gelmf->NSides()-1);
            TPZGeoElSide neighgeoside = fElementVec[i].Reference();
            thisgeoside.SideTransform3(neighgeoside, LocalTransf);
            TPZGeoElSide highdim(neighgeoside.Element(), neighgeoside.Element()->NSides()-1);
            trVec[i] = neighgeoside.SideToSideTransform(highdim).Multiply(LocalTransf);
            
        }
    }
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
*/
template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    if (var >= 99) {
        TPZCompEl::Solution(qsi, var, sol);
        return;
    }
    
    auto *material =
        dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(this->Material());
    if(!material){
        sol.Resize(0);
        return;
    }

    TPZGeoEl *Ref2D = this->Reference();
    int matid = Ref2D->MaterialId();
    auto *mat2d =
        dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(this->Mesh()->FindMaterial(matid));
    TPZMaterialDataT<STATE> data2d;
    
    TPZManVector<TPZTransform<> > trvec;
    AffineTransform(trvec);
    
    TPZManVector<REAL,3> myqsi(qsi);
    myqsi.resize(qsi.size());
    
    int64_t nref = fElementVec.size();
    TPZManVector<TPZMaterialDataT<STATE>,2> datavec;
    datavec.resize(nref);
    InitMaterialData(datavec);
    material->FillDataRequirements(datavec);
    for (int64_t iref = 0; iref<nref; iref++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        trvec[iref].Apply(qsi, myqsi);
        datavec[iref].p = msp->MaxOrder();
        
        TPZMaterialData::MShapeFunctionType shapetype = datavec[iref].fShapeType;
        msp->ComputeRequiredData(datavec[iref], myqsi);
        constexpr bool hasPhi{true};
        datavec[iref].xParametric = myqsi;
        msp->ComputeSolution(myqsi, datavec[iref],hasPhi);
        
        datavec[iref].x.Resize(3);
        msp->Reference()->X(myqsi, datavec[iref].x);
    }
    
    material->Solution(datavec, var, sol);
    for (int64_t iref = 0; iref<nref; iref++)
    {
        
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
        if(!msp) continue;
        msp->CleanupMaterialData(datavec[iref]);
    }
}

#include "pzaxestools.h"

template <class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::InitMaterialData(TPZVec<TPZMaterialDataT<STATE> > &dataVec, TPZVec<int64_t> *indices)
{
    int64_t nref = this->fElementVec.size();
    
#ifdef PZDEBUG
    if (nref != dataVec.size()) {
        PZError << "Error at " << __PRETTY_FUNCTION__ << " The number of materials can not be different from the size of the fElementVec !\n";
        DebugStop();
    }
#endif
    if(indices){
        int64_t nindices = indices->size();
        TPZVec<int> nshape(nindices);
        for (int64_t iref = 0; iref <nindices; iref++) {
            int64_t indiciref = indices->operator[](iref);
            if(fElementVec[indiciref])
            {
                dataVec[indiciref].gelElId = fElementVec[indiciref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[indiciref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[indiciref]);
        }
    }else{
        TPZVec<int> nshape(nref);
        for (int64_t iref = 0; iref < nref; iref++)
        {
            if(fElementVec[iref])
            {
                dataVec[iref].gelElId = fElementVec[iref].Element()->Reference()->Id();
            }
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(fElementVec[iref].Element());
            if (!msp) {
                continue;
            }
            // precisa comentar essa parte se for calcular os vetores no pontos de integracao.
            msp->InitMaterialData(dataVec[iref]);
        }
    }
    
    int n_active_approx_spaces = TPZMultiphysicsElement::fActiveApproxSpace.size();
    if (n_active_approx_spaces == 0) { /// it preserves the integrity for old version of multiphycis codes.
        TPZMultiphysicsElement::fActiveApproxSpace.Resize(nref, ETestTrial);
    }
    
    for (int64_t iref = 0; iref < nref; iref++) {
        dataVec[iref].fActiveApproxSpace = TPZMultiphysicsElement::fActiveApproxSpace[iref];
    }
    auto * mat =
        dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(this->Material());
    mat->FillDataRequirements(dataVec);
    
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::ComputeRequiredData(TPZVec<TPZMaterialDataT<STATE>> &dataVec,TPZVec<REAL> &qsi)
{
    if(dataVec.size())
    {
        for (int i = 0; i < 2; i++)        
        {
            auto sbfemvol = dynamic_cast<TPZInterpolationSpace *>(fElementVec[i].Element());
            if (!sbfemvol)
            {
                DebugStop();
            }
            
            sbfemvol->Reference()->X(qsi, dataVec[i].x);        

            TPZVec<TPZTransform<> > trvec;
            AffineTransform(trvec);
            TPZManVector<REAL> locpt(sbfemvol->Reference()->Dimension());
            trvec[i].Apply(qsi, locpt);

            sbfemvol->ComputeRequiredData(dataVec[i],locpt);
        }
    }
    else
        DebugStop();
}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::EvaluateError(TPZVec<REAL> &errors,bool store_error)
{
    auto *nullmat = dynamic_cast<TPZNullMaterialCS<STATE> *>(this->Material());
    if(nullmat) return;

    auto *mat = dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(this->Material());
    auto *matError = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);

    if (!matError || !(matError->HasExactSol()))
    {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" the material has no associated exact solution\n";
        PZError<<"Aborting...";
        DebugStop();
    }
    if (dynamic_cast<TPZBndCond *> (matError)) return;

    int problemdimension = this->Mesh()->Dimension();
    TPZGeoEl *ref = this->Reference();
    if (ref->Dimension() < problemdimension) return;

    int NErrors = matError->NEvalErrors();
    errors.Resize(NErrors);
    errors.Fill(0.);

    int dim = Dimension();
    TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides() - 1, 5);
    int maxIntOrder = intrule->GetMaxOrder();
    TPZManVector<int, 3> prevorder(dim), maxorder(dim, maxIntOrder);
    intrule->GetOrder(prevorder);
    intrule->SetOrder(maxorder);

    int ndof = this->Material()->NStateVariables();
    TPZManVector<STATE, 10> u_exact(ndof);
    TPZFNMatrix<90, STATE> du_exact(dim, ndof);

    TPZManVector<REAL, 10> intpoint(problemdimension), values(NErrors);
    values.Fill(0.0);
    REAL weight;
    TPZManVector<STATE, 9> flux_el(0, 0.);

    TPZManVector<TPZMaterialDataT<STATE>,2> datavec(2);
    InitMaterialData(datavec);

    const int64_t nref = fElementVec.size();
    int iactive = -1;
    for (unsigned int i = 0; i < nref; ++i) {
        if (datavec[i].fShapeType != TPZMaterialData::EEmpty) {
            datavec[i].fNeedsSol = true;
            iactive = i;
        }
    }
    if (iactive == -1) DebugStop();

    int nintpoints = intrule->NPoints();
    TPZFNMatrix<9, REAL> jac, axe, jacInv;
    REAL detJac;

    for (int nint = 0; nint < nintpoints; nint++)
    {
        intrule->Point(nint, intpoint, weight);

        ref->Jacobian(intpoint, jac, axe, detJac, jacInv);
        ComputeRequiredData(datavec, intpoint);
        
        weight *= fabs(detJac);
        matError->Errors(datavec, values);

        for (int ier = 0; ier < NErrors; ier++)
        {
            errors[ier] += values[ier] * weight;
        }
    }
    //Norma sobre o elemento
    for (int ier = 0; ier < NErrors; ier++) {
        errors[ier] = sqrt(errors[ier]);
    }//for ier

    if (store_error) {
        int64_t index = TPZMultiphysicsElement::Index();
        TPZFMatrix<STATE> &elvals = this->Mesh()->ElementSolution();
        if (elvals.Cols() < NErrors) {
            std::cout << "The element solution of the mesh should be resized before EvaluateError\n";
            DebugStop();
        }
        for (int ier = 0; ier < NErrors; ier++) {
            elvals(index, ier) = errors[ier];
        }
    }
    intrule->SetOrder(prevorder);

}

template<class TGeometry>
void TPZSBFemVolumeMultiphysics<TGeometry>::InitializeIntegrationRule()
{
    if (fIntRule) {
        DebugStop();
    }
    int nsides = this->Reference()->NSides();
    fIntRule = this->Reference()->CreateSideIntegrationRule(nsides-1, 1);
}

template<class TGeometry>
int TPZSBFemVolumeMultiphysics<TGeometry>::Dimension() const
{
    return this->Reference()->Dimension();
}

template<class TGeometry>
int TPZSBFemVolumeMultiphysics<TGeometry>::NShapeF() const
{
    int nc = fElementVec1D[6]->NConnects();
    int nshape = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = fElementVec1D[6]->Connect(ic);
        nshape += c.NShape();
    }
    return nshape;
}

TPZCompEl * CreateSBFemMultiphysicsLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh)
{
    return new TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoLinear>(mesh, gel);
}

TPZCompEl * CreateSBFemMultiphysicsQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh)
{
    return new TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad>(mesh, gel);
}

TPZCompEl * CreateSBFemMultiphysicsCubeEl(TPZGeoEl *gel, TPZCompMesh &mesh)
{
    return new TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoCube>(mesh, gel);
}

TPZCompEl * CreateSBFemMultiphysicsPrismaEl(TPZGeoEl *gel, TPZCompMesh &mesh)
{
    return new TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoPrism>(mesh, gel);    
}

template class TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoPoint>;
template class TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoLinear>;
template class TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoQuad>;
template class TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoCube>;
template class TPZSBFemVolumeMultiphysics<pzgeom::TPZGeoPrism>;
