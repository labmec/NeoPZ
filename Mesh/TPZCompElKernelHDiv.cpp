/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDiv.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel2D.h"
#include "pzcmesh.h"
#include "TPZMatSingleSpace.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template <class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh, TPZGeoEl *gel) : TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId), TPZCompElH1<TSHAPE>(mesh, gel, 1)
{
    gel->SetReference(this);
    int matid = gel->MaterialId();
#ifdef PZDEBUG
    if (mesh.FindMaterial(matid) == 0)
    {
        DebugStop();
    }
#endif

    this->fPreferredOrder = mesh.GetDefaultOrder();
    for (int i = 0; i < TSHAPE::NSides; i++)
    {
        fConnectIndexes[i] = this->CreateMidSideConnect(i);
        mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
    }

    // Updates the number of shape functions and also the integration rule
    for (int icon = 0; icon < this->NConnects(); icon++)
    {
        TPZConnect &c = this->Connect(icon);
        int nShapeF = NConnectShapeF(icon, c.Order());
        if (c.NShape() != nShapeF)
            DebugStop();
            // c.SetNShape(nShapeF);

#ifdef PZDEBUG
        if (c.HasDependency())
        {
            TPZConnect::TPZDepend<STATE> *dep = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(c.FirstDepend());
            while (dep)
            {
                if (dep->fDepMatrix.Rows() != c.NShape())
                {
                    DebugStop();
                }
                dep = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(dep->fNext);
            }
        }
#endif
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial *mat = this->Material();
        if (mat)
            nvar = mat->NStateVariables();
        int blsize = this->Mesh()->Block().Size(seqnum);
        if (blsize != nShapeF * nvar)
            DebugStop();
        // this->Mesh()->Block().Set(seqnum, nvar * nShapeF);
    }
    this->AdjustIntegrationRule();
}

template <class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
    data.gelElId = this->Reference()->Id();
    auto *mat =
        dynamic_cast<TPZMatSingleSpace *>(this->Material());
#ifdef PZDEBUG
    if (!mat)
    {
        DebugStop();
    }
#endif

    TPZManVector<int64_t, TSHAPE::NCornerNodes> ids(TSHAPE::NCornerNodes);
    TPZManVector<int, TSHAPE::NSides> orders(TSHAPE::NSides - TSHAPE::NCornerNodes);
    TPZManVector<int, TSHAPE::NFacets> sideorient(TSHAPE::NFacets, 0);
    TPZGeoEl *gel = this->Reference();
    for (int i = 0; i < TSHAPE::NCornerNodes; i++)
        ids[i] = gel->Node(i).Id();
    for (int i = TSHAPE::NCornerNodes; i < TSHAPE::NSides; i++)
        orders[i - TSHAPE::NCornerNodes] = this->Connect(i).Order();
    if (TSHAPE::Type() == ETriangle)
    {
        orders[TSHAPE::NSides - TSHAPE::NCornerNodes - 1]++;
    }
    TPZShapeData &shapedata = data;

    TPZShapeH1<TSHAPE>::Initialize(ids, orders, shapedata);

    mat->FillDataRequirements(data);
    const int dim = this->Dimension();
    const int nshape = data.fH1.fPhi.Rows();
    const int nstate = this->Material()->NStateVariables();
    data.fShapeType = TPZMaterialData::EScalarShape;
    data.phi.Redim(nshape, 1);
    data.fDeformedDirections.Redim(3, nshape);
    //    data.dphi.Redim(dim,nshape);
    data.dphix.Redim(dim, nshape);
    data.axes.Redim(dim, 3);
    data.jacobian.Redim(dim, dim);
    data.jacinv.Redim(dim, dim);
    data.x.Resize(3);
    if (data.fNeedsSol)
    {
        uint64_t ulen, durow, ducol;
        mat->GetSolDimensions(ulen, durow, ducol);
        data.SetSolSizes(nstate, ulen, durow, ducol);
    }
    if (data.fNeedsNeighborCenter)
    {
        TPZGeoElSide gelside(gel);
        gelside.CenterX(data.XCenter);
    }

    // TPZCompElH1<TSHAPE>::InitMaterialData(data);

    // int nshape = this->NShapeF();

    data.divphi.Resize(nshape, 1);
    data.divphi.Zero();

    data.fVecShapeIndex.Resize(nshape);
    for (int i = 0; i < nshape; i++)
    {
        data.fVecShapeIndex[i] = std::make_pair(i, 1);
    }
}

template <class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data)
{

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;
    TPZCompElH1<TSHAPE>::ComputeShape(qsi, data);
    auto dim = TSHAPE::Dimension;

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    TPZFMatrix<REAL> auxPhi(dim, nshape);
    auxPhi.Zero();

    TPZShapeHDivKernel2D<TSHAPE>::Shape(qsi, shapedata, auxPhi, data.divphi);

    switch (dim)
    {
    case 1:
        if (data.phi.Rows() > 1)
        {
            for (int i = 0; i < data.phi.Rows(); i++)
            {
                data.phi(i, 0) = -auxPhi(0, i) / data.detjac;    // Used to assemble the load vector
                data.fDeformedDirections(0, i) = data.phi(i, 0); // Used to compute the solution
            }
        }
        break;
    case 2:
    {
        TPZFMatrix<REAL> gradx;
        this->Reference()->GradX(qsi, gradx);
        const auto nCorner = TSHAPE::NCornerNodes;
        gradx.MultAdd(auxPhi, data.phi, data.fDeformedDirections, 1 / data.detjac, 0);
        data.phi.Resize(auxPhi.Cols(), 3);
        data.phi = 1.;
    }
    break;

    default:
        DebugStop();
        break;
    }

} // void

template <class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Solution(TPZVec<REAL> &qsi, int var, TPZVec<STATE> &sol)
{

    TPZMaterialDataT<STATE> data;
    constexpr bool hasPhi{false};
    this->ComputeSolution(qsi, data, hasPhi);

    sol.Resize(3);

    // REAL Sol = data.sol[0];
    // data.sol.resize(3);
    // data.sol[0] = Sol;
    sol[0] = data.sol[0][0];
    data.sol[0].Resize(3);
    // sol = std::move(data.sol[0]);
}

template <class TSHAPE>
template <class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{

    const int dim = 3;
    const int nstate = this->Material()->NStateVariables();
    const int ncon = this->NConnects();

    TPZFMatrix<TVar> &MeshSol = this->Mesh()->Solution();

    int64_t numbersol = MeshSol.Cols();

    if (numbersol != 1)
    {
        DebugStop();
    }
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    data.divsol.Resize(numbersol);

    for (int64_t is = 0; is < numbersol; is++)
    {
        data.sol[is].Resize(dim * nstate);
        data.sol[is].Fill(0);
        data.dsol[is].Redim(dim * nstate, dim);
        data.divsol[is].Resize(nstate);
        data.divsol[is].Fill(0.);
    }
    TPZFNMatrix<220, REAL> dphix(3, data.dphix.Cols());
    TPZFMatrix<REAL> &dphi = data.dphix;
    ;

    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    TPZBlock &block = this->Mesh()->Block();
    int ishape = 0, ivec = 0, counter = 0;

    for (int in = 0; in < ncon; in++)
    {
        TPZConnect *df = &this->Connect(in);
        int64_t dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
        // pos : position of the block in the solution matrix
        int64_t pos = block.Position(dfseq);

        /// ish loops of the number of shape functions associated with the block
        for (int ish = 0; ish < dfvar / nstate; ish++)
        {
            ishape = data.fVecShapeIndex[counter].first;
            for (int idf = 0; idf < nstate; idf++)
            {
                TVar meshsol = MeshSol(pos + ish * nstate + idf, 0);
                // Compute rotated flux
                data.sol[0][dim * idf + 0] -= meshsol * data.fDeformedDirections(0, ishape);
                data.sol[0][dim * idf + 1] -= meshsol * data.fDeformedDirections(1, ishape);
                data.divsol[0][idf] += meshsol * data.fH1.fPhi(ishape, 0);
            }
            counter++;
        }
    }
}

template <class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NConnectShapeF(int connect, int order) const
{
    if (TSHAPE::Type() == ETriangle && connect == 6)
        order++; // For triangles, the internal shape function is one degree higher
    return TPZCompElH1<TSHAPE>::NConnectShapeF(connect, order);
}

#include "pzshapepoint.h"
#include "pzshapelinear.h"

using namespace pzshape;

template <class TSIDESHAPE>
static void SideShape(const TPZVec<REAL> &point, TPZVec<int64_t> &ids, TPZVec<int> &ord, TPZShapeData &data)
{
    TPZShapeH1<TSIDESHAPE>::Initialize(ids, ord, data);
    TPZShapeH1<TSIDESHAPE>::Shape(point, data);
}

/**compute the values of the shape function of the side*/
template <class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SideShapeFunction(int side, TPZVec<REAL> &point, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{

    TPZGeoEl *ref = this->Reference();
    if (side == TSHAPE::NSides - 1)
    {
        DebugStop();
        return;
    }
    int nc = TSHAPE::NContainedSides(side);
    int nn = TSHAPE::NSideNodes(side);
    TPZManVector<int64_t, 27> id(nn);
    TPZManVector<int, 27> order(nc - nn);
    int n, c;
    for (n = 0; n < nn; n++)
    {
        int nodloc = TSHAPE::SideNodeLocId(side, n);
        id[n] = ref->NodePtr(nodloc)->Id();
    }
    for (c = nn; c < nc; c++)
    {
        int conloc = TSHAPE::ContainedSideLocId(side, c);
        order[c - nn] = this->Connect(conloc).Order();
    }

    TPZShapeData data;

    switch (TSHAPE::Type(side))
    {
    case EPoint:
        SideShape<TPZShapePoint>(point, id, order, data);
        break;
    case EOned:
        SideShape<TPZShapeLinear>(point, id, order, data);
        break;
    default:
        DebugStop();
        break;
    }
    phi = data.fH1.fPhi;
    dphi = data.fH1.fDPhi;
}

#include "pzshapepoint.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"

using namespace pzshape;

template class TPZCompElKernelHDiv<TPZShapeLinear>;
template class TPZCompElKernelHDiv<TPZShapeTriang>;
template class TPZCompElKernelHDiv<TPZShapeQuad>;

// BC
TPZCompEl *CreateHDivKernelBoundLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HDivFamily hdivfam)
{
    return new TPZCompElKernelHDiv<TPZShapeLinear>(mesh, gel);
}

// Domain
TPZCompEl *CreateHDivKernelQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HDivFamily hdivfam)
{
    return new TPZCompElKernelHDiv<TPZShapeQuad>(mesh, gel);
}
TPZCompEl *CreateHDivKernelTriangleEl(TPZGeoEl *gel, TPZCompMesh &mesh, const HDivFamily hdivfam)
{
    return new TPZCompElKernelHDiv<TPZShapeTriang>(mesh, gel);
}