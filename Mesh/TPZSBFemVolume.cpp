//
//  TPZSBFemVolume.cpp
//  PZ
//
//  Created by Philippe Devloo on 4/4/16.
//
//

#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "pzmultiphysicselement.h"
#include "TPZMaterialT.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "TPZMaterial.h"
#include "pzelmat.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "TPZBndCond.h"
#include "pzvec_extras.h"
#include "pzcmesh.h"
#include "pzaxestools.h"

#include <algorithm>

#ifdef LOG4CXX
static LoggerPtr loggercoefmatrices(Logger::getLogger("pz.mesh.sbfemvolume.coefmatrices"));
static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
static LoggerPtr loggerLBF(Logger::getLogger("pz.mesh.sbfemvolume.bodyloads"));
static LoggerPtr loggerEvaluateError(Logger::getLogger("pz.mesh.sbfemvolume.error"));
#endif

TPZSBFemVolume::TPZSBFemVolume(TPZCompMesh &mesh, TPZGeoEl *gel) : TPZInterpolationSpace(mesh, gel), fElementGroupIndex(-1), fSkeleton(-1), fDensity(1.) {

}



/// Compute the K matrices

void TPZSBFemVolume::ComputeKMatrices(TPZElementMatrixT<STATE> &E0, TPZElementMatrixT<STATE> &E1, TPZElementMatrixT<STATE> &E2, TPZElementMatrixT<STATE> &M0)
{
    // do all the computations here

    TPZElementMatrixT<STATE> efmat(Mesh(), TPZElementMatrixT<STATE>::EF);

    TPZGeoEl *Ref2D = Reference();
    TPZGeoMesh *gmesh = Ref2D->Mesh();

    TPZCompMesh *cmesh = Mesh();

    TPZMultiphysicsElement *mpSkel = dynamic_cast<TPZMultiphysicsElement *> (cmesh->Element(fSkeleton));
    TPZInterpolatedElement *CSkeleton = SkeletonElement();

    CSkeleton->InitializeElementMatrix(E0, efmat);
    CSkeleton->InitializeElementMatrix(E1, efmat);
    CSkeleton->InitializeElementMatrix(E2, efmat);
    CSkeleton->InitializeElementMatrix(M0, efmat);

    TPZGeoEl *Ref1D = CSkeleton->Reference();
    int dim1 = Ref1D->Dimension();

    int matid = Ref2D->MaterialId();
    int dim2 = Ref2D->Dimension();

    // find the first face side
    int nsides = Ref2D->NSides();
    int is;
    for (is = 0; is < nsides; is++) {
        if (Ref2D->SideDimension(is) == dim1) {
            break;
        }
    }
    int faceside = is;

    TPZGeoElSide thisside(Ref2D, faceside);

    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    auto *mat2dSingle = dynamic_cast<TPZMatSingleSpaceT<STATE>*>(mat2d);
    auto *mat2dMult = dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(mat2d);
    
    if (!mat2dSingle && !mat2dMult) {
        DebugStop();
    }

    int nstate = mat2d->NStateVariables();

    TPZGeoElSide SkeletonSide(Ref1D, Ref1D->NSides() - 1);

    TPZTransform<REAL> tr(dim2, dim1);
    tr = SkeletonSide.NeighbourSideTransform(thisside);
    TPZTransform<REAL> t2 = Ref2D->SideToSideTransform(thisside.Side(), Ref2D->NSides() - 1);
    tr = t2.Multiply(tr);
    // create a one-d integration rule
    TPZIntPoints &intpoints = CSkeleton->GetIntegrationRule();

    TPZMaterialDataT<STATE> data1d;
    TPZManVector<TPZMaterialDataT<STATE>,2> datavec(2);
    TPZMaterialDataT<STATE> &data2d = datavec[0];
    CSkeleton->InitMaterialData(data1d);
    CSkeleton->InitMaterialData(data2d);
    int nshape = data2d.fH1.fPhi.Rows();
    data2d.fH1.fPhi.Redim(nshape * 2, 1);
    data2d.fH1.fDPhi.Redim(dim2, 2 * nshape);
    data2d.dphix.Redim(dim2, 2 * nshape);
    data2d.dsol[0].Redim(dim2, nstate);

    TPZFNMatrix<200, STATE> ek(nshape * nstate * 2, nshape * nstate * 2, 0.), ef(nshape * nstate * 2, 1, 0.);
    int npoint = intpoints.NPoints();
    for (int ip = 0; ip < npoint; ip++) {
        TPZManVector<REAL, 3> xi(dim1), xiquad(dim2), xivol(dim2);
        REAL weight;
        intpoints.Point(ip, xi, weight);
        tr.Apply(xi, xiquad);
        xivol = xiquad;
        xivol[dim2 - 1] = -0.5;
        TPZFNMatrix<9, REAL> jacobian(dim1, dim1), axes(dim1, 3), jacinv(dim1, dim1);
        REAL detjac;
        Ref1D->Jacobian(xi, jacobian, axes, detjac, jacinv);
        Ref2D->Jacobian(xiquad, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
        Ref2D->X(xivol, data2d.x);
        CSkeleton->ComputeRequiredData(data1d, xi);
#ifdef PZDEBUG
        // if the dimension of the problem is 2, we assume that the 1D axes corresponds to the first axis of the 2D problem
        if (dim2 == 2) {
            REAL norm = 0;
            for (int i = 0; i < 3; i++) {
                norm += (axes(0, i) - data2d.axes(0, i))*(axes(0, i) - data2d.axes(0, i));
            }
            norm = sqrt(norm);
            if (norm > 1.e-8) {
                Ref1D->Print();
                Ref2D->Print();
                axes.Print("axes 1D", std::cout);
                data2d.axes.Print("axes 2D", std::cout);
                {
                    std::ofstream out("gmesh.txt");
                    gmesh->Print(out);
                }
                DebugStop();
            }
        }
#endif
        // adjust the axes of the 3D element to match the axes of the side element
        if (dim2 == 3) {
            //            TPZFNMatrix<9,REAL> jacorig(data2d.jacobian);
            AdjustAxes3D(axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
#ifdef LOG4CXX2
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "x 2d " << data1d.x << std::endl;
                data1d.axes.Print("axes 2D", sout);
                data2d.axes.Print("axes 3D", sout);
                data2d.jacobian.Print("jacobian", sout);
                //                jacorig.Print("jacobian original",sout);
                sout << "detjac = " << data2d.detjac << std::endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
        }
#ifdef LOG4CXX2
        if (logger->isDebugEnabled()) {
            std::stringstream sout;
            TPZFNMatrix<9> axest, gradx(3, dim2);
            data2d.axes.Transpose(&axest);
            axest.Multiply(data2d.jacobian, gradx);
            gradx.Print("gradx ", sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        ExtendShapeFunctions(data1d, data2d);

        weight *= fabs(data2d.detjac)*2.;
        for (int i = 0; i < nshape; i++) {
            for (int j = 0; j < nshape; j++) {
                for (int st = 0; st < nstate; st++) {
                    M0.fMat(i * nstate + st, j * nstate + st) += weight * data1d.fH1.fPhi(i, 0) * data1d.fH1.fPhi(j, 0) * fDensity;
                }
            }
        }
        // compute the contributions to K11 K12 and K22
        if (mat2dSingle) {
            mat2dSingle->Contribute(data2d, weight, ek, ef);
        } else if(mat2dMult) {
            if(Norm(data2d.axes)<1.e-6) {
                DebugStop();
            }
            // data2d.axes.Print("axes = ",std::cout);
            TPZAxesTools<REAL>::Axes2XYZ(data2d.dphix, data2d.fDeformedDirections, data2d.axes);
            int nshape = data2d.dphix.Cols();
            for(int ish = 0; ish<nshape; ish++) {
                REAL tmp = data2d.fDeformedDirections(0,ish);
                data2d.fDeformedDirections(0,ish) = -data2d.fDeformedDirections(1,ish);
                data2d.fDeformedDirections(1,ish) = tmp;
            }
            // data2d.dphix.Print("dphix = ",std::cout);
            // data2d.fDeformedDirections.Print("DeformedDirections = ",std::cout);
            // datavec[0].fDeformedDirections.Print("DeformedDirections = ",std::cout);
            mat2dMult->Contribute(datavec, weight, ek, ef);
        } else {
            DebugStop();
        }
    }
    for (int i = 0; i < nstate * nshape; i++) {
        for (int j = 0; j < nstate * nshape; j++) {
            E0.fMat(i, j) = ek(i, j);
            E1.fMat(j, i) = ek(i, j + nstate * nshape);
            E2.fMat(i, j) = ek(i + nstate*nshape, j + nstate * nshape);
        }
    }
#ifdef LOG4CXX
    if (loggercoefmatrices->isDebugEnabled()) {
        std::stringstream sout;
        E0.fMat.Print("E0 = ",sout);
        E1.fMat.Print("E1 = ",sout);
        E2.fMat.Print("E2 = ",sout);
        LOGPZ_DEBUG(loggercoefmatrices, sout.str())
    }
#endif
}

/// adjust the axes and jacobian of the 3D element
void TPZSBFemVolume::AdjustAxes3D(const TPZFMatrix<REAL> &axes2D, TPZFMatrix<REAL> &axes3D, TPZFMatrix<REAL> &jac3D, TPZFMatrix<REAL> &jacinv3D, REAL detjac)
{
    TPZManVector<REAL, 3> ax1(3), ax2(3), ax3(3);
    for (int i = 0; i < 3; i++) {
        ax1[i] = axes2D.g(0, i);
        ax2[i] = axes2D.g(1, i);
        Cross(ax1, ax2, ax3);
    }
    for (int i = 0; i < 3; i++) {
        axes3D(0, i) = ax1[i];
        axes3D(1, i) = ax2[i];
        axes3D(2, i) = ax3[i];
        if (detjac < 0.) {
            axes3D(2, i) *= -1.;
        }
    }
    TPZFNMatrix<9, REAL> jacnew(3, 3), axest(3, 3), jacinv(3, 3);
    axes3D.Transpose(&axest);
    axes3D.Multiply(jac3D, jacnew);
    jacinv3D.Multiply(axest, jacinv);
    jac3D = jacnew;
    jacinv3D = jacinv;
#ifdef PZDEBUG
    // check whether the axes are orthogonal and whether the jacobian is still the inverse of jacinv
    {
        TPZFNMatrix<9, REAL> ident1(3, 3, 0.), ident2(3, 3, 0.), identity(3, 3);
        identity.Identity();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    ident1(i, j) += axes3D(i, k) * axes3D(j, k);
                    ident2(i, j) += jac3D(i, k) * jacinv3D(k, j);
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (fabs(ident1(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
                if (fabs(ident2(i, j) - identity(i, j)) > 1.e-6) {
                    DebugStop();
                }
            }
        }
    }
#endif
}

/// extend the border shape functions for SBFem computations
void TPZSBFemVolume::ExtendShapeFunctions(TPZMaterialDataT<STATE> &data1d, TPZMaterialDataT<STATE> &data2d)
{
    int dim = Reference()->Dimension();
    int64_t nshape = data2d.fH1.fPhi.Rows() / 2;
    for (int ish = 0; ish < nshape; ish++) {
        data2d.fH1.fPhi(ish + nshape, 0) = data1d.fH1.fPhi(ish, 0);
        for (int d = 0; d < dim - 1; d++) {
            data2d.fH1.fDPhi(d, ish + nshape) = data1d.fH1.fDPhi(d, ish);
            data2d.fH1.fDPhi(d, ish) = 0.;
        }
        data2d.fH1.fDPhi(dim - 1, ish) = -data1d.fH1.fPhi(ish) / 2.;
        data2d.fH1.fDPhi(dim - 1, ish + nshape) = 0.;
    }
    TPZInterpolationSpace::Convert2Axes(data2d.fH1.fDPhi, data2d.jacinv, data2d.dphix);
}

TPZCompEl * CreateSBFemCompEl(TPZGeoEl *gel, TPZCompMesh &mesh)
{
    return new TPZSBFemVolume(mesh, gel);
}

/// initialize the data structures of the eigenvectors and eigenvalues associated with this volume element
void TPZSBFemVolume::SetPhiEigVal(TPZFMatrix<std::complex<double> > &phi, TPZManVector<std::complex<double> > &eigval)
{
    fEigenvalues = eigval;
    int nrow = fLocalIndices.size();
    fPhi.Resize(nrow, phi.Cols());
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < phi.Cols(); j++) {
            fPhi(i, j) = phi(fLocalIndices[i], j);
        }
    }
}

void TPZSBFemVolume::SetCoefNonHomogeneous(TPZFNMatrix<100,std::complex<double>> &phi, TPZManVector<std::complex<double> > &eigval, TPZFNMatrix<100,std::complex<double> > &phiinv, TPZFNMatrix<100,std::complex<double> > &rot)
{
    fEigenvaluesBubble = eigval;
    int nrow = fLocalIndices.size();
    fPhiBubble.Resize(nrow, phi.Cols());
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < phi.Cols(); j++) {
            fPhiBubble(i, j) = phi(fLocalIndices[i], j);
        }
    }
    fPhiInv = phiinv;
    fPhiInvBubbles = rot;
}

/** @brief Loads the solution within the internal data structure of the element */

/**
 * Is used to initialize the solution of connect objects with dependency. \n
 * Is also used to load the solution within SuperElements
 */
void TPZSBFemVolume::LoadCoef(TPZFMatrix<std::complex<double> > &coef)
{
    fCoeficients = coef;
}

/// Return the Computational Skeleton element
TPZInterpolatedElement *TPZSBFemVolume::SkeletonElement() {
    TPZCompMesh *cmesh = Mesh();
    TPZMultiphysicsElement *mpSkel = dynamic_cast<TPZMultiphysicsElement *> (cmesh->Element(fSkeleton));
    TPZInterpolatedElement *CSkeleton = 0;
    if(mpSkel) {
        CSkeleton = dynamic_cast<TPZInterpolatedElement *> (mpSkel->Element(0));
    } else {
        CSkeleton = dynamic_cast<TPZInterpolatedElement *> (cmesh->Element(fSkeleton));
    }
    if(!CSkeleton) DebugStop();
    return CSkeleton;
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZSBFemVolume::ReallyComputeSolution(TPZMaterialDataT<STATE>& data)
{
    TPZVec<REAL> &qsi = data.xParametric;
    TPZSolVec<STATE> &sol = data.sol;
    TPZGradSolVec<STATE> &dsol = data.dsol;
    TPZFMatrix<REAL> &axes = data.axes;

    TPZCompMesh *cmesh = Mesh();
    sol.Resize(fCoeficients.Cols());
    dsol.Resize(fCoeficients.Cols());
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);

    int dim = Ref2D->Dimension();
    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }
    TPZInterpolatedElement *CSkeleton = SkeletonElement();

    TPZMaterialDataT<STATE> data1d, data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);
    CSkeleton->InitMaterialData(data1d);
    TPZGeoEl *Ref1D = CSkeleton->Reference();

    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    if (dim == 3) {

        AdjustAxes3D(data1d.axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
    }
    axes = data2d.axes;
    CSkeleton->ComputeRequiredData(data1d, qsilow);

    int nshape = data1d.fH1.fPhi.Rows();
    int nstate = mat2d->NStateVariables();
#ifdef PZDEBUG
    if (fPhi.Cols() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        TPZManVector<std::complex<double> > coefcol(fCoeficients.Rows());
        for (int i = 0; i < fCoeficients.Rows(); i++) {
            coefcol[i] = fCoeficients(i, 0);
        }
        std::stringstream sout;
        sout << "coefficients " << coefcol << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    for (int s = 0; s < sol.size(); s++) {
        TPZManVector<REAL, 10> uh_xi(fPhi.Rows(), 0.), Duh_xi(fPhi.Rows(), 0.);
        int nphixi = fPhi.Rows();
        int numeig = fPhi.Cols();
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += (fCoeficients(c, s) * xiexp * fPhi(i, c)).real();
                Duh_xi[i] += (-fCoeficients(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c)).real();
            }
        }
#ifdef LOG4CXX2
        if (s == 0 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.fH1.fPhi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        sol[s].Resize(nstate);
        sol[s].Fill(0.);
        TPZFNMatrix<9, STATE> dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                sol[s][istate] += data1d.fH1.fPhi(ishape) * uh_xi[ishape * nstate + istate];
                dsolxi[istate] += data1d.fH1.fPhi(ishape) * Duh_xi[ishape * nstate + istate];
                for (int d = 0; d < dim - 1; d++) {
                    dsolxieta(d, istate) += data1d.fH1.fDPhi(d, ishape) * uh_xi[ishape * nstate + istate];
                }
                dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
            }
        }
        dsol[s].Resize(dim, nstate);
        dsol[s].Zero();
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    dsol[s](d1, istate) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }
    TPZMatCombinedSpacesT<STATE> *matcomb = dynamic_cast<TPZMatCombinedSpacesT<STATE> *> (mat2d);
    if(matcomb) {
        data.divsol = sol;
        data.sol[0].resize(3);
        data.sol[0].Fill(0.);
        TPZFNMatrix<3,STATE> dsolLoc(3,1,0.);
        TPZAxesTools<STATE>::Axes2XYZ(dsol[0], dsolLoc, data.axes);
        data.sol[0][0] = -dsolLoc(1,0);
        data.sol[0][1] = dsolLoc(0,0);
    }
}

void TPZSBFemVolume::ComputeSolutionWithBubbles(TPZVec<REAL> &qsi,
                                       TPZSolVec<STATE> &sol, TPZGradSolVec<STATE> &dsol, TPZFMatrix<REAL> &axes)
{    
    TPZCompMesh *cmesh = Mesh();
    int nrows = fCoeficients.Rows();
    if (fCoeficients.Rows() ==0)
    {
        TPZSBFemElementGroup * sbgr = dynamic_cast<TPZSBFemElementGroup *> (fElementGroup);
        sbgr->LoadSolution();
    }
    
    sol.Resize(fCoeficients.Cols());
    dsol.Resize(fCoeficients.Cols());
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    
    int dim = Ref2D->Dimension();
    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }
    TPZInterpolatedElement *CSkeleton = SkeletonElement();
    TPZMaterialDataT<STATE> data1d, data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);
    CSkeleton->InitMaterialData(data1d);
    TPZGeoEl *Ref1D = CSkeleton->Reference();
    
    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    if (dim == 3)
    {
        AdjustAxes3D(data1d.axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
    }
    axes = data2d.axes;
    CSkeleton->ComputeRequiredData(data1d, qsilow);
    CSkeleton->Shape(qsilow, data1d.fH1.fPhi, data1d.fH1.fDPhi);
    
    int nshape = data1d.fH1.fPhi.Rows();
    int nstate = mat2d->NStateVariables();
    int numeig = fEigenvalues.size();
#ifdef PZDEBUG
    if (fPhiInv.Rows() + fPhiInvBubbles.Rows() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        TPZManVector<std::complex<double> > coefcol(fCoeficients.Rows());
        for (int i = 0; i < fCoeficients.Rows(); i++) {
            coefcol[i] = fCoeficients(i, 0);
        }
        std::stringstream sout;
        sout << "coefficients " << coefcol << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    for (int s = 0; s < sol.size(); s++)
    {
        int nphixi = fPhi.Rows();
        TPZManVector<std::complex<double>, 10> uh_xi(nphixi, 0.), Duh_xi(nphixi, 0.);

        TPZFNMatrix<200,std::complex<double>> umat0(nphixi, numeig, 0), Dumat0(nphixi, numeig, 0);
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += (fCoeficients(c, s) * xiexp * fPhi(i, c));
                Duh_xi[i] += (-fCoeficients(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c));
            }
        }
        int numeigbubbles = fEigenvaluesBubble.size();
        TPZFNMatrix<200,std::complex<double>> umat(nphixi, numeigbubbles, 0), Dumat(nphixi, numeigbubbles, 0);
        for (int c = 0; c < numeigbubbles; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvaluesBubble[c])) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvaluesBubble[c] + 1.)) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvaluesBubble[c]);
                xiexpm1 = pow(sbfemparam, -fEigenvaluesBubble[c] - 1.);
            }
            for (int i = 0; i < nphixi; i++) {
                // umat(i,c) += xiexp * fPhiBubble(i, c);
                // Dumat(i,c) += (fEigenvaluesBubble[c] + xiexpm1 * fPhiBubble(i, c));
                uh_xi[i] += (fCoeficients(c + numeig, s) * xiexp * fPhiBubble(i, c));
                Duh_xi[i] += (-fCoeficients(c + numeig, s)*(fEigenvaluesBubble[c]) * xiexpm1 * fPhiBubble(i, c));
            }
        }
        
#ifdef LOG4CXX
        if (s == 0 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.fH1.fPhi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        sol[s].Resize(nstate);
        sol[s].Fill(0.);
        TPZFNMatrix<9, STATE> dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                sol[s][istate] += data1d.fH1.fPhi(ishape) * uh_xi[ishape * nstate + istate].real();
                dsolxi[istate] += data1d.fH1.fPhi(ishape) * Duh_xi[ishape * nstate + istate].real();
                for (int d = 0; d < dim - 1; d++) {
                    dsolxieta(d, istate) += data1d.fH1.fDPhi(d, ishape) * uh_xi[ishape * nstate + istate].real();
                }
                dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
            }
        }
        dsol[s].Resize(dim, nstate);
        dsol[s].Zero();
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    dsol[s](d1, istate) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }
    
    
#ifdef LOG4CXX
   if (loggerEvaluateError->isDebugEnabled()) {
        TPZFNMatrix<200, std::complex<double> > sbfemmat(numeig,numeig,0.);
        for (int i=0; i<numeig; i++) {
            sbfemmat(i,i) = pow(sbfemparam, -Eigenvalues()[i].real());
        }
        TPZFMatrix<std::complex<double> > res;
        sbfemmat.Multiply(fPhiInv, res);
        TPZFNMatrix<200, std::complex<REAL>> umatrix;
        fPhi.Multiply(res, umatrix);
        std::stringstream sout;
        umatrix.Print("umatrix = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(loggerEvaluateError, sout.str());
   }
#endif
    
}

/**
 * @brief Computes the shape function set at the point x.
 * @param qsi point in master element coordinates
 * @param phi vector of values of shapefunctions, dimension (numshape,1)
 * @param dphi matrix of derivatives of shapefunctions in master element coordinates, dimension (dim,numshape)
 */

/**
 * This method uses the order of interpolation
 * of the element along the sides to compute the number of shapefunctions
 */
void TPZSBFemVolume::Shape(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphidxi)
{
    TPZCompMesh *cmesh = Mesh();
    if(fElementGroupIndex < 0)
    {
        DebugStop();
    }
    TPZCompEl *celgroup = cmesh->Element(fElementGroupIndex);
    TPZSBFemElementGroup *elgr = dynamic_cast<TPZSBFemElementGroup *> (celgroup);
    TPZFMatrix<std::complex<double> > &CoefficientLoc = elgr->PhiInverse();
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        CoefficientLoc.Print("Coefficients = ", sout, EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);
    int dim = Ref2D->Dimension();
    int nstate = mat2d->NStateVariables();

    phi.Redim(CoefficientLoc.Cols() * nstate, 1);
    dphidxi.Redim(dim*nstate, CoefficientLoc.Cols());

    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }
    TPZInterpolatedElement *CSkeleton = SkeletonElement();
    TPZMaterialDataT<STATE> data1d, data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);
    CSkeleton->InitMaterialData(data1d);
    TPZGeoEl *Ref1D = CSkeleton->Reference();

    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    if (dim == 3) {

        AdjustAxes3D(data1d.axes, data2d.axes, data2d.jacobian, data2d.jacinv, data2d.detjac);
    }
    CSkeleton->ComputeRequiredData(data1d, qsilow);

    int nshape = data1d.fH1.fPhi.Rows();
#ifdef PZDEBUG
    if (fPhi.Cols() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif
    phi.Zero();
#ifdef LOG4CXX2
    if (logger->isDebugEnabled()) {
        int eq = 1;
        TPZManVector<std::complex<double> > coefcol(CoefficientLoc.Rows());
        for (int i = 0; i < CoefficientLoc.Rows(); i++) {
            coefcol[i] = CoefficientLoc(i, eq);
        }
        std::stringstream sout;
        sout << "coefficients " << coefcol << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    for (int s = 0; s < CoefficientLoc.Cols(); s++) {
        TPZManVector<std::complex<double>, 10> uh_xi(fPhi.Rows(), 0.), Duh_xi(fPhi.Rows(), 0.);
        int nphixi = fPhi.Rows();
        int numeig = fPhi.Cols();
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += CoefficientLoc(c, s) * xiexp * fPhi(i, c);
                Duh_xi[i] += -CoefficientLoc(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c);
            }
        }
#ifdef LOG4CXX2
        if (s == 1 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.fH1.fPhi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZFNMatrix<9, STATE> dsollow(dim - 1, nstate, 0.), dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                phi(s * nstate + istate, 0) += data1d.fH1.fPhi(ishape) * uh_xi[ishape * nstate + istate].real();
                dsolxi[istate] += data1d.fH1.fPhi(ishape) * Duh_xi[ishape * nstate + istate].real();
                for (int d = 0; d < dim - 1; d++) {
                    dsollow(d, istate) += data1d.fH1.fDPhi(d, ishape) * uh_xi[ishape * nstate + istate].real();
                }
            }
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d = 0; d < dim - 1; d++) {
                dsolxieta(d, istate) = dsollow(d, istate);
            }
            dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    dphidxi(istate * nstate + d1, s) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }

}

/**
 * @brief Calculates the solution - sol - for the variable var
 * at point qsi, where qsi is expressed in terms of the
 * master element coordinates
 * @param qsi master element coordinate
 * @param var variable name
 * @param sol vetor for the solution
 */
void TPZSBFemVolume::Solution(TPZVec<REAL> &qsi, int var, TPZVec<STATE> &sol) {
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    TPZCompMesh *cmesh = Mesh();

    auto *mat2d =
        dynamic_cast<TPZMatSingleSpaceT<STATE>*>(cmesh->FindMaterial(matid));
    auto *mat2dCS = dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(cmesh->FindMaterial(matid));
    TPZManVector<TPZMaterialDataT<STATE>,2> datavec(2);
    TPZMaterialDataT<STATE> &data2d = datavec[0];

    if(TPZSBFemElementGroup::gDefaultPolynomialOrder == 0)
    {
        data2d.xParametric = qsi;
        ReallyComputeSolution(data2d);
    }
    else
    {
        ComputeSolutionWithBubbles(qsi, data2d.sol, data2d.dsol, data2d.axes);
    }
    data2d.x.Resize(3, 0.);
    Reference()->X(qsi, data2d.x);
    if(mat2d) {
        mat2d->Solution(data2d, var, sol);
    } else if(mat2dCS) {
        mat2dCS->Solution(datavec, var, sol);
    } else {
        DebugStop();
    }

}

void TPZSBFemVolume::CreateGraphicalElement(TPZGraphMesh &graphmesh, int dimension) {

    TPZGeoEl *ref = Reference();
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

void TPZSBFemVolume::EvaluateError(TPZVec<REAL> &errors,bool store_error)
{
    auto *material = this->Material();
	auto* matError =
        dynamic_cast<TPZMatErrorSingleSpace<STATE> *>(this->Material());
    auto *matErrorCS =
        dynamic_cast<TPZMatErrorCombinedSpaces<STATE> *>(this->Material());
    if (!matError && !matErrorCS)
    {
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" the material has no error interface\n";
        PZError<<"Aborting...";
        DebugStop();
    }
    if (matError && !matError->HasExactSol()) {
        std::cout << "Exiting EvaluateError - null error - no exact solution.";
        DebugStop();
    }
    if (matErrorCS && !matErrorCS->HasExactSol()) {
        std::cout << "Exiting EvaluateError - null error - no exact solution.";
        DebugStop();
    }
    if (dynamic_cast<TPZBndCond *> (matError))
    {
        std::cout << "Exiting EvaluateError - null error - boundary condition material.";
        DebugStop();
    }
    int NErrors = 0;
    if(matError)
    {
        NErrors = matError->NEvalErrors();
    }
    else
    {
        NErrors = matErrorCS->NEvalErrors();
    }
    errors.Resize(NErrors);
    errors.Fill(0.);
    int problemdimension = Mesh()->Dimension();
    TPZGeoEl *ref = Reference();
    if (ref->Dimension() < problemdimension) return;

    int dim = Dimension();
    TPZAutoPointer<TPZIntPoints> intrule = ref->CreateSideIntegrationRule(ref->NSides() - 1, 5);
    int maxIntOrder = intrule->GetMaxOrder();
    
    TPZManVector<int, 3> prevorder(dim), maxorder(dim, maxIntOrder);
    intrule->GetOrder(prevorder);
    intrule->SetOrder(maxorder);
    
    int ndof = material->NStateVariables();
    TPZManVector<STATE, 10> u_exact(ndof);
    TPZFNMatrix<90, STATE> du_exact(dim, ndof);
    TPZManVector<REAL, 10> intpoint(problemdimension), values(NErrors);
    values.Fill(0.0);
    REAL weight;
    TPZManVector<STATE, 9> flux_el(0, 0.);

    TPZManVector<TPZMaterialDataT<STATE>,2> datavec(2);
    TPZMaterialDataT<STATE> &data = datavec[0];
    data.x.Resize(3);
    int nintpoints = intrule->NPoints();

    for (int nint = 0; nint < nintpoints; nint++) {

        intrule->Point(nint, intpoint, weight);

        ref->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);

        weight *= fabs(data.detjac);
        
        if(TPZSBFemElementGroup::gDefaultPolynomialOrder == 0)
        {
            data.xParametric = intpoint;
            ReallyComputeSolution(data);
        }
        else
        {
            ComputeSolutionWithBubbles(intpoint, data.sol, data.dsol, data.axes);
        }
        //contribuicoes dos erros
        ref->X(intpoint, data.x);
        if(matError)
        {
            matError->Errors(data, values);
        }
        else
        {
            matErrorCS->Errors(datavec, values);
        }

        for (int ier = 0; ier < NErrors; ier++)
            errors[ier] += values[ier] * weight;

    }//fim for : integration rule
    //Norma sobre o elemento
    for (int ier = 0; ier < NErrors; ier++) {
        errors[ier] = sqrt(errors[ier]);
    }//for ier

    // intrule->SetOrder(prevorder);
    intrule->SetOrder(maxorder);

}

void TPZSBFemVolume::SetIntegrationRule(int order) 
{
    if (!fIntRule) {
        InitializeIntegrationRule();
    }
    int dim = Reference()->Dimension();
    TPZManVector<int, 3> ordervec(dim, order);
    ordervec[dim - 1] = 10;
    if (10 < order + 6) ordervec[dim - 1] = order + 6;
    fIntRule->SetOrder(ordervec);
}

void TPZSBFemVolume::SetElementGroupIndex(int64_t index)
{
    fElementGroupIndex = index;
    std::map<int64_t, int> globtolocal;
    TPZCompEl *celgr = Mesh()->Element(index);
    fElementGroup = celgr;
    int nc = celgr->NConnects();
    TPZManVector<int, 10> firsteq(nc + 1, 0);
    for (int ic = 0; ic < nc; ic++) {
        globtolocal[celgr->ConnectIndex(ic)] = ic;
        TPZConnect &c = celgr->Connect(ic);
        firsteq[ic + 1] = firsteq[ic] + c.NShape() * c.NState();
    }
    int neq = 0;
    TPZCompEl *celskeleton = Mesh()->Element(fSkeleton);
    nc = celskeleton->NConnects();
    for (int ic = 0; ic < nc; ic++) {
        TPZConnect &c = celskeleton->Connect(ic);
        neq += c.NShape() * c.NState();
    }
    fLocalIndices.Resize(neq);
    int count = 0;
    for (int ic = 0; ic < nc; ic++) {
        int64_t cindex = celskeleton->ConnectIndex(ic);
#ifdef PZDEBUG
        if (globtolocal.find(cindex) == globtolocal.end()) {
            DebugStop();
        }
#endif
        TPZConnect &c = celskeleton->Connect(ic);
        int neq = c.NShape() * c.NState();
        int locfirst = firsteq[globtolocal[cindex]];
        for (int eq = 0; eq < neq; eq++) {
            fLocalIndices[count++] = locfirst + eq;
        }
    }
#ifdef PZDEBUG
    if (count != neq) DebugStop();
#endif
}

void TPZSBFemVolume::InitMaterialData(TPZMaterialData &data)
{
    data.fShapeType = TPZMaterialData::EVecShape;
    data.gelElId = this->Reference()->Id();
    TPZManVector<TPZMaterialDataT<STATE>,2> datavec;
    auto *mat =
        dynamic_cast<TPZMatSingleSpaceT<STATE>*>(Material());
    auto *matMF = dynamic_cast<TPZMatCombinedSpacesT<STATE>*>(Material());
#ifdef PZDEBUG
    if (!mat && !matMF) {
        DebugStop();
    }
#endif
    if (mat) {
        mat->FillDataRequirements(data);
    } else {
        datavec.Resize(2);
        matMF->FillDataRequirements(datavec);
        data = datavec[0];
    }
    const int dim = this->Dimension();
    const int nshape = this->NShapeF();
    const int nstate = this->Material()->NStateVariables();
    data.fH1.fPhi.Redim(nstate * nshape*dim, 1);
    data.fH1.fDPhi.Redim(dim*nstate, nshape * nstate);
    data.dphix.Redim(dim*nstate, nshape * nstate);
    data.axes.Redim(dim, 3);
    data.jacobian.Redim(dim, dim);
    data.jacinv.Redim(dim, dim);
    data.x.Resize(3);
    if (mat && data.fNeedsSol)
    {
        uint64_t ulen,durow,ducol;
        mat->GetSolDimensions(ulen,durow,ducol);
        data.SetSolSizes(nstate, ulen, durow, ducol);
	}
}//void

void TPZSBFemVolume::ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X,
        TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,
        REAL &detjac, TPZFMatrix<REAL> &jacinv,
        TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<REAL> &dphidx) {
    TPZGeoEl * ref = this->Reference();
    if (!ref) {
        PZError << "\nERROR AT " << __PRETTY_FUNCTION__ << " - this->Reference() == NULL\n";
        return;
    }//if

    ref->Jacobian(intpoint, jacobian, axes, detjac, jacinv);
    this->Shape(intpoint, phi, dphidx);

}

void TPZSBFemVolume::BuildCornerConnectList(std::set<int64_t>& connectindexes) const {
    if (fSkeleton == -1) {
        DebugStop();
    }
    Mesh()->Element(fSkeleton)->BuildCornerConnectList(connectindexes);
}

void TPZSBFemVolume::PRefine(int order) {
    TPZCompEl *cel = Mesh()->Element(fSkeleton);
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
    intel->PRefine(order);
}

void TPZSBFemVolume::SetPreferredOrder(int order) {
    fPreferredOrder = order;
    TPZCompEl *cel = Mesh()->Element(fSkeleton);
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
    intel->SetPreferredOrder(order);
}

void TPZSBFemVolume::SetSkeleton(int64_t skeleton) {
#ifdef PZDEBUG
    if (fSkeleton != -1) {
        DebugStop();
    }
    if (fLocalIndices.size()) {
        DebugStop();
    }
#endif
    fSkeleton = skeleton;
    TPZCompEl *cel = Mesh()->Element(fSkeleton);
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *> (cel);
    TPZMultiphysicsElement *mpcel = dynamic_cast<TPZMultiphysicsElement *> (cel);
    if(mpcel)
    {
        intel = dynamic_cast<TPZInterpolationSpace *>(mpcel->Element(0));
    }
    if(!intel) DebugStop();
    int order = intel->GetPreferredOrder();
    SetIntegrationRule(2 * order);
}

void TPZSBFemVolume::LocalBodyForces(TPZFNMatrix<100,std::complex<double>> &f, TPZFNMatrix<100,std::complex<double>> &fbubble, TPZManVector<std::complex<double>> &eigval, TPZManVector<std::complex<double>> &eigvalbubbles, int icon)
{
    TPZCompMesh *cmesh = Mesh();
    TPZGeoEl *Ref2D = Reference();
    int matid = Ref2D->MaterialId();
    auto *mat2d = dynamic_cast<TPZMaterialT<STATE>*>(cmesh->FindMaterial(matid));    
    int dim2 = Ref2D->Dimension();
    
    TPZMaterial * material = Material();
    if (!material) {
        PZError << "TPZSBFemVolume::LocalBodyForces : no material for this element\n";
        Print(PZError);
        return;
    }
    if (dynamic_cast<TPZBndCond *> (material)) {
        std::cout << "TPZSBFemVolume::LocalBodyForces : null error - boundary condition material.";
        DebugStop();
    }
    int problemdimension = Mesh()->Dimension();
    if (Ref2D->Dimension() < problemdimension) return;
    
    TPZAutoPointer<TPZIntPoints> intrule = Ref2D->CreateSideIntegrationRule(Ref2D->NSides() - 1, 7);
    int maxIntOrder = intrule->GetMaxOrder();
    TPZManVector<int, 3> maxorder(Dimension(), maxIntOrder);
    intrule->SetOrder(maxorder);
    
    int nstate = material->NStateVariables();
    TPZManVector<REAL, 10> intpoint(problemdimension);
    REAL weight;
    
    TPZMaterialDataT<STATE> data;
    TPZMaterialDataT<STATE> data1d;
    this->InitMaterialData(data);
    int npts = intrule->NPoints();
    
    TPZManVector<REAL> bodyforce(nstate, 0.);
    TPZFMatrix<REAL> dbodyforce(nstate, dim2, 0.);
    
    TPZInterpolatedElement *CSkeleton = SkeletonElement();
    CSkeleton->InitMaterialData(data1d);
    
    TPZGeoEl *Ref1D = CSkeleton->Reference();
    
    int nphixi = fPhi.Rows();
    int numeig = fPhi.Cols();
    int numbubbles = fPhiBubble.Cols();
    
    TPZFNMatrix<200,std::complex<double>> eflocal(nphixi, numeig,0);
    TPZFNMatrix<200,std::complex<double>> eflocalbubble(nphixi, numbubbles, 0);
    
    for (int ipts=0; ipts<npts; ipts++) {
        intrule->Point(ipts, intpoint, weight);
        
        REAL sbfemparam = (1. - intpoint[dim2 - 1]) / 2.;
        if (sbfemparam < 0.) {
            std::cout << "sbfemparam " << sbfemparam << std::endl;
            sbfemparam = 0.;
        }
        if (IsZero(sbfemparam)) {
            for (int i = 0; i < dim2 - 1; i++) {
                intpoint[i] = 0.;
            }
            if (dim2 == 2) {
                sbfemparam = 1.e-6;
                intpoint[dim2 - 1] = 1. - 2.e-6;
            } else {
                sbfemparam = 1.e-4;
                intpoint[dim2 - 1] = 1. - 2.e-4;
            }
        }
        
        Ref2D->Jacobian(intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
        weight *= fabs(data.detjac);
        TPZManVector<REAL,3> intskel(problemdimension-1);
        for (int i=0; i<problemdimension-1; i++) {
            intskel[i] = intpoint[i];
        }
        CSkeleton->Shape(intskel, data1d.fH1.fPhi, data1d.fH1.fDPhi);
        
        Ref2D->X(intpoint, data.x);
        if (mat2d->HasForcingFunction())
        {
            mat2d->ForcingFunction()(data.x,bodyforce);
        }

        TPZManVector<std::complex<double>> eigvalb(fEigenvaluesBubble);
        
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            if (IsZero(-eigval[c] - 0.5*(dim2-2))){
                xiexp = 1;
            } else {
                xiexp = pow(sbfemparam, -eigval[c] - 0.5*(dim2-2) );
            }
            for (int i = 0; i < data1d.fH1.fPhi.Rows(); i++) {
                for (int istate = 0; istate < nstate; istate++) {
                    eflocal(i*nstate + istate, c) += xiexp * data1d.fH1.fPhi(i,0) * bodyforce[istate] * weight;
                }
            }
        }
        for (int c = 0; c < numbubbles; c++) {
            std::complex<double> xiexp;
            if (IsZero(eigvalbubbles[c])) {
                xiexp = 1;
            } else {
                xiexp = pow(sbfemparam, -eigvalbubbles[c]);
            }
            for (int i = 0; i < data1d.fH1.fPhi.Rows(); i++) {
                for (int istate = 0; istate < nstate; istate++) {
                    eflocalbubble(i*nstate + istate, c) += xiexp * data1d.fH1.fPhi(i,0) * bodyforce[istate] * weight;
                }
            }
        }
    }
    
    // f = \int \xi^0.5(d-1)*xi^\lambda
    for (int c = 0; c < numeig; c++) {
        for (int i = 0; i < nphixi; i++) {
            f(c,0) += eflocal(i,c) * fPhi(i,c);
        }
    }
    for (int c = 0; c < numbubbles; c++) {
        for (int i = 0; i < nphixi; i++) {
            fbubble(c,0) += eflocalbubble(i,c) * fPhiBubble(i,c);
        }
    }
    
// #ifdef PZ_LOG
//     if (loggerLBF->isDebugEnabled()) {
//         std::stringstream sout;
//         eflocal.Print("eflocal = ", sout, EMathematicaInput);
//         f.Print("f = ", sout, EMathematicaInput);
//         eflocalbubble.Print("eflocalbubble = ", sout, EMathematicaInput);
//         fbubble.Print("fbubble = ", sout, EMathematicaInput);
//         LOGPZ_DEBUG(loggerLBF, sout.str())
//     }
// #endif
}
