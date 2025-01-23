//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZMixedDarcyFlow.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
#include "TPZLapack.h"

#define USEBLAS

TPZMixedDarcyFlow::TPZMixedDarcyFlow() : TPZRegisterClassId(&TPZMixedDarcyFlow::ClassId),
                                         TBase(), fDim(-1) {}

[[maybe_unused]] TPZMixedDarcyFlow::TPZMixedDarcyFlow(int id, int dim) : TPZRegisterClassId(&TPZMixedDarcyFlow::ClassId),
                                                        TBase(id), fDim(dim)
{
}

/**
         copy constructor
 */
TPZMixedDarcyFlow::TPZMixedDarcyFlow(const TPZMixedDarcyFlow &copy) : TBase(copy), fDim(copy.fDim)
{
    
}
/**
         copy constructor
 */
TPZMixedDarcyFlow &TPZMixedDarcyFlow::operator=(const TPZMixedDarcyFlow &copy)
{
    TBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}


void TPZMixedDarcyFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {

    STATE force = 0;
    if (fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x, res);
        force = res[0];
    }
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<9, REAL> dphiPXY(3, dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);

    REAL &faceSize = datavec[0].HSize;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();

    int nactive = 0;
    for (const auto &i : datavec) {
        if (i.fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
    if (nactive >= 4) {
        if(nactive%2 != 0) DebugStop();
        int numavg = (nactive-2)/2;
        for(int iavg = 0; iavg<numavg; iavg++)
        {
            if(datavec[2+2*iavg].phi.Rows() != 1) DebugStop();
            if(datavec[2+2*iavg+1].phi.Rows() != 1) DebugStop();
        }
        if (phrp + phrq + 2*numavg != ek.Rows()) {
            DebugStop();
        }
    } else {
        if (phrp + phrq != ek.Rows()) {
            DebugStop();
        }
    }
#endif

#if defined(USEBLAS) && defined(USING_LAPACK)
    TPZFNMatrix<3, REAL> ivec(3, phrq, 0.);
    for (int iq = 0; iq < phrq; iq++){
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        for (int id = 0; id < 3; id++) {
            ivec(id, iq) = datavec[0].fDeformedDirections(id, ivecind);
        }
    }


    /**
       We cannot make the BLAS calls directly since
       the dimensions dont match:
       we are filling ek "blockwise".
       any ideas on how to use TPZFMatrix::MultAdd?
     **/
    {
        const double alpha = weight*inv_perm;
        constexpr double beta = 1.0;
        double *A, *B, *C;
        int m,n,k;
        m = phrq;
        n = phrq;
        k = 3;
        int LDA,LDB,LDC;
        LDA = 3;
        LDB = 3;
        LDC = ek.Rows();
        A = &ivec(0,0);
        B = &ivec(0,0);
        C = &ek(0,0);
        
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
        // ek = (weight/K) * (ivec^T ivec) + ek
    }
    
    //contribution matrix B
    if(phrp>0){
        double *A, *B, *C;
        const double alpha = - weight;
        constexpr double beta = 1.0;
        int m,n,k;
        m = phrq;
        n = phrp;
        k = 1;
        int LDA,LDB,LDC;
        LDA = phrq;
        LDB = phrp;
        LDC = ek.Rows();
        A = &divQ(0,0);
        B = &phip(0,0);
        C = &ek(0,phrq);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
        //ek = -weight* divQ phip^T + ek
    }
    //contribution matrix B^t
    if(phrp>0){
        double *A, *B, *C;
        const double alpha = - weight;
        constexpr double beta = 1.0;
        int m,n,k;
        m = phrp;
        n = phrq;
        k = 1;
        int LDA,LDB,LDC;
        LDA = phrp;
        LDB = phrq;
        LDC = ek.Rows();
        A = &phip(0,0);
        B = &divQ(0,0);
        C = &ek(phrq,0);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    m, n, k,alpha , A, LDA, B, LDB, beta, C, LDC);
    }

#else
    
    //Calculate the matrix contribution for flux. Matrix A
    for (int iq = 0; iq < phrq; iq++) {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = datavec[0].fDeformedDirections(id, ivecind);
        }

        TPZFNMatrix<3, REAL> ivecZ(3, 1, 0.);
        TPZFNMatrix<3, REAL> jvecZ(3, 1, 0.);
        for (int jq = 0; jq < phrq; jq++) {
            TPZFNMatrix<3, REAL> jvec(3, 1, 0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;

            for (int id = 0; id < 3; id++) {
                jvec(id, 0) = datavec[0].fDeformedDirections(id, jvecind);
            }

            //dot product between Kinv[u]v
            jvecZ.Zero();
            for (int id = 0; id < 3; id++) {
                jvecZ(id, 0) += inv_perm * jvec(id, 0);
            }
            REAL prod1 = ivec(0, 0) * jvecZ(0, 0) + ivec(1, 0) * jvecZ(1, 0) + ivec(2, 0) * jvecZ(2, 0);
            ek(iq, jq) += weight * phiQ(ishapeind, 0) * phiQ(jshapeind, 0) * prod1;
        }
    }

    // Coupling terms between flux and pressure. Matrix B
    for (int iq = 0; iq < phrq; iq++) {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;

        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = datavec[0].fDeformedDirections(id, ivecind);
            //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
            //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        }
        TPZFNMatrix<3, REAL> axesvec(3, 1, 0.);
        datavec[0].axes.Multiply(ivec, axesvec);

        REAL divwq = datavec[0].divphi(ivecind);

        for (int jp = 0; jp < phrp; jp++) {

            REAL fact = (-1.) * weight * phip(jp, 0) * divwq;
            // Matrix B
            ek(iq, phrq + jp) += fact;

            // Matrix B^T
            ek(phrq + jp, iq) += fact;

        }
    }
    
#endif

    // source term related to the pressure equation
    for (int ip = 0; ip < phrp; ip++) {
        ef(phrq + ip, 0) += (-1.) * weight * force * phip(ip, 0);
    }

    if (nactive >= 4) {
        if(nactive%2 != 0) DebugStop();
        int numav = (nactive-2)/2;
        for(int iav = 0; iav < numav; iav++)
        {
            for (int ip = 0; ip < phrp; ip++) {
                ek(phrq + ip, phrq + phrp + 2*iav) += phip(ip, 0) * weight;
                ek(phrq + phrp + 2*iav, phrq + ip) += phip(ip, 0) * weight;
            }
            ek(phrp + phrq + 2*iav + 1, phrq + phrp + 2*iav) += -weight;
            ek(phrq + phrp + 2*iav, phrp + phrq + 2*iav + 1) += -weight;
        }
    }
}

void TPZMixedDarcyFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else if (bc.Type() == 5) {
            v2 = res[0];
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;

        case 1 :            // Neumann condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += TPZMaterial::fBigNumber * v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {

                    ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            break;

        case 2 :            // mixed condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += weight / v1 * phiQ(iq, 0) * phiQ(jq, 0);
                }
            }
            break;

        case 4:
            //this case implemented the general Robin boundary condition
            // sigma.n = Km(u-u_D)+g
            //val1(0,0) = Km
            //val2(1,0) = g
            if (IsZero(bc.Val1()(0, 0))) {

                for (int iq = 0; iq < phrq; iq++) {
                    ef(iq, 0) += TPZMaterial::fBigNumber * normflux * phiQ(iq, 0) * weight;
                    for (int jq = 0; jq < phrq; jq++) {
                        ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    }
                }

            } else {

                REAL InvKm = 1. / bc.Val1()(0, 0);
                REAL g = normflux;
                for (int in = 0; in < phiQ.Rows(); in++) {
                    //<(InvKm g - u_D)*(v.n)
                    ef(in, 0) += (STATE) (InvKm * g - u_D) * phiQ(in, 0) * weight;
                    for (int jn = 0; jn < phiQ.Rows(); jn++) {
                        //InvKm(sigma.n)(v.n)
                        ek(in, jn) += (STATE) (InvKm * phiQ(in, 0) * phiQ(jn, 0) * weight);
                    }
                }
            }

            break;

        case 5:
            TPZFMatrix<REAL> &phi = datavec[0].fH1.fPhi;
            for (int in = 0; in < phi.Rows(); in++) {
                //<(InvKm g - u_D)*(v.n)
                ef(in, 0) += TPZMaterial::fBigNumber * v2 * phi(in, 0) * weight;
                for (int jn = 0; jn < phi.Rows(); jn++) {
                    //InvKm(sigma.n)(v.n)
                    ek(in, jn) += TPZMaterial::fBigNumber * phi(in, 0) * phi(jn, 0) * weight;
                }
            }

    }
}

void TPZMixedDarcyFlow::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) {
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    TPZManVector<STATE, 10> SolP, SolQ;
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;

    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    if(SolP.size() == 0) SolP.Resize(1,0.);

    if (var == 1) { //function (state variable Q)
        for (int i = 0; i < fDim; i++) {
            solOut[i] = datavec[0].sol[0][i];

        }
        return;
    }

    if (var == 2) {
        solOut[0] = SolP[0];//function (state variable p)
        return;
    }

    if (var == 3) {
        solOut[0] = datavec[0].dsol[0](0, 0);
        solOut[1] = datavec[0].dsol[0](1, 0);
        solOut[2] = datavec[0].dsol[0](2, 0);
        return;
    }

    if (var == 4) {
        solOut[0] = datavec[0].dsol[0](0, 1);
        solOut[1] = datavec[0].dsol[0](1, 1);
        solOut[2] = datavec[0].dsol[0](2, 1);
        return;
    }

    if (var == 5) {
        solOut[0] = datavec[0].dsol[0](0, 0) + datavec[0].dsol[0](1, 1);
        return;
    }

    // Exact solution
    if (var == 6) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, flux);
        }
        solOut[0] = exactSol[0];
        return;
    } // var6

    if (var == 7) {

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, gradu);
        }

        for (int i = 0; i < fDim; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    } // var7

    if (var == 8) {
        solOut[0] = datavec[1].p;
        return;
    }

    if (var == 9) {

        if(datavec[1].fShapeType == TPZMaterialData::EEmpty) return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1.,0.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1,0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < fDim; i++) {
            solOut[i] = dsoldx(i, 0);
        }

        return;
    }

    if (var == 10) {
        solOut[0] = 0.;
        // solOut[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        for (int j = 0; j < fDim; j++) {
            solOut[0] += datavec[0].dsol[0](j, j);
        }
        return;
    }

    if (var == 11) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        fExactSol(datavec[0].x, exactSol, flux);
        solOut[0] = flux(2, 0);
        return;
    }
    if (var == 12) {
        for (int i = 0; i < fDim; i++) {
            solOut[i] = 0.;
        }
        for (int i = 0; i < fDim; i++) {
            solOut[i] -= inv_perm * datavec[0].sol[0][i];
        }
        return;
    }
    if (var == 13) {
        solOut[0] = perm;
        return;
    }

    if (datavec.size() == 4) {
        if (var == 14) {
            solOut[0] = datavec[2].sol[0][0];
            return;
        }
        if (var == 15) {
            solOut[0] = datavec[3].sol[0][0];
            return;
        }

    }

    if (var == 16) { //ExactFluxShiftedOrigin
        // Solution EArcTan returns NAN for (x,y) == (0,0). Replacing data.x by
        // inf solves this problem.
        STATE infinitesimal = 0.0000000001;
        TPZManVector<REAL, 3> inf = {infinitesimal, infinitesimal, infinitesimal};

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            if (datavec[0].x[0] == 0. && datavec[0].x[1] == 0.) {
                fExactSol(inf, exactSol, gradu);
            } else {
                fExactSol(datavec[0].x, exactSol, gradu);
            }
        }
        for (int i = 0; i < 3; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    }

    if (var == 17) {
        TPZVec<STATE> divsigma(1, 0.);
        if (fForcingFunction) {
            fForcingFunction(datavec[0].x, divsigma);
        }
        solOut[0] = divsigma[0];
        return;
    }
}

void TPZMixedDarcyFlow::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

    /**
     * datavec[0]= Flux
     * datavec[1]= Pressure
     *
     * Errors:
     * [0] L2 for pressure
     * [1] L2 for flux
     * [2] L2 for div(flux)
     * [3] Grad pressure (Semi H1)
     * [4] Hdiv norm
    **/
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1,0);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZManVector<STATE,1> divsigma(1,0.);

    TPZManVector<STATE,1> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    if(data[1].sol[0].size())
        pressurefem[0] = data[1].sol[0][0];

    const STATE perm = GetPermeability(data[0].x);
    const STATE inv_perm = 1 / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexact(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexact[i] = -perm * gradpressure[i];
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] - fluxexact[i]) * inv_perm * (fluxfem[i] - fluxexact[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
}

int TPZMixedDarcyFlow::VariableIndex(const std::string &name) const {
    if (!strcmp("Flux", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 2;
    if (!strcmp("GradFluxX", name.c_str())) return 3;
    if (!strcmp("GradFluxY", name.c_str())) return 4;
    if (!strcmp("DivFlux", name.c_str())) return 5;
    if (!strcmp("ExactPressure", name.c_str())) return 6;
    if (!strcmp("ExactFlux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("GradPressure", name.c_str())) return 9;
    if (!strcmp("Divergence", name.c_str())) return 10;
    if (!strcmp("ExactDiv", name.c_str())) return 11;
    if (!strcmp("Derivative", name.c_str())) return 12;
    if (!strcmp("Permeability", name.c_str())) return 13;
    if (!strcmp("g_average", name.c_str())) return 14;
    if (!strcmp("u_average", name.c_str())) return 15;
    if (!strcmp("ExactFluxShiftedOrigin", name.c_str())) return 16;
    if (!strcmp("EstimatedError", name.c_str())) return 100;
    if (!strcmp("TrueError", name.c_str())) return 101;
    if (!strcmp("EffectivityIndex", name.c_str())) return 102;
    if (!strcmp("ExactDivSigma", name.c_str())) return 17;
    DebugStop();
    return -1;
}

int TPZMixedDarcyFlow::NSolutionVariables(int var) const {
    if (var == 1) return 3;
    if (var == 2) return 1;
    if (var == 3) return 3;
    if (var == 4) return 3;
    if (var == 5) return 1;
    if (var == 6) return 1;
    if (var == 7) return 3;
    if (var == 8) return 1;
    if (var == 9) return 3;
    if (var == 10 || var == 11) return 1;
    if (var == 12) return 3;
    if (var == 13) return 1;
    if (var == 14) return 1;
    if (var == 15) return 1;
    if (var == 16) return 3;
    if (var == 17) return 1;
    if (var == 100) return 1;
    if (var == 101) return 1;
    if (var == 102) return 1;

    DebugStop();
    return -1;
}

void TPZMixedDarcyFlow::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

int TPZMixedDarcyFlow::ClassId() const {
    return Hash("TPZMixedDarcyFlow") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZMixedDarcyFlow::NewMaterial() const {
    return new TPZMixedDarcyFlow(*this);
}

void TPZMixedDarcyFlow::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZMixedDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
    }
}

void
TPZMixedDarcyFlow::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    // default is no specific data requirements
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++) {
        datavec[iref].SetAllRequirements(false);
        datavec[iref].fNeedsSol = false;
    }
    datavec[0].fNeedsNormal = true;
    if (type == 50) {
        for (int iref = 0; iref < nref; iref++) {
            datavec[iref].fNeedsSol = false;
        }
    }
}
#undef USEBLAS
