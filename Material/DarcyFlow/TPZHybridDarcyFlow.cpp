//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZHybridDarcyFlow.h"
#include "pzaxestools.h"
#include "TPZLapack.h"

TPZHybridDarcyFlow::TPZHybridDarcyFlow() : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                               TPZMatCombinedSpacesT<STATE>(), TPZDarcyFlow() {}

TPZHybridDarcyFlow::TPZHybridDarcyFlow(int id, int dim) : TPZRegisterClassId(&TPZHybridDarcyFlow::ClassId),
                                        TPZDarcyFlow(id,dim) {}



/*
int TPZHybridDarcyFlow::VariableIndex(const std::string &name) const {

    if (!strcmp("Solution", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 1;
    if (!strcmp("Derivative", name.c_str())) return 2;
    if (!strcmp("GradU", name.c_str())) return 2;
    if (!strcmp("KDuDx", name.c_str())) return 3;
    if (!strcmp("KDuDy", name.c_str())) return 4;
    if (!strcmp("KDuDz", name.c_str())) return 5;
    if (!strcmp("NormKDu", name.c_str())) return 6;
    if (!strcmp("MinusKGradU", name.c_str())) return 7;
    if (!strcmp("Flux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("ExactPressure", name.c_str())) return 9;
    if (!strcmp("ExactSolution", name.c_str())) return 9;
    if (!strcmp("ExactFlux", name.c_str())) return 10;
    if (!strcmp("Div", name.c_str())) return 11;
    if (!strcmp("Divergence", name.c_str())) return 11;
    if (!strcmp("ExactDiv", name.c_str())) return 12;
    if (!strcmp("ExactDivergence", name.c_str())) return 12;
    if (!strcmp("FluxL2", name.c_str())) return 13;

    return TPZDarcyFlow::VariableIndex(name);
}

int TPZHybridDarcyFlow::NSolutionVariables(int var) const {

    if (var == 1) return 1;      // Solution/Pressure
    if (var == 2) return fDim;   // Derivative/GradU
    if (var == 3) return 1;      // KDuDx;
    if (var == 4) return 1;      // KDuDy;
    if (var == 5) return 1;      // KDuDz;
    if (var == 6) return 1;      // NormKDu;
    if (var == 7) return fDim;   // MinusKGradU/Flux;
    if (var == 8) return 1;      // POrder
    if (var == 9) return 1;      // ExactPressure/ExactSolution
    if (var == 10) return fDim;  // ExactFlux
    if (var == 11) return 1;     // Div/Divergence
    if (var == 12) return 1;     // ExactDiv/ExactDivergence
    if (var == 13) return fDim;  // FluxL2

    return TPZDarcyFlow::NSolutionVariables(var);
}
*/


int TPZHybridDarcyFlow::ClassId() const {
    return Hash("TPZHybridDarcyFlow") ^ TPZDarcyFlow::ClassId() << 1;
}

TPZMaterial *TPZHybridDarcyFlow::NewMaterial() const {
    return new TPZHybridDarcyFlow(*this);
}

void TPZHybridDarcyFlow::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << TPZDarcyFlow::Id() << "\n";
    out << "Dimension: " << TPZDarcyFlow::Dimension() << "\n\n";
}

/** @name Contribute */
/** @{ */
/**
 * @brief It computes a contribution to the stiffness matrix
 * and load vector at one integration point.
 * @param[in] datavec stores all input data
 * @param[in] weight is the weight of the integration rule
 * @param[out] ek is the element matrix
 * @param[out] ef is the rhs vector
 */
void TPZHybridDarcyFlow::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                        REAL weight,TPZFMatrix<STATE> &ek,
                        TPZFMatrix<STATE> &ef)
{
    /**
    datavec[1] L2 mesh (phi's)
    datavec[0] Hdiv mesh,
    datavec[2] Interface Mesh
    datavec[3] Interface Mesh

    Implement the matrix
    |Sk Ck^T |  = |f1|
    |Ck  0   |    |f2|
    Sk = int_K K graduk.gradv dx = int_K K gradphi_i.gradphi_j dx
    CK = int_partialK lambda_k*uk dx = int_K phi_i dx
    f1 = int_K f*v dx = int_K f*phi_j dx
    ck = int_partialK phi_i*mu_j dx
    f2 = int_partialK g*mu_j dx
    **/

    TPZFMatrix<REAL>  &phi = datavec[1].phi;
    TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    TPZVec<REAL>  &x = datavec[1].x;

    int phr = phi.Rows();

    STATE fXfLoc = 0;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        fForcingFunction(x, res);
        fXfLoc = res[0];
    }

    STATE KPerm = GetPermeability(datavec[0].x);

#if defined(USING_MKL)
    {
        double *A, *B, *C;
        int m, n, k;
        double alpha, beta;
        m = phr, k = fDim, n = phr;
        alpha = weight*KPerm;
        beta = 1.0;
        int LDA,LDB,LDC;
        LDA = dphi.Rows();
        LDB = dphi.Rows();
        LDC = ek.Rows();
        if(LDC != phr+2 && datavec.size() ==4)
            DebugStop();
        A = &dphi(0,0);
        B = &dphi(0,0);
        C = &ek(0,0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                       phr, phr, LDA,alpha , A, LDA, B, LDB, beta, C, LDC);
    }
    {
        //saxpy implementation
        int N = phr;
        double alpha = weight*fXfLoc;
        double *X = &phi(0,0);
        int incX = 1;
        double *Y = &ef(0,0);
        int incY = 1;
        cblas_daxpy(N,alpha,X,incX, Y,incY);
    }
    if(datavec.size() >=4){
        int N = phr;
        double alpha = weight;
        double *X = &phi(0,0);
        int incX = 1;
        double *Y = &ek(0,phr);
        int incY = 1;
        cblas_daxpy(N,alpha,X,incX, Y,incY);
        Y = &ek(phr,0);
        incY = ek.Rows();
        cblas_daxpy(N,alpha,X,incX, Y,incY);
    }
#else
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);

        //matrix Sk
        for( int jn = 0; jn < phr; jn++ ) {
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    if(datavec.size() >=4){
        for (int in =0; in < phr; in++) {
            ek(phr,in) += weight*phi(in,0);//lambda*phi
            ek(in,phr) += weight*phi(in,0);
        }
    }
#endif
    //equacoes de restricao de pressao media
    if(datavec.size() >=4) {
        ek(phr, phr + 1) -= weight;
        ek(phr + 1, phr) -= weight;
    }
}

void TPZHybridDarcyFlow::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCondT<STATE> &bc)
{
    TPZFMatrix<REAL>  &phi_u = datavec[1].phi;
    TPZFMatrix<REAL>  &phi_flux = datavec[0].phi;
    //    TPZFMatrix<REAL> &axes = data.axes;
    int phr_primal = phi_u.Rows();
    int phr_hybrid = phi_flux.Rows();
    bool primal = true;   /// weather pressure or flux is hybridized
    TPZManVector<REAL,3> x(3);
    if(phr_hybrid)
    {
        primal = false;
        x = datavec[0].x;
    }
    else
    {
        x = datavec[1].x;
    }
    short in,jn;
    STATE v2[1];
    v2[0] = bc.Val2()[0];

    if(bc.HasForcingFunctionBC()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
        TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
        bc.ForcingFunctionBC()(x, res, dres);
        v2[0] = res[0];
    }

    if(primal)
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += (STATE)(fBigNumber* phi_u(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_primal; jn++) {
                        ek(in,jn) += fBigNumber * phi_u(in,0) * phi_u(jn,0) * weight;
                    }
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in,0) += v2[0] * (STATE)(phi_u(in,0) * weight);
                }
                break;
            case 2 :        // mixed condition
                for(in = 0 ; in < phr_primal; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_u(in, 0) * weight);
                    for (jn = 0 ; jn < phi_u.Rows(); jn++) {
                        ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi_u(in,0) * phi_u(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
            default:
                DebugStop();
        }
    } else
    {
        switch (bc.Type()) {
            case 0 :            // Dirichlet condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in,0) += v2[0] * (STATE)(phi_flux(in,0) * weight);
                }
                break;
            case 1 :            // Neumann condition
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in,0) += (STATE)(fBigNumber* phi_flux(in,0) * weight) * v2[0];
                    for (jn = 0 ; jn < phr_hybrid; jn++) {
                        ek(in,jn) += fBigNumber * phi_flux(in,0) * phi_flux(jn,0) * weight;
                    }
                }
                break;
            case 2 :        // mixed condition
                DebugStop();
                for(in = 0 ; in < phr_hybrid; in++) {
                    ef(in, 0) += v2[0] * (STATE)(phi_flux(in, 0) * weight);
                    for (jn = 0 ; jn < phi_flux.Rows(); jn++) {
                        ek(in,jn) += 1./bc.Val1()(0,0) * (STATE)(phi_flux(in,0) * phi_flux(jn,0) * weight);     // peso de contorno => integral de contorno
                    }
                }
                break;
        }
    }
}

void TPZHybridDarcyFlow::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    if(!fExactSol) return;

    TPZDarcyFlow::Errors(data[1],errors);
    return;
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE> u_exact(1);
    TPZFNMatrix<9,STATE> du_exact;


    if(this->fExactSol){

        this->fExactSol(data[1].x,u_exact,du_exact);
    }

    REAL pressure = data[1].sol[0][0];

    // errors[0] norm L2 || u ||_l2

    errors[0] = (pressure-u_exact[0])*(pressure-u_exact[0]);//exact error pressure

    // errors[1] Semi norm H1 || grad u ||_l2

    TPZManVector<STATE,3> sol(1),dsol(3,0.);

    TPZFMatrix<REAL> &dsolaxes = data[1].dsol[0];
    TPZFNMatrix<9,REAL> flux(3,0);
    TPZAxesTools<REAL>::Axes2XYZ(dsolaxes, flux, data[1].axes);

    for(int id=0; id<fDim; id++) {
        REAL diff = fabs(flux(id,0) - du_exact(id,0));
        errors[1]  += diff*diff;
    }

    // error[2] H1 norm

    errors[2] = errors[0] +errors[1];

    // error[3] Energy norm || u ||_e = a(u,u)= int_K K gradu.gradu dx

    STATE KPerm = GetPermeability(data[0].x);


    TPZFNMatrix<9,REAL> gradpressure(fDim,1),Kgradu(fDim,1);
    for (int i=0; i<fDim; i++) {
        gradpressure(i,0) = du_exact(i,0);
        Kgradu(i,0) = gradpressure(0)*KPerm;
    }

    REAL energy = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
            double cperm =0.;
            if(i==j)
                cperm = KPerm;
            energy += cperm*fabs(flux(j,0) - du_exact(j,0))*fabs(flux(i,0) - du_exact(i,0));
        }
    }

    errors[3] = energy;
}
/**@}*/
