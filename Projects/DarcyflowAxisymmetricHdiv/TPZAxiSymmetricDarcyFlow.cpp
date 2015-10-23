/*
 *  TPZAxiSymmetricDarcyFlow.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZAxiSymmetricDarcyFlow.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzfmatrix.h"

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow() : TPZDiscontinuousGalerkin()
{
    
    fepsilon = 10.0e-7;
    
    fSimulationData=NULL;
    fReservoirdata=NULL;
    fPetrophysicdata=NULL;
    fluid_alpha=NULL;
    fluid_beta=NULL;
    fluid_gamma=NULL;
    
    fnstate_vars = 0;
    fBulkVelocity.Resize(3,0.0);
    fAveragePressure = 0.0;
    fWaterSaturation = 0.0;
    fOilSaturation = 0.0;
    
    fWaterVelocity.Resize(3,0.0);
    fOilVelocity.Resize(3,0.0);
    fGasVelocity.Resize(3,0.0);
    
    fWaterPressure.Resize(4,0.0);
    fOilPressure.Resize(4,0.0);
    fGasPressure.Resize(4,0.0);
    
    fWaterDensity.Resize(4,0.0);
    fOilDensity.Resize(4,0.0);
    fGasDensity.Resize(4,0.0);
    
    fFWater.Resize(4,0.0);
    fFOil.Resize(4,0.0);
    fFGas.Resize(4,0.0);
    
    fWaterMobility.Resize(4,0.0);
    fOilMobility.Resize(4,0.0);
    fGasMobility.Resize(4,0.0);
    
    fTotalMobility.Resize(4,0.0);
    fTotalMobility.Resize(4,0.0);
    
}

TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(int matid) : TPZDiscontinuousGalerkin(matid)
{
    fepsilon = 10.0e-7;
    
    fSimulationData=NULL;
    fReservoirdata=NULL;
    fPetrophysicdata=NULL;
    
    fBulkVelocity.Resize(3,0.0);
    fAveragePressure = 0.0;
    fWaterSaturation = 0.0;
    fOilSaturation = 0.0;
    
    fWaterVelocity.Resize(3,0.0);
    fOilVelocity.Resize(3,0.0);
    fGasVelocity.Resize(3,0.0);
    
    fWaterPressure.Resize(4,0.0);
    fOilPressure.Resize(4,0.0);
    fGasPressure.Resize(4,0.0);
    
    fWaterDensity.Resize(4,0.0);
    fOilDensity.Resize(4,0.0);
    fGasDensity.Resize(4,0.0);
    
    fFWater.Resize(4,0.0);
    fFOil.Resize(4,0.0);
    fFGas.Resize(4,0.0);
    
    fWaterMobility.Resize(4,0.0);
    fOilMobility.Resize(4,0.0);
    fGasMobility.Resize(4,0.0);
    
    fTotalMobility.Resize(4,0.0);
    fTotalMobility.Resize(4,0.0);
}


TPZAxiSymmetricDarcyFlow::TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat) : TPZDiscontinuousGalerkin(mat)
{
    fSimulationData     = mat.fSimulationData;
    fReservoirdata      = mat.fReservoirdata;
    fPetrophysicdata    = mat.fPetrophysicdata;
    
    fBulkVelocity       = mat.fBulkVelocity;
    fAveragePressure    = mat.fAveragePressure;
    fWaterSaturation    = mat.fWaterSaturation;
    fOilSaturation      = mat.fOilSaturation;
    
    fWaterVelocity      = mat.fWaterVelocity;
    fOilVelocity        = mat.fOilVelocity;
    fGasVelocity        = mat.fGasVelocity;
    
    fWaterPressure      = mat.fWaterPressure;
    fOilPressure        = mat.fOilPressure;
    fGasPressure        = mat.fGasPressure;
    
    fWaterDensity       = mat.fWaterDensity;
    fOilDensity         = mat.fOilDensity;
    fGasDensity         = mat.fGasDensity;
    
    fFWater             = mat.fFWater;
    fFOil               = mat.fFOil;
    fFGas               = mat.fFGas;
    
    fWaterMobility      = mat.fWaterMobility;
    fOilMobility        = mat.fOilMobility;
    fGasMobility        = mat.fGasMobility;
    
    fTotalMobility      = mat.fTotalMobility;
    fTotalDensity       = mat.fTotalDensity;
}

TPZAxiSymmetricDarcyFlow::~TPZAxiSymmetricDarcyFlow()
{
    
}

void TPZAxiSymmetricDarcyFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

void TPZAxiSymmetricDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

void TPZAxiSymmetricDarcyFlow::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

int TPZAxiSymmetricDarcyFlow::VariableIndex(const std::string &name) {
    if (!strcmp("WeightedPressure", name.c_str())) return 0;
    if (!strcmp("BulkVelocity", name.c_str())) return 1;
    if (!strcmp("Salpha", name.c_str())) return 2;
    if (!strcmp("Sbeta", name.c_str())) return 3;
    if (!strcmp("Sgamma", name.c_str())) return 4;
    if (!strcmp("Rhoalpha", name.c_str())) return 5;
    if (!strcmp("Rhobeta", name.c_str())) return 6;
    if (!strcmp("Rhogamma", name.c_str())) return 7;
    if (!strcmp("Porosity", name.c_str())) return 8;
    if (!strcmp("DivOfBulkVeclocity", name.c_str())) return 9;
    if (!strcmp("ExactSaturation", name.c_str())) return 10;
    if (!strcmp("ExactP", name.c_str())) return 11;
    if (!strcmp("Rhs", name.c_str())) return 12;
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

int TPZAxiSymmetricDarcyFlow::NSolutionVariables(int var) {
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 3; // Vector
        case 2:
            return 1; // Scalar
        case 3:
            return 1; // Scalar
        case 4:
            return 1; // Scalar
        case 5:
            return 1; // Scalar
        case 6:
            return 1; // Scalar
        case 7:
            return 1; // Scalar
        case 8:
            return 1; // Scalar
        case 9:
            return 1; // Scalar
        case 10:
            return 1; // Scalar
        case 11:
            return 1; // Scalar
        case 12:
            return 1; // Scalar
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

void TPZAxiSymmetricDarcyFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    int Qblock = 0;
    int Pblock = 1;
    int Swblock = 2;
    int Soblock = 3;
    
    TPZManVector<REAL,3> Q = datavec[Qblock].sol[0];
    REAL P = datavec[Pblock].sol[0][0];
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    TPZManVector<REAL> state_vars(2,0.0);
    state_vars[1] = P;
    
    REAL S_alpha = 0.0;
    REAL S_beta = 0.0;
    REAL S_gamma = 0.0;
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    
    if (fSimulationData->IsTwoPhaseQ()) {
        state_vars.Resize(3);
        
        fluid_beta->Density(rho_beta, state_vars);
        S_alpha = datavec[Swblock].sol[0][0];
        state_vars[2] = S_alpha;
        S_beta = 1.0 - S_alpha;
    }
    
    if (fSimulationData->IsThreePhaseQ()) {
        state_vars.Resize(4);
        
        fluid_beta->Density(rho_beta, state_vars);
        fluid_gamma->Density(rho_gamma, state_vars);
        S_alpha = datavec[Swblock].sol[0][0];
        S_beta = datavec[Soblock].sol[0][0];
        state_vars[2] = S_alpha;
        state_vars[3] = S_beta;
    }
    
    
    S_gamma = 1.0 - S_alpha - S_beta;
    
    
    
    
    REAL time;
    TPZVec<STATE> Saturation(1,0.0);
    TPZFMatrix<STATE> GradS(1,0);
    
    // Petrophysics data
    
    //    REAL Po, Pw, Pg;
    //    REAL pcow, pcgo, pcgw;
    //    REAL dPcowdSw, dPcgodSo, dPcgwdSw;
    
    
    //    this->fPetrophysicdata->Pcwo(Sw, pcow, dPcowdSw);
    //    this->fPetrophysicdata->Pcog(So, pcgo, dPcgodSo);
    //    this->fPetrophysicdata->Pcgo(Sw, pcgw, dPcgwdSw);
    
    // Recovering phase pressures
    //    Po = P + Sw * pcow - Sg * pcgo;
    //    Pw = P - So * pcow - Sg * pcgw;
    //    Pg = P + So * pcgo + So * pcgw;
    
    // Rock and fluids parameters
    
    //    REAL rockporosity, waterdensity, oildensity;
    //    REAL drockporositydP, dwaterdensitydPw, doildensitydPo;
    
    REAL rockporosity;
    REAL drockporositydP;
    this->fReservoirdata->Porosity(P, rockporosity, drockporositydP);
    //    this->fFluidmodeldata->WaterDensity(Pw, waterdensity, dwaterdensitydPw);
    //    this->fFluidmodeldata->OilDensity(Po, oildensity, doildensitydPo);
    
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        case 0:
        {
            Solout[0] = P;
        }
            break;
        case 1:
        {
            Solout[0] = Q[0]; // Bulk mass velocity
            Solout[1] = Q[1]; // Bulk mass velocity
        }
            break;
        case 2:
        {
            Solout[0] = S_alpha;
        }
            break;
        case 3:
        {
            Solout[0] = S_beta;
        }
            break;
        case 4:
        {
            Solout[0] = S_gamma;
        }
            break;
        case 5:
        {
            Solout[0] = rho_alpha[0];
        }
            break;
        case 6:
        {
            Solout[0] = rho_beta[0];
        }
            break;
        case 7:
        {
            Solout[0] = rho_gamma[0];
        }
            break;
        case 8:
        {
            Solout[0] = rockporosity;
        }
            break;
        case 9:
        {
            Solout[0] = dQdx(0,0) + dQdx(1,1) + dQdx(2,2);
        }
            break;
        case 10:
        {
            time=this->fSimulationData->GetTime();
            fTimedependentFunctionExact->Execute(datavec[Qblock].x, time, Saturation,GradS);
            Solout[0] = Saturation[0];
            
        }
            break;
        case 11:
        {
            //            REAL epsilon = 0.01;
            //            REAL xc = 0.5;
            //            REAL yc = 0.5;
            REAL x = datavec[Pblock].x[0];
            REAL y = datavec[Pblock].x[1];
            //            REAL Pressure = exp(-((x-xc)*(x-xc)+(y-yc)*(y-yc))/epsilon);
            
            REAL Pressure = x*y*(1-x)*(1-y)*(sin(M_PI*x))*(cos(M_PI*y));
            Solout[0] = Pressure;
            
        }
            break;
        case 12:
        {
            TPZManVector<STATE,1> fvalue(1,0.0);
            TPZFMatrix<REAL> Grad;
            if(fTimeDependentForcingFunction)
            {
                fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
            }
            Solout[0] = fvalue[0];
        }
            break;
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

// Divergence on deformed element
void TPZAxiSymmetricDarcyFlow::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<REAL> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    
    TPZFNMatrix<660,REAL> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<REAL> Qaxes = datavec[ublock].axes;
    TPZFMatrix<REAL> QaxesT;
    TPZFMatrix<REAL> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<REAL> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<REAL> GradOfX;
    TPZFMatrix<REAL> GradOfXInverse;
    TPZFMatrix<REAL> VectorOnMaster;
    TPZFMatrix<REAL> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola)
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fNormalVec(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fNormalVec(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fNormalVec(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            DivergenceofPhi(iq,0) =  (1.0/JacobianDet) * ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) + dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) );
        }
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) + datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex);
        }
    }
    
    return;
    
}

// Jacobian contribution
void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    this->ContributeDarcy(datavec, weight, ek, ef);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeAlpha(datavec, weight, ek, ef);
    }
    
    return;
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2        = datavec[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2        = datavec[Soblock].phi; // For L2  test functions So
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions
    
    
    TPZFMatrix<REAL> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSwL2    = phiSwL2.Rows();                                   // For L2   Sw
    int nphiSoL2    = phiSoL2.Rows();                                   // For L2   So
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSw   = nphiPL2       + iniP;
    int iniSo   = nphiSwL2      + iniSw;
    
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL Sw             = datavec[Swblock].sol[0][0];
    REAL So             = datavec[Soblock].sol[0][0];
    //    REAL Sg             = 1.0 - So - Sw;
    
    
    TPZFMatrix<REAL> Graduaxes = datavec[ublock].dsol[0]; // Piola divengence may works, needed set piola computation on the solution elchiv method!!!
    TPZFMatrix<REAL> GradPaxes = datavec[Pblock].dsol[0];
    
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    REAL deltat = fSimulationData->GetDeltaT();
    this->UpdateStateVariables(u, P, Sw, So);
    this->PhaseFractionalFlows();
    
    // Rock and fluids parameters
    TPZFMatrix<REAL> KInverse = fReservoirdata->KabsoluteInv();
    REAL rockporosity, drockporositydP;
    this->fReservoirdata->Porosity(fAveragePressure, rockporosity, drockporositydP);
    
    
    // Defining local variables
    TPZFMatrix<REAL> oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> oneoverlambda_Kinv_jphiuHdiv(2,1);
    TPZFMatrix<REAL> Gravity(2,1);
    
    oneoverlambda_Kinv_u(0,0) = (1.0/fTotalMobility[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    oneoverlambda_Kinv_u(1,0) = (1.0/fTotalMobility[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    Gravity = fSimulationData->GetGravity();
    
    
    //    REAL divu = 0.0;
    TPZFMatrix<REAL> iphiuHdiv(2,1);
    TPZFMatrix<REAL> jphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    int jshapeindex;
    int jvectorindex;
    
    REAL dgmcdP = (fFOil[0] * fOilDensity[1] + fFWater[1]* fWaterDensity[0] + fFOil[1] * fOilDensity[0] + fFWater[0]* fWaterDensity[1]);
    REAL dgmcdSw = (fFOil[0] * fOilDensity[2] + fFWater[2]* fWaterDensity[0] + fFOil[2] * fOilDensity[0] + fFWater[0]* fWaterDensity[2]);
    REAL dgmcdSo = (fFOil[0] * fOilDensity[3] + fFWater[3]* fWaterDensity[0] + fFOil[3] * fOilDensity[0] + fFWater[0]* fWaterDensity[3]);
    
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
            
            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
            
            oneoverlambda_Kinv_jphiuHdiv(0,0) = (1.0/fTotalMobility[0]) * (KInverse(0,0)*jphiuHdiv(0,0) + KInverse(0,1)*jphiuHdiv(1,0));
            oneoverlambda_Kinv_jphiuHdiv(1,0) = (1.0/fTotalMobility[0]) * (KInverse(1,0)*jphiuHdiv(0,0) + KInverse(1,1)*jphiuHdiv(1,0));
            
            ek(iq + iniu,jq + iniu) +=  weight * ((oneoverlambda_Kinv_jphiuHdiv(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_jphiuHdiv(1,0)*iphiuHdiv(1,0)));
        }
        
        // du/dalphap terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iq + iniu,jp + iniP) += weight * ((-(fTotalMobility[1]/fTotalMobility[0]))*((oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))) - DivergenceOnDeformed(iq,0) - dgmcdP *(Gravity(0,0)*iphiuHdiv(0,0) + Gravity(1,0)*iphiuHdiv(1,0)) )* phiPL2(jp,0) ;
        }
        
        // du/dalphaSw terms
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(iq + iniu, jsw + iniSw) += weight * (-(fTotalMobility[2]/fTotalMobility[0])* phiSwL2(jsw,0) *((oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))) - dgmcdSw *(Gravity(0,0)*iphiuHdiv(0,0) + Gravity(1,0)*iphiuHdiv(1,0)) * phiSwL2(jsw,0));
        }
        
        // du/dalphaSo terms
        for (int jso = 0; jso < nphiSoL2; jso++)
        {
            ek(iq + iniu, jso + iniSo) += weight * (-(fTotalMobility[3]/fTotalMobility[0])* phiSoL2(jso,0) *((oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0)) ) - dgmcdSo *(Gravity(0,0)*iphiuHdiv(0,0) + Gravity(1,0)*iphiuHdiv(1,0)) * phiSoL2(jso,0) );
        }
    }
    
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            ek(ip + iniP, jq + iniu) += -1.0 * weight * DivergenceOnDeformed(jq,0) * phiPL2(ip,0);
        }
        
        //- (1.0/deltat) * rockporosity * (Sw * fWaterDensity[0] +  So * fOilDensity[0] )
        
        // dp/dalphap terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(ip + iniP,jp + iniP) += weight * ( (-1.0/deltat) * rockporosity * (Sw * fWaterDensity[1] +  So * fOilDensity[1] ) )* phiPL2(jp,0) * phiPL2(ip,0) ;
        }
        
        // dp/dalphaSw terms
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(ip + iniP, jsw  + iniSw) += weight * ( (-1.0/deltat) * rockporosity * (fWaterDensity[0])) * phiSwL2(jsw,0) * phiPL2(ip,0);
        }
        
        // dp/dalphaSo terms
        for (int jso = 0; jso < nphiSoL2; jso++)
        {
            ek(ip + iniP, jso  + iniSo) += weight * ( (-1.0/deltat) * rockporosity * (fOilDensity[0])) * phiSwL2(jso,0) * phiPL2(ip,0);
        }
        
    }
    
    // Transport equations
    //  n+1 state computations
    for (int isw = 0; isw < nphiSwL2; isw++)
    {
        
        //      ef(isw  + iniSw ) += weight * (1.0/deltat) * rockporosity * Sw * fWaterDensity[0] * phiSwL2(isw,0);
        
        // dWaterDensity/dalphaP terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(isw  + iniSw, jp  + iniP) += weight * (1.0/deltat) * rockporosity * Sw * fWaterDensity[1] * phiPL2(jp,0) * phiSwL2(isw,0);
        }
        
        // dSw/dalphaSw terms
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(isw  + iniSw, jsw  + iniSw) += weight * (1.0/deltat) * rockporosity * phiSwL2(jsw,0) * fWaterDensity[0] * phiSwL2(isw,0);
        }
        
    }
    
    for (int iso = 0; iso < nphiSoL2; iso++)
    {
        
        //      ef(iso  + iniSo ) += weight * (1.0/deltat) * rockporosity * So * fOilDensity[0] * phiSoL2(iso,0);
        
        //        // dWaterDensity/dalphaP terms
        //        for (int jp = 0; jp < nphiPL2; jp++)
        //        {
        //            ek(iso  + iniSo, jp  + iniP) += weight * (1.0/deltat) * rockporosity * So * fOilDensity[1] * phiPL2(jp,0) * phiSoL2(iso,0);
        //        }
        
        // dSo/dalphaSo terms
        for (int jso = 0; jso < nphiSoL2; jso++)
        {
            ek(iso  + iniSo, jso  + iniSo) += weight * (1.0/deltat) * rockporosity * phiSoL2(jso,0) * fOilDensity[0] * phiSoL2(iso,0);
        }
    }
    
    
    this->Contribute(datavec,weight,ef);
    
}

// Residual contribution
void TPZAxiSymmetricDarcyFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    this->ContributeDarcy(datavec, weight, ef);
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeAlpha(datavec, weight, ef);
    }
    
    return;
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2        = datavec[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2        = datavec[Soblock].phi; // For L2  test functions So
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    TPZFMatrix<REAL> dphiSwL2   = datavec[Swblock].dphix; // Derivative For L2  test functions Sw
    TPZFMatrix<REAL> dphiSoL2   = datavec[Soblock].dphix; // Derivative For L2  test functions So
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSwL2    = phiSwL2.Rows();                                   // For L2   Sw
    int nphiSoL2    = phiSoL2.Rows();                                   // For L2   So
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSw   = nphiPL2       + iniP;
    int iniSo   = nphiSwL2      + iniSw;
    
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL Sw             = datavec[Swblock].sol[0][0];
    REAL So             = datavec[Soblock].sol[0][0];
    //    REAL Sg             = 1.0 - So - Sw;
    
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0]; // Piola divengence may works, needed set piola computation on the solution elchiv method!!!
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFMatrix<STATE> GradSwaxes = datavec[Swblock].dsol[0];
    TPZFMatrix<STATE> GradSoaxes = datavec[Soblock].dsol[0];
    
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    TPZFNMatrix<660,REAL> GradSw;
    TPZAxesTools<REAL>::Axes2XYZ(GradSwaxes, GradSw, datavec[Swblock].axes);
    TPZFNMatrix<660,REAL> GradSo;
    TPZAxesTools<REAL>::Axes2XYZ(GradSoaxes, GradSo, datavec[Soblock].axes);
    
    
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    REAL deltat = fSimulationData->GetDeltaT();
    this->UpdateStateVariables(u, P, Sw, So);
    this->PhaseFractionalFlows();
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    REAL rockporosity, drockporositydP;
    this->fReservoirdata->Porosity(fAveragePressure, rockporosity, drockporositydP);
    
    
    // Defining local variables
    TPZFMatrix<STATE> oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    
    oneoverlambda_Kinv_u(0,0) = (1.0/fTotalMobility[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    oneoverlambda_Kinv_u(1,0) = (1.0/fTotalMobility[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    Gravity = fSimulationData->GetGravity();
    
    gm(0,0) = (fFOil[0] * fOilDensity[0] + fFWater[0]* fWaterDensity[0]) * Gravity(0,0);
    gm(1,0) = (fFOil[0] * fOilDensity[0] + fFWater[0]* fWaterDensity[0]) * Gravity(1,0);
    
    REAL divu = 0.0;
    TPZFMatrix<STATE> iphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    
    /////////////////////////////////
    // Last State n
    if (fSimulationData->IsnStep()) {
        
        //  n state computations
        for (int ip = 0; ip < nphiPL2; ip++)
        {
            ef(ip + iniP) += 1.0 * weight * (1.0/deltat) * rockporosity * (Sw * fWaterDensity[0] +  So * fOilDensity[0] )*  phiPL2(ip,0);
        }
        
        for (int isw = 0; isw < nphiSwL2; isw++)
        {
            ef(isw  + iniSw ) += -1.0 * weight * (1.0/deltat) * rockporosity * Sw * fWaterDensity[0]  * phiSwL2(isw,0);
        }
        
        for (int iso = 0; iso < nphiSoL2; iso++)
        {
            ef(iso  + iniSo ) += -1.0 * weight * (1.0/deltat) * rockporosity * So * fOilDensity[0] * phiSoL2(iso,0);
        }
        
        return;
    }
    // Last State n
    /////////////////////////////////
    
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * ((oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))
                                   - P * DivergenceOnDeformed(iq,0)
                                   - (gm(0,0)*iphiuHdiv(0,0) + gm(1,0)*iphiuHdiv(1,0)) );
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    if(fForcingFunction)
    {
        fForcingFunction->Execute(datavec[Pblock].x,fvalue);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // uses this for constant jacobian elements
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        ef(ip + iniP) += weight * (- divu + fvalue[0] - (1.0/deltat) * rockporosity * (Sw * fWaterDensity[0] +  So * fOilDensity[0]) ) * phiPL2(ip,0);
        
    }
    
    
    // Transport equations
    
    for (int isw = 0; isw < nphiSwL2; isw++)
    {
        ef(isw  + iniSw ) += weight * (1.0/deltat) * rockporosity * Sw * fWaterDensity[0] * phiSwL2(isw,0);
        //        ef(isw  + iniSw ) += weight * (-1.0 * fFWater[2] *(u[0]*GradSw[0]) + (u[1]*GradSw[1])) * phiSwL2(isw,0);
        //        ef(isw  + iniSw ) += weight * (-1.0 * fFWater[0] *(u[0]*dphiSwL2(0,isw)) + (u[1]*dphiSwL2(1,isw)));
    }
    
    for (int iso = 0; iso < nphiSoL2; iso++)
    {
        ef(iso  + iniSo ) += weight * (1.0/deltat) * rockporosity * So * fOilDensity[0] * phiSoL2(iso,0);
        //        ef(iso  + iniSo ) += weight * (-1.0 * fFOil[2] *(u[0]*GradSo[0]) + (u[1]*GradSo[1])) * phiSoL2(iso,0);
        //        ef(iso  + iniSo ) += weight * (-1.0 * fFOil[0] *(u[0]*dphiSoL2(0,iso)) + (u[1]*dphiSoL2(1,iso)));
    }
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeInterfaceAlpha(data, datavecleft, datavecright, weight, ek, ef);
    }
    
    return;
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2L        = datavecleft[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2L        = datavecleft[Soblock].phi; // For L2  test functions So
    TPZFMatrix<REAL> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2R        = datavecright[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2R        = datavecright[Soblock].phi; // For L2  test functions So
    TPZFMatrix<REAL> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSwL2L    = phiSwL2L.Rows();                                   // For L2   Sw
    int nphiSoL2L    = phiSoL2L.Rows();                                   // For L2   So
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSwL   = nphiPL2L       + iniPL;
    int iniSoL   = nphiSwL2L      + iniSwL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSwL2R    = phiSwL2R.Rows();                                   // For L2   Sw
    int nphiSoL2R    = phiSoL2R.Rows();                                   // For L2   So
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSwR   = nphiPL2R       + iniPR;
    int iniSoR   = nphiSwL2R      + iniSwR;
    
    int iblock = iniSoL + nphiSoL2L;
    int jblock = iniSoL + nphiSoL2L;
    
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SwL             = datavecleft[Swblock].sol[0][0];
    REAL SoL             = datavecleft[Soblock].sol[0][0];
    //    REAL SgL             = 1.0 - SoL - SwL;
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL PR              = datavecright[Pblock].sol[0][0];
    REAL SwR             = datavecright[Swblock].sol[0][0];
    REAL SoR             = datavecright[Soblock].sol[0][0];
    //    REAL SgR             = 1.0 - SoR - SwR;
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Upwind scheme
    // set for upwind derivatives
    
    TPZFMatrix<REAL> phiuH1;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2; // For L2  test functions So
    TPZFMatrix<STATE> dphiuH1; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2; // Derivative For L2  test functions
    
    int nphiuHdiv   = 0;       // For Hdiv u
    int nphiPL2     = 0;                                    // For L2   P
    int nphiSwL2    = 0;                                   // For L2   Sw
    int nphiSoL2    = 0;                                   // For L2   So
    int iniu    = 0;
    int iniP    = 0;
    int iniSw   = 0;
    int iniSo   = 0;
    
    int iblockt = 0;
    int jblockt = 0;
    
    if (uLn >= 0.0 )
    {
        this->UpdateStateVariables(uL, PL, SwL, SoL);
        this->PhaseFractionalFlows();
        
        phiuH1 = phiuH1L;  // For H1  test functions u
        phiPL2 = phiPL2L;  // For L2  test functions P
        phiSwL2 = phiSwL2L; // For L2  test functions Sw
        phiSoL2 = phiSoL2L; // For L2  test functions So
        dphiuH1 = dphiuH1L; // Derivative For H1  test functions
        dphiPL2 = dphiPL2L; // Derivative For L2  test functions
        
        nphiuHdiv   = nphiuHdivL;       // For Hdiv u
        nphiPL2     = nphiPL2L;                                    // For L2   P
        nphiSwL2    = nphiSwL2L;                                   // For L2   Sw
        nphiSoL2    = nphiSoL2L;                                   // For L2   So
        iniu    = iniuL;
        iniP    = iniPL;
        iniSw   = iniSwL;
        iniSo   = iniSoL;
        
        iblockt = 0;
        jblockt = 0;
        
    }
    else
    {
        this->UpdateStateVariables(uR, PR, SwR, SoR);
        this->PhaseFractionalFlows();
        
        phiuH1 = phiuH1R;  // For H1  test functions u
        phiPL2 = phiPL2R;  // For L2  test functions P
        phiSwL2 = phiSwL2R; // For L2  test functions Sw
        phiSoL2 = phiSoL2R; // For L2  test functions So
        dphiuH1 = dphiuH1R; // Derivative For H1  test functions
        dphiPL2 = dphiPL2R; // Derivative For L2  test functions
        
        nphiuHdiv   = nphiuHdivR;       // For Hdiv u
        nphiPL2     = nphiPL2R;                                    // For L2   P
        nphiSwL2    = nphiSwL2R;                                   // For L2   Sw
        nphiSoL2    = nphiSoL2R;                                   // For L2   So
        iniu    = iniuR;
        iniP    = iniPR;
        iniSw   = iniSwR;
        iniSo   = iniSoR;
        
        iblockt = iblock;
        jblockt = jblock;
        
    }
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    
    TPZFMatrix<REAL> jphiuHdiv(3,1,0.0);
    
    
    for (int isw = 0; isw < nphiSwL2L; isw++)
    {
        
        // dP/dalphaP
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(isw + iniSwL, jp + iniP + jblockt) += 1.0 * weight * fFWater[1] * phiPL2(jp,0) * phiSwL2L(isw,0) * uLn;
        }
        
        // dSw/dalphaSw
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(isw + iniSwL, jsw + iniSw + jblockt) += 1.0 * weight * fFWater[2] * phiSwL2(jsw,0) * phiSwL2L(isw,0) * uLn;
        }
        
    }
    
    for (int isw = 0; isw < nphiSwL2R; isw++)
    {
        
        // dP/dalphaP
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iblock + isw + iniSwR, jp + iniP + jblockt) += -1.0 * weight * fFWater[1] * phiPL2(jp,0) * phiSwL2R(isw,0) * uLn;
        }
        
        // dSw/dalphaSw
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(iblock + isw + iniSwR, jsw + iniSw + jblockt) += -1.0 * weight * fFWater[2] * phiSwL2(jsw,0) * phiSwL2R(isw,0) * uLn;
        }
        
    }
    
    
    
    
    for (int iso = 0; iso < nphiSoL2L; iso++)
    {
        //        // duL/dalphauL
        //        for (int jq = 0; jq < nphiuHdivL; jq++)
        //        {
        //            int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
        //            int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
        //            jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
        //            jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
        //            jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
        //
        //            REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
        //
        //            ek(iso + iniSoL, jq + iniuL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * jphiuHdivn;
        //        }
        //
        // dP/dalphaP
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iso + iniSoL, jp + iniP + jblockt) += 1.0 * weight * fFOil[1] * phiPL2(jp,0) * phiSoL2L(iso,0) * uLn;
        }
        
        // dSw/dalphaSw
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(iso + iniSoL, jsw + iniSw + jblockt) += 1.0 * weight * fFOil[2] * phiSwL2(jsw,0) * phiSoL2L(iso,0) * uLn;
        }
        
        // dSo/dalphaSo
        for (int jso = 0; jso < nphiSoL2; jso++)
        {
            ek(iso + iniSoL, jso + iniSo + jblockt) += 1.0 * weight * fFOil[3] * phiSoL2(jso,0) * phiSoL2L(iso,0) * uLn;
        }
    }
    
    for (int iso = 0; iso < nphiSoL2R; iso++)
    {
        //        // duL/dalphauL
        //        for (int jq = 0; jq < nphiuHdivL; jq++)
        //        {
        //            int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
        //            int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
        //            jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
        //            jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
        //            jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
        //
        //            REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
        //
        //            ek(iblock + iso + iniSoR, jq + iniuL) += -1.0 * weight * fFOil[0] * phiSoL2R(iso,0) * jphiuHdivn;
        //        }
        //
        // dP/dalphaP
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iblock + iso + iniSoR, jp + iniP + jblockt) += -1.0 * weight * fFOil[1] * phiPL2(jp,0) * phiSoL2R(iso,0) * uLn;
        }
        
        // dSw/dalphaSw
        for (int jsw = 0; jsw < nphiSwL2; jsw++)
        {
            ek(iblock + iso + iniSoR, jsw + iniSw + jblockt) += -1.0 * weight * fFOil[2] * phiSwL2(jsw,0) * phiSoL2R(iso,0) * uLn;
        }
        
        // dSo/dalphaSo
        for (int jso = 0; jso < nphiSoL2; jso++)
        {
            ek(iblock + iso + iniSoR, jso + iniSo + jblockt) += -1.0 * weight * fFOil[3] * phiSoL2(jso,0) * phiSoL2R(iso,0) * uLn;
        }
        
    }
    
    //    //Gravity
    //
    //    TPZFMatrix<STATE> Gravity(2,1);
    //    TPZFMatrix<STATE> KGravityL(2,1);
    //    TPZFMatrix<STATE> KGravityR(2,1);
    //
    //    Gravity = fSimulationData->GetGravity();
    //
    //    REAL epsilon = fepsilon;
    //
    //    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    //    REAL ndotG = Gravity(0,0)*n[0] + Gravity(1,0)*n[1];
    //
    //    // Computing Gravitational segregational function on the left side
    //
    //    this->UpdateStateVariables(uL, PL, SwL, SoL);
    //    this->PhaseFractionalFlows();
    //
    //    KGravityL(0,0) = K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0);
    //    KGravityL(1,0) = K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0);
    //    REAL lambdaDensitydiffL = fTotalMobility[0] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffLdP  = fTotalMobility[0] * (fWaterDensity[1] - fOilDensity[1]) + fTotalMobility[1] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffLdSw = fTotalMobility[0] * (fWaterDensity[2] - fOilDensity[2]) + fTotalMobility[2] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffLdSo = fTotalMobility[0] * (fWaterDensity[3] - fOilDensity[3]) + fTotalMobility[3] * (fWaterDensity[0] - fOilDensity[0]);
    //
    //    REAL fstrL = 0.0;
    //    REAL dfstrLdP  = 0.0;
    //    REAL dfstrLdSw = 0.0;
    //    REAL dfstrLdSo = 0.0;
    //
    //    if (ndotG <= 0.0) {
    //        // Expelling oil
    //        if (SoL < epsilon) {
    //            this->UpdateStateVariables(uL, PL, SwL, SoL);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //            dfstrLdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrLdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrLdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uL, PL, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //            dfstrLdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrLdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrLdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //    }
    //    else
    //    {
    //        // Expelling water
    //        if (SwL < 1.0 - epsilon) {
    //            this->UpdateStateVariables(uL, PL, SwL, SoL);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //            dfstrLdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrLdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrLdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uL, PL, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //            dfstrLdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrLdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrLdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //
    //    }
    //
    //    TPZFMatrix<STATE> qgL(2,1);
    //    TPZFMatrix<STATE> dqgLdP(2,1);
    //    TPZFMatrix<STATE> dqgLdSw(2,1);
    //    TPZFMatrix<STATE> dqgLdSo(2,1);
    //
    //    qgL(0,0) = fstrL * lambdaDensitydiffL * KGravityL(0,0);
    //    qgL(1,0) = fstrL * lambdaDensitydiffL * KGravityL(1,0);
    //
    //    dqgLdP(0,0) = (dfstrLdP * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdP) * KGravityL(0,0);
    //    dqgLdP(1,0) = (dfstrLdP * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdP) * KGravityL(1,0);
    //
    //    dqgLdSw(0,0) = (dfstrLdSw * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdSw) * KGravityL(0,0);
    //    dqgLdSw(1,0) = (dfstrLdSw * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdSw) * KGravityL(1,0);
    //
    //    dqgLdSo(0,0) = (dfstrLdSo * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdSo) * KGravityL(0,0);
    //    dqgLdSo(1,0) = (dfstrLdSo * lambdaDensitydiffL + fstrL * dlambdaDensitydiffLdSo) * KGravityL(1,0);
    //
    //    REAL qgLdotn = qgL(0,0)*n[0] + qgL(1,0)*n[1];
    //    REAL dqgLdotndP  = dqgLdP(0,0)*n[0] + dqgLdP(1,0)*n[1];
    //    REAL dqgLdotndSw = dqgLdSw(0,0)*n[0] + dqgLdSw(1,0)*n[1];
    //    REAL dqgLdotndSo = dqgLdSo(0,0)*n[0] + dqgLdSo(1,0)*n[1];
    //
    //
    //
    //    // Computing Gravitational segregational function on the right side
    //
    //    this->UpdateStateVariables(uR, PR, SwR, SoR);
    //    this->PhaseFractionalFlows();
    //
    //    KGravityR(0,0) = K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0);
    //    KGravityR(1,0) = K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0);
    //    REAL lambdaDensitydiffR = fTotalMobility[0] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffRdP  = fTotalMobility[0] * (fWaterDensity[1] - fOilDensity[1]) + fTotalMobility[1] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffRdSw = fTotalMobility[0] * (fWaterDensity[2] - fOilDensity[2]) + fTotalMobility[2] * (fWaterDensity[0] - fOilDensity[0]);
    //    REAL dlambdaDensitydiffRdSo = fTotalMobility[0] * (fWaterDensity[3] - fOilDensity[3]) + fTotalMobility[3] * (fWaterDensity[0] - fOilDensity[0]);
    //
    //    REAL fstrR = 0.0;
    //    REAL dfstrRdP  = 0.0;
    //    REAL dfstrRdSw = 0.0;
    //    REAL dfstrRdSo = 0.0;
    //    // std::cout<<" data.normal: "<< data.normal<<std::endl;
    //    if (ndotG >= 0.0) {
    //        // Expelling oil
    //        if (SoR < epsilon) {
    //            this->UpdateStateVariables(uR, PR, SwR, SoR);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //            dfstrRdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrRdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrRdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uR, PR, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //            dfstrRdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrRdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrRdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //    }
    //    else
    //    {
    //        // Expelling water
    //        if (SwR < 1.0 - epsilon) {
    //            this->UpdateStateVariables(uR, PR, SwR, SoR);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //            dfstrRdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrRdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrRdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uR, PR, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //            dfstrRdP  = fFOil[0] * fFWater[1] + fFOil[1] * fFWater[0];
    //            dfstrRdSw = fFOil[0] * fFWater[2] + fFOil[2] * fFWater[0];
    //            dfstrRdSo = fFOil[0] * fFWater[3] + fFOil[3] * fFWater[0];
    //        }
    //
    //    }
    //
    //
    //    TPZFMatrix<STATE> qgR(2,1);
    //    TPZFMatrix<STATE> dqgRdP(2,1);
    //    TPZFMatrix<STATE> dqgRdSw(2,1);
    //    TPZFMatrix<STATE> dqgRdSo(2,1);
    //
    //    qgR(0,0) = fstrR * lambdaDensitydiffR * KGravityR(0,0);
    //    qgR(1,0) = fstrR * lambdaDensitydiffR * KGravityR(1,0);
    //
    //    dqgRdP(0,0) = (dfstrRdP * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdP) * KGravityR(0,0);
    //    dqgRdP(1,0) = (dfstrRdP * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdP) * KGravityR(1,0);
    //
    //    dqgRdSw(0,0) = (dfstrRdSw * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdSw) * KGravityR(0,0);
    //    dqgRdSw(1,0) = (dfstrRdSw * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdSw) * KGravityR(1,0);
    //
    //    dqgRdSo(0,0) = (dfstrRdSo * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdSo) * KGravityR(0,0);
    //    dqgRdSo(1,0) = (dfstrRdSo * lambdaDensitydiffR + fstrR * dlambdaDensitydiffRdSo) * KGravityR(1,0);
    //
    //    REAL qgRdotn = qgR(0,0)*n[0] + qgR(1,0)*n[1];
    //    REAL dqgRdotndP  = dqgRdP(0,0)*n[0] + dqgRdP(1,0)*n[1];
    //    REAL dqgRdotndSw = dqgRdSw(0,0)*n[0] + dqgRdSw(1,0)*n[1];
    //    REAL dqgRdotndSo = dqgRdSo(0,0)*n[0] + dqgRdSo(1,0)*n[1];
    //
    //    // Computing the minimum flux at the edge
    //
    //    REAL gqdotn = 0.0;
    //    REAL dgqdotndP  = 0.0;
    //    REAL dgqdotndSw = 0.0;
    //    REAL dgqdotndSo = 0.0;
    //
    //    int nphiu = 0;
    //    int nphiP = 0;
    //    int nphiSw = 0;
    //    int nphiSo = 0;
    //    TPZFMatrix<REAL> phiP;
    //    TPZFMatrix<REAL> phiSw;
    //    TPZFMatrix<REAL> phiSo;
    //    int block = 0;
    //
    //    if (fabs(qgRdotn) >= fabs(qgLdotn)) {
    //        gqdotn = qgLdotn;
    //        dgqdotndP  = dqgLdotndP;
    //        dgqdotndSw = dqgLdotndSw;
    //        dgqdotndSo = dqgLdotndSo;
    //        nphiu = nphiuHdivL;
    //        nphiP = nphiPL2L;
    //        nphiSw = nphiSwL2L;
    //        nphiSo = nphiSoL2L;
    //        phiP  = phiPL2L;
    //        phiSw = phiSwL2L;
    //        phiSo = phiSoL2L;
    //        block = 0;
    //    }
    //    else
    //    {
    //        gqdotn = qgRdotn;
    //        dgqdotndP  = dqgRdotndP;
    //        dgqdotndSw = dqgRdotndSw;
    //        dgqdotndSo = dqgRdotndSo;
    //        nphiu = nphiuHdivR;
    //        nphiP = nphiPL2R;
    //        nphiSw = nphiSwL2R;
    //        nphiSo = nphiSoL2R;
    //        phiP  = phiPL2R;
    //        phiSw = phiSwL2R;
    //        phiSo = phiSoL2R;
    //        block = iblock;
    //    }
    //
    //
    //    for (int isw = 0; isw < nphiSwL2L; isw++)
    //    {
    //        for (int jp = 0; jp < nphiP; jp++)
    //        {
    //            ek(isw + iniSwL, jp + block + nphiu) += 1.0 * weight * dgqdotndP * phiP(jp,0) * phiSwL2L(isw,0);
    //        }
    //
    //        for (int jsw = 0; jsw < nphiSw; jsw++)
    //        {
    //            ek(isw + iniSwL, jsw + block + nphiP + nphiu) += 1.0 * weight * dgqdotndSw * phiSw(jsw,0) * phiSwL2L(isw,0);
    //        }
    //
    //        for (int jso = 0; jso < nphiSo; jso++)
    //        {
    //            ek(isw + iniSwL, jso + block + nphiSw + nphiP + nphiu) += 1.0 * weight * dgqdotndSo * phiSo(jso,0) * phiSwL2L(isw,0);
    //        }
    //
    //    }
    //    for (int isw = 0; isw < nphiSwL2R; isw++)
    //    {
    //        for (int jp = 0; jp < nphiP; jp++)
    //        {
    //            ek(iblock + isw + iniSwR, jp + block + nphiu) += -1.0 * weight * dgqdotndP * phiP(jp,0) * phiSwL2R(isw,0);
    //        }
    //
    //        for (int jsw = 0; jsw < nphiSw; jsw++)
    //        {
    //            ek(iblock + isw + iniSwR, jsw + block + nphiP + nphiu) += -1.0 * weight * dgqdotndSw * phiSw(jsw,0) * phiSwL2R(isw,0);
    //        }
    //
    //        for (int jso = 0; jso < nphiSo; jso++)
    //        {
    //            ek(iblock + isw + iniSwR, jso + block + nphiSw + nphiP + nphiu) += -1.0 * weight * dgqdotndSo * phiSo(jso,0) * phiSwL2R(isw,0);
    //        }
    //    }
    //
    //
    //    for (int iso = 0; iso < nphiSoL2L; iso++)
    //    {
    //        for (int jp = 0; jp < nphiP; jp++)
    //        {
    //            ek(iso + iniSoL, jp + block + nphiu) += -1.0 * weight * dgqdotndP * phiP(jp,0) * phiSoL2L(iso,0);
    //        }
    //
    //        for (int jsw = 0; jsw < nphiSw; jsw++)
    //        {
    //            ek(iso + iniSoL, jsw + block + nphiP + nphiu) += -1.0 * weight * dgqdotndSw * phiSw(jsw,0) * phiSoL2L(iso,0);
    //        }
    //
    //        for (int jso = 0; jso < nphiSo; jso++)
    //        {
    //            ek(iso + iniSoL, jso + block + nphiSw + nphiP + nphiu) += -1.0 * weight * dgqdotndSo * phiSo(jso,0) * phiSoL2L(iso,0);
    //        }
    //    }
    //    for (int iso = 0; iso < nphiSoL2R; iso++)
    //    {
    //        for (int jp = 0; jp < nphiP; jp++)
    //        {
    //            ek(iblock + iso + iniSoR, jp + block + nphiu) += 1.0 * weight * dgqdotndP * phiP(jp,0) * phiSoL2R(iso,0);
    //        }
    //
    //        for (int jsw = 0; jsw < nphiSw; jsw++)
    //        {
    //            ek(iblock + iso + iniSoR, jsw + block + nphiP + nphiu) += 1.0 * weight * dgqdotndSw * phiSw(jsw,0) * phiSoL2R(iso,0);
    //        }
    //
    //        for (int jso = 0; jso < nphiSo; jso++)
    //        {
    //            ek(iblock + iso + iniSoR, jso + block + nphiSw + nphiP + nphiu) += 1.0 * weight * dgqdotndSo * phiSo(jso,0) * phiSoL2R(iso,0);
    //        }
    //    }
    
    
    this->ContributeInterface(data, datavecleft, datavecright, weight, ef);
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeInterfaceAlpha(data, datavecleft, datavecright, weight, ef);
    }
    
    return;
    
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2L        = datavecleft[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2L        = datavecleft[Soblock].phi; // For L2  test functions So
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2R        = datavecright[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2R        = datavecright[Soblock].phi; // For L2  test functions So
    TPZFMatrix<STATE> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSwL2L    = phiSwL2L.Rows();                                   // For L2   Sw
    int nphiSoL2L    = phiSoL2L.Rows();                                   // For L2   So
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSwL   = nphiPL2L       + iniPL;
    int iniSoL   = nphiSwL2L      + iniSwL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSwL2R    = phiSwL2R.Rows();                                   // For L2   Sw
    int nphiSoL2R    = phiSoL2R.Rows();                                   // For L2   So
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSwR   = nphiPL2R       + iniPR;
    int iniSoR   = nphiSwL2R      + iniSwR;
    
    int iblock = iniSoL + nphiSoL2L;
    //    int jblock = iniSoL + nphiSoL2L;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SwL             = datavecleft[Swblock].sol[0][0];
    REAL SoL             = datavecleft[Soblock].sol[0][0];
    //    REAL SgL             = 1.0 - SoL - SwL;
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL PR              = datavecright[Pblock].sol[0][0];
    REAL SwR             = datavecright[Swblock].sol[0][0];
    REAL SoR             = datavecright[Soblock].sol[0][0];
    //    REAL SgR             = 1.0 - SoR - SwR;
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    if (uLn >= 0.0 ) {
        this->UpdateStateVariables(uL, PL, SwL, SoL);
        this->PhaseFractionalFlows();
    }
    else
    {
        this->UpdateStateVariables(uR, PR, SwR, SoR);
        this->PhaseFractionalFlows();
    }
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    
    // contribuicao para elemento esquerdo
    for (int isw = 0; isw < nphiSwL2L; isw++)
    {
        ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
    }
    // contribuicao para elemento direito
    for (int isw = 0; isw < nphiSwL2R; isw++)
    {
        ef(iblock + isw + iniSwR) += -1.0 * weight * fFWater[0] * phiSwL2R(isw,0) * uLn;
    }
    
    for (int iso = 0; iso < nphiSoL2L; iso++)
    {
        ef(iso + iniSoL) += 1.0 * weight * (1-fFWater[0]) * phiSoL2L(iso,0) * uLn;
    }
    for (int iso = 0; iso < nphiSoL2R; iso++)
    {
        ef(iblock + iso + iniSoR) += -1.0 * weight * (1-fFWater[0]) * phiSoL2R(iso,0) * uLn;
    }
    
    
    //    //Gravity
    //
    //    TPZFMatrix<STATE> Gravity(2,1);
    //    TPZFMatrix<STATE> KGravityL(2,1);
    //    TPZFMatrix<STATE> KGravityR(2,1);
    //
    //    Gravity=fSimulationData->GetGravity();
    //
    //    TPZFMatrix<STATE> qgL(2,1);
    //    TPZFMatrix<STATE> qgR(2,1);
    //
    //    REAL epsilon = fepsilon;
    //
    //    TPZFMatrix<STATE> K = fReservoirdata->Kabsolute();
    //    REAL ndotG = Gravity(0,0)*n[0] + Gravity(1,0)*n[1];
    //
    //
    //    // Computing Gravitational segregational function on the left side
    //
    //    this->UpdateStateVariables(uL, PL, SwL, SoL);
    //    this->PhaseFractionalFlows();
    //
    //    //std::cout<<"gelElId= " << datavecright[Pblock].gelElId << std::endl;
    //
    //    KGravityL(0,0) = K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0);
    //    KGravityL(1,0) = K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0);
    //    REAL lambdaDensitydiffL = fTotalMobility[0] * (fWaterDensity[0] - fOilDensity[0]);
    //
    //    REAL fstrL = 0.0;
    //
    //    if (ndotG <= 0.0) {
    //
    //        // Expelling oil
    //        if (SoL < epsilon) {
    //            this->UpdateStateVariables(uL, PL, SwL, SoL);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uL, PL, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //        }
    //    }
    //    else
    //    {
    //        // Expelling water
    //        if (SwL < 1.0 - epsilon) {
    //            this->UpdateStateVariables(uL, PL, SwL, SoL);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uL, PL, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrL = fFOil[0] * fFWater[0];
    //        }
    //
    //    }
    
    //    qgL(0,0) = fstrL * lambdaDensitydiffL * KGravityL(0,0);
    //    qgL(1,0) = fstrL * lambdaDensitydiffL * KGravityL(1,0);
    //    REAL qgLdotn = qgL(0,0)*n[0] + qgL(1,0)*n[1];
    //
    //    // Computing Gravitational segregational function on the right side
    //
    //    this->UpdateStateVariables(uR, PR, SwR, SoR);
    //    this->PhaseFractionalFlows();
    //
    //    KGravityR(0,0) = K(0,0)*Gravity(0,0) + K(0,1)*Gravity(1,0);
    //    KGravityR(1,0) = K(1,0)*Gravity(0,0) + K(1,1)*Gravity(1,0);
    //    REAL lambdaDensitydiffR = fTotalMobility[0] * (fWaterDensity[0] - fOilDensity[0]);
    //
    //    REAL fstrR = 0.0;
    //
    //    if (ndotG >= 0.0) {
    //
    //        // Expelling oil
    //        if (SoR < epsilon) {
    //            this->UpdateStateVariables(uR, PR, SwR, SoR);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uR, PR, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //        }
    //    }
    //    else
    //    {
    //        // Expelling water
    //        if (SwR < 1.0 - epsilon) {
    //            this->UpdateStateVariables(uR, PR, SwR, SoR);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //        }
    //        else
    //        {
    //            this->UpdateStateVariables(uR, PR, 1.0 - epsilon, epsilon);
    //            this->PhaseFractionalFlows();
    //            fstrR = fFOil[0] * fFWater[0];
    //        }
    //
    //    }
    //
    //    qgR(0,0) = fstrR * lambdaDensitydiffR * KGravityR(0,0);
    //    qgR(1,0) = fstrR * lambdaDensitydiffR * KGravityR(1,0);
    //    REAL qgRdotn = qgR(0,0)*n[0] + qgR(1,0)*n[1];
    //
    //    // Computing the minimum flux at the edge
    //
    //    REAL gqdotn = 0.0;
    //
    //    if (fabs(qgRdotn) >= fabs(qgLdotn)) {
    //        gqdotn = qgLdotn;
    //    }
    //    else
    //    {
    //        gqdotn = qgRdotn;
    //    }
    //
    //
    //
    //    for (int isw = 0; isw < nphiSwL2L; isw++)
    //    {
    //        ef(isw + iniSwL) += 1.0 * weight * gqdotn * phiSwL2L(isw,0);
    //    }
    //    for (int isw = 0; isw < nphiSwL2R; isw++)
    //    {
    //        ef(iblock + isw + iniSwR) += -1.0 * weight * gqdotn * phiSwL2R(isw,0);
    //    }
    //
    //    for (int iso = 0; iso < nphiSoL2L; iso++)
    //    {
    //        ef(iso + iniSoL) += -1.0 * weight * gqdotn * phiSoL2L(iso,0);
    //    }
    //    for (int iso = 0; iso < nphiSoL2R; iso++)
    //    {
    //        ef(iblock + iso + iniSoR) += 1.0 * weight * gqdotn * phiSoL2R(iso,0);
    //    }
    
}

void TPZAxiSymmetricDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeBCInterfaceAlpha(data, datavecleft, weight, ek, ef, bc);
    }
    
    return;
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2L        = datavecleft[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2L        = datavecleft[Soblock].phi; // For L2  test functions So
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSwL2L    = phiSwL2L.Rows();                                   // For L2   Sw
    int nphiSoL2L    = phiSoL2L.Rows();                                   // For L2   So
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSwL   = nphiPL2L       + iniPL;
    int iniSoL   = nphiSwL2L      + iniSwL;
    
    
    //    int iblock = iniSoL + nphiSoL2L;
    //    int jblock = iniSoL + nphiSoL2L;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SwL             = datavecleft[Swblock].sol[0][0];
    REAL SoL             = datavecleft[Soblock].sol[0][0];
    //    REAL SgL             = 1.0 - SoL - SwL;
    
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    REAL  Value, Sw, So;
    Value = bc.Val2()(0,0);         //  un ou P
    Sw = bc.Val2()(1,0);         //  Sw in
    So = bc.Val2()(2,0);         //  So in
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    TPZFMatrix<REAL> jphiuHdiv(3,1,0.0);
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            if (uLn >= 0.0 ) {
                // Outflow boundary condition
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            else
            {
                // Inflow boundary condition
                if (uLn < 0.0 && fabs(uLn) < 1.0e-24) { std::cout << "Warning: Inflow boundary condition detected in outflow boundary condition: uLn = " << uLn << "\n";}
                this->UpdateStateVariables(uL, Value, Sw, So);
                this->PhaseFractionalFlows();
            }
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                //ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                    ek(isw + iniSwL, jq + iniuL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * jphiuHdivn;
                //                }
                
                //                // dSw/dalphaSw
                //                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                //                {
                //                    ek(isw + iniSwL, jsw + iniSwL) += 1.0 * weight * fFWater[2] * phiSwL2L(jsw,0) * phiSwL2L(isw,0) * uLn;
                //                }
                
                //                // dSo/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(isw + iniSwL, jso + iniSoL) += 1.0 * weight * fFWater[3] * phiSoL2L(jso,0) * phiSwL2L(isw,0) * uLn;
                //                }
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                //ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                    ek(iso + iniSoL, jq + iniuL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * jphiuHdivn;
                //                }
                
                //                // dSw/dalphaSw
                //                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                //                {
                //                    ek(iso + iniSoL, jsw + iniSwL) += 1.0 * weight * fFOil[2] * phiSwL2L(jsw,0) * phiSoL2L(iso,0) * uLn;
                //                }
                
                //                // dSo/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(iso + iniSoL, jso + iniSoL) += 1.0 * weight * fFOil[3] * phiSoL2L(jso,0) * phiSoL2L(iso,0) * uLn;
                //                }
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            this->UpdateStateVariables(uL, PL, Sw, So);
            this->PhaseFractionalFlows();
            //            }
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                //ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                  int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                  int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                  jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                  jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                  jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                  REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                  ek(isw + iniSwL, jq + iniuL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * jphiuHdivn;
                //                }
                
                //                // dP/dalphaP
                //                for (int jp = 0; jp < nphiPL2L; jp++)
                //                {
                //                    ek(isw + iniSwL, jp + iniPL) += 1.0 * weight * fFWater[1] * phiPL2L(jp,0) * phiSwL2L(isw,0) * Value;
                //                }
                
                //                // dSw/dalphaSw
                //                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                //                {
                //                    ek(isw + iniSwL, jsw + iniSwL) += 1.0 * weight * fFWater[2] * phiSwL2L(jsw,0) * phiSwL2L(isw,0) * Value;
                //                }
                
                //                // dSo/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(isw + iniSwL, jso + iniSoL) += 1.0 * weight * fFWater[3] * phiSoL2L(jso,0) * phiSwL2L(isw,0) * Value;
                //                }
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                //ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                    ek(iso + iniSoL, jq + iniuL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * jphiuHdivn;
                //                }
                
                //                // dP/dalphaP
                //                for (int jp = 0; jp < nphiPL2L; jp++)
                //                {
                //                    ek(iso + iniSoL, jp + iniPL) += 1.0 * weight * fFOil[1] * phiPL2L(jp,0) * phiSoL2L(iso,0) * Value;
                //                }
                
                //                // dSo/dalphaSw
                //                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                //                {
                //                    ek(iso + iniSoL, jsw + iniSwL) += 1.0 * weight * fFOil[2] * phiSwL2L(jsw,0) * phiSoL2L(iso,0) * Value;
                //                }
                
                //                // dSo/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(iso + iniSoL, jso + iniSoL) += 1.0 * weight * fFOil[3] * phiSoL2L(jso,0) * phiSoL2L(iso,0) * Value;
                //                }
                
            }
            
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            
            if (uLn >= 0.0 ) {
                // Outflow boundary condition
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            else
            {
                // Inflow boundary condition
                if (uLn < 0.0 && fabs(uLn) < 1.0e-24) { std::cout << "Warning: Inflow boundary condition detected in outflow boundary condition: uLn = " << uLn << "\n";}
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                //ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                    ek(isw + iniSwL, jq + iniuL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * jphiuHdivn;
                //                }
                
                // dSw/dalphaSw
                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                {
                    ek(isw + iniSwL, jsw + iniSwL) += 1.0 * weight * fFWater[2] * phiSwL2L(jsw,0) * phiSwL2L(isw,0) * uLn;
                }
                
                //                // dSw/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(isw + iniSwL, jso + iniSoL) += 1.0 * weight * fFWater[3] * phiSoL2L(jso,0) * phiSwL2L(isw,0) * uLn;
                //                }
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                //ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
                //                // duL/dalphauL
                //                for (int jq = 0; jq < nphiuHdivL; jq++)
                //                {
                //                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                //                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                //                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                //                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                //                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                //
                //                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                //
                //                    ek(iso + iniSoL, jq + iniuL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * jphiuHdivn;
                //                }
                
                // dSo/dalphaSw
                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                {
                    ek(iso + iniSoL, jsw + iniSwL) += 1.0 * weight * fFOil[2] * phiSwL2L(jsw,0) * phiSoL2L(iso,0) * uLn;
                }
                
                //                // dSo/dalphaSo
                //                for (int jso = 0; jso < nphiSoL2L; jso++)
                //                {
                //                    ek(iso + iniSoL, jso + iniSoL) += 1.0 * weight * fFOil[3] * phiSoL2L(jso,0) * phiSoL2L(iso,0) * uLn;
                //                }
                
            }
            
            
        }
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            
            
            //            if (uLn >= 0.0 ) {
            // Outflow boundary condition
            this->UpdateStateVariables(uL, PL, SwL, SoL);
            this->PhaseFractionalFlows();
            //            }
            //            else
            //            {
            //                // Inflow boundary condition
            //                if (uLn < 0.0 && fabs(uLn) < 1.0e-24) { std::cout << "Boundary condition error: inflow detected in outflow boundary condition: uLn = " << uLn << "\n";}
            //                this->UpdateStateVariables(uL, PL, Sw, So);
            //                this->PhaseFractionalFlows();
            //            }
            // return;
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                //ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
                
                // duL/dalphauL
                for (int jq = 0; jq < nphiuHdivL; jq++)
                {
                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                    
                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                    
                    ek(isw + iniSwL, jq + iniuL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * jphiuHdivn;
                }
                
                
                // dP/dalphaP
                for (int jp = 0; jp < nphiPL2L; jp++)
                {
                    ek(isw + iniSwL, jp + iniPL) += 1.0 * weight * fFWater[1] * phiPL2L(jp,0) * phiSwL2L(isw,0) * Value;
                }
                
                // dSw/dalphaSw
                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                {
                    ek(isw + iniSwL, jsw + iniSwL) += 1.0 * weight * fFWater[2] * phiSwL2L(jsw,0) * phiSwL2L(isw,0) * Value;
                }
                
                // dSo/dalphaSo
                for (int jso = 0; jso < nphiSoL2L; jso++)
                {
                    ek(isw + iniSwL, jso + iniSoL) += 1.0 * weight * fFWater[3] * phiSoL2L(jso,0) * phiSwL2L(isw,0) * Value;
                }
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                //ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
                // duL/dalphauL
                for (int jq = 0; jq < nphiuHdivL; jq++)
                {
                    int jvectorindex    = datavecleft[ublock].fVecShapeIndex[jq].first;
                    int jshapeindex     = datavecleft[ublock].fVecShapeIndex[jq].second;
                    jphiuHdiv(0,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(0,jvectorindex);
                    jphiuHdiv(1,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(1,jvectorindex);
                    jphiuHdiv(2,0) = phiuH1L(jshapeindex,0)*datavecleft[ublock].fNormalVec(2,jvectorindex);
                    
                    REAL jphiuHdivn = jphiuHdiv(0,0)*n[0] + jphiuHdiv(1,0)*n[1] + jphiuHdiv(2,0)*n[2];
                    
                    ek(iso + iniSoL, jq + iniuL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * jphiuHdivn;
                }
                
                // dP/dalphaP
                for (int jp = 0; jp < nphiPL2L; jp++)
                {
                    ek(iso + iniSoL, jp + iniPL) += 1.0 * weight * fFOil[1] * phiPL2L(jp,0) * phiSoL2L(iso,0) * Value;
                }
                
                // dSw/dalphaSw
                for (int jsw = 0; jsw < nphiSwL2L; jsw++)
                {
                    ek(iso + iniSoL, jsw + iniSwL) += 1.0 * weight * fFOil[2] * phiSwL2L(jsw,0) * phiSoL2L(iso,0) * Value;
                }
                
                // dSo/dalphaSo
                for (int jso = 0; jso < nphiSoL2L; jso++)
                {
                    ek(iso + iniSoL, jso + iniSoL) += 1.0 * weight * fFOil[3] * phiSoL2L(jso,0) * phiSoL2L(iso,0) * Value;
                }
                
            }
            
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            DebugStop();
        }
            break;
    }
    
    this->ContributeBCInterface(data, datavecleft, weight, ef, bc);
    
    return;
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    if (fSimulationData->IsTwoPhaseQ()) {
        this->ContributeBCInterfaceAlpha(data, datavecleft, weight, ef, bc);
    }
    
    return;
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Swblock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    int Soblock = 3;        // So Oil Saturation needs L2 scalar functions      (phiSoL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSwL2L        = datavecleft[Swblock].phi; // For L2  test functions Sw
    TPZFMatrix<REAL> phiSoL2L        = datavecleft[Soblock].phi; // For L2  test functions So
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSwL2L    = phiSwL2L.Rows();                                   // For L2   Sw
    int nphiSoL2L    = phiSoL2L.Rows();                                   // For L2   So
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSwL   = nphiPL2L       + iniPL;
    int iniSoL   = nphiSwL2L      + iniSwL;
    
    
    //    int iblock = iniSoL + nphiSoL2L;
    //    int jblock = iniSoL + nphiSoL2L;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SwL             = datavecleft[Swblock].sol[0][0];
    REAL SoL             = datavecleft[Soblock].sol[0][0];
    
    //    REAL SgL             = 1.0 - SoL - SwL;
    
    TPZFMatrix<STATE> iphiuHdiv(3,1);
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    REAL  Value, Sw, So;
    Value = bc.Val2()(0,0);         //  un ou P
    Sw = bc.Val2()(1,0);         //  Sw in
    So = bc.Val2()(2,0);         //  So in
    
    
    //  Compute axuliar functions for the current values of time, u, P, Sw and So
    TPZFMatrix<REAL> jphiuHdiv(3,1,0.0);
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            if (uLn >= 0.0 ) {
                // Outflow boundary condition
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            else
            {
                // Inflow boundary condition
                if (uLn < 0.0 && fabs(uLn) < 1.0e-24)
                {
                    std::cout << "Warning: Inflow boundary condition detected in outflow boundary condition: uLn = " << uLn << "\n";
                }
                
                this->UpdateStateVariables(uL, Value, Sw, So);
                this->PhaseFractionalFlows();
            }
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            //            if (uLn >= 0.0 ) {
            //                // Outflow boundary condition
            //                this->UpdateStateVariables(uL, PL, SwL, SoL);
            //                this->PhaseFractionalFlows();
            //            }
            //            else
            //            {
            //                // Inflow boundary condition
            //                if (uLn < 0.0 && fabs(uLn) < 1.0e-24) { std::cout << "Boundary condition error: inflow detected in outflow boundary condition: uLn = " << uLn << "\n";}
            this->UpdateStateVariables(uL, PL, Sw, So);
            this->PhaseFractionalFlows();
            //            }
            
            for (int iq = 0; iq < nphiuHdivL; iq++)
            {
                
                int ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
                int ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
                
                iphiuHdiv(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                iphiuHdiv(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                
                REAL iphiuHdivdotn = iphiuHdiv(0,0)*n[0] + iphiuHdiv(1,0)*n[1];
                
                ef(iq + iniuL) += 1.0 * weight * PL * iphiuHdivdotn;
                
            }
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
            }
            
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            
            if (uLn >= 0.0 ) {
                // Outflow boundary condition
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            else
            {
                // Inflow boundary condition
                if (uLn < 0.0 && fabs(uLn) < 1.0e-24) { std::cout << "Warning: Inflow boundary condition detected in outflow boundary condition: uLn = " << uLn << "\n";}
                this->UpdateStateVariables(uL, Value, SwL, SoL);
                this->PhaseFractionalFlows();
            }
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * uLn;
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * uLn;
                
            }
            
            
        }
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            
            //            if (uLn >= 0.0 ) {
            // Outflow boundary condition
            this->UpdateStateVariables(uL, PL, SwL, SoL);
            this->PhaseFractionalFlows();
            //            }
            //            else
            //            {
            //                // Inflow boundary condition
            //                if (uLn < 0.0 && fabs(uLn) > 1.0e-8) { std::cout << "Boundary condition error: inflow detected in outflow boundary condition: uLn = " << uLn << "\n";}
            //                this->UpdateStateVariables(uL, PL, Sw, So);
            //                this->PhaseFractionalFlows();
            //            }
            
            for (int iq = 0; iq < nphiuHdivL; iq++)
            {
                
                int ivectorindex = datavecleft[ublock].fVecShapeIndex[iq].first;
                int ishapeindex = datavecleft[ublock].fVecShapeIndex[iq].second;
                
                iphiuHdiv(0,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(0,ivectorindex);
                iphiuHdiv(1,0) = phiuH1L(ishapeindex,0) * datavecleft[ublock].fNormalVec(1,ivectorindex);
                
                REAL iphiuHdivdotn = iphiuHdiv(0,0)*n[0] + iphiuHdiv(1,0)*n[1];
                
                ef(iq + iniuL) += 1.0 * weight * PL * iphiuHdivdotn;
                
            }
            
            
            for (int isw = 0; isw < nphiSwL2L; isw++)
            {
                ef(isw + iniSwL) += 1.0 * weight * fFWater[0] * phiSwL2L(isw,0) * Value;
                
            }
            
            for (int iso = 0; iso < nphiSoL2L; iso++)
            {
                ef(iso + iniSoL) += 1.0 * weight * fFOil[0] * phiSoL2L(iso,0) * Value;
                
            }
            break;
        }
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            DebugStop();
        }
            break;
    }
    
    return;
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight,TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    
    this->ContributeBCDarcy(datavec, weight, ek, ef, bc);
    return;
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFNMatrix<9,STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFNMatrix<9,STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = datavec[Qblock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    // Computing normal flux
    STATE Qn = Q[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = PhiH1.Rows();  // For Hdiv
    //    int nPhiL2   = WL2.Rows();                                  // For L2
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    
    
    REAL Value;
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * ( (Value ) * PhiH1(iq,0));
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            
            Value = bc.Val2()(0,0);         //  NormalFlux
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (Qn - Value) ) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (PhiH1(jq,0)) * PhiH1(iq,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ef(iq) += weight * ( (Value ) * PhiH1(iq,0));
                
            }
        }
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (Qn - Value) ) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (PhiH1(jq,0)) * PhiH1(iq,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    return;
    
}




void TPZAxiSymmetricDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdP(2,1);
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    
    Oneoverlambda_Kinv_u(0,0) = (1.0/lambda[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    Oneoverlambda_Kinv_u(1,0) = (1.0/lambda[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    gm(0,0) = rhof[0] * Gravity(0,0);
    gm(1,0) = rhof[0] * Gravity(1,0);
    
    dgmdP(0,0) = rhof[2] * Gravity(0,0);
    dgmdP(1,0) = rhof[2] * Gravity(1,0);
    
    REAL divu = 0.0;
    TPZFMatrix<STATE> iphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    TPZFMatrix<STATE> jphiuHdiv(2,1);
    int jshapeindex;
    int jvectorindex;
    
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * ((Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))
                                   - P * DivergenceOnDeformed(iq,0)
                                   - (gm(0,0)*iphiuHdiv(0,0) + gm(1,0)*iphiuHdiv(1,0)) );
        
        // du/dalphau terms
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            jvectorindex = datavec[ublock].fVecShapeIndex[jq].first;
            jshapeindex = datavec[ublock].fVecShapeIndex[jq].second;
            
            jphiuHdiv(0,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(0,jvectorindex);
            jphiuHdiv(1,0) = phiuH1(jshapeindex,0) * datavec[ublock].fNormalVec(1,jvectorindex);
            
            Oneoverlambda_Kinv_Phiu(0,0) = (1.0/lambda[0]) * (KInverse(0,0)*jphiuHdiv(0,0) + KInverse(0,1)*jphiuHdiv(1,0));
            Oneoverlambda_Kinv_Phiu(1,0) = (1.0/lambda[0]) * (KInverse(1,0)*jphiuHdiv(0,0) + KInverse(1,1)*jphiuHdiv(1,0));
            
            ek(iq + iniu,jq + iniu) +=  weight * ((Oneoverlambda_Kinv_Phiu(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_Phiu(1,0)*iphiuHdiv(1,0)));
        }
        
        // du/dalphap terms
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iq + iniu,jp + iniP) += weight * ( (-lambda[2]/lambda[0])*(Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0)) - DivergenceOnDeformed(iq,0) - (dgmdP(0,0)*iphiuHdiv(0,0) + dgmdP(1,0)*iphiuHdiv(1,0)) )* phiPL2(jp,0) ;
        }
        
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    TPZFMatrix<REAL> Grad;
    if(fTimeDependentForcingFunction)
    {
        fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // Note::  uses this for HdivPiola = 1 or constant Jacobian mappings
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        
        ef(ip + iniP) += weight * (- divu + fvalue[0] - (1.0/dt) * phi * rho[0]) * phiPL2(ip,0);
        
        for (int jq = 0; jq < nphiuHdiv; jq++)
        {
            ek(ip + iniP,jq + iniu) += weight * (- DivergenceOnDeformed(jq,0)) * phiPL2(ip,0);
        }
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(ip + iniP,jp + iniP) += weight * (- (1.0/dt) * phi * rho[2] * phiPL2(jp,0) ) * phiPL2(ip,0);
        }
        
    }
    
    
}


void TPZAxiSymmetricDarcyFlow::ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    TPZFMatrix<STATE> DivergenceOnDeformed;
    // Compute the divergence on deformed element by piola contravariant transformation
    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    
    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    TPZFNMatrix<660,REAL> GradP;
    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    
    
    Oneoverlambda_Kinv_u(0,0) = (1.0/lambda[0])* (KInverse(0,0)*u[0] + KInverse(0,1)*u[1]);
    Oneoverlambda_Kinv_u(1,0) = (1.0/lambda[0])* (KInverse(1,0)*u[0] + KInverse(1,1)*u[1]);
    
    gm(0,0) = rhof[0] * Gravity(0,0);
    gm(1,0) = rhof[0] * Gravity(1,0);
    
    REAL divu = 0.0;
    TPZFMatrix<STATE> iphiuHdiv(2,1);
    int ishapeindex;
    int ivectorindex;
    
    /////////////////////////////////
    // Last State n
    if (fSimulationData->IsnStep()) {
        
        //  n state computations
        for (int ip = 0; ip < nphiPL2; ip++)
        {
            ef(ip + iniP) += 1.0 * weight * (1.0/dt) * phi * rho[0]*  phiPL2(ip,0);
        }
        
        return;
    }
    // Last State n
    /////////////////////////////////
    
    
    for (int iq = 0; iq < nphiuHdiv; iq++)
    {
        
        /* $ \underset{\Omega_{e}}{\int}\left(K\lambda\right)^{-1}\mathbf{q}\cdot\mathbf{v}\;\partial\Omega_{e}-\underset{\Omega_{e}}{\int}P\; div\left(\mathbf{v}\right)\partial\Omega-\underset{\Omega_{e}}{\int}\nabla\left(\rho_{f}g\; z\right)\cdot\mathbf{v} $ */
        
        ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
        ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
        
        iphiuHdiv(0,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(0,ivectorindex);
        iphiuHdiv(1,0) = phiuH1(ishapeindex,0) * datavec[ublock].fNormalVec(1,ivectorindex);
        
        ef(iq + iniu) += weight * ((Oneoverlambda_Kinv_u(0,0)*iphiuHdiv(0,0) + Oneoverlambda_Kinv_u(1,0)*iphiuHdiv(1,0))
                                   - P * DivergenceOnDeformed(iq,0)
                                   - (gm(0,0)*iphiuHdiv(0,0) + gm(1,0)*iphiuHdiv(1,0)) );
        
    }
    
    TPZManVector<STATE,1> fvalue(1,0.0);
    TPZFMatrix<REAL> Grad;
    if(fTimeDependentForcingFunction)
    {
        fTimeDependentForcingFunction->Execute(datavec[Pblock].x,fSimulationData->GetTime(),fvalue,Grad);
    }
    
    
    divu = (Graduaxes(0,0) + Graduaxes(1,1) + Graduaxes(2,2)); // Note::  uses this for constant jacobian elements
    
    /* $ - \underset{\Omega}{\int}w\; div\left(\mathbf{q}\right)\partial\Omega $ */
    for (int ip = 0; ip < nphiPL2; ip++)
    {
        ef(ip + iniP) += weight * (- divu + fvalue[0] - (1.0/dt) * phi * rho[0]) * phiPL2(ip,0);
        
    }
    
    return;
    
}


void TPZAxiSymmetricDarcyFlow::ContributeBCDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // At each Integration Point.
    
    int Qblock = 0;
    int Pblock = 1;
    
    // Getting test and basis functions
    TPZFNMatrix<9,STATE> PhiH1 = datavec[Qblock].phi; // For H1   test functions
    TPZFNMatrix<9,STATE> WL2   = datavec[Pblock].phi; // For HL2  test functions
    
    TPZFMatrix<STATE> dPhiH1 = datavec[Qblock].dphix; // Derivative For H1   test functions
    TPZFMatrix<STATE> dWL2   = datavec[Pblock].dphix; // Derivative For HL2  test functions
    
    // Getting Linear combinations of basis functions
    TPZManVector<STATE> Q = datavec[Qblock].sol[0];
    TPZManVector<STATE> P = datavec[Pblock].sol[0];
    
    // Computing normal flux
    REAL Qn = Q[0];
    
    TPZFMatrix<STATE> dQdx = datavec[Qblock].dsol[0];
    TPZFMatrix<STATE> dPdx = datavec[Pblock].dsol[0];
    
    // Number of phis
    int nPhiHdiv = PhiH1.Rows();  // For Hdiv
    //    int nPhiL2   = WL2.Rows();                                  // For L2
    
    TPZFMatrix<STATE> iPhiHdiv(2,1);
    TPZFMatrix<STATE> jPhiHdiv(2,1);
    
    
    
    REAL Value;
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * ( (Value ) * PhiH1(iq,0));
            }
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            
            
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            //             TPZManVector<STATE,1> fvalue(1,0.0);
            //             TPZFMatrix<REAL> Grad;
            //             if(fTimedependentBCForcingFunction)
            //             {
            //                 fTimedependentBCForcingFunction->Execute(datavec[Qblock].x,fSimulationData->GetTime(),fvalue,Grad);
            //                 Value = fvalue[0];
            //             }
        
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (Qn - Value)) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (PhiH1(jq,0)) * PhiH1(iq,0);
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            
            Value = bc.Val2()(0,0);         //  Pressure
            //             TPZManVector<STATE,1> fvalue(1,0.0);
            //             TPZFMatrix<REAL> Grad;
            //             if(fTimedependentFunctionExact)
            //             {
            //                 fTimedependentFunctionExact->Execute(datavec[Qblock].x,fSimulationData->GetTime(),fvalue,Grad);
            //                 Value = fvalue[0];
            //             }
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                
                ef(iq) += weight * ( (Value ) * PhiH1(iq,0));
                
            }
        }
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            Value = bc.Val2()(0,0);         //  NormalFlux
            
            for (int iq = 0; iq < nPhiHdiv; iq++)
            {
                ef(iq) += weight * (gBigNumber * (Qn - Value) ) * PhiH1(iq,0);
                
                for (int jq = 0; jq < nPhiHdiv; jq++)
                {
                    ek(iq,jq) += gBigNumber * weight * (PhiH1(jq,0)) * PhiH1(iq,0);
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            DebugStop();
        }
            break;
    }
    
    return;
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Salphablock = 2;         // Salpha Saturation needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> Sa_phiPL2         = datavec[Salphablock].phi;  // For L2  test functions S_alpha
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    //    TPZFMatrix<STATE> DivergenceOnDeformed;
    //    // Compute the divergence on deformed element by piola contravariant transformation
    //    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSaL2     = Sa_phiPL2.Rows();                                    // For L2   S_alpha
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSa    = nphiPL2     + iniP;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL S_alpha              = datavec[Salphablock].sol[0][0];
    //
    //    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    //    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    //    TPZFNMatrix<660,REAL> GradP;
    //    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdP(2,1);
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    
    
    for (int isw = 0; isw < nphiSaL2; isw++)
    {
        
        ef(isw  + iniSa ) += 1.0 * weight * (1.0/dt) * phi * S_alpha * rho_alpha[0]  * Sa_phiPL2(isw,0);
        
        for (int jp = 0; jp < nphiSaL2; jp++)
        {
            ek(isw  + iniSa, jp  + iniP ) += 1.0 * weight * (1.0/dt) * phi * S_alpha * rho_alpha[2] * phiPL2(jp,0) * Sa_phiPL2(isw,0);
            
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(isw  + iniSa, jsw  + iniSa ) += 1.0 * weight * (1.0/dt) * phi * Sa_phiPL2(jsw,0) * rho_alpha[0]  * Sa_phiPL2(isw,0);
            
        }
    }
    
}

/**
 * It computes a contribution to the load vector at one integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    
    
    // Getting data for Mixed-Darcy flow problem
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Salphablock = 2;         // Salpha Saturation needs L2 scalar functions     (phiPL2)
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2         = datavec[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> Sa_phiPL2         = datavec[Salphablock].phi;  // For L2  test functions S_alpha
    
    TPZFMatrix<REAL> dphiuH1   = datavec[ublock].dphix; // Derivative For H1  test functions u
    TPZFMatrix<REAL> dphiPL2   = datavec[Pblock].dphix; // Derivative For L2  test functions P
    
    //    TPZFMatrix<STATE> DivergenceOnDeformed;
    //    // Compute the divergence on deformed element by piola contravariant transformation
    //    this->ComputeDivergenceOnDeformed(datavec, DivergenceOnDeformed);
    
    
    // Blocks dimensions and lengths
    int nphiuHdiv   = datavec[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2     = phiPL2.Rows();                                    // For L2   P
    int nphiSaL2     = Sa_phiPL2.Rows();                                    // For L2   S_alpha
    int iniu    = 0;
    int iniP    = nphiuHdiv     + iniu;
    int iniSa    = nphiPL2     + iniP;
    
    // Getting linear combinations from different approximation spaces
    TPZManVector<REAL,3> u      = datavec[ublock].sol[0];
    REAL P              = datavec[Pblock].sol[0][0];
    REAL S_alpha              = datavec[Salphablock].sol[0][0];
    //
    //    TPZFMatrix<STATE> Graduaxes = datavec[ublock].dsol[0];
    //    TPZFMatrix<STATE> GradPaxes = datavec[Pblock].dsol[0];
    //    TPZFNMatrix<660,REAL> GradP;
    //    TPZAxesTools<REAL>::Axes2XYZ(GradPaxes, GradP, datavec[Pblock].axes);
    
    
    // Rock and fluids parameters
    TPZFMatrix<STATE> KInverse = fReservoirdata->KabsoluteInv();
    TPZFMatrix<STATE> Oneoverlambda_Kinv_u(2,1);
    TPZFMatrix<REAL> Oneoverlambda_Kinv_Phiu(2,1);
    TPZFMatrix<STATE> Gravity(2,1);
    TPZFMatrix<STATE> gm(2,1);
    TPZFMatrix<STATE> dgmdP(2,1);
    
    Gravity = fSimulationData->GetGravity();
    
    REAL phi, dphidP;
    this->fReservoirdata->Porosity(P, phi, dphidP);
    
    // time
    REAL dt = fSimulationData->GetDeltaT();
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    this->ComputeProperties(datavec, props);
    TPZManVector<REAL> rho_alpha        = props[0];
    TPZManVector<REAL> rho_beta         = props[1];
    TPZManVector<REAL> rho_gamma        = props[2];
    TPZManVector<REAL> f_alpha          = props[3];
    TPZManVector<REAL> f_beta           = props[4];
    TPZManVector<REAL> f_gamma          = props[5];
    TPZManVector<REAL> lambda           = props[6];
    TPZManVector<REAL> rho              = props[7];
    TPZManVector<REAL> rhof             = props[8];
    
    /////////////////////////////////
    // Last State n
    if (fSimulationData->IsnStep()) {
        
        for (int isw = 0; isw < nphiSaL2; isw++)
        {
            ef(isw  + iniSa ) += -1.0 * weight * (1.0/dt) * phi * S_alpha * rho_alpha[0]  * Sa_phiPL2(isw,0);
        }
        
        return;
    }
    // Last State n
    /////////////////////////////////
    
    
    for (int isw = 0; isw < nphiSaL2; isw++)
    {
        ef(isw  + iniSa ) += 1.0 * weight * (1.0/dt) * phi * S_alpha * rho_alpha[0]  * Sa_phiPL2(isw,0);
    }
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SaL             = datavecleft[Sablock].sol[0][0];
    
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
//
//                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
//                    {
//                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                }
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * qn;
//                    }
//
//                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
//                    {
//                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * qn;
//                    }
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " <<  uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
                
            }
            
        }
            
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                
//                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
//                state_vars[0] = qn;                             //  qn
//                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
//                state_vars[2] = S_alpha;                        //  S_alpha
//                state_vars[3] = S_beta;                         //  S_beta
                
//                this->ComputeProperties(state_vars, props);
//                TPZManVector<REAL> f_alpha          = props[3];
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * qn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * qn;
                    }
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " << uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
//                    for (int jp = 0; jp < nphiPL2L; jp++)
//                    {
//                        ek(isw + iniSaL,jp + iniPL) += 1.0 * weight * f_alpha[2] * phiPL2L(jp,0)  * phiSaL2L(isw,0) * uLn;
//                    }
                    
                    for (int jsw = 0; jsw < nphiSaL2L; jsw++)
                    {
                        ek(isw + iniSaL,jsw + iniSaL) += 1.0 * weight * f_alpha[3] * phiSaL2L(jsw,0)  * phiSaL2L(isw,0) * uLn;
                    }
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            DebugStop();
        }
            break;
    }
    
    return;
}


/**
 * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{


    // Getting data from different approximation spaces

    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)

    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions


    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;


    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SaL             = datavecleft[Sablock].sol[0][0];


    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];


    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  PD inflow
        {
            
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN inflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                //                std::cout << "Outflow condition detected in inflow condition qn = " << uLn << std::endl;
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                    
                }
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD outflow
        {
            REAL PD         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
                
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = uLn;                            //  qn
                state_vars[1] = PD;                             //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                //                this->ComputeProperties(state_vars, props);
                //                TPZManVector<REAL> f_alpha          = props[3];
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " <<  uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            
            break;
            
        case 3 :    // Neumann BC  QN outflow
        {
            REAL qn         = bc.Val2()(0,0);
            REAL S_alpha    = bc.Val2()(1,0);
            REAL S_beta     = bc.Val2()(2,0);
            
            // Returning needed relationships
            TPZVec< TPZManVector<REAL>  > props;
            
            if ((uLn >= 0.0) && fabs(uLn) > fepsilon ) {
                
                TPZManVector<REAL> state_vars(fnstate_vars+2,0.0);
                state_vars[0] = qn;                             //  qn
                state_vars[1] = datavecleft[Pblock].sol[0][0];  //  P
                state_vars[2] = S_alpha;                        //  S_alpha
                state_vars[3] = S_beta;                         //  S_beta
                
                this->ComputeProperties(state_vars, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                //                this->ComputeProperties(datavecleft, props);
                //                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * qn;
                    
                }
            }
            else
            {
                std::cout << "Inflow  condition detected in outflow condition qn = " << uLn << std::endl;
                
                this->ComputeProperties(datavecleft, props);
                TPZManVector<REAL> f_alpha          = props[3];
                
                for (int isw = 0; isw < nphiSaL2L; isw++)
                {
                    ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
                    
                }
                
            }
            
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            DebugStop();
        }
            break;
    }

    return;


}


/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ek[out] is the stiffness matrix
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    // Full implicit case: there is no n state computations here
    if (fSimulationData->IsnStep()) {
        return;
    }
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSaR   = nphiPL2R       + iniPR;
    
    int iblock = iniSaL + nphiSaL2L;
    int iblockt = 0;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SaL             = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL PR              = datavecright[Pblock].sol[0][0];
    REAL SaR             = datavecright[Sablock].sol[0][0];
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1;
    TPZFMatrix<REAL> phiPL2;
    TPZFMatrix<REAL> phiSaL2;
    TPZFMatrix<STATE> dphiuH1;
    TPZFMatrix<STATE> dphiPL2;
    int nphiuHdiv   = 0;
    int nphiPL2     = 0;
    int nphiSaL2    = 0;
    int iniu    = 0;
    int iniP    = 0;
    int iniSa   = 0;
    
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    
    if ((uLn >= 0.0) && fabs(uLn) > fepsilon) {
        this->ComputeProperties(datavecleft, props);
        phiuH1      = phiuH1L;
        phiPL2      = phiPL2L;
        phiSaL2     = phiSaL2L;
        dphiuH1     = dphiuH1L;
        dphiPL2     = dphiPL2L;
        nphiuHdiv   = nphiuHdivL;
        nphiPL2     = nphiPL2L;
        nphiSaL2    = nphiSaL2L;
        iniu        = iniuL;
        iniP        = iniPL;
        iniSa       = iniSaL;
        iblockt      = 0;
    }
    else
    {
        this->ComputeProperties(datavecright, props);
        phiuH1      = phiuH1R;
        phiPL2      = phiPL2R;
        phiSaL2     = phiSaL2R;
        dphiuH1     = dphiuH1R;
        dphiPL2     = dphiPL2R;
        nphiuHdiv   = nphiuHdivR;
        nphiPL2     = nphiPL2R;
        nphiSaL2    = nphiSaL2R;
        iniu        = iniuR;
        iniP        = iniPR;
        iniSa       = iniSaR;
        iblockt     = iblock;
    }
    
    TPZManVector<REAL> f_alpha          = props[3];
    
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(isw + iniSaL,iblockt + jp + iniP) += 1.0 * weight * f_alpha[2] * phiPL2(jp,0)  * phiSaL2L(isw,0) * uLn;
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(isw + iniSaL,iblockt + jsw + iniSa) += 1.0 * weight * f_alpha[3] * phiSaL2(jsw,0)  * phiSaL2L(isw,0) * uLn;
        }
    }
    
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * f_alpha[0] * phiSaL2R(isw,0) * uLn;
        
        for (int jp = 0; jp < nphiPL2; jp++)
        {
            ek(iblock + isw + iniSaR, iblockt + jp + iniP ) += -1.0 * weight * f_alpha[2] * phiPL2(jp,0) * phiSaL2R(isw,0) * uLn;
        }
        
        for (int jsw = 0; jsw < nphiSaL2; jsw++)
        {
            ek(iblock + isw + iniSaR, iblockt + jsw + iniSa) += -1.0 * weight * f_alpha[3] * phiSaL2(jsw,0) * phiSaL2R(isw,0) * uLn;
        }
        
    }
    
    return;
    
}

/**
 * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
 * @param data[in] stores all input data
 * @param weight[in] is the weight of the integration rule
 * @param ef[out] is the load vector
 * @param bc[in] is the boundary condition material
 * @since April 16, 2007
 */
void TPZAxiSymmetricDarcyFlow::ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef){
    
    
    // Getting data from different approximation spaces
    
    int ublock = 0;         // u Bulk velocity needs H1 scalar functions        (phiuH1) for the construction of Hdiv basis functions phiuHdiv
    int Pblock = 1;         // P Average Pressure needs L2 scalar functions     (phiPL2)
    int Sablock = 2;        // Sw Water Saturation needs L2 scalar functions    (phiSwL2)
    
    // Getting test and basis functions for the left side
    TPZFMatrix<REAL> phiuH1L         = datavecleft[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2L         = datavecleft[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2L        = datavecleft[Sablock].phi; // For L2  test functions Sw
    TPZFMatrix<STATE> dphiuH1L   = datavecleft[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2L   = datavecleft[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Getting test and basis functions for the right side
    TPZFMatrix<REAL> phiuH1R         = datavecright[ublock].phi;  // For H1  test functions u
    TPZFMatrix<REAL> phiPL2R         = datavecright[Pblock].phi;  // For L2  test functions P
    TPZFMatrix<REAL> phiSaL2R        = datavecright[Sablock].phi; // For L2  test functions Sa
    TPZFMatrix<STATE> dphiuH1R   = datavecright[ublock].dphix; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiPL2R   = datavecright[Pblock].dphix; // Derivative For L2  test functions
    
    
    // Blocks dimensions and lengths for the left side
    int nphiuHdivL   = datavecleft[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2L     = phiPL2L.Rows();                                    // For L2   P
    int nphiSaL2L    = phiSaL2L.Rows();                                   // For L2   Sa
    int iniuL    = 0;
    int iniPL    = nphiuHdivL     + iniuL;
    int iniSaL   = nphiPL2L       + iniPL;
    
    // Blocks dimensions and lengths for the right side
    int nphiuHdivR   = datavecright[ublock].fVecShapeIndex.NElements();       // For Hdiv u
    int nphiPL2R     = phiPL2R.Rows();                                    // For L2   P
    int nphiSaL2R    = phiSaL2R.Rows();                                   // For L2   Sa
    int iniuR    = 0;
    int iniPR    = nphiuHdivR     + iniuR;
    int iniSaR   = nphiPL2R       + iniPR;
    
    int iblock = iniSaL + nphiSaL2L;
    
    // Getting linear combinations from different approximation spaces for the left side
    TPZManVector<REAL,3> uL      = datavecleft[ublock].sol[0];
    REAL PL              = datavecleft[Pblock].sol[0][0];
    REAL SaL             = datavecleft[Sablock].sol[0][0];
    
    // Getting linear combinations from different approximation spaces for the right side
    TPZManVector<REAL,3> uR      = datavecright[ublock].sol[0];
    REAL PR              = datavecright[Pblock].sol[0][0];
    REAL SaR             = datavecright[Sablock].sol[0][0];
    
    TPZManVector<REAL,3> n =  data.normal;
    REAL uLn = uL[0]*n[0] + uL[1]*n[1] + uL[2]*n[2];
    
    
    // Returning needed relationships
    TPZVec< TPZManVector<REAL>  > props;
    
    if ((uLn >= 0.0) && fabs(uLn) > fepsilon)
    {
        this->ComputeProperties(datavecleft, props);
    }
    else
    {
        this->ComputeProperties(datavecright, props);
    }
    
    TPZManVector<REAL> f_alpha          = props[3];
    
    
    for (int isw = 0; isw < nphiSaL2L; isw++)
    {
        ef(isw + iniSaL) += 1.0 * weight * f_alpha[0] * phiSaL2L(isw,0) * uLn;
    }
    
    for (int isw = 0; isw < nphiSaL2R; isw++)
    {
        ef(iblock + isw + iniSaR) += -1.0 * weight * f_alpha[0] * phiSaL2R(isw,0) * uLn;
    }
    
    
}



void TPZAxiSymmetricDarcyFlow::UpdateStateVariables(TPZManVector<REAL,3> u, REAL P, REAL Sw, REAL So)
{
    fBulkVelocity       = u;
    fAveragePressure    = P;
    fWaterSaturation    = Sw;
    fOilSaturation      = So;
}

// systematic methods computed on the current u, p, Sw and So values.



void TPZAxiSymmetricDarcyFlow::PhasePressures()
{
    // Petrophysics data
    
    REAL Sw, So, Sg, P;
    P = fAveragePressure;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    REAL pcow = 0, pcgo = 0, pcgw = 0;
    REAL dPcowdSw = 0, dPcgodSo = 0;
    
    //    this->fPetrophysicdata->Pcow(Sw, pcow, dPcowdSw);
    //    this->fPetrophysicdata->Pcgo(So, pcgo, dPcgodSo);
    //    this->fPetrophysicdata->Pcgw(Sw, pcgw, dPcgwdSw); // or Pcgo(So) + Pcow(Sw)
    
    // Recovering phase pressures and their derivatives
    
    fWaterPressure[0] = P - So * pcow - (1.0 - So - Sw) * pcgw;
    fWaterPressure[1] = 1.0;
    fWaterPressure[2] = pcgo - So * dPcowdSw + pcow - Sg * dPcowdSw;
    fWaterPressure[3] = pcgo - Sg * dPcgodSo;
    
    fOilPressure[0] = P + Sw * pcow - (1.0 - So - Sw) * pcgo;
    fOilPressure[1] = 1.0;
    fOilPressure[2] = Sw * dPcowdSw + pcow + pcgo;
    fOilPressure[3] = pcgo - Sg * dPcgodSo;
    
    // Computing fluid mobilities and derivatives for two-phase
    
    fGasPressure[0] = 0.0;
    fGasPressure[1] = 0.0;
    fGasPressure[2] = 0.0;
    fGasPressure[3] = 0.0;
    
    
    //    // Computing fluid mobilities and derivatives for three-phase
    //
    //    fGasPressure[0] = P + So * pcgo + So * pcgw;
    //    fGasPressure[1] = 1.0;
    //    fGasPressure[2] = pcgo + Sw * dPcowdSw + pcow;
    //    fGasPressure[3] = pcgo + So * dPcgodSo + Sw * dPcgodSo;
    
}

void TPZAxiSymmetricDarcyFlow::PhaseDensities()
{
    
    this->PhasePressures();
    
    //    fFluidmodeldata->WaterDensity(fWaterPressure[0], fWaterDensity[0], fWaterDensity[1]);
    //    fFluidmodeldata->OilDensity(fOilPressure[0], fOilDensity[0], fOilDensity[1]);
    //    fFluidmodeldata->GasDensity(fGasPressure[0], fGasDensity[0], fGasDensity[1]);
    
    // Chain Rule
    fWaterDensity[1]    *= fWaterPressure[1];
    fOilDensity[1]      *= fOilPressure[1];
    fGasDensity[1]      *= fGasPressure[1];
}

void TPZAxiSymmetricDarcyFlow::PhaseMobilities()
{
    REAL Sw, So, Sg;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    REAL waterdensity = 0, oildensity= 0, gasdensity= 0;
    REAL dwaterdensitydP= 0, doildensitydP= 0, dgasdensitydP= 0;
    
    REAL waterviscosity= 0, oilviscosity= 0;
    REAL dwaterviscositydPw= 0, doilviscositydPo= 0;
    
    REAL krw= 0, kro= 0;
    REAL dkrwdSw= 0, dkrodSo= 0;
    
    
    this->PhaseDensities();
    waterdensity        = fWaterDensity[0];
    dwaterdensitydP     = fWaterDensity[1];
    oildensity          = fOilDensity[0];
    doildensitydP       = fOilDensity[1];
    gasdensity          = fGasDensity[0];
    dgasdensitydP       = fGasDensity[1];
    
    //    this->fFluidmodeldata->WaterViscosity(fWaterPressure[0], waterviscosity, dwaterviscositydPw);
    //    this->fFluidmodeldata->OilViscosity(fOilPressure[0], oilviscosity, doilviscositydPo);
    //    this->fFluidmodeldata->GasViscosity(fGasPressure[0], gasviscosity, dgasviscositydPg);
    //
    //    this->fPetrophysicdata->Krw(Sw, krw, dkrwdSw);
    //    this->fPetrophysicdata->Kro(1.0-Sw, kro, dkrodSo); // here appears the two-phase dependence
    //    this->fPetrophysicdata->Krg(1.0-(1.0-Sw)-Sw, krg, dkrgdSg);   // here appears the two-phase dependence
    
    // Computing fluid mobilities and derivatives for two-phase
    
    fWaterMobility[0] = waterdensity * krw / waterviscosity;
    fWaterMobility[1] = krw*(dwaterdensitydP/waterviscosity) - krw*((waterdensity * dwaterviscositydPw)/(waterviscosity*waterviscosity));
    fWaterMobility[2] = 1.0*dkrwdSw*(waterdensity/waterviscosity);
    fWaterMobility[3] = (-0.0)*dkrwdSw*(waterdensity/waterviscosity);   // here appears the two-phase dependence
    
    fOilMobility[0] = oildensity * kro / oilviscosity;
    fOilMobility[1] = kro*(doildensitydP/oilviscosity) - kro*((oildensity * doilviscositydPo)/(oilviscosity*oilviscosity));
    fOilMobility[2] = (-1.0)*dkrodSo*(oildensity/oilviscosity);         // here appears the two-phase dependence
    fOilMobility[3] = 0.0*dkrodSo*(oildensity/oilviscosity);
    
    fGasMobility[0] = 0.0;
    fGasMobility[1] = 0.0;
    fGasMobility[2] = 0.0;
    fGasMobility[3] = 0.0;
    
    //    // Computing fluid mobilities and derivatives for three-phase
    //
    //    fWaterMobility[0] = waterdensity * krw / waterviscosity;
    //    fWaterMobility[1] = krw*(dwaterdensitydP/waterviscosity) - krw*((waterdensity * dwaterviscositydPw)/(waterviscosity*waterviscosity));
    //    fWaterMobility[2] = dkrwdSw*(waterdensity/waterviscosity);
    //    fWaterMobility[3] = 0.0;
    //
    //    fOilMobility[0] = oildensity * kro / oilviscosity;
    //    fOilMobility[1] = kro*(doildensitydP/oilviscosity) - kro*((oildensity * doilviscositydPo)/(oilviscosity*oilviscosity));
    //    fOilMobility[2] = 0.0;
    //    fOilMobility[3] = dkrodSo*(oildensity/oilviscosity);
    //
    //    fGasMobility[0] = gasdensity * krg / gasviscosity;
    //    fGasMobility[1] = krg*(dgasdensitydP/gasviscosity) - krw*((gasdensity * dgasviscositydPg)/(gasviscosity*gasviscosity));
    //    fGasMobility[2] = (dkrgdSg*(-1.0)*(gasdensity/gasviscosity));
    //    fGasMobility[3] = (dkrgdSg*(-1.0)*(gasdensity/gasviscosity));
    
}

void TPZAxiSymmetricDarcyFlow::TotalDensity()
{
    REAL Sw, So, Sg;
    Sw = fWaterSaturation;
    So = fOilSaturation;
    Sg = 1.0 - Sw - So;
    
    this->PhaseDensities();
    
    fTotalDensity[0] = fWaterDensity[0] * Sw + fOilDensity[0] * So + fGasDensity[0] * (1.0 - Sw - So);
    fTotalDensity[1] = fWaterDensity[1] * Sw + fOilDensity[1] * So + fGasDensity[1] * (1.0 - Sw - So);
    fTotalDensity[2] = fWaterDensity[0] - fGasDensity[0];
    fTotalDensity[3] = fOilDensity[0]   - fGasDensity[0];
    
}

void TPZAxiSymmetricDarcyFlow::TotalMobility()
{
    this->PhaseMobilities();
    
    fTotalMobility[0] = fWaterMobility[0] + fOilMobility[0] + fGasMobility[0];
    fTotalMobility[1] = fWaterMobility[1] + fOilMobility[1] + fGasMobility[1];
    fTotalMobility[2] = fWaterMobility[2] + fOilMobility[2] + fGasMobility[2];
    fTotalMobility[3] = fWaterMobility[3] + fOilMobility[3] + fGasMobility[3];
    
    
}

void TPZAxiSymmetricDarcyFlow::PhaseFractionalFlows()
{
    this->TotalMobility();
    
    fFWater[0] = (fWaterMobility[0] / fTotalMobility[0]);
    fFWater[1] = (fWaterMobility[1] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFWater[2] = (fWaterMobility[2] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFWater[3] = (fWaterMobility[3] / fTotalMobility[0]) - ((fWaterMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
    fFOil[0] = (fOilMobility[0] / fTotalMobility[0]);
    fFOil[1] = (fOilMobility[1] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFOil[2] = (fOilMobility[2] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFOil[3] = (fOilMobility[3] / fTotalMobility[0]) - ((fOilMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
    fFGas[0] = (fGasMobility[0] / fTotalMobility[0]);
    fFGas[1] = (fGasMobility[1] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[1])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFGas[2] = (fGasMobility[2] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[2])/(fTotalMobility[0]*fTotalMobility[0])) ;
    fFGas[3] = (fGasMobility[3] / fTotalMobility[0]) - ((fGasMobility[0] * fTotalMobility[3])/(fTotalMobility[0]*fTotalMobility[0])) ;
    
}



int TPZAxiSymmetricDarcyFlow::ClassId() const {
    return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    
}

// -------------------------------------------------------------------------------------------

void TPZAxiSymmetricDarcyFlow::Read(TPZStream &buf, void *context) {
    TPZDiscontinuousGalerkin::Read(buf, context);
    
}



// System Properties

void TPZAxiSymmetricDarcyFlow::ComputeProperties(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL> > & props){
    
    int nphases = fnstate_vars - 1;
    int n = fnstate_vars + 2;
    int p       = 1;
    int s_alpha = 2;
    //    int s_beta  = 3;
    
    TPZManVector<REAL> state_vars(fnstate_vars+1,0.0);
    props.Resize(9);
    
    TPZManVector<REAL> rho_alpha(n,0.0);
    TPZManVector<REAL> mu_alpha(n,0.0);
    TPZManVector<REAL> kr_alpha(n,0.0);
    TPZManVector<REAL> l_alpha(n,0.0);
    TPZManVector<REAL> f_alpha(n,0.0);
    TPZManVector<REAL> Pc_alpha(n,0.0);
    TPZManVector<REAL> vars_alpha(n-1,0.0);
    
    TPZManVector<REAL> rho_beta(n,0.0);
    TPZManVector<REAL> mu_beta(n,0.0);
    TPZManVector<REAL> kr_beta(n,0.0);
    TPZManVector<REAL> l_beta(n,0.0);
    TPZManVector<REAL> f_beta(n,0.0);
    TPZManVector<REAL> Pc_beta(n,0.0);
    TPZManVector<REAL> vars_beta(n-1,0.0);
    
    TPZManVector<REAL> rho_gamma(n,0.0);
    TPZManVector<REAL> mu_gamma(n,0.0);
    TPZManVector<REAL> kr_gamma(n,0.0);
    TPZManVector<REAL> l_gamma(n,0.0);
    TPZManVector<REAL> f_gamma(n,0.0);
    TPZManVector<REAL> Pc_gamma(n,0.0);
    TPZManVector<REAL> vars_gamma(n-1,0.0);
    
    // Strong formulation constitutive functions
    TPZManVector<REAL> rho(n,0.0);
    TPZManVector<REAL> rhof(n,0.0);
    TPZManVector<REAL> l(n,0.0);
    
    
    switch (nphases) {
        case 1:
        {
            // Getting state variables
            REAL P = datavec[p].sol[0][0];
            
            // Compute P_alpha
            state_vars[1] = P;
            REAL P_alpha = P;
            state_vars[1] = P_alpha;
            
            fluid_alpha->Density(rho_alpha, state_vars);
            fluid_alpha->Viscosity(mu_alpha, state_vars);
            
            // Rho computation
            rho = rho_alpha;
            
            // Rhof computation
            rhof = rho_alpha;
            
            // lambda computation
            l[0] = rho_alpha[0]/mu_alpha[0];
            l[2] = rho_alpha[2]/mu_alpha[0] - (rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]);
            
        }
            break;
        case 2:
        {
            
            // Getting state variables
            REAL P          = datavec[p].sol[0][0];
            REAL S_alpha    = datavec[s_alpha].sol[0][0];
            REAL S_beta     = 1.0 - S_alpha;
            
            vars_alpha[1] = P;
            vars_alpha[2] = S_alpha;

            vars_beta[1] = P;
            vars_beta[2] = S_beta;
            
            // Compute P_alpha  and  P_beta
            fluid_beta->Pc(Pc_beta, state_vars);
            vars_alpha[1]   = P - S_alpha   * Pc_beta[1];
            vars_beta[1]    = P + S_beta    * Pc_beta[1];
            
            
            
            fluid_alpha->Density(rho_alpha, vars_alpha);
            fluid_alpha->Viscosity(mu_alpha, vars_alpha);
            fluid_alpha->Kr(kr_alpha, vars_alpha);
            
            
            fluid_beta->Density(rho_beta, vars_beta);
            fluid_beta->Viscosity(mu_beta, vars_beta);
            fluid_beta->Kr(kr_beta, vars_beta);
            
            // Rho computation
            rho[0] = rho_alpha[0] * S_alpha + rho_beta[0] * S_beta;
            rho[1] = rho_alpha[1] * S_alpha + rho_beta[1] * S_beta;
            rho[2] = rho_alpha[2] * S_alpha + rho_beta[2] * S_beta;
            rho[3] = rho_alpha[3] * S_alpha + rho_beta[3] * S_beta;
            
            // Two-Phase restriction
            kr_beta[3] *= -1.0;
            
            // Fluid Mobilities
            
            l_alpha[0] = kr_alpha[0]*(rho_alpha[0]/mu_alpha[0]);
            l_alpha[1] = 0.0;
            l_alpha[2] = kr_alpha[0]*(rho_alpha[2]/mu_alpha[0]-(rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]));
            l_alpha[3] = kr_alpha[3]*(rho_alpha[0]/mu_alpha[0]);
            
            l_beta[0] = kr_beta[0]*(rho_beta[0]/mu_beta[0]);
            l_beta[1] = 0.0;
            l_beta[2] = kr_beta[0]*(rho_beta[2]/mu_beta[0]-(rho_beta[0]*mu_beta[2])/(mu_beta[0]*mu_beta[0]));
            l_beta[3] = kr_beta[3]*(rho_beta[0]/mu_beta[0]);
            
            // lambda computation
            l[0] = l_alpha[0] + l_beta[0];
            l[1] = l_alpha[1] + l_beta[1];
            l[2] = l_alpha[2] + l_beta[2];
            l[3] = l_alpha[3] + l_beta[3];
            
            // Fractional flows
            f_alpha[0] = l_alpha[0]/l[0];
            f_alpha[1] = l_alpha[1]/l[0]-(l_alpha[0]*l[1])/(l[0]*l[0]);
            f_alpha[2] = l_alpha[2]/l[0]-(l_alpha[0]*l[2])/(l[0]*l[0]);
            f_alpha[3] = l_alpha[3]/l[0]-(l_alpha[0]*l[3])/(l[0]*l[0]);
            
            f_beta[0] = l_beta[0]/l[0];
            f_beta[1] = l_beta[1]/l[0]-(l_beta[0]*l[1])/(l[0]*l[0]);
            f_beta[2] = l_beta[2]/l[0]-(l_beta[0]*l[2])/(l[0]*l[0]);
            f_beta[3] = l_beta[3]/l[0]-(l_beta[0]*l[3])/(l[0]*l[0]);
            
            // Rhof computation
            rhof[0] = f_alpha[0] * rho_alpha[0] +  f_beta[0] * rho_beta[0];
            rhof[1] = f_alpha[0] * rho_alpha[1] + f_alpha[1] * rho_alpha[0]  + f_beta[0] * rho_beta[1] + f_beta[1] * rho_beta[0];
            rhof[2] = f_alpha[0] * rho_alpha[2] + f_alpha[2] * rho_alpha[0]  + f_beta[0] * rho_beta[2] + f_beta[2] * rho_beta[0];
            rhof[3] = f_alpha[0] * rho_alpha[3] + f_alpha[3] * rho_alpha[0]  + f_beta[0] * rho_beta[3] + f_beta[3] * rho_beta[0];
            
        }
            break;
        case 3:
        {
            std::cout<< "3-phases system not implemented"<< std::endl;
            DebugStop();
            
        }
            break;
        default:
        {
            std::cout<< "4-phases system not implemented"<< std::endl;
            DebugStop();
        }
            break;
    }
    
    // Returning needed relationships
    
    props[0] = rho_alpha;
    props[1] = rho_beta;
    props[2] = rho_gamma;
    props[3] = f_alpha;
    props[4] = f_beta;
    props[5] = f_gamma;
    props[6] = l;
    props[7] = rho;
    props[8] = rhof;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::ComputeProperties(TPZManVector<REAL> &state_vars, TPZVec<TPZManVector<REAL> > & props){
    
    int nphases = fnstate_vars - 1;
    int n = fnstate_vars + 2;
    int p       = 1;
    int s_alpha = 2;
    //    int s_beta  = 3;
    
    props.Resize(9);
    
    TPZManVector<REAL> rho_alpha(n,0.0);
    TPZManVector<REAL> mu_alpha(n,0.0);
    TPZManVector<REAL> kr_alpha(n,0.0);
    TPZManVector<REAL> l_alpha(n,0.0);
    TPZManVector<REAL> f_alpha(n,0.0);
    TPZManVector<REAL> Pc_alpha(n,0.0);
    TPZManVector<REAL> vars_alpha = state_vars;
    
    TPZManVector<REAL> rho_beta(n,0.0);
    TPZManVector<REAL> mu_beta(n,0.0);
    TPZManVector<REAL> kr_beta(n,0.0);
    TPZManVector<REAL> l_beta(n,0.0);
    TPZManVector<REAL> f_beta(n,0.0);
    TPZManVector<REAL> Pc_beta(n,0.0);
    TPZManVector<REAL> vars_beta = state_vars;
    
    TPZManVector<REAL> rho_gamma(n,0.0);
    TPZManVector<REAL> mu_gamma(n,0.0);
    TPZManVector<REAL> kr_gamma(n,0.0);
    TPZManVector<REAL> l_gamma(n,0.0);
    TPZManVector<REAL> f_gamma(n,0.0);
    TPZManVector<REAL> Pc_gamma(n,0.0);
    TPZManVector<REAL> vars_gamma = state_vars;
    
    // Strong formulation constitutive functions
    TPZManVector<REAL> rho(n,0.0);
    TPZManVector<REAL> rhof(n,0.0);
    TPZManVector<REAL> l(n,0.0);
    
    
    switch (nphases) {
        case 1:
        {
            // Getting state variables
            REAL P = state_vars[1];
            
            // Compute P_alpha
            REAL P_alpha = P;
            state_vars[1] = P_alpha;
            
            fluid_alpha->Density(rho_alpha, state_vars);
            fluid_alpha->Viscosity(mu_alpha, state_vars);
            
            // Rho computation
            rho = rho_alpha;
            
            // Rhof computation
            rhof = rho_alpha;
            
            // lambda computation
            l[0] = rho_alpha[0]/mu_alpha[0];
            l[2] = rho_alpha[2]/mu_alpha[0] - (rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]);
            
        }
            break;
        case 2:
        {
            
            // Getting state variables
            REAL P          = state_vars[1];
            REAL S_alpha    = state_vars[2];
            REAL S_beta     = 1.0 - S_alpha;
            
            // Compute P_alpha  and  P_beta
            fluid_beta->Pc(Pc_beta, state_vars);
            vars_alpha[1]   = P - S_alpha   * Pc_beta[1];
            vars_beta[1]    = P + S_beta    * Pc_beta[1];
            
            vars_alpha[2] = S_alpha;
            vars_beta[2] = S_beta;
            
            fluid_alpha->Density(rho_alpha, vars_alpha);
            fluid_alpha->Viscosity(mu_alpha, vars_alpha);
            fluid_alpha->Kr(kr_alpha, vars_alpha);
            
            
            fluid_beta->Density(rho_beta, vars_beta);
            fluid_beta->Viscosity(mu_beta, vars_beta);
            fluid_beta->Kr(kr_beta, vars_beta);
            
            // Rho computation
            rho[0] = rho_alpha[0] * S_alpha + rho_beta[0] * S_beta;
            rho[1] = rho_alpha[1] * S_alpha + rho_beta[1] * S_beta;
            rho[2] = rho_alpha[2] * S_alpha + rho_beta[2] * S_beta;
            rho[3] = rho_alpha[3] * S_alpha + rho_beta[3] * S_beta;
            
            // Two-Phase restriction
            kr_beta[3] *= -1.0;
            
            // Fluid Mobilities
            
            l_alpha[0] = kr_alpha[0]*(rho_alpha[0]/mu_alpha[0]);
            l_alpha[1] = 0.0;
            l_alpha[2] = kr_alpha[0]*(rho_alpha[2]/mu_alpha[0]-(rho_alpha[0]*mu_alpha[2])/(mu_alpha[0]*mu_alpha[0]));
            l_alpha[3] = kr_alpha[3]*(rho_alpha[0]/mu_alpha[0]);
            
            l_beta[0] = kr_beta[0]*(rho_beta[0]/mu_beta[0]);
            l_beta[1] = 0.0;
            l_beta[2] = kr_beta[0]*(rho_beta[2]/mu_beta[0]-(rho_beta[0]*mu_beta[2])/(mu_beta[0]*mu_beta[0]));
            l_beta[3] = kr_beta[3]*(rho_beta[0]/mu_beta[0]);
            
            // lambda computation
            l[0] = l_alpha[0] + l_beta[0];
            l[1] = l_alpha[1] + l_beta[1];
            l[2] = l_alpha[2] + l_beta[2];
            l[3] = l_alpha[3] + l_beta[3];
            
            // Fractional flows
            f_alpha[0] = l_alpha[0]/l[0];
            f_alpha[1] = l_alpha[1]/l[0]-(l_alpha[0]*l[1])/(l[0]*l[0]);
            f_alpha[2] = l_alpha[2]/l[0]-(l_alpha[0]*l[2])/(l[0]*l[0]);
            f_alpha[3] = l_alpha[3]/l[0]-(l_alpha[0]*l[3])/(l[0]*l[0]);
            
            f_beta[0] = l_beta[0]/l[0];
            f_beta[1] = l_beta[1]/l[0]-(l_beta[0]*l[1])/(l[0]*l[0]);
            f_beta[2] = l_beta[2]/l[0]-(l_beta[0]*l[2])/(l[0]*l[0]);
            f_beta[3] = l_beta[3]/l[0]-(l_beta[0]*l[3])/(l[0]*l[0]);
            
            // Rhof computation
            rhof[0] = f_alpha[0] * rho_alpha[0] +  f_beta[0] * rho_beta[0];
            rhof[1] = f_alpha[0] * rho_alpha[1] + f_alpha[1] * rho_alpha[0]  + f_beta[0] * rho_beta[1] + f_beta[1] * rho_beta[0];
            rhof[2] = f_alpha[0] * rho_alpha[2] + f_alpha[2] * rho_alpha[0]  + f_beta[0] * rho_beta[2] + f_beta[2] * rho_beta[0];
            rhof[3] = f_alpha[0] * rho_alpha[3] + f_alpha[3] * rho_alpha[0]  + f_beta[0] * rho_beta[3] + f_beta[3] * rho_beta[0];
            
        }
            break;
        case 3:
        {
            std::cout<< "3-phases system not implemented"<< std::endl;
            DebugStop();
            
        }
            break;
        default:
        {
            std::cout<< "4-phases system not implemented"<< std::endl;
            DebugStop();
        }
            break;
    }
    
    // Returning needed relationships
    
    props[0] = rho_alpha;
    props[1] = rho_beta;
    props[2] = rho_gamma;
    props[3] = f_alpha;
    props[4] = f_beta;
    props[5] = f_gamma;
    props[6] = l;
    props[7] = rho;
    props[8] = rhof;
    
    return;
    
}

void TPZAxiSymmetricDarcyFlow::Rho_alpha(TPZVec<REAL> P_alpha){
    
}

void TPZAxiSymmetricDarcyFlow::Rho_beta(TPZVec<REAL> P_beta){
    
    
    
}

void TPZAxiSymmetricDarcyFlow::Rho_gamma(TPZVec<REAL> P_gamma){
    
    
    
}

