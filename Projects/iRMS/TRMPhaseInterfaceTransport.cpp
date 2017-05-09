//
//  TRMPhaseInterfaceTransport.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMPhaseInterfaceTransport.h"

/**
 * Empty Constructor
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport() : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>()
{
    
}

/** Creates a material object and inserts it in the vector of
 *  material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(int matid) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(matid)
{
    
}


/** Creates a material object based on the referred object and
 *  inserts it in the vector of material pointers of the mesh.
 */
TRMPhaseInterfaceTransport::TRMPhaseInterfaceTransport(const TRMPhaseInterfaceTransport &mat) : TPZMatWithMem<TRMPhaseInterfaceMemory, TPZDiscontinuousGalerkin>(mat)
{
    
}

/**
 * Destructor
 */
TRMPhaseInterfaceTransport::~TRMPhaseInterfaceTransport()
{

}

/** Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as necessary. Each derived class may optimize performance by selecting
 * only the necessary data.
 * @since April 10, 2007
 */
void TRMPhaseInterfaceTransport::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    DebugStop();
}

void TRMPhaseInterfaceTransport::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

/** print out the data associated with the material */
void TRMPhaseInterfaceTransport::Print(std::ostream &out)
{
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}


/** returns the variable index associated with the name */
int TRMPhaseInterfaceTransport::VariableIndex(const std::string &name)
{
    DebugStop();
    return 0;
}


/** returns the number of variables associated with the variable
 indexed by var.  var is obtained by calling VariableIndex */
int TRMPhaseInterfaceTransport::NSolutionVariables(int var)
{
    DebugStop();
    return 0;
}


/** Computes the divergence over the parametric space */
void TRMPhaseInterfaceTransport::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    DebugStop();
}


/** returns the solution associated with the var index based on
 * the finite element approximation */
void TRMPhaseInterfaceTransport::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout)
{
    DebugStop();
}

// Contribute Methods being used

void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{

    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            ContributeBCInterface_ab(data, datavecleft, weight, ek, ef, bc);
        }
            break;
        case 3:
        {
            ContributeBCInterface_abc(data, datavecleft, weight, ek, ef, bc);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    
}


void TRMPhaseInterfaceTransport::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    
    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            ContributeInterface_ab(data, datavecleft, datavecright, weight, ek, ef);
        }
            break;
        case 3:
        {
            ContributeInterface_abc(data, datavecleft, datavecright, weight, ek, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }

}


void TRMPhaseInterfaceTransport::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{

    switch (fSimulationData->SystemType().size()) {
        case 1:
        {
            DebugStop();
        }
            break;
        case 2:
        {
            ContributeInterface_ab(data, datavecleft, datavecright, weight, ef);
        }
            break;
        case 3:
        {
            ContributeInterface_abc(data, datavecleft, datavecright, weight, ef);
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }

}

// two phase case


void TRMPhaseInterfaceTransport::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    REAL sa = datavecleft[sb_a].sol[0][0];
    
    int nphis_a_l     = phi_ss_l.Rows();
    int firsts_a_l    = 0;
    
    int global_point_index = data.intGlobPtIndex;
    
    // Get the pressure at the integrations points
    TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  & material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond > & >(bc);
    TRMPhaseInterfaceMemory &point_memory = material_bc_mem.GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n_l();
    REAL sa_avg_n    = sa;//point_memory.sa_n_l(); @omar:: saturation is not updated at faces
    REAL un_l    = point_memory.un();
    
    TPZManVector<STATE,3> n = data.normal;
    REAL p_l                  = p_avg_n;
    REAL s_l                  = sa_avg_n;
    
    
    //  Average values p_a
    STATE p_a_l    = p_l;
    STATE s_a_l    = s_l;
    
    STATE beta = 0.0;
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1);
    
    
    REAL Value_m    = 0.0;
    REAL Value_s    = 0.0;
    if (bc.HasTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> f(2);
        TPZFMatrix<double> gradf;
        REAL time = fSimulationData->t();
        bc.TimedependentBCForcingFunction()->Execute(datavecleft[sb_a].x, time, f, gradf);
        Value_m = f[0];
        Value_s = f[1];
    }
    else{
        Value_m = bc.Val2()(0,0);
    }
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD outlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0])*phi_ss_l(is,0)*un_l;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
                }
                
                
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            // upwinding
            if (Value_m > 0) {
                beta = 1.0;
            }
            
            STATE un_N = un_l;//Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = s_a_l;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_N;
                }
                
                
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            
            // upwinding
            if (un_l < 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_l;
                
            }
            
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            // upwinding
            if (Value_m < 0) {
                beta = 1.0;
            }
            
            STATE un_N = un_l;//Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            // upwinding
            beta = 1.0;
            
            STATE un_N = 0.0;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_s;
            
            this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ss_l(is,0)*un_N;
                
                
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


void TRMPhaseInterfaceTransport::ContributeBCInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


void TRMPhaseInterfaceTransport::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    
    TPZFNMatrix<100,STATE> phi_ss_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ss_r       = datavecright[sb_a].phi;
    
    int nphis_a_l     = phi_ss_l.Rows();
    int firsts_a_l    = 0;
    
    int nphis_a_r     = phi_ss_r.Rows();
    int firsts_a_r    = firsts_a_l + nphis_a_l;
    
    TPZManVector<STATE,3> n = data.normal;
    
    REAL sa_avg_l                  = datavecleft[sb_a].sol[0][0];
    REAL sa_avg_r                  = datavecright[sb_a].sol[0][0];
    REAL p_avg_l = 0.0;
    REAL p_avg_r = 0.0;
    REAL un_l = 0.0;
    
    // Interface memory
    // Get the pressure at the integrations points
    long global_point_index = data.intGlobPtIndex;
    TRMPhaseInterfaceMemory &point_memory = GetMemory()[global_point_index];
    un_l = point_memory.un();
    p_avg_l = point_memory.p_avg_n_l();
    p_avg_r = point_memory.p_avg_n_r();
    
    //  Average values p_a
    STATE p_a_l    = p_avg_l;
    STATE s_a_l    = sa_avg_l;
    STATE p_a_r    = p_avg_r;
    STATE s_a_r    = sa_avg_r;
    
    TPZManVector<STATE, 10> fa_l,v_l(nvars+1),fa_r,v_r(nvars+1);
    TPZManVector<STATE, 10> rho_w_l,rho_w_r;
    TPZManVector<STATE, 10> rho_o_l,rho_o_r;
    TPZManVector<STATE, 10> f_l,f_r;
    TPZManVector<STATE, 10> l_l,l_r;
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    
    this->fSimulationData->PetroPhysics()->fa(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(fa_r, v_r);

    v_l[0] = point_memory.p_avg_n_l(); // last state saturations
    v_r[0] = point_memory.p_avg_n_r(); // last state saturations
    v_l[1] = point_memory.sa_l(); // last state saturations
    v_r[1] = point_memory.sa_r(); // last state saturations
    
    this->fSimulationData->PetroPhysics()->fa(f_l, v_l);
    this->fSimulationData->PetroPhysics()->fa(f_r, v_r);
    
    this->fSimulationData->PetroPhysics()->l(l_l, v_l);
    this->fSimulationData->PetroPhysics()->l(l_r, v_r);
    
    this->fSimulationData->AlphaProp()->Density(rho_w_l, v_l);
    this->fSimulationData->AlphaProp()->Density(rho_w_r, v_r);

    this->fSimulationData->BetaProp()->Density(rho_o_l, v_l);
    this->fSimulationData->BetaProp()->Density(rho_o_r, v_r);
    
    
    // Explicit Upstream Differencing for Gravity term
    // Upstream Differencing for Multiphase Flow in Reservoir Simulation
    // Yann Brenier and Jérôme Jaffré
    
    TPZFMatrix<REAL> & k_avg = point_memory.K_0();
    TPZManVector<REAL,3>  g = fSimulationData->Gravity();
    REAL n_dot_K_g = 0.0;
    for (int i = 0; i < 3; i++) {
        REAL K_g = 0.0;
        for (int j =0; j<3; j++) {
             K_g += k_avg(i,j)*g[j];
        }
        n_dot_K_g += n[i]*K_g;
    }
    
<<<<<<< Updated upstream
    n_dot_K_g *= 2.0;
=======
//    n_dot_K_g *= 0.0;
>>>>>>> Stashed changes
    
    REAL delta_rho_l =rho_w_l[0]-rho_o_l[0];
    REAL delta_rho_r =rho_w_r[0]-rho_o_r[0];

    REAL theta_1 = un_l + l_l[0]*(delta_rho_l)*(1.0-f_l[0])*n_dot_K_g;
    REAL theta_2 = un_l - l_r[0]*(delta_rho_r)*(f_r[0])*n_dot_K_g;
    
    bool take_l = false;
    bool take_r = false;

    take_l = (0.0 <= theta_1) && (theta_1 <= theta_2);
    take_r = (theta_1 <= 0.0) && (0.0 <= theta_2);
    
<<<<<<< Updated upstream
=======
    take_l = (0.0 <= theta_1) && (theta_1 <= theta_2);
    take_r = (theta_1 <= 0.0) && (0.0 <= theta_2);
    
//    if (un_l >= 0.0) {
//
//    }
//    else{
//        take_l = (0.0 >= theta_1) && (theta_1 >= theta_2);
//        take_r = (theta_1 >= 0.0) && (0.0 >= theta_2);
//    }
>>>>>>> Stashed changes
    
    REAL qgn_l = l_l[0]*(delta_rho_l)*(1.0-f_l[0])*n_dot_K_g; // last
    state saturations
    REAL qgn_r = l_r[0]*(delta_rho_r)*(1.0-f_r[0])*n_dot_K_g; // last state saturations
    
    if (take_l && take_r) {
        std::cout << "Jaffré:: both true! " << std::endl;
        DebugStop();
    }
    
//    if (!take_l && !take_r) {
//        std::cout << "Jaffré:: Just water up and down! " << std::endl;
//    }
    
//    // apply gravity segregation
//    if(take_l){
//        un_l += qgn_l;
//    }
//    
//    if(take_r){
//        un_l += qgn_r;
//    }
    
    REAL epsilon = 1.0e-20;
    REAL qgn_avg = 2.0*(qgn_l*qgn_r)/(qgn_l + qgn_r + epsilon);
    un_l += qgn_avg;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_l(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ss_r(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_r, js + firsts_a_l) += -1.0*weight * beta * fa_l[2] * phi_ss_l(js,0) * phi_ss_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fa_r[2] * phi_ss_r(js,0) * phi_ss_r(is,0)*un_l;
        }
        
    }
}


void TRMPhaseInterfaceTransport::ContributeInterface_ab(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    TPZFMatrix<STATE>  ek_fake(ef.Rows(),ef.Rows(),0.0);
    this->ContributeInterface_ab(data,datavecleft,datavecright, weight, ek_fake, ef);
    return;
}


// three phase case


void TRMPhaseInterfaceTransport::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    int sb_b    = 1;
    
    TPZFNMatrix<100,STATE> phi_ssa_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l       = datavecleft[sb_b].phi;
    REAL sa = datavecleft[sb_a].sol[0][0];
    REAL sb = datavecleft[sb_b].sol[0][0];
    
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firsts_a_l    = 0;
    int firsts_b_l    = firsts_a_l + nphis_a_l;
    
    int global_point_index = data.intGlobPtIndex;
    
    // Get the pressure at the integrations points
    TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond>  & material_bc_mem = dynamic_cast<TPZMatWithMem<TRMPhaseInterfaceMemory,TPZBndCond > & >(bc);
    TRMPhaseInterfaceMemory &point_memory = material_bc_mem.GetMemory()[global_point_index];
    REAL p_avg_n    = point_memory.p_avg_n_l();
    REAL sa_avg_n    = sa;//point_memory.sa_n_l(); @omar:: saturation is not updated at faces
    REAL sb_avg_n    = sb;//point_memory.sb_n_l(); @omar:: saturation is not updated at faces
    REAL un_l    = point_memory.un();
    
    TPZManVector<STATE,3> n = data.normal;
    REAL p_l                  = p_avg_n;
    REAL sa_l                 = sa_avg_n;
    REAL sb_l                 = sb_avg_n;
    
    
    //  Average values p_a
    STATE p_a_l    = p_l;
    STATE s_a_l    = sa_l;
    STATE s_b_l    = sb_l;
    
    STATE beta = 0.0;
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1);
    
    REAL Value_m    = 0.0;
    REAL Value_sa   = 0.0;
    REAL Value_sb   = 0.0;
    if (bc.HasTimedependentBCForcingFunction()) {
        TPZManVector<STATE,2> f(3);
        TPZFMatrix<double> gradf;
        REAL time = 0.0;
        bc.TimedependentBCForcingFunction()->Execute(datavecleft[sb_a].x, time, f, gradf);
        Value_m  = f[0];
        Value_sa = f[1];
        Value_sb = f[2];
    }
    else{
        Value_m  = bc.Val2()(0,0);
        Value_sa = bc.Val2()(1,0);
        Value_sb = bc.Val2()(2,0);
    }
    
    switch (bc.Type()) {
            
        case 0 :    // Dirichlet BC  PD outlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = s_a_l;
            v_l[2] = s_b_l;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0])*phi_ssa_l(is,0)*un_l;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_l;
                }
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0])*phi_ssb_l(is,0)*un_l;
                
                for (int js = 0; js < nphis_b_l; js++) {
                    ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_l;
                }
            }
            
        }
            break;
            
        case 1 :    // Neumann BC  QN outlet
        {
            
            // upwinding
            if (Value_m > 0) {
                beta = 1.0;
            }
            
            STATE un_N = un_l;//Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = s_a_l;
            v_l[2] = s_b_l;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
                
                for (int js = 0; js < nphis_a_l; js++) {
                    ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_N;
                }
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
                
                for (int js = 0; js < nphis_b_l; js++) {
                    ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_N;
                }
            }
            
        }
            break;
            
        case 2 :    // Dirichlet BC  PD inlet
        {
            
            // upwinding
            if (un_l > 0) {
                beta = 1.0;
            }
            
            STATE p_D = Value_m;
            
            v_l[0] = p_D;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_l;
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_l;
                
            }
            
        }
            break;
            
        case 3 :    // Neumann BC  QN inlet
        {
            
            // upwinding
            if (Value_m < 0) {
                beta = 1.0;
            }
            
            STATE un_N = un_l;//Value_m;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
            }
            
        }
            break;
            
        case 4 :    // Neumann BC  Impervious bc
        {
            
            // upwinding
            beta = 1.0;
            
            STATE un_N = 0.0;
            
            v_l[0] = p_a_l;
            v_l[1] = Value_sa;
            v_l[2] = Value_sb;
            
            this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
            this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
            
            for (int is = 0; is < nphis_a_l; is++) {
                
                ef(is + firsts_a_l) += +1.0*weight * beta*fa_l[0]*phi_ssa_l(is,0)*un_N;
                
            }
            
            for (int is = 0; is < nphis_b_l; is++) {
                
                ef(is + firsts_b_l) += +1.0*weight * beta*fb_l[0]*phi_ssb_l(is,0)*un_N;
                
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


void TRMPhaseInterfaceTransport::ContributeBCInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


void TRMPhaseInterfaceTransport::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    int sb_b    = 1;
    
    TPZFNMatrix<100,STATE> phi_ssa_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l       = datavecleft[sb_b].phi;
    
    TPZFNMatrix<100,STATE> phi_ssa_r       = datavecright[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_r       = datavecright[sb_b].phi;
    
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firsts_a_l    = 0;
    int firsts_b_l    = firsts_a_l + nphis_a_l;
    
    int nphis_a_r     = phi_ssa_r.Rows();
    int nphis_b_r     = phi_ssb_r.Rows();
    int firsts_a_r    = firsts_b_l + nphis_b_l;
    int firsts_b_r    = firsts_a_r + nphis_a_r;
    
    TPZManVector<STATE,3> n = data.normal;
    
    REAL sa_avg_l                  = datavecleft[sb_a].sol[0][0];
    REAL sb_avg_l                  = datavecleft[sb_b].sol[0][0];
    
    REAL sa_avg_r                  = datavecright[sb_a].sol[0][0];
    REAL sb_avg_r                  = datavecright[sb_b].sol[0][0];
    
    REAL p_avg_l = 0.0;
    REAL p_avg_r = 0.0;
    REAL un_l = 0.0;
    
    // Interface memory
    // Get the pressure at the integrations points
    long global_point_index = data.intGlobPtIndex;
    TRMPhaseInterfaceMemory &point_memory = GetMemory()[global_point_index];
    un_l = point_memory.un();
    p_avg_l = point_memory.p_avg_n_l();
    p_avg_r = point_memory.p_avg_n_r();
    //    sa_avg_l = point_memory.sa_n_l(); @omar:: saturation is not updated
    //    sa_avg_r = point_memory.sa_n_r();
    
    //  Average values p_a
    STATE p_a_l    = p_avg_l;
    STATE s_a_l    = sa_avg_l;
    STATE s_b_l    = sb_avg_l;
    STATE p_a_r    = p_avg_r;
    STATE s_a_r    = sa_avg_r;
    STATE s_b_r    = sb_avg_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1),fa_r,fb_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_l[2] = s_b_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    v_r[2] = s_b_r;
    
    this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa_3p(fa_r, v_r);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_r, v_r);
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_l(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_l, js + firsts_a_l) += +1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_l, js + firsts_a_r) += +1.0*weight * (1.0-beta) * fa_r[2] * phi_ssa_r(js,0) * phi_ssa_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_r(is,0)*un_l;
        
        for (int js = 0; js < nphis_a_l; js++) {
            ek(is + firsts_a_r, js + firsts_a_l) += -1.0*weight * beta * fa_l[2] * phi_ssa_l(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_a_r; js++) {
            ek(is + firsts_a_r, js + firsts_a_r) += -1.0*weight * (1.0-beta) * fa_r[2] * phi_ssa_r(js,0) * phi_ssa_r(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_b_l; is++) {
        
        ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_l(is,0)*un_l;
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_b_l, js + firsts_b_l) += +1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_r; js++) {
            ek(is + firsts_b_l, js + firsts_b_r) += +1.0*weight * (1.0-beta) * fb_r[3] * phi_ssb_r(js,0) * phi_ssb_l(is,0)*un_l;
        }
        
    }
    
    for (int is = 0; is < nphis_b_r; is++) {
        
        ef(is + firsts_b_r) += -1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_r(is,0)*un_l;
        
        for (int js = 0; js < nphis_b_l; js++) {
            ek(is + firsts_b_r, js + firsts_b_l) += -1.0*weight * beta * fb_l[3] * phi_ssb_l(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
        for (int js = 0; js < nphis_b_r; js++) {
            ek(is + firsts_b_r, js + firsts_b_r) += -1.0*weight * (1.0-beta) * fb_r[3] * phi_ssb_r(js,0) * phi_ssb_r(is,0)*un_l;
        }
        
    }
    
}


void TRMPhaseInterfaceTransport::ContributeInterface_abc(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef)
{
    if (!fSimulationData->IsCurrentStateQ()) {
        return;
    }
    
    int nvars = 4; // {p,sa,sb,t}
    int sb_a    = 0;
    int sb_b    = 1;
    
    TPZFNMatrix<100,STATE> phi_ssa_l       = datavecleft[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_l       = datavecleft[sb_b].phi;
    
    TPZFNMatrix<100,STATE> phi_ssa_r       = datavecright[sb_a].phi;
    TPZFNMatrix<100,STATE> phi_ssb_r       = datavecright[sb_b].phi;
    
    int nphis_a_l     = phi_ssa_l.Rows();
    int nphis_b_l     = phi_ssb_l.Rows();
    int firsts_a_l    = 0;
    int firsts_b_l    = firsts_a_l + nphis_a_l;
    
    int nphis_a_r     = phi_ssa_r.Rows();
    int nphis_b_r     = phi_ssb_r.Rows();
    int firsts_a_r    = firsts_b_l + nphis_b_l;
    int firsts_b_r    = firsts_a_r + nphis_a_r;
    
    TPZManVector<STATE,3> n = data.normal;
    
    REAL sa_avg_l                  = datavecleft[sb_a].sol[0][0];
    REAL sb_avg_l                  = datavecleft[sb_b].sol[0][0];
    
    REAL sa_avg_r                  = datavecright[sb_a].sol[0][0];
    REAL sb_avg_r                  = datavecright[sb_b].sol[0][0];
    
    REAL p_avg_l = 0.0;
    REAL p_avg_r = 0.0;
    REAL un_l = 0.0;
    
    // Interface memory
    // Get the pressure at the integrations points
    long global_point_index = data.intGlobPtIndex;
    TRMPhaseInterfaceMemory &point_memory = GetMemory()[global_point_index];
    un_l = point_memory.un();
    p_avg_l = point_memory.p_avg_n_l();
    p_avg_r = point_memory.p_avg_n_r();
    //    sa_avg_l = point_memory.sa_n_l(); @omar:: saturation is not updated
    //    sa_avg_r = point_memory.sa_n_r();
    
    //  Average values p_a
    STATE p_a_l    = p_avg_l;
    STATE s_a_l    = sa_avg_l;
    STATE s_b_l    = sb_avg_l;
    STATE p_a_r    = p_avg_r;
    STATE s_a_r    = sa_avg_r;
    STATE s_b_r    = sb_avg_r;
    
    STATE beta = 0.0;
    // upwinding
    if (un_l > 0.0) {
        beta = 1.0;
    }
    
    TPZManVector<STATE, 10> fa_l,fb_l,v_l(nvars+1),fa_r,fb_r,v_r(nvars+1);
    v_l[0] = p_a_l;
    v_l[1] = s_a_l;
    v_l[2] = s_b_l;
    v_r[0] = p_a_r;
    v_r[1] = s_a_r;
    v_r[2] = s_b_r;
    
    this->fSimulationData->PetroPhysics()->fa_3p(fa_l, v_l);
    this->fSimulationData->PetroPhysics()->fa_3p(fa_r, v_r);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_l, v_l);
    this->fSimulationData->PetroPhysics()->fb_3p(fb_r, v_r);
    
    for (int is = 0; is < nphis_a_l; is++) {
        
        ef(is + firsts_a_l) += +1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_l(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_a_r; is++) {
        
        ef(is + firsts_a_r) += -1.0*weight * (beta*fa_l[0] + (1.0-beta)*fa_r[0])*phi_ssa_r(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_b_l; is++) {
        
        ef(is + firsts_b_l) += +1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_l(is,0)*un_l;
        
    }
    
    for (int is = 0; is < nphis_b_r; is++) {
        
        ef(is + firsts_b_r) += -1.0*weight * (beta*fb_l[0] + (1.0-beta)*fb_r[0])*phi_ssb_r(is,0)*un_l;
        
    }
    
}


/**
 * Unique identifier for serialization purposes
 */
int TRMPhaseInterfaceTransport::ClassId() const
{
    DebugStop();
}


/**
 * Save the element data to a stream
 */
void TRMPhaseInterfaceTransport::Write(TPZStream &buf, int withclassid)
{
    DebugStop();
}


/**
 * Read the element data from a stream
 */
void TRMPhaseInterfaceTransport::Read(TPZStream &buf, void *context)
{
    DebugStop();
}



/// Copy the n+1 data to the n data
void TRMPhaseInterfaceTransport::UpdateMemory()
{
    DebugStop();
}

