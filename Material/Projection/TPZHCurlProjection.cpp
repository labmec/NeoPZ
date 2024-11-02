#include "TPZHCurlProjection.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"


/*
  @fran: TODO: these materials are inconsistent.
  They must be rewritten taken into account that
  HCurl basis functions are always 3d and their curl
  is a scalar in 1d/2d
**/
template<class TVar>
TPZHCurlProjection<TVar>::TPZHCurlProjection(int id, int dim) :
    TPZRegisterClassId(&TPZHCurlProjection::ClassId),
    TBase(id){
    SetDimension(dim);
}

template<class TVar>
TPZMaterial * TPZHCurlProjection<TVar>::NewMaterial() const{
	return new TPZHCurlProjection(*this);
}

template<class TVar>
void TPZHCurlProjection<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	/*
     * HCurl approximation spaces in NeoPZ are created 
     * from the combination of scalar functions with 
     * constant vector fields
     */
    const auto &phiHCurl = data.phi;
    const auto &curlPhi = data.curlphi;

    //last position of solLoc is the curl
    TPZManVector<TVar,6> solLoc(fDim + fCurlDim);
    if(!this->HasForcingFunction()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" a forcing function (solution to be projected)\n";
        PZError<<"has not been set! Aborting...\n";
        DebugStop();
    }
    this->fForcingFunction(data.x,solLoc);
    TPZManVector<TVar,3> curlSol(fCurlDim,0.);
    for(int i = 0; i < fCurlDim; i++){curlSol[i] = solLoc[fDim+i];}
    
    const auto nHCurlFunctions = phiHCurl.Rows();
    for (auto iVec = 0; iVec < nHCurlFunctions; iVec++) {
        TVar load = 0.;
        
        for(auto x = 0; x < fDim; x++)
            {load += phiHCurl.GetVal(iVec, x) * solLoc[x];}
        for(auto x = 0; x < fCurlDim; x++)
            {load += curlPhi.GetVal(x, iVec) * curlSol[x];}
        
        ef(iVec,0) += load * weight;
        
        for (auto jVec = 0; jVec < nHCurlFunctions; jVec++) {

            STATE phiIdotPhiJ = 0.;
            for(auto x = 0; x < fDim; x++)
                {phiIdotPhiJ += phiHCurl.GetVal(iVec, x) * phiHCurl.GetVal(jVec, x);}
            STATE curlIcurlJ = 0.;

            for(auto x = 0; x < fCurlDim; x++)
                {curlIcurlJ +=curlPhi.GetVal(x,iVec) * curlPhi.GetVal(x,jVec);}

            ek(iVec, jVec) += (phiIdotPhiJ + curlIcurlJ) * weight;
        }
    }
}

template<class TVar>
void TPZHCurlProjection<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
    const auto & phiHCurl = data.phi;
    const auto nHCurlFunctions = phiHCurl.Rows();
    const auto phiDim = phiHCurl.Cols();
    const auto &BIG = TPZMaterial::fBigNumber;

    const TVar v2 = [&bc = std::as_const(bc),
                     &data = std::as_const(data)](){
        if(bc.HasForcingFunctionBC()){
            TPZManVector<TVar> res(3);
            TPZFNMatrix<9,TVar> dummy(3,3,0.);
            bc.ForcingFunctionBC()(data.x,res,dummy);
            return res[0];
        }else {
            return bc.Val2()[0];
        }
    }();
    switch (bc.Type()) {
        case 0:
            for (int i = 0; i < nHCurlFunctions; i++) {
                TVar rhs = 0, stiff = 0;
                for(auto x = 0; x <phiDim; x++) rhs+= phiHCurl(i, x) * BIG * v2;
                ef(i, 0) += rhs * weight;
                for (int j = 0; j < nHCurlFunctions; j++) {
                    for(auto x = 0; x <phiDim; x++) stiff+= phiHCurl(i, x) * phiHCurl(j, x) * BIG;
                    ek(i, j) += stiff * weight;
                }
            }
            break;
        case 1:
            DebugStop();
            break;
        case 2:
            DebugStop();
            break;
    }
	
}
template<class TVar>
void TPZHCurlProjection<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=fDim;
    du_row=fCurlDim;
    du_col=1;
}


template<class TVar>
int TPZHCurlProjection<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Curl",name.c_str())) return ECurl;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZHCurlProjection<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return fDim;
  if (var == ECurl) return fCurlDim;
    
	
  return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
int TPZHCurlProjection<TVar>::NSolutionVariablesBC(int var) const{
    if(fDim ==3){
        if(var == ESolution) return 2;
        if(var == ECurl) return 1;
    }else{
        if(var == ESolution) return 1;
        //no idea on how to plot the curl in 1d bound els
    }
    
	
  return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZHCurlProjection<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                     int var, TPZVec<TVar> &solOut)
{
    const auto &sol = data.sol[0];
    const auto &curlsol = data.curlsol[0];
    switch (var) {
        case ESolution: // E
        {
            for(auto x = 0; x < fDim; x++) solOut[x] = sol[x];
        } break;
        case ECurl: // curlE
        {
            for(auto x = 0; x < fCurlDim; x++) solOut[x] = curlsol[x];
        } break;
        default:
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"\nCouldnt identify solution index. Aborting...\n";
            DebugStop();
    }
}

template<class TVar>
void TPZHCurlProjection<TVar>::SolutionBC(const TPZMaterialDataT<TVar> &data,
                                          int var, TPZVec<TVar> &solOut)
{
    const int dim = this->Dimension();
    this->SetDimension(dim-1);
    this->Solution(data,var,solOut);
    this->SetDimension(dim);
}

template<class TVar>
void TPZHCurlProjection<TVar>::Errors(const TPZMaterialDataT<TVar>&data,
                                      TPZVec<REAL> &values)
{
    const auto &x = data.x;
    const auto &u = data.sol[0];
    const auto &curlU = data.curlsol[0];
#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
#endif
    TPZManVector<TVar,3> u_exact(fDim);
    TPZFNMatrix<3,TVar> curlU_exact(fCurlDim,1);
    this->ExactSol()(x, u_exact, curlU_exact);
    values.Fill(0.0);

    // values[0] : E error using HCurl norm (values[1]+values[2])
    // values[1] : E error using L2 norm
    // values[2] : E error using HCurl semi-norm

    // values[1] : E error using L2 norm
    for (int id = 0; id < fDim; id++) {
        const TVar diffU = u[id]-u_exact[id];
        if constexpr (is_complex<TVar>::value){
            values[1] += std::norm(std::conj(diffU)*diffU);
        }else{
            values[1] += diffU*diffU;
        }
    }

    // values[2] : curlE error using L2 norm
    for (int x = 0; x < fCurlDim; x++) {
        const TVar diffCurl = curlU[x] - curlU_exact(x,0);
        if constexpr (is_complex<TVar>::value){
            values[2] += std::norm(std::conj(diffCurl)*diffCurl);
        }else{
            values[2] += diffCurl*diffCurl;
        }
    }

    // values[0] : E error using HCurl norm (values[1]+values[2])
    values[0] = values[1] + values[2];
}

template<class TVar>
int TPZHCurlProjection<TVar>::IntegrationRuleOrder(const int elPMaxOrder) const
{
    //for order k functions of order k+1
    return  TBase::IntegrationRuleOrder(elPMaxOrder+1);
}

template<class TVar>
int TPZHCurlProjection<TVar>::ClassId() const{
    return Hash("TPZHCurlProjection") ^ TBase::ClassId() << 1;
}


template class TPZHCurlProjection<STATE>;
template class TPZHCurlProjection<CSTATE>;