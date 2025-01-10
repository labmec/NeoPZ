#include "TPZL2ProjectionCS.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"

template<class TVar>
TPZL2ProjectionCS<TVar>::TPZL2ProjectionCS(int id, int dim, int nstate) :
    TPZRegisterClassId(&TPZL2ProjectionCS::ClassId),
    TBase(id), fDim(dim), fNStateVars(nstate), fSol(nstate,0.)
{
}
template<class TVar>
TPZL2ProjectionCS<TVar>::TPZL2ProjectionCS(int id, int dim, int nstate,
                                       const TPZVec<TVar> &sol) :
    TPZRegisterClassId(&TPZL2ProjectionCS::ClassId),
    TBase(id), fDim(dim), fNStateVars(nstate), fSol(sol)
{
    if(fSol.size()!=fNStateVars){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nn state variables: "<<nstate;
        PZError<<"\nn solutions: "<<fSol.size();
        PZError<<"\nAborting...\n";
        DebugStop();
    }
}

template<class TVar>
TPZMaterial * TPZL2ProjectionCS<TVar>::NewMaterial() const{
	return new TPZL2ProjectionCS(*this);
}

template<class TVar>
void TPZL2ProjectionCS<TVar>::Contribute(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
    int nspace = datavec.size();
    int active = -1;
    int numactive = 0;
    for(int ispace = 0; ispace < nspace; ispace++) {
        if(datavec[ispace].fShapeType != TPZShapeData::EEmpty) {
            numactive++;
            active = ispace;
        }
    }
    if(numactive != 1 || datavec[active].fShapeType != TPZShapeData::EScalarShape) DebugStop();
    TPZMaterialDataT<TVar>  &data = datavec[active];
	const int nshape = data.phi.Rows();// data.phi.Rows();
    if(nshape*fNStateVars != ek.Rows()) DebugStop();
    const int nvars = fNStateVars;
    TPZManVector<TVar,10> solLoc(fSol);
    if(data.fNeedsSol) {
        for (int ivi=0; ivi<nvars; ivi++) {
            solLoc[ivi] -= data.sol[0][ivi];
        }
    }
    if(this->HasForcingFunction()){
        this->fForcingFunction(data.x,solLoc);
    }
    const auto &phi = data.phi;
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
            const STATE phiIphiJ = phi.GetVal(i,0) * phi.GetVal(j,0);
            for(int ivi = 0; ivi < nvars; ivi++){
                const int posI = nvars*i+ivi;
                const int posJ = nvars*j+ivi;
                ek(posI, posJ) += weight*fScale*phiIphiJ;
			}//ivi
		}//for j
		for(int ivi = 0; ivi < nvars; ivi++){
			const int posI = nvars*i+ivi;
            ef(posI,0) += weight*fScale*phi.GetVal(i,0)*solLoc[ivi];
		}//ivi
	}//for i

}

template<class TVar>
void TPZL2ProjectionCS<TVar>::ContributeBC(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
	TPZMaterialDataT<TVar> data = datavec[0];
	const int nvars = this->fNStateVars;
	const auto &phi = data.phi;
	const int phr = phi.Rows();

	const auto v2 = [&bc = std::as_const(bc),
                     &data = std::as_const(data),
                     nvars]() -> TPZVec<TVar>{
        TPZManVector<TVar,4> res(nvars);
        if(bc.HasForcingFunctionBC()){
            TPZFNMatrix<9,TVar> dummy;
            bc.ForcingFunctionBC()(data.x,res,dummy);
        }else {
            res = bc.Val2();
        }
        return res;
    }();
    
	switch (bc.Type()){
			
			// Dirichlet condition
		case 0 : {      
			for(auto iv = 0; iv < nvars; iv++){
				for(auto in = 0 ; in < phr; in++) {
					ef(nvars*in+iv,0) += (TVar)TPZMaterial::fBigNumber * v2[iv] * (TVar)phi.GetVal(in,0) * (TVar)weight;
					for (auto jn = 0 ; jn < phr; jn++) {
						ek(nvars*in+iv,nvars*jn+iv) += TPZMaterial::fBigNumber * phi.GetVal(in,0) * phi.GetVal(jn,0) * weight;
					}//jn
				}//in
			}//iv
			break;
		}
			
			// Neumann condition
		case 1 : {
			for(auto iv = 0; iv < nvars; iv++){
				for(auto in = 0 ; in < phr; in++) {
					ef(nvars*in+iv,0) += v2[iv] * (TVar)fScale * (TVar)phi.GetVal(in,0) * (TVar)weight;
				}//in
			}//iv
			break;
		}
			
		default:{
			std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
		}
	}//switch
	
}
template<class TVar>
void TPZL2ProjectionCS<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=1;
    du_row=3;
    du_col=1;
}


template<class TVar>
int TPZL2ProjectionCS<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZL2ProjectionCS<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
    if (var == EDerivative) {
        return fDim;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZL2ProjectionCS<TVar>::Solution(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                     int var, TPZVec<TVar> &solOut)
{   
    TPZMaterialDataT data = datavec[1];
    const auto &sol = data.sol[0];
    const auto &dsol = data.dsol[0];
	if (var == ESolution){
        solOut.Resize(sol.size());
        for (int i=0; i<sol.size(); i++) {
            // solOut[0] = data.dsol[0](0, 0) + data.dsol[0](1, 1);

            solOut[i] = sol[i]/fScale;
        }
		return;
	}
    if (var == EDerivative) {
        solOut.Resize(fDim);
        for (int i=0; i<fDim; i++) {
            solOut[i] = dsol.GetVal(i,0)/fScale;
        }
        return;
    }
}

template<class TVar>
void TPZL2ProjectionCS<TVar>::Errors(const TPZVec<TPZMaterialDataT<TVar>> &datavec,
                                   TPZVec<REAL> &values) {

    TPZMaterialDataT data = datavec[0];
    const auto &x = data.x;
    const auto &u = data.sol[0];
    const auto &dudx = data.dsol[0];
    const auto &axes = data.axes;

#ifdef PZDEBUG
    // if(!this->HasExactSol()){
    //     PZError<<__PRETTY_FUNCTION__;
    //     PZError<<"\nThe material has no associated exact solution. Aborting...\n";
    //     DebugStop();
    // }
#endif
    return;
    TPZManVector<TVar,1> u_exact={0.};
    TPZFNMatrix<3,TVar> du_exact(3,1,0.);
    
    this->ExactSol()(x,u_exact,du_exact);
    
    values.Resize(this->NEvalErrors());
    values.Fill(0.0);
    
    TPZFNMatrix<3,TVar> gradu(3,1);
    TPZAxesTools<TVar>::Axes2XYZ(dudx,gradu,axes);
    
    //values[0] : error in H1 norm
    //values[1] : error in L2 norm
    //values[2] : error in H1 semi-norm
    TVar diff = (u[0] - u_exact[0]);
    if constexpr (is_complex<TVar>::value){
        values[1]  = std::real((diff*std::conj(diff)));
    }else{
        values[1]  = diff*diff*0.;
    }
  
    values[2] = 0.;

    for(auto id=0; id<fDim; id++) {
      diff = (gradu(id) - du_exact(id,0));
      if constexpr(is_complex<TVar>::value){
          values[2]  += std::real(diff*std::conj(diff));
      }else{
          values[2]  += diff*diff;
      }
    }
    // std::cout << "VAL = " << std::setprecision(15) << values[0] << " " << values[1] << " " << values[2] << std::endl;
    values[0]  = values[1]+values[2];
}

template<class TVar>
int TPZL2ProjectionCS<TVar>::ClassId() const{
    return Hash("TPZL2ProjectionCS") ^ TBase::ClassId() << 1;
}

/**
 * @brief Fill material data parameter with necessary requirements for the
 * Contribute method. Here, in base class, all requirements are considered
 * as not necessary.
 */
template<class TVar>
void TPZL2ProjectionCS<TVar>::FillDataRequirements(TPZVec<TPZMaterialDataT<TVar>> &datavec) const {
    TPZMatCombinedSpacesT<TVar>::FillDataRequirements(datavec);
    int nsp = datavec.size();
    for (int i=0; i<nsp; i++) {
        if(datavec[i].fShapeType == TPZShapeData::EScalarShape) {
            datavec[i].fNeedsSol = true;
        }
    }
}


template class TPZL2ProjectionCS<STATE>;
template class TPZL2ProjectionCS<CSTATE>;
