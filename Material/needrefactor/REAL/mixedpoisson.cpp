/**
 * @file
 * @brief Contains the methods of the TPZMixedPoisson class (multiphysics environment)
 * @author Agnaldo Farias
 * @date 2012/05/28
 */

#include "mixedpoisson.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef PZ_LOG
static TPZLogger logdata("pz.mixedpoisson.data");
static TPZLogger logerror("pz.mixedpoisson.error");
#endif

TPZMixedPoisson::TPZMixedPoisson(): TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMatPoisson3d() {
    fvisc = 1.;
	ff = 0.;
    fIsStabilized = false;
    fdelta1 = 0.;
    fdelta2 = 0.;
    fUseHdois = false;
    fh2 = 1.;

    fPermeabilityFunction = NULL;
}

TPZMixedPoisson::TPZMixedPoisson(int matid, int dim): TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMatPoisson3d(matid,dim) {
    if (dim < 1) {
        DebugStop();
    }
    fvisc = 1.;
	ff = 0.;
    fIsStabilized = false;
    fdelta1 = 0.;
    fdelta2 = 0.;
    fUseHdois = false;
    fh2 = 1.;

    fPermeabilityFunction = NULL;
}

TPZMixedPoisson::~TPZMixedPoisson() {
}

TPZMixedPoisson::TPZMixedPoisson(const TPZMixedPoisson &cp) :TPZRegisterClassId(&TPZMixedPoisson::ClassId), TPZMatPoisson3d(cp) {
    fvisc = cp.fvisc;
    ff = cp.ff;
    fIsStabilized = cp.fIsStabilized;
    fdelta1 = cp.fdelta1;
    fdelta2 = cp.fdelta2;
    fUseHdois = cp.fUseHdois;
    fh2 = cp.fh2;

    fPermeabilityFunction = cp.fPermeabilityFunction;
}

TPZMixedPoisson & TPZMixedPoisson::operator=(const TPZMixedPoisson &copy){
    TPZMatPoisson3d::operator=(copy);
    fvisc = copy.fvisc;
    ff = copy.ff;
    fIsStabilized = copy.fIsStabilized;
    fdelta1 = copy.fdelta1;
    fdelta2 = copy.fdelta2;
    fUseHdois = copy.fUseHdois;
    fh2 = copy.fh2;

    fPermeabilityFunction = copy.fPermeabilityFunction;
    return *this;
} 

int TPZMixedPoisson::NStateVariables() const {
	return 1;
}

void TPZMixedPoisson::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Base Class properties :";
	TPZMatPoisson3d::Print(out);
	out << "\n";
}

void TPZMixedPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {

//#ifdef PZDEBUG
//    int nref = datavec.size();
//
//    if (nref != 2) {
//        std::cout << " Error. The size of the datavec is different from 2." << std::endl;
//        DebugStop();
//    }
//#endif
    
    
    
    STATE force = ff;
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		force = res[0];
	}
    
    TPZFNMatrix<9,STATE> PermTensor;
    TPZFNMatrix<9,STATE> InvPermTensor;
    
    GetPermeabilities(datavec[1].x, PermTensor, InvPermTensor);
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9,REAL> dphiPXY(3,dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    
    
    REAL &faceSize = datavec[0].HSize;
    if(fUseHdois==true){
        fh2 = faceSize*faceSize;
    }else fh2 = 1.;
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
	
    int nactive = 0;
    for (int i=0; i<datavec.size(); i++) {
        if (datavec[i].fActiveApproxSpace == ETestTrial) {
            nactive++;
        } else if(datavec[i].fActiveApproxSpace != ENotActive) {
            DebugStop();
        }
    }
#ifdef PZDEBUG
    if(nactive == 4)
    {
        int phrgb = datavec[2].phi.Rows();
        int phrub = datavec[3].phi.Rows();
        if(phrp+phrq+phrgb+phrub != ek.Rows())
        {
            DebugStop();
        }
    }else
    {
        if(phrp+phrq != ek.Rows())
        {
            DebugStop();
        }
    }
#endif
	//Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }
        
        //Inserindo termo de estabilizacao no termo de fonte
        REAL divqi = 0.;
        if(fIsStabilized)
        {
            //calculando div(qi)
            TPZFNMatrix<3,REAL> axesvec(3,1,0.);
            datavec[0].axes.Multiply(ivec,axesvec);
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divqi += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            ef(iq, 0) += weight*(fdelta2*divqi*force);
        }
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
            
            
            //Inserindo termos de estabilizacao na matriz do fluxo
            if(fIsStabilized==true)
            {
                //termos de delta1
                //dot product between uKinv[v]
                ivecZ.Zero();
                for(int id=0; id<3; id++){
                    for(int jd=0; jd<3; jd++){
                        ivecZ(id,0) += InvPermTensor(id,jd)*ivec(jd,0);
                    }
                }
                //ivecZ.Print("mat2 = ");
                REAL prod2 = ivecZ(0,0)*jvec(0,0) + ivecZ(1,0)*jvec(1,0) + ivecZ(2,0)*jvec(2,0);
                ek(iq,jq) += (-1.)*weight*fh2*fdelta1*fvisc*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod2;
                
                
                //termos de delta2: //dot product between divQu.divQv
                REAL divqj = 0.;
                TPZFNMatrix<3,REAL> axesvec(3,1,0.);
                datavec[0].axes.Multiply(jvec,axesvec);
                //calculando div(qj)
                for(int jloc=0; jloc<fDim; jloc++)
                {
                    divqj += axesvec(jloc,0)*dphiQ(jloc,jshapeind);
                }
                ek(iq,jq) += weight*fdelta2*divqi*divqj;
            }
            
        }
    }
    
	
	// Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
            //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
            //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
        datavec[0].axes.Multiply(ivec,axesvec);
        
        REAL divwq = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
            
            
            //Inserindo termo de estabilizacao: delta1
            if(fIsStabilized==true)
            {
                //produto gardPu.Qv
                REAL dotVGradP = 0.;
    
                for(int k =0; k<fDim; k++)
                {
                    dotVGradP += ivec(k,0)*phiQ(ishapeind,0)*dphiP(k,jp);
                }
                
                REAL integration = (-1.)*weight*fh2*fdelta1*dotVGradP;
                
                // Estabilizacao delta1 na Matrix B
                ek(iq, phrq+jp) += integration;
                
                // Estabilizacao delta1 na Matrix BˆT
                ek(phrq+jp,iq) += integration;
            }
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
    //Contribution for estabilization delta1 for gradPu*gradPv. Matrix D
    if(fIsStabilized==true)
    {
        //produto KgradPu x KgradPv

        TPZFNMatrix<3,STATE> dphiPuZ(dphiP.Rows(),dphiP.Cols(),0.), dphiPXYState(3,1);
        for(int i=0; i<3; i++) dphiPXYState(i,0) = dphiPXY(i,0);
        PermTensor.Multiply(dphiPXYState, dphiPuZ);
        
        REAL umSobreVisc = 1./fvisc;
        for(int ip=0; ip<phrp; ip++)
        {
            for(int jp=0; jp<phrp; jp++)
            {
                for(int k =0; k<3; k++)
                {
                    ek(phrq+ip, phrq+jp) += (-1.)*weight*fh2*fdelta1*umSobreVisc*dphiPuZ(k,ip)*dphiPXY(k,jp);
                }
                
            }
        }
    }
    if(nactive == 4)
    {
        for(int ip=0; ip<phrp; ip++)
        {
            ek(phrq+ip,phrq+phrp) += phip(ip,0)*weight;
            ek(phrq+phrp,phrq+ip) += phip(ip,0)*weight;
        }
        ek(phrp+phrq+1,phrq+phrp) += -weight;
        ek(phrq+phrp,phrp+phrq+1) += -weight;
    }
    //
//    #ifdef PZ_LOG
//        if(logdata.isDebugEnabled())
//        {
//            std::stringstream sout;
//            sout<<"\n\n Matriz ek e vetor fk \n ";
//            ek.Print("ekmph = ",sout,EMathematicaInput);
//            ef.Print("efmph = ",sout,EMathematicaInput);
//            LOGPZ_DEBUG(logdata,sout.str());
//        }
//    #endif
    
}


//void TPZMixedPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
//    
//#ifdef PZDEBUG
//	int nref =  datavec.size();
//	if (nref != 2 ) {
//        std::cout << " Erro. The size of the datavec is different from 2 \n";
//		DebugStop();
//	}
//#endif
//    
//    STATE force = ff;
//    if(fForcingFunction) {            
//		TPZManVector<STATE> res(1);
//		fForcingFunction->Execute(datavec[1].x,res);   
//		force = res[0];
//	}
//    
//    // Setting the phis
//    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
//    TPZFMatrix<REAL>  &phip =  datavec[1].phi;
//	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
//    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
//    
//    REAL &faceSize = datavec[0].HSize;
//    if(fUseHdois==true){
//        fh2 = faceSize*faceSize;
//    }else fh2 = 1.;
//    
//    int phrq, phrp;
//    phrp = phip.Rows();
//    phrq = datavec[0].fVecShapeIndex.NElements();
//	
//	//Calculate the matrix contribution for flux. Matrix A
//    REAL InvK = fvisc/fk;
//    for(int iq=0; iq<phrq; iq++)
//    {
//        //ef(iq, 0) += 0.;
//        int ivecind = datavec[0].fVecShapeIndex[iq].first;
//        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
//        TPZFNMatrix<3> ivec(3,1);
//        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
//        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
//        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
//        
//        //Inserindo termo de estabilizacao no termo de fonte
//        REAL divqi = 0.;
//        if(fIsStabilized==true)
//        {
//            //calculando div(qi)
//            TPZFNMatrix<3> axesvec(3,1);
//            datavec[0].axes.Multiply(ivec,axesvec);
//            for(int iloc=0; iloc<fDim; iloc++)
//            {
//                divqi += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
//            }
//            ef(iq, 0) += weight*(fdelta2*divqi*force);
//        }
//        
//        for (int jq=0; jq<phrq; jq++)
//        {
//            TPZFNMatrix<3> jvec(3,1);
//            int jvecind = datavec[0].fVecShapeIndex[jq].first;
//            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
//            jvec(0,0) = datavec[0].fDeformedDirections(0,jvecind);
//            jvec(1,0) = datavec[0].fDeformedDirections(1,jvecind);
//            jvec(2,0) = datavec[0].fDeformedDirections(2,jvecind);
//            
//            //dot product between u and v
//            REAL prod = ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);          
//            ek(iq,jq) += InvK*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
//            
//            
//            //Inserindo termos de estabilizacao na matriz do fluxo
//            if(fIsStabilized==true)
//            {
//                //termos de delta1
//                ek(iq,jq) += (-1.)*weight*fh2*fdelta1*InvK*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
//                
//                
//                //termos de delta2
//                REAL divqj = 0.;
//                TPZFNMatrix<3> axesvec(3,1);
//                datavec[0].axes.Multiply(jvec,axesvec);
//                //calculando div(qj)
//                for(int jloc=0; jloc<fDim; jloc++)
//                {
//                    divqj += axesvec(jloc,0)*dphiQ(jloc,jshapeind);
//                }
//                ek(iq,jq) += weight*fdelta2*divqi*divqj;
//            }
//
//        }
//    }
//    
//	
//	// Coupling terms between flux and pressure. Matrix B
//    for(int iq=0; iq<phrq; iq++)
//    {
//        int ivecind = datavec[0].fVecShapeIndex[iq].first;
//        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
//        
//        TPZFNMatrix<3> ivec(3,1);
//        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
//        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
//        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
//        TPZFNMatrix<3> axesvec(3,1);
//        datavec[0].axes.Multiply(ivec,axesvec);
//        
//        REAL divwq = 0.;
//        for(int iloc=0; iloc<fDim; iloc++)
//        {
//            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
//        }
//        for (int jp=0; jp<phrp; jp++) {
//            
//            REAL fact = (-1.)*weight*phip(jp,0)*divwq;            
//            // Matrix B
//            ek(iq, phrq+jp) += fact;
//            
//            // Matrix B^T
//            ek(phrq+jp,iq) += fact;
//            
//            
//            //Inserindo termo de estabilizacao: delta1
//            if(fIsStabilized==true)
//            {
//                REAL dotVGradP = 0.;
//                for(int k =0; k<fDim; k++)
//                {
//                    dotVGradP += ivec(k,0)*phiQ(ishapeind,0)*dphiP(k,jp);
//                }
//                
//                REAL integration = (-1.)*weight*fh2*fdelta1*dotVGradP;
//                
//                // Estabilizacao delta1 na Matrix B
//                ek(iq, phrq+jp) += integration;
//                
//                // Estabilizacao delta1 na Matrix BˆT
//                ek(phrq+jp,iq) += integration;
//            }
//        }
//    }
//    
//    //termo fonte referente a equacao da pressao
//    for(int ip=0; ip<phrp; ip++){
//         ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
//    }
//    
//    //Contribution for estabilization delta1 for gradPu*gradPv. Matrix D
//    if(fIsStabilized==true)
//    {
//        REAL mK = fk/fvisc;
//        for(int ip=0; ip<phrp; ip++)
//        {
//             for(int jp=0; jp<phrp; jp++)
//             {
//                 for(int k =0; k<fDim; k++)
//                 {
//                     ek(phrq+ip, phrq+jp) += (-1.)*weight*fh2*fdelta1*mK*dphiP(k,ip)*dphiP(k,jp);
//                 }
//                 
//             }
//        }
//    }
////
////#ifdef PZ_LOG
////    if(logdata.isDebugEnabled())
////	{
////        std::stringstream sout;
////        sout<<"\n\n Matriz ek e vetor fk \n ";
////        ek.Print("ekmph = ",sout,EMathematicaInput);
////        ef.Print("efmph = ",sout,EMathematicaInput);
////        LOGPZ_DEBUG(logdata,sout.str());
////	}
////#endif
//    
//}

void TPZMixedPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
    int nref =  datavec.size();
//    if (nref != 2 ) {
//        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
//        DebugStop();
//    }
//	if (bc.Type() > 2 ) {
//        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
//		DebugStop();
//	}
#endif
    
    
    int dim = Dimension();
	
	TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
	int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()(0,0);
    REAL v1 = bc.Val1()(0,0);
    REAL u_D = 0;
    REAL normflux = 0.;
    
    if(bc.HasForcingFunction())
    {
		TPZManVector<STATE> res(3);
        TPZFNMatrix<9,STATE> gradu(dim,1);
        bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
        TPZFNMatrix<9,STATE> PermTensor, InvPermTensor;
        GetPermeabilities(datavec[0].x, PermTensor, InvPermTensor);
        
        
        for(int i=0; i<3; i++)
        {
            for(int j=0; j<dim; j++)
            {
                normflux += datavec[0].normal[i]*PermTensor(i,j)*gradu(j,0);
            }
        }
        
        
        if(bc.Type() == 0||bc.Type() == 4)
        {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        }
        else if(bc.Type() == 1 || bc.Type() == 2)
        {
            v2 = -normflux;
            if(bc.Type() ==2)
            {
                v2 = -res[0]+v2/v1;
            }
        }
        else
        {
            DebugStop();
        }
	}else
    {
        v2 = bc.Val2()(0,0);
    }

	switch (bc.Type()) {
		case 0 :		// Dirichlet condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
			
		case 1 :			// Neumann condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                }
            }
			break;
        
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
				ef(iq,0) += v2*phiQ(iq,0)*weight;
				for (int jq = 0; jq < phrq; jq++) {
					ek(iq,jq) += weight/v1*phiQ(iq,0)*phiQ(jq,0);
				}
			}
            break;
            
        case 4:
            //this case implemented the general Robin boundary condition
            // sigma.n = Km(u-u_D)+g
            //val1(0,0) = Km
            //val2(1,0) = g
            if(IsZero(bc.Val1()(0,0))){
                
                for(int iq=0; iq<phrq; iq++)
                {
                    ef(iq,0)+= gBigNumber*normflux*phiQ(iq,0)*weight;
                    for (int jq=0; jq<phrq; jq++) {
                        
                        ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight;
                    }
                }
                
            }
            else{
            
            REAL InvKm = 1./bc.Val1()(0,0);
            REAL g = normflux;
            for(int in = 0 ; in < phiQ.Rows(); in++) {
                //<(InvKm g - u_D)*(v.n)
                ef(in, 0) +=  (STATE)(InvKm*g -u_D)*phiQ(in,0)* weight;
                for (int jn = 0 ; jn < phiQ.Rows(); jn++) {
                    //InvKm(sigma.n)(v.n)
                    ek(in,jn) += (STATE)(InvKm*phiQ(in,0) * phiQ(jn,0) * weight);
                }
            }
            }
        
            break;
        
	}
    
}

/** Returns the variable index associated with the name */
int TPZMixedPoisson::VariableIndex(const std::string &name){
	if(!strcmp("Flux",name.c_str()))        return  31;
	if(!strcmp("Pressure",name.c_str()))    return  32;
    if(!strcmp("GradFluxX",name.c_str()))   return  33;
    if(!strcmp("GradFluxY",name.c_str()))   return  34;
    if(!strcmp("DivFlux",name.c_str()))   return  35;
    
    if(!strcmp("ExactPressure",name.c_str()))  return 36;
    if(!strcmp("ExactFlux",name.c_str()))  return 37;
    
    if(!strcmp("POrder",name.c_str()))        return  38;
    if(!strcmp("GradPressure",name.c_str()))        return  39;
    if(!strcmp("Divergence",name.c_str()))      return  40;
    if(!strcmp("ExactDiv",name.c_str()))        return  41;
    if (!strcmp("Derivative",name.c_str())) {
        return 42;
    }
    if (!strcmp("Permeability",name.c_str())) {
        return 43;
    }
    if(!strcmp("g_average",name.c_str()))        return  44;
    if(!strcmp("u_average",name.c_str()))        return  45;

    if(!strcmp("ExactFluxShiftedOrigin",name.c_str()))  return  46;
    return TPZMatPoisson3d::VariableIndex(name);
    
}

int TPZMixedPoisson::NSolutionVariables(int var){
	if(var == 31) return 3;
	if(var == 32 || var==8) return 1;
    if(var == 33) return 3;
    if(var == 34) return 3;
    if(var == 35) return 1;
    if(var == 36) return 1;
    if(var == 37) return fDim;
    if(var == 38) return 1;
    if(var == 39) return fDim;
    if(var == 40 || var == 41) return 1;
    if(var == 42) return 3;
    if(var == 43) return 1;
    if(var == 44) return 1;
    if(var == 45) return 1;
    if(var == 46) return fDim;
    return TPZMatPoisson3d::NSolutionVariables(var);
}

void TPZMixedPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> SolP, SolQ;
    
    TPZFNMatrix<9,STATE> PermTensor; this->GetPermeability(PermTensor);
    TPZFNMatrix<9,STATE> InvPermTensor; this ->GetInvPermeability(InvPermTensor);

    //int rtens = 2*fDim;
    if(fPermeabilityFunction){
        PermTensor.Redim(fDim,fDim);
        InvPermTensor.Redim(fDim,fDim);
        TPZFNMatrix<18,STATE> resultMat(2*fDim,fDim);
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        
        for(int id=0; id<fDim; id++){
            for(int jd=0; jd<fDim; jd++){
                
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+fDim,jd);
            }
        }
    }

   // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 31){ //function (state variable Q)
        for (int i=0; i<3; i++)
        {
            Solout[i] = datavec[0].sol[0][i];
            
        }
		return;
	}
    
    if(var == 32){
		Solout[0] = SolP[0];//function (state variable p)
		return;
	}
    
    if(var==33){
        Solout[0]=datavec[0].dsol[0](0,0);
        Solout[1]=datavec[0].dsol[0](1,0);
        Solout[2]=datavec[0].dsol[0](2,0);
        return;
    }

    if(var==34){
        Solout[0]=datavec[0].dsol[0](0,1);
        Solout[1]=datavec[0].dsol[0](1,1);
        Solout[2]=datavec[0].dsol[0](2,1);
        return;
    }
    
    if(var==35){
        Solout[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        return;
    }

    // Exact solution
    if (var == 36) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        if (fExactSol) {
            fExactSol->Execute(datavec[0].x, exactSol, flux);
        }
        Solout[0] = exactSol[0];
        return;
    } // var6

    if (var == 37) {

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            fExactSol->Execute(datavec[0].x, exactSol, gradu);
        }

        TPZFNMatrix<3, REAL> flux(3, 1);

        PermTensor.Multiply(gradu, flux);

        for (int i = 0; i < fDim; i++) {
            Solout[i] = -flux(i, 0);
        }

        return;
    } // var7

    if(var==38){
        Solout[0] = datavec[1].p;
        return;
    }
    
    if(var==39){

        TPZFMatrix<REAL> dsoldx(fDim,1.);
        TPZFMatrix<REAL> dsoldaxes(fDim,1);
        
        dsoldaxes = datavec[1].dsol[0];

        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        
        for (int i=0; i<fDim; i++)
        {
            Solout[i] = dsoldx(i,0);
        }
        
        
        return;
    }
    
    if(var==40){
        Solout[0]= 0.;
       // Solout[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        for (int j=0; j<fDim; j++) {
            Solout[0] += datavec[0].dsol[0](j,j);
        }
        return;
    }
    
    if(var==41){
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        fExactSol->Execute(datavec[0].x, exactSol, flux);
        Solout[0]=flux(2,0);
        return;
    }
    if(var == 42)
    {
        for(int i=0; i<fDim; i++)
        {
            Solout[i] = 0.;
        }
        for (int i=0; i<fDim; i++) {
            for (int j=0; j<fDim; j++) {
                Solout[i] -= InvPermTensor(i,j)*datavec[0].sol[0][i];
            }
        }
        return;
    }
    if(var ==43)
    {
        Solout[0] = PermTensor(0,0);
        return;
    }
    
    if(datavec.size() == 4)
    {
        if(var ==44)
        {
            Solout[0] = datavec[2].sol[0][0];
            return;
        }
        if(var ==45)
        {
            Solout[0] = datavec[3].sol[0][0];
            return;
        }
        
    }

    if(var == 46) { //ExactFluxShiftedOrigin
        // Solution EArcTan returns NAN for (x,y) == (0,0). Replacing data.x by
        // inf solves this problem.
        STATE infinitesimal = 0.0000000001;
        TPZManVector<REAL, 3> inf = {infinitesimal, infinitesimal, infinitesimal};

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            if (datavec[0].x[0] == 0. && datavec[0].x[1] == 0.) {
                fExactSol->Execute(inf, exactSol, gradu);
            } else {
                fExactSol->Execute(datavec[0].x, exactSol, gradu);
            }
        }
        TPZFNMatrix<3, REAL> flux(3, 1);

        PermTensor.Multiply(gradu, flux);

        for (int i = 0; i < fDim; i++) {
            Solout[i] = flux(i, 0);
        }

        return;
    }

    TPZMatPoisson3d::Solution(datavec, var, Solout);
}


void TPZMixedPoisson::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(false);
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
	}
}

void TPZMixedPoisson::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
/**
datavec[0]= Flux
datavec[1]= Pressure
 
 Errors
 [0] L2 for pressure
 [1] L2 for flux
 [2] L2 for div(flux)
 [3] Grad pressure (Semi H1)
 [4] Hdiv norm
**/

errors.Resize(NEvalErrors());
errors.Fill(0.0);

int dim = fDim;

TPZManVector<STATE,3> fluxfem(3),pressurefem(1);
fluxfem = data[0].sol[0];
STATE divsigmafem=data[0].divsol[0][0];
    
TPZGradSolVec &dsol=data[1].dsol;


TPZVec<STATE> divsigma(1);

if(this->fExactSol) {
    this->fExactSol->Execute(data[0].x, u_exact, du_exact);
}
if(this->fForcingFunction){
   this->fForcingFunction->Execute(data[0].x,divsigma);
}

REAL residual = 0.;
residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
pressurefem[0] = data[1].sol[0][0];

TPZFNMatrix<9,STATE> PermTensor; this->GetPermeability(PermTensor);
TPZFNMatrix<9,STATE> InvPermTensor; this->GetInvPermeability(InvPermTensor);
    
if(fPermeabilityFunction){
   PermTensor.Redim(3,3);
   InvPermTensor.Redim(3,3);
   TPZFNMatrix<3,STATE> resultMat;
   TPZManVector<STATE> res;
   fPermeabilityFunction->Execute(data[0].x,res,resultMat);
   for(int id=0; id<dim; id++){
       for(int jd=0; jd<dim; jd++){
           PermTensor(id,jd) = resultMat(id,jd);
           InvPermTensor(id,jd) = resultMat(id+dim,jd);
       }
   }
}

    TPZManVector<STATE,3> gradpressurefem(fDim,0.);
    this->Solution(data,VariableIndex("GradPressure"), gradpressurefem);

TPZFNMatrix<3,STATE> fluxexactneg;

//sigma=-K grad(u)

   TPZFNMatrix<9,STATE> gradpressure(3,1);
   for (int i=0; i<3; i++) {
       gradpressure(i,0) = du_exact[i];
   }
   PermTensor.Multiply(gradpressure,fluxexactneg);



REAL L2flux = 0., L2grad = 0.;
for (int i=0; i<fDim; i++) {
 //   std::cout<<"fluxo fem = "<<fluxfem[i]<<std::endl;
   for (int j=0; j<fDim; j++) {
       L2flux += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -
       
   }
    
   // std::cout<<"grad sol = "<<gradpressurefem[i]<<std::endl;
    L2grad += (gradpressure(i,0)-gradpressurefem[i])*(gradpressure(i,0)-gradpressurefem[i]);
}
errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//L2 error for pressure
errors[1] = L2flux;//L2 error for flux
errors[2] = residual;//L2 for div
errors[3] = L2grad;
errors[4] = L2flux + residual;
                             

}

void TPZMixedPoisson::Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors)
{
/**
datavec[0]= Flux
datavec[1]= Pressure

 Errors
 [0] L2 for pressure
 [1] L2 for flux
 [2] L2 for div(flux)
 [3] Grad pressure (Semi H1)
 [4] Hdiv norm
**/

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    int dim = fDim;

    TPZManVector<STATE,3> fluxfem(3),pressurefem(1);
    fluxfem = data[0].sol[0];
    STATE divsigmafem=data[0].divsol[0][0];

    TPZGradSolVec &dsol=data[1].dsol;


    TPZVec<STATE> divsigma(1);

    TPZVec<STATE> u_exact(1,0);
    TPZFMatrix<STATE> du_exact(3,1,0);
    if(this->fExactSol) {
        this->fExactSol->Execute(data[0].x, u_exact, du_exact);
    }
    if(this->fForcingFunction){
        this->fForcingFunction->Execute(data[0].x,divsigma);
    }

    REAL residual = 0.;
    residual = (divsigma[0] - divsigmafem)*(divsigma[0] - divsigmafem);
    pressurefem[0] = data[1].sol[0][0];

    TPZFNMatrix<9,STATE> PermTensor; this->GetPermeability(PermTensor);
    TPZFNMatrix<9,STATE> InvPermTensor; this->GetInvPermeability(InvPermTensor);

    if(fPermeabilityFunction){
        PermTensor.Redim(3,3);
        InvPermTensor.Redim(3,3);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(data[0].x,res,resultMat);
        for(int id=0; id<dim; id++){
            for(int jd=0; jd<dim; jd++){
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+dim,jd);
            }
        }
    }

    TPZManVector<STATE,3> gradpressurefem(fDim,0.);
    this->Solution(data,VariableIndex("GradPressure"), gradpressurefem);

    TPZFNMatrix<3,STATE> fluxexactneg;

//sigma=-K grad(u)

    TPZFNMatrix<9,STATE> gradpressure(3,1);
    for (int i=0; i<3; i++) {
        gradpressure(i,0) = du_exact[i];
    }
    PermTensor.Multiply(gradpressure,fluxexactneg);



    REAL L2flux = 0., L2grad = 0.;
    for (int i=0; i<fDim; i++) {
        //   std::cout<<"fluxo fem = "<<fluxfem[i]<<std::endl;
        for (int j=0; j<fDim; j++) {
            L2flux += (fluxfem[i]+fluxexactneg(i,0))*InvPermTensor(i,j)*(fluxfem[j]+fluxexactneg(j,0));//Pq esta somando: o fluxo fem esta + e o exato -

        }

        // std::cout<<"grad sol = "<<gradpressurefem[i]<<std::endl;
        L2grad += (gradpressure(i,0)-gradpressurefem[i])*(gradpressure(i,0)-gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0]-u_exact[0])*(pressurefem[0]-u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;


}

void TPZMixedPoisson::ErrorsBC(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors,TPZBndCond &bc){
    //Add on Robin part the term ||e||_Gamma_R = K_R (e.n)^2
    // with e = sigma_fem - sigma_ex

    if(bc.Type() == 4){

    REAL InvKm = 1/bc.Val1()(0,0);

    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    int dim = fDim;


    REAL normalflux_fem = 0;//u.n
    normalflux_fem = data[0].sol[0][0];



    if(this->fExactSol){

       this->fExactSol->Execute(data[0].x,u_exact,du_exact);
    }


    TPZFNMatrix<9,REAL> PermTensor; this->GetPermeability(PermTensor);



    if(fPermeabilityFunction){
       PermTensor.Redim(3,3);
       TPZFNMatrix<3,STATE> resultMat;
       TPZManVector<STATE> res;
       fPermeabilityFunction->Execute(data[0].x,res,resultMat);
       for(int id=0; id<dim; id++){
           for(int jd=0; jd<dim; jd++){
               PermTensor(id,jd) = resultMat(id,jd);
           }
       }
    }

        TPZFNMatrix<3,REAL> fluxexactneg;

        //sigma=-K grad(u)

           TPZFNMatrix<9,REAL> gradpressure(3,1);
           for (int i=0; i<3; i++) {
               gradpressure(i,0) = du_exact[i];
           }
           PermTensor.Multiply(gradpressure,fluxexactneg);

        REAL normalflux_ex=0.;
        for(int i=0;i<fDim;i++){

            normalflux_ex += fluxexactneg(i,0)*data[0].normal[i];

        }


    REAL robinBCterm = 0.;

  //  std::cout<<"normalflux fem = "<<normalflux_fem<<" normalflux ex = "<<normalflux_ex<<std::endl;

    robinBCterm = (normalflux_fem +normalflux_ex)*InvKm*(normalflux_fem +normalflux_ex);//Pq esta somando: o fluxo fem esta + e o exato -


    errors[1] = robinBCterm;//L2 error for flux
    //errors[0] = L2 error for pressure
    //errors[2] = L2 for div
    //errors[3] = L2 for grad;
    //errors[4] = L2flux + residual;

    }
}

void TPZMixedPoisson::ErrorsBC(TPZVec<TPZMaterialData> &data,  TPZVec<REAL> &errors,TPZBndCond &bc){
    //Add on Robin part the term ||e||_Gamma_R = K_R (e.n)^2
    // with e = sigma_fem - sigma_ex

    if(bc.Type() == 4){

        REAL InvKm = 1/bc.Val1()(0,0);

        errors.Resize(NEvalErrors());
        errors.Fill(0.0);

        int dim = fDim;


        REAL normalflux_fem = 0;//u.n
        normalflux_fem = data[0].sol[0][0];


        TPZVec<STATE> u_exact(1,0); TPZFMatrix<STATE> du_exact(3,1,0);
        if(this->fExactSol){

            this->fExactSol->Execute(data[0].x,u_exact,du_exact);
        }


        TPZFNMatrix<9,REAL> PermTensor; this->GetPermeability(PermTensor);



        if(fPermeabilityFunction){
            PermTensor.Redim(3,3);
            TPZFNMatrix<3,STATE> resultMat;
            TPZManVector<STATE> res;
            fPermeabilityFunction->Execute(data[0].x,res,resultMat);
            for(int id=0; id<dim; id++){
                for(int jd=0; jd<dim; jd++){
                    PermTensor(id,jd) = resultMat(id,jd);
                }
            }
        }

        TPZFNMatrix<3,REAL> fluxexactneg;

        //sigma=-K grad(u)

        TPZFNMatrix<9,REAL> gradpressure(3,1);
        for (int i=0; i<3; i++) {
            gradpressure(i,0) = du_exact[i];
        }
        PermTensor.Multiply(gradpressure,fluxexactneg);

        REAL normalflux_ex=0.;
        for(int i=0;i<fDim;i++){

            normalflux_ex += fluxexactneg(i,0)*data[0].normal[i];

        }


        REAL robinBCterm = 0.;

        //  std::cout<<"normalflux fem = "<<normalflux_fem<<" normalflux ex = "<<normalflux_ex<<std::endl;

        robinBCterm = (normalflux_fem +normalflux_ex)*InvKm*(normalflux_fem +normalflux_ex);//Pq esta somando: o fluxo fem esta + e o exato -


        errors[1] = robinBCterm;//L2 error for flux
        //errors[0] = L2 error for pressure
        //errors[2] = L2 for div
        //errors[3] = L2 for grad;
        //errors[4] = L2flux + residual;

    }
}



int TPZMixedPoisson::ClassId() const{
    return Hash("TPZMixedPoisson") ^ TPZMatPoisson3d::ClassId() << 1;
}

