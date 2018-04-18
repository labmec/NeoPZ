/*
 *  TPZMohrCoulomb.h
 *  FEMPZ
 *
 *  Created by Diogo Cecilio on 5/4/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


/* Generated by Together */// $Id: TPZMohrCoulomb.h,v 1.2 2010-06-11 22:12:14 diogo Exp $

#ifndef TPZMOHRCOULOMBNETO_H
#define TPZMOHRCOULOMBNETO_H

#include "pzlog.h"
#include "TPZTensor.h"
#include "pzvec_extras.h"
#include "TPZPlasticState.h"

#ifdef LOG4CXX
static LoggerPtr loggerMohrCoulomb(Logger::getLogger("pz.plasticity.mohrcoulombneto"));
#endif



class TPZMohrCoulombNeto
{

    REAL fYoung;
    REAL fPoisson;
    REAL fPhi;
    REAL fPsi;
    REAL coesion;
    
public:
    
    /// Internal structure to represent the plastic memory (plastic deformation and damage)
    struct TPlasticState
    {
        TPlasticState() : fEpsPlastic(), fEpsPlasticBar(0.)
        {
            
        }
        
        TPlasticState(const TPlasticState &copy) : fEpsPlastic(copy.fEpsPlastic), fEpsPlasticBar(copy.fEpsPlasticBar)
        {
            
        }
        
        TPlasticState &operator=(const TPlasticState &copy)
        {
            fEpsPlastic = copy.fEpsPlastic;
            fEpsPlasticBar = copy.fEpsPlasticBar;
            return *this;
        }
        
        void Print(std::ostream &out) const
        {
            out << "Plastic Deformation tensor ";
            fEpsPlastic.Print(out);
            out << "Acumulated plastic deformation " << fEpsPlasticBar << std::endl;
        }
        
        /// plastic deformation tensor
        TPZTensor<REAL> fEpsPlastic;
        
        /// accumulated damage
        REAL fEpsPlasticBar;
    };
    
    
    /// structure which contains the decision tree of the return map
    // we can only expect a consistent tangent matrix if the decision tree remains the same
    struct TComputeSequence
    {
        TComputeSequence() : fWhichPlane(ENoPlane), fGamma(0)
        {
            
        }
        
        TComputeSequence(const TComputeSequence &copy) : fWhichPlane(copy.fWhichPlane), fGamma(copy.fGamma)
        {
            
        }
        
        TComputeSequence &operator=(const TComputeSequence &copy)
        {
            fWhichPlane = copy.fWhichPlane;
            fGamma = copy.fGamma;
            return *this;
        }
        
        enum MPlane {ENoPlane, EElastic, EMainPlane, ERightEdge, ELeftEdge, EHydroStatic };
        
        MPlane fWhichPlane;
        
        TPZManVector<REAL> fGamma;
        

    };
    
protected:
    
    /// information of the plastic state of the material point
    TPZMohrCoulombNeto::TPlasticState fState;

    
public:
    
    TPZMohrCoulombNeto() : fYoung(20000.), fPoisson(0.), fPhi(M_PI/9.), fPsi(M_PI/9.) 
    {
        this->coesion = 9.35;
    }
    REAL Lambda()
    {
        return fPoisson*fYoung/(1.+fPoisson)*(1.-2.*fPoisson);
    }
    REAL Mu()
    {
        return fYoung/(2.*(1.+fPoisson));
    }
    REAL G()
    {
        return fYoung/(2.*(1+fPoisson));
    }
    REAL K()
    {
        return fYoung/(3.*(1.-2.*fPoisson));
    }
    
    void Print(std::ostream &out) const
    {
        out << "TPZMohrCoulombNeto\n";
        fState.Print(out);
    }
    
    /// the hardening function and its derivative
    template<class T>
    void PlasticityFunction(T epsp, T &sigmay, T &H) const
    {
        //sigmay = T(15.)+(T(2571.43)-T(2.95238e6)*epsp)*(T(-0.0035)+epsp);
       // H = T(12904.8)-T(5.90476e6)*epsp;
       //sigmay=20.;
       // H=0.;
            sigmay = T(15.)+(T(2571.43))*(T(-0.0035)) + T(12904.8)*epsp;
            H = T(12904.8);
    }
    
    /// a piecewise linear hardening function
    template<class T>
    void PieceWise(T epbar, T &m, T & fx)const
    {
        
        ifstream file("curvadehardening.txt");
        TPZFMatrix<REAL> mat;
        int sz;
        file >>sz;
        mat.Resize(sz,2);
        
        for(int i=0;i<sz-1;i++)
        {
            
            file >> mat(i,0);
            file >> mat(i,1);
        }
        
        //cout << "\n mat = "<< mat <<endl;
        for(int i=0;i<sz-1;i++)
        {
            if(epbar >= mat(i,0) && epbar <=  mat(i+1,0))
            {
                REAL x0 = mat(i,0);
                REAL y0 = mat(i,1);
                REAL x = mat(i+1,0);
                REAL y = mat(i+1,1);
                m = (y - y0)/(x - x0);
                fx=m*(epbar-x0)+y0;
                // return;
            }
            
        }
    }
    
    
    template<class T>
    TPZTensor<T> SigmaElast(const TPZTensor<T> &deform)
    {
        T trdeform = deform.I1();
        TPZTensor<T> result;
        result.Identity();
        result *= (Lambda()*trdeform);
        result.Add(deform,2.*Mu());
        return result;
    }
    
    template<class T>
    typename TPZTensor<T>::TPZDecomposed SigmaTrial(const TPZTensor<T> &epstotal)
    {
        TPZTensor<T> epslocal(epstotal);
        for(int i=0; i<6; i++)
        {
            epslocal[i] -= T(fState.fEpsPlastic[i]);
        }
        TPZTensor<T> sigma;
        sigma = SigmaElast(epslocal);
        typename TPZTensor<T>::TPZDecomposed sigma_trial;
        sigma.EigenSystem(sigma_trial);
        
#ifdef LOG4CXX
        if (loggerMohrCoulomb->isDebugEnabled()) {
            std::stringstream sout;
            sout << "Input stress tensor ";
            sigma.Print(sout);
            sout << "Input tensor in decomposed form\n";
            sigma_trial.Print(sout);
            LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
        }
#endif
        return sigma_trial;
    }
    
    void ComputeSigmaTangent(TPZTensor<REAL> &epstotal, TPZTensor<REAL> &sigma, TPZFNMatrix<36,REAL> &tangent, const TComputeSequence &memory)
    {
        typedef TFad<6,REAL> fadtype;
        TPZTensor<fadtype> epstotalFAD, sigmaElastFAD, sigmaFAD;
        for (int i=0; i<6; i++) {
            epstotalFAD[i].val() = epstotal[i];
            epstotalFAD[i].fastAccessDx(i) = 1.;
        }
        sigmaElastFAD = SigmaElast(epstotalFAD);
        
        switch (memory.fWhichPlane) {
            case TComputeSequence::ENoPlane:
                DebugStop();
                break;
            case TComputeSequence::EElastic:
                sigmaFAD = sigmaElastFAD;
                break;
            case TComputeSequence::EMainPlane:
            case TComputeSequence::ELeftEdge:
            case TComputeSequence::ERightEdge:
            {
                TPZTensor<fadtype>::TPZDecomposed sigma_trial = SigmaTrial(epstotalFAD);
                TPZTensor<fadtype>::TPZDecomposed sigma_projected;
                TComputeSequence locmem(memory);
                switch (memory.fWhichPlane) {
                    case TComputeSequence::EMainPlane:
                        ReturnMapPlane<fadtype>(sigma_trial, sigma_projected, locmem);                        
                        break;
                    case TComputeSequence::ELeftEdge:
                        ReturnMapLeftEdge<fadtype>(sigma_trial, sigma_projected, locmem);
                        break;
                    case TComputeSequence::ERightEdge:
                        ReturnMapRightEdge<fadtype>(sigma_trial, sigma_projected, locmem);
                        break;
                    default:
                        DebugStop();
                        break;
                }
                sigmaFAD = TPZTensor<fadtype>(sigma_projected);
            }   
            default:
                break;
        }
        for (int i=0; i<6; i++) {
            sigma[i] = sigmaFAD[i].val();
            for (int j=0; j<6; j++) {
                tangent(i,j) = sigmaFAD[i].fastAccessDx(j);
            }
        }
    }
    
    void CommitDeformation(TPZTensor<REAL> &epstotal, TComputeSequence &memory)
    {
        const REAL cosphi = cos(fPhi);
        TPZTensor<REAL> sigma;
        switch (memory.fWhichPlane) {
            case TComputeSequence::EElastic:
            {
                TPZTensor<REAL> epslocal(epstotal);
                epslocal -= fState.fEpsPlastic;
                sigma = SigmaElast(epslocal);
                //sigma = SigmaElast(epstotal);
            }
                break;
            case TComputeSequence::EMainPlane:
            case TComputeSequence::ELeftEdge:
            case TComputeSequence::ERightEdge:
            {
                TPZTensor<REAL>::TPZDecomposed sigma_trial = SigmaTrial(epstotal);
                TPZTensor<REAL>::TPZDecomposed sigma_projected;
                TComputeSequence locmem(memory);
                switch (memory.fWhichPlane) {
                    case TComputeSequence::EMainPlane:
                        ReturnMapPlane<REAL>(sigma_trial, sigma_projected, locmem);
                        fState.fEpsPlasticBar+=(locmem.fGamma[0]*2.*cosphi);
                        break;
                    case TComputeSequence::ELeftEdge:
                        ReturnMapLeftEdge<REAL>(sigma_trial, sigma_projected, locmem);
                        fState.fEpsPlasticBar+=(locmem.fGamma[0]+locmem.fGamma[1])*2.*cosphi;
                        break;
                    case TComputeSequence::ERightEdge:
                        ReturnMapRightEdge<REAL>(sigma_trial, sigma_projected, locmem);
                        fState.fEpsPlasticBar+=(locmem.fGamma[0]+locmem.fGamma[1])*2.*cosphi;
                        break;
                    default:
                        DebugStop();
                        break;
                }
                 sigma = TPZTensor<REAL>(sigma_projected);
                break;
            }
            default:
            {
                break;
            }
        }

        
		 	REAL tempval;
		 	TPZTensor<REAL> StressDeviatoric,P,I,epsplastic(epstotal),epselastic;
		
		 	P.XX()=1; I.XX()=1;
		 	P.YY()=1; I.YY()=1;
			P.ZZ()=1; I.ZZ()=1;
		
		
		 	tempval=(sigma.I1()/3.)*(1./(3.*K()));
		 	P.Multiply(tempval,1);
		
		
		 	cout << " \n P = "<< P <<endl;
		 	sigma.S(StressDeviatoric);
		
		 	StressDeviatoric*=1./(2.*G());
		 	cout << " \n S = "<< StressDeviatoric <<endl;
		
		 	StressDeviatoric.Add(P,1);
		 	epselastic=StressDeviatoric;
        
		 	epsplastic-=epselastic;
        
        
            fState.fEpsPlastic = epsplastic;
        

        
    }
    
    template<class T>
    TComputeSequence ComputeSigma(TPZTensor<T> &epstotal, TPZTensor<T> &sigma)
    {
        typename TPZTensor<T>::TPZDecomposed sigma_trial = SigmaTrial(epstotal);
        TComputeSequence memory;
        T phi = PhiPlane<T>(sigma_trial);
        if (TPZExtractVal::val(phi) <= 0.) {
            memory.fWhichPlane = TComputeSequence::EElastic;
            memory.fGamma.Resize(0);
            sigma = TPZTensor<T>(sigma_trial);
              //state.fEpsT = epstotal;
            return memory;
        }
        typename TPZTensor<T>::TPZDecomposed sigma_projected;
        memory.fGamma.Resize(1);
        memory.fGamma[0] = 0.;
        if (ReturnMapPlane<T>(sigma_trial, sigma_projected, memory)) {
            sigma = TPZTensor<T>(sigma_projected);
            memory.fWhichPlane = TComputeSequence::EMainPlane;
        }
        else {
            memory.fGamma.Resize(2);
            memory.fGamma[0] = 0.;
            memory.fGamma[1] = 0.;


            const REAL sinpsi = sin(fPsi);
            TPZManVector<T,3> &eigenvalues = sigma_trial.fEigenvalues;
            REAL val = (1-sinpsi)*TPZExtractVal::val(eigenvalues[0])-2.*TPZExtractVal::val(eigenvalues[2])+(1+sinpsi)*TPZExtractVal::val(eigenvalues[1]);
            if (val > 0.) {
                ReturnMapRightEdge<T>(sigma_trial, sigma_projected, memory);
                memory.fWhichPlane = TComputeSequence::ERightEdge;
            }
            else {
                ReturnMapLeftEdge<T>(sigma_trial, sigma_projected, memory);
                memory.fWhichPlane = TComputeSequence::ELeftEdge;
            }
#ifdef LOG4CXX
            {
                std::stringstream sout;
                sout << "After the map to the edge, sigma_projected :\n";
                sigma_projected.Print(sout);
                LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
            }
#endif
            sigma = TPZTensor<T>(sigma_projected);
        }
        CommitDeformation(epstotal,memory);
        return memory;
    }
    
    template<class T>
    T PhiPlane(typename TPZTensor<T>::TPZDecomposed &sigma) const
    {
        const REAL sinphi = sin(fPhi);
        const REAL cosphi = cos(fPhi);
        T sigmay,H;
        PlasticityFunction(T(fState.fEpsPlasticBar),sigmay, H);
        return sigma.fEigenvalues[0]-sigma.fEigenvalues[2]+(sigma.fEigenvalues[0]+sigma.fEigenvalues[2])*sinphi-2.*sigmay*cosphi;
    }

    template<class T>
    bool ReturnMapPlane(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected, 
                            TComputeSequence &memory)
    {
        sigma_projected = sigma_trial;
        TPZManVector<T,3> &eigenvalues = sigma_projected.fEigenvalues;
//        TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;
        const REAL sinphi = sin(fPhi);
        const REAL sinpsi = sin(fPsi);
        const REAL cosphi = cos(fPhi);
        const REAL sinphi2 = sinphi*sinphi;
        const REAL cosphi2 = 1.-sinphi2;
        const REAL constA = 4.* G() *(1.+ sinphi*sinpsi/3.) + 4.*K() * sinphi*sinpsi;
        T sigmay,H;
        T epsbar = T(fState.fEpsPlasticBar+memory.fGamma[0]*2.*cosphi);//diogo aqui
        PlasticityFunction(epsbar,sigmay, H);
        T phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi;
        T gamma = memory.fGamma[0];
        REAL phival = TPZExtractVal::val(phi);
        REAL tolerance = 1.e-8;
        do {
            T denom = -constA- T(4.*cosphi2)*H;
//            T d = T(-4.*G()*(1.+sinphi*sinpsi/3.)-4.*K()*sinphi*sinpsi)-T(4.*cosphi2)*H;
            T deriv_gamma = -phi/denom;
            gamma += deriv_gamma;
            epsbar = T(fState.fEpsPlasticBar)+gamma*T(2.*cosphi);///errado esta inicializando toda vez. diogo aqui
            PlasticityFunction(epsbar, sigmay, H);
            if (TPZExtractVal::val(H) < 0.) {
                DebugStop();
            }
            phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi-constA*gamma;
            phival = TPZExtractVal::val(phi);
            
        } while (abs(phival) > tolerance);
        

        memory.fGamma[0] = TPZExtractVal::val(gamma);
        eigenvalues[0] -= T(2.*G()*(1+sinpsi/3.)+2.*K()*sinpsi)*gamma;
        eigenvalues[1] += T((4.*G()/3. - K()*2.)*sinpsi)*gamma;
        eigenvalues[2] += T(2.*G()*(1-sinpsi/3.)-2.*K()*sinpsi)*gamma;
#ifdef PZDEBUG
        phi = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*sinphi-2.*sigmay*cosphi;
#endif
        return (TPZExtractVal::val(eigenvalues[0])>TPZExtractVal::val(eigenvalues[1]) && TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]));
    }
    
    template<class T>
    bool ReturnMapLeftEdge(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected,
                           TComputeSequence &memory)
    {
        
        sigma_projected = sigma_trial;
        TPZManVector<T,3> &eigenvalues = sigma_projected.fEigenvalues;
//        TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;        
        const REAL sinphi = sin(fPhi);
        const REAL sinpsi = sin(fPsi);
        const REAL cosphi = cos(fPhi);
        const REAL sinphi2 = sinphi*sinphi;
        const REAL cosphi2 = 1.-sinphi2;
        TPZManVector<T,2> gamma(2,0.),phi(2,0.),sigma_bar(2,0.),ab(2,0.);
        gamma[0] = memory.fGamma[0];
        gamma[1] = memory.fGamma[1];
        TPZManVector<REAL,2> phival(2,0.);
        TPZFNMatrix<4,T> d(2,2,0.), dinverse(2,2,0.);
        sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
        sigma_bar[1] = eigenvalues[1]-eigenvalues[2]+(eigenvalues[1]+eigenvalues[2])*T(sinphi);
        T sigmay,H;
        T epsbar = T(fState.fEpsPlasticBar) + (gamma[0]+gamma[1])*T(2.*cosphi);//diogo aqui
        PlasticityFunction(epsbar,sigmay, H);
        phi[0] = sigma_bar[0] - T(2.*cosphi)*sigmay;
        phi[1] = sigma_bar[1] - T(2.*cosphi)*sigmay;
        ab[0] = T(4.*G()*(1+sinphi*sinpsi/3.)+4.*K()*sinphi*sinpsi);
        ab[1] = T(2.*G()*(1.-sinphi-sinpsi-sinphi*sinpsi/3.)+4.*K()*sinphi*sinpsi);
        T residual =1;
        REAL tolerance = 1.e-8;
        do {
            d(0,0) = -ab[0]-T(4.*cosphi2)*H;
            d(1,0) = -ab[1]-T(4.*cosphi2)*H;
            d(0,1) = -ab[1]-T(4.*cosphi2)*H;
            d(1,1) = -ab[0]-T(4.*cosphi2)*H;
            T detd = d(0,0)*d(1,1)-d(0,1)*d(1,0);
            dinverse(0,0) = d(1,1)/detd;
            dinverse(1,0) = -d(1,0)/detd;
            dinverse(0,1) = -d(0,1)/detd;
            dinverse(1,1) = d(0,0)/detd;
            gamma[0] -= (dinverse(0,0)*phi[0]+dinverse(0,1)*phi[1]);
            gamma[1] -= (dinverse(1,0)*phi[0]+dinverse(1,1)*phi[1]);
          //T epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi); diogo aqui
            epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
            PlasticityFunction(epsbar, sigmay, H);
            phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T(2.*cosphi)*sigmay;
            phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[0] - T(2.*cosphi)*sigmay;
            phival[0] = TPZExtractVal::val(phi[0]);
            phival[1] = TPZExtractVal::val(phi[1]);
            residual=(fabs(phival[0])+fabs(phival[1]))/sigmay;//aqui diogo
        //} while (abs(phival[0]) > tolerance || abs(phival[1]) > tolerance);//aqui diogo
        }while (residual>tolerance);//aqui diogo
//        eigenvalues[0] -= T(2.*G()*(1+sinpsi/3.)+2.*K()*sinpsi)*gamma;
//        eigenvalues[1] += T((4.*G()/3. - K()*2.)*sinpsi)*gamma;
//        eigenvalues[2] += T(2.*G()*(1-sinpsi/3.)-2.*K()*sinpsi)*gamma;

        memory.fGamma[0] = TPZExtractVal::val(gamma[0]);
        memory.fGamma[1] = TPZExtractVal::val(gamma[1]);
        eigenvalues[0] -= T(2.*G()*(1+sinpsi/3.)+2.*K()*sinpsi)*gamma[0]+T((4.*G()/3.-2.*K())*sinpsi)*gamma[1];
        eigenvalues[1] += T((4.*G()/3.- K()*2.)*sinpsi)*gamma[0]-T(2.*G()*(1.+sinpsi/3.)+2.*K()*sinpsi)*gamma[1];
        eigenvalues[2] -= T(2.*G()*(1-sinpsi/3.)-2.*K()*sinpsi)*(gamma[0]+gamma[1]);
        return (TPZExtractVal::val(eigenvalues[0])>TPZExtractVal::val(eigenvalues[1]) && TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]));
    }
    
    template<class T>
    bool ReturnMapRightEdge(const typename TPZTensor<T>::TPZDecomposed &sigma_trial, typename TPZTensor<T>::TPZDecomposed &sigma_projected,
                            TComputeSequence &memory)
    {
        sigma_projected = sigma_trial;
        TPZManVector<T,3> &eigenvalues = sigma_projected.fEigenvalues;
//      TPZManVector<TPZTensor<T>,3> &eigenvectors = sigma_projected.fEigenvectors;
        const REAL sinphi = sin(fPhi);
        const REAL sinpsi = sin(fPsi);
        const REAL cosphi = cos(fPhi);
        const REAL sinphi2 = sinphi*sinphi;
        const REAL cosphi2 = 1.-sinphi2;
        const REAL KV = K();
        const REAL GV = G();
        TPZManVector<T,2> gamma(2,0.),phi(2,0.),sigma_bar(2,0.),ab(2,0.);
        gamma[0] = memory.fGamma[0];
        gamma[1] = memory.fGamma[1];
        TPZManVector<REAL,2> phival(2,0.);
        TPZFNMatrix<4,T> d(2,2,0.), dinverse(2,2,0.);
        sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
        sigma_bar[1] = eigenvalues[0]-eigenvalues[1]+(eigenvalues[0]+eigenvalues[1])*T(sinphi);
        T sigmay,H;
        T epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
        PlasticityFunction(epsbar,sigmay, H);
        phi[0] = sigma_bar[0] - T(2.*cosphi)*sigmay;
        phi[1] = sigma_bar[1] - T(2.*cosphi)*sigmay;
        ab[0] = T(4.*GV*(1+sinphi*sinpsi/3.)+4.*KV*sinphi*sinpsi);
        ab[1] = T(2.*GV*(1.+sinphi+sinpsi-sinphi*sinpsi/3.)+4.*KV*sinphi*sinpsi);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "phi = " << phi << std::endl;
            LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
        }
#endif
        
//#ifdef PZDEBUG
//        gamma[0] = 0.;
//        gamma[1] = 1.;
//        T v[3];
//        v[0] = -T(2.*GV*(1+sinpsi/3.)+2.*KV*sinpsi)*(gamma[0]+gamma[1]);
//        v[1] = T((4.*GV/3.- KV*2.)*sinpsi)*gamma[0]+T(2.*GV*(1.-sinpsi/3.)-2.*KV*sinpsi)*gamma[1];
//        v[2] = T(2.*GV*(1-sinpsi/3.)-2.*KV*sinpsi)*gamma[0]+T((4.*GV/3.-2.*KV)*sinpsi)*gamma[1];
//        eigenvalues[0] = sigma_trial.fEigenvalues[0]+v[0];
//        eigenvalues[1] = sigma_trial.fEigenvalues[1]+v[1];
//        eigenvalues[2] = sigma_trial.fEigenvalues[2]+v[2];
//        
//        T test1 = (v[0]-v[2])+(v[0]+v[2])*sinphi;
//        T test2 = (v[0]-v[1])+(v[0]+v[1])*sinphi;
//
//        sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
//        sigma_bar[1] = eigenvalues[0]-eigenvalues[1]+(eigenvalues[0]+eigenvalues[1])*T(sinphi);
//        T A = phi[0] - sigma_bar[0] + T(2.*cosphi)*sigmay;
//        T B = phi[1] - sigma_bar[1] + T(2.*cosphi)*sigmay;
//
//#endif
        REAL tolerance = 1.e-8;
        int iter = 0;
        T residual =1;
        do {
#ifdef LOG4CXX
            {
                std::stringstream sout;
                sout << "epsbar = " << epsbar << std::endl;
                sout << "sigmay = " << sigmay << std::endl;
                sout << "H = " << H << std::endl;
                LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
            }
#endif
            d(0,0) = -ab[0]-T(4.*cosphi2)*H;
            d(1,0) = -ab[1]-T(4.*cosphi2)*H;
            d(0,1) = -ab[1]-T(4.*cosphi2)*H;
            d(1,1) = -ab[0]-T(4.*cosphi2)*H;
            T detd = d(0,0)*d(1,1)-d(0,1)*d(1,0);
            dinverse(0,0) = d(1,1)/detd;
            dinverse(1,0) = -d(1,0)/detd;
            dinverse(0,1) = -d(0,1)/detd;
            dinverse(1,1) = d(0,0)/detd;
            gamma[0] -= (dinverse(0,0)*phi[0]+dinverse(0,1)*phi[1]);
            gamma[1] -= (dinverse(1,0)*phi[0]+dinverse(1,1)*phi[1]);
            epsbar = T(fState.fEpsPlasticBar)+(gamma[0]+gamma[1])*T(2.*cosphi);
            PlasticityFunction(epsbar, sigmay, H);
            if (TPZExtractVal::val(H) < 0.) {
                DebugStop();
            }
            iter++;
            phi[0] = sigma_bar[0] - ab[0]*gamma[0] - ab[1]*gamma[1] - T(2.*cosphi)*sigmay;
            phi[1] = sigma_bar[1] - ab[1]*gamma[0] - ab[0]*gamma[1] - T(2.*cosphi)*sigmay;
            phival[0] = TPZExtractVal::val(phi[0]);
            phival[1] = TPZExtractVal::val(phi[1]);
#ifdef LOG4CXX
            {
                std::stringstream sout;
                sout << "iter = " << iter << " phi = " << phival << std::endl;
                LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
            }
#endif
            residual=(fabs(phival[0])+fabs(phival[1]))/sigmay;//aqui diogo
            cout << "\n residula = "<< endl;
            //} while (abs(phival[0]) > tolerance || abs(phival[1]) > tolerance);//aqui diogo
        }while (residual>tolerance);//aqui diogo
        
//        eigenvalues[0] -= T(2.*GV*(1+sinpsi/3.)+2.*KV*sinpsi)*gamma;
//        eigenvalues[1] += T((4.*GV/3. - KV*2.)*sinpsi)*gamma;
//        eigenvalues[2] += T(2.*GV*(1-sinpsi/3.)-2.*KV*sinpsi)*gamma;
       // epsbar = T(state.fAlpha)+(gamma[0]+gamma[1])*T(2.*cosphi);
       
            
        memory.fGamma[0] = TPZExtractVal::val(gamma[0]);
        memory.fGamma[1] = TPZExtractVal::val(gamma[1]);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "gamma = " << gamma << std::endl;
            sout << "phival = " << phival << std::endl;
            sout << "ab = " << ab << std::endl;
            sout << "sigma_bar = " << sigma_bar << std::endl;
            d.Print("Jacobian",sout);
            dinverse.Print("Inverse Jacobian",sout);
            sout << "epsbar = " << epsbar << std::endl;
            LOGPZ_DEBUG(loggerMohrCoulomb, sout.str())
        }
#endif
        eigenvalues[0] -= T(2.*GV*(1+sinpsi/3.)+2.*KV*sinpsi)*(gamma[0]+gamma[1]);
        eigenvalues[1] += T((4.*GV/3.- KV*2.)*sinpsi)*gamma[0]+T(2.*GV*(1.-sinpsi/3.)-2.*KV*sinpsi)*gamma[1];
        eigenvalues[2] += T(2.*GV*(1-sinpsi/3.)-2.*KV*sinpsi)*gamma[0]+T((4.*GV/3.-2.*KV)*sinpsi)*gamma[1];
/*#ifdef PZDEBUG
        sigma_bar[0] = eigenvalues[0]-eigenvalues[2]+(eigenvalues[0]+eigenvalues[2])*T(sinphi);
        sigma_bar[1] = eigenvalues[0]-eigenvalues[1]+(eigenvalues[0]+eigenvalues[1])*T(sinphi);
        //epsbar = T(state.fAlpha)+(gamma[0]+gamma[1])*T(2.*cosphi);
        PlasticityFunction(epsbar, sigmay, H);

        phi[0] = sigma_bar[0] - T(2.*cosphi)*sigmay;
        phi[1] = sigma_bar[1] - T(2.*cosphi)*sigmay;
#endif*/
        return (TPZExtractVal::val(eigenvalues[0])>TPZExtractVal::val(eigenvalues[1]) && TPZExtractVal::val(eigenvalues[1]) > TPZExtractVal::val(eigenvalues[2]));        
    }
};


#endif //TPZMohrCoulomb_H
