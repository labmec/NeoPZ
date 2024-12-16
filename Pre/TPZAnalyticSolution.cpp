
#include "TPZAnalyticSolution.h"



#include "pzlog.h"

#include "fadType.h"


#ifndef STATE_COMPLEX
static Fad<STATE> atan2(Fad<STATE> y, Fad<STATE> x)
{
    int sz = x.size();
    STATE xval,yval;
    xval = x.val();
    yval = y.val();
    STATE angle = atan2(yval,xval);
    Fad<STATE> result(sz,angle);
    STATE r = x.val()*x.val()+y.val()*y.val();
    for (int i=0; i<sz; i++) {
        result.fastAccessDx(i) = (-y.val()*x.fastAccessDx(i) + x.val()*y.fastAccessDx(i))/r;
    }
    return result;
}
#else
static Fad<STATE> atan2(Fad<STATE> y, Fad<STATE> x)
{
    int sz = x.size();
    double xval,yval,xvalimag,yvalimag;
    xval = x.val().real();
    yval = y.val().real();
    xvalimag = x.val().imag();
    yvalimag = y.val().imag();
    if(abs(xvalimag) > 1.e-12 || abs(yvalimag)> 1.e-12)
    {
        DebugStop();
    }
    REAL angle = atan2(yval,xval);
    Fad<STATE> result(sz,angle);
    double r = xval*xval+yval*yval;
    for (int i=0; i<sz; i++) {
        result.fastAccessDx(i) = (-y.val()*x.fastAccessDx(i) + x.val()*y.fastAccessDx(i))/r;
    }
    return result;
}
#endif

static FADFADSTATE FADatan2(FADFADSTATE y, FADFADSTATE x)
{
    int sz = x.size();
    FADFADSTATE result(sz,atan2(y.val(),x.val()));
    Fad<STATE> r = x.val()*x.val()+y.val()*y.val();
    for (int i=0; i<sz; i++) {
        result.fastAccessDx(i) = (-y.val()*x.fastAccessDx(i) + x.val()*y.fastAccessDx(i))/r;
    }
    return result;
}

static FADFADSTATE FADsin(FADFADSTATE x)
{
    
    Fad<STATE> sinaval = sin(x.val());
    Fad<STATE> cosaval = cos(x.val());
    int sz = x.size();
    FADFADSTATE sina(sz,sinaval);
    for (int i=0; i<sz; i++) {
        sina.fastAccessDx(i) = cosaval*x.dx(i);
    }
    return sina;
}

static FADFADSTATE FADcos(FADFADSTATE x)
{
    Fad<STATE> sinaval = sin(x.val());
    Fad<STATE> cosaval = cos(x.val());
    int sz = x.size();
    FADFADSTATE cosa(sz,cosaval);
    for (int i=0; i<sz; i++) {
        cosa.fastAccessDx(i) = -sinaval*x.dx(i);
    }
    return cosa;
}

static FADFADSTATE FADexp(FADFADSTATE x)
{
    Fad<STATE> expaval = exp(x.val());
    int sz = x.size();
    FADFADSTATE expa(sz,expaval);
    for (int i=0; i<sz; i++) {
        expa.fastAccessDx(i) = expaval*x.dx(i);
    }
    return expa;
}

static FADFADSTATE FADsqrt(FADFADSTATE x)
{
    Fad<STATE> fadres = sqrt(x.val());
    int sz = x.size();
    FADFADSTATE resa(sz,fadres);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = REAL(0.5)/fadres*x.dx(i);
    }
    return resa;
}

static FADFADSTATE FADatan(FADFADSTATE x)
{
    Fad<STATE> fadres = atan(x.val());
    int sz = x.size();
    FADFADSTATE resa(sz,fadres);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = 1./(1.+x.val()*x.val())*x.dx(i);
    }
    return resa;
}

static FADFADSTATE FADcosh(FADFADSTATE x)
{
    Fad<STATE> coshval = cosh(x.val());
    Fad<STATE> sinhval = sinh(x.val());
    int sz = x.size();
    FADFADSTATE resa(sz,coshval);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = sinhval*x.dx(i);
    }
    return resa;
}

static FADFADSTATE FADsinh(FADFADSTATE x)
{
    Fad<STATE> coshval = cosh(x.val());
    Fad<STATE> sinhval = sinh(x.val());
    int sz = x.size();
    FADFADSTATE resa(sz,sinhval);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = coshval*x.dx(i);
    }
    return resa;
}

static const REAL FI = 1.;
static const REAL a = 0.5;
static const REAL b = 0.5;

TPZAnalyticSolution::TPZAnalyticSolution(const TPZAnalyticSolution &cp) :fSignConvention(cp.fSignConvention) {
    std::cout << "TPZAnalyticSolution::TPZAnalyticSolution(const TPZAnalyticSolution &cp): One should not invoke this copy constructor";
    DebugStop();
}

TPZAnalyticSolution & TPZAnalyticSolution::operator=(const TPZAnalyticSolution &copy) {
    std::cout << "TPZAnalyticSolution & TPZAnalyticSolution::operator=(const TPZAnalyticSolution &copy): One should not invoke this copy constructor";
    DebugStop();
    fSignConvention = copy.fSignConvention;
    return *this;
}

REAL TElasticity2DAnalytic::gE = 1.;

REAL TElasticity2DAnalytic::gPoisson = 0.3;

int TElasticity2DAnalytic::gOscilatoryElasticity = 0;

TElasticity2DAnalytic::TElasticity2DAnalytic(const TElasticity2DAnalytic &cp) : TPZAnalyticSolution(cp),fProblemType(cp.fProblemType) {
    std::cout << "TElasticity2DAnalytic::TElasticity2DAnalytic(const TElasticity2DAnalytic &cp): One should not invoke this copy constructor";
    DebugStop();
}

TElasticity2DAnalytic & TElasticity2DAnalytic::operator=(const TElasticity2DAnalytic &copy){
    fProblemType = copy.fProblemType;
    gE = copy.gE;
    gPoisson = copy.gPoisson;
    return *this;
}

template<>
void TElasticity2DAnalytic::uxy(const TPZVec<FADFADSTATE > &x, TPZVec<FADFADSTATE > &disp) const
{
    typedef FADFADSTATE TVar;
    if(fProblemType == Etest1)
    {   
        FADFADSTATE tmp = (FADFADSTATE)(1./27.)*x[0]*x[0]*x[1]*x[1];
        disp[0] = tmp*FADcos((FADFADSTATE)(6.*M_PI)*x[0])*FADsin((FADFADSTATE)(7.*M_PI)*x[1]);
        disp[1] = (FADFADSTATE)(0.2)*FADexp(x[1])*FADsin((FADFADSTATE)(4.*M_PI)*x[0]);
    
    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1.-x[0]*x[0])*(1.+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1.-x[1]*x[1])*(1.+x[0]*x[0]*x[0]*x[0]));
    }
      
    else if(fProblemType ==ERot)//rotation
    {      
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] =(FADFADSTATE)-x[1];
        disp[1] =(FADFADSTATE) x[0];
      
    }
    
    else if(fProblemType == EShear)//pure shear
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) x[1];
        disp[1] += (FADFADSTATE) 0. ;
    }
    else if(fProblemType == EStretchx)//strech x
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) x[0];
        disp[1] += (FADFADSTATE) 0.;
    }
    else if(fProblemType == EUniAxialx)
    {
        if (fPlaneStress == 0) {
            disp[0] = x[0]*(1.-gPoisson*gPoisson)/gE;
            disp[1] = -x[1]*(1.+gPoisson)*gPoisson/gE;
        }
        else
        {
            disp[0] = x[0]/gE;
            disp[1] = -x[1]*gPoisson/gE;
        }
    }
    else if(fProblemType ==EStretchy)//strech y
    {    
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) 0.;
        disp[1] += (FADFADSTATE) x[1];
    }
    else if(fProblemType==EDispx)
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[1] +=   0.;
    }
    else if(fProblemType==EDispy)
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) 0.;
        disp[1] += (FADFADSTATE) 1.;
    }
    else if (fProblemType == EThiago){
//        disp[0] = (x[0])*(x[0]-1)*(x[1])*(x[1]-1);
//        disp[1] = (x[0])*(x[0]-1)*(x[1])*(x[1]-1);
        disp[0] = FADcos(M_PI * x[0]) * FADsin(2 * M_PI * x[1]);
        disp[1] = FADcos(M_PI * x[1]) * FADsin(M_PI * x[0]);
 
    }
    
    else if (fProblemType==EHomoDir){
        disp[0] =  FADsin(2 * M_PI * x[1])*(-1.+FADcos(2*M_PI*x[0]))+FADsin( M_PI * x[1])*FADsin( M_PI * x[0]);
        disp[1] = FADsin(2 * M_PI * x[0])*(1.-FADcos(2*M_PI*x[1]))+FADsin( M_PI * x[1])*FADsin( M_PI * x[0]);
        
    }
    
    else if(fProblemType == EPoly){
        disp[0] = 0. * x[0]; //*x[0]*x[1];
        disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
    } else if(fProblemType==EBend)
    {
        typedef TVar FADFADSTATE;
        TVar poiss = gPoisson;
        TVar elast = gE;
        if(fPlaneStress == 0)
        {
            poiss = poiss/(1.-poiss);
            elast /= (1-gPoisson*gPoisson);
        }
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est,nust,G;
        REAL MI = 5, h = 1.;
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1-gPoisson);
//            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
//            nust = gPoisson/(1.+gPoisson);
        }
        else
        {
            Est = gE;
            nust = gPoisson;
        }
        disp[0] = MI*h*h*x[1]/(2.*G)+ MI*x[0]*x[0]*x[1]/(2.*Est)-MI *x[1]*x[1]*x[1]/(6.*G)+MI*nust*x[1]*x[1]*x[1]/(6.*Est);
        disp[1] = -MI*x[0]*x[0]*x[0]/(6.*Est)-MI*nust*x[0]*x[1]*x[1]/(2.*Est);
    }
    else if(fProblemType == ESquareRoot)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1.-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1.+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    else if(fProblemType == ESquareRootLower)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        if (shapeFAD::val(theta) > 0.) {
            theta -= (2.*M_PI);
        }
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1.-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1.+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    else if(fProblemType == ESquareRootUpper)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        if (shapeFAD::val(theta) < 0.) {
            theta += (2.*M_PI);
        }
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1.-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1.+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1./(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    
    else if (fProblemType==ELShape){
        
        TVar mu=1.;
        TVar lambda=5.;
        TVar alfa=0.6;
        TVar A= 31./9.;
        TVar theta = FADatan2(x[1],x[0]);
        
        auto thetaval = shapeFAD::val(theta);
        if (thetaval < (0.)) theta += 2. * M_PI;
        TVar eps = 1e-15;
        
    
        // Verification to avoid numerical errors when x > 0 and y = 0
        if (x[0] > TVar (0.) && x[1] < TVar (1e-15) &&  x[1] >  TVar (-1e-15)) {
           disp[0] = 0.;
            disp[1] = 0.;
        }
        else{
            
            TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
            TVar gamma= pow(r,alfa)/(2.*mu);
            
            TVar cosT = FADcos(alfa*theta) - FADcos((alfa-2.)*theta);
            TVar sinT = A*FADsin(alfa*theta) + FADsin((alfa-2.)*theta);
            
            disp[0] =gamma*(cosT);
            disp[1] =gamma*(sinT);
        }
        
        
    }else if (fProblemType==ECrack){
        
        TVar Kolosov = 3. - 4.*gPoisson;
        TVar mu = gE / (2.*(1.+gPoisson));
        TVar theta = FADatan2(x[1],x[0]);
        
        auto thetaval = shapeFAD::val(theta);
        if (thetaval < (0.)) theta += 2. * M_PI;
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        TVar a = 1.;
        TVar sinfty = 1.;
        TVar tinfty = 1.;
        TVar K1 = sinfty * FADsqrt(M_PI*a);//Modo de abertura
        TVar K2 = tinfty * FADsqrt(M_PI*a);//Modo de deslizamento


        TVar cosT2 = FADcos(0.5*theta);
        TVar sinT2 = FADsin(0.5*theta);

        TVar mult1 = FADsqrt(r) / (2. * mu * FADsqrt(2.*M_PI));
        TVar mult2 = Kolosov - 1. + 2.*sinT2*sinT2;
        TVar mult3 = Kolosov + 1. + 2.*cosT2*cosT2;
        TVar mult4 = Kolosov + 1. - 2.*cosT2*cosT2;
        TVar mult5 = Kolosov - 1. - 2.*sinT2*sinT2;
        
        disp[0] = K1*mult1*cosT2*mult2 + K2*mult1*sinT2*mult3;
        disp[1] = K1*mult1*sinT2*mult3 + K2*mult1*cosT2*mult4;
        
        
    } else if (fProblemType==EHarmonic){        
        disp[0] = FADsin(2.*M_PI*x[1])*(-1.+FADcos(2.*M_PI*x[0])) + FADsin(M_PI*x[0])*FADsin(M_PI*x[1]);;
        disp[1] = FADsin(2.*M_PI*x[0])*(1.-FADcos(2.*M_PI*x[1])) + FADsin(M_PI*x[0])*FADsin(M_PI*x[1]);;
    }
    else
    {
        DebugStop();
    }
}

template<typename TVar1, typename TVar2>
void TElasticity2DAnalytic::uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &disp) const {
    
    if (fProblemType == Etest1) {
        disp[0] = TVar2(1. / 27.) * x[0] * x[0] * x[1] * x[1] * cos(TVar2(6. * M_PI) * x[0]) * sin(TVar2(7. * M_PI) * x[1]);
        disp[1] = TVar2(0.2) * exp(x[1]) * sin(TVar2(4. * M_PI) * x[0]);
    } else if (fProblemType == Etest2) {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((TVar2(1) - x[0] * x[0])*(1. + x[1] * x[1] * x[1] * x[1]));
        disp[1] = ((TVar2(1) - x[1] * x[1])*(1. + x[0] * x[0] * x[0] * x[0]));
    }
    else if (fProblemType == ERot)//rotation
    {
        disp[0] = (TVar2) - x[1];
        disp[1] = (TVar2) x[0];
    }
    else if (fProblemType == EShear)//pure shear
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar2) x[1];
        disp[1] += (TVar2) 0.;
    } else if (fProblemType == EStretchx)//strech x
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar2) x[0];
        disp[1] += (TVar2) 0.;
    } else if (fProblemType == EUniAxialx) {
        if (fPlaneStress == 0) {
            disp[0] = x[0]*(1. - gPoisson * gPoisson) / gE;
            disp[1] = -x[1]*(1. + gPoisson) * gPoisson / gE;
        } else {
            disp[0] = x[0] / gE;
            disp[1] = -x[1] * gPoisson / gE;
        }
    } else if (fProblemType == EStretchy)//strech y
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar2) 0.;
        disp[1] += (TVar2) x[1];
    } else if (fProblemType == EDispx) {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += 1.;
        disp[1] += 0.;
    } else if (fProblemType == EDispy) {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar2) 0.;
        disp[1] += (TVar2) 1.;
    } else if (fProblemType == EThiago){
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;


        //disp[0] = x[0]*x[1]*(x[0]-1)*(x[1]-1);
        //disp[1] = x[0]*x[1]*(x[0]-1)*(x[1]-1);
        disp[0] = (TVar2) cos(M_PI * x[0])*(TVar2) sin(2 * M_PI * x[1]);
        disp[1] = (TVar2) cos(M_PI * x[1])*(TVar2) sin(M_PI * x[0]);
        

        
    } 
    else if (fProblemType==EHomoDir){
        disp[0] =  (TVar2) sin(2 * M_PI * x[1])*(-1.+ (TVar2) cos(2*M_PI*x[0]))+ (TVar2) sin( M_PI * x[1])* (TVar2) sin( M_PI * x[0]);
        disp[1] = (TVar2) sin(2 * M_PI * x[0])*(1.-(TVar2) cos(2*M_PI*x[1]))+(TVar2) sin( M_PI * x[1])*(TVar2) sin( M_PI * x[0]);
        
    }
    
    
    else if(fProblemType == EPoly){
        disp[0] = 0. * x[0]; //*x[0]*x[1];
        disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
    } else if (fProblemType == EBend) {
        TVar2 poiss = gPoisson;
        TVar2 elast = gE;
        if (fPlaneStress == 0) {
            poiss = poiss / (1. - poiss);
            elast /= (1 - gPoisson * gPoisson);
        }
        disp[0] = 5. * x[0] * x[1] / elast;
        disp[1] = (-poiss * 5. * x[1] * x[1] / 2. - 5. * x[0] * x[0] / 2.) / elast;
    } else if (fProblemType == ELoadedBeam) {
        TVar2 Est, nust, G;
        REAL MI = 5, h = 1.;
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
        } else {
            Est = gE;
            nust = gPoisson;
        }
        disp[0] = MI * h * h * x[1] / (2. * G) + MI * x[0] * x[0] * x[1] / (2. * Est) - MI * x[1] * x[1] * x[1] / (6. * G) + MI * nust * x[1] * x[1] * x[1] / (6. * Est);
        disp[1] = -MI * x[0] * x[0] * x[0] / (6. * Est) - MI * nust * x[0] * x[1] * x[1] / (2. * Est);
    } else if (fProblemType == ESquareRoot) {
#ifdef STATE_COMPLEX
        DebugStop();
#else
        TVar2 Est, nust, G, kappa;
        TVar2 theta = atan2(x[1], x[0]);
        TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
        }
        TVar2 costh = cos(theta / 2.);
        TVar2 sinth = sin(theta / 2.);
        disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
        disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
        //        std::cout << "SQ x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
    } else if (fProblemType == ESquareRootLower) {
#ifdef STATE_COMPLEX
        DebugStop();
#else
        TVar2 Est, nust, G, kappa;
        TVar2 theta = atan2(x[1], x[0]);
        if (shapeFAD::val(theta) > 0.) {
            theta -= (2. * M_PI);
        }
        TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
        }
        TVar2 costh = cos(theta / 2.);
        TVar2 sinth = sin(theta / 2.);
        disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
        disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
        //        std::cout << "SQL x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
    } else if (fProblemType == ESquareRootUpper) {
#ifdef STATE_COMPLEX
        DebugStop();
#else
        TVar2 Est, nust, G, kappa;
        TVar2 theta = atan2(x[1], x[0]);
        if (shapeFAD::val(theta) < 0.) {
            theta += (2. * M_PI);
        }
        TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
        }
        TVar2 costh = cos(theta / 2.);
        TVar2 sinth = sin(theta / 2.);
        disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
        disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
        //        std::cout << "SQU x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
    } 
    
    else if (fProblemType==ELShape){
        
        TVar2 mu=1.;
        TVar2 lambda=5.;
        TVar2 alfa=0.6;
        TVar2 A= 31./9.;
        TVar2 theta = atan2(x[1],x[0]);
        TVar2 eps = 1e-15;
        
        auto thetaval = shapeFAD::val(theta);
        if (thetaval < (0.)) theta += 2. * M_PI;
        
        // Verification to avoid numerical errors when x > 0 and y = 0
        if (x[0] > 0. && x[1] < eps && x[1] >  1e-15) {
           disp[0] = 0.;
           disp[1] = 0.;
        }
        
        
        else{
            
            TVar2 r = sqrt(x[0]*x[0]+x[1]*x[1]);
            TVar2 gamma= pow(r,alfa)/(2.*mu);
            
            TVar2 cosT = cos(alfa*theta) - cos((alfa-2.)*theta);
            TVar2 sinT = A*sin(alfa*theta) + sin((alfa-2.)*theta);
            
            disp[0] = gamma*(cosT);
            disp[1] =gamma*(sinT);
            
        }
        
        
        

        
        
    } else if (fProblemType==ECrack){
        
        TVar2 Kolosov = 3. - 4.*gPoisson;
        TVar2 mu = gE / (2.*(1.+gPoisson));
        TVar2 theta = atan2(x[1],x[0]);
        
        auto thetaval = shapeFAD::val(theta);
        if (thetaval < (0.)) theta += 2. * M_PI;
        TVar2 r = sqrt(x[0]*x[0]+x[1]*x[1]);
        TVar2 a = 1.;
        TVar2 sinfty = 1.;
        TVar2 tinfty = 1.;
        TVar2 K1 = sinfty * sqrt(M_PI*a);//Modo de abertura
        TVar2 K2 = tinfty * sqrt(M_PI*a);//Modo de deslizamento


        TVar2 cosT2 = cos(0.5*theta);
        TVar2 sinT2 = sin(0.5*theta);

        TVar2 mult1 = sqrt(r) / (2. * mu * sqrt(2.*M_PI));
        TVar2 mult2 = Kolosov - 1. + 2.*sinT2*sinT2;
        TVar2 mult3 = Kolosov + 1. + 2.*cosT2*cosT2;
        TVar2 mult4 = Kolosov + 1. - 2.*cosT2*cosT2;
        TVar2 mult5 = Kolosov - 1. - 2.*sinT2*sinT2;
        
        disp[0] = K1*mult1*cosT2*mult2 + K2*mult1*sinT2*mult3;
        disp[1] = K1*mult1*sinT2*mult3 + K2*mult1*cosT2*mult4;
        
        
    }
    else if (fProblemType==EHarmonic){        
        disp[0] = sin(2.*M_PI*x[1])*(-1.+cos(2.*M_PI*x[0])) + sin(M_PI*x[0])*sin(M_PI*x[1]);
        disp[1] = sin(2.*M_PI*x[0])*(1.-cos(2.*M_PI*x[1])) + sin(M_PI*x[0])*sin(M_PI*x[1]);
    }
    else {
        DebugStop();
    }
}

template
void TElasticity2DAnalytic::uxy<REAL,REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &disp) const;

#if (REAL != STATE)

template
void TElasticity2DAnalytic::uxy<REAL,STATE>(const TPZVec<REAL> &x, TPZVec<STATE> &disp) const;

template
void TElasticity2DAnalytic::uxy<STATE,STATE>(const TPZVec<STATE> &x, TPZVec<STATE> &disp) const;

#endif

template<class TVar>
void TElasticity2DAnalytic::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu)
{
    if(gOscilatoryElasticity != 0)
    {
        Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    }
    else
    {
        Elast = gE;
    }
    nu = TVar(gPoisson);
}

//template<>
//void TElasticity2DAnalytic::Elastic(const TPZVec<double> &x, double &Elast, double &nu)
//{
//  Elast = 1000.;
//    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
//    Elast = 1.;
//    nu = 0.;
//}

void TElasticity2DAnalytic::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    //TPZManVector<STATE> xstate(x.size());
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    
    STATE E = 1,nu = 0.3;
    Elastic(xstate,E,nu);
    result[0] = E;
    result[1] = nu;
}


template<>
void TElasticity2DAnalytic::graduxy(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &grad) const
{
    TPZManVector<Fad<Fad<STATE> >,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<Fad<STATE> > temp = Fad<Fad<STATE> >(3,Fad<STATE>(3,0.));
//      Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        temp.val()= x[i];
        Fad<STATE> temp3(3,0.);
        for(int j=0; j<3; j++)
        {
            temp.fastAccessDx(j) = temp3;
        }
        Fad<STATE> temp2(3,1.);
        temp.fastAccessDx(i) = temp2;
//      Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        xfad[i] = temp;    
//      xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<STATE> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(j,i) = result[i].d(j);
        }
    }
}

template
void TElasticity2DAnalytic::graduxy<REAL,REAL>(const TPZVec<REAL> &x, TPZFMatrix<REAL> &grad) const;

#if (REAL != STATE)

template
void TElasticity2DAnalytic::graduxy<REAL,STATE>(const TPZVec<REAL> &x, TPZFMatrix<STATE> &grad) const;

template
void TElasticity2DAnalytic::graduxy<STATE,STATE>(const TPZVec<STATE> &x, TPZFMatrix<STATE> &grad) const;

#endif

void TElasticity2DAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const {
    TPZManVector<STATE> xst(3);
    for(int i=0; i<3; i++) xst[i] = x[i];
    uxy(xst,u);
    graduxy(xst,gradu);
}

void TElasticity2DAnalytic::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<STATE>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<STATE>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,2);
    u[0] = result[0].val();
    u[1] = result[1].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

void TElasticity2DAnalytic::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
{
    TPZFNMatrix<4,STATE> grad;
    REAL E, nu;
    sigma.Resize(2,2);
    Elastic(x, E, nu);
    graduxy(x,grad);
    //uxy(x,u);
    if (fPlaneStress == 0)
    {
        STATE Fac = E/((STATE)1.+nu)/((STATE(1.)-STATE(2.)*nu));
        sigma(0,0) = Fac*((STATE(1.)-nu)*grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*((STATE(1.)-nu)*grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(STATE(2.)*(STATE(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
    else
    {
        STATE Fac = E/((STATE)1.-nu*nu);
        sigma(0,0) = Fac*(grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*(grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(STATE(2.)*(STATE(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
}

void TElasticity2DAnalytic::Sigma(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &sigma) const {
    TPZFNMatrix<4,Fad<STATE> > grad;
    sigma.Resize(2,2);
    Fad<STATE>  E, nu;
    Elastic(x, E, nu);
    graduxy(x,grad);
    if (fPlaneStress == 0)
    {
        Fad<STATE>  Fac = E/(Fad<STATE>(1.)+nu)/((Fad<STATE>(1.)-Fad<STATE>(2.)*nu));
        sigma(0,0) = Fac*((Fad<STATE>(1.)-nu)*grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*((Fad<STATE>(1.)-nu)*grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(Fad<STATE>(2.)*(Fad<STATE>(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
    else
    {
        typedef Fad<STATE> TVar;
        TVar Fac = E/((TVar)1.-nu*nu);
        sigma(0,0) = Fac*(grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*(grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);

    }
}

template<class TVar>
//void TElasticity2DAnalytic::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma) const
void TElasticity2DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<4, Fad<TVar> > sigma(2,2);
    Sigma(xfad,sigma);
    divsigma[0] = sigma(0,0).dx(0)+sigma(0,1).dx(1);
    divsigma[1] = sigma(1,0).dx(0)+sigma(1,1).dx(1);
}



template
void TElasticity2DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<STATE> &divsigma) const;


TElasticity3DAnalytic::TElasticity3DAnalytic(const TElasticity3DAnalytic &cp) : TPZAnalyticSolution(cp),fProblemType(cp.fProblemType),fE(cp.fE),fPoisson(cp.fPoisson) {
    std::cout << "TElasticity3DAnalytic::TElasticity3DAnalytic(const TElasticity3DAnalytic &cp): One should not invoke this copy constructor";
    DebugStop();
}

TElasticity3DAnalytic & TElasticity3DAnalytic::operator=(const TElasticity3DAnalytic &copy){
    fProblemType = copy.fProblemType;
    fE = copy.fE;
    fPoisson = copy.fPoisson;
    return *this;
}

template<class TVar>
void TElasticity3DAnalytic::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    
    if(fProblemType == Etest1)
    {
        disp[0] = TVar(1./27.)*x[0]*x[0]*x[1]*x[1]*cos(TVar(6.*M_PI)*x[0])*sin(TVar(7.*M_PI)*x[1]);
        disp[1] = TVar(0.2)*exp(x[1])*sin(TVar(4.*M_PI)*x[0]);
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1.-x[0]*x[0])*(1.+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1.-x[1]*x[1])*(1.+x[0]*x[0]*x[0]*x[0]));
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType ==ERot)//rotation
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = (TVar)-x[1];
        disp[1] = (TVar)x[0];
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType == EShear)//pure shear
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) x[1];
        disp[1] += (TVar) 0. ;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EStretchx)//strech x
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) x[0];
        disp[1] += (TVar) 0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EUniAxialx)//strech x
    {
        disp[0] += (TVar) x[0]/fE;
        disp[1] += (TVar) -x[1]*fPoisson/fE;
        disp[2] = -x[2]*fPoisson/fE;
    }
    else if(fProblemType ==EStretchy)//strech y
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) 0.;
        disp[1] += (TVar) x[1];
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType ==EStretchz)//strech z
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[2] = x[0]*0.;
        disp[0] += (TVar) 0.;
        disp[1] += x[0]*TVar(0.);
        disp[2] += (TVar) x[2];
    }
    else if(fProblemType==EDispx)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[0] +=   0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispy)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) 0.;
        disp[1] += (TVar) 1.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EBend)
    {
        TVar poiss = fPoisson;
        TVar elast = fE;
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        disp[2] = x[0]*0.;
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est = fE,nust = fPoisson,G;
        REAL MI = 5, h = 1.;
        G = fE/(2.*(1.+fPoisson));
        // returning the solution of plane strain
        Est = fE/((1.-fPoisson*fPoisson));
        nust = fPoisson/(1-fPoisson);
        int offset = 2;
        int i0 = offset;
        int i1 = (offset+1)%3;
        int i2 = (offset+2)%3;
        
        disp[i0] = MI*h*h*x[i1]/(2.*G)+MI * x[i0]*x[i0]*x[i1]/(2.*Est)-MI *x[i1]*x[i1]*x[i1]/(6.*G)+MI*nust*x[i1]*x[i1]*x[i1]/(6.*Est);
        disp[i1] = -MI*x[i0]*x[i0]*x[i0]/(6.*Est)-MI*nust*x[i0]*x[i1]*x[i1]/(2.*Est);
        disp[i2] = x[i2]*0.;
    }
    else if(fProblemType == ETestShearMoment)
    {
        disp[0] = -((TVar) FI/fE*fPoisson) *x[0]*x[1]*x[2];
        disp[1] = ((TVar) FI/fE) *((TVar)(fPoisson/2.)*(x[0]*x[0]-x[1]*x[1])*x[2]-((TVar)(1./6.)*x[2]*x[2]*x[2]));
        TVar disp2A = FI/fE*(1./2.* x[1]*(fPoisson*x[0]*x[0]+x[2]*x[2])+1./6.*fPoisson*x[1]*x[1]*x[1]+(1.+fPoisson)*(b*b*x[1]-1./3.*x[1]*x[1]*x[1])-1./3.*a*a*fPoisson*x[1]);
        TVar series = (TVar) 0.;
        TVar minusone = (TVar) -1;
        for (int i=1; i<5; i++) {
            series += (minusone/TVar(i*i*i)*cos(i*M_PI*x[0]/a)*sinh(i*M_PI*x[1]/a)/cosh(i*M_PI*b/a));
            minusone *= (TVar) -1.;
        }
        series *= (-4.*a*a*a*fPoisson/(M_PI*M_PI*M_PI));
        disp[2] = disp2A+series;
    }
    else if(fProblemType == ESphere)
    {
        TPZManVector<TVar,3> xc(x);
        TVar radius2 = xc[0]*xc[0]+xc[1]*xc[1]+xc[2]*xc[2];
        TVar radius = sqrt(radius2);
        TVar acube = 10.*10.*10.;
        TVar bcube = 50.*50.*50.;
        TVar pi = 6.;
        TVar ur = pi*acube/(2.*fE*(bcube-acube)*radius2)*(2.*(1.-2.*fPoisson)*radius2*radius+(1.+fPoisson)*bcube);
        TVar cosphi = xc[2]/radius;
        TVar xyradius = sqrt(xc[0]*xc[0]+xc[1]*xc[1]);
        TVar sinphi = xyradius/radius;
        TVar costheta = xc[0]/xyradius;
        TVar sintheta = xc[1]/xyradius;
        disp[2] = ur*cosphi;
        disp[0] = ur*sinphi*costheta;
        disp[1] = ur*sinphi*sintheta;
    }
    else if(fProblemType == EYotov)
	{
        TPZManVector<TVar,3> xc(x);
		const TVar x1 = x[0], x2 = x[1], x3 = x[2];
		const TVar Pi = M_PI;

        disp[0] = -0.1*(-1. + exp(x1))*sin(Pi*x1)*sin(Pi*x2);
        disp[1] = (1. - exp(x1))*(-0.5 - ((1. + sqrt(3))*(-0.5 + x2))/(2.*sqrt(2.)) + x2 + ((-1 + sqrt(3.))*(-0.5 + x3))/(2.*sqrt(2)));
        disp[2] = (1. - exp(x1))*(-0.5 - ((-1. + sqrt(3))*(-0.5 + x2))/(2.*sqrt(2.)) - ((1 + sqrt(3.))*(-0.5 + x3))/(2.*sqrt(2)) + x3);
	}
    else{
        DebugStop();
    }
}


template<>
void TElasticity3DAnalytic::uxy(const TPZVec<FADFADSTATE > &x, TPZVec<FADFADSTATE > &disp) const
{
    typedef FADFADSTATE TVar;
    if(fProblemType == Etest1)
    {
        FADFADSTATE tmp = (FADFADSTATE)(1./27.)*x[0]*x[0]*x[1]*x[1];
        disp[0] = tmp*FADcos((FADFADSTATE)(6.*M_PI)*x[0])*FADsin((FADFADSTATE)(7.*M_PI)*x[1]);
        disp[1] = (FADFADSTATE)(0.2)*FADexp(x[1])*FADsin((FADFADSTATE)(4.*M_PI)*x[0]);
        disp[2] = x[0]*TVar(0.);

    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1.-x[0]*x[0])*(1.+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1.-x[1]*x[1])*(1.+x[0]*x[0]*x[0]*x[0]));
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType ==ERot)//rotation
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] =(FADFADSTATE)-x[1];
        disp[1] =(FADFADSTATE) x[0];
        disp[2] = x[0]*TVar(0.);

    }
    
    else if(fProblemType == EShear)//pure shear
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) x[1];
        disp[1] += (FADFADSTATE) 0. ;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EStretchx)//strech x
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) x[0];
        disp[1] += (FADFADSTATE) 0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EUniAxialx)//strech x
    {
        disp[0] += (TVar) x[0]/fE;
        disp[1] += (TVar) -x[1]*fPoisson/fE;
        disp[2] = -x[2]*fPoisson/fE;
    }
    else if(fProblemType ==EStretchy)//strech y
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) 0.;
        disp[1] += (FADFADSTATE) x[1];
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType ==EStretchz)//strech z
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[2] = x[0]*0.;
        disp[0] += (FADFADSTATE) 0.;
        disp[1] += (FADFADSTATE) 0.;
        disp[2] += (FADFADSTATE) x[2];
    }
    else if(fProblemType==EDispx)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[1] +=   0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispy)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADSTATE) 0.;
        disp[1] += (FADFADSTATE) 1.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EBend)
    {
        TVar poiss = fPoisson;
        TVar elast = fE;
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        disp[2] = x[0]*0.;
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est = fE,nust = fPoisson,G;
        REAL MI = 5, h = 1.;
        G = fE/(2.*(1.+fPoisson));
        // returning the solution of plane strain
        Est = fE/((1.-fPoisson*fPoisson));
        nust = fPoisson/(1.-fPoisson);
        int offset = 2;
        int i0 = offset;
        int i1 = (offset+1)%3;
        int i2 = (offset+2)%3;

        disp[i0] = MI*h*h*x[i1]/(2.*G)+MI * x[i0]*x[i0]*x[i1]/(2.*Est)-MI *x[i1]*x[i1]*x[i1]/(6.*G)+MI*nust*x[i1]*x[i1]*x[i1]/(6.*Est);
        disp[i1] = -MI*x[i0]*x[i0]*x[i0]/(6.*Est)-MI*nust*x[i0]*x[i1]*x[i1]/(2.*Est);
        disp[i2] = x[i2]*0.;
    }
    else if(fProblemType == ETestShearMoment)
    {
        disp[0] = -((TVar) FI/fE*fPoisson) *x[0]*x[1]*x[2];
        disp[1] = ((TVar) FI/fE) *((TVar)(fPoisson/2.)*(x[0]*x[0]-x[1]*x[1])*x[2]-((TVar)(1./6.)*x[2]*x[2]*x[2]));
        TVar disp2A = FI/fE*(1./2.* x[1]*(fPoisson*x[0]*x[0]+x[2]*x[2])+1./6.*fPoisson*x[1]*x[1]*x[1]+(1.+fPoisson)*(b*b*x[1]-1./3.*x[1]*x[1]*x[1])-1./3.*a*a*fPoisson*x[1]);
        TVar series = (TVar) 0.;
        TVar minusone = (TVar) -1;
        for (int i=1; i<5; i++) {
            TVar ivar(i);
            series += (minusone/(TVar)(i*i*i)*FADcos(ivar*M_PI*x[0]/a)*FADsinh(ivar*M_PI*x[1]/a)/FADcosh(ivar*M_PI*b/a));
            minusone *= (TVar) -1.;
        }
        series *= (-4.*a*a*a*fPoisson/(M_PI*M_PI*M_PI));
        disp[2] = disp2A+series;
    }
    else if(fProblemType == ESphere)
    {
        TPZManVector<TVar,3> xc(x);
        TVar radius2 = xc[0]*xc[0]+xc[1]*xc[1]+xc[2]*xc[2];
        TVar radius = FADsqrt(radius2);
        TVar acube = 10.*10.*10.;
        TVar bcube = 50.*50.*50.;
        TVar pi = 6.;
        TVar ur = pi*acube/(2.*fE*(bcube-acube)*radius2)*(2.*(1.-2.*fPoisson)*radius2*radius+(1.+fPoisson)*bcube);
        TVar cosphi = xc[2]/radius;
        TVar xyradius = FADsqrt(xc[0]*xc[0]+xc[1]*xc[1]);
        TVar sinphi = xyradius/radius;
        TVar costheta = xc[0]/xyradius;
        TVar sintheta = xc[1]/xyradius;
        disp[2] = ur*cosphi;
        disp[0] = ur*sinphi*costheta;
        disp[1] = ur*sinphi*sintheta;
    }
    else if(fProblemType == EYotov)
	{
		TPZManVector<TVar,3> xc(x);
		const TVar x1 = x[0], x2 = x[1], x3 = x[2];
		const TVar Pi = M_PI;

		disp[0] = -0.1*(-1. + FADexp(x1))*FADsin(Pi*x1)*FADsin(Pi*x2);
		disp[1] = (1. - FADexp(x1))*(-0.5 - ((1. + FADsqrt(3))*(-0.5 + x2))/(2.*FADsqrt(2.)) + x2 + ((-1 + FADsqrt(3.))*(-0.5 + x3))/(2.*FADsqrt(2)));
		disp[2] = (1. - FADexp(x1))*(-0.5 - ((-1. + FADsqrt(3))*(-0.5 + x2))/(2.*FADsqrt(2.)) - ((1 + FADsqrt(3.))*(-0.5 + x3))/(2.*FADsqrt(2)) + x3);
	}
    else{
        DebugStop();
    }
}


template<class TVar>
void TElasticity3DAnalytic::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu) const
{
    //    Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    Elast = TVar(fE);
    nu = TVar(fPoisson);
}

//template<>
//void TElasticity3DAnalytic::Elastic(const TPZVec<double> &x, double &Elast, double &nu)
//{
//  Elast = 1000.;
//    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
//    Elast = 1.;
//    nu = 0.;
//}

void TElasticity3DAnalytic::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    STATE E = 1.,nu = 0.3;
    DebugStop();
//    Elastic(xstate,E,nu);
    result[0] = E;
    result[1] = nu;
}


template<class TVar>
void TElasticity3DAnalytic::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad) const
{
    int sz = x.size();
#ifdef PZDEBUG
    if(sz != 3)
    {
        DebugStop();
    }
#endif
    TPZManVector<Fad<TVar>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        Fad<TVar> temp = Fad<TVar>(sz,i,x[i]);
        xfad[i] = temp;
        
    }
    TPZManVector<Fad<TVar>,3> result(3);
    uxy(xfad,result);
    grad.Resize(sz,sz);
    for (int i=0; i<sz; i++) {
        for (int j=0; j<sz; j++)
        {
            grad(j,i) = result[i].d(j);
        }
    }
    
}

void TElasticity3DAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    int sz = x.size();
    TPZManVector<Fad<STATE>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        Fad<STATE> temp = Fad<STATE>(sz,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<STATE>,3> result(3);
    uxy(xfad,result);
    gradu.Redim(sz,sz);
    for (int i=0; i<sz; i++) {
        u[i] = result[i].val();
        for (int j=0; j<sz; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

/*
template<>
void TElasticity3DAnalytic::graduxy(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &grad)
{
    TPZManVector<Fad<Fad<STATE> >,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<Fad<STATE> > temp = Fad<Fad<STATE> >(3,Fad<STATE>(3,0.));
        //      Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        temp.val()= x[i];
        Fad<STATE> temp3(3,0.);
        for(int j=0; j<3; j++)
        {
            temp.fastAccessDx(j) = temp3;
        }
        Fad<STATE> temp2(3,1.);
        temp.fastAccessDx(i) = temp2;
        //      Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        xfad[i] = temp;
        //      xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<STATE> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
}
 */

template<class TVar>
void TElasticity3DAnalytic::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZFNMatrix<9,TVar> grad;
    TVar E, nu;
    Elastic(x, E, nu);
    TVar Fac = E/((TVar)1.+nu)/((TVar(1.)-TVar(2.)*nu));
    graduxy(x,grad);
    //uxy(x,u);
    sigma.Resize(3,3);
    sigma(0,0) = Fac*((TVar(1.)-nu)*grad(0,0)+nu*grad(1,1)+nu*grad(2,2));
    sigma(1,1) = Fac*((TVar(1.)-nu)*grad(1,1)+nu*grad(0,0)+nu*grad(2,2));
    sigma(2,2) = Fac*((TVar(1.)-nu)*grad(2,2)+nu*grad(0,0)+nu*grad(1,1));
    sigma(0,1) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    sigma(0,2) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,2)+grad(2,0));
    sigma(2,0) = sigma(0,2);
    sigma(1,2) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(1,2)+grad(2,1));
    sigma(2,1) = sigma(1,2);
}

void TElasticity3DAnalytic::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const
{
    TPZManVector<STATE,3> xloc(3,0.);
    for (int i=0; i<3; i++) {
        xloc[i] = x[i];
    }
    Sigma<STATE>(xloc,tensor);
}

/*
template<>
void TElasticity3DAnalytic::Sigma(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &sigma)
{
    TPZFNMatrix<4,Fad<STATE> > grad;
    Fad<STATE>  E, nu;
    Elastic(x, E, nu);
    Fad<STATE>  Fac = E/(Fad<STATE>(1.)+nu)/((Fad<STATE>(1.)-Fad<STATE>(2.)*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((Fad<STATE>(1.)-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((Fad<STATE>(1.)-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(Fad<STATE>(2.)*(Fad<STATE>(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    
}
*/
template
void TElasticity3DAnalytic::Sigma<REAL>(const TPZVec<REAL> &x, TPZFMatrix<REAL> &divsigma) const;

template<class TVar>
void TElasticity3DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const
{
    int sz = x.size();
    TPZManVector<Fad<TVar>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        xfad[i] = Fad<TVar>(sz,i,x[i]);
    }
    TPZFNMatrix<9, Fad<TVar> > sigma(3,3);
    Sigma(xfad,sigma);
    for (int i=0; i<3; i++) {
        divsigma[i] = sigma(i,0).dx(0)+sigma(i,1).dx(1)+sigma(i,2).dx(2);
    }
}



template
void TElasticity3DAnalytic::DivSigma<STATE>(const TPZVec<REAL> &x, TPZVec<STATE> &divsigma) const;
template
void TElasticity3DAnalytic::Sigma<Fad<STATE> >(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &sigma) const;

double TLaplaceExample1::gC = 1.;

TLaplaceExample1::TLaplaceExample1(const TLaplaceExample1 &cp) : TPZAnalyticSolution(cp),fExact(cp.fExact) {
    std::cout << "TLaplaceExample1::TLaplaceExample1(const TLaplaceExample1 &cp): One should not invoke this copy constructor";
    DebugStop();
}

TLaplaceExample1 & TLaplaceExample1::operator=(const TLaplaceExample1 &copy){
    fExact = copy.fExact;
    return *this;
}

template<class TVar>
void TLaplaceExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    TPZManVector<TVar,3> xloc(x);
    for (int i = 0; i<xloc.size(); i++) {
        xloc[i] -= fCenter[i];
    }
    TVar r2 = 0.;
    for (int i=0; i<fDimension; i++) {
        r2 += xloc[i]*xloc[i];
    }
    TVar r = sqrt(r2);
    disp[0] = xloc[0]*(TVar)(0.);
    
    REAL xval[2] = {shapeFAD::val(x[0]),shapeFAD::val(x[1])};
     
    switch (fExact) {
        case EConst:
            disp[0] += gC;
            break;
        case EX:
            disp[0] += xloc[0];
            break;
        case ESinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]);
        }
            break;
        case E10SinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]*10.);
        }
            break;
            
            case E2SinSin:
            {
                disp[0] += (TVar)(1.);
                for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]*2.);
            }
                break;
        case ESinDist:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]*1.1);
        }
            break;
        case ECosCos:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= cos((TVar)M_PI*xloc[i]/2);
        }
            break;
        case EArcTan://(1+0.3sin(10Pi x))*(1+0.5cos(10Pi r)*arctan(100*(r-0.5))
        {
            TVar atanco = (r-(TVar)0.5)*100.;
            TVar freq = 10.;
            TVar mult = (TVar(1)+TVar(0.3)*sin(TVar(M_PI)*xloc[0]*freq))*(TVar(1)+TVar(0.5)*cos(TVar(M_PI)*r*freq));
            disp[0] = atan(atanco)*mult;
        }
            break;
        case EArcTanSingular://5*(0.5*pi + arctang(20(0.25 - r^2))))
        {
            REAL B = 5.;
            if(fDimension==1)
                B *= 0.25;
            else if(fDimension==3)
                B *= 4;
            // Argument value (arc) to compute ArcTangent( arc )
            TVar RCircle = 0.5;
            TVar Force = 20.;
            TVar arc = Force*(RCircle*RCircle-r2);
            TVar Prod = 1.;
            for (int i=0; i<fDimension; i++) {
                Prod *= x[i]*(1.-x[i]);
            }
            TVar temp = 0.5*M_PI + atan(arc);
            disp[0] = B*temp;
        }
            break;

        case ESteepWave:
        {
            TPZManVector<TVar,3> x0(fDimension,TVar(0.));
            TVar r2Circle = TVar(0);
            TVar r0 = TVar(0.7);
            TVar alpha = TVar(100);

            for (int i=0; i<fDimension; i++) {
//                x0[i] = TVar(-0.05);
                r2Circle += (xloc[i]-x0[i])*(xloc[i]-x0[i]);
            }
            TVar rCircle = sqrt(r2Circle);
            disp[0] = atan(alpha*(rCircle-r0));
        }
            break;
            
            //----
        case ESinMark://(r^(2/3)-r^2)sin(20/3) para homogeneo dirichlet e r^(2/3)sin(20/3) para f=0
        {

            TVar theta = atan2(xloc[1], xloc[0]);//theta=arctan(y/x)
            auto thetaval = shapeFAD::val(theta);
            if (thetaval < (0.)) theta += 2. * M_PI;

            // Verification to avoid numerical errors when x > 0 and y = 0
            if (xval[0] > 0 && xval[1] < 1e-15 && xval[1] > -1.e-15) {
               disp[0] = 0.;
            }
            else {
                TVar factor = pow(r,TVar (2.)/TVar (3.));//pow(r,TVar (2.)/TVar (3.))-pow(r,TVar (2.));//
                disp[0] = factor * (sin((TVar) (2.) * theta / TVar(3.)));
            }


        }
            break;
            
            //--
     
        case ESinSinDirNonHom: //sin(2 pi x)sin(2pi y)+1/(x+y+1)
        {
            
            disp[0]=(sin((TVar)M_PI*2.*xloc[0]))*sin((TVar)M_PI*2.*xloc[1])+(TVar)1./(xloc[0]+xloc[1]+(TVar)(1.));
            
        }
            
            break;
            
        case ESteklovNonConst://Steklov function for eigenvalue lambda=0.53544094560246 and permeability Omega1=Omega=3, Omega2=Omega4=5
        {
            TVar coefs[] = {1., 0.44721359549995787, 2.3333333333333326,
                -0.7453559924999296, 0.5555555555555556,
                -0.9441175904999111, -0.48148148148148173,
                -2.4017026424997736};
            TVar lambda = 0.53544094560246;
            TVar t = atan2(xloc[1], xloc[0]);
            REAL tval = shapeFAD::val(t);
            if(tval < (0.)) t += 2.*M_PI;
            
            if((xval[0] >=(0.)) && (xval[1] >=(0.))){
               // std::cout<<"1o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                
                disp[0]=pow(r, lambda)*(TVar(coefs[0])*cos(lambda *t) + TVar(coefs[1])*sin(lambda*t) );
               // std::cout<<"valor da funcao no 1o. Q "<<disp[0]<<std::endl;
               // disp[0]=pow(r, lambda)*(cos(lambda *t)+TVar(-1.)*TVar(0.1)*sin(lambda*t));
                
            }
            
            if(( xval[0] <= (0.)) && (xval[1] >=(0.))){
               // std::cout<<"2o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[2])*cos(lambda*t) + TVar(coefs[3])*sin(lambda* t));
                
             //    std::cout<<"valor da funcao no 2o. Q "<<disp[0]<<std::endl;
            }
            
            if((xval[0] <(0.)) && ( xval[1] <= (0.))){
               // std::cout<<"3o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                disp[0]= pow(r, lambda)*(TVar(coefs[4] )*cos(lambda*t) + TVar(coefs[5])*sin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(0.882757 )*cos(lambda*t) + TVar(-1.)*TVar(0.480355)*sin(lambda* t));
              //   std::cout<<"valor da funcao no 3o. Q "<<disp[0]<<std::endl;
            }
            if(( xval[0] >= (0.)) && ( xval[1] < (0.))){
              //  std::cout<<"4o. Q "<<xloc<< " r " << r << " th " << t << std::endl;

                disp[0]= pow(r, lambda)*(TVar(coefs[6])*cos(lambda*t) +  TVar(coefs[7])*sin(lambda* t));
                
               // std::cout<<"valor da funcao no 4o. Q "<<disp[0]<<std::endl;
                
            }
            
            
        }
            break;
            
            case EGalvisNonConst:
        {
            
            TVar k1 = 2;
            TVar k2 = 5;
            
            
            if((xval[0] <= (0.)) && (xval[1] <= (0.))){
                
                TVar u1 = sin(M_PI*(xloc[0]+TVar(1.))/(k1+1.));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"3o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if(( xval[0] > (0.)) && (xval[1] < (0.))){
             
                
                TVar u1 = sin(M_PI*(k1*xloc[0]+TVar(1.))/(k1 + TVar(1.)));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"4o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if((xval[0] < (0.)) && ( xval[1] >= (0.))){
              
                TVar u1 = sin(M_PI*(xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"2o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            if(( xval[0] >= (0.)) && ( xval[1] >= (0.))){
               
                
                TVar u1 = sin(M_PI*(k1*xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"1o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
                
            }
            
            
        }
            break;
            
        case EBoundaryLayer:{
            TVar term1 = xloc[0]*xloc[1]*(1.-xloc[0])*(1.-xloc[1]);
            TVar term2 = exp(TVar(10.)*xloc[0])*exp(TVar(10.)*xloc[1]);
            TVar factor = TVar(537930);
            
            disp[0] = (term1*term2)/factor;
            
          //  std::cout<<"Pto "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
            
            
        }
            break;
            
        case EBubble:{
            
            disp[0] = xloc[0]*xloc[1]*xloc[2]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1])*(TVar(1.)-xloc[2]);
           // std::cout<<"pto "<<xloc <<" f(x) "<<disp<<std::endl;
        }
            break;
        case EBubble2D:{

            disp[0] = xloc[0]*xloc[1]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1]);
        }
            break;
        case ESinCosCircle:{
            
            TVar coef = pow(r, TVar(4.));
            TVar theta=atan2(xloc[1],xloc[0]);
            disp[0] = coef*sin(TVar(2.)*theta)*cos(TVar(2.)*theta);
            
            
        }
            break;
        case EHarmonic:
            disp[0] = exp(M_PI*x[0])*sin(M_PI*x[1]);
            break;

        case EHarmonic2:
        {
            TVar a1 = 1./4;
            TVar alpha = M_PI/2;
            disp[0] = x[0]*a1*cos(x[0]*alpha)*cosh(x[1]*alpha) + x[1]*a1*sin(x[0]*alpha)*sinh(x[1]*alpha);
        }
            break;

        case EHarmonicPoly:
        {
            disp[0] = - 3*x[0]*x[0]*x[1] - 3*x[1]*x[1]*x[0] + x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
        }
            break;
        case ESquareRoot:
        {
            TVar r = sqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = atan2(x[1],x[0]);
            disp[0] = pow(2.,1/4.)*sqrt(r)*cos(theta/2);
            //disp[0] = pow(2.,-1/4.)*sqrt(x[0] + sqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;

        case ESquareRootLower:
        {
            TVar r = sqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = atan2(x[1],x[0]);
            if (shapeFAD::val(theta) > 0.) {
                theta -= (2.*M_PI);
            }
            disp[0] = pow(2.,1/4.)*sqrt(r)*cos(theta/2);
            //disp[0] = pow(2.,-1/4.)*sqrt(x[0] + sqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;

        case ESquareRootUpper:
        {
            TVar r = sqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = atan2(x[1],x[0]);
            if (shapeFAD::val(theta) < 0.) {
                theta += (2.*M_PI);
            }
            disp[0] = pow(2.,1/4.)*sqrt(r)*cos(theta/2);
            //disp[0] = pow(2.,-1/4.)*sqrt(x[0] + sqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;

        case ELaplace2D:
        {
            TVar Ck;
            disp[0] = 1.;
            for(int k = 1 ; k < fmaxIter+1; k++){
                TVar kpi = k*M_PI;
                TVar ck = -2/tanh(kpi)*(0.1013211836423378 -0.1013211836423378*cos(kpi)-0.1591549430918954 *k* sin(kpi))/(kpi*k*k);
                disp[0]+=ck/cosh(kpi)*cos(kpi*x[0])*cosh(kpi*(x[1]-1));
            }
        }
            break;

        default:
            disp[0] = 0.;
            break;
    }
}

template<>
void TLaplaceExample1::uxy(const TPZVec<FADFADSTATE > &x, TPZVec<FADFADSTATE > &disp) const
{
#ifdef PZDEBUG
    
    if(fDimension == -1){
        DebugStop();
    }
    
#endif

    typedef FADFADSTATE TVar;
    TPZManVector<TVar,3> xloc(x);
    for (int i = 0; i<xloc.size(); i++) {
        xloc[i] -= fCenter[i];
    }
    TVar r2 = 0.;
    for (int i=0; i<fDimension; i++) {
        r2 += xloc[i]*xloc[i];
    }
    TVar r = FADsqrt(r2);
    disp[0] = (TVar)(0.)*xloc[0];
    switch (fExact) {
        case EConst:
            disp[0] += gC;
            break;
        case EX:
            disp[0] += xloc[0];
            break;
        case ESinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]);
        }
            break;
        case E10SinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]*10.);
        }
            break;
        case E2SinSin:
            {
                disp[0] += (TVar)(1.);
                for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]*2.);
            }
            break;
            
            
        case ESinDist:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]*1.1);
        }
            break;
        case ECosCos:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADcos((TVar)M_PI*xloc[i]/2);
        }
            break;
        case EArcTan:
        {
            TVar atanco = (r-(TVar)0.5)*100.;
            TVar freq = 10.;
            TVar mult = (TVar(1)+TVar(0.3)*FADsin(TVar(M_PI)*xloc[0]*freq))*(TVar(1)+TVar(0.5)*FADcos(TVar(M_PI)*r*freq));
            disp[0] = FADatan(atanco)*mult;
        }
            break;
        case EArcTanSingular:
        {
            REAL B = 5.;
            if(fDimension==1)
                B *= 0.25;
            else if(fDimension==3)
                B *= 4;
            // Argument value (arc) to compute ArcTangent( arc )
            TVar RCircle = 0.5;
            TVar Force = 20.;
            TVar arc = Force*(RCircle*RCircle-r2);
            TVar Prod = 1.;
            for (int i=0; i<fDimension; i++) {
                Prod *= x[i]*(1.-x[i]);
            }
            TVar temp = 0.5*M_PI + FADatan(arc);
            disp[0] = B*temp;
        }
            break;

        case ESteepWave:
        {
            TPZManVector<TVar,3> x0(fDimension,TVar(0.));
            TVar r2Circle = TVar(0);
            TVar r0 = TVar(0.7);
            TVar alpha = TVar(100);

            for (int i=0; i<fDimension; i++) {
//                x0[i] = TVar(-0.05);
                r2Circle += (xloc[i]-x0[i])*(xloc[i]-x0[i]);
            }
            TVar rCircle = FADsqrt(r2Circle);
            disp[0] = FADatan(alpha*(rCircle-r0));
        }
            break;

        case ESinSinDirNonHom: //sin(2pi x)sin(2pi y)+1/(x+y+1)
        {
            
            disp[0]=(FADsin((TVar)M_PI*2.*xloc[0]))*FADsin((TVar)M_PI*2.*xloc[1])+(TVar)1./(xloc[0]+xloc[1]+(TVar)(1.));
            
        }
            
            break;
        case ESinMark://(r^(2/3)-r^2)sin(20/3) para f=0 usar r^(2/3)sin(20/3)
        {
            
            TVar theta=FADatan2(xloc[1],xloc[0]);//theta=atan(y/x)
#ifdef STATE_COMPLEX
            if( theta.val().val().real() < 0.) theta += 2.*M_PI;
#else
            if( theta < TVar(0.)) theta += 2.*M_PI;
#endif
            
            // Verification to avoid numerical errors when x > 0 and y = 0
#ifdef STATE_COMPLEX
            if ((xloc[0].val().val().real() > 0.) && (xloc[1].val().val().real() <  (1e-15)) && (xloc[1].val().val().real() > (-1e-15))) {
               disp[0] = TVar(0.);
            }
#else
            if ((xloc[0] > TVar(0.)) && (xloc[1] < TVar (1e-15)) && (xloc[1] > TVar(-1e-15))) {
               disp[0] = TVar(0.);
            }
#endif
            else{
                
            TVar factor = pow(r,TVar (2.)/TVar (3.));//pow(r,TVar (2.)/TVar (3.))-pow(r,TVar (2.));
            disp[0] = factor*(FADsin((TVar)(2.)*theta/TVar(3.)));
        }
            
        }
            break;
        case ESteklovNonConst://Steklov function for eigenvalue lambda=0.126902 and permeability Omega1=Omega=3=5, Omega2=Omega4=1
        {
            
            TVar coefs[] = {1., 0.44721359549995787, 2.3333333333333326,
                -0.7453559924999296, 0.5555555555555556,
                -0.9441175904999111, -0.48148148148148173,
                -2.4017026424997736};
            TVar lambda = 0.53544094560246;
#ifdef STATE_COMPLEX
            double xr = xloc[0].val().val().real();
            double yr = xloc[1].val().val().real();
#else
            double xr = xloc[0].val().val();
            double yr = xloc[1].val().val();
#endif
            TVar t = FADatan2(xloc[1], xloc[0]);
            double tval = atan2(yr,xr);
            if(tval < (0.)) t += TVar(2.*M_PI);

            if((xr >= (0.)) && (yr >= (0.))){
              //  std::cout<<"1o. Q"<<xloc<<std::endl;
                
                disp[0]=pow(r, lambda)*(TVar(coefs[0])*FADcos(lambda *t) + TVar(coefs[1])*FADsin(lambda*t));
              //  std::cout<<"valor da funcao no 1o. Q "<<disp[0]<<std::endl;
                // disp[0]=pow(r, lambda)*(cos(lambda *t)+TVar(-1.)*TVar(0.1)*sin(lambda*t));
                
            }
            
            if(( xr < (0)) && (yr >(0.))){
               // std::cout<<"2o. Q"<<xloc<<std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[2])*FADcos(lambda*t) + TVar(coefs[3])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(2.9604)*cos(lambda*t) +TVar(-1.)* TVar(9.60396)*sin(lambda* t));
                //std::cout<<"valor da funcao no 2o. Q "<<disp[0]<<std::endl;
            }
            
            if((xr < (0.)) && ( yr <= (0.))){
             //   std::cout<<"3o. Q"<<xloc<<std::endl;
                disp[0]= pow(r, lambda)*(TVar(coefs[4] )*FADcos(lambda*t) + TVar(coefs[5])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(0.882757 )*cos(lambda*t) + TVar(-1.)*TVar(0.480355)*sin(lambda* t));
               // std::cout<<"valor da funcao no 3o. Q "<<disp[0]<<std::endl;
            }
            if(( xr >= (0.)) && ( yr < 0.)){
               // std::cout<<"4o. Q"<<xloc<<std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[6])*FADcos(lambda*t) +  TVar(coefs[7])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(6.45646)*cos(lambda*t) +  TVar(7.70156 )*sin(lambda* t));
              //  std::cout<<"valor da funcao no 4o. Q "<<disp[0]<<std::endl;
                
            }

            
            
        }
            break;
            
        case EGalvisNonConst:
        {
            
            TVar k1 = 2;
            TVar k2 = 5;
#ifdef STATE_COMPLEX
            double xr = xloc[0].val().val().real();
            double yr = xloc[1].val().val().real();
#else
            double xr = xloc[0].val().val();
            double yr = xloc[1].val().val();
#endif

            
            if((xr <= (0.)) && (yr <= (0.))){
                
                TVar u1 = FADsin(M_PI*(xloc[0]+TVar(1.))/(k1+1.));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"3o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if(( xr > (0.)) && (yr < (0.))){
                
                
                TVar u1 = FADsin(M_PI*(k1*xloc[0]+TVar(1.))/(k1 + TVar(1.)));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"4o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if((xr < (0.)) && ( yr >= (0.))){
                
                TVar u1 = FADsin(M_PI*(xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
              //  std::cout<<"2o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            if(( xr >= (0.)) && ( yr >= (0.))){
                
                
                TVar u1 = FADsin(M_PI*(k1*xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
                //std::cout<<"1o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
                
            }
        }
            break;
            
        case EBoundaryLayer:{
            TVar term1 = xloc[0]*xloc[1]*(1.-xloc[0])*(1.-xloc[1]);
            TVar term2 = FADexp(TVar(10.)*xloc[0])*FADexp(TVar(10.)*xloc[1]);
            TVar factor = TVar(537930);
            
            disp[0] = (term1*term2)/factor;

            
        }
            break;
            
        case EBubble:{
            
            disp[0] = xloc[0]*xloc[1]*xloc[2]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1])*(TVar(1.)-xloc[2]);
         //   std::cout<<"pto "<<xloc <<" f(x) "<<disp<<std::endl;
            
        }
            break;
        case EBubble2D:{

            disp[0] = xloc[0]*xloc[1]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1]);
            //   std::cout<<"pto "<<xloc <<" f(x) "<<disp<<std::endl;

        }
            break;
            
        case ESinCosCircle:{
            
            TVar coef = pow(r, TVar(4.));
             TVar theta=FADatan2(xloc[1],xloc[0]);
            disp[0] = coef*FADsin(TVar(2.)*theta)*FADcos(TVar(2.)*theta);
           
            
        }
            break;

        case EHarmonic:
            disp[0] = FADexp(M_PI*x[0])*FADsin(M_PI*x[1]);
            break;

        case EHarmonic2:
        {
            TVar a1 = 1./4;
            TVar alpha = M_PI/2;
            disp[0] = x[0]*a1*FADcos(x[0]*alpha)*FADcosh(x[1]*alpha) + x[1]*a1*FADsin(x[0]*alpha)*FADsinh(x[1]*alpha);
            break;
        }
        case EHarmonicPoly:
        {
            disp[0] = - 3*x[0]*x[0]*x[1] - 3*x[1]*x[1]*x[0] + x[0]*x[0]*x[0] + x[1]*x[1]*x[1];
            break;
        }
        case ESquareRoot:
        {
            TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = FADatan2(x[1],x[0]);
            disp[0] = pow(2.,1/4.)*FADsqrt(r)*FADcos(theta/2);
            //disp[0] = pow(2.,-1/4.)*FADsqrt(x[0] + FADsqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;

        case ESquareRootLower:
        {
            TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = FADatan2(x[1],x[0]);
            if (shapeFAD::val(theta) > 0.) {
                theta -= (2.*M_PI);
            }
            disp[0] = pow(2.,1/4.)*FADsqrt(r)*FADcos(theta/2);
            //disp[0] = pow(2.,-1/4.)*FADsqrt(x[0] + FADsqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;
            
        case ESquareRootUpper:
        {
            TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
            TVar theta = FADatan2(x[1],x[0]);
            if (shapeFAD::val(theta) < 0.) {
                theta += (2.*M_PI);
            }
            disp[0] = pow(2.,1/4.)*FADsqrt(r)*FADcos(theta/2);
            //disp[0] = pow(2.,-1/4.)*FADsqrt(x[0] + FADsqrt(x[0]*x[0] + x[1]*x[1]));
        }
            break;
        case ELaplace2D:
        {
            TVar Ck;
            disp[0] = 1.;
            for(int k = 1 ; k < fmaxIter+1; k++){
                TVar kpi = k*M_PI;
                TVar ck = -2/(kpi)*(0.1013211836423378 -0.1013211836423378*FADcos(kpi)-0.1591549430918954 *k* FADsin(kpi))/(kpi*k*k);
                disp[0]+=ck*(FADcosh(kpi)/FADsinh(kpi))*FADcos(kpi*x[0])*FADcosh(kpi*(x[1]-1));
            }
        }
            break;
            
        default:
            disp[0] = xloc[0]*0.;
            break;
    }

}

void TLaplaceExample1::setPermeabilyTensor(TPZFNMatrix<9,REAL> K, TPZFNMatrix<9,REAL> invK)
{
    fTensorPerm = K;
    fInvPerm = invK;
}

template<class TVar>
void TLaplaceExample1::Permeability(const TPZVec<TVar> &x, TVar &Perm)
{
    Perm = (TVar)(1.);
}

void TLaplaceExample1::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE,3> xloc(x.size());
    for (unsigned int i = 0; i<xloc.size(); ++i) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0,0) = Perm;
    deriv(1,1) = Perm;
    deriv(2,2) = Perm;
    deriv(3,0) = 1./Perm;
    deriv(4,1) = 1./Perm;
    deriv(5,2) = 1./Perm;
}

template<class TVar>
void TLaplaceExample1::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<TVar> temp = Fad<TVar>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(3);
    for (int i=0; i<3; i++)
    {
            grad[i] = result[0].d(i);
    }
}

template<>
void TLaplaceExample1::graduxy<Fad<STATE>>(const TPZVec<Fad<STATE>> &x, TPZVec<Fad<STATE>> &grad) const
{
    typedef Fad<STATE> TVar;
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<TVar> temp = Fad<TVar>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(3);
    for (int i=0; i<3; i++)
    {
            grad[i] = result[0].d(i);
    }
}

void TLaplaceExample1::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<STATE>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<STATE> temp = Fad<STATE>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<STATE>,3> result(1);
    uxy(xfad,result);
    gradu.Redim(3,1);
    
    
    u[0] = result[0].val();
    for (int i=0; i<3; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }

    
}

template<class TVar>
void TLaplaceExample1::SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZManVector<TVar,3> grad;
    graduxy(x,grad);
    sigma.Resize(3,1);

    //Sigma_i = 1
    for(int rowIndex =0; rowIndex < 3 ; rowIndex++) sigma(rowIndex) =0;

    //Sigma = -K*grad[u]
    TVar Perm;
    for(int rowIndex =0; rowIndex < 3 ; rowIndex++) for(int colIndex=0; colIndex < 3; colIndex++)
    {
        Perm = (TVar) fTensorPerm.GetVal(rowIndex,colIndex);
        sigma(rowIndex) -= Perm*grad[colIndex];
    }
}

template<class TVar>
void TLaplaceExample1::DivSigma(const TPZVec<REAL> &x, TVar &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        xfad[i] = Fad<TVar>(3,i,x[i]);
    }
    TPZFNMatrix<3,Fad<TVar> > sigma(3,1);
    SigmaLoc(xfad,sigma);
    divsigma = sigma(0).dx(0)+sigma(1).dx(1)+sigma(2).dx(2);
    
}

template
void TLaplaceExample1::SigmaLoc(const TPZVec<STATE> &x, TPZFMatrix<STATE> &sigma) const;

template
void TLaplaceExample1::DivSigma<STATE>(const TPZVec<REAL> &x, STATE &divsigma) const;

template
void TLaplaceExample1::graduxy<STATE>(const TPZVec<STATE> &x, TPZVec<STATE> &grad) const;



//ExactFunc *Exact();


template<class TVar>
void TLaplaceExampleTimeDependent::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    switch(fProblemType)
    {
        case ELinear:
            disp[0] = x[0];
            break;
        case ESin:
            disp[0] = sin(M_PI*x[0])*sin(M_PI*x[1])*exp(-2.*M_PI*M_PI*fTime);
            break;
        case ECos:
            disp[0] = cos(M_PI_2*x[0])*cos(M_PI_2*x[1])*exp(-M_PI*M_PI_2*fTime);
            break;
        default:
            DebugStop();
    }
}

template<>
void TLaplaceExampleTimeDependent::uxy(const TPZVec<FADFADSTATE > &x, TPZVec<FADFADSTATE > &disp) const
{
    switch(fProblemType)
    {
        case ELinear:
            disp[0] = x[0];
            break;
        case ESin:
            disp[0] = FADsin(M_PI*x[0])*FADsin(M_PI*x[1])*FADexp(-2.*M_PI*M_PI*fTime);
            break;
        case ECos:
            disp[0] = FADcos(M_PI_2*x[0])*FADcos(M_PI_2*x[1])*FADexp(-M_PI*M_PI_2*fTime);
            break;
        default:
            DebugStop();
    }

}

template<class TVar>
void TLaplaceExampleTimeDependent::Permeability(const TPZVec<TVar> &x, TVar &Perm) const
{
    Perm = (TVar)(fK);
}


template<class TVar>
void TLaplaceExampleTimeDependent::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(2);
    for (int i=0; i<2; i++)
    {
        grad[i] = result[0].d(i);
    }
}

void TLaplaceExampleTimeDependent::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<STATE>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<STATE> temp = Fad<STATE>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<STATE>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,1);
    u[0] = result[0].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

template<class TVar>
void TLaplaceExampleTimeDependent::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZManVector<TVar,3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x,grad);
    sigma.Resize(2,1);
    sigma(0) = -Perm*grad[0];
    sigma(1) = -Perm*grad[1];
    
}


void TLaplaceExampleTimeDependent::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
{
    typedef STATE TVar;
    TPZManVector<TVar,3> grad, xst(3);
    for(int i=0; i<3; i++) xst[i] = x[i];
    TVar Perm;
    Permeability(xst, Perm);
    graduxy(xst,grad);
    sigma.Resize(2,1);
    sigma(0) = -Perm*grad[0];
    sigma(1) = -Perm*grad[1];
    
}

template<class TVar>
void TLaplaceExampleTimeDependent::DivSigma(const TPZVec<REAL> &x, TVar &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<3, Fad<TVar> > sigma(2,1);
    Sigma(xfad,sigma);
    divsigma = sigma(0).dx(0)+sigma(1).dx(1);
    
}


template
void TLaplaceExampleTimeDependent::DivSigma(const TPZVec<REAL> &x, STATE &divsigma) const;


template<class TVar>
void TStokesAnalytic::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &flux) const
{
    TVar x1 = x[0];
    TVar x2 = x[1];
    TVar x3 = x[2];
    REAL Re = 0.;
    REAL lambda = 0.;
    REAL xs = 0.; 
    
    switch(fExactSol)
    {
        case ESinCos:
            flux[0] = -1.*sin(x1)*sin(x2);
            flux[1] = -1.*cos(x1)*cos(x2);
            break;
        case ENoFlow:
            flux[0] = x1*0.;
            flux[1] = x1*0.;
            break;
        case ESinCosBDS:
	    if(fvisco==0) xs = 0.;
	    xs = 1./exp(fcBrinkman/fvisco);
            flux[0] = -xs*sin(x1)*sin(x2)+(1.-xs)*sin(x1)*sin(x2);
            flux[1] = -xs*cos(x1)*cos(x2)-(1.-xs)*cos(x1)*cos(x2);
            break;
	case ECouplingSD:
	    if(x2<0.){
		flux[0] = (exp(-x2)-exp(x2))*cos(x1);
		flux[1] = -(exp(-x2)+exp(x2))*sin(x1);			
	    }else if(x2>=0.){
		flux[0] = (2./Pi)*sin(Pi*x2)*cos(Pi*x2)*cos(x1);
		flux[1] = ((1./(Pi*Pi))*sin(Pi*x2)*sin(Pi*x2)-2.)*sin(x1);	
	    }
	    break;
	case ECouplingNSD:
	    if(x2<1.){
		flux[0] = -(1./8.)*Pi*Pi*x2*sin(Pi*x1/2.)*sin(1.-x2);
		flux[1] = (1./4.)*Pi*cos(Pi*x1/2.)*(-x2*cos(1.-x2)+sin(1.-x2));			
	    }else if(x2>=1.){
	        flux[0] = cos(Pi*x2/2.)*cos(Pi*x2/2.)*sin(Pi*x1/2.);
	        flux[1] = -cos(Pi*x1/2.)*((1./4.)*sin(Pi*x2)+Pi*x2/4.);	
	    }
	    break;
	case ESinCosBDS3D:
	    xs = 1./exp(fcBrinkman/fvisco);
            flux[0] = -xs*sin(x1)*sin(x2)+(1.-xs)*sin(x1)*sin(x2);
            flux[1] = xs*(-cos(x1)*cos(x2)-sin(x2)*sin(x3))+(1.-xs)*(-cos(x1)*cos(x2)+sin(x2)*sin(x3));
            flux[2] = -xs*cos(x2)*cos(x3)+(1.-xs)*(-cos(x2)*cos(x3));
            break;
	case EGatica3D:
    	    flux[0] = cos(Pi*x1)*sin(Pi*x2)*sin(Pi*x3);
	    flux[1] = sin(Pi*x1)*cos(Pi*x2)*sin(Pi*x3);
	    flux[2] = -2.*sin(Pi*x1)*sin(Pi*x2)*cos(Pi*x3);
	break;
        case ESinCos3D:
            flux[0] = -sin(x1)*sin(x2);
            flux[1] = -cos(x1)*cos(x2)-sin(x2)*sin(x3);
            flux[2] = -cos(x2)*cos(x3);
            break;
        case EKovasznay:
        case EKovasznayCDG:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            flux[0] = 1. - exp(lambda*x1)*cos(2.*Pi*x2);
            flux[1] = (lambda/(2.*Pi))*exp(lambda*x1)*sin(2.*Pi*x2);
            break;
        case EPconst:
            flux[0] = x1;
            flux[1] = -x2;
            break;
        default:
            DebugStop();
    }
}


template<>
void TStokesAnalytic::uxy(const TPZVec<FADFADSTATE > &x, TPZVec<FADFADSTATE > &flux) const
{
    FADFADSTATE x1 = x[0];
    FADFADSTATE x2 = x[1];
    FADFADSTATE x3 = x[2];
    REAL Re = 0.;
    REAL lambda = 0.;    
    FADFADSTATE xs = 0.;
    
    switch(fExactSol)
    {
        case ESinCos:
            flux[0] = -1.*FADsin(x1)*FADsin(x2);
            flux[1] = -1.*FADcos(x1)*FADcos(x2);
            break;
        case ENoFlow:
            flux[0] = (FADFADSTATE) x1*0.;
            flux[1] = (FADFADSTATE) x1*0.;
            break;
        case ESinCosBDS:
	    if(fvisco==0) xs = 0.;
	    xs = 1./FADexp(fcBrinkman/fvisco);
            flux[0] = -xs*FADsin(x1)*FADsin(x2)+(1.-xs)*FADsin(x1)*FADsin(x2);
            flux[1] = -xs*FADcos(x1)*FADcos(x2)-(1.-xs)*FADcos(x1)*FADcos(x2);
            break;
	case ECouplingSD:
	    if(x2< (FADFADSTATE) 0.){
		flux[0] = (FADexp(-x2)-FADexp(x2))*FADcos(x1);
		flux[1] = -(FADexp(-x2)+FADexp(x2))*FADsin(x1);			
	    }else if(x2>=(FADFADSTATE) 0.){
		flux[0] = (2./Pi)*FADsin(Pi*x2)*FADcos(Pi*x2)*FADcos(x1);
		flux[1] = ((1./(Pi*Pi))*FADsin(Pi*x2)*FADsin(Pi*x2)-2.)*FADsin(x1);	
	    }
	    break;
	case ECouplingNSD:
	    if(x2< (FADFADSTATE) 1.){
		flux[0] = -(1./8.)*Pi*Pi*x2*FADsin(Pi*x1/2.)*FADsin(1.-x2);
		flux[1] = (1./4.)*Pi*FADcos(Pi*x1/2.)*(-x2*FADcos(1.-x2)+FADsin(1.-x2));			
	    }else if(x2>= (FADFADSTATE) 1.){
	        flux[0] = FADcos(Pi*x2/2.)*FADcos(Pi*x2/2.)*FADsin(Pi*x1/2.);
	        flux[1] = -FADcos(Pi*x1/2.)*((1./4.)*FADsin(Pi*x2)+Pi*x2/4.);	
	    }
	    break;
	case ESinCosBDS3D:
	    xs = 1./exp(fcBrinkman/fvisco);
            flux[0] = -xs*FADsin(x1)*FADsin(x2)+(1.-xs)*FADsin(x1)*FADsin(x2);
            flux[1] = xs*(-FADcos(x1)*FADcos(x2)-FADsin(x2)*FADsin(x3))+(1.-xs)*(-FADcos(x1)*FADcos(x2)+FADsin(x2)*FADsin(x3));
            flux[2] = -xs*FADcos(x2)*FADcos(x3)+(1.-xs)*(-FADcos(x2)*FADcos(x3));
            break;
	case EGatica3D:
    	    flux[0] = FADcos(Pi*x1)*FADsin(Pi*x2)*FADsin(Pi*x3);
	    flux[1] = FADsin(Pi*x1)*FADcos(Pi*x2)*FADsin(Pi*x3);
	    flux[2] = -2.*FADsin(Pi*x1)*FADsin(Pi*x2)*FADcos(Pi*x3);
        case ESinCos3D:
            flux[0] = -FADsin(x1)*FADsin(x2);
            flux[1] = -FADcos(x1)*FADcos(x2)-FADsin(x2)*FADsin(x3);
            flux[2] = -FADcos(x2)*FADcos(x3);
            break;
        case EKovasznay:
        case EKovasznayCDG:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            flux[0] = 1. - FADexp(lambda*x1)*FADcos(2.*Pi*x2);
            flux[1] = (lambda/(2.*Pi))*FADexp(lambda*x1)*FADsin(2.*Pi*x2);
            break;
        case EPconst:
            flux[0] = x1;
            flux[1] = -x2;
            break;
        default:
            DebugStop();
    }
    
}

template<class TVar>
void TStokesAnalytic::pressure(const TPZVec<TVar> &x, TVar &p) const
{
    TVar x1 = x[0];
    TVar x2 = x[1];
    TVar x3 = x[2];
    TPZVec<TVar> flux(3,0.);
    REAL Re = 0.;
    REAL lambda = 0.;    
    
    switch(fExactSol)
    {
        case ESinCos:
        case ESinCosBDS:
            p = cos(x1)*sin(x2);
            break;
        case ENoFlow:
            p = multRa*(x2*x2*x2-(x2*x2)/2.+x2-7./12.);
            break;
	case ECouplingSD:
	    if(x2<0.){
		p = (-exp(-x2)+exp(x2))*sin(x1);		
	    }else if(x2>=0.){
		p = sin(x1)*sin(x2);
	    }
	    break;
	case ECouplingNSD:
	    if(x2<1.){
		p = -(Pi*x2/4.)*cos(Pi*x1/2.)*sin(1.-x2);		
	    }else if(x2>=1.){
                p = (Pi/4.)*cos(Pi*x1/2.)*(x2-1.-cos(Pi*x2))*sin(x2-1.);
	    }
	    break;
        case ESinCos3D:
        case ESinCosBDS3D:
            p = cos(x1)*sin(x2)+cos(x2)*sin(x3);
            break;
        case EGatica3D:
            p = sin(Pi*x1)*sin(Pi*x2)*sin(Pi*x3);
            break;
        case EKovasznay:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            p = -(1./2.)*exp(2.*lambda*x1);
            break;
        case EKovasznayCDG:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            flux[0] = 1. - exp(lambda*x1)*cos(2.*Pi*x2);
            flux[1] = (lambda/(2.*Pi))*exp(lambda*x1)*sin(2.*Pi*x2);
            p = -(1./2.)*exp(2.*lambda*x1);
            p += (1./2.)*(flux[0]*flux[0]+flux[1]*flux[1]);
            break;
        case EPconst:
            p = 1;
            break;
        default:
            DebugStop();
    }
}

template<>
void TStokesAnalytic::pressure(const TPZVec<FADFADSTATE > &x, FADFADSTATE &p) const
{
    FADFADSTATE x1 = x[0];
    FADFADSTATE x2 = x[1];
    FADFADSTATE x3 = x[2];
    TPZVec<FADFADSTATE > flux(3,0.);
    REAL Re = 0.;
    REAL lambda = 0.;
    FADFADSTATE fadRa = multRa;    

    switch(fExactSol)
    {
        case ESinCos:
        case ESinCosBDS:
            p = FADcos(x1)*FADsin(x2);
            break;
        case ENoFlow:
            p = (FADFADSTATE) fadRa*(x2*x2*x2-(x2*x2)/2.+x2-7./12.);
            break;
	case ECouplingSD:
	    if(x2< (FADFADSTATE) 0.){
		p = (-FADexp(-x2)+FADexp(x2))*FADsin(x1);		
	    }else if(x2>= (FADFADSTATE) 0.){
		p = FADsin(x1)*FADsin(x2);
	    }
	    break;
	case ECouplingNSD:
	    if(x2< (FADFADSTATE)1.){
		p = -(Pi*x2/4.)*FADcos(Pi*x1/2.)*FADsin(1.-x2);		
	    }else if(x2>= (FADFADSTATE)1.){
                p = (Pi/4.)*FADcos(Pi*x1/2.)*(x2-1.-FADcos(Pi*x2))*FADsin(x2-1.);
	    }
	    break;
        case ESinCos3D:
        case ESinCosBDS3D:
            p = FADcos(x1)*FADsin(x2)+FADcos(x2)*FADsin(x3);
            break;
        case EGatica3D:
            p = FADsin(Pi*x1)*FADsin(Pi*x2)*FADsin(Pi*x3);
            break;
        case EKovasznay:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            p = -(1./2.)*FADexp(2.*lambda*x1);
            break;
        case EKovasznayCDG:
	    Re = 1./fvisco; //Reynolds number
            lambda = Re/2.- sqrt(Re*Re/4.+4.*Pi*Pi); // Parameter for Navier-Stokes solution
            flux[0] = 1. - FADexp(lambda*x1)*FADcos(2.*Pi*x2);
            flux[1] = (lambda/(2.*Pi))*FADexp(lambda*x1)*FADsin(2.*Pi*x2);
            p = -(1./2.)*FADexp(2.*lambda*x1);
            p += (1./2.)*(flux[0]*flux[0]+flux[1]*flux[1]);
            break;
        case EPconst:
            p = 1;
            break;
        default:
            DebugStop();
    }
    
}

template<class TVar>
void TStokesAnalytic::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &gradu) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<TVar> temp = Fad<TVar>(3,i,x[i]);
        xfad[i] = temp;
    }
    //xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(3);
    uxy(xfad,result);
    gradu.Redim(3,3);
    for (int i=0; i<fDimension; i++) {
        for (int j=0; j<fDimension; j++)
        {
            gradu(i,j) = result[i].d(j);
        }
    }
}

template<class TVar>
void TStokesAnalytic::Duxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &Du) const
{
    TPZFMatrix<TVar> grad(3,3,0.), gradT(3,3,0.);
    graduxy(x,grad);
    grad.Transpose(&gradT);
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            Du(i,j) = 0.5*(grad(i,j)+gradT(i,j));
        }
    }
}


template<class TVar>
void TStokesAnalytic::SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZFMatrix<TVar> Du, pIdentity(sigma.Rows(),sigma.Cols());
    TVar p=0.;
    Duxy(x,Du);
    
    for(int i=0; i<Du.Rows(); i++)
    {
        for(int j=0; j<Du.Cols(); j++)
        {
            Du(i,j) *= 2.*fvisco;
        }
    }
    pressure(x, p);
    for (int i=0; i< pIdentity.Rows(); i++) {
        pIdentity(i,i) = p;
    }
    sigma = Du-pIdentity;
}

template<class TVar>
void TStokesAnalytic::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZFMatrix<TVar> Du(sigma.Rows(),sigma.Cols()), pIdentity(sigma.Rows(),sigma.Cols());
    TVar p=0.;
    Duxy(x,Du);
    for(int i=0; i<Du.Rows(); i++)
    {
        for(int j=0; j<Du.Cols(); j++)
        {
            Du(i,j) *= 2.*fvisco;
        }
    }
    pressure(x, p);
    for (int i=0; i< pIdentity.Rows(); i++) {
        pIdentity(i,i) = p;
    }
    
    for(int i=0; i<sigma.Rows(); i++)
    {
        for(int j=0; j<sigma.Cols(); j++)
        {
            sigma(i,j) = Du(i,j)-pIdentity(i,j);
        }
    }

}


void TStokesAnalytic::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
{
    typedef STATE TVar;
    TPZFMatrix<TVar> Du, pIdentity(sigma.Rows(),sigma.Cols());
    TVar p=0.;
    TPZManVector<STATE,3> xstate(3);
    for(int i=0; i<3; i++) xstate[i] = x[i];
    Duxy(xstate,Du);
    for(int i=0; i<Du.Rows(); i++)
    {
        for(int j=0; j<Du.Cols(); j++)
        {
            Du(i,j) *= 2.*fvisco;
        }
    }
    pressure(xstate, p);
    for (int i=0; i< pIdentity.Rows(); i++) {
        pIdentity(i,i) = p;
    }
    sigma = Du-pIdentity;
    
}




template
void TStokesAnalytic::SigmaLoc(const TPZVec<STATE> &x, TPZFMatrix<STATE> &sigma) const;


template<class TVar>
void TStokesAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const
{
    int sz = x.size();
    TPZManVector<Fad<TVar>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        xfad[i] = Fad<TVar>(sz,i,x[i]);
    }
    TPZFNMatrix<9, Fad<TVar> > sigma(3,3);
    Sigma(xfad,sigma);
    for (int i=0; i<3; i++) {
        divsigma[i] = sigma(i,0).dx(0)+sigma(i,1).dx(1)+sigma(i,2).dx(2);
    }
}

template
void TStokesAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<STATE> &divsigma) const;

void TStokesAnalytic::Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
{

    TPZManVector<STATE,3> locforce(3,0.),beta(3,0.),gradU_beta(3,0.),gradUt_beta(3,0.),xst(3);
    TPZFMatrix<STATE> grad(3,3,0.);
    for(int i=0; i<3; i++) xst[i] = x[i];
    DivSigma(x, locforce);

    switch(fProblemType)
    {
        case EStokes:
            force[0] = -locforce[0];
            force[1] = -locforce[1];
            force[2] = -locforce[2];
            break;

        case EBrinkman:
    	    uxy(xst,beta);
            force[0] = -locforce[0]+fcBrinkman*beta[0];
            force[1] = -locforce[1]+fcBrinkman*beta[1];
            force[2] = -locforce[2]+fcBrinkman*beta[2];
	    graduxy(xst,grad);
            force[3] = grad(0,0)+grad(1,1)+grad(2,2); //Pressure block term
            break;
            
        case ENavierStokes:
        case EOseen:

            graduxy(xst,grad);
            uxy(xst,beta);
            
            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    gradU_beta[e] += grad(e,f)*beta[f];
                }
            }
            
            force[0] = -locforce[0]+gradU_beta[0];
            force[1] = -locforce[1]+gradU_beta[1];
            force[2] = -locforce[2]+gradU_beta[2];
            break;

        case ENavierStokesCDG:
        case EOseenCDG:

            graduxy(xst,grad);
            uxy(xst,beta);

            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    gradU_beta[e] += grad(e,f)*beta[f];
                }
            }

            for (int e=0; e<3; e++) {
                for (int f=0; f<3; f++) {
                    gradUt_beta[e] += grad(f,e)*beta[f];
                }
            }

            force[0] = -locforce[0]+gradU_beta[0]-gradUt_beta[0];
            force[1] = -locforce[1]+gradU_beta[1]-gradUt_beta[1];
            force[2] = -locforce[2]+gradU_beta[2]-gradUt_beta[2];
            break;
            
        default:
            DebugStop();
    }



}

void TStokesAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &gradsol) const
{
//    TPZManVector<STATE> xst(3);
//    for(int i=0; i<3; i++) xst[i] = x[i];
//    uxy(xst,u);
//    graduxy(xst,gradu);
    
    TPZManVector<Fad<STATE>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<STATE> temp = Fad<STATE>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<STATE>,3> u_result(3);
    uxy(xfad,u_result);
    gradsol.Redim(3,3);
    sol.resize(4);
    for (int i = 0; i < 3; i++) {
        sol[i] = u_result[i].val();
    }
    for(int i=0; i<fDimension; i++) {
        for (int j=0; j<fDimension; j++)
        {
              gradsol(i,j) = u_result[j].d(i);
        }
    }
    Fad<STATE> p_result = 0.;
    pressure(xfad, p_result);
    sol[3] = p_result.val();

}

