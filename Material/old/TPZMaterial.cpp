/**
 * @file
 * @brief Contains implementations of the TPZMaterial methods.
 */

#include "TPZMaterial.h"
#include "pzmaterialdata.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzbndcond.h"
#include "pzreal.h"
#include "pzadmchunk.h"
#include "tpzintpoints.h"
#include "TPZExactFunction.h"
#include "pzlog.h"
#include "TPZPersistenceManager.h"

#ifdef PZ_LOG
#ifdef PZDEBUG
#define DEBUG2
#endif
static TPZLogger logger("pz.material");
#endif

//TPZVec< void(*) (const TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

using namespace std;
REAL TPZMaterial::gBigNumber = 1.e12;


TPZMaterial::TPZMaterial() : fNumLoadCases(1), fPostProcIndex(0) {
	this->fId = -666;
	this->fForcingFunction = NULL;
  this->fExactSol = NULL;
  this->fTimeDependentForcingFunction = NULL;
  this->fTimedependentFunctionExact = NULL;
  this->fLinearContext = true;
}

TPZMaterial::TPZMaterial(int id) : fId(id), fNumLoadCases(1), fPostProcIndex(0) {
	this->SetId(id);
    this->fForcingFunction = NULL;
    this->fExactSol = NULL;
    this->fTimeDependentForcingFunction = NULL;
    this->fTimedependentFunctionExact = NULL;
    this->fLinearContext = true;

}

TPZMaterial::~TPZMaterial()
{
    this->fId = -999;
}


TPZMaterial::TPZMaterial(const TPZMaterial &material) {
	fId = material.fId;
    fForcingFunction = material.fForcingFunction;
    fExactSol = material.fExactSol;
    fTimeDependentForcingFunction = material.fTimeDependentForcingFunction;
    fTimedependentFunctionExact = material.fTimedependentFunctionExact;
    fLinearContext = material.fLinearContext;
    fNumLoadCases = material.fNumLoadCases;
    fPostProcIndex = material.fPostProcIndex;
}

TPZMaterial &TPZMaterial::operator=(const TPZMaterial &material)
{
    fId = material.fId;
    fForcingFunction = material.fForcingFunction;
    fExactSol = material.fExactSol;
    fTimeDependentForcingFunction = material.fTimeDependentForcingFunction;
    fTimedependentFunctionExact = material.fTimedependentFunctionExact;
    fLinearContext = material.fLinearContext;
    fNumLoadCases = material.fNumLoadCases;
    fPostProcIndex = material.fPostProcIndex;
    return *this;
}


void TPZMaterial::SetExactSol(std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)>f,int p)
{
    TPZAutoPointer<TPZFunction<STATE>> fp(
        new TPZExactFunction<STATE>(f,p));
    fExactSol = fp;
}
void TPZMaterial::GetExactSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col)
{
  static bool firstTime = true;
  if(firstTime)
    {
      firstTime = false;
      std::cout << "Using default implementation of " << std::endl;
      std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "If needed, override this method" << std::endl;
    }
  u_len = 1; du_row = 3; du_col = 1;
}

void TPZMaterial::SetLinearContext(bool IsLinear){
	fLinearContext = IsLinear;
}

void TPZMaterial::FillDataRequirements(TPZMaterialData &data)
{
	data.SetAllRequirements(true);
	data.fNeedsNeighborSol = false;
	data.fNeedsNeighborCenter = false;
	data.fNeedsNormal = false;

}

void TPZMaterial::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(true);
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = false;
	}
	
}

void TPZMaterial::FillDataRequirements(std::map<int, TPZMaterialData > &datavec)
{
    for(auto &it : datavec)
    {
        it.second.SetAllRequirements(true);
        it.second.fNeedsNeighborSol = false;
        it.second.fNeedsNeighborCenter = false;
        it.second.fNeedsNormal = false;
    }
    
}

void TPZMaterial::Print(std::ostream & out) {
    out << __PRETTY_FUNCTION__ << std::endl;
	out << "Material Id = " << fId << std::endl;
    out << "Linear context " << fLinearContext << std::endl;
    out << "Num loadcases " << fNumLoadCases << std::endl;
    out << "Big number " << gBigNumber << std::endl;
    
    if (!fForcingFunction) {
        out << "Has no forcing function\n";
    }
    else {
        out << "Forcing function\n";
        fForcingFunction->Print(out);
    }
    if (!fExactSol) {
        out << "Has no exact forcing function\n";
    }
    else {
        out << "Forcing function exact\n";
        fExactSol->Print(out);
    }
    if (!fTimeDependentForcingFunction) {
        out << "Has no time dependent forcing function\n";
    }
    else {
        out << "Time dependent forcing function\n";
        fTimeDependentForcingFunction->Print(out);
    }
    if (!fTimedependentFunctionExact) {
        out << "No time dependent forcing function exact\n";
    }
    else {
        out << "Time dependent forcing function exact\n";
        fTimedependentFunctionExact->Print(out);
    }
    
}

int TPZMaterial::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"state")) return 0;
	if(!strcmp(name.c_str(),"State")) return 0;
	if(!strcmp(name.c_str(),"Solution")) return 0;
    if(!strcmp(name.c_str(),"GradState")) return 1;
	if(!strcmp(name.c_str(),"POrder")) return 99;
	if(!strcmp(name.c_str(),"Error")) return 100;
	if(!strcmp(name.c_str(),"TrueError")) return 101;
	if(!strcmp(name.c_str(),"EffectivityIndex")) return 102;
	
	if(!strcmp(name.c_str(),"L2Error")) return 103;
	if(!strcmp(name.c_str(),"SemiH1Error")) return 104;
	if(!strcmp(name.c_str(),"H1Error")) return 105;
	
	if(!strcmp(name.c_str(),"L2ErrorPerArea")) return 106;
	if(!strcmp(name.c_str(),"SemiH1ErrorPerArea")) return 107;
	if(!strcmp(name.c_str(),"H1ErrorPerArea")) return 108;
	if(!strcmp(name.c_str(),"dudxErrorPerArea")) return 109;
	if(!strcmp(name.c_str(),"dudyErrorPerArea")) return 110;
	if(!strcmp(name.c_str(),"ContDisc")) return 111;
    if(!strcmp(name.c_str(),"MaterialId")) return 98;
	
	
//	std::cout << __PRETTY_FUNCTION__ << " Variable " << name << " not found\n";
	
#ifdef PZ_LOG2
	{
		std::stringstream sout;
		sout << "Variable " << name << " not found";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;
}

int TPZMaterial::NSolutionVariables(int index) {
#ifdef STATE_COMPLEX
	if(index == 0) return NStateVariables()*2;    
#else
	if(index == 0) return NStateVariables();
#endif
	if(index == 99 || index == 98) return 1;
	if(index == 100) return 1;
	if(index == 101) return 1;
	if(index == 102) return 1;
	if (index == 103) return 1;
	if (index == 104) return 1;
	if (index == 105) return 1;
	if (index == 106) return 1;
	if (index == 107) return 1;
	if (index == 108) return 1;
	if (index == 109) return 1;
	if (index == 110) return 1;
	if (index == 111) return 1;
	PZError << "TPZMaterial::NSolutionVariables called index = " << index << "\n";
    DebugStop();
	return 0;
}

void TPZMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
}

void TPZMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    int nvec = datavec.size();
    int numdata = 0;
    int dataindex = -1;
    for (int iv=0; iv<nvec; iv++) {
        if(datavec[iv].fShapeType != TPZMaterialData::EEmpty)
        {
            numdata++;
            dataindex = iv;
        }
    }
    if (numdata == 1) {
        Solution(datavec[dataindex], var, Solout);
        return;
    }
    DebugStop();
}

void TPZMaterial::Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout)
{
    DebugStop();
}

void TPZMaterial::Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl *left, TPZCompEl *right)
{
    DebugStop();
}

void TPZMaterial::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
						   TPZVec<STATE> &Solout){
    if(var == 98){
        Solout[0] = this->Id();
        return;
    }
#ifdef STATE_COMPLEX
    if(var == 0) 
    {
        Solout[0] = Sol[0].real();
        Solout[1] = Sol[0].imag();
    }
    else if(var == 99 || var == 100 || var == 101 || var == 102) {
        PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
        Solout[0] = Sol[0].real(); // = 0.;
    } 
	else 
    {
        Solout.Resize(0);
    }
#else
    if(var == 0) Solout = Sol;
    else if(var == 99 || var == 100 || var == 101 || var == 102) {
    PZError << "TPZMaterial var = "<< var << " the element should treat this case\n";
        Solout[0] = Sol[0]; // = 0.;
    } else if(var == 1)
    {
        Solout.resize(Sol.size()*3);
        int64_t nsol = Sol.size();
        Solout.Fill(0.);
        int64_t dim = axes.Rows();
        for (int64_t is=0; is<nsol; is++) {
            for (int64_t d=0; d<dim; d++) {
                for (int64_t jco=0; jco<3; jco++) {
                    Solout[jco+3*is] += axes(d,jco)*DSol(d,is);
                }
            }
        }
    } else
    {
        DebugStop();
        Solout.Resize(0);
    }
#endif
}

TPZBndCond *TPZMaterial::CreateBC(TPZMaterial * reference, int id, int typ, const TPZFMatrix<STATE> &val1, const TPZFMatrix<STATE> &val2) {
	return new TPZBndCond(reference,id,typ,val1,val2);
}

void TPZMaterial::SetData(std::istream &data) {
	PZError << "TPZMaterial::SetData is called.\n";
	data >> fId;
}

TPZMaterial * TPZMaterial::NewMaterial() {
	PZError << "TPZMaterial::NewMaterial is called.\n";
	return 0;
}

void TPZMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef){
    TPZFMatrix<REAL>  &phi = data.phi;
    const auto phr = phi.Rows();
    for( auto in = 0; in < phr; in++ ) {
        ef(in) = phi(in,0);
        for( auto jn = 0; jn < phr; jn++ ) {
            ek(in,jn) = phi(in,0) * phi(jn,0);
        }
    }
}

void TPZMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	int nref=datavec.size();
    int ndif = 0;
    int onemat = 0;
    for (int ir = 0; ir < nref; ir++) {
        int nphis=datavec[ir].phi.Rows();
        if (datavec[ir].phi.Rows()) {
            onemat = ir;
            ndif++;
        }
    }
	if (ndif == 1) {
		this->Contribute(datavec[onemat], weight, ek,ef);
	}
    else
    {
        DebugStop();
    }
}

void TPZMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
							   TPZFMatrix<STATE> &ef, TPZBndCond &bc){
	int nref=datavec.size();
	if (nref== 1) {
		this->ContributeBC(datavec[0], weight, ek,ef,bc);
	}
    else
    {
        DebugStop();
    }
}

void TPZMaterial::Clone(std::map<int, TPZMaterial * >&matvec) {
	int matid = Id();
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZMaterial * newmat = NewMaterial();
	matvec[matid] = newmat;
}

/** Get the order of the integration rule necessary to integrate an
 * element with polinomial order p */
int TPZMaterial::IntegrationRuleOrder(int elPMaxOrder) const
{
    int order = 0;
    if(fForcingFunction){
        order = fForcingFunction->PolynomialOrder();
    }
    
    int pmax = elPMaxOrder;
    int integrationorder = 2*pmax;
    if (pmax < order) {
        integrationorder = order+pmax;
    }
    return  integrationorder;
}

int TPZMaterial::IntegrationRuleOrder(TPZVec<int> &elPMaxOrder) const {
    int order = 0;
    if (fForcingFunction) {
        order = fForcingFunction->PolynomialOrder();
    }

    int pmax = 0;
    for (int ip = 0; ip < elPMaxOrder.size(); ip++) {
        if (elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }
    int integrationorder = 2 * pmax;
    if (pmax < order) {
        integrationorder = order + pmax;
    }

    return integrationorder;
}

int TPZMaterial::ClassId() const{
    return Hash("TPZMaterial");
}

/* Saves the element data to a stream */
void TPZMaterial::Write(TPZStream &buf, int withclassid) const {
    buf.Write(&fId, 1);
    buf.Write(&gBigNumber, 1);
    TPZPersistenceManager::WritePointer(fForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fExactSol.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimeDependentForcingFunction.operator ->(), &buf);
    TPZPersistenceManager::WritePointer(fTimedependentFunctionExact.operator ->(), &buf);
    buf.Write(fLinearContext);
    buf.Write(&fNumLoadCases);
    buf.Write(&fPostProcIndex);
}

/* Reads the element data from a stream */
void TPZMaterial::Read(TPZStream &buf, void *context) {
    buf.Read(&fId, 1);
    buf.Read(&gBigNumber, 1);
    fForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fExactSol = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimeDependentForcingFunction = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    fTimedependentFunctionExact = TPZAutoPointerDynamicCast<TPZFunction<STATE>>(TPZPersistenceManager::GetAutoPointer(&buf));
    buf.Read(fLinearContext);
    buf.Read(&fNumLoadCases);
    buf.Read(&fPostProcIndex);
}

void TPZMaterial::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                   REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
    DebugStop();
}

void TPZMaterial::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                   REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);
}

void TPZMaterial::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright,
                                                   REAL weight, TPZFMatrix<STATE> &ef){
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeInterface(data, dataleft, dataright, weight, fakeek, ef);
}

void TPZMaterial::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}
//void TPZMaterial::ContributeInterface(TPZVec<TPZMaterialData> &datavec, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec,
//                                 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
//    std::cout << __PRETTY_FUNCTION__ << " please implement me\n";
//    DebugStop();
//}

void TPZMaterial::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    DebugStop();
}
/**
 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point to multiphysics simulation
 * @param data [in]
 * @param dataleft [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @param bc [in] is the boundary condition object
 * @since February 21, 2013
 */
void TPZMaterial::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<STATE> ek(ef.Rows(),ef.Rows());
    this->ContributeBCInterface(data, dataleft, weight, ek, ef, bc);
}

void TPZMaterial::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    DebugStop();
}


void TPZMaterial::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    TPZFMatrix<STATE> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeBCInterface(data, dataleft, weight, fakeek, ef, bc);
}

int TPZMaterial::IsInterfaceConservative(){
    return 0;
}

void TPZMaterial::InterfaceJump(TPZVec<REAL> &x,
                                             TPZSolVec &leftu,
                                             TPZSolVec &rightu,
                                             TPZSolVec &jump){
    int numbersol = leftu.size();
    for (int is=0; is<numbersol; is++) {
        const int n = leftu[is].NElements();
        jump[is].Resize(n);
        for(int i = 0; i < n; i++){
            jump[is][i] = leftu[is][i] - rightu[is][i];
        }
    }
}


void TPZMaterial::BCInterfaceJump(TPZVec<REAL> &x,
                                               TPZSolVec &leftu,
                                               TPZBndCond &bc,
                                               TPZSolVec & jump){
    PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
    DebugStop();
}

/// return the integration order as a function of interpolation orders of the left and right elements
int TPZMaterial::GetIntegrationOrder(TPZVec<int> &porder_left, TPZVec<int> &porder_right) const
{
    int maxl = 0, maxr = 0;
    for (auto porder: porder_left) {
        maxl = std::max(maxl,porder);
    }
    for (auto porder: porder_right) {
        maxr = std::max(maxr,porder);
    }
    return maxl+maxr;
}
