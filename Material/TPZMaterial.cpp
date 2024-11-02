#include "TPZMaterial.h"
#include "TPZMaterialData.h"
#include "TPZStream.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatCombinedSpaces.h"

#include <string>

int TPZMaterial::VariableIndex(const std::string &name) const{
	const TPZMatSingleSpace * matsingle = dynamic_cast<const TPZMatSingleSpace *>(this);
	if(matsingle) return matsingle->VariableIndex(name);
	const TPZMatCombinedSpaces *matcomb = dynamic_cast<const TPZMatCombinedSpaces *>(this);
	if(matcomb) return matcomb->VariableIndex(name);
	std::cout << __PRETTY_FUNCTION__ << " Variable " << name << " not found\n";
	
#ifdef PZ_LOG2
	{
		std::stringstream sout;
		sout << "Variable " << name << " not found";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;
}



int TPZMaterial::NSolutionVariables(int index) const{
	const TPZMatSingleSpace * matsingle = dynamic_cast<const TPZMatSingleSpace *>(this);
	if(matsingle) return matsingle->NSolutionVariables(index);
	const TPZMatCombinedSpaces *matcomb = dynamic_cast<const TPZMatCombinedSpaces *>(this);
	if(matcomb) return matcomb->NSolutionVariables(index);
	PZError<<__PRETTY_FUNCTION__;
	PZError<<" Implement me in your material for post processing solutions.\n";
	DebugStop();
	return 0;
}


TPZMaterial* TPZMaterial::NewMaterial() const{
	PZError<<__PRETTY_FUNCTION__;
	PZError<<" should not be called. Implement me! Aborting...\n";
	DebugStop();
	return nullptr;
}
void TPZMaterial::Clone(std::map<int, TPZMaterial * >&matvec) {
	int matid = Id();
	std::map<int, TPZMaterial * >::iterator matit;
	matit = matvec.find(matid);
	if(matit != matvec.end()) return;
	TPZMaterial * newmat = NewMaterial();
	matvec[matid] = newmat;
}

int TPZMaterial::ClassId() const{
    return Hash("TPZMaterial");
}

void TPZMaterial::Print(std::ostream & out) const {
	auto printbool = [](const bool val){
		if (val) return "yes";
		else return "no";
	};
    out << __PRETTY_FUNCTION__ << '\n';
	out << "Material Id: " << fId << '\n';
	out << "Dimension: " << Dimension() << '\n';
	out << "Num state variables: " << NStateVariables() << '\n';
	out << "Has forcing function: " << printbool(HasForcingFunction()) << '\n';
}


void TPZMaterial::Read(TPZStream& buf, void* context) {
    buf.Read(&fId);
}

void TPZMaterial::Write(TPZStream& buf, int withclassid) const {
    buf.Write(&fId);
}