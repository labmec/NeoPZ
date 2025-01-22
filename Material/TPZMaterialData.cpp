#include "TPZMaterialData.h"
#include <string>

TPZMaterialData::TPZMaterialData() :
    TPZRegisterClassId(&TPZMaterialData::ClassId), 
    normal(3,0.),x(3,0.), fUserData(0)
{

}

std::string TPZMaterialData::ShapeFunctionType() const
{
    switch(fShapeType){
        case EEmpty:
            return "NotInitialized";
        case EScalarShape:
            return "Scalar";
        case EVecandShape:
            return "Vector combined with Scalar";
        case EVecShape:
            return "Vector shape";
        default:
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Invalid value\n";
            DebugStop();
    }
    return "All Wrong!\n";
}

void TPZMaterialData::SetAllRequirements(bool set)
{
    this->fNeedsSol = set;
    this->fNeedsNeighborSol = set;
    this->fNeedsHSize = set;
    this->fNeedsNeighborCenter = set;
    this->fNeedsNormal = set;
}

void TPZMaterialData::Print(std::ostream &out) const
{
    out << "Shape function type " << ShapeFunctionType() << std::endl;
    out << "Active Approximation Space " << fActiveApproxSpace << std::endl;
    phi.Print("phi",out);
    fH1.fDPhi.Print("H1 fDPhi",out);
    dphix.Print("dphix",out);
    out << "Number dual functions " << numberdualfunctions << std::endl;
    divphi.Print("div phi",out);
    curlphi.Print("curl phi",out);
    axes.Print("axes",out);
    jacobian.Print("jacobian",out);
    jacinv.Print("jacinv",out);
    out << "normal " << normal << std::endl;
    out << "x " << x << std::endl;
    out << "xParametric " << xParametric << std::endl;
    out << "p " << p << std::endl;
    out << "HSize " << HSize << std::endl;
    out << "detjac " << detjac << std::endl;
    out << "XCenter " << XCenter << std::endl;
    out << "HDiv fMasterDirections" << fHDiv.fMasterDirections << std::endl;
    out << "HCurl fMasterDirections" << fHCurl.fMasterDirections << std::endl;
    out << "fDeformedDirections" << fDeformedDirections << std::endl;
    if(fNeedsDeformedDirectionsFad){
        fDeformedDirectionsFad.Print(out);
    }
    else
    {
        out << "No need for directions FAD\n";
    }
    out << "gelElId " << gelElId << std::endl;
    if (fVecShapeIndex.size()) {
        out << "VecShapeIndex: ";
        for (int64_t i = 0; i < fVecShapeIndex.size(); i++) {
            out << fVecShapeIndex[i].first << '/' << fVecShapeIndex[i].second << ' ';
        }
        out << '\n';
    }
    out << "NintPts " << NintPts << std::endl;
    out << "intLocPtIndex " << intLocPtIndex << std::endl;
    out << "intGlobPtIndex " << intGlobPtIndex << std::endl;
    out << "NeedsSol " << fNeedsSol << std::endl;
    out << "fNeedsNeighborSol " << fNeedsNeighborSol << std::endl;
    out << "fNeedsHSize " << fNeedsHSize << std::endl;
    out << "fNeedsNeighborCenter " << fNeedsNeighborCenter << std::endl;
    out << "fNeedsDeformedDirectionsFad " << fNeedsDeformedDirectionsFad << std::endl;
    out << "fNeedsNormal " << fNeedsNormal << std::endl;
}

/** Print the data in a format suitable for Mathematica */
void TPZMaterialData::PrintMathematica(std::ostream &out) const
{
    phi.Print("phi = ",out,EMathematicaInput);
    fH1.fDPhi.Print("H1 fDPhi = ",out,EMathematicaInput);
    dphix.Print("dphix = ",out,EMathematicaInput);
    axes.Print("axes = ",out,EMathematicaInput);
    jacobian.Print("jacobian = ",out,EMathematicaInput);
    jacinv.Print("jacinv = ",out,EMathematicaInput);
    out << "normal = {" << normal << "};" << std::endl;
    out << "x = {" << x << "};" << std::endl;
    out << "p = " << p << ";" << std::endl;
    out << "HSize = " << HSize << ";" << std::endl;
    out << "detjac = " << detjac << ";" << std::endl;
    out << "XCenter = {" << XCenter << "};" << std::endl;
    out << "HDivfMasterDirections" << fHDiv.fMasterDirections << std::endl;
    out << "HCurlfMasterDirections" << fHCurl.fMasterDirections << std::endl;
    out << "intLocPtIndex = " << intLocPtIndex << ";" <<std::endl;
    out << "intGlobPtIndex = " << intGlobPtIndex << ";" <<std::endl;
    out << "NintPts = " << NintPts << ";" <<std::endl;
    out << "gelElId = " << gelElId << ";" <<std::endl;
}


int TPZMaterialData::ClassId() const
{
    return Hash("TPZMaterialData");
}
void TPZMaterialData::Write(TPZStream &buf, int withclassid) const
{
    buf.Write(fNeedsSol);
    buf.Write(fNeedsNeighborSol);
    buf.Write(fNeedsHSize);
    buf.Write(fNeedsNeighborCenter);
    buf.Write(fNeedsNormal);
    buf.Write(fNeedsDeformedDirectionsFad);
    buf.Write((int)fActiveApproxSpace);
}
    
void TPZMaterialData::Read(TPZStream &buf, void *context)
{
    buf.Read(fNeedsSol);
    buf.Read(fNeedsNeighborSol);
    buf.Read(fNeedsHSize);
    buf.Read(fNeedsNeighborCenter);
    buf.Read(fNeedsNormal);
    buf.Read(fNeedsDeformedDirectionsFad);
    int64_t active;
    buf.Read(&active);
    fActiveApproxSpace = (MSpaceConfig) active;
}
