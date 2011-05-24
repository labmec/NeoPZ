//
//  TPBrSteamSimulation.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/15/11.
//  Copyright 2011 UNICAMP. All rights reserved.
//

#include "TPBrSteamSimulation.h"
#include "tpbrsteammesh.h"
#include "ThermalMethodsTables.h"
#include "pzlog.h"

WaterDataInStateOfSaturation waterdata;
OilData oildata;

int main()
{
#ifdef LOG4CXX
	InitializePZLOG();
#endif
#ifdef _AUTODIFF
//    TPBrSteamMesh::TPBrSteamMesh(int numcells, REAL temperature, REAL pressure, REAL WellRadius, REAL ReservoirRadius, REAL oilsaturation)
    TPBrSteamMesh mesh(1,100.,1.e4,0.2,200.,0.);
    mesh.Print();
    int neq = mesh.NumEquations();

    TPZFMatrix tangent(neq,neq),residual(neq,1);
    mesh.ComputeTangent(tangent, residual);
#endif
    return 0;
}