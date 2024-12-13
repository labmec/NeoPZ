//
//  pzmultiphysicselement.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/25/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "pzmultiphysicselement.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzlog.h"
#include "pzinterpolationspace.h"
#include "pzcompelwithmem.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZMultiPhysicsElement");
#endif

TPZMultiphysicsElement::TPZMultiphysicsElement(TPZCompMesh &mesh, const TPZMultiphysicsElement &copy) : TPZCompEl(mesh,copy)
{
  
}

void TPZMultiphysicsElement::CreateInterfaces()
{
    //nao verifica-se caso o elemento de contorno
    //eh maior em tamanho que o interface associado
    //caso AdjustBoundaryElement nao for aplicado
    //a malha eh criada consistentemente
    TPZGeoEl *ref = Reference();
    int nsides = ref->NSides();
    int InterfaceDimension = Mesh()->Dimension() - 1;
    int side;
    nsides--;//last face
    for(side=nsides;side>=0;side--) {
        if(ref->SideDimension(side) != InterfaceDimension) continue;
        TPZCompElSide thisside(this,side);
        if(this->ExistsInterface(side)) {
            //      cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
            continue;
        }
        TPZStack<TPZCompElSide> highlist;
        thisside.HigherLevelElementList(highlist,0,1);
        //a interface se cria uma vez so quando existem ambos
        //elementos esquerdo e direito (compu tacionais)
        if(!highlist.NElements()) {
            this->CreateInterface(side);//s�tem iguais ou grande => pode criar a interface
        } else {
            int64_t ns = highlist.NElements();
            int64_t is;
            for(is=0; is<ns; is++) {//existem pequenos ligados ao lado atual
                const int higheldim = highlist[is].Reference().Dimension();
                if(higheldim != InterfaceDimension) continue;
                // 	TPZCompElDisc *del = dynamic_cast<TPZCompElDisc *> (highlist[is].Element());
                // 	if(!del) continue;
                
                TPZCompEl *del = highlist[is].Element();
                if(!del) continue;
                
                TPZCompElSide delside( del, highlist[is].Side() );
                TPZMultiphysicsElement * delSp = dynamic_cast<TPZMultiphysicsElement *>(del);
                if (!delSp){
                    PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
                    return;
                }
                if ( delSp->ExistsInterface(highlist[is].Side()) ) {
                    //          cout << "TPZCompElDisc::CreateInterface inconsistent: interface already exists\n";
                }
                else {
                    delSp->CreateInterface(highlist[is].Side());
                }
            }
        }
    }
}

int TPZMultiphysicsElement::ComputeIntegrationOrder() const {
    DebugStop();
    return 0;
}

TPZMultiphysicsInterfaceElement * TPZMultiphysicsElement::CreateInterface(int side)
{
	//  LOGPZ_INFO(logger, "Entering CreateInterface");
	TPZMultiphysicsInterfaceElement * newcreatedinterface = NULL;
	
	TPZGeoEl *ref = Reference();
	if(!ref) {
		LOGPZ_ERROR(logger, "Exiting CreateInterface Null reference reached - NULL interface returned");
		return newcreatedinterface;
	}
	
	TPZCompElSide thisside(this,side);
	TPZStack<TPZCompElSide> list;
	list.Resize(0);
	thisside.EqualLevelElementList(list,0,0);//retorna distinto ao atual ou nulo
	int64_t size = list.NElements();
	//espera-se ter os elementos computacionais esquerdo e direito
	//ja criados antes de criar o elemento interface
    // try to create an interface element between myself and an equal sized neighbour
    for (int64_t is=0; is<size; is++)
    {
		//Interface has the same material of the neighbour with lesser dimension.
		//It makes the interface have the same material of boundary conditions (TPZCompElDisc with interface dimension)
		int matid;
		int thisdim = this->Dimension();
		int neighbourdim = list[is].Element()->Dimension();
        
        if(thisdim != neighbourdim)
        {
            if (thisdim < neighbourdim)
            {
                // return the material id of boundary condition IF the associated material is derived from bndcond
                TPZMaterial *mat = this->Material();
                TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                if(bnd)
                {
                    matid = this->Material()->Id();
                }
                else
                {
                    matid = this->Mesh()->Reference()->InterfaceMaterial(this->Reference()->MaterialId(),list[0].Element()->Reference()->MaterialId());
                    continue;
                }
            }
			else
            {
                TPZMaterial *mat = list[is].Element()->Material();
                TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                if (bnd) {
                    matid = bnd->Id();
                }
                else {
                    matid = this->Mesh()->Reference()->InterfaceMaterial(this->Reference()->MaterialId(), list[0].Element()->Reference()->MaterialId());
                    continue;
                }
            }
        }else
        {
            matid = this->Mesh()->Reference()->InterfaceMaterial(this->Reference()->MaterialId(), list[0].Element()->Reference()->MaterialId());
            if(matid == GMESHNOMATERIAL)
            {
                continue;
            }
        }
		
		
		TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid); //isto acertou as vizinhanas da interface geometrica com o atual
		if(!gel){
			DebugStop();
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << "CreateBCGeoEl devolveu zero!@@@@";
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
        bool withmem = fMesh->ApproxSpace().NeedsMemory();
		if(Dimension() > list[is].Reference().Dimension()) {
			//o de volume eh o direito caso um deles seja BC
			//a normal aponta para fora do contorno
			TPZCompElSide thiscompelside(this, thisside.Side());
			TPZCompElSide neighcompelside(list[is]);
            if (!withmem) {
                newcreatedinterface = new TPZMultiphysicsInterfaceElement(*fMesh,gel,thiscompelside,neighcompelside);
            }
            else
            {
                newcreatedinterface = new TPZCompElWithMem<TPZMultiphysicsInterfaceElement>(*fMesh,gel,thiscompelside,neighcompelside);
            }
		} else {
			//caso contrario ou caso ambos sejam de volume
			TPZCompElSide thiscompelside(this, thisside.Side());
			TPZCompElSide neighcompelside(list[is]);
            if (!withmem) {
                newcreatedinterface = new TPZMultiphysicsInterfaceElement(*fMesh,gel,neighcompelside,thiscompelside);
            }
            else
            {
                newcreatedinterface = new TPZCompElWithMem<TPZMultiphysicsInterfaceElement>(*fMesh,gel,neighcompelside,thiscompelside);
            }
		}
		
		
		
		return newcreatedinterface;
	}
	
	//If there is no equal or lower level element, we try the lower elements.
	//Higher elements will not be considered by this method. In that case the interface must be created by the neighbour.
	TPZCompElSide lower = thisside.LowerLevelElementList(0);
	if(lower.Exists()){
		//Interface has the same material of the neighbour with lesser dimension.
		//It makes the interface has the same material of boundary conditions (TPZCompElDisc with interface dimension)
		int matid = GMESHNOMATERIAL;
		int thisdim = this->Dimension();
		int neighbourdim = lower.Element()->Dimension();
        matid = this->Mesh()->Reference()->InterfaceMaterial(this->Material()->Id(), lower.Element()->Material()->Id() );
		
		if (matid == GMESHNOMATERIAL && thisdim == neighbourdim){
			//      matid = this->Material()->Id();
            //break;
        }
        else if(matid == GMESHNOMATERIAL && thisdim != neighbourdim)
        {
			if (thisdim < neighbourdim) 
            {
                // return the material id of boundary condition IF the associated material is derived from bndcond
                TPZMaterial *mat = this->Material();
                TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                if(bnd)
                {
                    matid = this->Material()->Id();
                }
                else {
                    //continue;
                }
            }
			else 
            {
                TPZMaterial *mat = lower.Element()->Material();
                TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
                if (bnd) {
                    matid = bnd->Id();
                }
                else {
                    //continue;
                }
            }
		}
		
        // return zero
        if(matid == GMESHNOMATERIAL)
        {
            return newcreatedinterface;
        }
		
		TPZCompEl *lowcel = lower.Element();
		//int lowside = lower.Side();
		//existem esquerdo e direito: this e lower
		TPZGeoEl *gel = ref->CreateBCGeoEl(side,matid);
		
        bool withmem = fMesh->ApproxSpace().NeedsMemory();
        
		if(Dimension() > lowcel->Dimension()){
			//para que o elemento esquerdo seja de volume
			TPZCompElSide thiscompelside(this, thisside.Side());
			TPZCompElSide lowcelcompelside(lower);
            if (!withmem) {
                newcreatedinterface = new TPZMultiphysicsInterfaceElement(*fMesh,gel,thiscompelside,lowcelcompelside);
            }
            else
            {
                newcreatedinterface = new TPZCompElWithMem<TPZMultiphysicsInterfaceElement>(*fMesh,gel,thiscompelside,lowcelcompelside);
            }
		} else {
			TPZCompElSide thiscompelside(this, thisside.Side());
			TPZCompElSide lowcelcompelside(lower);
#ifdef PZ_LOG_KEEP
            if (logger.isDebugEnabled())
			{
				std::stringstream sout;
				sout << __PRETTY_FUNCTION__ << " left element";
				sout << lowcelcompelside << thiscompelside;
				sout << "Left Element ";
				lowcelcompelside.Element()->Print(sout);
				sout << "Right Element ";
				thiscompelside.Element()->Print(sout);
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
            if (!withmem)
            {
                newcreatedinterface = new TPZMultiphysicsInterfaceElement(*fMesh,gel,lowcelcompelside,thiscompelside);
            }
            else
            {
                newcreatedinterface = new TPZCompElWithMem<TPZMultiphysicsInterfaceElement>(*fMesh,gel,lowcelcompelside,thiscompelside);
            }
		}
		
		
		return newcreatedinterface;
	}
	return newcreatedinterface;
}

bool TPZMultiphysicsElement::ExistsInterface(int side){
	
    TPZGeoElSide geosd(Reference(),side);
	TPZGeoElSide  neighside = geosd.Neighbour();
	while(neighside.Element() && neighside.Element() != geosd.Element()){
		TPZCompElSide neighcompside = neighside.Reference();
		neighside = neighside.Neighbour();
        TPZCompEl *cel = neighcompside.Element();
        if (!cel) {
            continue;
        }
        TPZMultiphysicsInterfaceElement *mfintel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (mfintel) {
            TPZCompElSide left = mfintel->Left();
            TPZCompElSide right = mfintel->Right();
            if(left.Element() == this || right.Element() == this)
            {
                return true;
            }
        }
	}
	return false;
}

void TPZMultiphysicsElement::RemoveInterfaces(){
	
	int nsides = Reference()->NSides();
	if (!this->Material()){
		std::stringstream mess;
		mess << __PRETTY_FUNCTION__ << " - this->Material() == NULL, I can't RemoveInterfaces()";
		PZError << mess.str() << std::endl;
		LOGPZ_ERROR(logger, mess.str());
		return;
	}
	int InterfaceDimension = this->Material()->Dimension() - 1;
	int is;
	TPZStack<TPZCompElSide> list,equal;
	for(is=0;is<nsides;is++){
		TPZCompElSide thisside(this,is);
		if(thisside.Reference().Dimension() != InterfaceDimension) continue;
		// procurar na lista de elementos iguais
		list.Resize(0);// o lado atual �uma face
		//thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
		RemoveInterface(is);// chame remove interface do elemento atual (para o side atual)
		thisside.HigherLevelElementList(list,0,0);// procurar na lista de elementos menores (todos)
		int64_t size = list.NElements(),i;            // 'isto pode incluir elementos interfaces'
		//tirando os elementos de interface da lista
		for(i=0;i<size;i++){
			if(list[i].Element()->Type() == EInterface) {
#ifdef PZ_LOG
				LOGPZ_DEBUG(logger, "Removing interface element from the list of higher level elements");
#endif
				//This need to be done because otherwise list could be invalidated when an interface is removed.
				list[i] = TPZCompElSide();//tirando interface
			}
		}
		for(i=0;i<size;i++){// percorre os elementos menores
			if(!list[i].Element()) continue;
			TPZGeoElSide geolist = list[i].Reference();//TESTE
			if(geolist.Dimension() != InterfaceDimension) continue;
			equal.Resize(0);// para cada elemento menor e' preciso verificar a dimensao,
			list[i].EqualLevelElementList(equal,0,0);//montar a lista de elementos iguais (todos)
			equal.Push(list[i]);//n� �incorporado no m�odo anterior
			int neq = equal.NElements(),k=-1;
			while(++k < neq) if(equal[k].Element()->Type() != EInterface) break;//procurando elemento descont�uo cujo
			if(!neq || k == neq){                               //lado faz parte da parti� do lado side do this
				LOGPZ_FATAL(logger, " Inconsistency of data");
				DebugStop();//elemento descont�uo n� achado: ERRO
			}// chame removeinterface do elemento menor
			
			TPZMultiphysicsElement * equalkSp = dynamic_cast<TPZMultiphysicsElement*>(equal[k].Element());
			if (!equalkSp){
				PZError << "\nERROR AT " << __PRETTY_FUNCTION__ <<  " - CASE NOT AVAILABLE\n";
				return;
			}
			equalkSp->RemoveInterface(equal[k].Side());
		}
	}
	
}

void TPZMultiphysicsElement::RemoveInterface(int side) {
	
	TPZStack<TPZCompElSide> list;
	list.Resize(0);
	TPZCompElSide thisside(this,side);
	thisside.EqualLevelElementList(list,0,0);// monta a lista de elementos iguais
	int64_t size = list.NElements(), i=-1;
	while(++i < size) if(list[i].Element()->Type() == EInterface) break;// procura aquele que e derivado de TPZInterfaceEl
	if(!size || i == size){
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ << " no interface element found\n";
			Print(sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return;// nada a ser feito
	}
	// aqui existe a interface
	TPZCompEl *cel = list[i].Element();
#ifdef PZ_LOG
	TPZGeoEl *gel = cel->Reference();
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " element index " << Index() << " side " << std::endl;
		sout << "geometric element reference count " << gel->NumInterfaces();
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	delete cel;
}

template<class TVar>
void TPZMultiphysicsElement::ComputeRequiredDataT(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, std::map<int, TPZMaterialDataT<TVar>> &datavec   )
{
    for (auto &it : datavec) {
        int elindex = it.first;
        TPZCompEl *cel = Element(elindex);
#ifdef PZDEBUG
        if(!cel) DebugStop();
#endif
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        TPZGeoEl *gel = intel->Reference();
        TPZManVector<REAL> locpt(gel->Dimension());
        trvec[elindex].Apply(point, locpt);
        datavec[elindex].intGlobPtIndex = -1;
        intel->ComputeRequiredData(it.second, locpt);
    }
}


template<class TVar>
void TPZMultiphysicsElement::ComputeRequiredDataT(TPZVec<REAL> &intpointtemp, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialDataT<TVar>> &datavec)
{
    int64_t ElemVecSize = NMeshes();
    for (int64_t iref = 0; iref < ElemVecSize; iref++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(Element(iref));
        if (!msp) {
            continue;
        }
        
        TPZManVector<REAL,3> intpoint(msp->Reference()->Dimension(),0.);
        
        trvec[iref].Apply(intpointtemp, intpoint);
        
        msp->ComputeRequiredData(datavec[iref], intpoint);
    }
}//ComputeRequiredData


void TPZMultiphysicsElement::TransferMultiphysicsElementSolution()
{
    if(fMesh->GetSolType() == EReal){
        return TransferMultiphysicsElementSolutionT<STATE>();
    }
    if(fMesh->GetSolType() == EComplex){
        return TransferMultiphysicsElementSolutionT<CSTATE>();
    }
    PZError<<__PRETTY_FUNCTION__<<'\n';
    PZError<<"Invalid type! Aborting...\n";
    DebugStop();
}
template<class TVar>
void TPZMultiphysicsElement::TransferMultiphysicsElementSolutionT()
{
    int nmeshes = this->NMeshes();
    int icon = 0;
    int nload = this->Mesh()->Solution().Cols();
    for (int imesh = 0; imesh < nmeshes; imesh++) {
        TPZCompEl *cel = this->ReferredElement(imesh);
        if (!cel) {
            continue;
        }
        if(!this->IsActiveApproxSpaces(imesh)){
            continue;
        }
        int ncon = cel->NConnects();
        for (int iconloc = 0; iconloc < ncon; iconloc++, icon++) {
            TPZConnect &con = this->Connect(icon);
            TPZConnect &conloc = cel->Connect(iconloc);
            int64_t seq = con.SequenceNumber();
            int64_t seqloc = conloc.SequenceNumber();
            int blsz = this->Mesh()->Block().Size(seq);

#ifdef PZDEBUG
            int blszloc = cel->Mesh()->Block().Size(seqloc);
            if (blsz != blszloc) {
                DebugStop();
            }
#endif
            int pos = this->Mesh()->Block().Position(seq);
            int posloc = cel->Mesh()->Block().Position(seqloc);
            TPZFMatrix<TVar> &celSol = cel->Mesh()->Solution();
            TPZFMatrix<TVar> &meshSol = this->Mesh()->Solution();

            int nrows = meshSol.Rows();
            int ncols = meshSol.Cols();
//            celSol.Resize(nrows,ncols);
            for (int ibl = 0; ibl < blsz; ibl++) {
                for (int iload = 0; iload < nload; iload++) {
                    celSol(posloc+ibl,iload) = meshSol(pos+ibl,iload);
                }
            }
        }
    }
}

void TPZMultiphysicsElement::EvaluateError(TPZVec<REAL> &errors, bool store_errors){
    
    DebugStop(); // Should never enter here
    

}//method

void TPZMultiphysicsElement::GetConnectMeshPairs(TPZVec<std::pair<int64_t,int>> &connectpairs)
{
    int nconnect = NConnects();
    connectpairs.Resize(nconnect);
    int nmeshes = NMeshes();
    int count = 0;
    for (int imesh = 0; imesh < nmeshes; imesh++) {
        if(!fActiveApproxSpace[imesh]) continue;
        TPZCompEl *cel = Element(imesh);
        if(!cel) continue;
        int nc = cel->NConnects();
        for (int ic=0; ic < nc; ic++) {
            connectpairs[count+ic] = std::pair<int64_t, int>(ConnectIndex(count+ic),imesh);
        }
        count += nc;
    }
#ifdef PZDEBUG
    if(count != nconnect) DebugStop();
#endif
}

