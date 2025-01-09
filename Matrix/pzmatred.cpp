/**
 * @file
 * @brief Contains the implementation of the TPZMatRed methods.
 */

#include <stdlib.h>
#include <stdio.h>
#include <fstream>



#include "pzmatred.h"
#include "pzfmatrix.h"
#include "pzstepsolver.h"
#include "tpzverysparsematrix.h"

#include "TPZPersistenceManager.h"

#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzmatred");
#else
static int logger;
#endif

/*************************** Public ***************************/

/******************/
/*** Construtor ***/

template<class TVar, class TSideMatrix>
TPZMatRed<TVar,  TSideMatrix>::TPZMatRed () : 
TPZRegisterClassId(&TPZMatRed::ClassId),
TPZMatrix<TVar>( 0, 0 ), fK11(0,0),fK01(0,0),fK10(0,0),fF0(0,1),fF1(0,1), fMaxRigidBodyModes(0), fNumberRigidBodyModes(0)
{
	fDim0=0;
	fDim1=0;
	fK01IsComputed = 0;
	fIsReduced = 0;
    fF0IsComputed = false;
}

template<class TVar, class TSideMatrix>
TPZMatRed<TVar, TSideMatrix>::TPZMatRed( int64_t dim, int64_t dim00 ):
TPZRegisterClassId(&TPZMatRed::ClassId),
TPZMatrix<TVar>( dim,dim ), fK11(dim-dim00,dim-dim00,0.), fK01(dim00,dim-dim00,0.),
fK10(dim-dim00,dim00,0.), fF0(dim00,1,0.),fF1(dim-dim00,1,0.), fMaxRigidBodyModes(0), fNumberRigidBodyModes(0)
{
	if(dim<dim00) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"dim k00> dim");
	fDim0=dim00;
	fDim1=dim-dim00;
	fK01IsComputed = 0;
    fF0IsComputed = false;
	fIsReduced = 0;
}

template<class TVar, class TSideMatrix>
TPZMatRed<TVar, TSideMatrix >::~TPZMatRed(){
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::SimetrizeMatRed() {
	// considering fK00 is simetric, only half of the object is assembled.
	// this method simetrizes the matrix object
	
	if(!fK00 || this->GetSymmetry() == SymProp::NonSym) return;
	fK01.Transpose(&fK10);
	
	fK11.Symetrize();
	
	/*
	int64_t row,col;
	for(row=0; row<fDim1; row++) {
		for(col=row+1; col<fDim1; col++) {
			(fK11)(col,row) = (fK11)(row,col);
		}
	}
	*/
}

template<class TVar, class TSideMatrix>
int
TPZMatRed<TVar, TSideMatrix>::PutVal(const int64_t r,const int64_t c,const TVar& value ){
	int64_t row(r),col(c);
  auto val = value;
	if (this->GetSymmetry() != SymProp::NonSym && row > col ) {
    Swap( &row, &col );
    if constexpr (is_complex<TVar>::value){
      if (this->GetSymmetry() == SymProp::Herm){
        val = std::conj(val);
      }
    }
  }
	if (row<fDim0 &&  col<fDim0)  fK00->PutVal(row,col,val);
	if (row<fDim0 &&  col>=fDim0)  fK01.PutVal(row,col-fDim0,(TVar)val);
	if (row>=fDim0 &&  col<fDim0)  fK10.PutVal(row-fDim0,col,(TVar)val);
	if (row>=fDim0 &&  col>=fDim0)  fK11.PutVal(row-fDim0,col-fDim0,(TVar)val);

	return( 1 );
}

template<class TVar, class TSideMatrix>
const TVar
TPZMatRed<TVar,TSideMatrix>::GetVal(const int64_t r,const int64_t c ) const {
	int64_t row(r),col(c);
	const auto sp = this->GetSymmetry();
	if (sp != SymProp::NonSym && row > col ) Swap( &row, &col );

  TVar val = -1;
	if (row<fDim0 &&  col<fDim0)  val = fK00->GetVal(row,col);
	else if (row<fDim0 &&  col>=fDim0)  val = fK01.GetVal(row,col-fDim0);
	else if (row>=fDim0 &&  col<fDim0)  val = fK10.GetVal(row-fDim0,col);
	else val = fK11.GetVal(row-fDim0,col-fDim0) ;

  if constexpr(is_complex<TVar>::value){
    if(sp == SymProp::Herm){val = std::conj(val);}
  }
  return val;
}

template<class TVar, class TSideMatrix>
TVar& TPZMatRed<TVar,TSideMatrix>::s(const int64_t r,const int64_t c ) {
	int64_t row(r),col(c);
	const auto sp = this->GetSymmetry();
  if constexpr (is_complex<TVar>::value) {
    if(sp == SymProp::Herm && row > col){
      PZError<<__PRETTY_FUNCTION__
             <<"\nTrying to get reference to lower triang section"
             <<"\nof symmetric-storage hermitian matrix.\nAborting..."
             <<std::endl;
      DebugStop();
    }
  }
	if (r < fDim0 && sp != SymProp::NonSym && row > col ) Swap( &row, &col );
	if (row<fDim0 &&  col<fDim0)  return ( fK00->s(row,col) );
	if (row<fDim0 &&  col>=fDim0)  return ( (TVar &)fK01.s(row,col-fDim0) );
	if (row>=fDim0 &&  col<fDim0)  return ( (TVar &)(fK10.s(row-fDim0,col)) );
	return ((TVar &)(fK11.s(row-fDim0,col-fDim0)) );
	
}

template<class TVar, class TSideMatrix>
int64_t TPZMatRed<TVar,TSideMatrix>::Size() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: should not be called for TPZMatRed.\nAborting...\n";
  DebugStop();
  return -1;
}
template<class TVar, class TSideMatrix>
TVar *&TPZMatRed<TVar,TSideMatrix>::Elem()
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: should not be called for TPZMatRed.\nAborting...\n";
  DebugStop();
  static TVar* t{nullptr};
  return t;
}
template<class TVar, class TSideMatrix>
const TVar *TPZMatRed<TVar,TSideMatrix>::Elem() const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: should not be called for TPZMatRed.\nAborting...\n";
  DebugStop();
  return nullptr;
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar,TSideMatrix>::SetSolver(TPZAutoPointer<TPZMatrixSolver<TVar> > solver)
{
	fK00=solver->Matrix();
	fSolver = solver;
    this->fSymProp = fK00->GetSymmetry();
}


template<class TVar,class TSideMatrix>
void
TPZMatRed<TVar, TSideMatrix>::SetK00(TPZAutoPointer<TPZMatrix<TVar> > K00)
{
	fK00=K00;
  this->fSymProp = fK00->GetSymmetry();
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar,TSideMatrix>::SetF(const TPZFMatrix<TVar> & F)
{
    if(F.Rows() != fDim0+fDim1){
        DebugStop();
    }
	
	int64_t FCols=F.Cols(),c,r,r1;
	
	fF0.Redim(fDim0,FCols);
	fF1.Redim(fDim1,FCols);
    fF0IsComputed = false;
	
	for(c=0; c<FCols; c++){
		r1=0;
		for(r=0; r<fDim0; r++){
			fF0.PutVal( r,c,F.GetVal(r,c) ) ;
		}
		//aqui r=fDim0
		for( ;r<fDim0+fDim1; r++){
			fF1.PutVal( r1++,c,F.GetVal(r,c) );
		}
	}
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        F.Print("F Input",sout);
        fF0.Print("fF0 Initialized",sout);
        fF1.Print("fF1 Initialized",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::F1Red(TPZFMatrix<TVar> &F1Red)
{
	if (!fDim0)
    {
        F1Red = fF1;
        return;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "fF0 input " << std::endl;
        fF0.Print("fF0",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    F1Red.Resize(fK11.Rows(),fF0.Cols());
    if (!fF0IsComputed)
    {
        DecomposeK00();
        fSolver->Solve(fF0,fF0);
        fF0IsComputed = true;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "After computing F0Invert" << std::endl;
        fF0.Print("F0Invert",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Input fF1" << std::endl;
        fF1.Print("fF1",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	//make [F1]=[F1]-[K10][F0Invert]
	fK10.MultAdd((fF0),fF1,(F1Red),-1,1);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        F1Red.Print("F1 Reduced", sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
	return;
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar,TSideMatrix>::K11Reduced(TPZFMatrix<TVar> &K11, TPZFMatrix<TVar> &F1)
{
	
	if(!fK01IsComputed)
	{
        DecomposeK00();
		SimetrizeMatRed();
		fSolver->Solve(fK01,fK01);
        TPZStepSolver<TVar> *step = dynamic_cast<TPZStepSolver<TVar> *>(fSolver.operator->());
        if (step->Singular().size())
        {
            std::cout << "Address " << (void *) step << " Number of singular modes " << step->Singular().size() << std::endl;
        }
		fK01IsComputed = true;
	}
	fK10.MultAdd(fK01,fK11,(K11),-1.,1.);
    F1Red(F1);

	return;
}

#include "tpzverysparsematrix.h"

template<>
void TPZMatRed<double, TPZVerySparseMatrix<double> >::K11Reduced(TPZFMatrix<double> &K11, TPZFMatrix<double> &F1)
{
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";	
	static TPZFMatrix<double> temp;
	return;
}

template<>
void TPZMatRed<float, TPZVerySparseMatrix<float> >::K11Reduced(TPZFMatrix<float> &K11, TPZFMatrix<float> &F1)
{
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";	
	static TPZFMatrix<float> temp;
	return;
}

template<>
void TPZMatRed<long double, TPZVerySparseMatrix<long double> >::K11Reduced(TPZFMatrix<long double> &K11, TPZFMatrix<long double> &F1)
{
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";	
	static TPZFMatrix<long double> temp;
	return;
}

template<>
void TPZMatRed<std::complex<double>, TPZVerySparseMatrix<std::complex<double> > >::K11Reduced(TPZFMatrix<std::complex<double> > &K11, TPZFMatrix<std::complex<double> > &F1)
{
	std::cout << __PRETTY_FUNCTION__ << " should never be called\n";	
	static TPZFMatrix<std::complex<double> > temp;
	return;
}


template<class TVar ,class TSideMatrix>
void TPZMatRed<TVar,TSideMatrix>::U1(TPZFMatrix<TVar> & F)
{
	TPZFNMatrix<1000,TVar> K1Red(fDim1,fDim1), F1Red(fDim1,fF1.Cols());
	K11Reduced(K1Red, F1Red);
	F=(F1Red);
	K1Red.SolveDirect( F ,ELU);
	
	
}

template<>
void TPZMatRed<REAL, TPZVerySparseMatrix<REAL> >::UGlobal(const TPZFMatrix<REAL> & U1, TPZFMatrix<REAL> & result)
{
	//[u0]=[A00^-1][F0]-[A00^-1][A01]
    //compute [F0]=[A00^-1][F0]
    if (!fF0IsComputed)
    {
        DecomposeK00();
        fSolver->Solve(fF0,fF0);
        fF0IsComputed = true;
    }
	if(!fK01IsComputed)
	{
		TPZFMatrix<REAL> k01(fK01);
        DecomposeK00();
		fSolver->Solve(k01,k01);
		fK01 = k01;
		fK01IsComputed = 1;
	}
	
	//make [u0]=[F0]-[U1]
	TPZFMatrix<REAL> u0( fF0.Rows() , fF0.Cols() );
	fK01.MultAdd(U1,(fF0),u0,-1,0);
	
	result.Redim( fDim0+fDim1,fF0.Cols() );
	int64_t c,r,r1;
	
	for(c=0; c<fF0.Cols(); c++)
	{
		r1=0;
		for(r=0; r<fDim0; r++)
		{
			result.PutVal( r,c,u0.GetVal(r,c) ) ;
		}
		//aqui r=fDim0
		for( ;r<fDim0+fDim1; r++)
		{
			result.PutVal( r,c,U1.GetVal(r1++,c) );
		}
	}
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix >::UGlobal(const TPZFMatrix<TVar> & U1, TPZFMatrix<TVar> & result)
{
	TPZFMatrix<TVar> u0( fF0.Rows() , fF0.Cols() );
	
	if(fK01IsComputed)
	{
		//[u0]=[A00^-1][F0]-[A00^-1][A01]
        //compute [F0]=[A00^-1][F0]
        if(!fF0IsComputed)
        {
            DecomposeK00();
            fSolver->Solve(fF0,fF0);
            fF0IsComputed = true;
        }
		//make [u0]=[F0]-[U1]
		fK01.MultAdd(U1,(fF0),u0,-1,1);
	} else {
        TPZFMatrix<TVar> K01U1(fK01.Rows(),U1.Cols(),0.);
        fK01.MultAdd(U1,fF0,K01U1,-1.,1.);
        DecomposeK00();
        fSolver->Solve(K01U1, u0);
	}
	
	//compute result
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
        U1.Print("U1 = ",sout,EMathematicaInput);
		fF0.Print("fF0 ",sout,EMathematicaInput);
		u0.Print("u0 " ,sout,EMathematicaInput);
		LOGPZ_DEBUG(logger,sout.str())   
		
	}
#endif
	
	result.Redim( fDim0+fDim1,fF0.Cols() );
	int64_t c,r,r1;
	
	for(c=0; c<fF0.Cols(); c++)
	{
		r1=0;
		for(r=0; r<fDim0; r++)
		{
			result.PutVal( r,c,u0.GetVal(r,c) ) ;
		}
		//aqui r=fDim0
		for( ;r<fDim0+fDim1; r++)
		{
			result.PutVal( r,c,U1.GetVal(r1++,c) );
		}
	}
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::UGlobal2(TPZFMatrix<TVar> & U1, TPZFMatrix<TVar> & result)
{
	TPZFMatrix<TVar> u0( fF0.Rows() , fF0.Cols() );
	
	if(fK01IsComputed)
	{
		//[u0]=[A00^-1][F0]-[A00^-1][A01][u1]
			//compute [F0]=[A00^-1][F0]
        if(!fF0IsComputed)
        {
            DecomposeK00();
			fSolver->Solve(fF0,fF0);
            fF0IsComputed = true;
        }
		//make [u0]=[F0]-[U1]
		fK01.MultAdd(U1,(fF0),u0,-1,1);
	} else {
        TPZFMatrix<TVar> K01U1(fK01.Rows(),U1.Cols(),0.);
        fK01.MultAdd(U1,fF0,K01U1,-1.,1.);
        DecomposeK00();
        fSolver->Solve(K01U1, u0);
	}
	
	//compute result
#ifdef PZ_LOG
	if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		fF0.Print("fF0 ",sout);
		u0.Print("u0 " ,sout);
		LOGPZ_DEBUG(logger,sout.str())   
		
	}
#endif
	
	result.Redim( fDim0+fDim1,fF0.Cols() );
    int64_t nglob = fDim0+fDim1;
	
	int64_t c,r,r1;
	
	for(c=0; c<fF0.Cols(); c++)
	{
		r1=0;
		for(r=0; r<fDim0; r++)
		{
			result(r,c) = u0(r,c) ;
		}
		//aqui r=fDim0
		for( ;r<nglob; r++)
		{
			result(r,c) = U1(r1++,c);
		}
	}
}


template<class TVar,class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::Print(const char *name , std::ostream &out ,const MatrixOutputFormat form ) const
{	
	if(form != EInputFormat) {
		out << "Writing matrix 'TPZMatRed<TSideMatrix>::" << name;
		out << "' (" << this->Dim() << " x " << this->Dim() << "):\n";
        out << "fIsReduced " << this->fIsReduced << std::endl;
        out << "fF0IsComputed " << this->fF0IsComputed << std::endl;
        out << "fK01IsComputed " << this->fK01IsComputed << std::endl;
        out << "fNumberRigidBodyModes " << this->fNumberRigidBodyModes << std::endl;
		out << std::endl;
		fK00->Print("K00 =",out,form);
		fK01.Print("K01 = ",out,form);
		fK10.Print("K10 = ",out,form);
		fK11.Print("K11 = ",out,form);
		
		
		fF0.Print("F0 = ",out,form);
		fF1.Print("F1 = ",out,form);
		
        out << "Matrix norms K00 " << Norm(*fK00.operator->()) << " K01 " << Norm(fK01) << " K10 " << Norm(fK10) << " K11 " << Norm(fK11);
		out << "\n\n";
	} else {
		TPZMatrix<TVar>::Print(name,out,form);
	}	
}

template<class TVar, class TSideMatrix>
int TPZMatRed<TVar,TSideMatrix>::Redim(int64_t dim, int64_t dim00){
	if(dim<dim00) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"dim k00> dim");
	if(fK00) fK00->Redim(dim00,dim00);
	
	fDim0=dim00;
	fDim1=dim-dim00;
	
	fK01.Redim(fDim0,fDim1);
	fK10.Redim(fDim1,fDim0);
	fK11.Redim(fDim1,fDim1);
	
	fF0.Redim(fDim0,1);
	fF1.Redim(fDim1,1);
	this->fRow = dim;
	this->fCol = dim;
    fIsReduced = false;
    fK01IsComputed = false;
    fF0IsComputed = false;

	return 0;
}


template<class TVar, class TSideMatrix>
int TPZMatRed<TVar, TSideMatrix>::Zero(){
	if(fK00) fK00->Zero();
    fIsReduced = false;
    fK01IsComputed = false;
    fF0IsComputed = false;
	fK01.Zero();
	fK10.Zero();
	fK11.Zero();
	fF0.Zero();
	fF1.Zero();
    TPZMatrix<TVar>::Redim(fDim0+fDim1,fDim0+fDim1);
	return 0;
}

// z = alpha A^opt x + beta y

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::MultAdd(const TPZFMatrix<TVar> &x,
										   const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
										   const TVar alpha,const TVar beta,
										   const int opt) const
{
	// #warning Not functional yet. Still need to Identify all the variables	
	if(!fIsReduced)
	{
		LOGPZ_ERROR(logger,"TPZMatRed not reduced, expect trouble")
		TPZMatrix<TVar>::MultAdd(x,y,z,alpha,beta,opt);
		return;
	}
	this->MultAddChecks(x,y,z,alpha,beta,opt);
	this->PrepareZ(y,z,beta,opt);
	
	if(!opt)
	{
		if(fK01IsComputed)
		{
			DebugStop();
		}
		
		TPZFMatrix<TVar> l_Res(fK01.Rows(), x.Cols(), 0);
		fK01.Multiply(x,l_Res,0);
		//fSolver->Solve(l_Res,l_Res);
#ifdef PZ_LOG
		if(logger.isDebugEnabled())
		{
			std::stringstream sout;
			l_Res.Print("Internal solution",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		TPZFMatrix<TVar> l_ResFinal(fK11.Rows(), x.Cols(), 0);
		fK11.Multiply(x,l_ResFinal,0);
#ifdef PZ_LOG
		if(logger.isDebugEnabled())
		{
			std::stringstream sout;
			l_ResFinal.Print("Intermediate product l_ResFinal",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fK10.MultAdd(l_Res,l_ResFinal,z,-alpha,alpha,opt);
#ifdef PZ_LOG
		if(logger.isDebugEnabled())
		{
			std::stringstream sout;
			z.Print("Final result z ",sout);
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
	}
	else // not implemented for the transpose case
	{
		DebugStop();
	}
}

/** @brief Decompose K00 and adjust K01 and K10 to reflect rigid body modes */
template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::DecomposeK00()
{
    if(fK00->IsDecomposed())
    {
        return;
    }
    TPZStepSolver<TVar> *stepsolve = dynamic_cast<TPZStepSolver<TVar> *>(fSolver.operator->());
    TPZStepSolver<TVar> *directsolve(0);
    if(!stepsolve)
    {
        DebugStop();
    }
    if(stepsolve->Solver() == TPZMatrixSolver<TVar>::EDirect)
    {
        directsolve = stepsolve;
    }
    if(!directsolve)
    {
        TPZMatrixSolver<TVar> *presolve = stepsolve->PreConditioner();
        TPZStepSolver<TVar> *prestep = dynamic_cast<TPZStepSolver<TVar> *>(presolve);
        if(prestep->Solver() == TPZMatrixSolver<TVar>::EDirect)
        {
            prestep->UpdateFrom(stepsolve->Matrix());
            directsolve = prestep;
        }
    }
    if (directsolve)
    {
        directsolve->Decompose();
        std::list<int64_t> &singular = directsolve->Singular();
        std::list<int64_t>::iterator it;
        int nsing = singular.size();
        if(nsing > fMaxRigidBodyModes-fNumberRigidBodyModes)
        {
            std::cout << "Number of rigid body modes larger than provision\n";
            std::cout << "Number of singular modes " << nsing << std::endl;
            std::cout << "Number of rigid body modes reserved " << fMaxRigidBodyModes << std::endl;
            std::cout << "Rigid body modes ";
            for (it=singular.begin(); it != singular.end(); it++) {
                std::cout << " " << *it;
            }
            std::cout << std::endl;
            //DebugStop();
        }
        for (it=singular.begin(); it != singular.end(); it++) {
            if(fNumberRigidBodyModes < fMaxRigidBodyModes)
            {
                fK01(*it,fDim1-fMaxRigidBodyModes+fNumberRigidBodyModes) = -1.;
                fK10(fDim1-fMaxRigidBodyModes+fNumberRigidBodyModes,*it) = -1.;
                fK11(fDim1-fMaxRigidBodyModes+fNumberRigidBodyModes,fDim1-fMaxRigidBodyModes+fNumberRigidBodyModes) = 1.;
                if(stepsolve != directsolve)
                {
                    TVar diag = stepsolve->Matrix()->GetVal(*it, *it)+ (TVar)1.;
                    stepsolve->Matrix()->PutVal(*it, *it, diag);
                }
            }
            fNumberRigidBodyModes++;
        }
    }
    else
    {
        DebugStop();
    }
}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::Write(TPZStream &buf, int withclassid) const {
    TPZMatrix<TVar>::Write(buf, withclassid);
    {//Ints
        buf.Write(&this->fDim0, 1);
        buf.Write(&this->fDim1, 1);
    }
    {//chars
        buf.Write(this->fIsReduced);
        buf.Write(this->fK01IsComputed);
        buf.Write(&this->fMaxRigidBodyModes, 1);
        buf.Write(&this->fNumberRigidBodyModes, 1);
    }
    {//Aggregates
        this->fF0.Write(buf, 0);
        this->fF1.Write(buf, 0);
        TPZPersistenceManager::WritePointer(this->fK00.operator ->(), &buf);
        this->fK01.Write(buf, 0);
        this->fK10.Write(buf, 0);
        this->fK11.Write(buf, 0);
        if (fSolver) {
            if (fSolver->Matrix() != fK00) {
                std::cout << "Error\n";
            } else {
                TPZPersistenceManager::WritePointer(fSolver.operator ->(), &buf);
                //TODO Enviar o solver, atenção com a Matrix do Solver;
            }

        } else {
            int flag = -1;
            buf.Write(&flag, 1);
        }

    }

}

template<class TVar, class TSideMatrix>
void TPZMatRed<TVar, TSideMatrix>::Read(TPZStream &buf, void *context) {
    TPZMatrix<TVar>::Read(buf, context);
    {//Ints
        buf.Read(&this->fDim0, 1);
        buf.Read(&this->fDim1, 1);
    }
    {//chars
        buf.Read(this->fIsReduced);
        buf.Read(this->fK01IsComputed);
        buf.Read(&this->fMaxRigidBodyModes, 1);
        buf.Read(&this->fNumberRigidBodyModes, 1);
    }
    {//Aggregates
        this->fF0.Read(buf, 0);
        this->fF1.Read(buf, 0);
        TPZAutoPointer<TPZSavable> sav = TPZPersistenceManager::GetAutoPointer(&buf);
        TPZAutoPointer<TPZMatrix<TVar>> mat = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(sav);
        if (sav && !mat) {
            DebugStop();
        }
        fK00 = mat;
        this->fK01.Read(buf, 0);
        this->fK10.Read(buf, 0);
        this->fK11.Read(buf, 0);
        sav = TPZPersistenceManager::GetAutoPointer(&buf);
        TPZAutoPointer<TPZMatrixSolver<TVar>>matsolv = TPZAutoPointerDynamicCast<TPZMatrixSolver<TVar>> (sav);
        if (sav && !matsolv) {
            DebugStop();
        }
        if (matsolv) {
            fSolver = matsolv;
        }
    }
}

#include "tpzverysparsematrix.h"

template class TPZMatRed<float, TPZVerySparseMatrix<float> >;
template class TPZMatRed<float, TPZFMatrix<float> >;

template class TPZMatRed<double, TPZVerySparseMatrix<double> >;
template class TPZMatRed<double, TPZFMatrix<double> >;

template class TPZMatRed<long double, TPZVerySparseMatrix<long double> >;
template class TPZMatRed<long double, TPZFMatrix<long double> >;

template class TPZMatRed<std::complex<double>, TPZVerySparseMatrix<std::complex<double> > >;

template class TPZMatRed<std::complex<float>, TPZFMatrix<std::complex<float> > >;
template class TPZMatRed<std::complex<double>, TPZFMatrix<std::complex<double> > >;
template class TPZMatRed<std::complex<long double>, TPZFMatrix<std::complex<long double> > >;


#ifndef BORLAND
template class TPZRestoreClass<TPZMatRed<REAL,TPZVerySparseMatrix<REAL> > >;
template class TPZRestoreClass<TPZMatRed<REAL, TPZFMatrix<REAL> > >;
#endif
