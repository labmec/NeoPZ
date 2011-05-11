/* Generated by Together */
class TPZEqnArray;

#ifndef TPZFRONTSYM_H
#define TPZFRONTSYM_H

#ifdef USING_ATLAS
extern "C"{
     #include <cblas.h>
     };
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>


#include "pzmatrix.h"
#include "pzstack.h"
#include "pzvec.h"
#include "TPZFront.h"
#include "TPZFileEqnStorage.h"
#include "TPZStackEqnStorage.h"

/** 
 * The Front matrix itself. \n
 * It is controled by TPZFrontMatrix.\n
 * TPZFrontSym is a symmetrical matrix. \n
 * It uses a Cholesky decomposition scheme.
 * @ingroup frontal
 */
/** Jorge
 * @brief Abstract class implements storage and decomposition process of the frontal matrix involving simetry characteristics
 */

class TPZFrontSym : public TPZFront {
public:
     /**Returns its type*/
     std::string GetMatrixType();

    /** Static main used for testing */
#ifndef WIN32
	 static void main();
#endif
    /** Simple destructor */
    ~TPZFrontSym();
    /** Simple constructor */
    TPZFrontSym();
    
    TPZFrontSym(const TPZFrontSym &cp) : TPZFront(cp),
    fDecomposeType(cp.fDecomposeType)
    {
    }
    /** Constructor with a initial size parameter */
	TPZFrontSym(int GlobalSize);

    /**
     * Decompose these equations and put the result in eqnarray
     * Default decompose method is Cholesky
     */
    void DecomposeEquations(
			    int mineq //! Starting index of equations to be decomposed
			    , int maxeq //! Finishing index of equations to be decomposed
			    , TPZEqnArray & result //! Result of decomposition
			    );

    /**
     * Decompose these equations in a symbolic way and store freed indexes in fFree 
     */
    void SymbolicDecomposeEquations(
          int mineq //! Initial equation index
          , int maxeq //! Final equation index
          );
	
	/** Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
	void SymbolicAddKel(TPZVec < int > & destinationindex);

    /**
     * Compress data structure 
     */
    void Compress();

	/**
	* Expand the front matrix
	*/
	void Expand(int largefrontsize);

    /**
     * Returns ith, jth element of matrix.
     * @associates <{mat(sourceindex[i],sourceindex[j])}>
     * @semantics += 
     */
     //REAL & Element(int i, int j);
     REAL & Element(int i, int j){
	if(i>j){
		int i_temp=i;
		i=j;
		j=i_temp;
		//cout << "Changing row column indexes !"<<endl;
	}
	return fData[(j*(j+1))/2+i];
}
    /**
     * Returns ith, jth element of matrix.
     * @associates <{mat(sourceindex[i],sourceindex[j])}>
     * @semantics += 
     */
    //REAL & Element(int i, int j);
    const REAL & Element(int i, int j) const {
        if(i>j){
            int i_temp=i;
            i=j;
            j=i_temp;
            //cout << "Changing row column indexes !"<<endl;
        }
        return fData[(j*(j+1))/2+i];
    }
    /**Add a contribution of a stiffness matrix*/
    void AddKel(TPZFMatrix &elmat, TPZVec<int> &destinationindex);

    /**Add a contribution of a stiffness matrix*/
    void AddKel(TPZFMatrix &elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex);    

	/** Reorders the elements of the frontmatrix into the full matrix
	*/
	virtual void ExtractFrontMatrix(TPZFMatrix &front);

private:    

    /**
     * Decomposes ieq equation and add the result to EqnArray 
     */
    void DecomposeOneEquation(
			      int ieq //! index of equation to be decomposed
			      , TPZEqnArray &eqnarray //! EqnArray to store resulting members
			      );
    /**
     * Sets the global equation as freed, allowing the space 
     * used by this equation to be used 
     * by future assembly processes 
     */
    void FreeGlobal(int global);
    /**
     * return a local index corresponding to a global equation number 
     */
    int Local(int global);
public:
    /** Returns the number of free equations */
	virtual int NFree();
    /** Resets data structure */
	void Reset(int GlobalSize=0);
    /** Allocates data for Front */
	void AllocData();

     /**
      * It prints TPZFront data 
      */
     void Print(const char *name, std::ostream& out=std::cout) const;
     void PrintGlobal(const char *name, std::ostream& out = std::cout);

     /**Returns decomposition type*/
     DecomposeType GetDecomposeType() const;

private:

    /**
     * Used Decomposition method 
     */
    DecomposeType fDecomposeType;

    /** @link dependency */
    /*#  TPZFileEqnStorage lnkTPZFileEqnStorage; */

    /** @link dependency */
    /*#  TPZStackEqnStorage lnkTPZStackEqnStorage; */
};
#endif //TPZFRONTSYM_H
