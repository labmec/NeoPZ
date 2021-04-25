/**
 * @file
 * @brief Contains TPZAnalysis class which implements the sequence of actions to perform a finite element analysis.
 */

#ifndef ANALYSISH
#define ANALYSISH

#include <mutex>              // for mutex
#include <set>                // for set
#include <functional>         // for function
#include "TPZGuiInterface.h"  // for TPZGuiInterface
#include "pzerror.h"          // for DebugStop
#include "pzmatrix.h"         // for TPZFMatrix, TPZMatrix
#include "TPZSolutionMatrix.h"// for TPZSolutionMatrix
#include "pzreal.h"           // for STATE, REAL
#include "TPZRenumbering.h"    // for TPZRenumbering
#include "TPZStructMatrix.h"      // for TPZStructMatrix
#include "pzvec.h"            // for TPZVec
#include "tpzautopointer.h"   // for TPZAutoPointer
#include "pzadmchunk.h"       // for TPZAdmChunkVector
class TPZCompEl;
class TPZCompMesh;
class TPZConnect;
class TPZGeoMesh;
class TPZGraphMesh;
template <class TVar> class TPZMatrixSolver;

/**
 * @ingroup analysis
 * @brief Implements the sequence of actions to perform a finite element analysis. \ref analysis "Analysis"
 */
/** This class will renumerate the nodes upon construction
 */
class TPZAnalysis : public TPZSavable {
	
public:
	
	/** @brief  Preconditioners which can be created by objects of this class */
	enum EPrecond { EJacobi, EBlockJacobi, EElement, ENodeCentered };
	
	
protected:
	/** @brief Geometric Mesh */
	TPZGeoMesh *fGeoMesh;
	/** @brief Computational mesh */
	TPZCompMesh *fCompMesh;
	/** @brief Graphical mesh */
	TPZGraphMesh *fGraphMesh[3];
	/** @brief Load vector */
	TPZSolutionMatrix fRhs;
	/** @brief Solution vector */
	TPZSolutionMatrix fSolution;
	/** @brief Type of solver to be applied*/
	TPZSolver *fSolver;
	/** @brief Scalar variables names - to post process */
	TPZVec<std::string> fScalarNames[3];
	/** @brief Vectorial variables names - to post process */
	TPZVec<std::string> fVectorNames[3];
    /** @brief Tensorial variables names - to post process */
    TPZVec<std::string> fTensorNames[3];
	/** @brief Time step */
	int fStep;
	/** @brief Time variable which is used in dx output */
	REAL fTime;
	
    /** @brief Number of threads to be used for post-processing error */
    int fNthreadsError;
    
	/** @brief Structural matrix */
	TPZAutoPointer<TPZStructMatrix>  fStructMatrix;
	
	/** @brief Renumbering scheme */
	TPZAutoPointer<TPZRenumbering> fRenumber;
	
	/** @brief Pointer for gui interface object */
	TPZAutoPointer<TPZGuiInterface> fGuiInterface;
	
	/** @brief Datastructure which defines postprocessing for one dimensional meshes */
	class TTablePostProcess : public TPZSavable {
        public :
		TPZVec<int64_t> fGeoElId;
		TPZVec<TPZCompEl *> fCompElPtr;
		int fDimension;
		TPZVec<REAL> fLocations;
		TPZVec<std::string> fVariableNames;
		TTablePostProcess();
		~TTablePostProcess();
                
                int ClassId() const override;
                
                void Write(TPZStream &buf, int withclassid) const override;

                void Read(TPZStream &buf, void *context) override;
	};
	
        TTablePostProcess fTable;
  public:
	
	/** @brief Create an TPZAnalysis object from one mesh pointer */
	TPZAnalysis(TPZCompMesh *mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
    	
	/** @brief Create an TPZAnalysis object from one mesh auto pointer object */
	TPZAnalysis(TPZAutoPointer<TPZCompMesh> mesh, bool mustOptimizeBandwidth = true, std::ostream &out = std::cout);
    
	/** @brief Defines gui interface object */
	void SetGuiInterface(TPZAutoPointer<TPZGuiInterface> gui){
		fGuiInterface = gui;
	}
	
	/** @brief Gets gui interface object */
	TPZAutoPointer<TPZGuiInterface> GetGuiInterface() const{
		return fGuiInterface;
	}
	
	/** @brief Returns if the process was canceled through gui interface */
	bool AmIKilled(){
		if(fGuiInterface){
			return fGuiInterface->AmIKilled();
		}
		else return false;
	}
	
	/** @brief Set the computational mesh of the analysis. */
	virtual void SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth);
	
	/** @brief Create an empty TPZAnalysis object */
	TPZAnalysis();
	
	/** @brief Destructor: deletes all protected dynamic allocated objects */
	virtual ~TPZAnalysis(void);
    
    /// deletes all data structures
    void CleanUp();
    
    /// Change the renumbering scheme
    void SetRenumber(TPZAutoPointer<TPZRenumbering> renumber)
    {
        fRenumber = renumber;
    }
    
	/** @brief Sets the computer connection block number from the graphical connections block number otimization */
	void OptimizeBandwidth();
	
	/** @brief Returns the dimension of the material which has the highest dimension */
	int HighestDimension();
	
	/** @brief Recompute the node sequence */
	void Resequence(int firstel = -1);
    
    /** @brief Determine the number of load cases from the material objects and return its value */
    /**
     * this method will modify the material objects so that they have all the same number of load cases
     * the number of load cases is the maximum value of load cases of all material objects
     */
    int ComputeNumberofLoadCases();
    
	
	/** @brief Assemble the stiffness matrix and load vector */
	virtual  void Assemble();
	
	/** @brief Assemble the load vector */
	virtual void AssembleResidual();
	
	/** @brief Invert the stiffness matrix */
	virtual void Solve();
	
	/** @brief Returns the load vector */
	TPZSolutionMatrix &Rhs() { return fRhs;}

	/** @brief Returns the solution matrix */
	TPZSolutionMatrix &Solution() { return fSolution;}
	
	/** @brief Returns the pointer to the computational mesh */
	TPZCompMesh *Mesh()const { return fCompMesh;}
	/** @brief Returns a reference to the structural matrix */
	TPZAutoPointer<TPZStructMatrix> StructMatrix() {
        if(!fStructMatrix)
        {
            DebugStop();
        }
        return fStructMatrix;
    }
	
	/** @brief Define the type of preconditioner used */
	/** This method will create the stiffness matrix but without assembling */
    template<class TVar>
	TPZMatrixSolver<TVar> *BuildPreconditioner(EPrecond preconditioner, bool overlap);
	
    /** @brief ste the step for post processing */
    void SetStep(int step)
    {
        fStep = step;
    }
    
    void SetThreadsForError(int nthreads)
    {
        fNthreadsError = nthreads;
    }
    
    int GetStep()
    {
        return fStep;
    }
    
	/** @brief Sets time will be used in dx files */
	void SetTime(REAL time);
	/** @brief Gets time used in dx files */
	REAL GetTime();

	/** @brief Graphic of the solution as V3DGrap visualization */
	void ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices);
    /** @brief Graphic of the solution as V3DGrap visualization */
    void ShowShape(const std::string &plotfile, TPZVec<int64_t> &equationindices, int matid, const TPZVec<std::string> &varname);
    /** @brief Make assembling and clean the load and solution vectors */
	void LoadShape(double dx,double dy, int64_t numelem,TPZConnect* nod);
	
	/** @brief Calls the appropriate sequence of methods to build a solution or a time stepping sequence */
	virtual void Run(std::ostream &out = std::cout);
	/** @brief Define GrapMesh as V3D, DX, MV or VTK depending on extension of the file */
	virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile);
    /** @brief Define GrapMesh as VTK with tensorial names depending on extension of the file */
    virtual void DefineGraphMesh(int dimension, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile);
    /** @brief Define GrapMesh as V3D, DX, MV or VTK depending on extension of the file */
    virtual void DefineGraphMesh(int dimension, const std::set<int> & matids ,const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const std::string &plotfile);
    /** @brief Define GrapMesh as VTK with tensorial names depending on extension of the file */
    virtual void DefineGraphMesh(int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames, const std::string &plotfile);
    
	/** @brief Clean the GrapMesh vector */
	virtual void CloseGraphMesh();
	
	/** @brief Defines the postprocessing parameters for the graphical grid */
	TPZGraphMesh *GraphMesh(int dimension) {
		return fGraphMesh[dimension-1];
	}
	/** @brief Draw solution over mesh for all dimensions */
	virtual void PostProcess(int resolution);
	/** @brief Draw solution over mesh by dimension  */	
	virtual void PostProcess(int resolution, int dimension);
	
    /** @brief Fill mat ids with materials with provided dimension wich are not boundary conditinos or interface  */
    void IdentifyPostProcessingMatIds(int dimension, std::set<int> & matids);
    
    
	/**
	 * @name Related over data structure to post processing
	 * @{
	 */
	
	/** @brief Fill the computational element vector to post processing depending over geometric mesh defined */
	virtual void DefineElementTable(int dimension, TPZVec<int64_t> &GeoElIds, TPZVec<REAL> &points);
	/** @brief Sets the names of the variables into the data structure for post processing */	
	virtual void SetTableVariableNames(TPZVec<std::string> varnames);
	/** @brief Prepare data to print post processing and print coordinates */
	virtual void PrePostProcessTable(std::ostream &out_file);
	/** @brief Print the solution related with the computational element vector in post process */
	virtual void PostProcessTable(std::ostream &out_file);

    
    /** @brief Compute and print the local error over all elements in data structure of post process, also compute global errors in several norms */
	void PostProcessTable(  TPZFMatrix<REAL> &pos,std::ostream &out= std::cout );

	/** @} */

	/** @brief Load the solution into the computable grid */
	virtual void LoadSolution();
	/** @brief Load the solution into the computable mesh considering sol as Solution vector of the analysis */
	virtual void LoadSolution(const TPZFMatrix<STATE> &sol){
		fSolution = sol;
		this->LoadSolution();
	}

    //! Integrate the postprocessed variable name over the elements included in the set matids
    TPZVec<STATE> Integrate(const std::string &varname, const std::set<int> &matids);
    
	//! Sets an exact solution in all the materials of the associated mesh
    void SetExact(std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> f, int pOrder = 1);
	/** @brief Compute the local error over all elements and global errors in several norms and print out */
	virtual void PostProcess(TPZVec<REAL> &loc, std::ostream &out = std::cout);
    
    /**
     * @brief Compute the local error over all elements and global errors in several norms and print out
     * without calculating the errors of the variables for hdiv spaces.
     */
    virtual void PostProcessError(TPZVec<REAL> &, bool store_error = true, std::ostream &out = std::cout);
    
    virtual void PostProcessErrorSerial(TPZVec<REAL> &, bool store_error = true, std::ostream &out = std::cout);
    
    virtual void PostProcessErrorParallel(TPZVec<REAL> &, bool store_error = true, std::ostream &out = std::cout);
    
    void CreateListOfCompElsToComputeError(TPZAdmChunkVector<TPZCompEl *> &elvec);
	
	/** @brief Print connect and solution information */
	void Print( const std::string &name , std::ostream &out );
    
    /// Print the residual vector for those elements with entry above a given tolerance
    void PrintVectorByElement(std::ostream &out, TPZFMatrix<STATE> &vec, REAL tol = 1.e-10);

    inline TPZSolver * Solver() {return fSolver;}
	/** @brief Get the matrix solver */
    template<class TVar>
	TPZMatrixSolver<TVar> & MatrixSolver();
	/** @brief Run and print the solution step by step */
	void AnimateRun(int64_t num_iter, int steps,
					TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames, const std::string &plotfile);
	/** @brief Set solver matrix */
	void SetSolver(const TPZSolver &solver);
	/** @brief Set structural matrix as auto pointer for analysis */
	void SetStructuralMatrix(TPZAutoPointer<TPZStructMatrix> strmatrix);
	/** @brief Set structural matrix for analysis */	
	void SetStructuralMatrix(TPZStructMatrix &strmatrix);
  
    int ClassId() const override;
    
    void Write(TPZStream &buf, int withclassid) const override;

    void Read(TPZStream &buf, void *context) override;
    
    struct ThreadData {

      TPZAdmChunkVector<TPZCompEl *> fElvec;
    
      ThreadData(TPZAdmChunkVector<TPZCompEl *> &elvec, bool store_error);
    
    ~ThreadData();
    
    static void *ThreadWork(void *threaddata);
    
    int64_t fNextElement;
    
    int ftid;
      
      bool fStoreError = false;
    
    // Vector with errors. Assuming no more than a 100 threads
    TPZManVector<TPZManVector<REAL,10>,100> fvalues;
    
    protected:
    //mutex for assigning thread id
    std::mutex fMutexThreadId;
    //mutex for accessing next element on list
    std::mutex fMutexAccessEl;
    
  };
  
  friend struct ThreadData;

private:
  /** @brief Build a sequence solver based on the block graph and its colors */
  template <class TVar>
  TPZMatrixSolver<TVar> *
  BuildSequenceSolver(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex,
                      int64_t neq, int numcolors, TPZVec<int> &colors);

  template <class TVar>
  void
  SetExactInternal(std::function<void(
               const TPZVec<REAL> &loc, TPZVec<TVar> &result,
               TPZFMatrix<TVar> &deriv)>
               f,int pOrder = 1);
  template <class TVar> void AssembleInternal();
  template <class TVar> void SolveInternal();
  template <class TVar>
  void ShowShapeInternal(const TPZStack<std::string> &scalnames,
                         const TPZStack<std::string> &vecnames,
                         const std::string &plotfile,
                         TPZVec<int64_t> &equationindices);
  template <class TVar>
  void AnimateRunInternal(int64_t num_iter, int steps,
                          TPZVec<std::string> &scalnames,
                          TPZVec<std::string> &vecnames,
                          const std::string &plotfile);
};




inline void TPZAnalysis::SetTime(REAL time){
	this->fTime = time;
}

inline REAL TPZAnalysis::GetTime(){
	return this->fTime;
}

extern template
TPZMatrixSolver<STATE> *TPZAnalysis::BuildPreconditioner<STATE>(
    EPrecond preconditioner,bool overlap);
extern template
TPZMatrixSolver<STATE> &TPZAnalysis::MatrixSolver<STATE>();

#endif
