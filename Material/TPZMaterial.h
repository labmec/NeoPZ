/**
 * @file TPZMaterial.h
 * @brief Header file for abstract class TPZMaterial.\n
 * It implements a type-agnostic interface for implementing the weak statement of the differential equation within the PZ environment.
 */

#ifndef TPZMATERIAL_H
#define TPZMATERIAL_H
#include "TPZSavable.h"
#include "TPZMatTypes.h"

class TPZStream;
class TPZMaterialData;
template<class T>
class TPZVec;
template<class T>
class TPZFMatrix;

/**
 * @ingroup material
 * @brief This abstract class is the basis of the TPZMaterial hierarchy.
 * It contains the type-agnostic interface that materials should implement,
 * and it applies for materials using both single approximation space
 * and combined approximation spaces (multiphysics materials).
 * @note Actual materials should derive from TPZMatBase
 */
class TPZMaterial : public virtual TPZSavable{
public:
    /** @name Constructors */
    /** @{ */
    //! Default constructor.
    TPZMaterial() = default;
    //! Constructor taking material identifier.
    explicit TPZMaterial(int id) : fId(id){
    }
    //! Copy constructor.
    TPZMaterial(const TPZMaterial&) = default;
    //! Move constructor.
    TPZMaterial(TPZMaterial&&) = default;
    /** @} */
    //! Destructor.
    virtual ~TPZMaterial() = default;
    //! Copy assignment operator
    TPZMaterial& operator=(const TPZMaterial&) = default;
    //! Move assignment operator
    TPZMaterial& operator=(TPZMaterial&&) = default;
    /** @name Utilities*/
    /** @{ */
    //! Returns the name of the material.
    [[nodiscard]] virtual std::string Name() const { return "no_name";}
    //! Returns the integrable dimension of the material.
    [[nodiscard]] virtual int Dimension() const = 0;
    //! Returns the material identifier.
    [[nodiscard]] virtual int Id() const {return fId;}
    //! Sets the material identifier.
    void SetId(int id) {fId = id;}
    //! Whether a forcing function has been set.
    [[nodiscard]] virtual bool HasForcingFunction() const = 0;
    //! Prints data associated with the material.
    virtual void Print(std::ostream &out = std::cout) const;

     /** @brief Create another material of the same type
      */
    virtual TPZMaterial * NewMaterial() const;
    /**
        @brief Change the value of the penalty constant
     */
    void SetBigNumber(REAL bignumber){
        fBigNumber = bignumber;
    }
    /**
        @brief Access method to the penaly constant
     */
    REAL BigNumber() {
        return fBigNumber;
    }
    /** @} */

    /** @name PostProcess
     * @{
     */
    //! Returns the number of state variables associated with the material.
    [[nodiscard]] virtual int NStateVariables() const = 0;
    
    /** @brief Returns the variable index associated with a given name */
    [[nodiscard]] virtual int VariableIndex(const std::string &name) const;

    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    [[nodiscard]] virtual int NSolutionVariables(int var) const;



    /**
     * @}
     */
    
    /** @brief Creates a copy of the material object and put it in the vector which is passed on */
    virtual void Clone(std::map<int, TPZMaterial * > &matvec);

    [[nodiscard]] int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override;

    void Write(TPZStream& buf, int withclassid) const override;

protected:
    /*! 
      Dummy function so we can distinguish the classes 
      materials should inherit from.
    */

    /*
      Do NOT implement this function in your material. Instead,
      inherit your material from
      TPZMatBase
     */
    virtual bool IsMatImpl() = 0;
    /**
       @brief Big number for penalization method and Dirichlet BCs.
       The default value should be sufficiently big for most cases.
       However, you may adapt it to your needs
    */
    REAL fBigNumber{pow(10,std::numeric_limits<STATE>::max_digits10*2/3)};
    //using float because it is already big enough
private:
    //! Material identifier.
    int fId{-666};
};
#endif
