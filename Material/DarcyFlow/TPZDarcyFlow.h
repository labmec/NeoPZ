//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZDARCYFLOW_H
#define TPZDARCYFLOW_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZIsotropicPermeability.h"

/**
 * @ingroup material
 * @brief This class implements an H1-conforming approximation for the Darcy flow equation for isotropic materials.
 *
 * The Darcy flow equation is given by: \f[- \nabla \cdot K \nabla u = f,\f] where \f$u\f$ is the pressure field
 * to be solved, \f$K\f$ is the permeability tensor and \f$f\f$ is the source term.
 *
 * @see TPZMixedDarcyFlow For an approximation using the mixed method.
 * @see TPZIsotropicPermeability For setting the permeability field.
 */

class TPZDarcyFlow : public TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>,
        TPZMatErrorSingleSpace<STATE>, TPZIsotropicPermeability> {

    // type alias to improve constructor readability
    using TBase = TPZMatBase<STATE, TPZMatSingleSpaceT<STATE>, TPZMatErrorSingleSpace<STATE>, TPZIsotropicPermeability>;

public:
    /**
     * @brief Default constructor
     */
    TPZDarcyFlow();

    /**
	 * @brief Class constructor
	 * @param [in] id material id
	 * @param [in] dim problem dimension
	 */
    TPZDarcyFlow(int id, int dim);

            TPZDarcyFlow(const TPZDarcyFlow &copy);

            TPZDarcyFlow& operator=(const TPZDarcyFlow &copy);

    /**
	 * @brief Returns a 'std::string' with the name of the material
	 */
    [[nodiscard]] std::string Name() const override { return "TPZDarcyFlow"; }

    /**
	 * @brief Returns the problem dimension
	 */
    [[nodiscard]] int Dimension() const override { return this->fDim; }

    /**
	 * @brief Returns the number of state variables
	 */
    [[nodiscard]] int NStateVariables() const override { return 1; }

    /**
	 * @brief Returns the number of errors to be evaluated
     *
     * Returns the number of errors to be evaluated, that is, the number of error norms associated
     * with the problem.
     */
    int NEvalErrors() const override { return 3; }

    /**
     * @brief Sets problem dimension
     */
    virtual void SetDimension(int dim);

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     */
    void Contribute(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek,
                    TPZFMatrix<STATE> &ef) override;

    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the element matrix
     * @param[out] ef is the rhs vector
     * @param[in] bc is the boundary condition material
     */
    void ContributeBC(const TPZMaterialDataT<STATE> &data, STATE weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                      TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Returns an integer associated with a post-processing variable name
     * @param [in] name string containing the name of the post-processing variable. Ex: "Pressure".
     */
    [[nodiscard]] int VariableIndex(const std::string &name) const override;

    /**
     * @brief Returns an integer with the dimension of a post-processing variable
     * @param [in] var index of the post-processing variable, according to TPZDarcyFlow::VariableIndex method.
     */
    [[nodiscard]] int NSolutionVariables(int var) const override;

    /**
     * @brief Returns the solution associated with the var index based on the
     * finite element approximation at a point
     * @param [in] data material data associated with a given integration point
     * @param [in] var index of the variable to be calculated
     * @param [out] solOut vector to store the solution
     */
    void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &solOut) override;

    /**
     * @brief Get the dimensions of the solution for each state variable
     *
     * This will be used for initializing the corresponding TPZMaterialData
     * @param [out] u_len solution dimension
     * @param [out] du_row number of rows of the derivative
     * @param [out] du_col number of columns of the derivative
     */
    void GetSolDimensions(uint64_t &u_len, uint64_t &du_row, uint64_t &du_col) const override;

    /**
     * @brief Calculates the approximation error at a point
     * @param [in] data material data of the integration point
     * @param [out] errors calculated errors
     */
    void Errors(const TPZMaterialDataT<STATE> &data, TPZVec<REAL> &errors) override;

    /*
     * @brief fill requirements for volumetric contribute
     */
    void FillDataRequirements(TPZMaterialData &data) const override;

    /*
     * @brief fill requirements for boundary contribute
     */
    void FillBoundaryConditionDataRequirements(int type, TPZMaterialData &data) const override;

    /**
     * @brief Returns an unique class identifier
     */
    [[nodiscard]] int ClassId() const override;

    /**
     * @brief Creates another material of the same type
     */
    [[nodiscard]] TPZMaterial *NewMaterial() const override;

    /**
     * @brief Prints data associated with the material.
     */
    void Print(std::ostream & out) const override;

protected:
    /**
     * @brief Problem dimension
     */
    int fDim;
};

#endif //TPZDARCYFLOW_H
