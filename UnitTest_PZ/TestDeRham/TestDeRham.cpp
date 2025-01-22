/**
 * @file
 * @brief Contains Unit Tests for checking compatibility of approx spaces
 */

#include "TPZMatDeRhamH1.h"
#include "TPZMatDeRhamHCurl.h"
#include "TPZMatDeRhamHDiv.h"
#include "TPZMatDeRhamL2.h"
#include "TPZMatDeRhamH1HCurl.h"
#include "TPZMatDeRhamHCurlHDiv.h"
#include "TPZMatDeRhamHDivL2.h"
#include "TPZNullMaterial.h"
#include "TPZGeoMeshTools.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZLinearAnalysis.h"
#include "pzfstrmatrix.h"
#include "TPZMatrixSolver.h"
#include "TPZElementMatrixT.h"
#include "MMeshType.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzshapetetra.h"
#include "pzshapecube.h"
#include "TPZShapeHDivConstant.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_approx.hpp>

enum class ESpace
{
  H1,
  HCurl,
  HDiv,
  HDivConst,
  L2
};

static const std::map<ESpace, const char *> names({{ESpace::H1, "H1"},
                                                   {ESpace::HCurl, "HCurl"},
                                                   {ESpace::HDiv, "HDiv"},
                                                   {ESpace::HDivConst, "HDivConst"},
                                                   {ESpace::L2, "L2"}});

/** @brief Compares the dimensions of the span of the differential operator associated
    with the left space against the kernel of the right space.
    Let us call the operator op (h1-> grad, hcurl -> curl, hdiv -> div, l2 -> id).
    This function creates the matrix op(phi_i,phi_j) and compares if
    rank(M_left) = ker(M_right).
    @note For L2, the comparison is rank(M_left) = rank(M_right), and the operator
    associated with the L2 approximation space is the identity (op(phi_i) = phi_i).
    @param [in] kRight polynomial order of the rightmost approximation space (except for HDiv and HDivConst)*/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckRankKerDim(int kRight);

/** @brief Checks if the range of the differential operator of the left approximation space is contained in the functions of the approximation space on the right.
    Let us call the operator op (h1-> grad, hcurl -> curl, hdiv -> div, l2 -> id).
    Denoting by
    - phi_i : basis functions of the left space
    - op: operator of the left space
    - varphi_i : basis functions of the right space
    This function creates the matrix
    M = [ A B ] = [ op(phi_i).op(phi_j)  op(phi_i). varphi_j ]
        [ C D ]   [ varphi_i .op(phi_j)  varphi_i . varphi_j ]

    And then compare if rank(M) = rank(D).
    @param [in] kRight polynomial order of the rightmost approximation space*/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckInclusion(int kRight);

TEMPLATE_TEST_CASE("Dimension Compatibility", "[derham_tests]",
                   (typename std::integral_constant<int, 2>),
                   (typename std::integral_constant<int, 3>))
{
  // SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  constexpr int dim = TestType::value;

  ESpace leftSpace = GENERATE(ESpace::H1,
                              ESpace::HCurl,
                              ESpace::HDiv,
                              ESpace::HDivConst);

  SECTION(names.at(leftSpace))
  {
    int k = GENERATE(1, 2, 3);
    SECTION(std::to_string(k))
    {
      switch (leftSpace)
      {
      case ESpace::H1:
      {
        CheckRankKerDim<ESpace::H1, ESpace::HCurl, dim>(k);
        break;
      }
      case ESpace::HCurl:
      {
        if constexpr (dim == 2)
        {
          CheckRankKerDim<ESpace::HCurl, ESpace::L2, dim>(k);
        }
        else
        {
          CheckRankKerDim<ESpace::HCurl, ESpace::HDiv, dim>(k);
        }
        break;
      }
      case ESpace::HDiv:
      {
        CheckRankKerDim<ESpace::HDiv, ESpace::L2, dim>(k);
        break;
      }
      case ESpace::HDivConst:
      {
        CheckRankKerDim<ESpace::HDivConst, ESpace::L2, dim>(k);
        break;
      }
      case ESpace::L2:
        break;
      }
    }
  }
}

TEMPLATE_TEST_CASE("Inclusion", "[derham_tests]",
                   (typename std::integral_constant<int, 2>),
                   (typename std::integral_constant<int, 3>))
{
  // SVD requires LAPACK
#ifndef PZ_USING_LAPACK
  return;
#endif
  constexpr int dim = TestType::value;

  ESpace leftSpace = GENERATE(ESpace::H1,
                              ESpace::HCurl,
                              ESpace::HDiv,
                              ESpace::HDivConst);
  SECTION(names.at(leftSpace))
  {
    int k = GENERATE(1, 2, 3);
    SECTION(std::to_string(k))
    {
      switch (leftSpace)
      {
      case ESpace::H1:
      {
        CheckInclusion<ESpace::H1, ESpace::HCurl, dim>(k);
        break;
      }
      case ESpace::HCurl:
      {
        CheckInclusion<ESpace::HCurl, ESpace::HDiv, dim>(k);
        break;
      }
      case ESpace::HDiv:
      {
        CheckInclusion<ESpace::HDiv, ESpace::L2, dim>(k);
        break;
      }
      case ESpace::HDivConst:
      {
        CheckInclusion<ESpace::HDivConst, ESpace::L2, dim>(k);
        break;
      }
      default:
        break;
      }
    }
  }
}

/*******************************************************************************
 ***********************************AUX FUNCS***********************************
 *******************************************************************************/

/** @brief Creates appropriate computational meshes to study De Rham compatibility of FEM approximation spaces.
    The mesh will consist of elements of type elType (inside the function,
    the option singleElement decides whether to create a mesh with only one
    element or not).
    @param [in] kRight polynomial order of the rightmost approximation space
    @param [in] elType element type to be created
    @param [out] kLeft polynomial order of the leftmost approximation space
    @param [out] gmesh created geometric mesh
    @param [out] cmeshL computational mesh of leftmost approximaton space
    @param [out] cmeshR computational mesh of rightmost approximaton space*/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CreateMeshes(int &kRight, MMeshType elType,
                  int &kLeft, TPZAutoPointer<TPZGeoMesh> &gmesh,
                  TPZAutoPointer<TPZCompMesh> &cmeshL,
                  TPZAutoPointer<TPZCompMesh> &cmeshR);

/** @brief Creates a simple geometric mesh for testing purposes
    @param [in] dim dimension of the domain
    @param [in] elType element type
    @param [in] matId material identifier of all elements
*/
TPZAutoPointer<TPZGeoMesh> CreateGMesh(const int dim, const MMeshType elType,
                                       const int matId);
/** @brief Creates a simple computational mesh for testing purposes
    @param [in] gMesh geometric mesh
    @param [in] mat material to be inserted in the mesh
    @param [in] matId material identifier
    @param [in] k default polynomial order
    @param [in] space approximation space to be created
*/
TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gMesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space);

/*
  @brief Calculates the rank of an S matrix obtained by SVD
  @param [in] S matrix
  @param[in] tol arithmetic tolerance
*/
int CalcRank(const TPZFMatrix<STATE> &S, const STATE tol);

/*******************************************************************************
***********************************TEST FUNCS***********************************
*******************************************************************************/

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckRankKerDim(int kRight)
{

  /******************************************
  CHECK FUNCTION DECLARATION FOR MORE DETAILS
  *******************************************/
  // select element type
  auto elType = GENERATE(MMeshType::ETriangular,
                         MMeshType::EQuadrilateral,
                         MMeshType::ETetrahedral,
                         MMeshType::EHexahedral,
                         MMeshType::EPrismatic);

  SECTION(MMeshType_Name(elType))
  {
    if (elType == MMeshType::EPrismatic && leftSpace == ESpace::HDivConst)
    {
      return; // TODOFIX
    }
    const auto elDim = MMeshType_Dimension(elType);
    if (elType == MMeshType::EPrismatic && leftSpace == ESpace::HCurl)
    {
      return; // TODOFIX
    }

    // if dimension does not correspond, skip
    if (elDim != dim)
      return;

    int kLeft;
    TPZAutoPointer<TPZGeoMesh> gmesh{nullptr};
    TPZAutoPointer<TPZCompMesh> cmeshL{nullptr}, cmeshR{nullptr};
    // creates computational meshes
    CreateMeshes<leftSpace, rightSpace, dim>(kRight, elType,
                                             kLeft, gmesh, cmeshL, cmeshR);

    // for debuggin
    const auto nameLeft = names.at(leftSpace);
    const auto nameRight = names.at(rightSpace);
    const auto elName = MMeshType_Name(elType);
    CAPTURE(nameLeft, nameRight, elName);

    // stores the matrices
    TPZAutoPointer<TPZMatrix<STATE>> matPtrL = nullptr;
    TPZAutoPointer<TPZMatrix<STATE>> matPtrR = nullptr;

    {
      constexpr int nThreads{4};
      TPZLinearAnalysis anL(cmeshL, RenumType::ENone);
      TPZLinearAnalysis anR(cmeshR, RenumType::ENone);

      TPZFStructMatrix<STATE> strmtrxL(cmeshL);
      strmtrxL.SetNumThreads(nThreads);
      anL.SetStructuralMatrix(strmtrxL);

      TPZFStructMatrix<STATE> strmtrxR(cmeshR);
      strmtrxR.SetNumThreads(nThreads);
      anR.SetStructuralMatrix(strmtrxR);

      TPZStepSolver<STATE> step;
      step.SetDirect(ELU);

      anL.SetSolver(step);
      anR.SetSolver(step);

      anL.Assemble();
      anR.Assemble();

      matPtrL = anL.MatrixSolver<STATE>().Matrix();
      matPtrR = anR.MatrixSolver<STATE>().Matrix();
    }

    TPZFMatrix<STATE> matL(*matPtrL);
    TPZFMatrix<STATE> matR(*matPtrR);

    const int dimL = matL.Rows();
    const int dimR = matR.Rows();

    TPZFMatrix<STATE> SL, SR;

    {
      TPZFMatrix<STATE> Udummy, VTdummy;
      matL.SVD(Udummy, SL, VTdummy, 'N', 'N');
      matR.SVD(Udummy, SR, VTdummy, 'N', 'N');
    }
    static constexpr auto tol = std::numeric_limits<STATE>::epsilon() * 100000;

    TPZManVector<int, 27> connectorderL(cmeshL->NConnects(), 0), connectorderR(cmeshR->NConnects(), 0);
    for (int i = 0; i < cmeshL->NConnects(); i++)
    {
      connectorderL[i] = cmeshL->ConnectVec()[i].Order();
    }
    for (int i = 0; i < cmeshR->NConnects(); i++)
    {
      connectorderR[i] = cmeshR->ConnectVec()[i].Order();
    }
    TPZManVector<int, 27> nshapeL(cmeshL->NConnects(), 0), nshapeR(cmeshR->NConnects(), 0);
    for (int i = 0; i < cmeshL->NConnects(); i++)
    {
      nshapeL[i] = cmeshL->ConnectVec()[i].NShape();
    }
    for (int i = 0; i < cmeshR->NConnects(); i++)
    {
      nshapeR[i] = cmeshR->ConnectVec()[i].NShape();
    }

    // rank of op(phi_i,phi_j) of left space
    const int rankL = CalcRank(SL, tol);
    const int rankR = CalcRank(SR, tol);
    // dimension of the nullspace of op(phi_i,phi_j) of right space
    const int kerR = dimR - rankR;

    CAPTURE(kLeft, kRight, dimL, rankL, dimR, kerR, rankR);
    CAPTURE(connectorderL, connectorderR);
    CAPTURE(nshapeL, nshapeR);
    // CAPTURE(SL, SR);
    if constexpr (rightSpace == ESpace::L2)
    { // for L2 the operator is -> 0
      REQUIRE(rankL == rankR);
    }
    else
    {
      REQUIRE(rankL == kerR);
    }
  }
}

template <ESpace leftSpace, ESpace rightSpace, int dim>
void CheckInclusion(int kRight)
{
  /******************************************
  CHECK FUNCTION DECLARATION FOR MORE DETAILS
  *******************************************/

  constexpr int matId = 1;
  // select element type
  auto elType = GENERATE(MMeshType::ETriangular, MMeshType::EQuadrilateral,
                         MMeshType::ETetrahedral, MMeshType::EHexahedral,
                         MMeshType::EPrismatic);

  SECTION(MMeshType_Name(elType))
  {
    const auto elDim = MMeshType_Dimension(elType);

    if (elType == MMeshType::EPrismatic && leftSpace == ESpace::HDivConst)
    {
      return; // TODOFIX
    }

    if (elType == MMeshType::EPrismatic && leftSpace == ESpace::HCurl)
    {
      return; // TODOFIX
    }
    // skips if dimension is not compatible
    if (elDim != dim)
      return;

    int kLeft;
    TPZAutoPointer<TPZGeoMesh> gmesh{nullptr};
    TPZAutoPointer<TPZCompMesh> cmeshL{nullptr}, cmeshR{nullptr};
    // for debugging
    const auto nameLeft = names.at(leftSpace);
    const auto nameRight = names.at(rightSpace);
    const auto elName = MMeshType_Name(elType);
    CAPTURE(nameLeft, nameRight, elName);

    // creates computational meshes
    CreateMeshes<leftSpace, rightSpace, dim>(kRight, elType, kLeft,
                                             gmesh, cmeshL, cmeshR);

    // creates multiphysics mesh
    TPZMultiphysicsCompMesh *cmeshMF = new TPZMultiphysicsCompMesh(gmesh);

    /*
      Now, the appropriate material will be selected given their
      combination of left and right approximation spaces.
      This material will be responsible for creating the matrix M
     */
    auto *matMF = []() -> TPZMaterial *
    {
      if constexpr (rightSpace == ESpace::H1)
      {
        return nullptr;
      }
      else if constexpr (rightSpace == ESpace::HCurl)
      {
        return new TPZMatDeRhamH1HCurl(matId, dim);
      }
      else if constexpr (rightSpace == ESpace::HDiv)
      {
        if constexpr (dim == 3)
        {
          return new TPZMatDeRhamHCurlHDiv(matId, dim);
        }
        return nullptr;
      }
      else if constexpr (rightSpace == ESpace::L2)
      {
        return new TPZMatDeRhamHDivL2(matId, dim);
      }
    }();

    if (!matMF)
      return;

    cmeshMF->InsertMaterialObject(matMF);
    cmeshMF->SetDimModel(dim);
    cmeshMF->SetAllCreateFunctionsMultiphysicElem();

    TPZManVector<MSpaceConfig, 2> active_approx_spaces(2, ETestTrial);
    TPZManVector<TPZCompMesh *, 2> meshVec(2);
    meshVec[0] = cmeshL.operator->();
    meshVec[1] = cmeshR.operator->();

    cmeshMF->BuildMultiphysicsSpace(active_approx_spaces, meshVec);
    cmeshMF->CleanUpUnconnectedNodes();
    cmeshMF->LoadReferences();

    TPZLinearAnalysis analysis(cmeshMF, RenumType::ENone);
    TPZFStructMatrix<STATE> strmtrx(cmeshMF);
    analysis.SetStructuralMatrix(strmtrx);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELU);
    analysis.SetSolver(step);

    analysis.Assemble();
    // it is a full matrix
    auto mfMatrix =
        dynamic_cast<TPZFMatrix<STATE> &>(
            analysis.MatrixSolver<STATE>().Matrix().operator*());

    const int rLeftEqs = cmeshR->NEquations();
    const int lLeftEqs = cmeshL->NEquations();
    TPZFMatrix<STATE> rightMatrix(rLeftEqs, rLeftEqs);

    // rightMatrix is the matrix phi_i phi_j for the approx. space on the right
    mfMatrix.GetSub(lLeftEqs, lLeftEqs, rLeftEqs, rLeftEqs, rightMatrix);

    TPZFMatrix<STATE> Sfull, Sright;

    {
      TPZFMatrix<STATE> Udummy, VTdummy;
      mfMatrix.SVD(Udummy, Sfull, VTdummy, 'N', 'N');
      rightMatrix.SVD(Udummy, Sright, VTdummy, 'N', 'N');
    }

    /**
       We want to compare if rank(rightMatrix) == rank(mfMatrix).
       This will ensure us that op(phi_i).op(phi_j) of the left space
       are contained in the functions of the approximation space on the right.
     */
    static constexpr auto tol = std::numeric_limits<STATE>::epsilon() * 10000;
    auto fullDim = Sfull.Rows();
    auto fullRank = CalcRank(Sfull, tol);
    auto rightDim = Sright.Rows();
    auto rightRank = CalcRank(Sright, tol);

    CAPTURE(kLeft, kRight, fullDim, fullRank, rightDim, rightRank);
    REQUIRE(fullRank == rightRank);
  }
}

/*******************************************************************************
 ***********************************AUX FUNCS***********************************
 *******************************************************************************/
template <ESpace leftSpace, ESpace rightSpace, int dim>
void CreateMeshes(int &kRight, MMeshType elType,
                  int &kLeft, TPZAutoPointer<TPZGeoMesh> &gmesh,
                  TPZAutoPointer<TPZCompMesh> &cmeshL,
                  TPZAutoPointer<TPZCompMesh> &cmeshR)
{

  constexpr int matId = 1;

  const auto elDim = MMeshType_Dimension(elType);

  // for CheckInclusion test, the following materials will be ignored
  auto *matLeft = [&kLeft, &kRight, elType]()
      -> TPZMaterial *
  {
    if constexpr (leftSpace == ESpace::H1)
    {
      kLeft = kRight + 1;
      return new TPZMatDeRhamH1(matId, dim);
    }
    else if constexpr (leftSpace == ESpace::HCurl)
    {
      if (elType == MMeshType::EQuadrilateral || elType == MMeshType::EHexahedral)
      {
        kLeft = kRight;
      }
      else
      {
        kLeft = kRight + 1;
      }
      return new TPZMatDeRhamHCurl(matId, dim);
    }
    else if constexpr (leftSpace == ESpace::HDiv)
    {
      kLeft = kRight;
      return new TPZMatDeRhamHDiv(matId, dim);
    }
    else if constexpr (leftSpace == ESpace::HDivConst)
    {
      kLeft = kRight;
      kRight = 0;
      return new TPZMatDeRhamHDiv(matId, dim);
    }
    else if constexpr (leftSpace == ESpace::L2)
    {
      return nullptr;
    }
  }();

  auto *matRight = []()
      -> TPZMaterial *
  {
    if constexpr (rightSpace == ESpace::H1)
    {
      return nullptr;
    }
    else if constexpr (rightSpace == ESpace::HCurl)
    {
      return new TPZMatDeRhamHCurl(matId, dim);
    }
    else if constexpr (rightSpace == ESpace::HDiv)
    {
      return new TPZMatDeRhamHDiv(matId, dim);
    }
    else if constexpr (rightSpace == ESpace::HDivConst)
    {
      return new TPZMatDeRhamHDiv(matId, dim);
    }
    else if constexpr (rightSpace == ESpace::L2)
    {
      return new TPZMatDeRhamL2(matId, dim);
    }
  }();

  REQUIRE((matLeft != nullptr && matRight != nullptr));
  // generates mesh with a single element
  constexpr bool singleElement{true};

  if constexpr (singleElement)
  {
    constexpr bool createBoundEls{false};
    gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(elType, matId, createBoundEls);
  }
  else
  {
    gmesh = CreateGMesh(dim, elType, matId);
  }

  auto EnrichMesh = [](TPZAutoPointer<TPZCompMesh> cmesh, const int k)
  {
    for (auto cel : cmesh->ElementVec())
    {
      if (cel->Dimension() != dim)
        continue;
      auto intel =
          dynamic_cast<TPZInterpolatedElement *>(cel);
      const auto nsides = intel->Reference()->NSides();
      intel->SetSideOrder(nsides - 1, k + 1);
      intel->AdjustIntegrationRule();
    }
    cmesh->ExpandSolution();
  };

  cmeshL = CreateCMesh(gmesh, matLeft, matId, kLeft, leftSpace);

  // Enriching the internal functions for Hdiv or HdivConst
  if constexpr (leftSpace == ESpace::HDiv || leftSpace == ESpace::HDivConst)
  {
    int enrichLvl = GENERATE(-1, 0, 1);
    // SECTION("Enrichment Level: " + std::to_string(enrichLvl+1));
    EnrichMesh(cmeshL, kLeft + enrichLvl);
    if constexpr (leftSpace == ESpace::HDiv)
    { // In this case, the L2 space order has to be equal to the Hdiv enriched space
      kRight = kLeft + enrichLvl + 1;
      cmeshR = CreateCMesh(gmesh, matRight, matId, kRight, rightSpace);
    }
    else
    {
      cmeshR = CreateCMesh(gmesh, matRight, matId, kRight, rightSpace);
    }
  }
  else
  {
    cmeshR = CreateCMesh(gmesh, matRight, matId, kRight, rightSpace);
  }

  if (elType == MMeshType::ETetrahedral)
  {
    if constexpr (leftSpace == ESpace::HCurl)
    {
      EnrichMesh(cmeshL, kLeft);
    }
    else if constexpr (rightSpace == ESpace::HCurl)
    {
      EnrichMesh(cmeshR, kRight);
    }

    if constexpr (leftSpace == ESpace::H1)
    {
      EnrichMesh(cmeshL, kLeft);
    }
  }
}

TPZAutoPointer<TPZGeoMesh> CreateGMesh(const int dim, const MMeshType elType,
                                       const int matId)
{
  constexpr int ndiv{2};
  TPZManVector<int, 3> nDivVec(dim, ndiv);
  TPZManVector<REAL, 3> minX({0, 0, 0});
  TPZManVector<REAL, 3> maxX({1, 1, 1});
  if (dim == 2)
    maxX[2] = 0.0;
  constexpr bool createBoundEls{false};
  TPZAutoPointer<TPZGeoMesh> gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(
      dim, minX, maxX, {matId}, nDivVec, elType, createBoundEls);
  return gmesh;
}

TPZAutoPointer<TPZCompMesh> CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh,
                                        TPZMaterial *mat, const int matId,
                                        const int k, ESpace space)
{

  TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
  const int nel = cmesh->NElements();
  cmesh->SetDefaultOrder(k);
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->InsertMaterialObject(mat);
  switch (space)
  {
  case ESpace::H1:
    cmesh->SetAllCreateFunctionsContinuous();
    break;
  case ESpace::HCurl:
    cmesh->SetAllCreateFunctionsHCurl();
    break;
  case ESpace::HDiv:
    cmesh->SetAllCreateFunctionsHDiv();
    break;
  case ESpace::HDivConst:
    cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    cmesh->SetAllCreateFunctionsHDiv();
    break;
  case ESpace::L2:
    if (k == 0)
      cmesh->SetAllCreateFunctionsDiscontinuous();
    else
      cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    break;
  }
  cmesh->AutoBuild();
  cmesh->CleanUpUnconnectedNodes();
  return cmesh;
}

int CalcRank(const TPZFMatrix<STATE> &S, const STATE tol)
{
  int rank = 0;
  const int dimMat = S.Rows();
  for (int i = 0; i < dimMat; i++)
  {
    rank += S.GetVal(i, 0) > tol ? 1 : 0;
  }
  return rank;
};
