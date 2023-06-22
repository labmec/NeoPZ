#include "pzreal.h"
#include "TPZSimpleTimer.h"
#include "tpzautopointer.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMatrixSolver.h"

// #define PZGMRESDEBUG

template<class TVar>
void GeneratePlaneRotation(const TVar dx, const TVar dy, TVar &cs, TVar &sn, RTVar tol);

template<class TVar>
void ApplyPlaneRotation(TVar &dx, TVar &dy, const TVar cs, const TVar sn);

template<class TVar>
void Update(TPZFMatrix<TVar>&x, const int k, const TPZFMatrix<TVar>&H,
            const TPZVec<TVar> &s, const TPZVec<TPZFMatrix<TVar>>&v);

template<class TVar>
int GMRES(const TPZMatrix<TVar> &A, TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &b,
          TPZMatrixSolver<TVar> &M, TPZFMatrix<TVar> &H, const int krylovdim,
          int64_t &max_iter, RTVar &tol,
          TPZFMatrix<TVar> *residual, const int fromcurrent,
          bool print_res=true){

  constexpr RTVar ztol = 10*std::numeric_limits<RTVar>::epsilon();
  //for printing in full precision
  std::ios cout_state(nullptr);
  cout_state.copyfmt(std::cout);
    
  std::cout << std::setprecision(std::numeric_limits<STATE>::max_digits10);

  
  TPZSimpleTimer gmres("GMRES");
  //allocating structures
  TPZVec<TPZFMatrix<TVar>> v(krylovdim+1);
  TPZVec<TVar> cs(krylovdim+1), sn(krylovdim+1), s(krylovdim+1);

  //compute rhs norm
  TPZFMatrix<TVar> r;
  M.Solve(b,r);
  const RTVar normb = [&r]{
    RTVar norm_b = Norm(r);
    if(norm_b < ztol) return (RTVar)1;
    return norm_b;
  }();
	
  if(fromcurrent){//we actually need to compute the residual
    //r = b-A*x without dynamic allocation
    A.MultAdd(x,b,r,-1,1);
    M.Solve(r,r);
  }else{//residual is equal to rhs for initial sol == 0   
    x.Zero();
  }
	
  RTVar beta = Norm(r);
  std::cout<<"initial precond correction "<<beta<<std::endl;
  RTVar resid = beta;
  if (resid <= tol) {
    tol = resid;
    max_iter = 0;
    std::cout.copyfmt(cout_state);
    return 0;
  }
#ifdef PZGMRESDEBUG
  TPZFMatrix<TVar> xcp = x;
  TPZFMatrix<TVar> rcp = r;
#endif
  //we avoid allocating w at every step
  TPZFMatrix<TVar> w(r.Rows(),1,0), w1;
  int iter = 1;
  while (iter <= max_iter) {
    r.MultiplyByScalar((1.0/beta),v[0]);
    s = 0.0;
    s[0] = beta;

    int i = 0;
    for (; i < krylovdim && iter <= max_iter; i++, iter++) {
      A.Multiply(v[i],w1);
      M.Solve(w1,w);
      for (int k = 0; k <= i; k++) {
        const auto dotwv = Dot(w,v[k]);
        H(k, i) = dotwv;
        w -= dotwv * v[k];
      }
      const auto normw = Norm(w);
      if(normw < ztol){
        i++;
        break;
      }
      H(i+1, i) = normw;
      w *= ((TVar)1.0 / H(i+1, i));
      v[i+1] = w;

      for (int k = 0; k < i; k++){
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      }
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i], ztol);
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

      resid = std::abs(s[i+1])/normb;
      if (resid < tol) {
        i++;
        break;
      }
#ifdef PZGMRESDEBUG
      xcp = x;
      Update(xcp,i,H,s,v);
      const auto du = Norm(x-xcp);
      A.MultAdd(xcp,b,rcp,-1,1);
      M.Solve(rcp,rcp);
      auto nres = Norm(rcp);
      std::cout << iter
                <<"\n\t|s(i+1)|   :" << std::abs(s[i+1])
                <<"\n\t|res|      :" << nres
                <<"\n\t|res|/normb:" << nres/normb
                <<"\n\t|du|:" << du
                << std::endl;
#else      
      if(print_res){std::cout << iter << "\t" << resid << std::endl;}
#endif
    }
    Update(x, i - 1, H, s, v);
    //r = b-A*x without dynamic allocation
    A.MultAdd(x,b,r,-1,1);
    M.Solve(r,r);
    beta = Norm(r);
    resid = beta/normb;
    if(print_res){std::cout<<"beta "<<beta<<" resid "<<resid<<" normb "<<normb<<std::endl;}
    if (resid < tol) {
      break;
    }
  }
  max_iter = iter;
  tol = resid;
  std::cout.copyfmt(cout_state);
  return 1;
}

/** @ingroup util */
/** @brief Compute the values cs and sn parameters to rotation
    @note These rotations are based on lapack Xlartg routines. (X=s,d,c,z)
    They will result in a plane rotation with real cosine and complex sine*/
template<class TVar>
void GeneratePlaneRotation(TVar dx, TVar dy, TVar &cs, TVar &sn, RTVar tol)
{
  
  if(abs(dy) < tol){
    cs = 1.;
    sn = 0;
  }else if (abs(dx)<tol){
    cs = 0;
    sn = dy/std::abs(dy);
  }else{
    //note: std::norm is actually the squared magnitude
    const TVar sq = sqrt(std::norm(dx) + std::norm(dy));
    cs = abs(dx)/sq;
    if constexpr (is_complex<TVar>::value){
      sn = (dx/abs(dx))*std::conj(dy)/sq;
    }else{
      sn = (dx/abs(dx))*dy/sq;
    }
  }
  
}

/** @ingroup util */
/** @brief Makes rotation of the plane based on the cs and sn parameters
    @note These rotations are based on lapack Xlartg routines. (X=s,d,c,z)
    They will result in a plane rotation with real cosine and complex sine*/
template<class TVar>
void ApplyPlaneRotation(TVar &dx, TVar &dy, const TVar cs, const TVar sn)
{
  const TVar temp = cs * dx + sn * dy;
  if constexpr (is_complex<TVar>::value){
    dy = -std::conj(sn) * dx + std::conj(cs) * dy;
  }else{
    dy = -sn * dx + cs * dy;
  }
  dx = temp;
}

template<class TVar>
void Update(TPZFMatrix<TVar>&x, const int k, const TPZFMatrix<TVar>&H,
            const TPZVec<TVar> &s, const TPZVec<TPZFMatrix<TVar>>&v)
{
  TPZVec<TVar> y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= H(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= H(j,i) * y[i];
  }

  for (int j = 0; j <= k; j++)
    x += v[j] * y[j];
}