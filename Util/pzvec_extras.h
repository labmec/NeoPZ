/**
 * @file
 * @brief Extra utilities for TPZVec. Implementations of the saxpy, sscal, sdot, intercept, max and min functions.
 */

#ifndef PZVEC_EXTRAS_H
#define PZVEC_EXTRAS_H

#include <algorithm>
#include "pzstack.h"


/**
 * \addtogroup util
 * @{
 */


//! Append to TPZVec
template<class T, class R>//utility function
void AppendToVec(TPZVec<T> &vec, R t){
  static_assert(std::is_convertible_v<R, T>);
  const auto sz = vec.size();
  vec.Resize(sz+1);
  vec[sz] = (T)t;
}

//! Removes duplicates AND sort vector contents
template<class T>
void RemoveDuplicates(TPZVec<T> &v)
{
    std::sort( v.begin(), v.end() );
    auto newend = std::unique( v.begin(), v.end() );
    const auto p = newend - v.begin();
    v.Resize(p);
}

/**
 * @since Jan 9, 2002
 * @author Cantao!
 * @brief Performs a saxpy operation: x <- x + s * y.
 */
template< class T1, class T2, class Scalar >
void saxpy(TPZVec< T1 >& x, const TPZVec< T2 >& y, Scalar s) {
    int size = x.NElements();

#ifdef PZDEBUG
    if (size != y.NElements()) {
        PZError << "SAXPY error!" << std::endl
                << "Vectors with different sizes #x = " << size
                << ", #y = " << y.NElements() << std::endl;

        PZError.flush();
    }
#endif

    for (int ii = 0; ii < size; ii++) {
        x[ ii ] += s * y[ ii ];
    }
}

/**
 * @brief Performs a sscal operation: x <- x * s.
 * @since Mar 20, 2003
 * @author Tiago Forti
 * @author Erick Santos
 */
template< class T1, class Scalar >
void sscal(TPZVec< T1 > & x, const Scalar s) {

    int size = x.NElements();

    for (int ii = 0; ii < size; ii++) {
        x[ii] *= s;
    }
}

/**
 * @brief substracts two vectors
 */
template <class T>
TPZVec<T> operator-(const TPZVec<T> &a, const TPZVec<T> &b) {
    if (a.size() != b.size()) {
        DebugStop();
    }
    TPZVec<T> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

/**
 * @brief substracts two vectors
 */
template <class T>
TPZVec<T> &operator-=(TPZVec<T> &a, const TPZVec<T> &b) {
    if (a.size() != b.size()) {
        DebugStop();
    }
    for (int i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}

/**
 * @brief divides vector by scalar
 */
template <class T>
TPZVec<T> &operator/=(TPZVec<T> &a, const T& b) {
    
    for (int i = 0; i < a.size(); i++) {
        a[i] /= b;
    }
    return a;
}


/**
 * @brief Performs a sdot operation: dot <- Transpose[x] * y
 *
 * @since Mar 20, 2003
 * @author Tiago Forti
 * @author Erick Santos
 */
template< class T1 >
double sdot(TPZVec< T1 > & x, TPZVec< T1 > & y) {

    int size = x.NElements();
    double sum = 0.0;

#ifndef PZNODEBUG
    if (size != y.NElements()) {
        PZError << "SDOT error!" << std::endl
                << "Vectors with different sizes #x = " << size
                << ", #y = " << y.NElements() << std::endl;

        PZError.flush();
    }
#endif

    for (int ii = 0; ii < size; ii++) {
        sum += x[ii] * y[ii];
    }

    return sum;
}

template < class T1 >
REAL dist(TPZVec<T1> &vec1, TPZVec<T1> &vec2) {
#ifdef PZDEBUG
    if (vec1.size() != vec2.size()) {
        DebugStop();
    }
#endif
    REAL dist = 0.;
    for (int i = 0; i < vec1.size(); i++) {
        dist += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    dist = sqrt(dist);
    return dist;
}

//--| SORTING |-----------------------------------------------------------------

/** @brief Sorting the elements into v */
template< class T >
TPZVec< T >& Sort(TPZVec< T >& v) {
    std::sort(static_cast<T*> (v.begin()), v.begin() + v.NElements());

    return v;
}

/** @brief Finds if exists the element e into vector v */
template< class T >
int64_t Find(const TPZVec< T >& v, const T& e) {
    auto* found = std::find(v.begin(), v.end(), e);

    const int64_t dist = found - v.begin();

    return (dist == v.NElements() ? -1 : dist);
}

/** @brief Returns the minimum element into v */
template< class T >
T Min(TPZVec< T >& v) {
    int nel = v.NElements();

    T m = v[ 0 ];

    for (int ii = 1; ii < nel; ii++) {
        m = min(m, v[ ii ]);
    }

    return m;
}

/** @brief Returns the maximum element into v */
template< class T >
T Max(TPZVec< T >& v) {
    int nel = v.NElements();

    T m = v[ 0 ];

    for (int ii = 1; ii < nel; ii++) {
        m = max(m, v[ ii ]);
    }

    return m;
}

/** @brief Gets commom elements into the one and two vectors */
template< class T, int N >
void Intersect(const TPZVec< T > &one, const TPZVec< T > &two, TPZStack< T, N > &result) {
    int firstc, secondc, nfirst, nsecond;
    nfirst = one.NElements();
    nsecond = two.NElements();
    firstc = 0;
    secondc = 0;
    while (firstc < nfirst && secondc < nsecond) {
        while (firstc < nfirst && one[firstc] < two[secondc]) {
            firstc++;
        }
        if (firstc == nfirst) break;
        while (secondc < nsecond && two[secondc] < one[firstc]) {
            secondc++;
        }
        if (firstc < nfirst && secondc < nsecond && one[firstc] == two[secondc]) {
            result.Push(one[firstc]);
            firstc++;
            secondc++;
        }
    }

}

/** @brief Gets commom elements into the one, two and three vectors */
template< class T, int N >
void Intersect(const TPZVec< T > &one, const TPZVec< T > &two, const TPZVec< T > &three, TPZStack< T, N > &result) {
    int firstc, secondc, thirdc, nfirst, nsecond, nthird;
    nfirst = one.NElements();
    nsecond = two.NElements();
    nthird = three.NElements();
    firstc = 0;
    secondc = 0;
    thirdc = 0;
    while (firstc < nfirst && secondc < nsecond && thirdc < nthird) {
        while (firstc < nfirst && (one[firstc] < two[secondc] || one[firstc] < three[thirdc])) {
            firstc++;
        }
        if (firstc == nfirst)break;
        while (secondc < nsecond && (two[secondc] < one[firstc] || two[secondc] < three[thirdc])) {
            secondc++;
        }
        if (secondc == nsecond) break;
        while (thirdc < nthird && (three[thirdc] < one[firstc] || three[thirdc] < two[secondc])) {
            thirdc++;
        }
        if (firstc < nfirst && secondc < nsecond && thirdc < nthird && one[firstc] == two[secondc] && one[firstc] == three[thirdc]) {
            result.Push(one[firstc]);
            firstc++;
            secondc++;
            thirdc++;
        }
    }

}

/** @brief Gets commom elements into the one and two vectors */
template< class T>
RType(T) Norm(const TPZVec< T > &one) {
    RType(T) res = 0.;
    int size = one.NElements();
    if constexpr (is_complex<T>::value){
        T cres{0};
        for (int i = 0; i < size; i++) {
            cres += one[i] * std::conj(one[i]);
        }
        res = cres.real();
    }else{
        for (int i = 0; i < size; i++) {
            res += one[i]*one[i];
        }
    }
    
    return sqrt(res);
}

template<class T>
void Cross(const TPZVec <T> &x1, const TPZVec<T> &x2, TPZVec<T> &result) {
#ifdef PZDEBUG
    if (x1.size() != 3) {
        DebugStop();
    }
    if (x2.size() != 3) {
        DebugStop();
    }
    if (result.size() != 3) {
        DebugStop();
    }
#endif
    result[0] = x1[1] * x2[2] - x2[1] * x1[2];
    result[1] = x1[2] * x2[0] - x2[2] * x1[0];
    result[2] = x1[0] * x2[1] - x2[0] * x1[1];
}

template<class T>
T Dot(const TPZVec <T> &x1, const TPZVec<T> &x2) {
    T result = T(0.);
    int size = x1.NElements();

#ifdef PZDEBUG
    if (size != x2.NElements()) {
        PZError << "Dot error!" << std::endl
                << "Vectors with different sizes #x1 = " << size
                << ", #x2 = " << x2.NElements() << std::endl;

        PZError.flush();
    }
#endif

    for (int i = 0; i < size; ++i) {
        result = result + x1[i] * x2[i];
    }
    return result;
}

/** @} */

#endif //PZVEC_EXTRAS_H

