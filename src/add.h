/* Scythe_add.h
 *
 * This file provides my  additional Code to the
 * Scythe Statistical Library.
 *
 *
 * Copyright (C) 2005, S M MWalili
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  A copy of this license is included
 * with this library (LICENSE.GPL).
 *
 * This library utilizes code from a number of other open source
 * projects.  Specific copyright information is provided with the
 * applicable code.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * This  program uses the
 * Scythe Statistical Library
 * Copyright (C) 2003, Andrew D. Martin, Kevin M. Quinn, and Daniel
 * Pemstein.  All Rights Reserved.
 *
 * This code is written by:
 *
 * S M Mwalili
 * Biostatistical Centre
 * UZ St Rafael
 * Kapucynenvoer 35
 * B-3000 Leuven
 * Belgium
 *
 * Tel: +32-16-336887
 * Fax: +32-16-336900
 * Samuel.Mwalili@med.kuleuven.ac.be
 */


#ifndef SCYTHE_ADD_H
#define SCYTHE_ADD_H

#include "matrix.h"

namespace SCYTHE {

   // delete the r-th row and c-th column of  a matrix
   template <class T>
   Matrix<T> delrc(const Matrix<T> &A, const int &a, const int &b);

   // delete  c-th column of  a matrix
    template <class T>
    Matrix<T>  delc (const Matrix<T> &, const int &);

   // return a subset of matrix A[a:c,b:d]
    template <class T>
    Matrix<T> subrc (const Matrix<T> &, const int &, const int &,
                             const int &,const int &);

   // Diagonal matrix from a vector
   template <class T>
    Matrix<T> diagv(Matrix<T>&);

    // returns true if the matrix is positive definite, false otherwise
    template <class T>
    bool ispd (const Matrix<T> &);

    // rotation required by JACOBI function
    template <class T>
    void ROT(Matrix<T>& A, const double&, const double& , const int&, const int&,
                               const int& k,   const int&);

    // computing eigen value/ eigen vector for symmetric matrix
     template <class T>
     void JACOBI(const Matrix<T>& B, Matrix<T>& d, Matrix<T>& V);

    // evaluate the conditional variance of MVN
    template <class T>
    double cnd_var(const Matrix<T>&, const int &);

    // evaluate the conditional mean of MVN
    template <class T>
    double cnd_mean(const Matrix<T>&, const Matrix<T>&,
                    const Matrix<T>&, const int &);

    // return A %x% ones (K) = a K times stack of A
    template <class T>
    Matrix<T> xvec(const Matrix<T>& A, const int & K);

    // extract the unique elements of the correlation matrix
    template <class T>
    Matrix<T> UNIK(const Matrix<T> &A);

    // evaluate the conditional moments of MVN in efficient way
    template <class T>
	void cndnorm(const Matrix<T> &Z, const Matrix<T> &Xbeta,
             const Matrix<T> &invS, const int& p,
             const int& j, double& muj, double& varj);

} // end namespace SCYTHE
#if defined (__GNUG__) || defined (__MWERKS__) || defined (_MSC_VER) || \
    defined (EXPLICIT_TEMPLATE_INSTANTIATION)
  // Necessary for template instantiation with some compilers.
# include "add.cc"
#endif  /* EXPLICIT_TEMPLATE_INSTANTIATION */

#endif /* SCYTHE_ADD_H */

