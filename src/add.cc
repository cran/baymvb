/* Scythe_Add.cc
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

#ifndef SCYTHE_ADD_CC
#define SCYTHE_ADD_CC

#include <cmath>
#include <algorithm>
#include <set>
#include "distributions.h"
#include "error.h"
#include "util.h"
#include "ide.h"
#include "stat.h"
#include "la.h"
#include "add.h"

namespace SCYTHE {

   // delete the r-th row and c-th column of  a matrix
   template <class T>
   Matrix<T> delrc(const Matrix<T> &A, const int &a, const int &b)
   {
     Matrix<T> temp(A.rows()-1, A.cols()-1, false);
     int cnt = -1;
 
     for (int i = 0; i < A.rows(); ++i) {
      for (int j = 0; j < A.cols(); ++j) {
         if(i!=(a-1) & j!=(b-1))
             temp[++cnt] = A(i,j);
       }
     }
     return temp;
  }
 
   // delete  c-th column of  a matrix
   template <class T>
   Matrix<T>
   delc (const Matrix<T> &A, const int &a) 
   {
     Matrix<T> temp((A.rows() -1),1,false);
     int cnt = -1;
     for (int i = 0; i < A.rows(); ++i)
        if(i!=(a-1)) 
            temp[++cnt] = A[i];
      return temp;
   }
   
   // return a subset of matrix A[a:c,b:d] 
   template <class T>
   Matrix<T>
   subrc (const Matrix<T> &A, const int &a, const int &b,
                              const int &c,const int &d) 
   {
     int k = -1;
     Matrix<T> temp((c - a + 1), (d - b + 1), false);
     for (int i = a-1; i < c; ++i)
       for (int j = b-1; j < d; ++j)
         temp[++k] = A[i * A.cols() + j];
     return temp;
  }
  
  // write a diagonal matrix from  a vector A 
  template <class T>
  Matrix<T> diagv (Matrix<T>& A)
  {
      int n = A.size();
      Matrix<T> temp =  Matrix<T> (n,n);
      for(int j=0; j<n; ++j)
       temp(j,j) = A[j];
       return temp;
  }
  
  // true if the matrix is positive definite
  template <class T>
   bool
   ispd (const Matrix<T> &A){
 
    if (! A.isSquare()) {
      return false;
    }
    
    Matrix<T> temp (A.rows(), A.cols(), false);
    register T h;
    
    for (int i = 0; i < A.rows(); ++i) {
      for (int j = i; j < A.cols(); ++j) {
        h = A(i,j);
        for (int k = 0; k < i; ++k) {
          h -= temp(i, k) * temp(j, k);
        }
        if (i == j) {
          if (h <= (T) 0) {
            return false;
          }
          temp(i,i) = ::sqrt(h);
        } else {
          temp(j,i) = (((T) 1) / temp(i,i)) * h;
          temp(i,j) = (T) 0;
        }
      }
    }
  
    return true;
  }
  
    // rotation function required by JACOBI
    template <class T>
    void ROT(Matrix<T>& A, const double& s, const double& tau, const int& i, const int& j,  
                               const int& k,   const int& l)
     {
             double g, h;
             g = A(i,j);
             h = A(k,l);
             A(i,j) = g-s*(h+g*tau);
             A(k,l) = h+s*(g-h*tau);
      }  
         
    // evaluate the eigein values and vectors of a symmetric matrix
     template <class T>
     void JACOBI(Matrix<T>& A, Matrix<T>& d, Matrix<T>& V)
    //Computes all eigenvalues and eigenvectors of a real symmetric matrix A(1..n,1..n). On
    //output, elements of a above the diagonal are destroyed. d(1..n) returns the eigenvalues of a.
    //V(1..n,1..n) is a matrix whose columns contain, on output, the normalized eigenvectors of
    //A. nrot returns the number of Jacobi rotations that were required.
    {  
        int nrot;
        int n = A.rows();
        int j,iq,ip,i;
        double tresh,theta,tau,t,sm,s,h,g,c;
        Matrix<double>b=Matrix<double>(1,n);
        Matrix<double>z=Matrix<double>(1,n);
        for (ip=0;ip<n;ip++){//Initialize to the identity matrix.
             for (iq=0;iq<n;iq++) V(ip,iq)=0.0;
            V(ip,ip)=1.0;
        }
        for (ip=0;ip<n;ip++) { //Initialize b and d to the diagonal of A
         b[ip]=d[ip]=A(ip,ip);
         z[ip]=0.0; //This vector will accumulate terms
        }
    
        nrot=0;
        for (i=1;i<=50;i++) {
            sm=0.0;
            for (ip=0;ip<n-1;ip++) { //Sum off-diagonal elements.
                for (iq=ip+1;iq<n;iq++)
                    sm += abs(A(ip,iq));
            }
        if (sm == 0.0) { //The normal return, which relies on quadratic convergence to machine underflow.
        return;
        }
        if (i < 4)
        tresh=0.2*sm/(n*n); //...on the first three sweeps.
            else
                tresh=0.0;// ...thereafter.
        for (ip=0;ip<n-1;ip++) {
            for (iq=ip+1;iq<n;iq++) {
            g=100.0*abs(A(ip,iq));
            //After four sweeps, skip the rotation if the off-diagonal element is small.
            if (i > 4 && (double)(abs(d[ip])+g) == (double)abs(d[ip])
                && (double)(abs(d[iq])+g) == (double)abs(d[iq]))
                A(ip,iq)=0.0;
            else if (abs(A(ip,iq)) > tresh) {
                h=d[iq]-d[ip];
                if ((double)(abs(h)+g) == (double)abs(h))
                    t=(A(ip,iq))/h;
                    else {
                        theta=0.5*h/(A(ip,iq));
                        t=1.0/(abs(theta)+::sqrt(1.0+theta*theta));
                    if (theta < 0.0) t = -t;
                    }
                    c=1.0/::sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*A(ip,iq);
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    A(ip,iq)=0.0;
                    for (j=0;j<ip;j++) { 
                        ROT(A,s,tau,j,ip,j,iq);
                    }
                    for (j=ip+1;j<iq;j++) { 
                        ROT(A,s,tau,ip,j,j,iq);
                    }
                    for (j=iq+1;j<n;j++) {
                        ROT(A,s,tau,ip,j,iq,j);
                    }
                    for (j=0;j<n;j++) {
                        ROT(V,s,tau,j,ip,j,iq);
                    }
                    ++nrot;
                }
            }
        }
        for (ip=0;ip<n;ip++) {
        b[ip] += z[ip];
        d[ip]=b[ip]; //Update d
        z[ip]=0.0; // and reinitialize z.
        }
    }
     throw scythe_invalid_arg(__FILE__, __PRETTY_FUNCTION__,
                 __LINE__, "Too many iterations in routine jacobi");
    
    }
    
    // evaluate the conditional variance of MVN
    template <class T>
    double cnd_var(const Matrix<T>& A, const int & k)
    {
        Matrix<double> S12(A.rows()-1,1);
        int cnt = -1;
        int k0 = k-1;
        double S11 = A(k0,k0);
        for(int j=0; j < A.rows(); ++j)
        if(k0!=j) 
            S12[++cnt] = A(k0,j);
        Matrix<double> S22 = delrc(A,k,k);
        double var = S11-(t(S12)*invpd(S22)*S12)[0];
       return var;
    } 
    
    // evaluate the conditional mean of MVN
    template <class T>
    double cnd_mean(const Matrix<T>& A, const Matrix<T>& Z,
                    const Matrix<T>& mu, const int & k)
    {
        Matrix<double> S12(1,A.rows()-1);
        Matrix<double> zstar(A.rows()-1,1);
        Matrix<double> mustar(A.rows()-1,1);

        int cnt = -1;
        int cnz = -1;
        int k0 = k-1;
        double muj = mu[k0];
         
        //double S11 = A(k0,k0);
        for(int j=0; j < A.rows(); ++j)
        if(k0!=j){
            S12[++cnt]  = A(k0,j);
            zstar[++cnz]= Z[j];
            mustar[cnz]   = mu[j];
        }
        Matrix<double> S22 = delrc(A,k,k);
        double mean = muj + (S12*invpd(S22)*(zstar-mustar))[0];
       return mean;
    } 
    // eveluate both mean and standard deviation in efficient way
    template <class T>
	void cndnorm(const Matrix<T> &Z, const Matrix<T> &Xbeta, 
             const Matrix<T> &invS, const int& p, 
             const int& j, double& muj, double& sdj)
	{
	//	function to compute moments of Z[j] | Z[-j]  //

		int cnt,i,j0;
		double varj;
		j0=j-1;
		cnt = p*j0;
		varj = 1./invS[cnt+j0];

		muj = 0.0;
		for (i=0 ; i < p; ++i)
			{
			if (i != j0) 
				{muj +=  - varj*invS[cnt+i]*(Z[i]-Xbeta[i]);
				}
			}
		muj=Xbeta[j0] + muj;
		sdj=::sqrt(varj);
	}

    // return A %x% ones (K)
         template <class T>
         Matrix<T> xvec(const Matrix<T>& A, const int & K)
    {
            Matrix<T> temp = A%ones<T>(K,1);
            return temp; 
    } 

     // exract the p*(p-1)/2 unique elements of a correlation matrix
        template <class T>
        Matrix<T>
        UNIK(const Matrix<T> &A)
         {
         if (! A.isSquare()) {
            throw scythe_dimension_error(__FILE__, __PRETTY_FUNCTION__,
            __LINE__, "Matrix not square");
            }
        const int K = A.rows();
        const int K2 = K*(K-1)/2;
        int cnt = -1;
        Matrix<T> temp =   Matrix<T> (K2,1);
          for(int i = 0; i < K - 1; ++i) 
            for(int j = i; j < K; ++j) 
                 if(i < j){
                    temp[++cnt] = A(i,j);
                } 
                return temp;
    }
    
 
} // end namespace SCYTHE

#endif /* SCYTHE_ADD_CC */
