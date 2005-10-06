/* baymvb.cc
 *
 * This file provides implementations of functions for fitting
 * a Bayesian Multivariate Logistic Regression Model.
 *
 * MCMCmvlogit.cc
 * Copyright (C) 2005, Samuel M. Mwalili
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  A copy of this license is included
 * with this library (LICENSE.GPL).
 *
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
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
 * Samuel.Mwalili@med.kuleuven.be
 */

//inlcude header files

//inlcude header files
#include <iostream>
#include <math.h>
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "add.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
             
//define constant
#define ML_PI     3.14159265
#define ML_NU     7.30000000 
#define ML_S2     2.38853440 
#define ML_MIN(a,b) a<b?a:b

//namespapces
using namespace SCYTHE;
using namespace std;

    	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   	 //        function prototytpes
    	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
	//---------------------------------------------------------------------------
	//-                     posterior density of R                              -
	//---------------------------------------------------------------------------
           double logp_R(const Matrix<double>& R, const Matrix<double>& R0,const Matrix<double>& G0,
    			  const Matrix<double> *Xarr, const Matrix<double>& beta, 
    			  const Matrix<double> *Zarr,const Matrix<double>&PHI, 
    			  const double& zeta, const double& sd)
	  {   
		  // define constants
          	const int N = PHI.rows();
          	const int J = R.rows();
		const double q = zeta/(::sqrt(2)*sd);
	
 		// log eigen value contribution
		double log_ee = ::log(pnorm(q)-pnorm(-q));

		//log prior for R contribution
		double log_pi_R = -0.5*((t(UNIK(R)-UNIK(R0))*G0*(UNIK(R)-UNIK(R0)))[0]);

		//log mvn likelihood contribution
		double log_SSE = 0.0;
		const Matrix<double> invR  = invpd(R);
		const double lndR = log(abs(det(R)));
	
   		for(int i=0; i<N; i++){
		 	Matrix<double> invwR = (PHI[i]/ML_S2)*invR;
		 	Matrix<double> Xb = Xarr[i] * beta;
   		 	log_SSE -= static_cast<double>(J)*log(ML_S2/PHI[i]) +  lndR + (t(Zarr[i]-Xb)*invwR*(Zarr[i]-Xb))[0];
			}
		  return log_ee + log_pi_R + 0.5*log_SSE;
	  }

    extern "C"{
  
      void baymvb(double *sampledata, const int *samplerow,const int *samplecol, 
 		const double *Ydata,const int *Yrow, const int *Ycol, 
 		const double *Xdata,const int *Xrow, const int *Xcol, 
 		const int *burnin, 	const int *mcmc, const int *thin, 
 		const int *lecuyer, const int *seedarray, 
 		const int *lecuyerstream,
 		const double *betastartdata, const int *betastartrow, const int *betastartcol,
 		const double *Rstartdata, const int *Rstartrow, const int *Rstartcol, 
		const double *b0data, const int *b0row, const int *b0col, 
		const double *B0data, const int *B0row, const int *B0col, 
		const double *R0data, const int *R0row, const int *R0col, 
		const double *G0data, const int *G0row, const int *G0col, 
		const int* nn, const int* nvar, const int* np,  
		const int* distr, const double *std) {
		 
		// pull together Matrix objects
		const Matrix <double> Y = r2scythe(*Yrow, *Ycol, Ydata);
		const Matrix <double> X = r2scythe(*Xrow, *Xcol, Xdata);
		Matrix <double> beta = r2scythe(*betastartrow, *betastartcol, betastartdata);
		Matrix <double> R = r2scythe(*Rstartrow, *Rstartcol, Rstartdata);

		const Matrix <double> b0 = r2scythe(*b0row, *b0col, b0data);
		const Matrix <double> B0 = r2scythe(*B0row, *B0col, B0data);
    		const Matrix <double> R0 = r2scythe(*R0row, *R0col, R0data);
		const Matrix <double> G0 = r2scythe(*G0row, *G0col, G0data);

		// define constants
		const int tot_iter = *burnin + *mcmc;  // total number of mcmc iterations
		const int nstore = *mcmc / *thin;      // number of draws to store
		const int P = np[0];                   // # of regression coefficients
		const int J = nvar[0];                 // # of variables in the model
 		const int N = nn[0];                   // # of subjects in the model
		const int K = P + J*(J-1)/2;           // total # parameters to store
             	const double sd = std[0];
             	const int uniks       = J*(J-1)/2;
             
		// create arrays of matrices for data access
   		  Matrix<double>* Yarr = new Matrix<double>[N]; // Matrix<double> Yarr[Mn];
   		  Matrix<double>* Xarr = new Matrix<double>[N]; //Matrix<double> Xarr[Mn];
   		  for(int i = 0; i < N; ++i) {
	 			 Yarr[i] = Matrix<double>(J,1);
	   		     Xarr[i] = Matrix<double>(J, P, 0.0);
	   		     int start = i * J;
	   		     for(int j=0; j<J; ++j){
				     Yarr[i](j,0) = Y(start + j,0);
				     for(int m=0; m<P; ++m) {
			   		   Xarr[i](j,m) = X(start + j, m);    
				     }
	   		     }
   		 }
 
		// storage matrix or matrices
		Matrix<double> storemat(nstore, K);
		 
		// initialize rng stream
		 rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

		// starting values
		Matrix<double> PHI = ones<double>(N,1); 
		if(!(*distr)) PHI = ML_S2*PHI;
		 
		Matrix<double>* Zarr = new Matrix<double>[N];
 			for(int i = 0; i < N; ++i) {
	 		 	Zarr[i] = Matrix<double>(J,1);
  		     		for(int j = 0; j < J; ++j) {
		   			if (Yarr[i](j,0) == 1.0){
		 				Zarr[i](j,0) = ::fabs(stream->rnorm(0.0,1.0));
	  				}
	  				if (Yarr[i](j,0) == 0.0){
	    				Zarr[i](j,0) = -::fabs(stream->rnorm(0.0,1.0));
					}
					if (Yarr[i](j,0) != 1.0 && Yarr[i](j,0) != 0.0){
	    				Zarr[i](j,0) = stream->rnorm(0.0,1.0);
					}
     				}  			      
 		 	}
 			 
       
       		int count = 0;
       		int accepts = 0; 
       		
       		int progress = 1;
		int itemp =  static_cast<int> (floor((double)tot_iter/10));
       		
       		for (int iter = 0; iter < tot_iter; ++iter){
       							 
               	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          	//                    sample Z | X, R, phi
          	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          			
          	const double PHI_alpha = (ML_NU + J)/2.0;
       		Matrix<double> beta_var_sum(P,P);
               	Matrix<double> beta_mean_sum(P,1);
               	Matrix<double> invR = invpd(R);
                double z_mu, z_var;
       			for (int i=0; i<N; ++i){
       				const Matrix<double> Z_mean = Xarr[i] * beta;	
       			   	const Matrix<double> invwR = (PHI[i]/ML_S2)*invR;
       				for(int j = 0; j < J; ++j) {
       					cndnorm(Zarr[i],Z_mean,invwR, J, j + 1, z_mu, z_var);	 		
       			  		if (Yarr[i](j,0) == 1.0){
       					Zarr[i](j,0) = stream->rtbnorm_combo(z_mu, z_var, 0); 
       					if(isinf(Zarr[i](j,0)))
       						Zarr[i](j,0) = 1;
       					}
       					if (Yarr[i](j,0) == 0.0){
       					Zarr[i](j,0) = stream->rtanorm_combo(z_mu, z_var, 0); 
       					if(isinf(Zarr[i](j,0)))
       						Zarr[i](j,0) = -1;
       					}
       					if (Yarr[i](j,0) != 1.0 && Yarr[i](j,0) != 0.0)
       					Zarr[i](j,0) = stream->rnorm(z_mu, ::sqrt(z_var));
       				}//e.o.j
       				
       
       		//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                //            Sample     phi | Z, X, R
               	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               	 if(*distr){
       		 const Matrix<double> tZRZ  = t(Zarr[i] - Xarr[i] * beta)*invR*(Zarr[i] - Xarr[i] * beta);
       		 const double PHI_beta  = (ML_NU + (1/ML_S2)*tZRZ[0])/2.0;
       		  PHI[i] = stream->rgamma(PHI_alpha,PHI_beta); 
                  }
               	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               	//            Sample     beta | Z, X, R
               	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
               	beta_var_sum = beta_var_sum   + PHI[i]*(t(Xarr[i]) * invR * Xarr[i]);
               	beta_mean_sum = beta_mean_sum + PHI[i]*(t(Xarr[i]) * invR * Zarr[i]);
               		
               	}//e.o.i
               	
               	Matrix<double> beta_sim_var  = invpd(B0 + (1/ML_S2) * beta_var_sum);
               	Matrix<double> beta_sim_mean = beta_sim_var * (B0 * b0 + (1/ML_S2) * beta_mean_sum);
               	               beta          = beta_sim_mean + cholesky(beta_sim_var) * stream->rnorm(P,1);
        
  	      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	      //            Sample     R | Z, X
	      //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	      
	       //(i) sample standard normal variates
	    	Matrix<double> z = stream->rnorm(uniks,1);
	    	
	      //(ii) generate signed distance   
       		Matrix<double> zvalues =  Matrix<double>(J,1);
		Matrix<double> zvector =  eye<double>(J);
		     
		JACOBI(R,zvalues,zvector); // eigen values/vectors
		
		const double zeta    = min(zvalues);           // min eigen value
		const double d_lower = -zeta/::sqrt(2);
		const double d_upper =  zeta/::sqrt(2);
	       
		const double d = stream->rtnorm(0.0,sd*sd,d_lower,d_upper);
		     
	       //(iii) calculate H
	    	Matrix<double> H = 0*eye<double>(J);
	      	double zss = (t(z)*z)[0];
			int cnt = -1;      	
	      	for(int i = 0; i < J - 1; ++i) 
		      	for(int j = i; j < J; ++j) 
			     	 if(i < j){
				      	H(i,j) = d*z[++cnt]/::sqrt(zss); 
				      	H(j,i) = H(i,j); 
			      	}		      	
             
		// (iv) draw candidate
	      	const Matrix<double> R_can = R + H;
       		Matrix<double> zval_star =  Matrix<double>(J,1);
		Matrix<double> zvec_star =  eye<double>(J);

	        JACOBI(R_can,zval_star,zvec_star);  	// eigen values/vectors
		double zsta  = min(zval_star);  	// min eigen value
		    
              
	        if(ispd(R_can)){
	         const  double logp_R_cur = logp_R(R    ,R0, G0, Xarr,beta, Zarr, PHI, zeta, sd);
	         const  double logp_R_can = logp_R(R_can,R0, G0, Xarr,beta, Zarr, PHI, zsta, sd);
	      	 const double ratio = ML_MIN(::exp(logp_R_can - logp_R_cur),1.0); 
	      	 
	      	 
	   	if ((stream->runif() < ratio)){
			R = R_can;
		   	++accepts;
	     	 	}
	       }
	       
	       // save values
	       if (iter >= *burnin && (iter%*thin == 0)) {
	       		for(int j = 0; j < P; ++j) {
           			storemat(count,j) = beta[j];
            		}
	       		const Matrix<double> unik_R = UNIK(R);
	       		for(int j = 0; j < uniks; ++j) {
            			storemat(count,j+P) = unik_R(j,0);
	       		}
          		++count;
         	}

 	     if((iter + 1) >= itemp){ 
 		Rprintf("\tMCMC sampling %3i percent complete ... [MH-rate = %3.5f]\n", progress*10,
 		100.0*static_cast<double>(accepts) / static_cast<double>(iter + 1));
 		itemp +=  static_cast<int> (floor((double)tot_iter/10));
 		progress += 1;
 	    } 
			
         }//iter
	
	// print the the acceptance rate to the console 
	  Rprintf("\n------------------------------------------------------\n");
	  if(!(*distr))  Rprintf("... Posterior draws from multivariate probit model ...\n");
	  if(*distr)     Rprintf("... Posterior draws from multivariate t-link model ...\n");
	  Rprintf("... The Metropolis acceptance rate was %3.5f    ...", 
	   100.0*static_cast<double>(accepts) / static_cast<double>(tot_iter));
	  Rprintf("\n------------------------------------------------------\n");

	 R_CheckUserInterrupt(); // allow user interrupt  
	
	
	delete stream; // clean up random number stream
 	
	// return posterior denisty sample to R
  	int loop = samplerow[0] * samplecol[0];
  	for (int i=0; i<loop; ++i) {
     	 	sampledata[i] = storemat[i];
      	}
        
    	delete [] Zarr;
   	delete [] Xarr;
  	delete [] Yarr;
     } // baymvb
 
} //extern C
