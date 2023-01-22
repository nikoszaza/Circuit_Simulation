#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define THREADS_NUM 6
#define CHUNKSIZE 4000
#include "csparse.h"
#define MAX_ITOL 1e-22

double *ax, *r, *z, *m, *p, *q, *rt, *zt, *pt, *qt;


// a * b
double dot_prod(double a[], double b[],int size){
	double dot_prod = 0;
	int chunk = CHUNKSIZE, i;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(a, b, size, chunk, dot_prod) private(i)
	{

	  #pragma omp for schedule(dynamic, chunk) reduction(+:dot_prod)
	  for(i = 0; i < size; i++)
	  	dot_prod = dot_prod + a[i] * b[i];

	} /* Done */

	return dot_prod;
}

// y = Ax
void matrix_vector_mul(double **A, double *x, double *y, int size){
	int i,j, chunk = CHUNKSIZE;

	#pragma omp parallel num_threads(THREADS_NUM) shared(A,x,y,size) private(i,j)
	{
		#pragma omp for schedule(dynamic, chunk) collapse(1) 
		for(i=0; i<size; i++){
			y[i] = 0;
			for(j =0; j < size; j++)
				y[i] += A[i][j]*x[j];
		}
	}	
	
}

// y = A^Tx
void matrix_vector_mul_trans(double **A, double *x, double *y, int size){
	int i,j, chunk = CHUNKSIZE;

	#pragma omp parallel num_threads(THREADS_NUM) shared(A,x,y,size) private(i,j)
	{
		#pragma omp for schedule(dynamic, chunk) collapse(1) 
		for(i=0; i<size; i++){
			y[i] = 0;
			for(j =0; j < size; j++)
				y[i] += A[j][i]*x[j];
		}
	}	
	
}

// z = ax + by
void axpy(double* x, double a, double *y, double b, double* z, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(x,a,b,y,z,size) private(i)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++)
			z[i] = a*x[i] + b*y[i];
	}
		
}

// y = Ax + y
void gaxpy(double** A, double* x, double* y, int size){
	double* temp = (double*)malloc(sizeof(double)*size);
	
	//temp = Ax
	matrix_vector_mul(A,x,temp, size);
	
	// y = temp + y
	axpy(temp,1,y,1,y,size);
	
	free(temp);
}

// Z = aX + bY
double** add(double**X, double**Y, double a, double b, int size){
	double** Z = malloc(sizeof(double*)*size);
	for(int i = 0; i < size; i++){
		Z[i] = (double*)malloc(sizeof(double)*size);
	}
	
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			Z[i][j] = a*X[i][j] + b*Y[i][j];
		}
	}
	
	return Z;
}

// y = Ax for sparse A
int sparse_matrix_vector_mul(cs* A, double* x, double* y){
	int p, j, n, *Ap, *Ai;
    double *Ax ;
    if (!CS_CSC (A) || !x || !y) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
	for (j = 0 ; j < n ; j++){
		y[j] = 0.0;
	}
		
	for (j = 0 ; j < n ; j++){
	    for (p = Ap [j] ; p < Ap [j+1] ; p++){
	        y [Ai [p]] += Ax [p] * x [j];
	    }
	}
	
	return (1) ;
}


// y = A^Tx for sparse A
int sparse_matrix_vector_mul_trans(cs* A, double* x, double* y){
	int p, i, n, *Ap, *Ai, chunk = CHUNKSIZE ;
    double *Ax ;
    if (!CS_CSC (A) || !x || !y) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
    #pragma omp parallel num_threads(THREADS_NUM) private(i,p)
	{
	
		#pragma omp for schedule(dynamic, chunk) collapse(1)
	    for (i = 0 ; i < n ; i++)
	    {
	    	y[i] = 0.0;
	        for (p = Ap [i] ; p < Ap [i+1] ; p++)
	        {
	            y [i] += Ax [p] * x [Ai[p]];
	        }
	    }
	}
	
	return (1) ;
}

int sparse_prec(cs* A,double* m){
	int p, j, n, *Ap, *Ai, chunk = CHUNKSIZE ;
    double *Ax ;
    if (!CS_CSC (A) || !m) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
    #pragma omp parallel num_threads(THREADS_NUM) private(j,p)
    {
    	#pragma omp for schedule(dynamic, chunk) collapse(1)
	    for (j = 0 ; j < n ; j++)
	    {
	        for (p = Ap [j] ; p < Ap [j+1] ; p++)
	        {
	            if(Ai[p] == j){
	            	if( (Ax[p] < -1e-12) || (Ax[p] > 1e-12)){
	            		m[j] = 1/Ax[p];
					}
					else{
						m[j] = 1;
					}
	            	break;
				}
	        }
	        if(p == Ap[j+1]){
	    		m[j] = 1;
			}
	    } 
	}
	
	return(1);
}

// euclidian norm
double norm2(double* x, int size){
	int i, chunk = CHUNKSIZE;
	double norm = 0;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(x,size) private(i)
	{
		#pragma omp for schedule(dynamic, chunk) reduction(+:norm)
		for(i = 0; i<size; i++)
			norm = norm + x[i]*x[i];
	}
	
	return sqrt(norm);
}

// get diagonal as preconditioner
void prec(double **A, double* m, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(A, m)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++){
			if( (A[i][i] > 1e-15) || (A[i][i] < -1e-15))
				m[i] = 1/A[i][i];
			else
				m[i] = 1;
		}
	}
}

// solve Mr = z, where is M the inverse of a diagonal matrix
void diag_solve(double* m, double* r, double* z, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(m,r,z)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++)
			z[i] = m[i]*r[i];
	}
}

// cpy from x to y
void vec_cpy(double *x, double *y, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(x,y)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++)
			y[i] = x[i];
	}
}

//for cg and sparse cg, remember for sparse size = n
void cg_calloc(int size){
    ax = (double*)calloc(sizeof(double),size);
	r = (double*)calloc(sizeof(double),size);
	z = (double*)calloc(sizeof(double),size);
	m = (double*)calloc(sizeof(double),size);
	p = (double*)calloc(sizeof(double),size);
	q = (double*)calloc(sizeof(double),size);
}
//for cg and sparse cg, remember for sparse size = n
void free_cg(){
    free(ax);
	free(r);
	free(z);
	free(m);
	free(p);
	free(q);
}

// conjugate gradient
// returns the number of iterations
int cg(double **A, double *x, double *b, double itol, int size ){
	double rnorm, bnorm = norm2(b, size), rho, rho1, beta, alpha;
	int iter = 0;
	
	//init
	matrix_vector_mul(A, x, ax, size);
	prec(A,m,size);
	axpy(b, 1, ax, -1, r, size);

	
	if(bnorm < 1e-14)
		bnorm = 1;
			
	rnorm = norm2(r, size);
	//main flow control	
	while( (rnorm/bnorm > itol) && (iter < size)){
		iter++;
		diag_solve(m, r, z, size);
		rho = dot_prod(r, z, size);
	
		if(iter == 1){
			vec_cpy(z,p, size);
		}
		else{
			beta = rho/rho1;
			axpy(z, 1, p, beta, p, size);
		}
		rho1 = rho;
		matrix_vector_mul(A, p, q, size);
		alpha = rho / dot_prod(p,q,size);
		axpy(x,1,p,alpha,x,size);
		axpy(r,1,q,(-alpha),r,size);

		rnorm = norm2(r, size);		
	}
	
	
	return iter;
}


int sparse_cg(cs* A, double *x, double *b, double itol){
	if (!CS_CSC (A)){
		printf("Error: A is not CSC!\n");
		return -1;
	}
	
	int n = A->n;
	
	double rnorm, bnorm = norm2(b, n), rho, rho1, beta, alpha;
	int iter = 0;
	
	sparse_prec(A,m);
	//init
	sparse_matrix_vector_mul(A, x, ax);
	axpy(b, 1, ax, -1, r, n);
	
	if(bnorm < 1e-14)
		bnorm = 1;
			
	rnorm = norm2(r, n);
	//main flow control	
	while( (rnorm/bnorm > itol) && (iter < n)){
		iter++;
		diag_solve(m, r, z, n);
		rho = dot_prod(r, z, n);
	
		if(iter == 1){
			vec_cpy(z,p, n);
		}
		else{
			beta = rho/rho1;
			axpy(z, 1, p, beta, p, n);
		}
		rho1 = rho;
		sparse_matrix_vector_mul(A, p, q);
		alpha = rho / dot_prod(p,q,n);
		axpy(x,1,p,alpha,x,n);
		axpy(r,1,q,(-alpha),r,n);

		rnorm = norm2(r, n);		
	}
	
      
	
	return iter;	
}
//for bicg and sparse bicg, remember for sparse size = n
void bicg_calloc(int size){
    ax = (double*)calloc(sizeof(double),size);
	r = (double*)calloc(sizeof(double),size);
	rt = (double*)calloc(sizeof(double),size);
	z = (double*)calloc(sizeof(double),size);
	zt = (double*)calloc(sizeof(double),size);
	m = (double*)calloc(sizeof(double),size);
	p = (double*)calloc(sizeof(double),size);
	pt = (double*)calloc(sizeof(double),size);
	q = (double*)calloc(sizeof(double),size);
	qt = (double*)calloc(sizeof(double),size);
	
}
//for bicg and sparse bicg, remember for sparse size = n
void free_bicg(){
    free(ax);
	free(r);
	free(rt);
	free(zt);
	free(z);
	free(m);
	free(p);
	free(pt);
	free(q);
	free(qt);
}

int bicg(double **A, double *x, double *b, double itol, int size ){
	double rnorm, bnorm = norm2(b, size), rho, rho1, beta, alpha, omega;
	int iter = 0;
	
	prec(A,m, size);
	//init
	matrix_vector_mul(A, x, ax, size);
	
	axpy(b, 1, ax, -1, r, size);
	vec_cpy(r, rt, size);
	
	if(bnorm < 1e-14)
		bnorm = 1;
			
	rnorm = norm2(r, size);
	//main flow control
	while( (rnorm/bnorm > itol) && (iter < size)){
		iter++;
		diag_solve(m, r, z, size);
		diag_solve(m, rt, zt, size);
		rho = dot_prod(rt, z, size);
		
		if((rho < MAX_ITOL) && (rho > -MAX_ITOL)){
			//printf("Error(1): bi-cg cannot continue\n");
			return - 1;
		}
		
		if(iter == 1){
			vec_cpy(z,p, size);
			vec_cpy(zt,pt, size);
		}
		else{
			beta = rho/rho1;
			axpy(z, 1, p, beta, p, size);
			axpy(zt, 1, pt, beta, pt, size);
		}
		rho1 = rho;
		matrix_vector_mul(A, p, q, size);
		matrix_vector_mul_trans(A, pt, qt, size);
		omega = dot_prod(pt, q, size);
		if((omega < MAX_ITOL) && (omega > -MAX_ITOL)){
			//printf("Error(2): bi-cg cannot continue\n");
			return -1;
		}
		alpha = rho / omega;
		axpy(x,1,p,alpha,x,size);
		axpy(r,1,q,(-alpha),r,size);
		axpy(rt,1,qt,(-alpha),rt,size);
		
		rnorm = norm2(r, size);
		
	}
	
	
	return iter;
}


int sparse_bicg(cs* A, double *x, double *b, double itol ){
	if (!CS_CSC (A)){
		printf("Error: A is not CSC!\n");
		return -1;
	}
	int n = A->n;
	
	double rnorm, bnorm = norm2(b, n), rho, rho1, beta, alpha, omega;
	int iter = 0;
	
	sparse_prec(A,m);
	//init
	sparse_matrix_vector_mul(A, x, ax);
	
	axpy(b, 1, ax, -1, r, n);
	vec_cpy(r, rt, n);
	
	if(bnorm < 1e-18)
		bnorm = 1;
			
	rnorm = norm2(r, n);
	//main flow control	
	while( (rnorm/bnorm > itol) && (iter < n)){
		iter++;
		diag_solve(m, r, z, n);
		diag_solve(m, rt, zt, n);
		rho = dot_prod(rt, z, n);
		
		if((rho < MAX_ITOL) && (rho > -MAX_ITOL)){
			//printf("Error(1): bi-cg cannot continue\n");
			return -1;
		}
		
		if(iter == 1){
			vec_cpy(z,p, n);
			vec_cpy(zt,pt, n);
		}
		else{
			beta = rho/rho1;
			axpy(z, 1, p, beta, p, n);
			axpy(zt, 1, pt, beta, pt, n);
		}
		rho1 = rho;
		sparse_matrix_vector_mul(A, p, q);
		sparse_matrix_vector_mul_trans(A, pt, qt);
		omega = dot_prod(pt, q, n);
		if((omega < MAX_ITOL) && (omega > -MAX_ITOL)){
			//printf("Error(2): bi-cg cannot continue\n");
			return -1;
		}
		alpha = rho / omega;
		axpy(x,1,p,alpha,x,n);
		axpy(r,1,q,(-alpha),r,n);
		axpy(rt,1,qt,(-alpha),rt,n);
		
		rnorm = norm2(r, n);
		
	}
	
	
	
	return iter;
}
