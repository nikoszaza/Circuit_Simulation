#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define THREADS_NUM 1
#define CHUNKSIZE 4000
#include "csparse.h"
#include "CXSparse/Include/cs.h"
#define MAX_ITOL 1e-22

double *ax, *r, *z, *m, *p, *q, *rt, *zt, *pt, *qt;
double complex *ax_complex, *r_complex, *z_complex, *m_complex, *p_complex, *q_complex, *rt_complex, *zt_complex, *pt_complex, *qt_complex;

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


double complex dot_prod_complex(double complex a[], double complex b[],int size){
	double complex dot_prod = 0;
	int chunk = CHUNKSIZE, i;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(a, b, size, chunk, dot_prod) private(i)
	{

	  #pragma omp for schedule(dynamic, chunk) reduction(+:dot_prod)
	  for(i = 0; i < size; i++)
	  	dot_prod = dot_prod + conj(a[i])* b[i];

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


void matrix_vector_mul_complex(double complex **A, double complex *x, double complex *y, int size){
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

void matrix_vector_mul_trans_complex(double complex **A, double complex *x, double complex *y, int size){
	int i,j, chunk = CHUNKSIZE;

	#pragma omp parallel num_threads(THREADS_NUM) shared(A,x,y,size) private(i,j)
	{
		#pragma omp for schedule(dynamic, chunk) collapse(1) 
		for(i=0; i<size; i++){
			y[i] = 0;
			for(j =0; j < size; j++)
				y[i] += conj(A[j][i])*x[j];
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


void axpy_complex(double complex* x, double complex a, double complex *y, double complex b, double complex* z, int size){
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


void gaxpy_complex(double complex** A, double complex* x, double complex* y, int size){
	double complex* temp = (double complex*)malloc(sizeof(double complex)*size);
	
	//temp = Ax
	matrix_vector_mul_complex(A,x,temp, size);
	
	// y = temp + y
	axpy_complex(temp,1,y,1,y,size);
	
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


// Z = aX + bY
double complex** add_complex(double complex **X, double complex **Y, double complex a, double complex b, int size){
	double complex** Z = malloc(sizeof(double complex*)*size);
	for(int i = 0; i < size; i++){
		Z[i] = (double complex*)malloc(sizeof(double complex)*size);
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


int sparse_matrix_vector_mul_complex(cs_cl* A, double complex* x, double complex* y){
	int p, j, n;
	long int *Ap, *Ai;
    double complex *Ax ;
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

int sparse_matrix_vector_mul_trans_complex(cs_cl* A, double complex* x, double complex* y){
	int p, i, n, chunk = CHUNKSIZE ;
	long int *Ap, *Ai;
    double complex *Ax ;
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
    #pragma omp parallel num_threads(THREADS_NUM) private(i,p)
	{
	
		#pragma omp for schedule(dynamic, chunk) collapse(1)
	    for (i = 0 ; i < n ; i++)
	    {
	    	y[i] = 0.0;
	        for (p = Ap [i] ; p < Ap [i+1] ; p++)
	        {
	            y [i] += conj(Ax [p])* x [Ai[p]];
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


int sparse_prec_complex(cs_cl* A,double complex* m){
	int p, j, n, chunk = CHUNKSIZE ;
	long int *Ap, *Ai;
    double complex *Ax ;
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    
    #pragma omp parallel num_threads(THREADS_NUM) private(j,p)
    {
    	#pragma omp for schedule(dynamic, chunk) collapse(1)
	    for (j = 0 ; j < n ; j++)
	    {
	        for (p = Ap [j] ; p < Ap [j+1] ; p++)
	        {
	            if(Ai[p] == j){
	            	if( cabs(Ax[p]) > 1e-12 ){
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

double norm2_complex(double complex* x, int size){
	int i, chunk = CHUNKSIZE;
	double norm = 0;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(x,size) private(i)
	{
		#pragma omp for schedule(dynamic, chunk) reduction(+:norm)
		for(i = 0; i<size; i++)
			norm = norm + cabs(x[i])*cabs(x[i]);
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


void prec_complex(double complex **A, double complex* m_complex, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(A, m_complex)
	{
		for(i = 0; i<size; i++){
			if( cabs(A[i][i]) > 1e-12 )
				m_complex[i] = 1/A[i][i];
			else
				m_complex[i] = 1;
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


// solve Mr = z, where is M the inverse of a diagonal matrix
void diag_solve_complex(double complex* m_complex, double complex* r, double complex* z, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(m,r,z)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++)
			z[i] = m_complex[i]*r[i];
	}
}


void diag_solve_tran_complex(double complex* m, double complex* r, double complex* z, int size){
	int i, chunk = CHUNKSIZE;
	
	#pragma omp parallel num_threads(THREADS_NUM) shared(m,r,z)
	{
		#pragma omp for schedule(dynamic, chunk)
		for(i = 0; i<size; i++)
			z[i] = conj(m[i])*r[i];
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


void vec_cpy_complex(double complex *x, double complex *y, int size){
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

void bicg_calloc_complex(int size){
    ax_complex = (double complex*)calloc(sizeof(double complex),size);
	r_complex = (double complex*)calloc(sizeof(double complex),size);
	rt_complex = (double complex*)calloc(sizeof(double complex),size);
	z_complex = (double complex*)calloc(sizeof(double complex),size);;
	zt_complex = (double complex*)calloc(sizeof(double complex),size);
	m_complex = (double complex*)calloc(sizeof(double complex),size);;
	p_complex = (double complex*)calloc(sizeof(double complex),size);
	pt_complex = (double complex*)calloc(sizeof(double complex),size);
	q_complex = (double complex*)calloc(sizeof(double complex),size);
	qt_complex = (double complex*)calloc(sizeof(double complex),size);
	
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


void free_bicg_complex(){
	free(ax_complex);
	free(r_complex);
	free(rt_complex);
	free(zt_complex);
	free(z_complex);
	free(m_complex);
	free(p_complex);
	free(pt_complex);
	free(q_complex);
	free(qt_complex);
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


int bicg_complex(double complex **A, double complex *x, double complex *b, double itol, int size ){
	double rnorm, bnorm = norm2_complex(b, size);
	double complex rho, rho1, beta, alpha, omega;
	int iter = 0;
	
	prec_complex(A,m_complex, size);
	//init
	matrix_vector_mul_complex(A, x, ax_complex, size);
	
	axpy_complex(b, 1, ax_complex, -1, r_complex, size);
	vec_cpy_complex(r_complex, rt_complex, size);
	
	if(bnorm < 1e-14)
		bnorm = 1;
			
	rnorm = norm2_complex(r_complex, size);
	//main flow control
	while( (rnorm/bnorm > itol) && (iter < size)){
		iter++;
		diag_solve_complex(m_complex, r_complex, z_complex, size);
		diag_solve_tran_complex(m_complex, rt_complex, zt_complex, size);
		rho = dot_prod_complex(rt_complex, z_complex, size);
		
		if(cabs(rho) < MAX_ITOL){
			printf("Error(1): bi-cg cannot continue\n");
			return - 1;
		}
		//printf("Norm %f\n", rnorm/bnorm);
		
		if(iter == 1){
			vec_cpy_complex(z_complex,p_complex, size);
			vec_cpy_complex(zt_complex,pt_complex, size);
		}
		else{
			beta = rho/rho1;
			axpy_complex(z_complex, 1, p_complex, beta, p_complex, size);
			axpy_complex(zt_complex, 1, pt_complex, conj(beta), pt_complex, size);
		}
		rho1 = rho;
		matrix_vector_mul_complex(A, p_complex, q_complex, size);
		matrix_vector_mul_trans_complex(A, pt_complex, qt_complex, size);
		omega = dot_prod_complex(pt_complex, q_complex, size);
		if(cabs(omega)< MAX_ITOL){
			printf("Error(2): bi-cg cannot continue\n");
			return -1;
		}
		alpha = rho / omega;
		axpy_complex(x,1,p_complex,alpha,x,size);
		axpy_complex(r_complex,1,q_complex,(-alpha),r_complex,size);
		axpy_complex(rt_complex,1,qt_complex,(-conj(alpha)),rt_complex,size);
		
		rnorm = norm2_complex(r_complex, size);
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


int sparse_bicg_complex(cs_cl* A, double complex *x, double complex *b, double itol ){
	if (!CS_CSC (A)){
		printf("Error: A is not CSC!\n");
		return -1;
	}
	int n = A->n;
	double rnorm, bnorm = norm2_complex(b, n);
	double complex rho, rho1, beta, alpha, omega;
	int iter = 0;
	
	
	//init
	sparse_matrix_vector_mul_complex(A, x, ax_complex);
	sparse_prec_complex(A,m_complex);
	
	axpy_complex(b, 1, ax_complex, -1, r_complex, n);
	vec_cpy_complex(r_complex, rt_complex, n);
	
	if(bnorm < 1e-14)
		bnorm = 1;
			
	rnorm = norm2_complex(r_complex, n);
	//main flow control
	while( (rnorm/bnorm > itol) && (iter < n)){
		iter++;
		diag_solve_complex(m_complex, r_complex, z_complex, n);
		diag_solve_tran_complex(m_complex, rt_complex, zt_complex, n);
		rho = dot_prod_complex(rt_complex, z_complex, n);
		
		if(cabs(rho) < MAX_ITOL){
			printf("Error(1): bi-cg cannot continue\n");
			return - 1;
		}
		//printf("Norm %f\n", rnorm/bnorm);
		
		if(iter == 1){
			vec_cpy_complex(z_complex,p_complex, n);
			vec_cpy_complex(zt_complex,pt_complex, n);
		}
		else{
			beta = rho/rho1;
			axpy_complex(z_complex, 1, p_complex, beta, p_complex, n);
			axpy_complex(zt_complex, 1, pt_complex, conj(beta), pt_complex, n);
		}
		rho1 = rho;
		sparse_matrix_vector_mul_complex(A, p_complex, q_complex);
		sparse_matrix_vector_mul_trans_complex(A, pt_complex, qt_complex);
		omega = dot_prod_complex(pt_complex, q_complex, n);
		if(cabs(omega) < MAX_ITOL){
			printf("Error(2): bi-cg cannot continue\n");
			return -1;
		}
		alpha = rho / omega;
		axpy_complex(x,1,p_complex,alpha,x,n);
		axpy_complex(r_complex,1,q_complex,(-alpha),r_complex,n);
		axpy_complex(rt_complex,1,qt_complex,(-conj(alpha)),rt_complex,n);
		
		rnorm = norm2_complex(r_complex, n);
	}
	
	
	return iter;
}
