#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <complex.h>
#include "hash.h"
#include "list.h"
#include "direct_methods.h"
#include "blas.h"
#include "csparse.h"
#include "CXSparse/Include/cs.h"
#define INITIAL_SIZE 10000
#define IMPORT_LINE_SIZE 1024
#define MAX_LEN 36
#define MAX_LEN_SP 96

typedef struct {
	double** A; //array for i,v,l,r
	double** B; //array for l,c
	cs* B_sp; //sparse version of B
	cs* A_sp; //sparse version of A
	cs* C; //compressed column version of A_sp
	cs* D; //compressed column version of B_sp
	double* b; //dc operation point b
	double complex *b_complex; // complex b, used for AC
	double* x0; //dc solution/ initial tran solution
	int m1, m2; //m1: #nodes - 1 (gnd), m2: #i,l
}dc_mna_t;

//arrays for sparse direct methods
css *S = NULL;
csn *N = NULL;
cs_cls *S_complex = NULL;
cs_cln *N_complex = NULL;


//splits the src_str to two new strings str1, str2, based on the delimeter
void split_string(char* src_str, char* str1, char* str2, char delimeter){
	int i,j, str_len = strlen(src_str);
	
	for(i = 0; i < str_len; i++){
		if(src_str[i] == delimeter)
			break;
		
		str1[i]	= src_str[i];	
	}
	
	str1[i] = '\0';
	
	for(j = i+1; j < str_len; j++){
		str2[j - i - 1] = src_str[j];
	}
	
	str2[j- i - 1] = '\0';
}


dc_mna_t create_mna_dc(list_t* list, ht_t* ht, int m1, int m2, bool is_sparse){
	int m2_counter = 0, ka=0, kb = 0;
	dc_mna_t dc_mna;
	dc_mna.m1 = m1;
	dc_mna.m2 = m2;
	dc_mna.b = (double*)calloc((m1+m2), sizeof(double));
	dc_mna.b_complex = (double complex*)calloc((m1+m2), sizeof(double complex));
	dc_mna.x0 = (double*)calloc((m1+m2), sizeof(double));
	node_t* curr = list->head;
	
	printf("Creating dc mna...\n");
	if(!is_sparse){
		dc_mna.A = (double**)calloc((m1+m2), sizeof(double*));
		dc_mna.B = (double**)calloc((m1+m2), sizeof(double*));
		for(int i=0; i < (m1+m2); i++){
			dc_mna.A[i] = (double*)calloc((m1+m2), sizeof(double));
			dc_mna.B[i] = (double*)calloc((m1+m2), sizeof(double));
		}
		
		while(curr != NULL){
			int node1_code = lookup(ht, curr->node1);
			int node2_code = lookup(ht, curr->node2);
			switch(curr->el){
				case 'v':
					if(node1_code >= 0){
						if(node2_code >= 0){
							
							dc_mna.A[dc_mna.m1 + m2_counter][node1_code] += 1;
							dc_mna.A[dc_mna.m1 + m2_counter][node2_code] -= 1;
							dc_mna.A[node2_code][dc_mna.m1 + m2_counter] -= 1;
							dc_mna.A[node1_code][dc_mna.m1 + m2_counter] += 1;
							
						} else {
							dc_mna.A[dc_mna.m1 + m2_counter][node1_code] += 1;
							dc_mna.A[node1_code][dc_mna.m1 + m2_counter] += 1;
						}
					} else{
						if(node2_code >= 0){
							dc_mna.A[dc_mna.m1 + m2_counter][node2_code] -= 1;
							dc_mna.A[node2_code][dc_mna.m1 + m2_counter] -= 1;
						}
						//else add nothing
					}
					dc_mna.b[dc_mna.m1 + m2_counter] += curr->val;
					dc_mna.b_complex[dc_mna.m1 + m2_counter] += curr->val_complex;
					m2_counter++;
					break;
				case 'i':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							dc_mna.b[node1_code] -= curr->val;
							dc_mna.b[node2_code] += curr->val;
							dc_mna.b_complex[node1_code] -= curr->val_complex;
							dc_mna.b_complex[node2_code] += curr->val_complex;
						} else {
							dc_mna.b[node1_code] -= curr->val;
							dc_mna.b_complex[node1_code] -= curr->val_complex;
						}
					} else {
						if(node2_code >=0) {
							dc_mna.b[node2_code] += curr->val;
							dc_mna.b_complex[node2_code] += curr->val_complex;
						}
					}
					break;
				case 'r':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							dc_mna.A[node1_code][node1_code] += 1/curr->val;
							dc_mna.A[node1_code][node2_code] -= 1/curr->val;
							dc_mna.A[node2_code][node1_code] -= 1/curr->val;
							dc_mna.A[node2_code][node2_code] += 1/curr->val;
						} else {
							dc_mna.A[node1_code][node1_code] += 1/curr->val;
						}
					} else {
						if(node2_code >=0) {
							dc_mna.A[node2_code][node2_code] += 1/curr->val;
						}
					}
					break;
				case 'c':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							dc_mna.B[node1_code][node1_code] += curr->val;
							dc_mna.B[node1_code][node2_code] -= curr->val;
							dc_mna.B[node2_code][node1_code] -= curr->val;
							dc_mna.B[node2_code][node2_code] += curr->val;
						} else {
							dc_mna.B[node1_code][node1_code] += curr->val;
						}
					} else {
						if(node2_code >=0) {
							dc_mna.B[node2_code][node2_code] += curr->val;
						}
					}
					break;
				case 'l':
					if(node1_code >= 0){
						if(node2_code >= 0){
							dc_mna.A[dc_mna.m1 + m2_counter][node1_code] += 1;
							dc_mna.A[dc_mna.m1 + m2_counter][node2_code] -= 1;
							dc_mna.A[node2_code][dc_mna.m1 + m2_counter] -= 1;
							dc_mna.A[node1_code][dc_mna.m1 + m2_counter] += 1;
						} else {
							dc_mna.A[dc_mna.m1 + m2_counter][node1_code] += 1;
							dc_mna.A[node1_code][dc_mna.m1 + m2_counter] += 1;
						}
					} else{
						if(node2_code >= 0){
							dc_mna.A[dc_mna.m1 + m2_counter][node2_code] -= 1;
							dc_mna.A[node2_code][dc_mna.m1 + m2_counter] -= 1;
						}
						//else add nothing
					}
					dc_mna.B[dc_mna.m1 + m2_counter][dc_mna.m1 + m2_counter] = -curr->val;
					m2_counter++;
					break;
				default:
					break;
			}
			
			curr = curr->nxt;
		}
	}
	else {
		int nz = get_nz(*list,'A');
		dc_mna.A_sp=cs_spalloc(m1+m2,m1+m2,nz,1,1);
		dc_mna.A_sp->nz=nz;
		
		nz = get_nz(*list,'B');
		dc_mna.B_sp=cs_spalloc(m1+m2,m1+m2,nz,1,1);
		dc_mna.B_sp->nz=nz;
		
		while(curr != NULL){
			int node1_code = lookup(ht, curr->node1);
			int node2_code = lookup(ht, curr->node2);
			switch(curr->el){
				case 'v':
					if(node1_code >= 0){
						if(node2_code >= 0){
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = -1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = -1;
							ka++;		
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							
						} else {
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
						}
					} else{
						if(node2_code >= 0){
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = -1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = -1;
							ka++;
						}
						//else add nothing
					}
					dc_mna.b[dc_mna.m1 + m2_counter] += curr->val;
					dc_mna.b_complex[dc_mna.m1 + m2_counter] += curr->val_complex;
					m2_counter++;
					break;
				case 'i':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							dc_mna.b[node1_code] -= curr->val;
							dc_mna.b[node2_code] += curr->val;
							dc_mna.b_complex[node1_code] -= curr->val_complex;
							dc_mna.b_complex[node2_code] += curr->val_complex;
						} else {
							dc_mna.b[node1_code] -= curr->val;
							dc_mna.b_complex[node1_code] -= curr->val_complex;
						}
					} else {
						if(node2_code >=0) {
							dc_mna.b[node2_code] += curr->val;
							dc_mna.b_complex[node2_code] += curr->val_complex;
						}
					}
					break;
				case 'r':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1/curr->val;
							ka++;
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = -1/curr->val;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = -1/curr->val;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = 1/curr->val;
							ka++;
								
						} else {
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1/curr->val;
							ka++;
							
						}
					} else {
						if(node2_code >=0) {
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = 1/curr->val;
							ka++;
						}
					}
					break;
				case 'c':
					if(node1_code >= 0) {
						if(node2_code >= 0) {
							
							dc_mna.B_sp->i[kb] = node1_code;
							dc_mna.B_sp->p[kb] = node1_code;
							dc_mna.B_sp->x[kb] = curr->val;
							kb++;
							
							dc_mna.B_sp->i[kb] = node1_code;
							dc_mna.B_sp->p[kb] = node2_code;
							dc_mna.B_sp->x[kb] = -curr->val;
							kb++;
							
							dc_mna.B_sp->i[kb] = node2_code;
							dc_mna.B_sp->p[kb] = node1_code;
							dc_mna.B_sp->x[kb] = -curr->val;
							kb++;
							
							dc_mna.B_sp->i[kb] = node2_code;
							dc_mna.B_sp->p[kb] = node2_code;
							dc_mna.B_sp->x[kb] = curr->val;
							kb++;
								
						} else {
							dc_mna.B_sp->i[kb] = node1_code;
							dc_mna.B_sp->p[kb] = node1_code;
							dc_mna.B_sp->x[kb] = curr->val;
							kb++;
							
						}
					} else {
						if(node2_code >=0) {
							dc_mna.B_sp->i[kb] = node2_code;
							dc_mna.B_sp->p[kb] = node2_code;
							dc_mna.B_sp->x[kb] = curr->val;
							kb++;
						}
					}
					break;
				case 'l':
					if(node1_code >= 0){
						if(node2_code >= 0){
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = -1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = -1;
							ka++;		
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							
						} else {
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node1_code;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node1_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = 1;
							ka++;
						}
					} else{
						if(node2_code >= 0){
							dc_mna.A_sp->i[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->p[ka] = node2_code;
							dc_mna.A_sp->x[ka] = -1;
							ka++;
							
							dc_mna.A_sp->i[ka] = node2_code;
							dc_mna.A_sp->p[ka] = dc_mna.m1 + m2_counter;
							dc_mna.A_sp->x[ka] = -1;
							ka++;			
						}
						//else add nothing
					}
					dc_mna.B_sp->i[kb] = dc_mna.m1 + m2_counter;
					dc_mna.B_sp->p[kb] = dc_mna.m1 + m2_counter;
					dc_mna.B_sp->x[kb] = -curr->val;
					kb++;
					
					m2_counter++;
					break;
				default:
					break;
			}
			
			curr = curr->nxt;
		}
		
		dc_mna.C = cs_compress(dc_mna.A_sp);
		dc_mna.D = cs_compress(dc_mna.B_sp);
		cs_spfree(dc_mna.A_sp);
		cs_spfree(dc_mna.B_sp);
		if(!cs_dupl(dc_mna.C) || !cs_dupl(dc_mna.D)){
			printf("Error, cannot compress mna!\n");
		}
	}
	
	return dc_mna;
}

void print_mna_dc(dc_mna_t dc_mna){
	printf("\n-----MNA DC ANALYSIS----\n");
	printf("n-1 = %d\nm2 elements = %d\n", dc_mna.m1, dc_mna.m2);
	
	
	printf("----------A------------\n");
	
	for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
		for(int j=0; j < (dc_mna.m1+ dc_mna.m2); j++){
			printf("% 11.5f ", dc_mna.A[i][j]);
		}
		printf("\n");
	}
	
	printf("\n----------b------------\n");
	for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
		printf("% 11.5f ", dc_mna.b[i]);
	}
	printf("\n");
}

void delete_mna_dc(dc_mna_t dc_mna, bool is_sparse, bool is_iter){
	
	if(!is_sparse){
		for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
			free(dc_mna.A[i]);
		}
		free(dc_mna.A);	
	}
	else{
		cs_spfree(dc_mna.C);
	}
	
	free(dc_mna.b);
	free(dc_mna.b_complex);
	free(dc_mna.x0);
}


void sparse_lu_dec(dc_mna_t dc_mna){
	S=cs_sqr(2,dc_mna.C,0);
	if(S == NULL){
		printf("S is null\n!");
	}
	N=cs_lu(dc_mna.C,S,1);
	if(N == NULL){
		printf("N is null\n!");
	}
}

void sparse_lu_dec_complex(cs_cl* A){
	S_complex=cs_cl_sqr(2,A,0);
	if(S_complex == NULL){
		printf("S is null\n!");
	}
	N_complex=cs_cl_lu(A,S_complex,1);
	if(N_complex == NULL){
		printf("N is null\n!");
	}
}

void sparse_tran_lu_dec(cs* A_tran){
	S=cs_sqr(2,A_tran,0);
	if(S == NULL){
		printf("S is null\n!");
	}
	N=cs_lu(A_tran,S,1);
	if(N == NULL){
		printf("N is null\n!");
	}
}

void sparse_chol_dec(dc_mna_t dc_mna){
	S=cs_schol(1,dc_mna.C);
	if(S == NULL){
		printf("S is null\n!");
	}
	N=cs_chol(dc_mna.C,S);
	if(N == NULL){
		printf("N is null\n!");
	}
}

void sparse_tran_chol_dec(cs* A_tran){
	S=cs_schol(1,A_tran);
	if(S == NULL){
		printf("S is null\n!");
	}
	N=cs_lu(A_tran,S,1);
	if(N == NULL){
		printf("N is null\n!");
	}
}

void sparse_lu_solve(double* b, double* y, int n){
	cs_ipvec(N->pinv,b,y,n);
	cs_lsolve(N->L,y);
	cs_usolve(N->U,y);
	cs_ipvec(S->q,y,b,n);
}

void sparse_lu_solve_complex(double complex* b, double complex* y, int n){
	cs_cl_ipvec(N_complex->pinv,b,y,n);
	cs_cl_lsolve(N_complex->L,y);
	cs_cl_usolve(N_complex->U,y);
	cs_cl_ipvec(S_complex->q,y,b,n);
}

void sparse_chol_solve(double* b, double* y, int n){
	cs_ipvec(S->pinv,b,y,n);
	cs_lsolve(N->L,y);
	cs_ltsolve(N->L,y);
	cs_pvec(S->pinv,y,b,n);
}


void solve_direct_dcop(ht_t *ht, dc_mna_t dc_mna, bool is_spd, double** L, double** U, int* p, double* pb, double* x, double* y, int size, bool is_sparse){
	FILE* fp_dc = fopen("dc_op", "w+");
	htentry_t* node;
	if(!is_sparse){
		if(is_spd){
			printf("Solving with Cholesky decomposition!\n");
			if(!chol_dec(dc_mna.A, L, U, size)){
				printf("Warning!: The mna matrix is not SPD\n");
				lu_dec(dc_mna.A, L, U, p, size);		
			}
		}
		else{
			printf("Solving with LU decomposition!\n");
			lu_dec(dc_mna.A, L, U, p, size);			
		}
		
		// solving for dc op
		for(int i = 0; i < size; i++){
			pb[i] = dc_mna.b[p[i]];
		}
		solve_ld(L, pb, size, y);
		solve_ud(U, y, size, x);
		
		//copying initial solution to x0
		for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
			dc_mna.x0[i] = x[i];
		}
		
		for(int i = 0; i < ht->capacity; i++){
			node = ht->table[i];
			if(node)
				fprintf(fp_dc, "%s  %.12f\n", node->string, x[node->id]);
		}
		for(int i = dc_mna.m1; i < dc_mna.m1 + dc_mna.m2; i++){
			fprintf(fp_dc, "i%d  %.12f\n", i ,x[i]);
		}
	}
	else{
		double* b_sp = (double*)malloc(sizeof(double)*size);
		for(int i = 0; i < size; i++){
			b_sp[i] = dc_mna.b[i];
		}
		if(is_spd){
			printf("Solving with sparse Cholesky decomposition!\n");
			sparse_chol_dec(dc_mna);
			sparse_chol_solve(b_sp,y,size);
		}
		else{
			printf("Solving with sparse LU decomposition!\n");
			sparse_lu_dec(dc_mna);
			sparse_lu_solve(b_sp,y,size);
		}
		
		//copying initial solution to x0
		for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
			dc_mna.x0[i] = b_sp[i];
		}
		
		for(int i = 0; i < ht->capacity; i++){
			node = ht->table[i];
			if(node)
				fprintf(fp_dc, "%s  %.12f\n", node->string, b_sp[node->id]);
		}
		for(int i = dc_mna.m1; i < dc_mna.m1 + dc_mna.m2; i++){
			fprintf(fp_dc, "i%d  %.12f\n", i ,b_sp[i]);
		}
		free(b_sp);
	}
	fclose(fp_dc);
	
}

void solve_iter_dcop(ht_t* ht, dc_mna_t dc_mna, bool is_spd, double* x, double itol, int size, bool is_sparse){
	htentry_t* node;
	
	if(!is_sparse){
	
		if(is_spd){
			printf("Solving with Conjugate Gradients (itol=%.9f) algorithm!\n", itol);
			cg_calloc(size);
			cg(dc_mna.A,x,dc_mna.b,itol, size);
			free_cg();
		}
		else{
			printf("Solving with Bi-Conjugate Gradients (itol=%.9f) algorithm!\n", itol);
			bicg_calloc(size);
			bicg(dc_mna.A,x,dc_mna.b,itol, size);
			free_bicg();
		
		}
	}
	else{
		if(is_spd){
			printf("Solving with sparse Conjugate Gradients (itol=%.9f) algorithm!\n", itol);
			cg_calloc(size);
			sparse_cg(dc_mna.C,x,dc_mna.b,itol);
			free_cg();
		}
		else{
			printf("Solving with sparse Bi-Conjugate Gradients (itol=%.9f) algorithm!\n", itol);
			bicg_calloc(size);
			sparse_bicg(dc_mna.C,x,dc_mna.b,itol);	
            free_bicg();
		}
	}
	FILE* fp_dc = fopen("dc_op", "w+");
    
	//copying initial solution to x0
	for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
		dc_mna.x0[i] = x[i];
	}
	
	for(int i = 0; i < ht->capacity; i++){
		node = ht->table[i];
		if(node)
			fprintf(fp_dc, "%s  %.12f\n", node->string, x[node->id]);
	}
	for(int i = dc_mna.m1; i < dc_mna.m1 + dc_mna.m2; i++){
		fprintf(fp_dc, "i%d  %.12f\n", i ,x[i]);
	}
	fclose(fp_dc);
}

void plot_dcsweep(char* input_var, char* node_name, int dcsweep_counter){
	char filename[MAX_LEN_SP], temp_var[MAX_LEN], temp_name[MAX_LEN];
	int i;
	
	sprintf(filename, "dcsweep%d_%s_%s", dcsweep_counter, input_var, node_name);
	
	FILE *gnuplot = popen("gnuplot -persist", "w");
	
	for(i = 0; i < strlen(input_var); i++){
		if(input_var[i] == '_')
			temp_var[i] = '.';
		else
			temp_var[i] = input_var[i];	
	}
	temp_var[i] = '\0';
	
	for(i = 0; i < strlen(node_name); i++){
		if(node_name[i] == '_')
			temp_name[i] = '.';
		else
			temp_name[i] = node_name[i];	
	}
	temp_name[i] = '\0';
	
    if (!gnuplot) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    fprintf(gnuplot, "plot \"%s\" t 'dcsweep%d:%s:%s' w lp\n", filename, dcsweep_counter, temp_var, temp_name);
    fflush(gnuplot);

    pclose(gnuplot);
}

void plot_tran(char* node_name, int tran_counter){
	char filename[MAX_LEN_SP], temp_name[MAX_LEN];
	int i;
	
	sprintf(filename, "tran%d_%s", tran_counter, node_name);
	
	FILE *gnuplot = popen("gnuplot -persist", "w");
	
	for(i = 0; i < strlen(node_name); i++){
		if(node_name[i] == '_')
			temp_name[i] = '.';
		else
			temp_name[i] = node_name[i];	
	}
	temp_name[i] = '\0';
	
    if (!gnuplot) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    fprintf(gnuplot, "plot \"%s\" t 'tran%d:%s' w lp\n", filename, tran_counter, temp_name);
    fflush(gnuplot);

    pclose(gnuplot);
}



void plot_ac(char* node_name, int ac_counter, char* ac_type){
	char filename[MAX_LEN_SP], temp_name[MAX_LEN];
	int i;
	
	sprintf(filename, "ac%d_%s", ac_counter, node_name);
	
	FILE *gnuplot = popen("gnuplot -persist", "w");
	
	for(i = 0; i < strlen(node_name); i++){
		if(node_name[i] == '_')
			temp_name[i] = '.';
		else
			temp_name[i] = node_name[i];	
	}
	temp_name[i] = '\0';
	
    if (!gnuplot) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    if(!strcmp("log", ac_type)){
    	fprintf(gnuplot, "set logscale x\n");
	}
	
    fprintf(gnuplot, "plot \"%s\" using 1:2 t 'tran%d:%s.mag' w lp\n" , filename, ac_counter, temp_name);
    fflush(gnuplot);
    pclose(gnuplot);
    
    gnuplot = popen("gnuplot -persist", "w");
    if (!gnuplot) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    if(!strcmp("log", ac_type)){
    	fprintf(gnuplot, "set logscale x\n");
	}
    fprintf(gnuplot, "plot \"%s\" using 1:3 t 'tran%d:%s.phase' w lp\n",filename, ac_counter, temp_name );
    fflush(gnuplot);
    pclose(gnuplot);

    
}



void solve_direct_dcsweep(dc_mna_t dc_mna, ht_t* ht, list_t list, FILE* out_fp[], int num_prints, int nodes_print[], double** L, double** U, int* p, double* pb, double* bdc, double* x, double* y, char* input_var, double start_val, double end_val, double inc, bool is_sparse, bool is_spd){
	int pos_node, neg_node;
	double *b_sp = (double*)malloc(sizeof(double)*(get_ht_size(ht) + dc_mna.m2));
	node_t* el_ptr;
	
	for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
		bdc[i] = dc_mna.b[i];
	}
					
	// current source in dc sweep
	if(input_var[0] == 'i'){
		el_ptr = return_element(list, &input_var[1], 'i');
		pos_node = lookup(ht, el_ptr->node1);
		neg_node = lookup(ht, el_ptr->node2);
						
		// starting values for b
		if(pos_node >= 0)
			bdc[pos_node] = bdc[pos_node] + el_ptr->val - start_val;
						
		if(neg_node >= 0)
			bdc[neg_node] = bdc[neg_node] - el_ptr->val + start_val; 
						
		for(double val = start_val; val <= end_val; val += inc){
			if(!is_sparse){
				//init Pb
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					pb[i] = bdc[p[i]];
				}
				//print_vector(bdc, get_ht_size(ht) + m2);
				//solve Upper and Lower Triangular
				solve_ld(L, pb, get_ht_size(ht) + dc_mna.m2, y);
				solve_ud(U, y, get_ht_size(ht) + dc_mna.m2, x);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
			}
			else{
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					b_sp[i] = bdc[i];
				}
				if(!is_spd)
					sparse_lu_solve(b_sp, y, get_ht_size(ht) + dc_mna.m2);
				else
					sparse_chol_solve( b_sp, y, get_ht_size(ht) + dc_mna.m2);
					
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, b_sp[nodes_print[i]]);
				}
			}
			//print_vector(x, get_ht_size(ht) + m2);
			//increment for bdc
			if(pos_node >= 0)
				bdc[pos_node] -= inc;
								
			if(neg_node >= 0)	
				bdc[neg_node] += inc;
		}
	}
	else{ // voltage source	
		el_ptr = return_element(list, &input_var[1], 'v');
		pos_node = k_value(list, el_ptr->name) + get_ht_size(ht);
		//printf("k_value = %d\n", k_value(list, el_ptr->name) );
		bdc[pos_node] = bdc[pos_node] - el_ptr->val + start_val;
						
		for(double val = start_val; val <= end_val; val += inc){
			if(!is_sparse){
				//init Pb
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					pb[i] = bdc[p[i]];
				}
				//print_vector(bdc, get_ht_size(ht) + m2);
				//solve Upper and Lower Triangular
				solve_ld(L, pb, get_ht_size(ht) + dc_mna.m2, y);
				solve_ud(U, y, get_ht_size(ht) + dc_mna.m2, x);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
			}
			else{
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					b_sp[i] = bdc[i];
				}
				if(!is_spd)
					sparse_lu_solve(b_sp, y, get_ht_size(ht) + dc_mna.m2);
				else
					sparse_chol_solve(b_sp, y, get_ht_size(ht) + dc_mna.m2);
					
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, b_sp[nodes_print[i]]);
				}
			}
			//print_vector(x, get_ht_size(ht) + m2);
			//increment for bdc
			bdc[pos_node] += inc;
		}
						
	}
	if(is_sparse)
		free(b_sp);
	//closing out files
	for(int i = 0; i < num_prints; i++){
		fclose(out_fp[i]);
	}	
}


void solve_iter_dcsweep(dc_mna_t dc_mna, ht_t* ht, list_t list, FILE* out_fp[], int num_prints, int nodes_print[], double itol, double* bdc, double* x, char* input_var, double start_val, double end_val, double inc, bool is_sparse){
	int pos_node, neg_node;
	node_t* el_ptr;
	
	for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
		bdc[i] = dc_mna.b[i];
		x[i] = 0;
	}
					
	// current source in dc sweep
	if(input_var[0] == 'i'){
		el_ptr = return_element(list, &input_var[1], 'i');
		pos_node = lookup(ht, el_ptr->node1);
		neg_node = lookup(ht, el_ptr->node2);
						
		// starting values for b
		if(pos_node >= 0)
			bdc[pos_node] = bdc[pos_node] + el_ptr->val - start_val;
						
		if(neg_node >= 0)
			bdc[neg_node] = bdc[neg_node] - el_ptr->val + start_val; 
		if(!is_sparse){
		
		    bicg_calloc(get_ht_size(ht) + dc_mna.m2);
		
			for(double val = start_val; val <= end_val; val += inc){
				
				bicg(dc_mna.A,x,bdc,itol, get_ht_size(ht) + dc_mna.m2);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				if(pos_node >= 0)
					bdc[pos_node] -= inc;
									
				if(neg_node >= 0)	
					bdc[neg_node] += inc;
			}
			
			free_bicg();
		}
		else{
		
		    bicg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				
				sparse_bicg(dc_mna.C,x,bdc,itol);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				if(pos_node >= 0)
					bdc[pos_node] -= inc;
									
				if(neg_node >= 0)	
					bdc[neg_node] += inc;
			}
			
			free_bicg();
		}
	}
	else{ // voltage source	
		el_ptr = return_element(list, &input_var[1], 'v');
		pos_node = k_value(list, el_ptr->name) + get_ht_size(ht);
		//printf("k_value = %d\n", k_value(list, el_ptr->name) );
		bdc[pos_node] = bdc[pos_node] - el_ptr->val + start_val;
		
		if(!is_sparse){	
		
		    bicg_calloc(get_ht_size(ht) + dc_mna.m2);
		
			for(double val = start_val; val <= end_val; val += inc){
				
				bicg(dc_mna.A,x,bdc,itol, get_ht_size(ht) + dc_mna.m2);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				bdc[pos_node] += inc;
			}
			
			free_bicg();
		}
		else{
		
		    bicg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				
				sparse_bicg(dc_mna.C,x,bdc,itol);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				bdc[pos_node] += inc;
			}
			
			free_bicg();
		}
						
	}
					
	//closing out files
	for(int i = 0; i < num_prints; i++){
		fclose(out_fp[i]);
	}
}

void solve_iter_dcsweep_spd(dc_mna_t dc_mna, ht_t* ht, list_t list, FILE* out_fp[], int num_prints, int nodes_print[], double itol, double* bdc, double* x, char* input_var, double start_val, double end_val, double inc, bool is_sparse){
	int pos_node, neg_node;
	node_t* el_ptr;
	
	for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
		bdc[i] = dc_mna.b[i];
		x[i] = 0;
	}
					
	// current source in dc sweep
	if(input_var[0] == 'i'){
		el_ptr = return_element(list, &input_var[1], 'i');
		pos_node = lookup(ht, el_ptr->node1);
		neg_node = lookup(ht, el_ptr->node2);
						
		// starting values for b
		if(pos_node >= 0)
			bdc[pos_node] = bdc[pos_node] + el_ptr->val - start_val;
						
		if(neg_node >= 0)
			bdc[neg_node] = bdc[neg_node] - el_ptr->val + start_val; 
		
		if(!is_sparse){	
		    cg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				
				cg(dc_mna.A,x,bdc,itol, get_ht_size(ht) + dc_mna.m2);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				if(pos_node >= 0)
					bdc[pos_node] -= inc;
									
				if(neg_node >= 0)	
					bdc[neg_node] += inc;
			}
			
			free_cg();
		}
		else{
		
		    cg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				sparse_cg(dc_mna.C,x,bdc,itol);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				if(pos_node >= 0)
					bdc[pos_node] -= inc;
									
				if(neg_node >= 0)	
					bdc[neg_node] += inc;
			}
			
			free_cg();
		}
	}
	else{ // voltage source	
		el_ptr = return_element(list, &input_var[1], 'v');
		pos_node = k_value(list, el_ptr->name) + get_ht_size(ht);
		//printf("k_value = %d\n", k_value(list, el_ptr->name) );
		bdc[pos_node] = bdc[pos_node] - el_ptr->val + start_val;
		
		if(!is_sparse){	
		
		    cg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				
				cg(dc_mna.A,x,bdc,itol, get_ht_size(ht) + dc_mna.m2);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				bdc[pos_node] += inc;
			}
			
			free_cg();
		}
		else{
		
		    cg_calloc(get_ht_size(ht) + dc_mna.m2);
		        
			for(double val = start_val; val <= end_val; val += inc){
				sparse_cg(dc_mna.C,x,bdc,itol);
								
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", val, x[nodes_print[i]]);
				}
								
				//print_vector(x, get_ht_size(ht) + m2);
				//increment for bdc
				bdc[pos_node] += inc;
			}
			
			free_cg();
		}					
	}				
	//closing out files
	for(int i = 0; i < num_prints; i++){
		fclose(out_fp[i]);
	}
}


void solve_direct_tran(dc_mna_t dc_mna, ht_t* ht, list_t list, tran_list_t* tran_list, FILE* out_fp[], int num_prints, int nodes_print[], double end_val, double inc, bool is_sparse, bool is_spd, char method){
	double *x_cur = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double *x_prev = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double *b_tran = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double t_cur;
	double* e_cur = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double* e_prev = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double* y = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	int pos_node, neg_node;
	tran_node_t* curr;
	
	//initializes x_k-1 to x0 and e_prev to e_cur to b
	for(int i =0; i < dc_mna.m1 + dc_mna.m2; i++){
		x_prev[i] = dc_mna.x0[i];
		e_cur[i] = dc_mna.b[i];
	}
	
	if(is_sparse){
		cs* A_tran;
		cs* B_tran;
		if(method =='t'){
			A_tran = cs_add(dc_mna.C,dc_mna.D,1,2/inc);
			B_tran = cs_add(dc_mna.C,dc_mna.D,-1,2/inc);
			if(is_spd){
				sparse_tran_chol_dec(A_tran);
				
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_chol_solve(b_tran,y,dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, b_tran[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = b_tran[i];
					}
				}		
			}
			else{
				sparse_tran_lu_dec(A_tran);
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					//making e_tk
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_lu_solve(b_tran,y,dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, b_tran[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = b_tran[i];
					}
				}
			}
		}
		else{
		    A_tran = cs_add(dc_mna.C,dc_mna.D,1,1/inc);
		    B_tran = cs_add(dc_mna.C,dc_mna.D,0,1/inc);
			if(is_spd){
				sparse_tran_chol_dec(A_tran);
				
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_chol_solve(b_tran,y,dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, b_tran[nodes_print[i]]);
					}
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = b_tran[i];
					}
				}		
			}
			else{
				sparse_tran_lu_dec(A_tran);
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_lu_solve(b_tran,y,dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, b_tran[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = b_tran[i];
					}
				}
			}
		}
		cs_spfree(B_tran);
	}
	else{
		double** A_tran;
		double** B_tran;
		double** L_tran = (double**)malloc(sizeof(double*)*(dc_mna.m1+dc_mna.m2));
		double** U_tran = (double**)malloc(sizeof(double*)*(dc_mna.m1+dc_mna.m2));
		int* p_tran = (int*)malloc(sizeof(int)*(dc_mna.m1+dc_mna.m2));
		double* pb = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
		for(int i = 0; i < (dc_mna.m1+dc_mna.m2); i++){
			L_tran[i] = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
			U_tran[i] = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
		}
		
		if(method == 't'){
			A_tran = add(dc_mna.A, dc_mna.B, 1, 2/inc, (dc_mna.m1+dc_mna.m2));
			B_tran = add(dc_mna.A, dc_mna.B, -1, 2/inc, (dc_mna.m1+dc_mna.m2));
			if(is_spd){
				chol_dec(A_tran, L_tran, U_tran, (dc_mna.m1+dc_mna.m2));
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					//pb = b_tran for chol
					solve_ld(L_tran, b_tran, get_ht_size(ht) + dc_mna.m2, y);
					solve_ud(U_tran, y, get_ht_size(ht) + dc_mna.m2, x_cur);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
			}
			else{
				lu_dec(A_tran,L_tran,U_tran,p_tran,(dc_mna.m1+dc_mna.m2));
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					// line permutation for right hand side
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						pb[i] = b_tran[p_tran[i]];
					}
					solve_ld(L_tran, pb, get_ht_size(ht) + dc_mna.m2, y);
					solve_ud(U_tran, y, get_ht_size(ht) + dc_mna.m2, x_cur);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
			}
		}
		else{
			A_tran = add(dc_mna.A, dc_mna.B, 1, 1/inc, (dc_mna.m1+dc_mna.m2));
			B_tran = add(dc_mna.A, dc_mna.B, 0, 1/inc, (dc_mna.m1+dc_mna.m2));
			if(is_spd){
				chol_dec(A_tran, L_tran, U_tran, (dc_mna.m1+dc_mna.m2));
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					//pb = b_tran for chol
					solve_ld(L_tran, b_tran, get_ht_size(ht) + dc_mna.m2, y);
					solve_ud(U_tran, y, get_ht_size(ht) + dc_mna.m2, x_cur);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
			}
			else{
				lu_dec(A_tran,L_tran,U_tran,p_tran,(dc_mna.m1+dc_mna.m2));
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					// line permutation for right hand side
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						pb[i] = b_tran[p_tran[i]];
					}
					solve_ld(L_tran, pb, get_ht_size(ht) + dc_mna.m2, y);
					solve_ud(U_tran, y, get_ht_size(ht) + dc_mna.m2, x_cur);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
			}	
		}
		
		for(int i = 0; i <(dc_mna.m1+dc_mna.m2); i++){
			free(A_tran[i]);
			free(B_tran[i]);
			free(L_tran[i]);
			free(U_tran[i]);
		}
		free(A_tran);
		free(B_tran);
		free(L_tran);
		free(U_tran);
		free(p_tran);
	}
	
	free(y);
	free(x_cur);
	free(x_prev);
	free(b_tran);
	free(e_cur);
	free(e_prev);
	
	//closing out files
	for(int i = 0; i < num_prints; i++){
		fclose(out_fp[i]);
	}
}


void solve_iter_tran(dc_mna_t dc_mna, ht_t* ht, list_t list, tran_list_t* tran_list, FILE* out_fp[], int num_prints, int nodes_print[], double end_val, double inc, bool is_sparse, bool is_spd, char method, int itol){
	double *x_cur = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double *x_prev = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double *b_tran = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double t_cur;
	double* e_cur = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	double* e_prev = (double*)malloc(sizeof(double)*(dc_mna.m1+dc_mna.m2));
	int pos_node, neg_node;
	tran_node_t* curr;
	
	//initializes x_k-1 to x0 and e_prev to e_cur to b
	for(int i =0; i < dc_mna.m1 + dc_mna.m2; i++){
		x_prev[i] = dc_mna.x0[i];
		x_cur[i] = 0;
		e_cur[i] = dc_mna.b[i];
	}
	
	if(is_sparse){
		cs* A_tran;
		cs* B_tran;
		if(method =='t'){
			A_tran = cs_add(dc_mna.C,dc_mna.D,1,2/inc);
			B_tran = cs_add(dc_mna.C,dc_mna.D,-1,2/inc);
			if(is_spd){	
			        
			    cg_calloc(dc_mna.m1 + dc_mna.m2);
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_cg(A_tran, x_cur, b_tran, itol);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_cg();
			}
			else{
			        
			    bicg_calloc(dc_mna.m1 + dc_mna.m2);
			        
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					//making e_tk
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_bicg(A_tran, x_cur, b_tran, itol);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_bicg();
			}
		}
		else{
		    A_tran = cs_add(dc_mna.C,dc_mna.D,1,1/inc);
		    B_tran = cs_add(dc_mna.C,dc_mna.D,0,1/inc);
			if(is_spd){
			       
			    cg_calloc(dc_mna.m1 + dc_mna.m2);
			        
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_cg(A_tran, x_cur, b_tran, itol);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_cg();
			}
			else{
			        
			    bicg_calloc(dc_mna.m1 + dc_mna.m2);
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					if(!cs_gaxpy(B_tran,x_prev, b_tran)){
						printf("Error using gaxpy\n");
						exit(-1);
					}
					
					sparse_bicg(A_tran, x_cur, b_tran, itol);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, b_tran[nodes_print[i]]);
					}
				}
				
				free_bicg();
			}
		}
		cs_spfree(B_tran);
	}
	else{
		double** A_tran;
		double** B_tran;
		
		if(method == 't'){
			A_tran = add(dc_mna.A, dc_mna.B, 1, 2/inc, (dc_mna.m1+dc_mna.m2));
			B_tran = add(dc_mna.A, dc_mna.B, -1, 2/inc, (dc_mna.m1+dc_mna.m2));
			if(is_spd){
			
			    cg_calloc(dc_mna.m1+dc_mna.m2);
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					cg(A_tran, x_cur, b_tran, itol, dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_cg();
				
			}
			else{
			
			    bicg_calloc(dc_mna.m1 + dc_mna.m2);			
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_prev[i] = e_cur[i];
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i] + e_prev[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					bicg(A_tran, x_cur, b_tran, itol, dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_bicg();
			}
		}
		else{
			A_tran = add(dc_mna.A, dc_mna.B, 1, 1/inc, (dc_mna.m1+dc_mna.m2));
			B_tran = add(dc_mna.A, dc_mna.B, 0, 1/inc, (dc_mna.m1+dc_mna.m2));
			if(is_spd){
			
			    cg_calloc(dc_mna.m1+dc_mna.m2);
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					cg(A_tran, x_cur, b_tran, itol, dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_cg();
				
			}
			else{
			
			    bicg_calloc(dc_mna.m1 + dc_mna.m2);
			
				for(t_cur = inc; t_cur <= end_val; t_cur += inc){
					//making new b_tran
					for(int i = 0; i<dc_mna.m1 + dc_mna.m2; i++){
						e_cur[i] = dc_mna.b[i]; 
					}
					curr = tran_list->head;
					while(curr != NULL){
						//if curr tran element is i
						if(curr->tran_node->el == 'i'){
							pos_node = lookup(ht, curr->tran_node->node1);
							neg_node = lookup(ht, curr->tran_node->node2);
							if(pos_node >= 0){
								e_cur[pos_node] = e_cur[pos_node] + curr->tran_node->val - (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
							if(neg_node >= 0){
								e_cur[neg_node] = e_cur[neg_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
							}
								
						}
						//if curr tran element is v
						else{
							pos_node = k_value(list, curr->tran_node->name) + get_ht_size(ht);
							e_cur[pos_node] = e_cur[pos_node] - curr->tran_node->val + (*curr->tran_node->tran_spec->fun)(curr->tran_node->tran_spec->args, curr->tran_node->tran_spec->num_args, t_cur);
						}
						
						curr = curr->next;
					}
					// we will store e_cur and e_prev to b_tran
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						b_tran[i] = e_cur[i];
					}
					
					gaxpy(B_tran,x_prev, b_tran, dc_mna.m1 + dc_mna.m2 );
					
					bicg(A_tran, x_cur, b_tran, itol, dc_mna.m1 + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%16.9f \t\t %16.9f\n", t_cur, x_cur[nodes_print[i]]);
					}
					
					//updating x_prev
					for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
						x_prev[i] = x_cur[i];
					}
				}
				
				free_bicg();
			}	
		}	
	}
	
	free(x_cur);
	free(x_prev);
	free(b_tran);
	free(e_cur);
	free(e_prev);
	
	//closing out files
	for(int i = 0; i < num_prints; i++){
		fclose(out_fp[i]);
	}
}


double complex** copy_dc_mna(dc_mna_t dc_mna){
	//printf("Copying dc mna to AC matrix\n");
	double complex** A_complex = (double complex**)malloc(sizeof(double complex*)*(dc_mna.m1 + dc_mna.m2));
	for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
		A_complex[i] = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		for(int j = 0; j < dc_mna.m1 + dc_mna.m2; j++){
			A_complex[i][j] = dc_mna.A[i][j] + 0*I;
		}
	}
	
	return A_complex;
}


cs_cl* copy_dc_mna_sparse(dc_mna_t dc_mna){
	int nz = dc_mna.A_sp->nz;
	
	cs_cl* A = cs_cl_spalloc(dc_mna.m1+dc_mna.m2,dc_mna.m1+dc_mna.m2,nz,1,1);
	A->nz=nz;
	
	for(int k = 0; k < nz; k++){
		A->i[k] = dc_mna.A_sp->i[k];
		A->x[k] = dc_mna.A_sp->x[k];
		A->p[k] = dc_mna.A_sp->p[k];
	}
	
	return A;
}


void free_ac_matrices(dc_mna_t dc_mna, double complex **A_complex1, double complex **A_complex2){
	//printf("Freeing AC matrix\n");
	for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
		free(A_complex1[i]);
		free(A_complex2[i]);
	}
	
	free(A_complex1);
	free(A_complex2);
}


double complex** create_ac_matrix(dc_mna_t dc_mna, ht_t* ht, ac_list_t ac_list, double complex** A_complex_init, double freq){
	ac_node_t *cur_node = ac_list.head;
	double complex** A_complex_cur = (double complex**)malloc(sizeof(double complex*)*(dc_mna.m1 + dc_mna.m2));
	
	//printf("Creating matrix for freq %f\n", freq);
	for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
		A_complex_cur[i] = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		for(int j = 0; j < dc_mna.m1 + dc_mna.m2; j++){
			A_complex_cur[i][j] = A_complex_init[i][j];
		}
	}
	
	while(cur_node != NULL){
		int node1_code = lookup(ht, cur_node->ac_node->node1);
		int node2_code = lookup(ht, cur_node->ac_node->node2);
		if(cur_node->ac_node->el == 'c'){
			//printf("Node C with val %f\n", cur_node->ac_node->val);
			if(node1_code >= 0){
				if(node2_code >= 0){
					A_complex_cur[node1_code][node1_code] += I*2*PI*freq*cur_node->ac_node->val;
					A_complex_cur[node1_code][node2_code] -= I*2*PI*freq*cur_node->ac_node->val;
					A_complex_cur[node2_code][node1_code] -= I*2*PI*freq*cur_node->ac_node->val;
					A_complex_cur[node2_code][node2_code] += I*2*PI*freq*cur_node->ac_node->val;
				}
				else{
					A_complex_cur[node1_code][node1_code] += I*2*PI*freq*cur_node->ac_node->val;
				}
			}
			else{
				if(node2_code >= 0){
					A_complex_cur[node2_code][node2_code] += I*2*PI*freq*cur_node->ac_node->val;
				}
				// else add nothing
			}
		}
		else if(cur_node->ac_node->el == 'l'){
			//printf("Node L with val %f\n", cur_node->ac_node->val);
			A_complex_cur[dc_mna.m1 + cur_node->ac_node->m2_pos][dc_mna.m1 + cur_node->ac_node->m2_pos] -= I*2*PI*freq*cur_node->ac_node->val;
		}
		
		cur_node = cur_node->next;
	}
	
	return A_complex_cur;
}


cs_cl* create_ac_matrix_sparse(dc_mna_t dc_mna, ht_t* ht, list_t list, ac_list_t ac_list, double freq){
	node_t *cur_node = list.head;
	int nz = get_nz_complex(list), m2_counter = 0, k = 0;
	cs_cl *A_complex_cur;
	A_complex_cur=cs_cl_spalloc(dc_mna.m1+dc_mna.m2,dc_mna.m1+dc_mna.m2,nz,1,1);
	A_complex_cur->nz=nz;
	
	while(cur_node != NULL){
		int node1_code = lookup(ht, cur_node->node1);
		int node2_code = lookup(ht, cur_node->node2);
		if(cur_node->el == 'c'){
			if(node1_code >= 0){
				if(node2_code >= 0){
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = I*2*PI*freq*cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -I*2*PI*freq*cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = -I*2*PI*freq*cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = I*2*PI*freq*cur_node->val;
					k++;
				}
				else{
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = I*2*PI*freq*cur_node->val;
					k++;
				}
			}
			else{
				if(node2_code >= 0){
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = I*2*PI*freq*cur_node->val;
					k++;
				}
				// else add nothing
			}
		}
		else if(cur_node->el == 'l'){
			if(node1_code >= 0){
				if(node2_code >= 0){
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1;
					k++;
							
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -1;
					k++;
							
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = -1;
					k++;		
							
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = 1;
					k++;
				}
				else{
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1;
					k++;
					
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = 1;
					k++;
				}
			}
			else{
				if(node2_code >= 0){
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -1;
					k++;
							
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = -1;
					k++;
				}
				// else add nothing
			}
			A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
			A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
			A_complex_cur->x[k] = -2*PI*freq*cur_node->val;
			k++;
			m2_counter++;
		}
		else if(cur_node->el == 'r'){
			if(node1_code >= 0){
				if(node2_code >= 0){
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1/cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -1/cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = -1/cur_node->val;
					k++;
					
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = 1/cur_node->val;
					k++;
				}
				else{
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1/cur_node->val;
					k++;
				}
			}
			else{
				if(node2_code >= 0){
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = 1/cur_node->val;
					k++;
				}
				// else add nothing
			}
		}
		else if(cur_node->el == 'v'){
			if(node1_code >= 0){
				if(node2_code >= 0){
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1;
					k++;
							
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -1;
					k++;
							
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = -1;
					k++;		
							
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = 1;
					k++;
				}
				else{
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node1_code;
					A_complex_cur->x[k] = 1;
					k++;
					
					A_complex_cur->i[k] = node1_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = 1;
					k++;
				}
			}
			else{
				if(node2_code >= 0){
					A_complex_cur->i[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->p[k] = node2_code;
					A_complex_cur->x[k] = -1;
					k++;
							
					A_complex_cur->i[k] = node2_code;
					A_complex_cur->p[k] = dc_mna.m1 + m2_counter;
					A_complex_cur->x[k] = -1;
					k++;
				}
				// else add nothing
			}
			m2_counter++;
		}
		
		cur_node = cur_node->nxt;
	}
	
	
	cs_cl *A_com = cs_cl_compress(A_complex_cur);
	cs_cl_spfree(A_complex_cur);
	cs_cl_dupl(A_com);
	
	return A_com;
}


int get_dec_pow(double x){
	int i;
	
	for(i = -10; i <= 10; i++){
		if(x <= pow(10, i)- pow(10,-10))
			break;
	}
	
	return i - 1;
}


void solve_direct_ac(dc_mna_t dc_mna, ht_t* ht, list_t list, ac_list_t ac_list, FILE* out_fp[], int num_prints, int nodes_print[], char* ac_type, int num_points, double start_val, double end_val, bool is_sparse ){
	double inc;
	
	if(!is_sparse){
		double complex **A_complex_init = copy_dc_mna(dc_mna);
		double complex **A_complex_cur;
		double complex **L = (double complex**)malloc(sizeof(double complex*)*(dc_mna.m1 + dc_mna.m2));
		double complex **U = (double complex**)malloc(sizeof(double complex*)*(dc_mna.m1 + dc_mna.m2));
		int *p = (int*)malloc(sizeof(int)*(dc_mna.m1 + dc_mna.m2));
		double complex *pb = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		double complex *x = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		double complex *y= (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		for(int i=0; i<dc_mna.m1 + dc_mna.m2; i++){
			L[i] = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
			U[i] = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		}
		printf("Solve Direct AC\n");
		// if sweep is linear
		if(!strcmp("lin", ac_type)){
			inc = (end_val - start_val)/num_points;
			for(double cur_val = start_val; cur_val <= end_val; cur_val += inc){
				A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
				lu_dec_complex(A_complex_cur, L, U, p, dc_mna.m1+dc_mna.m2);
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					pb[i] = dc_mna.b_complex[p[i]];
				}
				//solve Upper and Lower Triangular
				solve_ld_complex(L, pb, get_ht_size(ht) + dc_mna.m2, y);
				solve_ud_complex(U, y, get_ht_size(ht) + dc_mna.m2, x);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f \t\t %16.9f\n", cur_val, cabs(x[nodes_print[i]]), 180*carg(x[nodes_print[i]])/PI);
				}
			}	
		}
		// if sweep is log
		else if(!strcmp("log", ac_type)){
			double log_inc = 10.0/num_points;
			double cur_val = start_val;
			double cur_dec_pow = get_dec_pow(cur_val);
			while(cur_val <= end_val){
				A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
				lu_dec_complex(A_complex_cur, L, U, p, dc_mna.m1+dc_mna.m2);
				for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
					pb[i] = dc_mna.b_complex[p[i]];
				}
				//solve Upper and Lower Triangular
				solve_ld_complex(L, pb, get_ht_size(ht) + dc_mna.m2, y);
				solve_ud_complex(U, y, get_ht_size(ht) + dc_mna.m2, x);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
				}
				
				cur_val += log_inc*pow(10, cur_dec_pow);
				
				//solve for end_val
				/*if(cur_val > end_val){
					cur_val = end_val;
					A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
					lu_dec_complex(A_complex_cur, L, U, p, dc_mna.m1+dc_mna.m2);
					for(int i = 0; i < get_ht_size(ht) + dc_mna.m2; i++){
						pb[i] = dc_mna.b_complex[p[i]];
					}
					//solve Upper and Lower Triangular
					solve_ld_complex(L, pb, get_ht_size(ht) + dc_mna.m2, y);
					solve_ud_complex(U, y, get_ht_size(ht) + dc_mna.m2, x);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
					}
					break;
				}*/
				cur_dec_pow = get_dec_pow(cur_val);
				
			}
		}
		free_ac_matrices(dc_mna, A_complex_init, A_complex_cur);
		free(pb);
		free(p);
		free(y);
		free(x);
		for(int i=0; i<dc_mna.m1 + dc_mna.m2; i++){
			free(L[i]);
			free(U[i]);
		}
		free(L);
		free(U);
	}
	else{
		double complex* b_sp = (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		double complex* y= (double complex*)malloc(sizeof(double complex)*(dc_mna.m1 + dc_mna.m2));
		cs_cl *A_complex_cur;
		printf("Solve Direct sparse AC\n");
		
		if(!strcmp("lin", ac_type)){
			inc = (end_val - start_val)/num_points;
			for(double cur_val = start_val; cur_val <= end_val; cur_val += inc){
				
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					b_sp[i] = dc_mna.b_complex[i];
				}
				
				A_complex_cur = create_ac_matrix_sparse(dc_mna, ht, list, ac_list, cur_val);
				sparse_lu_dec_complex(A_complex_cur);
				sparse_lu_solve_complex(b_sp,y,dc_mna.m1 + dc_mna.m2);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, cabs(b_sp[nodes_print[i]]), 180*carg(b_sp[nodes_print[i]])/PI);
				}
				
			}
		}
		else if(!strcmp("log", ac_type)){
			double log_inc = 10.0/num_points;
			double cur_val = start_val;
			double cur_dec_pow = get_dec_pow(cur_val);
			while(cur_val <= end_val){
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					b_sp[i] = dc_mna.b_complex[i];
				}
				A_complex_cur = create_ac_matrix_sparse(dc_mna, ht, list, ac_list, cur_val);
				sparse_lu_dec_complex(A_complex_cur);
				sparse_lu_solve_complex(b_sp,y,dc_mna.m1 + dc_mna.m2);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(b_sp[nodes_print[i]])), 180*carg(b_sp[nodes_print[i]])/PI);
				}
				cur_val += log_inc*pow(10, cur_dec_pow);
				cur_dec_pow = get_dec_pow(cur_val);
			}
		}
		cs_cl_spfree(A_complex_cur);
		free(y);
		free(b_sp);
	}
	
	
}

void solve_iter_ac(dc_mna_t dc_mna, ht_t* ht, list_t list, ac_list_t ac_list, FILE* out_fp[], int num_prints, int nodes_print[], char* ac_type, int num_points, double start_val, double end_val, bool is_sparse, double itol){
	double inc;
	
	bicg_calloc_complex(get_ht_size(ht) + dc_mna.m2);
	if(!is_sparse){
		double complex **A_complex_init = copy_dc_mna(dc_mna);
		double complex **A_complex_cur;
		double complex *x = (double complex*)calloc(sizeof(double complex),(dc_mna.m1 + dc_mna.m2));
		
		printf("Solve Iter AC\n");
		// if sweep is linear
		if(!strcmp("lin", ac_type)){
			inc = (end_val - start_val)/num_points;
			
			for(double cur_val = start_val; cur_val <= end_val; cur_val += inc){
				A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					x[i] = 1;
				}
				bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol, get_ht_size(ht) + dc_mna.m2);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%16.9f \t\t %16.9f \t\t %16.9f\n", cur_val, cabs(x[nodes_print[i]]), 180*carg(x[nodes_print[i]])/PI);
				}
			}
				
		}
		// if sweep is log
		else if(!strcmp("log", ac_type)){
			double log_inc = 10.0/num_points;
			double cur_val = start_val;
			double cur_dec_pow = get_dec_pow(cur_val);
			while(cur_val <= end_val){
				A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					x[i] = 1;
				}
				bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol, get_ht_size(ht) + dc_mna.m2);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
				}
				
				cur_val += log_inc*pow(10, cur_dec_pow);
				
				//solve for end_val
				/*if(cur_val > end_val){
					cur_val = end_val;
					A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
					
					bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol, get_ht_size(ht) + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
					}
					break;
				}*/
				cur_dec_pow = get_dec_pow(cur_val);
				
			}
		}
		free_ac_matrices(dc_mna, A_complex_init, A_complex_cur);
		free(x);
	}
	else{
		cs_cl *A_complex_cur;
		double complex *x = (double complex*)calloc(sizeof(double complex),(dc_mna.m1 + dc_mna.m2));
		printf("Solve Sparse Iter AC\n");
		// if sweep is linear
		if(!strcmp("lin", ac_type)){
			inc = (end_val - start_val)/num_points;
			
			for(double cur_val = start_val; cur_val <= end_val; cur_val += inc){
				A_complex_cur = create_ac_matrix_sparse(dc_mna, ht, list, ac_list, cur_val);
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					x[i] = 1;
				}
				sparse_bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, cabs(x[nodes_print[i]]), 180*carg(x[nodes_print[i]])/PI);
				}
			}
				
		}
		// if sweep is log
		else if(!strcmp("log", ac_type)){
			double log_inc = 10.0/num_points;
			double cur_val = start_val;
			double cur_dec_pow = get_dec_pow(cur_val);
			while(cur_val <= end_val){
				A_complex_cur = create_ac_matrix_sparse(dc_mna, ht, list, ac_list, cur_val);
				for(int i = 0; i < dc_mna.m1 + dc_mna.m2; i++){
					x[i] = 1;
				}
				sparse_bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol);
				
				//print in file
				for(int i = 0; i < num_prints; i++){
					fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
				}
				
				cur_val += log_inc*pow(10, cur_dec_pow);
				
				//solve for end_val
				/*if(cur_val > end_val){
					cur_val = end_val;
					A_complex_cur = create_ac_matrix(dc_mna, ht, ac_list, A_complex_init, cur_val);
					
					bicg_complex(A_complex_cur, x, dc_mna.b_complex, itol, get_ht_size(ht) + dc_mna.m2);
					
					//print in file
					for(int i = 0; i < num_prints; i++){
						fprintf(out_fp[i], "%20.9f \t\t %16.9f \t\t %16.9f\n", cur_val, 20*log10(cabs(x[nodes_print[i]])), 180*carg(x[nodes_print[i]])/PI);
					}
					break;
				}*/
				cur_dec_pow = get_dec_pow(cur_val);
				
			}
		}
		cs_cl_spfree(A_complex_cur);
		free(x);	
	}
	free_bicg_complex();
}


int main(int argc, char *argv[]) {
	FILE * fp, *out_fp[MAX_LEN];
	char line[IMPORT_LINE_SIZE], element, name[MAX_LEN], node1[MAX_LEN], node2[MAX_LEN], w_str[MAX_LEN], l_str[MAX_LEN];
	char node3[MAX_LEN], node4[MAX_LEN], val_str[MAX_LEN], model_name[MAX_LEN], area_str[MAX_LEN], node_names[MAX_LEN][MAX_LEN];
	char *token = NULL, method ='t', analysis='n', *line_cpy, *ac_specs, *ac_type;
	char input_var[MAX_LEN], str1[MAX_LEN], str2[MAX_LEN];
	char node_name[MAX_LEN], out_name[MAX_LEN_SP];
	int i, j, m2 = 0, dcsweep_counter = 0, tran_counter = 0, ac_counter = 0, num_points, num_prints=0, nodes_print[MAX_LEN];
	int num_args;
	double val, area, w, l, start_val, end_val, inc, itol = 1e-3, mag, phase;
	ht_t* ht = init_ht(INITIAL_SIZE);	
	list_t list;
	tran_list_t tran_list;
	ac_list_t ac_list;
	dc_mna_t dc_mna;
	bool is_spd = false, is_iter=false, is_sparse=false;
	double *x, *y, *pb, *bdc, args[MAX_LEN];
	double **L, **U;
	int *p;
	tran_spec_t* tran_spec;
	
	list_init(&list, &tran_list, &ac_list);

	if(argc == 2){
		fp = fopen(argv[1], "r");
	}
	else{
		printf("Invalid number of arguments!\n");
		return(-1);
	}
	
	if(!fp){
		printf("Invalid argument\n");
		return(-2);
	}
	
	printf("Parsing netlist...\n");
	while(1) {
		num_args = area = w = l = 0;
		tran_spec = NULL;
		mag = phase = 0;
	
		if(fgets(line, IMPORT_LINE_SIZE, fp) == NULL) { //read file line by line
			break;
		}
		//printf("line is: %s", line);
		
		//convert string to lowercase characters
		i = 0;
		while(line[i] != '\0') {
			line[i] = tolower(line[i]);
			i++;
		}
		//printf("%s", line);
		
		//replace all tabs with whitespaces
		i = 0;
		while(line[i] != '\0') {
			if(line[i] == '\t') {
				line[i] = ' ';
			}
			i++;
		}
		
		if(line[0] == '*') { //ignore lines starting with '*'
			continue;
		}
		else if(line[0] == 'v') {
			line[strlen(line) - 1] = '\0';
			line_cpy = strdup(line);
			element = 'v';
			
			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
			m2++;
			
			//check for ac
			//printf("Rem line: %s\n", &line_cpy[strlen(name)+strlen(node1)+strlen(node2)+strlen(val_str)+4]);
			ac_specs = strchr(&line_cpy[strlen(name)+strlen(node1)+strlen(node2)+strlen(val_str)+4], 'c');
			
			if (ac_specs){
				ac_specs = &ac_specs[2];
				//node3 has AC mag and node4 has AC phase
				split_string(ac_specs, node3, node4, ' ');
				mag = atof(node3);
				phase = atof(node4);
				//deg to rad
				phase = (phase * PI)/180;
			}
			token = strtok(NULL, "(");
			if(token != NULL){
				//if it is not ac
				if(token[0] != 'a'){
					// if there is transient spec
					// we will use node3 for fun name
					//printf("Element v%s has transient %s\n", name, token);
					strcpy(node3, token);
					
					//if is not pwl
					if(strcmp(token,"pwl")){
						while( (token = strtok(NULL, " ") ) ){
							args[num_args] = atof(token);
							//printf("Token: %f\n", args[num_args]);
							num_args++;
						}
					}
					else{	
						token = strtok(NULL, ")");
						split_string(token, l_str, w_str, ' ');
						args[num_args] = atof(l_str);
						args[num_args + 1] = atof(w_str);
						//printf("Token1: %f, Token2: %f\n", args[num_args], args[num_args+1]);
						num_args +=2;
					
						while((token = strtok(NULL, ")"))){
							split_string(&token[1], l_str, w_str, ' ');
							args[num_args] = atof(l_str);
							args[num_args + 1] = atof(w_str);
							//printf("Token1: %f, Token2: %f\n", args[num_args], args[num_args+1]);
							num_args +=2;
						}
					}
					tran_spec = tran_spec_create(node3, args, num_args);
					//printf("Value: %f\n", (*tran_spec->fun)(tran_spec->args, tran_spec->num_args, 46.5));
				}
			}
			
			
		}
		else if(line[0] == 'r') {
			element = 'r';

			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
		}
		else if(line[0] == 'i') {
			element = 'i';
			line_cpy = strdup(line);
			
			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
			
			//check for ac
			ac_specs = strchr(&line_cpy[strlen(name)+strlen(node1)+strlen(node2)+strlen(val_str)+4], 'c');
			
			if (ac_specs){
				ac_specs = &ac_specs[2];;
				//node3 has AC mag and node4 has AC phase
				split_string(ac_specs, node3, node4, ' ');
				mag = atof(node3);
				phase = atof(node4);
				phase = (phase * PI)/180;
			}
			
			token = strtok(NULL, "(");
			if(token != NULL){
				//if is not ac 
				if (token[0] != 'a'){
					// if there is transient spec
					// we will use node3 for fun name
					//printf("Element v%s has transient %s\n", name, token);
					strcpy(node3, token);
					
					//if is not pwl
					if(strcmp(token,"pwl")){
						while( (token = strtok(NULL, " ") ) ){
							args[num_args] = atof(token);
							//printf("Token: %f\n", args[num_args]);
							num_args++;
						}
					}
					else{	
						token = strtok(NULL, ")");
						split_string(token, l_str, w_str, ' ');
						args[num_args] = atof(l_str);
						args[num_args + 1] = atof(w_str);
						//printf("Token1: %f, Token2: %f\n", args[num_args], args[num_args+1]);
						num_args +=2;
					
						while((token = strtok(NULL, ")"))){
							split_string(&token[1], l_str, w_str, ' ');
							args[num_args] = atof(l_str);
							args[num_args + 1] = atof(w_str);
							//printf("Token1: %f, Token2: %f\n", args[num_args], args[num_args+1]);
							num_args +=2;
						}
					}
					tran_spec = tran_spec_create(node3, args, num_args);
					//printf("Value: %f\n", (*tran_spec->fun)(tran_spec->args, tran_spec->num_args, 46.5));
				}
			}
		}
		else if(line[0] == 'c') {
			element = 'c';

			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
		}
		else if(line[0] == 'l') {
			element = 'l';

			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
			m2++;
		}
		else if(line[0] == 'd') {
			element = 'd';
			
			token = strtok(&line[1], " ");
			strcpy(name, token);	
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(model_name, token);
			
			token = strtok(NULL, " ");
			
			if(token == NULL){
				// if there is no area then model_name has an extra \n at the end
				model_name[strlen(model_name) - 1] = '\0';
				area = 1;
			}
			else{
				strcpy(area_str, token);
				area = atof(area_str);
			}
		}
		else if(line[0] == 'm') {
			element = 'm';
			
			token = strtok(&line[1], " ");
			strcpy(name, token);	
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(node3, token);
			token = strtok(NULL, " ");
			strcpy(node4, token);
			token = strtok(NULL, " ");
			strcpy(model_name, token);
			token = strtok(NULL, " ");
			strcpy(l_str, &token[2]);
			token = strtok(NULL, " ");
			strcpy(w_str, &token[2]);
			
			l = atof(l_str);
			w = atof(w_str);
		}
		else if(line[0] == 'q') {
			element = 'q';
			
			token = strtok(&line[1], " ");
			strcpy(name, token);	
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(node3, token);
			token = strtok(NULL, " ");
			strcpy(model_name, token);
			
			token = strtok(NULL, " ");
			
			if(token == NULL){
				// if there is no area then model_name has an extra \n at the end
				model_name[strlen(model_name) - 1] = '\0';
				area = 1;
			}
			else{
				strcpy(area_str, token);
				area = atof(area_str);
			}
		}
		else if(line[0] == '.'){
			line[strlen(line)-1] = '\0';
			token = strtok(&line[1], " ");
			if(!strcmp("options" , token)) {
				while( (token = strtok(NULL, " ") ) ){	
					//printf("Token:%s\n", token);
					if(!strcmp("spd" , token)){
						is_spd = true;
					}
					else if(!strcmp("iter" , token)){
						is_iter = true;
					}
					else if(!strcmp("sparse" , token)){
						is_sparse = true;
					}
					else{ //itol or method
						split_string(token, str1, str2, '=');
						
						//itol for iter method
						if(!strcmp("itol" , str1)){
							itol = atof(str2);
						}
						else{ //method for trans
							if(!strcmp("tr" , str2)){
								method = 't';
							}
							else if(!strcmp("be" ,str2)){
								method = 'b';
							}	
						}
					}
				}
			}
			else if(!strcmp("dc" , token)){ //dc sweep
				dcsweep_counter++;
				token = strtok(NULL, " ");
				strcpy(input_var, token);
				token = strtok(NULL, " ");
				start_val = atof(token);
				token = strtok(NULL, " ");
				end_val =  atof(token);
				token = strtok(NULL, "\n");
				inc = atof(token);
				analysis = 'd';
				num_prints = 0;
				printf("--------DC SWEEP: %d---------\n", dcsweep_counter);
			}
			else if(!strcmp("tran", token)){ //tran 
				tran_counter++;
				start_val = 0; //assume starting time = 0
				token = strtok(NULL, " ");
				inc =  atof(token);
				token = strtok(NULL, "\n");
				end_val = atof(token);
				analysis = 't';
				num_prints = 0;
				printf("----------TRAN: %d----------\n", tran_counter);
				printf("Timestep: %f, final time:%f, method: %c\n", inc, end_val, method);
			}
			else if(!strcmp("ac", token)){ //ac
				ac_counter++;
				token = strtok(NULL, " ");
				ac_type = strdup(token);
				token = strtok(NULL, " ");
				num_points = atoi(token);
				token = strtok(NULL, " ");
				start_val = atof(token);
				token = strtok(NULL, "\n");
				end_val = atof(token);
				analysis = 'a';
				num_prints = 0;
				printf("------------AC: %d----------\n", ac_counter);
				printf("Type: %s, number of points %d, start frequency: %f, end frequency: %f\n", ac_type, num_points, start_val, end_val);
			}
			break; 
			
			
		}
		else {
			continue; //for an unrecognizable character
		}
		
		val = atof(val_str);
		
		if(strcmp(node1, "0"))
			insert_ht(ht, node1);
			
		if(strcmp(node2, "0"))
			insert_ht(ht, node2);
		
		list_insert(&list, &tran_list, &ac_list, element, name, node1, node2, node3, node4, model_name, l, w, area, val, tran_spec, mag, phase);	
	}
	//print_tran_list(tran_list);
	dc_mna = create_mna_dc(&list, ht, get_ht_size(ht), m2, is_sparse);
		
	x = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
	bdc = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
	if(!is_iter){
		if(!is_sparse){
			//matrices row
			L = (double**)calloc(sizeof(double*), get_ht_size(ht) + m2);
			U = (double**)calloc(sizeof(double*), get_ht_size(ht) + m2);
			//vectors
			p = (int*)calloc(sizeof(int), get_ht_size(ht) + m2);
			pb = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
			//metrices cells
			for(int i = 0; i < (get_ht_size(ht) + m2); i++) {
				L[i] = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
				U[i] = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
			}
		}
		
		y = (double*)calloc(sizeof(double), get_ht_size(ht) + m2);
		
		solve_direct_dcop(ht, dc_mna, is_spd, L, U, p, pb,x, y,  get_ht_size(ht) + m2, is_sparse);
	}
	else{
		solve_iter_dcop(ht, dc_mna, is_spd, x, itol, get_ht_size(ht) + m2, is_sparse);
	}
	
	
	//print_ht(ht);
	//print_list(list);
	//print_mna_dc(dc_mna);
	//print_vector(dc_mna.b, get_ht_size(ht) + m2);
	
	while(1){
		if(fgets(line, IMPORT_LINE_SIZE, fp) == NULL) { //read file line by line
			break;
		}
		
		//convert string to lowercase characters
		i = 0;
		while(line[i] != '\0') {
			line[i] = tolower(line[i]);
			i++;
		}
		
		//replace all tabs with whitespaces
		i = 0;
		while(line[i] != '\0') {
			if(line[i] == '\t') {
				line[i] = ' ';
			}
			i++;
		}
		//printf("%s\n", line);
		if(line[0] == '*') { //ignore lines starting with '*'
			continue;
		}
		else if(line[0] == '.'){
			
			token = strtok(&line[1], " ");
			if(!strcmp("dc" , token)){ //dc sweep
				//If there is a previous sweep before the current
				if(analysis == 'd'){
					if(!is_iter)
						solve_direct_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, L, U, p, pb, bdc, x, y, input_var, start_val, end_val, inc, is_sparse, is_spd);
					else{
						if(is_spd)
							solve_iter_dcsweep_spd(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
						else
							solve_iter_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
					}
					for(int i = 0; i < num_prints; i++){
						plot_dcsweep(input_var, node_names[i], dcsweep_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
							
				}
				//If there is a previous tran before the current dcsweep
				else if(analysis == 't'){
					if(!is_iter){
						solve_direct_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method);
					}
					else{
						solve_iter_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_tran(node_names[i], tran_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				//If there is a previous AC analysis
				else if(analysis == 'a'){
					if(!is_iter){
						solve_direct_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse);
					}
					else{
						solve_iter_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_ac(node_names[i], ac_counter, ac_type);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				
				analysis = 'd';
				token = strtok(NULL, " ");
				strcpy(input_var, token);
				token = strtok(NULL, " ");
				start_val = atof(token);
				token = strtok(NULL, " ");
				end_val =  atof(token);
				token = strtok(NULL, "\n");
				inc = atof(token);
				dcsweep_counter++;
				num_prints = 0;
				printf("--------DC SWEEP: %d---------\n", dcsweep_counter);	
			}
			else if(!strcmp("tran" , token)){
				//If there is a previous sweep before the current
				if(analysis == 'd'){
					if(!is_iter)
						solve_direct_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, L, U, p, pb, bdc, x, y, input_var, start_val, end_val, inc, is_sparse, is_spd);
					else{
						if(is_spd)
							solve_iter_dcsweep_spd(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
						else
							solve_iter_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
					}
					for(int i = 0; i < num_prints; i++){
						plot_dcsweep(input_var, node_names[i], dcsweep_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
							
				}
				//If there is a previous tran before the current tran
				else if(analysis == 't'){
					if(!is_iter){
						solve_direct_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method);
					}
					else{
						solve_iter_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_tran(node_names[i], tran_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				//If there is a previous AC analysis
				else if(analysis == 'a'){
					if(!is_iter){
						solve_direct_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse);
					}
					else{
						solve_iter_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_ac(node_names[i], ac_counter, ac_type);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				
				tran_counter++;
				start_val = 0; //assume starting time = 0
				token = strtok(NULL, " ");
				inc =  atof(token);
				token = strtok(NULL, "\n");
				end_val = atof(token);
				analysis = 't';
				num_prints = 0;
				printf("----------TRAN: %d----------\n", tran_counter);
				printf("Timestep: %f, final time:%f, method: %c\n", inc, end_val, method);
			}
			else if(!strcmp("ac", token)){ //ac
				//If there is a previous dc sweep before the current ac
				if(analysis == 'd'){
					if(!is_iter)
						solve_direct_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, L, U, p, pb, bdc, x, y, input_var, start_val, end_val, inc, is_sparse, is_spd);
					else{
						if(is_spd)
							solve_iter_dcsweep_spd(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
						else
							solve_iter_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
					}
					for(int i = 0; i < num_prints; i++){
						plot_dcsweep(input_var, node_names[i], dcsweep_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
							
				}
				//If there is a previous tran before the current ac
				else if(analysis == 't'){
					if(!is_iter){
						solve_direct_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method);
					}
					else{
						solve_iter_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_tran(node_names[i], tran_counter);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				//If there is a previous AC analysis
				else if(analysis == 'a'){
					if(!is_iter){
						solve_direct_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse);
					}
					else{
						solve_iter_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse, itol);
					}
					for(int i = 0; i < num_prints; i++){
						plot_ac(node_names[i], ac_counter, ac_type);		
					}
					fprintf(stdout, "Click Enter to quit...\n");
    				getchar();
				}
				
				ac_counter++;
				token = strtok(NULL, " ");
				ac_type = strdup(token);
				token = strtok(NULL, " ");
				num_points = atoi(token);
				token = strtok(NULL, " ");
				start_val = atof(token);
				token = strtok(NULL, "\n");
				end_val = atof(token);
				analysis = 'a';
				num_prints = 0;
				printf("------------AC: %d----------\n", ac_counter);
				printf("Type: %s, number of points %d, start frequency: %f, end frequency: %f\n", ac_type, num_points, start_val, end_val);
			}
			else if(!strcmp("print" , token) || !strcmp("plot", token)){
					
				//parsing nodes to print
				while((token = strtok(NULL, "v"))){
					//extracting the node name
					strcpy(node_name, &token[1]);
					j = 0;
					while(node_name[j] != ')'){
						j++;
					}
					node_name[j] = '\0';
						
					//printf("Nodes to print: %s -> Hash value: %d\n", node_name, lookup(ht, node_name));
					if(analysis == 'd')
						sprintf(out_name, "dcsweep%d_%s_%s", dcsweep_counter, input_var, node_name);
					else if(analysis == 't')
						sprintf(out_name, "tran%d_%s", tran_counter, node_name);
					else if(analysis == 'a')
						sprintf(out_name, "ac%d_%s", ac_counter, node_name);
					nodes_print[num_prints] = lookup(ht, node_name);
					strcpy(node_names[num_prints], node_name); 
					out_fp[num_prints] = fopen(out_name, "w+");
					num_prints++;
				}
			}
		}
		else{
			continue;
		}	
	}
	
	//If last dc sweep analysis
	if(analysis == 'd'){
		if(!is_iter)
			solve_direct_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, L, U, p, pb, bdc, x, y, input_var, start_val, end_val, inc, is_sparse, is_spd);
		else{
			if(is_spd)
				solve_iter_dcsweep_spd(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
			else
				solve_iter_dcsweep(dc_mna, ht, list, out_fp, num_prints, nodes_print, itol, bdc, x, input_var, start_val, end_val, inc, is_sparse);
		}
		for(int i = 0; i < num_prints; i++){
			plot_dcsweep(input_var, node_names[i], dcsweep_counter);
		}
		fprintf(stdout, "Click Enter to quit...\n");
    	getchar();	
	}
	//If last tran analysis
	else if(analysis == 't'){
		if(!is_iter){
			solve_direct_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method);
		}
		else{
			solve_iter_tran(dc_mna, ht, list, &tran_list, out_fp, num_prints, nodes_print, end_val, inc, is_sparse, is_spd, method, itol);
		}
		for(int i = 0; i < num_prints; i++){
			plot_tran(node_names[i], tran_counter);		
		}
		fprintf(stdout, "Click Enter to quit...\n");
    	getchar();
	}
	//If last AC analysis
	else if(analysis == 'a'){
		if(!is_iter){
			solve_direct_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse);
		}
		else{
			solve_iter_ac(dc_mna, ht, list, ac_list, out_fp, num_prints, nodes_print, ac_type, num_points, start_val, end_val, is_sparse, itol);
		}
		for(int i = 0; i < num_prints; i++){
			plot_ac(node_names[i], ac_counter, ac_type);		
		}
		fprintf(stdout, "Click Enter to quit...\n");
    	getchar();
	}
	
	//closing files and freeing dynamicly allocated memory
	fclose(fp);
	
	list_delete(&list, &tran_list, &ac_list);
	delete_mna_dc(dc_mna, is_sparse, is_iter);
	
	free(x);
	free(bdc);
	if(!is_iter){
		free(y);
		if(!is_sparse){
			free(pb);
			free(p);
			for(int i = 0; i < (get_ht_size(ht) + m2); i++){
				free(L[i]);
				free(U[i]);
			}
			free(L);
			free(U);
		}
		else{
			cs_sfree(S);
			cs_nfree(N);
		}
	}
	delete_ht(ht);

	return 0;
}



