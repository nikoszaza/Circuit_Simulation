#include <math.h>
#include <complex.h>

double abs_val(double x){
	if(x < 1.0E-15){
		return(-x);
	}
	return(x);
}


void print_array(double** A, int size){
	printf("\n");
	for(int i =0; i < size; i++){
		for(int j =0; j < size; j++){
			printf(" %f ", A[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_vector(double* a, int size){
	printf("\n");
	for(int i =0; i < size; i++){
		printf(" %f ", a[i]);
	}
	printf("\n");
}

void print_vector_int(int* a, int size){
	printf("\n");
	for(int i =0; i < size; i++){
		printf(" %d ", a[i]);
	}
	printf("\n");
}

void lu_dec(double** A, double** L, double** U, int* p, int size){
	double* temp_row;
	double x;
	int m, i, j, temp;
	
	// copying A to U before LU decomposition and adding 1s to L diag
	for(i = 0; i<size; i++){
		for(j = 0; j<size; j++){
			U[i][j] = A[i][j];
			L[i][j] = 0;
		}
		p[i] = i;
		L[i][i] = 1;
	}
	
	
	for(int k = 0; k < size; k++){
		// partial pivoting
		//printf("\nStep %d\n", k);
		//printf("Current array U");
		//print_array(U, 3);
		x = fabs(U[k][k]);
		m = k;
		for(i = k + 1; i < size; i++){
			if( fabs(U[i][k]) > x){
				x = fabs(U[i][k]);
				m = i;
			}
		}
		//printf("Found pivot %f at row %d\n", x, m);
		
		temp_row = U[k];
		U[k] = U[m];
		U[m] = temp_row;
		
		temp = p[k];
		p[k] = p[m];
		p[m] = temp;
		
		
		for(i = k+1; i < size; i++){
			L[i][k] = U[i][k] / U[k][k];
			for(j = k+1; j < size; j++){
				U[i][j] = U[i][j] - L[i][k]*U[k][j];
			}
		}
	}
}


void lu_dec_complex(double complex** A, double complex** L, double complex** U, int* p, int size){
	double complex* temp_row;
	double x;
	int m, i, j, temp;
	
	// copying A to U before LU decomposition and adding 1s to L diag
	for(i = 0; i<size; i++){
		for(j = 0; j<size; j++){
			U[i][j] = A[i][j];
			L[i][j] = 0;
		}
		p[i] = i;
		L[i][i] = 1;
	}
	
	
	for(int k = 0; k < size; k++){
		// partial pivoting
		//printf("\nStep %d\n", k);
		//printf("Current array U");
		//print_array(U, 3);
		x = cabs(U[k][k]);
		m = k;
		for(i = k + 1; i < size; i++){
			if( cabs(U[i][k]) > x){
				x = cabs(U[i][k]);
				m = i;
			}
		}
		//printf("Found pivot %f at row %d\n", x, m);
		
		temp_row = U[k];
		U[k] = U[m];
		U[m] = temp_row;
		
		temp = p[k];
		p[k] = p[m];
		p[m] = temp;
		
		
		for(i = k+1; i < size; i++){
			L[i][k] = U[i][k] / U[k][k];
			for(j = k+1; j < size; j++){
				U[i][j] = U[i][j] - L[i][k]*U[k][j];
			}
		}
	}
}



// Performing cholesky decomposition
// Returns false is A is NOT SPD 
bool chol_dec(double** A, double** L, double** LT, int size ){
	
	for(int i =0 ; i < size; i++){
		for(int j = 0; j < size; j++){
			L[i][j] = 0;
		}
	}
	
	for (int j = 0; j < size; j++) {
    	float sum = 0;
	    for (int k = 0; k < j; k++) 
	        sum += L[j][k] * L[j][k];
	    
	    if(A[j][j] - sum < 1e-15)
	    	return false;
	    	
	    L[j][j] = sqrt(A[j][j] - sum);
	
	    for (int i = j + 1; i < size; i++) {
	        sum = 0;
	        for (int k = 0; k < j; k++) {
	            sum += L[i][k] * L[j][k];
	        }
	        L[i][j] = (1.0 / L[j][j] * (A[i][j] - sum));
	    }
	}
	
	for (int i =0; i < size; i++){
		for(int j = 0; j < size; j++)
			LT[i][j] = L[j][i];
	}
	return true;
}

// Solves Ly = Pb
void solve_ld(double** L, double* pb, int size, double* y){
	double temp;
	for(int k = 0; k < size; k++){
		temp = pb[k];
		for(int j = 0; j < k ; j++){
			temp -= L[k][j]*y[j];
		}
		y[k] = temp / L[k][k];
	}
}

// Solves Ux = y
void solve_ud(double** U, double* y, int size, double* x){
	double temp;
	for(int k = size-1; k >= 0; k--){
		temp = y[k];
		for(int j = k + 1; j < size; j++){
			temp -= U[k][j]*x[j];
		}
		x[k] = temp/U[k][k];
	}
}

//complex versions
void solve_ld_complex(double complex** L, double complex* pb, int size, double complex* y){
	double temp;
	for(int k = 0; k < size; k++){
		temp = pb[k];
		for(int j = 0; j < k ; j++){
			temp -= L[k][j]*y[j];
		}
		y[k] = temp / L[k][k];
	}
}

void solve_ud_complex(double complex** U, double complex* y, int size, double complex* x){
	double temp;
	for(int k = size-1; k >= 0; k--){
		temp = y[k];
		for(int j = k + 1; j < size; j++){
			temp -= U[k][j]*x[j];
		}
		x[k] = temp/U[k][k];
	}
}
