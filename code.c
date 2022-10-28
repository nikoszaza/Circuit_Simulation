#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "hash.h"
#include "list.h"
#define IMPORT_LINE_SIZE 128
#define MAX_LEN 64

typedef struct {
	double** A;
	double* b;
	int m1, m2;
}dc_mna_t;

dc_mna_t create_mna_dc(list_t* list, ht_t* ht, int m1, int m2){
	int m2_counter = 0;
	dc_mna_t dc_mna;
	dc_mna.m1 = m1;
	dc_mna.m2 = m2;
	dc_mna.A = (double**)calloc((m1+m2), sizeof(double*));
	dc_mna.b = (double*)calloc((m1+m2), sizeof(double));
	node_t* curr = list->head;
	
	for(int i=0; i < (m1+m2); i++){
		dc_mna.A[i] = (double*)calloc((m1+m2), sizeof(double));
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
				m2_counter++;
				break;
			case 'i':
				if(node1_code >= 0) {
					if(node2_code >= 0) {
						dc_mna.b[node1_code] -= curr->val;
						dc_mna.b[node2_code] += curr->val;
					} else {
						dc_mna.b[node1_code] -= curr->val;	
					}
				} else {
					if(node2_code >=0) {
						dc_mna.b[node2_code] += curr->val;
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
						dc_mna.A[node2_code][node2_code] -= 1/curr->val;
					}
				}
				break;
			case 'c':
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
				m2_counter++;
				break;
			default:
				break;
		}
		
		curr = curr->nxt;
	}
	
	return dc_mna;
}

void print_mna_dc(dc_mna_t dc_mna){
	printf("-----MNA DC ANALYSIS----\n");
	printf("m1 elements = %d\nm2 elements = %d\n", dc_mna.m1, dc_mna.m2);
	
	
	printf("----------A------------\n");
	
	for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
		for(int j=0; j < (dc_mna.m1+ dc_mna.m2); j++){
			printf("%.6f ", dc_mna.A[i][j]);
		}
		printf("\n");
	}
	
	printf("\n----------B------------\n");
	for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
		printf("%.6f ", dc_mna.b[i]);
	}
	printf("\n");
}

void delete_mna_dc(dc_mna_t dc_mna){
	for(int i=0; i < (dc_mna.m1 + dc_mna.m2); i++){
		free(dc_mna.A[i]);
	}
	
	free(dc_mna.A);
	free(dc_mna.b);
}



int main(int argc, char *argv[]) {
	FILE * fp;
	char line[IMPORT_LINE_SIZE], element, name[MAX_LEN], node1[MAX_LEN], node2[MAX_LEN], w_str[MAX_LEN], l_str[MAX_LEN];
	char node3[MAX_LEN], node4[MAX_LEN], val_str[MAX_LEN], model_name[MAX_LEN], area_str[MAX_LEN];
	char *token = NULL;
	int i, m2 = 0;
	double val, area, w, l;
	ht_t* ht = init_ht(5);	
	list_t list;
	dc_mna_t dc_mna;
	
	list_init(&list);
	

	if(argc == 2){
		fp = fopen(argv[1], "r");
	}
	else{
		printf("Too few arguments!\n");
		return(-1);
	}
	
	if(!fp){
		printf("Invalid argument\n");
		return(-2);
	}
	
	while(1) {
		area = w = l = 0;
		
		
		if(fgets(line, IMPORT_LINE_SIZE, fp) == NULL) { //read file line by line
			break;
		}
		printf("line is: %s", line);
		
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

			token = strtok(&line[1], " ");
			strcpy(name, token);
			token = strtok(NULL, " ");
			strcpy(node1, token);
			token = strtok(NULL, " ");
			strcpy(node2, token);
			token = strtok(NULL, " ");
			strcpy(val_str, token);
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
		else {
			continue;
		}
		
		val = atof(val_str);
		
		if(strcmp(node1, "0"))
			insert_ht(ht, node1);
			
		if(strcmp(node2, "0"))
			insert_ht(ht, node2);
		
		list_insert(&list, element, name, node1, node2, node3, node4, model_name, l, w, area, val);	
	}
	dc_mna = create_mna_dc(&list, ht, get_ht_size(ht), m2);
	
	print_ht(ht);
	print_list(list);
	print_mna_dc(dc_mna);
	
	//closing file and freeing dynamicly allocated memory
	fclose(fp);
	delete_ht(ht);
	list_delete(&list);
	delete_mna_dc(dc_mna);
	
	return 0;
}
