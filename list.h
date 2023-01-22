#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#define PI 3.14159265358979323846


//args[0]:i1, args[1]: i2, args[2]: td1
//args[3]:tc1, args[4]: td2, args[5]:tc2
double exp_tran(double* args, int num_args, double t){
	if(t <= args[2])
		return args[0];
	else if(t <= args[4])
		return args[0] + (args[1]-args[0])*(1-exp(-(t-args[2])/args[3]));
	else 
		return args[0] + (args[1]-args[0])*(exp(-(t-args[4])/args[5])-exp(-(t-args[2])/args[3]));
}

//args[0]: i1, args[1]: ia, args[2]: fr
//args[3]: td, args[4]: df, args[5]: ph
double sin_tran(double* args, int num_args, double t){
	if(t <= args[3])
		return args[0] + args[1]*sin(2*PI*args[5]/360);
	else
		return args[0] + args[1]*sin(2*PI*args[2]*(t-args[3])+ 2*PI*args[5]/360)*exp(-(t-args[3])*args[4]);
}


//args[0]:i1, args[1]:i2, args[2]:td
//args[3]: tr, args[4]: tf, args[5]: pw, args[6]: per
double pulse_tran(double* args, int num_args, double t){
	double t_per = fmod(t,args[6]);
	
	if(t_per <= args[2])
		return args[0];
	else if(t_per <= args[2] + args[3])
		return (args[1]-args[0])/args[3]*(t_per-args[2]) + args[0];
	else if(t_per <= args[2] + args[3]+ args[5])
		return args[1];
	else if(t_per <= args[2] + args[3]+ args[5] + args[4])
		return (args[0]-args[1])/args[4]*(t_per-(args[2] + args[3]+ args[5]))+ args[1];
	else
		return args[0];
}

// args[0]: t1, args[1]: i1
// args[2]: t2, args[3]: i2 ...
double pwl_tran(double* args, int num_args, double t){
	
	for(int i = 0; i < num_args; i += 2){
		if(t < args[i])
			return (args[i+1] - args[i-1])/(args[i] - args[i-2])*(t-args[i-2])+args[i-1];
	}
	
	return 0;
}

//struct for transient specs
typedef struct{
	double* args;
	int num_args;
	double (*fun)(double*, int, double);
}tran_spec_t;


//struct for transistor
typedef struct Transistor{
	char* node3;
	char* node4;
	char* model_name;
	double l;
	double w;
	double area;
}trans_t;


typedef struct Element{
	char el;
	char* name;
	char* node1;
	char* node2;
	double val;
	trans_t* trans_data;
	tran_spec_t* tran_spec;
	struct Element* nxt;
}node_t;

typedef struct List{
	struct Element *head;
	struct Element *tail;
	int size;
}list_t;


// node for transient analysis list
typedef struct Tran_Node{
	node_t* tran_node;
	struct Tran_Node* next;
}tran_node_t;


//transient analysis list
typedef struct {
	tran_node_t* head;
	tran_node_t* tail;
	int size;
}tran_list_t;

//initializes both list and tran_list
void list_init(list_t* list, tran_list_t* tran_list){
	list->head = NULL;
	list->tail = NULL;
	list->size = 0;
	
	tran_list->head = NULL;
	tran_list->tail = NULL;
	tran_list->size = 0;
}


tran_spec_t* tran_spec_create(char* fun_name, double* args, int num_args){
	tran_spec_t* tran_spec = (tran_spec_t*)malloc(sizeof(tran_spec_t));
	//printf("Fun name: %s\n", fun_name);
	if(!strcmp(fun_name,"exp"))
		tran_spec->fun = &exp_tran;
	else if(!strcmp(fun_name,"sin"))
		tran_spec->fun = &sin_tran;
	else if(!strcmp(fun_name,"pulse")){
		tran_spec->fun = &pulse_tran;
		//printf("Hi\n");
	}
	else if(!strcmp(fun_name,"pwl"))
		tran_spec->fun = &pwl_tran;
	else{
		tran_spec->fun = NULL;
		//printf("Hello\n");
	}	
	tran_spec->num_args = num_args;
	tran_spec->args = (double*)malloc(sizeof(double)*num_args);
	
	for(int i =0; i<num_args; i++){
		tran_spec->args[i] = args[i];
	}
	return tran_spec;
}

//deletes both list and tran_list
void list_delete(list_t* list, tran_list_t* tran_list) {
	int i;
	struct Element *current, *next;
	
	tran_node_t* tran_cur, *tran_next;
	
	for(i = 0, current = list->head; i < list->size; i++, current = next) {
		free(current->name);
		free(current->node1);
		free(current->node2);
		if(current->trans_data != NULL) {
			if(current->el == 'd'){
				free(current->trans_data->model_name);
				free(current->trans_data);
			}
			else if(current->el == 'q'){
				free(current->trans_data->node3);
				free(current->trans_data->model_name);
				free(current->trans_data);
			}
			else if(current->el == 'm'){
				free(current->trans_data->node3);
				free(current->trans_data->node4);
				free(current->trans_data->model_name);
				free(current->trans_data);
			}
		}
		next = current->nxt;
		free(current);
	}
	
	tran_cur = tran_list->head;
	while(tran_cur != NULL){
		tran_next = tran_cur->next;
		free(tran_cur);
		tran_cur = tran_next;
	}
}

void tran_list_insert(tran_list_t* tran_list, node_t* node){
	tran_node_t* new_tran_node = malloc(sizeof(tran_node_t));
	new_tran_node->tran_node = node;
	new_tran_node->next = NULL;
	
	if(tran_list->size == 0){
		tran_list->head = tran_list->tail = new_tran_node;
		tran_list->size++;
	}
	else{
		tran_list->tail->next = new_tran_node;
		tran_list->tail = new_tran_node;
		tran_list->size++;
	}
}

void list_insert(list_t* list, tran_list_t* tran_list, char el, char* name, char* node1, char* node2, char* node3, char* node4, char* model_name, double l, double w, double area, double val, tran_spec_t* tran_spec){
	node_t* new_node = malloc(sizeof(node_t));
	new_node->el = el;
	
	new_node->name = malloc(sizeof(char)*(strlen(name) + 1));
	strcpy(new_node->name, name);
	
	new_node->node1 = malloc(sizeof(char)*(strlen(node1) + 1));
	strcpy(new_node->node1, node1);
	
	new_node->node2 = malloc(sizeof(char)*(strlen(node2) + 1));
	strcpy(new_node->node2, node2);
	
	new_node->tran_spec = tran_spec;
	if(tran_spec){
		tran_list_insert(tran_list, new_node);
	}
		
	
	if(el == 'd'){
		new_node->trans_data = malloc(sizeof(trans_t));
		new_node->trans_data->model_name = malloc(sizeof(char)*(strlen(model_name)+1));
		strcpy(new_node->trans_data->model_name, model_name);
			
		new_node->trans_data->area = area;
	}
	else if(el == 'q'){
		new_node->trans_data = malloc(sizeof(trans_t));	
		
		new_node->trans_data->node3 = malloc(sizeof(char)*(strlen(node3)+1));
		strcpy(new_node->trans_data->node3, node3);
		
		new_node->trans_data->model_name = malloc(sizeof(char)*(strlen(model_name)+1));
		strcpy(new_node->trans_data->model_name, model_name);
		
		new_node->trans_data->area = area;
	}
	else if(el == 'm'){
		new_node->trans_data = malloc(sizeof(trans_t));	
		
		new_node->trans_data->node3 = malloc(sizeof(char)*(strlen(node3)+1));
		strcpy(new_node->trans_data->node3, node3);
		
		new_node->trans_data->node4 = malloc(sizeof(char)*(strlen(node4)+1));
		strcpy(new_node->trans_data->node4, node4);
		
		new_node->trans_data->model_name = malloc(sizeof(char)*(strlen(model_name)+1));
		strcpy(new_node->trans_data->model_name, model_name);
		
		new_node->trans_data->l = l;
		new_node->trans_data->w = w;
	}
	else{
		new_node->trans_data = NULL;
	}
	
	new_node->val = val;
	new_node->nxt = NULL;
	
	if(list->size == 0){
		list->head = list->tail = new_node;
		list->size++;
	}
	else{
		list->tail->nxt = new_node;
		list->tail = new_node;
		list->size++;
	}
}

node_t* return_element(list_t list, char *var_name, char el){
	node_t* cur = list.head;
	
	while(cur != NULL){
		if((el == cur->el) && !strcmp(cur->name, var_name))
			return cur; 
		cur = cur->nxt;
	}
	
	return NULL;
}

int k_value(list_t list, char *name){
	int k = -1;
	node_t* cur = list.head;
	
	while(cur != NULL){
		if( ( (cur->el == 'v') && (strcmp(name, cur->name)) ) || (cur->el == 'l') )
			k++;
		else if((cur->el == 'v') && (!strcmp(name, cur->name))){
			k++;
			return k;
		}
		cur = cur->nxt;
	}
	return -1;
}


void print_list(list_t list){
	node_t* cur = list.head;
	
	printf("---------LIST-------------\n");
	while(cur != NULL){
		if(cur->el == 'd'){
			printf("d%s %s %s %s %f\n", cur->name, cur->node1, cur->node2, cur->trans_data->model_name, cur->trans_data->area);
		} else if(cur->el == 'q'){
			printf("q%s %s %s %s %s %f\n", cur->name, cur->node1, cur->node2, cur->trans_data->node3, cur->trans_data->model_name, cur->trans_data->area);
		} else if(cur->el == 'm'){
			printf("m%s %s %s %s %s %s L=%f W=%f\n", cur->name, cur->node1, cur->node2, cur->trans_data->node3, cur->trans_data->node4, cur->trans_data->model_name, cur->trans_data->l, cur->trans_data->w);
		} else {
			printf("%c%s %s %s %f\n", cur->el, cur->name, cur->node1, cur->node2, cur->val);	
		}
		
		cur = cur->nxt;
	}
	printf("----------END-------------\n");
}

void print_tran_list(tran_list_t tran_list){
	tran_node_t* cur = tran_list.head;
	
	printf("------TRANSIENT LIST--------\n");
	while(cur != NULL){
		printf("Node %c%s has transient spec\n", cur->tran_node->el, cur->tran_node->name);
		cur = cur->next;
	}
	printf("----------END-------------\n");
}

//returns number of non-zero for array
//array_type = A or B
//A is for i,l,v,r
//B is for c,l
int get_nz(list_t list, char array_type){
	int nz = 0;
	node_t* cur = list.head;
	
	if(array_type == 'A'){
		while(cur != NULL){
			if(cur->el == 'r'){
				if( !strcmp(cur->node1,"0") || !strcmp(cur->node2,"0"))
					nz++;
				else
					nz += 4;
			}
			else if(cur->el == 'v' || cur->el == 'l'){
				if( !strcmp(cur->node1,"0") || !strcmp(cur->node2,"0"))
					nz += 2;
				else
					nz += 4;
			}
			
			cur = cur->nxt;
		}
	}
	else{
		while(cur != NULL){
			if(cur->el == 'c'){
				if( !strcmp(cur->node1,"0") || !strcmp(cur->node2,"0"))
					nz++;
				else
					nz += 4;
			}
			else if(cur->el == 'l'){
				nz++;
			}
			
			cur = cur->nxt;
		}	
	}
	return nz;
}
