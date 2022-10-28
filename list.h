#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

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
	struct Element* nxt;
}node_t;

typedef struct List{
	struct Element *head;
	struct Element *tail;
	int size;
}list_t;


void list_init(list_t* list){
	list->head = NULL;
	list->tail = NULL;
	list->size = 0;
}

void list_delete(list_t* list) {
	int i;
	struct Element *current, *next;
	
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
}

void list_insert(list_t* list, char el, char* name, char* node1, char* node2, char* node3, char* node4, char* model_name, double l, double w, double area, double val){
	node_t* new_node = malloc(sizeof(node_t));
	new_node->el = el;
	
	new_node->name = malloc(sizeof(char)*(strlen(name) + 1));
	strcpy(new_node->name, name);
	
	new_node->node1 = malloc(sizeof(char)*(strlen(node1) + 1));
	strcpy(new_node->node1, node1);
	
	new_node->node2 = malloc(sizeof(char)*(strlen(node2) + 1));
	strcpy(new_node->node2, node2);
	
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
