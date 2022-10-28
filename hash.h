#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// hashtable for mapping
// works for insert, not delete/pop element

typedef struct HashTableEntry{
	char* string;
	int id; // the identifier for the string	
}htentry_t;


typedef struct HashTable{
	htentry_t** table;	
	int size;
	int capacity;
}ht_t;

ht_t* init_ht(int capacity){
	ht_t* ht = (ht_t*)malloc(sizeof(ht_t));
	ht->size = 0;
	ht->capacity = capacity;
	ht->table = malloc(sizeof(htentry_t*)*capacity); 	
	
	for(int i = 0; i<capacity; i++){
		(ht->table)[i] = NULL;
	}
	
	return ht;
}


int hash_fun(char* string){

	int hash_key = 0;
	for(int i = 0; i<strlen(string); i++){
		hash_key += string[i] << i ;
	}
	
	return hash_key;
}


// lookup: searches for string in hashtable and retuns its id, otherwise it returns -1
int lookup(ht_t* ht, char* str){
	int pos = hash_fun(str) % ht->capacity;
	
	while(1){
		// not in ht
		if( (ht->table)[pos] == NULL){
			return -1;
		}
		
		// found
		if( strcmp( (ht->table)[pos]->string , str ) == 0 ){
			return (ht->table)[pos]->id;
		}
		
		// not in this pos, check next pos
		pos = (pos + 1) % ht->capacity;
	}	
}


int get_ht_size(ht_t* ht){
	return ht->size;
}


// used by rehash
void insert_new_ht(ht_t* new_ht, htentry_t* htentry){
	int pos = hash_fun(htentry->string) % new_ht->capacity;
	
	while(1){
		if( (new_ht->table)[pos] == NULL ){
			(new_ht->table)[pos] = htentry;
			new_ht->size++;
			return;		
		}
		
		pos = (pos + 1) % new_ht->capacity;
	}
}

// doubles table capacity when size/capacity > threshold
// returns new hashtable
void rehash(ht_t* ht){
	ht_t* new_ht = init_ht(2*(ht->capacity));
	
	for(int i = 0; i < ht->capacity; i++){
		if( (ht->table)[i] != NULL ){
			insert_new_ht(new_ht, (ht->table)[i]);
		}
	}
	ht->capacity = 2*ht->capacity;
	ht->table = realloc(ht->table, sizeof(htentry_t*)*ht->capacity);
	for(int i = 0; i < new_ht->capacity; i++){
		(ht->table)[i] = (new_ht->table)[i];
	}
	
	free(new_ht->table);
	free(new_ht);
}

// returns true if insertion was succesful, otherwise false
bool insert_ht(ht_t* ht, char* str) {
	int pos;
	htentry_t* new_entry_ptr;
	
	if( (float)ht->size/(float)ht->capacity > 0.7 ){
		rehash(ht);
	}
	
	pos = hash_fun(str) % ht->capacity;
	while(1){
		
		// not in ht
		if( (ht->table)[pos] == NULL){
			
			new_entry_ptr = malloc(sizeof(htentry_t));
			new_entry_ptr->string = malloc(sizeof(char)*(strlen(str) + 1));
			new_entry_ptr->id = ht->size;
			strcpy(new_entry_ptr->string, str);
			(ht->table)[pos] = new_entry_ptr;
			(ht->size)++;
			//printf("Insert for %s in pos %d OK\n", str, pos);
			
			return true;
		}
		
		// found
		if( strcmp( (ht->table)[pos]->string , str ) == 0 ){
			// printf("Insert for %s in pos %d NOK\n", str, pos);
			return false;
		}
		
		// pos has other element, check next pos
		pos = (pos + 1) % ht->capacity;
	}		
}

void delete_ht(ht_t* ht){
	for(int i = 0; i < ht->capacity; i++ ){
		if( (ht->table)[i] != NULL ) {
			free(ht->table[i]->string);
			free(ht->table[i]);
		}
	}
	free(ht->table);
	free(ht);
}

void print_ht(ht_t* ht){
	printf("------------HashTable----------\n");
	printf("Size = %d, Capacity = %d\n", ht->size, ht->capacity);
	for(int i = 0; i < ht->capacity; i++ ){
		if( (ht->table)[i] != NULL ) {
			printf("Pos %d: %s -> %d\n", i, (ht->table)[i]->string, (ht->table)[i]->id);
		}
	}
}

/*int main(void){
	char str1[] = "1", str2[] = "22", str3[] = "3", str4[] = "412", str5[] = "a", str6[] = "gg";
	ht_t* ht = init_ht(3); 
	insert_ht(ht, str1);
	insert_ht(ht, str2);
	insert_ht(ht, str3);
	insert_ht(ht, str4);
	insert_ht(ht, str5);
	insert_ht(ht, str6);
	print_ht(ht);
	delete_ht(ht);
	return 0;
}*/
