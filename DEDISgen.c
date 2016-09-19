/* DEDISbench
 * (c) 2010 2010 U. Minho. Written by J. Paulo
 */

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <openssl/sha.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include "berk.c"

//max path size of a folder/file
#define MAXSIZEP 10000
#define FOLDER 1
#define DEVICE 2

//TODO for now this is static
#define READSIZE 5242880


//TODO define these new variables
int nr_sizes_proc=0;
int *sizes_proc;


//TODO in the future this will be an argument
#define PRINTDB "gendbs/printdb/"
#define DUPLICATEDB "gendbs/duplicatedb/"


//total blocks scanned
uint64_t *total_blocks;
//identical blocks (that could be eliminated)
uint64_t *eq;
//distinct blocks
uint64_t *dif;
//distinct blocks with duplicates
uint64_t *distinctdup;
//blocks that were appended with zeros due to their size
//uint64_t *zeroed_blocks;

//duplicated disk space
uint64_t *space;


//number files scaneed
uint64_t nfiles=0;

//This is the processed blocks of READSIZE so it is not an array
uint64_t processed_blocks=0;

//duplicates DB
DB ***dbporiginal; // DB structure handle
DB_ENV ***envporiginal;


//print file DB
DB ***dbprinter; // DB structure handle
DB_ENV ***envprinter;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX_SIZE 256


struct define_entropy {
	double ent;
	struct define_entropy_operations *buffer_to_populate;
};

struct define_entropy_operations {
	int (*fill)(struct define_entropy *, void *, unsigned int);
};

int fill_nothing_inside_block(struct define_entropy *ds, void *buffer_to_populate, unsigned int size);
int populate_DEDISbench_like_initial(struct define_entropy *ds, void *buffer_to_populate, unsigned int size);
int contigous_block_fill(struct define_entropy *ds, void *buffer_to_populate, unsigned int size);
int contigous_block_fill_random(struct define_entropy *ds, void *buffer_to_populate, unsigned int size);
int populating_buffer_with_actual_entropy(struct define_entropy *ds, void *buffer_to_populate, unsigned int size);
int random_shuffle(unsigned char array[], unsigned int size);
double calculate_entropy_of_prob_dist(double prob_dist[], unsigned int size);
double test_check_buffer_to_populate_output(void *buffer_to_populate, unsigned int size);
int binary_search(double index, double array[], unsigned int size);
int cum_dist_distribution_entropy_gen(double prob_dist[], unsigned int size, double cum_dist[]);
double secant_root(double (*func)(double, double, double), double lower_lim, double higher_lim, double begin, double end);
double transformation(double lower_lim, double higher_lim, double value_found);
double finding_root(double (*func)(double, double, double), double lower_lim, double higher_lim);
int prob_dist_distribution_entropy_gen(double prob_dist[], int size , double ent);


/*
	buffer_to_populatefer will not change in allocated memory, and will be filled with zeroes
*/
int fill_nothing_inside_block(struct define_entropy *ds, void *buffer_to_populate, unsigned int size){
	return 0;
}

/*
	This function is filling the buffer_to_populatefer with strings just like in DEDISbench with constant characters (can be changed to any character)
*/
int populate_DEDISbench_like_initial(struct define_entropy *ds, void *buffer_to_populate, unsigned int size){

	int i=0;
	for(i=0; i < size; i++){
		((char *)buffer_to_populate)[i] = 0;
	}
	return 0;
}

/*
	contigous_block_fill fills the buffer_to_populate with random data according to
	a prob_dist has the entropy specified in the define_entropy.

	It calculates the cum_dist then populates the buffer_to_populatefer by searching in the cum_dist for
    the random number.
*/

int contigous_block_fill(struct define_entropy *ds, void *buffer_to_populate, unsigned int size){

	int i = 0;
	double prob_dist[MAX_SIZE];
	double cum_dist[MAX_SIZE];
	unsigned char characters_available[MAX_SIZE];

	//Calculate the probability distribution according to the given entropy
	prob_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, ds-> ent);

	//Calculate cum_dist from the prob_dist
	cum_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, cum_dist);

	//initializing the character table available which is 256 for a maximum of 8 bits/byte
	for(i=0; i< MAX_SIZE; i++)
		characters_available[i] = (unsigned char)i;

	//shuffle the symbols table
	random_shuffle(characters_available, MAX_SIZE);

	for(i=0; i < size; i++){
		((unsigned char*)buffer_to_populate)[i] = characters_available[binary_search(rand()/(double)RAND_MAX, cum_dist, MAX_SIZE)];
	}
	return 0;
}

/*
	contigous_block_fill_random fills the buffer_to_populatefer with random data according to
	the prob_dist.

	It will generate contiguous segments of data in the buffer_to_populate to make different same
	size buffer_to_populate look different (shuffle the symbols table).
*/

int contigous_block_fill_random(struct define_entropy *ds, void *buffer_to_populate, unsigned int size){

	int i = 0;
	double prob_dist[MAX_SIZE];
	double cum_dist[MAX_SIZE];
	unsigned char characters_available[MAX_SIZE];

	//Calculate prob_dist according to the given entropy
	prob_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, ds-> ent);


	//Calculate cum_dist from the prob_dist
	cum_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, cum_dist);

	//character table initialization
	for(i=0; i< MAX_SIZE; i++)
		characters_available[i] = (unsigned char)i;

	//shuffle the character table
    random_shuffle(characters_available, MAX_SIZE);

	int k=0;
	int j=0;

	for(i=0; i< MAX_SIZE; i++){
		for(j=0; j < (int)(prob_dist[i]*size); j++){
			((unsigned char*)buffer_to_populate)[k] = characters_available[i];
			k++;
		}
	}
	//remaining elements
	for(; k < size; k++){
		((unsigned char *)buffer_to_populate)[k] = characters_available[binary_search(rand()/(double)RAND_MAX, cum_dist, MAX_SIZE)];
	}

	return 0;
}


/*allocating memory for buffer and filling with data from contigous_block_fill_random.

This is an example of how the data can populate a given buffer.
*/


int populating_buffer_with_actual_entropy(struct define_entropy *ds, void *buffer_to_populate, unsigned int size){
	void *tmp = malloc(size);
    int err = contigous_block_fill_random(ds, tmp, size);
    memcpy(buffer_to_populate, tmp, size);
    free(tmp);
    return err;
}


/*generates a random number between start and end for the permutation of the character table*/

inline int random_int(int start, int end){
	return start + rand()%(end - start);
}

/*swaps two values in the array*/

inline void swap(unsigned char array[], int i, int j){
	int temp = array[i];
	array[i] = array[j];
	array[j] = temp;
}

/*uses the above two functions to randomly shuffle the values in the array*/\

int random_shuffle(unsigned char array[], unsigned int size){

	int i,j;
	if(size <= 0) return -1;
	for(i=0; i < size; i++){
		j = random_int(i,size);
		swap(array, i, j);
	}
	return 0;
}

/*uses the Shannon entropy formula to caluclate the probability density and as a result the entropy */

double calculate_entropy_of_prob_dist(double prob_dist[], unsigned int size){

	double entropy =0;
	int i= 0;

	if(size <= 0)
		return -1;

	for(i=0; i<size; i++){
		if(prob_dist[i] > 0)
			entropy+= prob_dist[i]*log2(1.0/prob_dist[i]);
	}
	return entropy;
}

/*checking for whether the buffer entropy is equal to the user defined entropy or marginally close*/

double test_check_buffer_to_populate_output(void *buffer_to_populate, unsigned int size){

	int i=0;
	double prob_dist[MAX_SIZE];
	for(i=0; i < MAX_SIZE; i++)
		prob_dist[i] = 0.0;
	for(i=0; i < size; i++)
		prob_dist[((unsigned char*)buffer_to_populate)[i]]+= 1.0/size;
	return calculate_entropy_of_prob_dist(prob_dist, MAX_SIZE);
}

int binary_search(double index, double array[], unsigned int size){

	int start = 0;
	int end = size-1;
	int mid;

	if(size <= 0)
		return -2;
	if(index <= array[0])
		return 0;
	while(end-start > 1){
		mid = (end+start)/2;
		if(index == array[mid])
			return mid;
		if(index < array[mid])
			end = mid;
		else
			start = mid;
	}
	return end;
}

/*cumalative distribution function that is derived from the probability distribution and is a non-decreasing function*/

int cum_dist_distribution_entropy_gen(double prob_dist[], unsigned int size, double cum_dist[]){

	int i;
	if (size <= 0)
		return -1;

	cum_dist[0] = prob_dist[0];
	for(i=1; i < size; i++){
		cum_dist[i] = cum_dist[i-1] + prob_dist[i];
	}
	return 0;
}

/*method used to find roots for the equation that is dependent on one variable - eplison to generate the user-defined probability density, cumalitive density and character table to populate the buffer of a given size */

double secant_root(double (*func)(double, double, double), double lower_lim, double higher_lim, double begin, double end){
	int i =0;
	double a, b;
	for(i=0; i<20 ;i++){
		b = func(lower_lim, higher_lim, end);
		a = (end-begin) / (b - func(lower_lim, higher_lim, begin)) * b;
		if(fabs(a) < 5e-11)
			return end;
		begin = end;
		end = end - a;
	}
	return end;
}

/* the variable - eplison that is shifted to the other side of the Shannon entropy equation to and used to find the root of the equation*/

double transformation(double lower_lim, double higher_lim, double value_found){
	double prob = 1.0/lower_lim;
    double shift_func = higher_lim - ((lower_lim-2.0)/lower_lim)*log2(1.0/prob);
    return (prob+value_found)*log2(1.0/(prob+value_found)) + (prob-value_found)*log2(1.0/(prob-value_found)) - shift_func;
}

/*this function calls the above functions to find the root using a lower and higher limit of -1e-10+ 1.0/lower_lim and
1e-10 respectively */

double finding_root(double (*func)(double, double, double), double lower_lim, double higher_lim){
	return secant_root(func, lower_lim, higher_lim, 1e-10, -1e-10+ 1.0/lower_lim);
}

/* this function produces the user defined entropy by calculating the root from the secant method
and provides the manipulated probability distribution */

int prob_dist_distribution_entropy_gen(double prob_dist[], int size , double ent){

	int i = 0;
	if(size <= 0)
		return -1;
	for(i=0; i< size; i++)
		prob_dist[i]=0.0;
    double temp = (int) ceil(pow(2.0,ent));
    for(i = 0; i < (int)temp; i++){
        prob_dist[i] = 1.0/temp;
	}
    double reduction = finding_root(transformation, temp, ent);
    prob_dist[0] += reduction;
    prob_dist[1] -= reduction;
	return 0;
}



//calculate block checksum and see if already exists at the hash table
int check_duplicates(unsigned char* block,uint64_t bsize,int id_blocksize){

  //calculate sha1 hash
  unsigned char *result = malloc(sizeof(unsigned char)*20);
  result[0]=0;
  SHA1(block, bsize,result);

  //generate a visual (ASCII) representation for the hashes
    char *start=malloc(sizeof(char)*41);
    start[40]=0;
    int j =0;

    for(j=0;j<40;j+=2)
    {
	start[j]=(result[j/2]&0x0f)+'a';
	start[j+1]=((result[j/2]&0xf0)>>4)+'a';
    }

    free(result);

    //one more block scanned
    total_blocks[id_blocksize]++;

    //get hash from berkley db
    struct hash_value hvalue;

    //printf("before introducing in db\n");
    int ret = get_db(start,&hvalue,dbporiginal[id_blocksize],envporiginal[id_blocksize]);
    

    //if hash entry does not exist
    if(ret == DB_NOTFOUND){
    	// printf("not found\n");

       //unique blocks++
       dif[id_blocksize]++;

       //set counter to 1
       hvalue.cont=1;

       //insert into into the hashtable
       put_db(start,&hvalue,dbporiginal[id_blocksize],envporiginal[id_blocksize]);
    }
    else{

    	// printf("after search\n");

    	if(hvalue.cont==1){
    		//this is a distinct block with duplicate
    		distinctdup[id_blocksize]++;
    	}

    	//found duplicated block
        eq[id_blocksize]++;
    	//space that could be spared
    	space[id_blocksize]=space[id_blocksize]+bsize;

    	//increase counter
        hvalue.cont++;

        //insert counter in the right entry
        put_db(start,&hvalue,dbporiginal[id_blocksize],envporiginal[id_blocksize]);

    }

    free(start);

    processed_blocks++;
    if(processed_blocks%100000==0){
    	printf("processed %llu blocks\n",(long long unsigned int) processed_blocks);
    }

  return 0;
}


//given a file extract blocks and check for duplicates
int extract_blocks(char* filename){

	printf("Processing %s \n",filename);
    int fd = open(filename,O_RDONLY | O_LARGEFILE);
    uint64_t off=0;
    if(fd){

      //declare a block
      unsigned char block_read[READSIZE];
      bzero(block_read,sizeof(block_read));

      //read first block from file
      int aux = pread(fd,block_read,READSIZE,off);
      off+= READSIZE;

      //check if the file still has more blocks and if size <bzise discard
      //the last incomplete block
      while(aux>0){

      	//Process fixed size dups all sizes specified
         int curr_sizes_proc=0;
         while(curr_sizes_proc<nr_sizes_proc){

        	 int size_block=sizes_proc[curr_sizes_proc];
        	 int size_proced=0;
        	 unsigned char *block_proc;

        	 while(aux-size_proced>=size_block){

        		 block_proc=malloc(size_block);
        		 bzero(block_proc,size_block);

        		 memcpy(block_proc,&block_read[size_proced],size_block);

        		 //index the block and find duplicates
        		 check_duplicates(block_proc,size_block,curr_sizes_proc);

        		 free(block_proc);

        		 size_proced+=size_block;
        	 }

        	 curr_sizes_proc++;
       	 }

       	 //free this block from memory and process another
         //free(block_read);
         //block_read=malloc(sizeof(unsigned char)*READSIZE);
         aux = pread(fd,block_read,READSIZE,off);
         off+=READSIZE;

      }

      
    close(fd);
    }
    else{
      fprintf(stderr,"error opening file %s\n",filename);
      exit(1);

    }

  return 0;

}

//search a directory for files inside and return the number of files found and their path(nfiles,lfiles)
int search_dir(char* path){

 //directory information
 struct dirent *dp=malloc(sizeof(struct dirent));

 //opens the path and check the files inside
 DIR *dirp = opendir(path);
 if(dirp){
	//while opendir is not null
    while ((dp = readdir(dirp)) != NULL){
      //exclude . and ..
      if(strcmp(dp->d_name,".")!=0 && strcmp(dp->d_name,"..")!=0){

         //build the full path
         char newpath[MAXSIZEP];

         strcpy(newpath,path);
         strcat(newpath,dp->d_name);

         if(dp->d_type==DT_DIR){
               strcat(newpath,"/");
               // recursively process the files inside the diretory
               search_dir(newpath);
         }
         else{
            if(dp->d_type==DT_REG){
            	//If it is a regular file then start segmenting files and indexing blocks

            	extract_blocks(newpath);
            	nfiles++;

            }

         }
      }
    }
 }
 closedir(dirp);
 free(dp);

 //printf("return search dir\n");

 return 0;
}



int gen_output(DB **dbpor,DB_ENV **envpor,DB **dbprint,DB_ENV **envprint){

	//Iterate through original DB and insert in print DB
	int ret;

	DBT key, data;

	DBC *cursorp;

	(*dbpor)->cursor(*dbpor, NULL, &cursorp, 0);

	// Initialize our DBTs.
	memset(&key, 0, sizeof(DBT));
	memset(&data, 0, sizeof(DBT));
	// Iterate over the database, retrieving each record in turn.
	while ((ret = cursorp->get(cursorp, &key, &data, DB_NEXT)) == 0) {

	   //get hash from berkley db
	   struct hash_value hvalue;
	   uint64_t ndups = (unsigned long long int)((struct hash_value *)data.data)->cont;
	   ndups=ndups-1;
	   //char ndups[25];

	   //key
	   //sprintf(ndups,"%llu",(unsigned long long int)((struct hash_value *)data.data)->cont);

	   //see if entry already exists and
	   //Insert new value in hash for print number_of_duplicates->numberof blocks
	   int retprint = get_db_print(&ndups,&hvalue,dbprint,envprint);

	   //if hash entry does not exist
	   if(retprint == DB_NOTFOUND){

		  hvalue.cont=1;
		  //insert into into the hashtable
		  put_db_print(&ndups,&hvalue,dbprint,envprint);
       }
	   else{

		  //increase counter
		  hvalue.cont++;
		  //insert counter in the right entry
		  put_db_print(&ndups,&hvalue,dbprint,envprint);
       }

	}
	if (ret != DB_NOTFOUND) {
	    perror("failed while iterating");
	}


	if (cursorp != NULL)
	    cursorp->close(cursorp);

	return 0;
}


//Aux function to split the string with multiple sizes
void str_split(char* a_str)
{
    char* tmp = a_str;
   
    //starts with one because one element does not need the comma
    while (*tmp)
    {
        if (',' == *tmp)
        {
            nr_sizes_proc++;
        }
        tmp++;
    }

  	//increment one more due to the last str
  	nr_sizes_proc++;

    sizes_proc = malloc(sizeof(int)*nr_sizes_proc);

    char* token = strtok(a_str, ",");

    int i=0;
    while (token)
    {
        sizes_proc[i] = atoi(token);
        token = strtok(NULL, ",");
        i++;
    }

}

void usage(void)
{
	printf("Usage:\n");
	printf(" -f or -d\t(Find duplicates in folder -f or in a Disk Device -d)\n");
	printf(" -p<value>\t\t(Path for the folder or disk device)\n");
	printf(" -h\t\t\t(Help)\n");
	exit (8);
}

void help(void){

	printf(" Help:\n\n");
	printf(" -f or -d\t(Find duplicates in folder -f or in a Disk Device -d)\n");
	printf(" -p<value>\t\t(Path for the folder or disk device)\n");
	printf("\n Optional Parameters\n\n");
	printf(" -o<value>\t\t(Path for the output distribution file. If not specified this is not generated.)\n");
	printf(" -z<value>\t(Path for the folder where duplicates databases are created default: ./gendbs/duplicatedb/)\n");
	printf(" -b<value>\t( Size of blocks to analyse in bytes eg: -b1024,4096,8192 default: -b4096\n");
	exit (8);

}




int main (int argc, char *argv[]){

	//directory or disk device
	int devicetype=-1;
    //path to the device
	char devicepath[100];

	//path output log
	int outputfile=0;
	char outputpath[100];

	//path of databases folder
	int dbfolder=0;
	char dbfolderpath[100];

	
  	while ((argc > 1) && (argv[1][0] == '-'))
  	{
		switch (argv[1][1])
		{
			case 'f':
			//Test if -d is not being used also
			if(devicetype!=DEVICE)
				devicetype=FOLDER;
			else{
			   printf("Cannot use both -f and -d\n");
			   usage();
			}
			break;
			case 'd':
			//test if -f is not being used also
			if(devicetype!=FOLDER)
				devicetype=DEVICE;
			else{
			    printf("Cannot use both -f and -d\n\n");
			    usage();
			}
			break;
			case 'p':
				strcpy(devicepath,&argv[1][2]);
			break;
			case 'o':
				outputfile=1;
				strcpy(outputpath,&argv[1][2]);
				break;
			case 'b':
				str_split(&argv[1][2]);
				break;
			case 'z':
				dbfolder=1;
				strcpy(dbfolderpath,&argv[1][2]);
				break;
			case 'h':
				help();
				break;
			default:
				printf("Wrong Argument: %s\n", argv[1]);
				usage();
				exit(0);
				break;
			}

			++argv;
			--argc;
	}


	//test if iotype is defined
	if(devicetype!=FOLDER && devicetype!=DEVICE){
		printf("missing -f or -d\n\n");
		usage();
		exit(0);
	}
	//test if testype is defined
	if(strlen(devicepath)==0){
		printf("missing -p<value>\n\n");
		usage();
		exit(0);
	}
	//test if blocksize >0
	if(nr_sizes_proc==0){
		nr_sizes_proc=1;
		sizes_proc=malloc(sizeof(int));
		sizes_proc[0]=4096;
	}

	//Initialize variables
	total_blocks=malloc(sizeof(uint64_t)*nr_sizes_proc);
	eq=malloc(sizeof(uint64_t)*nr_sizes_proc);
	dif=malloc(sizeof(uint64_t)*nr_sizes_proc);
	distinctdup=malloc(sizeof(uint64_t)*nr_sizes_proc);
	//zeroed_blocks=malloc(sizeof(uint64_t)*nr_sizes_proc);
	space=malloc(sizeof(uint64_t)*nr_sizes_proc);

	dbporiginal=malloc(sizeof(DB**)*nr_sizes_proc);
	envporiginal=malloc(sizeof(DB_ENV**)*nr_sizes_proc);

	dbprinter=malloc(sizeof(DB**)*nr_sizes_proc);
	envprinter=malloc(sizeof(DB_ENV**)*nr_sizes_proc);


	int aux=0;
	for(aux=0;aux<nr_sizes_proc;aux++){

		dbporiginal[aux]=malloc(sizeof(DB *));
		envporiginal[aux]=malloc(sizeof(DB_ENV *));
		dbprinter[aux]=malloc(sizeof(DB *));
		envprinter[aux]=malloc(sizeof(DB_ENV *));

		char printdbpath[100];
		char duplicatedbpath[100];
		char sizeid[5];
		sprintf(sizeid,"%d",aux);


		//if a folder were specified for databases
		if(dbfolder==1){
			strcpy(printdbpath,PRINTDB);
			strcat(printdbpath,sizeid);
			strcpy(duplicatedbpath,dbfolderpath);
			strcat(duplicatedbpath,sizeid);
		}
		else{
			strcpy(printdbpath,PRINTDB);
			strcat(printdbpath,sizeid);
			strcpy(duplicatedbpath,DUPLICATEDB);
			strcat(duplicatedbpath,sizeid);
		}

		char mkcmd[200];
		sprintf(mkcmd, "mkdir -p %s", printdbpath);
		int ress = system(mkcmd);
		sprintf(mkcmd, "mkdir -p %s", duplicatedbpath);
		ress=system(mkcmd);
		if(ress<0)
	    	perror("Error creating folders for databases\n");


		printf("Removing old databases\n");
		//remove databases if exist
		remove_db(duplicatedbpath,dbporiginal[aux],envporiginal[aux]);
		remove_db(printdbpath,dbprinter[aux],envprinter[aux]);


		printf("Initing new database\n");
		init_db(duplicatedbpath,dbporiginal[aux],envporiginal[aux]);

		if(outputfile==1){
			init_db(printdbpath,dbprinter[aux],envprinter[aux]);
		}
	}

	//initialize analysis variables
	bzero(total_blocks,nr_sizes_proc*(sizeof(uint64_t)));
	//identical chunks (that could be eliminated)
	bzero(eq,nr_sizes_proc*(sizeof(uint64_t)));
	//distinct chunks
	bzero(dif,nr_sizes_proc*(sizeof(uint64_t)));
	//distinct chunks with duplicates
	bzero(distinctdup,nr_sizes_proc*(sizeof(uint64_t)));
	//chunks that were appended with zeros due to their size
	//bzero(zeroed_blocks,nr_sizes_proc*(sizeof(uint64_t)));
	//duplicated disk space
	bzero(space,nr_sizes_proc*(sizeof(uint64_t)));

		
	//check if it is a folder or device and start processing
	if(devicetype==FOLDER){
		printf("start processing folder %s\n",devicepath);
		search_dir(devicepath);
	}
	else{
		printf("start processing device %s\n",devicepath);
		extract_blocks(devicepath);
	}

	for(aux=0;aux<nr_sizes_proc;aux++){

		fprintf(stderr,"\n\n\nResults for %d\n",sizes_proc[aux]);
		fprintf(stderr,"files scanned %llu\n",(unsigned long long int)nfiles);
		fprintf(stderr,"total blocks scanned %llu\n",(unsigned long long int)total_blocks[aux]);
		//fprintf(stderr,"total blocks with zeros appended %llu\n",(unsigned long long int)zeroed_blocks[aux]);
		//blocks without any duplicate are the distinct block minus the distinct blocks with duplicates
		uint64_t zerodups=dif[aux]-distinctdup[aux];
		fprintf(stderr,"blocks without duplicates %llu\n",(unsigned long long int)zerodups);
		fprintf(stderr,"distinct blocks with duplicates %llu\n",(unsigned long long int)distinctdup[aux]);
		fprintf(stderr,"duplicate blocks %llu\n",(unsigned long long int)eq[aux]);
		fprintf(stderr,"space saved %llu Bytes\n",(unsigned long long int)space[aux]);

		//if outputdist was chosen and specified generate it
		if(outputfile==1){

			printf("before generating output dist\n");

			gen_output(dbporiginal[aux],envporiginal[aux],dbprinter[aux],envprinter[aux]);

			printf("before printing output dist\n");

			char outputfilename[100];
			char sizeid[10];

			sprintf(sizeid,"%d",sizes_proc[aux]);
			strcpy(outputfilename,outputpath);
			strcat(outputfilename,sizeid);
			FILE* fpp=fopen(outputfilename,"w");

			print_elements_print(dbprinter[aux], envprinter[aux],fpp);
			fclose(fpp);

			close_db(dbprinter[aux],envprinter[aux]);
			//TODO this is not removed now but in the future it can be...
			//remove_db(printdbpath,dbprinter[aux],envprinter[aux]);

		}	

		close_db(dbporiginal[aux],envporiginal[aux]);
		//TODO this is not removed to keep the database for dedisgen-utils
		//remove_db(duplicatedbpath,dbporiginal,envporiginal);
	}

	//free memory
	free(total_blocks);
	free(eq);
	free(dif);
	free(distinctdup);
	//free(zeroed_blocks);
	free(space);

	free(dbporiginal);
	free(envporiginal);
	free(dbprinter);
	free(envprinter);

return 0;

/*main is used for testing this code*/
/*
int main(int argc, char* argv[]){
	//int i=0;

	double prob_dist[MAX_SIZE];
	double cum_dist[MAX_SIZE];
	double index = .99;
	prob_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, 7.5);
	cum_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, cum_dist);
	int return_index = binary_search(index, cum_dist, MAX_SIZE);
	//test_check_prob_dist_output(cum_dist, MAX_SIZE);
	printf("%f is found at %d, which is between %f and %f\n", index, return_index, cum_dist[return_index-1], cum_dist[return_index]);
	for(int i=0; i< 801; i++){
		prob_dist_distribution_entropy_gen(prob_dist, MAX_SIZE, i/100.0);
		printf("Entropy requested: %f\n",i/100.0);
		//test_check_prob_dist_output(prob_dist, 5);
		printf("Entropy generated: %f\n\n", calculate_entropy_of_prob_dist(prob_dist, MAX_SIZE));
	}
	return 0;
}
*/
/*
int main(int argc, char **argv){

	FILE* pFile;
	struct define_entropy ds;
	unsigned int size = 256*1024;
	void* buffer_to_populate=malloc(size);
	ds.ent = 3.00;
	populating_buffer_to_populatefer_with_actual_entropy(&ds, buffer_to_populate, size);
	pFile = fopen("/home/rahulraju93/result33","w+");
    	if (pFile ){
        		fwrite(buffer_to_populate,1,sizeof(buffer_to_populate),pFile);
		   }
	fclose(pFile);
	printf("The entropy of buffer_to_populate: %f\n", test_check_buffer_to_populatefer_output(buffer_to_populate, size));
	free(buffer_to_populate);
	return 0;
}
*/


}






