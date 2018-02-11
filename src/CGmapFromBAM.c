/*
    cgmaptools - CGmapFromBAM.c

    Copyright (C) Weilong Guo
    Contact: Weilong Guo <guoweilong@126.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>




int RmOverlap = 0;


// ======== Option Parser for C ========== BEGIN
// Code revised from https://github.com/clibs/commander
// Max options that can be defined.
#define COMMANDER_MAX_OPTIONS 128
// Max arguments that can be passed.
#define COMMANDER_MAX_ARGS 128
// Command struct.
struct command;
// Option callback.
typedef void (* command_callback_t)(struct command *self);
// Command option.
typedef struct {
  int optional_arg;
  int required_arg;
  char *argname;
  char *large;
  const char *small;
  const char *large_with_arg;
  const char *description;
  command_callback_t cb;
} command_option_t;
// Command.
typedef struct command {
  void *data;
  const char *usage;
  const char *arg;
  const char *name;
  //const char *version;
  int option_count;
  command_option_t options[COMMANDER_MAX_OPTIONS];
  int argc;
  char *argv[COMMANDER_MAX_ARGS];
  char **nargv;
} command_t;
// prototypes
//void command_init(command_t *self, const char *name, const char *version);
void command_init(command_t *self, const char *name);
void command_free(command_t *self);
void command_help(command_t *self);
void command_option(command_t *self, const char *small, const char *large, const char *desc, command_callback_t cb);
void command_parse(command_t *self, int argc, char **argv);
// include
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
// Output error and exit.
static void error(char *msg) {
  fprintf(stderr, "%s\n", msg);
  exit(1);
}
// Output command version.
//static void command_version(command_t *self) {
//  printf("%s\n", self->version);
//  command_free(self);
//  exit(0);
//}
// Output command help.
void command_help(command_t *self) {
  printf("\n");
  printf("  Usage: %s %s\n", self->name, self->usage);
  printf("\n");
  printf("  Options:\n");
  printf("\n");
  int i;
  for (i = 0; i < self->option_count; ++i) {
    command_option_t *option = &self->options[i];
    printf("    %s, %-18s %s\n"
      , option->small
      , option->large_with_arg
      , option->description);
  }
  printf("\n");
  command_free(self);
  exit(0);
}
// Initialize with program `name` and `version`.
//void command_init(command_t *self, const char *name, const char *version) {
void command_init(command_t *self, const char *name) {
  self->arg = NULL;
  self->name = name;
  //self->version = version;
  self->option_count = self->argc = 0;
  self->usage = "[options]";
  self->nargv = NULL;
  //command_option(self, "-V", "--version", "output program version", command_version);
  command_option(self, "-h", "--help", "output help information", command_help);
}
// Free up commander after use.
void command_free(command_t *self) {
  int i;
  for (i = 0; i < self->option_count; ++i) {
    command_option_t *option = &self->options[i];
    free(option->argname);
    free(option->large);
  }
  if (self->nargv) {
    for (i = 0; self->nargv[i]; ++i) {
      free(self->nargv[i]);
    }
    free(self->nargv);
  }
}
// Parse argname from `str`. For example
// Take "--required <arg>" and populate `flag`
// with "--required" and `arg` with "<arg>".
static void parse_argname(const char *str, char *flag, char *arg) {
  int buffer = 0;
  size_t flagpos = 0;
  size_t argpos = 0;
  size_t len = strlen(str);
  size_t i;
  for (i = 0; i < len; ++i) {
    if (buffer || '[' == str[i] || '<' == str[i]) {
      buffer = 1;
      arg[argpos++] = str[i];
    } else {
      if (' ' == str[i]) continue;
      flag[flagpos++] = str[i];
    }
  }
  arg[argpos] = '\0';
  flag[flagpos] = '\0';
}
// Normalize the argument vector by exploding
// multiple options (if any). For example
// "foo -abc --scm git" -> "foo -a -b -c --scm git"
static char ** normalize_args(int *argc, char **argv) {
  int size = 0;
  int alloc = *argc + 1;
  char **nargv = malloc(alloc * sizeof(char *));
  int i;
  for (i = 0; argv[i]; ++i) {
    const char *arg = argv[i];
    size_t len = strlen(arg);
    // short flag
    if (len > 2 && '-' == arg[0] && !strchr(arg + 1, '-')) {
      alloc += len - 2;
      nargv = realloc(nargv, alloc * sizeof(char *));
      size_t  j;
      for (j = 1; j < len; ++j) {
        nargv[size] = malloc(3);
        sprintf(nargv[size], "-%c", arg[j]);
        size++;
      }
      continue;
    }
    // regular arg
    nargv[size] = malloc(len + 1);
    strcpy(nargv[size], arg);
    size++;
  }
  nargv[size] = NULL;
  *argc = size;
  return nargv;
}
// Define an option.
void command_option(command_t *self, const char *small, const char *large, const char *desc, command_callback_t cb) {
  if (self->option_count == COMMANDER_MAX_OPTIONS) {
    command_free(self);
    error("Maximum option definitions exceeded");
  }
  int n = self->option_count++;
  command_option_t *option = &self->options[n];
  option->cb = cb;
  option->small = small;
  option->description = desc;
  option->required_arg = option->optional_arg = 0;
  option->large_with_arg = large;
  option->argname = malloc(strlen(large) + 1);
  assert(option->argname);
  option->large = malloc(strlen(large) + 1);
  assert(option->large);
  parse_argname(large, option->large, option->argname);
  if ('[' == option->argname[0]) option->optional_arg = 1;
  if ('<' == option->argname[0]) option->required_arg = 1;
}
// Parse `argv` (internal).
// Input arguments should be normalized first
// see `normalize_args`.
static void command_parse_args(command_t *self, int argc, char **argv) {
  int literal = 0;
  int i, j;
  for (i = 1; i < argc; ++i) {
    const char *arg = argv[i];
    for (j = 0; j < self->option_count; ++j) {
      command_option_t *option = &self->options[j];
      // match flag
      if (!strcmp(arg, option->small) || !strcmp(arg, option->large)) {
        self->arg = NULL;
        // required
        if (option->required_arg) {
          arg = argv[++i];
          if (!arg || '-' == arg[0]) {
            fprintf(stderr, "%s %s argument required\n", option->large, option->argname);
            command_free(self);
            exit(1);
          }
          self->arg = arg;
        }
        // optional
        if (option->optional_arg) {
          if (argv[i + 1] && '-' != argv[i + 1][0]) {
            self->arg = argv[++i];
          }
        }
        // invoke callback
        option->cb(self);
        goto match;
      }
    }
    // --
    if ('-' == arg[0] && '-' == arg[1] && 0 == arg[2]) {
      literal = 1;
      goto match;
    }
    // unrecognized
    if ('-' == arg[0] && !literal) {
      fprintf(stderr, "unrecognized flag %s\n", arg);
      command_free(self);
      exit(1);
    }
    int n = self->argc++;
    if (n == COMMANDER_MAX_ARGS) {
      command_free(self);
      error("Maximum number of arguments exceeded");
    }
    self->argv[n] = (char *) arg;
    match:;
  }
}
// Parse `argv` (public).
void command_parse(command_t *self, int argc, char **argv) {
  self->nargv = normalize_args(&argc, argv);
  command_parse_args(self, argc, self->nargv);
  self->argv[self->argc] = NULL;
}
// ======== Option Parser for C ============= END








// ===== Dictionary Implementary for C ====== BEGIN
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#define DICT_MALLOC malloc
#define DICT_FREE free
// dict_t pair struct.
typedef struct dict_pair {
  struct dict_pair *prev;
  struct dict_pair *next;
  char *key;
  void *val;
} dict_pair_t;
// dict_t struct.
typedef struct {
  dict_pair_t *head;
  dict_pair_t *tail;
  void (*free)(char *key, void *val);
} dict_t;
// dict_t iterator struct.
typedef struct {
  dict_pair_t *next;
} dict_iterator_t;
// Node prototypes.
dict_pair_t *
dict_pair_new(char *key, void *val);
// dict_t prototypes.
dict_t * dict_new();
dict_pair_t * dict_set(dict_t *self, char *key, void *val);
dict_pair_t * dict_get(dict_t *self, char *key);
void dict_remove(dict_t *self, char *key);
void dict_destroy(dict_t *self);
// dict_t iterator prototypes.
dict_iterator_t * dict_iterator_new(dict_t *dict);
dict_pair_t * dict_iterator_next();
void dict_iterator_destroy(dict_iterator_t *self);
// Allocate a new dict_iterator_t. NULL on failure.
dict_iterator_t * dict_iterator_new(dict_t *dict) {
  dict_pair_t *pair = dict->head;
  dict_iterator_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_iterator_t))))
    return NULL;
  self->next = pair;
  return self;
}
// Return the next dict_pair_t or NULL when no more
// nodes remain in the dict.
dict_pair_t * dict_iterator_next(dict_iterator_t *self) {
  dict_pair_t *curr = self->next;
  if (curr)
    self->next = curr->next;
  return curr;
}
// Free the dict iterator.
void dict_iterator_destroy(dict_iterator_t *self) {
  DICT_FREE(self);
  self = NULL;
}
// dict.c
// Allocates a new dict_pair_t. NULL on failure.
dict_pair_t * dict_pair_new(char *key, void *val) {
  dict_pair_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_pair_t))))
    return NULL;
  //printf("malloc pair\n");
  self->prev = NULL;
  self->next = NULL;
  self->key = key;
  self->val = val;
  return self;
}
// Allocate a new dict_t. NULL on failure.
dict_t * dict_new() {
  dict_t *self;
  if (!(self = DICT_MALLOC(sizeof(dict_t))))
    return NULL;
  //printf("malloc dict\n");
  self->head = NULL;
  self->tail = NULL;
  self->free = NULL;
  return self;
}
// Free the dict.
void dict_destroy(dict_t *self) {
  dict_pair_t *next;
  dict_pair_t *curr = self->head;
  while (curr) {
    next = curr->next;
    if (self->free) {
        self->free(curr->key, curr->val);
        //printf("free 100\n");
    }
    DICT_FREE(curr);
    //printf("free pair\n");
    curr = next;
  }
  DICT_FREE(self);
  //printf("free dict\n");
}
// Add a key-value pair to dict.
dict_pair_t * dict_set(dict_t *self, char *key, void *val) {
  dict_pair_t *pair = dict_get(self, key);
  if (pair) {
    if (self->free) self->free(pair->key, pair->val);
    /*if (self->free) {
		free(pair->key);
		free(pair->val);
	}*/
    pair->val = val;
  } else {
    pair = dict_pair_new(key, val);
    if (self->head) {
      pair->prev = self->tail;
      pair->next = NULL;
      self->tail->next = pair;
      self->tail = pair;
    } else {
      self->head = self->tail = pair;
      pair->prev = pair->next = NULL;
    }
  }
  return pair;
}
// Return the node associated to val or NULL.
dict_pair_t * dict_get(dict_t *self, char *key) {
  dict_iterator_t *it = dict_iterator_new(self);
  dict_pair_t *pair;
  while ((pair = dict_iterator_next(it))) {
    //printf("DEBUG | %s\n", pair->key);
    if (strcmp(key, pair->key) == 0) {
      dict_iterator_destroy(it);
      return pair;
    }
  }
  dict_iterator_destroy(it);
  return NULL;
}
// Remove the given node from the dict, freeing it and it's value.
void dict_remove(dict_t *self, char *key) {
  dict_pair_t *pair = dict_get(self, key);
  if (!pair) return;
  pair->prev
    ? (pair->prev->next = pair->next)
    : (self->head = pair->next);
  pair->next
    ? (pair->next->prev = pair->prev)
    : (self->tail = pair->prev);
  if (self->free) self->free(pair->key, pair->val);
  DICT_FREE(pair);
}
// ===== Dictionary Implementary for C ====== END











// ============== Code body ================ BEGIN
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>
#include <ctype.h>
//#include "../include/zlib-1.2.8/zlib.h"
#include "../include/samtools-0.1.18/sam.h"
#include "../include/samtools-0.1.18/faidx.h"

typedef struct {
	char *current_chrom_seq;
	int current_tid;
	int current_chrom_length;
    char *genome_filename;
	//bool RmOverlap;
	//faidx_t *genome_idx;
	samfile_t *in;
	gzFile *ATCGmap;
	gzFile *CGmap;
	gzFile *wiggle;
} infoholder_t;

#define MAX_PATH 1024

int is_dir(const char *path) {
    struct stat s;
    if ( stat(path, &s) == 0 ) {
        if( s.st_mode & S_IFDIR ) {
            return 1;
        } else if( s.st_mode & S_IFREG ) {
            return 0;
        } else {
	        fprintf(stderr, "%s is neither a file nor a folder\n", path);
		    exit(1);
        }
    } else {
        fprintf(stderr, "Cannot access the directory: %s\n", path);
        exit(1);
    }
}

// Loads the sequence of a chromosome from the reference genome file.
void load_chromosome(int tid, infoholder_t *info){
	char *chrom_name = info->in->header->target_name[tid];
	fprintf(stderr, "# Processing reads from: %s\n", chrom_name);
	if (info->current_chrom_seq != NULL) {
		free(info->current_chrom_seq);
	}
    char genome_filename[MAX_PATH] = {0};
    if (is_dir(info->genome_filename)) {
        sprintf(genome_filename,"%s/%s.fa", info->genome_filename, chrom_name);
    } else {
        strcpy(genome_filename, info->genome_filename);
    }
	faidx_t *genome_idx = fai_load(genome_filename);
	if (genome_idx == NULL) {
		fprintf(stderr, "Fail to open index for genome file: %s\n", genome_filename);
		exit(1);
	}
	info->current_chrom_seq = fai_fetch(genome_idx, chrom_name, &info->current_chrom_length);
	info->current_tid = tid;
	if (info->current_chrom_seq == NULL) {
		fprintf(stderr, "Failed to load chromosome: %s\n", chrom_name);
		exit(1);
	}
	fai_destroy(genome_idx);
}

#define A 0
#define T 1
#define C 2
#define G 3
#define N 4

#define MAX_LINE 1024

//BAM nucleotides are encoded as 4-bit integers following this order: =ACMGRSVTWYHKDBN
//						  =  A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
int NT_BAM_TO_IDX[16] = { N, A, C, N, G, N, N, N, T, N, N, N, N, N, N, N};

// Returns the sequence of the read pointed by "b" as an array of chars.
char *get_seq(bam1_t *b, char *buf){
	uint8_t *seq = bam1_seq(b);
	int i;
	//printf("DEBUG | seqlen: %d\n", b->core.l_qseq);
	for (i = 0; i < b->core.l_qseq; i++) {
		buf[i] = bam_nt16_rev_table[bam1_seqi(seq, i)];
	}
	buf[i] = 0;
	return buf;
}

// calculates the reverse complement of the base "n"
static inline char rc(char n){
	return n == 'A'?'T':
		  (n == 'T'?'A':
		  (n == 'C'?'G':
		  (n == 'G'?'C':'N')));
}

// checks the context and stores it char *nuc, char* context, char *subcontext.
static inline void context_calling(const char *chrom_seq, int chrom_seq_length, int pos,
					 		char *nucleotide, char **context, char *subcontext) {
	char nuc = toupper(chrom_seq[pos]);
	char nuc2;
	*nucleotide = nuc;
	*context = "--";
	subcontext[0]= '-';
	subcontext[1]= '-';
    //subcontext = "--"; // This sentence cause error
//	printf("DEBUG | Enter context_call\n");
    if ((pos + 2) < chrom_seq_length && (pos - 2) >= 0){
//		printf("DEBUG | nuc = %c\n", nuc);
        if (nuc == 'C'){
            nuc2 = toupper(chrom_seq[pos + 1]);
//			printf("DEBUG | nuc2 = %c\n", nuc2);
            //subcontext[0] = nuc;
//			printf("DEBUG | == 0 ==\n");
            subcontext[0] = 'C';
//			printf("DEBUG | == 1 ==\n");
            subcontext[1] = nuc2;
//			printf("DEBUG | == 2 ==\n");
            if(nuc2 == 'G') {
                *context = "CG";
            } else if (nuc2 == 'A' || nuc2 == 'C' || nuc2 == 'T') {
//				printf("DEBUG | pos=%d; and chrom_seq length = %d", pos, strlen(chrom_seq) );
            	char nuc3 = toupper(chrom_seq[pos + 2]);
//				printf("DEBUG | nuc3 = %c\n", nuc3);
            	if (nuc3 == 'G') {
            		*context = "CHG";
            	} else if (nuc3 == 'A' || nuc3 == 'C' || nuc3 =='T') {
            		*context = "CHH";
            	}
            }
		} else if (nuc == 'G') {
            nuc2 = toupper(chrom_seq[pos - 1]);
            //subcontext[0] = rc(nuc);
			subcontext[0] = 'C';
            subcontext[1] = rc(nuc2);
			if (nuc2 == 'C'){
				*context = "CG";
			} else if (nuc2 == 'A' || nuc2 == 'G' || nuc2 == 'T') {
				char nuc3 = toupper(chrom_seq[pos - 2]);
				if (nuc3 == 'C') {
					*context = "CHG";
				} else if (nuc3 == 'A' || nuc3 == 'G' || nuc3 =='T') {
					*context = "CHH";
				}
            }
		}
	}
}

// Get the clean qname of reads, by removing mate infor
// "READXXXXCCAC#1" => "READXXXXCCAC"
static int QnameClean (char * qname, char * qname_clean) {
	size_t qname_len = strlen(qname);
	if (qname_len >= 3) {
	    char ch = qname[qname_len-2];
//		printf("DEBUG | Char=%c\n", ch);
		if (ch == '.' || ch == '#' || ch == ':') {
			strncpy(qname_clean, qname, qname_len-2);
			qname_clean[qname_len-2]='\0';
		} else {
			strcpy(qname_clean, qname);
		}
		return 1;
	}
	strcpy(qname_clean, qname);
	return 1;
}

// Seems not works good for linux
// small space may not return to system in real time
void FreeDictEle(char * key, void* value) {
	//free(key);
	free(value);
}

// Define a list of string point
// A pool of 10000 string for string sequence names
char * MEM[1000];
char ** MEM_pointer;

char ** InitStringList ( char ** p ) {
	int i;
	char ** q = p;
	for(i=0; i<1000; i++) {
		(*q) = DICT_MALLOC(50*sizeof(char));
		q++;
	}
	MEM_pointer = p;
	return p;
}

char * GetNextStringSpace (void) {
	//printf("debug 0\n");
	//printf("%d\n", MEM_pointer-MEM);
	if(MEM_pointer >= (MEM+999)) {
	    //printf("debug 1\n");
	    MEM_pointer = MEM;
	} else {
	    //printf("debug 2\n");
	    MEM_pointer++;
	    //printf("debug 3\n");
	}
	//printf("%d\n", MEM_pointer-MEM);
	return (*MEM_pointer);
}

// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pileups, void *vinfo)
{
	infoholder_t *info = (infoholder_t *) vinfo;
	int i;
	const bam_pileup1_t *pl;
	bam1_t *aln;
	int fwd_counts[5] = {0};
	int rev_counts[5] = {0};
	int current_tid = pileups->b->core.tid;
	char *chrom_name = info->in->header->target_name[current_tid];
	//char * aln_qname = NULL;
	if (current_tid != info->current_tid) {
		load_chromosome(current_tid, info);
		gzprintf(info->wiggle, "variableStep chrom=%s\n", chrom_name);
	}
	dict_t * Dict = dict_new();
	char * aln_qname ; // A pointer to dynamic space to store reads name
	int nucleotide;
	//Dict->free = FreeDictEle; // Point to function for free space of dict element
	//printf("DEBUG | pileup_func: before <for>\n");
	for (i = 0; i < n; i ++) {
        // printf("DEBUG | == %d ==\n", i);
		pl = pileups + i;
		aln = pl->b;
		if (RmOverlap) {
			//printf("DEBUG 0\n");
			//aln_qname = DICT_MALLOC(50*sizeof(char)); // Fixed Error here, should use dynamic space
			aln_qname = GetNextStringSpace(); // Find the space for a string
			//printf("GetSpace\n");
			//printf("%s\n", bam1_qname(aln));
			QnameClean(bam1_qname(aln), aln_qname);
			//printf("%s\n", aln_qname);
			//printf("DEBUG 1\n");
			if ( dict_get(Dict, aln_qname) == NULL ) {
			    //printf("DEBUG | First: %s\n", aln_qname);
				dict_set(Dict, aln_qname, NULL);
				// Redudent region A(1)
				if (!pl->indel){
					nucleotide = NT_BAM_TO_IDX[bam1_seqi(bam1_seq(aln), pl->qpos)];
					if (bam1_strand(aln)){ // negative strand
						rev_counts[nucleotide] ++;
					} else {    // positive strand
						fwd_counts[nucleotide] ++;
					}
				}
			} //else {
			    //printf("DEBUG | Second: %s\n", aln_qname);
			//}
			//printf("DEBUG 2\n");
	        //DICT_FREE(aln_qname); // avoid using too much memory
		} else {
            // Redudent region A(2)
			//printf("DEBUG 3\n");
            if (!pl->indel){
                nucleotide = NT_BAM_TO_IDX[bam1_seqi(bam1_seq(aln), pl->qpos)];
                if (bam1_strand(aln)){ // negative strand
                    rev_counts[nucleotide] ++;
                } else {    // positive strand
                    fwd_counts[nucleotide] ++;
                }
            }
		}
	}
    //
	//printf("DEBUG | pileup_func: after <for>\n");
    /*
    dict_pair_t *curr = Dict->head;
   	dict_pair_t *next;
  	while (curr) {
        	next = curr->next;
        	//if (Dict->free) Dict->free(curr->key, curr->val);
        	//printf("%s\n", curr->key);
        	//DICT_FREE(curr->key);
		curr->key = NULL;
		DICT_FREE(curr->val);
        	DICT_FREE(curr);
        	curr = next;
    }
    DICT_FREE(Dict);
	Dict = NULL;*/
	//
	dict_destroy(Dict);
	char *context;
	char subcontext[3] = {0};
	char nuc;
	context_calling(info->current_chrom_seq, info->current_chrom_length, pos,
                    &nuc, &context, subcontext);
    //	printf("DEBUG | pileup_func: after <context_call>\n");
	double meth_level = 0;
	int meth_level_is_available = 0;
	int meth_cytosines = 0;
	int unmeth_cytosines = 0;
	if (nuc == 'C'){
		// plus strand: take the ratio of C's to T's from reads that come from the forward strand
		meth_cytosines = fwd_counts[C];
		unmeth_cytosines = fwd_counts[T];
	} else if (nuc == 'G'){
		// minus strand: take the ratio of G's to A's from reads that come from the reverse strand
		meth_cytosines = rev_counts[G];
		unmeth_cytosines = rev_counts[A];
	}
	char meth_level_string[10] = {0};
	if (meth_cytosines + unmeth_cytosines > 0) {
		meth_level = ((double)meth_cytosines)/(meth_cytosines + unmeth_cytosines);
		meth_level_is_available = 1;
		sprintf(meth_level_string, "%.2lf", meth_level);
	} else {
		sprintf(meth_level_string, "na");
	}
    // output ATCGmap
	gzprintf(info->ATCGmap,
			"%s\t%c\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
			chrom_name,
			nuc,
			pos + 1,
			context,
			subcontext,
			fwd_counts[A], fwd_counts[T], fwd_counts[C], fwd_counts[G], fwd_counts[N],
			rev_counts[A], rev_counts[T], rev_counts[C], rev_counts[G], rev_counts[N],
			meth_level_string);
	if (meth_level_is_available){
	    // output wiggle
		gzprintf(info->wiggle, "%d\t%s\n", pos + 1, meth_level_string);
        // output CGmap
		gzprintf(info->CGmap,
				"%s\t%c\t%d\t%s\t%s\t%s\t%d\t%d\n",
				chrom_name,
				nuc,
				pos + 1,
				context,
				subcontext,
				meth_level_string,
				meth_cytosines,
				unmeth_cytosines + meth_cytosines);
	}
//	printf("DEBUG | pileup_func END\n");
	return 0;
}

// set zlib's internal buffer to 16MB
#define ZLIB_BUFFER 0x1000000

// the main method in the library
int call_methylation(char *sorted_bam_filename,
                     char *genome_filename,
                     char *ATCGmap_filename,
                     char *CGmap_filename,
                     char *wiggle_filename)
					//, bool RmOverlap)
	{
	infoholder_t info;
	info.current_chrom_seq = NULL;
	info.current_chrom_length = -1;
	info.current_tid = -1;
	//info.RmOverlap = RmOverlap;
	info.in = samopen(sorted_bam_filename, "rb", 0);
	if (info.in == 0) {
		fprintf(stderr, "Fail to open BAM file %s\n", sorted_bam_filename);
		exit(1);
	}
    info.ATCGmap = gzopen(ATCGmap_filename, "wb");
    gzbuffer(info.ATCGmap, ZLIB_BUFFER);
	if (info.ATCGmap == NULL) {
		fprintf(stderr, "Fail to open ATCGmap file: %s\n", ATCGmap_filename);
		exit(1);
	}
	info.CGmap = gzopen(CGmap_filename, "wb");
	gzbuffer(info.CGmap, ZLIB_BUFFER);
	if (info.CGmap == NULL) {
		fprintf(stderr, "Fail to open CGmap file: %s\n", CGmap_filename);
		exit(1);
	}
	info.wiggle = gzopen(wiggle_filename, "wb");
	gzbuffer(info.wiggle, ZLIB_BUFFER);
    if (info.wiggle == NULL) {
		fprintf(stderr, "Fail to open wiggle file: %s\n", wiggle_filename);
		exit(1);
	}
    info.genome_filename = genome_filename;
//	printf("DEBUG | Before sampileup\n");
	sampileup(info.in, -1, pileup_func, &info);
//	printf("DEBUG | After sampileup\n");
	samclose(info.in);
	if (info.current_chrom_seq != NULL) {
		free(info.current_chrom_seq);
	}
	gzclose(info.ATCGmap);
	gzclose(info.CGmap);
	gzclose(info.wiggle);
	return 0;
}
// ============== Code body ================ END





// ================  Main  ================= BEGIN
// Functions for call back to parameter
char BamFileName[1000];
char GenomeFileName[1000] = "";
char OutputPrefix[1000];

static void cmd_bam(command_t *self) {
	fprintf(stderr, "# input bam file: %s\n", self->arg);
	strcpy(BamFileName, self->arg);
}
static void cmd_genome(command_t *self) {
	fprintf(stderr, "# source genome file: %s\n", self->arg);
	strcpy(GenomeFileName, self->arg);
}
static void cmd_output(command_t *self) {
	fprintf(stderr, "# prefix for output CGmap: %s\n", self->arg);
	strcpy(OutputPrefix, self->arg);
}
static void cmd_RmOverlap(command_t *self) {
	RmOverlap = 1;
	fprintf(stderr, "# [selected mode] Remove overlap\n");
}
// main
// TODO: rm-CCGG for RRBS
// TODO: rm-XS:i:1 tag
int main(int argc, char **argv){
	command_t cmd;
	//command_init(&cmd, argv[0], "0.0.1");
	command_init(&cmd, "cgmaptools convert bam2cgmap");
	cmd.usage = "-b <BAM> -g <genome.fa> -o <prefix>\n"
	        "        (aka CGmapFromBAM)\n" \
			"  Description: Convert BAM file to CGmap/ATCGmap format.\n" \
			"  Notice: \n" \
			"    * For BS-Seeker2, CGmapTools seamlessly comparable with its BAM.\n" \
			"    * For Bismark, CGmapTools comparable with single end mode; \n" \
			"    *     and comparable with v0.8.2 and older version of bismarks for paired-end.\n" \
			"    *     Bismark changed FLAG strategy since v0.8.3, which raise unexpected count of nucleotides.\n" \
			"    * For BSmap, we can not guarantee generating right CGmap from our experience.\n" \
			"    * For more information, please contact us. Contribution to improve comparability is welcome!\n" \
			"  Contact: Guo, Weilong; guoweilong@126.com\n" \
			"  Last update: 2017-12-13";
	command_option(&cmd, "-b", "--bam <arg>", "input bam file, should be sorted first", cmd_bam);
	command_option(&cmd, "-g", "--genome <arg>", "genome file, fasta", cmd_genome);
	command_option(&cmd, "-O", "--rmOverlap", "Removed overlapped region for paired-end library if specified.", cmd_RmOverlap);
	command_option(&cmd, "-o", "--output [arg]", "prefix for output files", cmd_output);
	command_parse(&cmd, argc, argv);
	/*printf("Additional args:\n");
	int i;
	for (i = 0; i < cmd.argc; ++i) {
	printf("  - '%s'\n", cmd.argv[i]);
	}*/
	command_free(&cmd);
	/*int GWL[5] = {1,2,3};
	for (i = 0; i< 5; i++) {
	printf("%d\n", GWL[i]);
	}*/
	//
	/*
	char tmp_qname[100] = "SEQIDNNNNCAD1235#2";
	char tmp_qnameclean[100];
	QnameClean(tmp_qname, tmp_qnameclean);
	printf("Original: %s\tCleaned: %s\n", tmp_qname, tmp_qnameclean );
	*/
	// ===========================
	if (GenomeFileName[0] == 0) {
	    fprintf(stderr, "[Error] genome file is not specified.\n");
	    exit(1);
	}

	//printf("%s\t%s\t%s\n", BamFileName, GenomeFileName, OutputPrefix);
	char ATCGmap_filename[1000], CGmap_filename[1000], wiggle_filename[1000];
	sprintf(ATCGmap_filename, "%s.ATCGmap.gz", OutputPrefix);
	sprintf(CGmap_filename,   "%s.CGmap.gz",   OutputPrefix);
	sprintf(wiggle_filename,  "%s.wig.gz",     OutputPrefix);
	InitStringList(MEM);
	call_methylation(BamFileName, GenomeFileName, ATCGmap_filename, CGmap_filename, wiggle_filename);

	return 0;
}
// ================  Main  ================= END


// todo: fix the following error
/*
$ ../bin/CGmapFromBAM
open: No such file or directory
Fail to open BAM file
*/






