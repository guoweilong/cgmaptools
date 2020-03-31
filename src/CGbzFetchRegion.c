/*
    cgmaptools - CGbzFetchRegion.c

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
  const char *version;
  int option_count;
  command_option_t options[COMMANDER_MAX_OPTIONS];
  int argc;
  char *argv[COMMANDER_MAX_ARGS];
  char **nargv;
} command_t;
// prototypes
void command_init(command_t *self, const char *name, const char *version);
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
static void command_version(command_t *self) {
  printf("%s\n", self->version);
  command_free(self);
  exit(0);
}
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
  //command_free(self);
  exit(0);
}
// Initialize with program `name` and `version`.
void command_init(command_t *self, const char *name, const char *version) {
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








// ============== Code body ================ BEGIN
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>
#include <ctype.h>

#include "../include/samtools-0.1.18/razf.h"

#define MAX_PATH 1024

/*
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
        fprintf(stderr, "Cannot access: %s\n", path);
        exit(1);
    }
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


// calculates the reverse complement of the base "n"
inline char rc(char n){
	return n == 'A'?'T':
		  (n == 'T'?'A':
		  (n == 'C'?'G':
		  (n == 'G'?'C':'N')));
}

// set zlib's internal buffer to 16MB
#define ZLIB_BUFFER 0x1000000
*/


/*! @typedef
  @abstract Structure for a line of CGmap format.
  @field chr		the name of chromosome
  @field nuc		the nucleotide, 'C' or 'G'
  @field pos		position in string
  @field context	"CG", "CHG", "CHH" or "--"
  @field dinu		"CA", "CC", "CG", "CT" or "--"
  @field freq		string for methylation percentage
  @field nMC		string for number of methylated sytosimes
  @field nC			string for number of all cytosines

  @discussion		
 */
typedef struct {
	char chr[100];
	char nuc[5];
	char pos[20];
	char context[5];
	char dinuc[5];
	char freq[30];
	char nMC[10];
	char nC[10];
}CGmapT;


/*!
  @abstract		Get the next token from a string
  @param  p		SAM file handler
  @return		pointer to start of next token
 */
char * GetNextToken(char * p, char * dest, char delim){
	char * q = p;
	int i = 0;
	for(q = p; (*q) && (*q!=delim); q++, i++) ;
	strncpy(dest, p, i);
	dest[i] = '\0';
	return q;
}

int CGmapLineParser(char * buf, CGmapT * cgmap){
	if(strlen(buf)<1)
		return 0;
	char * p = buf;
	p = GetNextToken(p, cgmap->chr, '\t'); p++;
	p = GetNextToken(p, cgmap->nuc, '\t'); p++;
	p = GetNextToken(p, cgmap->pos, '\t'); p++;
	p = GetNextToken(p, cgmap->context, '\t'); p++;
	p = GetNextToken(p, cgmap->dinuc, '\t'); p++;
	p = GetNextToken(p, cgmap->freq, '\t'); p++;
	p = GetNextToken(p, cgmap->nMC, '\t'); p++;
	p = GetNextToken(p, cgmap->nC, '\t');
	return 1;
}


int CGmapT_to_String(CGmapT * cgmap, char * str){
	sprintf( str, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
		cgmap->chr, cgmap->nuc, cgmap->pos, cgmap->context, 
		cgmap->dinuc, cgmap->freq, cgmap->nMC, cgmap->nC);
	return 1;
}


int CGmapLinePrint(CGmapT * cgmap){
	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
		cgmap->chr, cgmap->nuc, cgmap->pos, cgmap->context, 
		cgmap->dinuc, cgmap->freq, cgmap->nMC, cgmap->nC);
	return 1;
}


typedef struct {
	// One byte: 8764321
	uint32_t pos;
	//
	uint32_t info;
	// 32 bits; from left most
	// 1 bit: 0 for + strand; 1 for - strand;
	// 2,3 bits: (Dinuc) 00=CA; 01=CC; 10=CT; 11=CG;
	// 4 bits: (Context) 0=CNH; 1=CNG;
	//    addition: 2,3,4: 111 is not 'CGG', but '-' (Unknown)
	// 5-16 bits (12 bits) unsigned short mMC;
	// 17-24 bits (12 bits) unsigned short NC;
	//char info;
	//unsigned short nMC;
	//unsigned short nC;
}CGbzT;

#define MAX_NC 16383

int CGmapT_to_CGbzT (CGmapT *cgmap, CGbzT * cgbz) {
	cgbz->pos = (uint32_t) atoi(cgmap->pos);
	uint32_t i = 0;
	// strand
	if ( cgmap->nuc[0] == 'G' ) i = 1;
	// context
	uint32_t context_i = 0;
	if ( cgmap->context[0] == '-' || cgmap->dinuc[0] == '-' ) {
		context_i = 7;
	} else {
		switch(cgmap->dinuc[1]){
		case 'A':
			if(cgmap->context[2]=='G') {
				context_i = 1; break; // CAG 001
			} else {
				context_i = 0; break; // CAH 000
			}
		case 'C':
			if(cgmap->context[2]=='G') {
				context_i = 3; break; // CCG 011
			} else {
				context_i = 2; break; // CCH 010
			}
		case 'T':
			if(cgmap->context[2]=='G') {
				context_i = 5; break; // CTG 101
			} else {
				context_i = 4; break; // CTH 100
			}
		case 'G':
			context_i = 6; break; // CG 110
		default :
			context_i = 7; break;	// '-', 'CN', other 111
		}
	}
	i = (i<<3 | context_i);
	//
	uint32_t nMC = (uint32_t) atoi(cgmap->nMC);
	uint32_t nC = (uint32_t) atoi(cgmap->nC);
	//printf("%d\t%d\n", nMC, nC);
	//
	while (nC > MAX_NC) {
		nC/=2;
		nMC/=2;
	}
	i = (i<<14 | nMC);
	i = (i<<14 | nC);
	//
	cgbz->info = i;
	//
	return 1;
}

int CGbzT_to_CGmapT (CGbzT * cgbz, CGmapT *cgmap, char* CHR) {
	strcpy(cgmap->chr, CHR);
	sprintf(cgmap->pos, "%d", cgbz->pos);
	uint32_t i = cgbz->info;
	// strand
	if(i>>31 & 1) {
		strcpy(cgmap->nuc, "G");
	} else {
		strcpy(cgmap->nuc, "C");
	}
	// dinuc
	uint32_t context_i = (i>>28 & 7);
	switch(context_i) {
		case 0: strcpy(cgmap->context, "CHH"); strcpy(cgmap->dinuc, "CA"); break;
		case 1: strcpy(cgmap->context, "CHG"); strcpy(cgmap->dinuc, "CA"); break;
		case 2: strcpy(cgmap->context, "CHH"); strcpy(cgmap->dinuc, "CC"); break;
		case 3: strcpy(cgmap->context, "CHG"); strcpy(cgmap->dinuc, "CC"); break;
		case 4: strcpy(cgmap->context, "CHH"); strcpy(cgmap->dinuc, "CT"); break;
		case 5: strcpy(cgmap->context, "CHG"); strcpy(cgmap->dinuc, "CT"); break;
		case 6: strcpy(cgmap->context, "CG"); strcpy(cgmap->dinuc, "CG"); break;
		case 7: strcpy(cgmap->context, "-"); strcpy(cgmap->dinuc, "-"); break;
	}
	//
	uint32_t nMC = ((i<<4) >> 18);
	sprintf(cgmap->nMC, "%d", nMC);
	uint32_t nC = ((i<<18) >> 18);
	sprintf(cgmap->nC, "%d", nC);
	sprintf(cgmap->freq, "%.2f", (double)nMC/nC);
	return 1;
}




typedef struct ChrInfo
{
	char CHR[118];
	uint32_t count;
} ChrInfo;

typedef struct ChrInfoChain
{
	ChrInfo data;
	struct ChrInfoChain *prev, *next;
} ChrInfoChain;


// Return N_chr : count of chromosomes
//
int IndexInfoFromCGmapFile ( char* CGmapFN, ChrInfoChain * pChainRoot ) {
    gzFile CGmap = gzopen(CGmapFN, "r");
    if (CGmap == NULL) {
            fprintf(stderr, "Fail to open CGmap file: %s\n", CGmapFN);
            exit(1);
    }
	char buf[1000];
	char pre_CHR[100] = "\0";
	int len = 1000;
	ChrInfoChain * curChain = NULL;
	uint32_t N_chr = 0;
	while ( gzgets(CGmap, buf, len) ) {
		CGmapT cgmap_tmp;
		CGmapLineParser(buf, &cgmap_tmp);
		if( strcmp(	cgmap_tmp.chr, pre_CHR) == 0 ) {
			curChain->data.count++;
		} else {
			strcpy( pre_CHR, cgmap_tmp.chr );
			if( N_chr == 0 ){ // First Node
				//ChainRoot = malloc( sizeof(ChrInfoChain) ); 
				//  The initial ChainRoot should point to an object
				curChain = pChainRoot;
				strcpy( pChainRoot->data.CHR, cgmap_tmp.chr );
				pChainRoot->data.count = 1;
				pChainRoot->prev = pChainRoot->next = NULL;
			} else { // other Nodes
				curChain->next = malloc( sizeof(ChrInfoChain) );
				strcpy( curChain->next->data.CHR, cgmap_tmp.chr );
				curChain->next->data.count = 1;
				curChain->next->prev = curChain;
				curChain->next->next = NULL;
				curChain = curChain->next;
			}
			N_chr++;
		}
	}
	fprintf(stderr, "# Size of uint32_t: %d\n", (int)sizeof(uint32_t) );
	fprintf(stderr, "# Size of ChrInfoChain: %d\n", (int)sizeof(ChrInfo) );
	fprintf(stderr, "# Total %d chrs\n", N_chr);
	for(curChain = pChainRoot; curChain != NULL; curChain = curChain->next) {
		fprintf(stderr, "# %s\t%d\n", curChain->data.CHR, curChain->data.count);
	}
	gzclose(CGmap);
	return N_chr;
}

int CGmapFile_To_CGbzFile ( char* CGmapFN, char* CGbzFN ) {
    gzFile CGmap = gzopen(CGmapFN, "r");
    if (CGmap == NULL) {
            fprintf(stderr, "Fail to open CGmap file: %s\n", CGmapFN);
            exit(1);
    }
	char buf[1000];
	int len = 1000;
	// Index part
	// -- N_chr
	ChrInfoChain ChainRoot;
	uint32_t N_chr;
	N_chr = IndexInfoFromCGmapFile(CGmapFN, &ChainRoot);
	//
	RAZF * razf_F = razf_open(CGbzFN, "w");
	// -- chromosome informaiton
	razf_write( razf_F, &N_chr, 1*sizeof(uint32_t) );
	ChrInfoChain * pChain;
	for (pChain = &ChainRoot; pChain!=NULL; pChain=pChain->next) {
		razf_write( razf_F, &(pChain->data), 1*sizeof(ChrInfo) );
	}
	// Content part
	while ( gzgets(CGmap, buf, len) ) {
		CGmapT cgmap_tmp;
		CGmapLineParser( buf, &cgmap_tmp );
		CGbzT cgbz_tmp;
		CGmapT_to_CGbzT( &cgmap_tmp, &cgbz_tmp );
		razf_write( razf_F, &cgbz_tmp, 1*sizeof(CGbzT) );
	}
	razf_close(razf_F);
	gzclose(CGmap);
	return 1;
}


int CGbzFile_To_CGmapFile ( char* CGbzFN, char* CGmapFN ) {
    gzFile CGmap = gzopen(CGmapFN, "wb9");
    if (CGmap == NULL) {
            fprintf(stderr, "Fail to open CGmap file: %s\n", CGmapFN);
            exit(1);
    }
	char buf[1000];
	int len = 1000;
	uint32_t N_chr;
	ChrInfoChain * pChainRoot = NULL;
	RAZF * razf_F = razf_open(CGbzFN, "r");
	// Index part
	// -- N_chr
	razf_read(razf_F, &N_chr, 1*sizeof(uint32_t));
	int i;
	pChainRoot = malloc( sizeof(ChrInfoChain) ); 
	ChrInfoChain * pChain;
	pChain = pChainRoot;
	razf_read(razf_F, &(pChain->data), 1*sizeof(ChrInfo));
	pChain->prev = NULL; pChain->next = NULL;
	for (i=1; i<N_chr; i++, pChain = pChain->next) {
		pChain->next = malloc( sizeof(ChrInfoChain) ); 
		razf_read(razf_F, &(pChain->next->data), 1*sizeof(ChrInfo));
		pChain->next->prev = pChain; pChain->next->next = NULL;
	}
	//
	CGmapT cgmap_tmp;
	CGbzT cgbz_tmp;
	pChain = pChainRoot;
	int count_site = 0;
	while ( pChain && razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT)) ) {
		count_site++;
		if ( count_site > pChain->data.count ) {
			pChain = pChain->next;
			count_site = 1;
		}
		CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, pChain->data.CHR);
		CGmapT_to_String(&cgmap_tmp, buf);
		gzwrite(CGmap, buf, strlen(buf));
	}
	razf_close(razf_F);
	gzclose(CGmap);
	return 1;
}


int CGbzFile_Seek_Region ( char* CGbzFN, char* CHR, int64_t R_left, int64_t R_right ) {
	char buf[1000];
	int len = 1000;
	uint32_t N_chr;
	ChrInfoChain * pChainRoot = NULL;
	RAZF * razf_F = razf_open(CGbzFN, "r");

	CGbzT cgbz_tmp;
	CGmapT cgmap_tmp;

	// Index part
	// -- N_chr
	razf_read(razf_F, &N_chr, 1*sizeof(uint32_t));
	// Search for CHR
	int i=0;
	size_t start_off = 1*sizeof(uint32_t) + N_chr*sizeof(ChrInfo);
	uint32_t count_in_chr;
	while ( i<N_chr ) {
		ChrInfo chrinfo;
		razf_read(razf_F, &chrinfo, 1*sizeof(ChrInfo));
		//printf("%s\t%s\n", chrinfo.CHR, CHR);
		if ( strcmp(chrinfo.CHR, CHR)==0 ) { // matched
			count_in_chr = chrinfo.count;
			break;
		} else { // not match
			start_off += chrinfo.count * sizeof(CGbzT);
			//printf("==> %ld\t%d\n", start_off, chrinfo.count);
		}
		i++;
	}
	if (i == N_chr)	return 0;
	//printf(" %d\n", count_in_chr);
	/*
	razf_seek(razf_F, start_off, SEEK_SET);
	razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT) );
	CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, CHR );
	CGmapLinePrint( &cgmap_tmp );
	*/
	//printf("Start to binary search.\n");
	// binary search
	int32_t start, mid, end;
	start = 1;
	end = count_in_chr;
	while (start <= end) {
        mid = (start + end) / 2;
		//printf("==> %d\t%d\t%d\n", start, mid, end);
		razf_seek(razf_F, (start_off+(mid-1)*sizeof(CGbzT)), SEEK_SET);
		razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT) );
		CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, CHR );
		//CGmapLinePrint( &cgmap_tmp );
		uint32_t mid_pos = cgbz_tmp.pos;
		//printf("%d\t%d\t%d\t%d:%ld-%ld\n", start, mid, end, mid_pos, R_left, R_right);
        if (mid_pos >= R_left){
            end = mid - 1;
        } else if (mid_pos < R_right) {
            start = mid + 1;
        } else {
            break;
        }
	}
	//printf("%ld\t%ld\t%ld\n", start, mid, end);
	//printf("Start to output sites.\n");
	// Output sites
	if(start < count_in_chr ) {
		razf_seek(razf_F, start_off+(start-1)*sizeof(CGbzT), SEEK_SET);
		razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT) );
		if (cgbz_tmp.pos <= R_right) {
			CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, CHR );
			CGmapLinePrint( &cgmap_tmp );
			start++;
		} else {
			razf_close(razf_F);
			return 1;
		}
	}
	while ( razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT)) && (start<count_in_chr) ) {
		if (cgbz_tmp.pos > R_right) {
			break;
		}
		CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, CHR );
		CGmapLinePrint( &cgmap_tmp );
		start++;
	}
	razf_close(razf_F);
	return 1;
}

int ReadFromNLine ( char* CGbzFN, int N_line ) {
	char buf[1000];
	int len = 1000;
	RAZF * razf_F = razf_open(CGbzFN, "r");
	// Index part
	// -- N_chr
	uint32_t N_chr;
	razf_read(razf_F, &N_chr, 1*sizeof(uint32_t));
	// Search for CHR
	//size_t start_off = 1*sizeof(uint32_t) + N_chr*sizeof(ChrInfo);
	uint32_t count_in_chr;
	int64_t fp_pos = razf_tell(razf_F);
	razf_seek(razf_F, fp_pos + N_chr*sizeof(ChrInfo) + N_line*sizeof(CGbzT), SEEK_SET);
	CGbzT cgbz_tmp;
	CGmapT cgmap_tmp;
	int i=0;
	while ( razf_read(razf_F, &cgbz_tmp, 1*sizeof(CGbzT)) ) {
		CGbzT_to_CGmapT( &cgbz_tmp, &cgmap_tmp, "chr1");
		CGmapLinePrint( &cgmap_tmp );
		i++;
		if(i==10) break;
	}
	razf_close(razf_F);
	return 1;
}


// ============== Code body ================ END





// ================  Main  ================= BEGIN
// Functions for call back to parameter 
char CGbzFN[1000];
char CHR[100] = "";
int64_t leftPos = 0;
int64_t rightPos = 0;

static void cmb_CGbz(command_t *self) {
	fprintf(stderr, "# genome file: %s\n", self->arg);
	strcpy(CGbzFN, self->arg);
}
static void cmb_CHR(command_t *self) {
	fprintf(stderr, "# chromosome name: %s\n", self->arg);
	strcpy(CHR, self->arg);
}
static void cmb_leftPos(command_t *self) {
	fprintf(stderr, "# start position: %s\n", self->arg);
	leftPos = atoi(self->arg);
}
static void cmb_rightPos(command_t *self) {
	fprintf(stderr, "# end position: %s\n", self->arg);
	rightPos = atoi(self->arg);
}

// main
int main(int argc, char **argv){
	command_t cmd;
	command_init(&cmd, "cgmaptools fetch cgbz", "");
	cmd.usage = "-b <CGbz> -C <CHR> -L <LeftPos> -R <RightPos>\n" \
			"         (aka CGvzFetchRegion)\n" \
			"  Description: Convert CGbz file to CGmap format.\n" \
			"  Contact: Guo, Weilong; guoweilong@126.com\n" \
			"  Last update: 2016-12-07";
	command_option(&cmd, "-b", "--CGbz <arg>", "output CGbz file", cmb_CGbz);
	command_option(&cmd, "-C", "--CHR <arg>", "specify the chromosome name", cmb_CHR);
	command_option(&cmd, "-L", "--leftPos <arg>", "the left position", cmb_leftPos);
	command_option(&cmd, "-R", "--rightPos <arg>", "the right position", cmb_rightPos);
	command_parse(&cmd, argc, argv);
	command_free(&cmd);
	//
	if ( strcmp(CGbzFN, "")==0 ) {
		command_help(&cmd);
		return 0;
	}
	//
	/* Test CGbzFile_Seek_Region */
	CGbzFile_Seek_Region ( CGbzFN, CHR, leftPos, rightPos);

	/* Test CGmapLineParser
	char CGmapLine[100] = "chr1\tG\t3000851\tCHH\tCC\t0.1\t1\t10";
	CGmapT cgmap;
	CGmapLineParser(CGmapLine, & cgmap);
	CGmapLinePrint(& cgmap); 
	*/
	/* Test two formats convert
	printf("Convert to CGbz.\n");
	CGbzT cgbz_tmp;
	CGmapT_to_CGbzT(&cgmap, &cgbz_tmp);
	printf("Convert back to CGmap.\n");
	CGmapT cgmap_tmp;
	CGmapT_to_CGbzT(&cgmap, &cgbz_tmp);
	CGmapLinePrint(&cgmap_tmp);
	*/

	//char CGmapFN[100] = "/home/guoweilong/Cluster/Methylomes/Stad2011/mES_Wt/merge.RmCx_RmChrM.CGmap.gz";
	//char CGmapFN[100] = "/home/guoweilong/Cluster/Methylomes/Smith2012/mICM_Rrbs/merge.RmCX.CGmap.gz";	
	//char CGmapFN[100] = "TestDS/prefix.CGmap.gz";
	//char CGbzFN[100] = "TestDS/prefix.CGbz";
	//char CGmapFN2[100] = "TestDS/2nd.CGmap.gz";

	/* Test two file formats conversion
	char CGmapFN[100] = "TestDS/MultChr.CGmap.gz";
	char CGbzFN[100] = "TestDS/MultChr.CGbz";
	char CGmapFN2[100] = "TestDS/MultChr_2.CGmap.gz";
	CGmapFile_To_CGbzFile(CGmapFN, CGbzFN);
	CGbzFile_To_CGmapFile(CGbzFN, CGmapFN2);
	*/

	/* Test ReadFromNLine
	char CGbzFN[100] = "TestDS/MultChr.CGbz";
	ReadFromNLine(CGbzFN, 1000);
	*/

	/* Test two structs conversion
	char CGmapLine[100] = "chr1\tC\t87488906\tCHH\tCA\t0.3\t53800\t157900";
	CGmapT cgmap_tmp, cgmap_tmp2;
	CGmapLineParser(CGmapLine, &cgmap_tmp);
	CGbzT cgbz_tmp;
	CGmapT_to_CGbzT(&cgmap_tmp, &cgbz_tmp);
	CGbzT_to_CGmapT(&cgbz_tmp, &cgmap_tmp2, "chr1");
	CGmapLinePrint(& cgmap_tmp);
	CGmapLinePrint(& cgmap_tmp2);
	*/	

	/* Size of types
	printf("sizeof int      : %ld\n", sizeof(int );
	printf("sizeof short    : %ld\n", sizeof(short);
	printf("sizeof size_t   : %ld\n", sizeof(size_t) );
	printf("sizeof int32_t  : %ld\n", sizeof(int32_t) );
	printf("sizeof int64_t  : %ld\n", sizeof(int64_t) );
	printf("sizeof uint32_t : %ld\n", sizeof(uint32_t) );
	printf("sizeof uint64_t : %ld\n", sizeof(uint64_t) );
	*/

	return 0;
}
// ================  Main  ================= END









