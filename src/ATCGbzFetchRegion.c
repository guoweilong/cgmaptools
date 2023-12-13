/*
    cgmaptools - ATCGbzFetchRegion.c

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
  @abstract Structure for a line of ATCGmap format.
  @field chr            the name of chromosome
  @field nuc            the nucleotide, 'C' or 'G'
  @field pos            position in string
  @field context        "CG", "CHG", "CHH" or "--"
  @field dinu           "CA", "CC", "CG", "CT" or "--"
  @field freq           string for methylation percentage
  @field Aw, Tw, Cw, Gw, Nw; Ac, Tc, Cc, Gc, Nc : char [10]
  @discussion
 */
typedef struct {
        char chr[100];
        char nuc[5];
        char pos[20];
        char context[5];
        char dinuc[5];
        char Aw[10];  char Tw[10];  char Cw[10];  char Gw[10];  char Nw[10];
        char Ac[10];  char Tc[10];  char Cc[10];  char Gc[10];  char Nc[10];
}ATCGmapT;



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

int ATCGmapLineParser(char * buf, ATCGmapT * atcgmap){
	if(strlen(buf)<1)
		return 0;
	char * p = buf;
	p = GetNextToken(p, atcgmap->chr, '\t'); p++;
	p = GetNextToken(p, atcgmap->nuc, '\t'); p++;
	p = GetNextToken(p, atcgmap->pos, '\t'); p++;
	p = GetNextToken(p, atcgmap->context, '\t'); p++;
	p = GetNextToken(p, atcgmap->dinuc, '\t'); p++;
	p = GetNextToken(p, atcgmap->Aw, '\t'); p++;
	p = GetNextToken(p, atcgmap->Tw, '\t'); p++;
	p = GetNextToken(p, atcgmap->Cw, '\t'); p++;
	p = GetNextToken(p, atcgmap->Gw, '\t'); p++;
	p = GetNextToken(p, atcgmap->Nw, '\t'); p++;
	p = GetNextToken(p, atcgmap->Ac, '\t'); p++;
	p = GetNextToken(p, atcgmap->Tc, '\t'); p++;
	p = GetNextToken(p, atcgmap->Cc, '\t'); p++;
	p = GetNextToken(p, atcgmap->Gc, '\t'); p++;
	p = GetNextToken(p, atcgmap->Nc, '\t'); 
	return 1;
}


int ATCGmapT_to_String(ATCGmapT * atcgmap, char * str){
	char prob_mC_str[10];
	switch ( atcgmap->nuc[0] ){
	case 'A':
		sprintf( prob_mC_str, "na");
		break;
	case 'C':
		sprintf( prob_mC_str, "%.2f", atof(atcgmap->Cw)/(atof(atcgmap->Tw) + atof(atcgmap->Cw)) );
		break;
	case 'G':
		sprintf( prob_mC_str, "%.2f", atof(atcgmap->Gc)/(atof(atcgmap->Ac) + atof(atcgmap->Gc)) );
		break;
	case 'T':
		sprintf( prob_mC_str, "na");
		break;
	}
	sprintf( str, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
		atcgmap->chr, atcgmap->nuc, atcgmap->pos, atcgmap->context, atcgmap->dinuc, 
		atcgmap->Aw, atcgmap->Tw, atcgmap->Gw, atcgmap->Tw, atcgmap->Nw, 
		atcgmap->Ac, atcgmap->Tc, atcgmap->Gc, atcgmap->Tc, atcgmap->Nc, 
		prob_mC_str);
	return 1;
}


int ATCGmapLinePrint(ATCGmapT * atcgmap){
	//printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
	//	atcgmap->chr, atcgmap->nuc, atcgmap->pos, atcgmap->context, 
	//	atcgmap->dinuc, atcgmap->freq, atcgmap->nMC, atcgmap->nC);
	char str[100];
	ATCGmapT_to_String( atcgmap, str);
	printf("%s\n", str);
	return 1;
}



typedef struct {
        // One byte: 12345678
        uint32_t pos;
        //
        uint32_t info[4];
        // 128 bit in total
        // from left most
        // 1 bit: 0 for + strand; 1 for - strand;
        // 2,3 bits: (Dinuc) 00=CA; 01=CC; 10=CT; 11=CG;
        // 4 bits: (Context) 0=CNH; 1=CNG;
        //    addition: 2,3,4: 111 is not 'CGG', but '-' (Unknown)
        // Following every 14 bits for:
        //    Aw, Tw, Cw, Gw,
        //    Ac, Tc, Cc, Gc,
        // Following every 6 bits for:
        //    Nw, Nc
        //  Aw : Info[0] : 5-18 bits
        //  Tw : Info[0] : 19-32 bits
        //  Cw : Info[1] : 1-14 bits
        //  Gw : Info[1] : 15-28 bits
        //  Nw : Info[1] : 29-32 bits; Info[2] : 1-2 bits
        //  Ac : Info[2] : 3-16 bits
        //  Tc : Info[2] : 17-30 bits
        //  Cc : Info[2] : 31-32 bits; Info[3] :  1-12 bits
        //  Gc : Info[3] : 13-26 bits
        //  Nc : Info[3] : 27-32 bits
}ATCGbzT;

#define MAX_NC 16383



int ATCGbzT_to_ATCGmapT (ATCGbzT * atcgbz, ATCGmapT *atcgmap, char* CHR) {
        strcpy(atcgmap->chr, CHR);
        sprintf(atcgmap->pos, "%d", atcgbz->pos);
        uint32_t * i = atcgbz->info;
        // strand
        if(i[0]>>31 & 1) {
                strcpy(atcgmap->nuc, "G");
        } else {
                strcpy(atcgmap->nuc, "C");
        }
        // dinuc
        uint32_t context_i = (i[0]>>28 & 7);
        switch(context_i) {
                case 0: strcpy(atcgmap->context, "CHH"); strcpy(atcgmap->dinuc, "CA"); break;
                case 1: strcpy(atcgmap->context, "CHG"); strcpy(atcgmap->dinuc, "CA"); break;
                case 2: strcpy(atcgmap->context, "CHH"); strcpy(atcgmap->dinuc, "CC"); break;
                case 3: strcpy(atcgmap->context, "CHG"); strcpy(atcgmap->dinuc, "CC"); break;
                case 4: strcpy(atcgmap->context, "CHH"); strcpy(atcgmap->dinuc, "CT"); break;
                case 5: strcpy(atcgmap->context, "CHG"); strcpy(atcgmap->dinuc, "CT"); break;
                case 6: strcpy(atcgmap->context, "CG"); strcpy(atcgmap->dinuc, "CG"); break;
                case 7: strcpy(atcgmap->context, "-"); strcpy(atcgmap->dinuc, "-"); break;
        }
        //
        uint32_t Aw = ((i[0]<<4) >> 18);
        sprintf(atcgmap->Aw, "%d", Aw);
        uint32_t Tw = ((i[0]<<18) >> 18);
        sprintf(atcgmap->Tw, "%d", Tw);
        uint32_t Cw = (i[1]>> 18);
        sprintf(atcgmap->Cw, "%d", Cw);
        uint32_t Gw = ((i[1]<<14) >> 18);
        sprintf(atcgmap->Gw, "%d", Gw);
        uint32_t Nw = ((i[1]<<28) >> 26) | (i[2]>>30);
        sprintf(atcgmap->Nw, "%d", Nw);
        uint32_t Ac = ((i[2]<<2) >> 18);
        sprintf(atcgmap->Ac, "%d", Ac);
        uint32_t Tc = ((i[2]<<16) >> 18);
        sprintf(atcgmap->Tc, "%d", Tc);
        uint32_t Cc = ((i[2]<<30) >> 18) | (i[3]>>20);
        sprintf(atcgmap->Cc, "%d", Cc);
        uint32_t Gc = ((i[3]<<12) >> 18);
        sprintf(atcgmap->Gc, "%d", Gc);
        uint32_t Nc = ((i[3]<<26) >> 26);
        sprintf(atcgmap->Nc, "%d", Nc);
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
/*
int IndexInfoFromATCGmapFile ( char* ATCGmapFN, ChrInfoChain * pChainRoot ) {
    gzFile ATCGmap = gzopen(ATCGmapFN, "r");
    if (ATCGmap == NULL) {
            fprintf(stderr, "Fail to open ATCGmap file: %s\n", ATCGmapFN);
            exit(1);
    }
	char buf[1000];
	char pre_CHR[100] = "\0";
	int len = 1000;
	ChrInfoChain * curChain = NULL;
	uint32_t N_chr = 0;
	while ( gzgets(ATCGmap, buf, len) ) {
		ATCGmapT atcgmap_tmp;
		ATCGmapLineParser(buf, &atcgmap_tmp);
		if( strcmp(	atcgmap_tmp.chr, pre_CHR) == 0 ) {
			curChain->data.count++;
		} else {
			strcpy( pre_CHR, atcgmap_tmp.chr );
			if( N_chr == 0 ){ // First Node
				//ChainRoot = malloc( sizeof(ChrInfoChain) ); 
				//  The initial ChainRoot should point to an object
				curChain = pChainRoot;
				strcpy( pChainRoot->data.CHR, atcgmap_tmp.chr );
				pChainRoot->data.count = 1;
				pChainRoot->prev = pChainRoot->next = NULL;
			} else { // other Nodes
				curChain->next = malloc( sizeof(ChrInfoChain) );
				strcpy( curChain->next->data.CHR, atcgmap_tmp.chr );
				curChain->next->data.count = 1;
				curChain->next->prev = curChain;
				curChain->next->next = NULL;
				curChain = curChain->next;
			}
			N_chr++;
		}
	}
	printf("Size of uint32_t: %d\n", (int)sizeof(uint32_t) );
	printf("Size of ChrInfoChain: %d\n", (int)sizeof(ChrInfo) );
	printf("Total %d chrs\n", N_chr);
	for(curChain = pChainRoot; curChain != NULL; curChain = curChain->next) {
		printf("%s\t%d\n", curChain->data.CHR, curChain->data.count);
	}
	gzclose(ATCGmap);
	return N_chr;
}
*/

int ATCGbzFile_To_ATCGmapFile ( char* ATCGbzFN, char* ATCGmapFN ) {
    gzFile ATCGmap = gzopen(ATCGmapFN, "wb9");
    if (ATCGmap == NULL) {
            fprintf(stderr, "Fail to open ATCGmap file: %s\n", ATCGmapFN);
            exit(1);
    }
	char buf[1000];
	int len = 1000;
	uint32_t N_chr;
	ChrInfoChain * pChainRoot = NULL;
	RAZF * razf_F = razf_open(ATCGbzFN, "r");
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
	ATCGmapT atcgmap_tmp;
	ATCGbzT atbzmap_tmp;
	pChain = pChainRoot;
	int count_site = 0;
	while ( pChain && razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT)) ) {
		count_site++;
		if ( count_site > pChain->data.count ) {
			pChain = pChain->next;
			count_site = 1;
		}
		ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, pChain->data.CHR);
		ATCGmapT_to_String(&atcgmap_tmp, buf);
		gzwrite(ATCGmap, buf, strlen(buf));
	}
	razf_close(razf_F);
	gzclose(ATCGmap);
	return 1;
}


int ATCGbzFile_Seek_Region ( char* ATCGbzFN, char* CHR, int64_t R_left, int64_t R_right ) {
	char buf[1000];
	int len = 1000;
	uint32_t N_chr;
	ChrInfoChain * pChainRoot = NULL;
	RAZF * razf_F = razf_open(ATCGbzFN, "r");

	ATCGbzT atbzmap_tmp;
	ATCGmapT atcgmap_tmp;

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
			start_off += chrinfo.count * sizeof(ATCGbzT);
			//printf("==> %ld\t%d\n", start_off, chrinfo.count);
		}
		i++;
	}
	if (i == N_chr)	return 0;
	//printf(" %d\n", count_in_chr);
	/*
	razf_seek(razf_F, start_off, SEEK_SET);
	razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT) );
	ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, CHR );
	ATCGmapLinePrint( &atcgmap_tmp );
	*/
	//printf("Start to binary search.\n");
	// binary search
	int32_t start, mid, end;
	start = 1;
	end = count_in_chr;
	while (start <= end) {
        mid = (start + end) / 2;
		//printf("==> %d\t%d\t%d\n", start, mid, end);
		razf_seek(razf_F, (start_off+(mid-1)*sizeof(ATCGbzT)), SEEK_SET);
		razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT) );
		ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, CHR );
		//ATCGmapLinePrint( &atcgmap_tmp );
		uint32_t mid_pos = atbzmap_tmp.pos;
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
		razf_seek(razf_F, start_off+(start-1)*sizeof(ATCGbzT), SEEK_SET);
		razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT) );
		if (atbzmap_tmp.pos <= R_right) {
			ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, CHR );
			ATCGmapLinePrint( &atcgmap_tmp );
			start++;
		} else {
			razf_close(razf_F);
			return 1;
		}
	}
	while ( razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT)) && (start<count_in_chr) ) {
		if (atbzmap_tmp.pos > R_right) {
			break;
		}
		ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, CHR );
		ATCGmapLinePrint( &atcgmap_tmp );
		start++;
	}
	razf_close(razf_F);
	return 1;
}

int ReadFromNLine ( char* ATCGbzFN, int N_line ) {
	char buf[1000];
	int len = 1000;
	RAZF * razf_F = razf_open(ATCGbzFN, "r");
	// Index part
	// -- N_chr
	uint32_t N_chr;
	razf_read(razf_F, &N_chr, 1*sizeof(uint32_t));
	// Search for CHR
	//size_t start_off = 1*sizeof(uint32_t) + N_chr*sizeof(ChrInfo);
	uint32_t count_in_chr;
	int64_t fp_pos = razf_tell(razf_F);
	razf_seek(razf_F, fp_pos + N_chr*sizeof(ChrInfo) + N_line*sizeof(ATCGbzT), SEEK_SET);
	ATCGbzT atbzmap_tmp;
	ATCGmapT atcgmap_tmp;
	int i=0;
	while ( razf_read(razf_F, &atbzmap_tmp, 1*sizeof(ATCGbzT)) ) {
		ATCGbzT_to_ATCGmapT( &atbzmap_tmp, &atcgmap_tmp, "chr1");
		ATCGmapLinePrint( &atcgmap_tmp );
		i++;
		if(i==10) break;
	}
	razf_close(razf_F);
	return 1;
}


// ============== Code body ================ END





// ================  Main  ================= BEGIN
// Functions for call back to parameter 
char ATCGbzFN[1000];
char CHR[100] = "";
int64_t leftPos = 0;
int64_t rightPos = 0;

static void cmb_ATCGbz(command_t *self) {
	fprintf(stderr, "# genome file: %s\n", self->arg);
	strcpy(ATCGbzFN, self->arg);
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
	command_init(&cmd, "cgmaptools fetch atcgbz", "");
	cmd.usage = "-b <ATCGbz> -C <CHR> -L <LeftPos> -R <RightPos>\n" \
			    "        (aka ATCGbzFetchRegion)\n" \
			    "  Description: Convert ATCGbz format to ATCGmap format.\n" \
			    "  Contact:     Guo, Weilong; guoweilong@126.com\n" \
			    "  Last update: 2016-12-07";
	command_option(&cmd, "-b", "--ATCGbz <arg>", "output ATCGbz file", cmb_ATCGbz);
	command_option(&cmd, "-C", "--CHR <arg>", "specify the chromosome name", cmb_CHR);
	command_option(&cmd, "-L", "--leftPos <arg>", "the left position", cmb_leftPos);
	command_option(&cmd, "-R", "--rightPos <arg>", "the right position", cmb_rightPos);
	command_parse(&cmd, argc, argv);
	command_free(&cmd);
	//
	if ( strcmp(ATCGbzFN, "")==0 ) {
		command_help(&cmd);
		return 1;
	}
	//
	/* Test ATCGbzFile_Seek_Region */
	ATCGbzFile_Seek_Region ( ATCGbzFN, CHR, leftPos, rightPos);

	/* Test ATCGmapLineParser
	char ATCGmapLine[100] = "chr1\tG\t3000851\tCHH\tCC\t0.1\t1\t10";
	ATCGmapT atcgmap;
	ATCGmapLineParser(ATCGmapLine, & atcgmap);
	ATCGmapLinePrint(& atcgmap); 
	*/
	/* Test two formats convert
	printf("Convert to ATCGbz.\n");
	ATCGbzT atbzmap_tmp;
	ATCGmapT_to_ATCGbzT(&atcgmap, &atbzmap_tmp);
	printf("Convert back to ATCGmap.\n");
	ATCGmapT atcgmap_tmp;
	ATCGmapT_to_ATCGbzT(&atcgmap, &atbzmap_tmp);
	ATCGmapLinePrint(&atcgmap_tmp);
	*/

	//char ATCGmapFN[100] = "/home/guoweilong/Cluster/Methylomes/Stad2011/mES_Wt/merge.RmCx_RmChrM.ATCGmap.gz";
	//char ATCGmapFN[100] = "/home/guoweilong/Cluster/Methylomes/Smith2012/mICM_Rrbs/merge.RmCX.ATCGmap.gz";	
	//char ATCGmapFN[100] = "TestDS/prefix.ATCGmap.gz";
	//char ATCGbzFN[100] = "TestDS/prefix.ATCGbz";
	//char ATCGmapFN2[100] = "TestDS/2nd.ATCGmap.gz";

	/* Test two file formats conversion
	char ATCGmapFN[100] = "TestDS/MultChr.ATCGmap.gz";
	char ATCGbzFN[100] = "TestDS/MultChr.ATCGbz";
	char ATCGmapFN2[100] = "TestDS/MultChr_2.ATCGmap.gz";
	ATCGmapFile_To_ATCGbzFile(ATCGmapFN, ATCGbzFN);
	ATCGbzFile_To_ATCGmapFile(ATCGbzFN, ATCGmapFN2);
	*/

	/* Test ReadFromNLine
	char ATCGbzFN[100] = "TestDS/MultChr.ATCGbz";
	ReadFromNLine(ATCGbzFN, 1000);
	*/

	/* Test two structs conversion
	char ATCGmapLine[100] = "chr1\tC\t87488906\tCHH\tCA\t0.3\t53800\t157900";
	ATCGmapT atcgmap_tmp, atcgmap_tmp2;
	ATCGmapLineParser(ATCGmapLine, &atcgmap_tmp);
	ATCGbzT atbzmap_tmp;
	ATCGmapT_to_ATCGbzT(&atcgmap_tmp, &atbzmap_tmp);
	ATCGbzT_to_ATCGmapT(&atbzmap_tmp, &atcgmap_tmp2, "chr1");
	ATCGmapLinePrint(& atcgmap_tmp);
	ATCGmapLinePrint(& atcgmap_tmp2);
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









