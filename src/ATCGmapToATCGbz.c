/*
    cgmaptools - ATCGmapToATCGbz.c

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
  //command_free(self);
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








// ============== Code body ================ BEGIN
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>
#include <ctype.h>

#include "../include/samtools-0.1.18/razf.h"

/*! @typedef
  @abstract Structure for a line of ATCGmap format.
  @field chr		the name of chromosome
  @field nuc		the nucleotide, 'C' or 'G'
  @field pos		position in string
  @field context	"CG", "CHG", "CHH" or "--"
  @field dinu		"CA", "CC", "CG", "CT" or "--"
  @field freq		string for methylation percentage
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
	char mC[10];
	float Tw, Cw, Ac, Gc;
	switch (atcgmap->nuc[0] ) {
	case 'C' :
		Tw = atof(atcgmap->Tw);
		Cw = atof(atcgmap->Cw);
		if ( (Tw+Cw)>0 ) { 
			sprintf(mC, "%.2f", Cw/(Cw+Tw) );
		} else {
			sprintf( mC, "na");
		}
		break;
	case 'G' :
		Ac = atof(atcgmap->Ac);
		Gc = atof(atcgmap->Gc);
		if ( (Ac+Gc)>0 ) { 
			sprintf(mC, "%.2f", Gc/(Ac+Gc) );
		} else {
			sprintf( mC, "na");
		}
		break;
	default :
		sprintf( mC, "na");
	}	
	sprintf( str, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
		atcgmap->chr, atcgmap->nuc, atcgmap->pos, atcgmap->context, 
		atcgmap->dinuc, 
		atcgmap->Aw, atcgmap->Tw, atcgmap->Cw, atcgmap->Gw, atcgmap->Nw,
		atcgmap->Ac, atcgmap->Tc, atcgmap->Cc, atcgmap->Gc, atcgmap->Nc,
		mC);
	return 1;
}


int ATCGmapLinePrint(ATCGmapT * atcgmap){
	char str[200];
	ATCGmapT_to_String(atcgmap, str);
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

/*
chr1  G 2154  CHG CA  0 0 0 1 0 0 0 0 0 0 na
chr1  T 2155  --  --  0 1 0 0 0 0 0 0 0 0 na
chr1  A 2156  --  --  1 0 0 0 0 0 0 0 0 0 na
chr1  C 2157  CHH CC  0 1 0 0 0 0 0 0 0 0 0.0
chr1  C 2158  CHH CC  0 1 0 0 0 0 0 0 0 0 0.0
chr1  C 2159  CHH CA  0 1 0 0 0 0 0 0 0 0 0.0
chr1  A 2160  --  --  1 0 0 0 0 0 0 0 0 0 na
*/

int ATCGmapT_to_ATCGbzT (ATCGmapT *cgmap, ATCGbzT * cgbz) {
	cgbz->pos = (uint32_t) atoi(cgmap->pos);
	uint32_t * i = cgbz->info;
	i[0] = 0;
	// strand
	if ( cgmap->nuc[0] == 'G' ) i[0] = 1;
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
	i[0] = (i[0]<<3 | context_i);
	//
	uint32_t Aw = (uint32_t) atoi(cgmap->Aw);
	uint32_t Tw = (uint32_t) atoi(cgmap->Tw);
	uint32_t Cw = (uint32_t) atoi(cgmap->Cw);
	uint32_t Gw = (uint32_t) atoi(cgmap->Gw);
	uint32_t Nw = (uint32_t) atoi(cgmap->Nw);
	//
	uint32_t Ac = (uint32_t) atoi(cgmap->Ac);
	uint32_t Tc = (uint32_t) atoi(cgmap->Tc);
	uint32_t Cc = (uint32_t) atoi(cgmap->Cc);
	uint32_t Gc = (uint32_t) atoi(cgmap->Gc);
	uint32_t Nc = (uint32_t) atoi(cgmap->Nc);
	//
	while (	Aw>MAX_NC || Tw>MAX_NC || Cw>MAX_NC || Gw>MAX_NC || Nw>MAX_NC
	     || Ac>MAX_NC || Tc>MAX_NC || Cc>MAX_NC || Gc>MAX_NC || Nc>MAX_NC ) {
		Aw/=2; Tw/=2; Cw/=2; Gw/=2; Nw/=2;
		Ac/=2; Tc/=2; Cc/=2; Gc/=2; Nc/=2;
	}
	if (Nw>=64){Nw=63;}
	if (Nc>=64){Nc=63;}
	i[0] = (i[0]<<14 | ((Aw<<18)>>18) );
	i[0] = (i[0]<<14 | ((Tw<<18)>>18) );
	i[1] = Cw;
	i[1] = (i[1]<<14 | ((Gw<<18)>>18) );
	i[1] = (i[1]<<4 | ((Nw<<26)>>28) );
	i[2] = Nw;
	i[2] = (i[2]<<14 | ((Ac<<18)>>18) );
	i[2] = (i[2]<<14 | ((Tc<<18)>>18) );
	i[2] = (i[2]<<2 | ((Cc<<18)>>30) );
	i[3] = Cc;
	i[3] = (i[3]<<14 | ((Gc<<18)>>18) );
	i[3] = (i[3]<<6 | ((Nc<<26)>>26) );
	//
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
		if( strcmp( atcgmap_tmp.chr, pre_CHR) == 0 ) {
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
	fprintf(stderr, "# Size of uint32_t: %d\n", (int)sizeof(uint32_t) );
	fprintf(stderr, "# Size of ChrInfoChain: %d\n", (int)sizeof(ChrInfo) );
	fprintf(stderr, "# Total %d chrs\n", N_chr);
	for(curChain = pChainRoot; curChain != NULL; curChain = curChain->next) {
		fprintf(stderr, "# %s\t%d\n", curChain->data.CHR, curChain->data.count);
	}
	gzclose(ATCGmap);
	return N_chr;
}

int ATCGmapFile_To_ATCGbzFile ( char* ATCGmapFN, char* ATCGbzFN ) {
    gzFile ATCGmap = gzopen(ATCGmapFN, "r");
    if (ATCGmap == NULL) {
            fprintf(stderr, "Fail to open ATCGmap file: %s\n", ATCGmapFN);
            exit(1);
    }
	char buf[1000];
	int len = 1000;
	// Index part
	// -- N_chr
	ChrInfoChain ChainRoot;
	uint32_t N_chr;
	N_chr = IndexInfoFromATCGmapFile(ATCGmapFN, &ChainRoot);
	//
	RAZF * razf_F = razf_open(ATCGbzFN, "w");
	// -- chromosome informaiton
	razf_write( razf_F, &N_chr, 1*sizeof(uint32_t) );
	ChrInfoChain * pChain;
	for (pChain = &ChainRoot; pChain!=NULL; pChain=pChain->next) {
		razf_write( razf_F, &(pChain->data), 1*sizeof(ChrInfo) );
	}
	// Content part
	while ( gzgets(ATCGmap, buf, len) ) {
		ATCGmapT atcgmap_tmp;
		ATCGmapLineParser( buf, &atcgmap_tmp );
		ATCGbzT atcgbz_tmp;
		ATCGmapT_to_ATCGbzT( &atcgmap_tmp, &atcgbz_tmp );
		razf_write( razf_F, &atcgbz_tmp, 1*sizeof(ATCGbzT) );
	}
	razf_close(razf_F);
	gzclose(ATCGmap);
	return 1;
}

// ============== Code body ================ END





// ================  Main  ================= BEGIN
// Functions for call back to parameter 
char ATCGmapFN[1000] = "";
char ATCGbzFN[1000];

static void cmd_ATCGmap(command_t *self) {
	fprintf(stderr, "# input ATCGmap file: %s\n", self->arg);
	strcpy(ATCGmapFN, self->arg);
}
static void cmb_ATCGbz(command_t *self) {
	fprintf(stderr, "# output ATCGbz file: %s\n", self->arg);
	strcpy(ATCGbzFN, self->arg);
}

// main
int main(int argc, char **argv){
	command_t cmd;
	//command_init(&cmd, argv[0], "0.0.1");
	command_init(&cmd, "cgmaptools convert atcgmap2atcgbz");
	cmd.usage = "-c <ATCGmap> -b <ATCGbz>\n" \
			"         (aka ATCGmapToATCGbz)\n" \
			"  Description: Convert ATCGmap format to ATCGbz format.\n" \
			"  Contact: Guo, Weilong; guoweilong@126.com\n" \
			"  Last update: 2016-12-07";
	command_option(&cmd, "-c", "--ATCGmap <arg>", "ATCGmap file (gzipped)", cmd_ATCGmap);
	command_option(&cmd, "-b", "--ATCGbz <arg>", "ATCGbz file", cmb_ATCGbz);
	command_parse(&cmd, argc, argv);
	command_free(&cmd);
	//
	if ( strcmp(ATCGmapFN, "")==0 ) {
		command_help(&cmd);
		return 0;
	}
	//
	ATCGmapFile_To_ATCGbzFile(ATCGmapFN, ATCGbzFN);
	//
	return 0;
}
// ================  Main  ================= END









