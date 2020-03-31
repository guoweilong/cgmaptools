/*
    cgmaptools - ATCGmapMerge.cpp

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

/*
 Guo, Weilong; guoweilong@gmail.com; 2014-09-02
 */


#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
//#include <zlib.h>
/*If changed to supporting zlib, function pointer will be needed*/
using namespace std;
//#include "math.h"

struct parameter
{
	string fn1;	// -1
	string fn2;	// -2
};

parameter param;


void exit_with_help( void )
{
	printf(
		"Usage:	cgmaptools merge2 atcgmap -1 <ATCGmap> -2 <ATCGmap>\n"
		"       (aka ATCGmapMerge)\n"
		"Contact:     Guo, Weilong; guoweilong@126.com;\n"
		"Last Update: 2016-12-07\n"
        "Options:\n"
		"  -1	Input, 1st ATCGmap file\n"
		"  -2	Input, 2nd ATCGmap file\n"
        "Output to STDOUT in ATCGmap format\n"
        "Tips: Two input files should have the same order of chromosomes\n"
		);
	exit(0);
}

void parse_command_line(int argc, char **argv)
{
	int i;

	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
		case '1':
			if(i == argc)	exit_with_help();
			param.fn1 = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case '2':
			if(i == argc)	exit_with_help();	
			param.fn2 = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.fn1.length() && !param.fn2.length() ){
		exit_with_help();
	}
}


vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
inline vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
	// Skip delimiters at beginning.
	string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	vector<string> result;
	result.clear();

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		//__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");

		//if (pos == lastPos) result.push_back("");
		result.push_back(str.substr(lastPos, pos - lastPos));

		if (pos == string::npos) break;
		if (pos == str.length() - 1) {
			if (!skip_empty) result.push_back("");
			break;
		}
		// Skip delimiters.  Note the "not_of"
		lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
	return result;
}

void init(){
	// default values
	param.fn1 = "";
	param.fn2 = "";
}

//   0,   1,   2,       3,     4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,     15
// chr, nuc, pos, pattern, dinuc, WA, WT, WC, WG, WN, CA, CT, CC, CG, CN, methyl


class ATCGmapLine{
public:
	string	chr;
	char	nuc;
	long	pos;
	string	pattern;
	string	dinuc;
	unsigned int	WA;
	unsigned int	WT;
	unsigned int	WC;
	unsigned int	WG;
	unsigned int	WN;
	unsigned int	CA;
	unsigned int	CT;
	unsigned int	CC;
	unsigned int	CG;
	unsigned int	CN;
	//float	methyl;
	ATCGmapLine(){
		chr=string();
	};
	ATCGmapLine(vector<string> & tokens){
		chr = tokens[0];
		nuc = tokens[1][0];
		pos = atoi( tokens[2].c_str() );
		pattern = tokens[3];
		dinuc = tokens[4];
		WA = atoi( tokens[5].c_str() );
		WT = atoi( tokens[6].c_str() );
		WC = atoi( tokens[7].c_str() );
		WG = atoi( tokens[8].c_str() );
		WN = atoi( tokens[9].c_str() );
		CA = atoi( tokens[10].c_str() );
		CT = atoi( tokens[11].c_str() );
		CC = atoi( tokens[12].c_str() );
		CG = atoi( tokens[13].c_str() );
		CN = atoi( tokens[14].c_str() );
		//methyl = atof( tokens[15].c_str() );
	};

	ATCGmapLine(string CHR, char NUC, long POS, string PATTERN, string DINUC,
		    unsigned int wa, unsigned int wt, unsigned int wc, unsigned int wg, unsigned int wn, 
		    unsigned int ca, unsigned int ct, unsigned int cc, unsigned int cg, unsigned int cn):
		    chr(CHR), nuc(NUC), pos(POS), pattern(PATTERN), dinuc(DINUC), 
		   WA(wa), WT(wt), WC(wc), WG(wg), WN(wn), CA(ca), CT(ct), CC(cc), CG(cg), CN(cn) {};

	const ATCGmapLine operator+(ATCGmapLine &x){
		return ATCGmapLine( chr, nuc, pos, pattern, dinuc, WA+x.WA, WT+x.WT, WC+x.WC, WG+x.WG, WN+x.WN,
		                    CA+x.CA, CT+x.CT, CC+x.CC, CG+x.CG, CN+x.CN );
	};

	friend ostream& operator<< (ostream &os, const ATCGmapLine &x)  
	{
		if (!x.chr.empty()) {
			os << setprecision(2);
			os << x.chr << "\t" << x.nuc << "\t" << x.pos << "\t" << x.pattern << "\t" << x.dinuc << "\t" 
			   << x.WA << "\t" << x.WT << "\t" << x.WC << "\t" << x.WG << "\t" << x.WN << "\t" << x.CA << 
		           "\t" << x.CT << "\t" << x.CC << "\t" << x.CG << "\t" << x.CN << "\t";
			if (x.nuc == 'C' && x.WC+x.WT>0) {
				os << float(x.WC)/(x.WC+x.WT) << endl;  
			} else if (x.nuc == 'G' && x.CG+x.CA>0) {
				os << float(x.CG)/(x.CG+x.CA) << endl;  
			} else {
				os << "na" << endl;
			}
		}
		return os;  
	}  
};

/*
ATCGmapLine GetNextATCGmapElement(ifstream &IN) {
	char buffer[1000+1];

	IN.getline(buffer, 1000);
	if(!strlen(buffer))
		return ATCGmapLine();
	if (buffer[strlen(buffer)-1] == '\r') {
		buffer[strlen(buffer)-1] = '\0';
	}
	string tmp = buffer;
	vector<string> tokens = string_tokenize(tmp);
	return ATCGmapLine(tokens);
}
*/

#include <zlib.h>

int ATCGmapMerge (string fn1, string fn2) {
    gzFile IN_1 = gzopen(fn1.c_str(), "r");
    if (IN_1 == NULL) {
            fprintf(stderr, "# Fail to open ATCGmap file: %s\n", fn1.c_str());
            exit(1);
    }
    gzFile IN_2 = gzopen(fn2.c_str(), "r");
    if (IN_2 == NULL) {
            fprintf(stderr, "# Fail to open ATCGmap file: %s\n", fn2.c_str());
            exit(1);
    }
	int k = 0;
	string chr_pre = "";
	char buf_1[1000];
	char buf_2[1000];
	int len = 1000;
	char * EOF1 = gzgets(IN_1, buf_1, len);
	char * EOF2 = gzgets(IN_2, buf_2, len);
	ATCGmapLine Ele_1, Ele_2;
	while ( EOF1 && EOF2 ) {
	    vector<string> tokens_1 = string_tokenize(string(buf_1));
        Ele_1 = ATCGmapLine( tokens_1 );
	    vector<string> tokens_2 = string_tokenize(string(buf_2));
        Ele_2 = ATCGmapLine( tokens_2 );
        if(Ele_1.chr != Ele_2.chr){
			if (Ele_1.chr == chr_pre) {
				cout << Ele_1;
				EOF1 = gzgets(IN_1, buf_1, len);
			} else if (Ele_2.chr == chr_pre) {
				cout << Ele_2;
				EOF2 = gzgets(IN_2, buf_2, len);
			} else {
				cerr << "Debug\n";
			}
		} else {
			chr_pre= Ele_1.chr;
			if (Ele_1.pos < Ele_2.pos) {
				cout << Ele_1;
				EOF1 = gzgets(IN_1, buf_1, len);
			} else if (Ele_1.pos > Ele_2.pos) {
				cout << Ele_2;
				EOF2 = gzgets(IN_2, buf_2, len);
			} else {
				//ATCGmapLine Ele_3 = Ele_1 + Ele_2;
				cout << (Ele_1 + Ele_2);
				EOF1 = gzgets(IN_1, buf_1, len);
				EOF2 = gzgets(IN_2, buf_2, len);
			}
		}
	}
	if (EOF2) {
		cout << Ele_1;
		while ( (EOF1 = gzgets(IN_1, buf_1, len)) ) {
			cout << buf_1 << endl;
		}
	} else if (EOF1) {
		cout << Ele_2;
		while ( (EOF2 = gzgets(IN_2, buf_2, len)) ) {
			cout << buf_2 << endl;
		}
	}
	return 1;
}
/*
int ATCGmapMerge (string fn1, string fn2) {
	ifstream IN_1(fn1.c_str());
	if(!IN_1) {
		cout << "cannot open input file \"" << fn1.c_str() << "\"\n";
		exit(1);
	}
	ifstream IN_2(fn2.c_str());
	if(!IN_2) {
		cout << "cannot open input file \"" << fn2.c_str() << "\"\n";
		exit(1);
	}
	int k = 0;
	string chr_pre = "";
	ATCGmapLine Ele_1 = GetNextATCGmapElement(IN_1);
	ATCGmapLine Ele_2 = GetNextATCGmapElement(IN_2);
	while (!IN_1.eof() && !IN_2.eof()) {
		if(Ele_1.chr != Ele_2.chr){
			if (Ele_1.chr == chr_pre) {
				cout << Ele_1;
				Ele_1 = GetNextATCGmapElement(IN_1);
			} else if (Ele_2.chr == chr_pre) {
				cout << Ele_2;
				Ele_2 = GetNextATCGmapElement(IN_2);
			} else {
				cerr << "Debug\n";
			}
		} else {
			chr_pre= Ele_1.chr;
			if (Ele_1.pos < Ele_2.pos) {
				cout << Ele_1;
				Ele_1 = GetNextATCGmapElement(IN_1);
			} else if (Ele_1.pos > Ele_2.pos) {
				cout << Ele_2;
				Ele_2 = GetNextATCGmapElement(IN_2);
			} else {
				//ATCGmapLine Ele_3 = Ele_1 + Ele_2;
				cout << (Ele_1 + Ele_2);
				Ele_1 = GetNextATCGmapElement(IN_1);
				Ele_2 = GetNextATCGmapElement(IN_2);
			}
		}
		//cout << "test" << endl;
	}
	string line;
	if (IN_2.eof()) {
		cout << Ele_1;
		while ( getline(IN_1, line) ) {
			cout << line << endl;
		} 
	} else if (IN_1.eof()) {
		cout << Ele_2;
		while ( getline(IN_2, line) ) {
			cout << line << endl;
		} 
	}
	return 1;
}
*/


int main(int argc, char* argv[])
{
	init();

	parse_command_line(argc,argv);

	ATCGmapMerge( param.fn1, param.fn2);
    
	return 0;
}

// todo: error will raise if chr_list not complete in two files

