/*
    cgmaptools - CGmapSelectByRegion.cpp

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


/* CGmapSelectByRegion.cpp
 * Fuction: Find all the sites in CGmap file are or not covered by any region in the region.bed file
 * Author: GUO Weilong @ 2015-08-17
 * *********************
 * Modification:
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;


struct parameter
{
	string infile;	    // -i
	string regionfile;	// -r
	bool exclude;		// -d
};

parameter param;

class Region
{
public:
	unsigned int left;
	unsigned int right;
	Region () { left = 0; right = 0; }
	Region (unsigned int l, unsigned int r):left(l),right(r){
	}
};

class Site {
public:
	string 	chr;
	unsigned int	pos;
	char 	strand;
	string	pattern;
	unsigned int MC;
	unsigned int ALLC;
};


bool operator<(const Region &x, const Region &y)
{
	if (x.left != y.left) 
		return x.left < y.left;
	if (x.right != y.right)
		return x.right < y.right;
      // The function should return under all the conditions
      return false;
}

// For the same site in different strand, it's considered to be the same
bool operator==(const Region &x, const Region &y)
{
        return (x.left == y.left && x.right == y.right);
}

void exit_with_help( void )
{
	printf(
		"Usage:	cgmaptools select region [-i <CGmap/ATCGmap>] -r <BED> [-R]\n"
		"      (aka CGmapSelectByRegion)\n"
		"Description: Lines in input CGmap/ATCGmap be selected/excluded by BED file.\n"
		"             Strand is NOT considered.\n"
		"             Output to STDOUT in same format with input.\n"
		"Contact:     Guo, Weilong; guoweilong@126.com\n"
		"Last Update: 2016-12-07\n"
		"Options:\n"
		"  -i  Input, CGmap/ATCGmap file; use STDIN if not specified\n"
		"      Please use \"gunzip -c <input>.gz \" and pipe as input for gzipped file.\n"
		"      Ex: chr12\tG\t19898796\t...\n"
		"  -r  Input, Region file, BED file to store regions\n"
		"      At least 3 columns are required\n"
		"      Ex: chr12 19898766 19898966 XX XXX XXX\n"
		"  -R  [optional] Reverse selection. Sites in region file will be excluded when specified\n"
		"  -h  help\n"
		"Tips: program will do binary search for each site in regions\n"
		);
	exit(1);
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper); 
}

void parse_command_line(int argc, char **argv)
{
	int i;
	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
		case 'i':
			if(i == argc)	exit_with_help();
			param.infile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'r':
			if(i == argc)	exit_with_help();	
			param.regionfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'R':
			param.exclude = 1;
			break;
		case 'h':
			exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if ( !param.regionfile.length() ){
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
	param.infile = "";
	param.regionfile = "";
	param.exclude = 0;
}

int ReadRegionFile (map<string, vector<Region> > & Genome) {
	ifstream Rfile(param.regionfile.c_str());
	if(!Rfile) {
		cout << "cannot open input file \"" << param.regionfile.c_str() << "\"\n";
		exit(1);
	}
	int k = 0;
	while (!Rfile.eof()) {
		char buffer[1000+1];
		Rfile.getline(buffer, 1000);
		if(!strlen(buffer))
			continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		string chr = tokens[0];
		map< string,vector<Region> >::iterator giter;
		if ( (giter = Genome.find(chr)) == Genome.end() ) {
			vector<Region> trv;
			Genome[chr] = trv;
		}
		unsigned int left = atoi(tokens[1].c_str());
		unsigned int right = atoi(tokens[2].c_str());
		Region region(left, right);
		Genome[chr].push_back(region);
	}
	map<string, vector<Region> >::iterator giter;
	for ( giter = Genome.begin(); giter != Genome.end(); giter++ ) {
		string chr = giter->first;
		sort(Genome[chr].begin(), Genome[chr].end());
            unique(Genome[chr].begin(), Genome[chr].end());
	}
	Rfile.close();
	return 0;
}


#define MAX_LINE_LENGTH 10000

int CGmapSelectByRegion (map<string, vector<Region> > & Genome) {
	//ifstream CGmapF(param.infile.c_str());
	//chr1    G   3000851 CHH CC  0.1 1   10
	//chr1    C   3001624 CHG CA  0.0 0   9
	//chr1    C   3001631 CG  CG  1.0 5   5
	istream *pCGmapF = &cin;
	ifstream CGmapF;
	if (param.infile != "") {
		CGmapF.open(param.infile.c_str());
		if(!CGmapF) {
			cout << "cannot open input file" << param.infile.c_str() << endl;
			return -1;
		}
		pCGmapF = &CGmapF;
	}
	char buffer[MAX_LINE_LENGTH+1];
	while ( pCGmapF->getline(buffer, MAX_LINE_LENGTH) ) {
		if(!strlen(buffer))
			continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		Site site;
		site.chr = tokens[0]; //site.chr = string("chr") + tokens[0];
		string chr = site.chr;
		map< string,vector<Region> >::iterator giter;
		if ( (giter = Genome.find(chr)) == Genome.end() )
			continue;

		site.pos = atoi(tokens[2].c_str());
		size_t size = Genome[chr].size();
		if (size == 0)    continue;

		/*Binary search*/
		int start = 0;
		int end = (int)size - 1;
		int mid;
		while (start <= end) {
			mid = (start + end) / 2;
			if (site.pos < Genome[chr][mid].left){
				end = mid - 1;
			} else if (site.pos > Genome[chr][mid].right) {
				      start = mid + 1;
			} else {
				break;
			}
		}
		unsigned int left = Genome[chr][mid].left;
		unsigned int right = Genome[chr][mid].right;

		if (site.pos >= left && site.pos <= right) {
		    if (param.exclude == 0){
		        cout << tmp << endl;
		    }
		} else {
		    if (param.exclude == 1){
			    cout << tmp << endl;
			}
		}
	}
	
	CGmapF.close();
	return 0;
}

int main(int argc, char* argv[])
{
	init();

	parse_command_line(argc,argv);

	map< string, vector<Region> > Genome; 
	//store all regions in each chromosome

	//cout << "Before Read Region file" <<endl;
	ReadRegionFile(Genome);
	//cout << "After Read Region file" << endl;
	CGmapSelectByRegion(Genome);

	return 1;
}


