/*
    cgmaptools - CGmapMethInBed.cpp

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

/* CGmapMethylInBed.cpp 2014-01-03
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;


struct parameter
{
	string	CGmapfile;
	string  BEDfile;
	bool	by_strand;
	int	min_coverage;
	int	max_coverage;
};

parameter param;

class Site {
public:
	//string	chr;
	unsigned int	pos;
	char	strand;
	float	ML; // Methylation level
	//string context;
};


void exit_with_help( void )
{
	printf(
		"Usage:	cgmaptools mbed [-i <CGmap>]  -b <regin.bed> [-c 5 -C 500 -s] \n"
		"      (aka CGmapMethylInBed)\n"
		"Description: Calculated bulk average methylation levels in given regions.\n"
		"Contact:     Guo, Weilong; guoweilong@126.com\n"
		"Last Update: 2017-01-20\n"
		"Options:\n"
		"   -i  String, CGmap file; use STDIN if not specified\n"
		"       Please use \"gunzip -c <input>.gz \" and pipe as input for gzipped file.\n"
		"       Ex: chr1\tG\t3000851\tCHH\tCC\t0.1\t1\t10\n"
		"   -b  String, BED file, should have at least 4 columns\n"
		"       Ex: chr1\t3000000\t3005000\t-\n"
		"   -c  Int, minimum Coverage [Default: 5] \n"
		"   -C  Int, maximum Coverage [Default: 500] \n"
		"   -s  Strands would be distinguished when specified\n"
		"   -h  help\n\n"
		"Output to STDOUT:\n"
		"    Title         Count    mean_mC\n"
		"    sense         34       0.2353\n"
		"    antisense     54       0.2778\n"
		"    total         88       0.2614\n"
		"Notice:\n"
        "    The overlapping of regions would not be checked. \n"
		"    A site might be considered multiple times.\n\n"
		);
	exit(0);
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
			param.CGmapfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'b':
			if(i == argc)	exit_with_help();	
			param.BEDfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 's':
			//if(i == argc)	exit_with_help();	
			param.by_strand = true;
			//if(++i>argc)	exit_with_help();
			break;
		case 'c':
			if(i == argc)	exit_with_help();	
			param.min_coverage = atoi(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'C':
			if(i == argc)	exit_with_help();	
			param.max_coverage = atoi(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'h':
			exit_with_help();
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}
	if (  !param.BEDfile.length() ){
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
	param.CGmapfile = "";
	param.BEDfile = "";
	param.min_coverage = 5;
	param.max_coverage = 500;
	param.by_strand = false;
}


#define MAX_LINE_LENGTH 10000

int ReadMethylomeFromCGmap ( map< string, vector<Site> > & Genome ) {
	//
	//chr1    G   3000851 CHH CC  0.1 1   10
	//chr1    C   3001624 CHG CA  0.0 0   9
	//chr1    C   3001631 CG  CG  1.0 5   5
	//
	//ifstream CGmapF(param.CGmapfile.c_str());
	istream *pCGmapF = &cin;
	ifstream CGmapF;
	if (param.CGmapfile != "") {
		CGmapF.open(param.CGmapfile.c_str());
		if(!CGmapF) {
			cout << "cannot open input file" << param.CGmapfile.c_str() << endl;
			return -1;
		}
		pCGmapF = &CGmapF;
	}
	//CGmapF.reset(new ifstream(param.CGmapfile.c_str()));
	//CGmapF.open(param.CGmapfile.c_str());
	//while ( !CGmapF.eof() ) {
	//	CGmapF.getline(buffer, MAX_LINE_LENGTH);

	char buffer[MAX_LINE_LENGTH+1];
	//string tmp;
	while ( pCGmapF->getline(buffer, MAX_LINE_LENGTH) ) {
		if(!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		Site site;
		string chr = tokens[0];
		char nuc = tokens[1][0];
		if ( nuc == 'C' ) {
			site.strand = '+';
		} else {
			site.strand = '-';
		}
		site.pos = atoi ( tokens[2].c_str() );
		//site.context = tokens[3];
		site.ML = atof ( tokens[6].c_str() ) / atof (tokens[7].c_str() ) ;
		if ( Genome.find(chr) == Genome.end() ) {
			vector<Site> tsv;
			Genome[chr] = tsv;
		}
		int coverage = atoi(tokens[7].c_str());
		if ( (param.min_coverage <= coverage) && (coverage <=param.max_coverage) ) {
			Genome[chr].push_back(site);
		}
	}
	CGmapF.close();
	return 1;
}

int OutputMethylome( map< string, vector<Site> > & Genome, string fn ) {
	ofstream of( fn.c_str() );
	of.precision(2);
	for ( map< string, vector<Site> >::iterator iter=Genome.begin(); iter!=Genome.end(); iter++ ) {
        	string chr = iter->first;
		for ( vector<Site>::iterator citer=Genome[chr].begin(); citer!=Genome[chr].end(); citer++ ) {
			of << chr << "\t" << citer->pos << "\t" << citer->strand << "\t" << citer->ML << endl;		
		}
	}
	of.close();
	return 1; 
}


struct MethylRegion {
	// Watson Strand
	double TotalMethyl_W;
	int TotalSite_W;
	// Crick Strand
	double TotalMethyl_C;
	int TotalSite_C;

} ;

// Find all the methylated sites in the region [Left, Right],
// and put them in the lists of WStrand (Watson) and CStrand (Click)
int RegionMethylation ( vector<Site> & Chr, unsigned int Leftpos, unsigned int Rightpos, MethylRegion & mr ) {
	int low, high, mid;
	low = 0; 
	if (Chr.size() == 0) return 0;
	high = Chr.size() - 1;
	mid = ( low + high ) / 2;
	while ( low <= high ) {
		mid = ( low + high ) / 2;
		if ( Leftpos  > Chr[mid].pos ) { 
		// Not >=, so that can get low as the lowest one which >= Left
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	int chr_size = Chr.size();
	for ( unsigned int i = low; i<chr_size && Chr[i].pos <= Rightpos; i++ )
		//region.push_back( Chr[i] );
		if ( Chr[i].strand == '+' ) {
			mr.TotalMethyl_W +=Chr[i].ML;
			mr.TotalSite_W ++;
		} else {
			mr.TotalMethyl_C +=Chr[i].ML;
			mr.TotalSite_C ++;
		}

	return 1;
}



struct MethylSummary {
	double 	TotalMethyl_sense;
	double 	TotalMethyl_anti;
	unsigned int 	TotalSite_sense;
	unsigned int 	TotalSite_anti;
} ;


int AverageMethylInBins ( map< string, vector<Site> > & Methylome ) {

	ifstream BedFBedF(param.BEDfile.c_str());
	if(!BedFBedF) {
		cout << "cannot open input file" << param.BEDfile.c_str() << endl;
		return -1;
	}
	//double 	TotalMethyl_sense = 0;
	//double 	TotalMethyl_anti = 0;
	//unsigned int 	TotalSite_sense = 0;
	//unsigned int 	TotalSite_anti = 0;
	map<string, MethylSummary> methylSummary;
	while ( !BedFBedF.eof() ) {
		char buffer[MAX_LINE_LENGTH+1];
		BedFBedF.getline(buffer, MAX_LINE_LENGTH);
		if (!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		string tmp = buffer;
		vector<string> tokens = string_tokenize(tmp);
		string chr = tokens[0];
		if (methylSummary.find(chr) == methylSummary.end()) {
			MethylSummary ms  = {0, 0, 0, 0};
			methylSummary[chr] = ms;
		}
		int Leftpos = atoi( tokens[1].c_str() );
		int Rightpos = atoi( tokens[2].c_str() );
		MethylRegion mr = { 0, 0, 0, 0 };
		RegionMethylation( Methylome[chr], Leftpos, Rightpos, mr );
		char strand = tokens[3][0];
		if ( strand == '+' ) {
			methylSummary[chr].TotalMethyl_sense += mr.TotalMethyl_W;
			methylSummary[chr].TotalMethyl_anti  += mr.TotalMethyl_C;
			methylSummary[chr].TotalSite_sense   += mr.TotalSite_W;
			methylSummary[chr].TotalSite_anti    += mr.TotalSite_C;
		} else {
			methylSummary[chr].TotalMethyl_sense += mr.TotalMethyl_C;
			methylSummary[chr].TotalMethyl_anti  += mr.TotalMethyl_W;
			methylSummary[chr].TotalSite_sense   += mr.TotalSite_C;
			methylSummary[chr].TotalSite_anti    += mr.TotalSite_W;
		}
	}
	BedFBedF.close();

	//cout << "Sense     strand\t" << TotalMethyl_sense <<  "\t" << TotalSite_sense << endl;
	//cout << "Antisense strand\t" << TotalMethyl_anti  <<  "\t" << TotalSite_anti << endl;

	//if ( param.by_strand ) {
	map<string, MethylSummary>::iterator miter;
	cout << "chr\tsense_Count\tsense_mC\tanti_Count\tanti_mC\tall_Count\tall_mC\n";
	for (miter = methylSummary.begin(); miter!=methylSummary.end(); miter++) {
		cout.precision(4);
		string chr = miter->first;
		cout << chr;
		if (methylSummary[chr].TotalSite_sense > 0)  {
			cout << "\t" << methylSummary[chr].TotalSite_sense << "\t"
			<< (methylSummary[chr].TotalMethyl_sense / methylSummary[chr].TotalSite_sense);
		} else {
			cout << "\t0\tNaN" ;
		}
		if (methylSummary[chr].TotalSite_anti > 0)  {
			cout << "\t" << methylSummary[chr].TotalSite_anti << "\t"
			<< (methylSummary[chr].TotalMethyl_anti / methylSummary[chr].TotalSite_anti);
		} else {
			cout << "\t0\tNaN";
		}
		double TotalMethyl = methylSummary[chr].TotalMethyl_sense + methylSummary[chr].TotalMethyl_anti;
		unsigned int TotalSite = methylSummary[chr].TotalSite_sense + methylSummary[chr].TotalSite_anti;
		if (TotalSite > 0) {
			cout << "\t" << TotalSite << "\t" << TotalMethyl / TotalSite;
		} else {
			cout << "\t0\tNaN";
		}
		cout << endl;
	}
	/*} else {
		cout.precision(4);
		double TotalMethyl = TotalMethyl_sense + TotalMethyl_anti;			
		unsigned int TotalSite = TotalSite_sense + TotalSite_anti;
		if (TotalSite > 0) {
			cout << "total\t" << TotalSite << "\t" << TotalMethyl / TotalSite << endl;
		} else {
			cout << "total\t0\tNaN" << endl;
		}
	}*/
	return 1;
}


int main(int argc, char* argv[])
{
	init ();
	parse_command_line ( argc, argv );
	
	if (param.BEDfile == "") {
		exit_with_help();
	}

	map< string, vector<Site> > Methylome;

	ReadMethylomeFromCGmap ( Methylome );
	//OutputMethylome ( Methylome, "output.tmp" );
	AverageMethylInBins ( Methylome );
	return 0;
}



