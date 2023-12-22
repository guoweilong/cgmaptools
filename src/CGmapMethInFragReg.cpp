/*
    cgmaptools - CGmapMethInFragReg.cpp

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

/* CGmapMethylInFragmentedRegion.cpp 2015-01-25
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
	string  Regionfile;
	int	min_coverage;
	int	max_coverage;
	string ctx;
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
		"Usage:	cgmaptools mfg [-i <CGmap>]  -r <region> [-c 5 -C 500] \n"
		"Description: Calculated methylation profile across fragmented regions.\n"
		"Contact:     Guo, Weilong; guoweilong@126.com\n"
		"Last Update: 2017-01-20\n"
		"Options:\n"
		"   -i  String, CGmap file; use STDIN if not specified\n"
		"       Please use \"gunzip -c <input>.gz \" and pipe as input for gzipped file.\n"
		"       chr1\tG\t851\tCHH\tCC\t0.1\t1\t10\n"
		"   -r  String, Region file, at least 4 columns\n"
		"       Format: chr\tstrand\tpos_1\tpos_2\tpos_3\t...\n"
		"       Regions would be considered as [pos_1, pos_2), [pos_2, pos_3)\n"
		"       Strand information will be used for distinguish sense/antisense strand\n"
		"       Ex:\n"
		"       #chr\tstrand\tU1\tR1\tR2\tD1\tEnd\n"
		"       chr1\t+\t600\t700\t800\t900\t950\n"
		"       chr1\t-\t1600\t1500\t1400\t1300\t1250\n"
		"   -c  Int, minimum Coverage [Default: 5] \n"
		"   -C  Int, maximum Coverage [Default: 500] \n"
		"       Sites exceed the coverage range will be discarded\n"
		"   -x  String, context [use all sites by default] \n"
		"       string can be CG, CH, CHG, CHH, CA, CC, CT, CW\n"
		"   -h  help\n"
		"Output to STDOUT:\n"
		"   Region_ID       U1      R1      R2      D1\n"
        "   sense_ave_mC    0.50    0.40    0.30    0.20\n"
        "   sense_sum_mC    5.0     4.0     3.0     2.0\n"
        "   sense_sum_NO    10      10      10      10\n"
        "   anti_ave_mC     0.40    0.20    0.10    NaN\n"
        "   anti_sum_mC     8.0     4.0     2.0     0.0\n"
        "   anti_sum_NO     20      20      20      0\n"
        "   total_ave_mC    0.43    0.27    0.17    0.2\n"
        "   total_sum_mC    13.0    8.0     5.0     2.0\n"
        "   total_sum_NO    30      30      30      10\n"
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
			param.CGmapfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
			break;
		case 'r':
			if(i == argc)	exit_with_help();	
			param.Regionfile = string(argv[i]);
			if(++i>argc)	exit_with_help();
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
		case 'x':
			if(i == argc)	exit_with_help();
			param.ctx = string(argv[i]);
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
	if (  !param.Regionfile.length() ){
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
	param.Regionfile = "";
	param.min_coverage = 5;
	param.max_coverage = 500;
	param.ctx = "";
}


#define MAX_LINE_LENGTH 10000

int ValideCtx( string ctx, string dinuc) {
    if (param.ctx == "") return 1;
    if (param.ctx == "CG") {
        if (ctx == "CG") {return 1;} else {return 0;}
    }
    if (param.ctx == "CA") {
        if (dinuc == "CA") {return 1;} else {return 0;}
    }
    if (param.ctx == "CC") {
        if (dinuc == "CC") {return 1;} else {return 0;}
    }
    if (param.ctx == "CT") {
        if (dinuc == "CT") {return 1;} else {return 0;}
    }

    if (param.ctx == "CH") {
        if (ctx=="CHG" || ctx=="CHH") {return 1;} else {return 0;}
    }

    if (param.ctx == "CHG") {
        if (ctx=="CHG") {return 1;} else {return 0;}
    }
    if (param.ctx == "CHH") {
        if (ctx=="CHH") {return 1;} else {return 0;}
    }
    if (param.ctx == "CW") {
        if (dinuc=="CA" || dinuc=="CT") {return 1;} else {return 0;}
    }
    return 0;
}

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
		    if (ValideCtx(tokens[3], tokens[4])) {
			    Genome[chr].push_back(site);
			}
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
/*
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
			mr.TotalMethyl_W += Chr[i].ML;
			mr.TotalSite_W ++;
		} else {
			mr.TotalMethyl_C += Chr[i].ML;
			mr.TotalSite_C ++;
		}

	return 1;
}

*/

struct MethylSummary {
	double 	TotalMethyl_sense;
	double 	TotalMethyl_anti;
	unsigned int 	TotalSite_sense;
	unsigned int 	TotalSite_anti;
} ;

int RegionsMethylation( vector<Site> & Chr, vector<int> & pos_list, vector<MethylRegion> & mr) {
	int n_pos = pos_list.size();
	int nregion = n_pos - 1;
	int low, high, mid;
	int LeftMostPos = pos_list[0];
	//for(int i=0; i<n_pos; i++) {
	//    cout << pos_list[i] << "\t";
	//}
	//cout << endl;
	//cerr << "LeftMost=" << LeftMostPos << endl;
	low = 0;
	if (Chr.size() == 0) return 0;
	high = Chr.size() - 1;
	mid = ( low + high ) / 2;
	while ( low <= high ) {
		mid = ( low + high ) / 2;
		//cerr << "low=" << low << "\t mid=" << mid  << "\t high=" << high << endl;
		//cerr << "low=" << Chr[low].pos << "\t mid=" << Chr[mid].pos  << "\t high=" << Chr[high].pos << "\tLeftMostPos=" << LeftMostPos  << endl;
		if ( LeftMostPos  > Chr[mid].pos ) {
		// Not >=, so that can get low as the lowest one which >= Left
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	//cerr << "low=" << low << "\t mid=" << mid  << "\t high=" << high << endl;
	//cerr << "low=" << Chr[low].pos << "\t mid=" << Chr[mid].pos  << "\t high=" << Chr[high].pos << "\tLeftMostPos=" << LeftMostPos  << endl;

	int chr_size = Chr.size();
	int RightMostPos = pos_list[n_pos-1];
	int RegionID = 0;
	//cerr << "Chr[low].pos=" << Chr[low].pos << endl;
	for ( unsigned int i = low; i<chr_size && Chr[i].pos < RightMostPos; i++ ) {
		while ( Chr[i].pos >= pos_list[RegionID+1] && RegionID<(n_pos-1) ) {
			RegionID++;
		}
		if (Chr[i].pos < pos_list[RegionID+1]) {
		    //cerr << "pos=" << Chr[i].pos << "\t" << Chr[i].ML << endl;
			if ( Chr[i].strand == '+' ) {
				mr[RegionID].TotalMethyl_W += Chr[i].ML;
				//cerr << "ML\t" << Chr[i].ML << endl;
				mr[RegionID].TotalSite_W ++;
			} else {
				mr[RegionID].TotalMethyl_C += Chr[i].ML;
				mr[RegionID].TotalSite_C ++;
			}
		}
	}


	return 1;
}

int AverageMethylInFragmentedRegion ( map< string, vector<Site> > & Methylome ) {
	// Read the first line of RegionFile
	ifstream RegionF(param.Regionfile.c_str());
	if(!RegionF) {
		cout << "cannot open input file" << param.Regionfile.c_str() << endl;
		return -1;
	}
	char buffer[MAX_LINE_LENGTH+1];
	RegionF.getline(buffer, MAX_LINE_LENGTH);
	if ( !strlen(buffer) ) {
		return 0;
	}
	if ( buffer[strlen(buffer)-1] == '\r' ) {
		buffer[strlen(buffer)-1] = '\0';
	}
	string tmp = buffer;
	vector<string> tokens = string_tokenize(tmp);
	vector<string> header_tokens = string_tokenize(tmp);
	int ncol = tokens.size();
	int nregion = ncol - 3;
	vector<MethylSummary> ms;
	if ( ncol > 3 ) {
	    //Create a vector and initiation
		//ms = new MethylSummary[nregion];
		for (int i = 0; i < nregion; i++) {
			MethylSummary ms_tmp = {0, 0, 0, 0};
		    ms.push_back(ms_tmp);
		}
//		cerr << "Detect " << ncol << " columns in the 1st row.\n"
//		     << "Region Number: " << nregion << "\n";
	} else {
		cerr << "Error : The 1st line should have no less than 4 columns\n";
		exit(-1);
	}
	RegionF.close();
	//
	//Start to count by re-read the region file
	RegionF.open(param.Regionfile.c_str(), std::ifstream::in);
	//map<string, MethylSummary> methylSummary;
	int lineno = 0;
	//cerr << "Before while\n";
    RegionF.getline(buffer, MAX_LINE_LENGTH);
    lineno++;
    if (!strlen(buffer)) return 0;
    if (buffer[strlen(buffer)-1] == '\r') {
        buffer[strlen(buffer)-1] = '\0';
    }
    tmp = buffer;
    header_tokens = string_tokenize(tmp);
	cout.precision(4);
	//    The Title line
	cout << "Region_ID";
	for (int i = 2; i < header_tokens.size()-1; i++ ) {
		cout << "\t" << header_tokens[i] ;
	}
	cout << endl;
	while ( !RegionF.eof() ) {
		RegionF.getline(buffer, MAX_LINE_LENGTH);
		lineno++;
		if (!strlen(buffer)) continue;
		if (buffer[strlen(buffer)-1] == '\r') {
			buffer[strlen(buffer)-1] = '\0';
		}
		tmp = buffer;
		//cerr << tmp << endl;
		tokens = string_tokenize(tmp);
		int ncol_cur = tokens.size();
		if (ncol_cur != ncol) {
			cerr << "Line\t" << lineno << " : has " << ncol << " columns, different from the 1st line.\n";
			continue;
		}
		//cerr << "1\n";
		string chr = tokens[0];
		char strand = tokens[1][0];
		vector<int> pos_list;
		//cerr << "2\n";
		//cerr << tokens.size() << endl;
		if (strand == '+') {
			for (int i = 2; i <= nregion + 2; i++) {
				//cerr << "i=" << i << "\n";
				//cerr << "nregion=" << nregion << "\n";
				pos_list.push_back( atoi(tokens[i].c_str()) );
			}
		} else {
			for (int i = nregion + 2; i >= 2; i--) {
				pos_list.push_back( atoi(tokens[i].c_str())+1 );
			}
		}
		//cerr << "3\n";
		vector<MethylRegion> mr;
		for (int i = 0; i < nregion; i++) {
			MethylRegion mr_tmp = {0, 0, 0, 0};
			mr.push_back(mr_tmp);
		}
		//cerr << "lineno\t:" << lineno << "\n";
		RegionsMethylation(Methylome[chr], pos_list, mr);
		if (strand == '+') {
			for (int i = 0; i < nregion; i++) {
				ms[i].TotalMethyl_sense += mr[i].TotalMethyl_W;
				ms[i].TotalMethyl_anti  += mr[i].TotalMethyl_C;
				ms[i].TotalSite_sense   += mr[i].TotalSite_W;
				ms[i].TotalSite_anti    += mr[i].TotalSite_C;
			}
		} else {
			for (int i = 0; i < nregion; i++) {
				ms[i].TotalMethyl_sense += mr[nregion-1-i].TotalMethyl_C;
				ms[i].TotalMethyl_anti  += mr[nregion-1-i].TotalMethyl_W;
				ms[i].TotalSite_sense   += mr[nregion-1-i].TotalSite_C;
				ms[i].TotalSite_anti    += mr[nregion-1-i].TotalSite_W;
			}
		}
	}
	RegionF.close();

	// Start to output summary

	//    The "sense_mC" line
	cout << "sense_ave_mC";
	for (int i = 0; i < nregion; i++ ) {
		if (ms[i].TotalSite_sense > 0)  {
			cout << "\t" << ms[i].TotalMethyl_sense/ms[i].TotalSite_sense;
		} else {
			cout << "\tNaN" ;
		}
	}
	cout << endl;
	cout << "sense_sum_mC"; for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalMethyl_sense; }	cout << endl;
	cout << "sense_sum_NO";for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalSite_sense; }	cout << endl;
	// The "Antisense_mC" line
	cout << "anti_ave_mC";
	for (int i = 0; i < nregion; i++ ) {
		if (ms[i].TotalSite_anti > 0)  {
			cout << "\t" << ms[i].TotalMethyl_anti/ms[i].TotalSite_anti;
		} else {
			cout << "\tNaN";
		}
	}
	cout << endl;
	cout << "anti_sum_mC"; for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalMethyl_anti; }	cout << endl;
	cout << "anti_sum_NO";for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalSite_anti; }	cout << endl;
	// The "total_mC" line
	cout << "total_ave_mC";
	for (int i = 0; i < nregion; i++ ) {
		double TotalMethyl = ms[i].TotalMethyl_sense + ms[i].TotalMethyl_anti;
		unsigned int TotalSite = ms[i].TotalSite_sense + ms[i].TotalSite_anti;
		if (TotalSite > 0) {
			cout << "\t" << TotalMethyl / TotalSite;
		} else {
			cout << "\tNaN";
		}
	}
	cout << endl;
	cout << "total_sum_mC"; for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalMethyl_sense + ms[i].TotalMethyl_anti; }	cout << endl;
	cout << "total_sum_NO";for (int i = 0; i < nregion; i++ ) { cout << "\t" << ms[i].TotalSite_sense + ms[i].TotalSite_anti; }	cout << endl;
	//cerr << "Finisheddd" << endl;
	return 1;
}


int main(int argc, char* argv[])
{
	init ();
	parse_command_line ( argc, argv );
	
	if (param.Regionfile == "") {
		exit_with_help();
	}

	map< string, vector<Site> > Methylome;

	ReadMethylomeFromCGmap ( Methylome );
	//OutputMethylome ( Methylome, "output.tmp" );
	AverageMethylInFragmentedRegion ( Methylome );
	//cerr << "Finished" << endl;
	return 1;
}



