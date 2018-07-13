#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Rtypes.h"
#include "TVectorT.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"
#include "UniDbRun.h"

using namespace std;

const int run_period = 7;

int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tfileQuality /path/to/digi/file\n";
		return -1;
	}

	// Load input file and get run number from file name
	TString file = argv[1];
	TFile * inFile = new TFile(file);

	string f = string(argv[1]);
	string r = f.substr( f.find(".")-9 , 4  );
	int run_number = atoi(r.c_str());

	// Load field value
	UniDbRun* pCurrentRun = UniDbRun::GetRun(run_period, run_number);
	if (pCurrentRun == 0) {
		cerr << "Run does not exist in database; cannot access current magnetic field value for cuts\n"
			<< "\tBailing...\n";
		return -1;
	}
	double map_current = 55.87;
	double * field_voltage = pCurrentRun->GetFieldVoltage();
	if (*field_voltage < 10){
		cerr << "Magnetic field not on, I haven't calibrated this!!\n"
			<< "\tBailing...\n";
		return -2;
	}
	else if( fabs( (*field_voltage) - 87) < 1){
		cout << "SP-41 at 1400A\n";
	}
	else if( fabs( (*field_voltage) - 107.7) < 1){
		cout << "SP-41 at 1800A\n";
	}
	else if( fabs( (*field_voltage) - 123.7) < 1){
		cout << "SP-41 at 2200A\n";
	}
	else{
		cerr << "Magnetic field at unknown value, I haven't calibrated this!!\n"
			<< "\tBailing...\n";
		return -3;
	}

		
	return 0;
}

