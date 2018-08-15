// Using BmnTOF1Detector class to store everything rather than in one large macro
// and then want to check against tofAnalysis.cpp to make sure same output

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>

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
#include "TRandom3.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"
#include "BmnTof1Digit.h"
#include "BmnTOF1Detector.h"
#include "UniDbRun.h"

const double pedBC1 = 69.2885;
const double pedBC2 = -11.7212;
const double pedBC3 = -25.4808;
const double pedBC4 = 126.067;
const int run_period = 7;

using namespace std;

struct TOFHit{
	double t,a,y;
	int st;
	bool hit;
	TOFHit(){
		t = 0;
		a = 0;
		y = 0;
		st = -1;
		hit = false;
	}
};

void findIdx( TClonesArray* data, int &index , double refT);
double randomStrips();
double fixWalk( double amp, double shift, double par0, double par1, double par2);
double randomStripDiff();
double randomPosDiff();

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tcalcPedestals /path/to/all/digi/files\n";
		return -1;
	}

	// Init BmnTOF1Detector and load calibration files and geometry
	TString name;
	int NDet = 20;
	BmnTOF1Detector * Plane[NDet];
	for (int i = 0; i < NDet; i++) {
		name = Form("Plane%d", i);
		cout << endl << "==============Init " << name << "=============" << endl;
		Plane[i] = new BmnTOF1Detector(i, 1, NULL);
			// Load calibration files
		string path = "/input/TOF400_LRcorr_RUN7_SRC.dat";
		Plane[i]->SetCorrLR( path ); // Temp fix for this function, because overlaps with old class function that I'm phasing out
		Plane[i]->SetStripShift("/input/TOF400_StripShift_RUN7_SRC.dat");
		Plane[i]->SetWalkFunc("/input/TOF400_TimeWalk_RUN7_SRC.dat");
		//Plane[i]->TestPrint(19); // Prints out calibration info from TXT files for given strip
			// Load geometry file
		path = "/macro/run/geofile_full_src_noZshift.root";
        	Plane[i]->SetGeoFile( path );
	}


	
	return 0;
}

void findIdx( TClonesArray* data, int &index , double refT){
	double minT = 1e4;
	for( int m = 0 ; m < data->GetEntriesFast() ; m++){
		BmnTrigWaveDigit * signal = (BmnTrigWaveDigit*)data->At(m);
		double time = fabs(signal->GetTime() - refT);
		if( time < minT){
			minT = time;
			index = m;
		}
	}
}


double fixWalk( double amp, double shift, double par0, double par1, double par2){

	if( par0 == -1 || par1 == -1 || par2 == -1 || shift == -1)
		return 0;

	return par0 + par1*exp( -(amp - shift) / par2 );
}

double randomStripDiff(){
	TRandom3 * rand = new TRandom3(0);
	double dist = 0.;
	while( dist == 0.){
		int st1 = rand->Rndm() * 47;
		int st2 = rand->Rndm() * st1; // lower range possible
		int st3 = rand->Rndm() * (47-st1) + st1; // upper range possible
		double choice = rand->Rndm(); // choice between the two
		int st4;
		if( choice < 0.5) st4 = st2;
		else	st4 = st3;

		if( st1 > st4) dist = (st1 - st4);
		else	dist = (st4 - st1);
	}
	delete rand;

	return dist;
}

double randomPosDiff(){
	TRandom3 * rand = new TRandom3(0);
	double dist = (rand->Rndm() * 30. - 15.) - (rand->Rndm() * 30. - 15.);

	delete rand;
	return dist;
}

