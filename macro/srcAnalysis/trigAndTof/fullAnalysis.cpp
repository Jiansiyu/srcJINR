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
#include "TRandom3.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"
#include "BmnTof1Digit.h"
#include "UniDbRun.h"

const double pedBC1 = 69.2885;
const double pedBC2 = -11.7212;
const double pedBC3 = -25.4808;
const double pedBC4 = 126.067;
const int run_period = 7;

using namespace std;

struct TOFHit{
	double t,a,y;
	bool hit;
	TOFHit(){
		t = 0;
		a = 0;
		y = 0;
		hit = false;
	}
};

void findIdx( TClonesArray* data, int &index , double refT);
double randomStrips();

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tcalcPedestals /path/to/all/digi/files\n";
		return -1;
	}

	// Try opening the LR corr file if it exists:
	ifstream f_corr;
	string dir = std::getenv("VMCWORKDIR");
	dir += "/input/TOF400_LRcorr_RUN7_SRC.dat";
	cout << "Attempting to open LR corr file: " << dir << "\n";
	f_corr.open(dir);
	char line[256];
	f_corr.getline(line, 256);
	f_corr.getline(line, 256);
	int Pl, St;
	double Temp;
	double corrLR[20][48] = {0.};
	if (f_corr.is_open() == true){
 		while (!f_corr.eof()) {
			f_corr >> Pl >> St >> Temp;
			corrLR[Pl][St] = Temp;
		}
		cout << "\tLoaded LR corr file\n";
	}
	else{
		cout << "\tFailed to find LR corr file, setting all corrections to 0...\n";
	} 
	
	// Try opening the strip shifts file if it exists:
	ifstream f_corr_shift;
	dir = std::getenv("VMCWORKDIR");
	dir += "/input/TOF400_StripShift_RUN7_SRC.dat";
	cout << "Attempting to open StripShift file: " << dir << "\n";
	f_corr_shift.open(dir);
	char line2[256];
	f_corr_shift.getline(line2, 256);
	f_corr_shift.getline(line2, 256);
	int Pl2, St2;
	double Temp2, Temp3;
	double stripShift[20][48] = {0.};
	bool killStrip[20][48] = { false };
	if (f_corr_shift.is_open() == true){
 		while (!f_corr_shift.eof()) {
			f_corr_shift >> Pl2 >> St2 >> Temp2 >> Temp3;
			stripShift[Pl2][St2] = Temp2;
			if( stripShift[Pl2][St2] == -1)
				killStrip[Pl2][St2] = true;
		}
		cout << "\tLoaded strip shift file\n";
	}
	else{
		cout << "\tFailed to find strip shift file, setting all corrections to 0...\n";
	} 


	// Histograms for looking at cuts
	TH1D * hMult = new TH1D("hMult","hMult",20,0,20);
	TH1D ** hHits = new TH1D*[20];
	TString hName;
	for( int pl = 0 ; pl < 20 ; pl++){
		hName = Form("hHits_%i",pl);
		hHits[pl] = new TH1D(hName,hName,20,0,20);
	}

	const int files = argc - 1;
	for( int fi = 0 ; fi < files ; ++fi){
		// Set up the input file
		TFile * infile = NULL;
		infile = new TFile(argv[fi+1]);
		if (infile->IsZombie())
		{
			cerr << "Could not open file " << argv[fi+1] <<"\n"
				<< "\tBailing out\n";
			return -2;
		}

		// Get the run number from input file:
		TString file = argv[fi+1];
		TString run_number( file(file.Index(".")-9,4) );
		
		// Grab magnetic field current for storing different histograms
		int runNo = atoi( run_number.Data() );
		UniDbRun* pCurrentRun = UniDbRun::GetRun(run_period,runNo);
		double field_voltage;
		if( pCurrentRun == 0 ){
			if( (runNo > 3474 && runNo < 3485) || (runNo > 3434 && runNo < 3443) )
				field_voltage = 87.0;
			else if( (runNo > 3484 && runNo < 3496) || (runNo > 3442 && runNo < 3453) || (runNo > 3513) )
				field_voltage = 107.7;
			else if( (runNo > 3495 && runNo < 3501) || (runNo > 3452 && runNo < 3458) )
				field_voltage = 123.7;
		}
		else{
			field_voltage = *(pCurrentRun->GetFieldVoltage());
		}

		// With this run number, find corresponding check file so that we can
		// pull up the carbon in cut
		TFile * qualityFile = NULL;
		TString path = std::getenv("VMCWORKDIR");
		path = path + "/build/bin/qualityCheck/checked_" + run_number + ".root";
		qualityFile = new TFile(path);
		if( qualityFile->IsZombie() ){
			cerr << "No checked file for this run number. You need to run calcPedestals first.\n"
				<< "\tBailing...\n";
			return -3;
		}

		TVectorT<double> * carbonIn	= (TVectorT<double>*)qualityFile->Get("carbonIn");		
		TVectorT<double> * carbonInWidth= (TVectorT<double>*)qualityFile->Get("carbonInWidth");	

		double BC1_CPeak = (*carbonIn)[0];
		double BC2_CPeak = (*carbonIn)[1];
		double BC1_CWidth= (*carbonInWidth)[0];
		double BC2_CWidth= (*carbonInWidth)[1];

		
		// Set up the tree
		TClonesArray * bc1Data 	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * bc2Data	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * t0Data	= new TClonesArray("BmnTrigDigit");
		TClonesArray * tofData  = new TClonesArray("BmnTof1Digit");

		TClonesArray * bc3Data 	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * bc4Data  = new TClonesArray("BmnTrigWaveDigit");


		TTree * intree = NULL;
		intree = (TTree*) infile->Get("cbmsim");
		if (!intree)
		{
			cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";
			return -3;
		}
		else
		{
			//cerr << "Successfully opened file " << argv[fi+1] << " and saved it to address " << infile << "\n";
			//cerr << "Successfully loaded tree at address " << intree << " with events " << intree->GetEntries() << "\n";
		}
		
		const int nEvents = intree->GetEntries();
		cout << "Working on file " << argv[fi+1] << " with " << nEvents << " events and field voltage " << field_voltage << "\n";
		// TQDC Branches
		intree->SetBranchAddress("TQDC_BC1"	,&bc1Data);
		intree->SetBranchAddress("TQDC_T0"	,&bc2Data);
		intree->SetBranchAddress("TQDC_BC3"	,&bc3Data);
		intree->SetBranchAddress("TQDC_BC4"	,&bc4Data);
		// TDC Branches
		intree->SetBranchAddress("T0"		,&t0Data);
		// ToF400 Branches
		intree->SetBranchAddress("TOF400"	,&tofData);


		// Loop over events
		for (int event=0 ; event<nEvents ; event++){

			double adcBC1 = 0;
			double adcBC2 = 0;
			double adcBC3 = 0;
			double adcBC4 = 0;
			t0Data->Clear();
			bc1Data->Clear();	
			bc2Data->Clear();	
			bc3Data->Clear();
			bc4Data->Clear();
			tofData->Clear();
		

			intree->GetEvent(event);

			// Kill any event that doesn't have T0 entry, or that has more than one TDC
			if( t0Data->GetEntriesFast() != 1) continue;

			BmnTrigDigit * t0Signal = (BmnTrigDigit*) t0Data->At(0);
			double t0Time = t0Signal->GetTime();
			double t0Amp  = t0Signal->GetAmp();
				
			
			if( bc1Data->GetEntriesFast() ){
				int bc1Idx;
				findIdx(bc1Data,bc1Idx,t0Time);
				BmnTrigWaveDigit * signal = (BmnTrigWaveDigit *) bc1Data->At(bc1Idx);
				adcBC1 = signal->GetPeak() - pedBC1;
			}

			if( bc2Data->GetEntriesFast() ){
				int bc2Idx;
				findIdx(bc2Data,bc2Idx,t0Time);
				BmnTrigWaveDigit * signal = (BmnTrigWaveDigit *) bc2Data->At(bc2Idx);
				adcBC2 = signal->GetPeak() - pedBC2;
			}
			
			if( fabs( adcBC1 - (BC1_CPeak-pedBC1) ) > 2*BC1_CWidth )
				continue;
			if( fabs( adcBC2 - (BC2_CPeak-pedBC2) ) > 2*BC2_CWidth )
				continue;

			TOFHit tofEvents[20][48];
			std::vector<double> hitInfoTime[20][48];
			std::vector<double> hitInfoAmps[20][48];

			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){ // Loop through the events once and store left side info
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				int plane = signal->GetPlane();
				int strip = signal->GetStrip();
				int side = signal->GetSide();
				double tofT = signal->GetTime();
				double tofAmp = signal->GetAmplitude();
				if( plane < 0 || plane > 19) continue;
				if( strip < 0 || strip > 47) continue;
				if( killStrip[plane][strip] == true) continue;
				if( side == 1) continue;				

				if( side == 0){ // Store info
					hitInfoTime[plane][strip].push_back( tofT );
					hitInfoAmps[plane][strip].push_back( tofAmp );
				}
			}

			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){// Loop through events again and look for only the right side to create hit
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				int plane = signal->GetPlane();
				int strip = signal->GetStrip();
				int side = signal->GetSide();
				double tofT = signal->GetTime();
				double tofAmp = signal->GetAmplitude();
				if( plane < 0 || plane > 19) continue;
				if( strip < 0 || strip > 47) continue;
				if( killStrip[plane][strip] == true) continue;
				if( side == 0) continue;

				if( side == 1){
					
					for( int hit = 0 ; hit < hitInfoTime[plane][strip].size() ; hit++){
						double tmpT = hitInfoTime[plane][strip].at(hit);
						double tmpA = hitInfoAmps[plane][strip].at(hit);

						if( fabs( (1./0.06)*(tmpT - tofT + corrLR[plane][strip]) ) < 32 ){ // Create a valid strip hit
							double meanTime = (tmpT + tofT + corrLR[plane][strip])*0.5;
							double sumAmps = sqrt(tmpA * tofAmp);
							double pos = 0.5*(1./0.06)*(tmpT - tofT + corrLR[plane][strip]); // +/- 15cm
						
							// If there was already a created hit for this strip, only take the earliest one
							// so there will only be one hit per strip allowed
							if( tofEvents[plane][strip].hit == false || (tofEvents[plane][strip].hit == true && meanTime < tofEvents[plane][strip].t) ){
								tofEvents[plane][strip].t = meanTime - t0Time + (-6.1 + 27./sqrt(t0Amp) ) - stripShift[plane][strip];
								tofEvents[plane][strip].a = sumAmps;
								tofEvents[plane][strip].y = pos;
								tofEvents[plane][strip].hit = true;
							}
						}
					}
				}
			}
		
			// First I need to understand the multiplicity of ToF400 in our normal data running
			int totcnt = 0;
			for( int pl = 0 ; pl < 20 ; pl++){
				int cnt = 0;
				std::vector<int> fired;
				for( int st = 0 ; st < 48 ; st++)
					if( tofEvents[pl][st].hit == true){
						cnt++;
					}
				if( cnt > 0 ) hHits[pl] -> Fill( cnt );
				if( cnt == 2){
				}
				if( cnt == 3){
				}
				if( cnt == 4){
				}
				totcnt+= cnt;
			}
			hMult->Fill( totcnt );
						

		} // End of loop over events in file

	} // End of loop over files	
	
	
	TFile * outFile = new TFile("realData.root","RECREATE");
	outFile->cd();
	hMult->Write();
	for( int pl = 0; pl < 20 ; pl++){
		hHits[pl]	-> Write();
	}
	outFile->Write();
	outFile->Close();


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

double randomStrips(){
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

		if( st1 > st4) dist = (st1 - st4)*12.5;
		else	dist = (st4 - st1)*12.5;
	}
	delete rand;
	
	return dist;
}
