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
#include "BmnTof1Digit.h"
#include "UniDbRun.h"

const double pedBC1 = 96.7212;
const double pedBC2 = 99.625;
const double pedBC3 = -15.8942;
const double pedBC4 = 146.843;
const int run_period = 7;

using namespace std;

void findIdx( TClonesArray* data, int &index , double refT);

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
	}
	else{
		cout << "Failed to find LR corr file, setting all corrections to 0...\n";
	} 


	// Histograms for looking at cuts
	TH1D *** hLR_time	= new TH1D**[20];
	TH1D *** hL_amp		= new TH1D**[20];
	TH1D *** hR_amp		= new TH1D**[20];

	TH2D ** hHits		= new TH2D*[20];
	
	TString hName;
	for( int pl = 0; pl < 20 ; pl++){

		hLR_time[pl] = new TH1D*[48];
                hL_amp[pl] = new TH1D*[48]; 
                hR_amp[pl] = new TH1D*[48]; 
		
		hName	= Form("hHits_%i",pl);
		hHits[pl]	= new TH2D(hName,hName,48,0,48,15,0,15);
		
		for( int st = 0 ; st < 48 ; st++){
			hName	= Form("hLR_time_%i_%i",pl,st);
			hLR_time[pl][st] = new TH1D(hName,hName,4000,-50,50);
		
			hName	= Form("hL_amp_%i_%i",pl,st);
			hL_amp[pl][st]	= new TH1D(hName,hName,500,0,50);
			hName	= Form("hR_amp_%i_%i",pl,st);
			hR_amp[pl][st]	= new TH1D(hName,hName,500,0,50);
		}

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
		// TDC Branches
		intree->SetBranchAddress("T0"		,&t0Data);
		// ToF400 Branches
		intree->SetBranchAddress("TOF400"	,&tofData);

		// Loop over events
		for (int event=0 ; event<nEvents ; event++){
			double adcBC1 = 0;
			double adcBC2 = 0;
			t0Data->Clear();
			bc1Data->Clear();	
			bc2Data->Clear();	
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

			int flagHit[20][48] = { 0 };	// Storage to tell if a strip had a hit or not
			int numHits[20][48][2] = { 0 };	// Storage for the total number of hits in a strip -- basically should be 1 for all because that's what we see in multiplicity
			double amps[20][48][2] = { 0. };	// Storage for amplitudes
			double times[20][48][2] = { 0. };	// Storage for times
		
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				int plane = signal->GetPlane();
				int strip = signal->GetStrip();
				int side = signal->GetSide();
				double tofT = signal->GetTime();
				double tofAmp = signal->GetAmplitude();
				if( plane < 0 || plane > 19) continue;
				if( strip < 0 || strip > 47) continue;
				
				amps[plane][strip][side] = tofAmp;
				times[plane][strip][side] = tofT;
				numHits[plane][strip][side]++;

				if( times[plane][strip][0] != 0. && times[plane][strip][1] != 0. && \
				     amps[plane][strip][0] != 0. &&  amps[plane][strip][1] != 0. && \
				    fabs( (1./0.06)*(times[plane][strip][0] - times[plane][strip][1] + corrLR[plane][strip]) ) < 32 ){
					
					if( numHits[plane][strip][0] == 1 && numHits[plane][strip][1] == 1 && flagHit[plane][strip] == 0){ // we have a matched hit
						flagHit[plane][strip] = 1;
						
						hLR_time[plane][strip] -> Fill( times[plane][strip][0] - times[plane][strip][1] + corrLR[plane][strip] ); 
						//hMeanTime[plane][strip]-> Fill( 0.5*(times[plane][strip][0] + times[plane][strip][1] + corrLR[plane][strip] ) );
						
						hL_amp[plane][strip]  -> Fill(  amps[plane][strip][0] );
						hR_amp[plane][strip]  -> Fill(  amps[plane][strip][1] );
						
						hHits[plane]->Fill( strip , 1 );
						
					}

					// The above will let me align strips, but this goes to a question of how often to I have events that are not
					// a clean single hit in L/R -- but either multi hits in one side, or actual multi signals in a strip (i.e. double hits)

					// Now for a single event if a strip has more than one matched hit:
					else if( numHits[plane][strip][0] == 1 && numHits[plane][strip][1] > 1)	{}	// These are super rare to happend
					else if( numHits[plane][strip][0] == 1 && numHits[plane][strip][1] == 0){}	// so i'm going to just ignore these
					else if( numHits[plane][strip][0] > 1 && numHits[plane][strip][1] == 1)	{}	// events completely
					else if( numHits[plane][strip][0] == 0 && numHits[plane][strip][1] == 1){}
					else if( numHits[plane][strip][0] > 1 && numHits[plane][strip][1] > 1)	{	// Only 4% of time this occurs
						if( numHits[plane][strip][0] == 2 && numHits[plane][strip][1] == 2)
							hHits[plane]->Fill( strip , 2 );
						else if( numHits[plane][strip][0] == 3 && numHits[plane][strip][1] == 3)
							hHits[plane]->Fill( strip , 3 );
						else if( numHits[plane][strip][0] == 4 && numHits[plane][strip][1] == 4)
							hHits[plane]->Fill( strip , 4 );
						else
							hHits[plane]->Fill( strip , 10 );
					}
				}
				// double tdiff = tofT - t0Time +  (-6.1 + 27./sqrt(t0Amp) );
				
			}
			

			
		}
	}	
	
	//TApplication theApp("App",&argc,argv);

	//theApp.Run();
	
	TFile * outFile = new TFile("alignStrips.root","RECREATE");
	outFile->cd();
	for( int pl = 0; pl < 20 ; pl++){
 		for( int st = 0 ; st < 48 ; st++){
			//hLR_time[pl][st]	-> Write();
			//hL_amp[pl][st]	-> Write();
			//hR_amp[pl][st]	-> Write();
		}
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

