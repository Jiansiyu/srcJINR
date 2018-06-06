#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;
void findIdx( TClonesArray* data, int &index , double refT);
void fillPedestal( TH1D* hist, TClonesArray* data, int index);


int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tidentify12C /path/to/digi/file\n";
		return -1;
	}

	TString file = argv[1];
	TString run_number( file(file.Index(".")-9,4) );
	
	// Attempt to find a matchin pedestal file:
	TFile * infile_ped = NULL;
	TString path = "/home/segarrae/software/srcJINR/build/bin/tqdcPedestals" + run_number + ".root";
	infile_ped = new TFile(path);
	if (!infile_ped){
		cerr << "Could not open file " << path <<"\n"
			<< "\tBailing out\n";
		return -2;
	}
	else{
		//cerr << "Successfully opened file " << path << " and saved it to address " << infile_ped << "\n";
	}


	// Read in pedestals
	TH1D * BC1_ped = (TH1D * ) infile_ped->Get("BC1_ped");
	TH1D * BC2_ped = (TH1D * ) infile_ped->Get("BC2_ped");
	TH1D * BC3_ped = (TH1D * ) infile_ped->Get("BC3_ped");
	TH1D * BC4_ped = (TH1D * ) infile_ped->Get("BC4_ped");
	double pedBC1, pedBC2, pedBC3, pedBC4;
	pedBC1 = BC1_ped->GetXaxis()->GetBinCenter(BC1_ped->GetMaximumBin());
	pedBC2 = BC2_ped->GetXaxis()->GetBinCenter(BC2_ped->GetMaximumBin());
	pedBC3 = BC3_ped->GetXaxis()->GetBinCenter(BC3_ped->GetMaximumBin());
	pedBC4 = BC4_ped->GetXaxis()->GetBinCenter(BC4_ped->GetMaximumBin());
	infile_ped->Close();
	
	// Now we can read in the digi file and look at BC1, BC2, BC3, BC4 correlation with pedestal subtracted!
	TFile * inFile = NULL;
	inFile = new TFile(argv[1]);
	if (!inFile){
		cerr << "Could not open file " << argv[1] <<"\n"
			<< "\tBailing out\n";
		return -2;
	}
	else{
		//cerr << "Successfully opened file " << argv[1] << " and saved it to address " << inFile << "\n";
	}
	TTree * intree = (TTree *) inFile->Get("cbmsim");

	// Setup the tree to read in the counters
	TClonesArray * bc1Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc2Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc3Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc4Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * t0Data	= new TClonesArray("BmnTrigDigit");
	
	const int nEvents = intree->GetEntries();
	// TQDC Branches
	intree->SetBranchAddress("TQDC_BC1"	,&bc1Data);
	intree->SetBranchAddress("TQDC_T0"	,&bc2Data);
	intree->SetBranchAddress("TQDC_BC3"	,&bc3Data);
	intree->SetBranchAddress("TQDC_BC4"	,&bc4Data);
	
	// TDC Branches
	intree->SetBranchAddress("T0"		,&t0Data);


	// Setup an output file and tree
	TFile * outFile = new TFile("carbonIn.root","RECREATE");
	TTree * outTree = new TTree("events","");
	double adcBC1, adcBC2, adcBC3, adcBC4;
	double tBC1, tBC2, tBC3, tBC4;
	outTree->Branch("adcBC1"	,&adcBC1	,"adcBC1/D");
	outTree->Branch("adcBC2"	,&adcBC2	,"adcBC2/D");
	outTree->Branch("adcBC3"	,&adcBC3	,"adcBC3/D");
	outTree->Branch("adcBC4"	,&adcBC4	,"adcBC4/D");
	outTree->Branch("tBC1"		,&tBC1		,"tBC1/D");
	outTree->Branch("tBC2"		,&tBC2		,"tBC2/D");
	outTree->Branch("tBC3"		,&tBC3		,"tBC3/D");
	outTree->Branch("tBC4"		,&tBC4		,"tBC4/D");


	// Loop over events
	int cntT0 = 0;
	int cnt1 = 0;
	int cnt2 = 0;
	int cnt3 = 0;
	int cnt4 = 0;

	for (int event=0 ; event<nEvents ; event++)
	{
		//if (event % 1000== 0)
			//cerr << "Working on event " << event << "\n";

		intree->GetEvent(event);
		
		// Kill any event that doesn't have T0 entry, or that
		// has more than 1 TDC in T0 -- cuts out 7% of data	
		if( t0Data->GetEntriesFast() != 1) continue;
		cntT0 ++;

		// Do a DAQ counter
		if( bc1Data->GetEntriesFast() > 0) cnt1 ++;
		if( bc2Data->GetEntriesFast() > 0) cnt2 ++;
		if( bc3Data->GetEntriesFast() > 0) cnt3 ++;
		if( bc4Data->GetEntriesFast() > 0) cnt4 ++;

		// Get the T0 Time
		BmnTrigDigit * t0Signal = (BmnTrigDigit*) t0Data->At(0);
		double t0Time = t0Signal->GetTime();

		// For all the possible ADC inside of BC1, T0, I want to take
		// the one that has the closest time to T0 from TDC, and use the
		// others as the pedestal subtraction
		if( bc1Data->GetEntriesFast() && bc2Data->GetEntriesFast() && bc3Data->GetEntriesFast() && bc4Data->GetEntriesFast() ){
			int bc1Idx, bc2Idx, bc3Idx, bc4Idx;			

				// Find the index for BC1
			findIdx(bc1Data,bc1Idx,t0Time);
			
				// Find the index for BC2
			findIdx(bc2Data,bc2Idx,t0Time);
	
				// Find the index for BC3
			findIdx(bc3Data,bc3Idx,t0Time);
	
				// Find the index for BC4
			findIdx(bc4Data,bc4Idx,t0Time);

			// Now we use bc1Idx for actual analysis and use our pedestals
			BmnTrigWaveDigit * signalBC1 = (BmnTrigWaveDigit*) bc1Data->At(bc1Idx);	
			BmnTrigWaveDigit * signalBC2 = (BmnTrigWaveDigit*) bc2Data->At(bc2Idx);	
			BmnTrigWaveDigit * signalBC3 = (BmnTrigWaveDigit*) bc3Data->At(bc3Idx);	
			BmnTrigWaveDigit * signalBC4 = (BmnTrigWaveDigit*) bc4Data->At(bc4Idx);	

			adcBC1 = signalBC1->GetPeak() - pedBC1;
			adcBC2 = signalBC2->GetPeak() - pedBC2;
			adcBC3 = signalBC3->GetPeak() - pedBC3;
			adcBC4 = signalBC4->GetPeak() - pedBC4;
			tBC1   = signalBC1->GetTime() - t0Time;			
			tBC2   = signalBC2->GetTime() - t0Time;			
			tBC3   = signalBC3->GetTime() - t0Time;			
			tBC4   = signalBC4->GetTime() - t0Time;			

			outTree->Fill();
		}
	}

	// Get all the peaks for quality-control-check to see for significant drifts over time
	TH1D * h = new TH1D("trash","trash",400,0,4000);
	outTree->Draw("adcBC1>>trash");
	double bc1Peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	outTree->Draw("adcBC2>>trash");
	double bc2Peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	outTree->Draw("adcBC3>>trash");
	double bc3Peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
	outTree->Draw("adcBC4>>trash");
	double bc4Peak = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());

	// Look at the BC1-BC2 carbon-in peak and the fraction of events in that peak vs total events
	outTree->Draw("sqrt( adcBC1 * adcBC2 ) >>trash");
	TFitResultPtr fit = h->Fit("gaus","QES",0,4000);
	double mean = fit->Parameter(1);
	double sig = fit->Parameter(2);
	double maxX = mean + 1.*sig;
	TFitResultPtr fit2 = h->Fit("gaus","QES","",0., maxX);
	
	// Let's look at beam quality:
	int bmin = h->GetXaxis()->FindBin( fit2->Parameter(1) - 2*fit2->Parameter(2) );
	int bmax = h->GetXaxis()->FindBin( fit2->Parameter(1) + 2*fit2->Parameter(2) );
	double integral = h->Integral(bmin,bmax);
	double beamQ = integral / nEvents;

	// Let's do BC3-BC4 carbon out -- only useful for no-target runs
	outTree->Draw("sqrt( adcBC3 * adcBC4 ) >>trash");
	TFitResultPtr fitO = h->Fit("gaus","QES",0,4000);
	mean = fitO->Parameter(1);
	sig = fitO->Parameter(2);
	maxX = mean + 1.*sig;
	TFitResultPtr fit2O = h->Fit("gaus","QES","",0., maxX);
	
	// Let's look at beam quality:
	bmin = h->GetXaxis()->FindBin( fit2O->Parameter(1) - 2*fit2O->Parameter(2) );
	bmax = h->GetXaxis()->FindBin( fit2O->Parameter(1) + 2*fit2O->Parameter(2) );
	integral = h->Integral(bmin,bmax);
	double beamQout = integral / nEvents;
	
	// DAQ counters
	double rateT0 = cntT0 / (double) nEvents;
	double rate1  = cnt1 / (double) nEvents;
	double rate2  = cnt2 / (double) nEvents;
	double rate3  = cnt3 / (double) nEvents;
	double rate4  = cnt4 / (double) nEvents;
	
	cout << "\n**********************************************\n";
	cout << "\tPeak\tPedestal\tRate\n";
	cout << "BC1:\t" << bc1Peak << "\t" << pedBC1 << "\t\t" << rate1 << "\n";
	cout << "BC2:\t" << bc2Peak << "\t" << pedBC2 << "\t\t" << rate2 << "\n";
	cout << "BC3:\t" << bc3Peak << "\t" << pedBC3 << "\t\t" << rate3 << "\n";
	cout << "BC4:\t" << bc4Peak << "\t" << pedBC4 << "\t\t" << rate4 << "\n\n";
	cout << "Carbon-In Peak (BC1-BC2):\t" << fit2->Parameter(1) << "\n\twidth of:\t\t" << fit2->Parameter(2) << "\n\tbeam quality:\t\t" << beamQ<< "\n\n";
	cout << "Carbon-Out Peak (BC3-BC4):\t" << fit2O->Parameter(1) << "\n\twidth of:\t\t" << fit2O->Parameter(2) << "\n\tbeam quality:\t\t" << beamQout<< "\n\n";
	cout << "T0 Rate:\t" << rateT0 << "\n";
	cout << "**********************************************\n\n";
	

	outTree->Write();
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
