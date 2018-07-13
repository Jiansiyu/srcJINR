#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>

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

	// Histograms for looking at cuts
	TH1D * hMultTot 	= new TH1D("MultTot","",50,0,50);
	TH1D ** hMultPlane	= new TH1D*[20];
	TH1D *** hMultStrip	= new TH1D**[20];
	TH1D **** hMultSide	= new TH1D***[20];
	TH1D ** hStripsFired	= new TH1D*[20];
	TH1D ** hSingleStripTime= new TH1D*[20];
	TH1D ** hManyStripTime	= new TH1D*[20];
	TH1D ** hManyStripFastTime= new TH1D*[20];


	TString hName;
	for( int pl = 0; pl < 20 ; pl++){
		hName	= Form("MultPlane_%i",pl);
		hMultPlane[pl]	= new TH1D(hName,"",50,0,50);
		
		hName	= Form("StripsFired_%i",pl);
		hStripsFired[pl]= new TH1D(hName,"",2400,-600,600);
		
		hName	= Form("ManyStripTime_%i",pl);
		hManyStripTime[pl] = new TH1D(hName,"",8000,200,400);
		hName	= Form("ManyStripFastTime_%i",pl);
		hManyStripFastTime[pl] = new TH1D(hName,"",8000,200,400);
		hName	= Form("SingleStripTime_%i",pl);
		hSingleStripTime[pl] = new TH1D(hName,"",8000,200,400);
		
		hMultStrip[pl] = new TH1D*[48];
		hMultSide[pl]  = new TH1D**[48];
		
		for( int st = 0 ; st < 48 ; st++){
			hName	= Form("MultStrip_%i_%i",pl,st);
			hMultStrip[pl][st] = new TH1D(hName,"",50,0,50);
		
			hMultSide[pl][st] = new TH1D*[2];
			for( int si = 0 ; si < 2 ; si++){
				hName	= Form("MultSide_%i_%i_%i",pl,st,si);
				hMultSide[pl][st][si] = new TH1D(hName,"",50,0,50);
			}	

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

			int multSide[20][48][2] = { 0 };
			int multStrip[20][48]	= { 0 };
			int multTotal = 0;
			std::vector<double> stripsFired[20];
			std::vector<double> dt[20];
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				int plane = signal->GetPlane();
				int strip = signal->GetStrip();
				int side = signal->GetSide();
				double tofT = signal->GetTime();
				double tofAmp = signal->GetAmplitude();

				//if( tofAmp > 15 || tofAmp < 6) continue;
				
				// Figure out multiplicities	
				if( plane > -1 && plane < 20 && strip > -1 && strip < 48){
					multStrip[plane][strip]++;		// for each strip, how many hits are there? -- 2
					multSide[plane][strip][side]++;		// for each side,  how many hits are there? -- 1
				
					double tdiff = tofT - t0Time +  (-6.1 + 27./sqrt(t0Amp) );
					// PROBLEM -- in order to use time, i really need to align strip sides
					// and take the average time, not just one side time. so first I need
					// to align the strips L-R and create single strip hits

					if( multSide[plane][strip][0] == 1 && multSide[plane][strip][1] == 1 ){ // If I have a fired strip in the plane, save that strip
						stripsFired[plane].push_back(strip);
						dt[plane].push_back(tdiff);
					}
				}
				
			}	
			
			// Now I want to know for a given event that has multiple strips firing,
			// are they close to each other (i.e. they are the same event) or are they
			// very far from each other?
			

			for( int pl = 0 ; pl < 20 ; pl++ ){
				if( stripsFired[pl].size() > 0){
					multTotal++;					// How many planes fire in a given trigger event
					hMultPlane[pl]->Fill( stripsFired[pl].size() );	// How many strips fired in a given plane
					if( stripsFired[pl].size() > 1){		// For more than one strip firing, how does the distance between hits look like?
						double mean = std::accumulate(stripsFired[pl].begin(),stripsFired[pl].end(),0.)/stripsFired[pl].size();
						double fastest = 1e4;
						for( int val = 0 ; val < stripsFired[pl].size() ; val++){
							hStripsFired[pl] -> Fill( (stripsFired[pl].at(val) - mean)*12.5 ); // in units of mm
							hManyStripTime[pl] -> Fill( dt[pl].at(val) );
							if( dt[pl].at(val) < fastest) fastest = dt[pl].at(val);
						}	
						hManyStripFastTime[pl]-> Fill( fastest );
					}
					else{
						hSingleStripTime[pl] -> Fill( dt[pl].at(0) );
					}
				}
				stripsFired[pl].clear();

				for( int st = 0 ; st < 48 ; st++ ){
					if(  multStrip[pl][st] != 0)
						hMultStrip[pl][st]->Fill( multStrip[pl][st] );
					for( int si = 0 ; si < 2 ; si++){
						if(  multSide[pl][st][si] != 0)
							hMultSide[pl][st][si]->Fill( multSide[pl][st][si] );
					}
				}
			}
			if( multTotal != 0)
				hMultTot->Fill( multTotal);
			
		}
	}	
	
	//TApplication theApp("App",&argc,argv);

	//theApp.Run();
	
	TFile * outFile = new TFile("multiplicities.root","RECREATE");
	outFile->cd();
	hMultTot->Write();
	for( int pl = 0; pl < 20 ; pl++){
		hMultPlane[pl] -> Write();
		hStripsFired[pl] ->Write();
		hManyStripTime[pl]->Write();
		hManyStripFastTime[pl]->Write();
		hSingleStripTime[pl]->Write();
		/*
 		for( int st = 0 ; st < 48 ; st++){
			hMultStrip[pl][st] -> Write();
			for( int si = 0 ; si < 2 ; si++){
				hMultSide[pl][st][si] -> Write();
			}	
		}
		*/
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

