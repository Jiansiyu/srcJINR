// Using BmnTOF1Detector class to store everything rather than in one large macro
// and then want to check against tofAnalysis.cpp to make sure same output

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <ctime>

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


// Some constants we will need for analysis
const double pedBC1 = 69.2885;
const double pedBC2 = -11.7212;
const double pedBC3 = -25.4808;
const double pedBC4 = 126.067;
const int run_period = 7;



using namespace std;

double GrabField(TString run_number, const int period);
void checkQC( TString run_number, std::vector<double> *cuts );
void loadT0( TClonesArray *t0, double &time, double &amp);
void skimForCarbon( TClonesArray *bc1Data, TClonesArray *bc2Data, double t0Time, std::vector<double> cuts, bool &pass );
void findIdx( TClonesArray* data, int &index , double refT);

int main(int argc, char ** argv)
{
	std::clock_t totalStart, classInitTime, voltageTime, qcTime, timer;
	double loadEvent = 0, loadT0time = 0, skimCarbon_time = 0, tofClear_time = 0, tofInitSkim_time = 0, tofCreateHit_time = 0;
	totalStart = std::clock();

	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tcalcPedestals /path/to/all/digi/files\n";
		return -1;
	}

	TString hName;
	TH1D ** hStripMult      = new TH1D*[20];
	for( int pl = 0 ; pl < 20 ; pl++){
		hName = Form("hStripMult_%i",pl);
		hStripMult[pl]          = new TH1D(hName,hName,40,0,20);
	}


	////////////////////////////////////////////////////////////////////////////
	// Init BmnTOF1Detector and load calibration files and geometry
	classInitTime = std::clock();
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
	cout << "Class init timer: " << ( std::clock() - classInitTime ) / (double) CLOCKS_PER_SEC << "\n"; 

	////////////////////////////////////////////////////////////////////////////
	// Loop over files in arguments provided
	const int files = argc - 1;
	for( int fi = 0 ; fi < files ; ++fi){

		////////////////////////////////////////////////////////////////////////////
		// Initalize input file
		TFile * infile = NULL;
		infile = new TFile(argv[fi+1]);
		if (infile->IsZombie()){
			cerr << "Could not open file " << argv[fi+1] <<"\n"
			        << "\tBailing out\n";
			return -2;
		}
		
		////////////////////////////////////////////////////////////////////////////
		// Get the run number from input file:
		TString file = argv[fi+1];
		TString run_number( file(file.Index(".")-9,4) );
		
		////////////////////////////////////////////////////////////////////////////
		// Grab magnetic field current for storing different histograms
		voltageTime = std::clock();
		double field_voltage = GrabField(run_number,run_period);
		cout << "Votlage grab timer: " << ( std::clock() - voltageTime ) / (double) CLOCKS_PER_SEC << "\n"; 
		

		////////////////////////////////////////////////////////////////////////////
		// Load QC file for this run to get carbon position in BC1, BC2
		qcTime = std::clock();
		std::vector<double> carbonInfo;
		checkQC( run_number, &carbonInfo);
		if( carbonInfo.size() == 1) return -3;
		cout << "QC Processor timer: " << ( std::clock() - qcTime ) / (double) CLOCKS_PER_SEC << "\n"; 

		////////////////////////////////////////////////////////////////////////////
		// Setup intree for analysis
		TTree * intree = NULL;
		intree = (TTree*) infile->Get("cbmsim");
		if (!intree){
			cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";
			return -4;
		}
		TClonesArray * bc1Data  = new TClonesArray("BmnTrigWaveDigit");
		intree->SetBranchAddress("TQDC_BC1"     ,&bc1Data);
		TClonesArray * bc2Data  = new TClonesArray("BmnTrigWaveDigit");
		intree->SetBranchAddress("TQDC_T0"     ,&bc2Data);
		TClonesArray * bc3Data  = new TClonesArray("BmnTrigWaveDigit");
		intree->SetBranchAddress("TQDC_BC3"     ,&bc3Data);
		TClonesArray * bc4Data  = new TClonesArray("BmnTrigWaveDigit");
		intree->SetBranchAddress("TQDC_BC4"     ,&bc4Data);

		TClonesArray * t0Data   = new TClonesArray("BmnTrigDigit");
		intree->SetBranchAddress("T0"           ,&t0Data);
		TClonesArray * tofData  = new TClonesArray("BmnTof1Digit");
		intree->SetBranchAddress("TOF400"       ,&tofData);
		
		const int nEvents = intree->GetEntries();
		cout << "Working on file " << argv[fi+1] << " with " << nEvents << " events and field voltage " << field_voltage << "\n";
		
		////////////////////////////////////////////////////////////////////////////
		// Loop over all events in file
		for (int event=0 ; event<nEvents ; event++){
			timer = std::clock();
			bc1Data->Clear();
			bc2Data->Clear();
			bc3Data->Clear();
			bc4Data->Clear();
			
			intree->GetEvent(event);
			loadEvent+= ( std::clock() - timer );
	
			////////////////////////////////////////////////////////////////////////////
			// Demand that event has only 1 T0 TDC digit, otherwise skip event
			if( t0Data->GetEntriesFast() != 1) continue;
			timer = std::clock();
			double t0Time, t0Amp;
			loadT0( t0Data, t0Time, t0Amp );
			loadT0time+= ( std::clock() - timer );

			////////////////////////////////////////////////////////////////////////////
			// Demand that BC1-BC2 is within Carbon Peak
			timer = std::clock();
			bool pass;
			skimForCarbon( bc1Data, bc2Data, t0Time, carbonInfo, pass );
			if (!pass) continue;
			skimCarbon_time+= ( std::clock() - timer );

			////////////////////////////////////////////////////////////////////////////
			// Now process ToF400 digits
				
				// Make sure we clear all previous info of ToF400 detector
			timer = std::clock();
			for( int pl = 0; pl < 20; pl++)
				Plane[pl]->ClearHits();
			tofClear_time+= ( std::clock() - timer );

				// Initial skim for only LH side of strip
			timer = std::clock();
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				Plane[ signal->GetPlane() ]->InitSkim( signal );
			}
			tofInitSkim_time+= ( std::clock() - timer );

				// Secondary skim to match LH with RH side of strip and create strip hit
			timer = std::clock();
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				Plane[ signal->GetPlane() ]->CreateStripHit( signal , t0Time , t0Amp );
			}
			tofCreateHit_time+= ( std::clock() - timer );
			
			for( int pl = 0 ; pl < 20 ; pl ++){
				int stripMult = Plane[pl]->GetStripMult();
				if( stripMult > 0)	 hStripMult[pl]->Fill(stripMult);
				Plane[pl]->ClusterHits( );
			}
	

		} // end of loop over events in file


	} // end of loop over files

	double totTime = ( std::clock() - totalStart ) / (double) CLOCKS_PER_SEC;
	loadEvent = loadEvent / (double) CLOCKS_PER_SEC;
	loadT0time = loadT0time / (double) CLOCKS_PER_SEC;
	skimCarbon_time = skimCarbon_time / (double) CLOCKS_PER_SEC;
	tofClear_time = tofClear_time / (double) CLOCKS_PER_SEC;
	tofInitSkim_time = tofInitSkim_time / (double) CLOCKS_PER_SEC;
	tofCreateHit_time = tofCreateHit_time / (double) CLOCKS_PER_SEC;

	cout << "Total Time: " << totTime << "\n"
		<< "\t" << loadEvent << " " << loadT0time << " " << skimCarbon_time << "\n\t" << 
		tofClear_time << " " << tofInitSkim_time << " " << tofCreateHit_time << "\n";
	

	TFile * outFile = new TFile("CLASStofevents.root","RECREATE");
	outFile->cd();
	for( int pl = 0; pl < 20; pl++)
		hStripMult[pl]->Write();
	outFile->Write();
	outFile->Close();
	

	return 0;
}

double GrabField(TString run_number, const int period){
	int runNo = atoi( run_number.Data() );
	UniDbRun* pCurrentRun = UniDbRun::GetRun(period,runNo);
	double volt;
	if( pCurrentRun == 0 ){
		if( (runNo > 3474 && runNo < 3485) || (runNo > 3434 && runNo < 3443) )
		        volt = 87.0;
		else if( (runNo > 3484 && runNo < 3496) || (runNo > 3442 && runNo < 3453) || (runNo > 3513) )
		        volt = 107.7;
		else if( (runNo > 3495 && runNo < 3501) || (runNo > 3452 && runNo < 3458) )
		        volt = 123.7;
	}
	else{
		volt = *(pCurrentRun->GetFieldVoltage());
	}	


	return volt;
}

void checkQC( TString run_number, std::vector<double> *cuts ){
	TFile * qualityFile = NULL;
	TString path = std::getenv("VMCWORKDIR");
	path = path + "/build/bin/qualityCheck/checked_" + run_number + ".root";
	qualityFile = new TFile(path);
	if( qualityFile->IsZombie() ){
		cerr << "No checked file for this run number. You need to run calcPedestals first.\n"
		        << "\tBailing...\n";
		cuts->push_back(-1.);
		return;
	}
	TVectorT<double> * carbonIn     = (TVectorT<double>*)qualityFile->Get("carbonIn");
	TVectorT<double> * carbonInWidth= (TVectorT<double>*)qualityFile->Get("carbonInWidth");
	
	cuts->push_back( (*carbonIn)[0] );
	cuts->push_back( (*carbonIn)[1] );
	cuts->push_back( (*carbonInWidth)[0] );
	cuts->push_back( (*carbonInWidth)[1] );

	return;
}

void loadT0( TClonesArray *t0, double &time, double &amp){
	BmnTrigDigit * signal = (BmnTrigDigit*) t0->At(0);
	time = signal->GetTime();
	amp = signal->GetAmp();
}

void skimForCarbon( TClonesArray *bc1Data, TClonesArray *bc2Data, double t0Time, std::vector<double> cuts, bool &pass ){
	double adcBC1 = -1., adcBC2 = -1.;

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

	double BC1_CPeak = cuts.at(0);
	double BC1_CWidth = cuts.at(2);
	double BC2_CPeak = cuts.at(1);
	double BC2_CWidth = cuts.at(3);
	if(  ( fabs( adcBC1 - (BC1_CPeak-pedBC1) ) < 2*BC1_CWidth ) && ( fabs( adcBC2 - (BC2_CPeak-pedBC2) ) > 2*BC2_CWidth )   )
		pass = true;
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

