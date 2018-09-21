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
#include "BmnTOF1Conteiner.h"
#include "BmnMwpcHitFinder.h"
#include "BmnMwpcSegment.h"
#include "BmnMwpcTrackFinder.h"
#include "FairTrackParam.h"
#include "BmnTrack.h"
#include "UniDbRun.h"


// Some constants we will need for analysis
const double pedBC1 = 69.2885;
const double pedBC2 = -11.7212;
const double pedBC3 = -25.4808;
const double pedBC4 = 126.067;
const int run_period = 7;
const double a_out = 0.00173144;
const double b_out = 0.0384856;
const double c_out = 0.000015362;
const double a_in = 0.020542;
const double b_in = 0.0305108;
const double c_in = 0.0000114953;


using namespace std;

double GrabField(TString run_number, const int period);
void checkQC( TString run_number, std::vector<double> *cuts );
void loadT0( TClonesArray *t0, double &time, double &amp);
void skimForCarbon( TClonesArray *bc1Data, TClonesArray *bc2Data, double t0Time, std::vector<double> cuts, bool &pass , double &z2);
void grabZ2( TClonesArray *bc3Data, TClonesArray *bc4Data, double t0Time, double &z2);
void findIdx( TClonesArray* data, int &index , double refT);


int main(int argc, char ** argv)
{

	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tcalcPedestals /path/to/all/digi/files\n";
		return -1;
	}

	////////////////////////////////////////////////////////////////////////////
	// Output file
	TString hName;
	TH1D * hXdiff	= new TH1D("hXdiff","hXdiff",4000,-20,20);
	TH1D * hYdiff	= new TH1D("hYdiff","hYdiff",4000,-20,20);
	TH1D * hZdiff	= new TH1D("hZdiff","hZdiff",4000,-20,20);
	TH1D * hChi2	= new TH1D("hChi2","hChi2",200,0,20);
	TH1D ** hStripMult      = new TH1D*[20];
	TH1D ** hClusterMult	= new TH1D*[20];
	for( int pl = 0 ; pl < 20 ; pl++){
		hName = Form("hStripMult_%i",pl);
		hStripMult[pl]          = new TH1D(hName,hName,40,0,20);
		hName = Form("hClusterMult_%i",pl);
		hClusterMult[pl]        = new TH1D(hName,hName,40,0,20);
	}

	TFile * outFile = new TFile("reconstruction.root","RECREATE");
	TTree * outTree = new TTree("sk","sk");
	TClonesArray * tofHits 		= new TClonesArray("BmnTOF1Conteiner");
	TClonesArray * mwpcSegs 	= new TClonesArray("BmnMwpcSegment");
	TClonesArray * mwpcTracks	= new TClonesArray("BmnTrack");
	double z2_in, z2_out;
	outTree->Branch("tof400"	,&tofHits);
	outTree->Branch("mwpc_sg"	,&mwpcSegs);
	outTree->Branch("mwpc_tr"	,&mwpcTracks);
	outTree->Branch("z2_in"		,&z2_in);
	outTree->Branch("z2_out"	,&z2_out);

	////////////////////////////////////////////////////////////////////////////
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
		double field_voltage = GrabField(run_number,run_period);

		////////////////////////////////////////////////////////////////////////////
		// Load QC file for this run to get carbon position in BC1, BC2
		std::vector<double> carbonInfo;
		checkQC( run_number, &carbonInfo);
		if( carbonInfo.size() == 1) return -3;

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
		TClonesArray * mwpcData = new TClonesArray("BmnMwpcDigit");
		intree->SetBranchAddress("MWPC"		,&mwpcData);
		
		////////////////////////////////////////////////////////////////////////////
		// MWPC Hit Finder		
		cout << "Initializing MWPC Hit Finder...\n";
		BmnMwpcHitFinder* mwpcHM = new BmnMwpcHitFinder(true, 7, atoi(run_number.Data() ) );
		mwpcHM->Init();
		cout << "Initializing MWPC Track Finder...\n";
		BmnMwpcTrackFinder* mwpcTF = new BmnMwpcTrackFinder(true, 7, atoi(run_number.Data() ) );
		mwpcTF->Init();



		const int nEvents = intree->GetEntries();
		cout << "Working on file " << argv[fi+1] << " with " << nEvents << " events and field voltage " << field_voltage << "\n";
		////////////////////////////////////////////////////////////////////////////
		// Loop over all events in file
		for (int event=0 ; event<nEvents ; event++){
			if( (event % 5000) == 0) cout<< "Working on event " << event << "\n";
				// Input branches
			bc1Data->Clear();
			bc2Data->Clear();
			bc3Data->Clear();
			bc4Data->Clear();
			t0Data->Clear();
			tofData->Clear();
			mwpcData->Clear();
				// Output branches
			tofHits->Clear();
			mwpcSegs->Clear();
			mwpcTracks->Clear();
			
			intree->GetEvent(event);
	
			////////////////////////////////////////////////////////////////////////////
			// Demand that event has only 1 T0 TDC digit, otherwise skip event
			if( t0Data->GetEntriesFast() != 1) continue;
			double t0Time, t0Amp;
			loadT0( t0Data, t0Time, t0Amp );

			////////////////////////////////////////////////////////////////////////////
			// Demand that BC1-BC2 is within Carbon Peak
			bool pass = false;
			skimForCarbon( bc1Data, bc2Data, t0Time, carbonInfo, pass , z2_in);
			if (!pass) continue;

			////////////////////////////////////////////////////////////////////////////
			// Look at BC3-BC4 distribution
			grabZ2( bc3Data, bc4Data, t0Time, z2_out);
	

			////////////////////////////////////////////////////////////////////////////
			// Now process ToF400 digits
				// Make sure we clear all previous info of ToF400 detector
			for( int pl = 0; pl < 20; pl++)
				Plane[pl]->ClearHits();

				// Initial skim for only LH side of strip
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				Plane[ signal->GetPlane() ]->InitSkim( signal );
			}

				// Secondary skim to match LH with RH side of strip and create strip hit
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				Plane[ signal->GetPlane() ]->CreateStripHit( signal , t0Time , t0Amp );
			}

			int hitMult = 0;
			for( int pl = 0 ; pl < 20 ; pl ++){
				int stripMult = Plane[pl]->GetStripMult();
				if( stripMult == 0) continue;
				hStripMult[pl]->Fill(stripMult);
					// Cluster the strips in a plane
				Plane[pl]->ClusterHits( );
				int clustMult = Plane[pl]->GetNClusters();
				hitMult += clustMult;
				hClusterMult[pl]->Fill( clustMult );
					// Output hit information to a TClonesArray
				Plane[pl]->OutputToTree( tofHits );
			}
			/*
			for( int hits = 0 ; hits < tofHits->GetEntriesFast() ; hits++){
				BmnTOF1Conteiner * entry = (BmnTOF1Conteiner *)tofHits->At(hits);
				cout << entry->GetXGlobal() << " " << entry->GetYGlobal() << " " << entry->GetZGlobal() << "\n";
				cout << ntry->GetXLocal() << " " << entry->GetYLocal() << " " << entry->GetZLocal() << "\n";
				cout << entry->GetStrip() << " " << entry->GetPlane() << " " << entry->GetTime() << " " << entry->GetAmp() << "\n";
				cout << entry->GetdL() << "\n\n";
			}
			*/

			////////////////////////////////////////////////////////////////////////////
			// Now process MWPC digits
			mwpcHM->Exec("",mwpcData,mwpcSegs, event);
			mwpcTF->Exec("",mwpcSegs,mwpcTracks);
			for( int tr = 0 ; tr < mwpcTracks->GetEntriesFast() ; tr++){
				BmnMwpcTrack * thisTr = (BmnMwpcTrack *) mwpcTracks->At(tr);
				FairTrackParam * par = (FairTrackParam *) thisTr->GetParamFirst();
				
				//cout << par->GetTx() << " " << par->GetTy() << "\n";
			}





				// Fill output tree with the TClonesArray info from ToF400
			outTree->Fill();

		} // end of loop over events in file
	} // end of loop over files


	outFile->cd();
	hXdiff->Write();
	hYdiff->Write();
	hZdiff->Write();
	hChi2->Write();
	for( int pl = 0; pl < 20; pl++){
		//hStripMult[pl]->Write();
		//hClusterMult[pl]->Write();
	}
	outTree->Write();
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

void skimForCarbon( TClonesArray *bc1Data, TClonesArray *bc2Data, double t0Time, std::vector<double> cuts, bool &pass , double &z2){
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
		
	double x =  sqrt( adcBC1 * adcBC2 );
	z2 = a_in + b_in*x + c_in*x*x;

	double BC1_CPeak = cuts.at(0);
	double BC1_CWidth = cuts.at(2);
	double BC2_CPeak = cuts.at(1);
	double BC2_CWidth = cuts.at(3);

	

	if(  ( fabs( adcBC1 - (BC1_CPeak-pedBC1) ) < 2*BC1_CWidth ) && ( fabs( adcBC2 - (BC2_CPeak-pedBC2) ) < 2*BC2_CWidth )   )
		pass = true;
	
}

void grabZ2( TClonesArray *bc3Data, TClonesArray *bc4Data, double t0Time, double &z2){
	double x2 = -1;
	if( bc3Data->GetEntriesFast() && bc4Data->GetEntriesFast() ){
		int bc3Idx;
		findIdx(bc3Data,bc3Idx,t0Time);
		int bc4Idx;
		findIdx(bc4Data,bc4Idx,t0Time);

		BmnTrigWaveDigit * signalBC3 = (BmnTrigWaveDigit*) bc3Data->At(bc3Idx);
		BmnTrigWaveDigit * signalBC4 = (BmnTrigWaveDigit*) bc4Data->At(bc4Idx);

		x2 = sqrt( (signalBC3->GetPeak() - pedBC3)*(signalBC4->GetPeak() - pedBC4) );
		z2 = a_out + b_out*x2 + c_out*x2*x2;
	}
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

