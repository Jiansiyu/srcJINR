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

#include "BmnNewFieldMap.h"
#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"
#include "BmnTof1Digit.h"
#include "BmnTOF1Detector.h"
#include "BmnTOF1Conteiner.h"
#include "BmnMwpcGeometrySRC.h"
#include "BmnMwpcHitFinder.h"
#include "BmnMwpcSegment.h"
#include "BmnMwpcTrackFinder.h"
#include "FairTrackParam.h"
#include "BmnTrack.h"
#include "UniDbRun.h"
#include "BmnGemStripHit.h"
#include "BmnGemStripConfiguration.h"
#include "BmnGemStripHitMaker.h"


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
double randomStripDiff();

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
	TH2D ** hWireEvent	= new TH2D*[4];
	TH2D ** hWireEventFast	= new TH2D*[4];

	for( int st = 0 ; st < 4 ; st++){
		hName = Form("hWireEvent_%i",st);
		hWireEvent[st]		= new TH2D(hName,hName,2000,0.5,2000.5,600,0.5,600.5);
		hName = Form("hWireEventFast_%i",st);
		hWireEventFast[st]	= new TH2D(hName,hName,2000,0.5,2000.5,600,0.5,600.5);
	}

	TFile * outFile = new TFile("reconstruction.root","RECREATE");
	TTree * outTree = new TTree("sk","sk");
	TClonesArray * tofHits 		= new TClonesArray("BmnTOF1Conteiner");
	TClonesArray * mwpcSegs 	= new TClonesArray("BmnMwpcSegment"); // all segments
	TClonesArray * mwpc1Tracks	= new TClonesArray("BmnTrack"); // for chamber 1
	TClonesArray * mwpc2Tracks	= new TClonesArray("BmnTrack"); // for chamber 2
	TClonesArray * mwpc3Tracks	= new TClonesArray("BmnTrack"); // for chamber 3
	TClonesArray * mwpc4Tracks	= new TClonesArray("BmnTrack"); // for chamber 4
	TClonesArray * mwpcTracks	= new TClonesArray("BmnTrack"); // all tracks
	TClonesArray * gemHits		= new TClonesArray("BmnGemStripHit");

	double z2_in, z2_out;
	outTree->Branch("tof400"	,&tofHits);
	outTree->Branch("mwpc_sg"	,&mwpcSegs); // All segments
	outTree->Branch("mwpc_tr1"	,&mwpc1Tracks); // Aligned segments in ch 1
	outTree->Branch("mwpc_tr2"	,&mwpc2Tracks); // Aligned in ch 2
	outTree->Branch("mwpc_tr3"	,&mwpc3Tracks); // etc..
	outTree->Branch("mwpc_tr4"	,&mwpc4Tracks);
	outTree->Branch("mwpc_tr"	,&mwpcTracks); // Real 'tracks'
	outTree->Branch("z2_in"		,&z2_in);
	outTree->Branch("z2_out"	,&z2_out);
	outTree->Branch("gem"		,&gemHits); // For gems

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
	// Init MWPC geometry
	BmnMwpcGeometrySRC * mwpcGeo = new BmnMwpcGeometrySRC(7, 3430);
	TVector3 * ChCent = new TVector3[4];
	double ZCh_cent[4] = { 0. };
	double mwpc_shift[4][4] = { 0. };
	for( int i = 0 ; i < 4 ; i++){
		ChCent[i] = mwpcGeo -> GetChamberCenter(i);
		ZCh_cent[i] = ChCent[i].Z();
		mwpc_shift[i][0] = mwpcGeo -> GetTx(i);
		mwpc_shift[i][1] = ChCent[i].X();
		mwpc_shift[i][2] = mwpcGeo -> GetTy(i); 
		mwpc_shift[i][3] = ChCent[i].Y();
	}
	const double dw = mwpcGeo -> GetWireStep();
	const double midPl = 47.25;

	////////////////////////////////////////////////////////////////////////////
	// Init GEM geometry configuration for SRC run
	BmnGemStripConfiguration::GEM_CONFIG gem_config = BmnGemStripConfiguration::RunSRCSpring2018;


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
		// Grab magnetic field
		double field_voltage = GrabField(run_number,run_period);
		BmnNewFieldMap * magField = new BmnNewFieldMap("field_sp41v4_ascii_Extrap.root");
		double map_current = 55.87;
		double field_scale;
		bool isField = true;
		if( field_voltage < 10 ){
			field_scale = 0.;
			isField = false;
		}
		else{ field_scale = (field_voltage) / map_current; }
		magField->SetScale( field_scale );
		magField->Init();

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
		TClonesArray * gemData = new TClonesArray("BmnGemStripDigit");
		intree->SetBranchAddress("GEM"		,&gemData);
		
		////////////////////////////////////////////////////////////////////////////
		// MWPC Initialization for this Run		
		cout << "Initializing MWPC Hit Finder...\n";
		BmnMwpcHitFinder* mwpcHM = new BmnMwpcHitFinder(true, 7, atoi(run_number.Data() ) );
		mwpcHM->Init();
		cout << "Initializing MWPC Track Finder...\n";
		BmnMwpcTrackFinder* mwpcTF = new BmnMwpcTrackFinder(true, 7, atoi(run_number.Data() ) );
		mwpcTF->Init();

		////////////////////////////////////////////////////////////////////////////
		// GEM Initalization for this Run
		cout << "Initializing GEM Hit Maker...\n";
		BmnGemStripHitMaker * gemHM = new BmnGemStripHitMaker(7, true);
		gemHM->SetCurrentConfig( gem_config );
		gemHM->SetHitMatching( true );
		gemHM->Init( magField );


		const int nEvents = intree->GetEntries();
		cout << "Working on file " << argv[fi+1] << " with " << nEvents << " events and field voltage " << field_voltage << "\n";
		////////////////////////////////////////////////////////////////////////////
		// Loop over all events in file
		for (int event=0 ; event<nEvents ; event++){
			if( (event % 10000) == 0) cout<< "Working on event " << event << "\n";
				// Input branches
			bc1Data->Clear();
			bc2Data->Clear();
			bc3Data->Clear();
			bc4Data->Clear();
			t0Data->Clear();
			tofData->Clear();
			mwpcData->Clear();
			gemData->Clear();
				// Output branches
			tofHits->Clear();
			mwpcSegs->Clear();
			mwpc1Tracks->Clear();
			mwpc2Tracks->Clear();
			mwpc3Tracks->Clear();
			mwpc4Tracks->Clear();
			mwpcTracks->Clear();
			gemHits->Clear();
			
			
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
				//hStripMult[pl]->Fill(stripMult);
					// Cluster the strips in a plane
				Plane[pl]->ClusterHits( );
				int clustMult = Plane[pl]->GetNClusters();
				hitMult += clustMult;
				//hClusterMult[pl]->Fill( clustMult );
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

			
			if( z2_out <= 0 ) continue;
			if( event > 2000) break;
			////////////////////////////////////////////////////////////////////////////
			// Now process MWPC digits
			int flagHit[4][6][96] = { 0 };
			int nWiresPl[4][6] = { 0 };
			int nWiresPlBeam[4][6] = { 0 };
			int nWiresPlOutBeam[4][6] = { 0 };
			std::vector< double > coords[4][6];
			std::vector< int > wires[4][6];
			unsigned int times[4][6][96] = { 9999 };
			for( int en = 0 ; en < mwpcData -> GetEntriesFast() ; en++){
				BmnMwpcDigit * signal = (BmnMwpcDigit*)mwpcData->At(en);
				short st = signal->GetStation();
				int wire = signal->GetWireNumber();
				short pl = signal->GetPlane();
				unsigned int t = t0Time - signal->GetTime();

				// If wire fires multiple times, take the earliest timestamp
				if( t < times[st][pl][wire] )
					times[st][pl][wire] = t;
				
				// Cut physically around the beam position -- TODO: fix for all planes
				if( fabs( wire - 35 ) > 10 ) continue; // only valid for st 0, pl 0 right now
	
				wires[st][pl].push_back( wire );
				// Build multiplicity for how many wires fire in a single plane
				// and save which wires fired
				if( flagHit[st][pl][wire] == 0 ){
					nWiresPl[st][pl]++;
				}
				// Flag that a wire fired
				flagHit[st][pl][wire]++;

			}

			
			for( int st = 0 ; st < 4 ; st++){
				for( int pl = 0 ; pl < 6 ; pl++){

					// Cut on multiplicity of a plane:
					//if( nWiresPl[st][pl] > 7 ) continue;

					// Now let's look at how these events look like:
					std::vector<int> earlyWire;
					unsigned int earlyT = 9999;
					for( int wr = 0 ; wr < wires[st][pl].size() ; wr ++ ){
						hWireEvent[st] -> Fill( event, wires[st][pl].at(wr) + 96*pl );
						if( earlyT  > times[st][pl][ wires[st][pl].at(wr) ] ){
							earlyT = times[st][pl][ wires[st][pl].at(wr)];
							earlyWire.push_back( wires[st][pl].at(wr) );
						}
						else if( earlyT = times[st][pl][ wires[st][pl].at(wr) ] ){
							earlyWire.push_back( wires[st][pl].at(wr) );
						}
					}
					for( int wr = 0 ; wr < earlyWire.size() ; wr ++ )
						hWireEventFast[st] -> Fill( event, earlyWire.at(wr) + 96*pl );

				}// End pl loop
			}// End st loop
				
			
			//mwpcHM->Exec("",mwpcData,mwpcSegs, event);
				// Since I'm debugging segment matching, create "tracks" that are just segments in each chamber and align them and save to skim tree	

			////////////////////////////////////////////////////////////////////////////
			// Now process MWPC tracking -- the shitty one...
			//mwpcTF->Exec("",mwpcSegs,mwpcTracks);


			////////////////////////////////////////////////////////////////////////////
			// Process GEM digits
			// Issue currently with SRC XML file with no module 1 for stations 0-3
			//gemHM->Exec( "" , gemData  , gemHits , NULL , magField );
			




				// Fill output tree with the TClonesArray info from ToF400
			outTree->Fill();

		} // end of loop over events in file
	} // end of loop over files

	cout << "\tClosing and writing histograms/trees...\n";

	outFile->cd();
	for( int st = 0 ; st < 4 ; st++){
		hWireEvent[st]	-> Write();
		hWireEventFast[st]	-> Write();
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
double randomStripDiff(){
	TRandom3 * rand = new TRandom3(0);
	double dist = 0.;
	while( dist == 0.){
		int st1 = rand->Rndm() * 96;
		int st2 = rand->Rndm() * st1; // lower range possible
		int st3 = rand->Rndm() * (96-st1) + st1; // upper range possible
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

