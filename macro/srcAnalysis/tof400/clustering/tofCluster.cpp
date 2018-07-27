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
double randomStripDiff();
double randomPosDiff();
double fixWalk( double amp, double shift, double par0, double par1, double par2);

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

	// Try opening the timewalk file if it exists:
	ifstream f_corr_walk;
	dir = std::getenv("VMCWORKDIR");
	dir += "/input/TOF400_TimeWalk_RUN7_SRC.dat";
	cout << "Attempting to open TimeWalk file: " << dir << "\n";
	f_corr_walk.open(dir);
	char line3[256];
	f_corr_walk.getline(line3, 256);
	f_corr_walk.getline(line3, 256);
	int tmp_Pl, tmp_St, tmp_Pt;
	double tmp_Sh, tmp_p0, tmp_p1, tmp_p2, tmp_p0e, tmp_p1e, tmp_p2e;
	double walkFunc[20][48][4] = {0.};
	if (f_corr_walk.is_open() == true){
 		while (!f_corr_walk.eof()) {
			f_corr_walk >> tmp_Pl >> tmp_St >> tmp_Pt \
				>> tmp_Sh >> tmp_p0 >> tmp_p1 >> tmp_p2 \
				>> tmp_p0e >> tmp_p1e >> tmp_p2e;
			walkFunc[tmp_Pl][tmp_St][0] = tmp_Sh;
			walkFunc[tmp_Pl][tmp_St][1] = tmp_p0;
			walkFunc[tmp_Pl][tmp_St][2] = tmp_p1;
			walkFunc[tmp_Pl][tmp_St][3] = tmp_p2;
			if( tmp_Sh == -1)
				killStrip[tmp_Pl][tmp_St] = true;
		}
		cout << "\tLoaded time walk file\n";
	}
	else{
		cout << "\tFailed to find time walk file, setting all corrections to 0...\n";
	} 
	
	// Histograms for looking at cuts
	TH1D *** hSingle_meanTime 	= new TH1D**[20];	// single hits, mean time
	TH1D *** hMany_meanTime 	= new TH1D**[20];	// single hits, mean time
	TH1D *** hSingle_amp		= new TH1D**[20];	// single hits, amplitude
	TH1D *** hMany_amp		= new TH1D**[20];	// single hits, amplitude

	TH1D ** hMany_meanTime_diff	= new TH1D*[20];	// many strips, mean time diff
	TH1D ** hMany_amp_diff		= new TH1D*[20];	// many strips, amp diff
	TH1D ** hMany_StrSpread_diff	= new TH1D*[20];	// many strips, strip # diff
	TH1D ** hMany_StrSpread_diff_MC	= new TH1D*[20];	// many strips, strip # diff MC random
	TH1D ** hMany_YSpread_diff	= new TH1D*[20];	// many strips, strip y spread
	TH1D ** hMany_YSpread_diff_MC	= new TH1D*[20];	// many strips, strip y spread MC random

	TH2D ** hMany_2DSpread		= new TH2D*[20];
	TH2D ** hMany_2D_ampTime	= new TH2D*[20];
	TH2D ** hMany_2D_ampStrDiff	= new TH2D*[20];
	TH2D ** hMany_2D_ampYDiff	= new TH2D*[20];
	TH2D ** hMany_2D_timeStrDiff	= new TH2D*[20];
	TH2D ** hMany_2D_timeYDiff	= new TH2D*[20];
	
	TString hName;
	for( int pl = 0; pl < 20 ; pl++){
		hName	= Form("hMany_meanTime_diff_%i",pl);
		hMany_meanTime_diff[pl]		= new TH1D(hName,hName,8000,-100,100);

		hName	= Form("hMany_amp_diff_%i",pl);
		hMany_amp_diff[pl]		= new TH1D(hName,hName,4000,-100,100);

		hName	= Form("hMany_StrSpread_diff_%i",pl);
		hMany_StrSpread_diff[pl]	= new TH1D(hName,hName,96,0,96);

		hName	= Form("hMany_StrSpread_diff_MC_%i",pl);
		hMany_StrSpread_diff_MC[pl]	= new TH1D(hName,hName,96,0,96);
		
		hName	= Form("hMany_YSpread_diff_%i",pl);
		hMany_YSpread_diff[pl]		= new TH1D(hName,hName,2400,-60,60);

		hName	= Form("hMany_YSpread_diff_MC_%i",pl);
		hMany_YSpread_diff_MC[pl]	= new TH1D(hName,hName,2400,-60,60);
		
		hName	= Form("hMany_2DSpread_%i",pl);
		hMany_2DSpread[pl]		= new TH2D(hName,hName,96,0,96,2400,-60,60);

		hName	= Form("hMany_2D_ampTime_%i",pl);
		hMany_2D_ampTime[pl]		= new TH2D(hName,hName,4000,-100,100,8000,-100,100);
		
		hName	= Form("hMany_2D_ampStrDiff_%i",pl);
		hMany_2D_ampStrDiff[pl]		= new TH2D(hName,hName,4000,-100,100,96,0,96);

		hName	= Form("hMany_2D_ampYDiff_%i",pl);
		hMany_2D_ampYDiff[pl]		= new TH2D(hName,hName,4000,-100,100,2400,-60,60);

		hName	= Form("hMany_2D_timeStrDiff_%i",pl);
		hMany_2D_timeStrDiff[pl]	= new TH2D(hName,hName,8000,-100,100,96,0,96);
		
		hName	= Form("hMany_2D_timeYDiff_%i",pl);
		hMany_2D_timeYDiff[pl]		= new TH2D(hName,hName,8000,-100,100,2400,-60,60);

		hSingle_meanTime[pl]		= new TH1D*[48];
		hMany_meanTime[pl]		= new TH1D*[48];
		hSingle_amp[pl]			= new TH1D*[48];
		hMany_amp[pl]			= new TH1D*[48];
		
		for( int st = 0 ; st < 48 ; st++){
			hName	= Form("hSingle_meanTime_%i_%i",pl,st);
			hSingle_meanTime[pl][st]	= new TH1D(hName,hName,4000,-15,85);
			hName	= Form("hMany_meanTime_%i_%i",pl,st);
			hMany_meanTime[pl][st]		= new TH1D(hName,hName,4000,-15,85);
			hName	= Form("hSingle_amp_%i_%i",pl,st);
			hSingle_amp[pl][st]		= new TH1D(hName,hName,4000,-100,100);
			hName	= Form("hMany_amp_%i_%i",pl,st);
			hMany_amp[pl][st]		= new TH1D(hName,hName,4000,-100,100);
		}
	}
	// Output tree
	TFile * outFile = new TFile("clusteringPlots.root","RECREATE");


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

			TOFHit tofEvents[20][48];
			std::vector<double> hitInfoTime[20][48];
			std::vector<double> hitInfoAmps[20][48];
			std::vector<double> hitInfoTimeR[20][48];
			std::vector<double> hitInfoAmpsR[20][48];

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

					hitInfoTimeR[plane][strip].push_back( tofT );
                                        hitInfoAmpsR[plane][strip].push_back( tofAmp );
					
					for( int hit = 0 ; hit < hitInfoTime[plane][strip].size() ; hit++){
						double tmpT = hitInfoTime[plane][strip].at(hit);
						double tmpA = hitInfoAmps[plane][strip].at(hit);

						if( fabs( (1./0.06)*(tmpT - tofT + corrLR[plane][strip]) ) < 32 ){ // Create a valid strip hit
							double meanTime = (tmpT + tofT + corrLR[plane][strip])*0.5;
							double sumAmps = sqrt(tmpA * tofAmp);
							double pos = 0.5*(1./0.06)*(tmpT - tofT + corrLR[plane][strip]); // +/- 15cm
		
							if( sumAmps < 9 ) continue; // I don't want hits with low amplitude					
	
							// If there was already a created hit for this strip, only take the earliest one
							// so there will only be one hit per strip allowed
							if( tofEvents[plane][strip].hit == false || (tofEvents[plane][strip].hit == true && meanTime < tofEvents[plane][strip].t) ){
								double par0 = walkFunc[plane][strip][1];
								double par1 = walkFunc[plane][strip][2];
								double par2 = walkFunc[plane][strip][3];
								double shift = walkFunc[plane][strip][0];
								tofEvents[plane][strip].t = meanTime - t0Time + (-6.1 + 27./sqrt(t0Amp) ) - stripShift[plane][strip] - fixWalk(sumAmps, shift, par0, par1, par2);
	
								tofEvents[plane][strip].a = sumAmps;
								tofEvents[plane][strip].y = pos;
								tofEvents[plane][strip].hit = true;
							}
						}
					}
				}
			}
	
			// Now for each event, I've collected at most 1 hit per strip in the ToF400. Of course not all strips fired, and not all planes fired.			

			// Now fill histograms that I have created valid strip hits.
			// 	Look at how many strips fired in a trigger event
			//
			for( int pl = 0 ; pl < 20 ; pl++){
				int cnt = 0;
				std::vector<int> fired;
				for( int st = 0 ; st < 48 ; st++)
					if( tofEvents[pl][st].hit == true){
						fired.push_back(st);
						cnt++;
				}

			
				if( cnt == 0) {} // 0 strip hits in a plane
				if( cnt == 1) { // 1 strip hits in a plane
					hSingle_meanTime[pl][fired.at(0)]	-> Fill( tofEvents[pl][fired.at(0)].t	);
					hSingle_amp[pl][fired.at(0)]		-> Fill( tofEvents[pl][fired.at(0)].a	);
		
				}
				if( cnt == 2) { // 2 strip hits in a plane
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y) < 10. && fabs(fired.at(0) - fired.at(1) ) < 15. ){
						hMany_meanTime[pl][fired.at(0)]		-> Fill( tofEvents[pl][fired.at(0)].t	);
						hMany_meanTime[pl][fired.at(1)]		-> Fill( tofEvents[pl][fired.at(1)].t	);
						hMany_amp[pl][fired.at(0)]		-> Fill( tofEvents[pl][fired.at(0)].a   );
						hMany_amp[pl][fired.at(1)]		-> Fill( tofEvents[pl][fired.at(1)].a   );

						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(1) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(1) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a ,-fired.at(0) + fired.at(1) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ,-fired.at(0) + fired.at(1) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						/*
						//if( (fabs( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ) > 5 && fabs( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ) < 8) || (fabs( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ) > 12) ){
						if( (fabs(tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a - 6) < 2 )){
							cout 	<< pl << " " << fired.at(0) << " " << fired.at(1) << " "
								<< tofEvents[pl][fired.at(0)].t << " " << tofEvents[pl][fired.at(1)].t << " " 
								<< tofEvents[pl][fired.at(0)].a << " " << tofEvents[pl][fired.at(1)].a << " " 
								<< tofEvents[pl][fired.at(0)].y << " " << tofEvents[pl][fired.at(1)].y << "\n";
							// 1st strip characteristics:
							cout << fired.at(0) << "\n";
							for( int hitL = 0 ; hitL < hitInfoTime[pl][fired.at(0)].size() ; hitL++){
								for( int hitR = 0 ; hitR < hitInfoTimeR[pl][fired.at(0)].size() ; hitR++) {
									cout 	<< hitL << " " << hitR << " " 
										<< hitInfoTime[pl][fired.at(0)].at(hitL) << " " << hitInfoTimeR[pl][fired.at(0)].at(hitR) << " "
										<< corrLR[pl][fired.at(0)] << " " 
										<< 0.5*(hitInfoTime[pl][fired.at(0)].at(hitL) + hitInfoTimeR[pl][fired.at(0)].at(hitR) + corrLR[pl][fired.at(0)]) << " " 
										<< hitInfoAmps[pl][fired.at(0)].at(hitL) << " " << hitInfoAmpsR[pl][fired.at(0)].at(hitR) << " " 
										<< sqrt( hitInfoAmps[pl][fired.at(0)].at(hitL) * hitInfoAmpsR[pl][fired.at(0)].at(hitR) ) << "\n";
								}
							}
							cout << fired.at(1) << "\n";
							for( int hitL = 0 ; hitL < hitInfoTime[pl][fired.at(1)].size() ; hitL++){
								for( int hitR = 0 ; hitR < hitInfoTimeR[pl][fired.at(1)].size() ; hitR++) {
									cout 	<< hitL << " " << hitR << " " 
										<< hitInfoTime[pl][fired.at(1)].at(hitL) << " " << hitInfoTimeR[pl][fired.at(1)].at(hitR) << " "
										<< corrLR[pl][fired.at(1)] << " " 
										<< 0.5*(hitInfoTime[pl][fired.at(1)].at(hitL) + hitInfoTimeR[pl][fired.at(1)].at(hitR) + corrLR[pl][fired.at(1)]) << " " 
										<< hitInfoAmps[pl][fired.at(1)].at(hitL) << " " << hitInfoAmpsR[pl][fired.at(1)].at(hitR) << " " 
										<< sqrt( hitInfoAmps[pl][fired.at(1)].at(hitL) * hitInfoAmpsR[pl][fired.at(1)].at(hitR) ) << "\n";
								}
							}
							cout << "\n";
						}
						*/
					//
					
				}
				if( cnt == 3) { // 3 strip hits in a plane
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y) < 10. && fabs(fired.at(0) - fired.at(1) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(1) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(1) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a ,-fired.at(0) + fired.at(1) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t , -fired.at(0) + fired.at(1) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y) < 10. && fabs(fired.at(0) - fired.at(2) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(2) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(2) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a ,-fired.at(0) + fired.at(2) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t , -fired.at(0) + fired.at(2) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y) < 10. && fabs(fired.at(1) - fired.at(2) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(1) + fired.at(2) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(1) + fired.at(2) , tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a ,-fired.at(1) + fired.at(2) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t ,-fired.at(1) + fired.at(2) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
					//}
					
	
				}
				if( cnt == 4) { // 4 strip hits in a plane
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y) < 10. && fabs(fired.at(0) - fired.at(1) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(1) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(1) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a ,-fired.at(0) + fired.at(1) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(1)].a ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ,-fired.at(0) + fired.at(1) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(1)].t ,tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(1)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y) < 10. && fabs(fired.at(0) - fired.at(2) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(2) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(2) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a ,-fired.at(0) + fired.at(2) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t ,-fired.at(0) + fired.at(2) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(2)].t , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(2)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(3)].y) < 10. && fabs(fired.at(0) - fired.at(3) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(3)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(3)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(0) + fired.at(3) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(3)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(0) + fired.at(3) , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(3)].a , tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(3)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(3)].a ,-fired.at(0) + fired.at(3) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].a - tofEvents[pl][fired.at(3)].a , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(3)].t ,-fired.at(0) + fired.at(3) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(0)].t - tofEvents[pl][fired.at(3)].t , tofEvents[pl][fired.at(0)].y - tofEvents[pl][fired.at(3)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y) < 10. && fabs(fired.at(1) - fired.at(2) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(1) + fired.at(2) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(1) + fired.at(2) , tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a , tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a ,-fired.at(1) + fired.at(2) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(2)].a ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t ,-fired.at(1) + fired.at(2) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(2)].t ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(2)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(3)].y) < 10. && fabs(fired.at(1) - fired.at(3) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(3)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(3)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(1) + fired.at(3) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(3)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(1) + fired.at(3) , tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(3)].a , tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(3)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(3)].a ,-fired.at(1) + fired.at(3) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].a - tofEvents[pl][fired.at(3)].a ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(3)].t ,-fired.at(1) + fired.at(3) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(1)].t - tofEvents[pl][fired.at(3)].t ,tofEvents[pl][fired.at(1)].y - tofEvents[pl][fired.at(3)].y );
					//}
					//if( fabs( tofEvents[pl][fired.at(2)].y - tofEvents[pl][fired.at(3)].y) < 10. && fabs(fired.at(2) - fired.at(3) ) < 15.){ // these are correlated hits to be clustered
						hMany_meanTime_diff[pl]			-> Fill( tofEvents[pl][fired.at(2)].t - tofEvents[pl][fired.at(3)].t );
						hMany_amp_diff[pl]			-> Fill( tofEvents[pl][fired.at(2)].a - tofEvents[pl][fired.at(3)].a );
						hMany_StrSpread_diff[pl]		-> Fill( -fired.at(2) + fired.at(3) );
						hMany_StrSpread_diff_MC[pl]		-> Fill( randomStripDiff() );
						hMany_YSpread_diff[pl]			-> Fill( tofEvents[pl][fired.at(2)].y - tofEvents[pl][fired.at(3)].y );
						hMany_YSpread_diff_MC[pl]		-> Fill( randomPosDiff() );
						hMany_2DSpread[pl]	-> Fill( -fired.at(2) + fired.at(3) , tofEvents[pl][fired.at(2)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_ampTime[pl]	-> Fill( tofEvents[pl][fired.at(2)].a - tofEvents[pl][fired.at(3)].a , tofEvents[pl][fired.at(2)].t - tofEvents[pl][fired.at(3)].t );
						hMany_2D_ampStrDiff[pl]	-> Fill( tofEvents[pl][fired.at(2)].a - tofEvents[pl][fired.at(3)].a ,-fired.at(2) + fired.at(3) );
						hMany_2D_ampYDiff[pl]	-> Fill( tofEvents[pl][fired.at(2)].a - tofEvents[pl][fired.at(3)].a ,tofEvents[pl][fired.at(2)].y - tofEvents[pl][fired.at(3)].y );
						hMany_2D_timeStrDiff[pl]-> Fill( tofEvents[pl][fired.at(2)].t - tofEvents[pl][fired.at(3)].t ,-fired.at(2) + fired.at(3) );
						hMany_2D_timeYDiff[pl]	-> Fill( tofEvents[pl][fired.at(2)].t - tofEvents[pl][fired.at(3)].t ,tofEvents[pl][fired.at(2)].y - tofEvents[pl][fired.at(3)].y );
					//}
					
				}

			}
						

		} // End of loop over events in file

	} // End of loop over files	
	
	
	outFile->cd();
	for( int pl = 0; pl < 20 ; pl++){
		hMany_meanTime_diff[pl]		-> Write();
		hMany_amp_diff[pl]		-> Write();
		hMany_StrSpread_diff[pl]	-> Write();
		hMany_StrSpread_diff_MC[pl]	-> Write();
		hMany_YSpread_diff[pl]		-> Write();
		hMany_YSpread_diff_MC[pl]	-> Write();
		
		hMany_2D_ampStrDiff[pl]->GetXaxis()->SetTitle("Amplitude Difference");
		hMany_2D_ampStrDiff[pl]->GetYaxis()->SetTitle("X Difference [strips apart]");
		hMany_2D_ampStrDiff[pl]		-> Write();	

		hMany_2D_ampYDiff[pl]->GetXaxis()->SetTitle("Amplitude Difference");
		hMany_2D_ampYDiff[pl]->GetYaxis()->SetTitle("Y Difference [cm]");
		hMany_2D_ampYDiff[pl]		-> Write();	

		hMany_2D_timeStrDiff[pl]->GetXaxis()->SetTitle("Time Difference [ns]");
		hMany_2D_timeStrDiff[pl]->GetYaxis()->SetTitle("X Difference [strips apart]");
		hMany_2D_timeStrDiff[pl]	-> Write();	

		hMany_2D_timeYDiff[pl]->GetXaxis()->SetTitle("Time Difference [ns]");
		hMany_2D_timeYDiff[pl]->GetYaxis()->SetTitle("Y Difference [cm]");
		hMany_2D_timeYDiff[pl]		-> Write();	

		hMany_2DSpread[pl]->GetXaxis()->SetTitle("X Difference [strips apart]");
		hMany_2DSpread[pl]->GetYaxis()->SetTitle("Y Difference [cm]");
		hMany_2DSpread[pl]		-> Write();

		hMany_2D_ampTime[pl]->GetXaxis()->SetTitle("Amplitude Difference");
		hMany_2D_ampTime[pl]->GetYaxis()->SetTitle("Time Difference [ns]");
		hMany_2D_ampTime[pl]		-> Write();

		for( int st = 0 ; st < 48 ; st++){
			hSingle_meanTime[pl][st]	-> Write();
			hMany_meanTime[pl][st] 		-> Write();
			hSingle_amp[pl][st]		-> Write();
			hMany_amp[pl][st]		-> Write();
		}
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


double fixWalk( double amp, double shift, double par0, double par1, double par2){
	
	if( par0 == -1 || par1 == -1 || par2 == -1 || shift == -1)
		return 0;
	
	return par0 + par1*exp( -(amp - shift) / par2 );
}

