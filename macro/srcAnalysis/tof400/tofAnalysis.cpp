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

	// Init BmnTOF1Detector
	BmnTOF1Detector * Plane[20];
	stripCenters[20][48];
	TString name;
	for (Int_t i = 0; i < 20; i++) {
		name = Form("Plane%d", i);
		cout << endl << "==============Init " << name << "=============" << endl;
		Plane[i] = new BmnTOF1Detector(i, 1, NULL);
		Plane[i]->SetGeoFile("geofile_full_src_noZshift.root");
	}//*/

	/*
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
	double walkFunc[20][48][4] = {0.};
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

	// Set up histograms:
	TString hName;
	TH1D ** hStripMult 	= new TH1D*[20];	// number of strips fired in plane
	TH1D ** hClusterMult	= new TH1D*[20];	// number of clusters ID'd in plane
	TH1D ** hStripCluster	= new TH1D*[20];	// number of strips per cluster in plane

	TH1D ** hxDiff		= new TH1D*[20];
	TH1D ** hyDiff		= new TH1D*[20];
	TH1D ** hxDiffMC	= new TH1D*[20];
	TH1D ** hyDiffMC	= new TH1D*[20];
	TH1D ** htDiff		= new TH1D*[20];
	TH2D ** hyxDiff		= new TH2D*[20];
	TH2D ** hxtDiff		= new TH2D*[20];
	TH2D ** hytDiff		= new TH2D*[20];

	TH1D *** hToF_SingleEvents = new TH1D**[20];
	TH1D *** hToF_ClusteredEvents = new TH1D**[20];

	for( int pl = 0 ; pl < 20 ; pl++){
		hToF_SingleEvents[pl] = new TH1D*[48];
		hToF_ClusteredEvents[pl] = new TH1D*[48];		

		hName = Form("hStripMult_%i",pl);
		hStripMult[pl]		= new TH1D(hName,hName,40,0,20);
		hName = Form("hClusterMult_%i",pl);
		hClusterMult[pl]	= new TH1D(hName,hName,40,0,20);
		hName = Form("hStripCluster_%i",pl);
		hStripCluster[pl]	= new TH1D(hName,hName,40,0,20);

		hName = Form("hyDiff_%i",pl);
		hyDiff[pl]		= new TH1D(hName,hName,600,-30,30);
		hName = Form("hyDiffMC_%i",pl);
		hyDiffMC[pl]		= new TH1D(hName,hName,600,-30,30);
		hName = Form("hxDiff_%i",pl);
		hxDiff[pl]		= new TH1D(hName,hName,192,-48,0);
		hName = Form("hxDiffMC_%i",pl);
		hxDiffMC[pl]		= new TH1D(hName,hName,192,-48,0);
		hName = Form("htDiff_%i",pl);
		htDiff[pl]		= new TH1D(hName,hName,1000,-50,50);

		hName = Form("hyxDiff_%i",pl);
		hyxDiff[pl]		= new TH2D(hName,hName,600,-30,30,192,-96,96);
		hName = Form("hytDiff_%i",pl);
		hytDiff[pl]		= new TH2D(hName,hName,600,-30,30,1000,-50,50);
		hName = Form("hxtDiff_%i",pl);
		hxtDiff[pl]		= new TH2D(hName,hName,192,-96,96,1000,-50,50);

		for( int st = 0 ; st < 48 ; st++){
			hName = Form("hToF_SingleEvents_%i_%i",pl,st);
			hToF_SingleEvents[pl][st] = new TH1D(hName,hName,4000,-15,85);
			hName = Form("hToF_ClusteredEvents_%i_%i",pl,st);
			hToF_ClusteredEvents[pl][st] = new TH1D(hName,hName,4000,-15,85);
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
		TClonesArray * bc3Data 	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * bc4Data	= new TClonesArray("BmnTrigWaveDigit");
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
			std::vector<int> stripsFired[20];
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

							if( sumAmps < 9 ) continue; // I don't want hits with low amplitude					

							// If there was already a created hit for this strip, only take the earliest one
							// so there will only be one hit per strip allowed
							if( (tofEvents[plane][strip].hit == false) || (tofEvents[plane][strip].hit == true && meanTime < tofEvents[plane][strip].t) ){
								double par0 = walkFunc[plane][strip][1];
								double par1 = walkFunc[plane][strip][2];
								double par2 = walkFunc[plane][strip][3];
								double shift = walkFunc[plane][strip][0];

								tofEvents[plane][strip].t = meanTime - t0Time + (-6.1 + 27./sqrt(t0Amp) ) - stripShift[plane][strip] - fixWalk(sumAmps, shift, par0, par1, par2);
								tofEvents[plane][strip].a = sumAmps;
								tofEvents[plane][strip].y = pos;
								tofEvents[plane][strip].hit == true;
								tofEvents[plane][strip].st = strip;

								if( std::find( stripsFired[plane].begin() , stripsFired[plane].end() , strip) == stripsFired[plane].end())
									stripsFired[plane].push_back( strip );

							}
						}
					}
				}
			}

			std::vector<TOFHit> tofCluster[20];
			std::vector<int> clusterList[20][48];
			// Now we need to cluster the ToF Events if the size in a plane is > 1
			for( int pl = 0 ; pl < 20 ; pl++){
				if( stripsFired[pl].size() == 0) continue;
				//cout << "Number of strips fired in plane: " << stripsFired[pl].size() << "\n";
				hStripMult[pl]->Fill( stripsFired[pl].size() );				

				std::sort( stripsFired[pl].begin() , stripsFired[pl].end() );
				//cout << "Need to consider clustering for # strips: " << stripsFired[pl].size() << "\n\tstrips: ";
				//for( int i = 0 ; i < stripsFired[pl].size() ; i++) cout << stripsFired[pl].at(i) << " ";
				//cout << "\n";

				// if there is only 1 strip fired, save that
				if( stripsFired[pl].size() == 1){
					TOFHit selection = tofEvents[pl][stripsFired[pl].at(0)];
					tofCluster[pl].push_back( selection );
					//cout << "\tNumber of clusters: " << tofCluster[pl].size() << "\n";
					hStripCluster[pl]->Fill( 1 );
				}
				else{
					for( int i = 0 ; i < stripsFired[pl].size() ; i++){
						// if there is more than 1 strip fired, cluster them
						clusterList[pl][stripsFired[pl].at(i) ].push_back( stripsFired[pl].at(i) );
						for( int j = 0 ; j < stripsFired[pl].size() ; j++){
							if( j <= i ) continue;
							int stOne = stripsFired[pl].at(i);
							int stTwo = stripsFired[pl].at(j);

							TOFHit one = tofEvents[pl][stOne];
							TOFHit two = tofEvents[pl][stTwo];

							// if delta T > 2ns between two strips, they are separate hits
							// if hits are 4 strips apart, they are separate hits
							// if hits are more than 10cm away in y, they are separate hits
							hxDiff[pl]->Fill( one.st - two.st);
							hxDiffMC[pl]->Fill( -randomStripDiff() );
							hyDiff[pl]->Fill( one.y - two.y);
							hyDiffMC[pl]->Fill( randomPosDiff() );
							htDiff[pl]->Fill( one.t - two.t);
							hyxDiff[pl]->Fill( one.y - two.y , one.st - two.st );
							hxtDiff[pl]->Fill( one.st - two.st , one.t - two.t );
							hytDiff[pl]->Fill( one.y - two.y , one.t - two.t );

							if( (fabs( one.t - two.t ) > 2) || (fabs( one.y - two.y ) > 5) || (fabs( one.st - two.st ) > 10) ){
							}
							else{
								clusterList[pl][one.st].push_back( two.st );
							}
						}
						//cout << "\tcluster list of strip: " << stripsFired[pl].at(i) << "\n\t\t";
						//for( int j = 0 ; j < clusterList[pl][stripsFired[pl].at(i) ].size() ; j++)
						//	cout << clusterList[pl][stripsFired[pl].at(i) ].at(j) << " ";
						//cout << "\n";
					}

					std::vector< std::vector<int> > result;
					// For all strips, take the cluster list
					for( int st = 0 ; st < 48 ; st++){
						//cout << "Working on plane, strip " << pl << " " << st << "\n";
						// For all the groups I already have, find out if any intersection with this cluster list, and if so, add it to group
						bool insert = true;
						//cout << "number of current clusters size: " << result.size() << "\n";
						for( int group = 0 ; group < result.size() ; group++){
							std::vector<int> tmp;
							std::vector<int> newGrp;
							std::set_intersection( clusterList[pl][st].begin() , clusterList[pl][st].end(), \
									result.at(group).begin() , result.at(group).end() , back_inserter(tmp) ) ;
							//cout << "intersection of how many elements?: " << tmp.size() << "\n";
							if( tmp.size() ){
								// Union two sets into the group and exit
								std::set_union( clusterList[pl][st].begin() , clusterList[pl][st].end(), \
										result.at(group).begin() , result.at(group).end() , back_inserter( newGrp ) );
								insert = false;
								result.at(group).swap(newGrp);
								break;
							}
						}
						// Create new group if no group wiht an intersection
						if( insert ){
							if( clusterList[pl][st].size() > 0)
								result.push_back(  clusterList[pl][st] );
						}

					}

					for( int group = 0 ; group < result.size() ; group++){
					}

					//cout << "\tnumber of clusters: " << result.size() << "\n";
					hClusterMult[pl]->Fill( result.size() );
					int TEST = 0;
					for( int group = 0 ; group < result.size() ; group++){
						//cout << "\t\tcluster " << group << ": ";
						hStripCluster[pl]->Fill( result.at(group).size() );
						//cout << "\tNumber of strips in cluster: " << result.at(group).size() << "\n";
						TOFHit selectStrip;
						//for( int hit = 0 ; hit < result.at(group).size() ; hit++){
						//cout << result.at(group).at(hit) << " ";
							if( tofEvents[pl][strip].a > selectStrip.a )
								selectStrip = tofEvents[pl][strip];
						}
							if( tofEvents[pl][strip].a > selectStrip.a )
								selectStrip = tofEvents[pl][strip];
						}
						tofCluster[pl].push_back( selectStrip );
						//cout << "\n";
					}
				}

				// Now we have all tof400 hits properly clustered for each plane so we can look at these hits.
				// We can have multiple hits for a given plane, and multiple planes firing, and we can look at all of these
				// hits and fill the ToF spectrum
				for( int hit = 0 ; hit < tofCluster[pl].size() ; hit++){
					int strip = tofCluster[pl].at(hit).st;
					if( tofCluster[pl].size() == 1){
						hToF_SingleEvents[pl][strip] -> Fill( tofCluster[pl].at(hit).t );
					}
					else{
						hToF_ClusteredEvents[pl][strip] -> Fill( tofCluster[pl].at(hit).t );
					}
				}
			}


		} // End of loop over events in file

	} // End of loop over files	


	TFile * outFile = new TFile("tofevents.root","RECREATE");
	outFile->cd();
	for( int pl = 0 ; pl < 20 ; pl++){
		hStripMult[pl]->Write();	
		hClusterMult[pl]->Write();	
		hStripCluster[pl]->Write();	

		hxDiff[pl]->Write();
		hyDiff[pl]->Write();
		hxDiffMC[pl]->Write();
		hyDiffMC[pl]->Write();
		htDiff[pl]->Write();
		hyxDiff[pl]->Write();
		hxtDiff[pl]->Write();
		hytDiff[pl]->Write();

		for( int st = 0 ; st < 48 ; st++){
			hToF_SingleEvents[pl][st]->Write();
			hToF_ClusteredEvents[pl][st]->Write();
		}
	}
	outFile->Write();
	outFile->Close();
	*/
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

