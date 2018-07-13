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
	TH1D * hCarbonIn_1400 	= new TH1D("hCarbonIn_1400","",3000,0,3000);
	TH1D * hWhateverIn_1400	= new TH1D("hWhateverIn_1400","",3000,0,3000);
	TH1D * hCarbonIn_1800 	= new TH1D("hCarbonIn_1800","",3000,0,3000);
	TH1D * hWhateverIn_1800	= new TH1D("hWhateverIn_1800","",3000,0,3000);
	TH1D * hCarbonIn_2200 	= new TH1D("hCarbonIn_2200","",3000,0,3000);
	TH1D * hWhateverIn_2200	= new TH1D("hWhateverIn_2200","",3000,0,3000);
	TH1D * hTOFPlanes	= new TH1D("hTOFPlanes","",20,0,20);
	TH1D * hTOFStrips	= new TH1D("hTOFStrips","",50,0,50);
	TH1D * hTOFSides	= new TH1D("hTOFSides","",2,0,2);
	TH2D * hTOFdt0_t0Amp	= new TH2D("hTOFdt0_t0Amp","",200,10,30,200,200,400);
	TH2D * hTOFdt1_t0Amp	= new TH2D("hTOFdt1_t0Amp","",200,10,30,200,200,400);
	TH2D * hTOFdt0_tofAmp	= new TH2D("hTOFdt0_tofAmp","",1500,5,35,1000,200,400);
	TH2D * hTOFdt1_tofAmp	= new TH2D("hTOFdt1_tofAmp","",1500,5,35,1000,200,400);
	TH1D * hTOFEntry0	= new TH1D("hTOFEntry0","",10,0,10);
	TH1D * hTOFEntry1	= new TH1D("hTOFEntry1","",10,0,10);
	TH2D * hTOFStrips0_Amp	= new TH2D("hTOFStrips0_Amp","",50,0,50,3000,5,35);
	TH2D * hTOFStrips1_Amp	= new TH2D("hTOFStrips1_Amp","",50,0,50,3000,5,35);
	TH1D ** hTOFdt0 	= new TH1D*[24];
	TH1D ** hTOFdt1 	= new TH1D*[24];
	TH1D ** hTOFdt0_veto 	= new TH1D*[24];
	TH1D ** hTOFdt1_veto 	= new TH1D*[24];
	TH1D * hTOFdt0_all 	= new TH1D("hTOFdt0_all","",8000,200,400);
	TH1D * hTOFdt1_all 	= new TH1D("hTOFdt1_all","",8000,200,400);
	TH1D * hTOFdt0_veto_all	= new TH1D("hTOFdt0_veto_all","",8000,200,400);
	TH1D * hTOFdt1_veto_all = new TH1D("hTOFdt1_veto_all","",8000,200,400);
	for( int st = 0; st < 24 ; st++){
		TString hName 	= Form("hTOFdt0_%i", st );
		hTOFdt0[st] 	= new TH1D(hName,hName,8000,200,400);
		hName 		= Form("hTOFdt1_%i", st );
		hTOFdt1[st] 	= new TH1D(hName,hName,8000,200,400);
		hName 		= Form("hTOFdt0_veto_%i", st );
		hTOFdt0_veto[st]= new TH1D(hName,hName,8000,200,400);
		hName 		= Form("hTOFdt1_veto_%i", st );
		hTOFdt1_veto[st]= new TH1D(hName,hName,8000,200,400);
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

		TClonesArray * x1lData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * x1rData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * x2lData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * x2rData	= new TClonesArray("BmnTrigWaveDigit");

		TClonesArray * y1lData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * y1rData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * y2lData	= new TClonesArray("BmnTrigWaveDigit");
		TClonesArray * y2rData	= new TClonesArray("BmnTrigWaveDigit");


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
		// Arm branches
		intree->SetBranchAddress("TQDC_X1_Left"	,&x1lData);
		intree->SetBranchAddress("TQDC_X1_Right",&x1rData);
		intree->SetBranchAddress("TQDC_X2_Left"	,&x2lData);
		intree->SetBranchAddress("TQDC_X2_Right",&x2rData);
		intree->SetBranchAddress("TQDC_Y1_Left"	,&y1lData);
		intree->SetBranchAddress("TQDC_Y1_Right",&y1rData);
		intree->SetBranchAddress("TQDC_Y2_Left"	,&y2lData);
		intree->SetBranchAddress("TQDC_Y2_Right",&y2rData);

		// Loop over events
		for (int event=0 ; event<nEvents ; event++){
			double adcBC1 = 0;
			double adcBC2 = 0;
			t0Data->Clear();
			bc1Data->Clear();	
			bc2Data->Clear();	
			bc3Data->Clear();	
			bc4Data->Clear();	
			tofData->Clear();
			x1lData->Clear();
			x1rData->Clear();
			x2lData->Clear();
			x2rData->Clear();
			y1lData->Clear();
			y1rData->Clear();
			y2lData->Clear();
			y2rData->Clear();

			intree->GetEvent(event);
			
			bool veto = 	(x1lData->GetEntries() == 0) & (x1rData->GetEntries() == 0) & (x2lData->GetEntries() ==0) & (x2rData->GetEntries() == 0) & \
					(y1lData->GetEntries() == 0) & (y1rData->GetEntries() == 0) & (y2lData->GetEntries() ==0) & (y2rData->GetEntries() == 0);


			// Kill any event that doesn't have T0 entry, or that
			// has more than 1 TDC in T0 -- cuts out 7% of data	
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
			
		
			if( adcBC1 != 0 && adcBC2 != 0){
				if( fabs( (field_voltage) - 87 ) < 1)
					hWhateverIn_1400->Fill( sqrt( adcBC1 * adcBC2 ) );
				else if( fabs( (field_voltage) - 107.7) < 1)
					hWhateverIn_1800->Fill( sqrt( adcBC1 * adcBC2 ) );
				else if( fabs( (field_voltage) - 123.7) < 1)
					hWhateverIn_2200->Fill( sqrt( adcBC1 * adcBC2 ) );
			}

			if( fabs( adcBC1 - (BC1_CPeak-pedBC1) ) > 2*BC1_CWidth )
				continue;
			if( fabs( adcBC2 - (BC2_CPeak-pedBC2) ) > 2*BC2_CWidth )
				continue;
			
			if( fabs( (field_voltage) - 87 ) < 1)
				hCarbonIn_1400->Fill( sqrt( adcBC1 * adcBC2 ) );
			else if( fabs( (field_voltage) - 107.7) < 1)
				hCarbonIn_1800->Fill( sqrt( adcBC1 * adcBC2 ) );
			else if( fabs( (field_voltage) - 123.7) < 1)
				hCarbonIn_2200->Fill( sqrt( adcBC1 * adcBC2 ) );

			int cnt0 = 0;
			int cnt1 = 0;
	
			for( int en = 0 ; en < tofData->GetEntriesFast() ; en++ ){
				BmnTof1Digit * signal = (BmnTof1Digit*) tofData->At(en);
				hTOFPlanes->Fill( signal->GetPlane() );
				
				if( signal->GetPlane() == 6){
					hTOFStrips->Fill( signal->GetStrip() );
					
					if( signal->GetStrip() == 23){
						hTOFSides->Fill( signal->GetSide() );

						if( signal->GetSide()==0 && tofData->GetEntries() ==2){
							cnt0++;
							if( cnt0 == 1){
								hTOFdt0_t0Amp->Fill( t0Amp , 	signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
								hTOFdt0_tofAmp->Fill( signal->GetAmplitude() , 	signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
							}
						}				

						if( signal->GetSide()==1 && tofData->GetEntries() == 2){
							cnt1++;
							if( cnt1 == 1){
								hTOFdt1_t0Amp->Fill( t0Amp , 	signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
								hTOFdt1_tofAmp->Fill( signal->GetAmplitude() , 	signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
							}
						}
						
					}
				}



				if( signal->GetPlane() == 6 && signal->GetSide() == 0 && signal->GetStrip() > 10 && signal->GetStrip() < 35 && tofData->GetEntries() == 2){
					hTOFStrips0_Amp->Fill( signal->GetStrip() , signal->GetAmplitude() );
					if( signal->GetAmplitude() > 6 && signal->GetAmplitude() < 15){
						hTOFdt0[signal->GetStrip() - 11] -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						hTOFdt0_all -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						if( veto ){
							hTOFdt0_veto[signal->GetStrip() - 11] -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
							hTOFdt0_veto_all -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						}
					}
				}
				if( signal->GetPlane() == 6 && signal->GetSide() == 1 && signal->GetStrip() > 10 && signal->GetStrip() < 35 && tofData->GetEntries() == 2){
					hTOFStrips1_Amp->Fill( signal->GetStrip() , signal->GetAmplitude() );
					if( signal->GetAmplitude() > 6 && signal->GetAmplitude() < 15){
						hTOFdt1[signal->GetStrip() - 11] -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						hTOFdt1_all -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						if( veto ){
							hTOFdt1_veto[signal->GetStrip() - 11] -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
							hTOFdt1_veto_all -> Fill( signal->GetTime() - t0Time + (-6.1 + 27./sqrt(t0Amp) ) );
						}
					}
				}
				
			

			}	
			if( cnt0 != 0)
				hTOFEntry0->Fill( cnt0 );
			if( cnt1 != 0)
				hTOFEntry1->Fill( cnt1 );
			
		}
	}	
	
	//TApplication theApp("App",&argc,argv);
	TCanvas * c1 = new TCanvas("c1");
	c1->Divide(3,1);
	c1->cd(1);
	hWhateverIn_1400->SetLineColor(2);
	hWhateverIn_1400->SetTitle("BC1-BC2 on 1400A");
	hWhateverIn_1400->Draw();
	hCarbonIn_1400->Draw("same");
	c1->cd(2);
	hWhateverIn_1800->SetLineColor(2);
	hWhateverIn_1800->SetTitle("BC1-BC2 on 1800A");
	hWhateverIn_1800->Draw();
	hCarbonIn_1800->Draw("same");
	c1->cd(3);
	hWhateverIn_2200->SetLineColor(2);
	hWhateverIn_2200->SetTitle("BC1-BC2 on 2200A");
	hWhateverIn_2200->Draw();
	hCarbonIn_2200->Draw("same");
	c1->Update();
	
	TCanvas * c2 = new TCanvas("c2");
	hTOFPlanes->Draw();
	c2->Update();

	TCanvas * c3 = new TCanvas("c3");
	hTOFStrips->Draw();
	c3->Update();

	TCanvas * c4 = new TCanvas("c4");
	hTOFSides->Draw();
	c4->Update();

	TCanvas * c5 = new TCanvas("c5");
	c5->Divide(2,1);
	c5->cd(1);
	hTOFdt0_t0Amp->Draw("colz");
	c5->cd(2);
	hTOFdt1_t0Amp->Draw("colz");
	c5->Update();


	TCanvas * c6 = new TCanvas("c6");
	c6->Divide(2,1);
	c6->cd(1);
	hTOFdt0_tofAmp->Draw("colz");
	c6->cd(2);
	hTOFdt1_tofAmp->Draw("colz");
	c6->Update();

	TCanvas * c7 = new TCanvas("c7");
	hTOFEntry0->SetLineColor(2);
	hTOFEntry0->Draw();
	hTOFEntry1->Draw("same");
	c7->Update();

	TCanvas * c8 = new TCanvas("c8");
	c8->Divide(2,1);
	c8->cd(1);
	hTOFStrips0_Amp->Draw("colz");
	c8->cd(2);
	hTOFStrips1_Amp->Draw("colz");
	c8->Update();

	TCanvas * c9 = new TCanvas("c9");
	c9->Divide(2,1);
	c9->cd(1);
	for( int st = 0; st < 24 ; st++){
		if( st == 0)
			hTOFdt0[st]->Draw();
		else{
			hTOFdt0[st]->SetLineColor(st+1);
			hTOFdt0[st]->Draw("same");
		}
	}
	c9->cd(2);
	for( int st = 0; st < 24 ; st++){
		if( st == 0)
			hTOFdt1[st]->Draw();
		else{
			hTOFdt1[st]->SetLineColor(st+1);
			hTOFdt1[st]->Draw("same");
		}
	}
	c9->Update();

	TCanvas * c10 = new TCanvas("c10");
	c10->Divide(2,1);
	c10->cd(1);
	hTOFdt0_all->Draw();
	c10->cd(2);
	hTOFdt1_all->Draw();
	c10->Update();
	

	//theApp.Run();
	
	TFile * outFile = new TFile("test.root","RECREATE");
	outFile->cd();
	hCarbonIn_1400  ->Write();	
	hWhateverIn_1400->Write();
	hCarbonIn_1800  ->Write();	
	hWhateverIn_1800->Write();
	hCarbonIn_2200  ->Write();	
	hWhateverIn_2200->Write();
	hTOFPlanes      ->Write(); 
	hTOFStrips      ->Write(); 
	hTOFSides       ->Write(); 
	hTOFdt0_t0Amp   ->Write(); 
	hTOFdt1_t0Amp   ->Write(); 
	hTOFdt0_tofAmp  ->Write(); 
	hTOFdt1_tofAmp  ->Write(); 
	hTOFEntry0      ->Write(); 
	hTOFEntry1      ->Write(); 
	hTOFStrips0_Amp ->Write();	
	hTOFStrips1_Amp ->Write();	
	hTOFdt0_all     ->Write(); 
	hTOFdt1_all     ->Write(); 
	hTOFdt0_veto_all->Write(); 
	hTOFdt1_veto_all->Write(); 
	for( int st = 0; st < 24 ; st++){
		hTOFdt0[st] -> Write();
		hTOFdt1[st] -> Write();
		hTOFdt0_veto[st] -> Write();
		hTOFdt1_veto[st] -> Write();
	}
	hTOFdt0_all->Write();
	hTOFdt1_all->Write();
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

