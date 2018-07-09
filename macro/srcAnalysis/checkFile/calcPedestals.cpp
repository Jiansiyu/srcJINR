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
#include "TVectorT.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;

void findIdx( TClonesArray* data, int &index , double refT);
void fillPedestal( TH1D* pedHist, TH1D* sigHist, TH1D* tdiffPed, TH1D* tdiffSig, TClonesArray* data, int index, double startT);

int main(int argc, char ** argv)
{
	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tcalcPedestals /path/to/digi/file\n";
		return -1;
	}


	// Set up the input file
	TFile * infile = NULL;
	infile = new TFile(argv[1]);
	if (!infile)
	{
		cerr << "Could not open file " << argv[1] <<"\n"
			<< "\tBailing out\n";
		return -2;
	}
	else
	{
		//cerr << "Successfully opened file " << argv[1] << " and saved it to address " << infile << "\n";
	}

	// Get the run number from input file:
	TString file = argv[1];
	TString run_number( file(file.Index(".")-9,4) );
	// Setup output file
	
	TFile * outFile = new TFile("qualityCheck/checked_"+run_number+".root","RECREATE");
	TH1D * BC1_ped = new TH1D("BC1_ped","BC1_ped",4500,-500,4000);
	TH1D * BC2_ped = new TH1D("BC2_ped","BC2_ped",4500,-500,4000);
	TH1D * BC3_ped = new TH1D("BC3_ped","BC3_ped",4500,-500,4000);
	TH1D * BC4_ped = new TH1D("BC4_ped","BC4_ped",4500,-500,4000);
	
	TH1D * BC1_sig = new TH1D("BC1_sig","BC1_sig",4500,-500,4000);
	TH1D * BC2_sig = new TH1D("BC2_sig","BC2_sig",4500,-500,4000);
	TH1D * BC3_sig = new TH1D("BC3_sig","BC3_sig",4500,-500,4000);
	TH1D * BC4_sig = new TH1D("BC4_sig","BC4_sig",4500,-500,4000);

	TH1D * BC1_entry = new TH1D("BC1_entry","BC1_entry",15,0,15);
	TH1D * BC2_entry = new TH1D("BC2_entry","BC2_entry",15,0,15);
	TH1D * BC3_entry = new TH1D("BC3_entry","BC3_entry",15,0,15);
	TH1D * BC4_entry = new TH1D("BC4_entry","BC4_entry",15,0,15);

	TH1D * BC1_tdiff_sig = new TH1D("BC1_tdiff_sig","BC1_tdiff_sig",2000,-50,50);
	TH1D * BC2_tdiff_sig = new TH1D("BC2_tdiff_sig","BC2_tdiff_sig",2000,-50,50);
	TH1D * BC3_tdiff_sig = new TH1D("BC3_tdiff_sig","BC3_tdiff_sig",2000,-50,50);
	TH1D * BC4_tdiff_sig = new TH1D("BC4_tdiff_sig","BC4_tdiff_sig",2000,-50,50);

	TH1D * BC1_tdiff_ped = new TH1D("BC1_tdiff_ped","BC1_tdiff_ped",2000,-50,50);
	TH1D * BC2_tdiff_ped = new TH1D("BC2_tdiff_ped","BC2_tdiff_ped",2000,-50,50);
	TH1D * BC3_tdiff_ped = new TH1D("BC3_tdiff_ped","BC3_tdiff_ped",2000,-50,50);
	TH1D * BC4_tdiff_ped = new TH1D("BC4_tdiff_ped","BC4_tdiff_ped",2000,-50,50);


	// Set up the tree
	TClonesArray * bc1Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc2Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc3Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc4Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * t0Data	= new TClonesArray("BmnTrigDigit");

	TTree * intree = NULL;
	intree = (TTree*) infile->Get("cbmsim");
	if (!intree)
	{
		cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";
		return -3;
	}
	else
	{
		//cerr << "Successfully loaded tree at address " << intree << "\n";
	}
	
	const int nEvents = intree->GetEntries();
	// TQDC Branches
	intree->SetBranchAddress("TQDC_BC1"	,&bc1Data);
	intree->SetBranchAddress("TQDC_T0"	,&bc2Data);
	intree->SetBranchAddress("TQDC_BC3"	,&bc3Data);
	intree->SetBranchAddress("TQDC_BC4"	,&bc4Data);
	
	// TDC Branches
	intree->SetBranchAddress("T0"		,&t0Data);

	// Loop over events
	double cnt = 0;
	double cntBC1 = 0, cntBC2 = 0, cntBC3 = 0, cntBC4 = 0;
	double cntBC1_1 = 0, cntBC2_1 = 0, cntBC3_1 = 0, cntBC4_1 = 0;
	double pedBC1, pedBC2, pedBC3, pedBC4;

	for (int event=0 ; event<nEvents ; event++)
	{
		t0Data->Clear();
		bc1Data->Clear();	
		bc2Data->Clear();	
		bc3Data->Clear();	
		bc4Data->Clear();	

		//if (event % 1000== 0)
		//	cerr << "Working on event " << event << "\n";

		intree->GetEvent(event);
		
		// Kill any event that doesn't have T0 entry, or that
		// has more than 1 TDC in T0 -- cuts out 7% of data	
		if( t0Data->GetEntriesFast() != 1) continue;
		cnt ++;

		BC1_entry->Fill( bc1Data->GetEntriesFast() );
		BC2_entry->Fill( bc2Data->GetEntriesFast() );
		BC3_entry->Fill( bc3Data->GetEntriesFast() );
		BC4_entry->Fill( bc4Data->GetEntriesFast() );
	
		BmnTrigDigit * t0Signal = (BmnTrigDigit*) t0Data->At(0);
		double t0Time = t0Signal->GetTime();

		if( bc1Data->GetEntriesFast() == 1) cntBC1_1++;
		if( bc1Data->GetEntriesFast() ){
			cntBC1++;
			int bc1Idx;
			findIdx(bc1Data,bc1Idx,t0Time);
			fillPedestal( BC1_ped , BC1_sig , BC1_tdiff_ped , BC1_tdiff_sig , bc1Data , bc1Idx , t0Time);
		}

		if( bc2Data->GetEntriesFast() == 1) cntBC2_1++;
		if( bc2Data->GetEntriesFast() ){
			cntBC2++;
			int bc2Idx;
			findIdx(bc2Data,bc2Idx,t0Time);
			fillPedestal( BC2_ped , BC2_sig , BC2_tdiff_ped , BC2_tdiff_sig , bc2Data , bc2Idx , t0Time);
		}
		
		if( bc3Data->GetEntriesFast() == 1) cntBC3_1++;
		if( bc3Data->GetEntriesFast() ){
			cntBC3++;
			int bc3Idx;
			findIdx(bc3Data,bc3Idx,t0Time);
			fillPedestal( BC3_ped , BC3_sig , BC3_tdiff_ped , BC3_tdiff_sig , bc3Data , bc3Idx , t0Time);
		}
		
		if( bc4Data->GetEntriesFast() == 1) cntBC4_1++;
		if( bc4Data->GetEntriesFast() ){
			cntBC4++;
			int bc4Idx;
			findIdx(bc4Data,bc4Idx,t0Time);
			fillPedestal( BC4_ped , BC4_sig , BC4_tdiff_ped , BC4_tdiff_sig , bc4Data , bc4Idx , t0Time);

		}

	}


	
	pedBC1 = BC1_ped->GetXaxis()->GetBinCenter(BC1_ped->GetMaximumBin());
	pedBC2 = BC2_ped->GetXaxis()->GetBinCenter(BC2_ped->GetMaximumBin());
	pedBC3 = BC3_ped->GetXaxis()->GetBinCenter(BC3_ped->GetMaximumBin());
	pedBC4 = BC4_ped->GetXaxis()->GetBinCenter(BC4_ped->GetMaximumBin());

	TFitResultPtr init_bc1 = BC1_sig->Fit("gaus","QES","",0,4000);
	TFitResultPtr fin_bc1  = BC1_sig->Fit("gaus","QES","",0, init_bc1->Parameter(1) + init_bc1->Parameter(2) );
	TFitResultPtr init_bc2 = BC2_sig->Fit("gaus","QES","",0,4000);
	TFitResultPtr fin_bc2  = BC2_sig->Fit("gaus","QES","",0, init_bc2->Parameter(1) + init_bc2->Parameter(2) );
	
	TFitResultPtr tBC1 = BC1_tdiff_sig->Fit("gaus","QES");
	TFitResultPtr tBC2 = BC2_tdiff_sig->Fit("gaus","QES");
	TFitResultPtr tBC3 = BC3_tdiff_sig->Fit("gaus","QES");
	TFitResultPtr tBC4 = BC4_tdiff_sig->Fit("gaus","QES");

	cout << "\n********************************************************************************************\n";
	cout << "T0 Rate:\t" << std::setprecision(4) << cnt/intree->GetEntries() << "\n\n";
	cout << "\tPedestal\tTime Sig Peak\tTime Sig Res\tTime Ped Int\t# TDC = 1\n";
	cout <<  "BC1:\t" << pedBC1 << "\t\t" << tBC1->Parameter(1) << "\t\t" << tBC1->Parameter(2) << "\t\t" <<BC1_tdiff_ped->Integral() << "\t\t" << cntBC1_1/cntBC1 << "\n";
	cout <<  "BC2:\t" << pedBC2 << "\t\t" << tBC2->Parameter(1) << "\t\t" << tBC2->Parameter(2) << "\t\t" <<BC2_tdiff_ped->Integral() << "\t\t" << cntBC2_1/cntBC2 << "\n";
	cout <<  "BC3:\t" << pedBC3 << "\t\t" << tBC3->Parameter(1) << "\t\t" << tBC3->Parameter(2) << "\t\t" <<BC3_tdiff_ped->Integral() << "\t\t" << cntBC3_1/cntBC3 << "\n";
	cout <<  "BC4:\t" << pedBC4 << "\t\t" << tBC4->Parameter(1) << "\t\t" << tBC4->Parameter(2) << "\t\t" <<BC4_tdiff_ped->Integral() << "\t\t" << cntBC4_1/cntBC4 << "\n\n";
	cout << "Carbon-In Peak (BC1-BC2):\t" << fin_bc1->Parameter(1) << " , " << fin_bc2->Parameter(1) << "\n\twidth of:\t\t" << fin_bc1->Parameter(2) << " , " << fin_bc2->Parameter(2) << "\n\n";
	cout << "********************************************************************************************\n\n";
	
	
	TVectorT<double> ped(4);
	ped[0] = pedBC1;
	ped[1] = pedBC2;	
	ped[2] = pedBC3;	
	ped[3] = pedBC4;	

	TVectorT<double> tdcCnt(4);
	tdcCnt[0] = cntBC1_1/cntBC1;
	tdcCnt[1] = cntBC2_1/cntBC2;
	tdcCnt[2] = cntBC3_1/cntBC3;
	tdcCnt[3] = cntBC4_1/cntBC4;
	
	TVectorT<double> t0Rate(1);
	t0Rate[0] = cnt/intree->GetEntries();

	TVectorT<double> tDiff(4);
	tDiff[0] = tBC1->Parameter(1);
	tDiff[1] = tBC2->Parameter(1);
	tDiff[2] = tBC3->Parameter(1);
	tDiff[3] = tBC4->Parameter(1);

	TVectorT<double> tRes(4);
	tRes[0] = tBC1->Parameter(2);
	tRes[1] = tBC2->Parameter(2);
	tRes[2] = tBC3->Parameter(2);
	tRes[3] = tBC4->Parameter(2);

	TVectorT<double> pedTime(4);
	pedTime[0] = BC1_tdiff_ped->Integral();
	pedTime[1] = BC2_tdiff_ped->Integral();
	pedTime[2] = BC3_tdiff_ped->Integral();
	pedTime[3] = BC4_tdiff_ped->Integral();

	TVectorT<double> carbonIn(2);
	carbonIn[0] = fin_bc1->Parameter(1);
	carbonIn[1] = fin_bc2->Parameter(1);

	TVectorT<double> carbonInWidth(2);
	carbonInWidth[0] = fin_bc1->Parameter(2);
	carbonInWidth[1] = fin_bc2->Parameter(2);
	
	infile->Close();
	

	outFile->cd();
	ped.Write("ped");
	tdcCnt.Write("tdcCnt");
	t0Rate.Write("t0Rate");
	tDiff.Write("tDiff");
	tRes.Write("tRes");
	pedTime.Write("pedTime");
	carbonIn.Write("carbonIn");
	carbonInWidth.Write("carbonInWidth");

	BC1_ped->Write();
	BC2_ped->Write();
	BC3_ped->Write();
	BC4_ped->Write();

	BC1_sig->Write();
	BC2_sig->Write();
	BC3_sig->Write();
	BC4_sig->Write();

	BC1_entry->Write();
	BC2_entry->Write();
	BC3_entry->Write();
	BC4_entry->Write();

	BC1_tdiff_sig->Write();
	BC2_tdiff_sig->Write();
	BC3_tdiff_sig->Write();
	BC4_tdiff_sig->Write();

	BC1_tdiff_ped->Write();
	BC2_tdiff_ped->Write();
	BC3_tdiff_ped->Write();
	BC4_tdiff_ped->Write();

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

void fillPedestal( TH1D* pedHist, TH1D* sigHist, TH1D* tdiffPed, TH1D* tdiffSig, TClonesArray* data, int index, double startT){
	for( int m = 0 ; m < data->GetEntriesFast() ; m++){
		if( m == index){
			BmnTrigWaveDigit * signal = (BmnTrigWaveDigit*) data->At(m);
			sigHist->Fill ( signal->GetPeak() );
			tdiffSig->Fill( signal->GetTime() - startT );
		}
		else{
			BmnTrigWaveDigit * signal = (BmnTrigWaveDigit*) data->At(m);
			pedHist->Fill ( signal->GetPeak() );
			tdiffPed->Fill( signal->GetTime() - startT );
		}
	}

}

