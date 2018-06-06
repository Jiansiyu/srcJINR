#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"

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
		cerr << "Successfully opened file " << argv[1] << " and saved it to address " << infile << "\n";
	}

	// Get the run number from input file:
	TString file = argv[1];
	TString run_number( file(file.Index(".")-9,4) );
	// Setup output file
	
	TFile * outFile = new TFile("tqdcPedestals"+run_number+".root","RECREATE");
	TH1D * BC1_ped = new TH1D("BC1_ped","BC1_ped",4500,-500,4000);
	TH1D * BC2_ped = new TH1D("BC2_ped","BC2_ped",4500,-500,4000);
	TH1D * BC3_ped = new TH1D("BC3_ped","BC3_ped",4500,-500,4000);
	TH1D * BC4_ped = new TH1D("BC4_ped","BC4_ped",4500,-500,4000);

	// Set up the tree
	TClonesArray * bc1Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc2Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc3Data 	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc4Data	= new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * t0Data	= new TClonesArray("BmnTrigDigit");
	TClonesArray * mcpData	= new TClonesArray("BmnTrigDigit");
	TClonesArray * t02Data	= new TClonesArray("BmnTrigDigit");	

	TTree * intree = NULL;
	intree = (TTree*) infile->Get("cbmsim");
	if (!intree)
	{
		cerr << "Could not find cbmsim tree. Perhaps the wrong type of input file. Bailing out.\n";
		return -3;
	}
	else
	{
		cerr << "Successfully loaded tree at address " << intree << "\n";
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
	int cnt = 0;
	double pedBC1, pedBC2, pedBC3, pedBC4;

	for (int event=0 ; event<nEvents ; event++)
	{
		if (event % 1000== 0)
			cerr << "Working on event " << event << "\n";

		intree->GetEvent(event);
		
		// Kill any event that doesn't have T0 entry, or that
		// has more than 1 TDC in T0 -- cuts out 7% of data	
		if( t0Data->GetEntriesFast() != 1) continue;
		cnt ++;

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

			// Now we use bc1Idx for actual analysis, and any other entries are
			// used for pedestal subtraction
			fillPedestal( BC1_ped , bc1Data, bc1Idx);
			fillPedestal( BC2_ped , bc2Data, bc2Idx);
			fillPedestal( BC3_ped , bc3Data, bc3Idx);
			fillPedestal( BC4_ped , bc4Data, bc4Idx);
			
		}
	}

	infile->Close();
	
	outFile->cd();
	BC1_ped->Write();
	BC2_ped->Write();
	BC3_ped->Write();
	BC4_ped->Write();
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

void fillPedestal( TH1D* hist, TClonesArray* data, int index){
	for( int m = 0 ; m < data->GetEntriesFast() ; m++){
		if( m == index) continue;
		BmnTrigWaveDigit * signal = (BmnTrigWaveDigit*) data->At(m);
		hist->Fill(	signal->GetPeak()	);
	}

}

