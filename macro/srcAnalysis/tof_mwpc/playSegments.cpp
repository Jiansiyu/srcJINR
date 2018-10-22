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


#include "BmnTOF1Conteiner.h"
#include "BmnMwpcSegment.h"
#include "FairTrackParam.h"
#include "BmnTrack.h"

using namespace std;

int main( int argc, char ** argv){
	
	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tplaySegments /path/to/skim/file\n";
		return -1;
	}

	TH1D * hDAx1 = new TH1D("hDAx1","hDAx1",1000,-1,1);
	TH1D * hDAy1 = new TH1D("hDAy1","hDAy1",1000,-1,1);
	TH1D * hDAx2 = new TH1D("hDAx2","hDAx2",1000,-1,1);
	TH1D * hDAy2 = new TH1D("hDAy2","hDAy2",1000,-1,1);
	TH1D * hDAx1_fake = new TH1D("hDAx1_fake","hDAx1_fake",1000,-1,1);
	TH1D * hDAy1_fake = new TH1D("hDAy1_fake","hDAy1_fake",1000,-1,1);
	TH1D * hDAx2_fake = new TH1D("hDAx2_fake","hDAx2_fake",1000,-1,1);
	TH1D * hDAy2_fake = new TH1D("hDAy2_fake","hDAy2_fake",1000,-1,1);


	TFile * inFile = new TFile(argv[1]);
	TTree * inTree = (TTree*) inFile->Get("sk");
	TClonesArray * mwpc1Tracks = new TClonesArray("BmnTrack");
	TClonesArray * mwpc2Tracks = new TClonesArray("BmnTrack");
	TClonesArray * mwpc3Tracks = new TClonesArray("BmnTrack");
	TClonesArray * mwpc4Tracks = new TClonesArray("BmnTrack");
	inTree->SetBranchAddress("mwpc_tr1"	,&mwpc1Tracks);
	inTree->SetBranchAddress("mwpc_tr2"	,&mwpc2Tracks);
	inTree->SetBranchAddress("mwpc_tr3"	,&mwpc3Tracks);
	inTree->SetBranchAddress("mwpc_tr4"	,&mwpc4Tracks);

	const int nEv = inTree->GetEntries();

	for(int event = 0 ; event < nEv ; event++){
		
		if( (event % 10000) == 0) cout<< "Working on event " << event << "\n";

		mwpc1Tracks->Clear();
		mwpc2Tracks->Clear();
		mwpc3Tracks->Clear();
		mwpc4Tracks->Clear();

		inTree->GetEvent(event);

		double x1 = 0, y1 = 0, z1 = 0, tx1 = 0, ty1 = 0;
		// Only 1 segment so these should be 'good' matches
		if( mwpc1Tracks->GetEntriesFast() == 1 && mwpc2Tracks->GetEntriesFast() == 1){
			BmnMwpcSegment * mwpc1Sg = (BmnMwpcSegment*) mwpc1Tracks->At(0);
			FairTrackParam * par1 = (FairTrackParam *) mwpc1Sg->GetParamFirst();
			x1 = par1->GetX();
			y1 = par1->GetY();
			z1 = par1->GetZ();
			tx1 = par1->GetTx();
			ty1 = par1->GetTy();

			BmnMwpcSegment * mwpc2Sg = (BmnMwpcSegment*) mwpc2Tracks->At(0);
			FairTrackParam * par2 = (FairTrackParam *) mwpc2Sg->GetParamFirst();
			
			// Create a track between two segments using x1,x2,y1,y2
			// x(t) = x' * t + x0
			// y(t) = y' * t + y0
			// z(t) = t
			double Ax = ( par2->GetX() - par1->GetX() ) / fabs( par2->GetZ() - par1->GetZ() );
			double Ay = ( par2->GetY() - par1->GetY() ) / fabs( par2->GetZ() - par1->GetZ() );
			double delAx1 = par1->GetTx() - Ax;
			double delAx2 = par2->GetTx() - Ax;
			double delAy1 = par1->GetTy() - Ay;
			double delAy2 = par2->GetTy() - Ay;
			
			hDAx1-> Fill( delAx1 );
			hDAy1-> Fill( delAy1 );
			hDAx2-> Fill( delAx2 );
			hDAy2-> Fill( delAy2 );
		}

		// Now look for mismatch events to see how bkgrd should be
		int it = 1;
		if( event + it == nEv - 1) break;
		while( true ){
			mwpc2Tracks->Clear();
			inTree->GetEvent(event + it);
			if( mwpc2Tracks->GetEntries() == 1){
				BmnMwpcSegment * mwpc2Sg = (BmnMwpcSegment*) mwpc2Tracks->At(0);
				FairTrackParam * par2 = (FairTrackParam *) mwpc2Sg->GetParamFirst();
				double Ax = ( par2->GetX() - x1 ) / fabs( par2->GetZ() - z1 );
				double Ay = ( par2->GetY() - y1 ) / fabs( par2->GetZ() - z1 );
				double delAx1 = tx1 - Ax;
				double delAx2 = par2->GetTx() - Ax;
				double delAy1 = ty1 - Ay;
				double delAy2 = par2->GetTx() - Ay;

				hDAx1_fake-> Fill( delAx1 );
				hDAy1_fake-> Fill( delAy1 );
				hDAx2_fake-> Fill( delAx2 );
				hDAy2_fake-> Fill( delAy2 );
				break;
			}
			it++;
		}
		
		
		
	}

	TApplication theApp("App",&argc,argv);

	TCanvas * c1 = new TCanvas("c1");
	c1->Divide(2,2);
	c1->cd(1);
	hDAx1->SetLineColor(2);
	hDAx1->Draw();
	hDAx1_fake->Draw("same");
	c1->cd(2);
	hDAy1->SetLineColor(2);
	hDAy1->Draw();
	hDAy1_fake->Draw("same");
	c1->cd(3);
	hDAx2->SetLineColor(2);
	hDAx2->Draw();
	hDAx2_fake->Draw("same");
	c1->cd(4);
	hDAy2->SetLineColor(2);
	hDAy2->Draw();
	hDAy2_fake->Draw("same");
	c1->Update();


	theApp.Run();
	 return 0;
}
