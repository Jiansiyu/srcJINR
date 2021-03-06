// @(#)bmnroot/mwpc:$Id$
// Author: Vasilisa Lenivenko <vasilisa@jinr.ru> 2018-07-18

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcHitFinder                                                       //
//                                                                            //
// Implementation of an algorithm developed by                                //
// Vasilisa Lenivenko  and Vladimir Palchik                                   //
// to the BmnRoot software                                                    //
//                                                                            //
// The algorithm serves for searching for hits                                //
// in the Multi Wire Prop. Chambers of the BM@N experiment                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <Rtypes.h>
#include <climits>
#include <vector>
#include "TCanvas.h"

#include "BmnMwpcHitFinder.h"
#include "BmnTrack.h"
#include "BmnMwpcTrack.h"
#include "BmnMwpcSegment.h"
#include <BmnEventHeader.h>
#include <algorithm>


static Float_t workTime = 0.0;

using namespace std;
using namespace TMath;

struct segments {
	Int_t    Nhits     = 0;
	Double_t Chi2      = 0.;
	Double_t coord[6]  = {-999., -999., -999., -999., -999., -999.};
	Double_t clust[6]  = {0.,    0.,    0.,    0.,    0.,    0.};
	Double_t param[4]  = { 999.,  999., 999.,  999.};
};

bool compareSegments(const segments &a, const segments &b){
	return a.Chi2 < b.Chi2;
}

BmnMwpcHitFinder::BmnMwpcHitFinder(Bool_t isExp, Int_t runPeriod, Int_t runNumber) :
	fEventNo(0),
	expData(isExp){
		fRunPeriod        = runPeriod;
		fRunNumber        = runNumber;
		fInputBranchName  = "MWPC";
		fBmnEventHeaderBranchName = "EventHeader";
		fOutputBranchName = "BmnMwpcSegment";
		nInputDigits      = 3;
		nTimeSamples      = 3;
		kBig              = 100;
		kCh_min           = 0;
		if(fRunPeriod == 6 || (fRunPeriod == 7 && fRunNumber > 3588) ){
			kNumPairs       = 1;
			kCh_max         = 2;
		}else if(fRunPeriod == 7 && fRunNumber <= 3588){
			kNumPairs       = 2;
			kCh_max         = 4;
		}
		fBmnEvQualityBranchName = "BmnEventQuality";

	}
BmnMwpcHitFinder::~BmnMwpcHitFinder() {
}

InitStatus BmnMwpcHitFinder::Init() {
	if (!expData) return kERROR;
	if (fDebug) cout << " BmnMwpcHitFinder::Init() " << endl;



	//FairRootManager* ioman = FairRootManager::Instance();
	//fBmnMwpcDigitArray = (TClonesArray*) ioman->GetObject(fInputBranchName);
	//if (!fBmnMwpcDigitArray){
	//	cout<<"BmnMwpcHitFinder::Init(): branch "<<fInputBranchName<<" not found! Task will be deactivated"<<endl;
	//	SetActive(kFALSE);
	//	return kERROR;
	//}

	//fBmnMwpcEventHeader = (TClonesArray*) ioman->GetObject(fBmnEventHeaderBranchName);

	//fBmnMwpcSegmentsArray = new TClonesArray(fOutputBranchName);
	//ioman->Register(fOutputBranchName.Data(), "MWPC", fBmnMwpcSegmentsArray, kTRUE);

	fMwpcGeometrySRC = new BmnMwpcGeometrySRC(fRunPeriod, fRunNumber); 
	kNChambers = fMwpcGeometrySRC -> GetNChambers();
	kNPlanes   = fMwpcGeometrySRC -> GetNPlanes();
	kNWires    = fMwpcGeometrySRC -> GetNWires();
	if (fDebug) printf("C-P-W: %d %d %d\n", kNChambers, kNPlanes, kNWires);
	// cout<<" MWPC runPeriod "<<fRunPeriod<<" fRunNumber "<<fRunNumber<<endl;

	ChCent= new TVector3[kNChambers];
	ChZ =   new Float_t[kNChambers];
	Zmid = new Float_t[kNChambers];

	for (int i=0; i < kNChambers; ++i){ 
		TH1D *h;
		h = new TH1D(Form("Np_best_Ch%d", i), Form("Np_best_Ch%d", i), 6, 1.0, 7.0); fList.Add(h); hNp_best_Ch.push_back(h);
		h = new TH1D(Form("Nbest_Ch%d",   i), Form("Nbest_Ch%d",   i), 6, 0.0, 6.0); fList.Add(h); hNbest_Ch.push_back(h);

		h = new TH1D(Form("occupancyXp%d", i), Form("occupancyXp%d",i), 250, -10.,10.);fList.Add(h); hoccupancyXp.push_back(h);
		h = new TH1D(Form("occupancyUp%d", i), Form("occupancyUp%d",i), 250, -10.,10.);fList.Add(h); hoccupancyUp.push_back(h);
		h = new TH1D(Form("occupancyVp%d", i), Form("occupancyVp%d",i), 250, -10.,10.);fList.Add(h); hoccupancyVp.push_back(h);
		h = new TH1D(Form("occupancyXm%d", i), Form("occupancyXm%d",i), 250, -10.,10.);fList.Add(h); hoccupancyXm.push_back(h);
		h = new TH1D(Form("occupancyUm%d", i), Form("occupancyUm%d",i), 250, -10.,10.);fList.Add(h); hoccupancyUm.push_back(h);
		h = new TH1D(Form("occupancyVm%d", i), Form("occupancyVm%d",i), 250, -10.,10.);fList.Add(h); hoccupancyVm.push_back(h);

		h = new TH1D(Form("WireClust_Ch%d", i), Form("WireClust_Ch%d",i), 100, 0., 100.);    fList.Add(h); hWireClust.push_back(h);
		h = new TH1D(Form("ClusterSize_Ch%d", i), Form("ClusterSize_Ch%d", i), 15, 0.,15.);  fList.Add(h); hClusterSize.push_back(h);

		ChCent[i] = fMwpcGeometrySRC->GetChamberCenter(i);
		ChZ[i]    = ChCent[i].Z();
	}

	for (int i=0; i < 24; ++i){
		TH1D *h;
		h = new TH1D(Form("Time%d", i), Form("Time%d", i), 500, 0., 500.);            fList.Add(h); hTime.push_back(h);
		h = new TH1D(Form("Occupancy%d", i), Form("Occupancy%d", i), 100, 0., 100.);    fList.Add(h); hOccupancy.push_back(h);
	}

	// Segment finding cuts
	kMinHits  = 5;//5;// for alignment kMinHits  = 6;
	kChi2_Max = 20.;
	kmaxSeg = 5;//10;
	kChMaxAllWires = 200;
	// Constants
	dw        = fMwpcGeometrySRC -> GetWireStep(); //0.25; // [cm] // wires step
	dw_half   = 0.5 * dw;
	sq3       = sqrt(3.);
	sq12      = sqrt(12.);
	sigma     = dw / sq12;
	kMiddlePl = 47.25; // Center of wires plane
	// Matrices
	matrA  = new Double_t*[4];
	matrb  = new Double_t*[4];
	// Arrays
	kPln             = new Int_t*[kNChambers];
	iw_Ch            = new Int_t*[kNChambers];
	wire_Ch          = new Int_t**[kNChambers];
	xuv_Ch           = new Float_t**[kNChambers];
	Wires_Ch         = new Int_t**[kNChambers];    
	clust_Ch         = new Int_t**[kNChambers];   
	XVU_Ch           = new Float_t**[kNChambers];  
	Nhits_Ch         = new Int_t*[kNChambers];     
	Nseg_Ch          = new Int_t[kNChambers];      
	Nbest_Ch         = new Int_t[kNChambers];      
	ind_best_Ch      = new Int_t*[kNChambers];
	best_Ch_gl       = new Int_t*[kNChambers];
	Chi2_ndf_Ch      = new Double_t*[kNChambers];
	Chi2_ndf_best_Ch = new Double_t*[kNChambers];
	par_ab_Ch        = new Double_t**[kNChambers]; 
	XVU              = new Float_t*[kNChambers];   
	XVU_cl           = new Float_t*[kNChambers];  
	kZ_loc           = new Float_t*[kNChambers];
	z_gl             = new Float_t*[kNChambers];
	sigm2            = new Float_t[kNPlanes]; 
	ipl              = new Int_t[kNPlanes];
	z2               = new Float_t[kNPlanes];
	DigitsArray      = new Double_t**[kNChambers];
	ClusterSize      = new Int_t**[kNChambers];
	Nclust           = new Int_t*[kNChambers];
	Coord_wire       = new Double_t**[kNChambers];
	Coord_xuv        = new Double_t**[kNChambers];
	XVU_coord        = new Double_t**[kNChambers];
	Cluster_coord    = new Double_t**[kNChambers];
	Nhits_seg        = new Int_t*[kNChambers];
	Chi2_ndf_seg     = new Double_t*[kNChambers];
	Coor_seg         = new Double_t**[kNChambers];
	Cluster_seg      = new Double_t**[kNChambers];
	par_ab_seg       = new Double_t**[kNChambers];
	Nbest_seg        = new Int_t[kNChambers]; 


	for(Int_t i = 0; i < kNChambers; ++i){

		if (i== 0 || i== 2) { Zmid[i] = (ChZ[i]     - ChZ[i + 1]) *  0.5;}
		if (i== 1 || i== 3) { Zmid[i] = (ChZ[i - 1] - ChZ[i])     * -0.5;}
		if (fDebug) printf("Chamber %d Z: %f Zmid: %f\n", i, ChZ[i], Zmid[i]);

		kPln[i]             = new Int_t[kNPlanes];
		iw_Ch[i]            = new Int_t[kNPlanes]; 
		kZ_loc[i]           = new Float_t[kNPlanes];
		z_gl[i]             = new Float_t[kNPlanes];
		Nhits_Ch[i]         = new Int_t[kBig];
		wire_Ch[i]          = new Int_t*[kNWires];
		xuv_Ch[i]           = new Float_t*[kNWires];
		Wires_Ch[i]         = new Int_t*[kNPlanes];
		clust_Ch[i]         = new Int_t*[kNPlanes];
		XVU_Ch[i]           = new Float_t*[kNPlanes];
		par_ab_Ch[i]        = new Double_t*[4];
		XVU[i]              = new Float_t[kNPlanes];
		XVU_cl[i]           = new Float_t[kNPlanes];
		ind_best_Ch[i]      = new Int_t[kmaxSeg];
		best_Ch_gl[i]       = new Int_t[kmaxSeg];
		Chi2_ndf_Ch[i]      = new Double_t[kBig];
		Chi2_ndf_best_Ch[i] = new Double_t[kmaxSeg];
		DigitsArray[i]      = new Double_t*[kNPlanes];
		ClusterSize[i]      = new Int_t*[kNPlanes];
		Nclust[i]           = new Int_t[kNPlanes];
		Coord_wire[i]       = new Double_t*[kNPlanes];
		Coord_xuv[i]        = new Double_t*[kNPlanes];
		XVU_coord[i]        = new Double_t*[kNPlanes];
		Cluster_coord[i]    = new Double_t*[kNPlanes];
		Nhits_seg[i]        = new Int_t[kBig];
		Chi2_ndf_seg[i]     = new Double_t[kBig];
		Coor_seg[i]         = new Double_t*[kNPlanes];
		Cluster_seg[i]      = new Double_t*[kNPlanes];
		par_ab_seg[i]       = new Double_t*[4];

		for(int iWire = 0; iWire < kNWires; iWire++){
			wire_Ch[i][iWire] = new Int_t[kNPlanes];
			xuv_Ch[i][iWire]  = new Float_t[kNPlanes];
		}
		for(int iPla = 0; iPla < kNPlanes; ++iPla){
			Wires_Ch[i][iPla]     = new Int_t[kBig];
			clust_Ch[i][iPla]     = new Int_t[kBig];
			XVU_Ch[i][iPla]       = new Float_t[kBig];
			DigitsArray[i][iPla]  = new Double_t[kNWires];
			ClusterSize[i][iPla]  = new Int_t[kBig];
			Coord_wire[i][iPla]   = new Double_t[kBig];
			Coord_xuv[i][iPla]    = new Double_t[kBig];
			XVU_coord[i][iPla]    = new Double_t[kBig];
			Cluster_coord[i][iPla]= new Double_t[kBig];
			Coor_seg[i][iPla]     = new Double_t[kBig];
			Cluster_seg[i][iPla] = new Double_t[kBig];


			if (fRunPeriod == 6 || (fRunPeriod == 7 && fRunNumber > 3588) ){

				kPln[i][0] = 4;
				kPln[i][1] = 5;
				kPln[i][2] = 0;
				kPln[i][3] = 1;
				kPln[i][4] = 2;
				kPln[i][5] = 3;//{4,5,0,1,2,3,  7,11,6,10,9,8,  0,0,0,0,0,0,  0,0,0,0,0,0};                                                                                    
				kZ_loc[i][iPla] = -0.5 + iPla;
				if(iPla == 4) { kZ_loc[i][iPla] = -2.5;}
				if(iPla == 5) { kZ_loc[i][iPla] = -1.5;}
			}
			if (fRunPeriod == 7 && fRunNumber <= 3588 ){//SRC

				if ( i == 0 ){
					kPln[i][0] = 5;  kZ_loc[i][0] = -1.5;
					kPln[i][1] = 0;  kZ_loc[i][1] = -0.5;
					kPln[i][2] = 1;  kZ_loc[i][2] =  0.5;
					kPln[i][3] = 2;  kZ_loc[i][3] =  1.5;
					kPln[i][4] = 3;  kZ_loc[i][4] =  2.5;
					kPln[i][5] = 4;  kZ_loc[i][5] = -2.5;
				}
				else if(i == 1){
					kPln[i][0] = 1;  kZ_loc[1][0] = -1.5;
					kPln[i][1] = 0;  kZ_loc[i][1] = -2.5;
					kPln[i][2] = 5;  kZ_loc[i][2] =  2.5;
					kPln[i][3] = 4;  kZ_loc[i][3] =  1.5;
					kPln[i][4] = 3;  kZ_loc[i][4] =  0.5;
					kPln[i][5] = 2;  kZ_loc[i][5] = -0.5;
				}
				else if ( i == 2 || i == 3){
					kPln[i][0] = 4;
					kPln[i][1] = 5;
					kPln[i][2] = 0;
					kPln[i][3] = 1;
					kPln[i][4] = 2;
					kPln[i][5] = 3;                                                                                                                                                                       
					kZ_loc[i][iPla] = -0.5 + iPla;
					if(iPla == 4) { kZ_loc[i][iPla] = -2.5;}
					if(iPla == 5) { kZ_loc[i][iPla] = -1.5;}

				}// i 2  3
			}//7 run

		}//iPla

		for(int ii = 0; ii < 4; ++ii){ // 4 parameters: tan(x), tan(y), x ,y 
			par_ab_Ch[i][ii] = new Double_t[kBig];
			par_ab_seg[i][ii] = new Double_t[kBig];
			matrA[ii]        = new Double_t[4];
			matrb[ii]        = new Double_t[4];         
		}//4
	}//kChamber

	if (fRunPeriod == 6){
		//  kPln[1][0] = 1;  kZ_loc[1][0] =-2.5;//    kPln[3][0] = 1;//run6-II
		//    kPln[1][3] = 4;  kZ_loc[1][3] = 0.5;//    kPln[3][3] = 4;//
	}


	if (fDebug) printf("Chamber  Plane  kZ_loc   z_gl\n");
	for(Int_t i = 0; i < kNChambers; ++i){
		for(int ii = 0; ii < kNPlanes; ii++){
			z_gl[i][ii] =  Zmid[i] + kZ_loc[i][ii];
			if (fDebug) printf("%5d  %5d %8.4f  %8.4f\n", i, ii, kZ_loc[i][ii], z_gl[i][ii]);
		}
	}

	//fBmnEvQuality = (TClonesArray*) ioman->GetObject(fBmnEvQualityBranchName);
	cout << "\tSUCCESS!\n";
	return kSUCCESS;
}


void BmnMwpcHitFinder::Exec(Option_t* opt, TClonesArray* fBmnMwpcDigitArray, TClonesArray* fBmnMwpcSegmentsArray, UInt_t evHead) {
	if (!IsActive()) return;
	//if ( fBmnEvQuality) {
	//	BmnEventQuality* evQual = (BmnEventQuality*) fBmnEvQuality->UncheckedAt(0);
	//	if (!evQual->GetIsGoodEvent())
	//		return;
	//}

	clock_t tStart = clock();
	PrepareArraysToProcessEvent();


	//if (fDebug) cout << "\n======================== MWPC hit finder exec started =====================\n" << endl;
	//printf("Event number: %d\n", fEventNo++);

	//BmnEventHeader* evHeader = (BmnEventHeader*) fBmnMwpcEventHeader->UncheckedAt(0);
	//if (!evHeader) return;

	Short_t st, wire, pn, pl;
	UInt_t ts;
	//UInt_t evHead = evHeader->GetEventId();

	for (Int_t iDigit = 0; iDigit < fBmnMwpcDigitArray -> GetEntries(); iDigit++) {
		BmnMwpcDigit* digit = (BmnMwpcDigit*) fBmnMwpcDigitArray ->At (iDigit);

		st = digit   -> GetStation();
		wire = digit -> GetWireNumber();
		pl = digit   -> GetPlane();
		ts = digit   -> GetTime(); //ns

		Int_t ind = st*6 + pl;
		hTime.at(ind) -> Fill(ts);
		hOccupancy.at(ind) -> Fill(wire);
		if ( ts < 80 || ts > 280 ) continue;
		pn = kPln[st][pl];// made for the canonical sequence / x- v- u+ x+ v+ u-/

		//  Loop over repeated wires -- TODO: this isn't working properly because it only looks for repeated 0 wire??
		Bool_t repeat_wire = 0; 
		if (iw_Ch[st][pn] > 0) {
			for (Int_t ix = 0; ix < iw_Ch[st][pn]; ix++) {
				if (wire == wire_Ch[st][ ix ][pn]  ) {
					repeat_wire = 1;
					break;
				}
			}//ix
		}
		if (repeat_wire){ 
			continue;
		}

		if (iw_Ch[st][pn] >= 80){
			continue;
		}

		DigitsArray[st][pn][wire] = 1.;		// basically just 0,1 for a wire that has been hit

		iw_Ch[st][pn]++; 			// counts how many wires are in a given station in a given plane
							// and keeps the correct order of the planes in actual geometry

	}// iDigit


	for (Int_t iChamber = 0; iChamber < kNChambers; iChamber++) {

		Int_t counter = 0; Int_t counter_pl[6];

		for (Int_t iplane = 0; iplane < kNPlanes; iplane++) {
			counter_pl[iplane] = 0;
			counter += iw_Ch[iChamber][iplane];
			counter_pl[iplane] += iw_Ch[iChamber][iplane];   
			hWireClust.at(iChamber)->Fill(counter_pl[iplane]);   
		}//iplane

		if (counter < kMinHits || counter > kChMaxAllWires ) continue;

		Clustering(iChamber, ClusterSize,  DigitsArray, Coord_wire,  Coord_xuv, Nclust);

		for(Int_t iCase= 1; iCase < 9; iCase ++){
			SegmentFinder(iChamber, Nclust, Coord_xuv,  ClusterSize, Nseg_Ch, XVU_coord, Cluster_coord, Nhits_Ch, kMinHits, iCase, kBig);
		}
		if (fDebug) cout<<"  SegmentFinder: Nseg_["<<iChamber<<"]= "<<Nseg_Ch[iChamber]<<endl;

		ProcessSegments(iChamber, Nseg_Ch, XVU_coord, Cluster_coord, Nhits_Ch, kZ_loc, kMinHits, sigma, kChi2_Max, 
				Nhits_seg ,Chi2_ndf_seg, Coor_seg, Cluster_seg, par_ab_seg, Nbest_seg);
		if (fDebug)  cout<<"ProcessSegments: Nbest["<<iChamber<<"] "<<Nbest_seg[iChamber]<<endl;

		hNbest_Ch.at(iChamber) -> Fill(Nbest_Ch[iChamber]);
	}//iChamber

	vector<Double_t>vtmpCoor;
	vector<Double_t>vtmpClust;


	for (Int_t iChamber = 0; iChamber < kNChambers; iChamber++) { // For each chamber of MWPC
		for (Int_t ise = 0; ise < Nbest_seg[iChamber]; ise++) {
			if (Nhits_seg[iChamber][ise]  > 3) { // If I had more than 3 hits per segment,

				BmnMwpcSegment *pSeg = new ((*fBmnMwpcSegmentsArray)[fBmnMwpcSegmentsArray->GetEntriesFast()]) BmnMwpcSegment();
				//	BmnMwpcTrack *pSeg = new ((*fBmnMwpcSegmentsArray)[fBmnMwpcSegmentsArray->GetEntriesFast()]) BmnMwpcTrack();
				pSeg->SetChi2(Chi2_ndf_seg[iChamber][ise]);
				pSeg->SetNHits(Nhits_seg[iChamber][ise]);
				pSeg->SetFlag(ise);
				pSeg->SetLength(evHead);

				vtmpCoor.clear();
				vtmpClust.clear();

				for(Int_t i1 = 0 ; i1 < 6; i1++){
					vtmpCoor.push_back(Coor_seg[iChamber][i1][ise]);
					vtmpClust.push_back(Cluster_seg[iChamber][i1][ise]);
				}
				pSeg -> SetClust(vtmpClust);
				pSeg -> SetCoord(vtmpCoor);

				FairTrackParam pSegParams;
				pSegParams.SetPosition(TVector3(par_ab_seg[iChamber][1][ise], par_ab_seg[iChamber][3][ise],ChZ[iChamber])); 
				pSegParams.SetTx(par_ab_seg[iChamber][0][ise]);
				pSegParams.SetTy(par_ab_seg[iChamber][2][ise]);
				pSeg->SetParamFirst(pSegParams);
			}//if
		}//ise
	}//[iChamber]


	clock_t tFinish = clock();
	workTime += ((Float_t) (tFinish - tStart)) / CLOCKS_PER_SEC;
	if (fDebug)  cout << "\n======================== MWPC hit finder exec finished ====================" << endl;
}

void BmnMwpcHitFinder::Clustering(Int_t chNum, Int_t*** ClusterSize_, Double_t*** DigitsArray_, Double_t*** Coord_wire_, Double_t*** Coord_xuv_, Int_t **Nclust_){

	Int_t Nfirst[kCh_max][kNPlanes], Nlast[kCh_max][kNPlanes];

	for (Int_t ipll = 0; ipll < kNPlanes; ipll++) {
		Nfirst[chNum][ipll] = -1;
		Nlast[chNum][ipll] = -1;

		for (Int_t iwire = 0; iwire < 95; iwire++){
			ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]] = 0;
			// if (DigitsArray_[chNum][ipll][iwire] > 0) cout<<"  DigitsArray_["<<chNum<<"]["<<ipll<<"]["<<iwire<<"] "<<DigitsArray_[chNum][ipll][iwire]<<endl;

			if( Nfirst[chNum][ipll] < 0 &&  DigitsArray_[chNum][ipll][iwire] == 0.)   continue;
			if( Nfirst[chNum][ipll] < 0 &&  DigitsArray_[chNum][ipll][iwire] > 0. )   Nfirst[chNum][ipll] = iwire;  
			if( Nfirst[chNum][ipll] >=0 &&  DigitsArray_[chNum][ipll][iwire+1] == 0.) Nlast[chNum][ipll] = iwire;

			if (Nfirst[chNum][ipll] >= 0 && Nlast[chNum][ipll] > 0){
				ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]] = Nlast[chNum][ipll] - Nfirst[chNum][ipll] + 1;

				//   cout<<" Nfirst "<<Nfirst[chNum][ipll]<<" Nlast "<<Nlast[chNum][ipll]<<" ClusterSize_["<<chNum<<"]["<<ipll<<"] "<<ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]]<<endl;
				hClusterSize.at(chNum)->Fill(ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]]);

				if ( ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]] > 15.) continue;

				Coord_wire_[chNum][ipll][Nclust_[chNum][ipll]] = Nfirst[chNum][ipll]+ 0.5*ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]] - 0.5;
				Coord_xuv_[chNum][ipll][Nclust_[chNum][ipll]] = (Coord_wire_[chNum][ipll][Nclust_[chNum][ipll]]- kMiddlePl)* dw;

				// made for the canonical sequence / x- v- u+ x+ v+ u-/ 0 1 5 have back reading
				if (ipll == 0 || ipll == 1 || ipll == 5 ) Coord_xuv_[chNum][ipll][Nclust_[chNum][ipll]] = - Coord_xuv_[chNum][ipll][Nclust_[chNum][ipll]];

				hoccupancyXm.at(chNum)->Fill(Coord_xuv_[chNum][0][Nclust_[chNum][ipll]]);
				hoccupancyVm.at(chNum)->Fill(Coord_xuv_[chNum][1][Nclust_[chNum][ipll]]);
				hoccupancyUp.at(chNum)->Fill(Coord_xuv_[chNum][2][Nclust_[chNum][ipll]]);
				hoccupancyXp.at(chNum)->Fill(Coord_xuv_[chNum][3][Nclust_[chNum][ipll]]);
				hoccupancyVp.at(chNum)->Fill(Coord_xuv_[chNum][4][Nclust_[chNum][ipll]]);
				hoccupancyUm.at(chNum)->Fill(Coord_xuv_[chNum][5][Nclust_[chNum][ipll]]);

				//cout<<" Coord_xuv_["<< chNum<<"]["<<ipll<<"]["<<Nclust_[chNum][ipll]<<"] "<<Coord_xuv_[chNum][ipll][Nclust_[chNum][ipll]]<<endl;
				//  cout<<"  Nfirst "<< Nfirst[chNum][ipll]<<" Nlast "<<Nlast[chNum][ipll]<<" ClusterSize_ "<<ClusterSize_[chNum][ipll][Nclust_[chNum][ipll]]<<" Coord_wire_ "<<Coord_wire_[chNum][ipll][Nclust_[chNum][ipll]]<<endl;
				//  cout<<" Coord_xuv_ "<<Coord_xuv_[chNum][ipll][Nclust_[chNum][ipll]]<<endl;

				Nclust_[chNum][ipll]++;

				Nfirst[chNum][ipll] = -1;
				Nlast[chNum][ipll] = -1;

			}//if (Nfirst[chNum][ipll] >= 0		
		}//iwire 
	}// ipll

}//Clustering


void BmnMwpcHitFinder::SegmentFinder(Int_t chNum, Int_t **Nclust_, Double_t ***Coord_xuv_,  Int_t ***ClusterSize_,
		Int_t *Nsegm, Double_t ***XVU_coor, Double_t ***Cluster_coor, Int_t **Nhits_Ch_,Int_t minHits, Short_t code, Int_t kBig_ ) {

	//Coord_xuv_      - coordinates of all clusters 
	Int_t minHits4_5   = minHits;
	Double_t min_for_triples   = 5.;  // minimum delta wires for tree planes 
	Double_t min_for_conjugate = 3.; // minimum delta wires for conjugate plane 

	// code : first triples is 
	// 1  {X-, V-, U+}
	// 2  {X-, V+, U+}
	// 3  {X-, V-, U-}
	// 7  {X+, V-, U+}
	// 5  {X+, V+, U+}
	// 6  {X+, V-, U-}
	// 4  {X-, V+, U-}
	// 8  {X+, V+, U-}

	Int_t x = 0, v = 1, u = 2 , x1 = 3, v1 = 4, u1 = 5;//MK

		// For given 'code' this tells me which plane is 
		// x,v,u or x1,v1,u1. but i have 4 stations with 6 planes, so
		// why do i have 8 code values??

	switch (code) {
		case 1: x = 0; v = 1; u = 2; x1 = 3; v1 = 4; u1 = 5; break;
		case 2: x = 0; v = 4; u = 2; x1 = 3; v1 = 1; u1 = 5; break;
		case 3: x = 0; v = 1; u = 5; x1 = 3; v1 = 4; u1 = 2; break;
		case 7: x = 3; v = 1; u = 2; x1 = 0; v1 = 4; u1 = 5; break;
		case 5: x = 3; v = 4; u = 2; x1 = 0; v1 = 1; u1 = 5; break;
		case 6: x = 3; v = 1; u = 5; x1 = 0; v1 = 4; u1 = 2; break;
		case 4: x = 0; v = 4; u = 5; x1 = 3; v1 = 1; u1 = 2; break;
		case 8: x = 3; v = 4; u = 5; x1 = 0; v1 = 1; u1 = 2; break;
	}


	if (Nsegm[chNum] > kBig_ - 2) return;// MP


	for (Int_t ix = 0; ix < Nclust_[chNum][x]; ix++) { 

		for (Int_t iv = 0; iv < Nclust_[chNum][v]; iv++) {

			for (Int_t iu = 0; iu < Nclust_[chNum][u]; iu++) {

				if (Nsegm[chNum] > kBig_ - 2) return;

				// --for repeated triples--
				Bool_t it_was = 0;
				if (Nsegm[chNum] > 0) {
					for (Int_t iseg = 0; iseg < Nsegm[chNum]; iseg++) {
						Bool_t it_was_x = 0;
						Bool_t it_was_v = 0;
						Bool_t it_was_u = 0;

						if (XVU_coor[chNum][x][iseg] == Coord_xuv_[chNum][x][ix] ) it_was_x = 1; 	   
						if (XVU_coor[chNum][v][iseg] == Coord_xuv_[chNum][v][iv]) it_was_v = 1; 	   
						if (XVU_coor[chNum][u][iseg] == Coord_xuv_[chNum][u][iu]) it_was_u = 1; 

						it_was = it_was_x * it_was_v * it_was_u;

						if (it_was) {break; }
					}  // iseg
				}  //  Nsegm[chNum] > 0

				if (it_was) continue;
				// --for repeated triples.

				//--- main equation---// u + v - x = delta

				if (fabs(Coord_xuv_[chNum][u][iu] + Coord_xuv_[chNum][v][iv] - Coord_xuv_[chNum][x][ix]) < min_for_triples*dw ) {

					//  3p-candidate new Nsegm
					XVU_coor[chNum][x][Nsegm[chNum]]     = Coord_xuv_[chNum][x][ix];  
					XVU_coor[chNum][v][Nsegm[chNum]]     = Coord_xuv_[chNum][v][iv]; 
					XVU_coor[chNum][u][Nsegm[chNum]]     = Coord_xuv_[chNum][u][iu]; 
					Cluster_coor[chNum][x][Nsegm[chNum]] = ClusterSize_[chNum][x][ix];   
					Cluster_coor[chNum][v][Nsegm[chNum]] = ClusterSize_[chNum][v][iv];	
					Cluster_coor[chNum][u][Nsegm[chNum]] = ClusterSize_[chNum][u][iu];	

					XVU_coor[chNum][x1][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][v1][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][u1][Nsegm[chNum]]    = -999.;
					Cluster_coor[chNum][x1][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][v1][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][u1][Nsegm[chNum]]= -1;

					Nhits_Ch_[chNum][Nsegm[chNum]] = 3;
					// cout<<" Nhits 3 = "<<Nhits_Ch_[chNum][Nsegm[chNum]]<<endl;

					//-- if 3p-candidate was look for conjugate coord  
					if ( XVU_coor[chNum][x][Nsegm[chNum]] != -999.){
						for (Int_t ix2 = 0; ix2 < Nclust_[chNum][x1]; ix2++) { 

							if(abs(  XVU_coor[chNum][x][Nsegm[chNum]]  -  Coord_xuv_[chNum][x1][ix2] ) <  min_for_conjugate *dw ) {
								XVU_coor[chNum][x1][Nsegm[chNum]]     = Coord_xuv_[chNum][x1][ix2]; 
								Cluster_coor[chNum][x1][Nsegm[chNum]] = ClusterSize_[chNum][x1][ix2]; 
								Nhits_Ch_[chNum][Nsegm[chNum]]++;// 4 points   
							}//abs( 
						}//ix2
					} //if it was 

					//  cout<<" Nhits 4 = "<<Nhits_Ch_[chNum][Nsegm[chNum]]<<endl;

					if ( XVU_coor[chNum][v][Nsegm[chNum]] != -999.){
						for (Int_t iv2 = 0; iv2 < Nclust_[chNum][v1]; iv2++) { 

							if(abs( XVU_coor[chNum][v][Nsegm[chNum]]  -  Coord_xuv_[chNum][v1][iv2] ) <  min_for_conjugate *dw ) {
								XVU_coor[chNum][v1][Nsegm[chNum]]     = Coord_xuv_[chNum][v1][iv2]; 
								Cluster_coor[chNum][v1][Nsegm[chNum]] = ClusterSize_[chNum][v1][iv2];
								Nhits_Ch_[chNum][Nsegm[chNum]]++;// 5 points   
							}//abs( 
						}//iv2
					} //if it was 
					//  cout<<" Nhits 5 = "<<Nhits_Ch_[chNum][Nsegm[chNum]]<<endl;

					if ( XVU_coor[chNum][u][Nsegm[chNum]] != -999.){
						for (Int_t iu2 = 0; iu2 < Nclust_[chNum][u1]; iu2++) { 

							if(abs( XVU_coor[chNum][u][Nsegm[chNum]]  -  Coord_xuv_[chNum][u1][iu2] ) <  min_for_conjugate *dw ) {
								XVU_coor[chNum][u1][Nsegm[chNum]]     = Coord_xuv_[chNum][u1][iu2]; 
								Cluster_coor[chNum][u1][Nsegm[chNum]] = ClusterSize_[chNum][u1][iu2];
								Nhits_Ch_[chNum][Nsegm[chNum]]++;// 6 points   
							}//abs( 

						}//iu2
					} //if it was 

					// cout<<" Nhits 6 = "<<Nhits_Ch_[chNum][Nsegm[chNum]]<<endl;

					if (Nsegm[chNum] > 15) minHits4_5=5;
					if (Nsegm[chNum] > 30) minHits4_5=6;
					//  cout<<endl;
					// cout<<" Nsegm "<< Nsegm[chNum] <<" Nhits_Ch_["<<chNum<<"]["<<Nsegm[chNum]<<"] "<<Nhits_Ch_[chNum][Nsegm[chNum]]<<"  minHits4_5  "<< minHits4_5 <<endl;

					if (Nhits_Ch_[chNum][Nsegm[chNum]] >= minHits4_5) {

						Nsegm[chNum]++; 
					}           
					if (Nhits_Ch_[chNum][Nsegm[chNum]] < minHits4_5) {  // }else{
						//--segment deleting 
						Nhits_Ch_[chNum][Nsegm[chNum]] = 0;

					XVU_coor[chNum][x ][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][v ][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][u ][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][x1][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][v1][Nsegm[chNum]]    = -999.;
					XVU_coor[chNum][u1][Nsegm[chNum]]    = -999.;
					Cluster_coor[chNum][x ][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][v ][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][u ][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][x1][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][v1][Nsegm[chNum]]= -1;
					Cluster_coor[chNum][u1][Nsegm[chNum]]= -1;
				}//else
				}// u + v - x = delta

				//*/
				if (Nsegm[chNum] > kBig_ - 2) break;
			}//iu
			if (Nsegm[chNum] > kBig_ - 2) break;
		}//iv
		if (Nsegm[chNum] > kBig_ - 2) break;
	}//ix
	if (Nsegm[chNum] > kBig_ - 2)return;
	// cout<<" Nsegm["<<chNum<<"] "<<Nsegm[chNum]<<endl;

}//SegmentFinder1


void BmnMwpcHitFinder::ProcessSegments( Int_t chNum,  Int_t *Nsegm, Double_t ***XVU_coor, Double_t ***Cluster_coor, Int_t **Nhits_Ch_, Float_t **z_loc, Int_t Min_hits, Double_t sigma_, Double_t kChi2_Max_,
		Int_t **Nhits_seg_ , Double_t **Chi2_ndf_seg_, Double_t ***coor_seg, Double_t ***cluster_seg, Double_t ***Par_ab_seg, Int_t *Nbest_seg_) {

	segments seg[kNChambers][kBig];

	Double_t par_ab[kNChambers][kNPlanes][kBig];
	Double_t dx_[kNPlanes];
	Double_t Chi2[kNChambers][kBig]; 
	Double_t Chi2_ndf[kNChambers][kBig]; 
	Double_t sigm[kNPlanes], sigm2_[kNPlanes];
	Int_t h1[kNPlanes];

	Double_t delta =  3*dw;
	Double_t Chi2_Super_min = 0.1;
	Int_t Min_hits6 = Min_hits;

	// TEMPORARY OUT SEGMENTS ARRY
	int      OutSegCount = 0;
	segments OutSegArray[kmaxSeg];

	if (Nsegm[chNum] > kBig - 2) return;
	for (Int_t Nhitm = kNPlanes; Nhitm > Min_hits - 1; Nhitm--) {// -- first view 6 points 
		Bool_t ifNhitm = 0;

		if (Nhitm < Min_hits6) break;

		for (Int_t iseg = 0; iseg < Nsegm[chNum]; iseg++) {  //---cycle by segments
			if (Nbest_seg_[chNum] > kBig - 2) return;
			ifNhitm = 1;

			if (Nhits_Ch_[chNum][iseg] != Nhitm) continue;
			// cout<<" iseg "<<iseg<<" Nhits "<<Nhits_Ch_[chNum][iseg]<<endl; 

			for(Int_t i = 0; i < 6; i++){
				sigm[i]= 0.; 
				h1[i] = 0;
				if ( XVU_coor[chNum][i][iseg]  > -999.) {
					h1[i] = 1;
					sigm[i] = ( Cluster_coor[chNum][i][iseg]*dw)/sq12;
					sigm2_[i] = sigm[i]*sigm[i];
				}//if coord was
				if (fDebug) cout<<" chNum "<<chNum<<" i "<<i<<" iseg "<<iseg<<" z_loc "<<z_loc[chNum][i]<<" coor "<<XVU_coor[chNum][i][iseg]<<endl;
				//	cout<<" h "<<h1[i]<<" Cluster_coor["<<chNum<<"]["<<i<<"]["<<iseg<<"]= "<<Cluster_coor[chNum][i][iseg]<<" sigm "<<sigm[i]<<endl;
			}
			if (fDebug) cout<<endl;


			Amatr = new Double_t*[4]; bmatr = new Double_t*[4];
			for(Int_t ii=0; ii<4; ii++){
				Amatr[ii] = new Double_t[4]; bmatr[ii] = new Double_t[4];
			}
			for(Int_t im=0; im<4; im++){
				for(Int_t ii=0; ii<4; ii++){  
					Amatr[im][ii] = 0.;	bmatr[im][ii] = 0.;
				}
			}
			Double_t matrF[4] = {0,0,0,0};//free coef 

			FillFitMatrix(chNum, Amatr, z_loc, sigm2_, h1); 
			FillFreeCoefVector(chNum, matrF, XVU_coor, iseg, z_loc , sigm2_, h1);

			//Gaussian algorithm for 4x4 matrix inversion 
			Double_t A0matr[4][4];	 			
			for (Int_t i1 = 0; i1 < 4; i1++){
				for (Int_t j1 = 0; j1 < 4; j1++){
					A0matr[i1][j1] = Amatr[i1][j1];
				}
			}

			InverseMatrix(Amatr,bmatr);			 			  

			Double_t sum;		  
			Double_t A1[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};  
			//	  cout<<" A1 "<<endl;

			for (Int_t i1 = 0; i1 < 4; ++i1) 
				for (Int_t j1 = 0; j1 < 4; ++j1) {
					sum = 0; 
					for (Int_t k1 = 0; k1 < 4; ++k1) {
						Double_t a0 = A0matr[i1][k1];
						Double_t b0 = bmatr[k1][j1];        
						sum += a0 * b0;       
						A1[i1][j1] = sum;
					}
				}


			for(Int_t i1 = 0 ; i1 < 4; i1++){
				par_ab[chNum][i1][iseg] = 0;
				for(Int_t j1 = 0; j1 < 4; j1++){
					par_ab[chNum][i1][iseg] += bmatr[i1][j1]*matrF[j1];
					//if (fDebug) 	cout<<" i1 "<<i1<<" bmatr "<<bmatr[i1][j1]<<" F "<<matrF[j1] <<endl;
				}    
			} // i1

			//---Chi2 calculating---
			Chi2[chNum][iseg] =  0.;
			Chi2_ndf[chNum][iseg] =  0.;


			for(Int_t i1 = 0 ; i1 < 6; i1++){
				dx_[i1] = 0.;

				if ( XVU_coor[chNum][i1][iseg]  > -999.){
					if(i1==0 || i1==3) dx_[i1] = XVU_coor[chNum][i1][iseg] - par_ab[chNum][0][iseg]*z_loc[chNum][i1]-par_ab[chNum][1][iseg];
					if(i1==2 || i1==5) dx_[i1] = XVU_coor[chNum][i1][iseg]-0.5*(par_ab[chNum][0][iseg]+sq3*par_ab[chNum][2][iseg])*z_loc[chNum][i1]-0.5*(par_ab[chNum][1][iseg]+sq3*par_ab[chNum][3][iseg]);
					if(i1==1 || i1==4) dx_[i1] = XVU_coor[chNum][i1][iseg]-0.5*(par_ab[chNum][0][iseg]-sq3*par_ab[chNum][2][iseg])*z_loc[chNum][i1]-0.5*(par_ab[chNum][1][iseg]-sq3*par_ab[chNum][3][iseg]);
					Chi2[chNum][iseg] += dx_[i1]*dx_[i1]/(sigm2_[i1]);

					// cout<<"iseg "<<iseg <<" i "<<i1<<" dx_ "<<dx_[i1]<<" coor "<<XVU_coor[chNum][i1][iseg]<<" Chi2 "<<Chi2[chNum][iseg]<<" z "<<z_loc[chNum][i1]<<endl;
				}
			}//i1
			//---Chi2 calculating.---
			//  cout<<" Nhits["<<chNum<<"]["<<iseg<<"] "<<Nhits_Ch_[chNum][iseg]<<" Chi2[chNum][iseg] "<< Chi2[chNum][iseg]<<endl;
			if (Nhits_Ch_[chNum][iseg] > 4) {
				Chi2_ndf[chNum][iseg] = Chi2[chNum][iseg]/(Nhits_Ch_[chNum][iseg]-4);
				//	cout<<" Chi2_ndf[chNum][iseg] "<< Chi2_ndf[chNum][iseg]<<endl;
			}

			if (Chi2_ndf[chNum][iseg] > kChi2_Max_) { // --if bad chi2--

				if (Nhits_Ch_[chNum][iseg] <= Min_hits){// --if no enough points--
					Nhits_Ch_[chNum][iseg] = 0; 
					continue;
				}
				else  { //--reject most distant point--	
					for (Int_t i1 = 0; i1 < kNPlanes; i1++){
						if (XVU_coor[chNum][i1][iseg] > -999. && fabs(dx_[i1]) > delta){
							XVU_coor[chNum][i1][iseg] = -999.; 
							Nhits_Ch_[chNum][iseg]--;// --reject bad point
							continue;
						}//if point was
					}//i1
				}//else

			}// if bad chi2.

			//--if Chi2 too small
			if ( Chi2_ndf[chNum][iseg] < Chi2_Super_min) {
				Nhits_Ch_[chNum][iseg] = 0; 
				continue;
			}//--if Chi2 too small.


			if ( Chi2_ndf[chNum][iseg] > Chi2_Super_min && Chi2_ndf[chNum][iseg] < kChi2_Max_ && Nhits_Ch_[chNum][iseg] != 0){
				Nhits_Ch_[chNum][iseg] = 0;
				for (Int_t i1 = 0; i1 < kNPlanes; i1++){
					if (XVU_coor[chNum][i1][iseg] > -999.) Nhits_Ch_[chNum][iseg]++;
				}
			}


		}//--iseg------------------------------------------------------------------------------------
		if (!ifNhitm) continue;

		//--choice of min chi2 
		Double_t Chi2_best = 999.0;
		Int_t iseg_best = -1;
		for (Int_t iseg = 0; iseg < Nsegm[chNum]; iseg++) {
			if (Nhits_Ch_[chNum][iseg] != Nhitm) continue; 
			if (Chi2_ndf[chNum][iseg] >= Chi2_best) continue;
			Chi2_best = Chi2_ndf[chNum][iseg];
			iseg_best = iseg;
		} // iseg

		if (iseg_best == -1) continue; 

		//-- reject(common points)
		for (Int_t iseg = 0; iseg < Nsegm[chNum]; iseg++) {
			if (iseg == iseg_best)continue;
			for (Int_t i1 = 0; i1 < kNPlanes; i1++) {
				if ( XVU_coor[chNum][i1][iseg] > -999.) { 
					if( fabs(XVU_coor[chNum][i1][iseg] - XVU_coor[chNum][i1][iseg_best]) < 3*dw_half ) 
						Nhits_Ch_[chNum][iseg] = 0; 
				}
			}//i1
		}// iseg


		vector<segments> vtmpSeg;
		segments tmpSeg;

		for (int itSeg = 0; itSeg < Nsegm[chNum]; ++itSeg){
			if (Nhits_Ch_[chNum][itSeg] == Nhitm  && Chi2_ndf[chNum][itSeg] < kChi2_Max_){
				tmpSeg.Nhits = Nhits_Ch_[chNum][itSeg];
				tmpSeg.Chi2  = Chi2_ndf[chNum][itSeg];
				for(int ipla = 0; ipla < kNPlanes; ipla++){
					tmpSeg.coord[ipla] = XVU_coor[chNum][ipla][itSeg];
					tmpSeg.clust[ipla] = Cluster_coor[chNum][ipla][itSeg];
				}
				for(int ipar = 0; ipar < 4; ipar++){
					tmpSeg.param[ipar] = par_ab[chNum][ipar][itSeg];
				}

				vtmpSeg.push_back(tmpSeg);
			}
		}

		// vector sorting
		sort(vtmpSeg.begin(), vtmpSeg.end(), compareSegments);

		// storing
		for (int iterOut = 0; iterOut < vtmpSeg.size(); iterOut++){
			if (OutSegCount < kmaxSeg){
				OutSegArray[OutSegCount] = vtmpSeg.at(iterOut);
				OutSegCount++;
			}
		}

		// !!! vector clear for next Nhitm
		if ( Nhitm == 5) vtmpSeg.clear(); //clear for 4p-segments


	}//Nhitm--

	for (int iterOut = 0; iterOut < kmaxSeg; iterOut++){
		if(OutSegArray[iterOut].Chi2 > 0){
			Nbest_seg_[chNum]++;
			Nhits_seg_[chNum][iterOut] = OutSegArray[iterOut].Nhits;
			Chi2_ndf_seg_[chNum][iterOut] = OutSegArray[iterOut].Chi2;

			for (int iterCoord = 0; iterCoord < kNPlanes; iterCoord++) {
				coor_seg[chNum][iterCoord][iterOut] = OutSegArray[iterOut].coord[iterCoord];
				cluster_seg[chNum][iterCoord][iterOut] = OutSegArray[iterOut].clust[iterCoord];
			}
			for(int ipar = 0; ipar < 4; ipar++){
				Par_ab_seg[chNum][ipar][iterOut] = OutSegArray[iterOut].param[ipar];
			}
		}
	}

	// printf(">>>TSTING\n");
	// for (int iterOut = 0; iterOut < kmaxSeg; iterOut++){
	//   if(OutSegArray[iterOut].Chi2 > 0){
	//     printf("> Ch: %d, Seg %d, Nhits: %d, Chi2: %8.4f\n", chNum, iterOut, OutSegArray[iterOut].Nhits, OutSegArray[iterOut].Chi2);
	//     for (int iterCoord = 0; iterCoord < kNPlanes; iterCoord++) 
	// 	printf("> %8.4f  -  %8.4f\n", OutSegArray[iterOut].coord[iterCoord], OutSegArray[iterOut].clust[iterCoord]);
	//   }
	// }

}//ProcessSegments

void BmnMwpcHitFinder::PrepareArraysToProcessEvent(){
	//fBmnMwpcSegmentsArray->Clear();

	for(Int_t iCh = 0; iCh < kNChambers; iCh++){
		Nseg_Ch[iCh]  = 0;
		Nbest_Ch[iCh] = 0;
		Nbest_seg[iCh] = 0;	
		for(Int_t iPl = 0; iPl < kNPlanes; iPl++){
			iw_Ch[iCh][iPl]  = 0;
			XVU[iCh][iPl]    = 0;
			XVU_cl[iCh][iPl] = 0;
			Nclust[iCh][iPl] = 0;

			for(Int_t iWire=0; iWire<kNWires; iWire++){
				DigitsArray[iCh][iPl][iWire] = 0.;
			}

			for(Int_t iBig=0; iBig<kBig; iBig++){
				Wires_Ch[iCh][iPl][iBig]     = -1;
				clust_Ch[iCh][iPl][iBig]     = -1;
				XVU_Ch[iCh][iPl][iBig]       = -999.;      
				Coord_wire[iCh][iPl][iBig]   = -999.;  
				Coord_xuv[iCh][iPl][iBig]    = -999.; 
				ClusterSize[iCh][iPl][iBig]  = 0;
				XVU_coord[iCh][iPl][iBig]    = -999.; 
				Cluster_coord[iCh][iPl][iBig]= -1;
				Coor_seg[iCh][iPl][iBig]    = -999.; 
				Cluster_seg[iCh][iPl][iBig]= -1;  

			}
			for(Int_t iWire=0; iWire<kNWires; iWire++){
				wire_Ch[iCh][iWire][iPl] = 0;
				xuv_Ch[iCh][iWire][iPl]  = 0.;	
			}
		}

		for(Int_t ii = 0; ii < 4; ii++){
			for(Int_t jj=0; jj < kBig; jj++){
				par_ab_Ch[iCh][ii][jj] = 999.;
				par_ab_seg[iCh][ii][jj] = 999.;
			}
		}

		for(Int_t iBig=0; iBig<kBig; iBig++){     
			Nhits_Ch[iCh][iBig]    = 0;
			Nhits_seg[iCh][iBig]   = 0;
			Chi2_ndf_Ch[iCh][iBig] = 0;
		}

		for(Int_t i=0; i < kmaxSeg; i++){
			ind_best_Ch[iCh][i]      = 0;    
			best_Ch_gl[iCh][i]       = -1;    
			Chi2_ndf_best_Ch[iCh][i] = -999.; 
			Chi2_ndf_seg[iCh][i] = -999.;

		}
	}//iCh

	for(Int_t iPl=0; iPl<kNPlanes; iPl++){     
		sigm2[iPl] = sigma*sigma;  
		ipl[iPl] = 6;
		z2[iPl] = 0;
	}

	for(Int_t ii  = 0; ii < 4; ii++){
		for(Int_t jj  = 0; jj < 4; jj++){
			matrA[ii][jj] = 0.;
			matrb[ii][jj] = 0.;
		}
	}
}//PrepareArraysToProcessEvent


void BmnMwpcHitFinder::FillFitMatrix(Int_t chN, Double_t** AA, Float_t** z, Double_t* sigm2_, Int_t* h_) {
	// AA - matrix to be filledlayers)
	// sigm2 - square of sigma
	// h_ - array to include/exclude planes (h_[i] = 0 or 1)
	// Float_t z2_[nPlanes];
	Float_t z2_[6] = {z[chN][0] * z[chN][0], z[chN][1] * z[chN][1], z[chN][2] * z[chN][2], z[chN][3] * z[chN][3], z[chN][4] * z[chN][4], z[chN][5] * z[chN][5]}; //cm

	AA[0][0] += 2 * z2_[0] * h_[0] / sigm2_[0]  +      z2_[2] * h_[2] / (2 * sigm2_[2])   +      z2_[1] * h_[1] / (2 * sigm2_[1])  +  2 * z2_[3] * h_[3] / sigm2_[3] +      z2_[5] * h_[5] / (2 * sigm2_[5])  +      z2_[4] * h_[4] / (2 * sigm2_[4]); //Ax

	AA[0][1] += 2 * z[chN][0] * h_[0] / sigm2_[0]  +      z[chN][2] * h_[2] / (2 * sigm2_[2])  +      z[chN][1] * h_[1] / (2 * sigm2_[1])  +  2 * z[chN][3] * h_[3] / sigm2_[3]  +      z[chN][5] * h_[5] / (2 * sigm2_[5])  +      z[chN][4] * h_[4] / (2 * sigm2_[4]); //Bx

	AA[0][2] += sq3 * (z2_[2] * h_[2] / (2 * sigm2_[2])  -         z2_[1] * h_[1] / (2 * sigm2_[1])  +         z2_[5] * h_[5] / (2 * sigm2_[5])  -         z2_[4] * h_[4] / (2 * sigm2_[4])); //Ay

	AA[0][3] += sq3 * (z[chN][2] * h_[2] / (2 * sigm2_[2])   -         z[chN][1] * h_[1] / (2 * sigm2_[1])  +         z[chN][5] * h_[5] / (2 * sigm2_[5])  -         z[chN][4] * h_[4] / (2 * sigm2_[4])); //By

	AA[1][0] = AA[0][1];

	AA[1][1] +=   2 * h_[0] / sigm2_[0]  +  0.5 * h_[2] / sigm2_[2] + 0.5 * h_[1] / sigm2_[1]  +    2 * h_[3] / sigm2_[3] + 0.5 * h_[5] / sigm2_[5]  +  0.5 * h_[4] / sigm2_[4];

	AA[1][2] += sq3 * (z[chN][2] * h_[2] / sigm2_[2]  - z[chN][1] * h_[1] / sigm2_[1]  + z[chN][5] * h_[5] / sigm2_[5]  - z[chN][4] * h_[4] / sigm2_[4]) * 0.5;

	AA[1][3] += sq3 * (h_[2] / sigm2_[2]   -  h_[1] / sigm2_[1]  +   h_[5] / sigm2_[5]  -  h_[4] / sigm2_[4]) * 0.5;

	AA[2][0] = AA[0][2];

	AA[2][1] = AA[1][2];

	AA[2][2] += 3.0 * (z2_[2] * h_[2] / sigm2_[2]   +   z2_[1] * h_[1] / sigm2_[1]  +  z2_[5] * h_[5] / sigm2_[5]   +   z2_[4] * h_[4] / sigm2_[4]) * 0.5;

	AA[2][3] += 3.0 * (z[chN][2] * h_[2] / sigm2_[2]  +    z[chN][1] * h_[1] / sigm2_[1]   +  z[chN][5] * h_[5] / sigm2_[5]  +  z[chN][4] * h_[4] / sigm2_[4])   * 0.5;

	AA[3][0] = AA[0][3];
	AA[3][1] = AA[1][3];
	AA[3][2] = AA[2][3];
	AA[3][3] += 3.0 * (0.5 * h_[2] / sigm2_[2]  +  0.5 *  h_[1] / sigm2_[1]   +  0.5 *  h_[5] / sigm2_[5]  +  0.5 *  h_[4] / sigm2_[4]);
}


void BmnMwpcHitFinder::FillFreeCoefVector(Int_t ichNum, Double_t* F, Double_t*** XVU_, Int_t ise, Float_t** z, Double_t* sigmm2, Int_t* h_) {
	// F - vector to be filled
	// XVU_ - coordinates of segment in chamber (Is it correct definition?)
	// segIdx - index of current segment
	// z - local z-positions of planes(layers)
	// sigmm2 - square of sigma
	// h_ - array to include/exclude planes (h_[i] = 0 or 1)


	F[0] += 2 * XVU_[ichNum][0][ise] * z[ichNum][0] * h_[0] / sigmm2[0]  + XVU_[ichNum][1][ise]  * z[ichNum][1] * h_[1] / sigmm2[1] +  XVU_[ichNum][2][ise] * z[ichNum][2] * h_[2] / sigmm2[2] + 2 * XVU_[ichNum][3][ise] * z[ichNum][3] * h_[3] / sigmm2[3] + XVU_[ichNum][4][ise] * z[ichNum][4] * h_[4] / sigmm2[4] 
		+  XVU_[ichNum][5][ise] * z[ichNum][5] * h_[5] / sigmm2[5];

	F[1] += 2 * XVU_[ichNum][0][ise] * h_[0] / sigmm2[0] + XVU_[ichNum][1][ise] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * h_[2] / sigmm2[2] + 2 * XVU_[ichNum][3][ise] * h_[3] / sigmm2[3] + XVU_[ichNum][4][ise] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * h_[5] / sigmm2[5];

	F[2] += sq3*(-XVU_[ichNum][1][ise] * z[ichNum][1] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * z[ichNum][2] * h_[2] / sigmm2[2] - XVU_[ichNum][4][ise] * z[ichNum][4] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * z[ichNum][5] * h_[5] / sigmm2[5]);

	F[3] +=  sq3*(-XVU_[ichNum][1][ise] * h_[1] / sigmm2[1] + XVU_[ichNum][2][ise] * h_[2] / sigmm2[2] - XVU_[ichNum][4][ise] * h_[4] / sigmm2[4] + XVU_[ichNum][5][ise] * h_[5] / sigmm2[5]);

}

void BmnMwpcHitFinder::FillFreeCoefVectorXUV(Int_t ichNum ,  Double_t* F, Float_t** XVU_, Float_t** z, Float_t* sigm2_, Int_t* h_) {
	// F - vector to be filled
	// XVU_ - coordinates of segment in chamber (Is it correct definition?)
	// segIdx - index of current segment
	// z - local z-positions of planes(layers)
	// sigm2_ - square of sigma
	// h_ - array to include/exclude planes (h_[i] = 0 or 1)
	F[0] +=  2 * XVU_[ichNum][0]  * z[ichNum][0] * h_[0] / sigm2_[0]  +  XVU_[ichNum][1]   * z[ichNum][1] * h_[1] / sigm2_[1]  + XVU_[ichNum][2]  * z[ichNum][2] * h_[2] / sigm2_[2] + 2 * XVU_[ichNum][3]  * z[ichNum][3] * h_[3] / sigm2_[3]  +  XVU_[ichNum][4]  * z[ichNum][4] * h_[4] / sigm2_[4]  + XVU_[ichNum][5]  * z[ichNum][5] * h_[5] / sigm2_[5];
	F[1] += 2 * XVU_[ichNum][0]  * h_[0] / sigm2_[0] + XVU_[ichNum][1]  * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * h_[2] / sigm2_[2] + 2 * XVU_[ichNum][3]  * h_[3] / sigm2_[3] + XVU_[ichNum][4]  * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * h_[5] / sigm2_[5];
	F[2] += (-XVU_[ichNum][1]  * z[ichNum][1] * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * z[ichNum][2] * h_[2] / sigm2_[2] - XVU_[ichNum][4]  * z[ichNum][4] * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * z[ichNum][5] * h_[5] / sigm2_[5]);
	F[3] +=  (-XVU_[ichNum][1]  * h_[1] / sigm2_[1] + XVU_[ichNum][2]  * h_[2] / sigm2_[2] - XVU_[ichNum][4]  * h_[4] / sigm2_[4] + XVU_[ichNum][5]  * h_[5] / sigm2_[5]);

	F[2]=F[2]*sq3;
	F[3]=F[3]*sq3;
}


void BmnMwpcHitFinder::InverseMatrix(Double_t** AA, Double_t** bb) {
	// Gaussian algorithm for 4x4 matrix inversion 
	Double_t factor;
	Double_t temp[4];
	// Set b to I
	for (Int_t i1 = 0; i1 < 4; i1++){
		for (Int_t j1 = 0; j1 < 4; j1++){
			if (i1 == j1) bb[i1][j1] = 1.0;
			else bb[i1][j1] = 0.0;
		}
	}
	for (Int_t i1 = 0; i1 < 4; i1++) {
		for (Int_t j1 = i1 + 1; j1 < 4; j1++) {
			if (fabs(AA[i1][i1]) < fabs(AA[j1][i1])) {
				for (Int_t l1 = 0; l1 < 4; l1++) temp[l1] = AA[i1][l1];
				for (Int_t l1 = 0; l1 < 4; l1++) AA[i1][l1] = AA[j1][l1];
				for (Int_t l1 = 0; l1 < 4; l1++) AA[j1][l1] = temp[l1];
				for (Int_t l1 = 0; l1 < 4; l1++) temp[l1] = bb[i1][l1];
				for (Int_t l1 = 0; l1 < 4; l1++) bb[i1][l1] = bb[j1][l1];
				for (Int_t l1 = 0; l1 < 4; l1++) bb[j1][l1] = temp[l1];
			}
		}
		factor = AA[i1][i1];
		for (Int_t j1 = 4 - 1; j1>-1; j1--) {
			bb[i1][j1] /= factor;
			AA[i1][j1] /= factor;
		}
		for (Int_t j1 = i1 + 1; j1 < 4; j1++) {
			factor = -AA[j1][i1];
			for (Int_t k1 = 0; k1 < 4; k1++) {
				AA[j1][k1] += AA[i1][k1] * factor;
				bb[j1][k1] += bb[i1][k1] * factor;
			}
		}
	} // i1
	for (Int_t i1 = 3; i1 > 0; i1--) {
		for (Int_t j1 = i1 - 1; j1>-1; j1--) {
			factor = -AA[j1][i1];
			for (Int_t k1 = 0; k1 < 4; k1++) {
				AA[j1][k1] += AA[i1][k1] * factor;
				bb[j1][k1] += bb[i1][k1] * factor;
			}
		}
	} // i1
}    //end inverse

void BmnMwpcHitFinder::Finish() {

	if (fDebug){
		printf("MWPC hit finder: write hists to file... ");
		fOutputFileName = Form("hMWPChits_p%d_run%d.root", fRunPeriod, fRunNumber);
		cout<< fOutputFileName <<endl;
		TFile file(fOutputFileName, "RECREATE");
		if(fDoTest) fList.Write();
		file.Close();
	}

	if (fDebug) printf("done\n");
	delete fMwpcGeometrySRC;

	cout << "Work time of the MWPC hit finder: " << workTime << " s" << endl;
}

ClassImp(BmnMwpcHitFinder)


