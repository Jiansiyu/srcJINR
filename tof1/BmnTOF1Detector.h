/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BmnTOF1Detector.h
 * Author: mikhail
 *
 * Created on March 13, 2017, 4:52 PM
 */

#ifndef BMNTOF1DETECTOR_H
#define BMNTOF1DETECTOR_H 1

#include "TString.h"
#include "TSystem.h"
#include "BmnEnums.h"
#include "BmnTof1Digit.h"
#include "BmnEventHeader.h"
#include "BmnTrigDigit.h"
#include "BmnRunHeader.h"
#include "BmnTOF1Point.h"
#include "BmnTof1GeoUtils.h"
#include "BmnTofHit.h"
#include "BmnTOF1Conteiner.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "TVectorT.h"
#include "TDirectory.h"
#include <TGeoManager.h>
#include <TKey.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <map>
#include <deque>
#include <iterator>

using namespace std;

class BmnTOF1Detector {
	private:
		//////////////////////////////////////////////////////////////////////////////////
		// New holders for class
		static const Int_t fNStr = 48; // number of strips in ToF400
		TString fName; // name of class, it's just Plane_#
		int fNPlane;   // same as the nPlane from contructor
		int fFillHist; // what level to fill histograms
		double fStripLength, fSignalVelocity;		// strip length and velocity
		double fMaxDelta;				// how much leeway we allow for each strip in terms of length
		double fCorrLR[fNStr];
		double fStripShift[fNStr];
		double fWalkFunc[fNStr][4];
		bool fKilled[fNStr];
		double fGammaOffset[fNStr];
		TVector3 fCenterStrip[fNStr];

		// For initial processing
		double tempHitTime[fNStr][5];
		double tempHitAmps[fNStr][5];
		int tempCounter[fNStr];
		
		// For hit matching -- after this, there is
		// one hit for every strip (although most strips
		// don't fire of course)
		bool fFlagHit[fNStr];
		double fTof[fNStr];
		double fAmp[fNStr];
		double fYPos[fNStr];
		int fXPos[fNStr];
		
		// Storage for how many strips fired after hit matching
		int stripsFired[fNStr];
		int numStripsFired;

		// Holder for all possible clusters -- could be up to 48 of them
		int 	numClusters = 0;
		double 	final_Tof[fNStr];
		double 	final_Amp[fNStr];
		double 	final_YPos[fNStr];
		int 	final_XPos[fNStr];
		


		// Old holders that I'm phasing out
		TVector3 fCrossPoint[fNStr], fVectorTemp;
		BmnTrigDigit *fT0;

		Double_t CalculateDt(Int_t Str){ return 0.;};
		Bool_t GetCrossPoint(Int_t NStrip);
		void AddHit(Int_t Str, TClonesArray *TofHit){return;};
		void AddConteiner(Int_t Str, TClonesArray *TofHit){ return;};


	public:
		//////////////////////////////////////////////////////////////////////////////////
		// Constructors:
		BmnTOF1Detector(); //empty constructor
		BmnTOF1Detector(int NPlane, int FillHist, TTree *tree);

		//////////////////////////////////////////////////////////////////////////////////
		// New functions that I'm using:
		void SetCorrLR( string pathToFile );		// For single plane of ToF400, sets LR correction
		void SetStripShift( string pathToFile );	// For single plane of ToF400, sets strip shift
		void SetWalkFunc( string pathToFile );		// For single plane of ToF400, sets walk functions
		void SetGammaOffset( string pathToFile );	// For single plane of ToF400, moves gamma peak to correct location
		void TestPrint( int strip );			// For single plane of ToF400, prints out correction files loaded for given strip
		void SetGeoFile( string pathToFile );		// For single plane of ToF400, loads geometry file and sets strip centers in global coords

		void ClearHits();

		void InitSkim( BmnTof1Digit* tofDigit );
		void CreateStripHit( BmnTof1Digit* tofDigit , double t0Time , double t0Amp );
		void ClusterHits();
		void OutputToTree();
	
		double GetStripMult(){ return numStripsFired; };

		// OLD FUNCTIONS THAT I HAVE PHASED OUT
		Bool_t SetCorrLR( TString NameFile ){return kTRUE;};
		Bool_t SetCorrLR(Double_t *Mass){return kTRUE;};
		Bool_t SetCorrSlewing( TString NameFile ){return kTRUE;};
		Bool_t SetCorrTimeShift( TString NameFile ){return kTRUE;};
		Bool_t SetGeoFile(TString NameFile){return kTRUE;};
		Bool_t SaveHistToFile(TString NameFile){return kTRUE;};
		Int_t GetFillHistLevel(){return fFillHist;};
		void KillStrip(Int_t NumberOfStrip){return;};
		Bool_t SetGeo(BmnTof1GeoUtils *pGeoUtils){return kTRUE;};
		Bool_t GetXYZTime(Int_t Str, TVector3 *XYZ, Double_t *ToF){return kTRUE;};
		Double_t GetWidth(Int_t Str){return 0.;};
		Double_t GetTime(Int_t Str){return 0.;};
		Int_t FindHits(BmnTrigDigit *T0){return 0;};
		Bool_t SetDigit(BmnTof1Digit *TofDigit){return kTRUE;};
		void Clear(){return;}; // clear all holders for ToF400
		Int_t FindHits(BmnTrigDigit *T0, TClonesArray *TofHit) {return 0;};


		virtual ~BmnTOF1Detector() {}; // deconstructor

		ClassDef(BmnTOF1Detector, 4);

};

#endif
