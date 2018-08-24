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

#ifndef BMNMWPCDETECTOR_H
#define BMNMWPCDETECTOR_H 1

#include "TString.h"
#include "TSystem.h"
#include "BmnEnums.h"
#include "BmnEventHeader.h"
#include "BmnTrigDigit.h"
#include "BmnRunHeader.h"
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

class BmnMWPCDetector {
	private:
	public:
		BmnMWPCDetector(); //empty constructor
		BmnMWPCDetector(int NPlane, int FillHist, TTree *tree);
		virtual ~BmnMWPCDetector() {}; // deconstructor

		ClassDef(BmnMWPCDetector, 1);

};

#endif
