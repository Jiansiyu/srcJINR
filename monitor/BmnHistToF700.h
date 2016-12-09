/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BmnHistToF700.h
 * Author: ilnur
 *
 * Created on November 29, 2016, 1:24 PM
 */

#ifndef BMNHISTTOF700_H
#define BMNHISTTOF700_H 1
#include <TNamed.h>
#include <exception>
#include <stdlib.h>
#include <vector>
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TFolder.h"
#include "THttpServer.h"
#include "BmnEventHeader.h"

#include "BmnTof2Digit.h"
#include "BmnTof2Raw2Digit.h"

class BmnHistToF700 : public TNamed {
public:
    BmnHistToF700(TString title);
    virtual ~BmnHistToF700();
    void Register(THttpServer *serv);
    void SetDir(TFile *outFile = NULL, TTree *recoTree = NULL);
    void Reset();
    
    TClonesArray *Events = NULL;
    TH1D *histLeadingTime = NULL;
    TH1D *histLeadingTimeSpecific = NULL;
    TH1D *histAmp = NULL;
    TH1D *histAmpSpecific = NULL;
    TH1I *histStrip = NULL;
    TH1I *histStripSimult = NULL;
    TH2F *histState = NULL;
    
    void FillFromDigi(TClonesArray * ToF4Digits, BmnEventHeader * head, Int_t iEv);
    void SetSelection(Int_t Plane, Int_t Strip);
    void SaveHists();
    
    void  SetPlane(Int_t v){ fSelectedPlane = v;}
    void  SetSide (Int_t v){ fSelectedSide  = v;}
    void  SetStrip(Int_t v){ fSelectedStrip = v;}
    Int_t GetPlane() { return fSelectedPlane;}
    Int_t GetSide () { return fSelectedSide ;}
    Int_t GetStrip() { return fSelectedStrip;}    
private:
    THttpServer *fServer;
    TBranch * fEventsBranch;
    Int_t fSelectedPlane;
    Int_t fSelectedStrip;
    Int_t fSelectedSide;
    TTree *frecoTree;
    TH1I *histL = new TH1I("", "", TOF2_MAX_STRIPS_IN_CHAMBER, 0, TOF2_MAX_STRIPS_IN_CHAMBER);
    TH1I *histR = new TH1I("", "", TOF2_MAX_STRIPS_IN_CHAMBER, 0, TOF2_MAX_STRIPS_IN_CHAMBER);
    TH1I histSimultaneous;
    TString pathToImg = "/home/ilnur/Documents/BmnMonJS/public_html/img/";
    
void SaveHist(TH1 *hist, TString path) {
    TCanvas *c0 = new TCanvas("c0", hist->GetTitle());
    hist->Draw();
    c0->SaveAs(TString(path + TString(hist->GetTitle()) + ".png").Data());
    delete c0;
}
    
    ClassDef(BmnHistToF700, 1)
};
#endif /* BMNHISTTOF700_H */

