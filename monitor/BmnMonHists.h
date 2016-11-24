/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   BmnMonHists.h
 * Author: ilnur
 *
 * Created on November 15, 2016, 12:57 PM
 */
#ifndef BMNMONHISTS_H
#define BMNMONHISTS_H 1
#include <TNamed.h>
#include <exception>
#include <stdlib.h>
#include <vector>
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
        
#define MAX_STATIONS 40
#define MAX_MODULES 8
#define MAX_LAYERS 4
using namespace std;

class BmnMonHists : public TNamed {
public:
    
    BmnMonHists();
    BmnMonHists(const BmnMonHists& orig);
    virtual ~BmnMonHists();
    void Clear();
    void ClearTriggers();
    void ClearGEM();
    void ClearToF400();
    void ClearToF700();
    vector<vector<vector<TH1F*> > > histGemStrip;
//    TH1F *histGemStrip[MAX_STATIONS][MAX_MODULES][MAX_LAYERS];
    TH1D *histToF400LeadingTime;
    TH1D *histToF400Amp;
    TH1I *histToF400Strip;
    TH1I *histToF400StripSimult;
    TH2F *histToF400State;
    TH1D *histBC1TimeLen;
    TH1D *histBC2TimeLen;
    TH1D *histSDTimeLen;
    TH1D *histVDTimeLen;
    TH1D *histFDTimeLen;
    TH1D *histBDTimeLen;
    TH1I *histTriggers;
        
private:
ClassDef(BmnMonHists, 1)
};

#endif /* BMNMONHISTS_H */

