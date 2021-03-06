#ifndef BMNGEMSTRIPHITMAKER_H
#define BMNGEMSTRIPHITMAKER_H 1

#include <iostream>
#include <fstream>
#include <sstream>

#include "Rtypes.h"
#include "TClonesArray.h"
#include "TRegexp.h"
#include "TString.h"

#include "FairTask.h"
#include "FairMCPoint.h"

#include <UniDbDetectorParameter.h>
#include <UniDbRun.h>
#include "FairField.h"
#include "BmnGemStripDigit.h"
#include "BmnGemStripHit.h"
#include "BmnGemStripStationSet.h"
#include "BmnGemStripConfiguration.h"
#include "BmnGemStripTransform.h"
#include "BmnGemAlignmentCorrections.h"
#include "BmnGemAlignCorrections.h"
#include <BmnEventQuality.h>

using namespace std;

class BmnGemStripHitMaker : public FairTask {
public:

    BmnGemStripHitMaker();
    BmnGemStripHitMaker(Int_t, Bool_t);

    virtual ~BmnGemStripHitMaker();

    //void SetVerbosity(Bool_t verbose);

    virtual InitStatus Init( FairField *fField );

    virtual void Exec(Option_t* opt, TClonesArray* fBmnGemStripDigitsArray, TClonesArray* fBmnGemStripHitsArray, TClonesArray*  fBmnGemStripHitMatchesArray, FairField* fField);

    virtual void Finish();

	void ProcessDigits( TClonesArray* fBmnGemStripDigitsArray, TClonesArray* fBmnGemStripHitsArray, TClonesArray* fBmnGemStripHitMatchesArray, FairField* fField);

    void SetHitMatching(Bool_t opt = kTRUE) {
        fHitMatching = opt;
    }

    void SetCurrentConfig(BmnGemStripConfiguration::GEM_CONFIG config) {
        fCurrentConfig = config;
    }

    void SetAlignmentCorrectionsFileName(TString filename) {
        fAlignCorrFileName = filename;
    }

    void SetAlignmentCorrectionsFileName(Int_t file_number) {
        fRunId = file_number;
        fAlignCorrFileName = "alignment_GEM.root";
        UniDbDetectorParameter::ReadRootFile(fPeriodId, file_number, "BM@N", "alignment", (Char_t*) fAlignCorrFileName.Data());
    }

    inline Double_t GetLorentzByField(Double_t By, Int_t station) {
        // Appropriate Lorentz-corrections to the GEM-hits for RUN5, 6 as a function of voltage and so on
        // have been moved to the UniDb
        return lorCorrsCoeff[station][0] + lorCorrsCoeff[station][1] * By + lorCorrsCoeff[station][2] * By * By;
    }

    void SetSigmaX(Double_t sigX) {
        fSigmaX = sigX;
    }

    void SetSigmaY(Double_t sigY) {
        fSigmaY = sigY;
    }

    void SetSigmaZ(Double_t sigZ) {
        fSigmaZ = sigZ;
    }

    void SetSigmaXYZ(Double_t sigX, Double_t sigY, Double_t sigZ) {
        fSigmaX = sigX;
        fSigmaY = sigY;
        fSigmaZ = sigZ;
    }

private:

    TString fInputPointsBranchName;
    TString fInputDigitsBranchName;
    TString fInputDigitMatchesBranchName;

    TString fOutputHitsBranchName;
    TString fOutputHitMatchesBranchName;

    /** Input array of Gem Points **/
    TClonesArray* fBmnGemStripPointsArray;
    //TClonesArray* fBmnGemStripDigitsArray;
    TClonesArray* fBmnGemStripDigitMatchesArray;

    /** Output array of Gem Hits **/
    //TClonesArray* fBmnGemStripHitsArray;

    /** Output array of GEM Hit Matches **/
    //TClonesArray* fBmnGemStripHitMatchesArray;

    Bool_t fHitMatching;
    Int_t fRunId;
    Int_t fPeriodId;
    Bool_t fIsExp; // Specify type of input data (MC or real data)

    BmnGemStripConfiguration::GEM_CONFIG fCurrentConfig;

    BmnGemStripStationSet *StationSet; //Entire GEM detector

    BmnGemStripTransform *TransfSet; //Transformations for each module of the detector

    TString fAlignCorrFileName; // a file with geometry corrections
    void ReadAlignCorrFile(TString, Double_t***); // read corrections from the file
    Double_t*** corr; // array to store the corrections

    Double_t*** misAlign; // an array to introduce remain misalignment
    //FairField* fField;
    Double_t** lorCorrsCoeff;

    TString fBmnEvQualityBranchName;
    TClonesArray* fBmnEvQuality;

    Double_t fSigmaX;
    Double_t fSigmaY;
    Double_t fSigmaZ;

    ClassDef(BmnGemStripHitMaker, 1);
};

#endif
