// @(#)bmnroot/mwpc:$Id$


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// BmnMwpcTrackFinder                                                         //
//                                                                            //
//                                                                            //
// The algorithm serves for searching for track segments                      //
// in the MWPC of the BM@N experiment                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef BMNMWPCTRACKFINDER_H
#define BMNMWPCTRACKFINDER_H 


#include <TMath.h>
#include <TNamed.h>
#include <TClonesArray.h>
#include <TString.h>
#include "FairTask.h"
#include "BmnMwpcTrackToDC.h"
#include "BmnMwpcTrack.h"
#include "BmnMwpcHit.h"
#include "BmnMwpcGeometrySRC.h"
#include "BmnEnums.h"
#include "BmnMath.h"
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <sstream>

class TH1D;
class TH2D;

using namespace std;
using namespace TMath;

class BmnMwpcTrackFinder : public FairTask {
public:

    BmnMwpcTrackFinder() {};
    BmnMwpcTrackFinder(Bool_t, Int_t, Int_t);
    virtual ~BmnMwpcTrackFinder();

    virtual InitStatus Init();

    virtual void Exec(Option_t* opt, TClonesArray* fBmnMwpcSegmentsArray, TClonesArray* fBmnMwpcTracksArray);

    virtual void Finish();

private:
    Bool_t expData;
    Bool_t fDebug = 0;
    UInt_t fEventNo; // event counter
    Int_t fRunPeriod;
    Int_t fRunNumber;

    TString fInputBranchName;
    TString fOutputBranchName;
    TString fOutputBranchNameToDC;
    TString fOutputFileName;
   
    /** Input array of MWPC hits **/
    //TClonesArray* fBmnMwpcSegmentsArray;
    
    /** Output array of MWPC tracks **/
    //TClonesArray* fBmnMwpcTracksArray; 
    //TClonesArray* fBmnMwpcTracksArrayToDC; 
    
    BmnMwpcGeometrySRC* fMwpcGeo;

    vector<TH1D*> hpar_Ax_Ch, hpar_Bx_Ch, hpar_Ay_Ch, hpar_By_Ch, hChi2_match_pair, hChi2_ndf_Ch, hpar_Ax_pair, hpar_Bx_pair, hpar_Ay_pair, hpar_By_pair, hpar_theta_pair, hpar_phi_pair, Nomin_Ch, Denom_Ch, Eff_Ch, hX_in_target_pair, hY_in_target_pair, hdX_Zmid_pair, hdY_Zmid_pair, hdAx_Zmid_pair, hdAy_Zmid_pair, hdU_zmid_pair, hdV_zmid_pair, hdX_zmid_pair, hoccupancyXp, hoccupancyUp, hoccupancyVp,hoccupancyXm, hoccupancyUm, hoccupancyVm;
    vector<TH2D*> hAx_bx_in_target_pair, hAy_by_in_target_pair, hY_X_in_target_pair ;
    TH1D *hdX_target, *hdY_target, *hX_in_target, *hY_in_target, *hdAx_target, *hdAy_target, *hChi2_m_target, *hdX_pair01_inZpair1, *hdY_pair01_inZpair1, *hpar_Ax_3ch, *hpar_Bx_3ch, *hpar_Ay_3ch, *hpar_By_3ch;
    TH2D *hAx_bx_in_target, *hAy_by_in_target, *hY_X_in_target, *hdX_pair01_vs_x1, *hdY_pair01_vs_y1, *htheta_p1vsp0;

    Short_t kNChambers;
    Short_t kNPlanes;
    Int_t kBig;
    Int_t kCh_min;
    Int_t kCh_max;
    Int_t kNumPairs;
    TVector3 *ChCent;

    Float_t *kZmid;
    Float_t *ZCh;
    Float_t **kZ_loc;
    Float_t *kZ_midle_pair;
 
    Float_t **z_gl;
    Float_t **shift; 
    Float_t **shift_pair;
    
    Float_t kZ_to_pole;
    Float_t kZ_target;
    Float_t kZ_DC;
    Int_t kMinHits;
    Int_t kmaxPairs;
    Double_t kChi2_Max;

    Float_t dw;
    Float_t dw_half;
    Double_t sq3;
    Double_t sq12;
    Double_t sigma;
    Short_t kMiddlePl;
    Int_t *Nbest_3ch;
  
    Int_t **kPln;
    Float_t *sigm2;
    Float_t **sigma_delta;
    Int_t *ipl;
    Double_t **matrA;
    Double_t **matrb;
    Double_t **Amatr;
    Double_t **bmatr;

    Double_t ***par_ab_Ch;
    Double_t ***par_ab_pair;
    Double_t ***XVU_Ch;
    Double_t ***Clust_Ch;
    Double_t **Chi2_match_pair;
    Double_t **Chi2_ndf_pair;
    Double_t **Chi2_ndf_Ch;
    Int_t **ind_best_pair;
    Int_t **ind_best_Ch;
    Int_t **Nhits_match;
    Int_t *Nbest_pair;
    Int_t *Nbest_Ch;
    Double_t **par_ab_3ch;
    Int_t **Nhits_Ch;
    Int_t **Nhits_pair;

    TList fList;

    void PrepareArraysToProcessEvent();
    void SegmentParamAlignment(Int_t, Int_t *,  Double_t ***, Float_t **);
    
    void SegmentMatching(Int_t,Int_t *, Double_t ***, Float_t *, Int_t **,  Int_t *,Double_t **,Double_t ***, Int_t **, Int_t **);
    void SegmentFit(Int_t, Float_t **, Float_t *,Int_t *, Int_t **, Double_t ***,  Double_t **, Double_t ***, Double_t ***, Int_t **,Int_t **, Int_t **);
    void FillFitMatrix(Int_t, Double_t** , Float_t** , Double_t* , Int_t*);
    void FillFreeCoefVector(Int_t, Double_t*, Double_t***, Int_t, Float_t**, Double_t*, Int_t*);

    void InverseMatrix(Double_t**, Double_t**);
    void PairMatching( Int_t *, Double_t ***,  Float_t *);
    void FillEfficiency(Int_t, Double_t ***, Int_t **, Int_t, Int_t, Float_t, Float_t);
    void StraightByTwoPoints(Double_t ***, Int_t *, Double_t ***, Int_t *, Float_t *, Double_t **, Int_t *);
         
    ClassDef(BmnMwpcTrackFinder, 1)
};

#endif
