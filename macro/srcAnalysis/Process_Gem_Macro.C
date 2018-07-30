/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Process_Gem_Macro.C
 * Author: mikhailr
 *
 * Created on July 10, 2018, 12:54 PM
 */

#include <cstdlib>
#include <TAttMarker.h>
#include <TH1.h>
//#include "../../run/bmnloadlibs.C"

using namespace TMath;
using namespace std;

const Double_t cos15 = cos(15 * TMath::Pi() / 180.); //const Double_t cos15 = cos(15.4*TMath::Pi()/180.);
const Double_t sin15 = sin(15 * TMath::Pi() / 180.); //const Double_t sin15 = sin(15.4*TMath::Pi()/180.);

const Int_t nstations = 4;
const Int_t kBig = 200;
const Int_t nstrips = 935; //
const Double_t deltaX = 0.08; //cm
const Double_t deltaXp = deltaX; //cm

Double_t StripLenght[nstations] = {41., 41., 41., 41.}; //cm

const Int_t NDet = 20;

Double_t FindClusterCenter(Double_t* Ampl, Int_t nElements, Double_t &SumAmpll) {

    Double_t XQ = 0.;
    Double_t Q = 0.;

    for (Int_t i = 0; i < nElements; ++i) {
        Q += Ampl[i];
        //  cout<<" Ampl["<<i<<"] ="<<Ampl[i]<<endl;
        XQ += Double_t(i + 0.5) * Ampl[i]; //if strip start 1 //(i+0.5)*Ampl[i] -> if strip start 0
    }

    SumAmpll = Q;
    Double_t CoG = XQ / Q;
    //cout << "ClasterSize = " << nElements << "; Q = " << Q << "; Q/size = " << Q/nElements << endl;

    return (CoG);
}

void GetXYspatial(Int_t* NClX, Int_t* NClXp, Double_t** XCoor, Double_t** XpCoor, Double_t Xsp[nstations][kBig], Double_t Ysp[nstations][kBig], Int_t NXYsp[nstations]) {


    Double_t YCoor_cand = 0.;

    for (Int_t Iz = 0; Iz < nstations; Iz++) {
        // cout<<" NClX["<<Iz<<"] "<<NClX[Iz]<<" NClXp["<<Iz<<"] "<<NClXp[Iz]<<endl;

        for (Int_t clx = 0; clx < NClX[Iz]; ++clx) {
            for (Int_t clxp = 0; clxp < NClXp[Iz]; ++clxp) {
                if (XCoor[Iz][clx] > -900. && XpCoor[Iz][clxp] > -900.) {

                    //YCoor_cand = XpCoor[Iz][clxp] * sin15 + cos15 / sin15 * XCoor[Iz][clx]; //+15
                    YCoor_cand = (XpCoor[Iz][clxp] - cos15 * XCoor[Iz][clx]) / (sin15); //+15

                    //   cout<<"   XCoor["<<Iz<<"]["<<clx<<"]= "<<XCoor[Iz][clx]<<" XpCoor["<<Iz<<"]["<<clx<<"]= "<<XpCoor[Iz][clxp]<<" Ycand "<<YCoor_cand<<endl;


                    if (YCoor_cand >= 0 && YCoor_cand <= StripLenght[Iz]) {

                        Xsp[Iz][NXYsp[Iz]] = XCoor[Iz][clx];
                        Ysp[Iz][NXYsp[Iz]] = YCoor_cand;

			//	cout<<" zone "<< Iz <<" Xsp "<<Xsp[Iz][NXYsp[Iz]]<<" Ysp "<<Ysp[Iz][NXYsp[Iz]]<<endl;
                        NXYsp[Iz]++;

                    }

                    //   YCoor[iSt][iCl] = (   YCoord[iSt][iCl]  - cos(15*TMath::Pi()/180)* X_param[iSt] )/ sin(15*TMath::Pi() / 180);// +15 
                }
            }//NClXp
        }//NClX
        // cout<<" NXYsp= "<<NXYsp[Iz]<<endl;
    }// Iz


}// GetYortFromXSegment

void rotateAxis(Double_t *XYZ, Double_t *angleXYZ, Double_t *outXYZ) {


    Double_t X, Y, Z, XX, YY, ZZ, AngleXRad, AngleYRad, AngleZRad;
    AngleXRad = angleXYZ[0] * TMath::DegToRad();
    AngleYRad = angleXYZ[1] * TMath::DegToRad();
    AngleZRad = angleXYZ[2] * TMath::DegToRad();
    X = XYZ[0];
    Y = XYZ[1];
    Z = XYZ[2];
    //z- beam, y - vertical
    //    cout << "Before rotation \n";
    //    cout << "XYZ:" << XX << " : " << YY << " : " << ZZ << "\n";
    //cout << "angX:" << AngleXRad << " = " << angleXYZ[1] << " * " << TMath::DegToRad() << " : " << TMath::Sin(AngleXRad) << "\n";
    //cout << "angY:" << AngleYRad << " = " << angleXYZ[0] << " * " << TMath::DegToRad() << " : " << TMath::Cos(AngleYRad) << "\n";

    //Rotate around X
    XX = X;
    YY = Y * TMath::Cos(AngleXRad) + Z * TMath::Sin(AngleXRad);
    ZZ = -Y * TMath::Sin(AngleXRad) + Z * TMath::Cos(AngleXRad);
    //cout << "After Xrotation \n";
    //cout << "XYZ:" << XX << " : " << YY << " : " << ZZ << "\n";
    X = XX;
    Y = YY;
    Z = ZZ;

    //Rotate around Y
    XX = X * TMath::Cos(AngleYRad) - Z * TMath::Sin(AngleYRad);
    YY = Y;
    ZZ = X * TMath::Sin(AngleYRad) + Z * TMath::Cos(AngleYRad);
    //cout << "After Yrotation \n";
    //cout << "XYZ:" << XX << " : " << YY << " : " << ZZ << "\n";
    X = XX;
    Y = YY;
    Z = ZZ;

    //Rotate around Z
    XX = X * TMath::Cos(AngleZRad) + Y * TMath::Sin(AngleZRad);
    YY = -X * TMath::Sin(AngleZRad) + Y * TMath::Cos(AngleZRad);
    ZZ = Z;

    //cout << "After Zrotation \n";
    //cout << "XYZ:" << XX << " : " << YY << " : " << ZZ << "\n";

    outXYZ[0] = XX;
    outXYZ[1] = YY;
    outXYZ[2] = ZZ;

}

BmnTOF1Detector * Plane[NDet];

int Process_Gem_Macro(int run, int nEvForRead = 0) {
  //  bmnloadlibs();

    // Init BmnTOF1Detector
    TString name;
    for (Int_t i = 0; i < NDet; i++) {
        name = Form("Plane%d", i);
        cout << endl << "==============Init " << name << "=============" << endl;
        Plane[i] = new BmnTOF1Detector(i, 1, NULL);
        Plane[i]->SetCorrLR("TOF400_LRcorr_RUN7.dat");
        //Plane[i]->SetGeoFile("geofile_full_src.root");
        Plane[i]->SetGeoFile("geofile_full_src_noZshift.root");
    }//*/

    // ==============read input file ===========================================
    TString inFile = Form("bmn_run%d_digi.root", run);

    TFile *f_in1 = new TFile(inFile, "READ");
    TTree *t_in1 = (TTree *) f_in1->Get("cbmsim");

    TClonesArray *ToF400Digits = new TClonesArray("BmnTof1Digit");
    t_in1->SetBranchAddress("TOF400", &ToF400Digits);

    TClonesArray *T0Digits = new TClonesArray("BmnTrigDigit");
    t_in1->SetBranchAddress("BC2", &T0Digits);

    TClonesArray* eventHeader = NULL;
    t_in1->SetBranchAddress("EventHeader", &eventHeader);

    TClonesArray* gem = NULL;
    t_in1->SetBranchAddress("GEM", &gem);


    //============== creat output file =========================================
    TString FileOutName = Form("GemTof_run%04d.root", run);
    TFile *f_out = new TFile(FileOutName.Data(), "RECREATE");
    TTree *treeHitOut = new TTree("cbmsim", "cbmsim");

    vector <double> *XTof = new vector<double>;
    vector <double> *YTof = new vector<double>;
    vector <double> *ZTof = new vector<double>;
    vector <double> *AmpTof = new vector<double>;
    vector <double> *TimeTof = new vector<double>;
    vector <int> *PlaneTof = new vector<int>;
    vector <int> *StripTof = new vector<int>;
    vector <double> *AmpT0 = new vector<double>;
    vector <double> *TimeT0 = new vector<double>;
    vector <double> *XGem = new vector<double>;
    vector <double> *YGem = new vector<double>;
    vector <double> *ZGem = new vector<double>;
    vector <int> *ZoneGem = new vector<int>;
    UInt_t EventHeader_EventID;

    treeHitOut->Branch("EventID", &EventHeader_EventID);
    treeHitOut->Branch("XTof", &XTof);
    treeHitOut->Branch("YTof", &YTof);
    treeHitOut->Branch("ZTof", &ZTof);
    treeHitOut->Branch("AmpTof", &AmpTof);
    treeHitOut->Branch("TimeTof", &TimeTof);
    treeHitOut->Branch("StripTof", &StripTof);
    treeHitOut->Branch("PlaneTof", &PlaneTof);
    treeHitOut->Branch("AmpT0", &AmpT0);
    treeHitOut->Branch("TimeT0", &TimeT0);
    treeHitOut->Branch("XGem", &XGem);
    treeHitOut->Branch("YGem", &YGem);
    treeHitOut->Branch("ZGem", &ZGem);
    treeHitOut->Branch("ZoneGem", &ZoneGem);

    TList *hList = new TList();
    TH1I * hAmp[nstations][2];
    TH2I * hQ[nstations][2];
    TH2I * hGT[3];
    TH1D* Hist_occupancyX[nstations];
    TH1D* Hist_occupancyXp[nstations];
    TH1D* hSpatialXGem[nstations];
    TH1D* hSpatialYGem[nstations];
    TH2D* hSpatialY_XGem[nstations];

    for (Int_t i = 0; i < nstations; i++)
        for (Int_t l = 0; l < 2; l++) {
            hAmp[i][l] = new TH1I(Form("Amp_zone_%d_layer_%d", i, l), Form("Amp_zone_%d_layer_%d", i, l), 250, 0, 500);
            hList->Add(hAmp[i][l]);
            hQ[i][l] = new TH2I(Form("Q_zone_%d_layer_%d", i, l), Form("Q_zone_%d_layer_%d", i, l), 10, 0, 10, 250, 0, 500);
            hList->Add(hQ[i][l]);
        }

    for (Int_t i = 0; i < 3; i++) {
        hGT[i] = new TH2I(Form("GemHit_vs_ToFHit_%d", i), Form("GemHit_vs_ToFHit_%d", i), 10, 0, 10, 10, 0, 10);
        hGT[i]->GetXaxis()->SetTitle("Gem Xits");
        hGT[i]->GetYaxis()->SetTitle("Tof Xits");
        hList->Add(hGT[i]);
    }

    for ( Int_t i = 0; i < nstations; i++){
      Hist_occupancyX[i] = new TH1D(Form("Hist_occupancyX_%d", i), Form("occupancyX_%d", i),  nstrips, 0, nstrips); 
      hList->Add(Hist_occupancyX[i]);
      Hist_occupancyXp[i]= new TH1D(Form("Hist_occupancyXp_%d", i), Form("occupancyXp_%d", i),  nstrips, 0, nstrips); 
      hList->Add(Hist_occupancyXp[i]);
    }
    for ( Int_t i = 0; i < nstations; i++){
      hSpatialXGem[i] = new TH1D(Form("SpatialXGem_%d", i), Form("SpatialXGem_%d", i), 300,-150,150);
      hList->Add(hSpatialXGem[i]);
      hSpatialYGem[i] = new TH1D(Form("SpatialYGem_%d", i), Form("SpatialYGem_%d", i), 100,-50, 50);
      hList->Add(hSpatialYGem[i]);
      hSpatialY_XGem[i] = new TH2D(Form("SpatialY_XGem_%d", i), Form("SpatialY_XGem_%d", i), 300,-150,150, 100,-50, 50);
      hList->Add(hSpatialY_XGem[i]);
    }

    
    //-----Initialisation------
    Int_t NhitsXY[nstations], N_all_hits;
    Int_t NfirstX, NlastX, NfirstXp, NlastXp, ClusterSizeX, ClusterSizeXp;
    Int_t NclustX[nstations], NclustXp[nstations];
    Double_t SumAmpl[nstations];
    Double_t SumAmplXp;

    //--Arrays---

    Double_t Xspatial[nstations][kBig], Yspatial[nstations][kBig], Zspatial[nstations][kBig];

    Double_t** XCoord = new Double_t*[nstations];
    Double_t** XpCoord = new Double_t*[nstations];

    for (Int_t i = 0; i < nstations; i++) {
        XCoord[i] = new Double_t[kBig];
        XpCoord[i] = new Double_t[kBig];
    }

    Double_t DigitsArrayX[nstations][nstrips], DigitsArrayXp[nstations][nstrips];
    Double_t XClust[nstations][kBig], XpClust[nstations][kBig];
    Double_t Ampl_strX[nstrips], Ampl_strXp[nstrips];

    //                                 0        1      2       3
    Double_t shiftX[nstations]   = {-31.863, -32.84, -32.674, -32.01 };
    Double_t shiftY[nstations]   = {-41.858,  0.811, -41.609,  0.898};
    Double_t Zstation[nstations] = { 234.05, 236.15,  234.65, 236.75 };

    //                                0     1     2      3
    Double_t shiftXX[nstations] =  {-.22, -.58, -.46,  -2.45};
    Double_t shiftYY[nstations] =  {-.82, -.42+.18, -1.88, -1.24};

    /* Double_t XYZSAngleTune[nstations][3] = {
        {-0.17731, -0.584688, -0.0382435},
        {-0.143239, -0.654572, 0.137676},
        {-0.276742, 0.684091, -0.0771824},
        {-0.415814, 0.668466, -0.00626}
    }; */

    // /*
      Double_t XYZSAngleTune[nstations][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };//*/

    Double_t XYZSAngleRot[nstations][3] = {
        {0., -30., 0.},
        {0., -30., 0.},
        {0, 30, 0},
        {0, 30, 0}
    };

    Double_t Z0_SRC_target = -645.191;

    //-----------------------------------------------------start---------------------

    Long64_t nEvents = t_in1->GetEntries();
    if (nEvForRead == 0 || nEvForRead > nEvents) nEvForRead = nEvents;
    cout << " Process run= " << run << endl;
    cout << " events will be processed  = " << nEvForRead << endl;

    //read tree and fill array of amplitudes
    for (Int_t iEv = 0; iEv < nEvForRead; ++iEv) {//GEM

        ToF400Digits->Delete();
        T0Digits->Delete();
        eventHeader->Delete();
        gem->Delete();

        XTof->clear();
        YTof->clear();
        ZTof->clear();
        AmpTof->clear();
        TimeTof->clear();
        StripTof->clear();
        PlaneTof->clear();
        AmpT0->clear();
        TimeT0->clear();
        XGem->clear();
        YGem->clear();
        ZGem->clear();
        ZoneGem->clear();

        t_in1->GetEntry(iEv);

        BmnEventHeader* digHeader = (BmnEventHeader*) eventHeader->At(0);
        EventHeader_EventID = digHeader->GetEventId();
	// if (iEv % 1 == 0) {
	if (iEv % 10000 == 0) { 
            cout << "\n\n" << "====== Event number " << iEv << " ============== " << "\n";
            cout << "#EVENTHeader: " << EventHeader_EventID << endl;
        }

        for (Int_t i = 0; i < nstations; i++) {
            for (Int_t ii = 0; ii < nstrips; ii++) {
                DigitsArrayX[i][ii] = 0.;
                DigitsArrayXp[i][ii] = 0.;
            }
        }

        Int_t GemHit = 0;
        Int_t GemHitL = 0;
        Int_t GemHitR = 0;
        //loop over gem digits
        for (Int_t i = 0; i < gem->GetEntriesFast(); ++i) {

            BmnGemStripDigit* digit = (BmnGemStripDigit*) gem->At(i);

            Int_t stat = digit->GetStation();
            Int_t mod = digit->GetModule();
            Int_t layer = digit->GetStripLayer();
            Int_t strip = digit->GetStripNumber();
            Double_t Ampl = digit->GetStripSignal();

            Int_t NumStat = -1;
            if (stat == 0 && mod == 0) NumStat = 1; // X>0; Y>0  
            if (stat == 0 && mod == 1) NumStat = 3; // X<0; Y>0  
            if (stat == 1 && mod == 0) NumStat = 2; // X<0; Y<0  
            if (stat == 1 && mod == 1) NumStat = 0; // X>0; Y<0  
            //zone 0 X>0; Y<0 
            if (layer < 2 && NumStat != -1) hAmp[NumStat][layer]->Fill(Ampl);
            if (NumStat == 0 && Ampl < 40) Ampl = 0.;
            if (NumStat != 0 && Ampl < 80) Ampl = 0.;


            if (NumStat >= 0 && NumStat < 4) {
                // cout<<" stn= "<<stat<<" mod= "<<mod<<" layer= "<<layer<<" strip= "<<strip<<" Ampl= "<<Ampl<<endl;

                if (layer == 0) {		  
		          if (NumStat == 3 && (strip == 396 || strip == 402)) continue; //noisy strips

                    DigitsArrayX[NumStat][strip] = Ampl; Hist_occupancyX[NumStat]->Fill(strip);
                }
                if (layer == 1) {
		          if (NumStat == 1 && strip >= 807 && strip < 825)  continue; //noisy strips
		          if (NumStat == 2 && strip >= 808 && strip < 824) continue; //noisy strips
		          if (NumStat == 3 && (strip == 383 || strip == 390)) continue; //noisy strips
		          if (NumStat == 3 && strip >= 810 && strip < 824) continue; //noisy strips
		  
                    DigitsArrayXp[NumStat][strip] = Ampl; Hist_occupancyXp[NumStat]->Fill(strip);
                }

            }// if ( stat 
        }//end loop over gems digits

        // ---- per ev --- reset arrays
        for (Int_t izone = 0; izone < nstations; izone++) {
            NclustXp[izone] = 0;
            NclustX[izone] = 0;

            for (Int_t ii = 0; ii < kBig; ii++) {
                XClust[izone][ii] = 0.;
                XpClust[izone][ii] = 0.;
                XCoord[izone][ii] = -900.;
                XpCoord[izone][ii] = -900.;

                Xspatial[izone][ii] = -900.;
                Yspatial[izone][ii] = -900.;
                Zspatial[izone][ii] = -900.;

            }
        }


        for (Int_t istr = 0; istr < nstrips; istr++) {
            Ampl_strX[istr] = 0.;
            Ampl_strXp[istr] = 0.;

        }

        //---clustering -----------------------------------------------------------------   
        Double_t CoG;

        //loop over zones
        for (Int_t izone = 0; izone < nstations; izone++) {
            NfirstX = -1;
            NlastX = -1;
            ClusterSizeX = 0;
            SumAmpl[izone] = 0.;
            //----x-----
            for (Int_t istr = 1; istr < nstrips; istr++) {

                // looking for start and end of cluster
                if (NfirstX < 0 && DigitsArrayX[izone][istr] == 0.) continue;
                if (NfirstX < 0 && DigitsArrayX[izone][istr] > 0.) NfirstX = istr;
                if (NfirstX >= 0 && DigitsArrayX[izone][istr + 1] == 0.) NlastX = istr;

                // skip cluster width 1
                if (NlastX - NfirstX == 0) {
                    NfirstX = -1;
                    NlastX = -1;
                    continue;
                } // Misha's lines 

                // skip big clusters
                if (NlastX - NfirstX > 5) {//5
                    NfirstX = -1;
                    NlastX = -1;
                    continue;
                } // Misha's lines

                if (NfirstX >= 0 && NlastX > 0) {
                    ClusterSizeX = NlastX - NfirstX + 1;

                    for (Int_t is = NfirstX; is < NlastX + 1; is++)
                        Ampl_strX[is - NfirstX] = DigitsArrayX[izone][is];

                    //calculate center of cluster
                    CoG = FindClusterCenter(Ampl_strX, ClusterSizeX, SumAmpl[izone]);
                    hQ[izone][0]->Fill(ClusterSizeX, SumAmpl[izone]);
                    XClust[izone][NclustX[izone]] = NfirstX + CoG;
                    NclustX[izone]++;

                    //break in case clusters more then 20
                    if (NclustX[izone] > 20)
                        break;

                    NfirstX = -1;
                    NlastX = -1;
                }// if (NfirstX >= 0 && NlastX > 0){
            }//end loop over X


            //----xp-----
            if (NclustX[izone] == 0) continue;

            NfirstXp = -1;
            NlastXp = -1;
            ClusterSizeXp = 0;
            SumAmplXp = 0.;

            for (Int_t istr = 1; istr < nstrips; istr++) {

                // looking for start and end of cluster
                if (NfirstXp < 0 && DigitsArrayXp[izone][istr] == 0.) continue;
                if (NfirstXp < 0 && DigitsArrayXp[izone][istr] > 0.) NfirstXp = istr;
                if (NfirstXp >= 0 && DigitsArrayXp[izone][istr + 1] == 0.) NlastXp = istr;

                // skip cluster width 1
                if (NlastXp - NfirstXp == 0) {
                    NfirstXp = -1;
                    NlastXp = -1;
                    continue;
                } // Misha's lines 

                // skip big clusters
                if (NlastXp - NfirstXp > 5) {//5
                    NfirstXp = -1;
                    NlastXp = -1;
                    continue;
                }// Misha's lines

                if (NfirstXp >= 0 && NlastXp > 0) {
                    ClusterSizeXp = NlastXp - NfirstXp + 1;

                    for (Int_t is = NfirstXp; is < NlastXp + 1; is++)
                        Ampl_strXp[is - NfirstXp] = DigitsArrayXp[izone][is];

                    //calculate center of cluster
                    CoG = FindClusterCenter(Ampl_strXp, ClusterSizeXp, SumAmplXp);
                    hQ[izone][1]->Fill(ClusterSizeXp, SumAmpl[izone]);
                    XpClust[izone][NclustXp[izone]] = NfirstXp + CoG;
                    NclustXp[izone]++;

                    //break in case clusters more then 20
                    if (NclustXp[izone] > 20)
                        break;

                    NfirstXp = -1;
                    NlastXp = -1;
                }// if (NfirstXp >= 0 && NlastXp > 0)
            }//end loop over Xp

        }// end loop over zone


        // ======coordinate calculation=====================================================
        //move X and Xp coordinates on frame size
        for (Int_t ist = 0; ist < nstations; ist++) {
            NhitsXY[ist] = 0;

            for (Int_t cl = 0; cl < NclustX[ist]; ++cl){
	      XCoord[ist][cl] = deltaX * XClust[ist][cl]; //cm
	      //  cout<<" XCoord["<<ist<<"]["<<cl<<"]= "<<XCoord[ist][cl]<<endl;
	    }
            for (Int_t cl = 0; cl < NclustXp[ist]; ++cl){
	      XpCoord[ist][cl] = deltaXp * XpClust[ist][cl]; //cm 
	      // cout<<" XpCoord["<<ist<<"]["<<cl<<"]= "<<XpCoord[ist][cl]<<endl;
	    }
        }//ist

	//cout<<endl;

        GetXYspatial(NclustX, NclustXp, XCoord, XpCoord, Xspatial, Yspatial, NhitsXY);

        N_all_hits = 0;
        for (Int_t ist = 0; ist < nstations; ist++) {

            if (ist == 0 || ist == 1) GemHitL += NhitsXY[ist];
            if (ist == 2 || ist == 3) GemHitR += NhitsXY[ist];
            GemHit += NhitsXY[ist];

            for (Int_t cl = 0; cl < NhitsXY[ist]; ++cl) {

                // rotate 1Gem and 3Gem around Z on 180 deg
                if (ist == 1 || ist == 3) {
                    Xspatial[ist][cl] = -Xspatial[ist][cl] + 66.;
                    Yspatial[ist][cl] = -Yspatial[ist][cl] + 41.;
                }

                Xspatial[ist][cl] += shiftX[ist] + shiftXX[ist];
                Yspatial[ist][cl] += shiftY[ist] + shiftYY[ist];

                Double_t Out[3] = {-200, -200, -1000};
                Double_t XYZSpatial[3] = {Xspatial[ist][cl], Yspatial[ist][cl], 0};
                // local rotation 
                rotateAxis(XYZSpatial, XYZSAngleTune[ist], Out);
                Out[2] += Zstation[ist];

                //Global rotation
                rotateAxis(Out, XYZSAngleRot[ist], Out);
		            Out[2] +=-645.;// Out[3] +=-645.;
		            hSpatialXGem[ist]->Fill(Out[0]);
		            hSpatialYGem[ist]->Fill(Out[1]);
		            hSpatialY_XGem[ist]->Fill(Out[0],Out[1]);

                XGem->push_back(Out[0]);
                YGem->push_back(Out[1]);
                ZGem->push_back(Out[2]);
                ZoneGem->push_back(ist);

            }
        }

        //----------------------------ToF-------------------------
        Int_t TofHit = 0;
        Int_t TofHitL = 0;
        Int_t TofHitR = 0;
        Int_t CountMod0 = 0, CountMod1 = 0;
        for (Int_t i = 0; i < T0Digits->GetEntriesFast(); i++) {
            BmnTrigDigit* digT0 = (BmnTrigDigit*) T0Digits->At(i);
            if (digT0->GetMod() == 0) CountMod0++;
            if (digT0->GetMod() == 1) CountMod1++;
        }

        if (CountMod0 == 1 && CountMod1 >= 0) // check T0
        {
            //            cout << "Looking for right T0 digit" << "; nDigits = " << T0Digits->GetEntriesFast() << endl;
            BmnTrigDigit* digT0 = NULL;
            for (Int_t i = 0; i < T0Digits->GetEntriesFast(); i++) {
                digT0 = (BmnTrigDigit*) T0Digits->At(i);
                if (digT0->GetMod() == 0) break; // take first T0 digit. needed for ToF calculation.
            }

            if (digT0->GetAmp() >= 19.26 && digT0->GetAmp() <= 22.06) {

                //--------------------------- RPC --------------------------------------------------
                for (Int_t i = 0; i < NDet; i++)
                    Plane[i]->Clear();

                Bool_t FlagHit;

                //               cout << "Read TOF400 Digits" << "; NDigits = " << ToF400Digits->GetEntriesFast() << endl;
                for (Int_t iDig = 0; iDig < ToF400Digits->GetEntriesFast(); ++iDig) {
                    FlagHit = kFALSE;
                    BmnTof1Digit* digTof = (BmnTof1Digit*) ToF400Digits->At(iDig);
                    FlagHit = Plane[digTof->GetPlane()]->SetDigit(digTof);
                }

                AmpT0->push_back(digT0->GetAmp());
                TimeT0->push_back(digT0->GetTime());
                //                cout << "Fiil vectors" << endl;
                TVector3 XYZ;
                XYZ.SetXYZ(0., 0., 0.);
                Double_t ToF = 0;
                FlagHit = kFALSE;
                for (Int_t i = 0; i < NDet; i++) {// loop over TOF400 Planes
                    Int_t nHits = Plane[i] -> FindHits(digT0);
                    if (i < 10) TofHitL += nHits;
                    if (i >= 10) TofHitR += nHits;
                    TofHit += nHits;
                    //cout << "nHits in Plane " << i << " = " << nHits << endl;
                    if (nHits > 0) {
                        for (Int_t s = 0; s < 47; s++) {//loop over Plane's strips
                            XYZ.SetXYZ(0., 0., 0.);
                            ToF = 0;
                            FlagHit = kFALSE;
                            FlagHit = Plane[i]->GetXYZTime(s, &XYZ, &ToF);
                            if (FlagHit == kTRUE) {
 //cout << "Plane = " << i << "; Strip = " << s << "; X=" << XYZ.X() << "; Y=" << XYZ.Y() << "; Z=" << XYZ.Z() << "; Time=" << ToF << endl;
                                XTof->push_back(XYZ.X());
                                YTof->push_back(XYZ.Y());
                                //ZTof->push_back(XYZ.Z());
                                ZTof->push_back(XYZ.Z()-645.);
                                AmpTof->push_back(Plane[i]->GetWidth(s));
				// TimeTof->push_back(Plane[i]->GetTime(s));
				// TimeTof->push_back(ToF);
                                StripTof->push_back(s);
                                PlaneTof->push_back(i);
                            }
                        } // end for (Int_t s = 0; s < 47; s++)
                    }
                } // end for (Int_t i = 0; i < NDet; i++)

            }// end if (digT0->GetAmp() >= 17.3 && digT0->GetAmp() <= 19.2)
        } // end if (CountMod0 == 1 && CountMod1 <= 0)

        //                cout << "Fill output tree" << endl;
        if (XTof->size() > 0 && XGem->size() > 0) {
            treeHitOut->Fill();

            if (GemHitL != 0 && TofHitL != 0)
                hGT[0]->Fill(GemHitL, TofHitL);
            if (GemHitR != 0 && TofHitR != 0)
                hGT[1]->Fill(GemHitR, TofHitR);
            hGT[2]->Fill(GemHit, TofHit);
        }


    }//for (Int_t iEv = 0; iEv < nEvForRead; ++iEv)

    f_in1->Close();

    f_out->cd();
    treeHitOut->Write();
    hList->Write();
    f_out->Close();

    cout << "\n\n" << "=========Finished===========" << "\n";

    return 0;
}//Process

