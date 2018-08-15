#include "BmnTOF1Detector.h"
#include <iomanip>

ClassImp(BmnTOF1Detector)

	BmnTOF1Detector::BmnTOF1Detector() {
	}

//----------------------------------------------------------------------------------------

BmnTOF1Detector::BmnTOF1Detector(Int_t NPlane, Int_t FillHist = 0, TTree *tree = NULL) {
	Clear();
	memset(fCorrTimeShift, 0, sizeof (fCorrTimeShift));
	fNEvents = 0;
	//KillStrip(0);
	//KillStrip(47);
	fFillHist = FillHist;
	fNPlane = NPlane;
	fStripLength = 30; // cm
	fSignalVelosity = 0.06; // 0.06 ns/cm
	fMaxDelta = (fStripLength * 0.5 + 2.0) * fSignalVelosity; // + 20 mm on the strip edge

	for (Int_t i = 0; i < fNStr; i++)
		gSlew[i] = NULL;

	fName = Form("Plane_%d", NPlane);
	TString Name;

	if (fFillHist > 0) {
		fHistListStat = new TList();
		fHistListdt = new TList();

		Name.Clear();
		Name = Form("Hist_HitByCh_%s", fName.Data());
		hHitByCh = new TH1I(Name, Name, fNStr + 1, -0.5, fNStr + 0.5);
		fHistListStat->Add(hHitByCh);

		Name.Clear();
		Name = Form("Hist_HitPerEv_%s", fName.Data());
		hHitPerEv = new TH1I(Name, Name, fNStr + 1, -0.5, fNStr + 0.5);
		fHistListStat->Add(hHitPerEv);

		Name.Clear();
		Name = Form("Hist_HitLR_%s", fName.Data());
		hHitLR = new TH2I(Name, Name, fNStr, 0, fNStr, fNStr, 0, fNStr);
		fHistListStat->Add(hHitLR);

		Name.Clear();
		Name = Form("Hist_XY_%s", fName.Data());
		hXY = new TH2I(Name, Name, 240, -150, 150, 120, -75, 75);
		fHistListStat->Add(hXY);

		hDy_near = new TH1S(Form("hDy_near_%s", fName.Data()), Form("hDy_near_%s", fName.Data()), 400, -20, 20);
		hDtime_near = new TH1S(Form("hDtime_near_%s", fName.Data()), Form("hDtime_near_%s", fName.Data()), 400, -10, 10);
		hDWidth_near = new TH1S(Form("hDWidth_near_%s", fName.Data()), Form("hDWidth_near_%s", fName.Data()), 256, -28., 28.);
		hTempDtimeDy_near = new TH2S(Form("hTempDtimeDy_near_%s", fName.Data()), Form("hTempDtimeDy_near_%s", fName.Data()), 400, -10, 10, 200, -10, 10);
		hDy_acros = new TH1S(Form("hDy_acros_%s", fName.Data()), Form("hDy_acros_%s", fName.Data()), 400, -20, 20);
		hDtime_acros = new TH1S(Form("hDtime_acros_%s", fName.Data()), Form("hDtime_acros_%s", fName.Data()), 400, -10, 10);
		hDWidth_acros = new TH1S(Form("hDWidth_acros_%s", fName.Data()), Form("hDWidth_acros_%s", fName.Data()), 256, -28., 28.);
		hTempDtimeDy_acros = new TH2S(Form("hTempDtimeDy_acros_%s", fName.Data()), Form("hTempDtimeDy_acros_%s", fName.Data()), 400, -10, 10, 200, -10, 10);
		fHistListStat->Add(hDy_near);
		fHistListStat->Add(hDtime_near);
		fHistListStat->Add(hDWidth_near);
		fHistListStat->Add(hTempDtimeDy_near);
		fHistListStat->Add(hDy_acros);
		fHistListStat->Add(hDtime_acros);
		fHistListStat->Add(hDWidth_acros);
		fHistListStat->Add(hTempDtimeDy_acros);

		for (Int_t i = 0; i < fNStr + 1; i++){
			hdT_vs_WidthDet[i] = new TH2I (Form("dt_vs_WidthDet_str_%d_%s", i, fName.Data()), Form("dt_vs_WidthDet_str_%d_%s", i, fName.Data()), 1024, 0, 50, 1024, -4, 20);
			fHistListdt->Add(hdT_vs_WidthDet[i]);
		}
		for (Int_t i = 0; i < fNStr + 1; i++){
			hdT_vs_WidthT0[i] = new TH2I (Form("dt_vs_WidthT0_str_%d_%s", i, fName.Data()), Form("dt_vs_WidthT0_str_%d_%s", i, fName.Data()), 1024, 0, 50, 1024, -4, 20);
			fHistListdt->Add(hdT_vs_WidthT0[i]);
		}
		for (Int_t i = 0; i < fNStr + 1; i++){
			hdT[i] = new TH1I (Form("dt_str_%d_%s", i, fName.Data()), Form("dt_str_%d_%s", i, fName.Data()), 1024, -12, 12);
			fHistListdt->Add(hdT[i]);
		}

	} else {

		hHitByCh = NULL;
		hHitPerEv = NULL;
		hHitLR = NULL;
		hXY = NULL;

		hDy_near = NULL;
		hDtime_near = NULL;
		hDWidth_near = NULL;
		hTempDtimeDy_near = NULL;
		hDy_acros = NULL;
		hDtime_acros = NULL;
		hDWidth_acros = NULL;
		hTempDtimeDy_acros = NULL;

		for (Int_t i = 0; i < fNStr + 1; i++){
			hdT_vs_WidthDet[i] = NULL;
			hdT_vs_WidthT0[i] = NULL;        
			hdT[i] = NULL;
		}

	}

}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::Clear() {
	memset(fStripShift, 0, sizeof(fStripShift));
	memset(fWalkFunc, 0, sizeof(fWalkFunc));
	memset(fTimeL, 0, sizeof (fTimeL));
	memset(fTimeR, 0, sizeof (fTimeR));
	memset(fTimeLtemp, 0, sizeof (fTimeLtemp));
	memset(fTimeRtemp, 0, sizeof (fTimeRtemp));
	memset(fTime, 0, sizeof (fTime));
	memset(fWidthL, 0, sizeof (fWidthL));
	memset(fWidthR, 0, sizeof (fWidthR));
	memset(fWidthLtemp, 0, sizeof (fWidthLtemp));
	memset(fWidthRtemp, 0, sizeof (fWidthRtemp));
	memset(fWidth, 0, sizeof (fWidth));
	memset(fFlagHit, 0, sizeof (fFlagHit));
	memset(fTof, 0, sizeof (fTof));
	memset(fDigitL, 0, sizeof (fDigitL));
	memset(fDigitR, 0, sizeof (fDigitR));
	memset(fHit, 0, sizeof (fHit));
	memset(fKilled, 0, sizeof (fKilled));
	memset(fCorrLR, 0, sizeof (fCorrLR));
	fHit_Per_Ev = 0;

	for (Int_t i = 0; i < fNStr; i++)
		fCrossPoint[i].SetXYZ(0., 0., 0.);
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetDigit(BmnTof1Digit * TofDigit) {
	fStrip = TofDigit->GetStrip();
	if (fStrip < 0 || fStrip > fNStr) return kFALSE;
	//cout << " Plane = " << TofDigit->GetPlane() << "; Strip " << TofDigit->GetStrip() << "; Side " << TofDigit->GetSide() << "; Time " << TofDigit->GetTime() << "; Amp " << TofDigit->GetAmplitude() << endl;
	if (TofDigit->GetSide() == 0 && fFlagHit[fStrip] == kFALSE && fKilled[fStrip] == kFALSE) {
		fTimeLtemp[fStrip] = TofDigit->GetTime() - 2. * fCorrLR[fStrip];
		//cout << "Setting Shift: strip # " << fStrip << " shift val " << CorrLR[fStrip] << " curr timeL " << TofDigit->GetTime() << " shifted timeL " << fTimeLtemp[fStrip] <<  "\n";
		fWidthLtemp[fStrip] = TofDigit->GetAmplitude();
		fDigitL[fStrip]++;
	}
	if (TofDigit->GetSide() == 1 && fFlagHit[fStrip] == kFALSE && fKilled[fStrip] == kFALSE) {
		fTimeRtemp[fStrip] = TofDigit->GetTime();
		//cout << "Setting Shift: strip # " << fStrip << " shift val " << CorrLR[fStrip] << " curr timeR " << TofDigit->GetTime() << " shifted timeR " << fTimeRtemp[fStrip] <<  "\n";
		fWidthRtemp[fStrip] = TofDigit->GetAmplitude();
		fDigitR[fStrip]++;
	}
	if (
			fTimeRtemp[fStrip] != 0 && fTimeLtemp[fStrip] != 0
			&& TMath::Abs((fTimeLtemp[fStrip] - fTimeRtemp[fStrip]) * 0.5) <= fMaxDelta // cat for length of strip  
			//        && TMath::Abs((fWidthLtemp[fStrip] - fWidthRtemp[fStrip]) * 0.5) <= 1.5 // cat for Amplitude correlation
			//&& fFlagHit[fStrip] == kFALSE
	   )
		if (fFlagHit[fStrip] == kFALSE) {
			//cout << "Before set variable: " << fTimeL[fStrip] << " " << fTimeR[fStrip] << "\n";
			fTimeL[fStrip] = fTimeLtemp[fStrip];
			fTimeR[fStrip] = fTimeRtemp[fStrip];
			fWidthL[fStrip] = fWidthLtemp[fStrip];
			fWidthR[fStrip] = fWidthRtemp[fStrip];
			fFlagHit[fStrip] = kTRUE;
			fHit[fStrip]++;
			//cout << "After set variable: " << fTimeL[fStrip] << " " << fTimeR[fStrip] << "\n";
		} else
			fHit[fStrip]++;

	return fFlagHit[fStrip];
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::KillStrip(Int_t NumberOfStrip) {
	fKilled[NumberOfStrip] = kTRUE;
}

//----------------------------------------------------------------------------------------

Int_t BmnTOF1Detector::FindHits(BmnTrigDigit *T0) {
	return fHit_Per_Ev;
}

//----------------------------------------------------------------------------------------

Int_t BmnTOF1Detector::FindHits(BmnTrigDigit *T0, TClonesArray *TofHit) {
	fT0 = T0;
	fNEvents++;
	Bool_t flag;
	for (Int_t i = 0; i < fNStr; i++)
		if (
				fWidthL[i] != 0 && fWidthR[i] != 0
				//&& fFlagHit[fStrip] == kTRUE
		   ) {
			fHit_Per_Ev++;
			fWidth[i] = fWidthL[i] + fWidthR[i];
			fTime[i] = (fTimeL[i] + fTimeR[i]) * 0.5;
			flag = GetCrossPoint(i);
			if (fT0 == NULL) continue;
			if (fT0 != NULL) fTof[i] = CalculateDt(i);
			if (fFillHist > 0) {
				hdT_vs_WidthDet[i] -> Fill(fWidth[i], fTof[i]);
				hdT_vs_WidthT0[i] -> Fill(fT0->GetAmp(), fTof[i]);
				hdT[i] -> Fill(fTof[i]);
				hdT_vs_WidthDet[fNStr] -> Fill(fWidth[i], fTof[i]);
				hdT_vs_WidthT0[fNStr] -> Fill(fT0->GetAmp(), fTof[i]);
				hdT[fNStr] -> Fill(fTof[i]);
				if (i > 3) {
					if (fFlagHit[i - 1] == kTRUE) {
						hDy_near->Fill(fCrossPoint[i].Y() - fCrossPoint[i - 1].Y());
						hDtime_near->Fill(fTof[i] - fTof[i - 1]);
						hDWidth_near->Fill(fWidth[i] - fWidth[i - 1]);
						hTempDtimeDy_near->Fill(fTof[i] - fTof[i - 1], fCrossPoint[i].Y() - fCrossPoint[i - 1].Y());
					}
					if (fFlagHit[i - 2] == kTRUE) {
						hDy_acros->Fill(fCrossPoint[i].Y() - fCrossPoint[i - 2].Y());
						hDtime_acros->Fill(fTof[i] - fTof[i - 2]);
						hDWidth_acros->Fill(fWidth[i] - fWidth[i - 2]);
						hTempDtimeDy_acros->Fill(fTof[i] - fTof[i - 2], fCrossPoint[i].Y() - fCrossPoint[i - 2].Y());
					}
				}
			}
			TString Name = TofHit->GetClass()->GetName();
			if (Name == "BmnTofHit"){
				//            cout << " Fill BmnTofHit" << endl;
				AddHit(i, TofHit);
			}
			else if (Name == "BmnTOF1Conteiner"){
				//           cout << " Fill BmnTOF1Conteiner" << endl;
				AddConteiner(i, TofHit);
			}
		}

	if (fFillHist > 0)
		FillHist();

	return fHit_Per_Ev;
}

//------------------------------------------------------------------------------------------------------------------------

void BmnTOF1Detector::AddHit(Int_t Str, TClonesArray *TofHit) {

	fVectorTemp.SetXYZ(0.5, 0.36, 1.); // error for point dx = 0.5 cm; dy = 1.25/SQRT(12) = 0.36 cm; dz = 1(?)cm
	Int_t UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, Str + 1); // strip [0,47] -> [1, 48]
	BmnTofHit *pHit = new ((*TofHit)[TofHit->GetEntriesFast()]) BmnTofHit(UID, fCrossPoint[Str], fVectorTemp, -1);

	pHit->SetTimeStamp(fTof[Str]);
	pHit->AddLink(FairLink(0x1, -1));
	pHit->AddLink(FairLink(0x2, -1));
	pHit->AddLink(FairLink(0x4, UID));
}

//------------------------------------------------------------------------------------------------------------------------

void BmnTOF1Detector::AddConteiner(Int_t Str, TClonesArray *TofHit) {

	new((*TofHit)[TofHit->GetEntriesFast()]) BmnTOF1Conteiner(fNPlane, Str, fTimeL[Str], fTimeR[Str], fTime[Str], fWidthL[Str], fWidthR[Str], fWidth[Str], fCrossPoint[Str].x(), fCrossPoint[Str].y(), fCrossPoint[Str].z(), fT0->GetTime(), fT0->GetAmp());

}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::FillHist() {
	hHitPerEv->Fill(fHit_Per_Ev);
	for (Int_t i = 0; i < fNStr; i++)
		for (Int_t j = 0; j < fNStr; j++) {
			if (fWidthL[i] != 0 && fWidthR[j] != 0) {
				hHitLR->Fill(i, j);
				if (i == j) {
					hHitByCh->Fill(i);
					hXY->Fill(fCrossPoint[i].x(), fCrossPoint[i].y());
				}
			}
		}

}

//----------------------------------------------------------------------------------------

Double_t BmnTOF1Detector::CalculateDt(Int_t Str = 0) {
	Double_t dt = 0;
	Double_t T0Amp;
	dt = fTime[Str] - fT0->GetTime() - 270.; // RUN7 SRC
	//dt = fTime[Str] - fT0->GetTime();

	T0Amp = fT0->GetAmp();

	/*dt = dt - (-6.794 + 2.11 * T0Amp
	  - 0.1706 * T0Amp * T0Amp
	  + 0.004286 * T0Amp * T0Amp * T0Amp); //RUN7 SRC By BC4-BC2 distribution (from run 3463) not work!!!*/

	/*dt = dt - (1.947 - 0.5363 * T0Amp
	  + 0.03428 * T0Amp * T0Amp
	  - 0.0005853 * T0Amp * T0Amp * T0Amp);// RUN6 */

	/*dt = dt - (18.7191 - 3.17596 * T0Amp
	  + 0.172444 * T0Amp * T0Amp
	  - 0.00292435 * T0Amp * T0Amp * T0Amp) // RUN7 SRC from dt_vs_WidthT0 Plane7 str2. for line -3.15016 + 0.195976*w 
	  - (0.7491 - 0.0219 * T0Amp); // iteration 2*/

	/*dt = dt - (1.564 + 0.1065 * T0Amp
	  + 0.0 * T0Amp * T0Amp
	  - 0.0 * T0Amp * T0Amp * T0Amp);//RUN7 BM@N preliminarily*/

	if (gSlew[Str] != NULL) dt = dt - gSlew[Str]->Eval(fWidth[Str]) + fCorrTimeShift[Str]; // CorrTimeShift is ToF for Gamma
	//cout << dt << endl;
	return dt;
}

//----------------------------------------------------------------------------------------

TList* BmnTOF1Detector::GetList(Int_t n = 0) {
	if (fFillHist > 0) {
		if (n == 0) return fHistListStat;
	} else return NULL;
}

//----------------------------------------------------------------------------------------

TString BmnTOF1Detector::GetName() {
	return fName;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetCorrLR(Double_t* Mass) {
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetCorrLR( string pathToFile ) {
	//dir += "/input/TOF400_LRcorr_RUN7_SRC.dat"; 

	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	cout << "TOF400-Setup: Attempting to open LR corr file: " << dir << "\n"; 

	ifstream f_corr;
	f_corr.open(dir); 
	char line[256]; 
	f_corr.getline(line, 256); 
	f_corr.getline(line, 256); 
	int Pl, St; 
	double Temp; 
	if (f_corr.is_open() == true){ 
		while (!f_corr.eof()) { 
			f_corr >> Pl >> St >> Temp; 
			if (Pl == fNPlane){
				fCorrLR[St] = Temp; 
			}
		} 
		cout << "\tLoaded LR corr file\n"; 
	} 
	else{ 
		cout << "\tFailed to find LR corr file, setting all corrections to 0...\n"; 
	}  
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetStripShift( string pathToFile ) {
	//dir += "/input/TOF400_StripShift_RUN7_SRC.dat";
	
	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	cout << "TOF400-Setup: Attempting to open strip shift file: " << dir << "\n"; 
	
	ifstream f_corr;
	f_corr.open(dir); 
	char line[256]; 
	f_corr.getline(line, 256); 
	f_corr.getline(line, 256); 
	int Pl, St; 
	double Temp2, Temp3; 

	if (f_corr.is_open() == true){ 
		while (!f_corr.eof()) { 
			f_corr >> Pl >> St >> Temp2 >> Temp3;
			if (Pl == fNPlane){
				fStripShift[St] = Temp2;
				if( fStripShift[St] == -1)
					fKilled[St] = true;
			}
		}
		cout << "\tLoaded strip shift file\n";
	}
	else{
		cout << "\tFailed to find strip shift file, setting all corrections to 0...\n";
	} 
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetWalkFunc( string pathToFile ) {
        //dir += "/input/TOF400_TimeWalk_RUN7_SRC.dat";
	
	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	cout << "TOF400-Setup: Attempting to open strip walk-function file: " << dir << "\n"; 
	
	ifstream f_corr;
	f_corr.open(dir); 
	char line[256]; 
	f_corr.getline(line, 256); 
	f_corr.getline(line, 256); 
	int Pl, St;
        double tmp_Pt, tmp_Sh, tmp_p0, tmp_p1, tmp_p2, tmp_p0e, tmp_p1e, tmp_p2e;
        if (f_corr.is_open() == true){
                while (!f_corr.eof()) {
                        f_corr >> Pl >> St >> tmp_Pt \
                                >> tmp_Sh >> tmp_p0 >> tmp_p1 >> tmp_p2 \
                                >> tmp_p0e >> tmp_p1e >> tmp_p2e;
			if( Pl == fNPlane ){
				fWalkFunc[St][0] = tmp_Sh;
				fWalkFunc[St][1] = tmp_p0;
				fWalkFunc[St][2] = tmp_p1;
				fWalkFunc[St][3] = tmp_p2;
				if( tmp_Sh == -1)
					fKilled[St] = true;
			}
                }
                cout << "\tLoaded time walk file\n";
        }
        else{
                cout << "\tFailed to find time walk file, setting all corrections to 0...\n";
        } 

}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::TestPrint( int strip ){
	cout << fCorrLR[strip] << " " << fStripShift[strip] << " " << fWalkFunc[strip][0] << " " << fWalkFunc[strip][1] << " "
		<< fWalkFunc[strip][2] << " " << fWalkFunc[strip][3] << "\n";
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetCorrTimeShift( string pathToFile ){
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::GetCrossPoint(Int_t NStrip = 0) {

	fVectorTemp.SetXYZ(0., 0., 0.);
	if (TMath::Abs((fTimeL[NStrip] - fTimeR[NStrip]) * 0.5) >= fMaxDelta)
		return kFALSE; // estimated position out the strip edge.
	double dL = (fTimeL[NStrip] - fTimeR[NStrip]) * 0.5 / fSignalVelosity;
	fVectorTemp(0) = 0;
	fVectorTemp(1) = dL;
	fVectorTemp(2) = 0; //TMP ALIGMENT CORRECTIONS
	fCrossPoint[NStrip] = fCenterStrip[NStrip] + fVectorTemp;
	//    cout << "Z = " << fCrossPoint[NStrip].Z() << endl;
	return kTRUE;
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetGeoFile( string pathToFile ) {
	
	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	cout << "TOF400-Setup: Attempting to open geometry file: " << dir << "\n"; 

	// get gGeoManager from ROOT file 
	TFile* geoFile = new TFile( dir.c_str() , "READ");
	if( geoFile->IsZombie() ){
		cout << "\tError: could not open root geometry file: " + dir << endl;
		return;
	}
	TList* keyList = geoFile->GetListOfKeys();
	TIter next(keyList);
	TKey* key = (TKey*) next();
	TString className(key->GetClassName());
	if (className.BeginsWith("TGeoManager"))
		key->ReadObj();
	else {
		cout << "\tError: TGeoManager isn't top element in geometry file " + dir << endl;
		return;
	}

	// Found TGeoManager, so load geometry from file:
	BmnTof1GeoUtils *pGeoUtils = new BmnTof1GeoUtils(); 	// empty constructor
	pGeoUtils->ParseTGeoManager(false, NULL, true);		// read in TGeo file
								// setting global coords
								// for each strip
	Int_t UID;
	for (Int_t i = 0; i < fNStr; i++) {
		UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, i + 1); // strip [0,47] -> [1, 48] -- gives unique ID for a given strip in the entire system
		const LStrip1 *pStrip = pGeoUtils->FindStrip(UID);	 // finds this strip based on it's unique ID for fastest lookup
		fCenterStrip[i] = pStrip->center;			 // saves the strip center in global coordinates, which is center of rectangle of strip
	}
	geoFile->Close();
	pGeoUtils->~BmnTof1GeoUtils();

}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetGeo(BmnTof1GeoUtils *pGeoUtils) {
	Int_t UID;
	for (Int_t i = 0; i < fNStr; i++) {
		UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, i + 1); // strip [0,47] -> [1, 48]
		const LStrip1 *pStrip = pGeoUtils->FindStrip(UID);
		fCenterStrip[i] = pStrip->center;
		//        if (fNPlane >= 5) fCentrStrip[i].SetX(fCentrStrip[i].X()+5.5); // for field run only
		//        else fCentrStrip[i].SetX(fCentrStrip[i].X()+2.5);
	}
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::GetXYZTime(Int_t Str, TVector3 *XYZ, Double_t *ToF) {

	if (fTof[Str] == 0) return kFALSE;
	if (NULL == XYZ && NULL == ToF) return kFALSE;
	XYZ->SetXYZ(fCrossPoint[Str].x(), fCrossPoint[Str].y(), fCrossPoint[Str].z());
	*ToF = fTof[Str];
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Double_t BmnTOF1Detector::GetWidth(Int_t Str = 1) {
	return fWidth[Str];
}

//----------------------------------------------------------------------------------------

Double_t BmnTOF1Detector::GetTime(Int_t Str = 1) {
	return fTime[Str];
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SaveHistToFile(TString NameFile) {

	if (fFillHist > 0) {
		TFile *fileout = new TFile(NameFile.Data(), "UPDATE");
		Int_t ResWrite;

		TDirectory *Dir;
		TString Name;
		Name = Form("Tof400_%s", fName.Data());
		Dir = fileout->mkdir(Name.Data());
		Dir->cd();
		//Dir->pwd();

		TDirectory * DirStat;
		DirStat = Dir->mkdir("Statistic");
		DirStat -> cd();
		//DirStat->pwd();
		ResWrite = 0;
		ResWrite = fHistListStat->Write();
		//cout << "Res write = " << ResWrite << endl;

		TDirectory * DirdT;
		DirdT = Dir->mkdir("dt");
		DirdT -> cd();
		//DirStat->pwd();
		ResWrite = 0;
		ResWrite = fHistListdt->Write();
		//cout << "Res write = " << ResWrite << endl;

		fileout->Close(); //*/
		return kTRUE; //*/

	} else
		return kFALSE;

}

//----------------------------------------------------------------------------------------

