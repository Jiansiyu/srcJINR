#include "BmnTOF1Detector.h"
#include <iomanip>

ClassImp(BmnTOF1Detector)
BmnTOF1Detector::BmnTOF1Detector() {
}

//----------------------------------------------------------------------------------------

BmnTOF1Detector::BmnTOF1Detector(Int_t NPlane, Int_t FillHist = 0, TTree *tree = NULL) {
	Clear();

	fNPlane = NPlane;
	fStripLength = 30; // cm
	fSignalVelosity = 0.06; // 0.06 ns/cm
	fMaxDelta = (fStripLength * 0.5 + 2.0) * fSignalVelosity; // + 20 mm on the strip edge

}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::Clear() {
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
	memset(fCorrTimeShift, 0, sizeof (fCorrTimeShift));
	fHit_Per_Ev = 0;

	for (Int_t i = 0; i < fNStr; i++)
		fCrossPoint[i].SetXYZ(0., 0., 0.);
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::KillStrip(Int_t NumberOfStrip) {
	fKilled[NumberOfStrip] = kTRUE;
}

//----------------------------------------------------------------------------------------

Int_t BmnTOF1Detector::FindHits(BmnTrigDigit *T0, TClonesArray *TofHit) {
	fT0 = T0;
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
	for (Int_t i = 0; i < 48; i++)
		fCorrLR[i] = Mass[i];
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetCorrLR(TString NameFile) {
	char line[256];
	Int_t Pl, St;
	Double_t Temp;
	ifstream f_corr;
	TString dir = Form("%s%s%s", getenv("VMCWORKDIR"), "/input/", NameFile.Data());
	f_corr.open(dir);
	f_corr.getline(line, 256);
	f_corr.getline(line, 256);
	if (f_corr.is_open() == kTRUE) {
		while (!f_corr.eof()) {
			f_corr >> Pl >> St >> Temp;
			if (Pl == fNPlane) {
				fCorrLR[St] = Temp;
				f_corr >> Temp;
				// If diff between my shift and old shift is greater than the actual cable, throw 
				// strip out
				// if (TMath::Abs(Temp - fCorrLR[St]) > 2.) fCorrLR[St] = -11.9766;
				// cout << Pl << " " << St << " " << CorrLR[St] << " " << Temp << "\n";
			} else
				f_corr >> Temp;
		}
	} else {
		cout << "File " << NameFile.Data() << " for LR correction is not found" << endl;
		cout << "Check " << dir.Data() << " folder for file" << endl;
		return kFALSE;
	}
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetCorrSlewing(TString NameFile) {
	TString PathToFile = Form("%s%s%s", getenv("VMCWORKDIR"), "/input/", NameFile.Data());
	TString name, dirname;
	TDirectory *Dir;
	TFile *f_corr = new TFile(PathToFile.Data(), "READ");
	if (f_corr->IsOpen()) {
		dirname = Form("Plane_%d", fNPlane);
		f_corr->cd(dirname.Data());
		Dir = f_corr-> CurrentDirectory();
		for (Int_t i = 0; i < fNStr; i++) {
			name = Form("Graph_TA_Plane%d_str%d", fNPlane, i);
			gSlew[i] = (TGraphErrors*) Dir->Get(name.Data());
		}
	} else {
		cout << "File " << NameFile.Data() << " for Slewing correction is not found" << endl;
		cout << "Check " << PathToFile.Data() << " folder for file" << endl;
		return kFALSE;
	}
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetCorrTimeShift(TString NameFile) {
	char line[256];
	Int_t Pl, St;
	Double_t Temp;
	ifstream f_corr;
	TString dir = Form("%s%s%s", getenv("VMCWORKDIR"), "/input/", NameFile.Data());
	f_corr.open(dir);
	f_corr.getline(line, 256);
	f_corr.getline(line, 256);
	if (f_corr.is_open() == kTRUE) {
		while (!f_corr.eof()) {
			f_corr >> Pl >> St >> Temp;
			if (Pl == fNPlane) {
				fCorrTimeShift[St] = Temp;
				//cout << Pl << " " << St << " " << CorrTimeShift[St] << "\n";
			}
		}
	} else {
		cout << "File " << NameFile.Data() << " for TimeShift correction is not found" << endl;
		cout << "Check " << dir.Data() << " folder for file" << endl;
		return kFALSE;
	}
	return kTRUE;
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
	fCrossPoint[NStrip] = fCentrStrip[NStrip] + fVectorTemp;
	//    cout << "Z = " << fCrossPoint[NStrip].Z() << endl;
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetGeoFile(TString NameFile) {

	// get gGeoManager from ROOT file 
	TString PathToFile = Form("%s%s%s", getenv("VMCWORKDIR"), "/macro/run/", NameFile.Data());
	TFile* geoFile = new TFile(PathToFile, "READ");
	if (!geoFile->IsOpen()) {
		cout << "Error: could not open ROOT file with geometry: " + NameFile << endl;
		return kFALSE;
	}
	TList* keyList = geoFile->GetListOfKeys();
	TIter next(keyList);
	TKey* key = (TKey*) next();
	TString className(key->GetClassName());
	if (className.BeginsWith("TGeoManager"))
		key->ReadObj();
	else {
		cout << "Error: TGeoManager isn't top element in geometry file " + NameFile << endl;
		return kFALSE;
	}

	BmnTof1GeoUtils *pGeoUtils = new BmnTof1GeoUtils();
	pGeoUtils->ParseTGeoManager(false, NULL, true);

	Int_t UID;
	for (Int_t i = 0; i < fNStr; i++) {
		UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, i + 1); // strip [0,47] -> [1, 48]
		const LStrip1 *pStrip = pGeoUtils->FindStrip(UID);
		fCentrStrip[i] = pStrip->center;
		//        cout << "Strip = " << i << "; Centr XYZ = " << fCentrStrip[i].x() << "  " << fCentrStrip[i].y() << "  " << fCentrStrip[i].z() << endl;
		//        if (fNPlane >= 5) fCentrStrip[i].SetX(fCentrStrip[i].X() + 5.5); // for field run only
		//        else fCentrStrip[i].SetX(fCentrStrip[i].X() + 2.5);
	}
	geoFile->Close();
	pGeoUtils->~BmnTof1GeoUtils(); // deconstruct 
	return kTRUE;
}

//----------------------------------------------------------------------------------------

Bool_t BmnTOF1Detector::SetGeo(BmnTof1GeoUtils *pGeoUtils) {
	Int_t UID;
	for (Int_t i = 0; i < fNStr; i++) {
		UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, i + 1); // strip [0,47] -> [1, 48]
		const LStrip1 *pStrip = pGeoUtils->FindStrip(UID);
		fCentrStrip[i] = pStrip->center;
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

