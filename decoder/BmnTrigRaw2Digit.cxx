#include "BmnTrigRaw2Digit.h"
#include <climits>
#include <sstream>

const UShort_t kNCHANNELS = 8; // number of channels in one HPTDC

BmnTrigRaw2Digit::BmnTrigRaw2Digit( 	TString TrigPlaceMapFile, 
		TString TrigDetMapFile	,
		TString TrigINLTQDC1File,
		TString TrigINLTQDC2File,
		TString TrigINLTDC1File	,
		TString TrigINLTDC2File	,
		TTree *digiTree ){

	readMap(TrigDetMapFile);

	for (BmnTrigMapping &record : fMap) {
		TString detName = record.name;
		TString clsName = (detName.Contains("TQDC")) ? BmnTrigWaveDigit::Class_Name() : BmnTrigDigit::Class_Name();
		TBranch* br = digiTree->GetBranch(detName.Data());
		if (!br) {
			TClonesArray *ar = new TClonesArray(clsName.Data());
			ar->SetName(detName.Data());
			digiTree->Branch(detName.Data(), &ar);
			trigArrays.push_back(ar);
			record.branchRef = ar;
		} else // in case we have same names, but I don't have that for our trig map
		{
			for (auto tca : trigArrays)
				if (TString(tca->GetName()) == detName) {
					record.branchRef = tca;
					break;
				}
		}
	}
	readINLCorrections(TrigINLTDC1File , 72, 0);
	readINLCorrections(TrigINLTDC2File , 32, 0);
	readINLCorrections(TrigINLTQDC1File, 16, 1);
	readINLCorrections(TrigINLTQDC2File, 16, 2);

}

BmnStatus BmnTrigRaw2Digit::readMap(TString mappingFile) { // in mapping TQDC channels must appear earlier than TDC
	
	fMapFileName = TString(getenv("VMCWORKDIR")) + TString("/input/") + mappingFile;
	printf("Reading Triggers mapping file %s...\n", fMapFileName.Data());
	//========== read mapping file            ==========//
	fMapFile.open((fMapFileName).Data());
	if (!fMapFile.is_open()) {
		cout << "Error opening map-file (" << fMapFileName << ")!" << endl;
	}

	TString dummy;
	TString name;
	UInt_t ser;
	Short_t slot, ch, mod;

	// first two lines are empty for comments
	fMapFile >> dummy >> dummy >> dummy >> dummy >> dummy;
	fMapFile >> dummy;
	while (!fMapFile.eof()) {
		fMapFile >> name >> mod >> hex >> ser >> dec >> slot >> ch;
		if (!fMapFile.good()) break;
		BmnTrigMapping record;
		record.branchRef = NULL;
		record.name = name;
		record.serial = ser;
		record.module = mod;
		record.slot = slot;
		record.channel = ch;
		fMap.push_back(record);
	}
	fMapFile.close();
}

BmnStatus BmnTrigRaw2Digit::readINLCorrections(TString INLFile, int length, int iden) {

	for (int i = 0; i < length; i++)
		for (int j = 0; j < 1024; j++){
			if( length == 72) fINLTable72[i][j] = 0.;
			else if(length==32) fINLTable32[i][j] = 0.;
			else if( (length == 16) && (iden==1)) fINLTable16_1[i][j] = 0.;
			else if( (length == 16) && (iden==2)) fINLTable16_2[i][j] = 0.;
		}

	INLFile = TString(getenv("VMCWORKDIR")) + TString("/input/") + INLFile;
	fstream ff(INLFile, std::fstream::in);
	if(ff.fail()) {cerr << "Failed to open " << INLFile << endl; return kBMNSUCCESS;}

	//The format of the header seems to be [TDC-THESERIAL-inl_corr]
	ff.ignore(10, '-');
	int TDCSerial;
	ff >> std::hex >> TDCSerial >> std::dec;
	ff.ignore(1000, '\n');

	unsigned int chan_id = 0;
	unsigned int lines_num = 0;

	while(!ff.eof()) {
		string line; char dummy;

		std::getline(ff, line, '\n');
		if(ff.eof()) {break;}
		if(line == "") {continue;}
		istringstream ss(line);
		
		ss >> chan_id >> dummy;

		if(dummy != '='){
			cerr << "Wrong INL file format." << endl; ff.close(); 
			return kBMNERROR;
		}
		if(chan_id > length-1){
			cerr << "Wrong channel in in the INL file" << endl; ff.close(); 
			return kBMNERROR;
		}

		unsigned int i_bin = 0;
		while(ss.tellg() != -1) {
			if(i_bin > 1024) {
				cerr << "INL File contains too many bins in channel.\n";
				ff.close(); 
				return kBMNERROR;
			}
			if(ss.peek()==','){
				ss.ignore();
			}
			if( length == 72){ 
				ss >> fINLTable72[chan_id][i_bin]; 
				i_bin++;
			}
			else if(length == 32){
				ss>>fINLTable32[chan_id][i_bin]; 
				i_bin++;
			}
			else if( (length == 16) && (iden==1)){
				ss >> fINLTable16_1[chan_id][i_bin];
				i_bin++;
			}
			else if( (length == 16) && (iden==2)){
				ss >> fINLTable16_2[chan_id][i_bin];
				i_bin++;
			}
		}
		if(i_bin != 1024) {
			cout << "Warning: wrong number of bins in the INL file for channel " << chan_id << " (" << i_bin << ")" << endl;
		}
		lines_num++;
	}

	if( !((lines_num == length-1) || (lines_num == length) ) ) {
		cout << "Warning: wrong number of lines in the INL file (" << lines_num << endl;
	}


	return kBMNSUCCESS;
}

BmnStatus BmnTrigRaw2Digit::FillEvent(TClonesArray *tdc, TClonesArray *adc) {
	std::vector<double> times;
	std::vector<double> diff;
	for (Int_t iMap = 0; iMap < fMap.size(); ++iMap) {
		BmnTrigMapping tM = fMap[iMap];
		Short_t iMod = tM.module;
		TClonesArray *trigAr = tM.branchRef;

		// Matching of ADC/TDC based on the fact that (TDC_i - TDC_j) > 296ns correspond to different ADC
		// and (TDC_i - (ADC_i - Trig_i) ) ~ 0
		// 
		// However, it is possible that not all ADC's will have a TDC
		
		for (Int_t iAdc = 0; iAdc < adc->GetEntriesFast(); iAdc++) {
			times.clear();
			diff.clear();
			BmnTQDCADCDigit *adcDig = (BmnTQDCADCDigit*) adc->At(iAdc);
			if (adcDig->GetSerial() != tM.serial || adcDig->GetSlot() != tM.slot) continue;
			if (adcDig->GetChannel() != tM.channel) continue;
			Double_t adcTimestamp = adcDig->GetAdcTimestamp() * ADC_CLOCK_TQDC16VS;
			Double_t trgTimestamp = adcDig->GetTrigTimestamp() * ADC_CLOCK_TQDC16VS;

			for (Int_t iTdc = 0; iTdc < tdc->GetEntriesFast(); ++iTdc) {
				BmnTDCDigit* tdcDig = (BmnTDCDigit*) tdc->At(iTdc);
				if (tdcDig->GetSerial() != tM.serial || tdcDig->GetSlot() != tM.slot) continue;
				if (tdcDig->GetChannel() != tM.channel) continue;
				
				// Now of all the TDCs, find the TDC for which (TDC - trgTime + adcTime) ~ 0 within
				// some tolerance. The closest to 0 is our nominal TDC. And any nearby TDCs between 296ns
				// belong to the same ADC			
				double time = - 666.;
				if( tM.slot == 13)
					time = (tdcDig->GetValue() + fINLTable16_1[tM.channel][tdcDig->GetValue() % 1024]) * TDC_CLOCK / 1024;
				else if( tM.slot == 14)
					time = (tdcDig->GetValue() + fINLTable16_2[tM.channel][tdcDig->GetValue() % 1024]) * TDC_CLOCK / 1024;
				else{ printf("CRIT ERROR\n"); }				

				Double_t tdcTimestamp = tdcDig->GetTimestamp() * TDC_CLOCK;
				diff.push_back(fabs(time - (adcTimestamp - trgTimestamp)));
				times.push_back(time);

			}
			
			if( diff.size() > 0){
				std::vector<double>::iterator result = std::min_element( std::begin(diff) , std::end(diff) );
				int idx = std::distance(std::begin(diff) , result);
		
				// Found the match, so let's save that as and ADC and corresponding TDC
				double matchTime = times.at(idx);
				double minUsed = diff.at(idx);
				if( minUsed  > 300) matchTime = -666.;  // no TDC for this corresponding ADC
			
				new ((*trigAr)[trigAr->GetEntriesFast()]) BmnTrigWaveDigit(
						iMod,
						adcDig->GetShortValue(),
						adcDig->GetNSamples(),
						trgTimestamp,
						adcTimestamp,
						matchTime);
			}
			else{
				double matchTime = -666.;  // no TDC for this corresponding ADC
			
				new ((*trigAr)[trigAr->GetEntriesFast()]) BmnTrigWaveDigit(
						iMod,
						adcDig->GetShortValue(),
						adcDig->GetNSamples(),
						trgTimestamp,
						adcTimestamp,
						matchTime);
			
			}
		}

	}
	return kBMNSUCCESS;
}


// This is the easier of the two, just fill TDC information for trigger
BmnStatus BmnTrigRaw2Digit::FillEvent(TClonesArray *tdc) {
	for (Int_t iMap = 0; iMap < fMap.size(); ++iMap) {
		BmnTrigMapping tM = fMap[iMap];
		Short_t iMod = tM.module;
		TClonesArray *trigAr = tM.branchRef;

		for (Int_t iTdc = 0; iTdc < tdc->GetEntriesFast(); ++iTdc) {
			BmnTDCDigit* tdcDig1 = (BmnTDCDigit*) tdc->At(iTdc);
			
			if (tdcDig1->GetSerial() != tM.serial || tdcDig1->GetSlot() != tM.slot) continue;
			if (!tdcDig1->GetLeading()) continue; // use only leading digits
			
			UShort_t rChannel1 = tdcDig1->GetHptdcId() * kNCHANNELS + tdcDig1->GetChannel();
			if (rChannel1 != tM.channel) continue;
			BmnTDCDigit* nearestDig = NULL;
			UInt_t nearestTime = UINT_MAX;
			for (Int_t jTdc = 0; jTdc < tdc->GetEntriesFast(); ++jTdc) {
				if (iTdc == jTdc) continue;
				BmnTDCDigit* tdcDig2 = (BmnTDCDigit*) tdc->At(jTdc);
				if (tdcDig2->GetSerial() != tM.serial || tdcDig2->GetSlot() != tM.slot) continue;
				if (tdcDig2->GetLeading()) continue; // use only trailing digits as a pair to leading one
				UShort_t rChannel2 = tdcDig2->GetHptdcId() * kNCHANNELS + tdcDig2->GetChannel();
				if (rChannel1 != rChannel2) continue; // we need the same hptdc & channel to create pair
				Int_t dTime = tdcDig2->GetValue() - tdcDig1->GetValue();
				if (dTime < 0) continue; //time should be positive
				if (dTime < nearestTime) {
					nearestTime = dTime;
					nearestDig = tdcDig2;
				}
			}
			if (nearestDig != NULL) {
				Double_t tL, tT;
				// Now check slot ID to find out which INL correction table to use
				if (tM.slot == 10){
					tL = (tdcDig1->GetValue() + fINLTable72[rChannel1][tdcDig1->GetValue() % 1024]) * 24.0 / 1024;
					tT = (nearestDig->GetValue() + fINLTable72[rChannel1][nearestDig->GetValue() % 1024]) * 24.0 / 1024;
				}
				else if(tM. slot == 5){
					tL = (tdcDig1->GetValue() + fINLTable32[rChannel1][tdcDig1->GetValue() % 1024]) * 24.0 / 1024;
					tT = (nearestDig->GetValue() + fINLTable32[rChannel1][nearestDig->GetValue() % 1024]) * 24.0 / 1024;
				}
				new ((*trigAr)[trigAr->GetEntriesFast()]) BmnTrigDigit(iMod, tL, tT - tL);
			}
		}
	}
	return kBMNSUCCESS;
}

BmnStatus BmnTrigRaw2Digit::ClearArrays() {
	for (TClonesArray *ar : trigArrays)
		ar->Clear("C");
}

ClassImp(BmnTrigRaw2Digit)

