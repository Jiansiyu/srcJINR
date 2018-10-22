#include "BmnTOF1Detector.h"
#include <iomanip>

ClassImp(BmnTOF1Detector)

BmnTOF1Detector::BmnTOF1Detector(){
	// Empty constructor
}

//----------------------------------------------------------------------------------------

BmnTOF1Detector::BmnTOF1Detector(int NPlane, int FillHist = 0, TTree *tree = NULL) {
	fFillHist = FillHist;
	fNPlane = NPlane;
	fStripLength = 30; 		// cm
	fSignalVelocity = 0.06; 	// 0.06 ns/cm
	fMaxDelta = (fStripLength * 0.5 + 2.0) * fSignalVelocity; 
					// + 20 mm on the strip edge

	fName = Form("Plane_%d", NPlane);
	memset(fCorrLR, 0, sizeof (fCorrLR));
	memset(fStripShift, 0, sizeof(fStripShift));
	memset(fWalkFunc, 0, sizeof(fWalkFunc));
	memset(fGammaOffset, 0, sizeof (fGammaOffset));
	memset(fKilled, 0, sizeof (fKilled));
}


//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetCorrLR( string pathToFile ) {
	//dir += "/input/TOF400_LRcorr_RUN7_SRC.dat"; 

	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	//cout << "TOF400-Setup: Attempting to open LR corr file: " << dir << "\n"; 

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
		//cout << "\tLoaded LR corr file\n"; 
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
	//cout << "TOF400-Setup: Attempting to open strip shift file: " << dir << "\n"; 
	
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
		//cout << "\tLoaded strip shift file\n";
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
	//cout << "TOF400-Setup: Attempting to open strip walk-function file: " << dir << "\n"; 
	
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
                //cout << "\tLoaded time walk file\n";
        }
        else{
                cout << "\tFailed to find time walk file, setting all corrections to 0...\n";
        } 

}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetGammaOffset( string pathToFile ){}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::TestPrint( int strip ){
	cout << fCorrLR[strip] << " " << fStripShift[strip] << " " << fWalkFunc[strip][0] << " " << fWalkFunc[strip][1] << " "
		<< fWalkFunc[strip][2] << " " << fWalkFunc[strip][3] << "\n";
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::SetGeoFile( string pathToFile ) {
	
	string dir = std::getenv("VMCWORKDIR"); 
	dir += pathToFile;
	//cout << "TOF400-Setup: Attempting to open geometry file: " << dir << "\n"; 

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


void BmnTOF1Detector::SetGeo(BmnTof1GeoUtils *pGeoUtils) { // Used for run_reco_bmn.C
	int UID;
	for (int i = 0; i < fNStr; i++) {
		UID = BmnTOF1Point::GetVolumeUID(0, fNPlane + 1, i + 1); 
		const LStrip1 *pStrip = pGeoUtils->FindStrip(UID);	 
		fCenterStrip[i] = pStrip->center;			 
	}
}



//----------------------------------------------------------------------------------------

void BmnTOF1Detector::ClearHits(){
	memset( tempHitTime, 0, sizeof(tempHitTime));
	memset( tempHitAmps, 0, sizeof(tempHitAmps));
	memset( tempCounter, 0, sizeof(tempCounter));

	memset( fFlagHit, false, sizeof(fFlagHit));
	memset( fTof, 0, sizeof(fTof));
	memset( fAmp, 0, sizeof(fAmp));
	memset( fYPos, 0, sizeof(fYPos));
	memset( fXPos, 0, sizeof(fXPos));	

	memset( stripsFired, -1, sizeof(stripsFired));
	numStripsFired = 0;

	numClusters = 0;
	memset( final_Tof, 0, sizeof(final_Tof));
	memset( final_Amp, 0, sizeof(final_Amp));
	memset( final_YPos, 0, sizeof(final_YPos));
	memset( final_XPos, 0, sizeof(final_XPos));
	
	tmpVector.SetXYZ(0.,0.,0.);
	for (Int_t i = 0; i < fNStr; i++)
		final_Pos[i].SetXYZ(0.,0.,0.);
		
}

void BmnTOF1Detector::InitSkim( BmnTof1Digit* tofDigit ){
	int plane 	= tofDigit->GetPlane();
	int strip 	= tofDigit->GetStrip();
	int side 	= tofDigit->GetSide();
	double tofT 	= tofDigit->GetTime();
	double tofAmp 	= tofDigit->GetAmplitude();
	
	if( plane < 0 || plane > 19) 	return;
	if( strip < 0 || strip > 47) 	return;
	if( fKilled[strip] == true)	return;
	if( side == 1) 			return;

	if( side == 0){ // store info only for LH side 
		if(  tempCounter[strip] > 9 ){
			cerr << "Array not large enough, skipping this entry...\n";
			return;
		}
		tempHitTime[strip][tempCounter[strip]] = tofT;
		tempHitAmps[strip][tempCounter[strip]] = tofAmp;
		tempCounter[strip]++;
	}
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::CreateStripHit( BmnTof1Digit* tofDigit , double t0Time , double t0Amp){
	int plane 	= tofDigit->GetPlane();
	int strip 	= tofDigit->GetStrip();
	int side 	= tofDigit->GetSide();
	double tofT 	= tofDigit->GetTime();
	double tofAmp 	= tofDigit->GetAmplitude();
	
	if( plane < 0 || plane > 19) 	return;
	if( strip < 0 || strip > 47) 	return;
	if( fKilled[strip] == true)	return;
	if( side == 0) 			return;
	
	if( side == 1){ // store info for RH side and now look at LH side to create strip hit

		for( int hit = 0 ; hit < tempCounter[strip] ; hit++){	
			double tmpT = tempHitTime[strip][hit];
			double tmpA = tempHitAmps[strip][hit];

			//cout << tmpT << " " << tofT << "\n"
			//	<< tmpA << " " << tofAmp << "\n"
			//	<< strip << "\n";

			if( fabs( (1./0.06)*(tmpT - tofT + fCorrLR[strip]) ) < 32 ){
				double par0 = fWalkFunc[strip][1];
				double par1 = fWalkFunc[strip][2];
				double par2 = fWalkFunc[strip][3];
				double shift = fWalkFunc[strip][0];
				double meanTime = (tmpT + tofT + fCorrLR[strip])*0.5;
				double sumAmps = sqrt(tmpA * tofAmp);
				double pos = 0.5*(1./0.06)*(tmpT - tofT + fCorrLR[strip]); // +/- 15cm
		
				double fixWalk = par0 + par1*exp( -(sumAmps - shift) / par2 );
				if( par0 == -1 || par1 == -1 || par2 == -1 || shift == -1) fixWalk = 0;
				meanTime = meanTime - t0Time + (-6.1 + 27./sqrt(t0Amp) ) - fStripShift[strip] - fixWalk;
				//cout << "\tcompare to saved tof: " << meanTime << " " << fTof[strip] << " " << fFlagHit[strip] << "\n";

				
				if( sumAmps < 9 ) return; // I don't want hits with low amplitude

				if( ( fFlagHit[strip] == false) || ( fFlagHit[strip] == true && meanTime < fTof[strip] ) ){
					fTof[strip] = meanTime;
					fAmp[strip] = sumAmps;
					fYPos[strip] = pos;
					fXPos[strip] = strip;
					fFlagHit[strip] = true;
					//cout << "\tmatched: " << fTof[strip] << " " << fAmp[strip] << " " <<  fYPos[strip] << " " << fXPos[strip] << " " <<  fFlagHit[strip] << "\n";

					int *result = std::find(std::begin(stripsFired), std::end(stripsFired), strip);
					if( std::find( std::begin(stripsFired) , std::end(stripsFired) , strip) == std::end(stripsFired) ){
						stripsFired[numStripsFired] = strip;
						numStripsFired++;
					}
				}
				
			}

		}	

	}
	
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::ClusterHits(){
	if( numStripsFired == 0) return; // No strips fired in this plane

	std::vector<int> firedStrips;
	std::vector<int> clusterList[48];

	for( int i = 0 ; i < numStripsFired ; i++)
		firedStrips.push_back( stripsFired[i] );

		// Sort which which strips fired
	std::sort( firedStrips.begin() , firedStrips.end() );
	
		// Print option:
	//acout << "Need to consider clustering for # strips: " << numStripsFired << "\n\tstrips: ";
	//for( int i = 0 ; i < numStripsFired ; i++)
	//	cout << firedStrips.at(i) << " ";
	//cout << "\n";

		// if there is only 1 strip fired, save that
	if( numStripsFired == 1 ){
		//cout << "\tonly one strip fired, so save the information in tofCluster and skip to filling histograms\n";
		int strip = firedStrips.at(0);
		final_Tof[numClusters] = fTof[strip];
		final_Amp[numClusters] = fAmp[strip];
		final_YPos[numClusters] = fYPos[strip];
		final_XPos[numClusters] = fXPos[strip];
		
			// Create global position for the strip
		tmpVector.SetXYZ(0., fYPos[strip], 0.);
		final_Pos[numClusters] = fCenterStrip[strip] + tmpVector;
		
		numClusters++;
	}
	else{
		// Sort through multi-strip events to do clustering:
		//cout << "\tmore than one strip fired, so we need to try to cluster them together\n";
		for( int i = 0 ; i < firedStrips.size() ; i++){
				// Initializing cluster list for every strip that has itself in it
			clusterList[firedStrips.at(i)].push_back( firedStrips.at(i) );

			for( int j = 0 ; j < firedStrips.size() ; j++){
				if( j <= i ) continue;
				int stOne = firedStrips.at(i);
				int stTwo = firedStrips.at(j);

				double tDiff 	= fTof[stOne] - fTof[stTwo];
				double yDiff 	= fYPos[stOne] - fYPos[stTwo];
				int xDiff 	= fXPos[stOne] - fXPos[stTwo];
				//cout << fTof[stOne] << " " << fTof[stTwo] << "\n"
				//	<< fYPos[stOne] << " " << fYPos[stTwo] << "\n"
				//	<< fXPos[stOne] << " " << fXPos[stTwo] << "\n"
				//	<< fAmp[stOne] << " " << fAmp[stTwo] << "\n";
				if( (fabs( tDiff ) > 2) || (fabs( yDiff ) > 7) || (fabs( xDiff ) > 6) ){}
				else{
					clusterList[stOne].push_back( stTwo );
				}

			}
				// Print option:
			//cout << "\t\tresulting cluster list of strip: " << firedStrips.at(i) << "\n\t\t\t";
			//for( int j = 0 ; j < clusterList[firedStrips.at(i) ].size() ; j++)
			//	cout << clusterList[firedStrips.at(i) ].at(j) << " ";
			//cout << "\n";
			
		}
		        // Now with our clusterList, we need to do an intersection search and union clusters that are similar
		std::vector< std::vector<int> > result;
		//cout << "\tNow doing intersection/union for the cluster lists to get final result clusters\n";
			// For all strips, take the cluster list
		for( int idx = 0 ; idx <  firedStrips.size() ; idx++){
			int st = firedStrips.at(idx);
			//cout << "\t\tWorking on strip "  << st << "\n";
				// For all the groups I already have, find out if any intersection with this cluster list, and if so, add it to group
			bool insert = true;
			for( int group = 0 ; group < result.size() ; group++){
				std::vector<int> tmp;
				std::vector<int> newGrp;
				std::set_intersection( clusterList[st].begin() , clusterList[st].end(), \
                                                                        result.at(group).begin() , result.at(group).end() , back_inserter(tmp) ) ;
				
				//cout << "\t\t\tintersection of how many elements?: " << tmp.size() << "\n";
				if( tmp.size() ){
						// Union two sets into the group and exit
					std::set_union( clusterList[st].begin() , clusterList[st].end(), \
                                                                                result.at(group).begin() , result.at(group).end() , back_inserter( newGrp ) );	
					insert = false;
					result.at(group).swap(newGrp);
					break;
				}
			}
				// Create new group if no group wiht an intersection
			if( insert ){
				if( clusterList[st].size() > 0)
					result.push_back(  clusterList[st] );
			}
		}

		//cout << "\tfinal number of clusters: " << result.size() << "\n";
		for( int group = 0 ; group < result.size() ; group++){
				// For each cluster, choose a representative strip 
			//cout << "\t\tcluster " << group << ": ";
			//cout << "\t\tNumber of strips in cluster: " << result.at(group).size() << "\n";
			int strip = result.at(group).at(0);
			double sumCharge 	= fAmp[strip];
			double weightTime 	= fTof[strip] * fAmp[strip];
			double weightPos	= fYPos[strip] * fAmp[strip];
			
			int repStrip = strip;
			int repAmp = sumCharge;

			//cout << "\t\tCreating representative for cluster:\n";
			//cout << "\t\t\t" << fAmp[strip] << " " << fTof[strip] << " " << fYPos[strip] << "\n";
			for( int hit = 1 ; hit < result.at(group).size() ; hit++){
				strip = result.at(group).at(hit);
				//cout << "\t\t\t" << fAmp[strip] << " " << fTof[strip] << " " << fYPos[strip] << "\n";
				sumCharge	+= fAmp[strip];
				weightTime	+= fTof[strip] * fAmp[strip];
				weightPos	+= fYPos[strip] * fAmp[strip];

				if(fAmp[strip] > repAmp){
					repAmp = fAmp[strip];
					repStrip = strip;
				}
			}
			//cout << "\t\tResult: " << sumCharge / result.at(group).size() << " " << weightTime / sumCharge << " " << weightPos / sumCharge << "\n";
			final_Tof[numClusters] = weightTime / sumCharge;
			final_Amp[numClusters] = sumCharge / result.at(group).size();
			final_YPos[numClusters] = weightPos / sumCharge;
			final_XPos[numClusters] = repStrip;
			
				// Create global position for the strip
			tmpVector.SetXYZ(0., weightPos / sumCharge, 0.);
			final_Pos[numClusters] = fCenterStrip[repStrip] + tmpVector;
			
			numClusters++;
			
		}

	}

	return;
}

//----------------------------------------------------------------------------------------

void BmnTOF1Detector::OutputToTree(TClonesArray *TofHit){

	// For all the cluster strips we have, add them all to a branch
	for( int clust = 0 ; clust < numClusters ; clust++){
		new((*TofHit)[TofHit->GetEntriesFast()]) BmnTOF1Conteiner(
				fNPlane,
				final_XPos[clust],
				final_Tof[clust],
				final_Amp[clust],
				final_XPos[clust],
				final_YPos[clust],
				0,
				final_Pos[clust].x(),
				final_Pos[clust].y(),
				final_Pos[clust].z()
				);

	}
}


//------------------------------------------------------------------------------------------------------------------------
