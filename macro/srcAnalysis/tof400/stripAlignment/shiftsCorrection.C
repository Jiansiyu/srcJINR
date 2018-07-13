/*
   Written by Efrain Segarra
   Adapted version of ToF400 LR correction with some MASSIVE fixes
   When on real data, need to edit the minimum level based on grass of histogram
   */

void shiftsCorrection(TString file = "")
{

	TString inName;
	inName = file;
	cout << "Opening file " << inName << "\n";

	TFile *FileIn = new TFile(inName.Data(), "READ");
	if (FileIn->IsZombie())
	{
		cout << "Couldn't find file, bailing...\n";
		return -1;
	}

	TString hName;
	TH1D ***hDtLR = new TH1D**[20];
	
	
	ofstream f_call;
	TString NameCallFile = "TOF400_LRcorr_RUN7_SRC.dat";
	f_call.open(NameCallFile.Data());
	f_call << "Plane\tStrip\tShift" << endl << "=====================================================" << endl;
	
	for (int Plane = 0; Plane < 20; Plane++)
	{
		cout << "Working on plane " << Plane << "\n";
		hDtLR[Plane] = new TH1D*[48];
		for (int Str = 0; Str < 48; Str++)
		{	
			// For each strip, get the output of run_reco_src.C and grab the
			// output L-R histograms
			hName = Form("hLR_time_%d_%d", Plane, Str);
			hDtLR[Plane][Str] = (TH1D*) FileIn->Get(hName.Data());


			// For each strip, I want to search the histogram for the center bin that will give me the
			// largest integral. I know the strips are 30cm, with 0.06ns/cm signal velocity, so I expect 1.8ns
			// in length. However the distibution is really 2*L/v (3.6ns wide). I saved the histograms in 25 ps
			// bins, so I expect the length of a strip to be 144 bins wide.
			//
			// So I can take integral of histogram in 144 bins wide, and find the startBin that has the largest integral

			int nBins = hDtLR[Plane][Str]->GetNbinsX();
			double currMax = -1;
			int binSt = 0;
			int binEn = 0;
			for( int bi = 1 ; bi <= nBins ; bi++){
				double currInt = hDtLR[Plane][Str] -> Integral(bi,bi+144);	// Integral is [bin1,bin2] so I only want 144 bins
				if( currInt > currMax){
					currMax = currInt;
					binSt = bi;
					binEn = bi+144;
				}
			}
			
			double mean =  -(hDtLR[Plane][Str]->GetXaxis()->GetBinCenter(binEn) + hDtLR[Plane][Str]->GetXaxis()->GetBinCenter(binSt))/2.;
			f_call << Plane << "\t" << Str << "\t" << mean << endl;
			//f_call << Plane << "\t" << Str << "\t" << mean << "\t" << hDtLR[Plane][Str]->GetXaxis()->GetBinCenter(binEn) << "\t" << hDtLR[Plane][Str]->GetXaxis()->GetBinCenter(binSt) << "\t"  <<endl;
		}
	}



	f_call.close();
	



	FileIn->Close();

	
}
