void gammaPeak(TString file1){

	ofstream f_call;
	TString NameCallFile = "TOF400_StripShift_RUN7_SRC.dat";
	f_call.open(NameCallFile.Data());
	f_call << "Plane\tStrip\tShift\tResolution" << endl << "=====================================================" << endl;

	TFile * wallFile = new TFile(file1);
	TString hName;
	
	for( int plane = 0 ; plane < 20 ; plane++){
		cout << "Working on plane " << plane << "...\n";
		for( int strip = 0 ; strip < 48 ; strip++){
			hName = Form("hToF_All_%i_%i",plane,strip);
			TH1D * wall = (TH1D*)wallFile -> Get(hName);
			
			if( wall->GetEntries() < 1000){
				cout << "Not enough entries, plane or strip likely broken. Quitting on " << plane << " " << strip << "...\n";
				f_call << std::setprecision(6) << plane << "\t" << strip << "\t" << -1 << "\t" << -1 << "\n";
				delete wall;
				continue;
			}

			// Let's grab the peak of the histogram of Pb-wall
			double testMax = wall->GetXaxis()->GetBinCenter(wall->GetMaximumBin());
			TF1 * testFit = new TF1("testFit","gaus",-15,testMax+0.250);
			wall->Fit("testFit","QESRN");
			
			double init_par[3];
			double fin_par[3];
			testFit->GetParameters(&init_par[0]);

			TF1 * finFit = new TF1("finFit","gaus",-15,init_par[1]+0.5*init_par[2]);
			wall->Fit("finFit","QESRN");
			finFit->GetParameters(&fin_par[0]);

			int i = 1;
			while( fin_par[2] > init_par[2] ){
				finFit->SetRange( -15 , init_par[1] + (0.5+0.05*i)*init_par[2] );
				wall->Fit("finFit","QESRN");
				finFit->GetParameters(&fin_par[0]);
				if( 0.05*i >= 0.5){ 
					fin_par[2] = init_par[2];
					fin_par[1] = init_par[1];
					fin_par[0] = init_par[0];
					break;
				}
				i++;
			}
			if( fin_par[2] > 0.5 ){
				cout << "Likely error fitting plane " << plane << " strip " << strip << ": " << fin_par[1] << " " << fin_par[2] << "\n";		
				f_call << std::setprecision(6) << plane << "\t" << strip << "\t" << -1 << "\t" << -1 << "\n";
			}
			else{
				f_call << std::setprecision(6) << plane << "\t" << strip << "\t" << fin_par[1] << "\t" << fin_par[2] << "\n";
			}
			
			delete wall;
			delete testFit;
			delete finFit;		
		}
	}
}
