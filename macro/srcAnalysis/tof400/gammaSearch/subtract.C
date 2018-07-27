void subtract(TString file1 , TString file2, int plane){
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	TFile * wall_Fi = new TFile(file1);
	TFile * nowa_Fi = new TFile(file2);

	
	TString hName;
	hName = Form("hToF_All_%i_0",plane);
	TH1D * wall = (TH1D* ) wall_Fi->Get(hName);
	TH1D * nowall = (TH1D* ) nowa_Fi->Get(hName);

	hName = Form("hToF_Singles_%i_0",plane);
	TH1D * wall2 = (TH1D* ) wall_Fi->Get(hName);
	TH1D * nowall2 = (TH1D* ) nowa_Fi->Get(hName);
	
	wall->Add(wall2,-1);
	nowall->Add(nowall2,-1);

	for( int st = 1 ; st < 48 ; st++){
	
		hName = Form("hToF_Singles_%i_%i",plane,st);
		TH1D * data_wa1 = (TH1D *)wall_Fi->Get(hName);
		hName = Form("hToF_All_%i_%i",plane,st);
		TH1D * data_wa2 = (TH1D *)wall_Fi->Get(hName);
		
		data_wa2->Add(data_wa1,-1);		
		wall->Add(data_wa2,1);
		wall2->Add(data_wa1,1);

		hName = Form("hToF_Singles_%i_%i",plane,st);
		TH1D * data_n1 = (TH1D *)nowa_Fi->Get(hName);
		
		hName = Form("hToF_All_%i_%i",plane,st);
		TH1D * data_n2 = (TH1D *)nowa_Fi->Get(hName);
	
		data_n2->Add(data_n1,-1);
		nowall->Add(data_n2,1);
		nowall2->Add(data_n1,1);

		delete data_wa1, data_wa2, data_n1, data_n2;
	}

	TCanvas * c = new TCanvas("c");
	//wall->SetStats(0);
	//wall2->SetStats(0);
	//nowall->SetStats(0);
	//nowall2->SetStats(0);
	wall->GetXaxis()->SetRangeUser(-2,8);
	wall->SetTitle("");
	//wall->Draw();				// multi hits wall = blue
	wall2->SetLineColor(2);			// single hits wall = red
	//wall2->Draw("hist,same");
						// multi hits no wall = green
	nowall->Scale( wall->GetEntries() / nowall->GetEntries() );
	nowall->SetLineColor(8);
	//nowall->Draw("hist,same");

						// sing hits no wall = purple/pink
	nowall2->Scale( wall2->GetEntries() / nowall2->GetEntries() );
	nowall2->SetLineColor(6);
	//nowall2->Draw("hist,same");

	// Subtrack blue from green:
	wall->Add(nowall,-1);
	wall2->Add(nowall2,-1);
	//wall->Add(wall2,-1);
	//wall->Draw("hist");
	//wall2->Draw("same,hist");
	
	wall->Add(wall2,-1);
	wall->Rebin(2);
	wall->Draw("hist");
	
	TF1 * testFit = new TF1("testFit","gaus",-1.5,1.5);
	testFit->SetParameter(0,wall->GetMaximum() );
	testFit->SetParameter(1,0);
	testFit->SetParameter(2,0.2);
	wall->Fit("testFit","","QESR");
	double par[3];
	testFit->GetParameters(&par[0]);
	testFit->SetParameters(par);
	testFit->Draw("same");

	cout << sqrt( pow( par[2] - fabs(par[1]) , 2) - pow( 0.120 , 2) )*1000 << "\n";
	
	//wall->Draw("hist");
	//nowall->Draw("hist,same");
	wall->GetXaxis()->SetRangeUser(-2,8);
	c->Update();

}
