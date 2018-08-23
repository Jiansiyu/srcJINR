void subtract(TString file1 , TString file2, int plane){
	bool drawInd = true;

	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	TFile * wall_Fi = new TFile(file1);
	TFile * nowa_Fi = new TFile(file2);

	
	TString hName;
	hName = Form("hToF_ClusteredEvents_%i_0",plane);
	TH1D * multiClusterEvents_wall 		= (TH1D* ) wall_Fi->Get(hName);	// wall = clustered events with Pb-wall
	TH1D * multiClusterEvents_nowall 	= (TH1D* ) nowa_Fi->Get(hName);

	hName = Form("hToF_SingleEvents_%i_0",plane);
	TH1D * singleClusterEvents_wall 	= (TH1D* ) wall_Fi->Get(hName);
	TH1D * singleClusterEvents_nowall 	= (TH1D* ) nowa_Fi->Get(hName);	

	for( int st = 1 ; st < 48 ; st++){
		// Add all the strips for wall data file
		hName = Form("hToF_ClusteredEvents_%i_%i",plane,st);
		TH1D * data_multiCluster_wall 		= (TH1D *)wall_Fi->Get(hName);
		hName = Form("hToF_SingleEvents_%i_%i",plane,st);
		TH1D * data_singleCluster_wall 		= (TH1D *)wall_Fi->Get(hName);
		
		multiClusterEvents_wall->Add( data_multiCluster_wall );
		singleClusterEvents_wall->Add( data_singleCluster_wall );

		// Add all the strips for NO WALL data file
		hName = Form("hToF_ClusteredEvents_%i_%i",plane,st);
		TH1D * data_multiCluster_nowall 	= (TH1D *)nowa_Fi->Get(hName);
		hName = Form("hToF_SingleEvents_%i_%i",plane,st);
		TH1D * data_singleCluster_nowall 	= (TH1D *)nowa_Fi->Get(hName);
		
		multiClusterEvents_nowall->Add( data_multiCluster_nowall );
		singleClusterEvents_nowall->Add( data_singleCluster_nowall );

		delete data_multiCluster_wall, data_singleCluster_wall;
		delete data_multiCluster_nowall, data_singleCluster_nowall;
	}

	singleClusterEvents_wall->SetStats(0);
	multiClusterEvents_wall->SetStats(0);
	singleClusterEvents_nowall->SetStats(0);
	multiClusterEvents_nowall->SetStats(0);

	// Look at Pb and No-Pb wall data separately, comparing multi cluster and single cluster event times
	/*
	TCanvas * c = new TCanvas("Comparing Single and Multi Times - Pb Wall");
	singleClusterEvents_wall->SetTitle("Comparing Single and Multi Times - Pb Wall");
	singleClusterEvents_wall->Draw("hist");
	multiClusterEvents_wall->SetLineColor(2);
	multiClusterEvents_wall->Draw("hist,same");
	c->Update();
	
	TCanvas * c2 = new TCanvas("Comparing Single and Multi Times - No Pb Wall");
	singleClusterEvents_nowall->SetTitle("Comparing Single and Multi Times - No Pb Wall");
	singleClusterEvents_nowall->Draw("hist");
	multiClusterEvents_nowall->SetLineColor(2);
	multiClusterEvents_nowall->Draw("hist,same");
	c2->Update();
	*/

	// Now compare single cluster for Pb and No-Pb wall data
	singleClusterEvents_nowall->Scale( singleClusterEvents_wall->GetEntries() / singleClusterEvents_nowall->GetEntries() );
	multiClusterEvents_nowall->Scale( multiClusterEvents_wall->GetEntries() / multiClusterEvents_nowall->GetEntries() );
	

	TCanvas * c3 = new TCanvas("Comparing Single Cluster Times for Pb and No-Pb wall data");
	singleClusterEvents_wall->SetTitle("Comparing Single Cluster Times for Pb and No-Pb wall data");
	
	//Draw them individually or do subtraction of two and fit:
	if( drawInd){
		singleClusterEvents_wall->Draw("hist");
		singleClusterEvents_wall->GetXaxis()->SetRangeUser(-2,2);
		singleClusterEvents_nowall->SetLineColor(2);
		singleClusterEvents_nowall->Draw("hist,same");
	}
	else{
		singleClusterEvents_wall->Add( singleClusterEvents_nowall , -1 );
		singleClusterEvents_wall->GetXaxis()->SetRangeUser(-2,2);
			// Doing fit:
		TF1 * testFit = new TF1("testFit","gaus",-1.5,1.5);
		testFit->SetParameter(0, singleClusterEvents_wall->GetMaximum() );
		testFit->SetParameter(1,0);
		testFit->SetParameter(2,0.2);
		singleClusterEvents_wall->Fit("testFit","","QESR");
		double par[3];
		testFit->GetParameters(&par[0]);
		testFit->SetParameters(par);
		
		singleClusterEvents_wall->Draw("hist");
		testFit->Draw("same");
		cout << "Fit Results [in ps]:\n"
			<< "\tOffset: " << par[1]*1000. << "\n\tSigma: " << par[2]*1000. << "\n";
	}

	
	TCanvas * c4 = new TCanvas("Comparing Multiple Cluster Times for Pb and No-Pb wall data");
	multiClusterEvents_wall->SetTitle("Comparing Multiple Cluster Times for Pb and No-Pb wall data");
	
	if( drawInd){
		multiClusterEvents_wall->Draw("hist");
		multiClusterEvents_wall->GetXaxis()->SetRangeUser(-2,2);
		multiClusterEvents_nowall->SetLineColor(2);
		multiClusterEvents_nowall->Draw("hist,same");
	}
	else{
		multiClusterEvents_wall->Add( multiClusterEvents_nowall , -1 );
		multiClusterEvents_wall->GetXaxis()->SetRangeUser(-2,2);
			// Doing fit:
		TF1 * testFit = new TF1("testFit","gaus",-1.5,1.5);
		testFit->SetParameter(0, multiClusterEvents_wall->GetMaximum() );
		testFit->SetParameter(1,0);
		testFit->SetParameter(2,0.2);
		multiClusterEvents_wall->Fit("testFit","","QESR");
		double par[3];
		testFit->GetParameters(&par[0]);
		testFit->SetParameters(par);

		multiClusterEvents_wall->Draw("hist");
		testFit->Draw("same");
		cout << "Fit Results [in ps]:\n"
			<< "\tOffset: " << par[1]*1000. << "\n\tSigma: " << par[2]*1000. << "\n";
	}
	c4->Update();
	
	

}
