#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "Rtypes.h"
#include "TVectorT.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;

int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tfileQuality checkList.txt\n";
		return -1;
	}


	
	// Attempt to find a matchin pedestal file:
	string fileName = string(argv[1]);
	ifstream list (fileName);
	string run_number;

	std::vector<double> fileNo;
	std::vector<double> rate_t0;
	std::vector<double> pedestals_bc1;
	std::vector<double> pedestals_bc2;
	std::vector<double> pedestals_bc3;
	std::vector<double> pedestals_bc4;
	std::vector<double> tdiff_bc1;
	std::vector<double> tdiff_bc2;
	std::vector<double> tdiff_bc3;
	std::vector<double> tdiff_bc4;
	std::vector<double> tres_bc1;
	std::vector<double> tres_bc2;
	std::vector<double> tres_bc3;
	std::vector<double> tres_bc4;
	std::vector<double> tped_bc1;
	std::vector<double> tped_bc2;
	std::vector<double> tped_bc3;
	std::vector<double> tped_bc4;
	std::vector<double> tdc_bc1;
	std::vector<double> tdc_bc2;
	std::vector<double> tdc_bc3;
	std::vector<double> tdc_bc4;
	std::vector<double> carbon_bc1;
	std::vector<double> carbon_bc2;
	std::vector<double> carbonWi_bc1;
	std::vector<double> carbonWi_bc2;


	while( std::getline(list,run_number) ){
		cout << "Working on run: " << run_number << "\n";
		
		TFile * infile = NULL;
		
		TString path = std::getenv("VMCWORKDIR");
		path = path + "/build/bin/qualityCheck/checked_" + run_number + ".root";
		
		infile = new TFile(path);
		if (!infile){
			cerr << "Could not open file " << path <<"\n";
			continue;
		}

		TVectorT<double> * ped		= (TVectorT<double>*)	infile->Get("ped");			
		TVectorT<double> * tdcCnt	= (TVectorT<double>*)	infile->Get("tdcCnt");			
		TVectorT<double> * t0Rate	= (TVectorT<double>*)	infile->Get("t0Rate");			
		TVectorT<double> * tDiff	= (TVectorT<double>*)	infile->Get("tDiff");			
		TVectorT<double> * tRes		= (TVectorT<double>*)	infile->Get("tRes");			
		TVectorT<double> * pedTime	= (TVectorT<double>*)	infile->Get("pedTime");			
		TVectorT<double> * carbonIn	= (TVectorT<double>*)	infile->Get("carbonIn");		
		TVectorT<double> * carbonInWidth= (TVectorT<double>*)	infile->Get("carbonInWidth");	
	
		fileNo.push_back(  atof( run_number.c_str() ) );

		rate_t0.push_back( (*t0Rate)[0] );

		pedestals_bc1.push_back(	(*ped)[0]	);
		pedestals_bc2.push_back(	(*ped)[1]	);
		pedestals_bc3.push_back(	(*ped)[2]	);
		pedestals_bc4.push_back(	(*ped)[3]	);

		tdiff_bc1.push_back(	(*tDiff)[0]	);
		tdiff_bc2.push_back(	(*tDiff)[1]	);
		tdiff_bc3.push_back(	(*tDiff)[2]	);
		tdiff_bc4.push_back(	(*tDiff)[3]	);

		tres_bc1.push_back(	(*tRes)[0]	);
		tres_bc2.push_back(	(*tRes)[1]	);
		tres_bc3.push_back(	(*tRes)[2]	);
		tres_bc4.push_back(	(*tRes)[3]	);

		tped_bc1.push_back(	(*pedTime)[0]	);
		tped_bc2.push_back(	(*pedTime)[1]	);
		tped_bc3.push_back(	(*pedTime)[2]	);
		tped_bc4.push_back(	(*pedTime)[3]	);

		tdc_bc1.push_back(	(*tdcCnt)[0]	);
		tdc_bc2.push_back(	(*tdcCnt)[1]	);
		tdc_bc3.push_back(	(*tdcCnt)[2]	);
		tdc_bc4.push_back(	(*tdcCnt)[3]	);

		carbon_bc1.push_back(	(*carbonIn)[0]	);
		carbon_bc2.push_back(	(*carbonIn)[1]	);

		carbonWi_bc1.push_back(	(*carbonInWidth)[0]	);
		carbonWi_bc2.push_back(	(*carbonInWidth)[1]	);

		infile->Close();
	}
	
	int dim = fileNo.size();
	TGraph * hT0Rate = new TGraph(dim , &fileNo[0] , &rate_t0[0] );

	TGraph * hPedBC1 = new TGraph(dim , &fileNo[0] , &pedestals_bc1[0] );
	TGraph * hPedBC2 = new TGraph(dim , &fileNo[0] , &pedestals_bc2[0] );
	TGraph * hPedBC3 = new TGraph(dim , &fileNo[0] , &pedestals_bc3[0] );
	TGraph * hPedBC4 = new TGraph(dim , &fileNo[0] , &pedestals_bc4[0] );
	
	TGraph * hTdcBC1 = new TGraph(dim , &fileNo[0] , &tdc_bc1[0] );
	TGraph * hTdcBC2 = new TGraph(dim , &fileNo[0] , &tdc_bc2[0] );
	TGraph * hTdcBC3 = new TGraph(dim , &fileNo[0] , &tdc_bc3[0] );
	TGraph * hTdcBC4 = new TGraph(dim , &fileNo[0] , &tdc_bc4[0] );

	std::vector<double> trash;
	for( int i = 0 ; i < dim ; i++)
		trash.push_back(0);

	TGraphErrors * hTdiff_BC1 = new TGraphErrors(dim, &fileNo[0] , &tdiff_bc1[0] , &trash[0] , &tres_bc1[0] );
	TGraphErrors * hTdiff_BC2 = new TGraphErrors(dim, &fileNo[0] , &tdiff_bc2[0] , &trash[0] , &tres_bc2[0] );
	TGraphErrors * hTdiff_BC3 = new TGraphErrors(dim, &fileNo[0] , &tdiff_bc3[0] , &trash[0] , &tres_bc3[0] );
	TGraphErrors * hTdiff_BC4 = new TGraphErrors(dim, &fileNo[0] , &tdiff_bc4[0] , &trash[0] , &tres_bc4[0] );
	
	TGraphErrors * hCarbonBC1 = new TGraphErrors(dim, &fileNo[0] , &carbon_bc1[0] , &trash[0] , &carbonWi_bc1[0]);
	TGraphErrors * hCarbonBC2 = new TGraphErrors(dim, &fileNo[0] , &carbon_bc2[0] , &trash[0] , &carbonWi_bc2[0]);


	TApplication theApp("App",&argc,argv);
	
	TCanvas * c1 = new TCanvas("c1");
	hT0Rate->SetMarkerStyle(20);
	hT0Rate->Draw("AP");
	hT0Rate->GetXaxis()->SetTitle("Run Number");
	hT0Rate->GetYaxis()->SetTitleOffset(1.5);
	hT0Rate->GetYaxis()->SetRangeUser(0,1);
	hT0Rate->GetYaxis()->SetTitle("Fraction");
	hT0Rate->SetTitle("Fraction of Events where T0 has Single TDC");
	c1->Update();
	
	TCanvas * c2 = new TCanvas("c2");
	c2->Divide(2,2);
	c2->cd(1);
	hPedBC1->SetMarkerStyle(20);
	hPedBC1->Draw("AP");
	hPedBC1->GetXaxis()->SetTitle("Run Number");
	hPedBC1->GetYaxis()->SetTitleOffset(1.5);
	hPedBC1->GetYaxis()->SetTitle("ADC [a.u]");
	hPedBC1->SetTitle("Pedestal of BC1");
	hPedBC1->GetYaxis()->SetRangeUser(40,140);
	c2->cd(2);
	hPedBC2->SetMarkerStyle(20);
	hPedBC2->Draw("AP");
	hPedBC2->GetXaxis()->SetTitle("Run Number");
	hPedBC2->GetYaxis()->SetTitleOffset(1.5);
	hPedBC2->GetYaxis()->SetTitle("ADC [a.u]");
	hPedBC2->SetTitle("Pedestal of BC2");
	hPedBC2->GetYaxis()->SetRangeUser(0,150);
	c2->cd(3);
	hPedBC3->SetMarkerStyle(20);
	hPedBC3->Draw("AP");
	hPedBC3->GetXaxis()->SetTitle("Run Number");
	hPedBC3->GetYaxis()->SetTitleOffset(1.5);
	hPedBC3->GetYaxis()->SetTitle("ADC [a.u]");
	hPedBC3->SetTitle("Pedestal of BC3");
	hPedBC3->GetYaxis()->SetRangeUser(-65,35);
	c2->cd(4);
	hPedBC4->SetMarkerStyle(20);
	hPedBC4->Draw("AP");
	hPedBC4->GetXaxis()->SetTitle("Run Number");
	hPedBC4->GetYaxis()->SetTitleOffset(1.5);
	hPedBC4->GetYaxis()->SetTitle("ADC [a.u]");
	hPedBC4->SetTitle("Pedestal of BC4");
	hPedBC4->GetYaxis()->SetRangeUser(100,200);
	c2->Update();


	TCanvas * c3 = new TCanvas("c3");
	c3->Divide(2,2);
	c3->cd(1);
	hTdcBC1->SetMarkerStyle(20);
	hTdcBC1->Draw("AP");
	hTdcBC1->GetXaxis()->SetTitle("Run Number");
	hTdcBC1->GetYaxis()->SetTitleOffset(1.5);
	hTdcBC1->GetYaxis()->SetTitle("Fraction");
	hTdcBC1->SetTitle("BC1: Fraction of Events with only 1 TDC");
	hTdcBC1->GetYaxis()->SetRangeUser(0,1);
	c3->cd(2);
	hTdcBC2->SetMarkerStyle(20);
	hTdcBC2->Draw("AP");
	hTdcBC2->GetXaxis()->SetTitle("Run Number");
	hTdcBC2->GetYaxis()->SetTitleOffset(1.5);
	hTdcBC2->GetYaxis()->SetTitle("Fraction");
	hTdcBC2->SetTitle("BC2: Fraction of Events with only 1 TDC");
	hTdcBC2->GetYaxis()->SetRangeUser(0,1);
	c3->cd(3);
	hTdcBC3->SetMarkerStyle(20);
	hTdcBC3->Draw("AP");
	hTdcBC3->GetXaxis()->SetTitle("Run Number");
	hTdcBC3->GetYaxis()->SetTitleOffset(1.5);
	hTdcBC3->GetYaxis()->SetTitle("Fraction");
	hTdcBC3->SetTitle("BC3: Fraction of Events with only 1 TDC");
	hTdcBC3->GetYaxis()->SetRangeUser(0,1);
	c3->cd(4);
	hTdcBC4->SetMarkerStyle(20);
	hTdcBC4->Draw("AP");
	hTdcBC4->GetXaxis()->SetTitle("Run Number");
	hTdcBC4->GetYaxis()->SetTitleOffset(1.5);
	hTdcBC4->GetYaxis()->SetTitle("Fraction");
	hTdcBC4->SetTitle("BC4: Fraction of Events with only 1 TDC");
	hTdcBC4->GetYaxis()->SetRangeUser(0,1);
	c3->Update();



	TCanvas * c4 = new TCanvas("c4");
	c4->Divide(2,2);
	c4->cd(1);
	hTdiff_BC1->SetMarkerStyle(20);
	hTdiff_BC1->Draw("AP");
	hTdiff_BC1->GetXaxis()->SetTitle("Run Number");
	hTdiff_BC1->GetYaxis()->SetTitleOffset(1.5);
	hTdiff_BC1->GetYaxis()->SetTitle("T0 - BC1 Time [ns]");
	hTdiff_BC1->SetTitle("Time Difference of BC1 with T0 and Resolution");
	hTdiff_BC1->GetYaxis()->SetRangeUser(5.5,8.5);
	c4->cd(2);
	hTdiff_BC2->SetMarkerStyle(20);
	hTdiff_BC2->Draw("AP");
	hTdiff_BC2->GetXaxis()->SetTitle("Run Number");
	hTdiff_BC2->GetYaxis()->SetTitleOffset(1.5);
	hTdiff_BC2->GetYaxis()->SetTitle("T0 - BC2 Time [ns]");
	hTdiff_BC2->SetTitle("Time Difference of BC2 with T0 and Resolution");
	hTdiff_BC2->GetYaxis()->SetRangeUser(-31,-28);
	c4->cd(3);
	hTdiff_BC3->SetMarkerStyle(20);
	hTdiff_BC3->Draw("AP");
	hTdiff_BC3->GetXaxis()->SetTitle("Run Number");
	hTdiff_BC3->GetYaxis()->SetTitleOffset(1.5);
	hTdiff_BC3->GetYaxis()->SetTitle("T0 - BC3 Time [ns]");
	hTdiff_BC3->SetTitle("Time Difference of BC3 with T0 and Resolution");
	hTdiff_BC3->GetYaxis()->SetRangeUser(7.5,11.5);
	c4->cd(4);
	hTdiff_BC4->SetMarkerStyle(20);
	hTdiff_BC4->Draw("AP");
	hTdiff_BC4->GetXaxis()->SetTitle("Run Number");
	hTdiff_BC4->GetYaxis()->SetTitleOffset(1.5);
	hTdiff_BC4->GetYaxis()->SetTitle("T0 - BC4 Time [ns]");
	hTdiff_BC4->SetTitle("Time Difference of BC4 with T0 and Resolution");
	hTdiff_BC4->GetYaxis()->SetRangeUser(20,25);
	c4->Update();


	TCanvas * c5 = new TCanvas("c5");
	c5->Divide(2,1);
	c5->cd(1);
	hCarbonBC1->SetMarkerStyle(20);
	hCarbonBC1->Draw("AP");
	hCarbonBC1->GetXaxis()->SetTitle("Run Number");
	hCarbonBC1->GetYaxis()->SetTitleOffset(1.5);
	hCarbonBC1->GetYaxis()->SetTitle("ADC [a.u.]");
	hCarbonBC1->SetTitle("Carbon Peak of BC1");
	hCarbonBC1->GetYaxis()->SetRangeUser(500,1500);
	c5->cd(2);
	hCarbonBC2->SetMarkerStyle(20);
	hCarbonBC2->Draw("AP");
	hCarbonBC2->GetXaxis()->SetTitle("Run Number");
	hCarbonBC2->GetYaxis()->SetTitleOffset(1.5);
	hCarbonBC2->GetYaxis()->SetTitle("ADC [a.u.]");
	hCarbonBC2->SetTitle("Carbon Peak of BC2");
	hCarbonBC2->GetYaxis()->SetRangeUser(250,1250);
	c5->Update();


	cout << "BC1 Pedestal: " << accumulate(pedestals_bc1.begin() , pedestals_bc1.end() , 0.)/pedestals_bc1.size() << "\n";
	cout << "BC2 Pedestal: " << accumulate(pedestals_bc2.begin() , pedestals_bc2.end() , 0.)/pedestals_bc2.size() << "\n";
	cout << "BC3 Pedestal: " << accumulate(pedestals_bc3.begin() , pedestals_bc3.end() , 0.)/pedestals_bc3.size() << "\n";
	std::vector<double> sk_pedestals_bc4;
	for( int i = 0 ; i < pedestals_bc4.size() ; i++ )
		if( pedestals_bc4.at(i) < 160 )
			sk_pedestals_bc4.push_back(pedestals_bc4.at(i));

	cout << "BC4 Pedestal: " << accumulate(sk_pedestals_bc4.begin() , sk_pedestals_bc4.end() , 0.)/sk_pedestals_bc4.size() << "\n";


	theApp.Run();

	/*

	*/

	return 0;
}

