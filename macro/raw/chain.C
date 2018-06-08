#include <istream>
#include <fstream>

void chain(){
	TChain*c=new TChain("cbmsim");
	for( int i = 3400; i < 3500 ; i++){//for(int i=2227; i<2258; i++){
		ifstream f_file;
		TString name = Form("bmn_run%d_digi.root",i);
		f_file.open(name);
		if (f_file.is_open() == kTRUE)
			c->Add(Form("bmn_run%d_digi.root",i));
		f_file.close();
	}
	TFile*f=new TFile("bmn_runCHAIN_digi.root","recreate");
	c->Write();
	f->Close();
}
