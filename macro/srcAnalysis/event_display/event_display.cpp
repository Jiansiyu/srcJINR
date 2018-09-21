#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>

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
#include "TRandom3.h"
#include "TSystem.h"
#include "TGFrame.h"
#include "TPad.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGMenu.h"
#include "TRootEmbeddedCanvas.h"
#include "TGDockableFrame.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TGedFrame.h"
#include "TGLObject.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLEmbeddedViewer.h"
#include "TGLViewerBase.h"
#include "TGLCamera.h"
#include "TGLScenePad.h"
#include "TEveGeoNode.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEvePointSet.h"
#include "TEveLine.h"
#include "TEveElement.h"
#include "TGButton.h"
#include "TExec.h"
#include "TGFrame.h"
#include "TGTextEntry.h"
#include "TGTab.h"
#include "TGIcon.h"
#include "TASImage.h"
#include "TGWindow.h"
#include "TLatex.h"
#include "TError.h"
#include "TVector3.h"
#include "TGClient.h"
#include "TFrame.h"
#include "TGStatusBar.h"
#include "TAxis.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"
#include "BmnTof1Digit.h"
#include "BmnTOF1Detector.h"
#include "BmnTOF1Conteiner.h"
#include "BmnMwpcHitFinder.h"
#include "BmnMwpcSegment.h"
#include "BmnMwpcTrackFinder.h"
#include "FairTrackParam.h"
#include "BmnTrack.h"
#include "UniDbRun.h"
#include "SRCEventDisplay.h"


// Some constants we will need for analysis
const double pedBC1 = 69.2885;
const double pedBC2 = -11.7212;
const double pedBC3 = -25.4808;
const double pedBC4 = 126.067;
const int run_period = 7;
const double a_out = 0.00173144;
const double b_out = 0.0384856;
const double c_out = 0.000015362;
const double a_in = 0.020542;
const double b_in = 0.0305108;
const double c_in = 0.0000114953;


using namespace std;


int main(int argc, char ** argv)
{

	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tevent_display /path/to/all/digi/files\n";
		return -1;
	}

	TApplication theApp("runGUI", &argc, argv);
	new EventDisplay(gClient->GetRoot(),800,700);
	theApp.Run();

	return 0;
}

