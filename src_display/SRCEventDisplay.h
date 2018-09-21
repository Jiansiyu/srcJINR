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

class EventDisplay : public TGMainFrame {
	private:
		TRootEmbeddedCanvas  *fEcan;
		TGStatusBar          *fStatusBar;
	public:
		EventDisplay(const TGWindow *p, UInt_t w, UInt_t h);
		virtual ~EventDisplay();
		int graph;
		int divide;
		void DoExit();
		void SetStatusText(const char *txt, Int_t pi);
		void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);
		ClassDef(EventDisplay, 0)
};





