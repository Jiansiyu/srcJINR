#include "SRCEventDisplay.h"


void EventDisplay::DoExit()
{
	printf("Exit application...\n");
	gApplication->Terminate(0);
}
void EventDisplay::SetStatusText(const char *txt, Int_t pi)
{
	// Set text in status bar.
	fStatusBar->SetText(txt,pi);
}
void EventDisplay::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
	const char *text0, *text1, *text3;
	char text2[50];
	text0 = selected->GetTitle();
	SetStatusText(text0,0);
	text1 = selected->GetName();
	SetStatusText(text1,1);
	if (event == kKeyPress)
		sprintf(text2, "%c", (char) px);
	else
		sprintf(text2, "%d,%d", px, py);
	SetStatusText(text2,2);
	text3 = selected->GetObjectInfo(px,py);
	SetStatusText(text3,3);
}
EventDisplay::EventDisplay(const TGWindow *p, UInt_t w, UInt_t h) :
	TGMainFrame(p, w, h)
{
	tab = new TGCompositeFrame(this,1250,1000);
	
	fEcan = new TGLEmbeddedViewer(this);
	fEcan->TGLViewer::UseLightColorSet();

	/*
	graph = -1;
	divide = -1;
	fEcan = new TRootEmbeddedCanvas(0,this,1250,1000);
	Int_t wid = fEcan->GetCanvasWindowId();
	TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
	fEcan->AdoptCanvas(myc);
	myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","EventDisplay",this,
			"EventInfo(Int_t,Int_t,Int_t,TObject*)");
	*/
	//AddFrame(fEcan, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
	Int_t parts[] = {45, 15, 10, 30};
	fStatusBar = new TGStatusBar(this, 50, 10, kVerticalFrame);
	fStatusBar->SetParts(parts, 4);
	fStatusBar->Draw3DCorner(kFALSE);
	AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

	TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 200, 40);
	TGTextButton *exit = new TGTextButton(hframe, "&Exit GUI");
	exit->Connect("Pressed()", "EventDisplay", this, "DoExit()");
	hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

	AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));
	
	SetWindowName("Embedded Canvas Status Info");
	MapSubwindows();


	MapWindow();
}
EventDisplay::~EventDisplay()
{
	Cleanup();
	delete fEcan;
}
