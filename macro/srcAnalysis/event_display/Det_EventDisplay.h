#ifndef __DET_EVENTDISPLAY__
#define __DET_EVENTDISPLAY__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "orawtree.h"
#include "slowctrl.h"
#include "chef.h"
#include "eventdisplaytree.h"
#include <iostream>
#include <string>
#include <kftree.h>
#include "lumigemtree.h"
#include "mwpctree.h"

class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TPad;
class TGraph;
class TGLEmbeddedViewer;
class TGLViewer;
class TGLSAViewer;
class TGLScenePad;
class TGeoManager;
class TEveGeoTopNode;
class TGeoNode;
class TGeoVolume;
class TGTextButton;
class TGSplitButton;
class TGPopupMenu;
class TGWindow;
class TGTextEntry;
class TGIcon;
class TASImage;
class TGLabel;
class TGVerticalFrame;
class TString;
class TCanvas;
class TGMainFrame;
class TEveLine;
class TEvePointSet;

// Classes for handling visualization points
class visPoint
{
	public:
	double x, y, z, att;
};

class visSet
{
	public:

	visSet();
	virtual ~visSet();

	std::vector<visPoint> points;
	TEvePointSet * look;
	double att;
};

class Det_EventDisplay:public Plugin
{

 private:

	// ROOT GUI Stuff
   //TGWindow *viewwin;
	TGCompositeFrame *tab;
 	//TPad *p,  *p_1,*p_2, *p_3, *p_4;
 	//TGraph *g,*wh;
   TGVerticalFrame *vframe;
   TGVerticalFrame *keyframe;
	//TRootEmbeddedCanvas *c;
   //TCanvas *LVcan;
   //TGWindow *LVcan;
	TGLEmbeddedViewer *V;
   TGMainFrame *LVF;
   TGLEmbeddedViewer *LV1, *LV2, *LV3, *LV4;
	//TGeoManager *myGeoManager;
	TEveGeoTopNode *etn;
   TGTextButton *torButton;
   TGTextButton *wcspButton;
   TGTextButton *dwButton;
   TGTextButton *saveButton;
   TGTextEntry *saveName;
   TGSplitButton *camButton;
   TGPopupMenu *camMenu;
   TGSplitButton *wcButton;
   TGPopupMenu *wcMenu;
   TGSplitButton *tdButton;
   TGPopupMenu *tdMenu;
   TGTextButton *perspButton;
   TGTextButton *orthotopButton;
   TGTextButton *orthorearButton;
   TGTextButton *LVButton;
	TGTextButton *RDButton;
   TGTextButton *LEButton;
   // TGTextButton *gtButton;
   TGTextButton *lmButton;
   TGTextButton *symbButton;
	TGTextButton *pbblockButton;
   TGTextButton *mwpcButton;
   TGTextButton *sipmButton;
   TGTextButton *mwspButton;
   TGTextButton *beamButton;
   TGTextButton *targetButton;
   TGTextButton *twelvedegButton;
   TGTextButton *tofButton;
   TGTextButton *wcframeButton;
   TGTextButton *wcallButton;
   TGTextButton *wcwiresButton;
   TGTextButton *tricButton;
   TGTextButton *atButton;
   TGIcon *olymp;
   TGLabel *runl;
	TGLabel *mcl;
   TGLabel *toroidCurrent;
   TGLabel *beamCurrent;
   TGLabel *beamEnergy;
   TGLabel *targetFlow;
   TGLabel *species;
   TGLabel *trigger;
   TGLabel *runNumber;
	TGLabel *key0;
	TGLabel *key1;
	TGLabel *key2;
	TGLabel *key3;
	TGLabel *key4;
	TGLabel *key5;
	TGLabel *key6;
	TGLabel *key7;
	TGLabel *key8;
	TGLabel *key9;
	TGLabel *key10;


   // GDML node containers
   TGeoNode *tofnodes[36];
   TGeoNode *tornodes[8];
   // TGeoNode *gtnodes[2];
   TGeoNode *lmnodes[6];
   TGeoNode *symbnodes[20];
   TGeoNode *mwpcnodes[6];
   TGeoNode *sipmnodes[4];
   //TGeoNode *lcalnodes[40];
   TGeoNode *beamnodes[18];
   TGeoNode *targetnodes[8];
   TGeoNode *wcsecnodes[2];
   TGeoNode *wcframenodes[6];
   TGeoNode *wcgasnodes[6];
   TGeoNode *wcchamnodes[6];
   TGeoNode *wcspnodes[6];
   TGeoNode *wcslnodes[12];
   TGeoNode *wccellnodes[318];
	TGeoNode *wcwinnodes[14];
   TGeoNode *wclayernodes[36];
   TGeoNode *wirenodes[954];
   TGeoNode *mwpcbeamnodes[4];
	TGeoNode *mwpcwinnodes[24];
   TGeoNode *mwpcframenodes[6];
   TGeoNode *mwpcgasnodes[6];
   TGeoNode *mwpcelectnodes[12];
	TGeoNode *pbglassnodes[6];
	TGeoNode *pbblocknodes[8];
	TGeoNode *exitnode;

   // Various utilities for button, processes, etc.
   int torcyc;
   int wcspcyc;
   // int upw[24];
   // int gtcyc;
   int lmcyc;
   int symbcyc;
   int mwpccyc;
	int mwpcwincyc;
	bool mwpcwin;
   int sipmcyc;
   int mwspcyc;
   int beamcyc;
   int targetcyc;
   int twelvedegcyc;
   int tofcyc;
   int wcframecyc;
	int nframes;
	int nwin;
   int wcwirecyc;
	int wcwincyc;
   int wcallcyc;
	int pbblockcyc;
	int pbglasscyc;
	bool pbblock;
	bool pbglass;
   int nsymb; // Handles what SYMB geometry heirarachy is being used
   int curcam;
   bool multi;
   bool noat;
   bool deadw;
   bool awire;
   bool wflag;
   bool wirestate[954];
   bool atwire[954];
   bool deadID[954];
   bool lwinstate;
   bool lwinstarted;
   bool lumionly;
   bool kf;
   bool geoCheck(int gs);
   TString fname;
   TString fnamedef;
   ORTRunInfo * runinfo;
   ORTEventInfo * eventinfo;
   char ext[5];
   bool choosegeo;
   bool stoplogo;
   const char* geofile;
	bool mcmode;
   double * mccur;

	// Point set for visualization
	std::vector<TEvePointSet*> look;

   // Kalman Filter display objects
   KalmanFilterOutput *kftracks;
   std::vector<TEveLine*> kflines;

	// Vector of visualization points and associated utilities
	std::vector<visSet> viss;
	int natts;

   // LumiGEM and MWPC Hits objects
   LumiGEM * lumi;
   MWPC * mwpc;
   std::vector<TEveLine*> llines;
   TEvePointSet *lpset;

   // Slow Control information/functions
   slowctrl::manager* scmanager;
   slowctrl::datum *dtorcur;
   slowctrl::datum *dbeamcur;
   slowctrl::datum *dbeamE;
   slowctrl::datum *dtarflow;
   //slowctrl::datum *vex;

   // Chef Access
   Chef *chef;

 public:
  
	Det_EventDisplay(TTree *in, TTree *out, TFile *inf_, TFile *outf_, TObject *p);
	virtual ~Det_EventDisplay();

   // Geometry handling functions
   void GenerateDet();
   Long_t findnodes();

   // Button and keyboard control functions
   void toroidView(bool rdraw);
   void wcspView(bool rdraw);
   void deadWires();
   void deadWiresBut(bool rdraw);
   void saveImage();
   void Clearandblk();
   // void gtView(bool rdraw);
   void lmView(bool rdraw);
   void symbView(bool rdraw);
	void mwpcwinView(bool rdraw);
   void mwpcView(bool rdraw);
   void sipmView(bool rdraw);
   void mwspView(bool rdraw);
   void beamView(bool rdraw);
   void targetView(bool rdraw);
   void twelvedegView(bool rdraw);
   void tofView(bool rdraw);
   void wcframeView(bool rdraw);
	void wcwinView(bool rdraw);
   void wcallView(bool rdraw);
   void wcwireView(bool rdraw);
	void pbblockView(bool rdraw);
	void pbglassView(bool rdraw);

   // Camera functions
   void Persp();
   void OrthoTop();
   void OrthoRear();
   void TriColor();
	void RedrawReq();
   void atWires();
   void LumiView();
   void LumiViewStart();
   void ResetCameras();
   void camDefaults();

   // Menu handling
   void camSignals(int id);
   void wcSignals(int id);
   void tdSignals(int id);

   // Utilities
   void cursorsOff();
   void lumiViewOff();
   void LumiOnlyCom();
   void drawKFTracks();
   void drawLumiHits();
   void localToGlobalWC(unsigned int id, double * local, double * global);

   // Initialization and point visualization functions
   Long_t idDeadWire(int id);
   Long_t idLiveWire(int id);
	Long_t deadList();
	Long_t addVisPoint(int id,double x, double y, double z,double att = 0);
	void setPoints(std::vector<visPoint> visp,TEvePointSet * look);
	int pointColor(int att);

   // Command line/processing options
   Long_t useGeo(const char* cg);
   Long_t noLogo(bool yn);
   Long_t LumiOnly(bool leo);
   void LumiOnlyCyc();

	// Required functions for cooker/visco
  	Long_t defineHistograms();
  	Long_t startup();
   Long_t startup_kf();
   Long_t startup_mc();
  	Long_t execute();
   Long_t prepare();
 	Long_t finalize();
 	Long_t viscoEvent();
 	virtual Long_t cmdline(char * cmd);

 	ClassDef(Det_EventDisplay,1);

};

#endif

