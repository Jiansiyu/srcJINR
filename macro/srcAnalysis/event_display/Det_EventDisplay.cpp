//
// The OLYMPUS 3D Event Display
//
// Det_EventDisplay.cpp
// Brian S. Henderson (bhender1@mit.edu)
// Massachusetts Institute of Technology
// Last Updated: August 20, 2013
//
// This is the main Event Display plugin, used for visualizing data and tracks,
// as well as events created by the new Monte Carlo propagator.
// Tracks and hit points for the wire chambers, ToFs, and the Monte Carlo
// propagator are produced and added to the common geometry managers in the WC2,
// Det_ToF, and MCPropGeant4 plugins respectively, so please see the
// documentation for the Event Display processes in those files.  Tracks from
// the Geant4e track-fitter are produced by the DrawTrack function in
// TrackFit.cpp, while tracks from Kalman Filter are produced by this program
// (which reads the Kalman Filter plugin output).  LumiGEM hits are also
// produced internally to this program.  Please see descriptions in the
// functions below in the code.
//
// This documentation is also available in a perhaps easier to read format on
// OLYMPUS Wookieepedia:
//
// https://olympus-docu.hiskp.uni-bonn.de/dokuwiki/doku.php?id=plugins:eventdisplay
//
// If you would like anything added for your subsystem, something made more user
// friendly, etc, please send me an e-mail.
//
// PREREQUISITES FOR RUNNING THE EVENT DISPLAY:
//
// To run the Event Display, you must configure ROOT with the "--enable-gdml"
// flag appended when you run the "./configure" script before compiling.
//
// Please note the versions of ROOT more recent than 5.30 have been exhibiting
// some bad behavior with Event Display.  The program may seg fault when two
// element of the plugin try to change the geometry simultaneously, giving a 
// seg fault in ROOT's "ESigHandler".  I am working on a workaround and/or patch
// to ROOT that will help alleviate this, or perhaps it will be fixed in the
// next release of ROOT.  As far as I have tested, the Event Display works with
// minimal fuss with ROOT 5.30/06 (the version installed on ocompile as of
// December 2012).
//
// RECIPES TO RUN FOR THE EVENT DISPLAY:
//
// Run the Event Display with Geant4e tracking:
//       /path/to/bin recipes/EventDisplay.xml /path/to/yourinputROOTfile.root
//
// Run the Event Display with Kalman Filter tracking:
//       /path/to/bin recipes/EventDisplayKF.xml /path/to/yourinputROOTfile.root
//
// Run the Event Display to visualize new Monte Carlo propagated tracks:
//       /path/to/bin/visco recipes/MC/vis_propagate.xml /path/to/yourgeneratedfile.root
//
// Note, the Event Display can be run "as normal" (with either kind of tracking) on
// digitized MC data without any changes.
//
// The Event Display does not require an output ROOT file, but you may specify
// one after the input file if you desire.
//
// When using the Geant4e tracker, WC tracking may be turned off in the interest
// of speed by using the command line argument:
//
// -c WC2:useTracking:false
//
// NOTE: At this time, only one of the tracking methods may be used since
//       each instantiates Geant4 in a different way, and a workaround to make
//       them work simultaneously has not yet been implemented.  When this is
//       done, I will make it so that both can be used at once with a
//       "push-button" switch in the user interface.
//
// NOTE: Do not alter the Event Display recipes unless you are very confident
//       that you know what you are doing, or after you have consulted with me.
//       The ordering of processes is important to ensure that each plugin has
//       the information it needs when it starts.
//
// EVENT DISPLAY INTERFACE:
//
// The Event Display is fully controllable through a GUI incorporated into the
// Visual Cooker ("visco") framework.  The Event Display tab should be
// automatically opened when you run an Event Display recipe, but you may switch
// to tabs produced by other plugins at any time.  Note that by clicking the
// double vertical bar at the far left hand side, you can detach the Event
// Display into a separate window in order to simultaneously use another tab or
// to make the display truly full screen without the visco control window.
//
// A guide to the various controls implemented for the Event Display can be
// found on the OLYMPUS Wookieepedia version of this documentation:
//
// https://olympus-docu.hiskp.uni-bonn.de/dokuwiki/doku.php?id=plugins:eventdisplay
//
// Hits in the various detectors are visualized in the following ways:
//
//    Wire Chambers:
//       Hit Wire:      Entire planar area swept by wire colored red (unless
//                      tri-coloring is on, see below)
//       Tracks:        Drawn as trajectories produced by the tracking codes
//                      (with hits in "sensitive planes" used by the Geant4e
//                      tracker drawn as points along the tracks).  The colors
//                      used for the tracks are as follows:
//                     
//                      Blue:  Protons
//                      Red:   Electrons
//                      Gold:  Positrons
//                      Green: Other or Kalman Filter track (KF tracks currently
//                             don't produce a particle ID)
//
//    ToFs:
//       Hit ToF:       Entire ToF colored red, with a cyan horizontal line with
//                      a central point indicating the vertical position of the
//                      hit determined by the PMT time difference.
//
//    12 Degree Telescopes:
//       LumiGEM Hits:  2D hits are drawn as cyan points, 1D hits are drawn as
//                      cyan lines parallel to the appropriate axis in the GEM
//                      sensitive volumes.
//       MWPC Hits:     2D hits are drawn as cyan points, 1D will be implemented
//                      in the near future
//       Tracks:        Same as Wire Chamber tracks; see above.
//
//    SYMB:
//       NO CURRENT IMPLEMENTATION.  INPUT FROM THE SYMB GROUP AS TO WHAT THEY
//       WOULD LIKE WOULD BE MUCH APPRECIATED.
//
//		Monte Carlo Mode:
//    	Tracks:     	Drawn as trajectories produced produced by querying the
//								Geant4 step points.  The colors used for the tracks are
//								as follows:
//                  	
//                   	Blue:   Protons
//                   	Red:    Electrons
//                   	Gold:   Positrons
//                   	Green:  Photon
//                   	Violet: Other
//
// Down the left-hand column, below the OLYMPUS logo (unless it is turned off,
// see the option below) various parameters of the run are listed, including the
// run number, beam species, beam energy, beam current, target flow, and toroid
// current.  These values are updated every event.  In the very near future, I
// will update this to include the trigger fired on each event.
//
// Below the Run/Event Information is an array of buttons which allows you to
// turn each of the different elements of the detector geometry off in the
// visualization.  This is good for making nice pictures that highlight specific
// elements of the detector or for getting stuff you don't need out of the way
// (e.g. the toroid is off by default).  Multi-component detectors have menu
// buttons that allow you to turn off all or certain parts of the detectors.
//
// Additionally, this column of buttons contains two special options for the
// Wire Chambers.  The "WC Hit Tri-Coloring" button activates a feature in which
// hits in the wire chamber are colored based on how many wires in a cell fired.
// The button will change to indicate when this is feature is on.  The various
// colors signify:
//
//    Green:  All three wires in cell hit
//    Blue:   Any two wires in the cell hit
//    Red:    Only one wire in the cell hit
//
// The second button, "WC All Times", causes the display to show wires which
// fired in the event but registered times outside the accepted hit window by
// coloring them violet.  This button also changes to indicated when it is on.
//
// Along the bottom (above the usual visco controls), are several more Event
// Display controls.  The "Change Camera" menu allows you to switch between a
// free perspective view of the detector and preset top and rear orthographic
// views.  Note that the Event Display will remember any zooming, translations,
// rotations, etc. you have done manually for each of these three individual
// views when you switch between them.  You can reset the views to their 
// defaults by choosing "Reset Main Window Cameras" from this menu.
//
// Next to the camera menu, a button labeled "12 Degree Views" opens a new
// window that will eventually contain dedicated top and side orthographic views
// of the 12 degree arms, but this is a work in progress.  Opening it won't hurt
// anything, but it doesn't do much right now.
//
// In the bottom right, is the Event Display image saver.  By entering a
// filename in the text box and pressing return or clicking "Save Image",
// exactly what you see in the main viewer window at that time will be saved to
// the file you specify.  Your filename must include an extension (to tell ROOT
// what format you want), and must be one of the following 6 types supported by
// ROOT: ".png", ".pdf", ".eps", ".jpg", ".gif", or ".gif+", the latter of which
// appends the image to a GIF file to make animated GIFs.  By default, the Event
// Display automatically generates filenames with the run and event number using
// the PNG format, but if you save one image with a different format it will
// recognize this and change the default to the one you used.  By default, the
// image is saved in the directory that the Event Display command was called in,
// but you can always specify a path relative to this location in your filename.
//
// EVENT DISPLAY COMMAND LINE OPTIONS:
//
// The following command line options may be added when calling the Event
// Display by appending these after the ROOT files in the command line call to
// run the Event Display.
//
// CHANGE THE GEOMETRY FILE:
//
// The Event Display by default uses the most recent survey geometry, which are
// the "survey2013" family of GDML files (as of August 2013).  A variety of GDML
// files are available in src/gdml of your OLYMPUS repository, including the
// defaults.  Note that due to the different ways ROOT and Geant4 read GDML
// files, slightly different GDML files are produced (from the same source data)
// for each program.  The Event Display with work with any of the files named
// in the format "Root*.gdml" available in the src/gdml directory, but will
// likely fail for other geometries with a segmentation fault when it tries to
// access an unavailable geometry element.  The ROOT geometry files provided
// should suffice for almost all conceivable purposes, but if you must visualize
// a different geometry you can speak to me and I can show you other options or
// build in a workaround if absolutely necessary.
//
// To specify a non-default geometry file, use the command line option:
//
// -c EventDisplay:useGeo:\"filename\"
//
// where filename is the filename WITHOUT the .gdml extension.  By default, the
// program looks in your .olympus/shared/gdml directory, so either put any file
// you are using there or specify a relative path to that in the file name. For
// example, if you wanted the nominal position geometry, you would call:
//
// -c EventDisplay:useGeo:\"Root_nominal\"
//
// SAVE SPACE BY DISABLING LOGO AND COLOR KEY:
// 
// By default, the mighty OLYMPUS logo is displayed in the top-left hand corner
// of the Event Display window to make things look super-awesome when showing
// off our data to your friends, family, and colleagues.  Additionally, a color
// key is displayed to tell you what the tracks represent in each mode. If you
// have a small screen, however, these objects may push some of the option
// buttons off your screen making them inaccessible.  If this happens to you,
// you can disable the logo with the command line option:
//
// -c EventDisplay:noLogo:true
//



////////////////////////////////////////////////////////////////////////////////
// Notes on works in progress and other information for the program that      //
// serve as personal notes to myself and might be interesting for experts     //
// looking to see how things are implemented for their subsystem.             //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Dead wires now handled in Event Display init file, see instructions there  //
// (Although, I should migrate them to the WC chamber plugin, next time I am  //
//  in need of a mini-project)                                                //
////////////////////////////////////////////////////////////////////////////////

// V ***NOT YET IMPLEMENTED, THIS IS JUST A NOTE TO MYSELF*** V
// This is a list of the detector nodes used for cross plugin communication
// through the shared loaded geometry.
//
// tofnodes[0]: If lumi event only mode, linewidth = 2; else linewidth = 1;
// wirenodes[0]: If lumi event only mode, linewidth = 2; else linewidth = 1; 
// ^ ***NOT YET IMPLEMENTED, THIS IS JUST A NOTE TO MYSELF*** ^

// BEGIN OF CODE

// Associated header file
#include <Det_EventDisplay.h>

// C++ libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>

// OLYMPUS/cooker framework headers
#include "orawtree.h"
#include "kftree.h"
#include "slowctrl.h"

// ROOT openGL and EVE framework headers
#include "TSystem.h"
#include "TGFrame.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMarker.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGMenu.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGDockableFrame.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TGedFrame.h"
#include "TGLObject.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLEmbeddedViewer.h"
#include "TRootEmbeddedCanvas.h"
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

using namespace std;

// Constructor/Destructor
Det_EventDisplay::Det_EventDisplay(TTree *in, TTree *out, TFile *inf_, TFile *outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
   // Initialize the dead wire array
   for (int i=0; i<954; i++) deadID[i] = false;

   // By default, there is no geometry file input (use the default)  
   choosegeo = false;

   // By default, don't have lumi event only processing on
   lumionly = false;

   // By default, include the logo
   stoplogo = false;

	// By default, we aren't in MC mode
	mcmode = false;

   // By default, we'll use the slow control current rather than a memory object
   mccur = NULL;

   // Normally no Kalman Filter (overridden by startup_kf())
	// DON'T CHANGE THIS HERE!  USE THE KF RECIPE!
	kf = false;

	// Default number of point attribute sets (added as they are found in initialization)
	natts = 0;

};

Det_EventDisplay::~Det_EventDisplay()
{
};

// Instantiate the geometry
void Det_EventDisplay::GenerateDet()
{

   ///////////////////////////////////////////////////////////////////////////////
   // Dead wires now handled in Event Display init file, see instructions there //
   ///////////////////////////////////////////////////////////////////////////////

   // Import the default geometry file or use the command line geometry
	char buffer[1024];
   if (!choosegeo) geofile = "Root_survey2013"; // Default geometry file

	// Write out the file name
   snprintf(buffer,1000,"%s/.olympus/shared/gdml/%s.gdml",getenv("COOKERHOME"),geofile);

	//	TGeoManager::Import(buffer);
	gGeoManager=gEve->GetGeometry(buffer);
	gGeoManager->DefaultColors();
   gGeoManager->CloseGeometry();

   // Wrap the geometry for EVE, making sure to incorporate all needed levels
	etn = new TEveGeoTopNode(gGeoManager,gGeoManager->GetTopNode(),0,7);
	gEve->AddElement(etn);

   std::cout<<"\n\nConfiguring initial visualization settings.\n\nIgnore TGeoManager and TGLViewerbase info/errors...\n\n";

   // Get the pointers to all the nodes you might want to edit
   Long_t geostatus = findnodes();
   
   // Check to see if the geometry reading succeeded completely, abort if not
   bool geoOK = geoCheck((int)geostatus);
   if (!geoOK) abort();

   // Set default colors, visibility, etc.

   // Color the ToFs
   for (int i=0;i<36;i++)
   {
      tofnodes[i]->VisibleDaughters(0);
      tofnodes[i]->GetVolume()->SetLineColor(kGray+3);
   }
   tofView(false);

   // Set the coils to be invisible by default
   toroidView(false);

   // WC Sensitive planes
   wcspView(false);

   // GEM Trackers (on by default)
   // gtView(false);

   // SYMBs (on by default)
   symbView(false);

   // Beam Pipe (on by default)
   beamView(false);

   // Full twelve degree line (on by default, except the MWPC supports/frames)
   twelvedegView(false);
   mwspView(false);

   // Target Chamber, Cell, Collimator, Windows (on by default)
   targetView(false);

	// Wire chamber windows (off by default)
	wcwinView(false);

	// Lead shielding on by default
	if (pbblock) pbblockView(false);

   // Wire chamber frames (off by default)
   wcframeView(true);

   // Make the wires initially invisible execpt dead wires, but set them up for
   // first event
   wcwireView(false);
   awire = true; wcwirecyc++;
   deadWires();

	// Add the visualization points
   if (natts>0)
   {
      TEveElementList * holder =  new TEveElementList();
	   for (int j=0; j<natts; j++)
	   {
		   viss[j].look->SetMarkerSize(0.5);
		   viss[j].look->SetMarkerColor(pointColor(viss[j].att));
		   // cout<<(viss[j].look->GetMainColor())<<"\n";
		   setPoints(viss[j].points,viss[j].look);
		   holder->AddElement(viss[j].look);
	   }
      cout<<"\nAdding display points...\n";
      gEve->AddElement(holder);
      // delete holder;
   }

   // Make the world volume invisible
   gGeoManager->GetTopVolume()->SetVisibility(0);
   gGeoManager->GetTopNode()->SetVisibility(0);

   // Get ROOT to shut the front door
   gGeoManager->SetVerboseLevel(-1);

   // Draw the detector
   curcam = 0; // Start with perspective view (could change this to change default)
	gEve->FullRedraw3D(kTRUE); // Get good visual settings
   camDefaults(); // Set the default views
   gEve->FullRedraw3D(kFALSE); // Do the final redraw

};

// Function which checks if the geometry is OK and produces an appropriate error
// if not
bool Det_EventDisplay::geoCheck(int gs)
{
   // Return true (OK) immediately to avoid unnecessary checks
   if (gs == 0)
   {
      printf("\n\nThe geometry file was successfully indexed!\n\n");
      return true;
   }

   // Reference for the bits representing failures for different geometry
   // elements, fill into chars for later printing
   std::vector<const char*> failbits;
   failbits.push_back ("ToFs");                       // 0:  ToFs
   failbits.push_back ("Toroid Coils");               // 1:  Toroid coils
   failbits.push_back ("WC Sector Volumes");          // 2:  WC sector volumes
   failbits.push_back ("WC Gas Volumes");             // 3:  WC gas volumes
   failbits.push_back ("WC Frames");                  // 4:  WC frame volumes
   failbits.push_back ("WC Chamber Volumes");         // 5:  WC chamber volumes
   failbits.push_back ("WC Superlayers");             // 6:  WC superlayer volumes
   failbits.push_back ("WC Wires");                   // 7:  WC wires
   failbits.push_back ("WC Sensitive Planes");        // 8:  WC detector planes
   failbits.push_back ("GEM Trackers");               // 9:  GEM Trackers
   failbits.push_back ("LumiGEMs");                   // 10: LumiGEMs
   failbits.push_back ("Beamline Elements");          // 11: Beamline components
   failbits.push_back ("SYMB Detectors");             // 12: SYMBs
   failbits.push_back ("MWPC Support Beams");         // 13: MWPC beams
   failbits.push_back ("MWPC Main Volumes");          // 14: MWPC main volumes
   failbits.push_back ("MWPC Frames");                // 15: MWPC frames
   failbits.push_back ("MWPC Gas Volumes");           // 16: MWPC gas volumes
   failbits.push_back ("MWPC Electronics Boxes");     // 17: MWPC electronics
   failbits.push_back ("SiPMs");                      // 18: SiPMs
   failbits.push_back ("Target Chamber and Cell");    // 19: Target/Cell
	failbits.push_back ("WC Windows"); // 20: WC Windows

   printf("\n\n*************************************************************************\n\n");
   printf("  Error in Det_EventDisplay::GenerateDet().  The provided geometry file:\n\n");
   printf("  %s/.olympus/shared/gdml/%s.gdml\n\n",getenv("COOKERHOME"),geofile);
   printf("  does not include all necessary elements for the Event Display.\n\n");
	printf("  The geometry indexer (Det_EventDisplay::findnodes()) failed in the\n");
	printf("  following sections before aborting:\n\n");

   // Find the failure bit by looping over all possible failure
   for (unsigned int p = 0; p < failbits.size(); p++)
   {
      if ((gs & (1ull<<p)) != 0)
      {
			printf("     %s\n",failbits[p]);
      }
   }

   printf("\n  If this occurred for a standard olympus_release 'Root' GDML file,\n");
   printf("  please contact the maintainer of the OLYMPUS Event Display (see the\n");
   printf("  README file in the Event Display source directory for contact\n");
   printf("  information) with the following error code and your geometry file:\n\n");
	printf("  Error Code: %d\n\n",gs);
   printf("************************************************************************\n\n");

   // After printing the error, return false (geometry not OK)
   return false;
   
}

// Command line option for choosing a non-default geometry
Long_t Det_EventDisplay::useGeo(const char* cg)
{
   choosegeo = true;
   geofile = cg;   
   return 0;
}

// Command line option for suppressing the logo to save space in the viewer
Long_t Det_EventDisplay::noLogo(bool yn)
{
   stoplogo = yn;
   return 0;
}

// Methods required for cooker

Long_t Det_EventDisplay::defineHistograms()
{
  // Define histograms (none for now, can't think of any, but if you want
  // anything just let me know)
  return 0;
};

// Main startup function
Long_t Det_EventDisplay::startup()
{
	// Suppress the ROOT warnings that are meaningless
	gErrorIgnoreLevel = 2001;

   // Default save image extension
   sprintf(ext,"png");

   // 12 Degree View off to start
   lwinstate = false; lwinstarted = false;

   // Initialize the slow control manager and channels
   scmanager = (slowctrl::manager*)getMemoryObject("SlowCtrl Manager");
   if (!scmanager)
   {
      printf("\n\n*****************************************************************************\n\n");
      printf("  ERROR: The Event Display plugin requires the slow control plugin to be called\n");
		printf("         earlier in the recipe.");
      printf("\n\n*****************************************************************************\n\n");
      return -1;	
   }
   // Read in each of the slow control data for the event
   dtorcur = scmanager->getLastValidByName("TOR:CurrentIn");
   dbeamcur = scmanager->getLastValidByName("DORIS:Current");
   dbeamE = scmanager->getLastValidByName("DORIS:EDipole");
   dtarflow = scmanager->getLastValidByName("TGT:MFC:LOWER:in");

   // Get the run info
   runinfo = (ORTRunInfo*)getFileObject("RunInfo");

   // Get the LumiGEM hit info and MWPC hit info
   lumi = NULL;
   mwpc = NULL;
   getOutBranchObject("LumiGEMhits", (TObject**)&lumi);
   getOutBranchObject("MWPChits", (TObject**)&mwpc);
   if (!lumi)
   {
      debug(0, "\n\nDet_Event Display::Startup(): No LumiGEM hits found, so they won't be drawn.\n");
      debug(0, "Has the LumiGEM startup() routine been run?\n\n");
   }
   if (!mwpc)
   {
      debug(0, "\n\nDet_Event Display::Startup(): No MWPC hits found, so they won't be drawn.\n");
      debug(0, "Has the MWPC startup() routine been run?\n\n");
   }
   // Make the point set for 2D hits
   lpset = new TEvePointSet();
   //lpset->SetMarkerColor(kCyan);  // Determines color of ToF hits
   lpset->SetMarkerSize(1.5);     // Larger markers than default
   
   // Set cycle counters to defaults and other constant stuff
   torcyc = 3;
   wcspcyc = 1;
   // gtcyc = 0;
   symbcyc = 1;
   beamcyc = 1;
   twelvedegcyc = 1;
   targetcyc = 1;
   wcwirecyc = 0;
   wcallcyc = 0;
   tofcyc = 1;
   wcframecyc = 0;
	wcwincyc = 0;
	mwpcwincyc = 0;
	pbglasscyc = 0;
	pbblockcyc = 1;
   deadw = true;
   multi = false;
   noat = true;
   fnamedef = "Enter filename with extension (e.g. image.pdf)";

   // Add a tab in visco
   tab=addTab("Event Display");
  
   // As long as the tab is on, go nuts with laying things out
   if (tab)
   {

      // Layout frames/viewer
      TGHorizontalFrame *rframe = new TGHorizontalFrame(tab,0,0);
      vframe = new TGVerticalFrame(rframe,0,0);
      TGVerticalFrame *logo = new TGVerticalFrame(vframe,0,0);
      TGHorizontalFrame *wcopt = new TGHorizontalFrame(vframe,0,0);
      //TGVerticalFrame *wcopt1 = new TGVerticalFrame(wcopt,0,0);   // Frame layouts from previous versions
      //TGVerticalFrame *wcopt2 = new TGVerticalFrame(wcopt,0,0);
      //TGHorizontalFrame *tdopt = new TGHorizontalFrame(vframe,0,0);
      //TGVerticalFrame *tdopt1 = new TGVerticalFrame(tdopt,0,0);
      //TGVerticalFrame *tdopt2 = new TGVerticalFrame(tdopt,0,0);
      TGVerticalFrame *bot = new TGVerticalFrame(tab,0,0);
      TGHorizontalFrame *camframe = new TGHorizontalFrame(bot,0,0);
      V = new TGLEmbeddedViewer(rframe);

      // Make the separator for the bottom frame
      TGHorizontal3DLine *sep = new TGHorizontal3DLine(bot);

      // White background
      V->TGLViewer::UseLightColorSet();
      // Return keyboard actions to viewer from other GUI ojects when clicked
      //V->Connect("Clicked()","TGLViewer",V,"HandleEvent()"); // Not yet working

      // Labels
      TGLabel *viewl = new TGLabel(vframe,"Detector Display Options");
      viewl->SetTextJustify(kTextCenterX | kTextCenterY);
      viewl->SetWidth(600);

      runl = new TGLabel(vframe,"Run/Event Information");
      runl->SetTextJustify(kTextCenterX | kTextCenterY);
      runl->SetWidth(600);

		// Mode Label
		if (mcmode)
		{
			mcl = new TGLabel(vframe,"Monte Carlo Mode");
      	mcl->SetTextJustify(kTextCenterX | kTextCenterY);
		   mcl->SetWidth(600);
			mcl->SetTextColor(0xcc0000);
		}
		else if (kf)
		{
			mcl = new TGLabel(vframe,"Kalman Filter Raw Data Mode");
      	mcl->SetTextJustify(kTextCenterX | kTextCenterY);
		   mcl->SetWidth(600);
			mcl->SetTextColor(0xcc0000);
		}
		else
		{
			mcl = new TGLabel(vframe,"Raw Data Mode");
      	mcl->SetTextJustify(kTextCenterX | kTextCenterY);
		   mcl->SetWidth(600);
			mcl->SetTextColor(0xcc0000);
		}

		// Create the particle track color key
		bool key = false;
		if (mcmode && !stoplogo)
		{
	      keyframe = new TGVerticalFrame(vframe,0,0);
			TGHorizontalFrame * top = new TGHorizontalFrame(keyframe,0,0);
			TGHorizontalFrame * bot = new TGHorizontalFrame(keyframe,0,0);
			TGHorizontalFrame * bot2 = new TGHorizontalFrame(keyframe,0,0);

			key0 = new TGLabel(keyframe,"Track Color Key");
			key1 = new TGLabel(top,"Proton");
			key2 = new TGLabel(top,"Positron");
			key3 = new TGLabel(top,"Electron");
			key4 = new TGLabel(bot,"Photon");
			key6 = new TGLabel(bot,"Neutron");
			key7 = new TGLabel(bot,"Mu+");
			key8 = new TGLabel(bot2,"Mu-");
			key9 = new TGLabel(bot2,"Pi+");
			key10 = new TGLabel(bot2,"Pi-");
			key5 = new TGLabel(bot2,"Other");

      	key0->SetTextJustify(kTextCenterX | kTextTop);
      	key0->SetWidth(100);
      	key1->SetTextJustify(kTextCenterX | kTextTop);
      	key1->SetWidth(100);
			key2->SetTextJustify(kTextCenterX | kTextTop);
      	key2->SetWidth(100);
			key3->SetTextJustify(kTextCenterX | kTextTop);
      	key3->SetWidth(100);
			key4->SetTextJustify(kTextCenterX | kTextTop);
      	key4->SetWidth(100);
			key5->SetTextJustify(kTextCenterX | kTextTop);
      	key5->SetWidth(100);
			key6->SetTextJustify(kTextCenterX | kTextTop);
      	key6->SetWidth(100);
			key7->SetTextJustify(kTextCenterX | kTextTop);
      	key7->SetWidth(100);
			key8->SetTextJustify(kTextCenterX | kTextTop);
      	key8->SetWidth(100);
			key9->SetTextJustify(kTextCenterX | kTextTop);
      	key9->SetWidth(100);
			key10->SetTextJustify(kTextCenterX | kTextTop);
      	key10->SetWidth(100);


			key0->SetTextColor((Pixel_t)0x000000);
			key1->SetTextColor(0x0000ff);
			key2->SetTextColor(0x999933);
			key3->SetTextColor(0xff0000);
			key4->SetTextColor(0x00ff00);
			key5->SetTextColor(0xcc00ff);
			key6->SetTextColor(0x00cdcd);
			key7->SetTextColor(0x858585);
			key8->SetTextColor((Pixel_t)0x000000);
			key9->SetTextColor(0x008b8b);
			key10->SetTextColor(0xffa500);

			top->AddFrame(key1, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			top->AddFrame(key2, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			top->AddFrame(key3, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot->AddFrame(key4, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot->AddFrame(key6, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot->AddFrame(key7, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot2->AddFrame(key8, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot2->AddFrame(key9, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot2->AddFrame(key10, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot2->AddFrame(key5, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));

			keyframe->AddFrame(key0, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(top, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(bot, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(bot2, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));

			key = true;
		}
		else if (kf && !stoplogo)
		{
	      keyframe = new TGVerticalFrame(vframe,0,0);

			key0 = new TGLabel(keyframe,"Track Color Key");
			key1 = new TGLabel(keyframe,"All Kalman Filter Tracks");

      	key0->SetTextJustify(kTextCenterX | kTextTop);
      	key0->SetWidth(100);
      	key1->SetTextJustify(kTextCenterX | kTextTop);
      	key1->SetWidth(100);

			key0->SetTextColor((Pixel_t)0x000000);
			key1->SetTextColor(0x00ff00);

			keyframe->AddFrame(key0, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(key1, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));

			key = true;
		}
		else if (!stoplogo)
		{
	      keyframe = new TGVerticalFrame(vframe,0,0);
			TGHorizontalFrame * top = new TGHorizontalFrame(keyframe,0,0);
			TGHorizontalFrame * bot = new TGHorizontalFrame(keyframe,0,0);

			key0 = new TGLabel(keyframe,"Track Color Key");
			key1 = new TGLabel(top,"Proton");
			key2 = new TGLabel(top,"Positron");
			key3 = new TGLabel(bot,"Electron");
			key4 = new TGLabel(bot,"Other");

      	key0->SetTextJustify(kTextCenterX | kTextTop);
      	key0->SetWidth(100);
      	key1->SetTextJustify(kTextCenterX | kTextTop);
      	key1->SetWidth(100);
			key2->SetTextJustify(kTextCenterX | kTextTop);
      	key2->SetWidth(100);
			key3->SetTextJustify(kTextCenterX | kTextTop);
      	key3->SetWidth(100);
			key4->SetTextJustify(kTextCenterX | kTextTop);
      	key4->SetWidth(100);

			key0->SetTextColor((Pixel_t)0x000000);
			key1->SetTextColor(0x0000ff);
			key2->SetTextColor(0x999933);
			key3->SetTextColor(0xff0000);
			key4->SetTextColor(0x00ff00);

			top->AddFrame(key1, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			top->AddFrame(key2, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot->AddFrame(key3, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));
			bot->AddFrame(key4, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,0));

			keyframe->AddFrame(key0, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(top, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));
			keyframe->AddFrame(bot, new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1));

			key = true;
		}

      toroidCurrent = new TGLabel(vframe," ");
      toroidCurrent->SetTextJustify(kTextLeft | kTextTop);
      toroidCurrent->SetWidth(100);
      beamCurrent = new TGLabel(vframe,"(Process events to load data)");
      beamCurrent->SetTextJustify(kTextCenterX | kTextTop);
      beamCurrent->SetWidth(100);
      beamEnergy = new TGLabel(vframe," ");
      beamEnergy->SetTextJustify(kTextLeft | kTextTop);
      beamEnergy->SetWidth(100);
      targetFlow = new TGLabel(vframe," ");
      targetFlow->SetTextJustify(kTextLeft | kTextTop);
      targetFlow->SetWidth(100);
      species = new TGLabel(vframe," ");
      species->SetTextJustify(kTextLeft | kTextTop);
      species->SetWidth(100);
      char bufn[100];
      if (runinfo) sprintf(bufn,"Run Number: %d",runinfo->runNumber);
      else sprintf(bufn,"Run Number: N/A");
      runNumber = new TGLabel(vframe,bufn);
      runNumber->SetTextJustify(kTextLeft | kTextTop);
      runNumber->SetWidth(100);
      // trigger = new TGLabel(vframe," ");
      // trigger->SetTextJustify(kTextLeft | kTextTop);
      // trigger->SetWidth(100);

      toroidCurrent->SetMaxHeight(16);
      beamCurrent->SetMaxHeight(16);
      targetFlow->SetMaxHeight(16);
      beamEnergy->SetMaxHeight(16);
      species->SetMaxHeight(16);
      runNumber->SetMaxHeight(16);
      // trigger->SetMaxHeight(16);

      toroidCurrent->ChangeOptions(toroidCurrent->GetOptions() | kFixedHeight);
      beamCurrent->ChangeOptions(beamCurrent->GetOptions() | kFixedHeight);
      targetFlow->ChangeOptions(targetFlow->GetOptions() | kFixedHeight);
      beamEnergy->ChangeOptions(beamEnergy->GetOptions() | kFixedHeight);
      species->ChangeOptions(species->GetOptions() | kFixedHeight);
      runNumber->ChangeOptions(runNumber->GetOptions() | kFixedHeight);
      // trigger->ChangeOptions(runNumber->GetOptions() | kFixedHeight);

      /*TGLabel *caml = new TGLabel(camframe,"Cameras");
      caml->SetTextJustify(kTextCenterX | kTextCenterY);
      caml->SetWidth(600);*/
      
      // Detector Display buttons
      saveButton = new TGTextButton(camframe, "Save Image",1);
      saveButton->Connect("Clicked()","Det_EventDisplay",this,"saveImage()");
      saveName = new TGTextEntry(fnamedef,camframe);
      saveName->SetDefaultSize(370,23);
      saveName->Connect("DoubleClicked()","Det_EventDisplay",this,"Clearandblk()");
      saveName->Connect("ReturnPressed()","Det_EventDisplay",this,"saveImage()");

      torButton = new TGTextButton(vframe, "Toroid",1);
      torButton->Connect("Clicked()","Det_EventDisplay",this,"toroidView(=true)");
      //torButton->Resize(torButton->GetSize().fWidth,(torButton->GetSize().fHeight)*0.6);
      //torButton->Layout();

      symbButton = new TGTextButton(vframe, "SYMBs",1);
      symbButton->Connect("Clicked()","Det_EventDisplay",this,"symbView(=true)");

      // gtButton = new TGTextButton(vframe, "GEM Tracker",1);
      // gtButton->Connect("Clicked()","Det_EventDisplay",this,"gtView(=true)");

      beamButton = new TGTextButton(vframe, "Beam Pipe",1);
      beamButton->Connect("Clicked()","Det_EventDisplay",this,"beamView(=true)");

      tofButton = new TGTextButton(vframe, "ToFs",1);
      tofButton->Connect("Clicked()","Det_EventDisplay",this,"tofView(=true)");

      targetButton = new TGTextButton(vframe, "Target Chamber and Cell",1);
      targetButton->Connect("Clicked()","Det_EventDisplay",this,"targetView(=true)");

      // Consolidate the WC options into a split button
      wcMenu = new TGPopupMenu(gClient->GetRoot());
      wcMenu->AddEntry("All WC Elements",1);
      wcMenu->AddEntry("Frames",2);
		wcMenu->AddEntry("Windows",6);
      wcMenu->AddEntry("Chamber Planes",3);
      wcMenu->AddEntry("Dead Wires",4);
      wcMenu->AddEntry("All Wires",5);
      wcMenu->ChangeOptions(wcMenu->GetOptions() | kFixedSize);
      wcMenu->SetWidth(targetButton->GetWidth());
      wcButton = new TGSplitButton(wcopt,new TGHotString("            WC Elements    "),wcMenu,false);
      wcButton->Connect("ItemClicked(Int_t)","Det_EventDisplay",this,"wcSignals(Int_t)");
      wcButton->ChangeOptions(wcButton->GetOptions() | kFixedSize);
      wcButton->SetWidth(targetButton->GetWidth());
      //wcButton->SetStyle("classic");
      wcButton->MapWindow();
      //wcButton->SetLeftMargin(100);

      tricButton = new TGTextButton(vframe,"WC Hit Tri-Coloring (OFF)",1);
      tricButton->Connect("Clicked()","Det_EventDisplay",this,"TriColor()");

      atButton = new TGTextButton(vframe,"WC All Times (OFF)",1);
      atButton->Connect("Clicked()","Det_EventDisplay",this,"atWires()");

      // Consolidate the WC options into a split button
      tdMenu = new TGPopupMenu(gClient->GetRoot());
      tdMenu->AddEntry("All 12 Degree Elements",1);
      tdMenu->AddEntry("GEMs",2);
      tdMenu->AddEntry("MWPCs",3);
      tdMenu->AddEntry("SiPMs",4);
      tdMenu->AddEntry("Supports/Electronics",5);
      tdMenu->SetWidth(targetButton->GetWidth());
      tdButton = new TGSplitButton(vframe,new TGHotString("      12 Degree Elements"),tdMenu,false);
      tdButton->Connect("ItemClicked(Int_t)","Det_EventDisplay",this,"tdSignals(Int_t)");
      tdButton->Resize(targetButton->GetWidth(),targetButton->GetHeight());
      //tdButton->SetLeftMargin(100);

      // Camera buttons (consolidated to a single menu button)
      camMenu = new TGPopupMenu(gClient->GetRoot());
      camMenu->AddEntry("Perspective (Free)",1);
      camMenu->AddEntry("Orthographic Top",2);
      camMenu->AddEntry("Orthographic Rear",3);
      camMenu->AddSeparator();
      camMenu->AddEntry("Reset Main Window Cameras",4);
      camButton = new TGSplitButton(camframe,new TGHotString("Change &Camera"),camMenu,false);
      camButton->Connect("ItemClicked(Int_t)","Det_EventDisplay",this,"camSignals(Int_t)");

      // Lumi Window Button
      LVButton = new TGTextButton(camframe,"12 Degree Views",1);
      LVButton->Connect("Clicked()","Det_EventDisplay",this,"LumiView()");

		// Redraw Button
		RDButton = new TGTextButton(camframe,"Redraw Scene",1);
		RDButton->Connect("Clicked()","Det_EventDisplay",this,"RedrawReq()");

      // Old Camera Buttons (kept here for the time being)
      /*perspButton = new TGTextButton(camframe,"Perspective (Free)",1);
      orthotopButton = new TGTextButton(camframe,"Orthographic Top",1);
      orthorearButton = new TGTextButton(camframe,"Orthographic Rear",1);*/
      /*perspButton->Connect("Clicked()","Det_EventDisplay",this,"Persp()");
      orthotopButton->Connect("Clicked()","Det_EventDisplay",this,"OrthoTop()");
      orthorearButton->Connect("Clicked()","Det_EventDisplay",this,"OrthoRear()");*/
      
      // Lumi Events Only Button (currently not fully functional, left out)
      //LEButton = new TGTextButton(camframe,"Lumi Events Only (OFF)",1);
      //LEButton->Connect("Clicked()","Det_EventDisplay",this,"LumiOnlyCyc()");

      // Layout hints  (Padding order for reference L,R,T,B)
      TGLayoutHints *stackleft = new TGLayoutHints(kLHintsCenterY | kLHintsLeft,0,2,0,0);
      TGLayoutHints *stackright = new TGLayoutHints(kLHintsCenterY | kLHintsRight,2,0,0,0);
      TGLayoutHints *stacktop = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsTop,0,0,0,3);
      TGLayoutHints *stacktoplast = new TGLayoutHints(kLHintsCenterX | kLHintsExpandX | kLHintsTop,0,0,0,15);
      TGLayoutHints *stacktight = new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1);
      TGLayoutHints *stacktopover = new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,1);
      TGLayoutHints *stackbottom = new TGLayoutHints(kLHintsCenterX | kLHintsBottom,0,0,2,0);
      TGLayoutHints *submenu = new TGLayoutHints(kLHintsExpandX | kLHintsCenterY | kLHintsLeft, 0,2,0,0);
      TGLayoutHints *fill = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0,0,0,0);

      // Layout the display
      wcopt->AddFrame(wcButton, fill);

      vframe->AddFrame(logo, stacktop);
		vframe->AddFrame(mcl, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 5));
      vframe->AddFrame(runl, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 5));
      vframe->AddFrame(runNumber, stacktight);
      vframe->AddFrame(species, stacktight);
      //vframe->AddFrame(trigger, stacktight);
      vframe->AddFrame(beamEnergy, stacktight);
      vframe->AddFrame(beamCurrent, stacktight);
      vframe->AddFrame(targetFlow, stacktight);
      vframe->AddFrame(toroidCurrent, stacktight);
		if (key) vframe->AddFrame(keyframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 5, 0));
      vframe->AddFrame(viewl, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 5, 5));
      vframe->AddFrame(torButton, stacktop);
      vframe->AddFrame(beamButton, stacktop);
      vframe->AddFrame(targetButton, stacktop);
      vframe->AddFrame(tofButton, stacktop);
      vframe->AddFrame(wcopt, stacktop);
      vframe->AddFrame(tricButton, stacktop);
      vframe->AddFrame(atButton, stacktop);
      vframe->AddFrame(tdButton, stacktop);
      vframe->AddFrame(symbButton, stacktop);
      // vframe->AddFrame(gtButton, stacktoplast);
      vframe->MapSubwindows();
      vframe->MapWindow();

      // Make the almighty and holy OLYMPUS logo (unless told not to at command line)
      if (!stoplogo)
      {
         char filenameBuffer[1024];
         strcpy(filenameBuffer,getenv("COOKERHOME"));
         strcat(filenameBuffer,"/.olympus/shared/EventDisplay/olymp.png");
         int width = (atButton->GetWidth())*1.45;      
         const TGPicture *ipic =(TGPicture *)gClient->GetPicture(filenameBuffer,width,width*1417/2126);
         TGIcon *olymp = new TGIcon(logo,ipic,width,width*1417/2126);
         logo->AddFrame(olymp,new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 0, 0));
      }

      // Layout the bottom bar

      //camframe->AddFrame(caml, new TGLayoutHints(kLHintsCenterY |kLHintsLeft,2,2,0,6));
      camframe->AddFrame(camButton, stackleft);
      /*camframe->AddFrame(perspButton, stackleft);
      camframe->AddFrame(orthotopButton, stackleft);
      camframe->AddFrame(orthorearButton, stackleft);*/
      camframe->AddFrame(LVButton, stackleft);
		camframe->AddFrame(RDButton, stackleft);
		
      //camframe->AddFrame(LEButton, stackleft);
      camframe->AddFrame(saveButton, stackright);
      camframe->AddFrame(saveName, stackright);

      // Make the bottom bar with the separator
      bot->AddFrame(sep, new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,3,3));
      bot->AddFrame(camframe, new TGLayoutHints(kLHintsTop | kLHintsExpandX,5,0,0,0));

      // Place the main containers
      tab->AddFrame(rframe, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,2,0,0));
      tab->AddFrame(bot, new TGLayoutHints(kLHintsExpandX | kLHintsBottom,0,0,0,3));
      rframe->AddFrame(vframe, new TGLayoutHints( kLHintsLeft,5,5,0,0));
      rframe->AddFrame(V->GetFrame(), new TGLayoutHints(kLHintsLeft | kLHintsExpandY|kLHintsExpandX,0,0,0,0));
      
      // Instantiate the necessary EVE classes for track/hit displays
      TEveManager::Create(kFALSE);
      TEveViewer *eve_v=new TEveViewer("Eve Viewer");
      eve_v->SetGLViewer(V,V->GetFrame());
      eve_v->IncDenyDestroy();
      eve_v->AddScene(gEve->GetEventScene());
      gEve->GetViewers()->AddElement(eve_v);

      // Call the detector geometry drawer
      GenerateDet();

		// If MWPC windows were found add them to the 12 deg options
		if (mwpcwin) tdMenu->AddEntry("Windows/Foil",6);
		// Same for the lead glass
		if (pbglass) tdMenu->AddEntry("Lead Glass",7);

		// Add button for lead shielding if present
		if (pbblock)
		{
			pbblockButton = new TGTextButton(vframe, "Lead Shielding",1);
      	pbblockButton->Connect("Clicked()","Det_EventDisplay",this,"pbblockView(=true)");
      	vframe->AddFrame(pbblockButton, stacktop);
		}

      // For MC mode, we'll get the toroid current from a memory object in case
      // it has been manually set to ignore slow control
      if (mcmode)
         mccur = (double *)getMemoryObject("MC Current");

      // Make the Event Display visco tab current
      ((TGTab *) getMemoryObject("Tab Widget"))->SetTab(3);
		if (mcmode) ((TGTab *) getMemoryObject("Tab Widget"))->SetTab(2);

      // Call up the tlumi viewer by default
      // LumiView();

   }

	// Reset the ROOT warning level
	gErrorIgnoreLevel = 0;

//	// Some code to dump out position information

//   int slfw[12] = {0,54,111,189,270,372,477,531,588,666,747,849};

//   cout<<"\n\n";
//   for (unsigned int j=0; j<12; j++)
//   {
//      double loc0[3] = {0,0,0};
//      double locx[3] = {1,0,0};
//      double locy[3] = {0,1,0};
//      double locz[3] = {0,0,1};
//      double globo[3];
//      double glob2[3];
//      localToGlobalWC(slfw[j],loc0,globo);
//      cout<<j<<" "<<globo[0]<<" "<<globo[1]<<" "<<globo[2]<<"\n";
//      localToGlobalWC(slfw[j],locx,glob2);
//      cout<<j<<" "<<glob2[0]-globo[0]<<" "<<glob2[1]-globo[1]<<" "<<glob2[2]-globo[2]<<"\n";
//      localToGlobalWC(slfw[j],locy,glob2);
//      cout<<j<<" "<<glob2[0]-globo[0]<<" "<<glob2[1]-globo[1]<<" "<<glob2[2]-globo[2]<<"\n";
//      localToGlobalWC(slfw[j],locz,glob2);
//      cout<<j<<" "<<glob2[0]-globo[0]<<" "<<glob2[1]-globo[1]<<" "<<glob2[2]-globo[2]<<"\n";
//   }
//   cout<<"\n\n";
//   for (unsigned int j=0; j<954; j++)
//   {
//      double loc0[3] = {0,0,0};
//      double globo[3];
//      localToGlobalWC(j,loc0,globo);
//      cout<<j<<" "<<globo[0]<<" "<<globo[1]<<" "<<globo[2]<<"\n";
//   }
//   exit(0);

//	double loco[3] = {0,0,0};
//	ofstream pos;
//	pos.open("wirepos.txt");
//	for (int j=0; j<954; j++)
//	{
//      for (double y=0; y<100; y+=0.2)
//      {
//         for (int q=0; q<3; q+=2) loco[q] = 0;  
//         loco[1] = y;

//		   double globo[3];
//		   localToGlobalWC(j,loco,globo);

//         double xo0 = globo[0];
//         double yo0 = globo[1];
//		   pos<<j<<" "<<globo[0]<<" "<<globo[1]<<" "<<globo[2]<<"\n";

//         for (int q=0; q<3; q+=2) loco[q] = 0;  
//         loco[1] = -y;

//         localToGlobalWC(j,loco,globo);

//         double xo1 = globo[0];
//         double yo1 = globo[1];
//		   pos<<j<<" "<<globo[0]<<" "<<globo[1]<<" "<<globo[2]<<"\n";

//         double phi0 = atan2(yo0,xo0)*180/M_PI;
//         double phi1 = atan2(yo1,xo1)*180/M_PI;

//         if (phi0<-90) phi0+=360;
//         if (phi1<-90) phi1+=360;

//         if ((phi0<-10) || (phi0>10 && phi0<170) || (phi0>190)) break;
//         if ((phi1<-10) || (phi1>10 && phi1<170) || (phi1>190)) break;
//      }
//	}
//	pos.close();
//  
//   exit(0);

   return 0;
};

// MC Mode Startup Function
Long_t Det_EventDisplay::startup_mc()
{
	// Set MC mode
	mcmode = true;
	
	// Call regular startup
	startup();
};

// Menu handling

// Camera switching menu
void Det_EventDisplay::camSignals(int id)
{

   //std::cout<<"\n\nI am in the menu handler!\n\n";

   switch (id)
   {
      case 1:
         Persp();
         break;
      case 2:
         OrthoTop();
         break;
      case 3:
         OrthoRear();
         break;
      case 4:
         ResetCameras();
         break;
   }
}

// WC elements on/off menu
void Det_EventDisplay::wcSignals(int id)
{
   switch (id)
   {
      case 1:
         wcallView(true);
         break;
      case 2:
         wcframeView(true);
         break;
      case 3:
         wcspView(true);
         break;
      case 4:
         deadWiresBut(true);
         break;
      case 5:
         wcwireView(true);
         break;
		case 6:
			wcwinView(true);
			break;
   }
}

// Twelve degree elements menu
void Det_EventDisplay::tdSignals(int id)
{
   switch (id)
   {
      case 1:
         twelvedegView(true);
         break;
      case 2:
         lmView(true);
         break;
      case 3:
         mwpcView(true);
         break;
      case 4:
         sipmView(true);
         break;
      case 5:
         mwspView(true);
         break;
		case 6:
			mwpcwinView(true);
			break;
		case 7:
			pbglassView(true);
			break;
   }
}

// Start-up routine for Kalman filter recipe
Long_t Det_EventDisplay::startup_kf()
{
   kf = true; // Flag the Kalman Filter on
   startup(); // Regular startup

   // Open the KF tree to get the traces
   kftracks = NULL;
   getOutBranchObject("Tracks",(TObject**)&kftracks);
   if (!kftracks)
   {
      debug(0, "\n\nDet_EventDisplay::startup_kf(): Failed to find the Kalman Filter plugin output.\n");
      debug(0, "Has the Kalman Filter startup routine been run?\n\n");
      return -5;   
   }
   
   return 0;
}

// Cursor control function
void Det_EventDisplay::cursorsOff()
{
   saveName->SetState(false);
};

// Camera button functions

// Standard rotatable perspective view
void Det_EventDisplay::Persp()
{
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
   V->TGLViewer::SetStyle(TGLRnrCtx::kFill);
   curcam = 0;
};

// Orthographic view looking from above
void Det_EventDisplay::OrthoTop() 
{
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   V->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);
   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);
   curcam = 1;
};

// Orthographic view looking down the beam line
void Det_EventDisplay::OrthoRear()
{
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
   V->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);
   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);
   curcam = 2;
};

void Det_EventDisplay::RedrawReq()
{
   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);	
};

// Reset cameras to default views
void Det_EventDisplay::ResetCameras()
{
   // Reset the cameras and then run the initial settings
   V->TGLViewer::ResetCameras();
   camDefaults();

   // Go back to user's initial camera before reset
   switch (curcam)
   {
      case 0:
         V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
         V->TGLViewer::CurrentCamera().RotateRad(-3*M_PI/16,M_PI/4);
         break;
      case 1:
         V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
         break;
      case 2:
         V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
         break;   
   }
   
   // Redraw
   gEve->FullRedraw3D(kFALSE);
}

// Function to set default camera views
void Det_EventDisplay::camDefaults()
{
   // Rear view
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
   V->TGLViewer::ResetCurrentCamera();
   V->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoXOZ,1.2,0,0,M_PI/2,M_PI);

   // Top view
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   V->TGLViewer::ResetCurrentCamera();
   V->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoZOY,0.55,0,0,-M_PI/2,0);

   // Perspective view
   V->TGLViewer::SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
   V->TGLViewer::ResetCurrentCamera();
   V->TGLViewer::CurrentCamera().RotateRad(-3*M_PI/16,M_PI/4); // Rotate to a good view
}

// Functions to set Event Display to lumi events only
Long_t Det_EventDisplay::LumiOnly(bool leo)
{
   lumionly = leo;
   LumiOnlyCom();
   return 0;
}

void Det_EventDisplay::LumiOnlyCyc()
{
   lumionly = !lumionly;
   LumiOnlyCom();  
}

void Det_EventDisplay::LumiOnlyCom()
{
   if (lumionly)
   {
      tofnodes[0]->GetVolume()->SetLineWidth(2);
      wirenodes[0]->GetVolume()->SetLineWidth(2);
      //LEButton->SetTextColor(0x009900);
      //LEButton->SetText("Lumi Events Only (ON)");
   }
   else
   {
      tofnodes[0]->GetVolume()->SetLineWidth(1);
      wirenodes[0]->GetVolume()->SetLineWidth(1);
      //LEButton->SetTextColor(0x000000);
      //LEButton->SetText("Lumi Events Only (OFF)");
   }
}

// Functions to create/close the 12 degree line viewer window

// Build the window, viewers, etc.
void Det_EventDisplay::LumiViewStart()
{
   // Create the window
   LVF = new TGMainFrame(gClient->GetRoot(),800,700);
   //TRootEmbeddedCanvas *LVcan = new TRootEmbeddedCanvas("LVcan",LVF,800,600);
   //LVcan->GetCanvas()->Divide(2,2);
   LVF->Connect("CloseWindow()","Det_EventDisplay",this,"lumiViewOff()");  // Take care of a system close signal ("X" button or Alt-F4)
   //LVF->AddFrame(LVcan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   LVF->SetWindowName("12 Degree Telescope Viewer"); LVF->SetIconName("12 Degree Telescope Viewer");
   
   // Make the layout structures
   TGHorizontalFrame *Lframe = new TGHorizontalFrame(LVF,0,0);
   TGHorizontalFrame *Rframe = new TGHorizontalFrame(LVF,0,0);
   TGHorizontalFrame *Mframe = new TGHorizontalFrame(LVF,0,0);
   TGVerticalFrame *L1 = new TGVerticalFrame(Lframe,0,0);
   TGVerticalFrame *L2 = new TGVerticalFrame(Lframe,0,0);
   TGVerticalFrame *R1 = new TGVerticalFrame(Rframe,0,0);
   TGVerticalFrame *R2 = new TGVerticalFrame(Rframe,0,0);

   // Create lables and geometry viewers for each pad
   LV1 = new TGLEmbeddedViewer(L1); LV1->TGLViewer::UseLightColorSet();
   LV2 = new TGLEmbeddedViewer(L2); LV2->TGLViewer::UseLightColorSet();
   LV3 = new TGLEmbeddedViewer(R1); LV3->TGLViewer::UseLightColorSet();
   LV4 = new TGLEmbeddedViewer(R2); LV4->TGLViewer::UseLightColorSet();

   // Link to the EVE scenes
   TEveViewer *eve_1 = new TEveViewer("Lumi L Top");
   TEveViewer *eve_2 = new TEveViewer("Lumi L Side");
   TEveViewer *eve_3 = new TEveViewer("Lumi R Top");
   TEveViewer *eve_4 = new TEveViewer("Lumi R Side");
   eve_1->SetGLViewer(LV1,LV1->GetFrame());
   eve_2->SetGLViewer(LV2,LV2->GetFrame());
   eve_3->SetGLViewer(LV3,LV3->GetFrame());
   eve_4->SetGLViewer(LV4,LV4->GetFrame());
   eve_1->AddScene(gEve->GetGlobalScene());
   eve_2->AddScene(gEve->GetGlobalScene());
   eve_3->AddScene(gEve->GetGlobalScene());
   eve_4->AddScene(gEve->GetGlobalScene());
   eve_1->AddScene(gEve->GetEventScene());
   eve_2->AddScene(gEve->GetEventScene());
   eve_3->AddScene(gEve->GetEventScene());
   eve_4->AddScene(gEve->GetEventScene());
   gEve->GetViewers()->AddElement(eve_1);
   gEve->GetViewers()->AddElement(eve_2);
   gEve->GetViewers()->AddElement(eve_3);
   gEve->GetViewers()->AddElement(eve_4);

   // Make labels
   TGLabel *ltop = new TGLabel(L1,"Left (Top View, XZ Projection)");
   ltop->SetTextJustify(kTextCenterX | kTextTop);
   ltop->SetWidth(100);
   TGLabel *lside = new TGLabel(L2,"Left (Side View, Beam --> )");
   lside->SetTextJustify(kTextCenterX | kTextTop);
   lside->SetWidth(100);
   TGLabel *rtop = new TGLabel(R1,"Right (Top View, XZ Projection)");
   rtop->SetTextJustify(kTextCenterX | kTextTop);
   rtop->SetWidth(100);
   TGLabel *rside = new TGLabel(R2,"Right (Side View, Beam --> )");
   rside->SetTextJustify(kTextCenterX | kTextTop);
   rside->SetWidth(100);
   ltop->Layout(); lside->Layout(); rtop->Layout(); rside->Layout();

   /*TGLabel *ftest = new TGLabel(Mframe,"!@#$%^&*()<>?:{}[]|/.,");
   ftest->SetTextJustify(kTextCenterX | kTextTop);
   ftest->SetTextFont("symbol");
   Mframe->AddFrame(ftest, new TGLayoutHints(kLHintsCenterY | kLHintsExpandX));*/

   // Layout the window
   TGLayoutHints *stackleft = new TGLayoutHints(kLHintsCenterY | kLHintsLeft,0,2,0,0);
   TGLayoutHints *stacktop = new TGLayoutHints(kLHintsExpandX | kLHintsTop,0,0,0,3);
   L1->AddFrame(LV1->GetFrame(),stacktop); L2->AddFrame(LV2->GetFrame(),stacktop);
   R1->AddFrame(LV3->GetFrame(),stacktop); R2->AddFrame(LV4->GetFrame(),stacktop);
   L1->AddFrame(ltop,stacktop); L2->AddFrame(lside,stacktop);
   R1->AddFrame(rtop,stacktop); R2->AddFrame(rside,stacktop);
   Lframe->AddFrame(L1,stackleft); Lframe->AddFrame(L2,stackleft);
   Rframe->AddFrame(R1,stackleft); Rframe->AddFrame(R2,stackleft);
   L1->Layout(); L2->Layout(); R1->Layout(); R2->Layout();
   LVF->AddFrame(Lframe,stacktop);
   LVF->AddFrame(Rframe,stacktop);
   //LVF->AddFrame(Mframe,stacktop);
   LVF->Layout();

   // Some numbers for positioning the cameras
   static Double_t sx = 354.28; static Double_t sz = 1666.76;
   static Double_t mx = 564.18; static Double_t mz = 2626.51;
   static Double_t centl[3] = {0.,0.,50.}; // {(sx+mx)/2,0.,(mz+sz)/2};
   static Double_t centr[3] = {0.,0.,50.}; //{-(sx+mx)/2,0.,(mz+sz)/2};
   Double_t *clp; clp = &centl[0];
   Double_t *crp; crp = &centr[0];

   // Setup up the cameras for each view
   LV1->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   LV1->TGLViewer::ResetCurrentCamera();
   LV1->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoZOY,1,0,clp,-M_PI/2,-2*M_PI*12/260);
   LV1->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   LV1->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);
   LV3->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   LV3->TGLViewer::ResetCurrentCamera();
   LV3->TGLViewer::SetOrthoCamera(TGLViewer::kCameraOrthoZOY,1,0,crp,-M_PI/2,2*M_PI*12/260);
   LV3->TGLViewer::SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
   LV3->TGLViewer::SetStyle(TGLRnrCtx::kWireFrame);

   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);

}

// Handle the on/off of the lumi/window
void Det_EventDisplay::LumiView()
{

   // If this is the first time the lumi viewer has been turned on, build it
   if (!lwinstarted)
   {
      LumiViewStart(); // Call the startup routine
      lwinstarted = true; // From here on, the viewers, etc already exist
   }

   // If the lumi viewer window is not on, make it on
   if (!lwinstate)
   {
      // Draw the window
      LVF->MapSubwindows();
      LVF->MapWindow();
      lwinstate = true; // Set the window exists flag
   }
   // If it is on, close it
   else if (lwinstate)
   {
      LVF->UnmapWindow(); // Don't actually delete the window, causes problems
      lwinstate = false;
   } 
};

// Function to handle lumi window killed by system button
void Det_EventDisplay::lumiViewOff()
{
   LVF->DontCallClose();  // Block the system close
   LVF->UnmapWindow(); // Hide the window
   lwinstate = false;
};

void Det_EventDisplay::TriColor()
{

   // Immediately do the recoloring based on the previous state
   if (!multi)
   {
      for (int c = 0; c<318; c++)
      {
         if (((wirestate[c*3] && wirestate[c*3+1]) || (wirestate[c*3] && wirestate[c*3+2])) || (wirestate[c*3+2] && wirestate[c*3+1]))
         {
            wirenodes[c*3]->GetVolume()->SetLineColor(kBlue);
            wirenodes[c*3+1]->GetVolume()->SetLineColor(kBlue);
            wirenodes[c*3+2]->GetVolume()->SetLineColor(kBlue);
         }

         if ((wirestate[c*3] && wirestate[c*3+1]) && (wirestate[c*3+2]))
         {
            wirenodes[c*3]->GetVolume()->SetLineColor(kGreen+1);
            wirenodes[c*3+1]->GetVolume()->SetLineColor(kGreen+1);
            wirenodes[c*3+2]->GetVolume()->SetLineColor(kGreen+1);
         }
         
      }

   tricButton->SetTextColor(0x009900);
   tricButton->SetText("WC Hit Tri-Coloring (ON)");
   }
   else
   {
      for (int c = 0; c<954; c++)
      {
         if (wirestate[c]) wirenodes[c]->GetVolume()->SetLineColor(kRed);
      }
   tricButton->SetTextColor(0x000000);
   tricButton->SetText("WC Hit Tri-Coloring (OFF)");
   }   

   // Flip the flag and redraw scene
   deadWires(); // Keep dead wires drawn properly
   multi = !multi;
   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);
};

void Det_EventDisplay::atWires()
{

   // Immediately do the recoloring based on the previous state
   if (noat)
   {
      for (int c = 0; c<954; c++)
      {
         if (atwire[c]) 
         {
            wirenodes[c]->GetVolume()->SetVisibility(1);
         }
			if (wirenodes[c]->GetColour() == kViolet)
			{
				wirenodes[c]->GetVolume()->SetVisibility(1);
				// cout<<"Setting "<<c<<" as visible\n";
			}
		}   
   	atButton->SetTextColor(0x009900);
   	atButton->SetText("WC All Times (ON)");
   }
   else
   {
      for (int c = 0; c<954; c++)
      {
         if (atwire[c]) wirenodes[c]->GetVolume()->SetVisibility(0);
      }
   	atButton->SetTextColor(0x000000);
   	atButton->SetText("WC All Times (OFF)");
   }   

   // Flip the flag and redraw scene
   deadWires(); // Keep dead wires drawn properly
   noat = !noat;
   gEve->GetEventScene()->Changed();
   gEve->FullRedraw3D(kFALSE);
   //gEve->FullRedraw3D(kTRUE);
};

// Image saving function
void Det_EventDisplay::saveImage()
{

      TString fnametemp = saveName->GetText();

      bool gif = fnametemp.EndsWith(".gif");
      bool gifp = fnametemp.EndsWith(".gif+");
      bool jpg = fnametemp.EndsWith(".jpg");
      bool png = fnametemp.EndsWith(".png");
      bool eps = fnametemp.EndsWith(".eps");
      bool pdf = fnametemp.EndsWith(".pdf");
      bool overwrite = fnametemp.EqualTo("File already exists! Click save again to overwrite.");

      if ((((((!gif && !gifp) && !jpg) && !png) && !eps) && !pdf) && !overwrite)
      {
         saveName->SetTextColor(0xFF0000);
         saveName->SetText("Filename must have .gif, .gif+, .jpg, .png, .eps, or .pdf extension");
      }
      else
      {

         // Load filename if not preserving for overwrite
         if (!overwrite) fname = fnametemp;

         // Check if this file name already exists
         char tdir[400];
         sprintf(tdir,"%s/%s",gSystem->WorkingDirectory(),fname.Data());
         Long_t *id=0,*size=0,*flags=0,*mt=0;
         int fexists = gSystem->GetPathInfo(tdir,id,size,flags,mt);
         
         if (fexists!=0 || overwrite) // Doesn't exist or overwrite requested
         {
            V->TGLViewer::SavePicture(fname);
            saveName->SetText("Image saved in current directory.");
            saveName->SetTextColor(kBlack);

            // Reset the default extension to the user's expressed liking
            if (gif) sprintf(ext,"gif"); //ext = "gif";
            else if (gifp) sprintf(ext,"gif+"); //ext ="gif+";
            else if (jpg) sprintf(ext,"jpg"); //ext ="jpg";
            else if (png) sprintf(ext,"png"); //ext ="png";
            else if (eps) sprintf(ext,"eps"); //ext ="eps";
            else if (pdf) sprintf(ext,"pdf"); //ext ="pdf";
         }
         else
         {
            saveName->SetText("File already exists! Click save again to overwrite.");
            saveName->SetTextColor(0xFF0000);
         }
      }

};

// Helper function for image saving text box
void Det_EventDisplay::Clearandblk()
{
   saveName->Clear();
   saveName->SetTextColor(kBlack);
   if (fname) {saveName->SetText(fname);}
};

//
// Display on/off functions for various components
//

// Full twelve degree line as a unit
void Det_EventDisplay::twelvedegView(bool rdraw)
{

   // Set the cycle counters in unison then call individual units (redraw once)
   mwspcyc = twelvedegcyc;
   lmcyc = twelvedegcyc;
   mwpccyc = twelvedegcyc;
   sipmcyc = twelvedegcyc;
	pbglasscyc = twelvedegcyc;
	mwpcwincyc = twelvedegcyc;
   mwspView(false);
   lmView(false);
   mwpcView(false);
   sipmView(false);
	if (mwpcwin) mwpcwinView(false);
	if (pbglass) pbglassView(false);

   twelvedegcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::wcallView(bool rdraw)
{

   // Set the cycle counters in unison then call individual units (redraw once)
   wcwirecyc = wcallcyc;
   wcframecyc = wcallcyc;
	wcwincyc = wcallcyc;
   wcspcyc = wcallcyc;
   if (wcallcyc%2 == 0) {deadw = true;}
   else {deadw = false;}
   wcwireView(false);
   wcframeView(false);
   wcwinView(false);
   wcspView(false);
   deadWiresBut(false);

   wcallcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

// WC wires
void Det_EventDisplay::wcwireView(bool rdraw)
{
   switch (wcwirecyc%2)
   {
      case 0:
         // Set the wires to be invisible
         for (int i=0;i<954;i++) {wirenodes[i]->SetVisibility(0);}
         awire = false;
         break;
      case 1:
         // Set the desired wires to be visible
         for (int i=0;i<954;i++)
         {
            if (wirestate[i]) {wirenodes[i]->SetVisibility(1);}
         }
         deadWires();
         awire = true;
         break;
   }
   wcwirecyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

// WC frames
void Det_EventDisplay::wcframeView(bool rdraw)
{

	// cout<<"\nI think the number of frames is "<<nframes<<"\n";

   switch (wcframecyc%2)
   {
      case 0:
         // Set the frames to be invisible
         for (int i=0;i<nframes;i++) {wcframenodes[i]->SetVisibility(0);} // wcframenodes[i]->VisibleDaughters(1);}
         break;
      case 1:
         // Set the frames to be visible
         for (int i=0;i<nframes;i++) {wcframenodes[i]->GetVolume()->SetVisibility(1); wcframenodes[i]->SetVisibility(1);} // wcframenodes[i]->VisibleDaughters(1);}
         break;
   }
   wcframecyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

// WC windows
void Det_EventDisplay::wcwinView(bool rdraw)
{

	// cout<<"\nI think the number of windows is "<<nwin<<"\n";
   switch (wcwincyc%2)
   {
      case 0:
         // Set the frames to be invisible
         for (int i=0;i<nwin;i++) {wcwinnodes[i]->SetVisibility(0); wcwinnodes[i]->VisibleDaughters(1);}
         break;
      case 1:
         // Set the frames to be visible
         for (int i=0;i<nwin;i++) {wcwinnodes[i]->GetVolume()->SetVisibility(1); wcwinnodes[i]->SetVisibility(1); wcwinnodes[i]->VisibleDaughters(1);}
         break;
   }
   wcwincyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::mwpcView(bool rdraw)
{
   switch (mwpccyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<6;i++) 
         {
            mwpcgasnodes[i]->SetVisibility(0); mwpcgasnodes[i]->VisibleDaughters(0);
            //mwpcframenodes[i]->SetVisibility(0); mwpcframenodes[i]->VisibleDaughters(0);
         }
         break;
      case 1:
         // Set the trackerss to be visible
         for (int i=0;i<6;i++)
         {
            mwpcgasnodes[i]->SetVisibility(1); mwpcgasnodes[i]->VisibleDaughters(1);
            //mwpcframenodes[i]->SetVisibility(1); mwpcframenodes[i]->VisibleDaughters(1);
         }
         break;
   }
   mwpccyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::mwspView(bool rdraw)
{
   switch (mwspcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<12;i++) {mwpcelectnodes[i]->SetVisibility(0); mwpcelectnodes[i]->VisibleDaughters(0);}
         for (int i=0;i<4;i++) {mwpcbeamnodes[i]->SetVisibility(0); mwpcbeamnodes[i]->VisibleDaughters(0);}
         for (int i=0;i<6;i++) {mwpcframenodes[i]->SetVisibility(0); mwpcframenodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackerss to be visible
         for (int i=0;i<12;i++) {mwpcelectnodes[i]->SetVisibility(1); mwpcelectnodes[i]->VisibleDaughters(1);}
         for (int i=0;i<4;i++) {mwpcbeamnodes[i]->SetVisibility(1); mwpcbeamnodes[i]->VisibleDaughters(1);}
         for (int i=0;i<6;i++) {mwpcframenodes[i]->SetVisibility(1); mwpcframenodes[i]->VisibleDaughters(1);}
         break;
   }
   mwspcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::mwpcwinView(bool rdraw)
{
   switch (mwpcwincyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<24;i++) {mwpcwinnodes[i]->SetVisibility(0); mwpcwinnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackerss to be visible
         for (int i=0;i<24;i++) {mwpcwinnodes[i]->SetVisibility(1); mwpcwinnodes[i]->VisibleDaughters(1);}
         break;
   }
   mwpcwincyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::pbglassView(bool rdraw)
{
   switch (pbglasscyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<6;i++) {pbglassnodes[i]->SetVisibility(0);}
         break;
      case 1:
         // Set the trackers to be visible
         for (int i=0;i<6;i++) {pbglassnodes[i]->SetVisibility(1);}
         break;
   }
   pbglasscyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::sipmView(bool rdraw)
{
   switch (sipmcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<4;i++) {sipmnodes[i]->SetVisibility(0); sipmnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackers to be visible
         for (int i=0;i<4;i++) {sipmnodes[i]->SetVisibility(1); sipmnodes[i]->VisibleDaughters(1);}
         break;
   }
   sipmcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }


}

void Det_EventDisplay::pbblockView(bool rdraw)
{
   switch (pbblockcyc%2)
   {
      case 0:
         // Set the blocks to be invisible
         for (int i=0;i<8;i++) {pbblocknodes[i]->SetVisibility(0);}
         break;
      case 1:
         // Set the blocks to be visible
         for (int i=0;i<8;i++) {pbblocknodes[i]->SetVisibility(1);}
         break;
   }
   pbblockcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::symbView(bool rdraw)
{
   switch (symbcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0; i<nsymb; i++) {symbnodes[i]->SetVisibility(0); symbnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackers to be visible
         for (int i=0; i<nsymb; i++) {symbnodes[i]->SetVisibility(1); symbnodes[i]->VisibleDaughters(1);}
         break;
   }
   symbcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

void Det_EventDisplay::beamView(bool rdraw)
{
   switch (beamcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<18;i++) {beamnodes[i]->SetVisibility(0); beamnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackerss to be visible
         for (int i=0;i<18;i++) {beamnodes[i]->SetVisibility(1); beamnodes[i]->VisibleDaughters(1);}
         break;
   }
   beamcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }
}

void Det_EventDisplay::lmView(bool rdraw)
{
   switch (lmcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<6;i++) {lmnodes[i]->SetVisibility(0); lmnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackerss to be visible
         for (int i=0;i<6;i++) {lmnodes[i]->SetVisibility(1); lmnodes[i]->VisibleDaughters(1);}
         break;
   }
   lmcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

/* void Det_EventDisplay::gtView(bool rdraw)
{
   switch (gtcyc%2)
   {
      case 0:
         // Set the trackers to be invisible
         for (int i=0;i<2;i++) {gtnodes[i]->SetVisibility(0); gtnodes[i]->VisibleDaughters(0);}
         break;
      case 1:
         // Set the trackerss to be visible
         //cout<<"\n\n\n"<<"***WARNING***WARNING***WARNING***"<<"\n"<<"THE GEM TRACKERS DON'T EXIST"<<"\n\n\n";;
         for (int i=0;i<2;i++) {gtnodes[i]->SetVisibility(1); gtnodes[i]->VisibleDaughters(1);}
         break;
   }
   gtcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }


}*/

// Wire chamber sense planes
void Det_EventDisplay::wcspView(bool rdraw)
{
   switch (wcspcyc%2)
   {
      case 0:
         // Set the planes to be invisible
         for (int i=0;i<6;i++) {wcspnodes[i]->SetVisibility(0);}
         break;
      case 1:
         // Set the planes to be visible
         for (int i=0;i<6;i++) {wcspnodes[i]->SetVisibility(1); wcspnodes[i]->GetVolume()->SetLineColor(kGray+1);}
         break;
   }
   wcspcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }


}

// Dead wires (list set at beginning of file, will eventually be implemented )
void Det_EventDisplay::deadWiresBut(bool rdraw)
{

   if (deadw)
   {
      for (int u=0; u < 954; u++)
      {
         if (deadID[u]) {wirenodes[u]->SetVisibility(0);}
      }
      deadw = false;
   }
   else
   {
      for (int u=0; u < 954; u++)
      {
         if (deadID[u])
         {
            wirenodes[u]->SetVisibility(1);
            wirenodes[u]->GetVolume()->SetLineColor(14);
         }
      }
      deadw = true;
   }

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

// Function for the persistency of the dead wires across events
void Det_EventDisplay::deadWires()
{
   if (deadw)
   {
      for (int u=0; u < 954; u++)
      {
         if (deadID[u])
         {
            wirenodes[u]->SetVisibility(1);
            wirenodes[u]->GetVolume()->SetLineColor(14);
         }
      }
   }
   else
   {
      for (int u=0; u < 954; u++)
      {
         if (deadID[u]) {wirenodes[u]->SetVisibility(0);}
      }
   }  
}

// Cycle the toroid visibility (bottom only, all on, all off)
void Det_EventDisplay::toroidView(bool rdraw)
{

   switch (torcyc%4)
   {
      case 0:
         // Set the bottom coils to be visible
         for (int i=4;i<8;i++) {tornodes[i]->SetVisibility(1);}
         break;
      case 1:
         // Set the all coils to be visible, with very transparent top coils
         for (int i=0;i<8;i++)
			{
				tornodes[i]->SetVisibility(1);
				if (i<4)
				{
					tornodes[i]->GetVolume()->SetTransparency(75);
				}
			}
         break;
      case 2:
         // Set all coils to be visible
         for (int i=0;i<8;i++)
			{
				tornodes[i]->SetVisibility(1);
				tornodes[i]->GetVolume()->SetTransparency(0);
			}
         break;
      case 3:
         // Set all invisible
         for (int i=0;i<8;i++) {tornodes[i]->SetVisibility(0);}
         break;     
   }
   torcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }



}

// Cycle the target visibility (everything, just cell, nothing)
void Det_EventDisplay::targetView(bool rdraw)
{

   switch (targetcyc%4)
   {
      case 0:
         // Set all invisible
         for (int i=0;i<8;i++) {if (targetnodes[i]) targetnodes[i]->SetVisibility(0);}
         break; 
      case 1:
         // Set all elements to be visible
         for (int i=0;i<8;i++) {if (targetnodes[i]) targetnodes[i]->SetVisibility(1);}
         break;
      case 2:
         // Set just the cell, collimator, and WS to be visible
         for (int i=0;i<8;i++) {if (targetnodes[i]) targetnodes[i]->SetVisibility(0);}
			if (targetnodes[0]) targetnodes[0]->SetVisibility(1);
         if (targetnodes[1]) targetnodes[1]->SetVisibility(1);
			if (targetnodes[7]) targetnodes[7]->SetVisibility(1);
         break;
      case 3:
         // Set just the cell to be visible
         for (int i=0;i<8;i++) {if (targetnodes[i]) targetnodes[i]->SetVisibility(0);}
         if (targetnodes[1]) targetnodes[1]->SetVisibility(1);
         break;
   }
   targetcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

// Cycle the ToF visibility (left, right, all on, all off)
void Det_EventDisplay::tofView(bool rdraw)
{

   switch (tofcyc%2)
   {
      case 0:
         // Set all ToFs to be invisible
         for (int i=0;i<36;i++) {tofnodes[i]->SetVisibility(0);}
         break;
      case 1:
         // Set all ToFs to be visible
         for (int i=0;i<36;i++) {tofnodes[i]->SetVisibility(1);}
         break;
      /*case 2:
         // Set left ToFs visible
         for (int i=18;i<36;i++) {tofnodes[i]->SetVisibility(0);}
         break;
      case 3:
         // Set right ToFs visible
         for (int i=18;i<36;i++) {tofnodes[i]->SetVisibility(1);}
         for (int i=0;i<18;i++) {tofnodes[i]->SetVisibility(0);}
         break;*/
      
   }
   tofcyc++;

   if (rdraw)
   {
      gEve->GetEventScene()->Changed();
      gEve->FullRedraw3D(kFALSE);
   }

}

Long_t Det_EventDisplay::finalize()
{

   // Delete all the pointers
   delete etn;
   delete saveButton;
   delete saveName;
   delete torButton;
   delete symbButton;
   // delete gtButton;
   delete beamButton;
   delete tofButton;
   delete targetButton;
   delete tdButton;
   delete camButton;
   //delete camMenu;
   //delete wcMenu;
   delete wcButton;
   //delete perspButton;
   //delete orthotopButton;
   //delete orthorearButton;

   return 0;
};

Long_t Det_EventDisplay::viscoEvent()
{
  return 0;
};

// Event functions

Long_t Det_EventDisplay::prepare()
{
   //gEve->RemoveElement(trackp);
   return 0;
};

Long_t Det_EventDisplay::execute()
{

	// cout<<"\n";

   // Clear orphaned visualization objects
   gEve->ClearOrphanage();

	// Suppress the ROOT warnings that are meaningless
	gErrorIgnoreLevel = 2001;

   // Assert the current status of using lumi only mode to other plugins
   LumiOnlyCom();

   // Get the event information
   getBranchObject("EventInfo",(TObject**)&eventinfo);
   //eventinfo = (ORTEventInfo*)getFileObject("EventInfo");
   //int evn = eventinfo->eventNumber;
   int trig = eventinfo->trigFired;
   
   // Update the image filename
   char fnbuf[100];
   sprintf(fnbuf,"run_%u_event_%d.%s",(runinfo->runNumber),(eventinfo->eventNumber)-1,ext);
   saveName->SetText(fnbuf);

   // Parse the slow control data
   char bufa[100]; char bufb[100]; char bufc[100]; char bufd[100]; char bufe[100];
   char buff[100];
   if (!mccur)
      sprintf(bufa,"Toroid Current:  %+.1f A\n",dtorcur->value);
   else
      sprintf(bufa,"Toroid Current:  %+.1f A\n",*mccur);
   sprintf(bufb,"Target Flow:  %.3f sccm\n",dtarflow->value);
   sprintf(bufc,"Beam Current:  %.2f mA\n",dbeamcur->value);
   sprintf(bufd,"Beam Energy:  %.2f GeV\n",dbeamE->value);
   if (abs(dbeamcur->value)<0.1) sprintf(bufe,"Beam Species: Cosmic Rays");
   else if ((dbeamE->value)>0) sprintf(bufe,"Beam Species: Positrons");
   else
   {
      sprintf(bufe,"Beam Species: Electrons");
      sprintf(bufd,"Beam Energy:  %.2f GeV\n",-(dbeamE->value));
   }


   // Put together the trigger fired information
   // t# trigger pattern bit reference
   /* 0) kinematic trigger (0x1)
      1) Left SiPM and right TOF (0x2)
      2) Right SiPM & left TOF (0x4)
      3) Left Leadglass and right TOF (0x8)
      4) Right leadglass and left TOF (0x10)
      5) Left SiPM (0x20)
      6) Right SiPM (0x40)
      7) Left TOFs (0x80)
      8) Right TOFs (0x100)
      9) Left TOF OR & Right TOF coinc (0x200)
      10) Right TOF OR & Left TOF coinc (0x400) */
   bool trigbits[11];
   for (int t=0; t<11; t++)
   {
      if (trig & (int)pow(2,(double)t)) {trigbits[t] = true;}// cout<<t<<" fired\n";}
      else trigbits[t] = false;
   }

   // Update the run info display
   toroidCurrent->SetText(bufa);
   beamCurrent->SetText(bufc);
   beamCurrent->SetTextJustify(kTextLeft | kTextTop);
   targetFlow->SetText(bufb);
   beamEnergy->SetText(bufd);
   species->SetText(bufe);
   // trigger->SetText(buff);
   // trigger->Layout();
   toroidCurrent->Layout();
   beamCurrent->Layout();
   targetFlow->Layout();
   beamEnergy->Layout();
   species->Layout();
   vframe->Layout();

   deadWires(); // Maintain persistency of dead wire elements (on or off)

   // Collect lists of detector elements that are "hit"
   for (int w = 0; w < 954; w++)
   {
      if (wirenodes[w]->GetColour() == kRed) {wirestate[w] = true;}
      else {wirestate[w] = false;}

      if (wirenodes[w]->GetColour() == kViolet) {atwire[w] = true;}
      else {atwire[w] = false;}

      if (!awire) {wirenodes[w]->SetVisibility(0);} // Wire view persistence
   }

   // Set the anytime wires to be visible if desired
   if (noat)
   {
      for (int w = 0; w < 954; w++)
      {
         if (atwire[w]) {wirenodes[w]->SetVisibility(0);}
      }
   }

   // Do wire chamber hit recoloring when desired
   if (multi)
   {
      for (int c = 0; c<318; c++)
      {
         if (((wirestate[c*3] && wirestate[c*3+1]) || (wirestate[c*3] && wirestate[c*3+2])) || (wirestate[c*3+2] && wirestate[c*3+1]))
         {
            wirenodes[c*3]->GetVolume()->SetLineColor(kBlue);
            wirenodes[c*3+1]->GetVolume()->SetLineColor(kBlue);
            wirenodes[c*3+2]->GetVolume()->SetLineColor(kBlue);
         }

         if ((wirestate[c*3] && wirestate[c*3+1]) && (wirestate[c*3+2]))
         {
            wirenodes[c*3]->GetVolume()->SetLineColor(kGreen+1);
            wirenodes[c*3+1]->GetVolume()->SetLineColor(kGreen+1);
            wirenodes[c*3+2]->GetVolume()->SetLineColor(kGreen+1);
         }
         
      }
   }

   // Extract and draw any Kalman Filter tracks if present
   if (kf) drawKFTracks();

   // Draw any lumi hits if present (function does check)
   drawLumiHits();

   /*cout<<"\n"<<(((0x2 | 0x4) & trig))<<"\n";

   bool dontdraw= false;

   // If lumi only mode skip if the trigger isn't right show in display
   if ((((0x2 | 0x4) & trig) == 0) && lumionly)
   {
      cout<<"\n"<<"In the check"<<"\n";
      dontdraw = true;
   }
   
   cout<<"\n"<<"After the check"<<"\n";*/

   //if (!dontdraw)
   //{
      gEve->GetEventScene()->Changed(); // Get all changed scene elements
      gEve->FullRedraw3D(); // Redraw the scene with changes
   //}

   //std::cout<<"End of event.\n\n";   

	// Reset the ROOT warning level
	// gErrorIgnoreLevel = 0;

   return 0;
};

// Function to draw tracks produced by the Kalman Filter
void Det_EventDisplay::drawKFTracks()
{
   // Clear any previously drawn tracks
   for (std::vector<TEveLine* >::iterator iter=kflines.begin(); iter != kflines.end(); iter++)
	{
	   (*iter)->Reset();
	   delete (*iter);
	}
   kflines.clear();
   
   int test = 0;
   // Create enough point sets for all the tracks
   for (unsigned int n = 0; n < ((kftracks->tracks).size()); n++)
	{
	   kflines.push_back (new TEveLine());
      test++;
	}
   //std::cout<<"There might be "<<test<<" KF tracks!\n";

   int count = 0;
   // Loop through the Kalman Filter tracks for this event
   for (std::vector<KalmanFilterTrack>::iterator iter=(kftracks->tracks).begin(); iter != (kftracks->tracks).end(); iter++)
   {
      // Loop over the position vectors listed in the trace
      for (std::vector<TVector3>::iterator posi=(*iter).trace.position.begin(); posi != ((*iter).trace.position).end(); posi++)
      {
         kflines[count]->SetNextPoint((*posi).x()/10,(*posi).y()/10,(*posi).z()/10);
      }
      kflines[count]->SetMainColor(kGreen);
      kflines[count]->SetSmooth(true);
      count++;
   }
   //std::cout<<"There are "<<count<<" KF tracks!\n";

   // Add any filled lines to the display
   TEveElementList * holder = new TEveElementList();
   for (std::vector<TEveLine* >::reverse_iterator iter=kflines.rbegin();iter!=kflines.rend();iter++)
	{
	   holder->AddElement(*iter);
	}
   gEve->AddElement(holder);

   delete holder;
}

void Det_EventDisplay::drawLumiHits()
{
   // Clear the point set and lines from the previous event
   for (std::vector<TEveLine* >::iterator iter=llines.begin(); iter != llines.end(); iter++)
   {
      (*iter)->Reset();
      delete (*iter);
   }
   llines.clear();
   lpset->Reset();
   lpset->SetMarkerColor(kCyan);
   lpset->SetMarkerSize(1.5);

   if (lumi)
   {
      // Extract the hits from the LumiGEM output
      std::vector <LumiGEMhit> gemhits = lumi->hits;

      int odc = 0;  // 1D hit counter

      // Loop through the hits
      for (std::vector<LumiGEMhit>::iterator iter = gemhits.begin(); iter != gemhits.end(); iter++)
      {
         // Which hit is the gem in?
         int gid = (*iter).GEMid;

         // Check to see if the hit is 1D or 2D
         if (((int)((*iter).xl) == -1) || ((int)((*iter).yl) == -1)) // Hit is 1D
         {
            double locs[3]; double loce[3];  // GEM local coordinates

            // If yl = -1, it's a vertical hit
            if ((int)((*iter).yl) == -1)
            {
               // Set the end points of the line that will symbolize the hit in the
               // local GEM frame
               locs[0] = (*iter).xl*0.04-5.0; locs[1] =  5.0; locs[2] = 0;
               loce[0] = (*iter).xl*0.04-5.0; loce[1] = -5.0; loce[2] = 0;
            }
            else
            {
               // Set the end points of the line that will symbolize the hit in the
               // local GEM frame
               locs[0] =  5.0; locs[1] = (*iter).yl*0.04-5.0; locs[2] = 0;
               loce[0] = -5.0; loce[1] = (*iter).yl*0.04-5.0; loce[2] = 0;
            }

            // Convert to the master frame
            double mass[3]; double mase[3];
            lmnodes[gid]->LocalToMaster(locs,mass);
            lmnodes[gid]->LocalToMaster(loce,mase);

            // Make a new line to show the hit
            llines.push_back (new TEveLine());
            llines[odc]->SetNextPoint(mass[0],mass[1],mass[2]);
            llines[odc]->SetNextPoint(mase[0],mase[1],mase[2]);
            llines[odc]->SetMainColor(kCyan+3);
            llines[odc]->SetLineWidth(2);
            gEve->AddElement(llines[odc]);  // Add to the display
            odc++;  // Add to the 1D hit count
            
         }
         else  // Hit is 2D
         {
            // Coordinates in the local hit frame accounting for the numbering of
            // of the pitch coordinates xl and yl
            double loc[3] = {(*iter).xl*0.04-5.0, (*iter).yl*0.04-5.0, 0};
            // double loco[3] = {0,0,0};

            //cout<<"\nPitch: "<<(*iter).xl<<" "<<(*iter).yl<<"  Calc: "<<((*iter).xl*0.04-0.5)<<" "<<((*iter).yl*0.04-0.5)<<"\n";
            
            // Convert to global coordinates for the proper GEM
            double glob[3]; // double globo[3];
            lmnodes[gid]->LocalToMaster(loc,glob);
            // lmnodes[gid]->LocalToMaster(loco,globo);

            // Add point to the point set
            lpset->SetNextPoint(glob[0],glob[1],glob[2]);
            // lpset->SetNextPoint(globo[0],globo[1],globo[2]);
         }
      }
   }

   if (mwpc)
   {
      // Extract the hits from the LumiGEM output
      std::vector <MWPC_cluster> mwpchits = mwpc->hits;

      // Loop through the hits
      for (std::vector<MWPC_cluster>::iterator iter = mwpchits.begin(); iter != mwpchits.end(); iter++)
      {
         // Which hit is the gem in?
         int gid = (*iter).cham;

         // Coordinates in the local hit frame accounting for the numbering of
         // of the pitch coordinates xl and yl
         double loc[3] = {(*iter).xx/10, (*iter).yy/10, 0};
         // double loco[3] = {0,0,0};

         //cout<<"\nPitch: "<<(*iter).xl<<" "<<(*iter).yl<<"  Calc: "<<((*iter).xl*0.04-0.5)<<" "<<((*iter).yl*0.04-0.5)<<"\n";
         
         // Convert to global coordinates for the proper GEM
         double glob[3]; // double globo[3];
         mwpcnodes[gid]->LocalToMaster(loc,glob);
         // lmnodes[gid]->LocalToMaster(loco,globo);

         // Add point to the point set
         lpset->SetNextPoint(glob[0],glob[1],glob[2]);
         // lpset->SetNextPoint(globo[0],globo[1],globo[2]);
      }
   }

   // Make sure the point set is added to the EVE system
   gEve->AddElement(lpset);

};

void Det_EventDisplay::localToGlobalWC(unsigned int id, double * local, double * global)
{
   // Set all the appropriate index numbers for the given ID
   int sectn = 0;
   if (id>476) sectn = 1;

   int sln = 0;
   if (id>53)  sln = 1;
   if (id>110) sln = 2;
   if (id>188) sln = 3;
   if (id>269) sln = 4;
   if (id>371) sln = 5;
   if (id>476) sln = 6;
   if (id>530) sln = 7;
   if (id>587) sln = 8;
   if (id>665) sln = 9;
   if (id>746) sln = 10;
   if (id>848) sln = 11;

   int framen = sln/2;
   int layern = id%3 + sln*3;
   
   // Run through the heirarchy to set the output
   wirenodes[id]->LocalToMaster(local,global);
   wclayernodes[layern]->LocalToMaster(global,local);
   wcslnodes[sln]->LocalToMaster(local,global);
   wcchamnodes[framen]->LocalToMaster(global,local);
   wcgasnodes[framen]->LocalToMaster(local,global);
   wcframenodes[framen]->LocalToMaster(global,local);
   wcsecnodes[sectn]->LocalToMaster(local,global);

};

Long_t Det_EventDisplay::findnodes()
{

   int all = 0;  // Counter to keep track of if the whole geometry is found
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // ToFs   
   int ind = 0;
   TObjArray *l=gGeoManager->GetTopVolume()->GetNodes();

   for (int i=0;i<l->GetEntriesFast() && ind<36;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("TF_log")) 
         tofnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);//->GetVolume();

   // Check for all ToFs
   if (ind < 36) {all = all + 0x1;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
	
	// Lead glass calorimeters
	pbglass = false; // Old geometries have no lead glass
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<6;i++)
	{
      if (TString((*l)[i]->GetName()).BeginsWith("LG")) 
		{
      	pbglassnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
			pbglass = true;
		}
	}

	// Lead block shielding
	pbblock = false; // Old geometries have no lead shielding
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<8;i++)
	{
      if (TString((*l)[i]->GetName()).BeginsWith("Lead_Block")) 
		{
      	pbblocknodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
			pbblock = true;
		}
	}

   // Toroid Coils
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<8;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("Toroid")) 
      {
         tornodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);//->GetVolume();
      }
   // Check for all coils
   if (ind < 8) {all = all + 0x2;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Wire Chamber Components
   std::cout<<"Loading the drift chamber geometry...\n";
   int wcs = 0; // Counter for the sectors as they are found
   int wcf= 0; // Counter for the frames as they are found
	int wcw = 0; // Counter for the windows as they are found
   int wcg = 0; // Counter for the gas volumes as they are found
   int ch = 0; // Counter for the chambers as they are found
   int sl = 0; // Counter for the superlayers as they are found
   int lay = 0; // Counter for the layers as they are found
   int cells = 0; // Counter for the cells as they are found
   int wires = 0; // Counter for wires as they are found
   int sp = 0; // Counter for sensitive planes as they are found
   TObjArray *secs[2]; // Node list holders
   TObjArray *gas[6];
	TObjArray *frames[6];
   TObjArray *cham[6];
   TObjArray *superlayer[12];
   TObjArray *layer[36];
   TObjArray *clist[318];

   // Now, extract the sector volumes from this list
   for (int i=0; i<(l->GetEntriesFast()) && (wcs<2) ; i++)
   { 

      if ((TString((*l)[i]->GetName()).BeginsWith("WC")))
      {
         wcsecnodes[wcs] = ( (TGeoNodeMatrix *) (*l)[i]);
	 		wcsecnodes[wcs]->SetVisibility(0);
         secs[wcs++] = ( (TGeoNodeMatrix *) (*l)[i])->GetNodes();
      }
   }
   // Check for both sectors
   if (wcs < 2) {all = all + 0x4; return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Get the gas, window (including light blocker, etc.), and frame volumes for each sector
   for (int g=0; g<2; g++)
   {
      for (int k=0; ((k<(secs[g]->GetEntriesFast())) && (((wcg<2) || (wcf<6)) || (wcw<15))) ; k++)
      {
         if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_gas")))
         {
            wcgasnodes[wcg] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
            wcgasnodes[wcg]->SetVisibility(0);
            gas[wcg++] = ( (TGeoNodeMatrix *) (*secs[g])[k])->GetNodes();
         }
         if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_frame")))
         {
            wcframenodes[wcf] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
            frames[wcf++] = ( (TGeoNodeMatrix *) (*secs[g])[k])->GetNodes();
         }
			if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_window")))
			{
				wcwinnodes[wcw++] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
			}
			if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_lightguard")))
			{
				wcwinnodes[wcw++] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
			}
			if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_shield")))
			{
				wcwinnodes[wcw++] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
			}
			if ((TString((*secs[g])[k]->GetName()).BeginsWith("WC_foil")))
			{
				wcwinnodes[wcw++] = ( (TGeoNodeMatrix *) (*secs[g])[k]);
			}
      }
   }

	// for (int j=0; j<6; j++) cout<<(wcframenodes[j]->GetVolume()->GetName())<<"\n";
	// for (int j=0; j<14; j++) cout<<(wcwinnodes[j]->GetVolume()->GetName())<<"\n";
	// abort();

   // Check for both gas volumes and frames, and windows
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   if (wcf < 2) {all = all + 0x10;} // return (Long_t)all;}
	if (wcw < 10) {all = all + 0x100000;}; // cout<<wcw; return (Long_t)all;} // return (Long_t)all;}
	nframes = wcf;
	nwin = wcw;
   // std::cout<<wcw<<"\n\n";

	// Handle new frame hierarchy when present
   if (wcg == 0 && nframes == 6)
	{
		cout<<"\nSurvey2013 wire chamber geometry identified.  Processing...\n\n";
		for (int g=0; g<6; g++)
		{
			// cout<<(frames[g]->GetEntriesFast())<<" ";
			for (int k=0; k<(frames[g]->GetEntriesFast()) && (wcg<6) ; k++)
			{
				if ((TString((*frames[g])[k]->GetName()).BeginsWith("WC_gas")))
         	{
            	wcgasnodes[wcg] = ( (TGeoNodeMatrix *) (*frames[g])[k]);
            	wcgasnodes[wcg]->SetVisibility(0);
					// cout<<wcgasnodes[wcg]->CountDaughters()<<"\n";
            	gas[wcg++] = ( (TGeoNodeMatrix *) (*frames[g])[k])->GetNodes();
					// gas[wcg] = wcgasnodes[wcg]->GetNodes();
					// wcg++;
         	}
			}
		}		
	}
	else if (wcg != 2) {all = all + 0x8; return (Long_t)all;}

   // std::cout<<"I found "<<wcg<<" wire chamber gas volumes!\n";

   // Get the chamber volumes for each sector
   for (int s=0; s<nframes; s++)
   {
   	for (int k=0; k<(gas[s]->GetEntriesFast()) && (ch<6) ; k++)
   	{
         if ((TString((*gas[s])[k]->GetName()).BeginsWith("WC_chamber")))
         {
            wcchamnodes[ch] = ( (TGeoNodeMatrix *) (*gas[s])[k]);
            wcchamnodes[ch]->SetVisibility(0);
            cham[ch++] = ( (TGeoNodeMatrix *) (*gas[s])[k])->GetNodes();
         }
   	}
   }
   // Check for all chambers
   if (ch < 6) {all = all + 0x20; return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   std::cout<<"I found "<<ch<<" chamber volumes!\n";

   // Get the superlayer volumes for each sector
   for (int g=0; g<6; g++)
   {
      for (int k=0; k<(cham[g]->GetEntriesFast()) && (sl<12) ; k++)
      {
            if ((TString((*cham[g])[k]->GetName()).BeginsWith("WC_superlayer")))
            {
               wcslnodes[sl] = ( (TGeoNodeMatrix *) (*cham[g])[k]);
               wcslnodes[sl]->SetVisibility(0);
               superlayer[sl++] = ( (TGeoNodeMatrix *) (*cham[g])[k])->GetNodes();
            }
      }
   }
   // Check for all superlayers
   if (sl < 12) {all = all + 0x40; return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   std::cout<<"I found "<<sl<<" superlayer volumes!\n";

   // Extract the cell volumes from each superlayer
   for (int g=0; g<12; g++)
   {
      for (int j=0; j<(superlayer[g]->GetEntriesFast()) && (cells<318); j++)
      {
            if ((TString((*superlayer[g])[j]->GetName()).BeginsWith("WC_cell")))
            {
               wccellnodes[cells] = ( (TGeoNodeMatrix *) (*superlayer[g])[j]);
               wccellnodes[cells]->SetVisibility(0);
               clist[cells++] = ( (TGeoNodeMatrix *) (*superlayer[g])[j])->GetNodes();
            }
      }
   }

   bool nocells = false;
   std::cout<<"I found "<<cells<<" cells!\n";
   if (cells==0) {
      std::cout<<"\nThat's OK, you must be using a non-cells geometry.\n";
      std::cout<<"I'll try a different approach.\n\n";
      nocells = true;
   }

   // In the case that the cells were found get the wires that way
   if (!nocells)
   {
      // Extract the wires from each cell
      for (int c=0; c<318; c++)
      {
           for (int h=0; h<(clist[c]->GetEntriesFast()) && (wires<954); h++)
           {
              if ((TString((*clist[c])[h]->GetName()).BeginsWith("WC_wire0")))
              {
                  wirenodes[wires++]=( (TGeoNodeMatrix *) (*clist[c])[h]);
              }
              else if ((TString((*clist[c])[h]->GetName()).BeginsWith("WC_wire1")))
              {
                  wirenodes[wires++]=( (TGeoNodeMatrix *) (*clist[c])[h]);
                                 
              }
              else if ((TString((*clist[c])[h]->GetName()).BeginsWith("WC_wire2")))
              {
                  wirenodes[wires++]=( (TGeoNodeMatrix *) (*clist[c])[h]);
                                 
              }
           }
     }
      std::cout<<"I found "<<wires<<" wires!\n";
   }
   else
   {
      // Extract layers from the superlayers and then wires from layers instead
      for (int g=0; g<12; g++)
      {
         for (int j=0; j<(superlayer[g]->GetEntriesFast()) && (lay<36); j++)
         {
            if ((TString((*superlayer[g])[j]->GetName()).BeginsWith("WC_layer")))
            {
               wclayernodes[lay] = ( (TGeoNodeMatrix *) (*superlayer[g])[j]);
               wclayernodes[lay]->SetVisibility(0);
               layer[lay++] = ( (TGeoNodeMatrix *) (*superlayer[g])[j])->GetNodes();
            }
         }
      }

      std::cout<<"I found "<<lay<<" layers!\n";
      int sls[12];
      sls[0] = 0; sls[1] = 54; sls[2] = 111; sls[3] = 189;
      sls[4] = 270; sls[5] = 372; sls[6] = 477; sls[7] = 531;      
      sls[8] = 588; sls[9] = 666; sls[10] = 747; sls[11] = 849;
      for (int y=0; y<36; y++)
      {
         // Because of the index scheme relative to this geometry, some bojangling
         // has to be done to work around to get the right indices
         int lcount = 0;
         int loff = y%3;
         int tind = floor(y/3);
         int slstart = sls[tind];
         for (int j=0; j<(layer[y]->GetEntriesFast()) && (wires<954); j++)
         {
            if ((TString((*layer[y])[j]->GetName()).BeginsWith("WC_wire")))
            {
               wirenodes[slstart+loff+3*lcount] = ( (TGeoNodeMatrix *) (*layer[y])[j]);
               wires++;
               lcount++;
            }
         }
      }

      std::cout<<"I found "<<wires<<" wires!\n";
    
   } // End of no cells case

   // Check for all wires
   if (wires < 954) {all = all + 0x80;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   
   // Extract the sensitive planes from each sector
   for (int g=0; g<6; g++)
   {
      for (int j=0; j<(cham[g]->GetEntriesFast()) && (sp<6); j++)
      {
         if ((TString((*cham[g])[j]->GetName()).BeginsWith("WC_plane")))
         {
            wcspnodes[sp++] = ( (TGeoNodeMatrix *) (*cham[g])[j]);
         }
      }      
   }
   // Check for all SD planes
   if (sp<6) {all = all + 0x100;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   std::cout<<"I found "<<sp<<" SD planes!\n";

   // Gem Trackers (REMOVED BECAUSE THEY DON'T EXIST)
   /* ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<2;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("GT")) 
      {
         gtnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
      }
   // Check for all GTs
   if (ind < 2) {all = all + 0x200;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";*/

   // Lumi Gems
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<6;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("LM_log")) 
      {
         lmnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
      }
   // Check for all LumiGEMs
   if (ind < 6) {all = all + 0x400;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Beam line
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<18;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("BP")) 
      {
         beamnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
      }
   // Check for all beam line components
   if (ind < 18) {all = all + 0x800;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // SYMB
   ind = 0;
	TObjArray * ssy[2];
	TGeoNode * sleft, * sright;
   for (int i=0;i<l->GetEntriesFast() && ind<20;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("SM")) 
      {
         symbnodes[ind]=( (TGeoNodeMatrix *) (*l)[i]);
			ssy[ind++] = ( (TGeoNodeMatrix *) (*l)[i])->GetNodes();
      }
   nsymb = ind;
   // Check for all SYMB components
   if ((nsymb != 2) && (nsymb != 20)) {all = all + 0x1000;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   //cout<<"\n\n"<<nsymb<<"\n\n";

//	for (int m=0; m<2; m++)
//		for (int k=0; k<(ssy[m]->GetEntriesFast()); k++)
//		{
//		   if ((TString((*ssy[m])[k]->GetName()).BeginsWith("SM_log0")))
//		   {
//		      sleft = ( (TGeoNodeMatrix *) (*ssy[m])[k]);
//				cout<<"Found left!\n";
//		   }
//		   if ((TString((*ssy[m])[k]->GetName()).BeginsWith("SM_log10")))
//		   {
//		      sright = ( (TGeoNodeMatrix *) (*ssy[m])[k]);
//				cout<<"Found right!\n";
//		   }
//		}

//	for (int j=0; j<500; j++)
//	{
//		double theta = j/(2*3.1415);
//		double lx = 1.0225*cos(theta) - 5.57575;
//		double ly = 1.0225*sin(theta);

//		double pointa[3] = {lx,ly,-4.91};
//		double pointb[3];

//		sleft->LocalToMaster(pointa,pointb);
//		symbnodes[0]->LocalToMaster(pointb,pointa);
//		cout<<pointa[0]<<" "<<pointa[1]<<" "<<pointa[2]<<"\n";
//	}
//	for (int j=0; j<500; j++)
//	{
//		double theta = j/(2*3.1415);
//		double lx = 1.0225*cos(theta) - 5.57575;
//		double ly = 1.0225*sin(theta);

//		double pointa[3] = {lx,ly,4.88};
//		double pointb[3];

//		sleft->LocalToMaster(pointa,pointb);
//		symbnodes[0]->LocalToMaster(pointb,pointa);
//		cout<<pointa[0]<<" "<<pointa[1]<<" "<<pointa[2]<<"\n";
//	}
//	for (int j=0; j<500; j++)
//	{
//		double theta = j/(2*3.1415);
//		double lx = 1.0225*cos(theta) + 5.57575;
//		double ly = 1.0225*sin(theta);

//		double pointa[3] = {lx,ly,-4.91};
//		double pointb[3];

//		sright->LocalToMaster(pointa,pointb);
//		symbnodes[1]->LocalToMaster(pointb,pointa);
//		cout<<pointa[0]<<" "<<pointa[1]<<" "<<pointa[2]<<"\n";
//	}
//	for (int j=0; j<500; j++)
//	{
//		double theta = j/(2*3.1415);
//		double lx = 1.0225*cos(theta) + 5.57575;
//		double ly = 1.0225*sin(theta);

//		double pointa[3] = {lx,ly,4.88};
//		double pointb[3];

//		sright->LocalToMaster(pointa,pointb);
//		symbnodes[1]->LocalToMaster(pointb,pointa);
//		cout<<pointa[0]<<" "<<pointa[1]<<" "<<pointa[2]<<"\n";
//	}

//	abort();

	// K-Beam sections for lead glass (simply to make them invisible)
   for (int i=0;i<l->GetEntriesFast();i++)
      if (TString((*l)[i]->GetName()).BeginsWith("Kbeam")) 
      {
         ( (TGeoNodeMatrix *) (*l)[i])->SetVisibility(0);
      }

   // MWPCs
   int mwu = 0;
   TObjArray *mws[6];

   // Support beams
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<4;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("MW_Beam")) 
      {
         mwpcbeamnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
      }
   // Check for all MWPC support components
   if (ind < 4) {all = all + 0x2000;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Get main containers
   for (int i=0; i<(l->GetEntriesFast()) && (mwu<6) ; i++)
   { 
      if ((TString((*l)[i]->GetName()).BeginsWith("MW_log")))
      {

         mws[mwu++] = ( (TGeoNodeMatrix *) (*l)[i])->GetNodes();
         mwpcnodes[mwu-1] = ( (TGeoNodeMatrix *) (*l)[i]);
         mwpcnodes[mwu-1]->SetVisibility(0); mwpcnodes[mwu-1]->VisibleDaughters(1);
      }
   }
   // Check for all MWPCs
   if (mwu < 6) {all = all + 0x4000; return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Get the elements from the main containers
   int fn = 0;
   int gn = 0;
   int en = 0;
	int wn = 0;
	mwpcwin = false; // Old geometries have no MWPC windows
   for (int m=0; m<6; m++)
   {
      for (int k=0; k<(mws[m]->GetEntriesFast()); k++)
      {
         if ((TString((*mws[m])[k]->GetName()).BeginsWith("MW_Frame")))
         {
            mwpcframenodes[fn++] = ( (TGeoNodeMatrix *) (*mws[m])[k]);
         }

         if ((TString((*mws[m])[k]->GetName()).BeginsWith("MW_Gas")))
         {
            mwpcgasnodes[gn++] = ( (TGeoNodeMatrix *) (*mws[m])[k]);
         }

         if ((TString((*mws[m])[k]->GetName()).BeginsWith("MW_Elect")))
         {
            mwpcelectnodes[en++] = ( (TGeoNodeMatrix *) (*mws[m])[k]);
         }

			if ((TString((*mws[m])[k]->GetName()).BeginsWith("MW_Window")) || (TString((*mws[m])[k]->GetName()).BeginsWith("MW_Foil")))
			{
            mwpcwinnodes[wn++] = ( (TGeoNodeMatrix *) (*mws[m])[k]);
				mwpcwin = true; // We do have windows in this geometry
         }

      }
   }
   // Check for all MWPC components
   if (fn < 6) {all = all + 0x8000;} // return (Long_t)all;}  // Frames
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   if (gn < 6) {all = all + 0x10000;} // return (Long_t)all;}  // Gas volumes
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";
   if (en < 12) {all = all + 0x20000;} // return (Long_t)all;}  // Electronics Boxes
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Lumi ECals
   /*ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<40;i++)
      if (TString((*l)[i]->GetName()).BeginsWith("LC")) 
      {
         lcalnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
      }*/

   // SiPMs
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<4;i++)
	{
      if (TString((*l)[i]->GetName()).BeginsWith("SI")) 
      {
         sipmnodes[ind++] = ( (TGeoNodeMatrix *) (*l)[i]);
      }
	}
   // Check for all SiPMs
   if (ind < 4) {all = all + 0x40000;} // return (Long_t)all;}
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

	// Colton's exit box (may or may not be present, set to invisible when present)
   ind = 0;
   for (int i=0;i<l->GetEntriesFast() && ind<1;i++)
	{
      if (TString((*l)[i]->GetName()).BeginsWith("Exit")) 
      {
         exitnode = ( (TGeoNodeMatrix *) (*l)[i]);
			ind++;
      }
	}
	if (ind>0) exitnode->SetVisibility(0);

   // Target Chamber, etc.
   // Windows and window frames
	for (int j=0; j<8; j++) targetnodes[j] = NULL;
   ind = 2;
   int tcount = 0;  // Counter for the target objects for check
   for (int i=0;i<l->GetEntriesFast() && ind<6;i++)
   {
      if (TString((*l)[i]->GetName()).BeginsWith("Win")) 
      {
         targetnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
         tcount++;
      }
   }
   for (int i=0;i<l->GetEntriesFast() && ind<6;i++)
   {
      if (TString((*l)[i]->GetName()).BeginsWith("Fr")) 
      {
         targetnodes[ind++]=( (TGeoNodeMatrix *) (*l)[i]);
         tcount++;
      }
   }

	// Check for wakefield suppressor in world log (OLD HIERARCHY)
	bool fnd = false;
   for (int i=0;i<l->GetEntriesFast() && !fnd;i++)
   {
      if (TString((*l)[i]->GetName()).BeginsWith("WS")) 
      {
         targetnodes[7]=( (TGeoNodeMatrix *) (*l)[i]);
         tcount++;
      }
   }

   // Chamber, Cell, and Collimator
   TObjArray *TC = NULL;
   // Now, extract the target chamber volume
   bool tfound = false;
   for (int i=0; i<(l->GetEntriesFast()) && (!tfound) ; i++)
   {
      if ((TString((*l)[i]->GetName()).BeginsWith("TC")))
      {
         TC = ( (TGeoNodeMatrix *) (*l)[i])->GetNodes();
         ( (TGeoNodeMatrix *) (*l)[i])->SetVisibility(0);
         tfound = true;
      }
   }

   //int test = (TC->GetEntriesFast());

   //std::cout<<"Looking through "<<test<<" nodes\n";   

   for (int k=0; k<(TC->GetEntriesFast()); k++)
   {
   //std::cout<<"I am looking for the cell and coll.\n";
      if ((TString((*TC)[k]->GetName()).BeginsWith("Cell")))
      {
         targetnodes[1] = ( (TGeoNodeMatrix *) (*TC)[k]);
         //std::cout<<"I found the cell\n";
         tcount++;
      }
      else if ((TString((*TC)[k]->GetName()).BeginsWith("Coll")))
      {
         targetnodes[0] = ( (TGeoNodeMatrix *) (*TC)[k]);
         tcount++;
      }
      else if ((TString((*TC)[k]->GetName()).BeginsWith("WS")))
      {
         targetnodes[7] = ( (TGeoNodeMatrix *) (*TC)[k]);
         tcount++;
      }
      else if ((TString((*TC)[k]->GetName()).BeginsWith("TC_Frame")))
      {
         targetnodes[6] = ( (TGeoNodeMatrix *) (*TC)[k]);
         //std::cout<<"I found the target frame\n";
         tcount++;
      }
   }
   // Check for all target components
   if (tcount < 6) {all = all + 0x80000;} // return (Long_t)all;}  // This one
   //std::cout<<"\\n\\nall = "<<all<<"\\n\\n";

   // Return the bit count (0 if no failure)
   return (Long_t)all;

}

// Startup function to get the dead wire list from the WC2 plugin
Long_t Det_EventDisplay::deadList()
{
	for (int w=0; w<954; w++)
	{
		if (wirenodes[w]->GetColour() == 14) idDeadWire(w);
		else idLiveWire(w);
	}

	return 0;
};

// Function to set dead wire list
Long_t Det_EventDisplay::idDeadWire(int id)
{
   deadID[id] = true;
   return 0;
};

// Function to set a wire as live
Long_t Det_EventDisplay::idLiveWire(int id)
{
   deadID[id] = false;
   return 0;
};

// Functions to add a point (in initialization) to the visualization (give
// global coordinates in mm, can add local functionality if somebody needs it)
Long_t Det_EventDisplay::addVisPoint(int id,double x, double y, double z,double att)
{
	// Temporary visualization point
	visPoint temp;

	// cout<<"\n"<<id<<" "<<x<<" "<<y<<" "<<z;

	// Divide coordinates by 10 to go to cm (which is the EVE unit)
	temp.x = x/10;
	temp.y = y/10;
	temp.z = z/10;
	temp.att = att;

	// Add to the appropriate visualization set
	if (natts == 0)
	{
		// Create a temp set, fill it, and push it
		visSet tset;
		tset.att = att;
		tset.points.push_back(temp);
		viss.push_back(tset);
		natts++;
	}
	else
	{
		// Next check if this set already exists
		bool exists = false;
		int ind;
		for (int j = 0; j<natts; j++)
		{
			if (viss[j].att == att) {exists = true; ind = j;}
		}
		
		// Push on the set found or create a new set
		if (exists)
		{
			viss[ind].points.push_back(temp);
		}
		else
		{
			// Create a temp set, fill it, and push it
			// for (int q=0;q<100;q++) cout<<"IMMA CHARGIN MA LAZER!\n";
			visSet tset;
			tset.att = att;
			tset.points.push_back(temp);
			viss.push_back(tset);
			natts++;
		}
	}

	return 0;
};

void Det_EventDisplay::setPoints(std::vector<visPoint> visp,TEvePointSet * look)
{
	// Iterate through the vector of points and add them to the EVE point set
	for (std::vector<visPoint>::iterator iter=visp.begin(); iter != visp.end(); iter++)
	{
		// cout<<(*iter).x<<" "<<(*iter).y<<" "<<(*iter).z<<"\n";
		look->SetNextPoint((*iter).x,(*iter).y,(*iter).z);
	}
};

int Det_EventDisplay::pointColor(int catt)
{

	if (natts == 1) return kGreen+3;

	// Compute where on the full color scale we are
	int clevel =  catt;
	// cout<<catt<<" "<<natts<<" "<<clevel<<"\n\n";
	if (clevel == 0) return kRed+4;
	if (clevel == 1) return kRed+3;
	if (clevel == 2) return kRed+2;
	if (clevel == 3) return kRed+1;
	if (clevel == 4) return kRed;
	if (clevel == 5) return kRed-4;
	if (clevel == 6) return kRed-7;
	if (clevel == 7) return kRed-9;
	if (clevel == 8) return kRed-10;
	if (clevel == 9) return kBlue-10;
	if (clevel == 10) return kBlue-9;
	if (clevel == 11) return kBlue-7;
	if (clevel == 12) return kBlue-4;
	if (clevel == 13) return kBlue;
	if (clevel == 14) return kBlue+1;
	if (clevel == 15) return kBlue+2;
	if (clevel == 16) return kBlue+3;
	if (clevel == 17) return kBlue+4;

	return 0; // Point is white if offscale
}

// Methods required for cooker
Long_t Det_EventDisplay::cmdline(char * cmd)
{
  //add cmdline handling here

  return 0;
};

extern "C"{
  Plugin *factory(TTree *in, TTree *out, TFile *inf_, TFile *outf_, TObject *p)
  {
    return (Plugin*) new Det_EventDisplay(in,out,inf_,outf_,p);
  }
}

ClassImp(Det_EventDisplay);

// Constructor/destructor for point visualization container class

visSet::visSet()
{
	look = new TEvePointSet();
};

visSet::~visSet()
{
};

