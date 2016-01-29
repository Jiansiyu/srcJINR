/**
 * \file BmnTrackingQaReport.cxx
 * \author Semen Lebedev <s.lebedev@gsi.de> - Original author. First version of code for CBM experiment.
 * \author Sergey Merts <Sergey.Merts@gmail.com> - Modification for BMN experiment.
 * \date 2011-2014
 */
#include "BmnTrackingQaReport.h"
#include "report/BmnReportElement.h"
#include "report/BmnHistManager.h"
#include "report/BmnDrawHist.h"
#include "BmnUtils.h"
#include "TH1.h"
#include "map"
#include "TCanvas.h"
#include "TLine.h"
#include "BmnTrackingQa.h"
#include <boost/assign/list_of.hpp>
#include <vector>
#include <set>
#include "TString.h"

using lit::NumberToString;
using lit::FindAndReplace;
using lit::Split;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::set;
using std::endl;
using std::make_pair;
using std::pair;
using boost::assign::list_of;

string DefaultEfficiencyLabelFormatter(const string& histName, Float_t efficiency) {
    vector<string> split = Split(histName, '_');
    return split[1] + ":" + split[3] + "(" + NumberToString<Float_t > (efficiency, 1) + ")";
}

string ElectronIdEfficiencyLabelFormatter(const string& histName, Float_t efficiency) {
    vector<string> split = Split(histName, '_');
    return FindAndReplace(split[1], "Gem", "") + " (" + NumberToString<Float_t > (efficiency, 1) + ")";
}

string DefaultPionSuppressionLabelFormatter(const string& histName, Float_t efficiency) {
    vector<string> split = Split(histName, '_');
    return split[1] + " (" + NumberToString<Float_t > (efficiency, 1) + ")";
}

BmnTrackingQaReport::BmnTrackingQaReport() :
BmnSimulationReport(),
fGlobalTrackVariants() {
    SetReportName("tracking_qa");
}

BmnTrackingQaReport::BmnTrackingQaReport(vector<string> header) :
BmnSimulationReport(),
fGlobalTrackVariants(),
fHeader(header) {
    SetReportName("tracking_qa");
}

BmnTrackingQaReport::~BmnTrackingQaReport() {
}

void BmnTrackingQaReport::Create() {
    Out() << R()->DocumentBegin();
    Out() << R()->Title(0, GetTitle());
    Out() << PrintEventInfo();
    //Out() << PrintNofObjects();
    //Out() << PrintTrackHits();
    //Out() << PrintNofGhosts();
    Out().precision(3);
    Out() << PrintTrackingEfficiency(false);
    PrintCanvases();
    Out() << R()->DocumentEnd();
}

string BmnTrackingQaReport::PrintEventInfo() {
    Out() << "<h2>Event generator: QGSM</h2>" << endl;
    Out() << "<h2>Energy: 4 GeV/n</h2>" << endl;
    if (GetOnlyPrimes()) Out() << "<h2>Results only for primaries presented</h2>" << endl;
    Out() << "<h2>Number of events: " << HM()->H1("hen_EventNo_TrackingQa")->GetEntries() << "</h2>" << endl;
//    Out() << "<h2>Mean impact parameter: " << HM()->H1("Impact parameter")->GetMean() << "</h2>" << endl;
    Out() << "<h2>Mean multiplicity: " << HM()->H1("Multiplicity")->GetMean() << "</h2>" << endl;
    Out() << "<hr>" << endl;
    Out() << "<h3><font color=\"red\">Reconstructable</font> MC-track:</h3>" << "Monte Carlo track with at least <font color=\"red\">4</font> Monte Carlo points in GEM" << endl;
    Out() << "<h3><font color=\"red\">Good</font> track:</h3>" << "Reconstructed track with at least <font color=\"red\">4</font> hits in GEM and <font color=\"red\">60%</font> of them corresponded the same MC-track" << endl;
    Out() << "<h3><font color=\"red\">Clone</font> tracks:</h3>";
    Out() << "Two or more reconstructed tracks with reference to the same MC-track." << endl;
    Out() << "The number of clones is subtracted from number of good tracks before efficiency calculation." << endl;
    return "<hr>";
}

string BmnTrackingQaReport::PrintNofObjects() const {
    vector<TH1*> histos = HM()->H1Vector("hno_NofObjects_.+");
    Int_t nofHistos = histos.size();
    string str = R()->TableBegin("Average number of objects per event", list_of("Name")("Value"));
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        string cellName = Split(histos[iHist]->GetName(), '_')[2];
        str += R()->TableRow(list_of(cellName)(NumberToString<Int_t > (histos[iHist]->GetMean())));
    }
    str += R()->TableEnd();
    return str;
}

string BmnTrackingQaReport::PrintTrackHits() const {
    vector<TH1*> histos = HM()->H1Vector("hth_.+_TrackHits_All");
    Int_t nofHistos = histos.size();
    string str = R()->TableBegin("Average number of all/true/fake hits in tracks",
            list_of("")("all")("true")("fake")("true/all")("fake/all"));
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        string name = histos[iHist]->GetName();
        string cellName = Split(name, '_')[1];
        string all = NumberToString<Float_t > (histos[iHist]->GetMean(), 2);
        string trueh = NumberToString<Float_t > (HM()->H1(FindAndReplace(name, "_All", "_True"))->GetMean(), 2);
        string fakeh = NumberToString<Float_t > (HM()->H1(FindAndReplace(name, "_All", "_Fake"))->GetMean(), 2);
        string toa = NumberToString<Float_t > (HM()->H1(FindAndReplace(name, "_All", "_TrueOverAll"))->GetMean(), 2);
        string foa = NumberToString<Float_t > (HM()->H1(FindAndReplace(name, "_All", "_FakeOverAll"))->GetMean(), 2);
        str += R()->TableRow(list_of(cellName)(all) (trueh) (fakeh) (toa) (foa));
    }
    str += R()->TableEnd();
    return str;
}

string BmnTrackingQaReport::PrintNofGhosts() const {
    Float_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    vector<TH1*> histos = HM()->H1Vector("hng_NofGhosts_.+");
    Int_t nofHistos = histos.size();
    string str = R()->TableBegin("Average number of ghosts per event", list_of("Name")("Value"));
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        string cellName = Split(histos[iHist]->GetName(), '_')[2];
        str += R()->TableRow(list_of(cellName)(NumberToString<Float_t > (histos[iHist]->GetEntries() / nofEvents, 2)));
    }
    str += R()->TableEnd();
    return str;
}

string BmnTrackingQaReport::PrintTrackingEfficiency(Bool_t isPidEfficiency) const {
    string effRegex = "";
    if (isPidEfficiency) effRegex = "hte_(.)*_Eff_p";

    vector<TH1*> histos = HM()->H1Vector(effRegex);
    Int_t nofHistos = histos.size();
    if (nofHistos == 0) return "";

    // Find track and ring categories from the histogram names
    map<string, Int_t> catToCell;
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        string effName = histos[iHist]->GetName();
        pair<string, Int_t> tmp = make_pair(Split(effName, '_')[3], catToCell.size());
        catToCell.insert(tmp);
    }
    Int_t nofCats = catToCell.size();
    Int_t nofRows = nofHistos / nofCats;

    cout << "nofHistos = " << nofHistos << endl;
    cout << "nofCats = " << nofCats << endl;
    cout << "nofRows = " << nofRows << endl;

    vector<string> cat(nofCats);
    map<string, Int_t>::const_iterator it;
    for (it = catToCell.begin(); it != catToCell.end(); it++) {
        cat[(*it).second] = (*it).first;
    }

    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    string str = R()->TableBegin("Tracking efficiency", list_of(string("")).range(cat));
    Int_t histCounter = 0;
    for (Int_t iRow = 0; iRow < nofRows; iRow++) {
        vector<string> cells(nofCats);
        string rowName;
        for (Int_t iCat = 0; iCat < nofCats; iCat++) {
            string effName = histos[histCounter]->GetName();
            string accName = FindAndReplace(effName, "_Eff_", "_Acc_");
            string recName = FindAndReplace(effName, "_Eff_", "_Rec_");
            Float_t acc = HM()->H1(accName)->GetEntries() / nofEvents;
            Float_t rec = HM()->H1(recName)->GetEntries() / nofEvents;
            Float_t eff = (acc != 0.) ? 100. * rec / acc : 0.;
            string accStr = NumberToString<Float_t > (acc);
            string recStr = NumberToString<Float_t > (rec);
            string effStr = NumberToString<Float_t > (eff);
            vector<string> split = Split(effName, '_');
            cells[catToCell[split[3]]] = effStr + "(" + recStr + "/" + accStr + ")";
            histCounter++;
            rowName = split[1] + " (" + split[2] + ")";
        }
        str += R()->TableRow(list_of(rowName).range(cells));
    }
    str += R()->TableEnd();
    return str;
}

void BmnTrackingQaReport::Draw() {
    DrawEventsInfo("Distribution of impact parameter and multiplicity");
//    CalculateEfficiencyHistos();
//    FillGlobalTrackVariants();
    SetDefaultDrawStyle();
    //    DrawEfficiencyHistos();
    DrawEffGem("Distribution of MC-, reco- and fake-tracks vs P_{sim} per event for GEM TRACKS");
    DrawEffEtaGem("Distribution of MC-, reco- and fake-tracks vs Pseudorapidity per event for GEM TRACKS");
    DrawEffThetaGem("Distribution of MC-, reco- and fake-tracks vs theta per event for GEM TRACKS");
    //DrawEffGlob("Distribution of MC-, reco- and fake-tracks vs P_{sim} per event for GLOBAL TRACKS");
    DrawNhitsEtaGem("Distribution of GEM reconstructable MC-tracks (left) and MC-tracks corresponded to reconstructed tracks (right) vs number of MC-points and Pseudorapidity");
    DrawNhitsPGem("Distribution of GEM reconstructable MC-tracks (left) and MC-tracks corresponded to reconstructed tracks (right) vs number of MC-points and Momentum");
    DrawNhitsGem("Distribution of GEM RECO-tracks vs number of hits per track");
    //DrawNhitsGlob("Distribution of GLOBAL RECO-tracks vs number of hits per track");
    //    DrawEffGhostSeed("Distribution of MC-, reco- and fake-tracks vs P_{sim} per event for SEEDS only");
    //    DrawEffGhostGem("Distribution of MC-, reco- and fake-tracks vs P_{sim} per event for GEM TRACKS");
    //    DrawEffGhostGlob("Distribution of MC-, reco- and fake-tracks vs P_{sim} per event for GLOBAL TRACKS");
    // DrawYPtHistos();
    DrawEtaP("Distribution of MC-tracks, GEM-tracks and Global tracks in Pseudorapidity and Momentum");
    DrawPsimPrec("Reco vs MC for GEM-tracks and Global tracks");
    DrawPtSimPtRec("Pt Reco vs MC for GEM-tracks and Global tracks");
    DrawEtaSimEtaRec("Reco vs MC for Pseudorapidity and Momentum");
    DrawTxSimTxRec("Tx Reco vs MC for GEM-tracks and Global tracks");
    DrawTySimTyRec("Ty Reco vs MC for GEM-tracks and Global tracks");
    DrawPsimPrecComponentsGem("Reco vs MC for X-, Y- and Z-component of Momentum for GEM-tracks");
    //DrawPsimPrecComponentsGlob("Reco vs MC for X-, Y- and Z-component of Momentum for Global-tracks");
    DrawMomResGem("Momentum resolution for GEM-tracks", "momRes_2D_gem", "momRes_1D_gem");
//        DrawMomResGem("Momentum resolution for GEM-tracks");
    //DrawMomResGlob("Momentum resolution for Global-tracks");
    //    DrawHitsHistos();
}

void BmnTrackingQaReport::DrawEfficiencyHistos() {
    // Draw global tracking efficiency
    //for (UInt_t i = 0; i < fGlobalTrackVariants.size(); i++) {
    //string variant = fGlobalTrackVariants[i];
    //cout << "variant = " << variant << endl;
    //DrawEfficiency("Global tracking efficiency vs momentum (DETECTOR: " + variant + ")", "hte_Gem.*_" + variant + "_All_Eff_p", DefaultEfficiencyLabelFormatter);
    //DrawEfficiency("Global tracking efficiency vs transverse momentum (DETECTOR: " + variant + ")", "hte_Gem.*_" + variant + "_All_Eff_pt", DefaultEfficiencyLabelFormatter);
    //DrawEfficiency("Global tracking efficiency vs rapidity (DETECTOR: " + variant + ")", "hte_Gem.*_" + variant + "_All_Eff_y", DefaultEfficiencyLabelFormatter);
    //DrawEfficiency("Global tracking efficiency vs polar angle (DETECTOR: " + variant + ")", "hte_Gem.*_" + variant + "_All_Eff_Angle", DefaultEfficiencyLabelFormatter);
    //}

    // Draw local tracking efficiency
    //    vector<string> localTrackVariants = list_of("Gem")("Tof");
    //    vector<string> localTrackVariants;
    //    for (UInt_t i = 0; i < localTrackVariants.size(); i++) {
    //        string variant = localTrackVariants[i];
    //        string re = (variant == "Gem") ? "hte_Gem_Gem_All_Eff" : "hte_" + variant + "_.*_All_Eff";
    //        DrawEfficiency("tracking_qa_local_tracking_efficiency_" + variant + "_p", re + "_p", DefaultEfficiencyLabelFormatter);
    //        DrawEfficiency("tracking_qa_local_tracking_efficiency_" + variant + "_pt", re + "_pt", DefaultEfficiencyLabelFormatter);
    //        DrawEfficiency("tracking_qa_local_tracking_efficiency_" + variant + "_y", "hte_" + variant + "_.*" + variant + ".*_(All|Electron)_Eff_y", DefaultEfficiencyLabelFormatter);
    //
    //        string re2 = (variant == "Gem") ? "hte_Gem_Gem_All_Eff" : "hte_" + variant + "_.*_All_Eff";
    //        DrawEfficiency("tracking_qa_local_tracking_efficiency_" + variant + "_angle", re2 + "_Angle", DefaultEfficiencyLabelFormatter);
    //    }

    // Draw local accepted and reconstructed tracks vs number of points
    HM()->ShrinkEmptyBinsH1ByPattern("hte_.+_.+_.+_.+_Np");
    vector<string> accRecTracks = list_of("Gem")("Tof");
    for (UInt_t i = 0; i < accRecTracks.size(); i++) {
        string variant = accRecTracks[i];

        string re = (variant == "Gem") ? "hte_Gem_Gem_All_(Acc|Rec)_Np" : "hte_" + variant + "_.*_All_(Acc|Rec)_Np";
        DrawAccAndRec("Number of simulated and reconstructed tracks in vs Number of hits", re);

        //re = (variant == "Gem") ? "hte_Gem_Gem_All_(Acc|Rec)_p" : "hte_" + variant + "_.*_All_(Acc|Rec)_p";
        //DrawAccAndRec("Number of simulated and reconstructed tracks in " + variant + " vs momentum", re);
    }

}

void BmnTrackingQaReport::DrawEfficiency(const string& canvasName, const string& histNamePattern, string(*labelFormatter)(const string&, Float_t)) {
    vector<TH1*> histos = HM()->H1Vector(histNamePattern);
    if (histos.size() == 0) return;

    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 600, 500);
    canvas->SetGrid();
    canvas->cd();

    Int_t nofHistos = histos.size();
    vector<string> labels(nofHistos);
    vector<Float_t> efficiencies(nofHistos);
    for (UInt_t iHist = 0; iHist < nofHistos; iHist++) {
        string name = histos[iHist]->GetName();
        efficiencies[iHist] = CalcEfficiency(HM()->H1(FindAndReplace(name, "_Eff_", "_Rec_")), HM()->H1(FindAndReplace(name, "_Eff_", "_Acc_")), 100.);
        labels[iHist] = labelFormatter(name, efficiencies[iHist]);
    }

    DrawH1(histos, labels, kLinear, kLinear, true, 0.6, 0.9, 0.9, 1.0, "PE1");
    DrawMeanEfficiencyLines(histos, efficiencies);
}

void BmnTrackingQaReport::DrawEffGhostSeed(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("allSeedDistr")->Sumw2();
    HM()->H1("allSeedDistr")->Scale(1. / nofEvents);
    HM()->H1("recoSeedDistr")->Sumw2();
    HM()->H1("recoSeedDistr")->Scale(1. / nofEvents);
    HM()->H1("wellSeedDistr")->Sumw2();
    HM()->H1("wellSeedDistr")->Scale(1. / nofEvents);
    HM()->H1("ghostSeedDistr")->Sumw2();
    HM()->H1("ghostSeedDistr")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("allSeedDistr"));
    histos1.push_back(HM()->H1("recoSeedDistr"));
    histos1.push_back(HM()->H1("ghostSeedDistr"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Reco tracks");
    labels1.push_back("Ghost tracks");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    //    HM()->H1("EffSeedDistr")->Divide(HM()->H1("recoSeedDistr"), HM()->H1("allSeedDistr"), 1., 1., "B");
    HM()->H1("EffSeedDistr")->Divide(HM()->H1("wellSeedDistr"), HM()->H1("allSeedDistr"), 1., 1., "B");
    HM()->H1("EffSeedDistr")->Scale(100.0);
    HM()->H1("FakeSeedDistr")->Divide(HM()->H1("ghostSeedDistr"), HM()->H1("recoSeedDistr"), 1., 1., "B");
    HM()->H1("FakeSeedDistr")->Scale(100.0);
    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("EffSeedDistr"));
    histos2.push_back(HM()->H1("FakeSeedDistr"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawEffGhostGem(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("allGemDistr")->Sumw2();
    HM()->H1("allGemDistr")->Scale(1. / nofEvents);
    HM()->H1("recoGemDistr")->Sumw2();
    HM()->H1("recoGemDistr")->Scale(1. / nofEvents);
    HM()->H1("wellGemDistr")->Sumw2();
    HM()->H1("wellGemDistr")->Scale(1. / nofEvents);
    HM()->H1("ghostGemDistr")->Sumw2();
    HM()->H1("ghostGemDistr")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("allGemDistr"));
    histos1.push_back(HM()->H1("recoGemDistr"));
    histos1.push_back(HM()->H1("ghostGemDistr"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Reco tracks");
    labels1.push_back("Ghost tracks");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    //    HM()->H1("EffGemDistr")->Divide(HM()->H1("recoGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("EffGemDistr")->Divide(HM()->H1("wellGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("EffGemDistr")->Scale(100.0);
    HM()->H1("FakeGemDistr")->Divide(HM()->H1("ghostGemDistr"), HM()->H1("recoGemDistr"), 1., 1., "B");
    HM()->H1("FakeGemDistr")->Scale(100.0);
    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("EffGemDistr"));
    histos2.push_back(HM()->H1("FakeGemDistr"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawEffGlob(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("Sim_vs_P_glob")->Sumw2();
    HM()->H1("Sim_vs_P_glob")->Scale(1. / nofEvents);
    HM()->H1("Rec_vs_P_glob")->Sumw2();
    HM()->H1("Rec_vs_P_glob")->Scale(1. / nofEvents);
    HM()->H1("Well_vs_P_glob")->Sumw2();
    HM()->H1("Well_vs_P_glob")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_P_glob")->Sumw2();
    HM()->H1("Ghost_vs_P_glob")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Sim_vs_P_glob"));
    histos1.push_back(HM()->H1("Well_vs_P_glob"));
    histos1.push_back(HM()->H1("Ghost_vs_P_glob"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Reco tracks");
    labels1.push_back("Ghost tracks");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.8, 1.0, 0.99, "PE1", kFALSE);

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");

    //    HM()->H1("EffGemDistr")->Divide(HM()->H1("recoGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("Eff_vs_P_glob")->Divide(HM()->H1("Well_vs_P_glob"), HM()->H1("Sim_vs_P_glob"), 1., 1., "B");
    HM()->H1("Eff_vs_P_glob")->Scale(100.0);
    HM()->H1("Fake_vs_P_glob")->Divide(HM()->H1("Ghost_vs_P_glob"), HM()->H1("Rec_vs_P_glob"), 1., 1., "B");
    HM()->H1("Fake_vs_P_glob")->Scale(100.0);

    for (Int_t i = 0; i < HM()->H1("Eff_vs_P_glob")->GetNbinsX(); ++i) {
        if (HM()->H1("Eff_vs_P_glob")->GetBinContent(i) > 100.0)
            HM()->H1("Eff_vs_P_glob")->SetBinContent(i, 100.0);
        if (HM()->H1("Fake_vs_P_glob")->GetBinContent(i) > 100.0)
            HM()->H1("Fake_vs_P_glob")->SetBinContent(i, 100.0);
    }

    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("Eff_vs_P_glob"));
    histos2.push_back(HM()->H1("Fake_vs_P_glob"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawEffGem(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("Sim_vs_P_gem")->Sumw2();
    HM()->H1("Sim_vs_P_gem")->Scale(1. / nofEvents);
    HM()->H1("Rec_vs_P_gem")->Sumw2();
    HM()->H1("Rec_vs_P_gem")->Scale(1. / nofEvents);
    HM()->H1("Well_vs_P_gem")->Sumw2();
    HM()->H1("Well_vs_P_gem")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_P_gem")->Sumw2();
    HM()->H1("Ghost_vs_P_gem")->Scale(1. / nofEvents);
    HM()->H1("Split_vs_P_gem")->Sumw2();
    HM()->H1("Split_vs_P_gem")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Sim_vs_P_gem"));
    histos1.push_back(HM()->H1("Well_vs_P_gem"));
    histos1.push_back(HM()->H1("Ghost_vs_P_gem"));
    histos1.push_back(HM()->H1("Split_vs_P_gem"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Good tracks");
    labels1.push_back("Ghost tracks");
    labels1.push_back("Clones");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.8, 1.0, 0.99, "PE1", kFALSE);

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    labels2.push_back("Clones");

    //    HM()->H1("EffGemDistr")->Divide(HM()->H1("recoGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("Eff_vs_P_gem")->Divide(HM()->H1("Well_vs_P_gem"), HM()->H1("Sim_vs_P_gem"), 1., 1., "B");
    HM()->H1("Eff_vs_P_gem")->Scale(100.0);
    HM()->H1("Fake_vs_P_gem")->Divide(HM()->H1("Ghost_vs_P_gem"), HM()->H1("Rec_vs_P_gem"), 1., 1., "B");
    HM()->H1("Fake_vs_P_gem")->Scale(100.0);
    HM()->H1("SplitEff_vs_P_gem")->Divide(HM()->H1("Split_vs_P_gem"), HM()->H1("Rec_vs_P_gem"), 1., 1., "B");
    HM()->H1("SplitEff_vs_P_gem")->Scale(100.0);

    // Boundary checking.
    // These cases shouldn't happen, but they happen sometimes...
    for (Int_t i = 0; i < HM()->H1("Eff_vs_P_gem")->GetNbinsX(); ++i) {
        if (HM()->H1("Eff_vs_P_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Eff_vs_P_gem")->SetBinContent(i, 100.0);
            HM()->H1("Eff_vs_P_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_P_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Fake_vs_P_gem")->SetBinContent(i, 100.0);
            HM()->H1("Fake_vs_P_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_P_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("SplitEff_vs_P_gem")->SetBinContent(i, 100.0);
            HM()->H1("SplitEff_vs_P_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Eff_vs_P_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Eff_vs_P_gem")->SetBinContent(i, 0.0);
            HM()->H1("Eff_vs_P_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_P_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Fake_vs_P_gem")->SetBinContent(i, 0.0);
            HM()->H1("Fake_vs_P_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_P_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("SplitEff_vs_P_gem")->SetBinContent(i, 0.0);
            HM()->H1("SplitEff_vs_P_gem")->SetBinError(i, 0.0);
        }
    }

    HM()->H1("Eff_vs_P_gem")->SetMaximum(105.0);
    HM()->H1("Fake_vs_P_gem")->SetMaximum(105.0);
    HM()->H1("SplitEff_vs_P_gem")->SetMaximum(105.0);

    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("Eff_vs_P_gem"));
    histos2.push_back(HM()->H1("Fake_vs_P_gem"));
    histos2.push_back(HM()->H1("SplitEff_vs_P_gem"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1", kFALSE);
}

void BmnTrackingQaReport::DrawEffEtaGem(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("Sim_vs_Eta_gem")->Sumw2();
    HM()->H1("Sim_vs_Eta_gem")->Scale(1. / nofEvents);
    HM()->H1("Rec_vs_Eta_gem")->Sumw2();
    HM()->H1("Rec_vs_Eta_gem")->Scale(1. / nofEvents);
    HM()->H1("Well_vs_Eta_gem")->Sumw2();
    HM()->H1("Well_vs_Eta_gem")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_Eta_gem")->Sumw2();
    HM()->H1("Ghost_vs_Eta_gem")->Scale(1. / nofEvents);
    HM()->H1("Split_vs_Eta_gem")->Sumw2();
    HM()->H1("Split_vs_Eta_gem")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Sim_vs_Eta_gem"));
    histos1.push_back(HM()->H1("Well_vs_Eta_gem"));
    histos1.push_back(HM()->H1("Ghost_vs_Eta_gem"));
    histos1.push_back(HM()->H1("Split_vs_Eta_gem"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Reco tracks");
    labels1.push_back("Ghosts");
    labels1.push_back("Clones");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.8, 1.0, 0.99, "PE1", kFALSE);

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    labels2.push_back("Clones");

    //    HM()->H1("EffGemDistr")->Divide(HM()->H1("recoGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("Eff_vs_Eta_gem")->Divide(HM()->H1("Well_vs_Eta_gem"), HM()->H1("Sim_vs_Eta_gem"), 1., 1., "B");
    HM()->H1("Eff_vs_Eta_gem")->Scale(100.0);
    HM()->H1("Fake_vs_Eta_gem")->Divide(HM()->H1("Ghost_vs_Eta_gem"), HM()->H1("Rec_vs_Eta_gem"), 1., 1., "B");
    HM()->H1("Fake_vs_Eta_gem")->Scale(100.0);
    HM()->H1("SplitEff_vs_Eta_gem")->Divide(HM()->H1("Split_vs_Eta_gem"), HM()->H1("Rec_vs_Eta_gem"), 1., 1., "B");
    HM()->H1("SplitEff_vs_Eta_gem")->Scale(100.0);
    
    // Boundary checking.
    // These cases shouldn't happen, but they happen sometimes...
    for (Int_t i = 0; i < HM()->H1("Eff_vs_Eta_gem")->GetNbinsX(); ++i) {
        if (HM()->H1("Eff_vs_Eta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Eff_vs_Eta_gem")->SetBinContent(i, 100.0);
            HM()->H1("Eff_vs_Eta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_Eta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Fake_vs_Eta_gem")->SetBinContent(i, 100.0);
            HM()->H1("Fake_vs_Eta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_Eta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("SplitEff_vs_Eta_gem")->SetBinContent(i, 100.0);
            HM()->H1("SplitEff_vs_Eta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Eff_vs_Eta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Eff_vs_Eta_gem")->SetBinContent(i, 0.0);
            HM()->H1("Eff_vs_Eta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_Eta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Fake_vs_Eta_gem")->SetBinContent(i, 0.0);
            HM()->H1("Fake_vs_Eta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_Eta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("SplitEff_vs_Eta_gem")->SetBinContent(i, 0.0);
            HM()->H1("SplitEff_vs_Eta_gem")->SetBinError(i, 0.0);
        }
    }

    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("Eff_vs_Eta_gem"));
    histos2.push_back(HM()->H1("Fake_vs_Eta_gem"));
    histos2.push_back(HM()->H1("SplitEff_vs_Eta_gem"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawEffThetaGem(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("Sim_vs_Theta_gem")->Sumw2();
    HM()->H1("Sim_vs_Theta_gem")->Scale(1. / nofEvents);
    HM()->H1("Rec_vs_Theta_gem")->Sumw2();
    HM()->H1("Rec_vs_Theta_gem")->Scale(1. / nofEvents);
    HM()->H1("Well_vs_Theta_gem")->Sumw2();
    HM()->H1("Well_vs_Theta_gem")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_Theta_gem")->Sumw2();
    HM()->H1("Ghost_vs_Theta_gem")->Scale(1. / nofEvents);
    HM()->H1("Split_vs_Theta_gem")->Sumw2();
    HM()->H1("Split_vs_Theta_gem")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Sim_vs_Theta_gem"));
    histos1.push_back(HM()->H1("Well_vs_Theta_gem"));
    histos1.push_back(HM()->H1("Ghost_vs_Theta_gem"));
    histos1.push_back(HM()->H1("Split_vs_Theta_gem"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Reco tracks");
    labels1.push_back("Ghosts");
    labels1.push_back("Clones");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.8, 1.0, 0.99, "PE1", kFALSE);

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    labels2.push_back("Clones");

    //    HM()->H1("EffGemDistr")->Divide(HM()->H1("recoGemDistr"), HM()->H1("allGemDistr"), 1., 1., "B");
    HM()->H1("Eff_vs_Theta_gem")->Divide(HM()->H1("Well_vs_Theta_gem"), HM()->H1("Sim_vs_Theta_gem"), 1., 1., "B");
    HM()->H1("Eff_vs_Theta_gem")->Scale(100.0);
    HM()->H1("Fake_vs_Theta_gem")->Divide(HM()->H1("Ghost_vs_Theta_gem"), HM()->H1("Rec_vs_Theta_gem"), 1., 1., "B");
    HM()->H1("Fake_vs_Theta_gem")->Scale(100.0);
    HM()->H1("SplitEff_vs_Theta_gem")->Divide(HM()->H1("Split_vs_Theta_gem"), HM()->H1("Rec_vs_Theta_gem"), 1., 1., "B");
    HM()->H1("SplitEff_vs_Theta_gem")->Scale(100.0);

    // Boundary checking.
    // These cases shouldn't happen, but they happen sometimes...
    for (Int_t i = 0; i < HM()->H1("Eff_vs_Theta_gem")->GetNbinsX(); ++i) {
        if (HM()->H1("Eff_vs_Theta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Eff_vs_Theta_gem")->SetBinContent(i, 100.0);
            HM()->H1("Eff_vs_Theta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_Theta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("Fake_vs_Theta_gem")->SetBinContent(i, 100.0);
            HM()->H1("Fake_vs_Theta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_Theta_gem")->GetBinContent(i) > 100.0) {
            HM()->H1("SplitEff_vs_Theta_gem")->SetBinContent(i, 100.0);
            HM()->H1("SplitEff_vs_Theta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Eff_vs_Theta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Eff_vs_Theta_gem")->SetBinContent(i, 0.0);
            HM()->H1("Eff_vs_Theta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("Fake_vs_Theta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("Fake_vs_Theta_gem")->SetBinContent(i, 0.0);
            HM()->H1("Fake_vs_Theta_gem")->SetBinError(i, 0.0);
        }
        if (HM()->H1("SplitEff_vs_Theta_gem")->GetBinContent(i) < 0.0) {
            HM()->H1("SplitEff_vs_Theta_gem")->SetBinContent(i, 0.0);
            HM()->H1("SplitEff_vs_Theta_gem")->SetBinError(i, 0.0);
        }
    }

    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("Eff_vs_Theta_gem"));
    histos2.push_back(HM()->H1("Fake_vs_Theta_gem"));
    histos2.push_back(HM()->H1("SplitEff_vs_Theta_gem"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawNhitsEtaGem(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Nh_sim_Eta_sim_gem"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    DrawH2(HM()->H2("Nh_rec_Eta_rec_gem"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawNhitsPGem(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Nh_sim_P_sim_gem"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    DrawH2(HM()->H2("Nh_rec_P_rec_gem"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawNhitsGem(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 600, 600);
    canvas->SetGrid();
    HM()->H1("Well_vs_Nh_gem")->Sumw2();
    HM()->H1("Well_vs_Nh_gem")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_Nh_gem")->Sumw2();
    HM()->H1("Ghost_vs_Nh_gem")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Well_vs_Nh_gem"));
    histos1.push_back(HM()->H1("Ghost_vs_Nh_gem"));
    vector<string> labels1;
    labels1.push_back("Good tracks");
    labels1.push_back("Ghosts");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1", kFALSE);
}

void BmnTrackingQaReport::DrawNhitsGlob(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 600, 600);
    canvas->SetGrid();
    HM()->H1("Well_vs_Nh_glob")->Sumw2();
    HM()->H1("Well_vs_Nh_glob")->Scale(1. / nofEvents);
    HM()->H1("Ghost_vs_Nh_glob")->Sumw2();
    HM()->H1("Ghost_vs_Nh_glob")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("Well_vs_Nh_glob"));
    histos1.push_back(HM()->H1("Ghost_vs_Nh_glob"));
    vector<string> labels1;
    labels1.push_back("Good tracks");
    labels1.push_back("Ghosts");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.5, 0.9, 1.0, 0.99, "PE1", kFALSE);
}

void BmnTrackingQaReport::DrawEffGhostGlob(const string& canvasName) {
    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    HM()->H1("allGlobDistr")->Sumw2();
    HM()->H1("allGlobDistr")->Scale(1. / nofEvents);
    HM()->H1("recoGlobDistr")->Sumw2();
    HM()->H1("recoGlobDistr")->Scale(1. / nofEvents);
    HM()->H1("ghostGlobDistr")->Sumw2();
    HM()->H1("ghostGlobDistr")->Scale(1. / nofEvents);
    vector<TH1*> histos1;
    histos1.push_back(HM()->H1("allGlobDistr"));
    histos1.push_back(HM()->H1("recoGlobDistr"));
    histos1.push_back(HM()->H1("ghostGlobDistr"));
    vector<string> labels1;
    labels1.push_back("MC tracks");
    labels1.push_back("Good tracks");
    labels1.push_back("Ghosts");
    DrawH1(histos1, labels1, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");

    canvas->cd(2);
    vector<string> labels2;
    labels2.push_back("Efficiency");
    labels2.push_back("Ghosts");
    HM()->H1("EffGlobDistr")->Divide(HM()->H1("recoGlobDistr"), HM()->H1("allGlobDistr"), 1., 1., "B");
    HM()->H1("EffGlobDistr")->Scale(100.0);
    HM()->H1("FakeGlobDistr")->Divide(HM()->H1("ghostGlobDistr"), HM()->H1("recoGlobDistr"), 1., 1., "B");
    HM()->H1("FakeGlobDistr")->Scale(100.0);
    vector<TH1*> histos2;
    histos2.push_back(HM()->H1("EffGlobDistr"));
    histos2.push_back(HM()->H1("FakeGlobDistr"));
    DrawH1(histos2, labels2, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99, "PE1");
}

void BmnTrackingQaReport::DrawPsimPrec(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("P_rec_P_sim_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(2);
    //DrawH2(HM()->H2("P_rec_P_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawPtSimPtRec(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Pt_rec_Pt_sim_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(2);
    //DrawH2(HM()->H2("Pt_rec_Pt_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawTxSimTxRec(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Tx_rec_Tx_sim_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(2);
    //DrawH2(HM()->H2("Tx_rec_Tx_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawTySimTyRec(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Ty_rec_Ty_sim_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(2);
    //DrawH2(HM()->H2("Ty_rec_Ty_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawEtaSimEtaRec(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Eta_rec_Eta_sim_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(2);
    //DrawH2(HM()->H2("Eta_rec_Eta_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawPsimPrecComponentsGem(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 500);
    canvas->SetGrid();
    canvas->Divide(3, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Px_rec_Px_sim_gem"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    DrawH2(HM()->H2("Py_rec_Py_sim_gem"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(3);
    DrawH2(HM()->H2("Pz_rec_Pz_sim_gem"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawPsimPrecComponentsGlob(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 500);
    canvas->SetGrid();
    canvas->Divide(3, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("Px_rec_Px_sim_glob"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    DrawH2(HM()->H2("Py_rec_Py_sim_glob"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(3);
    DrawH2(HM()->H2("Pz_rec_Pz_sim_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawEventsInfo(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 500);
    canvas->SetGrid();
    canvas->Divide(3, 1);
    canvas->cd(1);
    DrawH1(HM()->H1("Impact parameter"), kLinear, kLinear, "", kRed, 2, 1, 1.1, 20, 33);
    canvas->cd(2);
    DrawH1(HM()->H1("Multiplicity"), kLinear, kLinear, "", kRed, 2, 1, 1.1, 20, 33);
    canvas->cd(3);
    DrawH2(HM()->H2("Impact_Mult"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawMomResGlob(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2("momRes_2D_glob"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    //    HM()->H1("momRes_1D_glob")->SetMaximum(50.0);
    HM()->H1("momRes_1D_glob")->SetMinimum(0.0);
    DrawH1(HM()->H1("momRes_1D_glob"), kLinear, kLinear, "PE1", kRed, 0.7, 0.75, 1.1, 20);
}

void BmnTrackingQaReport::DrawMomResGem(const string& canvasName, TString name2d, TString name1d) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1000, 500);
    canvas->SetGrid();
    canvas->Divide(2, 1);
    canvas->cd(1);
    DrawH2(HM()->H2(name2d.Data()), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    //    HM()->H1("momRes_1D_gem")->SetMaximum(50.0);
    HM()->H1(name1d.Data())->SetMinimum(0.0);
    DrawH1(HM()->H1(name1d.Data()), kLinear, kLinear, "PE1", kRed, 0.7, 0.75, 1.1, 20);
}

void BmnTrackingQaReport::DrawEtaP(const string& canvasName) {
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1500, 500);
    canvas->SetGrid();
    canvas->Divide(3, 1);
    canvas->cd(1);
    canvas->GetPad(1)->SetTitle("MC-tracks");
    DrawH2(HM()->H2("EtaP_sim"), kLinear, kLinear, kLinear, "colz");
    canvas->cd(2);
    canvas->GetPad(2)->SetTitle("GEM tracks only");
    DrawH2(HM()->H2("EtaP_rec_gem"), kLinear, kLinear, kLinear, "colz");
    //canvas->cd(3);
    //canvas->GetPad(3)->SetTitle("Global tracks");
    //DrawH2(HM()->H2("EtaP_rec_glob"), kLinear, kLinear, kLinear, "colz");
}

void BmnTrackingQaReport::DrawMeanEfficiencyLines(
        const vector<TH1*>& histos,
        const vector<Float_t>& efficiencies) {
    assert(histos.size() != 0 && efficiencies.size() == histos.size());

    Float_t minX = histos[0]->GetXaxis()->GetXmin();
    Float_t maxX = histos[0]->GetXaxis()->GetXmax();
    Int_t nofHistos = histos.size();
    for (UInt_t iHist = 0; iHist < nofHistos; iHist++) {
        TLine* line = new TLine(minX, efficiencies[iHist], maxX, efficiencies[iHist]);
        line->SetLineWidth(1);
        line->SetLineColor(histos[iHist]->GetLineColor());
        line->Draw();
    }
}

//void BmnTrackingQaReport:://DrawMeanLine(TH1* hist) {
//
//    Float_t minX = hist->GetXaxis()->GetXmin();
//    Float_t maxX = hist->GetXaxis()->GetXmax();
//    Int_t nonZeroBins = 0;
//    for (Int_t i = 0; i < hist->GetNbinsX(); ++i) {
//        if (hist->GetBinContent(i) != 0.0) nonZeroBins++;
//    }
//    TLine* line = new TLine(minX, hist->Integral() / nonZeroBins, maxX, hist->Integral() / nonZeroBins);
//    line->SetLineWidth(2);
//    line->SetLineColor(hist->GetLineColor());
//    line->Draw();
//}

void BmnTrackingQaReport::DrawAccAndRec(const string& canvasName, const string& histNamePattern) {
    vector<TH1*> histos = HM()->H1Vector(histNamePattern);
    if (histos.size() == 0) return;

    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 600, 500);
    canvas->SetGrid();
    canvas->cd();

    Int_t nofEvents = HM()->H1("hen_EventNo_TrackingQa")->GetEntries();
    Int_t nofHistos = histos.size();
    vector<string> labels(nofHistos);
    for (UInt_t iHist = 0; iHist < nofHistos; iHist++) {
        TH1* hist = histos[iHist];
        hist->Scale(1. / nofEvents);
        string name = hist->GetName();
        vector<string> split = Split(name, '_');
        labels[iHist] = split[4] + ":" + split[3] + "(" + NumberToString<Float_t > (hist->GetEntries() / nofEvents, 1) + ")";
    }

    DrawH1(histos, labels, kLinear, kLinear, true, 0.7, 0.75, 1.0, 0.99);
}

void BmnTrackingQaReport::DrawYPtHistos() {
    // Draw global tracking efficiency
    for (UInt_t i = 0; i < fGlobalTrackVariants.size(); i++) {
        string variant = fGlobalTrackVariants[i];
        string effHistName = "hte_" + variant + "_" + variant;
        DrawYPt("Pt and Y distribution of all MC-tracks (1), MC-tracks linked to RECO-tracks (2) and efficiency of reconstruction (3) for " + variant, effHistName + "_All_Eff_YPt");
    }
}

void BmnTrackingQaReport::DrawYPt(const string& canvasName, const string& effHistName, Bool_t drawOnlyEfficiency) {
    string accHistName = FindAndReplace(effHistName, "_Eff_", "_Acc_");
    string recHistName = FindAndReplace(effHistName, "_Eff_", "_Rec_");

    if (!(HM()->Exists(effHistName) && HM()->Exists(accHistName) && HM()->Exists(recHistName))) return;

    TH2* accHist = HM()->H2(accHistName);
    TH2* recHist = HM()->H2(recHistName);
    TH2* effHist = HM()->H2(effHistName);

    if (drawOnlyEfficiency) {
        TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 800, 800);
        //canvas->SetGrid();
        DrawH2(effHist);
    } else {
        TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 600);
        canvas->Divide(3, 1);
        canvas->SetGrid();
        canvas->cd(1);
        DrawH2(accHist);

        canvas->cd(2);
        DrawH2(recHist);

        canvas->cd(3);
        DrawH2(effHist);
    }
}

void BmnTrackingQaReport::DrawHitsHistos() {
    HM()->ShrinkEmptyBinsH1ByPattern("hth_.*(_All|_True|_Fake)");
    DrawHitsHistos("Number of all-real-fake tracks vs Number of hit", "hth_Gem_TrackHits");
}

void BmnTrackingQaReport::DrawHitsHistos(const string& canvasName, const string& hist) {
    if (!(HM()->Exists(hist + "_All") && HM()->Exists(hist + "_True") &&
            HM()->Exists(hist + "_Fake") && HM()->Exists(hist + "_TrueOverAll") &&
            HM()->Exists(hist + "_FakeOverAll"))) return;
    TCanvas* canvas = CreateCanvas(canvasName.c_str(), canvasName.c_str(), 1200, 600);
    canvas->Divide(2, 1);
    canvas->SetGrid();

    canvas->cd(1);
    TH1* hAll = HM()->H1(hist + "_All");
    TH1* hTrue = HM()->H1(hist + "_True");
    TH1* hFake = HM()->H1(hist + "_Fake");
    DrawH1(list_of(hAll)(hTrue) (hFake),
            list_of("all: " + NumberToString<Float_t > (hAll->GetMean(), 1))
            ("true: " + NumberToString<Float_t > (hTrue->GetMean(), 1))
            ("fake: " + NumberToString<Float_t > (hFake->GetMean(), 1)),
            kLinear, kLog, true, 0.25, 0.99, 0.55, 0.75);

    canvas->cd(2);
    TH1* hTrueOverAll = HM()->H1(hist + "_TrueOverAll");
    TH1* hFakeOverAll = HM()->H1(hist + "_FakeOverAll");
    DrawH1(list_of(hTrueOverAll)(hFakeOverAll),
            list_of("true/all: " + NumberToString<Float_t > (hTrueOverAll->GetMean()))
            ("fake/all: " + NumberToString<Float_t > (hFakeOverAll->GetMean())),
            kLinear, kLog, true, 0.25, 0.99, 0.55, 0.75);
}

Float_t BmnTrackingQaReport::CalcEfficiency(
        const TH1* histRec,
        const TH1* histAcc,
        Float_t scale) const {
    if (histAcc->Integral() == 0 || histRec->Integral() == 0) {
        return 0.;
    } else {
        return scale * Float_t(histRec->Integral()) / Float_t(histAcc->Integral());
    }
}

void BmnTrackingQaReport::FillGlobalTrackVariants() {
    fGlobalTrackVariants.clear();
    vector<TH1*> histos = HM()->H1Vector("hte_.*_Eff_p");
    Int_t nofHistos = histos.size();
    set<string> variants;
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        string effName = histos[iHist]->GetName();
        variants.insert(Split(effName, '_')[2]);
    }
    fGlobalTrackVariants.assign(variants.begin(), variants.end());
}

void BmnTrackingQaReport::DivideHistos(
        TH1* histo1,
        TH1* histo2,
        TH1* histo3,
        Float_t scale) {
    histo1->Sumw2();
    histo2->Sumw2();
    histo3->Sumw2();
    histo3->Divide(histo1, histo2, 1., 1., "B");
    histo3->Scale(scale);
}

void BmnTrackingQaReport::CalculateEfficiencyHistos() {
    vector<TH1*> effHistos = HM()->H1Vector("(hte|hpe)_.+_Eff_.+");
    Int_t nofEffHistos = effHistos.size();
    for (Int_t iHist = 0; iHist < nofEffHistos; iHist++) {
        TH1* effHist = effHistos[iHist];
        string effHistName = effHist->GetName();
        string accHistName = FindAndReplace(effHistName, "_Eff_", "_Acc_");
        string recHistName = FindAndReplace(effHistName, "_Eff_", "_Rec_");
        DivideHistos(HM()->H1(recHistName), HM()->H1(accHistName), effHist, 100.);
        effHist->SetMinimum(0.);
        effHist->SetMaximum(100.);
    }
}

void BmnTrackingQaReport::CalculatePionSuppressionHistos() {
    vector<TH1*> histos = HM()->H1Vector("hps_.+_PionSup_.+");
    Int_t nofHistos = histos.size();
    for (Int_t iHist = 0; iHist < nofHistos; iHist++) {
        TH1* psHist = histos[iHist];
        string psHistName = psHist->GetName();
        string recHistName = FindAndReplace(psHistName, "_PionSup_", "_Rec_");
        string pionRecHistName = FindAndReplace(psHistName, "_PionSup_", "_RecPions_");
        DivideHistos(HM()->H1(pionRecHistName), HM()->H1(recHistName), psHist, 1.);
        //  psHist->SetMinimum(1.);
        //  psHist->SetMaximum(20000.);
    }
}

ClassImp(BmnTrackingQaReport)
