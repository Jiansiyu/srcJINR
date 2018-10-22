// MWPC
// Silicon
//

// construct detector geometry
void geometry(FairRunSim *fRun)
{
    // Set Material file Name
    fRun->SetMaterials("media.geo");

    // -----   Create passive volumes   -------------------------
    FairModule* cave = new FairCave("CAVE");
    cave->SetGeometryFileName("cave.geo");
    fRun->AddModule(cave);

    FairModule* magnet = new FairMagnet("MAGNET");
    magnet->SetGeometryFileName("magnet_modified.root");
    fRun->AddModule(magnet);

    // -----   Create detectors        -------------------------
    FairDetector* tof = new BmnTOF("TOF", kTRUE);
    tof->SetGeometryFileName("tof700_run7_with_support.root");
    fRun->AddModule(tof);

    //CbmPsdv1* psd = new CbmPsdv1("PSD", kTRUE);
    //psd->SetXshift(45.);
    //psd->SetZposition(1000.);
    //psd->SetHole(1); // 0 for no hole
    //fRun->AddModule(psd);

    FairDetector* sts = new CbmSts("STS", kTRUE);
    sts->SetGeometryFileName("GEMS_RunSRCSpring2018.root");
    fRun->AddModule(sts);

    FairDetector* tof1 = new BmnTOF1("TOF1", kTRUE);
    tof1->SetGeometryFileName("TOF400_RUN7_SRC_new.root");
    fRun->AddModule(tof1);

    FairDetector* dch = new BmnDch("DCH", kTRUE);
    dch->SetGeometryFileName("DCH_RunWinter2016.root");
    fRun->AddModule(dch);

    FairDetector* mwpc = new BmnMwpc("MWPC", kTRUE);
    mwpc->SetGeometryFileName("MWPC_SRC2018.root");
    fRun->AddModule(mwpc);

    //FairDetector* emc = new BmnEcal("EMC", kTRUE);
    //emc->SetGeometryFileName("ecal_v1_0.geo");
    //fRun->AddModule(emc);

    FairDetector* silicon = new BmnSilicon("SILICON", kTRUE);
    silicon->SetGeometryFileName("Silicon_RunSRCSpring2018.root");
    fRun->AddModule(silicon);

}
