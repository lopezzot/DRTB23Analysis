//**************************************************
// \file TB23pion.C
// \brief: Analysis of 20223 TB pion runs
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
// 	    @lopezzot
// \start date: 10 July 2023
//**************************************************

// Usage: root -l TB23pion.C

#include "PhysicsEvent.h"
#include "utils.h"

#include <array>
#include <iostream>
#include <string>

ClassImp(EventOut);

struct AttenuatedEnergies
{
    double TotSiPMSene;
    double TotSiPMCene;
    double ErTotSiPMSene;
    double ErTotSiPMCene;
};

void DoAnalysis(TTree* tree, const int& runno, const int& energy);

std::array<double, 2> GetDWCoffset(const int& DWCIdx)
{

  std::array<double, 2> DWCoffset{};
  if (DWCIdx == 1) {  // DWC 1
    DWCoffset[0] = +1.351;
    DWCoffset[1] = -4.451;
  }
  else if (DWCIdx == 2) {  // DWC 2
    DWCoffset[0] = +4.613;
    DWCoffset[1] = +0.1348;
  }
  else {
    std::cout << "Wrong DWC index (1 or 2), going to std::abort()." << std::endl;
    std::abort();
  }

  return DWCoffset;
}

void TB23pion()
{
  // run numbers for pi+ runs
  // using vertical angle 0.0 deg
  // using orizonthal angle 0.0 deg
  //
  const std::array<int, 4> runNumbers{247, 246, 245, 244};
  // corresponding energies (GeV)
  const std::array<int, 4> energies{20, 60, 100, 180};
  // path to files + filename
  const std::string path("recoNtuple/physics_sps2023_run");

  // loop through runs
  for (int runno = 0; runno < 4; runno++) {
    std::string pathFile(path + std::to_string(runNumbers[runno]) + ".root");
    std::cout << "Using file: " << pathFile << std::endl;
    TFile file{pathFile.c_str()};
    TTree* tree = static_cast<TTree*>(file.Get("Ftree"));
    DoAnalysis(tree, runNumbers[runno], energies[runno]);
  }

  /*TGraphErrors GrAttenuationsipms{5, realdistances.data(), sipmsenergyatdist.data(), erdist.data(),
                                  ersipmsenergyatdist.data()};
  //.data() converts from std::array to array
  GrAttenuationsipms.SetTitle(
    "Attenuation SiPMSenergy;Distance to readout (mm);SiPM_S signal (GeV)");
  GrAttenuationsipms.SetName("AttSSiPM");
  TGraphErrors GrAttenuationsipmc{5, realdistances.data(), sipmcenergyatdist.data(), erdist.data(),
                                  ersipmcenergyatdist.data()};
  //.data() converts from std::array to array
  GrAttenuationsipmc.SetTitle(
    "Attenuation SiPMCenergy;Distance to readout (mm);SiPM_C signal (GeV)");
  GrAttenuationsipmc.SetName("AttCSiPM");

  TFile* FinalFile(TFile::Open("Attenuation.root", "RECREATE"));
  GrAttenuationsipms.Write();
  GrAttenuationsipmc.Write();
  FinalFile->Close();*/
}

void DoAnalysis(TTree* tree, const int& runno, const int& energy)
{
  TFile* analysisFile(
    TFile::Open(("AnalysisRun_" + std::to_string(runno) + ".root").c_str(), "RECREATE"));
  // SiPM histos
  TH1F H1SiPMSene{"H1SiPMSene", "H1SiPMSene", 100, 0., static_cast<double>(energy)};
  TH1F H1SiPMCene{"H1SiPMCene", "H1SiPMCene", 100, 0., static_cast<double>(energy)};
  TH2F H2SiPMSCene{"H2SiPMSCene", "H2SiPMSCene", 100, 0., static_cast<double>(energy), 100, 0., static_cast<double>(energy)};
  // PMT histos
  TH1F H1PMTSene{"H1PMTSene", "H1PMTSene", 100, 0., static_cast<double>(energy)};
  // Energy ratio histos
  TH1F H1RatioS{"H1RatioS", "H1RatioS", 100, 0., 1.3};
  // Geometry histos
  TH2F H2DWC1{"DWC1", "DWC1", 400, -20., 20., 400, -20, 20};
  TH2F H2DWC2{"DWC2", "DWC2", 400, -20., 20., 400, -20, 20};
  TH2F H2SiPMSbar{"SiPMSbar", "SiPMSbar", 400, -20., 40., 400, -20, 40};
  TH2F H2SiPMCbar{"SiPMCbar", "SiPMcbar", 400, -20., 40., 400, -20, 40};

  EventOut evtout{};
  EventOut* pevtout = &evtout;
  tree->SetBranchAddress("Events", &pevtout);

  int evtcounter = 0;
  // loop over events
  for (int evtno = 0; evtno < tree->GetEntries(); evtno++) {
    tree->GetEntry(evtno);
    // filter events
    if (!(utils::IsPionPsMu(pevtout->PShower, pevtout->MCounter)))
      continue;  // select pions
    double DWC1pos[2] = {pevtout->XDWC1, pevtout->YDWC1};
    double DWC2pos[2] = {pevtout->XDWC2, pevtout->YDWC2};
    if (!(utils::IsDWCradius(DWC1pos, 10.0, GetDWCoffset(1)))) continue;  // cut out-of-radius DWC1
    if (!(utils::IsDWCradius(DWC2pos, 10.0, GetDWCoffset(2)))) continue;  // cut out-of-radius DWC2
    // cut events with non-interacting pions
    if (!(utils::HasPionInteracted(pevtout->totSiPMSene, static_cast<double>(energy)))) continue;

    evtcounter++;
    // fill histos of this run
    // energy histos
    H1SiPMSene.Fill(pevtout->totSiPMSene);
    H1SiPMCene.Fill(pevtout->totSiPMCene);
    H2SiPMSCene.Fill(pevtout->totSiPMSene, pevtout->totSiPMCene);
    H1PMTSene.Fill(pevtout->SPMTenergy);
    // Energy ratio histos
    H1RatioS.Fill(pevtout->totSiPMSene/(pevtout->totSiPMSene+pevtout->SPMTenergy));
    // geometry histos
    H2DWC1.Fill(pevtout->XDWC1+GetDWCoffset(1)[0], pevtout->YDWC1+GetDWCoffset(1)[1]);
    H2DWC2.Fill(pevtout->XDWC2+GetDWCoffset(2)[0], pevtout->YDWC2+GetDWCoffset(2)[1]);
    double* sbar = utils::GetScinbar(pevtout->SiPMPheS);
    double* cbar = utils::GetCherbar(pevtout->SiPMPheC);
    H2SiPMSbar.Fill(sbar[0], sbar[1]);
    H2SiPMCbar.Fill(cbar[0], cbar[1]);
  }

  double SiPMSene = H1SiPMSene.GetMean();
  double SiPMCene = H1SiPMCene.GetMean();
  double ErSiPMSene = H1SiPMSene.GetMeanError();
  double ErSiPMCene = H1SiPMCene.GetMeanError();
  std::cout << "--->" << tree->GetEntries() << " evts, " << evtcounter
            << " selected, average SiPMene S: " << SiPMSene << " C: " << SiPMCene << " GeV"
            << std::endl;

  H1SiPMSene.Write();
  H1SiPMCene.Write();
  H2SiPMSCene.Write();
  H1PMTSene.Write();
  H1RatioS.Write();
  H2DWC1.Write();
  H2DWC2.Write();
  H2SiPMSbar.Write();
  H2SiPMCbar.Write();
  analysisFile->Close();

  //AttenuatedEnergies attene{SiPMSene, SiPMCene, ErSiPMSene, ErSiPMCene};
  //return attene;
}

//**************************************************
