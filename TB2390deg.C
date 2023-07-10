//**************************************************
// \file TB2390deg.C
// \brief: Analysis of 20223 TB 90 degree runs
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
// 	    @lopezzot
// \start date: 6 July 2023
//**************************************************

// Usage: root -l TB2390deg.C

#include "PhysicsEvent.h"
#include "utils.h"

#include <array>
#include <iostream>
#include <string>
#include <cstdlib>

ClassImp(EventOut);

struct AttenuatedEnergies
{
    double TotSiPMSene;
    double TotSiPMCene;
    double ErTotSiPMSene;
    double ErTotSiPMCene;
};

std::array<double, 2> GetDWCoffset(const int& DWCIdx)
{

  std::array<double, 2> DWCoffset{};
  if (DWCIdx == 1) {  // DWC 1
    DWCoffset[0] = 1.17;
    DWCoffset[1] = -4.3;
  }
  else if (DWCIdx == 2) {  // DWC 2
    DWCoffset[0] = 4.3;
    DWCoffset[1] = 0.0;
  }
  else {
    std::cout << "Wrong DWC index (1 or 2), going to std::abort()." << std::endl;
    std::abort();
  }

  return DWCoffset;
}

AttenuatedEnergies DoAnalysis(TTree* tree, const int& runno);

void TB2390deg()
{
  // run numbers for 90 deg (40 gev) e+ runs
  // using vertical angle 0.0 deg
  const std::array<int, 5> runNumbers{277, 276, 285, 279, 281};
  // distances (mm)
  const std::array<double, 5> distances{-352.9, -162.8, 33.8, 225.3, 416.0};  // mm (as on Twiki)
  std::array<double, 5> realdistances{};
  int distcounter = 0;
  // 33.8 mm corresponds to the center of the calo, calo is 1 m long
  for (const auto& d : distances) {
    realdistances[distcounter] = 500. + 33.8 - distances[distcounter];
    distcounter++;
  }
  // path to files + filename
  const std::string path("recoNtuple/physics_sps2023_run");
  // array with average reconstructed energies (initialized to 0)
  std::array<double, 5> sipmsenergyatdist{};
  std::array<double, 5> sipmcenergyatdist{};
  std::array<double, 5> ersipmsenergyatdist{};
  std::array<double, 5> ersipmcenergyatdist{};
  std::array<double, 5> erdist{};

  // loop through runs
  for (int runno = 0; runno < 5; runno++) {
    std::string pathFile(path + std::to_string(runNumbers[runno]) + ".root");
    std::cout << "Using file: " << pathFile << std::endl;
    TFile file{pathFile.c_str()};
    TTree* tree = static_cast<TTree*>(file.Get("Ftree"));
    AttenuatedEnergies attEne = DoAnalysis(tree, runNumbers[runno]);
    sipmsenergyatdist[runno] = attEne.TotSiPMSene;
    sipmcenergyatdist[runno] = attEne.TotSiPMCene;
    ersipmsenergyatdist[runno] = attEne.ErTotSiPMSene;
    ersipmcenergyatdist[runno] = attEne.ErTotSiPMCene;
  }

  TGraphErrors GrAttenuationsipms{5, realdistances.data(), sipmsenergyatdist.data(), erdist.data(),
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
  FinalFile->Close();
}

AttenuatedEnergies DoAnalysis(TTree* tree, const int& runno)
{
  TFile* analysisFile(
    TFile::Open(("AnalysisRun_" + std::to_string(runno) + ".root").c_str(), "RECREATE"));
  TH1F H1SiPMSene{"H1SiPMSene", "H1SiPMSene", 100, 0., 15.};
  TH1F H1SiPMCene{"H1SiPMCene", "H1SiPMCene", 100, 0., 15.};
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
    if (!(utils::IsPositronPsMu(pevtout->PShower, pevtout->MCounter)))
      continue;  // select
                 // positrons
    double DWC1pos[2] = {pevtout->XDWC1, pevtout->YDWC1};
    double DWC2pos[2] = {pevtout->XDWC2, pevtout->YDWC2};
    if (!(utils::IsDWCradius(DWC1pos, 10.0, GetDWCoffset(1)))) continue;  // cut out-of-radius DWC1
    if (!(utils::IsDWCradius(DWC2pos, 10.0, GetDWCoffset(2)))) continue;  // cut out-of-radius DWC2
    if (!(utils::IsSiPMSabovecut(pevtout->SiPMPheS, 0.1)))
      continue;  // cut SiPMSene < 0.1 GeV events

    evtcounter++;
    // fill histos of this run
    H1SiPMSene.Fill(pevtout->totSiPMSene);
    H1SiPMCene.Fill(pevtout->totSiPMCene);
    H2DWC1.Fill(pevtout->XDWC1, pevtout->YDWC1);
    H2DWC2.Fill(pevtout->XDWC2, pevtout->YDWC2);
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
  H2DWC1.Write();
  H2DWC2.Write();
  H2SiPMSbar.Write();
  H2SiPMCbar.Write();
  analysisFile->Close();

  AttenuatedEnergies attene{SiPMSene, SiPMCene, ErSiPMSene, ErSiPMCene};
  return attene;
}

//**************************************************
