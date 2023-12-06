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

#include <algorithm>
#include <array>
#include <iostream>
#include <string>

ClassImp(EventOut);

struct outputAnalysis
{
    double RatioSenergy;
    double ErRatioSenergy;
    double RatioCenergy;
    double ErRatioCenergy;
    double Senergy;
    double ErSenergy;
    double Cenergy;
    double ErCenergy;
};

outputAnalysis DoAnalysis(TTree* tree, const int& runno, const int& energy,
                          const bool isMC = false);
void CompareHistos(TFile* CompareFile);

void TB23pion()
{
  // run numbers for pi+ runs
  // using vertical angle 0.0 deg
  // using orizonthal angle 0.0 deg
  //
  const std::array<int, 4> runNumbers{247, 246, 245, 244};
  // corresponding energies (GeV)
  const std::array<int, 4> energies{20, 60, 100, 180};
  // corresponding energies as doubles
  std::array<double, 4> Denergies{};
  std::transform(energies.begin(), energies.end(), Denergies.begin(),
                 [](int value) { return static_cast<double>(value); });
  // energies errors (0)
  const std::array<double, 4> Erenergies{};
  // path to files + filename
  const std::string path("recoNtuple/physics_sps2023_run");
  // ratios of S energy central tower to all towers
  std::array<double, 4> ratioS{};
  // error on ratios of S energy central tower to all towers
  std::array<double, 4> erratioS{};
  // ratios of C energy central tower to all towers
  std::array<double, 4> ratioC{};
  // error on ratios of C energy central tower to all towers
  std::array<double, 4> erratioC{};
  // S energy all towers
  std::array<double, 4> Senergy{};
  // error S energy all towers
  std::array<double, 4> ErSenergy{};
  // C energy all towers
  std::array<double, 4> Cenergy{};
  // error C energy all towers
  std::array<double, 4> ErCenergy{};

  // loop through runs (test-beam data only)
  for (int runno = 0; runno < 4; runno++) {
    std::string pathFile(path + std::to_string(runNumbers[runno]) + ".root");
    std::cout << "Using file: " << pathFile << std::endl;
    TFile file{pathFile.c_str()};
    TTree* tree = static_cast<TTree*>(file.Get("Ftree"));
    outputAnalysis output = DoAnalysis(tree, runNumbers[runno], energies[runno]);
    ratioS[runno] = output.RatioSenergy;
    erratioS[runno] = output.ErRatioSenergy;
    ratioC[runno] = output.RatioCenergy;
    erratioC[runno] = output.ErRatioCenergy;
    Senergy[runno] = output.Senergy;
    ErSenergy[runno] = output.ErSenergy;
    Cenergy[runno] = output.Cenergy;
    ErCenergy[runno] = output.ErCenergy;
  }

  TGraphErrors GrRatioS{4, Denergies.data(), ratioS.data(), Erenergies.data(), erratioS.data()};
  //.data() converts from std::array to array
  GrRatioS.SetTitle("Central Tower S signal ratio;Beam energy (GeV);S ratio");
  GrRatioS.SetName("RatioS");

  TGraphErrors GrRatioC{4, Denergies.data(), ratioC.data(), Erenergies.data(), erratioC.data()};
  //.data() converts from std::array to array
  GrRatioC.SetTitle("Central Tower C signal ratio;Beam energy (GeV);C ratio");
  GrRatioC.SetName("RatioC");

  TGraphErrors GrSenergy{4, Denergies.data(), Senergy.data(), Erenergies.data(), ErSenergy.data()};
  GrSenergy.SetTitle("Scintillation;Beam energy (GeV);Signal S (GeV)");
  GrSenergy.SetName("Senergy");

  TGraphErrors GrCenergy{4, Denergies.data(), Cenergy.data(), Erenergies.data(), ErCenergy.data()};
  GrCenergy.SetTitle("Cherenkov;Beam energy (GeV);Signal C (GeV)");
  GrCenergy.SetName("Cenergy");

  TFile* FinalFile(TFile::Open("Pions.root", "RECREATE"));
  GrRatioS.Write();
  GrRatioC.Write();
  GrSenergy.Write();
  GrCenergy.Write();

  // run numbers for pi+ Geant4 runs
  // using vertical angle 0.0 deg
  // using orizonthal angle -0.5 deg
  //
  // const std::array<int, 2> G4runNumbers{84, 85}; // horiz=-0.5deg, x0y0 mm
  // const std::array<int, 2> G4runNumbers{87, 86}; // with noise, horiz=-0.5deg, x0y0 mm
  // const std::array<int, 2> G4runNumbers{90, 89}; //with noise, horiz=-0.5 deg, x=-5.0y0 mm
  // const std::array<int, 2> G4runNumbers{91, 92}; //with noise, horiz=-0.7 deg, x=0y0 mm
  // const std::array<int, 2> G4runNumbers{94, 93}; //with noise, horiz=-0.3 deg, x=0y0 mm
  // const std::array<int, 2> G4runNumbers{96, 95}; //with noise, horiz=-0.5 deg, x=3.0y0 mm
  // new converter
  // const std::array<int, 2> G4runNumbers{97, 98}; //with noise, horiz=-0.5 deg, x=-2.0y0 mm
  // const std::array<int, 2> G4runNumbers{99, 100}; //with noise, horiz=-0.5 deg, x=-5.0y0 mm
  // const std::array<int, 2> G4runNumbers{101, 102}; //with noise, horiz=-0.5 deg, x=-5.0y2.0 mm
  // const std::array<int, 2> G4runNumbers{103, 104}; //with noise, horiz=-0.5 deg, x=-5.0y2.0 mm,
  // kbirk*3
  const std::array<int, 4> G4runNumbers{105, 106, 107,
                                        108};  // with noise, horiz=-0.5 deg, x=-5.0y2.0 mm
  // corresponding energies (GeV)
  const std::array<int, 4> G4energies{20, 60, 100, 180};
  // corresponding energies as doubles
  std::array<double, 4> G4Denergies{};
  std::transform(G4energies.begin(), G4energies.end(), G4Denergies.begin(),
                 [](int value) { return static_cast<double>(value); });
  // energies errors (0)
  const std::array<double, 4> G4Erenergies{};
  // path to files + filename
  const std::string G4path("pionsimNtuple/physics_sps2023_run");
  // ratios of S energy central tower to all towers
  std::array<double, 4> G4ratioS{};
  // error on ratios of S energy central tower to all towers
  std::array<double, 4> G4erratioS{};
  // ratios of C energy central tower to all towers
  std::array<double, 4> G4ratioC{};
  // error on ratios of C energy central tower to all towers
  std::array<double, 4> G4erratioC{};
  // S energy all towers
  std::array<double, 4> G4Senergy{};
  // error S energy all towers
  std::array<double, 4> G4ErSenergy{};
  // C energy all towers
  std::array<double, 4> G4Cenergy{};
  // error C energy all towers
  std::array<double, 4> G4ErCenergy{};

  // loop through runs (Geant4 data only)
  for (int runno = 0; runno < G4energies.size(); runno++) {
    std::string pathFile(G4path + std::to_string(G4runNumbers[runno]) + ".root");
    std::cout << "Using file: " << pathFile << std::endl;
    TFile file{pathFile.c_str()};
    TTree* tree = static_cast<TTree*>(file.Get("Ftree"));
    bool isMC = true;
    outputAnalysis G4output = DoAnalysis(tree, G4runNumbers[runno], G4energies[runno], isMC);
    G4ratioS[runno] = G4output.RatioSenergy;
    G4erratioS[runno] = G4output.ErRatioSenergy;
    G4ratioC[runno] = G4output.RatioCenergy;
    G4erratioC[runno] = G4output.ErRatioCenergy;
    G4Senergy[runno] = G4output.Senergy;
    G4ErSenergy[runno] = G4output.ErSenergy;
    G4Cenergy[runno] = G4output.Cenergy;
    G4ErCenergy[runno] = G4output.ErCenergy;
  }

  TGraphErrors G4GrRatioS{4, G4Denergies.data(), G4ratioS.data(), G4Erenergies.data(),
                          G4erratioS.data()};
  //.data() converts from std::array to array
  G4GrRatioS.SetTitle("Central Tower S signal ratio;Beam energy (GeV);S ratio");
  G4GrRatioS.SetName("G4RatioS");

  TGraphErrors G4GrRatioC{4, G4Denergies.data(), G4ratioC.data(), G4Erenergies.data(),
                          G4erratioC.data()};
  //.data() converts from std::array to array
  G4GrRatioC.SetTitle("Central Tower C signal ratio;Beam energy (GeV);C ratio");
  G4GrRatioC.SetName("G4RatioC");

  TGraphErrors G4GrSenergy{4, G4Denergies.data(), G4Senergy.data(), G4Erenergies.data(),
                           G4ErSenergy.data()};
  G4GrSenergy.SetTitle("Scintillation;Beam energy (GeV);Signal S (GeV)");
  G4GrSenergy.SetName("G4Senergy");

  TGraphErrors G4GrCenergy{4, G4Denergies.data(), G4Cenergy.data(), G4Erenergies.data(),
                           G4ErCenergy.data()};
  G4GrCenergy.SetTitle("Cherenkov;Beam energy (GeV);Signal C (GeV)");
  G4GrCenergy.SetName("G4Cenergy");

  TCanvas* RatioCanvas = new TCanvas("RatioCanvas");
  GrRatioS.SetMarkerColor(kRed);
  GrRatioS.SetMarkerStyle(24);
  GrRatioS.SetLineColor(kRed);
  GrRatioS.GetHistogram()->SetMinimum(0.2);
  GrRatioS.GetHistogram()->SetMaximum(0.7);
  GrRatioS.Draw("APL");
  G4GrRatioS.SetMarkerColor(kRed);
  G4GrRatioS.SetMarkerStyle(25);
  G4GrRatioS.SetLineStyle(2);
  G4GrRatioS.SetLineColor(kRed);
  G4GrRatioS.Draw("SAME PL");
  GrRatioC.SetMarkerColor(kBlue);
  GrRatioC.SetMarkerStyle(24);
  GrRatioC.SetLineColor(kBlue);
  GrRatioC.Draw("SAME PL");
  G4GrRatioC.SetMarkerColor(kBlue);
  G4GrRatioC.SetMarkerStyle(25);
  G4GrRatioC.SetLineStyle(2);
  G4GrRatioC.Draw("SAME PL");

  FinalFile->cd();
  RatioCanvas->Write();
  G4GrRatioS.Write();
  G4GrRatioC.Write();
  G4GrSenergy.Write();
  G4GrCenergy.Write();

  CompareHistos(FinalFile);

  FinalFile->Close();
}

void CompareHistos(TFile* ComparisonFile)
{
  // TFile* ComparisonFile(TFile::Open("ComparisonPions.root", "RECREATE"));
  const std::array<int, 4> runNumbers{247, 246, 245, 244};  // experimental
  const std::array<int, 4> G4runNumbers{105, 106, 107, 108};  // simulation
  const std::array<int, 4> energies{20, 60, 100, 180};
  const std::string filename{"AnalysispionRun_"};
  const std::string G4filename{"G4AnalysispionRun_"};

  for (int i = 0; i < G4runNumbers.size(); i++) {
    TFile* expfile =
      TFile::Open((filename + std::to_string(runNumbers[i]) + ".root").c_str(), "READ");
    TFile* G4file =
      TFile::Open((G4filename + std::to_string(G4runNumbers[i]) + ".root").c_str(), "READ");

    // S-C tot histo
    auto H1Sene = (TH1F*)expfile->Get("H1Sene");
    H1Sene->Scale(1. / H1Sene->Integral());
    auto G4H1Sene = (TH1F*)G4file->Get("H1Sene");
    G4H1Sene->Scale(1. / G4H1Sene->Integral());
    TCanvas* c1 = new TCanvas((std::to_string(energies[i]) + "C1_" + std::to_string(i)).c_str());
    H1Sene->SetFillColor(kRed);
    H1Sene->SetLineColor(kRed);
    H1Sene->GetYaxis()->SetRangeUser(0, 0.07);
    H1Sene->Draw("HIST LF2");
    c1->Update();
    TPaveStats* stats = (TPaveStats*)c1->GetPrimitive("stats");
    stats->SetName("h1stats");
    stats->SetY1NDC(.8);
    stats->SetY2NDC(.9);
    stats->SetTextColor(kRed);
    G4H1Sene->SetLineColor(kBlue);
    G4H1Sene->Draw("SAMES");
    c1->Update();
    TPaveStats* stats2 = (TPaveStats*)c1->GetPrimitive("stats");
    stats2->SetName("h1stats2");
    stats2->SetY1NDC(.7);
    stats2->SetY2NDC(.8);
    stats2->SetTextColor(kBlue);
    c1->Update();
    auto legend = new TLegend(0.12, 0.8, 0.32, 0.9);
    legend->AddEntry(H1Sene, "Test-beam 2023", "f");
    legend->AddEntry(G4H1Sene, "Geant4", "f");
    legend->Draw();
    c1->Update();
    ComparisonFile->cd();
    c1->Write();
    // H1Cene
    auto H1Cene = (TH1F*)expfile->Get("H1Cene");
    H1Cene->Scale(1. / H1Cene->Integral());
    auto G4H1Cene = (TH1F*)G4file->Get("H1Cene");
    G4H1Cene->Scale(1. / G4H1Cene->Integral());
    TCanvas* c1C = new TCanvas((std::to_string(energies[i]) + "C1C_" + std::to_string(i)).c_str());
    H1Cene->SetFillColor(kRed);
    H1Cene->SetLineColor(kRed);
    H1Cene->GetYaxis()->SetRangeUser(0, 0.07);
    H1Cene->Draw("HIST LF2");
    c1C->Update();
    TPaveStats* statsC = (TPaveStats*)c1C->GetPrimitive("stats");
    statsC->SetName("h1stats");
    statsC->SetY1NDC(.8);
    statsC->SetY2NDC(.9);
    statsC->SetTextColor(kRed);
    G4H1Cene->SetLineColor(kBlue);
    G4H1Cene->Draw("SAMES");
    c1C->Update();
    TPaveStats* stats2C = (TPaveStats*)c1C->GetPrimitive("stats");
    stats2C->SetName("h1stats2");
    stats2C->SetY1NDC(.7);
    stats2C->SetY2NDC(.8);
    stats2C->SetTextColor(kBlue);
    c1C->Update();
    auto legendC = new TLegend(0.12, 0.8, 0.32, 0.9);
    legendC->AddEntry(H1Sene, "Test-beam 2023", "f");
    legendC->AddEntry(G4H1Sene, "Geant4", "f");
    legendC->Draw();
    c1C->Update();
    ComparisonFile->cd();
    c1C->Write();

    // S-C pmt histo
    auto H1PMTSene = (TH1F*)expfile->Get("H1PMTSene");
    H1PMTSene->Scale(1. / H1PMTSene->Integral());
    auto G4H1PMTSene = (TH1F*)G4file->Get("H1PMTSene");
    G4H1PMTSene->Scale(1. / G4H1PMTSene->Integral());
    TCanvas* c1PMT =
      new TCanvas((std::to_string(energies[i]) + "C1PMT_" + std::to_string(i)).c_str());
    H1PMTSene->SetFillColor(kRed);
    H1PMTSene->SetLineColor(kRed);
    H1PMTSene->GetYaxis()->SetRangeUser(0, 0.07);
    H1PMTSene->Draw("HIST LF2");
    c1PMT->Update();
    TPaveStats* statsPMT = (TPaveStats*)c1PMT->GetPrimitive("stats");
    statsPMT->SetName("h1stats");
    statsPMT->SetY1NDC(.8);
    statsPMT->SetY2NDC(.9);
    statsPMT->SetTextColor(kRed);
    G4H1PMTSene->SetLineColor(kBlue);
    G4H1PMTSene->Draw("SAMES");
    c1PMT->Update();
    TPaveStats* stats2PMT = (TPaveStats*)c1PMT->GetPrimitive("stats");
    stats2PMT->SetName("h1stats2");
    stats2PMT->SetY1NDC(.7);
    stats2PMT->SetY2NDC(.8);
    stats2PMT->SetTextColor(kBlue);
    c1PMT->Update();
    auto legendPMT = new TLegend(0.12, 0.8, 0.32, 0.9);
    legendPMT->AddEntry(H1PMTSene, "Test-beam 2023", "f");
    legendPMT->AddEntry(G4H1PMTSene, "Geant4", "f");
    legendPMT->Draw();
    c1PMT->Update();
    ComparisonFile->cd();
    c1PMT->Write();
    // H1CPMTene
    auto H1PMTCene = (TH1F*)expfile->Get("H1PMTCene");
    H1PMTCene->Scale(1. / H1PMTCene->Integral());
    auto G4H1PMTCene = (TH1F*)G4file->Get("H1PMTCene");
    G4H1PMTCene->Scale(1. / G4H1PMTCene->Integral());
    TCanvas* c1PMTC =
      new TCanvas((std::to_string(energies[i]) + "C1PMTC_" + std::to_string(i)).c_str());
    H1PMTCene->SetFillColor(kRed);
    H1PMTCene->SetLineColor(kRed);
    H1PMTCene->GetYaxis()->SetRangeUser(0, 0.07);
    H1PMTCene->Draw("HIST LF2");
    c1PMTC->Update();
    TPaveStats* statsPMTC = (TPaveStats*)c1PMTC->GetPrimitive("stats");
    statsPMTC->SetName("h1stats");
    statsPMTC->SetY1NDC(.8);
    statsPMTC->SetY2NDC(.9);
    statsPMTC->SetTextColor(kRed);
    G4H1PMTCene->SetLineColor(kBlue);
    G4H1PMTCene->Draw("SAMES");
    c1PMTC->Update();
    TPaveStats* stats2PMTC = (TPaveStats*)c1PMTC->GetPrimitive("stats");
    stats2PMTC->SetName("h1stats2");
    stats2PMTC->SetY1NDC(.7);
    stats2PMTC->SetY2NDC(.8);
    stats2PMTC->SetTextColor(kBlue);
    c1PMTC->Update();
    auto legendPMTC = new TLegend(0.12, 0.8, 0.32, 0.9);
    legendPMTC->AddEntry(H1Sene, "Test-beam 2023", "f");
    legendPMTC->AddEntry(G4H1Sene, "Geant4", "f");
    legendPMTC->Draw();
    c1PMTC->Update();
    ComparisonFile->cd();
    c1PMTC->Write();

    // S-C sipm histo
    auto H1sipmSene = (TH1F*)expfile->Get("H1SiPMSene");
    H1sipmSene->Scale(1. / H1sipmSene->Integral());
    auto G4H1sipmSene = (TH1F*)G4file->Get("H1SiPMSene");
    G4H1sipmSene->Scale(1. / G4H1sipmSene->Integral());
    TCanvas* c1sipm =
      new TCanvas((std::to_string(energies[i]) + "C1sipm_" + std::to_string(i)).c_str());
    H1sipmSene->SetFillColor(kRed);
    H1sipmSene->SetLineColor(kRed);
    H1sipmSene->GetYaxis()->SetRangeUser(0, 0.07);
    H1sipmSene->Draw("HIST LF2");
    c1sipm->Update();
    TPaveStats* statssipm = (TPaveStats*)c1sipm->GetPrimitive("stats");
    statssipm->SetName("h1stats");
    statssipm->SetY1NDC(.8);
    statssipm->SetY2NDC(.9);
    statssipm->SetTextColor(kRed);
    G4H1sipmSene->SetLineColor(kBlue);
    G4H1sipmSene->Draw("SAMES");
    c1sipm->Update();
    TPaveStats* stats2sipm = (TPaveStats*)c1sipm->GetPrimitive("stats");
    stats2sipm->SetName("h1stats2");
    stats2sipm->SetY1NDC(.7);
    stats2sipm->SetY2NDC(.8);
    stats2sipm->SetTextColor(kBlue);
    c1sipm->Update();
    auto legendsipm = new TLegend(0.12, 0.8, 0.32, 0.9);
    legendsipm->AddEntry(H1sipmSene, "Test-beam 2023", "f");
    legendsipm->AddEntry(G4H1sipmSene, "Geant4", "f");
    legendsipm->Draw();
    c1sipm->Update();
    ComparisonFile->cd();
    c1sipm->Write();
    // H1sipmCene
    auto H1sipmCene = (TH1F*)expfile->Get("H1SiPMCene");
    H1sipmCene->Scale(1. / H1sipmCene->Integral());
    auto G4H1sipmCene = (TH1F*)G4file->Get("H1SiPMCene");
    G4H1sipmCene->Scale(1. / G4H1sipmCene->Integral());
    TCanvas* c1sipmC =
      new TCanvas((std::to_string(energies[i]) + "C1SiPMC_" + std::to_string(i)).c_str());
    H1sipmCene->SetFillColor(kRed);
    H1sipmCene->SetLineColor(kRed);
    H1sipmCene->GetYaxis()->SetRangeUser(0, 0.07);
    H1sipmCene->Draw("HIST LF2");
    c1sipmC->Update();
    TPaveStats* statssipmC = (TPaveStats*)c1sipmC->GetPrimitive("stats");
    statssipmC->SetName("h1stats");
    statssipmC->SetY1NDC(.8);
    statssipmC->SetY2NDC(.9);
    statssipmC->SetTextColor(kRed);
    G4H1sipmCene->SetLineColor(kBlue);
    G4H1sipmCene->Draw("SAMES");
    c1sipmC->Update();
    TPaveStats* stats2sipmC = (TPaveStats*)c1sipmC->GetPrimitive("stats");
    stats2sipmC->SetName("h1stats2");
    stats2sipmC->SetY1NDC(.7);
    stats2sipmC->SetY2NDC(.8);
    stats2sipmC->SetTextColor(kBlue);
    c1sipmC->Update();
    auto legendsipmC = new TLegend(0.12, 0.8, 0.32, 0.9);
    legendsipmC->AddEntry(H1Sene, "Test-beam 2023", "f");
    legendsipmC->AddEntry(G4H1Sene, "Geant4", "f");
    legendsipmC->Draw();
    c1sipmC->Update();
    ComparisonFile->cd();
    c1sipmC->Write();
  }
}

outputAnalysis DoAnalysis(TTree* tree, const int& runno, const int& energy, const bool isMC = false)
{
  TFile* analysisFile;
  if (isMC == false) {
    analysisFile =
      TFile::Open(("AnalysispionRun_" + std::to_string(runno) + ".root").c_str(), "RECREATE");
  }
  else {
    analysisFile =
      TFile::Open(("G4AnalysispionRun_" + std::to_string(runno) + ".root").c_str(), "RECREATE");
  }
  // SiPM histos
  TH1F H1SiPMSene{"H1SiPMSene", "H1SiPMSene", 100, 0., static_cast<double>(energy)};
  H1SiPMSene.GetXaxis()->SetTitle("Signal (GeV)");
  TH1F H1SiPMCene{"H1SiPMCene", "H1SiPMCene", 100, 0., static_cast<double>(energy)};
  H1SiPMCene.GetXaxis()->SetTitle("Signal (GeV)");
  TH2F H2SiPMSCene{"H2SiPMSCene",
                   "H2SiPMSCene",
                   100,
                   0.,
                   static_cast<double>(energy),
                   100,
                   0.,
                   static_cast<double>(energy)};
  // PMT histos
  TH1F H1PMTSene{"H1PMTSene", "H1PMTSene", 100, 0., static_cast<double>(energy)};
  H1PMTSene.GetXaxis()->SetTitle("Signal (GeV)");
  TH1F H1PMTCene{"H1PMTCene", "H1PMTCene", 100, 0., static_cast<double>(energy)};
  H1PMTCene.GetXaxis()->SetTitle("Signal (GeV)");
  TH1F H1PMT1S{"H1PMT1S", "H1PMT1S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT2S{"H1PMT2S", "H1PMT2S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT3S{"H1PMT3S", "H1PMT3S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT4S{"H1PMT4S", "H1PMT4S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT5S{"H1PMT5S", "H1PMT5S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT6S{"H1PMT6S", "H1PMT6S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT7S{"H1PMT7S", "H1PMT7S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT8S{"H1PMT8S", "H1PMT8S", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT1C{"H1PMT1C", "H1PMT1C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT2C{"H1PMT2C", "H1PMT2C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT3C{"H1PMT3C", "H1PMT3C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT4C{"H1PMT4C", "H1PMT4C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT5C{"H1PMT5C", "H1PMT5C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT6C{"H1PMT6C", "H1PMT6C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT7C{"H1PMT7C", "H1PMT7C", 100, 0., static_cast<double>(energy) / 2.};
  TH1F H1PMT8C{"H1PMT8C", "H1PMT8C", 100, 0., static_cast<double>(energy) / 2.};
  // PMT+SiPM histos
  TH1F H1Sene{"H1Sene", "Scintillation Energy", 100, 0., 1.4 * static_cast<double>(energy)};
  H1Sene.GetXaxis()->SetTitle("Signal (GeV)");
  TH1F H1Cene{"H1Cene", "Cherenkov Energy", 100, 0., 1.4 * static_cast<double>(energy)};
  H1Cene.GetXaxis()->SetTitle("Signal (GeV)");
  TH2F H2SCene{"H2SCene",
               "H2SCene",
               100,
               0.,
               1.4 * static_cast<double>(energy),
               100,
               0.,
               1.4 * static_cast<double>(energy)};
  // Energy ratio histos
  TH1F H1RatioS{"H1RatioS", "H1RatioS", 100, 0., 1.3};
  TH1F H1RatioC{"H1RatioC", "H1RatioC", 100, 0., 1.3};
  // Geometry histos
  TH2F H2DWC1{"DWC1", "DWC1", 400, -20., 20., 400, -20, 20};
  TH2F H2DWC2{"DWC2", "DWC2", 400, -20., 20., 400, -20, 20};
  TH2F H2SiPMSbar{"SiPMSbar", "SiPMSbar", 400, -20., 40., 400, -20, 40};
  TH2F H2SiPMCbar{"SiPMCbar", "SiPMcbar", 400, -20., 40., 400, -20, 40};
  TH2F H2SXSiPMbarXDWC1{"XSiPMSbarXDWC1", "XSiPMSbarXDWC1", 400, -20., 40., 400, -20, 40};
  TH2F H2SXSiPMbarPMTS4{"XSiPMSbarPMTS4",           "XSiPMSbarPMTS4", 400, -20., 40., 100, 0.,
                        static_cast<double>(energy)};

  EventOut evtout{};
  EventOut* pevtout = &evtout;
  tree->SetBranchAddress("Events", &pevtout);

  int evtcounter = 0;
  // loop over events
  for (int evtno = 0; evtno < tree->GetEntries(); evtno++) {
    tree->GetEntry(evtno);
    // filter events
    if (isMC == false) {
      if (pevtout->TriggerMask == 6) continue;  // remove pedestal events
      if (!(utils::IsPionPsMu(pevtout->PShower, pevtout->MCounter))) continue;  // select pions
    }
    double DWC1pos[2] = {pevtout->XDWC1, pevtout->YDWC1};
    double DWC2pos[2] = {pevtout->XDWC2, pevtout->YDWC2};
    // if (!(utils::IsDWCradius(DWC1pos, 10.0, GetDWCoffset(1)))) continue;  // cut out-of-radius
    // DWC1 if (!(utils::IsDWCradius(DWC2pos, 10.0, GetDWCoffset(2)))) continue;  // cut
    // out-of-radius DWC2
    if (!(utils::IsDWCradius(DWC1pos, 10.0))) continue;
    if (!(utils::IsDWCradius(DWC2pos, 10.0))) continue;
    // cut events with non-interacting pions
    if (!(utils::HasPionInteracted(pevtout->totSiPMSene + pevtout->SPMTenergy, 4.0))) continue;

    evtcounter++;
    // fill histos of this run
    // energy histos
    H1SiPMSene.Fill(pevtout->totSiPMSene);
    H1SiPMCene.Fill(pevtout->totSiPMCene);
    H1PMT1S.Fill(pevtout->SPMT1);
    H1PMT2S.Fill(pevtout->SPMT2);
    H1PMT3S.Fill(pevtout->SPMT3);
    H1PMT4S.Fill(pevtout->SPMT4);
    H1PMT5S.Fill(pevtout->SPMT5);
    H1PMT6S.Fill(pevtout->SPMT6);
    H1PMT7S.Fill(pevtout->SPMT7);
    H1PMT8S.Fill(pevtout->SPMT8);
    H1PMT1C.Fill(pevtout->CPMT1);
    H1PMT2C.Fill(pevtout->CPMT2);
    H1PMT3C.Fill(pevtout->CPMT3);
    H1PMT4C.Fill(pevtout->CPMT4);
    H1PMT5C.Fill(pevtout->CPMT5);
    H1PMT6C.Fill(pevtout->CPMT6);
    H1PMT7C.Fill(pevtout->CPMT7);
    H1PMT8C.Fill(pevtout->CPMT8);
    H2SiPMSCene.Fill(pevtout->totSiPMSene, pevtout->totSiPMCene);
    H1PMTSene.Fill(pevtout->SPMTenergy);
    H1PMTCene.Fill(pevtout->CPMTenergy);
    H1Sene.Fill(pevtout->totSiPMSene + pevtout->SPMTenergy);
    H1Cene.Fill(pevtout->totSiPMCene + pevtout->CPMTenergy);
    H2SCene.Fill(pevtout->totSiPMSene + pevtout->SPMTenergy,
                 pevtout->totSiPMCene + pevtout->CPMTenergy);
    // Energy ratio histos
    H1RatioS.Fill(pevtout->totSiPMSene / (pevtout->totSiPMSene + pevtout->SPMTenergy));
    H1RatioC.Fill(pevtout->totSiPMCene / (pevtout->totSiPMCene + pevtout->CPMTenergy));
    // geometry histos
    H2DWC1.Fill(DWC1pos[0], DWC1pos[1]);
    H2DWC2.Fill(DWC2pos[0], DWC2pos[1]);
    double* sbar = utils::GetScinbar(pevtout->SiPMPheS);
    double* cbar = utils::GetCherbar(pevtout->SiPMPheC);
    H2SiPMSbar.Fill(sbar[0], sbar[1]);
    H2SiPMCbar.Fill(cbar[0], cbar[1]);
    H2SXSiPMbarXDWC1.Fill(DWC1pos[0], sbar[0]);
    H2SXSiPMbarPMTS4.Fill(DWC1pos[0], pevtout->SPMT4);
  }

  std::cout << "--->" << tree->GetEntries() << " evts, " << evtcounter
            << " selected, average SiPMene S: " << H1SiPMSene.GetMean()
            << " C: " << H1SiPMCene.GetMean() << " GeV" << std::endl;

  H1SiPMSene.Write();
  H1SiPMCene.Write();
  H1PMT1S.Write();
  H1PMT2S.Write();
  H1PMT3S.Write();
  H1PMT4S.Write();
  H1PMT5S.Write();
  H1PMT6S.Write();
  H1PMT7S.Write();
  H1PMT8S.Write();
  H1PMT1C.Write();
  H1PMT2C.Write();
  H1PMT3C.Write();
  H1PMT4C.Write();
  H1PMT5C.Write();
  H1PMT6C.Write();
  H1PMT7C.Write();
  H1PMT8C.Write();
  TCanvas SCanvas("SCanvas");
  SCanvas.Divide(3, 3);
  SCanvas.cd(9);
  H1PMT1S.Draw();
  SCanvas.cd(8);
  H1PMT2S.Draw();
  SCanvas.cd(7);
  H1PMT3S.Draw();
  SCanvas.cd(6);
  H1PMT4S.Draw();
  SCanvas.cd(5);
  H1SiPMSene.Draw();
  SCanvas.cd(4);
  H1PMT5S.Draw();
  SCanvas.cd(3);
  H1PMT6S.Draw();
  SCanvas.cd(2);
  H1PMT7S.Draw();
  SCanvas.cd(1);
  H1PMT8S.Draw();
  SCanvas.Write();
  TCanvas CCanvas("CCanvas");
  CCanvas.Divide(3, 3);
  CCanvas.cd(9);
  H1PMT1C.Draw();
  CCanvas.cd(8);
  H1PMT2C.Draw();
  CCanvas.cd(7);
  H1PMT3C.Draw();
  CCanvas.cd(6);
  H1PMT4C.Draw();
  CCanvas.cd(5);
  H1SiPMCene.Draw();
  CCanvas.cd(4);
  H1PMT5C.Draw();
  CCanvas.cd(3);
  H1PMT6C.Draw();
  CCanvas.cd(2);
  H1PMT7C.Draw();
  CCanvas.cd(1);
  H1PMT8C.Draw();
  CCanvas.Write();
  H2SiPMSCene.Write();
  H1PMTSene.Write();
  H1PMTCene.Write();
  H1Sene.Write();
  H1Cene.Write();
  H2SCene.Write();
  H1RatioS.Write();
  H1RatioC.Write();
  H2DWC1.Write();
  H2DWC2.Write();
  H2SiPMSbar.Write();
  H2SiPMCbar.Write();
  H2SXSiPMbarXDWC1.Write();
  H2SXSiPMbarPMTS4.Write();
  analysisFile->Close();

  outputAnalysis out{H1RatioS.GetMean(),      H1RatioS.GetMeanError(), H1RatioC.GetMean(),
                     H1RatioC.GetMeanError(), H1Sene.GetMean(),        H1Sene.GetMeanError(),
                     H1Cene.GetMean(),        H1Cene.GetMeanError()};
  return out;
}

//**************************************************
