//*************************************************
// \file TB23Shape.C
// \brief: em shower shape analysis
// \author: Lorenzo Pezzotti (CERN EP-SFT-sim)
//          @lopezzot
// \start date: 24 July 2023
//**************************************************

// Usage: root TB23Shape.C
#include "PhysicsEvent.h"
#include "utils.h"
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>
#include <TVector2.h>
#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

ClassImp(EventOut);

void DoAnalysis(const string RunNo, const double& cutradius, const int& beamene, const bool isMC=false);

void TB23Shape()
{
  const double cutradius = 5.0;  // mm
  //Test-beam data
  //

  // Energy scan e+ 2.5deg vertical, 1.5deg orizoonthal, PShower in
  /*DoAnalysis("191", cutradius, 10); // x=-4.1 mm, y=1170 mm
  DoAnalysis("185", cutradius, 20);
  DoAnalysis("193", cutradius, 30);
  DoAnalysis("186", cutradius, 40);
  DoAnalysis("192", cutradius, 50);
  DoAnalysis("187", cutradius, 60);
  DoAnalysis("188", cutradius, 80);
  DoAnalysis("189", cutradius, 100);
  DoAnalysis("190", cutradius, 120);*/

  // Angular scan 20 GeV e+ 2.5deg vertical, y=-1170 mm
  /*DoAnalysis("114", cutradius, 20); // 4.5 deg 
  DoAnalysis("112", cutradius, 20);   // 3.5 deg
  DoAnalysis("119", cutradius, 20);   // 2.0 deg
  DoAnalysis("120", cutradius, 20);   // 1.5 deg
  DoAnalysis("122", cutradius, 20);   // 1.0 deg
  DoAnalysis("123", cutradius, 20);   // 0.5 deg
  DoAnalysis("126", cutradius, 20);   // 0.0 deg
  DoAnalysis("127", cutradius, 20);   // -0.5 deg
  DoAnalysis("129", cutradius, 20);   // -1.0 deg
  DoAnalysis("130", cutradius, 20);   // -2.5 deg*/

  // Angular scan 20 GeV e+ 0.0deg vertical y=1176 mm
  /*DoAnalysis("255", cutradius, 20);  // 5.0deg, x=-40.7 mm
  DoAnalysis("254", cutradius, 20);  // 4.0deg, x=-32.2 mm
  DoAnalysis("253", cutradius, 20);  // 3.0deg, x=-23.4 mm
  DoAnalysis("252", cutradius, 20);  // 2.0deg, x=-14.7 mm
  DoAnalysis("251", cutradius, 20);  // 1.5deg, x=-10.2 mm
  DoAnalysis("256", cutradius, 20);  // 1.0deg, x=-8.1 mm
  DoAnalysis("257", cutradius, 20);  // 0.5deg, x=-4.1 mm
  DoAnalysis("258", cutradius, 20);  // 0.0deg, x=0.1 mm
  DoAnalysis("266", cutradius, 20);  // -2.0deg, x=10.1 mm
  DoAnalysis("269", cutradius, 20);  // -5.0deg, x=29.1 mm*/

  //Comparison with and without preshower
  //vertical 2.5 deg, horiz 2.5 deg
  //(remember to remove preshower cut in analysis)
  /*DoAnalysis("102", cutradius, 20);  // preshower in
  DoAnalysis("163", cutradius, 20);  // preshower off 
  DoAnalysis("105", cutradius, 40);  // preshower in
  DoAnalysis("164", cutradius, 40);  // preshower off
  DoAnalysis("107", cutradius, 60);  // preshower in
  DoAnalysis("165", cutradius, 60);  // prehoser off*/

  // Simulated data
  //
  const bool isMC = true;
  //preliminary tests at vert angle 0.0
  /*DoAnalysis("0", cutradius, 20, isMC); // 0.0deg
  DoAnalysis("1", cutradius, 20, isMC); // 0.5deg
  DoAnalysis("2", cutradius, 20, isMC); // 1.0deg
  DoAnalysis("3", cutradius, 20, isMC); // 1.4deg
  DoAnalysis("4", cutradius, 20, isMC); // 1.8deg
  DoAnalysis("5", cutradius, 20, isMC); // -1.0deg
  DoAnalysis("6", cutradius, 20, isMC); // -1.8deg
  DoAnalysis("7", cutradius, 20, isMC); // -1.8deg, new seeds
  DoAnalysis("8", cutradius, 20, isMC); // -1.8deg, no attenuation
  DoAnalysis("9", cutradius, 20, isMC); // -1.8deg, 1.0 beam radius
  DoAnalysis("10", cutradius, 20, isMC);// -1.8deg, x-8mmy0, 1.0 beam radius*/
 
  //runcards_angscan, vert ang 0.0
  //DoAnalysis("11", cutradius, 20, isMC); //card1, 0.0deg, x0.1mmy0
  DoAnalysis("12", cutradius, 20, isMC); //card2, 0.5deg, x0.1mmy0
  //DoAnalysis("13", cutradius, 20, isMC); //card3, 0.5deg, x-0.1mmy0
  //DoAnalysis("14", cutradius, 20, isMC); //card4, -0.5deg, x0mmy0
  //DoAnalysis("15", cutradius, 20, isMC); //card5, -0.5deg, x-1.7mmy-1.6mm
  //DoAnalysis("16", cutradius, 20, isMC); //card6, -0.5deg, x-2.0mmy1.6mm
  //DoAnalysis("17", cutradius, 20, isMC); //card7, -0.7deg, x-2.0mmy1.6mm
  //DoAnalysis("18", cutradius, 20, isMC); //card8, -1.0deg, x-3.0mmy1.6mm
  //DoAnalysis("19", cutradius, 20, isMC); //card9, -1.0deg, x-4.5mmy1.6mm
}

void DoAnalysis(const string RunNo, const double& cutradius, const int& beamene, const bool isMC=false)
{
  // Input file
  std::string pathFile{};
  if(!isMC){ pathFile = "recoNtuple/physics_sps2023_run" + RunNo + ".root";}
  else { pathFile = "simNtuple/physics_sps2023_run" + RunNo + ".root";}
  std::cout << "Using file: " << pathFile << std::endl;

  TFile file{pathFile.c_str()};
  TTree* tree = static_cast<TTree*>(file.Get("Ftree"));
  EventOut evtout{};
  EventOut* pevtout = &evtout;
  tree->SetBranchAddress("Events", &pevtout);

  // Output file
  TFile* analysisFile = TFile::Open(
    ("AnalysisShape_" + RunNo + "_" + std::to_string(beamene) + "GeV.root").c_str(), "RECREATE");

  const int points = 25;  // em shape samplings
  int pitch = 1.;  // bin distance (mm)

  // Prepare for analysis (Scintillation)
  double radialprof[points]{};  // to be filled with radial profile
  double radialprofer[points]{};  // to be filled with rad prof error
  TH1F radprofh1[points];  // one TH1 per each point
  for (auto& n : radprofh1) {
    n.SetTitle("scinprof");
    n.SetName("scinprof");
    n.SetBins(1000, 0., 1.0);
  }
  double fibers[points]{};
  TH1F fibersh1[points];
  for (auto& n : fibersh1) {
    n.SetTitle("scinfibers");
    n.SetName("scinfibers");
    n.SetBins(100, 0., 100.);
  }
  double cumulativeprof[points]{};  // to be filled with cumul prof
  double radii[points]{};  // radii for x-axis
  double radiier[points]{};
  std::fill(radiier, radiier + points, 0.5);  // 0.5 uncertanty for x-bin
  TH2F Slateralh2{"Slateral", "Slateral", 25, 0., 25., 1000, 0., 0.5};

  // Prepare for analysis (Cherenkov)
  double cradialprof[points]{};
  double cradialprofer[points]{};
  TH1F cradprofh1[points];
  for (auto& n : cradprofh1) {
    n.SetTitle("cherprof");
    n.SetName("cherprof");
    n.SetBins(1000, 0., 1.0);
  }
  double cfibers[points]{};
  TH1F cfibersh1[points];
  for (auto& n : cfibersh1) {
    n.SetTitle("cherfibers");
    n.SetName("cherfibers");
    n.SetBins(100, 0., 100.);
  }
  double ccumulativeprof[points]{};
  TH2F Clateralh2{"Clateral", "Clateral", 25, 0., 25., 1000, 0., 0.5};

  int cutentries = 0;  // evets after cuts
  double totS = 0.;  // total event SiPM S signal
  double totC = 0.;  // total event SiPM C signal

  // SiPMs barycenter plots
  TH2F H2SiPMSbar{"SiPMSbar", "SiPMSbar", 200, -20., 40., 200, -20, 40};
  TH2F H2SiPMCbar{"SiPMCbar", "SiPMcbar", 200, -20., 40., 200, -20, 40};
  // DWC plots
  TH2F H2DWC1{"DWC1", "DWC1", 400, -20., 20., 400, -20, 20};
  TH2F H2DWC2{"DWC2", "DWC2", 400, -20., 20., 400, -20, 20};

  // Loop over events
  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    totS = 0.;
    totC = 0.;

    // Select positrons (avoid if MC data)
    if(!isMC){
      if (!(utils::IsPositronPsMu(pevtout->PShower, pevtout->MCounter))) continue;
    }
    if(isMC){
      if (!(utils::IsPsAboveCut(pevtout->PShower, 5))) continue;
    }

    // Get DWC info and cut over DWC radius
    double DWC1pos[2] = {pevtout->XDWC1, pevtout->YDWC1};
    double DWC2pos[2] = {pevtout->XDWC2, pevtout->YDWC2};
    if(!(utils::IsDWCradius(DWC1pos,cutradius))) continue;
    if(!(utils::IsDWCradius(DWC2pos,cutradius))) continue;

    cutentries += 1;

    H2DWC1.Fill(pevtout->XDWC1, pevtout->YDWC1);
    H2DWC2.Fill(pevtout->XDWC2, pevtout->YDWC2);

    for (const auto& n : pevtout->SiPMPheS) {
      totS += n;
    }
    for (const auto& n : pevtout->SiPMPheC) {
      totC += n;
    }
    auto sbar = utils::GetScinbar(pevtout->SiPMPheS);
    auto cbar = utils::GetCherbar(pevtout->SiPMPheC);
    H2SiPMSbar.Fill(sbar[0], sbar[1]);
    H2SiPMCbar.Fill(cbar[0], cbar[1]);

    // Loop over SiPMs
    for (unsigned int index = 0; index < 160; index++) {
      auto r = utils::Getdist(utils::ScinSiPMmap(index), sbar);
      // cout<<"entry "<<i<<" Sipm "<<index<<" distance "<<r<<" round "<<(int)r/pitch<<endl;
      auto cr = utils::Getdist(utils::CherSiPMmap(index), cbar);
      int newindex = (int)r / pitch;  // r-index scintillation
      int cnewindex = (int)cr / pitch;  // r-index Cherenkov
      Slateralh2.Fill(r, pevtout->SiPMPheS[index] / totS);  // Fill S scatterplot
      Clateralh2.Fill(cr, pevtout->SiPMPheC[index] / totC);  // Fill C scatterplot
      if (newindex < points) {
        radialprof[newindex] += pevtout->SiPMPheS[index] / totS;
        fibers[newindex] += 1;
      }
      if (cnewindex < points) {
        cradialprof[cnewindex] += pevtout->SiPMPheC[index] / totC;
        cfibers[cnewindex] += 1;
      }
    }  // end SiPMs loop

    for (unsigned int i = 0; i < points; i++) {
      radprofh1[i].Fill(radialprof[i]);
    }
    for (unsigned int i = 0; i < points; i++) {
      radialprof[i] = 0.;
    }
    for (unsigned int i = 0; i < points; i++) {
      fibersh1[i].Fill(fibers[i]);
    }
    for (unsigned int i = 0; i < points; i++) {
      fibers[i] = 0.;
    }
    for (unsigned int i = 0; i < points; i++) {
      cfibersh1[i].Fill(cfibers[i]);
    }
    for (unsigned int i = 0; i < points; i++) {
      cfibers[i] = 0.;
    }
    for (unsigned int i = 0; i < points; i++) {
      cradprofh1[i].Fill(cradialprof[i]);
    }
    for (unsigned int i = 0; i < points; i++) {
      cradialprof[i] = 0.;
    }
  }  // end event loop

  cout << "-->entries: " << tree->GetEntries() << " used: " << cutentries << endl;

  analysisFile->cd();

  // Write auxiliary detectors plots
  H2SiPMSbar.Write();
  H2SiPMCbar.Write();
  H2DWC1.Write();
  H2DWC2.Write();

  // Finalize H2 scatter plots
  Slateralh2.GetXaxis()->SetTitle("Distance from shower axis [mm]");
  Slateralh2.GetYaxis()->SetTitle("Percentage of total S SiPM signal in fiber");
  Clateralh2.GetXaxis()->SetTitle("Distance from shower axis [mm]");
  Clateralh2.GetYaxis()->SetTitle("Percentage of total C SiPM signal in fiber");
  Slateralh2.Write();
  Clateralh2.Write();
  // Extract lateral profiles from scatter plots
  auto sprof = Slateralh2.ProfileX();
  auto cprof = Clateralh2.ProfileX();
  sprof->SetTitle("lateralprof");
  sprof->GetXaxis()->SetTitle("Distance from shower axis [mm]");
  sprof->GetYaxis()->SetTitle("Percentage of total SiPM signal in fiber");
  sprof->SetName("lateralprof");
  sprof->SetMarkerStyle(20);
  sprof->SetMarkerColor(2);
  sprof->SetLineColor(2);
  sprof->Write();
  cprof->SetTitle("cherlateralprof");
  cprof->SetName("cherlateralprof");
  cprof->GetXaxis()->SetTitle("Distance from shower axis [mm]");
  cprof->GetYaxis()->SetTitle("Percentage of total SiPM signal in fiber");
  cprof->SetMarkerStyle(23);
  cprof->SetMarkerColor(4);
  cprof->SetLineColor(4);
  cprof->Write();

  // Finalize radial profiles
  for (unsigned int i = 0; i < points; i++) {
    radprofh1[i].Write();
    cradprofh1[i].Write();
    fibersh1[i].Write();

    radialprof[i] = radprofh1[i].GetMean();
    radialprofer[i] = radprofh1[i].GetMeanError();
    fibers[i] = fibersh1[i].GetMean();

    cradialprof[i] = cradprofh1[i].GetMean();
    cradialprofer[i] = cradprofh1[i].GetMeanError();
    cfibers[i] = cfibersh1[i].GetMean();
    radii[i] = pitch / 2. + pitch * i;
  }
  // Finalize cumulative profiles
  double counter = 0;
  int index = 0;
  for (auto& n : radialprof) {
    counter += n;
    cumulativeprof[index] = counter;
    index += 1;
  }
  double ccounter = 0;
  int cindex = 0;
  for (auto& n : cradialprof) {
    ccounter += n;
    ccumulativeprof[cindex] = ccounter;
    cindex += 1;
  }

  // Scintillation radial and cumulative profile graph
  TGraphErrors Gr2{points, radii, radialprof, radiier, radialprofer};
  Gr2.SetTitle("radialprof");
  Gr2.SetName("radialprof");
  Gr2.GetXaxis()->SetTitle("Distance from shower axis [mm]");
  Gr2.GetYaxis()->SetTitle("Percentage of SiPM signal in 1 mm thick radial shell");
  Gr2.SetMarkerStyle(20);
  Gr2.SetMarkerColor(2);
  Gr2.SetLineColor(2);
  Gr2.Write();
  TGraph Gr3{points, radii, cumulativeprof};
  Gr3.SetTitle("cumulativeprof");
  Gr3.SetName("cumulativeprof");
  Gr3.GetXaxis()->SetTitle("Radius of cylinder around shower axis [mm]");
  Gr3.GetYaxis()->SetTitle("Percentage of SiPM signal");
  Gr3.SetMarkerStyle(20);
  Gr3.SetMarkerColor(2);
  Gr3.SetLineColor(2);
  Gr3.Write();

  // Cherenkov radial and cumulative profile graph
  TGraphErrors CGr2{points, radii, cradialprof, radiier, cradialprofer};
  CGr2.SetTitle("cherradialprof");
  CGr2.SetName("cherradialprof");
  CGr2.GetXaxis()->SetTitle("Distance from shower axis [mm]");
  CGr2.GetYaxis()->SetTitle("Percentage of SiPM signal in 1 mm thick radial shell");
  CGr2.SetMarkerStyle(23);
  CGr2.SetMarkerColor(4);
  CGr2.SetLineColor(4);
  CGr2.Write();
  TGraph CGr3{points, radii, ccumulativeprof};
  CGr3.SetTitle("chercumulativeprof");
  CGr3.SetName("chercumulativeprof");
  CGr3.GetXaxis()->SetTitle("Radius of cylinder around shower axis [mm]");
  CGr3.GetYaxis()->SetTitle("Percentage of SiPM signal");
  CGr3.SetMarkerStyle(23);
  CGr3.SetMarkerColor(4);
  CGr3.SetLineColor(4);
  CGr3.Write();

  // Canvas with lateral profiles (C and S)
  TCanvas C1laterals{"lateralprofs", "lateralprofs", 600, 600};
  sprof->GetYaxis()->SetRangeUser(0., 0.09);
  sprof->SetTitle("");
  sprof->SetStats(0.);
  sprof->Draw("");
  cprof->Draw("same P");
  TLegend C1lateralslegend{1. - 0.18, 0.7, 1. - 0.61, 0.89};
  C1lateralslegend.AddEntry(
    sprof, ("CERN SPS: Scintillation " + to_string(beamene) + " GeV e+, Run " + RunNo).c_str(),
    "ep");
  C1lateralslegend.AddEntry(
    cprof, ("CERN SPS: Cherenkov " + to_string(beamene) + " GeV e+, Run " + RunNo).c_str(), "ep");
  C1lateralslegend.Draw("same");
  C1laterals.SetLeftMargin(0.15);
  C1laterals.Write();

  // Canvas with radial profiles (C and S)
  TCanvas C1radials{"radialprofs", "radialprofs", 600, 600};
  Gr2.GetHistogram()->SetMinimum(0.);
  Gr2.GetHistogram()->SetMaximum(0.18);
  Gr2.SetTitle("");
  Gr2.Draw("APL");
  CGr2.Draw("same PL");
  TLegend C1radialslegend{1. - 0.18, 0.7, 1. - 0.61, 0.89};
  C1radialslegend.AddEntry(
    &Gr2, ("CERN SPS: Scintillation " + to_string(beamene) + " GeV e+, Run " + RunNo).c_str(),
    "ep");
  C1radialslegend.AddEntry(
    &CGr2, ("CERN SPS: Cherenkov " + to_string(beamene) + " GeV e+, Run " + RunNo).c_str(), "ep");
  C1radialslegend.Draw("same");
  C1radials.SetLeftMargin(0.15);
  C1radials.Write();

  // Canvas with cumulative profiles (C and S)
  TCanvas C1cumulatives{"cumulativeprofs", "cumulativeprofs", 600, 600};
  Gr3.GetHistogram()->SetMinimum(0.);
  Gr3.GetHistogram()->SetMaximum(1.1);
  Gr3.SetTitle("");
  Gr3.Draw("APL");
  CGr3.Draw("same PL");
  TLegend C1cumulativelegend{1. - 0.18, 0.7, 1. - 0.61, 0.89};
  C1cumulativelegend.AddEntry(
    &Gr3, ("CERN SPS: Scintillation " + std::to_string(beamene) + " GeV e+, Run" + RunNo).c_str(),
    "ep");
  C1cumulativelegend.AddEntry(
    &CGr3, ("CERN SPS: Cherenkov " + std::to_string(beamene) + " GeV e+, Run " + RunNo).c_str(),
    "ep");
  C1cumulativelegend.Draw("same");
  C1cumulatives.SetLeftMargin(0.15);
  C1cumulatives.Write();

  analysisFile->Close();
  exit(0);
}

//**************************************************
