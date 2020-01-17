// pionChargeSep.C
// Similar studies as before but separating pions by charge
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Systs/XSecSysts.h"
#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/common_fit_definitions.h"
// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"

#include "Utilities/rootlogon.C"

void setHistAttr(TH1 *h) 
{
  h->GetXaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(1.1);
  h->SetMarkerStyle(8);
  h->SetLineWidth(2);
  //  h->SetOption("hist");
}

void setHistAttr(TH2 *h2) 
{
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->GetXaxis()->SetTitleOffset(1.1);
  h2->GetYaxis()->SetTitleOffset(1.1);
  h2->GetZaxis()->SetTitleOffset(1.1);
  h2->SetMarkerStyle(8);
  h2->SetOption("colz");
}

// Sets bin content to zero if true MC stats are low
void suppressLowStats(TH1 *h, const double inputPOT, const double truePOT)
{
  for (int i=1; i<h->GetNbinsX()+1; i++) {
    if (h->GetBinContent(i) * (truePOT/inputPOT) < 100.) {
      // Set bin content to 0
      h->SetBinContent(i, 0);
      h->SetBinError(i, 0);
    }
  }
}

TH1D* ratioSuppressLowStats(TH1 *hgenie, TH1* hnuwro, const char* name, const char* title, 
			    const double inputPOT, const double truePOT)
{
  assert(hgenie->GetNbinsX() == hnuwro->GetNbinsX());
  TH1D *hout = new TH1D(name, title, hgenie->GetNbinsX(), 0., 3.);
  for (int i=1; i<hgenie->GetNbinsX()+1; i++) {
    if (hgenie->GetBinContent(i) * (truePOT/inputPOT) < 100.) {
      // Set bin content to 0
      hgenie->SetBinContent(i, 0);
      hgenie->SetBinError(i, 0);
      hnuwro->SetBinContent(i, 0);
      hnuwro->SetBinError(i, 0);
    }
  }
  hout->Divide(hnuwro, hgenie, 1., 1.);
  setHistAttr(hout);
  return hout;
}

using namespace ana;

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// True MC files POT
const double fhcPOT = 1.9342e20;
const double rhcPOT = 3.9302e20;

// Vars
const Var kRecoPipl  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var kRecoPimin = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi0   = SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi        = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kP             = SIMPLEVAR(dune.nP);
const Var kRecoP         = SIMPLEVAR(dune.gastpc_pro_mult);
const Var kRecoOtherHad  = SIMPLEVAR(dune.gastpc_other_had_mult); 
const Var kPi   = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);
const Var kFHC   = SIMPLEVAR(dune.isFHC);
const double mmu = 0.10566; // GeV/c^2

// Reco Q2
const Var kRecoQ2({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		  [](const caf::StandardRecord* sr) {
		    double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		    double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		    return q2;
		  });
// Reco W
const Var kRecoW({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		 [](const caf::StandardRecord* sr) {
		   double w = 0;
		   double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		   double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		   w = TMath::Sqrt(-q2 + 2 * 0.939 * (sr->dune.Ev_reco-sr->dune.Elep_reco) + 0.939*0.939);
		   return w;
		 });

const Binning simpleBins   = Binning::Simple(30, 0, 3);
const Binning binsQ2 = Binning::Simple(15, 0, 3);
const Binning binsW  = Binning::Simple(15, 0, 3);
const Binning binsParticle = Binning::Simple(6, -0.5, 5.5);
const Binning binsTheta = Binning::Simple(30, 0, 1.6);

const HistAxis axRecoQ2("Q^{2}_{reco} / GeV^{2}", binsQ2, kRecoQ2);
const HistAxis axRecoW("W_{reco} / GeV", binsW, kRecoW);

void pionChargeSep(const char *outfile, const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v8/") 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();
 
  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = {};
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }

  assert(systlist.size()==1);
  SystShifts nuwroShift(systlist.at(0), 1);

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", garDir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", garDir), kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  // FHC
  NoOscPredictionGenerator genQ2_1Pipl(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPipl==1 && kRecoPi==1);
  NoOscPredictionGenerator genQ2_1Pimin(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPimin==1 && kRecoPi==1);
  NoOscPredictionGenerator genQ2_1Pi0(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi0==1 && kRecoPi==1);
  PredictionInterp predQ2_1Pipl(systlist, 0, genQ2_1Pipl, loadersGArFHC);
  PredictionInterp predQ2_1Pimin(systlist, 0, genQ2_1Pimin, loadersGArFHC);
  PredictionInterp predQ2_1Pi0(systlist, 0, genQ2_1Pi0, loadersGArFHC);

  NoOscPredictionGenerator genW_1Pipl(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPipl==1 && kRecoPi==1);
  NoOscPredictionGenerator genW_1Pimin(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPimin==1 && kRecoPi==1);
  NoOscPredictionGenerator genW_1Pi0(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi0==1 && kRecoPi==1);
  PredictionInterp predW_1Pipl(systlist, 0, genW_1Pipl, loadersGArFHC);
  PredictionInterp predW_1Pimin(systlist, 0, genW_1Pimin, loadersGArFHC);
  PredictionInterp predW_1Pi0(systlist, 0, genW_1Pi0, loadersGArFHC);
  // RHC
  NoOscPredictionGenerator genQ2_1Pipl_rhc(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPipl==1 && kRecoPi==1);
  NoOscPredictionGenerator genQ2_1Pimin_rhc(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPimin==1 && kRecoPi==1);
  NoOscPredictionGenerator genQ2_1Pi0_rhc(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi0==1 && kRecoPi==1);
  PredictionInterp predQ2_1Pipl_rhc(systlist, 0, genQ2_1Pipl, loadersGArRHC);
  PredictionInterp predQ2_1Pimin_rhc(systlist, 0, genQ2_1Pimin, loadersGArRHC);
  PredictionInterp predQ2_1Pi0_rhc(systlist, 0, genQ2_1Pi0, loadersGArRHC);

  NoOscPredictionGenerator genW_1Pipl_rhc(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPipl==1 && kRecoPi==1);
  NoOscPredictionGenerator genW_1Pimin_rhc(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPimin==1 && kRecoPi==1);
  NoOscPredictionGenerator genW_1Pi0_rhc(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi0==1 && kRecoPi==1);
  PredictionInterp predW_1Pipl_rhc(systlist, 0, genW_1Pipl, loadersGArRHC);
  PredictionInterp predW_1Pimin_rhc(systlist, 0, genW_1Pimin, loadersGArRHC);
  PredictionInterp predW_1Pi0_rhc(systlist, 0, genW_1Pi0, loadersGArRHC);

  loadersGArFHC.Go();
  loadersGArRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  // FHC
  // Q2
  TH1 *hQ2_1Pipl_genie = predQ2_1Pipl.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pipl_nuwro = predQ2_1Pipl.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pipl_genie);
  setHistAttr(hQ2_1Pipl_nuwro);
  hQ2_1Pipl_genie->SetTitle("1#pi^{+} (GENIE)");
  hQ2_1Pipl_nuwro->SetTitle("1#pi^{+} (NuWro)");
  hQ2_1Pipl_genie->Write("hQ2_1Pipl_genie");
  hQ2_1Pipl_nuwro->Write("hQ2_1Pipl_nuwro");
  TH1 *hrQ2_1Pipl = ratioSuppressLowStats(hQ2_1Pipl_genie, hQ2_1Pipl_nuwro, "hrQ2_1Pipl", Form("NuWro/GENIE for HPgTPC, 1#pi^{+}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrQ2_1Pipl->Write();

  TH1 *hQ2_1Pimin_genie = predQ2_1Pimin.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pimin_nuwro = predQ2_1Pimin.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pimin_genie);
  setHistAttr(hQ2_1Pimin_nuwro);
  hQ2_1Pimin_genie->SetTitle("1#pi^{-} (GENIE)");
  hQ2_1Pimin_nuwro->SetTitle("1#pi^{-} (NuWro)");
  hQ2_1Pimin_genie->Write("hQ2_1Pimin_genie");
  hQ2_1Pimin_nuwro->Write("hQ2_1Pimin_nuwro");
  TH1 *hrQ2_1Pimin = ratioSuppressLowStats(hQ2_1Pimin_genie, hQ2_1Pimin_nuwro, "hrQ2_1Pimin", Form("NuWro/GENIE for HPgTPC, 1#pi^{-}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrQ2_1Pimin->Write();

  TH1 *hQ2_1Pi0_genie = predQ2_1Pi0.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pi0_nuwro = predQ2_1Pi0.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pi0_genie);
  setHistAttr(hQ2_1Pi0_nuwro);
  hQ2_1Pi0_genie->SetTitle("1#pi^{0} (GENIE)");
  hQ2_1Pi0_nuwro->SetTitle("1#pi^{0} (NuWro)");
  hQ2_1Pi0_genie->Write("hQ2_1Pi0_genie");
  hQ2_1Pi0_nuwro->Write("hQ2_1Pi0_nuwro");
  TH1 *hrQ2_1Pi0 = ratioSuppressLowStats(hQ2_1Pi0_genie, hQ2_1Pi0_nuwro, "hrQ2_1Pi0", Form("NuWro/GENIE for HPgTPC, 1#pi^{0}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrQ2_1Pi0->Write();

  // W
  TH1 *hW_1Pipl_genie = predW_1Pipl.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pipl_nuwro = predW_1Pipl.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pipl_genie);
  setHistAttr(hW_1Pipl_nuwro);
  hW_1Pipl_genie->SetTitle("1#pi^{+} (GENIE)");
  hW_1Pipl_nuwro->SetTitle("1#pi^{+} (NuWro)");
  hW_1Pipl_genie->Write("hW_1Pipl_genie");
  hW_1Pipl_nuwro->Write("hW_1Pipl_nuwro");
  TH1 *hrW_1Pipl = ratioSuppressLowStats(hW_1Pipl_genie, hW_1Pipl_nuwro, "hrW_1Pipl", Form("NuWro/GENIE for HPgTPC, 1#pi^{+}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrW_1Pipl->Write();

  TH1 *hW_1Pimin_genie = predW_1Pimin.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pimin_nuwro = predW_1Pimin.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pimin_genie);
  setHistAttr(hW_1Pimin_nuwro);
  suppressLowStats(hW_1Pimin_genie, pot_nd, fhcPOT);
  suppressLowStats(hW_1Pimin_nuwro, pot_nd, fhcPOT);
  hW_1Pimin_genie->SetTitle("1#pi^{-} (GENIE)");
  hW_1Pimin_nuwro->SetTitle("1#pi^{-} (NuWro)");
  hW_1Pimin_genie->Write("hW_1Pimin_genie");
  hW_1Pimin_nuwro->Write("hW_1Pimin_nuwro");
  TH1 *hrW_1Pimin = ratioSuppressLowStats(hW_1Pimin_genie, hW_1Pimin_nuwro, "hrW_1Pimin", Form("NuWro/GENIE for HPgTPC, 1#pi^{-}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrW_1Pimin->Write();

  TH1 *hW_1Pi0_genie = predW_1Pi0.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pi0_nuwro = predW_1Pi0.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pi0_genie);
  setHistAttr(hW_1Pi0_nuwro);
  suppressLowStats(hW_1Pi0_genie, pot_nd, fhcPOT);
  suppressLowStats(hW_1Pi0_nuwro, pot_nd, fhcPOT);
  hW_1Pi0_genie->SetTitle("1#pi^{0} (GENIE)");
  hW_1Pi0_nuwro->SetTitle("1#pi^{0} (NuWro)");
  hW_1Pi0_genie->Write("hW_1Pi0_genie");
  hW_1Pi0_nuwro->Write("hW_1Pi0_nuwro");
  TH1 *hrW_1Pi0 = ratioSuppressLowStats(hW_1Pi0_genie, hW_1Pi0_nuwro, "hrW_1Pi0", Form("NuWro/GENIE for HPgTPC, 1#pi^{0}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, fhcPOT);
  hrW_1Pi0->Write();

  THStack *hsQ2 = new THStack("hsQ2", Form("Reconstructed Q^{2} for HPgTPC: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro / GENIE", years));
  hrQ2_1Pipl->SetLineColor(kBlue);
  hrQ2_1Pipl->SetMarkerColor(kBlue);
  hrQ2_1Pi0->SetLineColor(kBlack);
  hrQ2_1Pi0->SetMarkerColor(kBlack);
  hrQ2_1Pimin->SetLineColor(kRed);
  hrQ2_1Pimin->SetMarkerColor(kRed);
  hsQ2->Add(hrQ2_1Pipl);
  hsQ2->Add(hrQ2_1Pimin);
  hsQ2->Add(hrQ2_1Pi0);
  hsQ2->Write();

  THStack *hsW = new THStack("hsW", Form("Reconstructed W for HPgTPC: %.2g years POT; W_{reco} / GeV; NuWro / GENIE", years));
  hrW_1Pipl->SetLineColor(kBlue);
  hrW_1Pipl->SetMarkerColor(kBlue);
  hrW_1Pi0->SetLineColor(kBlack);
  hrW_1Pi0->SetMarkerColor(kBlack);
  hrW_1Pimin->SetLineColor(kRed);
  hrW_1Pimin->SetMarkerColor(kRed);
  hsW->Add(hrW_1Pipl);
  hsW->Add(hrW_1Pimin);
  hsW->Add(hrW_1Pi0);
  hsW->Write();

  // RHC
  // Q2
  TH1 *hQ2_1Pipl_genie_rhc = predQ2_1Pipl_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pipl_nuwro_rhc = predQ2_1Pipl_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pipl_genie_rhc);
  setHistAttr(hQ2_1Pipl_nuwro_rhc);
  hQ2_1Pipl_genie_rhc->SetTitle("1#pi^{+} (GENIE)");
  hQ2_1Pipl_nuwro_rhc->SetTitle("1#pi^{+} (NuWro)");
  hQ2_1Pipl_genie_rhc->Write("hQ2_1Pipl_genie_rhc");
  hQ2_1Pipl_nuwro_rhc->Write("hQ2_1Pipl_nuwro_rhc");
  TH1 *hrQ2_1Pipl_rhc = ratioSuppressLowStats(hQ2_1Pipl_genie_rhc, hQ2_1Pipl_nuwro_rhc, "hrQ2_1Pipl_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{+}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, rhcPOT);
  hrQ2_1Pipl_rhc->Write();

  TH1 *hQ2_1Pimin_genie_rhc = predQ2_1Pimin_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pimin_nuwro_rhc = predQ2_1Pimin_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pimin_genie_rhc);
  setHistAttr(hQ2_1Pimin_nuwro_rhc);
  hQ2_1Pimin_genie_rhc->SetTitle("1#pi^{-} (GENIE)");
  hQ2_1Pimin_nuwro_rhc->SetTitle("1#pi^{-} (NuWro)");
  hQ2_1Pimin_genie_rhc->Write("hQ2_1Pimin_genie_rhc");
  hQ2_1Pimin_nuwro_rhc->Write("hQ2_1Pimin_nuwro_rhc");
  TH1 *hrQ2_1Pimin_rhc = ratioSuppressLowStats(hQ2_1Pimin_genie_rhc, hQ2_1Pimin_nuwro_rhc, "hrQ2_1Pimin_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{-}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, rhcPOT);
  hrQ2_1Pimin_rhc->Write();

  TH1 *hQ2_1Pi0_genie_rhc = predQ2_1Pi0_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hQ2_1Pi0_nuwro_rhc = predQ2_1Pi0_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hQ2_1Pi0_genie_rhc);
  setHistAttr(hQ2_1Pi0_nuwro_rhc);
  hQ2_1Pi0_genie_rhc->SetTitle("1#pi^{0} (GENIE)");
  hQ2_1Pi0_nuwro_rhc->SetTitle("1#pi^{0} (NuWro)");
  hQ2_1Pi0_genie_rhc->Write("hQ2_1Pi0_genie_rhc");
  hQ2_1Pi0_nuwro_rhc->Write("hQ2_1Pi0_nuwro_rhc");
  TH1 *hrQ2_1Pi0_rhc = ratioSuppressLowStats(hQ2_1Pi0_genie_rhc, hQ2_1Pi0_nuwro_rhc, "hrQ2_1Pi0_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{0}: %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE", years), pot_nd, rhcPOT); 
  hrQ2_1Pi0_rhc->Write();

  // W
  TH1 *hW_1Pipl_genie_rhc = predW_1Pipl_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pipl_nuwro_rhc = predW_1Pipl_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pipl_genie_rhc);
  setHistAttr(hW_1Pipl_nuwro_rhc);
  hW_1Pipl_genie_rhc->SetTitle("1#pi^{+} (GENIE)");
  hW_1Pipl_nuwro_rhc->SetTitle("1#pi^{+} (NuWro)");
  hW_1Pipl_genie_rhc->Write("hW_1Pipl_genie_rhc");
  hW_1Pipl_nuwro_rhc->Write("hW_1Pipl_nuwro_rhc");
  TH1 *hrW_1Pipl_rhc = ratioSuppressLowStats(hW_1Pipl_genie_rhc, hW_1Pipl_nuwro_rhc, "hrW_1Pipl_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{+}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, rhcPOT);
  hrW_1Pipl_rhc->Write();

  TH1 *hW_1Pimin_genie_rhc = predW_1Pimin_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pimin_nuwro_rhc = predW_1Pimin_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pimin_genie_rhc);
  setHistAttr(hW_1Pimin_nuwro_rhc);
  hW_1Pimin_genie_rhc->SetTitle("1#pi^{-} (GENIE)");
  hW_1Pimin_nuwro_rhc->SetTitle("1#pi^{-} (NuWro)");
  hW_1Pimin_genie_rhc->Write("hW_1Pimin_genie_rhc");
  hW_1Pimin_nuwro_rhc->Write("hW_1Pimin_nuwro_rhc");
  TH1 *hrW_1Pimin_rhc = ratioSuppressLowStats(hW_1Pimin_genie_rhc, hW_1Pimin_nuwro_rhc, "hrW_1Pimin_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{-}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, rhcPOT);
  hrW_1Pimin_rhc->Write();

  TH1 *hW_1Pi0_genie_rhc = predW_1Pi0_rhc.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hW_1Pi0_nuwro_rhc = predW_1Pi0_rhc.PredictSyst(0, nuwroShift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hW_1Pi0_genie_rhc);
  setHistAttr(hW_1Pi0_nuwro_rhc);
  hW_1Pi0_genie_rhc->SetTitle("1#pi^{0} (GENIE)");
  hW_1Pi0_nuwro_rhc->SetTitle("1#pi^{0} (NuWro)");
  hW_1Pi0_genie_rhc->Write("hW_1Pi0_genie_rhc");
  hW_1Pi0_nuwro_rhc->Write("hW_1Pi0_nuwro_rhc");
  TH1 *hrW_1Pi0_rhc = ratioSuppressLowStats(hW_1Pi0_genie_rhc, hW_1Pi0_nuwro_rhc, "hrW_1Pi0_rhc", Form("NuWro/GENIE for HPgTPC, 1#pi^{0}: %.2g years POT; W_{reco} / GeV; NuWro/GENIE", years), pot_nd, rhcPOT);
  hrW_1Pi0_rhc->Write();

  THStack *hsQ2_rhc = new THStack("hsQ2_rhc", Form("Reconstructed Q^{2} for HPgTPC (RHC): %.2g years POT; Q^{2}_{reco} / GeV^{2}; NuWro / GENIE", years));
  hrQ2_1Pipl_rhc->SetLineColor(kBlue);
  hrQ2_1Pipl_rhc->SetMarkerColor(kBlue);
  hrQ2_1Pi0_rhc->SetLineColor(kBlack);
  hrQ2_1Pi0_rhc->SetMarkerColor(kBlack);
  hrQ2_1Pimin_rhc->SetLineColor(kRed);
  hrQ2_1Pimin_rhc->SetMarkerColor(kRed);
  hsQ2_rhc->Add(hrQ2_1Pipl_rhc);
  hsQ2_rhc->Add(hrQ2_1Pimin_rhc);
  hsQ2_rhc->Add(hrQ2_1Pi0_rhc);
  hsQ2_rhc->Write();

  THStack *hsW_rhc = new THStack("hsW_rhc", Form("Reconstructed W for HPgTPC (RHC): %.2g years POT; W_{reco} / GeV; NuWro / GENIE", years));
  hrW_1Pipl_rhc->SetLineColor(kBlue);
  hrW_1Pipl_rhc->SetMarkerColor(kBlue);
  hrW_1Pi0_rhc->SetLineColor(kBlack);
  hrW_1Pi0_rhc->SetMarkerColor(kBlack);
  hrW_1Pimin_rhc->SetLineColor(kRed);
  hrW_1Pimin_rhc->SetMarkerColor(kRed);
  hsW_rhc->Add(hrW_1Pipl_rhc);
  hsW_rhc->Add(hrW_1Pimin_rhc);
  hsW_rhc->Add(hrW_1Pi0_rhc);
  hsW_rhc->Write();

  fout->Close();
  delete fout;
} // pionChargeSep
