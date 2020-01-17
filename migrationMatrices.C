// migrationMatrices.C
// Makes migration/confusion matrices for HPgTPC study
// I.e. true vs. reconstructed number of pions
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
  h->SetLineWidth(2);
  h->SetOption("hist");
}

void setHistAttr(TH2 *h2) 
{
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->SetOption("colz");
}

std::string catName(const int cat)
{
  std::string name;
  if (cat==1) name="0#pi";
  else if (cat==2) name="1#pi^{#pm}";
  else if (cat==3) name="1#pi^{0}";
  else if (cat==4) name="2#pi";
  else if (cat==5) name="Other";
  return name;
}

using namespace ana;

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// Vars
// Reconstructed particle multiplicities
const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var kPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi        = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kRecoP         = SIMPLEVAR(dune.gastpc_pro_mult);
const Var kRecoOtherHad  = SIMPLEVAR(dune.gastpc_other_had_mult); 
const Var kRecoHad       = SIMPLEVAR(dune.gastpc_pi_pl_mult)+SIMPLEVAR(dune.gastpc_pi_min_mult)+SIMPLEVAR(dune.gastpc_pi_0_mult)+SIMPLEVAR(dune.gastpc_pro_mult)+SIMPLEVAR(dune.gastpc_other_had_mult);

// True particle multiplicities
const Var kPipl      = SIMPLEVAR(dune.nipip);
const Var kPimin     = SIMPLEVAR(dune.nipim);
const Var kChargedPi = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip);
const Var kPi        = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);
const Var kP         = SIMPLEVAR(dune.nP);
const Var kOther     = SIMPLEVAR(dune.niother) + SIMPLEVAR(dune.nNucleus);
const Var kHad       = SIMPLEVAR(dune.nipip)+SIMPLEVAR(dune.nipim)+SIMPLEVAR(dune.nipi0)+SIMPLEVAR(dune.nikp)+SIMPLEVAR(dune.nikm)+SIMPLEVAR(dune.nik0)+SIMPLEVAR(dune.nP)+SIMPLEVAR(dune.nN)/*+SIMPLEVAR(dune.niother)+SIMPLEVAR(dune.nNucleus)*/;

// CM's categories
const Var kRecoCategory({"dune.gastpc_pi_pl_mult", "dune.gastpc_pi_min_mult", "dune.gastpc_pi_0_mult"},
			[](const caf::StandardRecord* sr) {
			  int cat = -1;
			  int nChargedPi = sr->dune.gastpc_pi_pl_mult + sr->dune.gastpc_pi_min_mult;
			  int nPi = sr->dune.gastpc_pi_pl_mult + sr->dune.gastpc_pi_min_mult + sr->dune.gastpc_pi_0_mult;
			  if (nPi==0) cat=1;
			  else if (nChargedPi==1 && nPi==1) cat=2;
			  else if (sr->dune.gastpc_pi_0_mult==1 && nPi==1) cat=3;
			  else if (nPi==2) cat=4;
			  else cat=5;
			  return cat;
			});

const Var kTrueCategory({"dune.nipip", "dune.nipim", "dune.nipi0"},
			[](const caf::StandardRecord* sr) {
			  int cat = -1;
			  int nChargedPi = sr->dune.nipip + sr->dune.nipim;
			  int nPi = sr->dune.nipip + sr->dune.nipim + sr->dune.nipi0;
			  if (nPi==0) cat=1;
			  else if (nChargedPi==1 && nPi==1) cat=2;
			  else if (sr->dune.nipi0==1 && nPi==1) cat=3;
			  else if (nPi==2) cat=4;
			  else cat=5;
			  return cat;
			});

// HC
const Var kFHC = SIMPLEVAR(dune.isFHC);

const Binning binsHadron = Binning::Simple(10, -0.5, 9.5);
const HistAxis axTrueRecoHad("N_{had, true}", binsHadron, kHad,
			     "N_{had, reco}", binsHadron, kRecoHad);

const Binning binsParticle = Binning::Simple(6, -0.5, 5.5);
const HistAxis axTrueRecoPi("N_{#pi, true}", binsParticle, kPi,
			    "N_{#pi, reco}", binsParticle, kRecoPi);
const HistAxis axTrueRecoP ("N_{P, true}", binsParticle, kP,
			    "N_{P, reco}", binsParticle, kRecoP);

const Binning binsCategory = Binning::Simple(5, 0.5, 5.5);
const HistAxis axCategory("True category", binsCategory, kTrueCategory, 
			  "Reco category", binsCategory, kRecoCategory);

void migrationMatrices(const char *outfile, 
		       const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v8/") 
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

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", garDir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", garDir), kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  // True vs reco particle multiplicities
  NoOscPredictionGenerator genFhcPi(axTrueRecoPi, kPassND_FHC_NUMU && kIsTrueGasFV);
  PredictionInterp predFhcPi(systlist, 0, genFhcPi, loadersGArFHC);
  NoOscPredictionGenerator genFhcP(axTrueRecoP, kPassND_FHC_NUMU && kIsTrueGasFV);
  PredictionInterp predFhcP(systlist, 0, genFhcP, loadersGArFHC);
  // Reco vs true hadron multiplicities
  NoOscPredictionGenerator genFhcHad(axTrueRecoHad, kPassND_FHC_NUMU && kIsTrueGasFV);
  PredictionInterp predFhcHad(systlist, 0, genFhcHad, loadersGArFHC);
  // Reco vs true categories
  NoOscPredictionGenerator genFhcCat(axCategory, kPassND_FHC_NUMU && kIsTrueGasFV);
  PredictionInterp predFhcCat(systlist, 0, genFhcCat, loadersGArFHC);

  loadersGArFHC.Go();
  // loadersGArRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  // True vs. reco pion multiplicity
  TH2 *h2Pi       = predFhcPi.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Pi_nuwro = predFhcPi.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  h2Pi->Scale(1./h2Pi->Integral());
  h2Pi_nuwro->Scale(1./h2Pi_nuwro->Integral());
  setHistAttr(h2Pi);
  setHistAttr(h2Pi_nuwro);
  h2Pi->SetTitle("True vs. reco pion multiplicity in HPgTPC; N_{#pi, true}; N_{#pi, reco}; Fraction of events");
  h2Pi_nuwro->SetTitle("True vs. reco pion multiplicity in HPgTPC (NuWro shifts); N_{#pi, true}; N_{#pi, reco}; Fraction of events");
  h2Pi->Write("h2Pi");
  h2Pi_nuwro->Write("h2Pi_nuwro");
  // True vs. reco proton multiplicity
  TH2 *h2P       = predFhcP.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2P_nuwro = predFhcP.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  h2P->Scale(1./h2P->Integral());
  h2P_nuwro->Scale(1./h2P_nuwro->Integral());
  setHistAttr(h2P);
  setHistAttr(h2P_nuwro);
  h2P->SetTitle("True vs. reco proton multiplicity in HPgTPC; N_{P, true}; N_{P, reco}; Fraction of events");
  h2P_nuwro->SetTitle("True vs. reco proton multiplicity in HPgTPC (NuWro shifts); N_{P, true}; N_{P, reco}; Fraction of events");
  h2P->Write("h2P");
  h2P_nuwro->Write("h2P_nuwro");

  // Same as the above but we are going to use a different normalisation
  // Normalise by the total number of true events
  TH2 *h2Pi_confus       = predFhcPi.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Pi_nuwro_confus = predFhcPi.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  for (int binX=1; binX<h2Pi_confus->GetNbinsX()+1; binX++) {
    double colInt       = h2Pi_confus->Integral(binX, binX, 0, h2Pi_confus->GetNbinsY()+1);
    double colInt_nuwro = h2Pi_nuwro_confus->Integral(binX, binX, 0, h2Pi_nuwro_confus->GetNbinsY()+1);
    for (int binY=1; binY<h2Pi_confus->GetNbinsY()+1; binY++) {
      h2Pi_confus->SetBinContent(binX, binY, h2Pi_confus->GetBinContent(binX, binY) / colInt);
      h2Pi_nuwro_confus->SetBinContent(binX, binY, h2Pi_nuwro_confus->GetBinContent(binX, binY) / colInt_nuwro);
    }
  }

  TH2 *h2P_confus       = predFhcP.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2P_nuwro_confus = predFhcP.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  for (int binX=1; binX<h2P_confus->GetNbinsX()+1; binX++) {
    double colInt       = h2P_confus->Integral(binX, binX, 0, h2P_confus->GetNbinsY()+1);
    double colInt_nuwro = h2P_nuwro_confus->Integral(binX, binX, 0, h2P_nuwro_confus->GetNbinsY()+1);
    for (int binY=1; binY<h2P_confus->GetNbinsY()+1; binY++) {
      h2P_confus->SetBinContent(binX, binY, h2P_confus->GetBinContent(binX, binY) / colInt);
      h2P_nuwro_confus->SetBinContent(binX, binY, h2P_nuwro_confus->GetBinContent(binX, binY) / colInt_nuwro);
    }
  }
  setHistAttr(h2Pi_confus);
  setHistAttr(h2Pi_nuwro_confus);
  setHistAttr(h2P_confus);
  setHistAttr(h2P_nuwro_confus);
  h2Pi_confus->SetTitle("#pi confusion matrix in HPgTPC; N_{#pi, true}; N_{#pi, reco}; Probability");
  h2Pi_nuwro_confus->SetTitle("#pi confusion matrix in HPgTPC (NuWro shifts); N_{#pi, true}; N_{#pi, reco}; Probability");
  h2P_confus->SetTitle("Proton confusion matrix in HPgTPC; N_{P, true}; N_{P, reco}; Probability");
  h2P_nuwro_confus->SetTitle("Proton confusion matrix in HPgTPC (NuWro shifts); N_{P, true}; N_{P, reco}; Probability");
  h2Pi_confus->Write("h2Pi_confus");
  h2Pi_nuwro_confus->Write("h2Pi_nuwro_confus");
  h2P_confus->Write("h2P_confus");
  h2P_nuwro_confus->Write("h2P_nuwro_confus");

  // Reco vs. true categories
  TH2 *h2Cat       = predFhcCat.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Cat_nuwro = predFhcCat.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  h2Cat->Scale(1. / h2Cat->Integral());
  h2Cat_nuwro->Scale(1. / h2Cat_nuwro->Integral());
  setHistAttr(h2Cat);
  setHistAttr(h2Cat_nuwro);
  h2Cat->SetTitle("True vs. reco interaction category in HPgTPC; True category; Reco category; Fraction of events");
  h2Cat_nuwro->SetTitle("True vs. reco interaction category in HPgTPC (NuWro shifts); True category; Reco category; Fraction of events");
  for (int binX=1; binX<h2Cat->GetNbinsX()+1; binX++) {
    h2Cat->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat_nuwro->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat_nuwro->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
  }
  h2Cat->Write("h2Cat");
  h2Cat_nuwro->Write("h2Cat_nuwro");
  // Confusion matrices
  TH2 *h2Cat_confus       = predFhcCat.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Cat_nuwro_confus = predFhcCat.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  for (int binX=1; binX<h2Cat_confus->GetNbinsX()+1; binX++) {
    double colInt       = h2Cat_confus->Integral(binX, binX, 0, h2Cat_confus->GetNbinsY()+1);
    double colInt_nuwro = h2Cat_nuwro_confus->Integral(binX, binX,0, h2Cat_nuwro_confus->GetNbinsY()+1);
    h2Cat_confus->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat_nuwro_confus->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    for (int binY=1; binY<h2Cat_confus->GetNbinsY()+1; binY++) {
      h2Cat_confus->SetBinContent(binX, binY, h2Cat_confus->GetBinContent(binX, binY) / colInt);
      h2Cat_nuwro_confus->SetBinContent(binX, binY, h2Cat_nuwro_confus->GetBinContent(binX, binY) / colInt_nuwro);
      h2Cat_confus->GetYaxis()->SetBinLabel(binY, catName(binY).c_str());
      h2Cat_nuwro_confus->GetYaxis()->SetBinLabel(binY, catName(binY).c_str());
    }
  }
  h2Cat_confus->SetTitle("True vs. reco interaction category in HPgTPC; True category; Reco category; Probability");
  h2Cat_nuwro_confus->SetTitle("True vs. reco interaction category in HPgTPC (NuWro shifts); True category; Reco category; Probability");
  setHistAttr(h2Cat_confus);
  setHistAttr(h2Cat_nuwro_confus);
  h2Cat_confus->Write("h2Cat_confus");
  h2Cat_nuwro_confus->Write("h2Cat_nuwro_confus");

  // True vs. reco hadron multiplicity
  TH2 *h2Had       = predFhcHad.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Had_nuwro = predFhcHad.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  h2Had->SetTitle("True vs. reco hadron multiplicity in HPgTPC; N_{had, true}; N_{had, reco}; Fraction of events");
  h2Had_nuwro->SetTitle("True vs. reco hadron multiplicity in HPgTPC (NuWro shifts); N_{had, true}; N_{had, reco}; Fraction of events");
  setHistAttr(h2Had);
  setHistAttr(h2Had_nuwro);
  h2Had->Scale(1. / h2Had->Integral());
  h2Had_nuwro->Scale(1. / h2Had_nuwro->Integral());
  h2Had->Write("h2Had");
  h2Had_nuwro->Write("h2Had_nuwro");
  // Confusion matrices
  TH2 *h2Had_confus       = predFhcHad.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Had_nuwro_confus = predFhcHad.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  for (int binX=1; binX<h2Had_confus->GetNbinsX()+1; binX++) {
    double colInt       = h2Had_confus->Integral(binX, binX, 0, h2Had_confus->GetNbinsY()+1);
    double colInt_nuwro = h2Had_nuwro_confus->Integral(binX, binX,0, h2Had_nuwro_confus->GetNbinsY()+1);
    for (int binY=1; binY<h2Had_confus->GetNbinsY()+1; binY++) {
      h2Had_confus->SetBinContent(binX, binY, h2Had_confus->GetBinContent(binX, binY) / colInt);
      h2Had_nuwro_confus->SetBinContent(binX, binY, h2Had_nuwro_confus->GetBinContent(binX, binY) / colInt_nuwro);
    }
  }
  setHistAttr(h2Had_confus);
  setHistAttr(h2Had_nuwro_confus);
  h2Had_confus->SetTitle("Hadron multiplicity confusion matrix in HPgTPC; N_{had, true}; N_{had, reco}; Probability");
  h2Had_nuwro_confus->SetTitle("Hadron multiplicity confusion matrix in HPgTPC (NuWro shifts); N_{had, true}; N_{had, reco}; Probability");
  h2Had_confus->Write("h2Had_confus");
  h2Had_nuwro_confus->Write("h2Had_nuwro_confus");

  fout->Close();
} // migrationMatrices
