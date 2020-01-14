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

using namespace ana;

// POT for n years
const double years = 3.5; // of POT
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

// True particle multiplicities
const Var kPipl      = SIMPLEVAR(dune.nipip);
const Var kPimin     = SIMPLEVAR(dune.nipim);
const Var kChargedPi = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip);
const Var kPi        = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);
const Var kP         = SIMPLEVAR(dune.nP);
const Var kOther     = SIMPLEVAR(dune.niother) + SIMPLEVAR(dune.nNucleus);

// HC
const Var kFHC = SIMPLEVAR(dune.isFHC);

const Binning binsParticle = Binning::Simple(6, -0.5, 5.5);

const HistAxis axTrueRecoPi("N_{#pi, true}", binsParticle, kPi,
			    "N_{#pi, reco}", binsParticle, kRecoPi);
const HistAxis axTrueRecoP("N_{P, true}", binsParticle, kP,
			   "N_{P, reco}", binsParticle, kRecoP);

void migrationMatrices(const char *outfile, 
		       const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v6/") 
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
  h2Pi->SetTitle("True vs. reco pion multiplicity in HPgTPC; N_{#pi, true}; N_{#pi, true}; Fraction of events");
  h2Pi_nuwro->SetTitle("True vs. reco pion multiplicity in HPgTPC (NuWro shifts); N_{#pi, true}; N_{#pi, true}; Fraction of events");
  h2Pi->Write("h2Pi");
  h2Pi_nuwro->Write("h2Pi_nuwro");
  // True vs. reco proton multiplicity
  TH2 *h2P       = predFhcP.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2P_nuwro = predFhcP.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  h2P->Scale(1./h2P->Integral());
  h2P_nuwro->Scale(1./h2P_nuwro->Integral());
  setHistAttr(h2P);
  setHistAttr(h2P_nuwro);
  h2P->SetTitle("True vs. reco proton multiplicity in HPgTPC; N_{P, true}; N_{P, true}; Fraction of events");
  h2P_nuwro->SetTitle("True vs. reco proton multiplicity in HPgTPC (NuWro shifts); N_{P, true}; N_{P, true}; Fraction of events");
  h2P->Write("h2P");
  h2P_nuwro->Write("h2P_nuwro");

  fout->Close();

} // migrationMatrices
