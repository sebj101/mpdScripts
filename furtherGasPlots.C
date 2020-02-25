// furtherGasPlots.C
// Plots showing things like the relationship between Q2 reco and true
// CAFAna includes
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
#include "TLine.h"

#include "Utilities/rootlogon.C"

using namespace ana;
const int nCats = 6;

void setHistAttr(TH1 *h) 
{
  h->GetXaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->SetMarkerStyle(8);
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
  h2->GetZaxis()->SetTitleOffset(1.);
  h2->SetOption("colz");
}

std::string catName(const int cat)
{
  std::string name;
  if (cat==1) name="0#pi";
  else if (cat==2) name="1#pi^{-}";
  else if (cat==3) name="1#pi^{+}";
  else if (cat==4) name="1#pi^{0}";
  else if (cat==5) name="2#pi";
  else if (cat==6) name=">2#pi";
  return name;
}

void doHistStuff(TH1* h, const char* name, const char* title)
{
  setHistAttr(h);
  h->SetTitle(title);
  h->SetName(name);
}

void doHistStuff(TH2* h, const char* name, const char* title)
{
  setHistAttr(h);
  h->SetTitle(title);
  h->SetName(name);
}

TH1D* ratioSuppressLowStats(TH1 *hgenie, TH1* hnuwro, const char* name, const char* title, 
			    const double inputPOT, const double truePOT)
{
  assert(hgenie->GetNbinsX() == hnuwro->GetNbinsX());
  TH1D *hout = new TH1D(name, title, hgenie->GetNbinsX(), 0., hgenie->GetXaxis()->GetBinLowEdge(hgenie->GetNbinsX()+1));
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

int colourFromCat(const int cat) 
{
  int col=0;
  if (cat<1 || cat>nCats) {
    std::cout<<"Invalid category!!!"<<std::endl;
  }
  if (cat==1) col = 632; // kRed
  else if (cat==2) col = 600; // kBlue
  else if (cat==3) col = 1;   // kBlack
  else if (cat==4) col = 417; // kGreen+1
  else if (cat==5) col = 618; // kMagenta+2
  else if (cat==6) col = 802; // kOrange+2
  //else col=1; // kBlack

  return col;
}

THStack* makeCatStack(const std::vector<TH1*> hvec, const char* name, const char* title)
{
  THStack *hs = new THStack(name, title);
  assert(hvec.size()==nCats);
  for (unsigned int i=0; i<hvec.size(); i++) {
    // Create clones of the histograms to preserve original formatting
    TH1D *h = (TH1D*)hvec.at(i)->Clone("h");
    h->SetLineColor(colourFromCat(i+1));
    h->SetMarkerColor(colourFromCat(i+1));
    hs->Add(h);
  }
  return hs;
}

TH2 *makeHistFromPred(PredictionInterp* pred, osc::IOscCalculatorAdjustable* calc, SystShifts shift,
		      const double pot, const char* name, const char* title)
{
  TH2 *h2 = pred->PredictSyst(calc, shift).FakeData(pot).ToTH2(pot);
  doHistStuff(h2, name, title);
  h2->Scale(1./h2->Integral();
  return h2;
}

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// True MC files POT
const double fhcPOT = 1.9342e20;
const double rhcPOT = 3.9302e20;

const double mmu = 0.10566; // GeV/c^2

// Vars
// Reconstructed particle multiplicities
const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var kPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi        = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kRecoP         = SIMPLEVAR(dune.gastpc_pro_mult);
const Var kRecoOtherHad  = SIMPLEVAR(dune.gastpc_other_had_mult); 
const Var kRecoHad       = SIMPLEVAR(dune.gastpc_pi_pl_mult)+SIMPLEVAR(dune.gastpc_pi_min_mult)+SIMPLEVAR(dune.gastpc_pro_mult)+SIMPLEVAR(dune.gastpc_other_had_mult);
// True particle multiplicities
const Var kPipl      = SIMPLEVAR(dune.nipip);
const Var kPimin     = SIMPLEVAR(dune.nipim);
const Var kChargedPi = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip);
const Var kPi        = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);
const Var kP         = SIMPLEVAR(dune.nP);
const Var kOther     = SIMPLEVAR(dune.niother) + SIMPLEVAR(dune.nNucleus);
const Var kHad       = SIMPLEVAR(dune.nipip)+SIMPLEVAR(dune.nipim)+SIMPLEVAR(dune.nipi0)+SIMPLEVAR(dune.nikp)+SIMPLEVAR(dune.nikm)+SIMPLEVAR(dune.nik0)+SIMPLEVAR(dune.nP)+SIMPLEVAR(dune.nN)/*+SIMPLEVAR(dune.niother)+SIMPLEVAR(dune.nNucleus)*/;
const Var kChargedHad = SIMPLEVAR(dune.nipip)+SIMPLEVAR(dune.nipim)+SIMPLEVAR(dune.nikp)+SIMPLEVAR(dune.nikm)+SIMPLEVAR(dune.nP);
const Var kRecoNumu  =  SIMPLEVAR(dune.reco_numu);
// True and reco lepton energy
const Var kLepE = SIMPLEVAR(dune.LepE);
const Var kLepEReco = SIMPLEVAR(dune.Elep_reco);
// True Q2
const Var kQ2 = SIMPLEVAR(dune.Q2);
// Reco Q2
const Var kRecoQ2({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		  [](const caf::StandardRecord* sr) {
		    double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		    double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		    return q2;
		  });
// Experimental Q2 with true variables
const Var kExpQ2({"dune.Ev", "dune.LepE", "dune.LepNuAngle"}, 
		 [](const caf::StandardRecord* sr) {
		   double pmu = sqrt(sr->dune.LepE*sr->dune.LepE - mmu*mmu);
		   double q2 = 2 * sr->dune.Ev * (sr->dune.LepE - pmu*TMath::Cos(sr->dune.LepNuAngle)) - mmu*mmu;
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

// CM's categories
const Var kRecoCategory({"dune.gastpc_pi_pl_mult", "dune.gastpc_pi_min_mult", "dune.gastpc_pi_0_mult"},
			[](const caf::StandardRecord* sr) {
			  int cat = -1;
			  int nChargedPi = sr->dune.gastpc_pi_pl_mult + sr->dune.gastpc_pi_min_mult;
			  int nPi = sr->dune.gastpc_pi_pl_mult + sr->dune.gastpc_pi_min_mult + sr->dune.gastpc_pi_0_mult;
			  if (nPi==0) cat=1;
			  else if (sr->dune.gastpc_pi_min_mult==1 && nPi==1) cat=2;
			  else if (sr->dune.gastpc_pi_pl_mult==1 && nPi==1) cat=3;
			  else if (sr->dune.gastpc_pi_0_mult==1 && nPi==1) cat=4;
			  else if (nPi==2) cat=5;
			  else if (nPi>2) cat=6;
			  return cat;
			});

const Var kTrueCategory({"dune.nipip", "dune.nipim", "dune.nipi0"},
			[](const caf::StandardRecord* sr) {
			  int cat = -1;
			  int nChargedPi = sr->dune.nipip + sr->dune.nipim;
			  int nPi = sr->dune.nipip + sr->dune.nipim + sr->dune.nipi0;
			  if (nPi==0) cat=1;
			  else if (sr->dune.gastpc_pi_min_mult==1 && nPi==1) cat=2;
			  else if (sr->dune.gastpc_pi_pl_mult==1 && nPi==1) cat=3;
			  else if (sr->dune.nipi0==1 && nPi==1) cat=4;
			  else if (nPi==2) cat=5;
			  else if (nPi>2) cat=6;
			  return cat;
			});

// HC
const Var kFHC = SIMPLEVAR(dune.isFHC);

const Binning binsCategory = Binning::Simple(nCats, 0.5, 0.5+nCats);

const Binning binsW  = Binning::Simple(20, 0., 3.);
const Binning binsQ2 = Binning::Simple(20, 0., 5.);

const HistAxis axRecoW ("W_{reco} / GeV", binsW, kRecoW);
const HistAxis axRecoQ2("Q^{2}_{reco} / (GeV)^{2}", binsQ2, kRecoQ2);
const HistAxis axExpQ2("Q^{2}_{exp} / (GeV)^{2}", binsQ2, kExpQ2);
// 2D axis for comparing reco to true Q2
const HistAxis axQ2Comp("Q^{2}_{true} / (GeV)^{2}", binsQ2, kQ2,
			"Q^{2}_{reco} / (GeV)^{2}", binsQ2, kRecoQ2);
const HistAxis axQ2CompTrueExp("Q^{2}_{true} / (GeV)^{2}", binsQ2, kQ2,
			       "Q^{2}_{exp} / (GeV)^{2}", binsQ2, kExpQ2);
const HistAxis axQ2CompRecoExp("Q^{2}_{reco} / (GeV)^{2}", binsQ2, kRecoQ2,
			       "Q^{2}_{exp} / (GeV)^{2}", binsQ2, kExpQ2);
const HistAxis axLepEComp("E_{#mu} / GeV", Binning::Simple(30, 0., 5.), kLepE,
			  "E_{#mu, reco} / GeV", Binning::Simple(30, 0., 5.), kLepEReco);
const HistAxis axLepEReco("E_{#mu, reco} / GeV", Binning::Simple(30, 0., 5.), kLepEReco);

// Cuts used in this analysis for various true and reco categories
// True selections
const Cut kPassCat1({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipip + sr->dune.nipim + sr->dune.nipi0==0
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat2({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipim==1 && 
			      sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==1
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat3({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipip==1 && 
			      sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==1
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat4({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipi0==1 && 
			      sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==1
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat5({},
		    [](const caf::StandardRecord*sr)
		    {
		      return (sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==2
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat6({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0>2 
			      && abs(sr->dune.LepPDG)==13);
		    });
// Reco selections
const Cut kPassRecoCat1({},
			[](const caf::StandardRecord* sr)
			{
			  return (sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==0);
			});
const Cut kPassRecoCat2({},
			[](const caf::StandardRecord* sr)
			{
			  return (sr->dune.gastpc_pi_min_mult==1 && 
				  sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==1);
			});
const Cut kPassRecoCat3({},
			[](const caf::StandardRecord* sr)
			{
			  return (sr->dune.gastpc_pi_pl_mult==1 && 
				  sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==1);
			});
const Cut kPassRecoCat4({},
			[](const caf::StandardRecord* sr)
			{
			  return (sr->dune.gastpc_pi_0_mult==1 && 
				  sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==1);
			});
const Cut kPassRecoCat5({},
			[](const caf::StandardRecord*sr)
			{
			  return (sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==2);
			});
const Cut kPassRecoCat6({},
			[](const caf::StandardRecord* sr)
			{
			  return (sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult>2);
				  
			});

void furtherGasPlots(const char *outfile, 
		     const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v8/") 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();
 
  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = {};
  std::vector<const ISyst*> fakedatavec = GetXSecSysts({"NuWroReweightFakeData"});
  SystShifts fakedata(fakedatavec.at(fakedatavec.size()-1), 1);
  systlist.insert(systlist.end(), fakedatavec.end()-1, fakedatavec.end());
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }

  assert(systlist.size()==1);

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;

  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", garDir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", garDir), kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  //-------------------------FHC----------------------------------------//
  // True vs. reco Q2
  NoOscPredictionGenerator genFhcQ2Comp(axQ2Comp, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcQ2CompCat1(axQ2Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcQ2CompCat2(axQ2Comp, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcQ2CompCat3(axQ2Comp, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcQ2CompCat4(axQ2Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcQ2CompCat5(axQ2Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcQ2CompCat6(axQ2Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcQ2Comp(systlist, 0, genFhcQ2Comp, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat1(systlist, 0, genFhcQ2CompCat1, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat2(systlist, 0, genFhcQ2CompCat2, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat3(systlist, 0, genFhcQ2CompCat3, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat4(systlist, 0, genFhcQ2CompCat4, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat5(systlist, 0, genFhcQ2CompCat5, loadersGArFHC);
  PredictionInterp predFhcQ2CompCat6(systlist, 0, genFhcQ2CompCat6, loadersGArFHC);
  NoOscPredictionGenerator genFhcLepEComp(axLepEComp, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcLepECompCat1(axLepEComp,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcLepECompCat2(axLepEComp,kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcLepECompCat3(axLepEComp,kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcLepECompCat4(axLepEComp,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcLepECompCat5(axLepEComp,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcLepECompCat6(axLepEComp,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcLepEComp(systlist, 0, genFhcLepEComp, loadersGArFHC);
  PredictionInterp predFhcLepECompCat1(systlist, 0, genFhcLepECompCat1, loadersGArFHC);
  PredictionInterp predFhcLepECompCat2(systlist, 0, genFhcLepECompCat2, loadersGArFHC);
  PredictionInterp predFhcLepECompCat3(systlist, 0, genFhcLepECompCat3, loadersGArFHC);
  PredictionInterp predFhcLepECompCat4(systlist, 0, genFhcLepECompCat4, loadersGArFHC);
  PredictionInterp predFhcLepECompCat5(systlist, 0, genFhcLepECompCat5, loadersGArFHC);
  PredictionInterp predFhcLepECompCat6(systlist, 0, genFhcLepECompCat6, loadersGArFHC);
  NoOscPredictionGenerator genFhcLepEReco(axLepEReco, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcLepERecoCat1(axLepEReco,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcLepERecoCat2(axLepEReco,kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcLepERecoCat3(axLepEReco,kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcLepERecoCat4(axLepEReco,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcLepERecoCat5(axLepEReco,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcLepERecoCat6(axLepEReco,kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcLepEReco(systlist, 0, genFhcLepEReco, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat1(systlist, 0, genFhcLepERecoCat1, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat2(systlist, 0, genFhcLepERecoCat2, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat3(systlist, 0, genFhcLepERecoCat3, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat4(systlist, 0, genFhcLepERecoCat4, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat5(systlist, 0, genFhcLepERecoCat5, loadersGArFHC);
  PredictionInterp predFhcLepERecoCat6(systlist, 0, genFhcLepERecoCat6, loadersGArFHC);
  // Comparisons of various types of Q2
  NoOscPredictionGenerator genFhcQ2CompTrueExp(axQ2CompTrueExp, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat1(axQ2CompTrueExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat2(axQ2CompTrueExp, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat3(axQ2CompTrueExp, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat4(axQ2CompTrueExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat5(axQ2CompTrueExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcQ2CompTrueExpCat6(axQ2CompTrueExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcQ2CompTrueExp(systlist, 0, genFhcQ2CompTrueExp, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat1(systlist, 0, genFhcQ2CompTrueExpCat1, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat2(systlist, 0, genFhcQ2CompTrueExpCat2, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat3(systlist, 0, genFhcQ2CompTrueExpCat3, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat4(systlist, 0, genFhcQ2CompTrueExpCat4, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat5(systlist, 0, genFhcQ2CompTrueExpCat5, loadersGArFHC);
  PredictionInterp predFhcQ2CompTrueExpCat6(systlist, 0, genFhcQ2CompTrueExpCat6, loadersGArFHC);
  NoOscPredictionGenerator genFhcQ2CompRecoExp(axQ2CompRecoExp, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat1(axQ2CompRecoExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat2(axQ2CompRecoExp, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat3(axQ2CompRecoExp, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat4(axQ2CompRecoExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat5(axQ2CompRecoExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcQ2CompRecoExpCat6(axQ2CompRecoExp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcQ2CompRecoExp(systlist, 0, genFhcQ2CompRecoExp, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat1(systlist, 0, genFhcQ2CompRecoExpCat1, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat2(systlist, 0, genFhcQ2CompRecoExpCat2, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat3(systlist, 0, genFhcQ2CompRecoExpCat3, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat4(systlist, 0, genFhcQ2CompRecoExpCat4, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat5(systlist, 0, genFhcQ2CompRecoExpCat5, loadersGArFHC);
  PredictionInterp predFhcQ2CompRecoExpCat6(systlist, 0, genFhcQ2CompRecoExpCat6, loadersGArFHC);

  //-------------------------RHC----------------------------------------//
  NoOscPredictionGenerator genRhcQ2Comp(axQ2Comp, kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genRhcQ2CompCat1(axQ2Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genRhcQ2CompCat2(axQ2Comp, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genRhcQ2CompCat3(axQ2Comp, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genRhcQ2CompCat4(axQ2Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genRhcQ2CompCat5(axQ2Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genRhcQ2CompCat6(axQ2Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predRhcQ2Comp(systlist, 0, genRhcQ2Comp, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat1(systlist, 0, genRhcQ2CompCat1, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat2(systlist, 0, genRhcQ2CompCat2, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat3(systlist, 0, genRhcQ2CompCat3, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat4(systlist, 0, genRhcQ2CompCat4, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat5(systlist, 0, genRhcQ2CompCat5, loadersGArRHC);
  PredictionInterp predRhcQ2CompCat6(systlist, 0, genRhcQ2CompCat6, loadersGArRHC);
  NoOscPredictionGenerator genRhcLepEComp(axLepEComp,kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genRhcLepECompCat1(axLepEComp,kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genRhcLepECompCat2(axLepEComp,kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genRhcLepECompCat3(axLepEComp,kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genRhcLepECompCat4(axLepEComp,kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genRhcLepECompCat5(axLepEComp,kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genRhcLepECompCat6(axLepEComp,kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predRhcLepEComp(systlist, 0, genRhcLepEComp, loadersGArRHC);
  PredictionInterp predRhcLepECompCat1(systlist, 0, genRhcLepECompCat1, loadersGArRHC);
  PredictionInterp predRhcLepECompCat2(systlist, 0, genRhcLepECompCat2, loadersGArRHC);
  PredictionInterp predRhcLepECompCat3(systlist, 0, genRhcLepECompCat3, loadersGArRHC);
  PredictionInterp predRhcLepECompCat4(systlist, 0, genRhcLepECompCat4, loadersGArRHC);
  PredictionInterp predRhcLepECompCat5(systlist, 0, genRhcLepECompCat5, loadersGArRHC);
  PredictionInterp predRhcLepECompCat6(systlist, 0, genRhcLepECompCat6, loadersGArRHC);

  loadersGArFHC.Go();
  // loadersGArRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  // Make the 2D hists
  // Q2 comp
  TH2 *h2FhcQ2Comp = makeHistFromPred(&predFhcQ2Comp, this_calc, kNoShift, pot_nd, "h2FhcQ2Comp", "CC inc. (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events");
  TH2 *h2FhcQ2CompCat1 = makeHistFromPred(&predFhcQ2CompCat1, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat1", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(1).c_str()));
  TH2 *h2FhcQ2CompCat2 = makeHistFromPred(&predFhcQ2CompCat2, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat2", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(2).c_str()));
  TH2 *h2FhcQ2CompCat3 = makeHistFromPred(&predFhcQ2CompCat3, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat3", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(3).c_str()));
  TH2 *h2FhcQ2CompCat4 = makeHistFromPred(&predFhcQ2CompCat4, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat4", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(4).c_str()));
  TH2 *h2FhcQ2CompCat5 = makeHistFromPred(&predFhcQ2CompCat5, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat5", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(5).c_str()));
  TH2 *h2FhcQ2CompCat6 = makeHistFromPred(&predFhcQ2CompCat6, this_calc, kNoShift, pot_nd, "h2FhcQ2CompCat6", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{reco} / (GeV)^{2}; Events", catName(6).c_str()));
  h2FhcQ2Comp->Write();
  h2FhcQ2CompCat1->Write();
  h2FhcQ2CompCat2->Write();
  h2FhcQ2CompCat3->Write();
  h2FhcQ2CompCat4->Write();
  h2FhcQ2CompCat5->Write();
  h2FhcQ2CompCat6->Write();

  TH2 *h2FhcQ2CompRecoExp = makeHistFromPred(&predFhcQ2CompRecoExp, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExp", "CC inc. (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events");
  TH2 *h2FhcQ2CompRecoExpCat1 = makeHistFromPred(&predFhcQ2CompRecoExpCat1, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat1", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(1).c_str()));
  TH2 *h2FhcQ2CompRecoExpCat2 = makeHistFromPred(&predFhcQ2CompRecoExpCat2, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat2", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(2).c_str()));
  TH2 *h2FhcQ2CompRecoExpCat3 = makeHistFromPred(&predFhcQ2CompRecoExpCat3, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat3", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(3).c_str()));
  TH2 *h2FhcQ2CompRecoExpCat4 = makeHistFromPred(&predFhcQ2CompRecoExpCat4, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat4", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(4).c_str()));
  TH2 *h2FhcQ2CompRecoExpCat5 = makeHistFromPred(&predFhcQ2CompRecoExpCat5, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat5", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(5).c_str()));
  TH2 *h2FhcQ2CompRecoExpCat6 = makeHistFromPred(&predFhcQ2CompRecoExpCat6, this_calc, kNoShift, pot_nd, "h2FhcQ2CompRecoExpCat6", Form("%s FHC (HPgTPC); Q^{2}_{reco} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(6).c_str()));
  h2FhcQ2CompRecoExp->Write();
  h2FhcQ2CompRecoExpCat1->Write();
  h2FhcQ2CompRecoExpCat2->Write();
  h2FhcQ2CompRecoExpCat3->Write();
  h2FhcQ2CompRecoExpCat4->Write();
  h2FhcQ2CompRecoExpCat5->Write();
  h2FhcQ2CompRecoExpCat6->Write();

  TH2 *h2FhcQ2CompTrueExp = makeHistFromPred(&predFhcQ2CompTrueExp, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExp", "CC inc. (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events");
  TH2 *h2FhcQ2CompTrueExpCat1 = makeHistFromPred(&predFhcQ2CompTrueExpCat1, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat1", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(1).c_str()));
  TH2 *h2FhcQ2CompTrueExpCat2 = makeHistFromPred(&predFhcQ2CompTrueExpCat2, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat2", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(2).c_str()));
  TH2 *h2FhcQ2CompTrueExpCat3 = makeHistFromPred(&predFhcQ2CompTrueExpCat3, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat3", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(3).c_str()));
  TH2 *h2FhcQ2CompTrueExpCat4 = makeHistFromPred(&predFhcQ2CompTrueExpCat4, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat4", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(4).c_str()));
  TH2 *h2FhcQ2CompTrueExpCat5 = makeHistFromPred(&predFhcQ2CompTrueExpCat5, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat5", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(5).c_str()));
  TH2 *h2FhcQ2CompTrueExpCat6 = makeHistFromPred(&predFhcQ2CompTrueExpCat6, this_calc, kNoShift, pot_nd, "h2FhcQ2CompTrueExpCat6", Form("%s FHC (HPgTPC); Q^{2}_{true} / (GeV)^{2}; Q^{2}_{exp} / (GeV)^{2}; Fraction of events", catName(6).c_str()));
  h2FhcQ2CompTrueExp->Write();
  h2FhcQ2CompTrueExpCat1->Write();
  h2FhcQ2CompTrueExpCat2->Write();
  h2FhcQ2CompTrueExpCat3->Write();
  h2FhcQ2CompTrueExpCat4->Write();
  h2FhcQ2CompTrueExpCat5->Write();
  h2FhcQ2CompTrueExpCat6->Write();

  // Lepton energy
  TH2 *h2FhcLepEComp = makeHistFromPred(&predFhcLepEComp, this_calc, kNoShift, pot_nd, "h2FhcLepEComp", "CC inc. (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events");
  TH2 *h2FhcLepECompCat1 = makeHistFromPred(&predFhcLepECompCat1, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat1", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(1).c_str()));
  TH2 *h2FhcLepECompCat2 = makeHistFromPred(&predFhcLepECompCat2, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat2", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(2).c_str()));
  TH2 *h2FhcLepECompCat3 = makeHistFromPred(&predFhcLepECompCat3, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat3", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(3).c_str()));
  TH2 *h2FhcLepECompCat4 = makeHistFromPred(&predFhcLepECompCat4, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat4", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(4).c_str()));
  TH2 *h2FhcLepECompCat5 = makeHistFromPred(&predFhcLepECompCat5, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat5", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(5).c_str()));
  TH2 *h2FhcLepECompCat6 = makeHistFromPred(&predFhcLepECompCat6, this_calc, kNoShift, pot_nd, "h2FhcLepECompCat6", Form("%s FHC (HPgTPC); E_{#mu, true} / GeV; E_{#mu, reco} / GeV; Events", catName(6).c_str()));
  h2FhcLepEComp->Write();
  h2FhcLepECompCat1->Write();
  h2FhcLepECompCat2->Write();
  h2FhcLepECompCat3->Write();
  h2FhcLepECompCat4->Write();
  h2FhcLepECompCat5->Write();
  h2FhcLepECompCat6->Write();
 
  TH1 *hFhcLepEReco     = predFhcLepEReco.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat1 = predFhcLepERecoCat1.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat2 = predFhcLepERecoCat2.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat3 = predFhcLepERecoCat3.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat4 = predFhcLepERecoCat4.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat5 = predFhcLepERecoCat5.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat6 = predFhcLepERecoCat6.PredictSyst(this_calc, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hFhcLepEReco,"hFhcLepEReco","CC inc.");
  doHistStuff(hFhcLepERecoCat1,"hFhcLepERecoCat1",Form("%s", catName(1).c_str()));
  doHistStuff(hFhcLepERecoCat2,"hFhcLepERecoCat2",Form("%s", catName(2).c_str()));
  doHistStuff(hFhcLepERecoCat3,"hFhcLepERecoCat3",Form("%s", catName(3).c_str()));
  doHistStuff(hFhcLepERecoCat4,"hFhcLepERecoCat4",Form("%s", catName(4).c_str()));
  doHistStuff(hFhcLepERecoCat5,"hFhcLepERecoCat5",Form("%s", catName(5).c_str()));
  doHistStuff(hFhcLepERecoCat6,"hFhcLepERecoCat6",Form("%s", catName(6).c_str()));
  hFhcLepEReco->Write();
  hFhcLepERecoCat1->Write();
  hFhcLepERecoCat2->Write();
  hFhcLepERecoCat3->Write();
  hFhcLepERecoCat4->Write();
  hFhcLepERecoCat5->Write();
  hFhcLepERecoCat6->Write();
  TH1 *hFhcLepEReco_n     = predFhcLepEReco.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat1_n = predFhcLepERecoCat1.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat2_n = predFhcLepERecoCat2.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat3_n = predFhcLepERecoCat3.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat4_n = predFhcLepERecoCat4.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat5_n = predFhcLepERecoCat5.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcLepERecoCat6_n = predFhcLepERecoCat6.PredictSyst(this_calc, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hFhcLepEReco_n,"hFhcLepEReco_n","CC inc.");
  doHistStuff(hFhcLepERecoCat1_n,"hFhcLepERecoCat1_n",Form("%s", catName(1).c_str()));
  doHistStuff(hFhcLepERecoCat2_n,"hFhcLepERecoCat2_n",Form("%s", catName(2).c_str()));
  doHistStuff(hFhcLepERecoCat3_n,"hFhcLepERecoCat3_n",Form("%s", catName(3).c_str()));
  doHistStuff(hFhcLepERecoCat4_n,"hFhcLepERecoCat4_n",Form("%s", catName(4).c_str()));
  doHistStuff(hFhcLepERecoCat5_n,"hFhcLepERecoCat5_n",Form("%s", catName(5).c_str()));
  doHistStuff(hFhcLepERecoCat6_n,"hFhcLepERecoCat6_n",Form("%s", catName(6).c_str()));
  hFhcLepEReco_n->Write();
  hFhcLepERecoCat1_n->Write();
  hFhcLepERecoCat2_n->Write();
  hFhcLepERecoCat3_n->Write();
  hFhcLepERecoCat4_n->Write();
  hFhcLepERecoCat5_n->Write();
  hFhcLepERecoCat6_n->Write();
  // Ratio
  TH1D *hrFhcLepEReco = ratioSuppressLowStats(hFhcLepEReco, hFhcLepEReco_n, "hrFhcLepEReco", "CC inc.; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat1 = ratioSuppressLowStats(hFhcLepERecoCat1, hFhcLepERecoCat1_n, "hrFhcLepERecoCat1", "0#pi; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat2 = ratioSuppressLowStats(hFhcLepERecoCat2, hFhcLepERecoCat2_n, "hrFhcLepERecoCat2", "1#pi^{-}; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat3 = ratioSuppressLowStats(hFhcLepERecoCat3, hFhcLepERecoCat3_n, "hrFhcLepERecoCat3", "1#pi^{+}; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat4 = ratioSuppressLowStats(hFhcLepERecoCat4, hFhcLepERecoCat4_n, "hrFhcLepERecoCat4", "1#pi^{0}; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat5 = ratioSuppressLowStats(hFhcLepERecoCat5, hFhcLepERecoCat5_n, "hrFhcLepERecoCat5", "2#pi; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcLepERecoCat6 = ratioSuppressLowStats(hFhcLepERecoCat6, hFhcLepERecoCat6_n, "hrFhcLepERecoCat6", ">2#pi; E_{#mu, reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  hrFhcLepEReco->Write();
  hrFhcLepERecoCat1->Write();
  hrFhcLepERecoCat2->Write();
  hrFhcLepERecoCat3->Write();
  hrFhcLepERecoCat4->Write();
  hrFhcLepERecoCat5->Write();
  hrFhcLepERecoCat6->Write();
  std::vector<TH1*> fhcLepERecoVec;
  fhcLepERecoVec.push_back(hrFhcLepERecoCat1);
  fhcLepERecoVec.push_back(hrFhcLepERecoCat2);
  fhcLepERecoVec.push_back(hrFhcLepERecoCat3);
  fhcLepERecoVec.push_back(hrFhcLepERecoCat4);
  fhcLepERecoVec.push_back(hrFhcLepERecoCat5);
  fhcLepERecoVec.push_back(hrFhcLepERecoCat6);
  THStack *hsFhcLepEReco = makeCatStack(fhcLepERecoVec, "hsFhcLepEReco", "NuWro/GENIE for various reconstructed final states; E_{#mu, reco} / GeV; NuWro/GENIE");
  hsFhcLepEReco->Write();

  fout->Close();
} // furtherGasPlots
