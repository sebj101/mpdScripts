// spreadTest.C
// Test out method for calculating error bars based upon confusion matrices
// Use CM's categories as a first past

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

void makeConfusionMatrix(TH2* h2)
{
  for (int binX=1; binX<h2->GetNbinsX()+1; binX++) {
    double colInt       = h2->Integral(binX, binX, 0, h2->GetNbinsY()+1);
    for (int binY=1; binY<h2->GetNbinsY()+1; binY++) {
      h2->SetBinContent(binX, binY, h2->GetBinContent(binX, binY) / colInt);
    }
  }
}

void doHistStuff (TH1* h, const char* name, const char* title)
{
  setHistAttr(h);
  h->SetTitle(title);
  h->Write(name);
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
  hout->Write();
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

TH1D* ratioWithSpread(const char* name, const char* title, 
		      const std::vector<TH1*> hNReco, const std::vector<TH1*> hGReco, 
		      const std::vector<TH1*> hNTrue, const std::vector<TH1*> hGTrue,
		      const TH2* matrix, const int cat, const double inputPOT, const double truePOT)
{
  if (cat<1 || cat>nCats) std::cout<<"Category is invalid!!!"<<std::endl;

  std::vector<TH1D*> rTvec;
  for (int i=1; i<=nCats; i++) {
    TH1D *tT = new TH1D(Form("tT%s%d",name,i), "", hNReco.at(i-1)->GetNbinsX(), 0., hNReco.at(i-1)->GetXaxis()->GetBinLowEdge(hNReco.at(i-1)->GetNbinsX()+1));
    tT->Divide(hNTrue.at(i-1), hGTrue.at(i-1));
    rTvec.push_back(tT);
  }
  TH1D *hr = new TH1D(name, title, hNReco.at(cat-1)->GetNbinsX(), 0., hNReco.at(cat-1)->GetXaxis()->GetBinLowEdge(hNReco.at(cat-1)->GetNbinsX()+1));
  
  for (int i=1; i<=hNReco.at(cat-1)->GetNbinsX(); i++) {
    double num = hNReco.at(cat-1)->GetBinContent(i);
    double denom = hGReco.at(cat-1)->GetBinContent(i);
    if (denom * (truePOT/inputPOT) < 100.) continue;
    // if (denom < 1000.) continue;
    // Get error from the spread in the smearing
    double true_cv = rTvec.at(cat-1)->GetBinContent(i);
    double spread = 0.;
    double norm = matrix->ProjectionX("", cat, cat)->Integral(1, nCats);
    for (int j=1; j<=nCats; j++) {
      double diff = rTvec.at(j-1)->GetBinContent(i) - true_cv;
      // Fraction of reco events from this true category
      double smear_wgt = matrix->GetBinContent(j, cat) / norm;
      spread += smear_wgt * pow(diff, 2);
    }
    spread = sqrt(spread);
    hr->SetBinContent(i, num/denom);
    hr->SetBinError(i, spread);
  }
  return hr;
}

using namespace ana;

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
			  else if (sr->dune.nipim==1 && nPi==1) cat=2;
			  else if (sr->dune.nipip==1 && nPi==1) cat=3;
			  else if (sr->dune.nipi0==1 && nPi==1) cat=4;
			  else if (nPi==2) cat=5;
			  else if (nPi>2) cat=6;
			  return cat;
			});

// HC
const Var kFHC = SIMPLEVAR(dune.isFHC);

const Binning binsCategory = Binning::Simple(nCats, 0.5, 0.5+nCats);
const HistAxis axCategory("True category", binsCategory, kTrueCategory, 
			  "Reco category", binsCategory, kRecoCategory);

const int nBinsKinematics = 20;
const Binning binsW  = Binning::Simple(20, 0., 5.);
const Binning binsQ2 = Binning::Simple(20, 0., 5.);
const HistAxis axRecoW ("W_{reco} / GeV", binsW, kRecoW);
const HistAxis axRecoQ2("Q^{2}_{reco} / (GeV)^{2}", binsQ2, kRecoQ2);

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

void spreadTest(const char *outfile, 
		const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v9/") 
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

  // Reco vs true categories
  NoOscPredictionGenerator genFhcCat(axCategory, kRecoNumu==1 && kIsTrueGasFV);
  PredictionInterp predFhcCat(systlist, 0, genFhcCat, loadersGArFHC);
  NoOscPredictionGenerator genRhcCat(axCategory, kRecoNumu==1 && kIsTrueGasFV);
  PredictionInterp predRhcCat(systlist, 0, genRhcCat, loadersGArRHC);
  //------------------------FHC---------------------------
  // True categories
  NoOscPredictionGenerator genFhcQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcQ2Cat1(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcQ2Cat2(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcQ2Cat3(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcQ2Cat4(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcQ2Cat5(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcQ2Cat6(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcQ2(systlist, 0, genFhcQ2, loadersGArFHC);
  PredictionInterp predFhcQ2Cat1(systlist, 0, genFhcQ2Cat1, loadersGArFHC);
  PredictionInterp predFhcQ2Cat2(systlist, 0, genFhcQ2Cat2, loadersGArFHC);
  PredictionInterp predFhcQ2Cat3(systlist, 0, genFhcQ2Cat3, loadersGArFHC);
  PredictionInterp predFhcQ2Cat4(systlist, 0, genFhcQ2Cat4, loadersGArFHC);
  PredictionInterp predFhcQ2Cat5(systlist, 0, genFhcQ2Cat5, loadersGArFHC);
  PredictionInterp predFhcQ2Cat6(systlist, 0, genFhcQ2Cat6, loadersGArFHC);
  NoOscPredictionGenerator genFhcW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genFhcWCat1(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genFhcWCat2(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genFhcWCat3(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genFhcWCat4(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genFhcWCat5(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genFhcWCat6(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predFhcW(systlist, 0, genFhcW, loadersGArFHC);
  PredictionInterp predFhcWCat1(systlist, 0, genFhcWCat1, loadersGArFHC);
  PredictionInterp predFhcWCat2(systlist, 0, genFhcWCat2, loadersGArFHC);
  PredictionInterp predFhcWCat3(systlist, 0, genFhcWCat3, loadersGArFHC);
  PredictionInterp predFhcWCat4(systlist, 0, genFhcWCat4, loadersGArFHC);
  PredictionInterp predFhcWCat5(systlist, 0, genFhcWCat5, loadersGArFHC);
  PredictionInterp predFhcWCat6(systlist, 0, genFhcWCat6, loadersGArFHC);
  // Reco categories
  NoOscPredictionGenerator genFhcQ2RecoCat1(axRecoQ2,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genFhcQ2RecoCat2(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genFhcQ2RecoCat3(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genFhcQ2RecoCat4(axRecoQ2,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genFhcQ2RecoCat5(axRecoQ2,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genFhcQ2RecoCat6(axRecoQ2,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  PredictionInterp predFhcQ2RecoCat1(systlist, 0, genFhcQ2RecoCat1, loadersGArFHC);
  PredictionInterp predFhcQ2RecoCat2(systlist, 0, genFhcQ2RecoCat2, loadersGArFHC);
  PredictionInterp predFhcQ2RecoCat3(systlist, 0, genFhcQ2RecoCat3, loadersGArFHC);
  PredictionInterp predFhcQ2RecoCat4(systlist, 0, genFhcQ2RecoCat4, loadersGArFHC);
  PredictionInterp predFhcQ2RecoCat5(systlist, 0, genFhcQ2RecoCat5, loadersGArFHC);
  PredictionInterp predFhcQ2RecoCat6(systlist, 0, genFhcQ2RecoCat6, loadersGArFHC);
  NoOscPredictionGenerator genFhcWRecoCat1(axRecoW,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genFhcWRecoCat2(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genFhcWRecoCat3(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genFhcWRecoCat4(axRecoW,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genFhcWRecoCat5(axRecoW,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genFhcWRecoCat6(axRecoW,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  PredictionInterp predFhcWRecoCat1(systlist, 0, genFhcWRecoCat1, loadersGArFHC);
  PredictionInterp predFhcWRecoCat2(systlist, 0, genFhcWRecoCat2, loadersGArFHC);
  PredictionInterp predFhcWRecoCat3(systlist, 0, genFhcWRecoCat3, loadersGArFHC);
  PredictionInterp predFhcWRecoCat4(systlist, 0, genFhcWRecoCat4, loadersGArFHC);
  PredictionInterp predFhcWRecoCat5(systlist, 0, genFhcWRecoCat5, loadersGArFHC);
  PredictionInterp predFhcWRecoCat6(systlist, 0, genFhcWRecoCat6, loadersGArFHC);
  //------------------------------------RHC--------------------------
  // True categories
  NoOscPredictionGenerator genRhcQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genRhcQ2Cat1(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genRhcQ2Cat2(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genRhcQ2Cat3(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genRhcQ2Cat4(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genRhcQ2Cat5(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genRhcQ2Cat6(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predRhcQ2(systlist, 0, genRhcQ2, loadersGArRHC);
  PredictionInterp predRhcQ2Cat1(systlist, 0, genRhcQ2Cat1, loadersGArRHC);
  PredictionInterp predRhcQ2Cat2(systlist, 0, genRhcQ2Cat2, loadersGArRHC);
  PredictionInterp predRhcQ2Cat3(systlist, 0, genRhcQ2Cat3, loadersGArRHC);
  PredictionInterp predRhcQ2Cat4(systlist, 0, genRhcQ2Cat4, loadersGArRHC);
  PredictionInterp predRhcQ2Cat5(systlist, 0, genRhcQ2Cat5, loadersGArRHC);
  PredictionInterp predRhcQ2Cat6(systlist, 0, genRhcQ2Cat6, loadersGArRHC);
  NoOscPredictionGenerator genRhcW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genRhcWCat1(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genRhcWCat2(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genRhcWCat3(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genRhcWCat4(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genRhcWCat5(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genRhcWCat6(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kPassCat6);
  PredictionInterp predRhcW(systlist, 0, genRhcW, loadersGArRHC);
  PredictionInterp predRhcWCat1(systlist, 0, genRhcWCat1, loadersGArRHC);
  PredictionInterp predRhcWCat2(systlist, 0, genRhcWCat2, loadersGArRHC);
  PredictionInterp predRhcWCat3(systlist, 0, genRhcWCat3, loadersGArRHC);
  PredictionInterp predRhcWCat4(systlist, 0, genRhcWCat4, loadersGArRHC);
  PredictionInterp predRhcWCat5(systlist, 0, genRhcWCat5, loadersGArRHC);
  PredictionInterp predRhcWCat6(systlist, 0, genRhcWCat6, loadersGArRHC);
  // Reco categories
  NoOscPredictionGenerator genRhcQ2RecoCat1(axRecoQ2,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genRhcQ2RecoCat2(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genRhcQ2RecoCat3(axRecoQ2, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genRhcQ2RecoCat4(axRecoQ2,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genRhcQ2RecoCat5(axRecoQ2,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genRhcQ2RecoCat6(axRecoQ2,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  PredictionInterp predRhcQ2RecoCat1(systlist, 0, genRhcQ2RecoCat1, loadersGArRHC);
  PredictionInterp predRhcQ2RecoCat2(systlist, 0, genRhcQ2RecoCat2, loadersGArRHC);
  PredictionInterp predRhcQ2RecoCat3(systlist, 0, genRhcQ2RecoCat3, loadersGArRHC);
  PredictionInterp predRhcQ2RecoCat4(systlist, 0, genRhcQ2RecoCat4, loadersGArRHC);
  PredictionInterp predRhcQ2RecoCat5(systlist, 0, genRhcQ2RecoCat5, loadersGArRHC);
  PredictionInterp predRhcQ2RecoCat6(systlist, 0, genRhcQ2RecoCat6, loadersGArRHC);
  NoOscPredictionGenerator genRhcWRecoCat1(axRecoW,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genRhcWRecoCat2(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genRhcWRecoCat3(axRecoW, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genRhcWRecoCat4(axRecoW,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genRhcWRecoCat5(axRecoW,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genRhcWRecoCat6(axRecoW,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  PredictionInterp predRhcWRecoCat1(systlist, 0, genRhcWRecoCat1, loadersGArRHC);
  PredictionInterp predRhcWRecoCat2(systlist, 0, genRhcWRecoCat2, loadersGArRHC);
  PredictionInterp predRhcWRecoCat3(systlist, 0, genRhcWRecoCat3, loadersGArRHC);
  PredictionInterp predRhcWRecoCat4(systlist, 0, genRhcWRecoCat4, loadersGArRHC);
  PredictionInterp predRhcWRecoCat5(systlist, 0, genRhcWRecoCat5, loadersGArRHC);
  PredictionInterp predRhcWRecoCat6(systlist, 0, genRhcWRecoCat6, loadersGArRHC);

  loadersGArFHC.Go();
  loadersGArRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  //---------------------------FHC------------------------------------
  // Reco vs. true categories
  TH2 *h2Cat       = predFhcCat.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Cat_nuwro = predFhcCat.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH2(pot_nd);
  h2Cat->Scale(1. / h2Cat->Integral());
  h2Cat_nuwro->Scale(1. / h2Cat_nuwro->Integral());
  setHistAttr(h2Cat);
  setHistAttr(h2Cat_nuwro);
  h2Cat->SetTitle("True vs. reco final state in HPgTPC; True category; Reco category; Fraction of events");
  h2Cat_nuwro->SetTitle("True vs. reco state in HPgTPC (NuWro shifts); True category; Reco category; Fraction of events");
  for (int binX=1; binX<h2Cat->GetNbinsX()+1; binX++) {
    h2Cat->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat_nuwro->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Cat_nuwro->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
  }
  h2Cat->Write("h2Cat");
  h2Cat_nuwro->Write("h2Cat_nuwro");
  makeConfusionMatrix(h2Cat);
  makeConfusionMatrix(h2Cat_nuwro);
  h2Cat->SetTitle("Final state confusion matrix in HPgTPC; True category; Reco category; Probability");
  h2Cat_nuwro->SetTitle("Final state confusion matrix in HPgTPC (NuWro shifts); True category; Reco category; Probability");
  h2Cat->Write("h2Cat_confus");
  h2Cat_nuwro->Write("h2Cat_nuwro_confus");

  // Now get the reco Q2 and reco W for the true and reco selections
  TH1 *hFhcQ2     = predFhcQ2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat1 = predFhcQ2Cat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat2 = predFhcQ2Cat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat3 = predFhcQ2Cat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat4 = predFhcQ2Cat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat5 = predFhcQ2Cat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat6 = predFhcQ2Cat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2_n     = predFhcQ2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat1_n = predFhcQ2Cat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat2_n = predFhcQ2Cat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat3_n = predFhcQ2Cat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat4_n = predFhcQ2Cat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat5_n = predFhcQ2Cat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2Cat6_n = predFhcQ2Cat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcW     = predFhcW.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat1 = predFhcWCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat2 = predFhcWCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat3 = predFhcWCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat4 = predFhcWCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat5 = predFhcWCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat6 = predFhcWCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcW_n     = predFhcW.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat1_n = predFhcWCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat2_n = predFhcWCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat3_n = predFhcWCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat4_n = predFhcWCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat5_n = predFhcWCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWCat6_n = predFhcWCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hFhcQ2,"hFhcQ2","CC inc.");
  doHistStuff(hFhcQ2Cat1,"hFhcQ2Cat1","0#pi");
  doHistStuff(hFhcQ2Cat2,"hFhcQ2Cat2","1#pi^{-}");
  doHistStuff(hFhcQ2Cat3,"hFhcQ2Cat3","1#pi^{+}");
  doHistStuff(hFhcQ2Cat4,"hFhcQ2Cat4","1#pi^{0}");
  doHistStuff(hFhcQ2Cat5,"hFhcQ2Cat5","2#pi");
  doHistStuff(hFhcQ2Cat6,"hFhcQ2Cat6",">2#pi");
  doHistStuff(hFhcQ2_n,"hFhcQ2_n","CC inc.");
  doHistStuff(hFhcQ2Cat1_n,"hFhcQ2Cat1_n","0#pi");
  doHistStuff(hFhcQ2Cat2_n,"hFhcQ2Cat2_n","1#pi^{-}");
  doHistStuff(hFhcQ2Cat3_n,"hFhcQ2Cat3_n","1#pi^{+}");
  doHistStuff(hFhcQ2Cat4_n,"hFhcQ2Cat4_n","1#pi^{0}");
  doHistStuff(hFhcQ2Cat5_n,"hFhcQ2Cat5_n","2#pi");
  doHistStuff(hFhcQ2Cat6_n,"hFhcQ2Cat6_n",">2#pi");
  doHistStuff(hFhcW,"hFhcW","CC inc.");
  doHistStuff(hFhcWCat1,"hFhcWCat1","0#pi");
  doHistStuff(hFhcWCat2,"hFhcWCat2","1#pi^{-}");
  doHistStuff(hFhcWCat3,"hFhcWCat3","1#pi^{+}");
  doHistStuff(hFhcWCat4,"hFhcWCat4","1#pi^{0}");
  doHistStuff(hFhcWCat5,"hFhcWCat5","2#pi");
  doHistStuff(hFhcWCat6,"hFhcWCat6",">2#pi");
  doHistStuff(hFhcW_n,"hFhcW_n","CC inc.");
  doHistStuff(hFhcWCat1_n,"hFhcWCat1_n","0#pi");
  doHistStuff(hFhcWCat2_n,"hFhcWCat2_n","1#pi^{-}");
  doHistStuff(hFhcWCat3_n,"hFhcWCat3_n","1#pi^{+}");
  doHistStuff(hFhcWCat4_n,"hFhcWCat4_n","1#pi^{0}");
  doHistStuff(hFhcWCat5_n,"hFhcWCat5_n","2#pi");
  doHistStuff(hFhcWCat6_n,"hFhcWCat6_n",">2#pi");
  TH1D *hrFhcQ2 = ratioSuppressLowStats(hFhcQ2, hFhcQ2_n, "hrFhcQ2", "CC inc.; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat1 = ratioSuppressLowStats(hFhcQ2Cat1, hFhcQ2Cat1_n, "hrFhcQ2Cat1", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat2 = ratioSuppressLowStats(hFhcQ2Cat2, hFhcQ2Cat2_n, "hrFhcQ2Cat2", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat3 = ratioSuppressLowStats(hFhcQ2Cat3, hFhcQ2Cat3_n, "hrFhcQ2Cat3", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat4 = ratioSuppressLowStats(hFhcQ2Cat4, hFhcQ2Cat4_n, "hrFhcQ2Cat4", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat5 = ratioSuppressLowStats(hFhcQ2Cat5, hFhcQ2Cat5_n, "hrFhcQ2Cat5", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2Cat6 = ratioSuppressLowStats(hFhcQ2Cat6, hFhcQ2Cat6_n, "hrFhcQ2Cat6", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  std::vector<TH1*> fhcQ2Vec;
  fhcQ2Vec.push_back(hrFhcQ2Cat1);
  fhcQ2Vec.push_back(hrFhcQ2Cat2);
  fhcQ2Vec.push_back(hrFhcQ2Cat3);
  fhcQ2Vec.push_back(hrFhcQ2Cat4);
  fhcQ2Vec.push_back(hrFhcQ2Cat5);
  fhcQ2Vec.push_back(hrFhcQ2Cat6);
  THStack *hsFhcQ2 = makeCatStack(fhcQ2Vec, "hsFhcQ2", "NuWro/GENIE for various true final states; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE");
  hsFhcQ2->Write();

  TH1D *hrFhcW = ratioSuppressLowStats(hFhcW, hFhcW_n, "hrFhcW", "CC inc.; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat1 = ratioSuppressLowStats(hFhcWCat1, hFhcWCat1_n, "hrFhcWCat1", "0#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat2 = ratioSuppressLowStats(hFhcWCat2, hFhcWCat2_n, "hrFhcWCat2", "1#pi^{-}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat3 = ratioSuppressLowStats(hFhcWCat3, hFhcWCat3_n, "hrFhcWCat3", "1#pi^{+}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat4 = ratioSuppressLowStats(hFhcWCat4, hFhcWCat4_n, "hrFhcWCat4", "1#pi^{0}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat5 = ratioSuppressLowStats(hFhcWCat5, hFhcWCat5_n, "hrFhcWCat5", "2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWCat6 = ratioSuppressLowStats(hFhcWCat6, hFhcWCat6_n, "hrFhcWCat6", ">2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  std::vector<TH1*> fhcWVec;
  fhcWVec.push_back(hrFhcWCat1);
  fhcWVec.push_back(hrFhcWCat2);
  fhcWVec.push_back(hrFhcWCat3);
  fhcWVec.push_back(hrFhcWCat4);
  fhcWVec.push_back(hrFhcWCat5);
  fhcWVec.push_back(hrFhcWCat6);
  THStack *hsFhcW = makeCatStack(fhcWVec, "hsFhcW", "NuWro/GENIE for various true final states; W_{reco} / GeV; NuWro/GENIE");
  hsFhcW->Write();
 
  // Clear these and put the non ratio plots in
  fhcWVec.clear();
  fhcWVec.push_back(hFhcWCat1);
  fhcWVec.push_back(hFhcWCat2);
  fhcWVec.push_back(hFhcWCat3);
  fhcWVec.push_back(hFhcWCat4);
  fhcWVec.push_back(hFhcWCat5);
  fhcWVec.push_back(hFhcWCat6);
  fhcQ2Vec.clear();
  fhcQ2Vec.push_back(hFhcQ2Cat1);
  fhcQ2Vec.push_back(hFhcQ2Cat2);
  fhcQ2Vec.push_back(hFhcQ2Cat3);
  fhcQ2Vec.push_back(hFhcQ2Cat4);
  fhcQ2Vec.push_back(hFhcQ2Cat5);
  fhcQ2Vec.push_back(hFhcQ2Cat6);

  std::vector<TH1*> fhcWVec_n;
  fhcWVec_n.push_back(hFhcWCat1_n);
  fhcWVec_n.push_back(hFhcWCat2_n);
  fhcWVec_n.push_back(hFhcWCat3_n);
  fhcWVec_n.push_back(hFhcWCat4_n);
  fhcWVec_n.push_back(hFhcWCat5_n);
  fhcWVec_n.push_back(hFhcWCat6_n);
  std::vector<TH1*> fhcQ2Vec_n;
  fhcQ2Vec_n.push_back(hFhcQ2Cat1_n);
  fhcQ2Vec_n.push_back(hFhcQ2Cat2_n);
  fhcQ2Vec_n.push_back(hFhcQ2Cat3_n);
  fhcQ2Vec_n.push_back(hFhcQ2Cat4_n);
  fhcQ2Vec_n.push_back(hFhcQ2Cat5_n);
  fhcQ2Vec_n.push_back(hFhcQ2Cat6_n);

  // Reco selections
  TH1 *hFhcQ2RecoCat1 = predFhcQ2RecoCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat2 = predFhcQ2RecoCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat3 = predFhcQ2RecoCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat4 = predFhcQ2RecoCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat5 = predFhcQ2RecoCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat6 = predFhcQ2RecoCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat1_n = predFhcQ2RecoCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat2_n = predFhcQ2RecoCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat3_n = predFhcQ2RecoCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat4_n = predFhcQ2RecoCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat5_n = predFhcQ2RecoCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcQ2RecoCat6_n = predFhcQ2RecoCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat1 = predFhcWRecoCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat2 = predFhcWRecoCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat3 = predFhcWRecoCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat4 = predFhcWRecoCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat5 = predFhcWRecoCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat6 = predFhcWRecoCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat1_n = predFhcWRecoCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat2_n = predFhcWRecoCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat3_n = predFhcWRecoCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat4_n = predFhcWRecoCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat5_n = predFhcWRecoCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hFhcWRecoCat6_n = predFhcWRecoCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hFhcQ2RecoCat1,"hFhcQ2RecoCat1","0#pi");
  doHistStuff(hFhcQ2RecoCat2,"hFhcQ2RecoCat2","1#pi^{-}");
  doHistStuff(hFhcQ2RecoCat3,"hFhcQ2RecoCat3","1#pi^{+}");
  doHistStuff(hFhcQ2RecoCat4,"hFhcQ2RecoCat4","1#pi^{0}");
  doHistStuff(hFhcQ2RecoCat5,"hFhcQ2RecoCat5","2#pi");
  doHistStuff(hFhcQ2RecoCat6,"hFhcQ2RecoCat6",">2#pi");
  doHistStuff(hFhcQ2RecoCat1_n,"hFhcQ2RecoCat1_n","0#pi");
  doHistStuff(hFhcQ2RecoCat2_n,"hFhcQ2RecoCat2_n","1#pi^{-}");
  doHistStuff(hFhcQ2RecoCat3_n,"hFhcQ2RecoCat3_n","1#pi^{+}");
  doHistStuff(hFhcQ2RecoCat4_n,"hFhcQ2RecoCat4_n","1#pi^{0}");
  doHistStuff(hFhcQ2RecoCat5_n,"hFhcQ2RecoCat5_n","2#pi");
  doHistStuff(hFhcQ2RecoCat6_n,"hFhcQ2RecoCat6_n",">2#pi");
  doHistStuff(hFhcWRecoCat1,"hFhcWRecoCat1","0#pi");
  doHistStuff(hFhcWRecoCat2,"hFhcWRecoCat2","1#pi^{-}");
  doHistStuff(hFhcWRecoCat3,"hFhcWRecoCat3","1#pi^{+}");
  doHistStuff(hFhcWRecoCat4,"hFhcWRecoCat4","1#pi^{0}");
  doHistStuff(hFhcWRecoCat5,"hFhcWRecoCat5","2#pi");
  doHistStuff(hFhcWRecoCat6,"hFhcWRecoCat6",">2#pi");
  doHistStuff(hFhcWRecoCat1_n,"hFhcWRecoCat1_n","0#pi");
  doHistStuff(hFhcWRecoCat2_n,"hFhcWRecoCat2_n","1#pi^{-}");
  doHistStuff(hFhcWRecoCat3_n,"hFhcWRecoCat3_n","1#pi^{+}");
  doHistStuff(hFhcWRecoCat4_n,"hFhcWRecoCat4_n","1#pi^{0}");
  doHistStuff(hFhcWRecoCat5_n,"hFhcWRecoCat5_n","2#pi");
  doHistStuff(hFhcWRecoCat5_n,"hFhcWRecoCat5_n",">2#pi");

  TH1D *hrFhcQ2RecoCat1 = ratioSuppressLowStats(hFhcQ2RecoCat1, hFhcQ2RecoCat1_n, "hrFhcQ2RecoCat1", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2RecoCat2 = ratioSuppressLowStats(hFhcQ2RecoCat2, hFhcQ2RecoCat2_n, "hrFhcQ2RecoCat2", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2RecoCat3 = ratioSuppressLowStats(hFhcQ2RecoCat3, hFhcQ2RecoCat3_n, "hrFhcQ2RecoCat3", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2RecoCat4 = ratioSuppressLowStats(hFhcQ2RecoCat4, hFhcQ2RecoCat4_n, "hrFhcQ2RecoCat4", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2RecoCat5 = ratioSuppressLowStats(hFhcQ2RecoCat5, hFhcQ2RecoCat5_n, "hrFhcQ2RecoCat5", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcQ2RecoCat6 = ratioSuppressLowStats(hFhcQ2RecoCat6, hFhcQ2RecoCat6_n, "hrFhcQ2RecoCat6", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, fhcPOT);
  std::vector<TH1*> fhcQ2RecoVec;
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat1);
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat2);
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat3);
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat4);
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat5);
  fhcQ2RecoVec.push_back(hrFhcQ2RecoCat6);
  THStack *hsFhcQ2Reco = makeCatStack(fhcQ2RecoVec, "hsFhcQ2Reco", "NuWro/GENIE for various reconstructed final states; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE");
  hsFhcQ2Reco->Write();

  TH1D *hrFhcWRecoCat1 = ratioSuppressLowStats(hFhcWRecoCat1, hFhcWRecoCat1_n, "hrFhcWRecoCat1", "0#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWRecoCat2 = ratioSuppressLowStats(hFhcWRecoCat2, hFhcWRecoCat2_n, "hrFhcWRecoCat2", "1#pi^{-}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWRecoCat3 = ratioSuppressLowStats(hFhcWRecoCat3, hFhcWRecoCat3_n, "hrFhcWRecoCat3", "1#pi^{+}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWRecoCat4 = ratioSuppressLowStats(hFhcWRecoCat4, hFhcWRecoCat4_n, "hrFhcWRecoCat4", "1#pi^{0}; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWRecoCat5 = ratioSuppressLowStats(hFhcWRecoCat5, hFhcWRecoCat5_n, "hrFhcWRecoCat5", "2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  TH1D *hrFhcWRecoCat6 = ratioSuppressLowStats(hFhcWRecoCat6, hFhcWRecoCat6_n, "hrFhcWRecoCat6", ">2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, fhcPOT);
  std::vector<TH1*> fhcWRecoVec;
  fhcWRecoVec.push_back(hrFhcWRecoCat1);
  fhcWRecoVec.push_back(hrFhcWRecoCat2);
  fhcWRecoVec.push_back(hrFhcWRecoCat3);
  fhcWRecoVec.push_back(hrFhcWRecoCat4);
  fhcWRecoVec.push_back(hrFhcWRecoCat5);
  fhcWRecoVec.push_back(hrFhcWRecoCat6);
  THStack *hsFhcWReco = makeCatStack(fhcWRecoVec, "hsFhcWReco", "NuWro/GENIE for various reconstructed final states; W_{reco} / GeV; NuWro/GENIE");
  hsFhcWReco->Write();

  fhcWRecoVec.clear();
  fhcWRecoVec.push_back(hFhcWRecoCat1);
  fhcWRecoVec.push_back(hFhcWRecoCat2);
  fhcWRecoVec.push_back(hFhcWRecoCat3);
  fhcWRecoVec.push_back(hFhcWRecoCat4);
  fhcWRecoVec.push_back(hFhcWRecoCat5);
  fhcWRecoVec.push_back(hFhcWRecoCat6);
  fhcQ2RecoVec.clear();
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat1);
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat2);
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat3);
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat4);
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat5);
  fhcQ2RecoVec.push_back(hFhcQ2RecoCat6);

  std::vector<TH1*> fhcWRecoVec_n;
  fhcWRecoVec_n.push_back(hFhcWRecoCat1_n);
  fhcWRecoVec_n.push_back(hFhcWRecoCat2_n);
  fhcWRecoVec_n.push_back(hFhcWRecoCat3_n);
  fhcWRecoVec_n.push_back(hFhcWRecoCat4_n);
  fhcWRecoVec_n.push_back(hFhcWRecoCat5_n);
  fhcWRecoVec_n.push_back(hFhcWRecoCat6_n);
  std::vector<TH1*> fhcQ2RecoVec_n;
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat1_n);
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat2_n);
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat3_n);
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat4_n);
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat5_n);
  fhcQ2RecoVec_n.push_back(hFhcQ2RecoCat6_n);

  // Make the same plots with the better errors

  TH1D* hrFhcQ2RecoCat1_spread = ratioWithSpread("hrFhcQ2RecoCat1_spread", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 1, pot_nd, fhcPOT);
  TH1D* hrFhcQ2RecoCat2_spread = ratioWithSpread("hrFhcQ2RecoCat2_spread", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 2, pot_nd, fhcPOT);
  TH1D* hrFhcQ2RecoCat3_spread = ratioWithSpread("hrFhcQ2RecoCat3_spread", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 3, pot_nd, fhcPOT);
  TH1D* hrFhcQ2RecoCat4_spread = ratioWithSpread("hrFhcQ2RecoCat4_spread", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 4, pot_nd, fhcPOT);
  TH1D* hrFhcQ2RecoCat5_spread = ratioWithSpread("hrFhcQ2RecoCat5_spread", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 5, pot_nd, fhcPOT);
  TH1D* hrFhcQ2RecoCat6_spread = ratioWithSpread("hrFhcQ2RecoCat6_spread", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", fhcQ2RecoVec_n, fhcQ2RecoVec, fhcQ2Vec_n, fhcQ2Vec, h2Cat, 6, pot_nd, fhcPOT);
  doHistStuff(hrFhcQ2RecoCat1_spread, "hrFhcQ2RecoCat1_spread", "0#pi");
  doHistStuff(hrFhcQ2RecoCat2_spread, "hrFhcQ2RecoCat2_spread", "1#pi^{-}");
  doHistStuff(hrFhcQ2RecoCat3_spread, "hrFhcQ2RecoCat3_spread", "1#pi^{+}");
  doHistStuff(hrFhcQ2RecoCat4_spread, "hrFhcQ2RecoCat4_spread", "1#pi^{0}");
  doHistStuff(hrFhcQ2RecoCat5_spread, "hrFhcQ2RecoCat5_spread", "2#pi");
  doHistStuff(hrFhcQ2RecoCat6_spread, "hrFhcQ2RecoCat6_spread", ">2#pi");
  std::vector<TH1*> fhcQ2RecoVec_spread;
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat1_spread);
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat2_spread);
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat3_spread);
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat4_spread);
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat5_spread);
  fhcQ2RecoVec_spread.push_back(hrFhcQ2RecoCat6_spread);
  THStack *hsFhcQ2Reco_spread = makeCatStack(fhcQ2RecoVec_spread, "hsFhcQ2Reco_spread", "NuWro/GENIE for various reconstructed final states; Q^{2}_{reco} / GeV; NuWro/GENIE");
  hsFhcQ2Reco_spread->Write();

  TH1D* hrFhcWRecoCat1_spread = ratioWithSpread("hrFhcWRecoCat1_spread", "0#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 1, pot_nd, fhcPOT);
  TH1D* hrFhcWRecoCat2_spread = ratioWithSpread("hrFhcWRecoCat2_spread", "1#pi^{-}; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 2, pot_nd, fhcPOT);
  TH1D* hrFhcWRecoCat3_spread = ratioWithSpread("hrFhcWRecoCat3_spread", "1#pi^{+}; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 3, pot_nd, fhcPOT);
  TH1D* hrFhcWRecoCat4_spread = ratioWithSpread("hrFhcWRecoCat4_spread", "1#pi^{0}; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 4, pot_nd, fhcPOT);
  TH1D* hrFhcWRecoCat5_spread = ratioWithSpread("hrFhcWRecoCat5_spread", "2#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 5, pot_nd, fhcPOT);
  TH1D* hrFhcWRecoCat6_spread = ratioWithSpread("hrFhcWRecoCat6_spread", ">2#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", fhcWRecoVec_n, fhcWRecoVec, fhcWVec_n, fhcWVec, h2Cat, 6, pot_nd, fhcPOT);
  doHistStuff(hrFhcWRecoCat1_spread, "hrFhcWRecoCat1_spread", "0#pi");
  doHistStuff(hrFhcWRecoCat2_spread, "hrFhcWRecoCat2_spread", "1#pi^{-}");
  doHistStuff(hrFhcWRecoCat3_spread, "hrFhcWRecoCat3_spread", "1#pi^{+}");
  doHistStuff(hrFhcWRecoCat4_spread, "hrFhcWRecoCat4_spread", "1#pi^{0}");
  doHistStuff(hrFhcWRecoCat5_spread, "hrFhcWRecoCat5_spread", "2#pi");
  doHistStuff(hrFhcWRecoCat6_spread, "hrFhcWRecoCat6_spread", ">2#pi");
  std::vector<TH1*> fhcWRecoVec_spread;
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat1_spread);
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat2_spread);
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat3_spread);
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat4_spread);
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat5_spread);
  fhcWRecoVec_spread.push_back(hrFhcWRecoCat6_spread);
  THStack *hsFhcWReco_spread = makeCatStack(fhcWRecoVec_spread, "hsFhcWReco_spread", "NuWro/GENIE for various reconstructed final states; W_{reco} / GeV; NuWro/GENIE");
  hsFhcWReco_spread->Write();

  //------------------------------------RHC------------------------------
  TH2 *h2Catrhc       = predRhcCat.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2Catrhc_nuwro = predRhcCat.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH2(pot_nd);
  h2Catrhc->Scale(1. / h2Catrhc->Integral());
  h2Catrhc_nuwro->Scale(1. / h2Catrhc_nuwro->Integral());
  setHistAttr(h2Catrhc);
  setHistAttr(h2Catrhc_nuwro);
  h2Catrhc->SetTitle("True vs. reco final state in HPgTPC; True category; Reco category; Fraction of events");
  h2Catrhc_nuwro->SetTitle("True vs. reco state in HPgTPC (NuWro shifts); True category; Reco category; Fraction of events");
  for (int binX=1; binX<h2Catrhc->GetNbinsX()+1; binX++) {
    h2Catrhc->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Catrhc_nuwro->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Catrhc->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2Catrhc_nuwro->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
  }
  h2Catrhc->Write("h2Catrhc");
  h2Catrhc_nuwro->Write("h2Catrhc_nuwro");
  makeConfusionMatrix(h2Catrhc);
  makeConfusionMatrix(h2Catrhc_nuwro);
  h2Catrhc->SetTitle("Final state confusion matrix in HPgTPC; True category; Reco category; Probability");
  h2Catrhc_nuwro->SetTitle("Final state confusion matrix in HPgTPC (NuWro shifts); True category; Reco category; Probability");
  h2Catrhc->Write("h2Catrhc_confus");
  h2Catrhc_nuwro->Write("h2Catrhc_nuwro_confus");

  // Now get the reco Q2 and reco W for the true and reco selections
  TH1 *hRhcQ2     = predRhcQ2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat1 = predRhcQ2Cat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat2 = predRhcQ2Cat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat3 = predRhcQ2Cat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat4 = predRhcQ2Cat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat5 = predRhcQ2Cat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat6 = predRhcQ2Cat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2_n     = predRhcQ2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat1_n = predRhcQ2Cat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat2_n = predRhcQ2Cat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat3_n = predRhcQ2Cat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat4_n = predRhcQ2Cat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat5_n = predRhcQ2Cat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2Cat6_n = predRhcQ2Cat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcW     = predRhcW.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat1 = predRhcWCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat2 = predRhcWCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat3 = predRhcWCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat4 = predRhcWCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat5 = predRhcWCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat6 = predRhcWCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcW_n     = predRhcW.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat1_n = predRhcWCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat2_n = predRhcWCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat3_n = predRhcWCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat4_n = predRhcWCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat5_n = predRhcWCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWCat6_n = predRhcWCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hRhcQ2,"hRhcQ2","CC inc.");
  doHistStuff(hRhcQ2Cat1,"hRhcQ2Cat1","0#pi");
  doHistStuff(hRhcQ2Cat2,"hRhcQ2Cat2","1#pi^{-}");
  doHistStuff(hRhcQ2Cat3,"hRhcQ2Cat3","1#pi^{+}");
  doHistStuff(hRhcQ2Cat4,"hRhcQ2Cat4","1#pi^{0}");
  doHistStuff(hRhcQ2Cat5,"hRhcQ2Cat5","2#pi");
  doHistStuff(hRhcQ2Cat6,"hRhcQ2Cat6",">2#pi");
  doHistStuff(hRhcQ2_n,"hRhcQ2_n","CC inc.");
  doHistStuff(hRhcQ2Cat1_n,"hRhcQ2Cat1_n","0#pi");
  doHistStuff(hRhcQ2Cat2_n,"hRhcQ2Cat2_n","1#pi^{-}");
  doHistStuff(hRhcQ2Cat3_n,"hRhcQ2Cat3_n","1#pi^{+}");
  doHistStuff(hRhcQ2Cat4_n,"hRhcQ2Cat4_n","1#pi^{0}");
  doHistStuff(hRhcQ2Cat5_n,"hRhcQ2Cat5_n","2#pi");
  doHistStuff(hRhcQ2Cat6_n,"hRhcQ2Cat6_n",">2#pi");
  doHistStuff(hRhcW,"hRhcW","CC inc.");
  doHistStuff(hRhcWCat1,"hRhcWCat1","0#pi");
  doHistStuff(hRhcWCat2,"hRhcWCat2","1#pi^{-}");
  doHistStuff(hRhcWCat3,"hRhcWCat3","1#pi^{+}");
  doHistStuff(hRhcWCat4,"hRhcWCat4","1#pi^{0}");
  doHistStuff(hRhcWCat5,"hRhcWCat5","2#pi");
  doHistStuff(hRhcWCat6,"hRhcWCat6",">2#pi");
  doHistStuff(hRhcW_n,"hRhcW_n","CC inc.");
  doHistStuff(hRhcWCat1_n,"hRhcWCat1_n","0#pi");
  doHistStuff(hRhcWCat2_n,"hRhcWCat2_n","1#pi^{-}");
  doHistStuff(hRhcWCat3_n,"hRhcWCat3_n","1#pi^{+}");
  doHistStuff(hRhcWCat4_n,"hRhcWCat4_n","1#pi^{0}");
  doHistStuff(hRhcWCat5_n,"hRhcWCat5_n","2#pi");
  doHistStuff(hRhcWCat6_n,"hRhcWCat6_n",">2#pi");
  TH1D *hrRhcQ2 = ratioSuppressLowStats(hRhcQ2, hRhcQ2_n, "hrRhcQ2", "CC inc.; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat1 = ratioSuppressLowStats(hRhcQ2Cat1, hRhcQ2Cat1_n, "hrRhcQ2Cat1", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat2 = ratioSuppressLowStats(hRhcQ2Cat2, hRhcQ2Cat2_n, "hrRhcQ2Cat2", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat3 = ratioSuppressLowStats(hRhcQ2Cat3, hRhcQ2Cat3_n, "hrRhcQ2Cat3", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat4 = ratioSuppressLowStats(hRhcQ2Cat4, hRhcQ2Cat4_n, "hrRhcQ2Cat4", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat5 = ratioSuppressLowStats(hRhcQ2Cat5, hRhcQ2Cat5_n, "hrRhcQ2Cat5", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2Cat6 = ratioSuppressLowStats(hRhcQ2Cat6, hRhcQ2Cat6_n, "hrRhcQ2Cat6", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  std::vector<TH1*> rhcQ2Vec;
  rhcQ2Vec.push_back(hrRhcQ2Cat1);
  rhcQ2Vec.push_back(hrRhcQ2Cat2);
  rhcQ2Vec.push_back(hrRhcQ2Cat3);
  rhcQ2Vec.push_back(hrRhcQ2Cat4);
  rhcQ2Vec.push_back(hrRhcQ2Cat5);
  rhcQ2Vec.push_back(hrRhcQ2Cat6);
  THStack *hsRhcQ2 = makeCatStack(rhcQ2Vec, "hsRhcQ2", "NuWro/GENIE for various true final states; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE");
  hsRhcQ2->Write();

  TH1D *hrRhcW = ratioSuppressLowStats(hRhcW, hRhcW_n, "hrRhcW", "CC inc.; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat1 = ratioSuppressLowStats(hRhcWCat1, hRhcWCat1_n, "hrRhcWCat1", "0#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat2 = ratioSuppressLowStats(hRhcWCat2, hRhcWCat2_n, "hrRhcWCat2", "1#pi^{-}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat3 = ratioSuppressLowStats(hRhcWCat3, hRhcWCat3_n, "hrRhcWCat3", "1#pi^{+}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat4 = ratioSuppressLowStats(hRhcWCat4, hRhcWCat4_n, "hrRhcWCat4", "1#pi^{0}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat5 = ratioSuppressLowStats(hRhcWCat5, hRhcWCat5_n, "hrRhcWCat5", "2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWCat6 = ratioSuppressLowStats(hRhcWCat6, hRhcWCat6_n, "hrRhcWCat6", ">2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  std::vector<TH1*> rhcWVec;
  rhcWVec.push_back(hrRhcWCat1);
  rhcWVec.push_back(hrRhcWCat2);
  rhcWVec.push_back(hrRhcWCat3);
  rhcWVec.push_back(hrRhcWCat4);
  rhcWVec.push_back(hrRhcWCat5);
  rhcWVec.push_back(hrRhcWCat6);
  THStack *hsRhcW = makeCatStack(rhcWVec, "hsRhcW", "NuWro/GENIE for various true final states; W_{reco} / GeV; NuWro/GENIE");
  hsRhcW->Write();
 
  // Clear these and put the non ratio plots in
  rhcWVec.clear();
  rhcWVec.push_back(hRhcWCat1);
  rhcWVec.push_back(hRhcWCat2);
  rhcWVec.push_back(hRhcWCat3);
  rhcWVec.push_back(hRhcWCat4);
  rhcWVec.push_back(hRhcWCat5);
  rhcWVec.push_back(hRhcWCat6);
  rhcQ2Vec.clear();
  rhcQ2Vec.push_back(hRhcQ2Cat1);
  rhcQ2Vec.push_back(hRhcQ2Cat2);
  rhcQ2Vec.push_back(hRhcQ2Cat3);
  rhcQ2Vec.push_back(hRhcQ2Cat4);
  rhcQ2Vec.push_back(hRhcQ2Cat5);
  rhcQ2Vec.push_back(hRhcQ2Cat6);

  std::vector<TH1*> rhcWVec_n;
  rhcWVec_n.push_back(hRhcWCat1_n);
  rhcWVec_n.push_back(hRhcWCat2_n);
  rhcWVec_n.push_back(hRhcWCat3_n);
  rhcWVec_n.push_back(hRhcWCat4_n);
  rhcWVec_n.push_back(hRhcWCat5_n);
  rhcWVec_n.push_back(hRhcWCat6_n);
  std::vector<TH1*> rhcQ2Vec_n;
  rhcQ2Vec_n.push_back(hRhcQ2Cat1_n);
  rhcQ2Vec_n.push_back(hRhcQ2Cat2_n);
  rhcQ2Vec_n.push_back(hRhcQ2Cat3_n);
  rhcQ2Vec_n.push_back(hRhcQ2Cat4_n);
  rhcQ2Vec_n.push_back(hRhcQ2Cat5_n);
  rhcQ2Vec_n.push_back(hRhcQ2Cat6_n);

  // Reco selections
  TH1 *hRhcQ2RecoCat1 = predRhcQ2RecoCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat2 = predRhcQ2RecoCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat3 = predRhcQ2RecoCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat4 = predRhcQ2RecoCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat5 = predRhcQ2RecoCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat6 = predRhcQ2RecoCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat1_n = predRhcQ2RecoCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat2_n = predRhcQ2RecoCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat3_n = predRhcQ2RecoCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat4_n = predRhcQ2RecoCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat5_n = predRhcQ2RecoCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcQ2RecoCat6_n = predRhcQ2RecoCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat1 = predRhcWRecoCat1.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat2 = predRhcWRecoCat2.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat3 = predRhcWRecoCat3.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat4 = predRhcWRecoCat4.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat5 = predRhcWRecoCat5.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat6 = predRhcWRecoCat6.PredictSyst(0, kNoShift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat1_n = predRhcWRecoCat1.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat2_n = predRhcWRecoCat2.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat3_n = predRhcWRecoCat3.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat4_n = predRhcWRecoCat4.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat5_n = predRhcWRecoCat5.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRhcWRecoCat6_n = predRhcWRecoCat6.PredictSyst(0, fakedata).FakeData(pot_nd).ToTH1(pot_nd);
  doHistStuff(hRhcQ2RecoCat1,"hRhcQ2RecoCat1","0#pi");
  doHistStuff(hRhcQ2RecoCat2,"hRhcQ2RecoCat2","1#pi^{-}");
  doHistStuff(hRhcQ2RecoCat3,"hRhcQ2RecoCat3","1#pi^{+}");
  doHistStuff(hRhcQ2RecoCat4,"hRhcQ2RecoCat4","1#pi^{0}");
  doHistStuff(hRhcQ2RecoCat5,"hRhcQ2RecoCat5","2#pi");
  doHistStuff(hRhcQ2RecoCat6,"hRhcQ2RecoCat6",">2#pi");
  doHistStuff(hRhcQ2RecoCat1_n,"hRhcQ2RecoCat1_n","0#pi");
  doHistStuff(hRhcQ2RecoCat2_n,"hRhcQ2RecoCat2_n","1#pi^{-}");
  doHistStuff(hRhcQ2RecoCat3_n,"hRhcQ2RecoCat3_n","1#pi^{+}");
  doHistStuff(hRhcQ2RecoCat4_n,"hRhcQ2RecoCat4_n","1#pi^{0}");
  doHistStuff(hRhcQ2RecoCat5_n,"hRhcQ2RecoCat5_n","2#pi");
  doHistStuff(hRhcQ2RecoCat6_n,"hRhcQ2RecoCat6_n",">2#pi");
  doHistStuff(hRhcWRecoCat1,"hRhcWRecoCat1","0#pi");
  doHistStuff(hRhcWRecoCat2,"hRhcWRecoCat2","1#pi^{-}");
  doHistStuff(hRhcWRecoCat3,"hRhcWRecoCat3","1#pi^{+}");
  doHistStuff(hRhcWRecoCat4,"hRhcWRecoCat4","1#pi^{0}");
  doHistStuff(hRhcWRecoCat5,"hRhcWRecoCat5","2#pi");
  doHistStuff(hRhcWRecoCat6,"hRhcWRecoCat6",">2#pi");
  doHistStuff(hRhcWRecoCat1_n,"hRhcWRecoCat1_n","0#pi");
  doHistStuff(hRhcWRecoCat2_n,"hRhcWRecoCat2_n","1#pi^{-}");
  doHistStuff(hRhcWRecoCat3_n,"hRhcWRecoCat3_n","1#pi^{+}");
  doHistStuff(hRhcWRecoCat4_n,"hRhcWRecoCat4_n","1#pi^{0}");
  doHistStuff(hRhcWRecoCat5_n,"hRhcWRecoCat5_n","2#pi");
  doHistStuff(hRhcWRecoCat5_n,"hRhcWRecoCat5_n",">2#pi");

  TH1D *hrRhcQ2RecoCat1 = ratioSuppressLowStats(hRhcQ2RecoCat1, hRhcQ2RecoCat1_n, "hrRhcQ2RecoCat1", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2RecoCat2 = ratioSuppressLowStats(hRhcQ2RecoCat2, hRhcQ2RecoCat2_n, "hrRhcQ2RecoCat2", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2RecoCat3 = ratioSuppressLowStats(hRhcQ2RecoCat3, hRhcQ2RecoCat3_n, "hrRhcQ2RecoCat3", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2RecoCat4 = ratioSuppressLowStats(hRhcQ2RecoCat4, hRhcQ2RecoCat4_n, "hrRhcQ2RecoCat4", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2RecoCat5 = ratioSuppressLowStats(hRhcQ2RecoCat5, hRhcQ2RecoCat5_n, "hrRhcQ2RecoCat5", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcQ2RecoCat6 = ratioSuppressLowStats(hRhcQ2RecoCat6, hRhcQ2RecoCat6_n, "hrRhcQ2RecoCat6", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", pot_nd, rhcPOT);
  std::vector<TH1*> rhcQ2RecoVec;
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat1);
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat2);
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat3);
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat4);
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat5);
  rhcQ2RecoVec.push_back(hrRhcQ2RecoCat6);
  THStack *hsRhcQ2Reco = makeCatStack(rhcQ2RecoVec, "hsRhcQ2Reco", "NuWro/GENIE for various reconstructed final states; Q^{2}_{reco} / GeV^{2}; NuWro/GENIE");
  hsRhcQ2Reco->Write();

  TH1D *hrRhcWRecoCat1 = ratioSuppressLowStats(hRhcWRecoCat1, hRhcWRecoCat1_n, "hrRhcWRecoCat1", "0#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWRecoCat2 = ratioSuppressLowStats(hRhcWRecoCat2, hRhcWRecoCat2_n, "hrRhcWRecoCat2", "1#pi^{-}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWRecoCat3 = ratioSuppressLowStats(hRhcWRecoCat3, hRhcWRecoCat3_n, "hrRhcWRecoCat3", "1#pi^{+}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWRecoCat4 = ratioSuppressLowStats(hRhcWRecoCat4, hRhcWRecoCat4_n, "hrRhcWRecoCat4", "1#pi^{0}; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWRecoCat5 = ratioSuppressLowStats(hRhcWRecoCat5, hRhcWRecoCat5_n, "hrRhcWRecoCat5", "2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  TH1D *hrRhcWRecoCat6 = ratioSuppressLowStats(hRhcWRecoCat6, hRhcWRecoCat6_n, "hrRhcWRecoCat6", ">2#pi; W_{reco} / GeV; NuWro/GENIE", pot_nd, rhcPOT);
  std::vector<TH1*> rhcWRecoVec;
  rhcWRecoVec.push_back(hrRhcWRecoCat1);
  rhcWRecoVec.push_back(hrRhcWRecoCat2);
  rhcWRecoVec.push_back(hrRhcWRecoCat3);
  rhcWRecoVec.push_back(hrRhcWRecoCat4);
  rhcWRecoVec.push_back(hrRhcWRecoCat5);
  rhcWRecoVec.push_back(hrRhcWRecoCat6);
  THStack *hsRhcWReco = makeCatStack(rhcWRecoVec, "hsRhcWReco", "NuWro/GENIE for various reconstructed final states; W_{reco} / GeV; NuWro/GENIE");
  hsRhcWReco->Write();

  rhcWRecoVec.clear();
  rhcWRecoVec.push_back(hRhcWRecoCat1);
  rhcWRecoVec.push_back(hRhcWRecoCat2);
  rhcWRecoVec.push_back(hRhcWRecoCat3);
  rhcWRecoVec.push_back(hRhcWRecoCat4);
  rhcWRecoVec.push_back(hRhcWRecoCat5);
  rhcWRecoVec.push_back(hRhcWRecoCat6);
  rhcQ2RecoVec.clear();
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat1);
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat2);
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat3);
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat4);
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat5);
  rhcQ2RecoVec.push_back(hRhcQ2RecoCat6);

  std::vector<TH1*> rhcWRecoVec_n;
  rhcWRecoVec_n.push_back(hRhcWRecoCat1_n);
  rhcWRecoVec_n.push_back(hRhcWRecoCat2_n);
  rhcWRecoVec_n.push_back(hRhcWRecoCat3_n);
  rhcWRecoVec_n.push_back(hRhcWRecoCat4_n);
  rhcWRecoVec_n.push_back(hRhcWRecoCat5_n);
  rhcWRecoVec_n.push_back(hRhcWRecoCat6_n);
  std::vector<TH1*> rhcQ2RecoVec_n;
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat1_n);
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat2_n);
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat3_n);
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat4_n);
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat5_n);
  rhcQ2RecoVec_n.push_back(hRhcQ2RecoCat6_n);

  // Make the same plots with the better errors
  TH1D* hrRhcQ2RecoCat1_spread = ratioWithSpread("hrRhcQ2RecoCat1_spread", "0#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 1, pot_nd, rhcPOT);
  TH1D* hrRhcQ2RecoCat2_spread = ratioWithSpread("hrRhcQ2RecoCat2_spread", "1#pi^{-}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 2, pot_nd, rhcPOT);
  TH1D* hrRhcQ2RecoCat3_spread = ratioWithSpread("hrRhcQ2RecoCat3_spread", "1#pi^{+}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 3, pot_nd, rhcPOT);
  TH1D* hrRhcQ2RecoCat4_spread = ratioWithSpread("hrRhcQ2RecoCat4_spread", "1#pi^{0}; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 4, pot_nd, rhcPOT);
  TH1D* hrRhcQ2RecoCat5_spread = ratioWithSpread("hrRhcQ2RecoCat5_spread", "2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 5, pot_nd, rhcPOT);
  TH1D* hrRhcQ2RecoCat6_spread = ratioWithSpread("hrRhcQ2RecoCat6_spread", ">2#pi; Q^{2}_{reco} / (GeV)^{2}; NuWro/GENIE", rhcQ2RecoVec_n, rhcQ2RecoVec, rhcQ2Vec_n, rhcQ2Vec, h2Catrhc, 6, pot_nd, rhcPOT);
  doHistStuff(hrRhcQ2RecoCat1_spread, "hrRhcQ2RecoCat1_spread", "0#pi");
  doHistStuff(hrRhcQ2RecoCat2_spread, "hrRhcQ2RecoCat2_spread", "1#pi^{-}");
  doHistStuff(hrRhcQ2RecoCat3_spread, "hrRhcQ2RecoCat3_spread", "1#pi^{+}");
  doHistStuff(hrRhcQ2RecoCat4_spread, "hrRhcQ2RecoCat4_spread", "1#pi^{0}");
  doHistStuff(hrRhcQ2RecoCat5_spread, "hrRhcQ2RecoCat5_spread", "2#pi");
  doHistStuff(hrRhcQ2RecoCat6_spread, "hrRhcQ2RecoCat6_spread", ">2#pi");
  std::vector<TH1*> rhcQ2RecoVec_spread;
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat1_spread);
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat2_spread);
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat3_spread);
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat4_spread);
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat5_spread);
  rhcQ2RecoVec_spread.push_back(hrRhcQ2RecoCat6_spread);
  THStack *hsRhcQ2Reco_spread = makeCatStack(rhcQ2RecoVec_spread, "hsRhcQ2Reco_spread", "NuWro/GENIE for various reconstructed final states; Q^{2}_{reco} / GeV; NuWro/GENIE");
  hsRhcQ2Reco_spread->Write();

  TH1D* hrRhcWRecoCat1_spread = ratioWithSpread("hrRhcWRecoCat1_spread", "0#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 1, pot_nd, rhcPOT);
  TH1D* hrRhcWRecoCat2_spread = ratioWithSpread("hrRhcWRecoCat2_spread", "1#pi^{-}; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 2, pot_nd, rhcPOT);
  TH1D* hrRhcWRecoCat3_spread = ratioWithSpread("hrRhcWRecoCat3_spread", "1#pi^{+}; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 3, pot_nd, rhcPOT);
  TH1D* hrRhcWRecoCat4_spread = ratioWithSpread("hrRhcWRecoCat4_spread", "1#pi^{0}; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 4, pot_nd, rhcPOT);
  TH1D* hrRhcWRecoCat5_spread = ratioWithSpread("hrRhcWRecoCat5_spread", "2#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 5, pot_nd, rhcPOT);
  TH1D* hrRhcWRecoCat6_spread = ratioWithSpread("hrRhcWRecoCat6_spread", ">2#pi; W_{reco} / (GeV)^{2}; NuWro/GENIE", rhcWRecoVec_n, rhcWRecoVec, rhcWVec_n, rhcWVec, h2Catrhc, 6, pot_nd, rhcPOT);
  doHistStuff(hrRhcWRecoCat1_spread, "hrRhcWRecoCat1_spread", "0#pi");
  doHistStuff(hrRhcWRecoCat2_spread, "hrRhcWRecoCat2_spread", "1#pi^{-}");
  doHistStuff(hrRhcWRecoCat3_spread, "hrRhcWRecoCat3_spread", "1#pi^{+}");
  doHistStuff(hrRhcWRecoCat4_spread, "hrRhcWRecoCat4_spread", "1#pi^{0}");
  doHistStuff(hrRhcWRecoCat5_spread, "hrRhcWRecoCat5_spread", "2#pi");
  doHistStuff(hrRhcWRecoCat6_spread, "hrRhcWRecoCat6_spread", ">2#pi");
  std::vector<TH1*> rhcWRecoVec_spread;
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat1_spread);
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat2_spread);
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat3_spread);
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat4_spread);
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat5_spread);
  rhcWRecoVec_spread.push_back(hrRhcWRecoCat6_spread);
  THStack *hsRhcWReco_spread = makeCatStack(rhcWRecoVec_spread, "hsRhcWReco_spread", "NuWro/GENIE for various reconstructed final states; W_{reco} / GeV; NuWro/GENIE");
  hsRhcWReco_spread->Write();

  TLine *l = new TLine(0., 1., 5., 1.);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->Write("l");

  fout->Close();
  delete fout;
} // spreadTest
