// statsCheck.C
// Check the statistical uncertainty on the MPD multi-pion channels
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
int nCats = 6;

void setHistAttr(TH2 *h2)
{
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->SetOption("colz");
}

void setHistAttr(TH1 *h)
{
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->SetOption("hist");
}

bool IsInFV(const double x, const double y, const double z)
{
  double r = sqrt( (y+75.)*(y+75.) + (z-955.)*(z-955.) );
  return ( r < 200. && abs(x) < 200. );
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
  return col;
}

double highQ2 = 3.;
int nQ2Bins   = 20;

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = 1. * POT120;
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
const Var kRecoNumu = SIMPLEVAR(dune.reco_numu);
// Reco Q2
const Var kRecoQ2({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		  [](const caf::StandardRecord* sr) {
		    double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		    double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		    return q2;
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

const Binning binsCategory = Binning::Simple(nCats, 0.5, 0.5+nCats);
const HistAxis axCatQ2("Q^{2}_{reco}", Binning::Simple(20, 0., 3.), kRecoQ2, 
		       "Reco category", binsCategory, kRecoCategory);
const HistAxis axRecoQ2("Q^{2}_{reco}", Binning::Simple(20, 0., 3.), kRecoQ2);

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

TH1 *statUnc(PredictionInterp* pred, const double pot, const char* name, const char* title)
{
  TH1 *h = pred->Predict(0).FakeData(pot).ToTH1(pot);
  for (int i=1; i<h->GetNbinsX()+1; i++) {
    double unc = sqrt(h->GetBinContent(i));
    h->SetBinContent(i, unc/h->GetBinContent(i));
  }
  setHistAttr(h);
  h->SetTitle(title);
  h->SetName(name);
  return h;
}

TH2 *statUnc2d(PredictionInterp* pred, const double pot, const char* name, const char* title)
{
  TH2 *h = pred->Predict(0).FakeData(pot).ToTH2(pot);
  for (int i=1; i<h->GetNbinsX()+1; i++) {
    for (int j=1; j<h->GetNbinsY()+1; j++) {
      double unc = sqrt(h->GetBinContent(i, j));
      h->SetBinContent(i, j, unc/h->GetBinContent(i, j));
    }
  }
  setHistAttr(h);
  h->SetTitle(title);
  h->SetName(name);
  return h;
}

void statsCheck(const char* outfile)
{
  rootlogon();

  TFile *fout = new TFile(outfile, "recreate");

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderGArFHC("/dune/data/users/sbjones/gasTpcCAF/v8/CAF_FHC.root", kBeam);
  SpectrumLoader loaderGArRHC("/dune/data/users/sbjones/gasTpcCAF/v8/CAF_RHC.root", kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  std::vector<const ISyst*> systlist = {};

  // Q2 vs category
  NoOscPredictionGenerator genFhcCat(axCatQ2, kPassND_FHC_NUMU && kIsTrueGasFV);
  PredictionInterp predFhcCat(systlist, 0, genFhcCat, loadersGArFHC);
  NoOscPredictionGenerator genRhcCat(axCatQ2, kPassND_RHC_NUMU && kIsTrueGasFV);
  PredictionInterp predRhcCat(systlist, 0, genRhcCat, loadersGArRHC);
  // 1D Q2 distributions
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

  loadersGArFHC.Go();
  loadersGArRHC.Go();

  fout->cd();
  // Make hists
  THStack *hsFhcQ2Cat1 = new THStack("hsFhcQ2Cat1", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(1).c_str()));
  THStack *hsFhcQ2Cat2 = new THStack("hsFhcQ2Cat2", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(2).c_str()));
  THStack *hsFhcQ2Cat3 = new THStack("hsFhcQ2Cat3", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(3).c_str()));
  THStack *hsFhcQ2Cat4 = new THStack("hsFhcQ2Cat4", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(4).c_str()));
  THStack *hsFhcQ2Cat5 = new THStack("hsFhcQ2Cat5", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(5).c_str()));
  THStack *hsFhcQ2Cat6 = new THStack("hsFhcQ2Cat6", Form("Statistical uncertainty for %s (FHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(6).c_str()));

  THStack *hsRhcQ2Cat1 = new THStack("hsRhcQ2Cat1", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(1).c_str()));
  THStack *hsRhcQ2Cat2 = new THStack("hsRhcQ2Cat2", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(2).c_str()));
  THStack *hsRhcQ2Cat3 = new THStack("hsRhcQ2Cat3", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(3).c_str()));
  THStack *hsRhcQ2Cat4 = new THStack("hsRhcQ2Cat4", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(4).c_str()));
  THStack *hsRhcQ2Cat5 = new THStack("hsRhcQ2Cat5", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(5).c_str()));
  THStack *hsRhcQ2Cat6 = new THStack("hsRhcQ2Cat6", Form("Statistical uncertainty for %s (RHC); Q^{2}_{reco} / (GeV)^{2}; Events", catName(6).c_str()));

  for (int y=1; y<=5; y++) {
    // 1D hists
    TH1 *hStatFhcQ2Cat1 = statUnc(&predFhcQ2RecoCat1, pot_nd*y, Form("hStatFhcQ2Cat1_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(1).c_str())); 
    TH1 *hStatFhcQ2Cat2 = statUnc(&predFhcQ2RecoCat2, pot_nd*y, Form("hStatFhcQ2Cat2_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(2).c_str()));
    TH1 *hStatFhcQ2Cat3 = statUnc(&predFhcQ2RecoCat3, pot_nd*y, Form("hStatFhcQ2Cat3_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(3).c_str()));
    TH1 *hStatFhcQ2Cat4 = statUnc(&predFhcQ2RecoCat4, pot_nd*y, Form("hStatFhcQ2Cat4_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(4).c_str()));
    TH1 *hStatFhcQ2Cat5 = statUnc(&predFhcQ2RecoCat5, pot_nd*y, Form("hStatFhcQ2Cat5_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(5).c_str()));
    TH1 *hStatFhcQ2Cat6 = statUnc(&predFhcQ2RecoCat6, pot_nd*y, Form("hStatFhcQ2Cat6_%dy", y), Form("Statistical uncertainty for %d years POT (FHC): %s", y, catName(6).c_str()));
    hStatFhcQ2Cat1->Write();
    hStatFhcQ2Cat2->Write();
    hStatFhcQ2Cat3->Write();
    hStatFhcQ2Cat4->Write();
    hStatFhcQ2Cat5->Write();
    hStatFhcQ2Cat6->Write();
    TH1 *hStatRhcQ2Cat1 = statUnc(&predRhcQ2RecoCat1, pot_nd*y, Form("hStatRhcQ2Cat1_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(1).c_str())); 
    TH1 *hStatRhcQ2Cat2 = statUnc(&predRhcQ2RecoCat2, pot_nd*y, Form("hStatRhcQ2Cat2_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(2).c_str()));
    TH1 *hStatRhcQ2Cat3 = statUnc(&predRhcQ2RecoCat3, pot_nd*y, Form("hStatRhcQ2Cat3_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(3).c_str()));
    TH1 *hStatRhcQ2Cat4 = statUnc(&predRhcQ2RecoCat4, pot_nd*y, Form("hStatRhcQ2Cat4_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(4).c_str()));
    TH1 *hStatRhcQ2Cat5 = statUnc(&predRhcQ2RecoCat5, pot_nd*y, Form("hStatRhcQ2Cat5_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(5).c_str()));
    TH1 *hStatRhcQ2Cat6 = statUnc(&predRhcQ2RecoCat6, pot_nd*y, Form("hStatRhcQ2Cat6_%dy", y), Form("Statistical uncertainty for %d years POT (RHC): %s", y, catName(6).c_str()));
    hStatRhcQ2Cat1->Write();
    hStatRhcQ2Cat2->Write();
    hStatRhcQ2Cat3->Write();
    hStatRhcQ2Cat4->Write();
    hStatRhcQ2Cat5->Write();
    hStatRhcQ2Cat6->Write();
    // 2D hists
    TH2 *h2StatFhcCat = statUnc2d(&predFhcCat, pot_nd*y, Form("h2StatFhcCat_%dy", y), Form("Statistical uncertainty for %d years POT (FHC)", y));
    TH2 *h2StatRhcCat = statUnc2d(&predRhcCat, pot_nd*y, Form("h2StatRhcCat_%dy", y), Form("Statistical uncertainty for %d years POT (RHC)", y));
    h2StatFhcCat->Write();
    h2StatRhcCat->Write();

    hStatFhcQ2Cat1->SetLineColor(52+y*4);
    hStatFhcQ2Cat2->SetLineColor(52+y*4);
    hStatFhcQ2Cat3->SetLineColor(52+y*4);
    hStatFhcQ2Cat4->SetLineColor(52+y*4);
    hStatFhcQ2Cat5->SetLineColor(52+y*4);
    hStatFhcQ2Cat6->SetLineColor(52+y*4);
    hsFhcQ2Cat1->Add(hStatFhcQ2Cat1);
    hsFhcQ2Cat2->Add(hStatFhcQ2Cat2);
    hsFhcQ2Cat3->Add(hStatFhcQ2Cat3);
    hsFhcQ2Cat4->Add(hStatFhcQ2Cat4);
    hsFhcQ2Cat5->Add(hStatFhcQ2Cat5);
    hsFhcQ2Cat6->Add(hStatFhcQ2Cat6);

    hStatRhcQ2Cat1->SetLineColor(52+y*4);
    hStatRhcQ2Cat2->SetLineColor(52+y*4);
    hStatRhcQ2Cat3->SetLineColor(52+y*4);
    hStatRhcQ2Cat4->SetLineColor(52+y*4);
    hStatRhcQ2Cat5->SetLineColor(52+y*4);
    hStatRhcQ2Cat6->SetLineColor(52+y*4);
    hsRhcQ2Cat1->Add(hStatRhcQ2Cat1);
    hsRhcQ2Cat2->Add(hStatRhcQ2Cat2);
    hsRhcQ2Cat3->Add(hStatRhcQ2Cat3);
    hsRhcQ2Cat4->Add(hStatRhcQ2Cat4);
    hsRhcQ2Cat5->Add(hStatRhcQ2Cat5);
    hsRhcQ2Cat6->Add(hStatRhcQ2Cat6);
  }

  hsFhcQ2Cat1->Write();
  hsFhcQ2Cat2->Write();
  hsFhcQ2Cat3->Write();
  hsFhcQ2Cat4->Write();
  hsFhcQ2Cat5->Write();
  hsFhcQ2Cat6->Write();

  hsRhcQ2Cat1->Write();
  hsRhcQ2Cat2->Write();
  hsRhcQ2Cat3->Write();
  hsRhcQ2Cat4->Write();
  hsRhcQ2Cat5->Write();
  hsRhcQ2Cat6->Write();


  fout->Close();
} // statsCheck
