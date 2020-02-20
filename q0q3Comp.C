// Looks at the 2D distributions of q0 and q3 and makes some reweightable hists
#include "StandardRecord/StandardRecord.h"
#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
// CAFAna includes
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Systs/XSecSysts.h"
#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/TDRLoaders.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/common_fit_definitions.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Experiment/ReactorExperiment.h"
// // ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"

#include "Utilities/rootlogon.C"

using namespace ana;

// Useful functions
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

TH2* makeHist2D(const char* name, const char* title, PredictionInterp *pred, const double pot,
		SystShifts shift = kNoShift)
{
  TH2 *h2 = pred->PredictSyst(0, shift).FakeData(pot).ToTH2(pot);
  h2->SetName(name);
  h2->SetTitle(title);
  setHistAttr(h2);
  return h2;
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
    h2->GetXaxis()->SetBinLabel(binX, catName(binX).c_str());
    h2->GetYaxis()->SetBinLabel(binX, catName(binX).c_str());
    for (int binY=1; binY<h2->GetNbinsY()+1; binY++) {
      h2->SetBinContent(binX, binY, h2->GetBinContent(binX, binY) / colInt);
    }
  }
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

TH2 * ratioSuppressLowStats(PredictionInterp *pred, const char* namestub, const char* title,
			    const double inputPOT, const double truePOT, SystShifts shift)
{
  TH2 *hgenie = pred->PredictSyst(0, kNoShift).FakeData(inputPOT).ToTH2(inputPOT);
  TH2 *hnuwro = pred->PredictSyst(0, shift).FakeData(inputPOT).ToTH2(inputPOT);
  hgenie->SetTitle(title);
  hnuwro->SetTitle(title);
  hgenie->Write(Form("%s_genie", namestub));
  hnuwro->Write(Form("%s_nuwro", namestub));

  for (int i=1; i<hgenie->GetNbinsX(); i++) {
    for (int j=1; j<hgenie->GetNbinsY(); j++) {
      if (hgenie->GetBinContent(i, j) * (truePOT/inputPOT) < 100.) {
	// Set bin content to 0
	hgenie->SetBinContent(i, j, 0);
	hgenie->SetBinError(i, j, 0);
	hnuwro->SetBinContent(i, j, 0);
	hnuwro->SetBinError(i, j, 0);
      }
    }
  }

  hnuwro->Divide(hgenie);
  hnuwro->SetTitle(Form("%s: NuWro/GENIE; q_{0, reco} / GeV; q_{3, reco} / GeV; NuWro/GENIE", title));
  hnuwro->SetName(Form("%s_ratio", namestub));
  return hnuwro;
}

// Vars
const Var kRecoEnergyND  = SIMPLEVAR(dune.Ev_reco);
const Var kTrueEnergy    = SIMPLEVAR(dune.Ev);
const Var kTrueLepEnergy = SIMPLEVAR(dune.LepE);
const Var kFHC = SIMPLEVAR(dune.isFHC);
const Var isCC = SIMPLEVAR(dune.isCC);
const Var kRecoNumu = SIMPLEVAR(dune.reco_numu);
// Reco multiplicities
const Var kRecoPipl      = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var kRecoPimin     = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi0       = SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi        = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
// True multiplicities
const Var kPi0       = SIMPLEVAR(dune.nipi0);
const Var kChargedPi = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip);
const Var kPi        = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);

const double mmu = 0.10566; // GeV/c^2

// True and reconstructed q0 and q3
const Var kq0 = SIMPLEVAR(dune.Ev) - SIMPLEVAR(dune.LepE);
const Var kq3({"dune.LepMomX", "dune.LepMomY", "dune.LepMomZ", 
      "dune.NuMomX", "dune.NuMomY", "dune.NuMomZ"},
	      [](const caf::StandardRecord* sr) {
		double q3 = 0.;
		if (sr->dune.LepMomZ != 0.) {
		  q3 = sqrt(pow(sr->dune.NuMomX-sr->dune.LepMomX, 2) + pow(sr->dune.NuMomY-sr->dune.LepMomY, 2) + pow(sr->dune.NuMomZ-sr->dune.LepMomZ, 2));
		}
		return q3;		
	      });

const Var kq0Reco = SIMPLEVAR(dune.Ev_reco) - SIMPLEVAR(dune.Elep_reco);
const Var kq3Reco({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		  [](const caf::StandardRecord* sr) {
		    double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		    //double q3 = pmu*pmu + sr->dune.Ev_reco*sr->dune.Ev_reco - 2*sr->dune.Ev_reco*pmu*TMath::Cos(sr->dune.theta_reco);
		    // Calculate Q2 and then use that to get q3
		    double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		    double q0 = sr->dune.Ev_reco - sr->dune.Elep_reco;
		    double q3 = sqrt(q2+q0*q0);
		    return q3;
		  });

// Energy bin edges
std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
				 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};

const std::vector<double> q0q3Edges = {0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 
				       1., 1.25, 1.5, 1.75,
				       2., 2.25, 2.5, 2.75,
				       3., 3.25, 3.5, 3.75,
				       4., 5., 6., 8.};
const Binning q0q3Bins = Binning::Custom(q0q3Edges);
const HistAxis axq0q3Reco("q_{0, reco} / GeV", q0q3Bins, kq0Reco, 
			  "q_{3, reco} / GeV", q0q3Bins, kq3Reco);
const HistAxis axq0q3("q_{0, true} / GeV", q0q3Bins, kq0, 
		      "q_{3, true} / GeV", q0q3Bins, kq3);
const HistAxis axq0Comp("q_{0, true} / GeV", q0q3Bins, kq0, 
			"q_{0, reco} / GeV", q0q3Bins, kq0Reco);
const HistAxis axq3Comp("q_{3, true} / GeV", q0q3Bins, kq3, 
			"q_{3, reco} / GeV", q0q3Bins, kq3Reco);

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

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// True MC files POT
const double fhcPOT = 1.9342e20;
const double rhcPOT = 3.9302e20;

void q0q3Comp(const char* outfile,
	      const char* cafs="/dune/data/users/sbjones/gasTpcCAF/v8/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  std::vector<const ISyst*> systlist = GetXSecSysts({"NuWroReweightFakeData"});
  std::vector<const ISyst*> onesyst = {};
  onesyst.push_back(systlist.at(systlist.size()-1));
  SystShifts fakedata(systlist.at(systlist.size()-1), 1);

  Loaders loadersFHC, loadersRHC;
  SpectrumLoader loaderFHC(Form("%s/CAF_FHC.root", cafs), kBeam);
  SpectrumLoader loaderRHC(Form("%s/CAF_RHC.root", cafs), kBeam);
  loadersFHC.AddLoader(&loaderFHC, caf::kNEARDET, Loaders::kMC);
  loadersRHC.AddLoader(&loaderRHC, caf::kNEARDET, Loaders::kMC);

  NoOscPredictionGenerator genq0CompFhc(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq3CompFhc(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3Fhc(axq0q3, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3RecoFhc(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV); 
  PredictionInterp predq0CompFhc(onesyst, 0, genq0CompFhc, loadersFHC);
  PredictionInterp predq3CompFhc(onesyst, 0, genq3CompFhc, loadersFHC);
  PredictionInterp predq0q3Fhc(onesyst, 0, genq0q3Fhc, loadersFHC);
  PredictionInterp predq0q3RecoFhc(onesyst, 0, genq0q3RecoFhc, loadersFHC);

  // Now do the same for various pion multiplicities
  NoOscPredictionGenerator genq0q3RecoFhcCat1(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3RecoFhcCat2(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3RecoFhcCat3(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3RecoFhcCat4(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3RecoFhcCat5(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3RecoFhcCat6(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);

  NoOscPredictionGenerator genq0q3FhcCat1(axq0q3, kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genq0q3FhcCat2(axq0q3, kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genq0q3FhcCat3(axq0q3, kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genq0q3FhcCat4(axq0q3, kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genq0q3FhcCat5(axq0q3, kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genq0q3FhcCat6(axq0q3, kIsTrueGasFV && kPassCat6);

  NoOscPredictionGenerator genq0q3RecoRhcCat1(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3RecoRhcCat2(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3RecoRhcCat3(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3RecoRhcCat4(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3RecoRhcCat5(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3RecoRhcCat6(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);

  NoOscPredictionGenerator genq0q3RhcCat1(axq0q3, kIsTrueGasFV && kPassCat1);
  NoOscPredictionGenerator genq0q3RhcCat2(axq0q3, kIsTrueGasFV && kPassCat2);
  NoOscPredictionGenerator genq0q3RhcCat3(axq0q3, kIsTrueGasFV && kPassCat3);
  NoOscPredictionGenerator genq0q3RhcCat4(axq0q3, kIsTrueGasFV && kPassCat4);
  NoOscPredictionGenerator genq0q3RhcCat5(axq0q3, kIsTrueGasFV && kPassCat5);
  NoOscPredictionGenerator genq0q3RhcCat6(axq0q3, kIsTrueGasFV && kPassCat6);

  PredictionInterp predq0q3RecoFhcCat1(onesyst, 0, genq0q3RecoFhcCat1, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat2(onesyst, 0, genq0q3RecoFhcCat2, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat3(onesyst, 0, genq0q3RecoFhcCat3, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat4(onesyst, 0, genq0q3RecoFhcCat4, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat5(onesyst, 0, genq0q3RecoFhcCat5, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat6(onesyst, 0, genq0q3RecoFhcCat6, loadersFHC);

  PredictionInterp predq0q3FhcCat1(onesyst, 0, genq0q3FhcCat1, loadersFHC);
  PredictionInterp predq0q3FhcCat2(onesyst, 0, genq0q3FhcCat2, loadersFHC);
  PredictionInterp predq0q3FhcCat3(onesyst, 0, genq0q3FhcCat3, loadersFHC);
  PredictionInterp predq0q3FhcCat4(onesyst, 0, genq0q3FhcCat4, loadersFHC);
  PredictionInterp predq0q3FhcCat5(onesyst, 0, genq0q3FhcCat5, loadersFHC);
  PredictionInterp predq0q3FhcCat6(onesyst, 0, genq0q3FhcCat6, loadersFHC);

  PredictionInterp predq0q3RecoRhcCat1(onesyst, 0, genq0q3RecoRhcCat1, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat2(onesyst, 0, genq0q3RecoRhcCat2, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat3(onesyst, 0, genq0q3RecoRhcCat3, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat4(onesyst, 0, genq0q3RecoRhcCat4, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat5(onesyst, 0, genq0q3RecoRhcCat5, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat6(onesyst, 0, genq0q3RecoRhcCat6, loadersRHC);

  PredictionInterp predq0q3RhcCat1(onesyst, 0, genq0q3RhcCat1, loadersRHC);
  PredictionInterp predq0q3RhcCat2(onesyst, 0, genq0q3RhcCat2, loadersRHC);
  PredictionInterp predq0q3RhcCat3(onesyst, 0, genq0q3RhcCat3, loadersRHC);
  PredictionInterp predq0q3RhcCat4(onesyst, 0, genq0q3RhcCat4, loadersRHC);
  PredictionInterp predq0q3RhcCat5(onesyst, 0, genq0q3RhcCat5, loadersRHC);
  PredictionInterp predq0q3RhcCat6(onesyst, 0, genq0q3RhcCat6, loadersRHC);

  loadersFHC.Go();
  loadersRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  TH2* h2q0CompFhc   = makeHist2D("h2q0CompFhc", "q_{0} comparison for HPgTPC", &predq0CompFhc, pot_nd);
  TH2* h2q3CompFhc   = makeHist2D("h2q3CompFhc", "q_{3} comparison for HPgTPC", &predq3CompFhc, pot_nd);
  TH2* h2q0q3Fhc     = makeHist2D("h2q0q3Fhc", "q_{0}, q_{3} in HPgTPC", &predq0q3Fhc, pot_nd);
  TH2* h2q0q3RecoFhc = makeHist2D("h2q0q3RecoFhc", "q_{0, reco}, q_{3, reco} in HPgTPC", &predq0q3RecoFhc, pot_nd);
  h2q0CompFhc->Write();
  h2q3CompFhc->Write();
  h2q0q3Fhc->Write();
  h2q0q3RecoFhc->Write();

  TH2 *hq0q3RecoFhcCat1 = ratioSuppressLowStats(&predq0q3RecoFhcCat1, "hq0q3RecoFhcCat1", Form("FHC, reco %s", catName(1).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat2 = ratioSuppressLowStats(&predq0q3RecoFhcCat2, "hq0q3RecoFhcCat2", Form("FHC, reco %s", catName(2).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat3 = ratioSuppressLowStats(&predq0q3RecoFhcCat3, "hq0q3RecoFhcCat3", Form("FHC, reco %s", catName(3).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat4 = ratioSuppressLowStats(&predq0q3RecoFhcCat4, "hq0q3RecoFhcCat4", Form("FHC, reco %s", catName(4).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat5 = ratioSuppressLowStats(&predq0q3RecoFhcCat5, "hq0q3RecoFhcCat5", Form("FHC, reco %s", catName(5).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat6 = ratioSuppressLowStats(&predq0q3RecoFhcCat6, "hq0q3RecoFhcCat6", Form("FHC, reco %s", catName(6).c_str()) , pot_nd, fhcPOT, fakedata);

  TH2 *hq0q3RecoRhcCat1 = ratioSuppressLowStats(&predq0q3RecoRhcCat1, "hq0q3RecoRhcCat1", Form("RHC, reco %s", catName(1).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat2 = ratioSuppressLowStats(&predq0q3RecoRhcCat2, "hq0q3RecoRhcCat2", Form("RHC, reco %s", catName(2).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat3 = ratioSuppressLowStats(&predq0q3RecoRhcCat3, "hq0q3RecoRhcCat3", Form("RHC, reco %s", catName(3).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat4 = ratioSuppressLowStats(&predq0q3RecoRhcCat4, "hq0q3RecoRhcCat4", Form("RHC, reco %s", catName(4).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat5 = ratioSuppressLowStats(&predq0q3RecoRhcCat5, "hq0q3RecoRhcCat5", Form("RHC, reco %s", catName(5).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat6 = ratioSuppressLowStats(&predq0q3RecoRhcCat6, "hq0q3RecoRhcCat6", Form("RHC, reco %s", catName(6).c_str()) , pot_nd, rhcPOT, fakedata);

  TH2 *hq0q3FhcCat1 = ratioSuppressLowStats(&predq0q3FhcCat1, "hq0q3FhcCat1", Form("FHC, true %s", catName(1).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat2 = ratioSuppressLowStats(&predq0q3FhcCat2, "hq0q3FhcCat2", Form("FHC, true %s", catName(2).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat3 = ratioSuppressLowStats(&predq0q3FhcCat3, "hq0q3FhcCat3", Form("FHC, true %s", catName(3).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat4 = ratioSuppressLowStats(&predq0q3FhcCat4, "hq0q3FhcCat4", Form("FHC, true %s", catName(4).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat5 = ratioSuppressLowStats(&predq0q3FhcCat5, "hq0q3FhcCat5", Form("FHC, true %s", catName(5).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat6 = ratioSuppressLowStats(&predq0q3FhcCat6, "hq0q3FhcCat6", Form("FHC, true %s", catName(6).c_str()) , pot_nd, fhcPOT, fakedata);

  TH2 *hq0q3RhcCat1 = ratioSuppressLowStats(&predq0q3RhcCat1, "hq0q3RhcCat1", Form("RHC, true %s", catName(1).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat2 = ratioSuppressLowStats(&predq0q3RhcCat2, "hq0q3RhcCat2", Form("RHC, true %s", catName(2).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat3 = ratioSuppressLowStats(&predq0q3RhcCat3, "hq0q3RhcCat3", Form("RHC, true %s", catName(3).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat4 = ratioSuppressLowStats(&predq0q3RhcCat4, "hq0q3RhcCat4", Form("RHC, true %s", catName(4).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat5 = ratioSuppressLowStats(&predq0q3RhcCat5, "hq0q3RhcCat5", Form("RHC, true %s", catName(5).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat6 = ratioSuppressLowStats(&predq0q3RhcCat6, "hq0q3RhcCat6", Form("RHC, true %s", catName(6).c_str()) , pot_nd, rhcPOT, fakedata);

  hq0q3RecoFhcCat1->Write();
  hq0q3RecoFhcCat2->Write();
  hq0q3RecoFhcCat3->Write();
  hq0q3RecoFhcCat4->Write();
  hq0q3RecoFhcCat5->Write();
  hq0q3RecoFhcCat6->Write();

  hq0q3RecoRhcCat1->Write();
  hq0q3RecoRhcCat2->Write();
  hq0q3RecoRhcCat3->Write();
  hq0q3RecoRhcCat4->Write();
  hq0q3RecoRhcCat5->Write();
  hq0q3RecoRhcCat6->Write();

  hq0q3FhcCat1->Write();
  hq0q3FhcCat2->Write();
  hq0q3FhcCat3->Write();
  hq0q3FhcCat4->Write();
  hq0q3FhcCat5->Write();
  hq0q3FhcCat6->Write();

  hq0q3RhcCat1->Write();
  hq0q3RhcCat2->Write();
  hq0q3RhcCat3->Write();
  hq0q3RhcCat4->Write();
  hq0q3RhcCat5->Write();
  hq0q3RhcCat6->Write();

  fout->Close();
  delete fout;
} // q0q3Comp
