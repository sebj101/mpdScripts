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
  h2->Scale(1./h2->Integral());
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
  else if (cat==7) name="0#pi-1p";
  else if (cat==8) name="0#pi->1p";
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

TH2* ratioSuppressLowStats(PredictionInterp *pred, const char* namestub, const char* title,
			    const double inputPOT, const double truePOT, SystShifts shift)
{
  TH2 *hgenie = pred->PredictSyst(0, kNoShift).FakeData(inputPOT).ToTH2(inputPOT);
  TH2 *hnuwro = pred->PredictSyst(0, shift).FakeData(inputPOT).ToTH2(inputPOT);
  hgenie->SetTitle(Form("%s GENIE; q_{3, reco} / GeV; q_{0, reco} / GeV; Events / (GeV)^{2}", title));
  hnuwro->SetTitle(Form("%s NuWro; q_{3, reco} / GeV; q_{0, reco} / GeV; Events / (GeV)^{2}", title));
  setHistAttr(hgenie);
  setHistAttr(hnuwro);
  hgenie->Scale(1, "width");
  hnuwro->Scale(1, "width");
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
  hnuwro->SetTitle(Form("%s: NuWro/GENIE; q_{3, reco} / GeV; q_{0, reco} / GeV; NuWro/GENIE", title));
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
const Var kRecoP         = SIMPLEVAR(dune.gastpc_pro_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi        = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
// True multiplicities
const Var kPi0       = SIMPLEVAR(dune.nipi0);
const Var kChargedPi = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip);
const Var kPi        = SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipi0);

const double mmu   = 0.10566; // GeV/c^2
const double mpipm = 0.13957;
const double mpi0  = 0.13498;

// True and reconstructed q0 and q3
const Var kq0({"dune.ePip", "dune.ePim", "dune.ePi0", "dune.eP", "dune.nipi0", "dune.nipim", "dune.nipip"},
	      [](const caf::StandardRecord* sr) {
		// Define as MINERvA available energy
		double q0 = sr->dune.ePi0 + sr->dune.nipi0*mpi0 + sr->dune.ePip + sr->dune.ePim + (sr->dune.nipip + sr->dune.nipim)*mpipm + sr->dune.eP;
		return q0;
	      });

const Var kq3({"dune.ePip", "dune.ePim", "dune.ePi0", "dune.eP", "dune.nipi0", "dune.nipim", "dune.nipip", "dune.LepMomZ", "dune.LepE", "dune.LepNuAngle"},
    //"dune.LepMomX", "dune.LepMomY", "dune.LepMomZ", "dune.NuMomX", "dune.NuMomY", "dune.NuMomZ"},
	      [](const caf::StandardRecord* sr) {
		double q3 = 0.;
		if (sr->dune.LepMomZ != 0.) {
		  //q3 = sqrt(pow(sr->dune.NuMomX-sr->dune.LepMomX, 2) + pow(sr->dune.NuMomY-sr->dune.LepMomY, 2) + pow(sr->dune.NuMomZ-sr->dune.LepMomZ, 2));
		  // Define as MINERvA available energy
		  double q0 = sr->dune.ePi0 + sr->dune.nipi0*mpi0 + sr->dune.ePip + sr->dune.ePim + (sr->dune.nipip + sr->dune.nipim)*mpipm + sr->dune.eP;
		  double pmu = sqrt( sr->dune.LepE*sr->dune.LepE - mmu*mmu );
		  double Q2 = 2*(q0+sr->dune.LepE) * (sr->dune.LepE - pmu*TMath::Cos(sr->dune.LepNuAngle)) - mmu*mmu;
		  q3 = sqrt( Q2 + q0*q0 );
		}
		return q3;		
	      });
// Reco kinematic quantities
const Var kq0Reco = SIMPLEVAR(dune.Ev_reco) - SIMPLEVAR(dune.Elep_reco);
const Var kq3Reco({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
                  [](const caf::StandardRecord* sr) {
		    double pmu = sqrt(sr->dune.Elep_reco*sr->dune.Elep_reco - mmu*mmu);
		    // Calculate Q2 and then use that to get q3
		    double q2 = 2 * sr->dune.Ev_reco * (sr->dune.Elep_reco - pmu*TMath::Cos(sr->dune.theta_reco)) - mmu*mmu;
		    double q0 = sr->dune.Ev_reco - sr->dune.Elep_reco;
		    double q3 = sqrt( q2+q0*q0 );
		    return q3;
		  });

// Energy bin edges
std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
				 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};

const std::vector<double> q0q3Edges = {0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 
				       1., 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875,
				       2., 2.25, 2.5, 2.75,
				       3., 3.25, 3.5, 3.75,
				       4., 5., 6.};
const Binning q0q3Bins = Binning::Custom(q0q3Edges);
const HistAxis axq0q3Reco("q_{3, reco} / GeV", q0q3Bins, kq3Reco,
			  "q_{0, reco} / GeV", q0q3Bins, kq0Reco);
const HistAxis axq0q3("q_{3, true} / GeV", q0q3Bins, kq3,
		      "q_{0, true} / GeV", q0q3Bins, kq0);
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
const Cut kPassCat7({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==0 && sr->dune.nP==1
			      && abs(sr->dune.LepPDG)==13);
		    });
const Cut kPassCat8({},
		    [](const caf::StandardRecord* sr)
		    {
		      return (sr->dune.nipip+sr->dune.nipim+sr->dune.nipi0==0 && sr->dune.nP>1
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
			[](const caf::StandardRecord*sr)
			{
			  return (sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult>2);
			});
const Cut kPassRecoCat7({},
			[](const caf::StandardRecord* sr)
			{
			  return ((sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==0) && sr->dune.gastpc_pro_mult==1);
			});
const Cut kPassRecoCat8({},
			[](const caf::StandardRecord* sr)
			{
			  return ((sr->dune.gastpc_pi_pl_mult+sr->dune.gastpc_pi_min_mult+sr->dune.gastpc_pi_0_mult==0) && sr->dune.gastpc_pro_mult>1);
			});

// POT for n years
const double years = 1.; // of POT
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// True MC files POT
const double fhcPOT = 1.9342e20;
const double rhcPOT = 3.9302e20;

void q0q3Comp(const char* outfile,
	      const char* cafs="/dune/data/users/sbjones/gasTpcCAF/v9/")
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

  // FHC
  NoOscPredictionGenerator genq0CompFhc(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq3CompFhc(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3Fhc(axq0q3, kPassND_FHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3RecoFhc(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV); 
  PredictionInterp predq0CompFhc(onesyst, 0, genq0CompFhc, loadersFHC);
  PredictionInterp predq3CompFhc(onesyst, 0, genq3CompFhc, loadersFHC);
  PredictionInterp predq0q3Fhc(onesyst, 0, genq0q3Fhc, loadersFHC);
  PredictionInterp predq0q3RecoFhc(onesyst, 0, genq0q3RecoFhc, loadersFHC);
  // RHC
  NoOscPredictionGenerator genq0CompRhc(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq3CompRhc(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3Rhc(axq0q3, kPassND_RHC_NUMU && kIsTrueGasFV); 
  NoOscPredictionGenerator genq0q3RecoRhc(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV); 
  PredictionInterp predq0CompRhc(onesyst, 0, genq0CompRhc, loadersRHC);
  PredictionInterp predq3CompRhc(onesyst, 0, genq3CompRhc, loadersRHC);
  PredictionInterp predq0q3Rhc(onesyst, 0, genq0q3Rhc, loadersRHC);
  PredictionInterp predq0q3RecoRhc(onesyst, 0, genq0q3RecoRhc, loadersRHC);

  // FHC true reco comps 
  NoOscPredictionGenerator genq0CompFhcCat1(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1); 
  NoOscPredictionGenerator genq0CompFhcCat2(axq0Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2); 
  NoOscPredictionGenerator genq0CompFhcCat3(axq0Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3); 
  NoOscPredictionGenerator genq0CompFhcCat4(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4); 
  NoOscPredictionGenerator genq0CompFhcCat5(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5); 
  NoOscPredictionGenerator genq0CompFhcCat6(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0CompFhcCat7(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat7);
  NoOscPredictionGenerator genq0CompFhcCat8(axq0Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat8);
  PredictionInterp predq0CompFhcCat1(onesyst, 0, genq0CompFhcCat1, loadersFHC);
  PredictionInterp predq0CompFhcCat2(onesyst, 0, genq0CompFhcCat2, loadersFHC);
  PredictionInterp predq0CompFhcCat3(onesyst, 0, genq0CompFhcCat3, loadersFHC);
  PredictionInterp predq0CompFhcCat4(onesyst, 0, genq0CompFhcCat4, loadersFHC);
  PredictionInterp predq0CompFhcCat5(onesyst, 0, genq0CompFhcCat5, loadersFHC);
  PredictionInterp predq0CompFhcCat6(onesyst, 0, genq0CompFhcCat6, loadersFHC);
  PredictionInterp predq0CompFhcCat7(onesyst, 0, genq0CompFhcCat7, loadersFHC);
  PredictionInterp predq0CompFhcCat8(onesyst, 0, genq0CompFhcCat8, loadersFHC);
  NoOscPredictionGenerator genq3CompFhcCat1(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1); 
  NoOscPredictionGenerator genq3CompFhcCat2(axq3Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2); 
  NoOscPredictionGenerator genq3CompFhcCat3(axq3Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3); 
  NoOscPredictionGenerator genq3CompFhcCat4(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq3CompFhcCat5(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5); 
  NoOscPredictionGenerator genq3CompFhcCat6(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6); 
  NoOscPredictionGenerator genq3CompFhcCat7(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat7); 
  NoOscPredictionGenerator genq3CompFhcCat8(axq3Comp, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat8); 
  PredictionInterp predq3CompFhcCat1(onesyst, 0, genq3CompFhcCat1, loadersFHC);
  PredictionInterp predq3CompFhcCat2(onesyst, 0, genq3CompFhcCat2, loadersFHC);
  PredictionInterp predq3CompFhcCat3(onesyst, 0, genq3CompFhcCat3, loadersFHC);
  PredictionInterp predq3CompFhcCat4(onesyst, 0, genq3CompFhcCat4, loadersFHC);
  PredictionInterp predq3CompFhcCat5(onesyst, 0, genq3CompFhcCat5, loadersFHC);
  PredictionInterp predq3CompFhcCat6(onesyst, 0, genq3CompFhcCat6, loadersFHC);
  PredictionInterp predq3CompFhcCat7(onesyst, 0, genq3CompFhcCat7, loadersFHC);
  PredictionInterp predq3CompFhcCat8(onesyst, 0, genq3CompFhcCat8, loadersFHC);

  // RHC true reco comps
  NoOscPredictionGenerator genq0CompRhcCat1(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1); 
  NoOscPredictionGenerator genq0CompRhcCat2(axq0Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2); 
  NoOscPredictionGenerator genq0CompRhcCat3(axq0Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3); 
  NoOscPredictionGenerator genq0CompRhcCat4(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4); 
  NoOscPredictionGenerator genq0CompRhcCat5(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5); 
  NoOscPredictionGenerator genq0CompRhcCat6(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0CompRhcCat7(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat7); 
  NoOscPredictionGenerator genq0CompRhcCat8(axq0Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat8);
  PredictionInterp predq0CompRhcCat1(onesyst, 0, genq0CompRhcCat1, loadersRHC);
  PredictionInterp predq0CompRhcCat2(onesyst, 0, genq0CompRhcCat2, loadersRHC);
  PredictionInterp predq0CompRhcCat3(onesyst, 0, genq0CompRhcCat3, loadersRHC);
  PredictionInterp predq0CompRhcCat4(onesyst, 0, genq0CompRhcCat4, loadersRHC);
  PredictionInterp predq0CompRhcCat5(onesyst, 0, genq0CompRhcCat5, loadersRHC);
  PredictionInterp predq0CompRhcCat6(onesyst, 0, genq0CompRhcCat6, loadersRHC);
  PredictionInterp predq0CompRhcCat7(onesyst, 0, genq0CompRhcCat7, loadersRHC);
  PredictionInterp predq0CompRhcCat8(onesyst, 0, genq0CompRhcCat8, loadersRHC);
  NoOscPredictionGenerator genq3CompRhcCat1(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1); 
  NoOscPredictionGenerator genq3CompRhcCat2(axq3Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2); 
  NoOscPredictionGenerator genq3CompRhcCat3(axq3Comp, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3); 
  NoOscPredictionGenerator genq3CompRhcCat4(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq3CompRhcCat5(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5); 
  NoOscPredictionGenerator genq3CompRhcCat6(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6); 
  NoOscPredictionGenerator genq3CompRhcCat7(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat7); 
  NoOscPredictionGenerator genq3CompRhcCat8(axq3Comp, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat8); 
  PredictionInterp predq3CompRhcCat1(onesyst, 0, genq3CompRhcCat1, loadersRHC);
  PredictionInterp predq3CompRhcCat2(onesyst, 0, genq3CompRhcCat2, loadersRHC);
  PredictionInterp predq3CompRhcCat3(onesyst, 0, genq3CompRhcCat3, loadersRHC);
  PredictionInterp predq3CompRhcCat4(onesyst, 0, genq3CompRhcCat4, loadersRHC);
  PredictionInterp predq3CompRhcCat5(onesyst, 0, genq3CompRhcCat5, loadersRHC);
  PredictionInterp predq3CompRhcCat6(onesyst, 0, genq3CompRhcCat6, loadersRHC);
  PredictionInterp predq3CompRhcCat7(onesyst, 0, genq3CompRhcCat7, loadersRHC);
  PredictionInterp predq3CompRhcCat8(onesyst, 0, genq3CompRhcCat8, loadersRHC);

  // Now do the same for various pion multiplicities
  NoOscPredictionGenerator genq0q3RecoFhcCat1(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3RecoFhcCat2(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3RecoFhcCat3(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3RecoFhcCat4(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3RecoFhcCat5(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3RecoFhcCat6(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0q3RecoFhcCat7(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat7);
  NoOscPredictionGenerator genq0q3RecoFhcCat8(axq0q3Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat8);

  NoOscPredictionGenerator genq0q3FhcCat1(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3FhcCat2(axq0q3,kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3FhcCat3(axq0q3,kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3FhcCat4(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3FhcCat5(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3FhcCat6(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0q3FhcCat7(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat7);
  NoOscPredictionGenerator genq0q3FhcCat8(axq0q3,kPassND_FHC_NUMU && kIsTrueGasFV && kPassRecoCat8);

  NoOscPredictionGenerator genq0q3RecoRhcCat1(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3RecoRhcCat2(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3RecoRhcCat3(axq0q3Reco, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3RecoRhcCat4(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3RecoRhcCat5(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3RecoRhcCat6(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0q3RecoRhcCat7(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat7);
  NoOscPredictionGenerator genq0q3RecoRhcCat8(axq0q3Reco, kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat8);

  NoOscPredictionGenerator genq0q3RhcCat1(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat1);
  NoOscPredictionGenerator genq0q3RhcCat2(axq0q3, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat2);
  NoOscPredictionGenerator genq0q3RhcCat3(axq0q3, kRecoNumu==1 && kIsTrueGasFV && kPassRecoCat3);
  NoOscPredictionGenerator genq0q3RhcCat4(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat4);
  NoOscPredictionGenerator genq0q3RhcCat5(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat5);
  NoOscPredictionGenerator genq0q3RhcCat6(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat6);
  NoOscPredictionGenerator genq0q3RhcCat7(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat7);
  NoOscPredictionGenerator genq0q3RhcCat8(axq0q3,kPassND_RHC_NUMU && kIsTrueGasFV && kPassRecoCat8);

  PredictionInterp predq0q3RecoFhcCat1(onesyst, 0, genq0q3RecoFhcCat1, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat2(onesyst, 0, genq0q3RecoFhcCat2, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat3(onesyst, 0, genq0q3RecoFhcCat3, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat4(onesyst, 0, genq0q3RecoFhcCat4, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat5(onesyst, 0, genq0q3RecoFhcCat5, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat6(onesyst, 0, genq0q3RecoFhcCat6, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat7(onesyst, 0, genq0q3RecoFhcCat7, loadersFHC);
  PredictionInterp predq0q3RecoFhcCat8(onesyst, 0, genq0q3RecoFhcCat8, loadersFHC);

  PredictionInterp predq0q3FhcCat1(onesyst, 0, genq0q3FhcCat1, loadersFHC);
  PredictionInterp predq0q3FhcCat2(onesyst, 0, genq0q3FhcCat2, loadersFHC);
  PredictionInterp predq0q3FhcCat3(onesyst, 0, genq0q3FhcCat3, loadersFHC);
  PredictionInterp predq0q3FhcCat4(onesyst, 0, genq0q3FhcCat4, loadersFHC);
  PredictionInterp predq0q3FhcCat5(onesyst, 0, genq0q3FhcCat5, loadersFHC);
  PredictionInterp predq0q3FhcCat6(onesyst, 0, genq0q3FhcCat6, loadersFHC);
  PredictionInterp predq0q3FhcCat7(onesyst, 0, genq0q3FhcCat7, loadersFHC);
  PredictionInterp predq0q3FhcCat8(onesyst, 0, genq0q3FhcCat8, loadersFHC);

  PredictionInterp predq0q3RecoRhcCat1(onesyst, 0, genq0q3RecoRhcCat1, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat2(onesyst, 0, genq0q3RecoRhcCat2, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat3(onesyst, 0, genq0q3RecoRhcCat3, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat4(onesyst, 0, genq0q3RecoRhcCat4, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat5(onesyst, 0, genq0q3RecoRhcCat5, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat6(onesyst, 0, genq0q3RecoRhcCat6, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat7(onesyst, 0, genq0q3RecoRhcCat7, loadersRHC);
  PredictionInterp predq0q3RecoRhcCat8(onesyst, 0, genq0q3RecoRhcCat8, loadersRHC);

  PredictionInterp predq0q3RhcCat1(onesyst, 0, genq0q3RhcCat1, loadersRHC);
  PredictionInterp predq0q3RhcCat2(onesyst, 0, genq0q3RhcCat2, loadersRHC);
  PredictionInterp predq0q3RhcCat3(onesyst, 0, genq0q3RhcCat3, loadersRHC);
  PredictionInterp predq0q3RhcCat4(onesyst, 0, genq0q3RhcCat4, loadersRHC);
  PredictionInterp predq0q3RhcCat5(onesyst, 0, genq0q3RhcCat5, loadersRHC);
  PredictionInterp predq0q3RhcCat6(onesyst, 0, genq0q3RhcCat6, loadersRHC);
  PredictionInterp predq0q3RhcCat7(onesyst, 0, genq0q3RhcCat7, loadersRHC);
  PredictionInterp predq0q3RhcCat8(onesyst, 0, genq0q3RhcCat8, loadersRHC);

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

  TH2* h2q0CompRhc   = makeHist2D("h2q0CompRhc", "q_{0} comparison for HPgTPC", &predq0CompRhc, pot_nd);
  TH2* h2q3CompRhc   = makeHist2D("h2q3CompRhc", "q_{3} comparison for HPgTPC", &predq3CompRhc, pot_nd);
  TH2* h2q0q3Rhc     = makeHist2D("h2q0q3Rhc", "q_{0}, q_{3} in HPgTPC", &predq0q3Rhc, pot_nd);
  TH2* h2q0q3RecoRhc = makeHist2D("h2q0q3RecoRhc", "q_{0, reco}, q_{3, reco} in HPgTPC", &predq0q3RecoRhc, pot_nd);
  h2q0CompRhc->Write();
  h2q3CompRhc->Write();
  h2q0q3Rhc->Write();
  h2q0q3RecoRhc->Write();

  // FHC
  TH2 *h2q0CompFhcCat1 = makeHist2D("h2q0CompFhcCat1", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(1).c_str()), &predq0CompFhcCat1, pot_nd);
  TH2 *h2q0CompFhcCat2 = makeHist2D("h2q0CompFhcCat2", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(2).c_str()), &predq0CompFhcCat2, pot_nd);
  TH2 *h2q0CompFhcCat3 = makeHist2D("h2q0CompFhcCat3", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(3).c_str()), &predq0CompFhcCat3, pot_nd);
  TH2 *h2q0CompFhcCat4 = makeHist2D("h2q0CompFhcCat4", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(4).c_str()), &predq0CompFhcCat4, pot_nd);
  TH2 *h2q0CompFhcCat5 = makeHist2D("h2q0CompFhcCat5", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(5).c_str()), &predq0CompFhcCat5, pot_nd);
  TH2 *h2q0CompFhcCat6 = makeHist2D("h2q0CompFhcCat6", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(6).c_str()), &predq0CompFhcCat6, pot_nd);
  TH2 *h2q0CompFhcCat7 = makeHist2D("h2q0CompFhcCat7", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(7).c_str()), &predq0CompFhcCat7, pot_nd);
  TH2 *h2q0CompFhcCat8 = makeHist2D("h2q0CompFhcCat8", Form("q_{0} comparison for HPgTPC, %s (FHC)", catName(8).c_str()), &predq0CompFhcCat8, pot_nd);

  TH2 *h2q3CompFhcCat1 = makeHist2D("h2q3CompFhcCat1", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(1).c_str()), &predq3CompFhcCat1, pot_nd);
  TH2 *h2q3CompFhcCat2 = makeHist2D("h2q3CompFhcCat2", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(2).c_str()), &predq3CompFhcCat2, pot_nd);
  TH2 *h2q3CompFhcCat3 = makeHist2D("h2q3CompFhcCat3", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(3).c_str()), &predq3CompFhcCat3, pot_nd);
  TH2 *h2q3CompFhcCat4 = makeHist2D("h2q3CompFhcCat4", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(4).c_str()), &predq3CompFhcCat4, pot_nd);
  TH2 *h2q3CompFhcCat5 = makeHist2D("h2q3CompFhcCat5", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(5).c_str()), &predq3CompFhcCat5, pot_nd);
  TH2 *h2q3CompFhcCat6 = makeHist2D("h2q3CompFhcCat6", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(6).c_str()), &predq3CompFhcCat6, pot_nd);
  TH2 *h2q3CompFhcCat7 = makeHist2D("h2q3CompFhcCat7", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(7).c_str()), &predq3CompFhcCat7, pot_nd);
  TH2 *h2q3CompFhcCat8 = makeHist2D("h2q3CompFhcCat8", Form("q_{3} comparison for HPgTPC, %s (FHC)", catName(8).c_str()), &predq3CompFhcCat8, pot_nd);
  // RHC
  TH2 *h2q0CompRhcCat1 = makeHist2D("h2q0CompRhcCat1", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(1).c_str()), &predq0CompRhcCat1, pot_nd);
  TH2 *h2q0CompRhcCat2 = makeHist2D("h2q0CompRhcCat2", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(2).c_str()), &predq0CompRhcCat2, pot_nd);
  TH2 *h2q0CompRhcCat3 = makeHist2D("h2q0CompRhcCat3", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(3).c_str()), &predq0CompRhcCat3, pot_nd);
  TH2 *h2q0CompRhcCat4 = makeHist2D("h2q0CompRhcCat4", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(4).c_str()), &predq0CompRhcCat4, pot_nd);
  TH2 *h2q0CompRhcCat5 = makeHist2D("h2q0CompRhcCat5", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(5).c_str()), &predq0CompRhcCat5, pot_nd);
  TH2 *h2q0CompRhcCat6 = makeHist2D("h2q0CompRhcCat6", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(6).c_str()), &predq0CompRhcCat6, pot_nd);
  TH2 *h2q0CompRhcCat7 = makeHist2D("h2q0CompRhcCat7", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(7).c_str()), &predq0CompRhcCat7, pot_nd);
  TH2 *h2q0CompRhcCat8 = makeHist2D("h2q0CompRhcCat8", Form("q_{0} comparison for HPgTPC, %s (RHC)", catName(8).c_str()), &predq0CompRhcCat8, pot_nd);

  TH2 *h2q3CompRhcCat1 = makeHist2D("h2q3CompRhcCat1", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(1).c_str()), &predq3CompRhcCat1, pot_nd);
  TH2 *h2q3CompRhcCat2 = makeHist2D("h2q3CompRhcCat2", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(2).c_str()), &predq3CompRhcCat2, pot_nd);
  TH2 *h2q3CompRhcCat3 = makeHist2D("h2q3CompRhcCat3", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(3).c_str()), &predq3CompRhcCat3, pot_nd);
  TH2 *h2q3CompRhcCat4 = makeHist2D("h2q3CompRhcCat4", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(4).c_str()), &predq3CompRhcCat4, pot_nd);
  TH2 *h2q3CompRhcCat5 = makeHist2D("h2q3CompRhcCat5", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(5).c_str()), &predq3CompRhcCat5, pot_nd);
  TH2 *h2q3CompRhcCat6 = makeHist2D("h2q3CompRhcCat6", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(6).c_str()), &predq3CompRhcCat6, pot_nd);
  TH2 *h2q3CompRhcCat7 = makeHist2D("h2q3CompRhcCat7", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(7).c_str()), &predq3CompRhcCat7, pot_nd);
  TH2 *h2q3CompRhcCat8 = makeHist2D("h2q3CompRhcCat8", Form("q_{3} comparison for HPgTPC, %s (RHC)", catName(8).c_str()), &predq3CompRhcCat8, pot_nd);

  TH2 *hq0q3RecoFhc = ratioSuppressLowStats(&predq0q3RecoFhc, "hq0q3RecoFhc", "FHC", pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoRhc = ratioSuppressLowStats(&predq0q3RecoRhc, "hq0q3RecoRhc", "RHC", pot_nd, rhcPOT, fakedata);
  hq0q3RecoFhc->Write();
  hq0q3RecoRhc->Write();

  TH2 *hq0q3RecoFhcCat1 = ratioSuppressLowStats(&predq0q3RecoFhcCat1, "hq0q3RecoFhcCat1", Form("FHC, reco %s", catName(1).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat2 = ratioSuppressLowStats(&predq0q3RecoFhcCat2, "hq0q3RecoFhcCat2", Form("FHC, reco %s", catName(2).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat3 = ratioSuppressLowStats(&predq0q3RecoFhcCat3, "hq0q3RecoFhcCat3", Form("FHC, reco %s", catName(3).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat4 = ratioSuppressLowStats(&predq0q3RecoFhcCat4, "hq0q3RecoFhcCat4", Form("FHC, reco %s", catName(4).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat5 = ratioSuppressLowStats(&predq0q3RecoFhcCat5, "hq0q3RecoFhcCat5", Form("FHC, reco %s", catName(5).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat6 = ratioSuppressLowStats(&predq0q3RecoFhcCat6, "hq0q3RecoFhcCat6", Form("FHC, reco %s", catName(6).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat7 = ratioSuppressLowStats(&predq0q3RecoFhcCat7, "hq0q3RecoFhcCat7", Form("FHC, reco %s", catName(7).c_str()) , pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3RecoFhcCat8 = ratioSuppressLowStats(&predq0q3RecoFhcCat8, "hq0q3RecoFhcCat8", Form("FHC, reco %s", catName(8).c_str()) , pot_nd, fhcPOT, fakedata);

  TH2 *hq0q3RecoRhcCat1 = ratioSuppressLowStats(&predq0q3RecoRhcCat1, "hq0q3RecoRhcCat1", Form("RHC, reco %s", catName(1).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat2 = ratioSuppressLowStats(&predq0q3RecoRhcCat2, "hq0q3RecoRhcCat2", Form("RHC, reco %s", catName(2).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat3 = ratioSuppressLowStats(&predq0q3RecoRhcCat3, "hq0q3RecoRhcCat3", Form("RHC, reco %s", catName(3).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat4 = ratioSuppressLowStats(&predq0q3RecoRhcCat4, "hq0q3RecoRhcCat4", Form("RHC, reco %s", catName(4).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat5 = ratioSuppressLowStats(&predq0q3RecoRhcCat5, "hq0q3RecoRhcCat5", Form("RHC, reco %s", catName(5).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat6 = ratioSuppressLowStats(&predq0q3RecoRhcCat6, "hq0q3RecoRhcCat6", Form("RHC, reco %s", catName(6).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat7 = ratioSuppressLowStats(&predq0q3RecoRhcCat7, "hq0q3RecoRhcCat7", Form("RHC, reco %s", catName(7).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RecoRhcCat8 = ratioSuppressLowStats(&predq0q3RecoRhcCat8, "hq0q3RecoRhcCat8", Form("RHC, reco %s", catName(8).c_str()) , pot_nd, rhcPOT, fakedata);

  TH2 *hq0q3Fhc     = ratioSuppressLowStats(&predq0q3Fhc, "hq0q3Fhc", "FHC, reco CC inc.; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat1 = ratioSuppressLowStats(&predq0q3FhcCat1, "hq0q3FhcCat1", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(1).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat2 = ratioSuppressLowStats(&predq0q3FhcCat2, "hq0q3FhcCat2", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(2).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat3 = ratioSuppressLowStats(&predq0q3FhcCat3, "hq0q3FhcCat3", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(3).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat4 = ratioSuppressLowStats(&predq0q3FhcCat4, "hq0q3FhcCat4", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(4).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat5 = ratioSuppressLowStats(&predq0q3FhcCat5, "hq0q3FhcCat5", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(5).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat6 = ratioSuppressLowStats(&predq0q3FhcCat6, "hq0q3FhcCat6", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(6).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat7 = ratioSuppressLowStats(&predq0q3FhcCat7, "hq0q3FhcCat7", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(7).c_str()), pot_nd, fhcPOT, fakedata);
  TH2 *hq0q3FhcCat8 = ratioSuppressLowStats(&predq0q3FhcCat8, "hq0q3FhcCat8", Form("FHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(8).c_str()), pot_nd, fhcPOT, fakedata);

  TH2 *hq0q3Rhc = ratioSuppressLowStats(&predq0q3Rhc, "hq0q3Rhc", "RHC, reco CC inc.; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat1 = ratioSuppressLowStats(&predq0q3RhcCat1, "hq0q3RhcCat1", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(1).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat2 = ratioSuppressLowStats(&predq0q3RhcCat2, "hq0q3RhcCat2", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(2).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat3 = ratioSuppressLowStats(&predq0q3RhcCat3, "hq0q3RhcCat3", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(3).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat4 = ratioSuppressLowStats(&predq0q3RhcCat4, "hq0q3RhcCat4", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(4).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat5 = ratioSuppressLowStats(&predq0q3RhcCat5, "hq0q3RhcCat5", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(5).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat6 = ratioSuppressLowStats(&predq0q3RhcCat6, "hq0q3RhcCat6", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(6).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat7 = ratioSuppressLowStats(&predq0q3RhcCat7, "hq0q3RhcCat7", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(7).c_str()) , pot_nd, rhcPOT, fakedata);
  TH2 *hq0q3RhcCat8 = ratioSuppressLowStats(&predq0q3RhcCat8, "hq0q3RhcCat8", Form("RHC, reco %s; p_{vis} / (GeV/c); E_{vis} / GeV; NuWro/GENIE", catName(8).c_str()) , pot_nd, rhcPOT, fakedata);

  // Make most of these again
  // FHC
  TH2 *q0q3FhcCat1 = makeHist2D("hq0q3FhcCat1", Form("FHC, true %s", catName(1).c_str()), &predq0q3FhcCat1, pot_nd);
  TH2 *q0q3FhcCat2 = makeHist2D("hq0q3FhcCat2", Form("FHC, true %s", catName(2).c_str()), &predq0q3FhcCat2, pot_nd);
  TH2 *q0q3FhcCat3 = makeHist2D("hq0q3FhcCat3", Form("FHC, true %s", catName(3).c_str()), &predq0q3FhcCat3, pot_nd);
  TH2 *q0q3FhcCat4 = makeHist2D("hq0q3FhcCat4", Form("FHC, true %s", catName(4).c_str()), &predq0q3FhcCat4, pot_nd);
  TH2 *q0q3FhcCat5 = makeHist2D("hq0q3FhcCat5", Form("FHC, true %s", catName(5).c_str()), &predq0q3FhcCat5, pot_nd);
  TH2 *q0q3FhcCat6 = makeHist2D("hq0q3FhcCat6", Form("FHC, true %s", catName(6).c_str()), &predq0q3FhcCat6, pot_nd);
  TH2 *q0q3FhcCat7 = makeHist2D("hq0q3FhcCat7", Form("FHC, true %s", catName(7).c_str()), &predq0q3FhcCat7, pot_nd);
  TH2 *q0q3FhcCat8 = makeHist2D("hq0q3FhcCat8", Form("FHC, true %s", catName(8).c_str()), &predq0q3FhcCat8, pot_nd);

  q0q3FhcCat1->Write();
  q0q3FhcCat2->Write();
  q0q3FhcCat3->Write();
  q0q3FhcCat4->Write();
  q0q3FhcCat5->Write();
  q0q3FhcCat6->Write();
  q0q3FhcCat7->Write();
  q0q3FhcCat8->Write();

  TH2 *q0q3RecoFhcCat1 = makeHist2D("hq0q3RecoFhcCat1", Form("FHC, reco %s", catName(1).c_str()), &predq0q3RecoFhcCat1, pot_nd);
  TH2 *q0q3RecoFhcCat2 = makeHist2D("hq0q3RecoFhcCat2", Form("FHC, reco %s", catName(2).c_str()), &predq0q3RecoFhcCat2, pot_nd);
  TH2 *q0q3RecoFhcCat3 = makeHist2D("hq0q3RecoFhcCat3", Form("FHC, reco %s", catName(3).c_str()), &predq0q3RecoFhcCat3, pot_nd);
  TH2 *q0q3RecoFhcCat4 = makeHist2D("hq0q3RecoFhcCat4", Form("FHC, reco %s", catName(4).c_str()), &predq0q3RecoFhcCat4, pot_nd);
  TH2 *q0q3RecoFhcCat5 = makeHist2D("hq0q3RecoFhcCat5", Form("FHC, reco %s", catName(5).c_str()), &predq0q3RecoFhcCat5, pot_nd);
  TH2 *q0q3RecoFhcCat6 = makeHist2D("hq0q3RecoFhcCat6", Form("FHC, reco %s", catName(6).c_str()), &predq0q3RecoFhcCat6, pot_nd);
  TH2 *q0q3RecoFhcCat7 = makeHist2D("hq0q3RecoFhcCat7", Form("FHC, reco %s", catName(7).c_str()), &predq0q3RecoFhcCat7, pot_nd);
  TH2 *q0q3RecoFhcCat8 = makeHist2D("hq0q3RecoFhcCat8", Form("FHC, reco %s", catName(8).c_str()), &predq0q3RecoFhcCat8, pot_nd);
  q0q3RecoFhcCat1->Write();
  q0q3RecoFhcCat2->Write();
  q0q3RecoFhcCat3->Write();
  q0q3RecoFhcCat4->Write();
  q0q3RecoFhcCat5->Write();
  q0q3RecoFhcCat6->Write();
  q0q3RecoFhcCat7->Write();
  q0q3RecoFhcCat8->Write();

  // RHC
  TH2 *q0q3RhcCat1 = makeHist2D("hq0q3RhcCat1", Form("RHC, true %s", catName(1).c_str()), &predq0q3RhcCat1, pot_nd);
  TH2 *q0q3RhcCat2 = makeHist2D("hq0q3RhcCat2", Form("RHC, true %s", catName(2).c_str()), &predq0q3RhcCat2, pot_nd);
  TH2 *q0q3RhcCat3 = makeHist2D("hq0q3RhcCat3", Form("RHC, true %s", catName(3).c_str()), &predq0q3RhcCat3, pot_nd);
  TH2 *q0q3RhcCat4 = makeHist2D("hq0q3RhcCat4", Form("RHC, true %s", catName(4).c_str()), &predq0q3RhcCat4, pot_nd);
  TH2 *q0q3RhcCat5 = makeHist2D("hq0q3RhcCat5", Form("RHC, true %s", catName(5).c_str()), &predq0q3RhcCat5, pot_nd);
  TH2 *q0q3RhcCat6 = makeHist2D("hq0q3RhcCat6", Form("RHC, true %s", catName(6).c_str()), &predq0q3RhcCat6, pot_nd);
  TH2 *q0q3RhcCat7 = makeHist2D("hq0q3RhcCat7", Form("RHC, true %s", catName(7).c_str()), &predq0q3RhcCat7, pot_nd);
  TH2 *q0q3RhcCat8 = makeHist2D("hq0q3RhcCat8", Form("RHC, true %s", catName(8).c_str()), &predq0q3RhcCat8, pot_nd);
  q0q3RhcCat1->Write();
  q0q3RhcCat2->Write();
  q0q3RhcCat3->Write();
  q0q3RhcCat4->Write();
  q0q3RhcCat5->Write();
  q0q3RhcCat6->Write();
  q0q3RhcCat7->Write();
  q0q3RhcCat8->Write();

  TH2 *q0q3RecoRhcCat1 = makeHist2D("hq0q3RecoRhcCat1", Form("RHC, reco %s", catName(1).c_str()), &predq0q3RecoRhcCat1, pot_nd);
  TH2 *q0q3RecoRhcCat2 = makeHist2D("hq0q3RecoRhcCat2", Form("RHC, reco %s", catName(2).c_str()), &predq0q3RecoRhcCat2, pot_nd);
  TH2 *q0q3RecoRhcCat3 = makeHist2D("hq0q3RecoRhcCat3", Form("RHC, reco %s", catName(3).c_str()), &predq0q3RecoRhcCat3, pot_nd);
  TH2 *q0q3RecoRhcCat4 = makeHist2D("hq0q3RecoRhcCat4", Form("RHC, reco %s", catName(4).c_str()), &predq0q3RecoRhcCat4, pot_nd);
  TH2 *q0q3RecoRhcCat5 = makeHist2D("hq0q3RecoRhcCat5", Form("RHC, reco %s", catName(5).c_str()), &predq0q3RecoRhcCat5, pot_nd);
  TH2 *q0q3RecoRhcCat6 = makeHist2D("hq0q3RecoRhcCat6", Form("RHC, reco %s", catName(6).c_str()), &predq0q3RecoRhcCat6, pot_nd);
  TH2 *q0q3RecoRhcCat7 = makeHist2D("hq0q3RecoRhcCat7", Form("RHC, reco %s", catName(7).c_str()), &predq0q3RecoRhcCat7, pot_nd);
  TH2 *q0q3RecoRhcCat8 = makeHist2D("hq0q3RecoRhcCat8", Form("RHC, reco %s", catName(8).c_str()), &predq0q3RecoRhcCat8, pot_nd);
  q0q3RecoRhcCat1->Write();
  q0q3RecoRhcCat2->Write();
  q0q3RecoRhcCat3->Write();
  q0q3RecoRhcCat4->Write();
  q0q3RecoRhcCat5->Write();
  q0q3RecoRhcCat6->Write();
  q0q3RecoRhcCat7->Write();
  q0q3RecoRhcCat8->Write();

  hq0q3RecoFhcCat1->Write();
  hq0q3RecoFhcCat2->Write();
  hq0q3RecoFhcCat3->Write();
  hq0q3RecoFhcCat4->Write();
  hq0q3RecoFhcCat5->Write();
  hq0q3RecoFhcCat6->Write();
  hq0q3RecoFhcCat7->Write();
  hq0q3RecoFhcCat8->Write();

  hq0q3RecoRhcCat1->Write();
  hq0q3RecoRhcCat2->Write();
  hq0q3RecoRhcCat3->Write();
  hq0q3RecoRhcCat4->Write();
  hq0q3RecoRhcCat5->Write();
  hq0q3RecoRhcCat6->Write();
  hq0q3RecoRhcCat7->Write();
  hq0q3RecoRhcCat8->Write();

  hq0q3Fhc->Write();
  hq0q3FhcCat1->Write();
  hq0q3FhcCat2->Write();
  hq0q3FhcCat3->Write();
  hq0q3FhcCat4->Write();
  hq0q3FhcCat5->Write();
  hq0q3FhcCat6->Write();
  hq0q3FhcCat7->Write();
  hq0q3FhcCat8->Write();

  hq0q3Rhc->Write();
  hq0q3RhcCat1->Write();
  hq0q3RhcCat2->Write();
  hq0q3RhcCat3->Write();
  hq0q3RhcCat4->Write();
  hq0q3RhcCat5->Write();
  hq0q3RhcCat6->Write();
  hq0q3RhcCat7->Write();
  hq0q3RhcCat8->Write();

  h2q0CompFhcCat1->Write();
  h2q0CompFhcCat2->Write();
  h2q0CompFhcCat3->Write();
  h2q0CompFhcCat4->Write();
  h2q0CompFhcCat5->Write();
  h2q0CompFhcCat6->Write();
  h2q0CompFhcCat7->Write();
  h2q0CompFhcCat8->Write();

  h2q3CompFhcCat1->Write();
  h2q3CompFhcCat2->Write();
  h2q3CompFhcCat3->Write();
  h2q3CompFhcCat4->Write();
  h2q3CompFhcCat5->Write();
  h2q3CompFhcCat6->Write();
  h2q3CompFhcCat7->Write();
  h2q3CompFhcCat8->Write();

  h2q0CompRhcCat1->Write();
  h2q0CompRhcCat2->Write();
  h2q0CompRhcCat3->Write();
  h2q0CompRhcCat4->Write();
  h2q0CompRhcCat5->Write();
  h2q0CompRhcCat6->Write();
  h2q0CompRhcCat7->Write();
  h2q0CompRhcCat8->Write();

  h2q3CompRhcCat1->Write();
  h2q3CompRhcCat2->Write();
  h2q3CompRhcCat3->Write();
  h2q3CompRhcCat4->Write();
  h2q3CompRhcCat5->Write();
  h2q3CompRhcCat6->Write();
  h2q3CompRhcCat7->Write();
  h2q3CompRhcCat8->Write();

  fout->Close();
  delete fout;
} // q0q3Comp
