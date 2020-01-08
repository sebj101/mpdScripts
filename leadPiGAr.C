// leadpiGAr.C
// Plot energy of leading and sub-leading pion for NuWro and nominal
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
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

/// Modes list:
/// * QE: 1
/// * Single Kaon: 2
/// * DIS: 3
/// * RES: 4
/// * COH: 5
/// * Diffractive: 6
/// * Nu-e El: 7
/// * IMD: 8
/// * AMnuGamma: 9
/// * MEC: 10
/// * COHEl: 11
/// * IBD: 12
/// * GlashowRES: 13
/// * IMDAnnihalation: 14

using namespace ana;

// POT for 3.5 years
const double pot_fd = 3.5 * POT120 * 40/1.13;
const double pot_nd = 3.5 * POT120;
// Vars
const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var KPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kEnu = SIMPLEVAR(dune.Ev);
const Var kRecoEnergyND  = SIMPLEVAR(dune.Ev_reco);
const Var kLeadPiE = SIMPLEVAR(dune.gastpc_lead_pi_E);
const Var kSubleadPiE = SIMPLEVAR(dune.gastpc_sublead_pi_E);
const Var kMode = SIMPLEVAR(dune.mode);
const Var kQ2 = SIMPLEVAR(dune.Q2);
const Var kFHC = SIMPLEVAR(dune.isFHC);
// Reco Q2
const Var kRecoQ2({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		  [](const caf::StandardRecord* sr) {
		    double q2 = 0;
		    q2 = 4 * sr->dune.Elep_reco *TMath::Sin(sr->dune.theta_reco/2.) * TMath::Sin(sr->dune.theta_reco/2.) *sr->dune.Ev_reco;
		    return q2;
		  });
// Reco W
const Var kRecoW({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		 [](const caf::StandardRecord* sr) {
		   double w = 0;
		   double q2 = 4 * sr->dune.Elep_reco *TMath::Sin(sr->dune.theta_reco/2.) * TMath::Sin(sr->dune.theta_reco/2.) *sr->dune.Ev_reco;
		   w = TMath::Sqrt(-q2 + 2 * 0.939 * (sr->dune.Ev_reco-sr->dune.Elep_reco) + 0.939*0.939);
		   return w;
		 });

std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
				 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
const Binning simpleBins = Binning::Simple(30, 0, 3);
const Binning simpleQ2Bins = Binning::Simple(30, 0, 6);
const Binning simpleWBins  = Binning::Simple(30, 0, 6);
const Binning binsEreco  = Binning::Custom(binEEdges);

const HistAxis axRecoW("W (GeV)", simpleWBins, kRecoW);
const HistAxis axND("E_{#nu, reco} / GeV", binsEreco, kRecoEnergyND);
const HistAxis axNDTrue("E_{#nu} / GeV", binsEreco, kEnu);
const HistAxis axQ2True("Q^{2} / GeV", simpleQ2Bins, kQ2);
const HistAxis axQ2Reco("Q^{2}_{reco} / GeV", simpleQ2Bins, kRecoQ2);
const HistAxis axLeadPi("E_{#pi, lead} (GeV)", simpleBins, kLeadPiE);
const HistAxis axSubleadPi("E_{#pi, sub-lead} (GeV)", simpleBins, kSubleadPiE);

void leadPiGAr(const char *outfile, const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v3/") 
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

  // Enu reco spectra
  NoOscPredictionGenerator genGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genGArNumuFHC0(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArNumuFHC1(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArNumuFHC2(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArNumuFHC3(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArNumuFHCHi(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
  PredictionInterp predGArNumuFHC(systlist, 0, genGArNumuFHC, loadersGArFHC);
  PredictionInterp predGArNumuFHC0(systlist, 0, genGArNumuFHC0, loadersGArFHC);
  PredictionInterp predGArNumuFHC1(systlist, 0, genGArNumuFHC1, loadersGArFHC);
  PredictionInterp predGArNumuFHC2(systlist, 0, genGArNumuFHC2, loadersGArFHC);
  PredictionInterp predGArNumuFHC3(systlist, 0, genGArNumuFHC3, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi(systlist, 0, genGArNumuFHCHi, loadersGArFHC);
  // QE
  NoOscPredictionGenerator genGArNumuFHC_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==1);
  NoOscPredictionGenerator genGArNumuFHC0_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0 && kMode==1);
  NoOscPredictionGenerator genGArNumuFHC1_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1 && kMode==1);
  NoOscPredictionGenerator genGArNumuFHC2_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2 && kMode==1);
  NoOscPredictionGenerator genGArNumuFHC3_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3 && kMode==1);
  NoOscPredictionGenerator genGArNumuFHCHi_qe(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3 && kMode==1);
  PredictionInterp predGArNumuFHC_qe(systlist, 0, genGArNumuFHC_qe, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_qe(systlist, 0, genGArNumuFHC0_qe, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_qe(systlist, 0, genGArNumuFHC1_qe, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_qe(systlist, 0, genGArNumuFHC2_qe, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_qe(systlist, 0, genGArNumuFHC3_qe, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_qe(systlist, 0, genGArNumuFHCHi_qe, loadersGArFHC);
  // DIS
  NoOscPredictionGenerator genGArNumuFHC_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==3);
  NoOscPredictionGenerator genGArNumuFHC0_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0 && kMode==3);
  NoOscPredictionGenerator genGArNumuFHC1_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1 && kMode==3);
  NoOscPredictionGenerator genGArNumuFHC2_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2 && kMode==3);
  NoOscPredictionGenerator genGArNumuFHC3_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3 && kMode==3);
  NoOscPredictionGenerator genGArNumuFHCHi_dis(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3 && kMode==3);
  PredictionInterp predGArNumuFHC_dis(systlist, 0, genGArNumuFHC_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_dis(systlist, 0, genGArNumuFHC0_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_dis(systlist, 0, genGArNumuFHC1_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_dis(systlist, 0, genGArNumuFHC2_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_dis(systlist, 0, genGArNumuFHC3_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_dis(systlist, 0, genGArNumuFHCHi_dis, loadersGArFHC);
  // RES
  NoOscPredictionGenerator genGArNumuFHC_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==4);
  NoOscPredictionGenerator genGArNumuFHC0_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0 && kMode==4);
  NoOscPredictionGenerator genGArNumuFHC1_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1 && kMode==4);
  NoOscPredictionGenerator genGArNumuFHC2_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2 && kMode==4);
  NoOscPredictionGenerator genGArNumuFHC3_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3 && kMode==4);
  NoOscPredictionGenerator genGArNumuFHCHi_res(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3 && kMode==4);
  PredictionInterp predGArNumuFHC_res(systlist, 0, genGArNumuFHC_res, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_res(systlist, 0, genGArNumuFHC0_res, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_res(systlist, 0, genGArNumuFHC1_res, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_res(systlist, 0, genGArNumuFHC2_res, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_res(systlist, 0, genGArNumuFHC3_res, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_res(systlist, 0, genGArNumuFHCHi_res, loadersGArFHC);
  // COH
  NoOscPredictionGenerator genGArNumuFHC_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==5);
  NoOscPredictionGenerator genGArNumuFHC0_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0 && kMode==5);
  NoOscPredictionGenerator genGArNumuFHC1_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1 && kMode==5);
  NoOscPredictionGenerator genGArNumuFHC2_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2 && kMode==5);
  NoOscPredictionGenerator genGArNumuFHC3_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3 && kMode==5);
  NoOscPredictionGenerator genGArNumuFHCHi_coh(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3 && kMode==5);
  PredictionInterp predGArNumuFHC_coh(systlist, 0, genGArNumuFHC_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_coh(systlist, 0, genGArNumuFHC0_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_coh(systlist, 0, genGArNumuFHC1_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_coh(systlist, 0, genGArNumuFHC2_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_coh(systlist, 0, genGArNumuFHC3_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_coh(systlist, 0, genGArNumuFHCHi_coh, loadersGArFHC);
  // MEC
  NoOscPredictionGenerator genGArNumuFHC_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==10);
  NoOscPredictionGenerator genGArNumuFHC0_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0 && kMode==10);
  NoOscPredictionGenerator genGArNumuFHC1_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1 && kMode==10);
  NoOscPredictionGenerator genGArNumuFHC2_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2 && kMode==10);
  NoOscPredictionGenerator genGArNumuFHC3_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3 && kMode==10);
  NoOscPredictionGenerator genGArNumuFHCHi_mec(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3 && kMode==10);
  PredictionInterp predGArNumuFHC_mec(systlist, 0, genGArNumuFHC_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_mec(systlist, 0, genGArNumuFHC0_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_mec(systlist, 0, genGArNumuFHC1_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_mec(systlist, 0, genGArNumuFHC2_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_mec(systlist, 0, genGArNumuFHC3_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_mec(systlist, 0, genGArNumuFHCHi_mec, loadersGArFHC);

  // True Q2 separated by interaction mode
  NoOscPredictionGenerator genGArNumuFHCQ2_mec(axQ2True, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==10);
  NoOscPredictionGenerator genGArNumuFHCQ2_coh(axQ2True, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==5);
  NoOscPredictionGenerator genGArNumuFHCQ2_res(axQ2True, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==4);
  NoOscPredictionGenerator genGArNumuFHCQ2_dis(axQ2True, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==3);
  NoOscPredictionGenerator genGArNumuFHCQ2_qe(axQ2True, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==1);
  PredictionInterp predGArNumuFHCQ2_mec(systlist, 0, genGArNumuFHCQ2_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2_coh(systlist, 0, genGArNumuFHCQ2_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2_res(systlist, 0, genGArNumuFHCQ2_res, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2_dis(systlist, 0, genGArNumuFHCQ2_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2_qe(systlist, 0, genGArNumuFHCQ2_qe, loadersGArFHC);

  // Reco Q2 separated by interaction mode
  NoOscPredictionGenerator genGArNumuFHCQ2Reco_mec(axQ2Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==10);
  NoOscPredictionGenerator genGArNumuFHCQ2Reco_coh(axQ2Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==5);
  NoOscPredictionGenerator genGArNumuFHCQ2Reco_res(axQ2Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==4);
  NoOscPredictionGenerator genGArNumuFHCQ2Reco_dis(axQ2Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==3);
  NoOscPredictionGenerator genGArNumuFHCQ2Reco_qe(axQ2Reco, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==1);
  PredictionInterp predGArNumuFHCQ2Reco_mec(systlist, 0, genGArNumuFHCQ2Reco_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2Reco_coh(systlist, 0, genGArNumuFHCQ2Reco_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2Reco_res(systlist, 0, genGArNumuFHCQ2Reco_res, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2Reco_dis(systlist, 0, genGArNumuFHCQ2Reco_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHCQ2Reco_qe(systlist, 0, genGArNumuFHCQ2Reco_qe, loadersGArFHC);

  // True energy separated by interaction mode
  NoOscPredictionGenerator genGArNumuFHCEnu_mec(axNDTrue, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==10);
  NoOscPredictionGenerator genGArNumuFHCEnu_coh(axNDTrue, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==5);
  NoOscPredictionGenerator genGArNumuFHCEnu_res(axNDTrue, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==4);
  NoOscPredictionGenerator genGArNumuFHCEnu_dis(axNDTrue, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==3);
  NoOscPredictionGenerator genGArNumuFHCEnu_qe(axNDTrue, kPassND_FHC_NUMU && kIsTrueGasFV && kMode==1);
  PredictionInterp predGArNumuFHCEnu_mec(systlist, 0, genGArNumuFHCEnu_mec, loadersGArFHC);
  PredictionInterp predGArNumuFHCEnu_coh(systlist, 0, genGArNumuFHCEnu_coh, loadersGArFHC);
  PredictionInterp predGArNumuFHCEnu_res(systlist, 0, genGArNumuFHCEnu_res, loadersGArFHC);
  PredictionInterp predGArNumuFHCEnu_dis(systlist, 0, genGArNumuFHCEnu_dis, loadersGArFHC);
  PredictionInterp predGArNumuFHCEnu_qe(systlist, 0, genGArNumuFHCEnu_qe, loadersGArFHC);

  // Reco W separated by pion multiplicity
  NoOscPredictionGenerator WGAr(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator WGArCC0Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator WGArCC1Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator WGArCC2Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator WGArCC3Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator WGArCCHiPi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);
  PredictionInterp predWGAr(systlist, 0, WGAr, loadersGArFHC);
  PredictionInterp predWGArCC0Pi(systlist, 0, WGArCC0Pi, loadersGArFHC);
  PredictionInterp predWGArCC1Pi(systlist, 0, WGArCC1Pi, loadersGArFHC);
  PredictionInterp predWGArCC2Pi(systlist, 0, WGArCC2Pi, loadersGArFHC);
  PredictionInterp predWGArCC3Pi(systlist, 0, WGArCC3Pi, loadersGArFHC);
  PredictionInterp predWGArCCHiPi(systlist, 0, WGArCCHiPi, loadersGArFHC);
  // Reco Q2 separated by pion multiplicity
  NoOscPredictionGenerator Q2GAr(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator Q2GArCC0Pi(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator Q2GArCC1Pi(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator Q2GArCC2Pi(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator Q2GArCC3Pi(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator Q2GArCCHiPi(axQ2Reco, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);
  PredictionInterp predQ2GAr(systlist, 0, Q2GAr, loadersGArFHC);
  PredictionInterp predQ2GArCC0Pi(systlist, 0, Q2GArCC0Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC1Pi(systlist, 0, Q2GArCC1Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC2Pi(systlist, 0, Q2GArCC2Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC3Pi(systlist, 0, Q2GArCC3Pi, loadersGArFHC);
  PredictionInterp predQ2GArCCHiPi(systlist, 0, Q2GArCCHiPi, loadersGArFHC);

  // Leading pi energies
  NoOscPredictionGenerator genGArNumuFHC_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genGArNumuFHC0_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArNumuFHC1_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArNumuFHC2_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArNumuFHC3_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArNumuFHCHi_lead(axLeadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
  PredictionInterp predGArNumuFHC_lead(systlist, 0, genGArNumuFHC_lead, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_lead(systlist, 0, genGArNumuFHC0_lead, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_lead(systlist, 0, genGArNumuFHC1_lead, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_lead(systlist, 0, genGArNumuFHC2_lead, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_lead(systlist, 0, genGArNumuFHC3_lead, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_lead(systlist, 0, genGArNumuFHCHi_lead, loadersGArFHC);

  NoOscPredictionGenerator genGArNumuFHC_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genGArNumuFHC0_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArNumuFHC1_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArNumuFHC2_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArNumuFHC3_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArNumuFHCHi_sublead(axSubleadPi, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);

  PredictionInterp predGArNumuFHC_sublead(systlist, 0, genGArNumuFHC_sublead, loadersGArFHC);
  PredictionInterp predGArNumuFHC0_sublead(systlist, 0, genGArNumuFHC0_sublead, loadersGArFHC);
  PredictionInterp predGArNumuFHC1_sublead(systlist, 0, genGArNumuFHC1_sublead, loadersGArFHC);
  PredictionInterp predGArNumuFHC2_sublead(systlist, 0, genGArNumuFHC2_sublead, loadersGArFHC);
  PredictionInterp predGArNumuFHC3_sublead(systlist, 0, genGArNumuFHC3_sublead, loadersGArFHC);
  PredictionInterp predGArNumuFHCHi_sublead(systlist, 0, genGArNumuFHCHi_sublead, loadersGArFHC);
  // RHC
  NoOscPredictionGenerator genGArNumuRHC_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genGArNumuRHC0_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArNumuRHC1_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArNumuRHC2_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArNumuRHC3_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArNumuRHCHi_lead(axLeadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);
  PredictionInterp predGArNumuRHC_lead(systlist, 0, genGArNumuRHC_lead, loadersGArRHC);
  PredictionInterp predGArNumuRHC0_lead(systlist, 0, genGArNumuRHC0_lead, loadersGArRHC);
  PredictionInterp predGArNumuRHC1_lead(systlist, 0, genGArNumuRHC1_lead, loadersGArRHC);
  PredictionInterp predGArNumuRHC2_lead(systlist, 0, genGArNumuRHC2_lead, loadersGArRHC);
  PredictionInterp predGArNumuRHC3_lead(systlist, 0, genGArNumuRHC3_lead, loadersGArRHC);
  PredictionInterp predGArNumuRHCHi_lead(systlist, 0, genGArNumuRHCHi_lead, loadersGArRHC);

  NoOscPredictionGenerator genGArNumuRHC_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV);
  NoOscPredictionGenerator genGArNumuRHC0_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArNumuRHC1_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArNumuRHC2_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArNumuRHC3_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArNumuRHCHi_sublead(axSubleadPi, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);

  PredictionInterp predGArNumuRHC_sublead(systlist, 0, genGArNumuRHC_sublead, loadersGArRHC);
  PredictionInterp predGArNumuRHC0_sublead(systlist, 0, genGArNumuRHC0_sublead, loadersGArRHC);
  PredictionInterp predGArNumuRHC1_sublead(systlist, 0, genGArNumuRHC1_sublead, loadersGArRHC);
  PredictionInterp predGArNumuRHC2_sublead(systlist, 0, genGArNumuRHC2_sublead, loadersGArRHC);
  PredictionInterp predGArNumuRHC3_sublead(systlist, 0, genGArNumuRHC3_sublead, loadersGArRHC);
  PredictionInterp predGArNumuRHCHi_sublead(systlist, 0, genGArNumuRHCHi_sublead, loadersGArRHC);

  loadersGArFHC.Go();
  loadersGArRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();
  // Nominal samples
  TH1 *hnumuFHC_lead   = predGArNumuFHC_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC0_lead  = predGArNumuFHC0_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC1_lead  = predGArNumuFHC1_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC2_lead  = predGArNumuFHC2_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC3_lead  = predGArNumuFHC3_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCHi_lead = predGArNumuFHCHi_lead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hnumuFHC_lead);
  setHistAttr(hnumuFHC0_lead);
  setHistAttr(hnumuFHC1_lead);
  setHistAttr(hnumuFHC2_lead);
  setHistAttr(hnumuFHC3_lead);
  setHistAttr(hnumuFHCHi_lead);
  hnumuFHC_lead->SetTitle("Lead #pi energy");
  hnumuFHC0_lead->SetTitle("Lead #pi energy, CC0#pi");
  hnumuFHC1_lead->SetTitle("Lead #pi energy, CC1#pi");
  hnumuFHC2_lead->SetTitle("Lead #pi energy, CC2#pi");
  hnumuFHC3_lead->SetTitle("Lead #pi energy, CC3#pi");
  hnumuFHCHi_lead->SetTitle("Lead #pi energy, CC>3#pi");
  hnumuFHC_lead->Write("hnumuFHC_lead");
  hnumuFHC0_lead->Write("hnumuFHC0_lead");
  hnumuFHC1_lead->Write("hnumuFHC1_lead");
  hnumuFHC2_lead->Write("hnumuFHC2_lead");
  hnumuFHC3_lead->Write("hnumuFHC3_lead");
  hnumuFHCHi_lead->Write("hnumuFHCHi_lead");

  TH1 *hnumuFHC_sublead   = predGArNumuFHC_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC0_sublead  = predGArNumuFHC0_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC1_sublead  = predGArNumuFHC1_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC2_sublead  = predGArNumuFHC2_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHC3_sublead  = predGArNumuFHC3_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCHi_sublead = predGArNumuFHCHi_sublead.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hnumuFHC_sublead);
  setHistAttr(hnumuFHC0_sublead);
  setHistAttr(hnumuFHC1_sublead);
  setHistAttr(hnumuFHC2_sublead);
  setHistAttr(hnumuFHC3_sublead);
  setHistAttr(hnumuFHCHi_sublead);
  hnumuFHC_sublead->SetTitle("Sub-leading #pi energy");
  hnumuFHC0_sublead->SetTitle("Sub-leading #pi energy, CC0#pi");
  hnumuFHC1_sublead->SetTitle("Sub-leading #pi energy, CC1#pi");
  hnumuFHC2_sublead->SetTitle("Sub-leading #pi energy, CC2#pi");
  hnumuFHC3_sublead->SetTitle("Sub-leading #pi energy, CC3#pi");
  hnumuFHCHi_sublead->SetTitle("Sub-leading #pi energy, CC>3#pi");
  hnumuFHC_sublead->Write("hnumuFHC_sublead");
  hnumuFHC0_sublead->Write("hnumuFHC0_sublead");
  hnumuFHC1_sublead->Write("hnumuFHC1_sublead");
  hnumuFHC2_sublead->Write("hnumuFHC2_sublead");
  hnumuFHC3_sublead->Write("hnumuFHC3_sublead");
  hnumuFHCHi_sublead->Write("hnumuFHCHi_sublead");

  // With NuWro shifts
  TH1 *hnumuFHCnuwro_lead   = predGArNumuFHC_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro0_lead  = predGArNumuFHC0_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro1_lead  = predGArNumuFHC1_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro2_lead  = predGArNumuFHC2_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro3_lead  = predGArNumuFHC3_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwroHi_lead = predGArNumuFHCHi_lead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hnumuFHCnuwro_lead);
  setHistAttr(hnumuFHCnuwro0_lead);
  setHistAttr(hnumuFHCnuwro1_lead);
  setHistAttr(hnumuFHCnuwro2_lead);
  setHistAttr(hnumuFHCnuwro3_lead);
  setHistAttr(hnumuFHCnuwroHi_lead);
  hnumuFHCnuwro_lead->SetTitle("Lead #pi energy");
  hnumuFHCnuwro0_lead->SetTitle("Lead #pi energy, CC0#pi (NuWro)");
  hnumuFHCnuwro1_lead->SetTitle("Lead #pi energy, CC1#pi (NuWro)");
  hnumuFHCnuwro2_lead->SetTitle("Lead #pi energy, CC2#pi (NuWro)");
  hnumuFHCnuwro3_lead->SetTitle("Lead #pi energy, CC3#pi (NuWro)");
  hnumuFHCnuwroHi_lead->SetTitle("Lead #pi energy, CC>3#pi (NuWro)");
  hnumuFHCnuwro_lead->Write("hnumuFHCnuwro_lead");
  hnumuFHCnuwro0_lead->Write("hnumuFHCnuwro0_lead");
  hnumuFHCnuwro1_lead->Write("hnumuFHCnuwro1_lead");
  hnumuFHCnuwro2_lead->Write("hnumuFHCnuwro2_lead");
  hnumuFHCnuwro3_lead->Write("hnumuFHCnuwro3_lead");
  hnumuFHCnuwroHi_lead->Write("hnumuFHCnuwroHi_lead");

  // Now do the ratios
  TH1D *hRatioFHC_lead = new TH1D("hRatioFHC_lead", "Lead #pi energy, GENIE/NuWro; E_{#pi, lead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC1_lead = new TH1D("hRatioFHC1_lead", "Lead #pi energy, CC1#pi, GENIE/NuWro; E_{#pi, lead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC2_lead = new TH1D("hRatioFHC2_lead", "Lead #pi energy, CC2#pi, GENIE/NuWro; E_{#pi, lead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC3_lead = new TH1D("hRatioFHC3_lead", "Lead #pi energy, CC3#pi, GENIE/NuWro; E_{#pi, lead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHCHi_lead = new TH1D("hRatioFHCHi_lead", "Lead #pi energy, CC>3#pi, GENIE/NuWro; E_{#pi, lead}; GENIE/NuWro", 30, 0., 3.);
  hRatioFHC_lead->Divide(hnumuFHC_lead, hnumuFHCnuwro_lead, 1., 1., "B");
  hRatioFHC1_lead->Divide(hnumuFHC1_lead, hnumuFHCnuwro1_lead, 1., 1., "B");
  hRatioFHC2_lead->Divide(hnumuFHC2_lead, hnumuFHCnuwro2_lead, 1., 1., "B");
  hRatioFHC3_lead->Divide(hnumuFHC3_lead, hnumuFHCnuwro3_lead, 1., 1., "B");
  hRatioFHCHi_lead->Divide(hnumuFHCHi_lead, hnumuFHCnuwroHi_lead, 1., 1., "B");
  hRatioFHC_lead->Write();
  hRatioFHC1_lead->Write();
  hRatioFHC2_lead->Write();
  hRatioFHC3_lead->Write();
  hRatioFHCHi_lead->Write();

  TH1 *hnumuFHCnuwro_sublead   = predGArNumuFHC_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro0_sublead  = predGArNumuFHC0_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro1_sublead  = predGArNumuFHC1_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro2_sublead  = predGArNumuFHC2_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwro3_sublead  = predGArNumuFHC3_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumuFHCnuwroHi_sublead = predGArNumuFHCHi_sublead.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hnumuFHCnuwro_sublead);
  setHistAttr(hnumuFHCnuwro0_sublead);
  setHistAttr(hnumuFHCnuwro1_sublead);
  setHistAttr(hnumuFHCnuwro2_sublead);
  setHistAttr(hnumuFHCnuwro3_sublead);
  setHistAttr(hnumuFHCnuwroHi_sublead);
  hnumuFHCnuwro_sublead->SetTitle("Sub-leading #pi energy");
  hnumuFHCnuwro0_sublead->SetTitle("Sub-leading #pi energy, CC0#pi (NuWro)");
  hnumuFHCnuwro1_sublead->SetTitle("Sub-leading #pi energy, CC1#pi (NuWro)");
  hnumuFHCnuwro2_sublead->SetTitle("Sub-leading #pi energy, CC2#pi (NuWro)");
  hnumuFHCnuwro3_sublead->SetTitle("Sub-leading #pi energy, CC3#pi (NuWro)");
  hnumuFHCnuwroHi_sublead->SetTitle("Sub-leading #pi energy, CC>3#pi (NuWro)");
  hnumuFHCnuwro_sublead->Write("hnumuFHCnuwro_sublead");
  hnumuFHCnuwro0_sublead->Write("hnumuFHCnuwro0_sublead");
  hnumuFHCnuwro1_sublead->Write("hnumuFHCnuwro1_sublead");
  hnumuFHCnuwro2_sublead->Write("hnumuFHCnuwro2_sublead");
  hnumuFHCnuwro3_sublead->Write("hnumuFHCnuwro3_sublead");
  hnumuFHCnuwroHi_sublead->Write("hnumuFHCnuwroHi_sublead");

  // Now do the ratios
  TH1D *hRatioFHC_sublead = new TH1D("hRatioFHC_sublead", "Sub-leading #pi energy, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC1_sublead  = new TH1D("hRatioFHC1_sublead", "Sub-leading #pi energy, CC1#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC2_sublead  = new TH1D("hRatioFHC2_sublead", "Sub-leading #pi energy, CC2#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHC3_sublead  = new TH1D("hRatioFHC3_sublead", "Sub-leading #pi energy, CC3#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro", 30, 0., 3.);
  TH1D *hRatioFHCHi_sublead  = new TH1D("hRatioFHCHi_sublead", "Sub-leading #pi energy, CC>3#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro", 30, 0., 3.);
  hRatioFHC_sublead->Divide(hnumuFHC_sublead, hnumuFHCnuwro_sublead, 1., 1., "B");
  hRatioFHC1_sublead->Divide(hnumuFHC1_sublead, hnumuFHCnuwro1_sublead, 1., 1., "B");
  hRatioFHC2_sublead->Divide(hnumuFHC2_sublead, hnumuFHCnuwro2_sublead, 1., 1., "B");
  hRatioFHC3_sublead->Divide(hnumuFHC3_sublead, hnumuFHCnuwro3_sublead, 1., 1., "B");
  hRatioFHCHi_sublead->Divide(hnumuFHCHi_sublead, hnumuFHCnuwroHi_sublead, 1., 1., "B"); 
  hRatioFHC1_sublead->SetTitle("Sub-leading #pi energy, CC1#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro");
  hRatioFHC2_sublead->SetTitle("Sub-leading #pi energy, CC2#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro");
  hRatioFHC3_sublead->SetTitle("Sub-leading #pi energy, CC3#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro");
  hRatioFHCHi_sublead->SetTitle("Sub-leading #pi energy, CC>3#pi, GENIE/NuWro; E_{#pi, sublead}; GENIE/NuWro");
  hRatioFHC_sublead->Write();
  hRatioFHC1_sublead->Write();
  hRatioFHC2_sublead->Write();
  hRatioFHC3_sublead->Write();
  hRatioFHCHi_sublead->Write();

  // Reconstructed kinematic histograms
  TH1 *hWGAr      = predWGAr.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArNuWro = predWGAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC0Pi      = predWGArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC0PiNuWro = predWGArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC1Pi      = predWGArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC1PiNuWro = predWGArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC2Pi      = predWGArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC2PiNuWro = predWGArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC3Pi      = predWGArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC3PiNuWro = predWGArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCCHiPi      = predWGArCCHiPi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCCHiPiNuWro = predWGArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  hWGArCC0Pi->SetTitle("0#pi");
  hWGArCC1Pi->SetTitle("1#pi");
  hWGArCC2Pi->SetTitle("2#pi");
  hWGArCC3Pi->SetTitle("3#pi");
  hWGArCCHiPi->SetTitle(">3#pi");
  hWGArCC0PiNuWro->SetTitle("0#pi");
  hWGArCC1PiNuWro->SetTitle("1#pi");
  hWGArCC2PiNuWro->SetTitle("2#pi");
  hWGArCC3PiNuWro->SetTitle("3#pi");
  hWGArCCHiPiNuWro->SetTitle(">3#pi");

  TH1 *hQ2GAr      = predQ2GAr.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArNuWro = predQ2GAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC0Pi      = predQ2GArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC0PiNuWro = predQ2GArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC1Pi      = predQ2GArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC1PiNuWro = predQ2GArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC2Pi      = predQ2GArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC2PiNuWro = predQ2GArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC3Pi      = predQ2GArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC3PiNuWro = predQ2GArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCCHiPi      = predQ2GArCCHiPi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCCHiPiNuWro = predQ2GArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  hQ2GArCC0Pi->SetTitle("0#pi");
  hQ2GArCC1Pi->SetTitle("1#pi");
  hQ2GArCC2Pi->SetTitle("2#pi");
  hQ2GArCC3Pi->SetTitle("3#pi");
  hQ2GArCCHiPi->SetTitle(">3#pi");
  hQ2GArCC0PiNuWro->SetTitle("0#pi");
  hQ2GArCC1PiNuWro->SetTitle("1#pi");
  hQ2GArCC2PiNuWro->SetTitle("2#pi");
  hQ2GArCC3PiNuWro->SetTitle("3#pi");
  hQ2GArCCHiPiNuWro->SetTitle(">3#pi");

  hWGArCC0Pi->Write("hWGArCC0Pi");
  hWGArCC0PiNuWro->Write("hWGArCC0PiNuWro");
  hWGArCC1Pi->Write("hWGArCC1Pi");
  hWGArCC1PiNuWro->Write("hWGArCC1PiNuWro");
  hWGArCC2Pi->Write("hWGArCC2Pi");
  hWGArCC2PiNuWro->Write("hWGArCC2PiNuWro");
  hWGArCC3Pi->Write("hWGArCC3Pi");
  hWGArCC3PiNuWro->Write("hWGArCC3PiNuWro");
  hWGArCCHiPi->Write("hWGArCCHiPi");
  hWGArCCHiPiNuWro->Write("hWGArCCHiPiNuWro");
  hQ2GArCC0Pi->Write("hQ2GArCC0Pi");
  hQ2GArCC0PiNuWro->Write("hQ2GArCC0PiNuWro");
  hQ2GArCC1Pi->Write("hQ2GArCC1Pi");
  hQ2GArCC1PiNuWro->Write("hQ2GArCC1PiNuWro");
  hQ2GArCC2Pi->Write("hQ2GArCC2Pi");
  hQ2GArCC2PiNuWro->Write("hQ2GArCC2PiNuWro");
  hQ2GArCC3Pi->Write("hQ2GArCC3Pi");
  hQ2GArCC3PiNuWro->Write("hQ2GArCC3PiNuWro");
  hQ2GArCCHiPi->Write("hQ2GArCCHiPi");
  hQ2GArCCHiPiNuWro->Write("hQ2GArCCHiPiNuWro");
  // Ratios
  THStack *hsWGArGenieNuWroRatios = new THStack("hsWGArGenieNuWroRatios", "GENIE/NuWro for W in GAr; W_{reco} / GeV; GENIE/NuWro");
  hWGArCC0Pi->Divide(hWGArCC0PiNuWro);
  hWGArCC1Pi->Divide(hWGArCC1PiNuWro);
  hWGArCC2Pi->Divide(hWGArCC2PiNuWro);
  hWGArCC3Pi->Divide(hWGArCC3PiNuWro);
  hWGArCCHiPi->Divide(hWGArCCHiPiNuWro);
  setHistAttr(hWGArCC0Pi);
  setHistAttr(hWGArCC1Pi);
  setHistAttr(hWGArCC2Pi);
  setHistAttr(hWGArCC3Pi);
  setHistAttr(hWGArCCHiPi);
  hWGArCC0Pi->SetLineColor(kBlack);
  hWGArCC1Pi->SetLineColor(kBlue);
  hWGArCC2Pi->SetLineColor(kRed);
  hWGArCC3Pi->SetLineColor(kGreen+2);
  hWGArCCHiPi->SetLineColor(kMagenta);
  hsWGArGenieNuWroRatios->Add(hWGArCC0Pi);
  hsWGArGenieNuWroRatios->Add(hWGArCC1Pi);
  hsWGArGenieNuWroRatios->Add(hWGArCC2Pi);
  hsWGArGenieNuWroRatios->Add(hWGArCC3Pi);
  hsWGArGenieNuWroRatios->Add(hWGArCCHiPi);
  hsWGArGenieNuWroRatios->Write();
  hWGArCC0Pi->Write("hWGArCC0PiRatio");
  hWGArCC1Pi->Write("hWGArCC1PiRatio");
  hWGArCC2Pi->Write("hWGArCC2PiRatio");
  hWGArCC3Pi->Write("hWGArCC3PiRatio");
  hWGArCCHiPi->Write("hWGArCCHiPiRatio");

  THStack *hsQ2GArGenieNuWroRatios = new THStack("hsQ2GArGenieNuWroRatios", "GENIE/NuWro for Q2 in GAr; Q^{2}_{reco} / GeV; GENIE/NuWro");
  hQ2GArCC0Pi->Divide(hQ2GArCC0PiNuWro);
  hQ2GArCC1Pi->Divide(hQ2GArCC1PiNuWro);
  hQ2GArCC2Pi->Divide(hQ2GArCC2PiNuWro);
  hQ2GArCC3Pi->Divide(hQ2GArCC3PiNuWro);
  hQ2GArCCHiPi->Divide(hQ2GArCCHiPiNuWro);
  setHistAttr(hQ2GArCC0Pi);
  setHistAttr(hQ2GArCC1Pi);
  setHistAttr(hQ2GArCC2Pi);
  setHistAttr(hQ2GArCC3Pi);
  setHistAttr(hQ2GArCCHiPi);
  hQ2GArCC0Pi->SetLineColor(kBlack);
  hQ2GArCC1Pi->SetLineColor(kBlue);
  hQ2GArCC2Pi->SetLineColor(kRed);
  hQ2GArCC3Pi->SetLineColor(kGreen+2);
  hQ2GArCCHiPi->SetLineColor(kMagenta);
  hsQ2GArGenieNuWroRatios->Add(hQ2GArCC0Pi);
  hsQ2GArGenieNuWroRatios->Add(hQ2GArCC1Pi);
  hsQ2GArGenieNuWroRatios->Add(hQ2GArCC2Pi);
  hsQ2GArGenieNuWroRatios->Add(hQ2GArCC3Pi);
  hsQ2GArGenieNuWroRatios->Add(hQ2GArCCHiPi);
  hsQ2GArGenieNuWroRatios->Write();
  hQ2GArCC0Pi->Write("hQ2GArCC0PiRatio");
  hQ2GArCC1Pi->Write("hQ2GArCC1PiRatio");
  hQ2GArCC2Pi->Write("hQ2GArCC2PiRatio");
  hQ2GArCC3Pi->Write("hQ2GArCC3PiRatio");
  hQ2GArCCHiPi->Write("hQ2GArCCHiPiRatio");

  // Enu reco samples
  TH1 *hnumufhc   = predGArNumuFHC.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  TH1 *hnumufhc0  = predGArNumuFHC0.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  TH1 *hnumufhc1  = predGArNumuFHC1.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  TH1 *hnumufhc2  = predGArNumuFHC2.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  TH1 *hnumufhc3  = predGArNumuFHC3.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  TH1 *hnumufhcHi = predGArNumuFHCHi.Predict(0).FakeData(pot_nd).ToTH1(pot_nd, kPOT, kBinDensity);
  setHistAttr(hnumufhc);
  setHistAttr(hnumufhc0);
  setHistAttr(hnumufhc1);
  setHistAttr(hnumufhc2);
  setHistAttr(hnumufhc3);
  setHistAttr(hnumufhcHi);
  hnumufhc->SetTitle("Inclusive");
  hnumufhc0->SetTitle("0#pi");
  hnumufhc1->SetTitle("1#pi");
  hnumufhc2->SetTitle("2#pi");
  hnumufhc3->SetTitle("3#pi");
  hnumufhcHi->SetTitle(">3#pi");
  hnumufhc->Write("hnumufhc");
  hnumufhc0->Write("hnumufhc0");
  hnumufhc1->Write("hnumufhc1");
  hnumufhc2->Write("hnumufhc2");
  hnumufhc3->Write("hnumufhc3");
  hnumufhcHi->Write("hnumufhcHi");
  THStack *hsEnuReco = new THStack("hsEnuReco", "HPgTPC #nu_{#mu} CC (FHC); E_{#nu, reco} / GeV; Events / GeV");
  hnumufhc0->SetLineColor(kBlack);
  hnumufhc0->SetFillColor(kBlack);
  hnumufhc1->SetLineColor(kBlue);
  hnumufhc1->SetFillColor(kBlue);
  hnumufhc2->SetLineColor(kRed);
  hnumufhc2->SetFillColor(kRed);
  hnumufhc3->SetLineColor(kGreen+2);
  hnumufhc3->SetFillColor(kGreen+2);
  hnumufhcHi->SetLineColor(kMagenta);
  hnumufhcHi->SetFillColor(kMagenta);
  hsEnuReco->Add(hnumufhc0);
  hsEnuReco->Add(hnumufhc1);
  hsEnuReco->Add(hnumufhc2);
  hsEnuReco->Add(hnumufhc3);
  hsEnuReco->Add(hnumufhcHi);
  hsEnuReco->Write();

  // Inclusive separated by interaction mode
  TH1* hnumufhc_qe  = predGArNumuFHC_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc_dis = predGArNumuFHC_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc_res = predGArNumuFHC_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc_coh = predGArNumuFHC_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc_mec = predGArNumuFHC_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhc_qe->SetTitle("Quasi-elastic");
  hnumufhc_dis->SetTitle("DIS");
  hnumufhc_res->SetTitle("Resonant");
  hnumufhc_coh->SetTitle("Coherent");
  hnumufhc_mec->SetTitle("MEC");
  hnumufhc_qe->Write("hnumufhc_qe");
  hnumufhc_dis->Write("hnumufhc_dis");
  hnumufhc_res->Write("hnumufhc_res");
  hnumufhc_coh->Write("hnumufhc_coh");
  hnumufhc_mec->Write("hnumufhc_mec");
  hnumufhc_qe->SetLineColor(kBlack);
  hnumufhc_qe->SetFillColor(kBlack);
  hnumufhc_dis->SetLineColor(kBlue);
  hnumufhc_dis->SetFillColor(kBlue);
  hnumufhc_res->SetLineColor(kRed);
  hnumufhc_res->SetFillColor(kRed);
  hnumufhc_coh->SetLineColor(kOrange+1);
  hnumufhc_coh->SetFillColor(kOrange+1);
  hnumufhc_mec->SetLineColor(kGreen+2);
  hnumufhc_mec->SetFillColor(kGreen+2);
  THStack *hsEnuReco_mode = new THStack("hsEnuReco_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsEnuReco_mode->Add(hnumufhc_coh);
  hsEnuReco_mode->Add(hnumufhc_qe);  
  hsEnuReco_mode->Add(hnumufhc_mec);
  hsEnuReco_mode->Add(hnumufhc_res);
  hsEnuReco_mode->Add(hnumufhc_dis);
  hsEnuReco_mode->Write();
  // 0pi separated by interaction mode
  TH1* hnumufhc0_qe  = predGArNumuFHC0_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc0_dis = predGArNumuFHC0_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc0_res = predGArNumuFHC0_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc0_coh = predGArNumuFHC0_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc0_mec = predGArNumuFHC0_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhc0_qe->SetTitle("Quasi-elastic");
  hnumufhc0_dis->SetTitle("DIS");
  hnumufhc0_res->SetTitle("Resonant");
  hnumufhc0_coh->SetTitle("Coherent");
  hnumufhc0_mec->SetTitle("MEC");
  hnumufhc0_qe->Write("hnumufhc0_qe");
  hnumufhc0_dis->Write("hnumufhc0_dis");
  hnumufhc0_res->Write("hnumufhc0_res");
  hnumufhc0_coh->Write("hnumufhc0_coh");
  hnumufhc0_mec->Write("hnumufhc0_mec");
  hnumufhc0_qe->SetLineColor(kBlack);
  hnumufhc0_qe->SetFillColor(kBlack);
  hnumufhc0_dis->SetLineColor(kBlue);
  hnumufhc0_dis->SetFillColor(kBlue);
  hnumufhc0_res->SetLineColor(kRed);
  hnumufhc0_res->SetFillColor(kRed);
  hnumufhc0_coh->SetLineColor(kOrange+1);
  hnumufhc0_coh->SetFillColor(kOrange+1);
  hnumufhc0_mec->SetLineColor(kGreen+2);
  hnumufhc0_mec->SetFillColor(kGreen+2);
  THStack *hsEnuReco0_mode = new THStack("hsEnuReco0_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsEnuReco0_mode->Add(hnumufhc0_coh);
  hsEnuReco0_mode->Add(hnumufhc0_mec);
  hsEnuReco0_mode->Add(hnumufhc0_qe);  
  hsEnuReco0_mode->Add(hnumufhc0_dis);
  hsEnuReco0_mode->Add(hnumufhc0_res);
  hsEnuReco0_mode->Write();
  // 1pi separated by interaction mode
  TH1* hnumufhc1_qe  = predGArNumuFHC1_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc1_dis = predGArNumuFHC1_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc1_res = predGArNumuFHC1_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc1_coh = predGArNumuFHC1_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc1_mec = predGArNumuFHC1_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhc1_qe->SetTitle("Quasi-elastic");
  hnumufhc1_dis->SetTitle("DIS");
  hnumufhc1_res->SetTitle("Resonant");
  hnumufhc1_coh->SetTitle("Coherent");
  hnumufhc1_mec->SetTitle("MEC");
  hnumufhc1_qe->Write("hnumufhc1_qe");
  hnumufhc1_dis->Write("hnumufhc1_dis");
  hnumufhc1_res->Write("hnumufhc1_res");
  hnumufhc1_coh->Write("hnumufhc1_coh");
  hnumufhc1_mec->Write("hnumufhc1_mec");
  hnumufhc1_qe->SetLineColor(kBlack);
  hnumufhc1_qe->SetFillColor(kBlack);
  hnumufhc1_dis->SetLineColor(kBlue);
  hnumufhc1_dis->SetFillColor(kBlue);
  hnumufhc1_res->SetLineColor(kRed);
  hnumufhc1_res->SetFillColor(kRed);
  hnumufhc1_coh->SetLineColor(kOrange+1);
  hnumufhc1_coh->SetFillColor(kOrange+1);
  hnumufhc1_mec->SetLineColor(kGreen+2);
  hnumufhc1_mec->SetFillColor(kGreen+2);
  THStack *hsEnuReco1_mode = new THStack("hsEnuReco1_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsEnuReco1_mode->Add(hnumufhc1_coh);
  hsEnuReco1_mode->Add(hnumufhc1_mec);
  hsEnuReco1_mode->Add(hnumufhc1_qe);  
  hsEnuReco1_mode->Add(hnumufhc1_dis);
  hsEnuReco1_mode->Add(hnumufhc1_res);
  hsEnuReco1_mode->Write();
  // 2pi
  TH1* hnumufhc2_qe  = predGArNumuFHC2_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc2_dis = predGArNumuFHC2_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc2_res = predGArNumuFHC2_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc2_coh = predGArNumuFHC2_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc2_mec = predGArNumuFHC2_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhc2_qe->SetTitle("Quasi-elastic");
  hnumufhc2_dis->SetTitle("DIS");
  hnumufhc2_res->SetTitle("Resonant");
  hnumufhc2_coh->SetTitle("Coherent");
  hnumufhc2_mec->SetTitle("MEC");
  hnumufhc2_qe->Write("hnumufhc2_qe");
  hnumufhc2_dis->Write("hnumufhc2_dis");
  hnumufhc2_res->Write("hnumufhc2_res");
  hnumufhc2_coh->Write("hnumufhc2_coh");
  hnumufhc2_mec->Write("hnumufhc2_mec");
  hnumufhc2_qe->SetLineColor(kBlack);
  hnumufhc2_qe->SetFillColor(kBlack);
  hnumufhc2_dis->SetLineColor(kBlue);
  hnumufhc2_dis->SetFillColor(kBlue);
  hnumufhc2_res->SetLineColor(kRed);
  hnumufhc2_res->SetFillColor(kRed);
  hnumufhc2_coh->SetLineColor(kOrange+1);
  hnumufhc2_coh->SetFillColor(kOrange+1);
  hnumufhc2_mec->SetLineColor(kGreen+2);
  hnumufhc2_mec->SetFillColor(kGreen+2);
  THStack *hsEnuReco2_mode = new THStack("hsEnuReco2_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsEnuReco2_mode->Add(hnumufhc2_coh);
  hsEnuReco2_mode->Add(hnumufhc2_mec);
  hsEnuReco2_mode->Add(hnumufhc2_qe);  
  hsEnuReco2_mode->Add(hnumufhc2_dis);
  hsEnuReco2_mode->Add(hnumufhc2_res);
  hsEnuReco2_mode->Write();
  // 3pi
  TH1* hnumufhc3_qe  = predGArNumuFHC3_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc3_dis = predGArNumuFHC3_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc3_res = predGArNumuFHC3_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc3_coh = predGArNumuFHC3_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhc3_mec = predGArNumuFHC3_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhc3_qe->SetTitle("Quasi-elastic");
  hnumufhc3_dis->SetTitle("DIS");
  hnumufhc3_res->SetTitle("Resonant");
  hnumufhc3_coh->SetTitle("Coherent");
  hnumufhc3_mec->SetTitle("MEC");
  hnumufhc3_qe->Write("hnumufhc3_qe");
  hnumufhc3_dis->Write("hnumufhc3_dis");
  hnumufhc3_res->Write("hnumufhc3_res");
  hnumufhc3_coh->Write("hnumufhc3_coh");
  hnumufhc3_mec->Write("hnumufhc3_mec");
  hnumufhc3_qe->SetLineColor(kBlack);
  hnumufhc3_qe->SetFillColor(kBlack);
  hnumufhc3_dis->SetLineColor(kBlue);
  hnumufhc3_dis->SetFillColor(kBlue);
  hnumufhc3_res->SetLineColor(kRed);
  hnumufhc3_res->SetFillColor(kRed);
  hnumufhc3_coh->SetLineColor(kOrange+1);
  hnumufhc3_coh->SetFillColor(kOrange+1);
  hnumufhc3_mec->SetLineColor(kGreen+2);
  hnumufhc3_mec->SetFillColor(kGreen+2);
  THStack *hsEnuReco3_mode = new THStack("hsEnuReco3_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsEnuReco3_mode->Add(hnumufhc3_coh);
  hsEnuReco3_mode->Add(hnumufhc3_mec);
  hsEnuReco3_mode->Add(hnumufhc3_qe);  
  hsEnuReco3_mode->Add(hnumufhc3_dis);
  hsEnuReco3_mode->Add(hnumufhc3_res);
  hsEnuReco3_mode->Write();
  // >3pi
  TH1* hnumufhcHi_qe  = predGArNumuFHCHi_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhcHi_dis = predGArNumuFHCHi_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhcHi_res = predGArNumuFHCHi_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhcHi_coh = predGArNumuFHCHi_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1* hnumufhcHi_mec = predGArNumuFHCHi_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhcHi_qe->SetTitle("Quasi-elastic");
  hnumufhcHi_dis->SetTitle("DIS");
  hnumufhcHi_res->SetTitle("Resonant");
  hnumufhcHi_coh->SetTitle("Coherent");
  hnumufhcHi_mec->SetTitle("MEC");
  hnumufhcHi_qe->Write("hnumufhcHi_qe");
  hnumufhcHi_dis->Write("hnumufhcHi_dis");
  hnumufhcHi_res->Write("hnumufhcHi_res");
  hnumufhcHi_coh->Write("hnumufhcHi_coh");
  hnumufhcHi_mec->Write("hnumufhcHi_mec");
  hnumufhcHi_qe->SetLineColor(kBlack);
  hnumufhcHi_qe->SetFillColor(kBlack);
  hnumufhcHi_dis->SetLineColor(kBlue);
  hnumufhcHi_dis->SetFillColor(kBlue);
  hnumufhcHi_res->SetLineColor(kRed);
  hnumufhcHi_res->SetFillColor(kRed);
  hnumufhcHi_coh->SetLineColor(kOrange+1);
  hnumufhcHi_coh->SetFillColor(kOrange+1);
  hnumufhcHi_mec->SetLineColor(kGreen+2);
  hnumufhcHi_mec->SetFillColor(kGreen+2);
  THStack *hsERecoHi_mode = new THStack("hsERecoHi_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu, reco} / GeV; Events / GeV");
  hsERecoHi_mode->Add(hnumufhcHi_coh);
  hsERecoHi_mode->Add(hnumufhcHi_mec);
  hsERecoHi_mode->Add(hnumufhcHi_qe);  
  hsERecoHi_mode->Add(hnumufhcHi_dis);
  hsERecoHi_mode->Add(hnumufhcHi_res);
  hsERecoHi_mode->Write();

  // True energy
  TH1 *hnumufhcEnu_qe = predGArNumuFHCEnu_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1 *hnumufhcEnu_dis = predGArNumuFHCEnu_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1 *hnumufhcEnu_res = predGArNumuFHCEnu_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1 *hnumufhcEnu_coh = predGArNumuFHCEnu_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  TH1 *hnumufhcEnu_mec = predGArNumuFHCEnu_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd,kPOT,kBinDensity);
  hnumufhcEnu_qe->SetTitle("Quasi-elastic");
  hnumufhcEnu_dis->SetTitle("DIS");
  hnumufhcEnu_res->SetTitle("Resonant");
  hnumufhcEnu_coh->SetTitle("Coherent");
  hnumufhcEnu_mec->SetTitle("MEC");
  hnumufhcEnu_qe->Write("hnumufhcEnu_qe");
  hnumufhcEnu_dis->Write("hnumufhcEnu_dis");
  hnumufhcEnu_res->Write("hnumufhcEnu_res");
  hnumufhcEnu_coh->Write("hnumufhcEnu_coh");
  hnumufhcEnu_mec->Write("hnumufhcEnu_mec");
  hnumufhcEnu_qe->SetLineColor(kBlack);
  hnumufhcEnu_qe->SetFillColor(kBlack);
  hnumufhcEnu_dis->SetLineColor(kBlue);
  hnumufhcEnu_dis->SetFillColor(kBlue);
  hnumufhcEnu_res->SetLineColor(kRed);
  hnumufhcEnu_res->SetFillColor(kRed);
  hnumufhcEnu_coh->SetLineColor(kOrange+1);
  hnumufhcEnu_coh->SetFillColor(kOrange+1);
  hnumufhcEnu_mec->SetLineColor(kGreen+2);
  hnumufhcEnu_mec->SetFillColor(kGreen+2);
  THStack *hsEnu_mode = new THStack("hsEnu_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; E_{#nu} / GeV; Events / GeV");
  hsEnu_mode->Add(hnumufhcEnu_coh);
  hsEnu_mode->Add(hnumufhcEnu_qe);  
  hsEnu_mode->Add(hnumufhcEnu_mec);
  hsEnu_mode->Add(hnumufhcEnu_res);
  hsEnu_mode->Add(hnumufhcEnu_dis);
  hsEnu_mode->Write();
 
  // Q2
  TH1 *hnumufhcQ2_qe = predGArNumuFHCQ2_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2_dis = predGArNumuFHCQ2_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2_res = predGArNumuFHCQ2_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2_coh = predGArNumuFHCQ2_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2_mec = predGArNumuFHCQ2_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  hnumufhcQ2_qe->SetTitle("Quasi-elastic");
  hnumufhcQ2_dis->SetTitle("DIS");
  hnumufhcQ2_res->SetTitle("Resonant");
  hnumufhcQ2_coh->SetTitle("Coherent");
  hnumufhcQ2_mec->SetTitle("MEC");
  hnumufhcQ2_qe->Write("hnumufhcQ2_qe");
  hnumufhcQ2_dis->Write("hnumufhcQ2_dis");
  hnumufhcQ2_res->Write("hnumufhcQ2_res");
  hnumufhcQ2_coh->Write("hnumufhcQ2_coh");
  hnumufhcQ2_mec->Write("hnumufhcQ2_mec");
  hnumufhcQ2_qe->SetLineColor(kBlack);
  hnumufhcQ2_qe->SetFillColor(kBlack);
  hnumufhcQ2_dis->SetLineColor(kBlue);
  hnumufhcQ2_dis->SetFillColor(kBlue);
  hnumufhcQ2_res->SetLineColor(kRed);
  hnumufhcQ2_res->SetFillColor(kRed);
  hnumufhcQ2_coh->SetLineColor(kOrange+1);
  hnumufhcQ2_coh->SetFillColor(kOrange+1);
  hnumufhcQ2_mec->SetLineColor(kGreen+2);
  hnumufhcQ2_mec->SetFillColor(kGreen+2);
  THStack *hsQ2_mode = new THStack("hsQ2_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; Q^{2} / GeV; Events");
  hsQ2_mode->Add(hnumufhcQ2_coh);
  hsQ2_mode->Add(hnumufhcQ2_qe);  
  hsQ2_mode->Add(hnumufhcQ2_mec);
  hsQ2_mode->Add(hnumufhcQ2_res);
  hsQ2_mode->Add(hnumufhcQ2_dis);
  hsQ2_mode->Write();

  // Reco Q2
  TH1 *hnumufhcQ2Reco_qe = predGArNumuFHCQ2Reco_qe.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2Reco_dis = predGArNumuFHCQ2Reco_dis.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2Reco_res = predGArNumuFHCQ2Reco_res.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2Reco_coh = predGArNumuFHCQ2Reco_coh.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hnumufhcQ2Reco_mec = predGArNumuFHCQ2Reco_mec.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  hnumufhcQ2Reco_qe->SetTitle("Quasi-elastic");
  hnumufhcQ2Reco_dis->SetTitle("DIS");
  hnumufhcQ2Reco_res->SetTitle("Resonant");
  hnumufhcQ2Reco_coh->SetTitle("Coherent");
  hnumufhcQ2Reco_mec->SetTitle("MEC");
  hnumufhcQ2Reco_qe->Write("hnumufhcQ2Reco_qe");
  hnumufhcQ2Reco_dis->Write("hnumufhcQ2Reco_dis");
  hnumufhcQ2Reco_res->Write("hnumufhcQ2Reco_res");
  hnumufhcQ2Reco_coh->Write("hnumufhcQ2Reco_coh");
  hnumufhcQ2Reco_mec->Write("hnumufhcQ2Reco_mec");
  hnumufhcQ2Reco_qe->SetLineColor(kBlack);
  hnumufhcQ2Reco_qe->SetFillColor(kBlack);
  hnumufhcQ2Reco_dis->SetLineColor(kBlue);
  hnumufhcQ2Reco_dis->SetFillColor(kBlue);
  hnumufhcQ2Reco_res->SetLineColor(kRed);
  hnumufhcQ2Reco_res->SetFillColor(kRed);
  hnumufhcQ2Reco_coh->SetLineColor(kOrange+1);
  hnumufhcQ2Reco_coh->SetFillColor(kOrange+1);
  hnumufhcQ2Reco_mec->SetLineColor(kGreen+2);
  hnumufhcQ2Reco_mec->SetFillColor(kGreen+2);
  THStack *hsQ2Reco_mode = new THStack("hsQ2Reco_mode", "HPgTPC #nu_{#mu} CC (FHC) separated by interaction mode; Q^{2}_{reco} / GeV; Events");
  hsQ2Reco_mode->Add(hnumufhcQ2Reco_coh);
  hsQ2Reco_mode->Add(hnumufhcQ2Reco_qe);  
  hsQ2Reco_mode->Add(hnumufhcQ2Reco_mec);
  hsQ2Reco_mode->Add(hnumufhcQ2Reco_res);
  hsQ2Reco_mode->Add(hnumufhcQ2Reco_dis);
  hsQ2Reco_mode->Write();

  fout->Close();
  delete fout;

} // leadPiGAr

