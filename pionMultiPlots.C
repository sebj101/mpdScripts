// pionMultiPlots.C
// Takes a look at the LAr and GAr samples energy reconstruction as function of pion multiplicity
// Also checks pion multiplicity as function of kinematic variables

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
#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Systs/XSecSysts.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/TDRLoaders.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
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

void pionMultiPlots(const char *outFile, bool doLAr=false, 
		    const char *lardir="/pnfs/dune/persistent/users/picker24/CAFv4/", 
		    const char *gardir="/dune/data/users/sbjones/gasTpcCAF/v2/") 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  // ND vars
  // Reco X
  const Var kRecoX({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
		   [](const caf::StandardRecord* sr) {
		     double x = 0;
		     x = (4 * sr->dune.Ev_reco * sr->dune.Elep_reco * TMath::Sin(sr->dune.theta_reco/2.) * TMath::Sin(sr->dune.theta_reco/2.)) / (2 * 0.939 * (sr->dune.Ev_reco - sr->dune.Elep_reco));
		     return x;
		   });
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
  const Var kRecoEnergyND  = SIMPLEVAR(dune.Ev_reco);
  // const Var kRecoYND       = (SIMPLEVAR(dune.Ev_reco)-SIMPLEVAR(dune.Elep_reco))/SIMPLEVAR(dune.Ev_reco);
  const Var kTrueYND       = (SIMPLEVAR(dune.Ev)-SIMPLEVAR(dune.LepE))/SIMPLEVAR(dune.Ev);
  const Var kTrueEnergy    = SIMPLEVAR(dune.Ev);
  const Var kTrueLepEnergy = SIMPLEVAR(dune.LepE);
  const Var kNPi = SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipi0);
  const Var kQ2  = SIMPLEVAR(dune.Q2);
  const Var kW   = SIMPLEVAR(dune.W);
  const Var kFHC = SIMPLEVAR(dune.isFHC);
  const Var Enu_reco_numu = SIMPLEVAR(dune.Ev_reco_numu);
  const Var Enu_reco_nue  = SIMPLEVAR(dune.Ev_reco_nue);
  const Var isCC = SIMPLEVAR(dune.isCC);
  const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
  const Var KPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
  const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
  const Var kRecoPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
  const Var kThetaReco = SIMPLEVAR(dune.theta_reco);
  // CV weighting
  //const Var kGENIEWeights = SIMPLEVAR(dune.total_xsSyst_cv_wgt);

  std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
				   2., 2.25, 2.5, 2.75, 
				   3., 3.25, 3.5, 3.75, 
				   4., 5., 6., 10.};
  std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
  const Binning binsEreco  = Binning::Custom(binEEdges);
  const Binning simpleBins   = Binning::Simple(20, 0, 5);
  const Binning simpleWBins  = Binning::Simple(30, 0, 6);
  const Binning simpleQ2Bins = Binning::Simple(30, 0, 6);
  const HistAxis axND("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND);
  const HistAxis axTrueRecoND("E_{#nu, true} (GeV)", simpleBins, kTrueEnergy,
			      "E_{#nu, reco} (GeV)", simpleBins, kRecoEnergyND);
  const HistAxis axQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kQ2);
  const HistAxis axW("W (GeV)", simpleWBins, kW);
  const HistAxis axRecoQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kRecoQ2);
  const HistAxis axRecoW("W (GeV)", simpleWBins, kRecoW);
  const HistAxis axisnumu("Reco #nu energy (GeV)", binsEreco, Enu_reco_numu);
  const HistAxis axisnue("Reco #nu energy (GeV)", binsEreco, Enu_reco_nue);
  // POT for 3.5 years
  const double pot_fd = 3.5 * POT120 * 40/1.13;
  const double pot_nd = 3.5 * POT120;
  // This is pretty annoying, but the above is for 7 years staged, which is 336 kT MW yr
  const double nom_exposure = 336.;

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  
  std::vector<const ISyst*> systlist = GetXSecSysts({"NuWroReweightFakeData"});

  TDRLoaders loadersFHC(TDRLoaders::kFHC);
  TDRLoaders loadersRHC(TDRLoaders::kRHC); 
  NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV); 
  NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV); 
  PredictionInterp predNumuFHCreco(systlist, this_calc, gennumureco, loadersFHC); 
  PredictionInterp predNumuRHCreco(systlist, this_calc, gennumureco, loadersRHC); 
  PredictionInterp predNueFHCreco(systlist, this_calc, gennuereco, loadersFHC); 
  PredictionInterp predNueRHCreco(systlist, this_calc, gennuereco, loadersRHC);
  
  Loaders loadersLArFHC;
  Loaders loadersLArRHC;
  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderLArFHC(Form("%s/ND_FHC_CAF.root", lardir), kBeam);
  SpectrumLoader loaderLArRHC(Form("%s/ND_RHC_CAF.root", lardir), kBeam);
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", gardir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", gardir), kBeam);
  
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);
  loadersLArFHC.AddLoader(&loaderLArFHC, caf::kNEARDET, Loaders::kMC);
  loadersLArRHC.AddLoader(&loaderLArRHC, caf::kNEARDET, Loaders::kMC);

  loadersFHC.Go();
  loadersRHC.Go();  

  TFile *fout = new TFile(outFile, "recreate");

  NoOscPredictionGenerator genLAr1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1);
  NoOscPredictionGenerator genLArCC0Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator genLArCC1Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator genLArCC2Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator genLArCC3Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);
  NoOscPredictionGenerator genLArCCHiPi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1 && kNPi > 3);

  NoOscPredictionGenerator genLAr(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1);
  NoOscPredictionGenerator genLArCC0Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator genLArCC1Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator genLArCC2Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator genLArCC3Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);
  NoOscPredictionGenerator genLArCCHiPi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1 && kNPi > 3);

  NoOscPredictionGenerator WLAr(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1);
  NoOscPredictionGenerator WLArCC0Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator WLArCC1Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator WLArCC2Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator WLArCC3Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);
  NoOscPredictionGenerator WLArCCHiPi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1 && kNPi > 3);

  NoOscPredictionGenerator Q2LAr(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1);
  NoOscPredictionGenerator Q2LArCC0Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator Q2LArCC1Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator Q2LArCC2Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator Q2LArCC3Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);
  NoOscPredictionGenerator Q2LArCCHiPi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && isCC==1 && kNPi > 3);

  // FHC
  PredictionInterp predLAr(systlist, this_calc, genLAr, loadersLArFHC);
  PredictionInterp predLArCC0Pi(systlist, this_calc, genLArCC0Pi, loadersLArFHC);
  PredictionInterp predLArCC1Pi(systlist, this_calc, genLArCC1Pi, loadersLArFHC);
  PredictionInterp predLArCC2Pi(systlist, this_calc, genLArCC2Pi, loadersLArFHC);
  PredictionInterp predLArCC3Pi(systlist, this_calc, genLArCC3Pi, loadersLArFHC);
  PredictionInterp predLArCCHiPi(systlist, this_calc, genLArCCHiPi, loadersLArFHC);

  PredictionInterp predLAr1d(systlist, this_calc, genLAr1d, loadersLArFHC);
  PredictionInterp predLArCC0Pi1d(systlist, this_calc, genLArCC0Pi1d, loadersLArFHC);
  PredictionInterp predLArCC1Pi1d(systlist, this_calc, genLArCC1Pi1d, loadersLArFHC);
  PredictionInterp predLArCC2Pi1d(systlist, this_calc, genLArCC2Pi1d, loadersLArFHC);
  PredictionInterp predLArCC3Pi1d(systlist, this_calc, genLArCC3Pi1d, loadersLArFHC);
  PredictionInterp predLArCCHiPi1d(systlist, this_calc, genLArCCHiPi1d, loadersLArFHC);

  PredictionInterp predWLAr(systlist, this_calc, WLAr, loadersLArFHC);
  PredictionInterp predWLArCC0Pi(systlist, this_calc, WLArCC0Pi, loadersLArFHC);
  PredictionInterp predWLArCC1Pi(systlist, this_calc, WLArCC1Pi, loadersLArFHC);
  PredictionInterp predWLArCC2Pi(systlist, this_calc, WLArCC2Pi, loadersLArFHC);
  PredictionInterp predWLArCC3Pi(systlist, this_calc, WLArCC3Pi, loadersLArFHC);
  PredictionInterp predWLArCCHiPi(systlist, this_calc, WLArCCHiPi, loadersLArFHC);

  PredictionInterp predQ2LAr(systlist, this_calc, Q2LAr, loadersLArFHC);
  PredictionInterp predQ2LArCC0Pi(systlist, this_calc, Q2LArCC0Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC1Pi(systlist, this_calc, Q2LArCC1Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC2Pi(systlist, this_calc, Q2LArCC2Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC3Pi(systlist, this_calc, Q2LArCC3Pi, loadersLArFHC);
  PredictionInterp predQ2LArCCHiPi(systlist, this_calc, Q2LArCCHiPi, loadersLArFHC);
  // RHC
  PredictionInterp predLArRHC(systlist, this_calc, genLAr, loadersLArRHC);
  PredictionInterp predLArCC0PiRHC(systlist, this_calc, genLArCC0Pi, loadersLArRHC);
  PredictionInterp predLArCC1PiRHC(systlist, this_calc, genLArCC1Pi, loadersLArRHC);
  PredictionInterp predLArCC2PiRHC(systlist, this_calc, genLArCC2Pi, loadersLArRHC);
  PredictionInterp predLArCC3PiRHC(systlist, this_calc, genLArCC3Pi, loadersLArRHC);

  PredictionInterp predLAr1dRHC(systlist, this_calc, genLAr1d, loadersLArRHC);
  PredictionInterp predLArCC0Pi1dRHC(systlist, this_calc, genLArCC0Pi1d, loadersLArRHC);
  PredictionInterp predLArCC1Pi1dRHC(systlist, this_calc, genLArCC1Pi1d, loadersLArRHC);
  PredictionInterp predLArCC2Pi1dRHC(systlist, this_calc, genLArCC2Pi1d, loadersLArRHC);
  PredictionInterp predLArCC3Pi1dRHC(systlist, this_calc, genLArCC3Pi1d, loadersLArRHC);

  PredictionInterp predWLArRHC(systlist, this_calc, WLAr, loadersLArRHC);
  PredictionInterp predWLArCC0PiRHC(systlist, this_calc, WLArCC0Pi, loadersLArRHC);
  PredictionInterp predWLArCC1PiRHC(systlist, this_calc, WLArCC1Pi, loadersLArRHC);
  PredictionInterp predWLArCC2PiRHC(systlist, this_calc, WLArCC2Pi, loadersLArRHC);
  PredictionInterp predWLArCC3PiRHC(systlist, this_calc, WLArCC3Pi, loadersLArRHC);

  PredictionInterp predQ2LArRHC(systlist, this_calc, Q2LAr, loadersLArRHC);
  PredictionInterp predQ2LArCC0PiRHC(systlist, this_calc, Q2LArCC0Pi, loadersLArRHC);
  PredictionInterp predQ2LArCC1PiRHC(systlist, this_calc, Q2LArCC1Pi, loadersLArRHC);
  PredictionInterp predQ2LArCC2PiRHC(systlist, this_calc, Q2LArCC2Pi, loadersLArRHC);
  PredictionInterp predQ2LArCC3PiRHC(systlist, this_calc, Q2LArCC3Pi, loadersLArRHC);

  // GAr ND samples
  NoOscPredictionGenerator genGAr(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator genGArCC0Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArCC1Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArCC2Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArCC3Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArCCHiPi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);

  NoOscPredictionGenerator genGAr1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator genGArCC0Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator genGArCC1Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator genGArCC2Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator genGArCC3Pi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator genGArCCHiPi1d(axND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);

  NoOscPredictionGenerator WGAr(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator WGArCC0Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator WGArCC1Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator WGArCC2Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator WGArCC3Pi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator WGArCCHiPi(axRecoW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);

  NoOscPredictionGenerator Q2GAr(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV);
  NoOscPredictionGenerator Q2GArCC0Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==0);
  NoOscPredictionGenerator Q2GArCC1Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==1);
  NoOscPredictionGenerator Q2GArCC2Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==2);
  NoOscPredictionGenerator Q2GArCC3Pi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi==3);
  NoOscPredictionGenerator Q2GArCCHiPi(axRecoQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueGasFV && kRecoPi>3);

  // FHC
  PredictionInterp predGAr(systlist, this_calc, genGAr, loadersGArFHC);
  PredictionInterp predGArCC0Pi(systlist, this_calc, genGArCC0Pi, loadersGArFHC);
  PredictionInterp predGArCC1Pi(systlist, this_calc, genGArCC1Pi, loadersGArFHC);
  PredictionInterp predGArCC2Pi(systlist, this_calc, genGArCC2Pi, loadersGArFHC);
  PredictionInterp predGArCC3Pi(systlist, this_calc, genGArCC3Pi, loadersGArFHC);
  PredictionInterp predGArCCHiPi(systlist, this_calc, genGArCCHiPi, loadersGArFHC);

  PredictionInterp predGAr1d(systlist, this_calc, genGAr1d, loadersGArFHC);
  PredictionInterp predGArCC0Pi1d(systlist, this_calc, genGArCC0Pi1d, loadersGArFHC);
  PredictionInterp predGArCC1Pi1d(systlist, this_calc, genGArCC1Pi1d, loadersGArFHC);
  PredictionInterp predGArCC2Pi1d(systlist, this_calc, genGArCC2Pi1d, loadersGArFHC);
  PredictionInterp predGArCC3Pi1d(systlist, this_calc, genGArCC3Pi1d, loadersGArFHC);
  PredictionInterp predGArCCHiPi1d(systlist, this_calc, genGArCCHiPi1d, loadersGArFHC);

  PredictionInterp predWGAr(systlist, this_calc, WGAr, loadersGArFHC);
  PredictionInterp predWGArCC0Pi(systlist, this_calc, WGArCC0Pi, loadersGArFHC);
  PredictionInterp predWGArCC1Pi(systlist, this_calc, WGArCC1Pi, loadersGArFHC);
  PredictionInterp predWGArCC2Pi(systlist, this_calc, WGArCC2Pi, loadersGArFHC);
  PredictionInterp predWGArCC3Pi(systlist, this_calc, WGArCC3Pi, loadersGArFHC);
  PredictionInterp predWGArCCHiPi(systlist, this_calc, WGArCCHiPi, loadersGArFHC);

  PredictionInterp predQ2GAr(systlist, this_calc, Q2GAr, loadersGArFHC);
  PredictionInterp predQ2GArCC0Pi(systlist, this_calc, Q2GArCC0Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC1Pi(systlist, this_calc, Q2GArCC1Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC2Pi(systlist, this_calc, Q2GArCC2Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC3Pi(systlist, this_calc, Q2GArCC3Pi, loadersGArFHC);
  PredictionInterp predQ2GArCCHiPi(systlist, this_calc, Q2GArCCHiPi, loadersGArFHC);
  // RHC
  PredictionInterp predGArRHC(systlist, this_calc, genGAr, loadersGArRHC);
  PredictionInterp predGArCC0PiRHC(systlist, this_calc, genGArCC0Pi, loadersGArRHC);
  PredictionInterp predGArCC1PiRHC(systlist, this_calc, genGArCC1Pi, loadersGArRHC);
  PredictionInterp predGArCC2PiRHC(systlist, this_calc, genGArCC2Pi, loadersGArRHC);
  PredictionInterp predGArCC3PiRHC(systlist, this_calc, genGArCC3Pi, loadersGArRHC);

  PredictionInterp predGAr1dRHC(systlist, this_calc, genGAr1d, loadersGArRHC);
  PredictionInterp predGArCC0Pi1dRHC(systlist, this_calc, genGArCC0Pi1d, loadersGArRHC);
  PredictionInterp predGArCC1Pi1dRHC(systlist, this_calc, genGArCC1Pi1d, loadersGArRHC);
  PredictionInterp predGArCC2Pi1dRHC(systlist, this_calc, genGArCC2Pi1d, loadersGArRHC);
  PredictionInterp predGArCC3Pi1dRHC(systlist, this_calc, genGArCC3Pi1d, loadersGArRHC);

  PredictionInterp predWGArRHC(systlist, this_calc, WGAr, loadersGArRHC);
  PredictionInterp predWGArCC0PiRHC(systlist, this_calc, WGArCC0Pi, loadersGArRHC);
  PredictionInterp predWGArCC1PiRHC(systlist, this_calc, WGArCC1Pi, loadersGArRHC);
  PredictionInterp predWGArCC2PiRHC(systlist, this_calc, WGArCC2Pi, loadersGArRHC);
  PredictionInterp predWGArCC3PiRHC(systlist, this_calc, WGArCC3Pi, loadersGArRHC);

  PredictionInterp predQ2GArRHC(systlist, this_calc, Q2GAr, loadersGArRHC);
  PredictionInterp predQ2GArCC0PiRHC(systlist, this_calc, Q2GArCC0Pi, loadersGArRHC);
  PredictionInterp predQ2GArCC1PiRHC(systlist, this_calc, Q2GArCC1Pi, loadersGArRHC);
  PredictionInterp predQ2GArCC2PiRHC(systlist, this_calc, Q2GArCC2Pi, loadersGArRHC);
  PredictionInterp predQ2GArCC3PiRHC(systlist, this_calc, Q2GArCC3Pi, loadersGArRHC);

  loadersGArFHC.Go();
  loadersGArRHC.Go();
  if (doLAr) {
    loadersLArFHC.Go();
    loadersLArRHC.Go();
  }
  std::cout<<systlist.size()<<std::endl;
  if (doLAr) {
    TH2 *h2LAr           = predLAr.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LAr);
    TH2 *h2LArNuWro      = predLAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArNuWro);
    TH2 *h2LArCC0Pi      = predLArCC0Pi.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC0Pi);
    TH2 *h2LArCC0PiNuWro = predLArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC0PiNuWro);
    TH2 *h2LArCC1Pi      = predLArCC1Pi.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC1Pi);
    TH2 *h2LArCC1PiNuWro = predLArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC1PiNuWro);
    TH2 *h2LArCC2Pi      = predLArCC2Pi.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC2Pi);
    TH2 *h2LArCC2PiNuWro = predLArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC2PiNuWro);
    TH2 *h2LArCC3Pi      = predLArCC3Pi.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC3Pi);
    TH2 *h2LArCC3PiNuWro = predLArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC3PiNuWro);
    TH2 *h2LArCCHiPi      = predLArCCHiPi.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCCHiPi);
    TH2 *h2LArCCHiPiNuWro = predLArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCCHiPiNuWro);

    TH1 *hLAr1d           = predLAr1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLAr1d);
    TH1 *hLArNuWro1d      = predLAr1d.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    setHistAttr(hLArNuWro1d);
    TH1 *hLArCC0Pi1d      = predLArCC0Pi1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLArCC0Pi1d);
    TH1 *hLArCC0PiNuWro1d = predLArCC0Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
    setHistAttr(hLArCC0PiNuWro1d);
    TH1 *hLArCC1Pi1d      = predLArCC1Pi1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLArCC1Pi1d);
    TH1 *hLArCC1PiNuWro1d = predLArCC1Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
    setHistAttr(hLArCC1PiNuWro1d);
    TH1 *hLArCC2Pi1d      = predLArCC2Pi1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLArCC2Pi1d);
    TH1 *hLArCC2PiNuWro1d = predLArCC2Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
    setHistAttr(hLArCC2PiNuWro1d);
    TH1 *hLArCC3Pi1d      = predLArCC3Pi1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLArCC3Pi1d);
    TH1 *hLArCC3PiNuWro1d = predLArCC3Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
    setHistAttr(hLArCC3PiNuWro1d);
    TH1 *hLArCCHiPi1d      = predLArCCHiPi1d.Predict(0).ToTH1(pot_nd);
    setHistAttr(hLArCCHiPi1d);
    TH1 *hLArCCHiPiNuWro1d = predLArCCHiPi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
    setHistAttr(hLArCCHiPiNuWro1d);

    TH1 *hWLAr      = predWLAr.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArNuWro = predWLAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC0Pi      = predWLArCC0Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC0PiNuWro = predWLArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC1Pi      = predWLArCC1Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC1PiNuWro = predWLArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC2Pi      = predWLArCC2Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC2PiNuWro = predWLArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC3Pi      = predWLArCC3Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC3PiNuWro = predWLArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCCHiPi      = predWLArCCHiPi.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCCHiPiNuWro = predWLArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

    TH1 *hQ2LAr      = predQ2LAr.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArNuWro = predQ2LAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC0Pi      = predQ2LArCC0Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC0PiNuWro = predQ2LArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC1Pi      = predQ2LArCC1Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC1PiNuWro = predQ2LArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC2Pi      = predQ2LArCC2Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC2PiNuWro = predQ2LArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC3Pi      = predQ2LArCC3Pi.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC3PiNuWro = predQ2LArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCCHiPi      = predQ2LArCCHiPi.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCCHiPiNuWro = predQ2LArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

    // LAr RHC
    TH2 *h2LArRHC           = predLArRHC.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArRHC);
    TH2 *h2LArNuWroRHC      = predLArRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArNuWroRHC);
    TH2 *h2LArCC0PiRHC      = predLArCC0PiRHC.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC0PiRHC); 
    TH2 *h2LArCC0PiNuWroRHC = predLArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC0PiNuWroRHC); 
    TH2 *h2LArCC1PiRHC      = predLArCC1PiRHC.Predict(0).ToTH2(pot_nd);
    setHistAttr(h2LArCC1PiRHC); 
    TH2 *h2LArCC1PiNuWroRHC = predLArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    setHistAttr(h2LArCC1PiNuWroRHC); 
    TH2 *h2LArCC2PiRHC      = predLArCC2PiRHC.Predict(0).ToTH2(pot_nd);
    TH2 *h2LArCC2PiNuWroRHC = predLArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
    TH2 *h2LArCC3PiRHC      = predLArCC3PiRHC.Predict(0).ToTH2(pot_nd);
    TH2 *h2LArCC3PiNuWroRHC = predLArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);

    TH1 *hWLArCC0PiRHC      = predWLArCC0PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC0PiNuWroRHC = predWLArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC1PiRHC      = predWLArCC1PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC1PiNuWroRHC = predWLArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC2PiRHC      = predWLArCC2PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC2PiNuWroRHC = predWLArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hWLArCC3PiRHC      = predWLArCC3PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hWLArCC3PiNuWroRHC = predWLArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

    TH1 *hQ2LArCC0PiRHC      = predQ2LArCC0PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC0PiNuWroRHC = predQ2LArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC1PiRHC      = predQ2LArCC1PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC1PiNuWroRHC = predQ2LArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC2PiRHC      = predQ2LArCC2PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC2PiNuWroRHC = predQ2LArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
    TH1 *hQ2LArCC3PiRHC      = predQ2LArCC3PiRHC.Predict(0).ToTH1(pot_nd);
    TH1 *hQ2LArCC3PiNuWroRHC = predQ2LArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

    // LAr FHC
    h2LAr->Write("h2LAr");
    h2LArNuWro->Write("h2LArNuWro");
    h2LArCC0Pi->Write("h2LArCC0Pi");
    h2LArCC0PiNuWro->Write("h2LArCC0PiNuWro");
    h2LArCC1Pi->Write("h2LArCC1Pi");
    h2LArCC1PiNuWro->Write("h2LArCC1PiNuWro");
    h2LArCC2Pi->Write("h2LArCC2Pi");
    h2LArCC2PiNuWro->Write("h2LArCC2PiNuWro");
    h2LArCC3Pi->Write("h2LArCC3Pi");
    h2LArCC3PiNuWro->Write("h2LArCC3PiNuWro");
    h2LArCCHiPi->Write("h2LArCCHiPi");
    h2LArCCHiPiNuWro->Write("h2LArCCHiPiNuWro");

    hLAr1d->Write("hLAr1d");
    hLArNuWro1d->Write("hLArNuWro1d");
    hLArCC0Pi1d->Write("hLArCC0Pi1d");
    hLArCC0PiNuWro1d->Write("hLArCC0PiNuWro1d");
    hLArCC1Pi1d->Write("hLArCC1Pi1d");
    hLArCC1PiNuWro1d->Write("hLArCC1PiNuWro1d");
    hLArCC2Pi1d->Write("hLArCC2Pi1d");
    hLArCC2PiNuWro1d->Write("hLArCC2PiNuWro1d");
    hLArCC3Pi1d->Write("hLArCC3Pi1d");
    hLArCC3PiNuWro1d->Write("hLArCC3PiNuWro1d");
    hLArCCHiPi1d->Write("hLArCCHiPi1d");
    hLArCCHiPiNuWro1d->Write("hLArCCHiPiNuWro1d");

    hWLArCC0Pi->Write("hWLArCC0Pi");
    hWLArCC0PiNuWro->Write("hWLArCC0PiNuWro");
    hWLArCC1Pi->Write("hWLArCC1Pi");
    hWLArCC1PiNuWro->Write("hWLArCC1PiNuWro");
    hWLArCC2Pi->Write("hWLArCC2Pi");
    hWLArCC2PiNuWro->Write("hWLArCC2PiNuWro");
    hWLArCC3Pi->Write("hWLArCC3Pi");
    hWLArCC3PiNuWro->Write("hWLArCC3PiNuWro");
    hWLArCCHiPi->Write("hWLArCCHiPi");
    hWLArCCHiPiNuWro->Write("hWLArCCHiPiNuWro");

    hQ2LArCC0Pi->Write("hQ2LArCC0Pi");
    hQ2LArCC0PiNuWro->Write("hQ2LArCC0PiNuWro");
    hQ2LArCC1Pi->Write("hQ2LArCC1Pi");
    hQ2LArCC1PiNuWro->Write("hQ2LArCC1PiNuWro");
    hQ2LArCC2Pi->Write("hQ2LArCC2Pi");
    hQ2LArCC2PiNuWro->Write("hQ2LArCC2PiNuWro");
    hQ2LArCC3Pi->Write("hQ2LArCC3Pi");
    hQ2LArCC3PiNuWro->Write("hQ2LArCC3PiNuWro");
    hQ2LArCCHiPi->Write("hQ2LArCCHiPi");
    hQ2LArCCHiPiNuWro->Write("hQ2LArCCHiPiNuWro");
    // LAr RHC
    h2LArRHC->Write("h2LArRHC");
    h2LArNuWroRHC->Write("h2LArNuWroRHC");
    h2LArCC0PiRHC->Write("h2LArCC0PiRHC");
    h2LArCC0PiNuWroRHC->Write("h2LArCC0PiNuWroRHC");
    h2LArCC1PiRHC->Write("h2LArCC1PiRHC");
    h2LArCC1PiNuWroRHC->Write("h2LArCC1PiNuWroRHC");
    h2LArCC2PiRHC->Write("h2LArCC2PiRHC");
    h2LArCC2PiNuWroRHC->Write("h2LArCC2PiNuWroRHC");
    h2LArCC3PiRHC->Write("h2LArCC3PiRHC");
    h2LArCC3PiNuWroRHC->Write("h2LArCC3PiNuWroRHC");

    hWLArCC0PiRHC->Write("hWLArCC0PiRHC");
    hWLArCC0PiNuWroRHC->Write("hWLArCC0PiNuWroRHC");
    hWLArCC1PiRHC->Write("hWLArCC1PiRHC");
    hWLArCC1PiNuWroRHC->Write("hWLArCC1PiNuWroRHC");
    hWLArCC2PiRHC->Write("hWLArCC2PiRHC");
    hWLArCC2PiNuWroRHC->Write("hWLArCC2PiNuWroRHC");
    hWLArCC3PiRHC->Write("hWLArCC3PiRHC");
    hWLArCC3PiNuWroRHC->Write("hWLArCC3PiNuWroRHC");

    hQ2LArCC0PiRHC->Write("hQ2LArCC0PiRHC");
    hQ2LArCC0PiNuWroRHC->Write("hQ2LArCC0PiNuWroRHC");
    hQ2LArCC1PiRHC->Write("hQ2LArCC1PiRHC");
    hQ2LArCC1PiNuWroRHC->Write("hQ2LArCC1PiNuWroRHC");
    hQ2LArCC2PiRHC->Write("hQ2LArCC2PiRHC");
    hQ2LArCC2PiNuWroRHC->Write("hQ2LArCC2PiNuWroRHC");
    hQ2LArCC3PiRHC->Write("hQ2LArCC3PiRHC");
    hQ2LArCC3PiNuWroRHC->Write("hQ2LArCC3PiNuWroRHC");

    THStack *hsWLArCC0Pi = new THStack("hsWLArCC0Pi", "W in LAr ND, CC0#pi; W; Events");
    THStack *hsWLArCC1Pi = new THStack("hsWLArCC1Pi", "W in LAr ND, CC1#pi; W; Events");
    THStack *hsWLArCC2Pi = new THStack("hsWLArCC2Pi", "W in LAr ND, CC2#pi; W; Events");
    THStack *hsWLArCC3Pi = new THStack("hsWLArCC3Pi", "W in LAr ND, CC3#pi; W; Events");
    hWLArCC0PiNuWro->SetLineColor(kRed);
    hWLArCC1PiNuWro->SetLineColor(kRed);
    hWLArCC2PiNuWro->SetLineColor(kRed);
    hWLArCC3PiNuWro->SetLineColor(kRed);
    hsWLArCC0Pi->Add(hWLArCC0Pi);
    hsWLArCC0Pi->Add(hWLArCC0PiNuWro);
    hsWLArCC1Pi->Add(hWLArCC1Pi);
    hsWLArCC1Pi->Add(hWLArCC1PiNuWro);
    hsWLArCC2Pi->Add(hWLArCC2Pi);
    hsWLArCC2Pi->Add(hWLArCC2PiNuWro);
    hsWLArCC3Pi->Add(hWLArCC3Pi);
    hsWLArCC3Pi->Add(hWLArCC3PiNuWro);
    hsWLArCC0Pi->Write("hsWLArCC0Pi");
    hsWLArCC1Pi->Write("hsWLArCC1Pi");
    hsWLArCC2Pi->Write("hsWLArCC2Pi");
    hsWLArCC3Pi->Write("hsWLArCC3Pi");

    THStack *hsQ2LArCC0Pi = new THStack("hsQ2LArCC0Pi", "Q^{2} in LAr ND, CC0#pi; Q^{2}; Events");
    THStack *hsQ2LArCC1Pi = new THStack("hsQ2LArCC1Pi", "Q^{2} in LAr ND, CC1#pi; Q^{2}; Events");
    THStack *hsQ2LArCC2Pi = new THStack("hsQ2LArCC2Pi", "Q^{2} in LAr ND, CC2#pi; Q^{2}; Events");
    THStack *hsQ2LArCC3Pi = new THStack("hsQ2LArCC3Pi", "Q^{2} in LAr ND, CC3#pi; Q^{2}; Events");
    hQ2LArCC0PiNuWro->SetLineColor(kRed);
    hQ2LArCC1PiNuWro->SetLineColor(kRed);
    hQ2LArCC2PiNuWro->SetLineColor(kRed);
    hQ2LArCC3PiNuWro->SetLineColor(kRed);
    hsQ2LArCC0Pi->Add(hQ2LArCC0Pi);
    hsQ2LArCC0Pi->Add(hQ2LArCC0PiNuWro);
    hsQ2LArCC1Pi->Add(hQ2LArCC1Pi);
    hsQ2LArCC1Pi->Add(hQ2LArCC1PiNuWro);
    hsQ2LArCC2Pi->Add(hQ2LArCC2Pi);
    hsQ2LArCC2Pi->Add(hQ2LArCC2PiNuWro);
    hsQ2LArCC3Pi->Add(hQ2LArCC3Pi);
    hsQ2LArCC3Pi->Add(hQ2LArCC3PiNuWro);
    hsQ2LArCC0Pi->Write("hsQ2LArCC0Pi");
    hsQ2LArCC1Pi->Write("hsQ2LArCC1Pi");
    hsQ2LArCC2Pi->Write("hsQ2LArCC2Pi");
    hsQ2LArCC3Pi->Write("hsQ2LArCC3Pi");

    THStack *hsQ2LAr = new THStack("hsQ2LAr", "LAr CC #nu_{#mu} events for various #pi multiplicities; Q^{2}; Events");
    hQ2LArCC0Pi->SetLineColor(kBlack);
    hQ2LArCC0Pi->SetFillColor(kBlack);
    hQ2LArCC1Pi->SetLineColor(kBlue);
    hQ2LArCC1Pi->SetFillColor(kBlue);
    hQ2LArCC2Pi->SetLineColor(kRed);
    hQ2LArCC2Pi->SetFillColor(kRed);
    hQ2LArCC3Pi->SetLineColor(kGreen+2);
    hQ2LArCC3Pi->SetFillColor(kGreen+2);
    hQ2LArCCHiPi->SetLineColor(kMagenta);
    hQ2LArCCHiPi->SetFillColor(kMagenta);
    hsQ2LAr->Add(hQ2LArCC0Pi);
    hsQ2LAr->Add(hQ2LArCC1Pi);
    hsQ2LAr->Add(hQ2LArCC2Pi);
    hsQ2LAr->Add(hQ2LArCC3Pi);
    hsQ2LAr->Add(hQ2LArCCHiPi);
    THStack *hsQ2LArRat = new THStack("hsQ2LArRat", "LAr ND CC #nu_{#mu} events; Q2; Fraction");
    hQ2LArCC0Pi->Divide(hQ2LAr);
    hQ2LArCC1Pi->Divide(hQ2LAr);
    hQ2LArCC2Pi->Divide(hQ2LAr);
    hQ2LArCC3Pi->Divide(hQ2LAr);
    hQ2LArCCHiPi->Divide(hQ2LAr);
    hsQ2LArRat->Add(hQ2LArCC0Pi);
    hsQ2LArRat->Add(hQ2LArCC1Pi);
    hsQ2LArRat->Add(hQ2LArCC2Pi);
    hsQ2LArRat->Add(hQ2LArCC3Pi); 
    hsQ2LArRat->Add(hQ2LArCCHiPi);
    TLegend *leg1 = new TLegend(0.5, 0.5, 0.7, 0.8);
    leg1->AddEntry(hQ2LArCC0Pi, "0#pi", "f");
    leg1->AddEntry(hQ2LArCC1Pi, "1#pi", "f");
    leg1->AddEntry(hQ2LArCC2Pi, "2#pi", "f");
    leg1->AddEntry(hQ2LArCC3Pi, "3#pi", "f");
    leg1->AddEntry(hQ2LArCCHiPi, ">3#pi", "f");
    leg1->Write("leg1"); 

    THStack *hsQ2LArNuWro = new THStack("hsQ2LArNuWro", "LAr CC #nu_{#mu} events for various #pi multiplicities (NuWro shifts); Q^{2}; Events");
    hQ2LArCC0PiNuWro->SetLineColor(kBlack);
    hQ2LArCC0PiNuWro->SetFillColor(kBlack);
    hQ2LArCC1PiNuWro->SetLineColor(kBlue);
    hQ2LArCC1PiNuWro->SetFillColor(kBlue);
    hQ2LArCC2PiNuWro->SetLineColor(kRed);
    hQ2LArCC2PiNuWro->SetFillColor(kRed);
    hQ2LArCC3PiNuWro->SetLineColor(kGreen+2);
    hQ2LArCC3PiNuWro->SetFillColor(kGreen+2);
    hQ2LArCCHiPiNuWro->SetLineColor(kMagenta);
    hQ2LArCCHiPiNuWro->SetFillColor(kMagenta);
    hsQ2LArNuWro->Add(hQ2LArCC0PiNuWro);
    hsQ2LArNuWro->Add(hQ2LArCC1PiNuWro);
    hsQ2LArNuWro->Add(hQ2LArCC2PiNuWro);
    hsQ2LArNuWro->Add(hQ2LArCC3PiNuWro);
    hsQ2LArNuWro->Add(hQ2LArCCHiPiNuWro);
    THStack *hsQ2LArNuWroRat = new THStack("hsQ2LArNuWroRat", "LAr ND CC #nu_{#mu} events (NuWro shifts); Q2; Fraction");
    hQ2LArCC0PiNuWro->Divide(hQ2LArNuWro);
    hQ2LArCC1PiNuWro->Divide(hQ2LArNuWro);
    hQ2LArCC2PiNuWro->Divide(hQ2LArNuWro);
    hQ2LArCC3PiNuWro->Divide(hQ2LArNuWro);
    hQ2LArCCHiPiNuWro->Divide(hQ2LArNuWro);
    hsQ2LArNuWroRat->Add(hQ2LArCC0PiNuWro);
    hsQ2LArNuWroRat->Add(hQ2LArCC1PiNuWro);
    hsQ2LArNuWroRat->Add(hQ2LArCC2PiNuWro);
    hsQ2LArNuWroRat->Add(hQ2LArCC3PiNuWro);
    hsQ2LArNuWroRat->Add(hQ2LArCCHiPiNuWro);

    THStack *hsWLAr = new THStack("hsWLAr", "LAr ND CC #nu_{#mu} events for various #pi multiplicities; W; Events");
    hWLArCC0Pi->SetLineColor(kBlack);
    hWLArCC0Pi->SetFillColor(kBlack);
    hWLArCC1Pi->SetLineColor(kBlue);
    hWLArCC1Pi->SetFillColor(kBlue);
    hWLArCC2Pi->SetLineColor(kRed);
    hWLArCC2Pi->SetFillColor(kRed);
    hWLArCC3Pi->SetLineColor(kGreen+2);
    hWLArCC3Pi->SetFillColor(kGreen+2);
    hWLArCCHiPi->SetLineColor(kMagenta);
    hWLArCCHiPi->SetFillColor(kMagenta);
    hsWLAr->Add(hWLArCC0Pi);
    hsWLAr->Add(hWLArCC1Pi);
    hsWLAr->Add(hWLArCC2Pi);
    hsWLAr->Add(hWLArCC3Pi);
    hsWLAr->Add(hWLArCCHiPi);
    THStack *hsWLArRat = new THStack("hsWLArRat", "LAr ND CC #nu_{#mu} events; W; Fraction");
    hWLArCC0Pi->Divide(hWLAr);
    hWLArCC1Pi->Divide(hWLAr);
    hWLArCC2Pi->Divide(hWLAr);
    hWLArCC3Pi->Divide(hWLAr);
    hWLArCCHiPi->Divide(hWLAr);
    hsWLArRat->Add(hWLArCC0Pi);
    hsWLArRat->Add(hWLArCC1Pi);
    hsWLArRat->Add(hWLArCC2Pi);
    hsWLArRat->Add(hWLArCC3Pi);
    hsWLArRat->Add(hWLArCCHiPi);

    THStack *hsWLArNuWro = new THStack("hsWLArNuWro", "LAr ND CC #nu_{#mu} events for various #pi multiplicities (NuWro shifts); W; Events");
    hWLArCC0PiNuWro->SetLineColor(kBlack);
    hWLArCC0PiNuWro->SetFillColor(kBlack);
    hWLArCC1PiNuWro->SetLineColor(kBlue);
    hWLArCC1PiNuWro->SetFillColor(kBlue);
    hWLArCC2PiNuWro->SetLineColor(kRed);
    hWLArCC2PiNuWro->SetFillColor(kRed);
    hWLArCC3PiNuWro->SetLineColor(kGreen+2);
    hWLArCC3PiNuWro->SetFillColor(kGreen+2);
    hWLArCCHiPiNuWro->SetLineColor(kMagenta);
    hWLArCCHiPiNuWro->SetFillColor(kMagenta);
    hsWLArNuWro->Add(hWLArCC0PiNuWro);
    hsWLArNuWro->Add(hWLArCC1PiNuWro);
    hsWLArNuWro->Add(hWLArCC2PiNuWro);
    hsWLArNuWro->Add(hWLArCC3PiNuWro);
    hsWLArNuWro->Add(hWLArCCHiPiNuWro);
    THStack *hsWLArNuWroRat = new THStack("hsWLArNuWroRat", "LAr ND CC #nu_{#mu} events (NuWro shifts); W; Fraction");
    hWLArCC0PiNuWro->Divide(hWLArNuWro);
    hWLArCC1PiNuWro->Divide(hWLArNuWro);
    hWLArCC2PiNuWro->Divide(hWLArNuWro);
    hWLArCC3PiNuWro->Divide(hWLArNuWro);
    hWLArCCHiPiNuWro->Divide(hWLArNuWro);
    hsWLArNuWroRat->Add(hWLArCC0PiNuWro);
    hsWLArNuWroRat->Add(hWLArCC1PiNuWro);
    hsWLArNuWroRat->Add(hWLArCC2PiNuWro);
    hsWLArNuWroRat->Add(hWLArCC3PiNuWro);
    hsWLArNuWroRat->Add(hWLArCCHiPiNuWro);

    hsWLAr->Write();
    hsWLArNuWro->Write();
    hsQ2LAr->Write();
    hsQ2LArNuWro->Write();
    hsWLArRat->Write();
    hsWLArNuWroRat->Write();
    hsQ2LArRat->Write();
    hsQ2LArNuWroRat->Write();

    // Now divide GENIE by NuWro for each of these samples
    hWLArCC0Pi->Divide(hWLArCC0PiNuWro);
    hWLArCC1Pi->Divide(hWLArCC1PiNuWro);
    hWLArCC2Pi->Divide(hWLArCC2PiNuWro);
    hWLArCC3Pi->Divide(hWLArCC3PiNuWro);
    hWLArCCHiPi->Divide(hWLArCCHiPiNuWro);
    hWLArCC0Pi->SetFillStyle(0);
    hWLArCC1Pi->SetFillStyle(0);
    hWLArCC2Pi->SetFillStyle(0);
    hWLArCC3Pi->SetFillStyle(0);
    hWLArCCHiPi->SetFillStyle(0);
    THStack *hsWLArGenieNuWroRatios = new THStack("hsWLArGenieNuWroRatios", "GENIE/NuWro for W in LAr; W; GENIE/NuWro");
    hsWLArGenieNuWroRatios->Add(hWLArCC0Pi);
    hsWLArGenieNuWroRatios->Add(hWLArCC1Pi);
    hsWLArGenieNuWroRatios->Add(hWLArCC2Pi);
    hsWLArGenieNuWroRatios->Add(hWLArCC3Pi);
    hsWLArGenieNuWroRatios->Add(hWLArCCHiPi);
    hsWLArGenieNuWroRatios->Write();
    hWLArCC0Pi->Write("hWLArCC0PiRatio");
    hWLArCC1Pi->Write("hWLArCC1PiRatio");
    hWLArCC2Pi->Write("hWLArCC2PiRatio");
    hWLArCC3Pi->Write("hWLArCC3PiRatio");
    hWLArCCHiPi->Write("hWLArCCHiPiRatio");

    hQ2LArCC0Pi->Divide(hQ2LArCC0PiNuWro);
    hQ2LArCC1Pi->Divide(hQ2LArCC1PiNuWro);
    hQ2LArCC2Pi->Divide(hQ2LArCC2PiNuWro);
    hQ2LArCC3Pi->Divide(hQ2LArCC3PiNuWro);
    hQ2LArCCHiPi->Divide(hQ2LArCCHiPiNuWro);
    hQ2LArCC0Pi->SetFillStyle(0);
    hQ2LArCC1Pi->SetFillStyle(0);
    hQ2LArCC2Pi->SetFillStyle(0);
    hQ2LArCC3Pi->SetFillStyle(0);
    hQ2LArCCHiPi->SetFillStyle(0);
    THStack *hsQ2LArGenieNuWroRatios = new THStack("hsQ2LArGenieNuWroRatios", "GENIE/NuWro for Q2 in LAr; Q2; GENIE/NuWro");
    hsQ2LArGenieNuWroRatios->Add(hQ2LArCC0Pi);
    hsQ2LArGenieNuWroRatios->Add(hQ2LArCC1Pi);
    hsQ2LArGenieNuWroRatios->Add(hQ2LArCC2Pi);
    hsQ2LArGenieNuWroRatios->Add(hQ2LArCC3Pi);
    hsQ2LArGenieNuWroRatios->Add(hQ2LArCCHiPi);
    hsQ2LArGenieNuWroRatios->Write();
    hQ2LArCC0Pi->Write("hQ2LArCC0PiRatio");
    hQ2LArCC1Pi->Write("hQ2LArCC1PiRatio");
    hQ2LArCC2Pi->Write("hQ2LArCC2PiRatio");
    hQ2LArCC3Pi->Write("hQ2LArCC3PiRatio");
    hQ2LArCCHiPi->Write("hQ2LArCCHiPiRatio");

    // FD plots
    TH1 *hNumuFHCreco = predNumuFHCreco.Predict(this_calc).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNumuFHCreco);
    TH1 *hNumuRHCreco = predNumuRHCreco.Predict(this_calc).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNumuRHCreco);
    TH1 *hNueFHCreco = predNueFHCreco.Predict(this_calc).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNueFHCreco);
    TH1 *hNueRHCreco = predNueRHCreco.Predict(this_calc).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNueRHCreco);
    hNumuFHCreco->Write("hNumuFHCreco");
    hNumuRHCreco->Write("hNumuRHCreco");
    hNueFHCreco->Write("hNueFHCreco");
    hNueRHCreco->Write("hNueRHCreco");
    // NuWro
    TH1 *hNumuFHCrecoNuwro = predNumuFHCreco.PredictSyst(this_calc, SystShifts(systlist.at(0), 1)).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNumuFHCrecoNuwro);
    TH1 *hNumuRHCrecoNuwro = predNumuRHCreco.PredictSyst(this_calc, SystShifts(systlist.at(0), 1)).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNumuRHCrecoNuwro);
    TH1 *hNueFHCrecoNuwro = predNueFHCreco.PredictSyst(this_calc, SystShifts(systlist.at(0), 1)).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNueFHCrecoNuwro);
    TH1 *hNueRHCrecoNuwro = predNueRHCreco.PredictSyst(this_calc, SystShifts(systlist.at(0), 1)).MockData(pot_fd).ToTH1(pot_fd);
    setHistAttr(hNueRHCrecoNuwro);
    hNumuFHCrecoNuwro->Write("hNumuFHCrecoNuwro");
    hNumuRHCrecoNuwro->Write("hNumuRHCrecoNuwro");
    hNueFHCrecoNuwro->Write("hNueFHCrecoNuwro");
    hNueRHCrecoNuwro->Write("hNueRHCrecoNuwro");
    hNumuFHCrecoNuwro->SetLineColor(kRed);
    hNumuRHCrecoNuwro->SetLineColor(kRed);
    hNueFHCrecoNuwro->SetLineColor(kRed);
    hNueRHCrecoNuwro->SetLineColor(kRed);
    THStack *hsNumuFHCreco = new THStack("hsNumuFHCreco", "#nu_{#mu} FHC; E_{#nu, reco} (GeV); Events / GeV");
    hNumuFHCreco->Scale(1.,"width");
    hNumuFHCrecoNuwro->Scale(1.,"width");
    hsNumuFHCreco->Add(hNumuFHCreco);
    hsNumuFHCreco->Add(hNumuFHCrecoNuwro);
    hsNumuFHCreco->Write();
    THStack *hsNumuRHCreco = new THStack("hsNumuRHCreco", "#nu_{#mu} RHC; E_{#nu, reco} (GeV); Events / GeV");
    hNumuRHCreco->Scale(1.,"width");
    hNumuRHCrecoNuwro->Scale(1.,"width");
    hsNumuRHCreco->Add(hNumuRHCreco);
    hsNumuRHCreco->Add(hNumuRHCrecoNuwro);
    hsNumuRHCreco->Write();
    THStack *hsNueFHCreco = new THStack("hsNueFHCreco", "#nu_{e} FHC; E_{#nu, reco} (GeV); Events / GeV");
    hNueFHCreco->Scale(1.,"width");
    hNueFHCrecoNuwro->Scale(1.,"width");
    hsNueFHCreco->Add(hNueFHCreco);
    hsNueFHCreco->Add(hNueFHCrecoNuwro);
    hsNueFHCreco->Write();
    THStack *hsNueRHCreco = new THStack("hsNueRHCreco", "#nu_{e} RHC; E_{#nu, reco} (GeV); Events / GeV");
    hNueRHCreco->Scale(1.,"width");
    hNueRHCrecoNuwro->Scale(1.,"width");
    hsNueRHCreco->Add(hNueRHCreco);
    hsNueRHCreco->Add(hNueRHCrecoNuwro);
    hsNueRHCreco->Write();
  }

  // GAr FHC
  TH2 *h2GAr           = predGAr.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GAr); 
  TH2 *h2GArNuWro      = predGAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArNuWro); 
  TH2 *h2GArCC0Pi      = predGArCC0Pi.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GArCC0Pi); 
  TH2 *h2GArCC0PiNuWro = predGArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArCC0PiNuWro); 
  TH2 *h2GArCC1Pi      = predGArCC1Pi.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GArCC1Pi); 
  TH2 *h2GArCC1PiNuWro = predGArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArCC1PiNuWro); 
  TH2 *h2GArCC2Pi      = predGArCC2Pi.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GArCC2Pi); 
  TH2 *h2GArCC2PiNuWro = predGArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArCC2PiNuWro); 
  TH2 *h2GArCC3Pi      = predGArCC3Pi.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GArCC3Pi); 
  TH2 *h2GArCC3PiNuWro = predGArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArCC3PiNuWro); 
  TH2 *h2GArCCHiPi      = predGArCCHiPi.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2GArCCHiPi); 
  TH2 *h2GArCCHiPiNuWro = predGArCCHiPi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  setHistAttr(h2GArCCHiPiNuWro); 

  TH1 *hGAr1d           = predGAr1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGAr1d);
  TH1 *hGArNuWro1d      = predGAr1d.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  setHistAttr(hGArNuWro1d);
  TH1 *hGArCC0Pi1d      = predGArCC0Pi1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGArCC0Pi1d);
  TH1 *hGArCC0PiNuWro1d = predGArCC0Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
  setHistAttr(hGArCC0PiNuWro1d);
  TH1 *hGArCC1Pi1d      = predGArCC1Pi1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGArCC1Pi1d);
  TH1 *hGArCC1PiNuWro1d = predGArCC1Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
  setHistAttr(hGArCC1PiNuWro1d);
  TH1 *hGArCC2Pi1d      = predGArCC2Pi1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGArCC2Pi1d);
  TH1 *hGArCC2PiNuWro1d = predGArCC2Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
  setHistAttr(hGArCC2PiNuWro1d);
  TH1 *hGArCC3Pi1d      = predGArCC3Pi1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGArCC3Pi1d);
  TH1 *hGArCC3PiNuWro1d = predGArCC3Pi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
  setHistAttr(hGArCC3PiNuWro1d);
  TH1 *hGArCCHiPi1d      = predGArCCHiPi1d.Predict(0).ToTH1(pot_nd);
  setHistAttr(hGArCCHiPi1d);
  TH1 *hGArCCHiPiNuWro1d = predGArCCHiPi1d.PredictSyst(0,SystShifts(systlist.at(0),1)).ToTH1(pot_nd);
  setHistAttr(hGArCCHiPiNuWro1d);

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
  
  // GAr RHC
  TH2 *h2GArRHC           = predGArRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArNuWroRHC      = predGArRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC0PiRHC      = predGArCC0PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC0PiNuWroRHC = predGArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC1PiRHC      = predGArCC1PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC1PiNuWroRHC = predGArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC2PiRHC      = predGArCC2PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC2PiNuWroRHC = predGArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC3PiRHC      = predGArCC3PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC3PiNuWroRHC = predGArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);

  TH1 *hWGArCC0PiRHC      = predWGArCC0PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC0PiNuWroRHC = predWGArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC1PiRHC      = predWGArCC1PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC1PiNuWroRHC = predWGArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC2PiRHC      = predWGArCC2PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC2PiNuWroRHC = predWGArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC3PiRHC      = predWGArCC3PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC3PiNuWroRHC = predWGArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  TH1 *hQ2GArCC0PiRHC      = predQ2GArCC0PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC0PiNuWroRHC = predQ2GArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC1PiRHC      = predQ2GArCC1PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC1PiNuWroRHC = predQ2GArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC2PiRHC      = predQ2GArCC2PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC2PiNuWroRHC = predQ2GArCC2PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC3PiRHC      = predQ2GArCC3PiRHC.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC3PiNuWroRHC = predQ2GArCC3PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  fout->cd();
  // GAr FHC  
  h2GAr->Write("h2GAr");
  h2GArNuWro->Write("h2GArNuWro");
  h2GArCC0Pi->Write("h2GArCC0Pi");
  h2GArCC0PiNuWro->Write("h2GArCC0PiNuWro");
  h2GArCC1Pi->Write("h2GArCC1Pi");
  h2GArCC1PiNuWro->Write("h2GArCC1PiNuWro");
  h2GArCC2Pi->Write("h2GArCC2Pi");
  h2GArCC2PiNuWro->Write("h2GArCC2PiNuWro");
  h2GArCC3Pi->Write("h2GArCC3Pi");
  h2GArCC3PiNuWro->Write("h2GArCC3PiNuWro");
  h2GArCCHiPi->Write("h2GArCCHiPi");
  h2GArCCHiPiNuWro->Write("h2GArCCHiPiNuWro");

  hGAr1d->Write("hGAr1d");
  hGArNuWro1d->Write("hGArNuWro1d");
  hGArCC0Pi1d->Write("hGArCC0Pi1d");
  hGArCC0PiNuWro1d->Write("hGArCC0PiNuWro1d");
  hGArCC1Pi1d->Write("hGArCC1Pi1d");
  hGArCC1PiNuWro1d->Write("hGArCC1PiNuWro1d");
  hGArCC2Pi1d->Write("hGArCC2Pi1d");
  hGArCC2PiNuWro1d->Write("hGArCC2PiNuWro1d");
  hGArCC3Pi1d->Write("hGArCC3Pi1d");
  hGArCC3PiNuWro1d->Write("hGArCC3PiNuWro1d");
  hGArCCHiPi1d->Write("hGArCCHiPi1d");
  hGArCCHiPiNuWro1d->Write("hGArCCHiPiNuWro1d");

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
  // GAr RHC
  h2GArRHC->Write("h2GArRHC");
  h2GArNuWroRHC->Write("h2GArNuWroRHC");
  h2GArCC0PiRHC->Write("h2GArCC0PiRHC");
  h2GArCC0PiNuWroRHC->Write("h2GArCC0PiNuWroRHC");
  h2GArCC1PiRHC->Write("h2GArCC1PiRHC");
  h2GArCC1PiNuWroRHC->Write("h2GArCC1PiNuWroRHC");
  h2GArCC2PiRHC->Write("h2GArCC2PiRHC");
  h2GArCC2PiNuWroRHC->Write("h2GArCC2PiNuWroRHC");
  h2GArCC3PiRHC->Write("h2GArCC3PiRHC");
  h2GArCC3PiNuWroRHC->Write("h2GArCC3PiNuWroRHC");

  hWGArCC0PiRHC->Write("hWGArCC0PiRHC");
  hWGArCC0PiNuWroRHC->Write("hWGArCC0PiNuWroRHC");
  hWGArCC1PiRHC->Write("hWGArCC1PiRHC");
  hWGArCC1PiNuWroRHC->Write("hWGArCC1PiNuWroRHC");
  hWGArCC2PiRHC->Write("hWGArCC2PiRHC");
  hWGArCC2PiNuWroRHC->Write("hWGArCC2PiNuWroRHC");
  hWGArCC3PiRHC->Write("hWGArCC3PiRHC");
  hWGArCC3PiNuWroRHC->Write("hWGArCC3PiNuWroRHC");

  hQ2GArCC0PiRHC->Write("hQ2GArCC0PiRHC");
  hQ2GArCC0PiNuWroRHC->Write("hQ2GArCC0PiNuWroRHC");
  hQ2GArCC1PiRHC->Write("hQ2GArCC1PiRHC");
  hQ2GArCC1PiNuWroRHC->Write("hQ2GArCC1PiNuWroRHC");
  hQ2GArCC2PiRHC->Write("hQ2GArCC2PiRHC");
  hQ2GArCC2PiNuWroRHC->Write("hQ2GArCC2PiNuWroRHC");
  hQ2GArCC3PiRHC->Write("hQ2GArCC3PiRHC");
  hQ2GArCC3PiNuWroRHC->Write("hQ2GArCC3PiNuWroRHC");

  THStack *hsWGArCC0Pi = new THStack("hsWGArCC0Pi", "W in GAr ND, CC0#pi; W; Events");
  THStack *hsWGArCC1Pi = new THStack("hsWGArCC1Pi", "W in GAr ND, CC1#pi; W; Events");
  THStack *hsWGArCC2Pi = new THStack("hsWGArCC2Pi", "W in GAr ND, CC2#pi; W; Events");
  THStack *hsWGArCC3Pi = new THStack("hsWGArCC3Pi", "W in GAr ND, CC3#pi; W; Events");
  hWGArCC0PiNuWro->SetLineColor(kRed);
  hWGArCC1PiNuWro->SetLineColor(kRed);
  hWGArCC2PiNuWro->SetLineColor(kRed);
  hWGArCC3PiNuWro->SetLineColor(kRed);
  hsWGArCC0Pi->Add(hWGArCC0Pi);
  hsWGArCC0Pi->Add(hWGArCC0PiNuWro);
  hsWGArCC1Pi->Add(hWGArCC1Pi);
  hsWGArCC1Pi->Add(hWGArCC1PiNuWro);
  hsWGArCC2Pi->Add(hWGArCC2Pi);
  hsWGArCC2Pi->Add(hWGArCC2PiNuWro);
  hsWGArCC3Pi->Add(hWGArCC3Pi);
  hsWGArCC3Pi->Add(hWGArCC3PiNuWro);
  hsWGArCC0Pi->Write("hsWGArCC0Pi");
  hsWGArCC1Pi->Write("hsWGArCC1Pi");
  hsWGArCC2Pi->Write("hsWGArCC2Pi");
  hsWGArCC3Pi->Write("hsWGArCC3Pi");

  THStack *hsQ2GArCC0Pi = new THStack("hsQ2GArCC0Pi", "Q^{2} in GAr ND, CC0#pi; Q^{2}; Events");
  THStack *hsQ2GArCC1Pi = new THStack("hsQ2GArCC1Pi", "Q^{2} in GAr ND, CC1#pi; Q^{2}; Events");
  THStack *hsQ2GArCC2Pi = new THStack("hsQ2GArCC2Pi", "Q^{2} in GAr ND, CC2#pi; Q^{2}; Events");
  THStack *hsQ2GArCC3Pi = new THStack("hsQ2GArCC3Pi", "Q^{2} in GAr ND, CC3#pi; Q^{2}; Events");
  hQ2GArCC0PiNuWro->SetLineColor(kRed);
  hQ2GArCC1PiNuWro->SetLineColor(kRed);
  hQ2GArCC2PiNuWro->SetLineColor(kRed);
  hQ2GArCC3PiNuWro->SetLineColor(kRed);
  hsQ2GArCC0Pi->Add(hQ2GArCC0Pi);
  hsQ2GArCC0Pi->Add(hQ2GArCC0PiNuWro);
  hsQ2GArCC1Pi->Add(hQ2GArCC1Pi);
  hsQ2GArCC1Pi->Add(hQ2GArCC1PiNuWro);
  hsQ2GArCC2Pi->Add(hQ2GArCC2Pi);
  hsQ2GArCC2Pi->Add(hQ2GArCC2PiNuWro);
  hsQ2GArCC3Pi->Add(hQ2GArCC3Pi);
  hsQ2GArCC3Pi->Add(hQ2GArCC3PiNuWro);
  hsQ2GArCC0Pi->Write("hsQ2GArCC0Pi");
  hsQ2GArCC1Pi->Write("hsQ2GArCC1Pi");
  hsQ2GArCC2Pi->Write("hsQ2GArCC2Pi");
  hsQ2GArCC3Pi->Write("hsQ2GArCC3Pi");

  THStack *hsQ2GAr = new THStack("hsQ2GAr", "GAr CC #nu_{#mu} events for various #pi multiplicities; Q^{2}; Events");
  hQ2GArCC0Pi->SetLineColor(kBlack);
  hQ2GArCC0Pi->SetFillColor(kBlack);
  hQ2GArCC1Pi->SetLineColor(kBlue);
  hQ2GArCC1Pi->SetFillColor(kBlue);
  hQ2GArCC2Pi->SetLineColor(kRed);
  hQ2GArCC2Pi->SetFillColor(kRed);
  hQ2GArCC3Pi->SetLineColor(kGreen+2);
  hQ2GArCC3Pi->SetFillColor(kGreen+2);
  hQ2GArCCHiPi->SetLineColor(kMagenta);
  hQ2GArCCHiPi->SetFillColor(kMagenta);
  hsQ2GAr->Add(hQ2GArCC0Pi);
  hsQ2GAr->Add(hQ2GArCC1Pi);
  hsQ2GAr->Add(hQ2GArCC2Pi);
  hsQ2GAr->Add(hQ2GArCC3Pi);
  hsQ2GAr->Add(hQ2GArCCHiPi);
  hsQ2GAr->Write();
  THStack *hsQ2GArRat = new THStack("hsQ2GArRat", "GAr ND CC #nu_{#mu} events; Q2; Fraction");
  hQ2GArCC0Pi->Divide(hQ2GAr);
  hQ2GArCC1Pi->Divide(hQ2GAr);
  hQ2GArCC2Pi->Divide(hQ2GAr);
  hQ2GArCC3Pi->Divide(hQ2GAr);
  hQ2GArCCHiPi->Divide(hQ2GAr);
  hsQ2GArRat->Add(hQ2GArCC0Pi);
  hsQ2GArRat->Add(hQ2GArCC1Pi);
  hsQ2GArRat->Add(hQ2GArCC2Pi);
  hsQ2GArRat->Add(hQ2GArCC3Pi); 
  hsQ2GArRat->Add(hQ2GArCCHiPi); 
  hsQ2GArRat->Write();

  THStack *hsQ2GArNuWro = new THStack("hsQ2GArNuWro", "GAr CC #nu_{#mu} events for various #pi multiplicities (NuWro shifts); Q^{2}; Events");
  hQ2GArCC0PiNuWro->SetLineColor(kBlack);
  hQ2GArCC0PiNuWro->SetFillColor(kBlack);
  hQ2GArCC1PiNuWro->SetLineColor(kBlue);
  hQ2GArCC1PiNuWro->SetFillColor(kBlue);
  hQ2GArCC2PiNuWro->SetLineColor(kRed);
  hQ2GArCC2PiNuWro->SetFillColor(kRed);
  hQ2GArCC3PiNuWro->SetLineColor(kGreen+2);
  hQ2GArCC3PiNuWro->SetFillColor(kGreen+2);
  hQ2GArCCHiPiNuWro->SetLineColor(kMagenta);
  hQ2GArCCHiPiNuWro->SetFillColor(kMagenta);
  hQ2GArCC0PiNuWro->SetTitle("Fraction of CC0#pi events");
  hQ2GArCC1PiNuWro->SetTitle("Fraction of CC1#pi events");
  hQ2GArCC2PiNuWro->SetTitle("Fraction of CC2#pi events");
  hQ2GArCC3PiNuWro->SetTitle("Fraction of CC3#pi events");
  hQ2GArCCHiPiNuWro->SetTitle("Fraction of CC events with >3#pi");
  hsQ2GArNuWro->Add(hQ2GArCC0PiNuWro);
  hsQ2GArNuWro->Add(hQ2GArCC1PiNuWro);
  hsQ2GArNuWro->Add(hQ2GArCC2PiNuWro);
  hsQ2GArNuWro->Add(hQ2GArCC3PiNuWro);
  hsQ2GArNuWro->Add(hQ2GArCCHiPiNuWro);
  hsQ2GArNuWro->Write();
  THStack *hsQ2GArNuWroRat = new THStack("hsQ2GArNuWroRat", "GAr ND CC #nu_{#mu} events (NuWro shifts); Q2; Fraction");
  hQ2GArCC0PiNuWro->Divide(hQ2GArNuWro);
  hQ2GArCC1PiNuWro->Divide(hQ2GArNuWro);
  hQ2GArCC2PiNuWro->Divide(hQ2GArNuWro);
  hQ2GArCC3PiNuWro->Divide(hQ2GArNuWro);
  hQ2GArCCHiPiNuWro->Divide(hQ2GArNuWro);
  hsQ2GArNuWroRat->Add(hQ2GArCC0PiNuWro);
  hsQ2GArNuWroRat->Add(hQ2GArCC1PiNuWro);
  hsQ2GArNuWroRat->Add(hQ2GArCC2PiNuWro);
  hsQ2GArNuWroRat->Add(hQ2GArCC3PiNuWro);
  hsQ2GArNuWroRat->Add(hQ2GArCCHiPiNuWro);
  hsQ2GArNuWroRat->Write();

  THStack *hsWGAr = new THStack("hsWGAr", "GAr ND CC #nu_{#mu} events for various #pi multiplicities; W; Events");
  hWGArCC0Pi->SetLineColor(kBlack);
  hWGArCC0Pi->SetFillColor(kBlack);
  hWGArCC1Pi->SetLineColor(kBlue);
  hWGArCC1Pi->SetFillColor(kBlue);
  hWGArCC2Pi->SetLineColor(kRed);
  hWGArCC2Pi->SetFillColor(kRed);
  hWGArCC3Pi->SetLineColor(kGreen+2);
  hWGArCC3Pi->SetFillColor(kGreen+2);
  hWGArCCHiPi->SetLineColor(kMagenta);
  hWGArCCHiPi->SetFillColor(kMagenta);
  hsWGAr->Add(hWGArCC0Pi);
  hsWGAr->Add(hWGArCC1Pi);
  hsWGAr->Add(hWGArCC2Pi);
  hsWGAr->Add(hWGArCC3Pi);
  hsWGAr->Add(hWGArCCHiPi);
  hsWGAr->Write();
  THStack *hsWGArRat = new THStack("hsWGArRat", "GAr ND CC #nu_{#mu} events; W; Fraction");
  hWGArCC0Pi->Divide(hWGAr);
  hWGArCC1Pi->Divide(hWGAr);
  hWGArCC2Pi->Divide(hWGAr);
  hWGArCC3Pi->Divide(hWGAr);
  hWGArCCHiPi->Divide(hWGAr);
  hsWGArRat->Add(hWGArCC0Pi);
  hsWGArRat->Add(hWGArCC1Pi);
  hsWGArRat->Add(hWGArCC2Pi);
  hsWGArRat->Add(hWGArCC3Pi);
  hsWGArRat->Add(hWGArCCHiPi);
  hsWGArRat->Write();

  THStack *hsWGArNuWro = new THStack("hsWGArNuWro", "GAr ND CC #nu_{#mu} events for various #pi multiplicities (NuWro shifts); W; Events");
  hWGArCC0PiNuWro->SetLineColor(kBlack);
  hWGArCC0PiNuWro->SetFillColor(kBlack);
  hWGArCC1PiNuWro->SetLineColor(kBlue);
  hWGArCC1PiNuWro->SetFillColor(kBlue);
  hWGArCC2PiNuWro->SetLineColor(kRed);
  hWGArCC2PiNuWro->SetFillColor(kRed);
  hWGArCC3PiNuWro->SetLineColor(kGreen+2);
  hWGArCC3PiNuWro->SetFillColor(kGreen+2);
  hWGArCCHiPiNuWro->SetLineColor(kMagenta);
  hWGArCCHiPiNuWro->SetFillColor(kMagenta);
  hsWGArNuWro->Add(hWGArCC0PiNuWro);
  hsWGArNuWro->Add(hWGArCC1PiNuWro);
  hsWGArNuWro->Add(hWGArCC2PiNuWro);
  hsWGArNuWro->Add(hWGArCC3PiNuWro);
  hsWGArNuWro->Add(hWGArCCHiPiNuWro);
  hsWGArNuWro->Write();
  THStack *hsWGArNuWroRat = new THStack("hsWGArNuWroRat", "GAr ND CC #nu_{#mu} events (NuWro shifts); W; Fraction");
  hWGArCC0PiNuWro->Divide(hWGArNuWro);
  hWGArCC1PiNuWro->Divide(hWGArNuWro);
  hWGArCC2PiNuWro->Divide(hWGArNuWro);
  hWGArCC3PiNuWro->Divide(hWGArNuWro);
  hWGArCCHiPiNuWro->Divide(hWGArNuWro);
  hsWGArNuWroRat->Add(hWGArCC0PiNuWro);
  hsWGArNuWroRat->Add(hWGArCC1PiNuWro);
  hsWGArNuWroRat->Add(hWGArCC2PiNuWro);
  hsWGArNuWroRat->Add(hWGArCC3PiNuWro);
  hsWGArNuWroRat->Add(hWGArCCHiPiNuWro);
  hsWGArNuWroRat->Write();

  hWGArCC0Pi->Divide(hWGArCC0PiNuWro);
  hWGArCC1Pi->Divide(hWGArCC1PiNuWro);
  hWGArCC2Pi->Divide(hWGArCC2PiNuWro);
  hWGArCC3Pi->Divide(hWGArCC3PiNuWro);
  hWGArCCHiPi->Divide(hWGArCCHiPiNuWro);
  hWGArCC0Pi->SetFillStyle(0);
  hWGArCC1Pi->SetFillStyle(0);
  hWGArCC2Pi->SetFillStyle(0);
  hWGArCC3Pi->SetFillStyle(0);
  hWGArCCHiPi->SetFillStyle(0);
  THStack *hsWGArGenieNuWroRatios = new THStack("hsWGArGenieNuWroRatios", "GENIE/NuWro for W in GAr; W; GENIE/NuWro");
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

  hQ2GArCC0Pi->Divide(hQ2GArCC0PiNuWro);
  hQ2GArCC1Pi->Divide(hQ2GArCC1PiNuWro);
  hQ2GArCC2Pi->Divide(hQ2GArCC2PiNuWro);
  hQ2GArCC3Pi->Divide(hQ2GArCC3PiNuWro);
  hQ2GArCCHiPi->Divide(hQ2GArCCHiPiNuWro);
  hQ2GArCC0Pi->SetFillStyle(0);
  hQ2GArCC1Pi->SetFillStyle(0);
  hQ2GArCC2Pi->SetFillStyle(0);
  hQ2GArCC3Pi->SetFillStyle(0);
  hQ2GArCCHiPi->SetFillStyle(0);
  THStack *hsQ2GArGenieNuWroRatios = new THStack("hsQ2GArGenieNuWroRatios", "GENIE/NuWro for Q2 in GAr; Q2; GENIE/NuWro");
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
  
  fout->Close();
  delete fout;

} // pionMultiPlots
