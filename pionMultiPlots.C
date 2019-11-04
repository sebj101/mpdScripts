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

void pionMultiPlots(const char *outFile, 
		    const char *lardir="/pnfs/dune/persistent/users/picker24/CAFv4/", 
		    const char *gardir="/dune/data/users/sbjones/gasTpcCAF/") 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  // ND vars
  const Var kRecoEnergyND  = SIMPLEVAR(dune.Ev_reco);
  const Var kRecoYND       = (SIMPLEVAR(dune.Ev_reco)-SIMPLEVAR(dune.Elep_reco))/SIMPLEVAR(dune.Ev_reco);
  const Var kTrueYND       = (SIMPLEVAR(dune.Ev)-SIMPLEVAR(dune.LepE))/SIMPLEVAR(dune.Ev);
  const Var kTrueEnergy    = SIMPLEVAR(dune.Ev);
  const Var kTrueLepEnergy = SIMPLEVAR(dune.LepE);
  const Var kNPi = SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipi0);
  const Var kQ2  = SIMPLEVAR(dune.Q2);
  const Var kW   = SIMPLEVAR(dune.W);
  const Var kFHC = SIMPLEVAR(dune.isFHC);
  const Var Enu_reco_numu = SIMPLEVAR(dune.Ev_reco_numu);
  const Var Enu_reco_nue  = SIMPLEVAR(dune.Ev_reco_nue);
  // CV weighting
  //const Var kGENIEWeights = SIMPLEVAR(dune.total_xsSyst_cv_wgt);

  std::vector<double> binEEdges = {0., 0.75, 1.25, 1.5, 1.75, 
				   2., 2.25, 2.5, 2.75, 
				   3., 3.25, 3.5, 3.75,
				   4., 5., 6., 10.};
  std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
  const Binning binsNDEreco  = Binning::Custom(binEEdges);
  const Binning simpleBins   = Binning::Simple(20, 0, 5);
  const Binning simpleWBins  = Binning::Simple(30, 0, 6);
  const Binning simpleQ2Bins = Binning::Simple(30, 0, 6);

  const HistAxis axTrueRecoND("E_{#nu, true} (GeV)", simpleBins, kTrueEnergy,
			      "E_{#nu, reco} (GeV)", simpleBins, kRecoEnergyND);
  const HistAxis axQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kQ2);
  const HistAxis axW("W (GeV)", simpleWBins, kW);
  const HistAxis axisnumu("Reco #nu energy (GeV)", Binning::Simple(40, 0, 10), Enu_reco_numu);
  const HistAxis axisnue("Reco #nu energy (GeV)", Binning::Simple(40, 0, 10), Enu_reco_nue);
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

  //loadersFHC.Go();
  //loadersRHC.Go();
  
  NoOscPredictionGenerator genLAr(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV);
  NoOscPredictionGenerator genLArCC0Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator genLArCC1Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator genLArCC2Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator genLArCC3Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC3Pi);

  NoOscPredictionGenerator WLAr(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV);
  NoOscPredictionGenerator WLArCC0Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator WLArCC1Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator WLArCC2Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator WLArCC3Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);

  NoOscPredictionGenerator Q2LAr(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV);
  NoOscPredictionGenerator Q2LArCC0Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator Q2LArCC1Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator Q2LArCC2Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator Q2LArCC3Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsTrueFV && kIsCC3Pi);

  // FHC
  PredictionInterp predLAr(systlist, this_calc, genLAr, loadersLArFHC);
  PredictionInterp predLArCC0Pi(systlist, this_calc, genLArCC0Pi, loadersLArFHC);
  PredictionInterp predLArCC1Pi(systlist, this_calc, genLArCC1Pi, loadersLArFHC);
  PredictionInterp predLArCC2Pi(systlist, this_calc, genLArCC2Pi, loadersLArFHC);
  PredictionInterp predLArCC3Pi(systlist, this_calc, genLArCC3Pi, loadersLArFHC);

  PredictionInterp predWLAr(systlist, this_calc, WLAr, loadersLArFHC);
  PredictionInterp predWLArCC0Pi(systlist, this_calc, WLArCC0Pi, loadersLArFHC);
  PredictionInterp predWLArCC1Pi(systlist, this_calc, WLArCC1Pi, loadersLArFHC);
  PredictionInterp predWLArCC2Pi(systlist, this_calc, WLArCC2Pi, loadersLArFHC);
  PredictionInterp predWLArCC3Pi(systlist, this_calc, WLArCC3Pi, loadersLArFHC);

  PredictionInterp predQ2LAr(systlist, this_calc, Q2LAr, loadersLArFHC);
  PredictionInterp predQ2LArCC0Pi(systlist, this_calc, Q2LArCC0Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC1Pi(systlist, this_calc, Q2LArCC1Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC2Pi(systlist, this_calc, Q2LArCC2Pi, loadersLArFHC);
  PredictionInterp predQ2LArCC3Pi(systlist, this_calc, Q2LArCC3Pi, loadersLArFHC);
  // RHC
  PredictionInterp predLArRHC(systlist, this_calc, genLAr, loadersLArRHC);
  PredictionInterp predLArCC0PiRHC(systlist, this_calc, genLArCC0Pi, loadersLArRHC);
  PredictionInterp predLArCC1PiRHC(systlist, this_calc, genLArCC1Pi, loadersLArRHC);
  PredictionInterp predLArCC2PiRHC(systlist, this_calc, genLArCC2Pi, loadersLArRHC);
  PredictionInterp predLArCC3PiRHC(systlist, this_calc, genLArCC3Pi, loadersLArRHC);

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
  NoOscPredictionGenerator genGAr(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)));
  NoOscPredictionGenerator genGArCC0Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator genGArCC1Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator genGArCC2Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator genGArCC3Pi(axTrueRecoND, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC3Pi);

  NoOscPredictionGenerator WGAr(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)));
  NoOscPredictionGenerator WGArCC0Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) &&kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator WGArCC1Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC1Pi);
  NoOscPredictionGenerator WGArCC2Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC2Pi);
  NoOscPredictionGenerator WGArCC3Pi(axW, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC3Pi);

  NoOscPredictionGenerator Q2GAr(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)));
  NoOscPredictionGenerator Q2GArCC0Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC0Pi);
  NoOscPredictionGenerator Q2GArCC1Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC1Pi);
  NoOscPredictionGenerator Q2GArCC2Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC2Pi);
  NoOscPredictionGenerator Q2GArCC3Pi(axQ2, ((kFHC==1 && kPassND_FHC_NUMU) || (kFHC!=1 && kPassND_RHC_NUMU)) && kIsCC3Pi);

  // FHC
  PredictionInterp predGAr(systlist, this_calc, genGAr, loadersGArFHC);
  PredictionInterp predGArCC0Pi(systlist, this_calc, genGArCC0Pi, loadersGArFHC);
  PredictionInterp predGArCC1Pi(systlist, this_calc, genGArCC1Pi, loadersGArFHC);
  PredictionInterp predGArCC2Pi(systlist, this_calc, genGArCC2Pi, loadersGArFHC);
  PredictionInterp predGArCC3Pi(systlist, this_calc, genGArCC3Pi, loadersGArFHC);

  PredictionInterp predWGAr(systlist, this_calc, WGAr, loadersGArFHC);
  PredictionInterp predWGArCC0Pi(systlist, this_calc, WGArCC0Pi, loadersGArFHC);
  PredictionInterp predWGArCC1Pi(systlist, this_calc, WGArCC1Pi, loadersGArFHC);
  PredictionInterp predWGArCC2Pi(systlist, this_calc, WGArCC2Pi, loadersGArFHC);
  PredictionInterp predWGArCC3Pi(systlist, this_calc, WGArCC3Pi, loadersGArFHC);

  PredictionInterp predQ2GAr(systlist, this_calc, Q2GAr, loadersGArFHC);
  PredictionInterp predQ2GArCC0Pi(systlist, this_calc, Q2GArCC0Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC1Pi(systlist, this_calc, Q2GArCC1Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC2Pi(systlist, this_calc, Q2GArCC2Pi, loadersGArFHC);
  PredictionInterp predQ2GArCC3Pi(systlist, this_calc, Q2GArCC3Pi, loadersGArFHC);
  // RHC
  PredictionInterp predGArRHC(systlist, this_calc, genGAr, loadersGArRHC);
  PredictionInterp predGArCC0PiRHC(systlist, this_calc, genGArCC0Pi, loadersGArRHC);
  PredictionInterp predGArCC1PiRHC(systlist, this_calc, genGArCC1Pi, loadersGArRHC);
  PredictionInterp predGArCC2PiRHC(systlist, this_calc, genGArCC2Pi, loadersGArRHC);
  PredictionInterp predGArCC3PiRHC(systlist, this_calc, genGArCC3Pi, loadersGArRHC);

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
  loadersLArFHC.Go();
  loadersLArRHC.Go();

  std::cout<<systlist.size()<<std::endl;

  TH2 *h2LAr           = predLAr.Predict(0).ToTH2(pot_nd);
  setHistAttr(h2LAr);
  TH2 *h2LArNuWro      = predLAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC0Pi      = predLArCC0Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC0PiNuWro = predLArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC1Pi      = predLArCC1Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC1PiNuWro = predLArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC2Pi      = predLArCC2Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC2PiNuWro = predLArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC3Pi      = predLArCC3Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC3PiNuWro = predLArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);

  TH1 *hWLArCC0Pi      = predWLArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC0PiNuWro = predWLArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC1Pi      = predWLArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC1PiNuWro = predWLArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC2Pi      = predWLArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC2PiNuWro = predWLArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC3Pi      = predWLArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC3PiNuWro = predWLArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  TH1 *hQ2LArCC0Pi      = predQ2LArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2LArCC0PiNuWro = predQ2LArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2LArCC1Pi      = predQ2LArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2LArCC1PiNuWro = predQ2LArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2LArCC2Pi      = predQ2LArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2LArCC2PiNuWro = predQ2LArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2LArCC3Pi      = predQ2LArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2LArCC3PiNuWro = predQ2LArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  // LAr RHC
  TH2 *h2LArRHC           = predLArRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArNuWroRHC      = predLArRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC0PiRHC      = predLArCC0PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC0PiNuWroRHC = predLArCC0PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2LArCC1PiRHC      = predLArCC1PiRHC.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArCC1PiNuWroRHC = predLArCC1PiRHC.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
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

  TH2 *h2GAr           = predGAr.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArNuWro      = predGAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC0Pi      = predGArCC0Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC0PiNuWro = predGArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC1Pi      = predGArCC1Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC1PiNuWro = predGArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC2Pi      = predGArCC2Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC2PiNuWro = predGArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);
  TH2 *h2GArCC3Pi      = predGArCC3Pi.Predict(0).ToTH2(pot_nd);
  TH2 *h2GArCC3PiNuWro = predGArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);

  TH1 *hWGArCC0Pi      = predWGArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC0PiNuWro = predWGArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC1Pi      = predWGArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC1PiNuWro = predWGArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC2Pi      = predWGArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC2PiNuWro = predWGArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWGArCC3Pi      = predWGArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWGArCC3PiNuWro = predWGArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  TH1 *hQ2GArCC0Pi      = predQ2GArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC0PiNuWro = predQ2GArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC1Pi      = predQ2GArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC1PiNuWro = predQ2GArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC2Pi      = predQ2GArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC2PiNuWro = predQ2GArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hQ2GArCC3Pi      = predQ2GArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hQ2GArCC3PiNuWro = predQ2GArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  
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

  TFile *fout = new TFile(outFile, "recreate");
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

  hWGArCC0Pi->Write("hWGArCC0Pi");
  hWGArCC0PiNuWro->Write("hWGArCC0PiNuWro");
  hWGArCC1Pi->Write("hWGArCC1Pi");
  hWGArCC1PiNuWro->Write("hWGArCC1PiNuWro");
  hWGArCC2Pi->Write("hWGArCC2Pi");
  hWGArCC2PiNuWro->Write("hWGArCC2PiNuWro");
  hWGArCC3Pi->Write("hWGArCC3Pi");
  hWGArCC3PiNuWro->Write("hWGArCC3PiNuWro");

  hQ2GArCC0Pi->Write("hQ2GArCC0Pi");
  hQ2GArCC0PiNuWro->Write("hQ2GArCC0PiNuWro");
  hQ2GArCC1Pi->Write("hQ2GArCC1Pi");
  hQ2GArCC1PiNuWro->Write("hQ2GArCC1PiNuWro");
  hQ2GArCC2Pi->Write("hQ2GArCC2Pi");
  hQ2GArCC2PiNuWro->Write("hQ2GArCC2PiNuWro");
  hQ2GArCC3Pi->Write("hQ2GArCC3Pi");
  hQ2GArCC3PiNuWro->Write("hQ2GArCC3PiNuWro");
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

  hWLArCC0Pi->Write("hWLArCC0Pi");
  hWLArCC0PiNuWro->Write("hWLArCC0PiNuWro");
  hWLArCC1Pi->Write("hWLArCC1Pi");
  hWLArCC1PiNuWro->Write("hWLArCC1PiNuWro");
  hWLArCC2Pi->Write("hWLArCC2Pi");
  hWLArCC2PiNuWro->Write("hWLArCC2PiNuWro");
  hWLArCC3Pi->Write("hWLArCC3Pi");
  hWLArCC3PiNuWro->Write("hWLArCC3PiNuWro");

  hQ2LArCC0Pi->Write("hQ2LArCC0Pi");
  hQ2LArCC0PiNuWro->Write("hQ2LArCC0PiNuWro");
  hQ2LArCC1Pi->Write("hQ2LArCC1Pi");
  hQ2LArCC1PiNuWro->Write("hQ2LArCC1PiNuWro");
  hQ2LArCC2Pi->Write("hQ2LArCC2Pi");
  hQ2LArCC2PiNuWro->Write("hQ2LArCC2PiNuWro");
  hQ2LArCC3Pi->Write("hQ2LArCC3Pi");
  hQ2LArCC3PiNuWro->Write("hQ2LArCC3PiNuWro");
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

  THStack *hsQ2GArCC0Pi = new THStack("hsQ2GArCC0Pi", "Q^2 in GAr ND, CC0#pi; Q^2; Events");
  THStack *hsQ2GArCC1Pi = new THStack("hsQ2GArCC1Pi", "Q^2 in GAr ND, CC1#pi; Q^2; Events");
  THStack *hsQ2GArCC2Pi = new THStack("hsQ2GArCC2Pi", "Q^2 in GAr ND, CC2#pi; Q^2; Events");
  THStack *hsQ2GArCC3Pi = new THStack("hsQ2GArCC3Pi", "Q^2 in GAr ND, CC3#pi; Q^2; Events");
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

  THStack *hsQ2LArCC0Pi = new THStack("hsQ2LArCC0Pi", "Q^2 in LAr ND, CC0#pi; Q^2; Events");
  THStack *hsQ2LArCC1Pi = new THStack("hsQ2LArCC1Pi", "Q^2 in LAr ND, CC1#pi; Q^2; Events");
  THStack *hsQ2LArCC2Pi = new THStack("hsQ2LArCC2Pi", "Q^2 in LAr ND, CC2#pi; Q^2; Events");
  THStack *hsQ2LArCC3Pi = new THStack("hsQ2LArCC3Pi", "Q^2 in LAr ND, CC3#pi; Q^2; Events");
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

  fout->Close();
  delete fout;

} // pionMultiPlots
