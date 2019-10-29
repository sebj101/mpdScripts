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
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/RefineSeeds.h"
#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "StandardRecord/StandardRecord.h"
//#include "CAFAna/Systs/GenieSysts.h"
#include "CAFAna/Systs/NuWroReweightFakeData.h"
#include "CAFAna/Systs/XSecSysts.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include <tuple>
#include "Utilities/rootlogon.C"

using namespace ana;

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

// CV weighting
//const Var kGENIEWeights = SIMPLEVAR(dune.total_cv_wgt);

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
// POT for 3.5 years
const double pot_fd = 3.5 * POT120 * 40/1.13;
const double pot_nd = 3.5 * POT120;
// This is pretty annoying, but the above is for 7 years staged, which is 336 kT MW yr
const double nom_exposure = 336.;

void pionMultiPlots(const char *lardir="/pnfs/dune/persistent/users/picker24/CAFv4/", 
		    const char *gardir="/dune/data/users/sbjones/gasTpcCAF/") 
{
  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  std::vector<const ISyst*> systlist = GetXSecSysts({"NuWroReweightFakeData"});

  Loaders loadersLAr;
  Loaders loadersGAr;
  SpectrumLoader loaderLArFHC(Form("%s/ND_FHC_CAF.root", lardir), kBeam, 0);
  SpectrumLoader loaderLArRHC(Form("%s/ND_RHC_CAF.root", lardir), kBeam, 0);
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", gardir), kBeam, 0);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", gardir), kBeam, 0);
  loadersGAr.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGAr.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);
  loadersLAr.AddLoader(&loaderLArFHC, caf::kNEARDET, Loaders::kMC);
  loadersLAr.AddLoader(&loaderLArRHC, caf::kNEARDET, Loaders::kMC);

  NoOscPredictionGenerator genLAr(axTrueRecoND, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				  kIsTrueFV);
  NoOscPredictionGenerator genLArCC0Pi(axTrueRecoND, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator genLArCC1Pi(axTrueRecoND, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator genLArCC2Pi(axTrueRecoND, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator genLArCC3Pi(axTrueRecoND, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC3Pi);
  std::vector<NoOscPredictionGenerator> genVec;
  genVec.push_back(genLAr);
  genVec.push_back(genLArCC0Pi);
  genVec.push_back(genLArCC1Pi);
  genVec.push_back(genLArCC2Pi);
  genVec.push_back(genLArCC3Pi);

  NoOscPredictionGenerator WLAr(axW, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				  kIsTrueFV);
  NoOscPredictionGenerator WLArCC0Pi(axW, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator WLArCC1Pi(axW, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator WLArCC2Pi(axW, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator WLArCC3Pi(axW, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC3Pi);
  std::vector<NoOscPredictionGenerator> WVec;
  genVec.push_back(WLAr);
  genVec.push_back(WLArCC0Pi);
  genVec.push_back(WLArCC1Pi);
  genVec.push_back(WLArCC2Pi);
  genVec.push_back(WLArCC3Pi);

  NoOscPredictionGenerator Q2LAr(axQ2, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				  kIsTrueFV);
  NoOscPredictionGenerator Q2LArCC0Pi(axQ2, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC0Pi);
  NoOscPredictionGenerator Q2LArCC1Pi(axQ2, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC1Pi);
  NoOscPredictionGenerator Q2LArCC2Pi(axQ2, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC2Pi);
  NoOscPredictionGenerator Q2LArCC3Pi(axQ2, ((kIsFHC && kPassND_FHC_NUMU) || (!kIsFHC && kPassND_FHC_NUMU)) && 
				       kIsTrueFV && kIsCC3Pi);
  std::vector<NoOscPredictionGenerator> Q2Vec;
  genVec.push_back(Q2LAr);
  genVec.push_back(Q2LArCC0Pi);
  genVec.push_back(Q2LArCC1Pi);
  genVec.push_back(Q2LArCC2Pi);
  genVec.push_back(Q2LArCC3Pi);


  PredictionInterp predLAr(systlist, this_calc, genLAr, loadersLAr);
  PredictionInterp predLArCC0Pi(systlist, this_calc, genLArCC0Pi, loadersLAr);
  PredictionInterp predLArCC1Pi(systlist, this_calc, genLArCC1Pi, loadersLAr);
  PredictionInterp predLArCC2Pi(systlist, this_calc, genLArCC2Pi, loadersLAr);
  PredictionInterp predLArCC3Pi(systlist, this_calc, genLArCC3Pi, loadersLAr);

  PredictionInterp predWLAr(systlist, this_calc, WLAr, loadersLAr);
  PredictionInterp predWLArCC0Pi(systlist, this_calc, WLArCC0Pi, loadersLAr);
  PredictionInterp predWLArCC1Pi(systlist, this_calc, WLArCC1Pi, loadersLAr);
  PredictionInterp predWLArCC2Pi(systlist, this_calc, WLArCC2Pi, loadersLAr);
  PredictionInterp predWLArCC3Pi(systlist, this_calc, WLArCC3Pi, loadersLAr);
  
  loadersLAr.Go();
  loadersGAr.Go();

  std::cout<<systlist.size()<<std::endl;

  TH2 *h2LAr = predLAr.Predict(0).ToTH2(pot_nd);
  TH2 *h2LArNuWro = predLAr.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH2(pot_nd);

  TH1 *hWLArCC0Pi      = predWLArCC0Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC0PiNuWro = predWLArCC0Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC1Pi      = predWLArCC1Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC1PiNuWro = predWLArCC1Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC2Pi      = predWLArCC2Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC2PiNuWro = predWLArCC2Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);
  TH1 *hWLArCC3Pi      = predWLArCC3Pi.Predict(0).ToTH1(pot_nd);
  TH1 *hWLArCC3PiNuWro = predWLArCC3Pi.PredictSyst(0, SystShifts(systlist.at(0), 1)).ToTH1(pot_nd);

  TFile *fout = new TFile("pionPlots.root", "recreate");
  fout->cd();
  h2LAr->Write("h2LAr");
  h2LArNuWro->Write("h2LArNuWro");
  hWLArCC0Pi->Write("hWLArCC0Pi");
  hWLArCC0PiNuWro->Write("hWLArCC0PiNuWro");
  hWLArCC1Pi->Write("hWLArCC1Pi");
  hWLArCC1PiNuWro->Write("hWLArCC1PiNuWro");
  hWLArCC2Pi->Write("hWLArCC2Pi");
  hWLArCC2PiNuWro->Write("hWLArCC2PiNuWro");
  hWLArCC3Pi->Write("hWLArCC3Pi");
  hWLArCC3PiNuWro->Write("hWLArCC3PiNuWro");

  fout->Close();
  delete fout;
} // pionMultiPlots
