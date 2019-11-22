// testNuWroFit.C
// Run an ND+FD and ND only fit to verify it gives the expected behaviour
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
#include "CAFAna/Analysis/common_fit_definitions.h"
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

// ND vars
// Reco X
// const Var kRecoX({"dune.Ev_reco", "dune.Elep_reco", "dune.theta_reco"},
// 		 [](const caf::StandardRecord* sr) {
// 		   double x = 0;
// 		   x = (4 * sr->dune.Ev_reco * sr->dune.Elep_reco * TMath::Sin(sr->dune.theta_reco/2.) * TMath::Sin(sr->dune.theta_reco/2.)) / (2 * 0.939 * (sr->dune.Ev_reco - sr->dune.Elep_reco));
// 		   return x;
// 		 });
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
const Var isCC = SIMPLEVAR(dune.isCC);
const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var KPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
const Var kThetaReco = SIMPLEVAR(dune.theta_reco);

std::vector<double> binEEdges = {0., 0.75, 1.25, 1.5, 1.75, 
				 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75,
				 4., 5., 6., 10.};
std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
const Binning binsEreco  = Binning::Custom(binEEdges);
const Binning simpleBins   = Binning::Simple(20, 0, 5);
const Binning simpleWBins  = Binning::Simple(20, 0, 6);
const Binning simpleQ2Bins = Binning::Simple(20, 0, 6);
const HistAxis axND1d("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND);
const HistAxis axND("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND,
		    "y_{reco}", binYEdges, kRecoYND);
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

void testNuWroFit(const char* outFile, const char* saveDir, bool makeInterps=false) 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, true, true,
						      false, false,
						      false, 20, 
						      true);
  std::vector<const ISyst*> fitsysts = GetListOfSysts(true, true, 
						      true, true, true,
						      false, false,
						      false, 20, 
						      true);
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  systlist.insert(systlist.end(), fakedata.begin(), fakedata.end());
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps of which "<<fitsysts.size()<<" are being fitted"<<std::endl;


  if (makeInterps || TFile(saveDir).IsZombie()) {
    TFile *stateFile = new TFile(saveDir, "recreate");

    TDRLoaders loadersFHC(TDRLoaders::kFHC);
    TDRLoaders loadersRHC(TDRLoaders::kRHC); 
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    PredictionInterp predNumuFHCreco(systlist, this_calc, gennumureco, loadersFHC); 
    PredictionInterp predNumuRHCreco(systlist, this_calc, gennumureco, loadersRHC); 
    PredictionInterp predNueFHCreco(systlist, this_calc, gennuereco, loadersFHC); 
    PredictionInterp predNueRHCreco(systlist, this_calc, gennuereco, loadersRHC);

    loadersFHC.Go();
    loadersRHC.Go();
    std::cout<<"Saving FD numu FHC"<<std::endl;
    predNumuFHCreco.SaveTo(stateFile->mkdir("fd_numu_fhc"));
    std::cout<<"Saving FD nue FHC"<<std::endl;
    predNueFHCreco.SaveTo(stateFile->mkdir("fd_nue_fhc"));
    std::cout<<"Saving FD numu RHC"<<std::endl;
    predNumuRHCreco.SaveTo(stateFile->mkdir("fd_numu_rhc"));
    std::cout<<"Saving FD nue RHC"<<std::endl;
    predNueRHCreco.SaveTo(stateFile->mkdir("fd_nue_rhc"));

    Loaders loadersLArFHC;
    Loaders loadersLArRHC;
    SpectrumLoader loaderLArFHC(Form("%s/ND_FHC_CAF.root", lardir), kBeam);
    SpectrumLoader loaderLArRHC(Form("%s/ND_RHC_CAF.root", lardir), kBeam);
    loadersLArFHC.AddLoader(&loaderLArFHC, caf::kNEARDET, Loaders::kMC);
    loadersLArRHC.AddLoader(&loaderLArRHC, caf::kNEARDET, Loaders::kMC);
    NoOscPredictionGenerator genNDNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueFV);
    PredictionInterp predNDNumuFHC(systlist, this_calc, genNDNumuFHC, loadersLArFHC);
    PredictionInterp predNDNumuRHC(systlist, this_calc, genNDNumuRHC, loadersLArRHC);
    loadersLArFHC.Go();
    loadersLArRHC.Go();
    std::cout<<"Saving ND numu FHC"<<std::endl;
    predNDNumuFHC.SaveTo(stateFile->mkdir("nd_lar_numu_fhc"));
    std::cout<<"Saving ND numu RHC"<<std::endl;
    predNDNumuRHC.SaveTo(stateFile->mkdir("nd_lar_numu_rhc"));

    Loaders loadersGArFHC;
    Loaders loadersGArRHC;
    SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", gardir), kBeam);
    SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", gardir), kBeam);
    loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
    loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);
    // Way too many samples but don't want to have to remake then again
    NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueFV);
    NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueFV);
    PredictionInterp predNDGArNumuFHC(systlist, this_calc, genNDGArNumuFHC, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC(systlist, this_calc, genNDGArNumuRHC, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d(systlist, this_calc, genNDGArNumuFHC1d, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d(systlist, this_calc, genNDGArNumuRHC1d, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2(systlist, this_calc, genNDGArNumuFHCQ2, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2(systlist, this_calc, genNDGArNumuRHCQ2, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW(systlist, this_calc, genNDGArNumuFHCW, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW(systlist, this_calc, genNDGArNumuRHCW, loadersGArRHC);
    // 0pi
    NoOscPredictionGenerator genNDGArNumuFHC_0pi(axND, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuRHC_0pi(axND, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuFHC1d_0pi(axND1d, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuRHC1d_0pi(axND1d, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuFHCQ2_0pi(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuRHCQ2_0pi(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuFHCW_0pi(axRecoW, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==0);
    NoOscPredictionGenerator genNDGArNumuRHCW_0pi(axRecoW, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==0);
    PredictionInterp predNDGArNumuFHC_0pi(systlist, this_calc, genNDGArNumuFHC_0pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC_0pi(systlist, this_calc, genNDGArNumuRHC_0pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d_0pi(systlist, this_calc, genNDGArNumuFHC1d_0pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d_0pi(systlist, this_calc, genNDGArNumuRHC1d_0pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2_0pi(systlist, this_calc, genNDGArNumuFHCQ2_0pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2_0pi(systlist, this_calc, genNDGArNumuRHCQ2_0pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW_0pi(systlist, this_calc, genNDGArNumuFHCW_0pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW_0pi(systlist, this_calc, genNDGArNumuRHCW_0pi, loadersGArRHC);
    // 1pi
    NoOscPredictionGenerator genNDGArNumuFHC_1pi(axND, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuRHC_1pi(axND, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuFHC1d_1pi(axND1d, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuRHC1d_1pi(axND1d, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuFHCQ2_1pi(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuRHCQ2_1pi(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuFHCW_1pi(axRecoW, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==1);
    NoOscPredictionGenerator genNDGArNumuRHCW_1pi(axRecoW, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==1);
    PredictionInterp predNDGArNumuFHC_1pi(systlist, this_calc, genNDGArNumuFHC_1pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC_1pi(systlist, this_calc, genNDGArNumuRHC_1pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d_1pi(systlist, this_calc, genNDGArNumuFHC1d_1pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d_1pi(systlist, this_calc, genNDGArNumuRHC1d_1pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2_1pi(systlist, this_calc, genNDGArNumuFHCQ2_1pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2_1pi(systlist, this_calc, genNDGArNumuRHCQ2_1pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW_1pi(systlist, this_calc, genNDGArNumuFHCW_1pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW_1pi(systlist, this_calc, genNDGArNumuRHCW_1pi, loadersGArRHC);
    // 2pi
    NoOscPredictionGenerator genNDGArNumuFHC_2pi(axND, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuRHC_2pi(axND, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuFHC1d_2pi(axND1d, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuRHC1d_2pi(axND1d, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuFHCQ2_2pi(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuRHCQ2_2pi(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuFHCW_2pi(axRecoW, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==2);
    NoOscPredictionGenerator genNDGArNumuRHCW_2pi(axRecoW, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==2);
    PredictionInterp predNDGArNumuFHC_2pi(systlist, this_calc, genNDGArNumuFHC_2pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC_2pi(systlist, this_calc, genNDGArNumuRHC_2pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d_2pi(systlist, this_calc, genNDGArNumuFHC1d_2pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d_2pi(systlist, this_calc, genNDGArNumuRHC1d_2pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2_2pi(systlist, this_calc, genNDGArNumuFHCQ2_2pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2_2pi(systlist, this_calc, genNDGArNumuRHCQ2_2pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW_2pi(systlist, this_calc, genNDGArNumuFHCW_2pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW_2pi(systlist, this_calc, genNDGArNumuRHCW_2pi, loadersGArRHC);
    // 3pi
    NoOscPredictionGenerator genNDGArNumuFHC_3pi(axND, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuRHC_3pi(axND, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuFHC1d_3pi(axND1d, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuRHC1d_3pi(axND1d, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuFHCQ2_3pi(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuRHCQ2_3pi(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuFHCW_3pi(axRecoW, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi==3);
    NoOscPredictionGenerator genNDGArNumuRHCW_3pi(axRecoW, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi==3);
    PredictionInterp predNDGArNumuFHC_3pi(systlist, this_calc, genNDGArNumuFHC_3pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC_3pi(systlist, this_calc, genNDGArNumuRHC_3pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d_3pi(systlist, this_calc, genNDGArNumuFHC1d_3pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d_3pi(systlist, this_calc, genNDGArNumuRHC1d_3pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2_3pi(systlist, this_calc, genNDGArNumuFHCQ2_3pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2_3pi(systlist, this_calc, genNDGArNumuRHCQ2_3pi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW_3pi(systlist, this_calc, genNDGArNumuFHCW_3pi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW_3pi(systlist, this_calc, genNDGArNumuRHCW_3pi, loadersGArRHC);
    loadersGArFHC.Go();
    loadersGArRHC.Go();

    predNDGArNumuFHC.SaveTo(stateFile->mkdir("nd_gar_numu_fhc"));
    predNDGArNumuRHC.SaveTo(stateFile->mkdir("nd_gar_numu_rhc"));
    predNDGArNumuFHC1d.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1d"));
    predNDGArNumuRHC1d.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1d"));
    predNDGArNumuFHCQ2.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_q2"));
    predNDGArNumuRHCQ2.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_q2"));
    predNDGArNumuFHCW.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_w"));
    predNDGArNumuRHCW.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_w"));

    predNDGArNumuFHC_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_0pi"));
    predNDGArNumuRHC_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_0pi"));
    predNDGArNumuFHC1d_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1d_0pi"));
    predNDGArNumuRHC1d_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1d_0pi"));
    predNDGArNumuFHCQ2_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_q2_0pi"));
    predNDGArNumuRHCQ2_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_q2_0pi"));
    predNDGArNumuFHCW_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_w_0pi"));
    predNDGArNumuRHCW_0pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_w_0pi"));

    predNDGArNumuFHC_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1pi"));
    predNDGArNumuRHC_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1pi"));
    predNDGArNumuFHC1d_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1d_1pi"));
    predNDGArNumuRHC1d_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1d_1pi"));
    predNDGArNumuFHCQ2_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_q2_1pi"));
    predNDGArNumuRHCQ2_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_q2_1pi"));
    predNDGArNumuFHCW_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_w_1pi"));
    predNDGArNumuRHCW_1pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_w_1pi"));

    predNDGArNumuFHC_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_2pi"));
    predNDGArNumuRHC_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_2pi"));
    predNDGArNumuFHC1d_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1d_2pi"));
    predNDGArNumuRHC1d_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1d_2pi"));
    predNDGArNumuFHCQ2_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_q2_2pi"));
    predNDGArNumuRHCQ2_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_q2_2pi"));
    predNDGArNumuFHCW_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_w_2pi"));
    predNDGArNumuRHCW_2pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_w_2pi"));

    predNDGArNumuFHC_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_3pi"));
    predNDGArNumuRHC_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_3pi"));
    predNDGArNumuFHC1d_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_1d_3pi"));
    predNDGArNumuRHC1d_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_1d_3pi"));
    predNDGArNumuFHCQ2_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_q2_3pi"));
    predNDGArNumuRHCQ2_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_q2_3pi"));
    predNDGArNumuFHCW_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_fhc_w_3pi"));
    predNDGArNumuRHCW_3pi.SaveTo(stateFile->mkdir("nd_gar_numu_rhc_w_3pi"));

    stateFile->Close();
    delete stateFile;
  }

  TFile *fout = new TFile(outFile, "recreate");

  fout->Close();
  delete fout;
} // testNuWroFit
