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
#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/TDRLoaders.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/common_fit_definitions.h"
#include "CAFAna/Analysis/Plots.h"
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

std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
const Binning binsEreco  = Binning::Custom(binEEdges);
const Binning binsYreco  = Binning::Custom(binYEdges);
const Binning simpleBins   = Binning::Simple(20, 0, 5);
const Binning simpleWBins  = Binning::Simple(20, 0, 6);
const Binning simpleQ2Bins = Binning::Simple(20, 0, 6);
const HistAxis axND1d("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND);
const HistAxis axND("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND,
		    "y_{reco}", binsYreco, kRecoYND);
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
// Oscillation variables to fit
std::vector<const IFitVar*> fitVars = {&kFitDmSq32NHScaled, &kFitTheta13, &kFitSinSqTheta23, &kFitDeltaInPiUnits, &kFitRho};

void testNuWroFit(const char* outFile, const char* saveDir, 
		  bool makeFDInterps=false, bool makeNDInterps=false,
		  bool makeNDGAr=false,
		  const char *lardir="/pnfs/dune/persistent/users/picker24/CAFv4/", 
		  const char *gardir="/dune/data/users/sbjones/gasTpcCAF/v2/")
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
  // for (int i=0; i<fakedata.size(); i++) {
  //   std::cout<<fakedata[i]->ShortName()<<std::endl;
  // }
  // std::cout<<fakedata.size()<<std::endl;
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps of which "<<fitsysts.size()<<" are being fitted"<<std::endl;


  if (makeFDInterps || TFile(Form("%s/state_FD_FHC.root", saveDir)).IsZombie() ||
      TFile(Form("%s/state_FD_FHC.root", saveDir)).IsZombie()) {
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
    TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", saveDir), "recreate");
    TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", saveDir), "recreate");
    std::cout<<"Saving FD numu FHC"<<std::endl;
    predNumuFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc"));
    std::cout<<"Saving FD nue FHC"<<std::endl;
    predNueFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc"));
    fFDfhc->Close();
    delete fFDfhc;
    std::cout<<"Saving FD numu RHC"<<std::endl;
    predNumuRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc"));
    std::cout<<"Saving FD nue RHC"<<std::endl;
    predNueRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc"));
    fFDrhc->Close();
    delete fFDrhc;
  }
  if (makeNDInterps || TFile(Form("%s/state_ND_FHC.root", saveDir)).IsZombie() ||
      TFile(Form("%s/state_ND_FHC.root", saveDir)).IsZombie()) {
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
    std::cout<<"Saving ND numu FHC to "<<Form("%s/state_ND_FHC.root", saveDir)<<std::endl;
    TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC.root", saveDir), "recreate");
    predNDNumuFHC.SaveTo(fNDfhc->mkdir("nd_interp_numu_fhc"));
    fNDfhc->Close();
    delete fNDfhc;
    std::cout<<"Saving ND numu RHC to "<<Form("%s/state_ND_RHC.root", saveDir)<<std::endl;
    TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC.root", saveDir), "recreate");
    predNDNumuRHC.SaveTo(fNDrhc->mkdir("nd_interp_numu_rhc"));
    fNDrhc->Close();
    delete fNDrhc;
  }

  if (makeNDGAr) {
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
    // >3pi
    NoOscPredictionGenerator genNDGArNumuFHC_hipi(axND, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuRHC_hipi(axND, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuFHC1d_hipi(axND1d, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuRHC1d_hipi(axND1d, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuFHCQ2_hipi(axRecoQ2, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuRHCQ2_hipi(axRecoQ2, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuFHCW_hipi(axRecoW, kPassND_FHC_NUMU && kIsTrueFV && kRecoPi>3);
    NoOscPredictionGenerator genNDGArNumuRHCW_hipi(axRecoW, kPassND_RHC_NUMU && kIsTrueFV && kRecoPi>3);
    PredictionInterp predNDGArNumuFHC_hipi(systlist, this_calc, genNDGArNumuFHC_hipi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC_hipi(systlist, this_calc, genNDGArNumuRHC_hipi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHC1d_hipi(systlist, this_calc, genNDGArNumuFHC1d_hipi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHC1d_hipi(systlist, this_calc, genNDGArNumuRHC1d_hipi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCQ2_hipi(systlist, this_calc, genNDGArNumuFHCQ2_hipi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCQ2_hipi(systlist, this_calc, genNDGArNumuRHCQ2_hipi, loadersGArRHC);
    PredictionInterp predNDGArNumuFHCW_hipi(systlist, this_calc, genNDGArNumuFHCW_hipi, loadersGArFHC);
    PredictionInterp predNDGArNumuRHCW_hipi(systlist, this_calc, genNDGArNumuRHCW_hipi, loadersGArRHC);

    loadersGArFHC.Go();
    loadersGArRHC.Go();

    TFile *fNDGArfhcAll = new TFile(Form("%s/state_ND_GAr_FHC_All.root", saveDir), "recreate");
    predNDGArNumuFHC.SaveTo(fNDGArfhcAll->mkdir("nd_gar_numu_fhc"));
    predNDGArNumuFHC1d.SaveTo(fNDGArfhcAll->mkdir("nd_gar_numu_fhc_1d"));
    predNDGArNumuFHCQ2.SaveTo(fNDGArfhcAll->mkdir("nd_gar_numu_fhc_q2"));
    predNDGArNumuFHCW.SaveTo(fNDGArfhcAll->mkdir("nd_gar_numu_fhc_w"));
    fNDGArfhcAll->Close();
    delete fNDGArfhcAll;
    std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    TFile *fNDGArrhcAll = new TFile(Form("%s/state_ND_GAr_RHC_All.root", saveDir), "recreate");
    predNDGArNumuRHC.SaveTo(fNDGArrhcAll->mkdir("nd_gar_numu_rhc"));    
    predNDGArNumuRHC1d.SaveTo(fNDGArrhcAll->mkdir("nd_gar_numu_rhc_1d"));
    predNDGArNumuRHCQ2.SaveTo(fNDGArrhcAll->mkdir("nd_gar_numu_rhc_q2"));
    predNDGArNumuRHCW.SaveTo(fNDGArrhcAll->mkdir("nd_gar_numu_rhc_w"));
    fNDGArrhcAll->Close();
    delete fNDGArrhcAll;
    std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;

    TFile *fNDGArfhc0pi = new TFile(Form("%s/state_ND_GAr_FHC_0pi.root", saveDir), "recreate");
    predNDGArNumuFHC_0pi.SaveTo(fNDGArfhc0pi->mkdir("nd_gar_numu_fhc_0pi"));
    predNDGArNumuFHC1d_0pi.SaveTo(fNDGArfhc0pi->mkdir("nd_gar_numu_fhc_1d_0pi"));
    predNDGArNumuFHCQ2_0pi.SaveTo(fNDGArfhc0pi->mkdir("nd_gar_numu_fhc_q2_0pi"));
    predNDGArNumuFHCW_0pi.SaveTo(fNDGArfhc0pi->mkdir("nd_gar_numu_fhc_w_0pi"));
    fNDGArfhc0pi->Close();
    delete fNDGArfhc0pi;
    std::cout<<"Saved GAr FHC 0pi samples"<<std::endl;
    TFile *fNDGArrhc0pi = new TFile(Form("%s/state_ND_GAr_RHC_0pi.root", saveDir), "recreate");
    predNDGArNumuRHC_0pi.SaveTo(fNDGArrhc0pi->mkdir("nd_gar_numu_rhc_0pi"));
    predNDGArNumuRHC1d_0pi.SaveTo(fNDGArrhc0pi->mkdir("nd_gar_numu_rhc_1d_0pi"));
    predNDGArNumuRHCQ2_0pi.SaveTo(fNDGArrhc0pi->mkdir("nd_gar_numu_rhc_q2_0pi"));
    predNDGArNumuRHCW_0pi.SaveTo(fNDGArrhc0pi->mkdir("nd_gar_numu_rhc_w_0pi"));
    fNDGArrhc0pi->Close();
    delete fNDGArrhc0pi;
    std::cout<<"Saved GAr RHC 0pi samples"<<std::endl;

    TFile *fNDGArfhc1pi = new TFile(Form("%s/state_ND_GAr_FHC_1pi.root", saveDir), "recreate");
    predNDGArNumuFHC_1pi.SaveTo(fNDGArfhc1pi->mkdir("nd_gar_numu_fhc_1pi"));
    predNDGArNumuFHC1d_1pi.SaveTo(fNDGArfhc1pi->mkdir("nd_gar_numu_fhc_1d_1pi"));
    predNDGArNumuFHCQ2_1pi.SaveTo(fNDGArfhc1pi->mkdir("nd_gar_numu_fhc_q2_1pi"));
    predNDGArNumuFHCW_1pi.SaveTo(fNDGArfhc1pi->mkdir("nd_gar_numu_fhc_w_1pi"));
    fNDGArfhc1pi->Close();
    delete fNDGArfhc1pi;
    std::cout<<"Saved GAr FHC 1pi samples"<<std::endl;
    TFile *fNDGArrhc1pi = new TFile(Form("%s/state_ND_GAr_RHC_1pi.root", saveDir), "recreate");
    predNDGArNumuRHC_1pi.SaveTo(fNDGArrhc1pi->mkdir("nd_gar_numu_rhc_1pi"));
    predNDGArNumuRHC1d_1pi.SaveTo(fNDGArrhc1pi->mkdir("nd_gar_numu_rhc_1d_1pi"));
    predNDGArNumuRHCQ2_1pi.SaveTo(fNDGArrhc1pi->mkdir("nd_gar_numu_rhc_q2_1pi"));
    predNDGArNumuRHCW_1pi.SaveTo(fNDGArrhc1pi->mkdir("nd_gar_numu_rhc_w_1pi"));
    fNDGArrhc1pi->Close();
    delete fNDGArrhc1pi;
    std::cout<<"Saved GAr RHC 1pi samples"<<std::endl;

    TFile *fNDGArfhc2pi = new TFile(Form("%s/state_ND_GAr_FHC_2pi.root", saveDir), "recreate");
    predNDGArNumuFHC_2pi.SaveTo(fNDGArfhc2pi->mkdir("nd_gar_numu_fhc_2pi"));
    predNDGArNumuFHC1d_2pi.SaveTo(fNDGArfhc2pi->mkdir("nd_gar_numu_fhc_1d_2pi"));
    predNDGArNumuFHCQ2_2pi.SaveTo(fNDGArfhc2pi->mkdir("nd_gar_numu_fhc_q2_2pi"));
    predNDGArNumuFHCW_2pi.SaveTo(fNDGArfhc2pi->mkdir("nd_gar_numu_fhc_w_2pi"));
    fNDGArfhc2pi->Close();
    delete fNDGArfhc2pi;
    std::cout<<"Saved GAr FHC 2pi samples"<<std::endl;
    TFile *fNDGArrhc2pi = new TFile(Form("%s/state_ND_GAr_RHC_2pi.root", saveDir), "recreate");
    predNDGArNumuRHC_2pi.SaveTo(fNDGArrhc2pi->mkdir("nd_gar_numu_rhc_2pi"));
    predNDGArNumuRHC1d_2pi.SaveTo(fNDGArrhc2pi->mkdir("nd_gar_numu_rhc_1d_2pi"));
    predNDGArNumuRHCQ2_2pi.SaveTo(fNDGArrhc2pi->mkdir("nd_gar_numu_rhc_q2_2pi"));
    predNDGArNumuRHCW_2pi.SaveTo(fNDGArrhc2pi->mkdir("nd_gar_numu_rhc_w_2pi"));
    fNDGArrhc2pi->Close();
    delete fNDGArrhc2pi;
    std::cout<<"Saved GAr RHC 2pi samples"<<std::endl;

    TFile *fNDGArfhc3pi = new TFile(Form("%s/state_ND_GAr_FHC_3pi.root", saveDir), "recreate");
    predNDGArNumuFHC_3pi.SaveTo(fNDGArfhc3pi->mkdir("nd_gar_numu_fhc_3pi"));
    predNDGArNumuFHC1d_3pi.SaveTo(fNDGArfhc3pi->mkdir("nd_gar_numu_fhc_1d_3pi"));
    predNDGArNumuFHCQ2_3pi.SaveTo(fNDGArfhc3pi->mkdir("nd_gar_numu_fhc_q2_3pi"));
    predNDGArNumuFHCW_3pi.SaveTo(fNDGArfhc3pi->mkdir("nd_gar_numu_fhc_w_3pi"));
    fNDGArfhc3pi->Close();
    delete fNDGArfhc3pi;
    std::cout<<"Saved GAr FHC 3pi samples"<<std::endl;
    TFile *fNDGArrhc3pi = new TFile(Form("%s/state_ND_GAr_RHC_3pi.root", saveDir), "recreate");
    predNDGArNumuRHC_3pi.SaveTo(fNDGArrhc3pi->mkdir("nd_gar_numu_rhc_3pi"));
    predNDGArNumuRHC1d_3pi.SaveTo(fNDGArrhc3pi->mkdir("nd_gar_numu_rhc_1d_3pi"));
    predNDGArNumuRHCQ2_3pi.SaveTo(fNDGArrhc3pi->mkdir("nd_gar_numu_rhc_q2_3pi"));
    predNDGArNumuRHCW_3pi.SaveTo(fNDGArrhc3pi->mkdir("nd_gar_numu_rhc_w_3pi"));
    fNDGArrhc3pi->Close();
    delete fNDGArrhc3pi;
    std::cout<<"Saved GAr RHC 3pi samples"<<std::endl;

    TFile *fNDGArfhchipi = new TFile(Form("%s/state_ND_GAr_FHC_hipi.root", saveDir), "recreate");
    predNDGArNumuFHC_hipi.SaveTo(fNDGArfhchipi->mkdir("nd_gar_numu_fhc_hipi"));
    predNDGArNumuFHC1d_hipi.SaveTo(fNDGArfhchipi->mkdir("nd_gar_numu_fhc_1d_hipi"));
    predNDGArNumuFHCQ2_hipi.SaveTo(fNDGArfhchipi->mkdir("nd_gar_numu_fhc_q2_hipi"));
    predNDGArNumuFHCW_hipi.SaveTo(fNDGArfhchipi->mkdir("nd_gar_numu_fhc_w_hipi"));
    fNDGArfhchipi->Close();
    delete fNDGArfhchipi;
    std::cout<<"Saved GAr FHC hipi samples"<<std::endl;
    TFile *fNDGArrhchipi = new TFile(Form("%s/state_ND_GAr_RHC_hipi.root", saveDir), "recreate");
    predNDGArNumuRHC_hipi.SaveTo(fNDGArrhchipi->mkdir("nd_gar_numu_rhc_hipi"));
    predNDGArNumuRHC1d_hipi.SaveTo(fNDGArrhchipi->mkdir("nd_gar_numu_rhc_1d_hipi"));
    predNDGArNumuRHCQ2_hipi.SaveTo(fNDGArrhchipi->mkdir("nd_gar_numu_rhc_q2_hipi"));
    predNDGArNumuRHCW_hipi.SaveTo(fNDGArrhchipi->mkdir("nd_gar_numu_rhc_w_hipi"));
    fNDGArrhchipi->Close();
    delete fNDGArrhchipi;
    std::cout<<"Saved GAr RHC hipi samples"<<std::endl;
    std::cout<<"All GAr PredictionInterps created"<<std::endl;
  }


  // Retrieve PredictionInterps
  TFile *fNDGArfhcAll = TFile::Open(Form("%s/state_ND_GAr_FHC_All.root", saveDir), "read");
  assert(fNDGArfhcAll);
  PredictionInterp& predNDGArNumuFHC = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_1d_0pi")).release();
  PredictionInterp& predNDGArNumuFHCQ2 = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhcAll->Close();
  std::cout<<"Loaded GAr FHC inclusive samples"<<std::endl;
  TFile *fNDGArrhcAll = TFile::Open(Form("%s/state_ND_GAr_RHC_All.root", saveDir), "read");
  assert(fNDGArrhcAll);
  PredictionInterp& predNDGArNumuRHC = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2 = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhcAll->Close();
  std::cout<<"Loaded GAr RHC inclusive samples"<<std::endl;

  TFile *fNDGArfhc0pi = TFile::Open(Form("%s/state_ND_GAr_FHC_0pi.root", saveDir), "read");
  assert(fNDGArfhc0pi);
  PredictionInterp& predNDGArNumuFHC_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_0pi")).release();
  PredictionInterp& predNDGArNumuFHC1d_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_1d_0pi")).release();
  PredictionInterp& predNDGArNumuFHCQ2_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_q2_0pi")).release();
  PredictionInterp& predNDGArNumuFHCW_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_w_0pi")).release();
  fNDGArfhc0pi->Close();
  std::cout<<"Loaded GAr FHC 0pi samples"<<std::endl;
  TFile *fNDGArrhc0pi = TFile::Open(Form("%s/state_ND_GAr_RHC_0pi.root", saveDir), "read");
  assert(fNDGArrhc0pi);
  PredictionInterp& predNDGArNumuRHC_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_0pi")).release();
  PredictionInterp& predNDGArNumuRHC1d_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_1d_0pi")).release();
  PredictionInterp& predNDGArNumuRHCQ2_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_q2_0pi")).release();
  PredictionInterp& predNDGArNumuRHCW_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_w_0pi")).release();
  fNDGArrhc0pi->Close();
  std::cout<<"Loaded GAr RHC 0pi samples"<<std::endl;

  TFile *fNDGArfhc1pi = TFile::Open(Form("%s/state_ND_GAr_FHC_1pi.root", saveDir), "read");
  assert(fNDGArfhc1pi);
  PredictionInterp& predNDGArNumuFHC_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_1pi")).release();
  PredictionInterp& predNDGArNumuFHC1d_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_1d_1pi")).release();
  PredictionInterp& predNDGArNumuFHCQ2_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_q2_1pi")).release();
  PredictionInterp& predNDGArNumuFHCW_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_w_1pi")).release();
  fNDGArfhc1pi->Close();
  std::cout<<"Loaded GAr FHC 1pi samples"<<std::endl;
  TFile *fNDGArrhc1pi = TFile::Open(Form("%s/state_ND_GAr_RHC_1pi.root", saveDir), "read");
  assert(fNDGArrhc1pi);
  PredictionInterp& predNDGArNumuRHC_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_1pi")).release();
  PredictionInterp& predNDGArNumuRHC1d_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_1d_1pi")).release();
  PredictionInterp& predNDGArNumuRHCQ2_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_q2_1pi")).release();
  PredictionInterp& predNDGArNumuRHCW_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_w_1pi")).release();
  fNDGArrhc1pi->Close();
  std::cout<<"Loaded GAr RHC 1pi samples"<<std::endl;

  TFile *fNDGArfhc2pi = TFile::Open(Form("%s/state_ND_GAr_FHC_2pi.root", saveDir), "read");
  assert(fNDGArfhc2pi);
  PredictionInterp& predNDGArNumuFHC_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_2pi")).release();
  PredictionInterp& predNDGArNumuFHC1d_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_1d_2pi")).release();
  PredictionInterp& predNDGArNumuFHCQ2_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_q2_2pi")).release();
  PredictionInterp& predNDGArNumuFHCW_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_w_2pi")).release();
  fNDGArfhc2pi->Close();
  std::cout<<"Loaded GAr FHC 2pi samples"<<std::endl;
  TFile *fNDGArrhc2pi = TFile::Open(Form("%s/state_ND_GAr_RHC_2pi.root", saveDir), "read");
  assert(fNDGArrhc2pi);
  PredictionInterp& predNDGArNumuRHC_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_2pi")).release();
  PredictionInterp& predNDGArNumuRHC1d_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_1d_2pi")).release();
  PredictionInterp& predNDGArNumuRHCQ2_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_q2_2pi")).release();
  PredictionInterp& predNDGArNumuRHCW_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_w_2pi")).release();
  fNDGArrhc2pi->Close();
  std::cout<<"Loaded GAr RHC 2pi samples"<<std::endl;

  TFile *fNDGArfhc3pi = TFile::Open(Form("%s/state_ND_GAr_FHC_3pi.root", saveDir), "read");
  assert(fNDGArfhc3pi);
  PredictionInterp& predNDGArNumuFHC_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_3pi")).release();
  PredictionInterp& predNDGArNumuFHC1d_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_1d_3pi")).release();
  PredictionInterp& predNDGArNumuFHCQ2_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_q2_3pi")).release();
  PredictionInterp& predNDGArNumuFHCW_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_w_3pi")).release();
  fNDGArfhc3pi->Close();
  std::cout<<"Loaded GAr FHC 3pi samples"<<std::endl;
  TFile *fNDGArrhc3pi = TFile::Open(Form("%s/state_ND_GAr_RHC_3pi.root", saveDir), "read");
  assert(fNDGArrhc3pi);
  PredictionInterp& predNDGArNumuRHC_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_3pi")).release();
  PredictionInterp& predNDGArNumuRHC1d_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_1d_3pi")).release();
  PredictionInterp& predNDGArNumuRHCQ2_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_q2_3pi")).release();
  PredictionInterp& predNDGArNumuRHCW_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_w_3pi")).release();
  fNDGArrhc3pi->Close();
  std::cout<<"Loaded GAr RHC 3pi samples"<<std::endl;

  TFile *fNDGArfhchipi = TFile::Open(Form("%s/state_ND_GAr_FHC_hipi.root", saveDir), "read");
  assert(fNDGArfhchipi);
  PredictionInterp& predNDGArNumuFHC_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_hipi")).release();
  PredictionInterp& predNDGArNumuFHC1d_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_1d_hipi")).release();
  PredictionInterp& predNDGArNumuFHCQ2_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_q2_hipi")).release();
  PredictionInterp& predNDGArNumuFHCW_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_w_hipi")).release();
  fNDGArfhchipi->Close();
  std::cout<<"Loaded GAr FHC hipi samples"<<std::endl;
  TFile *fNDGArrhchipi = TFile::Open(Form("%s/state_ND_GAr_RHC_hipi.root", saveDir), "read");
  assert(fNDGArrhchipi);
  PredictionInterp& predNDGArNumuRHC_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_hipi")).release();
  PredictionInterp& predNDGArNumuRHC1d_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_1d_hipi")).release();
  PredictionInterp& predNDGArNumuRHCQ2_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_q2_hipi")).release();
  PredictionInterp& predNDGArNumuRHCW_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_w_hipi")).release();
  fNDGArrhchipi->Close();
  std::cout<<"Loaded GAr RHC hipi samples"<<std::endl;

  // Get the readymade PredictionInterps for the LAr samples
  // Put the LAr samples in a vector 
  // Order is 0=nd numu fhc, 1=fd numu fhc, 2=fd nue fhc
  // 3=nd numu rhc, 4=fd numu rhc, 5=fd nue rhc
  TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC.root", saveDir), "read");
  TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC.root", saveDir), "read");
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", saveDir), "read");
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", saveDir), "read");

  std::vector<std::unique_ptr<ana::PredictionInterp>> predLAr;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fNDfhc->GetDirectory("nd_interp_numu_fhc")));
  std::cout<<"Loaded ND numu FHC"<<std::endl;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_numu_fhc")));
  std::cout<<"Loaded FD numu FHC"<<std::endl;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_nue_fhc")));
  std::cout<<"Loaded FD nue FHC"<<std::endl;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fNDrhc->GetDirectory("nd_interp_numu_rhc")));
  std::cout<<"Loaded ND numu RHC"<<std::endl;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_numu_rhc")));
  std::cout<<"Loaded FD numu RHC"<<std::endl;
  predLAr.emplace_back(LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_nue_rhc")));
  std::cout<<"Loaded FD nue RHC"<<std::endl;

  std::cout<<"Loaded all the LAr samples"<<std::endl;

  PredictionInterp &predNDLArNumuFHC = *predLAr[0].release();
  PredictionInterp &predFDNumuFHC    = *predLAr[1].release();
  PredictionInterp &predFDNueFHC     = *predLAr[2].release();
  PredictionInterp &predNDLArNumuRHC = *predLAr[3].release();
  PredictionInterp &predFDNumuRHC    = *predLAr[4].release();
  PredictionInterp &predFDNueRHC     = *predLAr[5].release();

  fNDfhc->Close();
  fNDrhc->Close();
  fFDfhc->Close();
  fFDrhc->Close();

  TFile *fout = new TFile(outFile, "recreate");
  fout->cd();

  cout<<this_calc->GetTh13()<<", "<<this_calc->GetTh23()<<endl;
  // Build experiments
  std::cout<<"Building SingleSampleExperiments"<<std::endl;
  SingleSampleExperiment nd_lar_fhc(&predNDLArNumuFHC, predNDLArNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));
  SingleSampleExperiment nd_lar_rhc(&predNDLArNumuRHC, predNDLArNumuRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));
  SingleSampleExperiment fd_numu_fhc(&predFDNumuFHC, predFDNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd));
  SingleSampleExperiment fd_numu_rhc(&predFDNumuRHC, predFDNumuRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd));
  SingleSampleExperiment fd_nue_fhc(&predFDNueFHC, predFDNueFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd));
  SingleSampleExperiment fd_nue_rhc(&predFDNueRHC, predFDNueRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd));
  std::cout<<"Built SingleSampleExperiments"<<std::endl;

  SystShifts junk = kNoShift;

  // Output data spectra
  TH1 *nd_lar_fhc_data = predNDLArNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd).ToTH1(pot_nd);
  nd_lar_fhc_data->SetName("nd_lar_fhc_data");
  nd_lar_fhc_data->Write();
  TH1 *pre_nd_lar_fhc = GetMCSystTotal(&predNDLArNumuFHC, this_calc, junk,                   
				       "prefit_nd_lar_fhc", pot_fd);
  pre_nd_lar_fhc->SetTitle(std::to_string(nd_lar_fhc.ChiSq(this_calc, junk)).c_str());
  pre_nd_lar_fhc->Write();
  TH1 *nd_lar_rhc_data = predNDLArNumuRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd).ToTH1(pot_nd);
  nd_lar_rhc_data->SetName("nd_lar_rhc_data");
  nd_lar_rhc_data->Write();
  TH1 *pre_nd_lar_rhc = GetMCSystTotal(&predNDLArNumuRHC, this_calc, junk,                   
				       "prefit_nd_lar_rhc", pot_fd);
  pre_nd_lar_rhc->SetTitle(std::to_string(nd_lar_rhc.ChiSq(this_calc, junk)).c_str());
  pre_nd_lar_rhc->Write();
  TH1 *fd_numu_fhc_data = predFDNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  fd_numu_fhc_data->SetName("fd_numu_fhc_data");
  fd_numu_fhc_data->Write();
  TH1 *pre_fd_numu_fhc = GetMCSystTotal(&predFDNumuFHC, this_calc, junk,               
					"prefit_fd_nue_fhc", pot_fd);
  pre_fd_numu_fhc->SetTitle(std::to_string(fd_numu_fhc.ChiSq(this_calc, junk)).c_str());
  pre_fd_numu_fhc->Write();
  TH1 *fd_numu_rhc_data = predFDNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  fd_numu_rhc_data->SetName("fd_numu_rhc_data");
  fd_numu_rhc_data->Write();
  TH1 *pre_fd_numu_rhc = GetMCSystTotal(&predFDNumuRHC, this_calc, junk,                 
					"prefit_fd_numu_rhc", pot_fd);
  pre_fd_numu_rhc->SetTitle(std::to_string(fd_numu_rhc.ChiSq(this_calc, junk)).c_str());
  pre_fd_numu_rhc->Write();
  TH1 *fd_nue_fhc_data = predFDNueFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  fd_nue_fhc_data->SetName("fd_nue_fhc_data");
  fd_nue_fhc_data->Write();
  TH1 *pre_fd_nue_fhc = GetMCSystTotal(&predFDNueFHC, this_calc, junk,                   
				       "prefit_fd_nue_fhc", pot_fd);
  pre_fd_nue_fhc->SetTitle(std::to_string(fd_nue_fhc.ChiSq(this_calc, junk)).c_str());
  pre_fd_nue_fhc->Write();
  TH1 *fd_nue_rhc_data = predFDNueFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  fd_nue_rhc_data->SetName("fd_nue_rhc_data");
  fd_nue_rhc_data->Write();
  TH1 *pre_fd_nue_rhc = GetMCSystTotal(&predFDNueRHC, this_calc, junk,                  
				       "prefit_fd_nue_rhc", pot_fd);
  pre_fd_nue_rhc->SetTitle(std::to_string(fd_nue_rhc.ChiSq(this_calc, junk)).c_str());
  pre_fd_nue_rhc->Write();

  MultiExperiment fd_only_expt, nd_fd_expt;
  fd_only_expt.Add(&fd_numu_fhc);
  fd_only_expt.Add(&fd_numu_rhc);
  fd_only_expt.Add(&fd_nue_fhc);
  fd_only_expt.Add(&fd_nue_rhc);

  nd_fd_expt.Add(&nd_lar_fhc);
  nd_fd_expt.Add(&nd_lar_rhc);
  nd_fd_expt.Add(&fd_numu_fhc);
  nd_fd_expt.Add(&fd_numu_rhc);
  nd_fd_expt.Add(&fd_nue_fhc);
  nd_fd_expt.Add(&fd_nue_rhc);

  cout<<this_calc->GetTh13()<<", "<<this_calc->GetTh23()<<endl;

  std::string covFileName = FindCAFAnaDir() + "/Systs/det_sys_cov.root";
  nd_fd_expt.AddCovarianceMatrix(covFileName, "nd_all_frac_cov", true, {0, 1});

  SystShifts fitsyst = kNoShift;

  auto start_fit = std::chrono::system_clock::now();
  // Now set up the fit itself                                  
  std::map<const IFitVar*, std::vector<double>> oscSeeds = {};                                        
  osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(1);
  std::cerr << "[INFO]: Beginning fit. " << BuildLogInfoString();
  Fitter fd_only_fit(&fd_only_expt, fitVars, fitsysts);
  double fd_only_chisq = fd_only_fit.Fit(testOsc, junk, oscSeeds, {}, Fitter::kVerbose);
  auto end_fit = std::chrono::system_clock::now();
  std::time_t end_fit_time = std::chrono::system_clock::to_time_t(end_fit);
  std::cerr << "[FIT]: Finished fit in "
            << std::chrono::duration_cast<std::chrono::seconds>(end_fit -
                                                                start_fit).count()
            << " s after " << fd_only_fit.GetNFCN() << " iterations "
            << BuildLogInfoString();

  TTree *fd_only_tree = new TTree("fd_only", "fd_only");
  std::vector<std::string> fParamNames;
  std::vector<double> fPreFitValues;
  std::vector<double> fPreFitErrors;
  std::vector<double> fPostFitValues;
  std::vector<double> fPostFitErrors;
  std::vector<double> fCentralValues;
  std::vector<double> fTrueValues;
  double fNFCN, fEDM, fChiSq;
  bool fIsValid;
  fd_only_tree->Branch("ParamNames", &fParamNames);
  fd_only_tree->Branch("PreFitValues", &fPreFitValues);
  fd_only_tree->Branch("PreFitErrors", &fPreFitErrors);
  fd_only_tree->Branch("PostFitValues", &fPostFitValues);
  fd_only_tree->Branch("PostFitErrors", &fPostFitErrors);
  fd_only_tree->Branch("CentralValues", &fCentralValues);
  fd_only_tree->Branch("TrueValues", &fTrueValues);
  fd_only_tree->Branch("NFCN", &fNFCN);
  fd_only_tree->Branch("EDM", &fEDM);
  fd_only_tree->Branch("ChiSq", &fChiSq);
  fd_only_tree->Branch("fIsValid", &fIsValid);

  fParamNames = fd_only_fit.GetParamNames();
  fPreFitValues = fd_only_fit.GetPreFitValues();
  fPreFitErrors = fd_only_fit.GetPreFitErrors();
  fPostFitValues = fd_only_fit.GetPostFitValues();
  fPostFitErrors = fd_only_fit.GetPostFitErrors();
  fCentralValues = fd_only_fit.GetCentralValues();
  fNFCN = fd_only_fit.GetNFCN();
  fEDM = fd_only_fit.GetEDM();
  fIsValid = fd_only_fit.GetIsValid();
  fChiSq = fd_only_chisq;

  fTrueValues.resize(fParamNames.size(), 0.);
  fTrueValues.at(0) = this_calc->GetDmsq32()*1000.;
  fTrueValues.at(1) = this_calc->GetTh13();
  fTrueValues.at(2) = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
  fTrueValues.at(3) = this_calc->GetdCP();
  fTrueValues.at(4) = this_calc->GetRho();

  fd_only_tree->Fill();
  fd_only_tree->Write();

  cout<<this_calc->GetTh13()<<", "<<this_calc->GetTh23()<<endl;

  TMatrixDSym *covar = (TMatrixDSym *)fd_only_fit.GetCovariance();
  TH2D hist_covar = TH2D(*covar);
  hist_covar.SetName("covar");                                                                         
  TH2D hist_corr = *make_corr_from_covar(&hist_covar);
  hist_covar.Write();
  hist_corr.Write();
  covar->Write("covar_mat");

  // Now do the joint fit
  Fitter nd_fd_fit(&nd_fd_expt, fitVars, fitsysts);
  double nd_fd_chisq = nd_fd_fit.Fit(testOsc, junk, oscSeeds, {}, Fitter::kVerbose);

  TTree *nd_fd_tree = new TTree("nd_fd", "nd_fd");
  nd_fd_tree->Branch("ParamNames", &fParamNames);
  nd_fd_tree->Branch("PreFitValues", &fPreFitValues);
  nd_fd_tree->Branch("PreFitErrors", &fPreFitErrors);
  nd_fd_tree->Branch("PostFitValues", &fPostFitValues);
  nd_fd_tree->Branch("PostFitErrors", &fPostFitErrors);
  nd_fd_tree->Branch("CentralValues", &fCentralValues);
  nd_fd_tree->Branch("TrueValues", &fTrueValues);
  nd_fd_tree->Branch("NFCN", &fNFCN);
  nd_fd_tree->Branch("EDM", &fEDM);
  nd_fd_tree->Branch("ChiSq", &fChiSq);
  nd_fd_tree->Branch("fIsValid", &fIsValid);

  fParamNames = nd_fd_fit.GetParamNames();
  fPreFitValues = nd_fd_fit.GetPreFitValues();
  fPreFitErrors = nd_fd_fit.GetPreFitErrors();
  fPostFitValues = nd_fd_fit.GetPostFitValues();
  fPostFitErrors = nd_fd_fit.GetPostFitErrors();
  fCentralValues = nd_fd_fit.GetCentralValues();
  fNFCN = nd_fd_fit.GetNFCN();
  fEDM = nd_fd_fit.GetEDM();
  fIsValid = nd_fd_fit.GetIsValid();
  fChiSq = nd_fd_chisq;
  nd_fd_tree->Fill();
  nd_fd_tree->Write();

  cout<<this_calc->GetTh13()<<", "<<this_calc->GetTh23()<<endl;

  fout->Close();
  delete fout;
} // testNuWroFit
