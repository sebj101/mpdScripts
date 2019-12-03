// makePredInterpsGas.C
// Makes the PredictionInterps for the necessary gas TPC samples
// Should hopefully make the times tolerable
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
// // ROOT includes
#include "TCanvas.h"
#include "TFile.h"

#include "Utilities/rootlogon.C"

using namespace ana;

// Define Vars for analysis
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
const Var kPiplmult  = SIMPLEVAR(dune.gastpc_pi_pl_mult);
const Var KPiminmult = SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoChargedPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult);
const Var kRecoPi = SIMPLEVAR(dune.gastpc_pi_pl_mult) + SIMPLEVAR(dune.gastpc_pi_min_mult) + SIMPLEVAR(dune.gastpc_pi_0_mult);
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
// Axes
const HistAxis axND1d("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND);
const HistAxis axND("E_{#nu, reco} (GeV)", binsEreco, kRecoEnergyND,
		    "y_{reco}", binsYreco, kRecoYND);
const HistAxis axTrueRecoND("E_{#nu, true} (GeV)", simpleBins, kTrueEnergy,
			    "E_{#nu, reco} (GeV)", simpleBins, kRecoEnergyND);
// const HistAxis axQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kQ2);
// const HistAxis axW("W (GeV)", simpleWBins, kW);
const HistAxis axRecoQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kRecoQ2);
const HistAxis axRecoW("W (GeV)", simpleWBins, kRecoW);

void makePredInterpsGas(const bool isFHC, const int multi, 
			const char* outDir="/pnfs/dune/persistent/users/sbjones/stateFiles/",
			const char *garDir="/dune/data/users/sbjones/gasTpcCAF/v2/") 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();
  // Check inputs
  if (isFHC) {
    std::cout<<"Forward horn current selected"<<std::endl;
  }
  else {
    std::cout<<"Reverse horn current selected"<<std::endl;
  }

  assert(multi >= -1 && multi <= 4);

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, true, true,
						      false, false,
						      false, 20, 
						      true);
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps"<<std::endl;

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", garDir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", garDir), kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  if (multi == -1) { // Inclusive sample
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_All.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_All.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
  else if (multi==0) {
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_0pi.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==0);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_0pi.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
  else if (multi==1) {
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_1pi.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==1);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_1pi.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
  else if (multi==2) {
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_2pi.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==2);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_2pi.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
  else if (multi==3) {
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_3pi.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi==3);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_3pi.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
  else if (multi==4) {
    if (isFHC) {
      NoOscPredictionGenerator genNDGArNumuFHC(axND, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuFHC1d(axND1d, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuFHCQ2(axRecoQ2, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuFHCW(axRecoW, kPassND_FHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      PredictionInterp predNDGArNumuFHC(systlist, 0, genNDGArNumuFHC, loadersGArFHC);
      PredictionInterp predNDGArNumuFHC1d(systlist, 0, genNDGArNumuFHC1d, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCQ2(systlist, 0, genNDGArNumuFHCQ2, loadersGArFHC);
      PredictionInterp predNDGArNumuFHCW(systlist, 0, genNDGArNumuFHCW, loadersGArFHC);
      loadersGArFHC.Go();
      TFile *fNDGArfhc = new TFile(Form("%s/state_ND_GAr_FHC_hipi.root", outDir), "recreate");
      predNDGArNumuFHC.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc"));
      predNDGArNumuFHC1d.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_1d"));
      predNDGArNumuFHCQ2.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_q2"));
      predNDGArNumuFHCW.SaveTo(fNDGArfhc->mkdir("nd_gar_numu_fhc_w"));
      fNDGArfhc->Close();
      delete fNDGArfhc;
      std::cout<<"Saved GAr FHC inclusive samples"<<std::endl;
    }
    else {
      NoOscPredictionGenerator genNDGArNumuRHC(axND, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuRHC1d(axND1d, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuRHCQ2(axRecoQ2, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      NoOscPredictionGenerator genNDGArNumuRHCW(axRecoW, kPassND_RHC_NUMU && kIsTrueGasFV && kRecoPi>3);
      PredictionInterp predNDGArNumuRHC(systlist, 0, genNDGArNumuRHC, loadersGArRHC);
      PredictionInterp predNDGArNumuRHC1d(systlist, 0, genNDGArNumuRHC1d, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCQ2(systlist, 0, genNDGArNumuRHCQ2, loadersGArRHC);
      PredictionInterp predNDGArNumuRHCW(systlist, 0, genNDGArNumuRHCW, loadersGArRHC);
      loadersGArRHC.Go();
      TFile *fNDGArrhc = new TFile(Form("%s/state_ND_GAr_RHC_hipi.root", outDir), "recreate");
      predNDGArNumuRHC.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc"));    
      predNDGArNumuRHC1d.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_1d"));
      predNDGArNumuRHCQ2.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_q2"));
      predNDGArNumuRHCW.SaveTo(fNDGArrhc->mkdir("nd_gar_numu_rhc_w"));
      fNDGArrhc->Close();
      delete fNDGArrhc;
      std::cout<<"Saved GAr RHC inclusive samples"<<std::endl;
    }
  }
} // makePredInterpsGas
