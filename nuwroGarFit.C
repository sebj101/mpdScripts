// nuwroGarFit.C
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

void nuwroGarFit(const char* outFile, 
		 const char* stateFileDir="/pnfs/dune/persistent/users/sbjones/stateFiles/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  this_calc->SetdCP(3*TMath::Pi()/2.);

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

  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps of which "<<fitsysts.size()<<" are being fitted"<<std::endl;

  // Retrieve PredictionInterps
  // GAr first
  TFile *fNDGArfhcAll = new TFile(Form("%s/state_ND_GAr_FHC_All.root", stateFileDir), "read");
  assert(fNDGArfhcAll);
  PredictionInterp& predNDGArNumuFHC = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2 = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW = *ana::LoadFrom<PredictionInterp>(fNDGArfhcAll->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhcAll->Close();
  std::cout<<"Loaded GAr FHC inclusive samples"<<std::endl;
  TFile *fNDGArrhcAll = new TFile(Form("%s/state_ND_GAr_RHC_All.root", stateFileDir), "read");
  assert(fNDGArrhcAll);
  PredictionInterp& predNDGArNumuRHC = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2 = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW = *ana::LoadFrom<PredictionInterp>(fNDGArrhcAll->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhcAll->Close();
  std::cout<<"Loaded GAr RHC inclusive samples"<<std::endl;

  TFile *fNDGArfhc0pi = new TFile(Form("%s/state_ND_GAr_FHC_0pi.root", stateFileDir), "read");
  assert(fNDGArfhc0pi);
  PredictionInterp& predNDGArNumuFHC_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc0pi->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhc0pi->Close();
  std::cout<<"Loaded GAr FHC 0pi samples"<<std::endl;
  TFile *fNDGArrhc0pi = new TFile(Form("%s/state_ND_GAr_RHC_0pi.root", stateFileDir), "read");
  assert(fNDGArrhc0pi);
  PredictionInterp& predNDGArNumuRHC_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW_0pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc0pi->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhc0pi->Close();
  std::cout<<"Loaded GAr RHC 0pi samples"<<std::endl;

  TFile *fNDGArfhc1pi = new TFile(Form("%s/state_ND_GAr_FHC_1pi.root", stateFileDir), "read");
  assert(fNDGArfhc1pi);
  PredictionInterp& predNDGArNumuFHC_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc1pi->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhc1pi->Close();
  std::cout<<"Loaded GAr FHC 1pi samples"<<std::endl;
  TFile *fNDGArrhc1pi = new TFile(Form("%s/state_ND_GAr_RHC_1pi.root", stateFileDir), "read");
  assert(fNDGArrhc1pi);
  PredictionInterp& predNDGArNumuRHC_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW_1pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc1pi->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhc1pi->Close();
  std::cout<<"Loaded GAr RHC 1pi samples"<<std::endl;

  TFile *fNDGArfhc2pi = new TFile(Form("%s/state_ND_GAr_FHC_2pi.root", stateFileDir), "read");
  assert(fNDGArfhc2pi);
  PredictionInterp& predNDGArNumuFHC_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc2pi->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhc2pi->Close();
  std::cout<<"Loaded GAr FHC 2pi samples"<<std::endl;
  TFile *fNDGArrhc2pi = new TFile(Form("%s/state_ND_GAr_RHC_2pi.root", stateFileDir), "read");
  assert(fNDGArrhc2pi);
  PredictionInterp& predNDGArNumuRHC_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW_2pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc2pi->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhc2pi->Close();
  std::cout<<"Loaded GAr RHC 2pi samples"<<std::endl;

  TFile *fNDGArfhc3pi = new TFile(Form("%s/state_ND_GAr_FHC_3pi.root", stateFileDir), "read");
  assert(fNDGArfhc3pi);
  PredictionInterp& predNDGArNumuFHC_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArfhc3pi->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhc3pi->Close();
  std::cout<<"Loaded GAr FHC 3pi samples"<<std::endl;
  TFile *fNDGArrhc3pi = new TFile(Form("%s/state_ND_GAr_RHC_3pi.root", stateFileDir), "read");
  assert(fNDGArrhc3pi);
  PredictionInterp& predNDGArNumuRHC_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW_3pi = *ana::LoadFrom<PredictionInterp>(fNDGArrhc3pi->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhc3pi->Close();
  std::cout<<"Loaded GAr RHC 3pi samples"<<std::endl;

  TFile *fNDGArfhchipi = new TFile(Form("%s/state_ND_GAr_FHC_hipi.root", stateFileDir), "read");
  assert(fNDGArfhchipi);
  PredictionInterp& predNDGArNumuFHC_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc")).release();
  PredictionInterp& predNDGArNumuFHC1d_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_1d")).release();
  PredictionInterp& predNDGArNumuFHCQ2_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_q2")).release();
  PredictionInterp& predNDGArNumuFHCW_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArfhchipi->GetDirectory("nd_gar_numu_fhc_w")).release();
  fNDGArfhchipi->Close();
  std::cout<<"Loaded GAr FHC hipi samples"<<std::endl;
  TFile *fNDGArrhchipi = new TFile(Form("%s/state_ND_GAr_RHC_hipi.root", stateFileDir), "read");
  assert(fNDGArrhchipi);
  PredictionInterp& predNDGArNumuRHC_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc")).release();
  PredictionInterp& predNDGArNumuRHC1d_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_1d")).release();
  PredictionInterp& predNDGArNumuRHCQ2_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_q2")).release();
  PredictionInterp& predNDGArNumuRHCW_hipi = *ana::LoadFrom<PredictionInterp>(fNDGArrhchipi->GetDirectory("nd_gar_numu_rhc_w")).release();
  fNDGArrhchipi->Close();
  std::cout<<"Loaded GAr RHC hipi samples"<<std::endl;

  // Get the readymade PredictionInterps for the LAr samples
  // Put the LAr samples in a vector 
  // Order is 0=nd numu fhc, 1=fd numu fhc, 2=fd nue fhc
  // 3=nd numu rhc, 4=fd numu rhc, 5=fd nue rhc
  TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC.root", stateFileDir), "read");
  TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC.root", stateFileDir), "read");
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "read");
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "read");

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

  // Build experiments
  std::cout<<"Building SingleSampleExperiments"<<std::endl;
  // LAr
  SingleSampleExperiment nd_lar_fhc(&predNDLArNumuFHC, predNDLArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  SingleSampleExperiment nd_lar_rhc(&predNDLArNumuRHC, predNDLArNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_lar_fhc_data    = predNDLArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_lar_fhc_noshift = predNDLArNumuFHC.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  TH2 *nd_lar_fhc_data_2d = predNDLArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *nd_lar_fhc_noshift_2d = predNDLArNumuFHC.Predict(this_calc).FakeData(pot_nd).ToTH2(pot_nd);
  nd_lar_fhc_data->SetName("nd_lar_fhc_data");
  nd_lar_fhc_noshift->SetName("nd_lar_fhc_noshift");
  nd_lar_fhc_data_2d->SetName("nd_lar_fhc_data_2d");
  nd_lar_fhc_noshift_2d->SetName("nd_lar_fhc_noshift_2d");
  nd_lar_fhc_data->Write();
  nd_lar_fhc_noshift->Write();
  nd_lar_fhc_data_2d->Write();
  nd_lar_fhc_noshift_2d->Write();
  TH1 *nd_lar_rhc_data    = predNDLArNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_lar_rhc_noshift = predNDLArNumuRHC.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  TH2 *nd_lar_rhc_data_2d = predNDLArNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *nd_lar_rhc_noshift_2d = predNDLArNumuRHC.Predict(this_calc).FakeData(pot_nd).ToTH2(pot_nd);
  nd_lar_rhc_data->SetName("nd_lar_rhc_data");
  nd_lar_rhc_noshift->SetName("nd_lar_rhc_noshift");
  nd_lar_rhc_data_2d->SetName("nd_lar_rhc_data_2d");
  nd_lar_rhc_noshift_2d->SetName("nd_lar_rhc_noshift_2d");
  nd_lar_rhc_data->Write();
  nd_lar_rhc_noshift->Write();
  nd_lar_rhc_data_2d->Write();
  nd_lar_rhc_noshift_2d->Write();
  SingleSampleExperiment fd_numu_fhc(&predFDNumuFHC, predFDNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  TH1 *fd_numu_fhc_data = predFDNumuFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *fd_numu_fhc_noshift = predFDNumuFHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  fd_numu_fhc_data->SetName("fd_numu_fhc_data");
  fd_numu_fhc_noshift->SetName("fd_numu_fhc_noshift");
  fd_numu_fhc_data->Write();
  fd_numu_fhc_noshift->Write();
  SingleSampleExperiment fd_numu_rhc(&predFDNumuRHC, predFDNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  TH1 *fd_numu_rhc_data = predFDNumuRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *fd_numu_rhc_noshift = predFDNumuRHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  fd_numu_rhc_data->SetName("fd_numu_rhc_data");
  fd_numu_rhc_noshift->SetName("fd_numu_rhc_noshift");
  fd_numu_rhc_data->Write();
  fd_numu_rhc_noshift->Write();
  SingleSampleExperiment fd_nue_fhc(&predFDNueFHC, predFDNueFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  TH1 *fd_nue_fhc_data = predFDNueFHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *fd_nue_fhc_noshift = predFDNueFHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  fd_nue_fhc_data->SetName("fd_nue_fhc_data");
  fd_nue_fhc_noshift->SetName("fd_nue_fhc_noshift");
  fd_nue_fhc_data->Write();
  fd_nue_fhc_noshift->Write();
  SingleSampleExperiment fd_nue_rhc(&predFDNueRHC, predFDNueRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  TH1 *fd_nue_rhc_data = predFDNueRHC.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *fd_nue_rhc_noshift = predFDNueRHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  fd_nue_rhc_data->SetName("fd_nue_rhc_data");
  fd_nue_rhc_noshift->SetName("fd_nue_rhc_noshift");
  fd_nue_rhc_data->Write();
  fd_nue_rhc_noshift->Write();
  // GAr
  // 2D Ereco vs Yrec
  SingleSampleExperiment nd_gar_fhc_all(&predNDGArNumuFHC, predNDGArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_all_data = predNDGArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_gar_fhc_all_noshift = predNDGArNumuFHC.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_all_data->SetName("nd_gar_fhc_all_data");
  nd_gar_fhc_all_noshift->SetName("nd_gar_fhc_all_noshift");
  nd_gar_fhc_all_data->Write();
  nd_gar_fhc_all_noshift->Write();
  SingleSampleExperiment nd_gar_rhc_all(&predNDGArNumuFHC, predNDGArNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_all_data = predNDGArNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_gar_rhc_all_noshift = predNDGArNumuRHC.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_all_data->SetName("nd_gar_rhc_all_data");
  nd_gar_rhc_all_noshift->SetName("nd_gar_rhc_all_noshift");
  nd_gar_rhc_all_data->Write();
  nd_gar_rhc_all_noshift->Write();

  SingleSampleExperiment nd_gar_fhc_0pi(&predNDGArNumuFHC_0pi, predNDGArNumuFHC_0pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_0pi_data = predNDGArNumuFHC_0pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_0pi_data->SetName("nd_gar_fhc_0pi_data");
  nd_gar_fhc_0pi_data->Write();
  SingleSampleExperiment nd_gar_rhc_0pi(&predNDGArNumuFHC_0pi, predNDGArNumuFHC_0pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_0pi_data = predNDGArNumuRHC_0pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_0pi_data->SetName("nd_gar_rhc_0pi_data");
  nd_gar_rhc_0pi_data->Write();

  SingleSampleExperiment nd_gar_fhc_1pi(&predNDGArNumuFHC_1pi, predNDGArNumuFHC_1pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_1pi_data = predNDGArNumuFHC_1pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_1pi_data->SetName("nd_gar_fhc_1pi_data");
  nd_gar_fhc_1pi_data->Write();
  SingleSampleExperiment nd_gar_rhc_1pi(&predNDGArNumuFHC_1pi, predNDGArNumuFHC_1pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_1pi_data = predNDGArNumuRHC_1pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_1pi_data->SetName("nd_gar_rhc_1pi_data");
  nd_gar_rhc_1pi_data->Write();

  SingleSampleExperiment nd_gar_fhc_2pi(&predNDGArNumuFHC_2pi, predNDGArNumuFHC_2pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_2pi_data = predNDGArNumuFHC_2pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_2pi_data->SetName("nd_gar_fhc_2pi_data");
  nd_gar_fhc_2pi_data->Write();
  SingleSampleExperiment nd_gar_rhc_2pi(&predNDGArNumuFHC_2pi, predNDGArNumuFHC_2pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_2pi_data = predNDGArNumuRHC_2pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_2pi_data->SetName("nd_gar_rhc_2pi_data");
  nd_gar_rhc_2pi_data->Write();

  SingleSampleExperiment nd_gar_fhc_3pi(&predNDGArNumuFHC_3pi, predNDGArNumuFHC_3pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_3pi_data = predNDGArNumuFHC_3pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_3pi_data->SetName("nd_gar_fhc_3pi_data");
  nd_gar_fhc_3pi_data->Write();
  SingleSampleExperiment nd_gar_rhc_3pi(&predNDGArNumuFHC_3pi, predNDGArNumuFHC_3pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_3pi_data = predNDGArNumuRHC_3pi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_3pi_data->SetName("nd_gar_rhc_3pi_data");
  nd_gar_rhc_3pi_data->Write();

  SingleSampleExperiment nd_gar_fhc_hipi(&predNDGArNumuFHC_hipi, predNDGArNumuFHC_hipi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_hipi_data = predNDGArNumuFHC_hipi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_hipi_data->SetName("nd_gar_fhc_hipi_data");
  nd_gar_fhc_hipi_data->Write();
  SingleSampleExperiment nd_gar_rhc_hipi(&predNDGArNumuFHC_hipi, predNDGArNumuFHC_hipi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_hipi_data = predNDGArNumuRHC_hipi.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_hipi_data->SetName("nd_gar_rhc_hipi_data");
  nd_gar_rhc_hipi_data->Write();
  std::cout<<"Built SingleSampleExperiments"<<std::endl;

  // 1D Erec
  SingleSampleExperiment nd_gar_fhc_1d_all(&predNDGArNumuFHC1d, predNDGArNumuFHC1d.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_fhc_1d_all_data = predNDGArNumuFHC1d.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_gar_fhc_1d_all_noshift = predNDGArNumuFHC1d.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_fhc_1d_all_data->SetName("nd_gar_fhc_1d_all_data");
  nd_gar_fhc_1d_all_noshift->SetName("nd_gar_fhc_1d_all_noshift");
  nd_gar_fhc_1d_all_data->Write();
  nd_gar_fhc_1d_all_noshift->Write();
  SingleSampleExperiment nd_gar_rhc_1d_all(&predNDGArNumuFHC1d, predNDGArNumuFHC1d.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd));
  TH1 *nd_gar_rhc_1d_all_data = predNDGArNumuRHC1d.PredictSyst(this_calc, fakedatashift).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *nd_gar_rhc_1d_all_noshift = predNDGArNumuRHC1d.Predict(this_calc).FakeData(pot_nd).ToTH1(pot_nd);
  nd_gar_rhc_1d_all_data->SetName("nd_gar_rhc_1d_all_data");
  nd_gar_rhc_1d_all_noshift->SetName("nd_gar_rhc_1d_all_noshift");
  nd_gar_rhc_1d_all_data->Write();
  nd_gar_rhc_1d_all_noshift->Write();

  fout->Close();
  delete fout;
} // nuwroGarFit
