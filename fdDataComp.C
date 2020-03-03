// fdDataComp.C
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

// Variables
const Var kTrueEnergy    = SIMPLEVAR(dune.Ev);
const Var kTrueLepEnergy = SIMPLEVAR(dune.LepE);
const Var kNPi = SIMPLEVAR(dune.nipip) + SIMPLEVAR(dune.nipim) + SIMPLEVAR(dune.nipi0);
const Var kQ2  = SIMPLEVAR(dune.Q2);
const Var kW   = SIMPLEVAR(dune.W);
const Var kFHC = SIMPLEVAR(dune.isFHC);
const Var Enu_reco_numu = SIMPLEVAR(dune.Ev_reco_numu);
const Var Enu_reco_nue  = SIMPLEVAR(dune.Ev_reco_nue);
const Var isCC = SIMPLEVAR(dune.isCC);
// Binnings 
const Binning simpleWBins  = Binning::Simple(20, 0., 3.);
const Binning simpleQ2Bins = Binning::Simple(20, 0., 3.);
std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
const Binning binsEreco  = Binning::Custom(binEEdges);
// Axes
const HistAxis axQ2("Q^{2} (GeV)^{2}", simpleQ2Bins, kQ2);
const HistAxis axW("W (GeV)", simpleWBins, kW);
const HistAxis axisnumu("Reco #nu energy (GeV)", binsEreco, Enu_reco_numu);
const HistAxis axisnue("Reco #nu energy (GeV)", binsEreco, Enu_reco_nue);
// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;

const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2/";

THStack *makeBigStacks(const char* name, const char* title, 		    
		       PredictionInterp* mcpred_wgt,
		       PredictionInterp* mcpred_lar_wgt, 
		       PredictionInterp *datapred, const SystShifts datasyst, 
		       osc::IOscCalculatorAdjustable* osc, const double pot)
{
  THStack *hs = new THStack(name, title);
  TH1 *hmc_wgt     = mcpred_wgt->Predict(osc).FakeData(pot).ToTH1(pot, kRed, kSolid, kPOT, kBinDensity);
  TH1 *hmc_lar_wgt = mcpred_lar_wgt->Predict(osc).FakeData(pot).ToTH1(pot, kBlue, kSolid, kPOT, kBinDensity);
  TH1 *hdata       = datapred->PredictSyst(osc, datasyst).FakeData(pot).ToTH1(pot, kBlack, kSolid, kPOT, kBinDensity);
  TH1 *hmc_no_wgt  = datapred->Predict(osc).FakeData(pot).ToTH1(pot, kBlack, kDashed, kPOT, kBinDensity);
  setHistAttr(hmc_wgt);
  setHistAttr(hmc_lar_wgt);
  setHistAttr(hdata);
  setHistAttr(hmc_no_wgt);
  hmc_wgt->SetTitle("Full weight");
  hmc_lar_wgt->SetTitle("Single weight");
  hdata->SetTitle("Data");
  hmc_no_wgt->SetTitle("Unweighted");
  hs->Add(hmc_wgt);
  hs->Add(hmc_lar_wgt);
  hs->Add(hdata);
  hs->Add(hmc_no_wgt);
  return hs;
}

THStack *makeDataMCComp(const char* name, const char* title, 		    
			PredictionInterp* mcpred_wgt,
			PredictionInterp *datapred, const SystShifts datasyst, 
			osc::IOscCalculatorAdjustable* osc, const double pot)
{
  THStack *hs = new THStack(name, title);
  TH1 *hmc_wgt     = mcpred_wgt->Predict(osc).FakeData(pot).ToTH1(pot, kRed, kSolid, kPOT, kBinDensity);
  TH1 *hdata       = datapred->PredictSyst(osc, datasyst).FakeData(pot).ToTH1(pot, kBlack, kDashed, kPOT, kBinDensity);
  TH1 *hmc_no_wgt  = datapred->Predict(osc).FakeData(pot).ToTH1(pot, kBlue, kSolid, kPOT, kBinDensity);
  setHistAttr(hmc_wgt);
  setHistAttr(hmc_no_wgt);
  setHistAttr(hdata);
  hmc_wgt->SetTitle("Weighted");
  hdata->SetTitle("Data");
  hmc_no_wgt->SetTitle("Unweighted MC");
  hs->Add(hmc_wgt);
  hs->Add(hmc_no_wgt);
  hs->Add(hdata);
  return hs;
}

THStack *makeDataQ0Q3Comp(const char* name, const char* title, 		    
			  PredictionInterp* mcpred_wgt,
			  PredictionInterp *datapred, const SystShifts datasyst, 
			  osc::IOscCalculatorAdjustable* osc, const double pot)
{
  THStack *hs = new THStack(name, title);
  TH1 *hmc_wgt     = mcpred_wgt->Predict(osc).FakeData(pot).ToTH1(pot, kRed, kSolid, kPOT, kBinDensity);
  TH1 *hdata       = datapred->PredictSyst(osc, datasyst).FakeData(pot).ToTH1(pot, kBlack, kSolid, kPOT, kBinDensity);
  TH1 *hmc_no_wgt  = datapred->Predict(osc).FakeData(pot).ToTH1(pot, kBlue, kDashed, kPOT, kBinDensity);
  setHistAttr(hmc_wgt);
  setHistAttr(hmc_no_wgt);
  setHistAttr(hdata);
  hmc_wgt->SetTitle("q_{0}, q_{3} weight");
  hdata->SetTitle("Data");
  hmc_no_wgt->SetTitle("Unweighted");
  hs->Add(hmc_wgt);
  hs->Add(hmc_no_wgt);
  hs->Add(hdata);
  return hs;
}

void fdDataComp(const char* outfile)
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true);

  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});

  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"fakedatashift is using "<<fakedata.at(fakedata.size()-1)->ShortName()<<std::endl;

  // Retrieve PredictionInterps
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "read");
  assert(fFDfhc);
  PredictionInterp& predNumuFhcReco = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_numu_fhc")).release();
  PredictionInterp& predNueFhcReco  = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_nue_fhc")).release();
  fFDfhc->Close();
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "read");
  assert(fFDrhc);
  PredictionInterp& predNumuRhcReco = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_numu_rhc")).release();
  PredictionInterp& predNueRhcReco  = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_nue_rhc")).release();
  fFDrhc->Close();

  std::cout<<"Loading weighted FHC samples"<<std::endl;
  TFile *fFDfhc_wgt = new TFile(Form("%s/state_FD_FHC_wgt.root", stateFileDir), "read");
  assert(fFDfhc_wgt);
  PredictionInterp& predNumuFhcReco_wgt = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt->GetDirectory("fd_interp_numu_fhc_wgt")).release();
  std::cout<<"Loaded numu"<<std::endl;
  PredictionInterp& predNueFhcReco_wgt  = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt->GetDirectory("fd_interp_nue_fhc_wgt")).release();
  std::cout<<"Loaded nue"<<std::endl;
  fFDfhc_wgt->Close();
  delete fFDfhc_wgt;
  std::cout<<"Loading weighted RHC samples"<<std::endl;
  TFile *fFDrhc_wgt = new TFile(Form("%s/state_FD_RHC_wgt.root", stateFileDir), "read");
  assert(fFDrhc_wgt);
  PredictionInterp& predNumuRhcReco_wgt = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt->GetDirectory("fd_interp_numu_rhc_wgt")).release();
  PredictionInterp& predNueRhcReco_wgt  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt->GetDirectory("fd_interp_nue_rhc_wgt")).release();
  fFDrhc_wgt->Close();
  delete fFDrhc_wgt;

  std::cout<<"Loading LAr weighted FHC samples"<<std::endl;
  TFile *fFDfhc_wgt_lar = new TFile(Form("%s/state_FD_FHC_wgt_lar.root", stateFileDir), "read");
  assert(fFDfhc_wgt_lar);
  PredictionInterp& predNumuFhcReco_wgt_lar = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_lar->GetDirectory("fd_interp_numu_fhc_wgt")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueFhcReco_wgt_lar  = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_lar->GetDirectory("fd_interp_nue_fhc_wgt")).release();
  fFDfhc_wgt_lar->Close();
  delete fFDfhc_wgt_lar;

  std::cout<<"Loading LAr weighted RHC samples"<<std::endl;
  TFile *fFDrhc_wgt_lar = new TFile(Form("%s/state_FD_RHC_wgt_lar.root", stateFileDir), "read");
  assert(fFDrhc_wgt_lar);
  PredictionInterp& predNumuRhcReco_wgt_lar = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_lar->GetDirectory("fd_interp_numu_rhc_wgt")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueRhcReco_wgt_lar  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_lar->GetDirectory("fd_interp_nue_rhc_wgt")).release();
  fFDrhc_wgt_lar->Close();
  delete fFDrhc_wgt_lar;
  std::cout<<"Loaded all samples"<<std::endl;
  // Q0Q3
  std::cout<<"Loading Q0q3 weighted FHC samples"<<std::endl;
  TFile *fFDfhc_wgt_q0q3 = new TFile(Form("%s/state_FD_FHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fFDfhc_wgt_q0q3);
  PredictionInterp& predNumuFhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_q0q3->GetDirectory("fd_interp_numu_fhc_wgt_q0q3")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueFhcReco_wgt_q0q3  = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_q0q3->GetDirectory("fd_interp_nue_fhc_wgt_q0q3")).release();
  fFDfhc_wgt_q0q3->Close();
  delete fFDfhc_wgt_q0q3;

  std::cout<<"Loading Q0q3 weighted RHC samples"<<std::endl;
  TFile *fFDrhc_wgt_q0q3 = new TFile(Form("%s/state_FD_RHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fFDrhc_wgt_q0q3);
  PredictionInterp& predNumuRhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_q0q3->GetDirectory("fd_interp_numu_rhc_wgt_q0q3")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueRhcReco_wgt_q0q3  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_q0q3->GetDirectory("fd_interp_nue_rhc_wgt_q0q3")).release();
  fFDrhc_wgt_q0q3->Close();
  delete fFDrhc_wgt_q0q3;
  std::cout<<"Loaded all samples"<<std::endl;

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  // Make comparisons of the various data
  THStack *hsnumufhc = makeBigStacks("hsnumufhc", "#nu_{#mu} FHC; E_{#nu, reco} / GeV; Events / GeV",
				     &predNumuFhcReco_wgt, &predNumuFhcReco_wgt_lar, &predNumuFhcReco, 
				     fakedatashift, this_calc, pot_fd);
  THStack *hsnumurhc = makeBigStacks("hsnumurhc", "#nu_{#mu} RHC; E_{#nu, reco} / GeV; Events / GeV",
				     &predNumuRhcReco_wgt, &predNumuRhcReco_wgt_lar, &predNumuRhcReco, 
				     fakedatashift, this_calc, pot_fd);
  THStack *hsnuefhc = makeBigStacks("hsnuefhc", "#nu_{e} FHC; E_{#nu, reco} / GeV; Events / GeV",
				     &predNueFhcReco_wgt, &predNueFhcReco_wgt_lar, &predNueFhcReco, 
				     fakedatashift, this_calc, pot_fd);
  THStack *hsnuerhc = makeBigStacks("hsnuerhc", "#nu_{e} RHC; E_{#nu, reco} / GeV; Events / GeV",
				     &predNueRhcReco_wgt, &predNueRhcReco_wgt_lar, &predNueRhcReco, 
				     fakedatashift, this_calc, pot_fd);

  THStack *hsnumufhc_q0q3 = makeDataQ0Q3Comp("hsnumufhc_q0q3", "#nu_{#mu} FHC; E_{#nu, reco} / GeV; Events / GeV", &predNumuFhcReco_wgt_q0q3, &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
  THStack *hsnumurhc_q0q3 = makeDataQ0Q3Comp("hsnumurhc_q0q3", "#nu_{#mu} RHC; E_{#nu, reco} / GeV; Events / GeV", &predNumuRhcReco_wgt_q0q3, &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
  THStack *hsnuefhc_q0q3 = makeDataQ0Q3Comp("hsnuefhc_q0q3", "#nu_{e} FHC; E_{#nu, reco} / GeV; Events / GeV", &predNueFhcReco_wgt_q0q3, &predNueFhcReco, fakedatashift, this_calc, pot_fd);
  THStack *hsnuerhc_q0q3 = makeDataQ0Q3Comp("hsnuerhc_q0q3", "#nu_{e} RHC; E_{#nu, reco} / GeV; Events / GeV", &predNueRhcReco_wgt_q0q3, &predNueRhcReco, fakedatashift, this_calc, pot_fd);

  hsnumufhc->Write();
  hsnumurhc->Write();
  hsnuefhc->Write();
  hsnuerhc->Write();

  hsnumufhc_q0q3->Write();
  hsnumurhc_q0q3->Write();
  hsnuefhc_q0q3->Write();
  hsnuerhc_q0q3->Write();

  for (int i=0; i<12; i++) {
    double dcp=(double)i/12. * 2. * TMath::Pi();
    this_calc->SetdCP(dcp);
    THStack *hsnumufhc_wgt = makeDataMCComp(Form("hsnumufhc_wgt_%d", i), "#nu_{#mu} FHC; E_{#nu, reco} / GeV; Events / GeV", &predNumuFhcReco_wgt, &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
    THStack *hsnumurhc_wgt = makeDataMCComp(Form("hsnumurhc_wgt_%d", i), "#nu_{#mu} RHC; E_{#nu, reco} / GeV; Events / GeV", &predNumuRhcReco_wgt, &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
    THStack *hsnuefhc_wgt = makeDataMCComp(Form("hsnuefhc_wgt_%d", i), "#nu_{e} FHC; E_{#nu, reco} / GeV; Events / GeV", &predNueFhcReco_wgt, &predNueFhcReco, fakedatashift, this_calc, pot_fd);
    THStack *hsnuerhc_wgt = makeDataMCComp(Form("hsnuerhc_wgt_%d", i), "#nu_{e} RHC; E_{#nu, reco} / GeV; Events / GeV", &predNueRhcReco_wgt, &predNueRhcReco, fakedatashift, this_calc, pot_fd);

    hsnumufhc_wgt->Write();
    hsnumurhc_wgt->Write();
    hsnuefhc_wgt->Write();
    hsnuerhc_wgt->Write();
  }

  fout->Close();
} // fdDataComp
