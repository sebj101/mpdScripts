// LArSamples.C
// Print out LAr spectra with and without the NuWro shifts
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
#include "CAFAna/Core/Ratio.h"
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

// Vars
const Var kRecoEnergyND  = SIMPLEVAR(dune.Ev_reco);
const Var kRecoYND       = (SIMPLEVAR(dune.Ev_reco)-SIMPLEVAR(dune.Elep_reco))/SIMPLEVAR(dune.Ev_reco);
const Var Enu_reco_numu  = SIMPLEVAR(dune.Ev_reco_numu);
const Var Enu_reco_nue   = SIMPLEVAR(dune.Ev_reco_nue);
// Axes 
std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};
const Binning binsEreco  = Binning::Custom(binEEdges);
const Binning binsYreco  = Binning::Custom(binYEdges);
const HistAxis axisnumu("Reco #nu energy (GeV)", binsEreco, Enu_reco_numu);
const HistAxis axisnue("Reco #nu energy (GeV)", binsEreco, Enu_reco_nue);
// POT for 3.5 years
const double pot_fd = 3.5 * POT120 * 40/1.13;
const double pot_nd = 3.5 * POT120;

void LArSamples(const char* outFile, 
		const char* stateFileDir="/pnfs/dune/persistent/users/sbjones/stateFiles/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, true, true,
						      false, false,
						      false, 20, 
						      true);

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  this_calc->SetdCP(3*TMath::Pi()/2.);
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"Using "<<fakedata.at(fakedata.size()-1)->ShortName()<<" as fake data shift"<<std::endl;

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
  delete fNDfhc;
  delete fNDrhc;
  delete fFDfhc;
  delete fFDrhc;

  TFile *fout = new TFile(outFile, "recreate");
  fout->cd();
  // Get the relevant spectra
  TH2 *h2_nd_lar_fhc_nuwro = predNDLArNumuFHC.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2_nd_lar_rhc_nuwro = predNDLArNumuRHC.PredictSyst(0, fakedatashift).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2_nd_lar_fhc_genie = predNDLArNumuFHC.Predict(0).FakeData(pot_nd).ToTH2(pot_nd);
  TH2 *h2_nd_lar_rhc_genie = predNDLArNumuRHC.Predict(0).FakeData(pot_nd).ToTH2(pot_nd);
  h2_nd_lar_fhc_nuwro->SetTitle("ND LAr #nu_{#mu} FHC (NuWro)");
  h2_nd_lar_rhc_nuwro->SetTitle("ND LAr #nu_{#mu} RHC (NuWro)");
  h2_nd_lar_fhc_genie->SetTitle("ND LAr #nu_{#mu} FHC");
  h2_nd_lar_rhc_genie->SetTitle("ND LAr #nu_{#mu} RHC");
  h2_nd_lar_fhc_nuwro->Write("h2_nd_lar_fhc_nuwro");
  h2_nd_lar_rhc_nuwro->Write("h2_nd_lar_rhc_nuwro");
  h2_nd_lar_fhc_genie->Write("h2_nd_lar_fhc_genie");
  h2_nd_lar_rhc_genie->Write("h2_nd_lar_rhc_genie");
  TH1D *h_nd_lar_fhc_nuwro = h2_nd_lar_fhc_nuwro->ProjectionX("h_nd_lar_fhc_nuwro", 1, h2_nd_lar_fhc_nuwro->GetNbinsY());
  TH1D *h_nd_lar_fhc_genie = h2_nd_lar_fhc_genie->ProjectionX("h_nd_lar_fhc_genie", 1, h2_nd_lar_fhc_genie->GetNbinsY());
  TH1D *h_nd_lar_rhc_nuwro = h2_nd_lar_rhc_nuwro->ProjectionX("h_nd_lar_rhc_nuwro", 1, h2_nd_lar_rhc_nuwro->GetNbinsY());
  TH1D *h_nd_lar_rhc_genie = h2_nd_lar_rhc_genie->ProjectionX("h_nd_lar_rhc_genie", 1, h2_nd_lar_rhc_genie->GetNbinsY());
  h_nd_lar_fhc_nuwro->Write();
  h_nd_lar_rhc_nuwro->Write();
  h_nd_lar_fhc_genie->Write();
  h_nd_lar_rhc_genie->Write();
  h_nd_lar_fhc_nuwro->SetLineColor(kRed);
  h_nd_lar_rhc_nuwro->SetLineColor(kRed);
  THStack *hs_nd_lar_fhc = new THStack("hs_nd_lar_fhc", "ND LAr #nu_{#mu} FHC; E_{#nu, reco} (GeV); Events");
  hs_nd_lar_fhc->Add(h_nd_lar_fhc_genie);
  hs_nd_lar_fhc->Add(h_nd_lar_fhc_nuwro);
  hs_nd_lar_fhc->Write();
  h_nd_lar_fhc_genie->Divide(h_nd_lar_fhc_nuwro);
  h_nd_lar_fhc_genie->SetTitle("GENIE/NuWro: ND LAr #nu_{#mu} FHC; E_{#nu, reco} (GeV); Ratio");
  h_nd_lar_fhc_genie->Write("h_nd_lar_fhc_ratio");
  THStack *hs_nd_lar_rhc = new THStack("hs_nd_lar_rhc", "ND LAr #nu_{#mu} RHC; E_{#nu, reco} (GeV); Events");
  hs_nd_lar_rhc->Add(h_nd_lar_rhc_genie);
  hs_nd_lar_rhc->Add(h_nd_lar_rhc_nuwro);
  hs_nd_lar_rhc->Write();
  h_nd_lar_rhc_genie->Divide(h_nd_lar_rhc_nuwro);
  h_nd_lar_rhc_genie->SetTitle("GENIE/NuWro: ND LAr #nu_{#mu} RHC; E_{#nu, reco} (GeV); Ratio");
  h_nd_lar_rhc_genie->Write("h_nd_lar_rhc_ratio");
  // 1D projections in energy
  for (int b=1; b<=h2_nd_lar_fhc_nuwro->GetNbinsY(); b++) {
    double lowerEdge = h2_nd_lar_fhc_nuwro->GetYaxis()->GetBinLowEdge(b);
    double upperEdge = lowerEdge + h2_nd_lar_fhc_nuwro->GetYaxis()->GetBinWidth(b);
    TH1D *h_nd_lar_fhc_nuwro_bin = h2_nd_lar_fhc_nuwro->ProjectionX(Form("h_nd_lar_fhc_nuwro_bin%d", b), b, b);
    TH1D *h_nd_lar_rhc_nuwro_bin = h2_nd_lar_rhc_nuwro->ProjectionX(Form("h_nd_lar_rhc_nuwro_bin%d", b), b, b);
    TH1D *h_nd_lar_fhc_genie_bin = h2_nd_lar_fhc_genie->ProjectionX(Form("h_nd_lar_fhc_genie_bin%d", b), b, b);
    TH1D *h_nd_lar_rhc_genie_bin = h2_nd_lar_rhc_genie->ProjectionX(Form("h_nd_lar_rhc_genie_bin%d", b), b, b);
    h_nd_lar_fhc_nuwro_bin->SetTitle(Form("ND LAr #nu_{#mu} FHC (NuWro) %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    h_nd_lar_rhc_nuwro_bin->SetTitle(Form("ND LAr #nu_{#mu} RHC (NuWro) %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    h_nd_lar_fhc_genie_bin->SetTitle(Form("ND LAr #nu_{#mu} FHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    h_nd_lar_rhc_genie_bin->SetTitle(Form("ND LAr #nu_{#mu} RHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    h_nd_lar_fhc_nuwro_bin->Write(Form("h_nd_lar_fhc_nuwro_bin%d", b));
    h_nd_lar_rhc_nuwro_bin->Write(Form("h_nd_lar_rhc_nuwro_bin%d", b));
    h_nd_lar_fhc_genie_bin->Write(Form("h_nd_lar_fhc_genie_bin%d", b));
    h_nd_lar_rhc_genie_bin->Write(Form("h_nd_lar_rhc_genie_bin%d", b));
    h_nd_lar_fhc_nuwro_bin->SetLineColor(kRed);
    h_nd_lar_rhc_nuwro_bin->SetLineColor(kRed);
    THStack *hs_nd_lar_fhc_bin = new THStack(Form("hs_nd_lar_fhc_bin%d", b), Form("ND LAr #nu_{#mu} FHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    hs_nd_lar_fhc_bin->Add(h_nd_lar_fhc_nuwro_bin);
    hs_nd_lar_fhc_bin->Add(h_nd_lar_fhc_genie_bin);
    hs_nd_lar_fhc_bin->Write();
    h_nd_lar_fhc_genie_bin->Divide(h_nd_lar_fhc_nuwro_bin);
    h_nd_lar_fhc_genie_bin->SetTitle(Form("GENIE/NuWro: ND LAr #nu_{#mu} FHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Ratio", lowerEdge, upperEdge));
    h_nd_lar_fhc_genie_bin->Write(Form("h_nd_lar_fhc_ratio_bin%d", b));
    THStack *hs_nd_lar_rhc_bin = new THStack(Form("hs_nd_lar_rhc_bin%d", b), Form("ND LAr #nu_{#mu} RHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Events", lowerEdge, upperEdge));
    hs_nd_lar_rhc_bin->Add(h_nd_lar_rhc_nuwro_bin);
    hs_nd_lar_rhc_bin->Add(h_nd_lar_rhc_genie_bin);
    hs_nd_lar_rhc_bin->Write();
    h_nd_lar_rhc_genie_bin->SetTitle(Form("GENIE/NuWro: ND LAr #nu_{#mu} RHC %.2g < y_{reco} < %.2g; E_{#nu, reco} (GeV); Ratio", lowerEdge, upperEdge));
    h_nd_lar_rhc_genie_bin->Divide(h_nd_lar_rhc_nuwro_bin);
    h_nd_lar_rhc_genie_bin->Write(Form("h_nd_lar_rhc_ratio_bin%d", b));
  }
  // FD numu
  TH1 *h_fd_numu_fhc_nuwro = predFDNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_numu_rhc_nuwro = predFDNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_numu_fhc_genie = predFDNumuFHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_numu_rhc_genie = predFDNumuRHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  h_fd_numu_fhc_nuwro->SetTitle("FD #nu_{#mu} FHC (NuWro)");
  h_fd_numu_rhc_nuwro->SetTitle("FD #nu_{#mu} RHC (NuWro)");
  h_fd_numu_fhc_genie->SetTitle("FD #nu_{#mu} FHC");
  h_fd_numu_rhc_genie->SetTitle("FD #nu_{#mu} RHC");
  h_fd_numu_fhc_nuwro->Write("h_fd_numu_fhc_nuwro");
  h_fd_numu_rhc_nuwro->Write("h_fd_numu_rhc_nuwro");
  h_fd_numu_fhc_genie->Write("h_fd_numu_fhc_genie");
  h_fd_numu_rhc_genie->Write("h_fd_numu_rhc_genie");
  h_fd_numu_fhc_nuwro->SetLineColor(kRed);
  h_fd_numu_rhc_nuwro->SetLineColor(kRed);
  THStack *hs_fd_numu_fhc = new THStack("hs_fd_numu_fhc", "FD #nu_{#mu} FHC; E_{#nu,reco} (GeV); Events");
  hs_fd_numu_fhc->Add(h_fd_numu_fhc_nuwro);
  hs_fd_numu_fhc->Add(h_fd_numu_fhc_genie);
  hs_fd_numu_fhc->Write();
  THStack *hs_fd_numu_rhc = new THStack("hs_fd_numu_rhc", "FD #nu_{#mu} RHC; E_{#nu,reco} (GeV); Events");
  hs_fd_numu_rhc->Add(h_fd_numu_rhc_nuwro);
  hs_fd_numu_rhc->Add(h_fd_numu_rhc_genie);
  hs_fd_numu_rhc->Write();
  // FD nue
  TH1 *h_fd_nue_fhc_nuwro = predFDNueFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_nue_rhc_nuwro = predFDNueRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_nue_fhc_genie = predFDNueFHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  TH1 *h_fd_nue_rhc_genie = predFDNueRHC.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd);
  h_fd_nue_fhc_nuwro->SetTitle("FD #nu_{e} FHC (NuWro)");
  h_fd_nue_rhc_nuwro->SetTitle("FD #nu_{e} RHC (NuWro)");
  h_fd_nue_fhc_genie->SetTitle("FD #nu_{e} FHC");
  h_fd_nue_rhc_genie->SetTitle("FD #nu_{e} RHC");
  h_fd_nue_fhc_nuwro->Write("h_fd_nue_fhc_nuwro");
  h_fd_nue_rhc_nuwro->Write("h_fd_nue_rhc_nuwro");
  h_fd_nue_fhc_genie->Write("h_fd_nue_fhc_genie");
  h_fd_nue_rhc_genie->Write("h_fd_nue_rhc_genie");
  h_fd_nue_fhc_nuwro->SetLineColor(kRed);
  h_fd_nue_rhc_nuwro->SetLineColor(kRed);
  THStack *hs_fd_nue_fhc = new THStack("hs_fd_nue_fhc", "FD #nu_{e} FHC; E_{#nu,reco} (GeV); Events");
  hs_fd_nue_fhc->Add(h_fd_nue_fhc_nuwro);
  hs_fd_nue_fhc->Add(h_fd_nue_fhc_genie);
  hs_fd_nue_fhc->Write();
  THStack *hs_fd_nue_rhc = new THStack("hs_fd_nue_rhc", "FD #nu_{e} RHC; E_{#nu,reco} (GeV); Events");
  hs_fd_nue_rhc->Add(h_fd_nue_rhc_nuwro);
  hs_fd_nue_rhc->Add(h_fd_nue_rhc_genie);
  hs_fd_nue_rhc->Write();

  // Ratios
  Ratio fd_numu_fhc_ratio(predFDNumuFHC.Predict(this_calc).FakeData(pot_fd), 
			  predFDNumuFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  Ratio fd_numu_rhc_ratio(predFDNumuRHC.Predict(this_calc).FakeData(pot_fd), 
			  predFDNumuRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  Ratio fd_nue_fhc_ratio(predFDNueFHC.Predict(this_calc).FakeData(pot_fd), 
			 predFDNueFHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  Ratio fd_nue_rhc_ratio(predFDNueRHC.Predict(this_calc).FakeData(pot_fd), 
			 predFDNueRHC.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
  TH1 *h_fd_numu_fhc_ratio = fd_numu_fhc_ratio.ToTH1();
  TH1 *h_fd_numu_rhc_ratio = fd_numu_rhc_ratio.ToTH1();
  TH1 *h_fd_nue_fhc_ratio  = fd_nue_fhc_ratio.ToTH1();
  TH1 *h_fd_nue_rhc_ratio  = fd_nue_rhc_ratio.ToTH1();
  h_fd_numu_fhc_ratio->SetTitle("FD #nu_{#mu} FHC GENIE/NuWro");
  h_fd_numu_rhc_ratio->SetTitle("FD #nu_{#mu} RHC GENIE/NuWro");
  h_fd_nue_fhc_ratio->SetTitle("FD #nu_{e} FHC GENIE/NuWro");
  h_fd_nue_rhc_ratio->SetTitle("FD #nu_{e} RHC GENIE/NuWro");
  h_fd_numu_fhc_ratio->Write("h_fd_numu_fhc_ratio");
  h_fd_numu_rhc_ratio->Write("h_fd_numu_rhc_ratio");
  h_fd_nue_fhc_ratio->Write("h_fd_nue_fhc_ratio");
  h_fd_nue_rhc_ratio->Write("h_fd_nue_rhc_ratio");

  std::cout<<"Wrote output file to "<<outFile<<std::endl;
  fout->Close();
  delete fout;
} // LArSamples
