// fdFitsWithNuWro.C
#include "StandardRecord/StandardRecord.h"
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
// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"

//#include "Utilities/rootlogon.C"

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
const Var kTrueEnergy = SIMPLEVAR(Ev);
const Var Enu_reco_numu = SIMPLEVAR(Ev_reco_numu);
const Var Enu_reco_nue  = SIMPLEVAR(Ev_reco_nue);
// Binnings
std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
const Binning binsEreco  = Binning::Custom(binEEdges);
// Axes
const HistAxis axisnumu("Reco #nu energy (GeV)", binsEreco, Enu_reco_numu);
const HistAxis axisnue("Reco #nu energy (GeV)", binsEreco, Enu_reco_nue);
// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
// Oscillation variables to fit
std::vector<const IFitVar*> fitVars = {&kFitDmSq32, &kFitTheta13, &kFitSinSqTheta23, 
				       &kFitDeltaInPiUnits, &kFitRho};

// Seeds for oscillation parameters
// Seed both octants for theta23
std::map<const IFitVar*, std::vector<double>> oscSeeds = {};
std::vector<double> theta23Seeds = {TMath::Sin(52.4*TMath::Pi()/180)*TMath::Sin(52.4*TMath::Pi()/180),
				    TMath::Sin(40.3*TMath::Pi()/180)*TMath::Sin(40.3*TMath::Pi()/180)};
// Seed both hierarchies
std::vector<double> dmsq32Seeds  = {2.525e-3 - 7.39e-5, -2.512e-3}; 
// Seed minimal and maximal delta
std::vector<double> deltaCPSeeds = {0., 1.5}; // Pi units

NuFitPenalizer penalty; // NuFit penalizer

void fdFitsWithNuWro(const char* outfile,
		     const int firstPoint=0, const int lastPoint=20, const int nPoints=20,
		     const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2")
{
  gROOT->SetBatch(kTRUE);
  // rootlogon();

  assert(firstPoint>=0 && firstPoint<=lastPoint && lastPoint<=nPoints);

  osc::IOscCalcAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true);

  std::vector<const ISyst*> fitsysts_norwt = GetListOfSysts(true, true, 
							    true, false, true,
							    false, false,
							    false, 20, 
							    true);
  RemoveSysts(fitsysts_norwt, {"NuWroReweightFakeData"});
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  // Little bit of printout to make sure everything is alright
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"fakedatashift is using "<<fakedata.at(fakedata.size()-1)->ShortName()<<std::endl;
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());

  std::cout<<"Loaded "<<fitsysts_norwt.size()<<" systs to be fitted"<<std::endl;
  std::cout<<"\nSysts to be fitted"<<std::endl;
  for (unsigned int i=0; i<fitsysts_norwt.size(); i++) {
    std::cout<<fitsysts_norwt[i]->ShortName()<<std::endl;
  }

  std::cout<<"Loading samples"<<std::endl;
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "read");
  assert(fFDfhc);
  PredictionInterp& predNumuFhcReco = *ana::LoadFrom<PredictionInterp>(fFDfhc, "fd_interp_numu_fhc").release();
  std::cout<<"Loaded numu FHC"<<std::endl;
  PredictionInterp& predNueFhcReco  = *ana::LoadFrom<PredictionInterp>(fFDfhc, "fd_interp_nue_fhc").release();
  fFDfhc->Close();
  std::cout<<"Loaded all FHC samples"<<std::endl;
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "read");
  assert(fFDrhc);
  PredictionInterp& predNumuRhcReco = *ana::LoadFrom<PredictionInterp>(fFDrhc, "fd_interp_numu_rhc").release();
  PredictionInterp& predNueRhcReco  = *ana::LoadFrom<PredictionInterp>(fFDrhc, "fd_interp_nue_rhc").release();
  std::cout<<"Loaded all FD samples"<<std::endl;
  fFDrhc->Close();

  TFile* fout = new TFile(outfile, "recreate");

  for (int i=firstPoint; i<=lastPoint; i++) {
    double dCP = ((double)i/(double)nPoints) * 2. * TMath::Pi();
    std::cout<<"\nRunning delta CP point "<<i<<" of "<<nPoints<<". dCP="<<dCP<<std::endl;
    this_calc->SetdCP(dCP);

    // Make the experiments
    std::cout<<"Building SingleSampleExperiments"<<std::endl;
    SingleSampleExperiment expNumuFhcReco(&predNumuFhcReco, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco(&predNueFhcReco, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco(&predNumuRhcReco, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco(&predNueRhcReco, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    expNumuFhcReco.SetMaskHist(0.5, 8.); 
    expNumuRhcReco.SetMaskHist(0.5, 8.);
    expNueFhcReco.SetMaskHist(0.5, 8.); 
    expNueRhcReco.SetMaskHist(0.5, 8.);
    std::cout<<"Built SingleSampleExperiments"<<std::endl;

    TH1* hNumuFhcReco_nuwro = predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueFhcReco_nuwro  = predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNumuRhcReco_nuwro = predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueRhcReco_nuwro  = predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    hNumuFhcReco_nuwro->SetTitle(Form("#nu_{#mu} FHC NuWro: #delta = %.2f", dCP));
    hNumuRhcReco_nuwro->SetTitle(Form("#nu_{#mu} RHC NuWro: #delta = %.2f", dCP));
    hNueFhcReco_nuwro->SetTitle(Form("#nu_{e} FHC NuWro: #delta = %.2f", dCP));
    hNueRhcReco_nuwro->SetTitle(Form("#nu_{e} RHC NuWro: #delta = %.2f", dCP));
    hNumuFhcReco_nuwro->SetLineColor(kBlack);
    hNumuRhcReco_nuwro->SetLineColor(kBlack);
    hNueFhcReco_nuwro->SetLineColor(kBlack);
    hNueRhcReco_nuwro->SetLineColor(kBlack);

    TH1* hNumuFhcReco_genie = predNumuFhcReco.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueFhcReco_genie  = predNueFhcReco.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNumuRhcReco_genie = predNumuRhcReco.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueRhcReco_genie  = predNueRhcReco.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    hNumuFhcReco_genie->SetTitle(Form("#nu_{#mu} FHC GENIE: #delta = %.2f", dCP));
    hNumuRhcReco_genie->SetTitle(Form("#nu_{#mu} RHC GENIE: #delta = %.2f", dCP));
    hNueFhcReco_genie->SetTitle(Form("#nu_{e} FHC GENIE: #delta = %.2f", dCP));
    hNueRhcReco_genie->SetTitle(Form("#nu_{e} RHC GENIE: #delta = %.2f", dCP));
    hNumuFhcReco_genie->SetLineColor(kBlue);
    hNumuRhcReco_genie->SetLineColor(kBlue);
    hNueFhcReco_genie->SetLineColor(kBlue);
    hNueRhcReco_genie->SetLineColor(kBlue);
    std::cout<<"Made data and nominal hists"<<std::endl;

    MultiExperiment multiExp({&expNumuFhcReco, &expNueFhcReco, &expNumuRhcReco, &expNueRhcReco});
    SystShifts testSysts = kNoShift;
    osc::IOscCalcAdjustable* testOsc = NuFitOscCalc(1);
    MinuitFitter fit(&multiExp, fitVars, fitsysts_norwt);
    double chisq  = fit.Fit(testOsc, testSysts, oscSeeds, {})->EvalMetricVal();

    std::cout<<"Now making postfit plots"<<std::endl;
    TH1* hNumuFhcReco_post = predNumuFhcReco.PredictSyst(testOsc, testSysts).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueFhcReco_post  = predNueFhcReco.PredictSyst(testOsc, testSysts).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNumuRhcReco_post = predNumuRhcReco.PredictSyst(testOsc, testSysts).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    TH1* hNueRhcReco_post  = predNueRhcReco.PredictSyst(testOsc, testSysts).FakeData(pot_fd).ToTH1(pot_fd, kPOT, kBinDensity);
    hNumuFhcReco_post->SetLineColor(kRed);
    hNumuRhcReco_post->SetLineColor(kRed);
    hNueFhcReco_post->SetLineColor(kRed);
    hNueRhcReco_post->SetLineColor(kRed);
    hNumuFhcReco_post->SetTitle(Form("#nu_{#mu} FHC postfit: #delta = %.2f", dCP));
    hNumuRhcReco_post->SetTitle(Form("#nu_{#mu} RHC postfit: #delta = %.2f", dCP));
    hNueFhcReco_post->SetTitle(Form("#nu_{e} FHC postfit: #delta = %.2f", dCP));
    hNueRhcReco_post->SetTitle(Form("#nu_{e} RHC postfit: #delta = %.2f", dCP));

    fout->cd();
    hNumuFhcReco_genie->Write(Form("hNumuFhcReco_genie_%d", i));
    hNumuRhcReco_genie->Write(Form("hNumuRhcReco_genie_%d", i));
    hNueFhcReco_genie->Write(Form("hNueFhcReco_genie_%d", i));
    hNueRhcReco_genie->Write(Form("hNueRhcReco_genie_%d", i));
    hNumuFhcReco_nuwro->Write(Form("hNumuFhcReco_nuwro_%d", i));
    hNumuRhcReco_nuwro->Write(Form("hNumuRhcReco_nuwro_%d", i));
    hNueFhcReco_nuwro->Write(Form("hNueFhcReco_nuwro_%d", i));
    hNueRhcReco_nuwro->Write(Form("hNueRhcReco_nuwro_%d", i));
    hNumuFhcReco_post->Write(Form("hNumuFhcReco_post_%d", i));
    hNumuRhcReco_post->Write(Form("hNumuRhcReco_post_%d", i));
    hNueFhcReco_post->Write(Form("hNueFhcReco_post_%d", i));
    hNueRhcReco_post->Write(Form("hNueRhcReco_post_%d", i));
  } // Loop over dcp points

  fout->Close();
  delete fout;

} // fdFitsWithNuWro

