// jointFitWithRwt.C
// Uses the new NuWro reweight syst (based upon the ND GAr data) to try and remove the delta CP bias
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

std::vector<double> fitPoint(MultiExperiment *exp, 
			     std::vector<const IFitVar*> oscVars,
			     std::vector<const ISyst*> systlist,
			     osc::IOscCalculatorAdjustable *fitOsc, SystShifts fitSyst,
			     ana::SeedList oscSeeds, osc::IOscCalculatorAdjustable *trueOsc)
{
  std::vector<double> params;
  params.resize(oscVars.size()+2, 0); // Fit values + bias and chi^2

  auto start_fit = std::chrono::system_clock::now();
  Fitter fit(exp, oscVars, systlist);
  double chisq = fit.Fit(fitOsc, fitSyst, oscSeeds, {}, Fitter::kVerbose);
  auto end_fit = std::chrono::system_clock::now();
  std::time_t end_fit_time = std::chrono::system_clock::to_time_t(end_fit);
  std::cerr << "[FIT]: Finished fit in "
	    << std::chrono::duration_cast<std::chrono::seconds>(end_fit -
								start_fit).count()                     
	    << " s after " << fit.GetNFCN() << " iterations "
	    << BuildLogInfoString();
  params.at(params.size()-1) = chisq;
  std::vector<double> fPostFitValues   = fit.GetPostFitValues();
  std::vector<std::string> fParamNames = fit.GetParamNames();
  for (unsigned int i=0; i<oscVars.size(); i++) {
    params.at(i) = fPostFitValues.at(i);
    if (fParamNames.at(i) == "delta(pi)") {
      double bias = (fPostFitValues.at(i)*TMath::Pi() - trueOsc->GetdCP())*(180./TMath::Pi());
      if (bias<-180.) bias+=360.;
      else if (bias>180.) bias-=360.;
      params.at(i) = fPostFitValues.at(i)*180.;
      params.at(params.size()-2) = bias;
    }
  }
  return params;
}

// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// Oscillation variables to fit
std::vector<const IFitVar*> fitVars = {&kFitDmSq32, &kFitTheta13, &kFitSinSqTheta23, 
				       &kFitDeltaInPiUnits/*, &kFitRho*/};

// Seeds for oscillation parameters
// Seed both octants for theta23
std::map<const IFitVar*, std::vector<double>> oscSeeds = {};
std::vector<double> theta23Seeds = {TMath::Sin(kNuFitTh23MaxNH)*TMath::Sin(kNuFitTh23MaxNH),
				    TMath::Sin(kNuFitTh23MinNH)*TMath::Sin(kNuFitTh23MinNH)};
// Seed both hierarchies
std::vector<double> dmsq32Seeds  = {kNuFitDmsq32CVNH, kNuFitDmsq32CVIH}; 
// Seed minimal and maximal delta
std::vector<double> deltaCPSeeds = {0., 1.5}; // Pi units

NuFitPenalizer penalty; // NuFit penalizer

void jointFitWithRwt(const char* outfile, 
		     const int firstPoint=0, const int lastPoint=20, const int nPoints=20,
		     const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();
  assert(firstPoint>=0 && firstPoint<=lastPoint && lastPoint<=nPoints);
  const ReactorExperiment* th13 = WorldReactorConstraint2017();

  // oscSeeds.insert({&kFitSinSqTheta23, theta23Seeds});
  // oscSeeds.insert({&kFitDmSq32, dmsq32Seeds});
  // oscSeeds.insert({&kFitDeltaInPiUnits, deltaCPSeeds});

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true);

  std::vector<const ISyst*> fitsysts = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true);
  RemoveSysts(fitsysts, {"NuWroRWT", "NuWroQ2", "NuWroW", "NuWroQ0Q3", "NuWroQ0Q3CC", "NuWroLAr", "NuWroLArCC", "NuWroRWTLar", "NuWroReweightFakeData"});
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});

  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"fakedatashift is using "<<fakedata.at(fakedata.size()-1)->ShortName()<<std::endl;
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());

  std::cout<<"Loaded "<<fitsysts.size()<<" systs to be fitted"<<std::endl;
  std::cout<<"\nSysts to be fitted"<<std::endl;
  for (unsigned int i=0; i<fitsysts.size(); i++) {
    std::cout<<fitsysts[i]->ShortName()<<std::endl;
  }

  // NuWroQ0Q3 weighted samples
  // FD
  std::cout<<"Loading Q0Q3 weighted samples"<<std::endl;
  TFile *fFDfhc_wgt_q0q3 = new TFile(Form("%s/state_FD_FHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fFDfhc_wgt_q0q3);
  PredictionInterp& predNumuFhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_q0q3->GetDirectory("fd_interp_numu_fhc_wgt_q0q3")).release();
  std::cout<<"Loaded FHC numu"<<std::endl;
  PredictionInterp& predNueFhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_q0q3->GetDirectory("fd_interp_nue_fhc_wgt_q0q3")).release();
  std::cout<<"Loaded FHC nue"<<std::endl;
  fFDfhc_wgt_q0q3->Close();
  delete fFDfhc_wgt_q0q3;
  TFile *fFDrhc_wgt_q0q3 = new TFile(Form("%s/state_FD_RHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fFDrhc_wgt_q0q3);
  PredictionInterp& predNumuRhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_q0q3->GetDirectory("fd_interp_numu_rhc_wgt_q0q3")).release();
  std::cout<<"Loaded RHC numu"<<std::endl;
  PredictionInterp& predNueRhcReco_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_q0q3->GetDirectory("fd_interp_nue_rhc_wgt_q0q3")).release();
  std::cout<<"Loaded RHC nue"<<std::endl;
  fFDrhc_wgt_q0q3->Close();
  delete fFDrhc_wgt_q0q3;
  // ND
  std::cout<<"Loading ND samples"<<std::endl;
  TFile *fNDfhc_wgt_q0q3 = new TFile(Form("%s/state_ND_FHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fNDfhc_wgt_q0q3);
  PredictionInterp& predNdFhc_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fNDfhc_wgt_q0q3->GetDirectory("nd_interp_numu_fhc_wgt")).release();
  fNDfhc_wgt_q0q3->Close();
  delete fNDfhc_wgt_q0q3;
  std::cout<<"Loaded FHC sample"<<std::endl;
  TFile *fNDrhc_wgt_q0q3 = new TFile(Form("%s/state_ND_RHC_wgt_q0q3.root", stateFileDir), "read");
  assert(fNDrhc_wgt_q0q3);
  PredictionInterp& predNdRhc_wgt_q0q3 = *ana::LoadFrom<PredictionInterp>(fNDrhc_wgt_q0q3->GetDirectory("nd_interp_numu_rhc_wgt")).release();
  fNDrhc_wgt_q0q3->Close();
  delete fNDrhc_wgt_q0q3;
  std::cout<<"Loaded RHC sample"<<std::endl;

  // Retrieve PredictionInterps
  // Unweighted samples
  // FD
  std::cout<<"Loading unweighted samples"<<std::endl;
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "read");
  assert(fFDfhc);
  PredictionInterp& predNumuFhcReco = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_numu_fhc")).release();
  PredictionInterp& predNueFhcReco  = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_nue_fhc")).release();
  fFDfhc->Close();
  std::cout<<"Loaded all FHC samples"<<std::endl;
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "read");
  assert(fFDrhc);
  PredictionInterp& predNumuRhcReco = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_numu_rhc")).release();
  PredictionInterp& predNueRhcReco  = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_nue_rhc")).release();
  fFDrhc->Close();
  std::cout<<"Loaded unweighted samples"<<std::endl;
  // ND
  std::cout<<"Loading ND samples"<<std::endl;
  TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC.root", stateFileDir), "read");
  assert(fNDfhc);
  PredictionInterp& predNdFhc = *ana::LoadFrom<PredictionInterp>(fNDfhc->GetDirectory("nd_interp_numu_fhc")).release();
  fNDfhc->Close();
  delete fNDfhc;
  std::cout<<"Loaded FHC sample"<<std::endl;
  TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC.root", stateFileDir), "read");
  assert(fNDrhc);
  PredictionInterp& predNdRhc = *ana::LoadFrom<PredictionInterp>(fNDrhc->GetDirectory("nd_interp_numu_rhc")).release();
  fNDrhc->Close();
  delete fNDrhc;
  std::cout<<"Loaded RHC sample"<<std::endl;

  // CC inc. Q0Q3 weighted samples
  // FD
  std::cout<<"Loading CC inc. weighted FHC samples"<<std::endl;
  TFile *fFDfhc_wgt_cc = new TFile(Form("%s/state_FD_FHC_wgt_cc.root", stateFileDir), "read");
  assert(fFDfhc_wgt_cc);
  PredictionInterp& predNumuFhcReco_wgt_cc = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_cc->GetDirectory("fd_interp_numu_fhc_wgt")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueFhcReco_wgt_cc  = *ana::LoadFrom<PredictionInterp>(fFDfhc_wgt_cc->GetDirectory("fd_interp_nue_fhc_wgt")).release();
  fFDfhc_wgt_cc->Close();
  delete fFDfhc_wgt_cc;
  std::cout<<"Loading CC inc. weighted RHC samples"<<std::endl;
  TFile *fFDrhc_wgt_cc = new TFile(Form("%s/state_FD_RHC_wgt_cc.root", stateFileDir), "read");
  assert(fFDrhc_wgt_cc);
  PredictionInterp& predNumuRhcReco_wgt_cc = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_cc->GetDirectory("fd_interp_numu_rhc_wgt")).release();
  std::cout<<"Loaded numu"<<std::endl; 
  PredictionInterp& predNueRhcReco_wgt_cc  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt_cc->GetDirectory("fd_interp_nue_rhc_wgt")).release();
  fFDrhc_wgt_cc->Close();
  delete fFDrhc_wgt_cc;
  // ND
  std::cout<<"Loading ND samples"<<std::endl;
  TFile *fNDfhc_wgt_cc = new TFile(Form("%s/state_ND_FHC_wgt_cc.root", stateFileDir), "read");
  assert(fNDfhc_wgt_cc);
  PredictionInterp& predNdFhc_wgt_cc = *ana::LoadFrom<PredictionInterp>(fNDfhc_wgt_cc->GetDirectory("nd_interp_numu_fhc_wgt")).release();
  fNDfhc_wgt_cc->Close();
  delete fNDfhc_wgt_cc;
  std::cout<<"Loaded FHC sample"<<std::endl;
  TFile *fNDrhc_wgt_cc = new TFile(Form("%s/state_ND_RHC_wgt_cc.root", stateFileDir), "read");
  assert(fNDrhc_wgt_cc);
  PredictionInterp& predNdRhc_wgt_cc = *ana::LoadFrom<PredictionInterp>(fNDrhc_wgt_cc->GetDirectory("nd_interp_numu_rhc_wgt")).release();
  fNDrhc_wgt_cc->Close();
  delete fNDrhc_wgt_cc;
  std::cout<<"Loaded RHC sample"<<std::endl;

  std::cout<<"Loaded all samples"<<std::endl;

  // Now do the fits at numerous delta CP points
  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  TTree *trFits = new TTree("fitResults", "fitResults");
  int ii;
  double deltaTrue, th13True, th23True, dmsq32True;
  trFits->Branch("ii", &ii);
  trFits->Branch("deltaTrue", &deltaTrue);
  trFits->Branch("th13True", &th13True);
  trFits->Branch("th23True", &th23True);
  trFits->Branch("dmsq32True", &dmsq32True);
  double deltaFit, deltaBias, th23Fit, th13Fit, dmsq32Fit, chi2;
  trFits->Branch("deltaFit", &deltaFit);
  trFits->Branch("deltaBias", &deltaBias);
  trFits->Branch("th13Fit", &th13Fit);
  trFits->Branch("th23Fit", &th23Fit);
  trFits->Branch("dmsq32Fit", &dmsq32Fit);
  trFits->Branch("chi2", &chi2);
  double deltaFit_wgt_q0q3, deltaBias_wgt_q0q3, th23Fit_wgt_q0q3, th13Fit_wgt_q0q3, dmsq32Fit_wgt_q0q3, chi2_wgt_q0q3;
  trFits->Branch("deltaFit_wgt_q0q3", &deltaFit_wgt_q0q3);
  trFits->Branch("deltaBias_wgt_q0q3", &deltaBias_wgt_q0q3);
  trFits->Branch("th13Fit_wgt_q0q3", &th13Fit_wgt_q0q3);
  trFits->Branch("th23Fit_wgt_q0q3", &th23Fit_wgt_q0q3);
  trFits->Branch("dmsq32Fit_wgt_q0q3", &dmsq32Fit_wgt_q0q3);
  trFits->Branch("chi2_wgt_q0q3", &chi2_wgt_q0q3);
  double deltaFit_wgt_cc, deltaBias_wgt_cc, th23Fit_wgt_cc, th13Fit_wgt_cc, dmsq32Fit_wgt_cc, chi2_wgt_cc;
  trFits->Branch("deltaFit_wgt_cc", &deltaFit_wgt_cc);
  trFits->Branch("deltaBias_wgt_cc", &deltaBias_wgt_cc);
  trFits->Branch("th13Fit_wgt_cc", &th13Fit_wgt_cc);
  trFits->Branch("th23Fit_wgt_cc", &th23Fit_wgt_cc);
  trFits->Branch("dmsq32Fit_wgt_cc", &dmsq32Fit_wgt_cc);
  trFits->Branch("chi2_wgt_cc", &chi2_wgt_cc);

  for (int i=firstPoint; i<=lastPoint; i++) {
    double dCP = ((double)i/(double)nPoints) * 2. * TMath::Pi();
    std::cout<<"\nRunning delta CP point "<<i<<" of "<<nPoints<<". dCP="<<dCP<<std::endl;
    this_calc->SetdCP(dCP);

    // True variables
    ii = i;
    deltaTrue  = dCP;
    th13True   = this_calc->GetTh13();
    th23True   = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
    dmsq32True = this_calc->GetDmsq32()*1000.; 

    // Make the experiments
    SingleSampleExperiment expNumuFhcReco(&predNumuFhcReco, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco(&predNueFhcReco, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco(&predNumuRhcReco, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco(&predNueRhcReco, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNdFhc(&predNdFhc, predNdFhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    SingleSampleExperiment expNdRhc(&predNdRhc, predNdRhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    expNumuFhcReco.SetMaskHist(0.5, 8.); 
    expNumuRhcReco.SetMaskHist(0.5, 8.);
    expNueFhcReco.SetMaskHist(0.5, 8.); 
    expNueRhcReco.SetMaskHist(0.5, 8.);
    expNdFhc.SetMaskHist(0.5, 10, 0, -1);
    expNdRhc.SetMaskHist(0.5, 10, 0, -1);

    // MultiExperiment multi({&expNumuFhcReco, &expNumuRhcReco, &expNueFhcReco, &expNueRhcReco, &expNdFhc, &expNdRhc});
    // SystShifts testSysts = kNoShift;
    // osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(1);
    // std::vector<double> out = fitPoint(&multi, fitVars, fitsysts, testOsc, testSysts, 
    // 				       oscSeeds, this_calc);
    chi2      = 0;//out.at(out.size()-1);
    deltaFit  = 0;//out.at(3);
    deltaBias = 0;//out.at(out.size()-2);
    th13Fit   = 0;//out.at(1);
    dmsq32Fit = 0;//out.at(0);
    th23Fit   = 0;//out.at(2);

    // Full Q0Q3 weighting
    SingleSampleExperiment expNumuFhcReco_wgt_q0q3(&predNumuFhcReco_wgt_q0q3, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_wgt_q0q3(&predNueFhcReco_wgt_q0q3, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_wgt_q0q3(&predNumuRhcReco_wgt_q0q3, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_wgt_q0q3(&predNueRhcReco_wgt_q0q3, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNdFhc_wgt_q0q3(&predNdFhc_wgt_q0q3, predNdFhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    SingleSampleExperiment expNdRhc_wgt_q0q3(&predNdRhc_wgt_q0q3, predNdRhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    expNumuFhcReco_wgt_q0q3.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_wgt_q0q3.SetMaskHist(0.5, 8.);
    expNueFhcReco_wgt_q0q3.SetMaskHist(0.5, 8.); 
    expNueRhcReco_wgt_q0q3.SetMaskHist(0.5, 8.);
    expNdFhc_wgt_q0q3.SetMaskHist(0.5, 10, 0, -1);
    expNdRhc_wgt_q0q3.SetMaskHist(0.5, 10, 0, -1);

    MultiExperiment multi_wgt_q0q3({&expNumuFhcReco_wgt_q0q3, &expNumuRhcReco_wgt_q0q3, &expNueFhcReco_wgt_q0q3, &expNueRhcReco_wgt_q0q3, &expNdFhc_wgt_q0q3, &expNdRhc_wgt_q0q3});
    SystShifts testSysts_wgt_q0q3 = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_wgt_q0q3 = NuFitOscCalc(1);
    std::vector<double> out_wgt_q0q3 = fitPoint(&multi_wgt_q0q3, fitVars, fitsysts, 
						testOsc_wgt_q0q3, testSysts_wgt_q0q3, 
						oscSeeds, this_calc);
    chi2_wgt_q0q3      = out_wgt_q0q3.at(out_wgt_q0q3.size()-1);
    deltaFit_wgt_q0q3  = out_wgt_q0q3.at(3);
    deltaBias_wgt_q0q3 = out_wgt_q0q3.at(out_wgt_q0q3.size()-2);
    th13Fit_wgt_q0q3   = out_wgt_q0q3.at(1);
    dmsq32Fit_wgt_q0q3 = out_wgt_q0q3.at(0);
    th23Fit_wgt_q0q3   = out_wgt_q0q3.at(2);

    // Q0Q3 CC inc. weighting
    SingleSampleExperiment expNumuFhcReco_wgt_cc(&predNumuFhcReco_wgt_cc, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_wgt_cc(&predNueFhcReco_wgt_cc, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_wgt_cc(&predNumuRhcReco_wgt_cc, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_wgt_cc(&predNueRhcReco_wgt_cc, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNdFhc_wgt_cc(&predNdFhc_wgt_cc, predNdFhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    SingleSampleExperiment expNdRhc_wgt_cc(&predNdRhc_wgt_cc, predNdRhc.PredictSyst(0, fakedatashift).FakeData(pot_nd));
    expNumuFhcReco_wgt_cc.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_wgt_cc.SetMaskHist(0.5, 8.);
    expNueFhcReco_wgt_cc.SetMaskHist(0.5, 8.); 
    expNueRhcReco_wgt_cc.SetMaskHist(0.5, 8.);
    expNdFhc_wgt_cc.SetMaskHist(0.5, 10, 0, -1);
    expNdRhc_wgt_cc.SetMaskHist(0.5, 10, 0, -1);

    MultiExperiment multi_wgt_cc({&expNumuFhcReco_wgt_cc, &expNumuRhcReco_wgt_cc, &expNueFhcReco_wgt_cc, &expNueRhcReco_wgt_cc, &expNdFhc_wgt_cc, &expNdRhc_wgt_cc});
    SystShifts testSysts_wgt_cc = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_wgt_cc = NuFitOscCalc(1);
    std::vector<double> out_wgt_cc = fitPoint(&multi_wgt_cc, fitVars, fitsysts, 
						testOsc_wgt_cc, testSysts_wgt_cc, 
						oscSeeds, this_calc);
    chi2_wgt_cc      = out_wgt_cc.at(out_wgt_cc.size()-1);
    deltaFit_wgt_cc  = out_wgt_cc.at(3);
    deltaBias_wgt_cc = out_wgt_cc.at(out_wgt_cc.size()-2);
    th13Fit_wgt_cc   = out_wgt_cc.at(1);
    dmsq32Fit_wgt_cc = out_wgt_cc.at(0);
    th23Fit_wgt_cc   = out_wgt_cc.at(2);

    trFits->Fill();
  } // Loop over true dcp values

  trFits->Write();
  fout->Close();
  delete fout;
} // jointFitWithRwt
