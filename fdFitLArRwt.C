// fdFitLArRwt.C
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

void fdFitLArRwt(const char* outfile, 
		 const int firstPoint=0, const int lastPoint=20, const int nPoints=20,
		 const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2")
		 
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  assert(firstPoint>=0 && firstPoint<=lastPoint && lastPoint<=nPoints);

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

  // Retrieve PredictionInterps
  // Unweighted samples
  // FD
  std::cout<<"Loading unweighted samples"<<std::endl;
  TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "read");
  assert(fFDfhc);
  PredictionInterp& predNumuFhcReco = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_numu_fhc")).release();
  PredictionInterp& predNueFhcReco  = *ana::LoadFrom<PredictionInterp>(fFDfhc->GetDirectory("fd_interp_nue_fhc")).release();
  fFDfhc->Close();
  delete fFDfhc;
  std::cout<<"Loaded all FHC samples"<<std::endl;
  TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "read");
  assert(fFDrhc);
  PredictionInterp& predNumuRhcReco = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_numu_rhc")).release();
  PredictionInterp& predNueRhcReco  = *ana::LoadFrom<PredictionInterp>(fFDrhc->GetDirectory("fd_interp_nue_rhc")).release();
  fFDrhc->Close();
  delete fFDrhc;
  std::cout<<"Loaded unweighted samples"<<std::endl;

  // Load LAr reweighted samples (separated by final state)
  std::cout<<"Loading weighted samples"<<std::endl;
  TFile *fFDfhc_lar_wgt = new TFile(Form("%s/state_FD_FHC_lar_wgt.root", stateFileDir), "read");
  assert(fFDfhc_lar_wgt);
  PredictionInterp& predNumuFhcReco_lar_wgt = *ana::LoadFrom<PredictionInterp>(fFDfhc_lar_wgt->GetDirectory("fd_interp_numu_fhc_lar_wgt")).release();
  PredictionInterp& predNueFhcReco_lar_wgt  = *ana::LoadFrom<PredictionInterp>(fFDfhc_lar_wgt->GetDirectory("fd_interp_nue_fhc_lar_wgt")).release();
  fFDfhc_lar_wgt->Close();
  delete fFDfhc_lar_wgt;
  std::cout<<"Loaded all FHC samples"<<std::endl;
  TFile *fFDrhc_lar_wgt = new TFile(Form("%s/state_FD_RHC_lar_wgt.root", stateFileDir), "read");
  assert(fFDrhc_lar_wgt);
  PredictionInterp& predNumuRhcReco_lar_wgt = *ana::LoadFrom<PredictionInterp>(fFDrhc_lar_wgt->GetDirectory("fd_interp_numu_rhc_lar_wgt")).release();
  PredictionInterp& predNueRhcReco_lar_wgt  = *ana::LoadFrom<PredictionInterp>(fFDrhc_lar_wgt->GetDirectory("fd_interp_nue_rhc_lar_wgt")).release();
  fFDrhc_lar_wgt->Close();
  delete fFDrhc_lar_wgt;
  std::cout<<"Loaded weighted samples"<<std::endl;

  // Load CC weighted samples
  std::cout<<"Loading CC weighted samples"<<std::endl;
  TFile *fFDfhc_lar_wgt_cc = new TFile(Form("%s/state_FD_FHC_lar_wgt_cc.root",stateFileDir),"read");
  assert(fFDfhc_lar_wgt_cc);
  PredictionInterp& predNumuFhcReco_lar_wgt_cc = *ana::LoadFrom<PredictionInterp>(fFDfhc_lar_wgt_cc->GetDirectory("fd_interp_numu_fhc_wgt")).release();
  PredictionInterp& predNueFhcReco_lar_wgt_cc  = *ana::LoadFrom<PredictionInterp>(fFDfhc_lar_wgt_cc->GetDirectory("fd_interp_nue_fhc_wgt")).release();
  fFDfhc_lar_wgt_cc->Close();
  delete fFDfhc_lar_wgt_cc;
  std::cout<<"Loaded all FHC samples"<<std::endl;
  TFile *fFDrhc_lar_wgt_cc = new TFile(Form("%s/state_FD_RHC_lar_wgt_cc.root",stateFileDir),"read");
  assert(fFDrhc_lar_wgt_cc);
  PredictionInterp& predNumuRhcReco_lar_wgt_cc = *ana::LoadFrom<PredictionInterp>(fFDrhc_lar_wgt_cc->GetDirectory("fd_interp_numu_rhc_lar_wgt")).release();
  std::cout<<"Got numu"<<std::endl;
  PredictionInterp& predNueRhcReco_lar_wgt_cc  = *ana::LoadFrom<PredictionInterp>(fFDrhc_lar_wgt_cc->GetDirectory("fd_interp_nue_rhc_lar_wgt")).release();
  std::cout<<"Got nue"<<std::endl;
  fFDrhc_lar_wgt_cc->Close();
  delete fFDrhc_lar_wgt_cc;
  std::cout<<"Loaded CC weighted samples"<<std::endl;

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
  double deltaFit_wgt, deltaBias_wgt, th23Fit_wgt, th13Fit_wgt, dmsq32Fit_wgt, chi2_wgt;
  trFits->Branch("deltaFit_wgt", &deltaFit_wgt);
  trFits->Branch("deltaBias_wgt", &deltaBias_wgt);
  trFits->Branch("th13Fit_wgt", &th13Fit_wgt);
  trFits->Branch("th23Fit_wgt", &th23Fit_wgt);
  trFits->Branch("dmsq32Fit_wgt", &dmsq32Fit_wgt);
  trFits->Branch("chi2_wgt", &chi2_wgt);
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
    expNumuFhcReco.SetMaskHist(0.5, 8.); 
    expNumuRhcReco.SetMaskHist(0.5, 8.);
    expNueFhcReco.SetMaskHist(0.5, 8.); 
    expNueRhcReco.SetMaskHist(0.5, 8.);

    MultiExperiment multi({&expNumuFhcReco, &expNumuRhcReco, &expNueFhcReco, &expNueRhcReco});
    SystShifts testSysts = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(1);
    std::vector<double> out = fitPoint(&multi, fitVars, fitsysts, testOsc, testSysts, 
    				       oscSeeds, this_calc);
    chi2      = out.at(out.size()-1);
    deltaFit  = out.at(3);
    deltaBias = out.at(out.size()-2);
    th13Fit   = out.at(1);
    dmsq32Fit = out.at(0);
    th23Fit   = out.at(2);

    // Full weighting
    SingleSampleExperiment expNumuFhcReco_wgt(&predNumuFhcReco_lar_wgt, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_wgt(&predNueFhcReco_lar_wgt, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_wgt(&predNumuRhcReco_lar_wgt, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_wgt(&predNueRhcReco_lar_wgt, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    expNumuFhcReco_wgt.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_wgt.SetMaskHist(0.5, 8.);
    expNueFhcReco_wgt.SetMaskHist(0.5, 8.); 
    expNueRhcReco_wgt.SetMaskHist(0.5, 8.);

    MultiExperiment multi_wgt({&expNumuFhcReco_wgt, &expNumuRhcReco_wgt, &expNueFhcReco_wgt, &expNueRhcReco_wgt});
    SystShifts testSysts_wgt = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_wgt = NuFitOscCalc(1);
    std::vector<double> out_wgt = fitPoint(&multi_wgt, fitVars, fitsysts, 
						testOsc_wgt, testSysts_wgt, 
						oscSeeds, this_calc);
    chi2_wgt      = out_wgt.at(out_wgt.size()-1);
    deltaFit_wgt  = out_wgt.at(3);
    deltaBias_wgt = out_wgt.at(out_wgt.size()-2);
    th13Fit_wgt   = out_wgt.at(1);
    dmsq32Fit_wgt = out_wgt.at(0);
    th23Fit_wgt   = out_wgt.at(2);

    // CC inc. weighting
    SingleSampleExperiment expNumuFhcReco_wgt_cc(&predNumuFhcReco_lar_wgt_cc, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_wgt_cc(&predNueFhcReco_lar_wgt_cc, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_wgt_cc(&predNumuRhcReco_lar_wgt_cc, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_wgt_cc(&predNueRhcReco_lar_wgt_cc, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    expNumuFhcReco_wgt_cc.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_wgt_cc.SetMaskHist(0.5, 8.);
    expNueFhcReco_wgt_cc.SetMaskHist(0.5, 8.); 
    expNueRhcReco_wgt_cc.SetMaskHist(0.5, 8.);

    MultiExperiment multi_wgt_cc({&expNumuFhcReco_wgt_cc, &expNumuRhcReco_wgt_cc, &expNueFhcReco_wgt_cc, &expNueRhcReco_wgt_cc});
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
  } // Loop over points

  fout->cd();
  trFits->Write();

  fout->Close();
  delete fout;
} // fdFitLArRwt
