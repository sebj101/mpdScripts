// fdFitWithRwgt.C
// Uses the new NuWro reweight syst (based upon the ND GAr data) to try and remove the delta CP bias
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

// Produces a vector of histograms with +/- 1, 2, 3 sigma shifts
std::vector<TH1*> makeSigmaHists(PredictionInterp* pred, const char* nameStub, const ISyst* syst, 
				 osc::IOscCalculatorAdjustable* osc, const double pot)
{
  std::vector<TH1*> hists;
  for (int i=-3; i<=3; i++) {
    SystShifts shift(syst, i);
    TH1 *h = pred->PredictSyst(osc, shift).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
    h->SetName(Form("%s_%s_%d", nameStub, syst->ShortName().c_str(), i));
    setHistAttr(h);
    hists.push_back(h);
    h->Write();
  }
  return hists;
}

THStack* makeSigmaStack(std::vector<TH1*> vec, const char* name, const char* title)
{
  THStack *hs = new THStack(name, title);
  for (unsigned int i=0; i<vec.size(); i++) {
    vec.at(i)->SetLineColor(52+i*3);
    hs->Add(vec.at(i));
  }
  return hs;
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
// Oscillation variables to fit
std::vector<const IFitVar*> fitVars = {&kFitDmSq32NHScaled, &kFitTheta13, &kFitSinSqTheta23, &kFitDeltaInPiUnits, &kFitRho};

void fdFitWithRwgt(const char* outfile, bool makeFDInterps=false,
		   const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/",
		   const char* cafs="/dune/data/users/marshalc/CAFs/nuwroRW/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

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
  // The same list of systs but don't include the reweight one
  std::vector<const ISyst*> fitsysts_norwt = GetListOfSysts(true, true, 
							    true, false, true,
							    false, false,
							    false, 20, 
							    true);
  RemoveSysts(fitsysts_norwt, {"NuWroRWT"});
  
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  // Sets up the systematics correctly

  // Little bit of printout to make sure everything is alright
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps of which "<<fitsysts.size()<<" are being fitted"<<std::endl;
  std::cout<<"Systs to make state files"<<std::endl;
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
  }
  
  std::cout<<"\nSysts to be fitted"<<std::endl;
  for (unsigned int i=0; i<fitsysts.size(); i++) {
    std::cout<<fitsysts[i]->ShortName()<<std::endl;
  }

  std::cout<<"\nSysts to be fitted no NuWro reweight"<<std::endl;
  for (unsigned int i=0; i<fitsysts_norwt.size(); i++) {
    std::cout<<fitsysts_norwt[i]->ShortName()<<std::endl;
  }

  if (makeFDInterps || TFile(Form("%s/state_FD_FHC.root", stateFileDir)).IsZombie()) {
    std::cout<<"(Re)making FHC PredictionInterps"<<std::endl;
    Loaders loadersFHC;
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    PredictionInterp predNumuFHCreco(systlist, this_calc, gennumureco, loadersFHC); 
    PredictionInterp predNueFHCreco(systlist, this_calc, gennuereco, loadersFHC); 
    loadersFHC.Go();

    std::cout<<"Saving FD numu FHC to "<<Form("%s/state_FD_FHC.root", stateFileDir)<<std::endl;
    TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "recreate");
    predNumuFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc"));
    std::cout<<"Saving FD nue FHC to "<<Form("%s/state_FD_FHC.root", stateFileDir)<<std::endl;
    predNueFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc"));
    fFDfhc->Close();
    delete fFDfhc;
  }
  if (makeFDInterps || TFile(Form("%s/state_FD_RHC.root", stateFileDir)).IsZombie()) {
    std::cout<<"(Re)making RHC PredictionInterps"<<std::endl;
    Loaders loadersRHC;
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    PredictionInterp predNumuRHCreco(systlist, this_calc, gennumureco, loadersRHC); 
    PredictionInterp predNueRHCreco(systlist, this_calc, gennuereco, loadersRHC);
    loadersRHC.Go();

    std::cout<<"Saving FD numu RHC to "<<Form("%s/state_FD_RHC.root", stateFileDir)<<std::endl;
    TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "recreate");
    predNumuRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc"));
    std::cout<<"Saving FD nue RHC to "<<Form("%s/state_FD_RHC.root", stateFileDir)<<std::endl;
    predNueRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc"));
    fFDrhc->Close();
    delete fFDrhc;
    std::cout<<"Made & saved all PredictionInterps"<<std::endl;
  } // make PredictionInterps

  // Retrieve PredictionInterps
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
  std::cout<<"Loaded all samples"<<std::endl;

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();
  // Make some +/- 1, 2, 3 sigma plots to check everything is in order
  std::vector<TH1*> numuFhcRecoVec = makeSigmaHists(&predNumuFhcReco, "numufhc", fitsysts.at(0/*fitsysts.size()-1*/), this_calc, pot_fd);
  std::vector<TH1*> nueFhcRecoVec  = makeSigmaHists(&predNueFhcReco, "nuefhc", fitsysts.at(0/*fitsysts.size()-1*/), this_calc, pot_fd);
  std::vector<TH1*> numuRhcRecoVec = makeSigmaHists(&predNumuRhcReco, "numurhc", fitsysts.at(0/*fitsysts.size()-1*/), this_calc, pot_fd);
  std::vector<TH1*> nueRhcRecoVec  = makeSigmaHists(&predNueRhcReco, "nuerhc", fitsysts.at(0/*fitsysts.size()-1*/), this_calc, pot_fd);
  THStack *hsNumuFhcReco = makeSigmaStack(numuFhcRecoVec, "hsNumuFhcReco", "#nu_{#mu} FHC; E / GeV; Events / GeV");
  THStack *hsNueFhcReco  = makeSigmaStack(nueFhcRecoVec, "hsNueFhcReco", "#nu_{e} FHC; E / GeV; Events / GeV");
  THStack *hsNumuRhcReco = makeSigmaStack(numuRhcRecoVec, "hsNumuRhcReco", "#nu_{#mu} RHC; E / GeV; Events / GeV");
  THStack *hsNueRhcReco  = makeSigmaStack(nueRhcRecoVec, "hsNueRhcReco", "#nu_{e} RHC; E / GeV; Events / GeV");
  hsNumuFhcReco->Write();
  hsNueFhcReco->Write();
  hsNumuRhcReco->Write();
  hsNueRhcReco->Write();

  // Make the experiments
  std::cout<<"Building SingleSampleExperiments"<<std::endl;
  SingleSampleExperiment expNumuFhcReco(&predNumuFhcReco, predNumuFhcReco.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));
  SingleSampleExperiment expNueFhcReco(&predNueFhcReco, predNueFhcReco.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));
  SingleSampleExperiment expNumuRhcReco(&predNumuRhcReco, predNumuRhcReco.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));
  SingleSampleExperiment expNueRhcReco(&predNueRhcReco, predNueRhcReco.PredictSyst(this_calc, SystShifts(fakedata.at(fakedata.size()-1), 1)).FakeData(pot_nd));

  std::cout<<"Built SingleSampleExperiments"<<std::endl;

  MultiExperiment multiExp;
  multiExp.Add(&expNumuFhcReco);
  multiExp.Add(&expNueFhcReco);
  multiExp.Add(&expNumuRhcReco);
  multiExp.Add(&expNueRhcReco);

  SystShifts fitsyst = kNoShift;
  SystShifts junk = kNoShift;

  auto start_fit = std::chrono::system_clock::now();
  // Now set up the fit itself                                  
  std::map<const IFitVar*, std::vector<double>> oscSeeds = {};                                        
  osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(1);
  std::cerr << "[INFO]: Beginning fit. " << BuildLogInfoString();
  Fitter fd_only_fit(&multiExp, fitVars, fitsysts);
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

  TMatrixDSym *covar = (TMatrixDSym *)fd_only_fit.GetCovariance();
  TH2D hist_covar = TH2D(*covar);
  hist_covar.SetName("covar");                                                                         
  TH2D hist_corr = *make_corr_from_covar(&hist_covar);
  hist_covar.Write();
  hist_corr.Write();
  covar->Write("covar_mat");

  // Do the same thing without the NuWro reweight syst as a check
  SystShifts fitsyst_norwt = kNoShift;
  SystShifts junk_norwt    = kNoShift;

  auto start_fit_norwt = std::chrono::system_clock::now();
  // Now set up the fit itself                                  
  osc::IOscCalculatorAdjustable* testOsc_norwt = NuFitOscCalc(1);
  std::cerr << "[INFO]: Beginning fit. " << BuildLogInfoString();
  Fitter fd_only_fit_norwt(&multiExp, fitVars, fitsysts_norwt);
  double fd_only_chisq_norwt = fd_only_fit_norwt.Fit(testOsc_norwt, junk_norwt, oscSeeds, {}, Fitter::kVerbose);
  auto end_fit_norwt = std::chrono::system_clock::now();
  std::time_t end_fit_time_norwt = std::chrono::system_clock::to_time_t(end_fit);
  std::cerr << "[FIT]: Finished fit in "
            << std::chrono::duration_cast<std::chrono::seconds>(end_fit_norwt -
                                                                start_fit_norwt).count()
            << " s after " << fd_only_fit_norwt.GetNFCN() << " iterations "
            << BuildLogInfoString();
  TTree *fd_only_tree_norwt = new TTree("fd_only_norwt", "fd_only_norwt");
  fd_only_tree_norwt->Branch("ParamNames", &fParamNames);
  fd_only_tree_norwt->Branch("PreFitValues", &fPreFitValues);
  fd_only_tree_norwt->Branch("PreFitErrors", &fPreFitErrors);
  fd_only_tree_norwt->Branch("PostFitValues", &fPostFitValues);
  fd_only_tree_norwt->Branch("PostFitErrors", &fPostFitErrors);
  fd_only_tree_norwt->Branch("CentralValues", &fCentralValues);
  fd_only_tree_norwt->Branch("TrueValues", &fTrueValues);
  fd_only_tree_norwt->Branch("NFCN", &fNFCN);
  fd_only_tree_norwt->Branch("EDM", &fEDM);
  fd_only_tree_norwt->Branch("ChiSq", &fChiSq);
  fd_only_tree_norwt->Branch("fIsValid", &fIsValid);
  fParamNames = fd_only_fit_norwt.GetParamNames();
  fPreFitValues = fd_only_fit_norwt.GetPreFitValues();
  fPreFitErrors = fd_only_fit_norwt.GetPreFitErrors();
  fPostFitValues = fd_only_fit_norwt.GetPostFitValues();
  fPostFitErrors = fd_only_fit_norwt.GetPostFitErrors();
  fCentralValues = fd_only_fit_norwt.GetCentralValues();
  fNFCN = fd_only_fit_norwt.GetNFCN();
  fEDM = fd_only_fit_norwt.GetEDM();
  fIsValid = fd_only_fit_norwt.GetIsValid();
  fChiSq = fd_only_chisq_norwt;

  fTrueValues.resize(fParamNames.size(), 0.);
  fTrueValues.at(0) = this_calc->GetDmsq32()*1000.;
  fTrueValues.at(1) = this_calc->GetTh13();
  fTrueValues.at(2) = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
  fTrueValues.at(3) = this_calc->GetdCP();
  fTrueValues.at(4) = this_calc->GetRho();

  fd_only_tree_norwt->Fill();
  fd_only_tree_norwt->Write();

  fout->Close();
} // fdFitWithRwgt
		  
