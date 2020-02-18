// fdFitWithRwgt.C
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

TH1* makePrefit(const char* name, PredictionInterp* pred, const SystShifts syst, 
		osc::IOscCalculatorAdjustable* osc, const double pot)
{
  TH1* h = pred->PredictSyst(osc, syst).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
  h->SetName(name);
  setHistAttr(h);
  h->SetTitle("Prefit");
  h->SetLineColor(kRed);
  return h;
}

TH1* makePostfit(const char* name, PredictionInterp* pred, const SystShifts syst, 
		 osc::IOscCalculatorAdjustable* osc, const double pot)
{
  TH1* h = pred->PredictSyst(osc, syst).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
  h->SetName(name);
  setHistAttr(h);
  h->SetTitle("Postfit");
  h->SetLineColor(kBlue);
  return h;
}

TH1* makeData(const char* name, PredictionInterp* pred, const SystShifts syst, 
	      osc::IOscCalculatorAdjustable* osc, const double pot)
{
  TH1* h = pred->PredictSyst(osc, syst).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
  h->SetName(name);
  setHistAttr(h);
  h->SetTitle("Data");
  h->SetLineColor(kBlack);
  return h;
}

THStack *makeFitStack(const char* name, const char* title, TH1* pre, TH1* data, TH1* post) 
{
  THStack *hs = new THStack(name, title);
  hs->Add(pre);
  hs->Add(data);
  hs->Add(post);
  return hs;
}

THStack *dataMCWgtComp(const char* name, const char* title, 
		       PredictionInterp* mcpred, 
		       PredictionInterp* mcpred_nowgt,
		       PredictionInterp *datapred, const SystShifts datasyst, 
		       osc::IOscCalculatorAdjustable* osc, const double pot)
{
  THStack *hs = new THStack(name, title);
  TH1 *hdata     = makeData(Form("%s_data", name), datapred, datasyst, osc, pot);
  TH1 *hmc       = mcpred->Predict(osc).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
  TH1 *hmc_nowgt = mcpred_nowgt->Predict(osc).FakeData(pot).ToTH1(pot, kPOT, kBinDensity);
  setHistAttr(hmc);
  setHistAttr(hmc_nowgt);
  hmc->SetName(Form("%s_mc", name));
  hmc->SetTitle("Weighted MC");
  hmc->SetLineColor(kRed);
  hmc_nowgt->SetName(Form("%s_mcnowgt", name));
  hmc_nowgt->SetTitle("Unweighted MC");
  hmc_nowgt->SetLineColor(kBlue);
  hs->Add(hmc);
  hs->Add(hmc_nowgt);
  hs->Add(hdata);
  return hs;
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

void fdFitWithRwgt(const char* outfile, 
		   const int firstPoint=0, const int lastPoint=20, const int nPoints=20,
		   bool makeFDInterps=false,
		   const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2",
		   const char* cafs="/dune/data/users/marshalc/CAFs/nuwroRW/")
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

  std::vector<const ISyst*> fitsysts_norwt = GetListOfSysts(true, true, 
							    true, false, true,
							    false, false,
							    false, 20, 
							    true);
  RemoveSysts(fitsysts_norwt, {"NuWroRWT", "NuWroReweightFakeData"});
  std::vector<const ISyst*> ndwgt = GetXSecSysts({"NuWroRWT"});
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  // Sets up the systematics correctly
  SystShifts ndshift(ndwgt.at(ndwgt.size()-1), 1); // ND derived weight
  // Little bit of printout to make sure everything is alright
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"fakedatashift is using "<<fakedata.at(fakedata.size()-1)->ShortName()<<std::endl;
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  std::cout<<"Loaded "<<systlist.size()<<" systs to make the PredictionInterps of which "<<fitsysts_norwt.size()<<" are being fitted"<<std::endl;
  std::cout<<"Systs to make state files"<<std::endl;
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist[i]->ShortName()<<std::endl;
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
  if (makeFDInterps || TFile(Form("%s/state_FD_FHC_wgt.root", stateFileDir)).IsZombie()) {
    std::cout<<"Making weighted FD FHC PredictionInterps"<<std::endl;
    Loaders loadersFHC;
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    PredictionInterp predNumuFHCreco_wgt(systlist, this_calc, gennumureco, loadersFHC, ndshift);
    PredictionInterp predNueFHCreco_wgt(systlist, this_calc, gennuereco, loadersFHC, ndshift); 
    loadersFHC.Go();

    TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC_wgt.root", stateFileDir), "recreate");
    predNumuFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc_wgt"));
    predNueFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc_wgt"));
    fFDfhc->Close();
    delete fFDfhc;
  }
  if (makeFDInterps || TFile(Form("%s/state_FD_RHC_wgt.root", stateFileDir)).IsZombie()) {
    std::cout<<"Remaking weighted FD RHC PredictionInterps"<<std::endl;
    Loaders loadersRHC;
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    PredictionInterp predNumuRHCreco_wgt(systlist, this_calc, gennumureco, loadersRHC, ndshift); 
    PredictionInterp predNueRHCreco_wgt(systlist, this_calc, gennuereco, loadersRHC, ndshift);
    loadersRHC.Go();

    std::cout<<"Saving weighted FD numu RHC to "<<Form("%s/state_FD_RHC_wgt.root", stateFileDir)<<std::endl;
    TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC_wgt.root", stateFileDir), "recreate");
    predNumuRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc_wgt"));
    std::cout<<"Saving weighted FD nue RHC to "<<Form("%s/state_FD_RHC_wgt.root", stateFileDir)<<std::endl;
    predNueRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc_wgt"));
    fFDrhc->Close();
    delete fFDrhc;
    std::cout<<"Made & saved all weighted FD RHC PredictionInterps"<<std::endl;
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
  std::cout<<"Loaded numu"<<std::endl;
  PredictionInterp& predNueRhcReco_wgt  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt->GetDirectory("fd_interp_nue_rhc_wgt")).release();
  std::cout<<"Loaded nue"<<std::endl;
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

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();

  TGraph *dCPTrue = new TGraph();
  TGraph *dCPBias = new TGraph();
  TGraph *dCPBias_nowgt = new TGraph();
  TGraph *dCPBias_noCorr = new TGraph();

  TGraph *theta13 = new TGraph();
  TGraph *theta23 = new TGraph();
  TGraph *theta13_nowgt = new TGraph();
  TGraph *theta23_nowgt = new TGraph();
  dCPBias_nowgt->SetLineColor(kRed);
  dCPBias_nowgt->SetMarkerColor(kRed);
  dCPBias_nowgt->SetLineWidth(2);
  dCPBias->SetLineColor(kBlack);
  dCPBias->SetLineWidth(2);

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
  double deltaFit_norwt, deltaBias_norwt, th23Fit_norwt, th13Fit_norwt, dmsq32Fit_norwt, chi2_norwt;
  trFits->Branch("deltaFit_norwt", &deltaFit_norwt);
  trFits->Branch("deltaBias_norwt", &deltaBias_norwt);
  trFits->Branch("th13Fit_norwt", &th13Fit_norwt);
  trFits->Branch("th23Fit_norwt", &th23Fit_norwt);
  trFits->Branch("dmsq32Fit_norwt", &dmsq32Fit_norwt);
  trFits->Branch("chi2_norwt", &chi2_norwt);
  double deltaFit_noTh13, deltaBias_noTh13, th23Fit_noTh13, th13Fit_noTh13, dmsq32Fit_noTh13, chi2_noTh13;
  trFits->Branch("deltaFit_noTh13", &deltaFit_noTh13);
  trFits->Branch("deltaBias_noTh13", &deltaBias_noTh13);
  trFits->Branch("th13Fit_noTh13", &th13Fit_noTh13);
  trFits->Branch("th23Fit_noTh13", &th23Fit_noTh13);
  trFits->Branch("dmsq32Fit_noTh13", &dmsq32Fit_noTh13);
  trFits->Branch("chi2_noTh13", &chi2_noTh13);
  double deltaFit_norwt_noTh13, deltaBias_norwt_noTh13, th23Fit_norwt_noTh13, th13Fit_norwt_noTh13, dmsq32Fit_norwt_noTh13, chi2_norwt_noTh13;
  trFits->Branch("deltaFit_norwt_noTh13", &deltaFit_norwt_noTh13);
  trFits->Branch("deltaBias_norwt_noTh13", &deltaBias_norwt_noTh13);
  trFits->Branch("th13Fit_norwt_noTh13", &th13Fit_norwt_noTh13);
  trFits->Branch("th23Fit_norwt_noTh13", &th23Fit_norwt_noTh13);
  trFits->Branch("dmsq32Fit_norwt_noTh13", &dmsq32Fit_norwt_noTh13);
  trFits->Branch("chi2_norwt_noTh13", &chi2_norwt_noTh13);
  double deltaFit_lar, deltaBias_lar, th23Fit_lar, th13Fit_lar, dmsq32Fit_lar, chi2_lar;
  trFits->Branch("deltaFit_lar", &deltaFit_lar);
  trFits->Branch("deltaBias_lar", &deltaBias_lar);
  trFits->Branch("th13Fit_lar", &th13Fit_lar);
  trFits->Branch("th23Fit_lar", &th23Fit_lar);
  trFits->Branch("dmsq32Fit_lar", &dmsq32Fit_lar);
  trFits->Branch("chi2_lar", &chi2_lar);

  for (int i=firstPoint; i<=lastPoint; i++) {
    double dCP = ((double)i/(double)nPoints) * 2. * TMath::Pi();
    std::cout<<"\nRunning delta CP point "<<i<<" of "<<nPoints<<". dCP="<<dCP<<std::endl;
    this_calc->SetdCP(dCP);
    dCPTrue->SetPoint(i, i, dCP);

    // True variables
    ii = i;
    deltaTrue  = dCP;
    th13True   = this_calc->GetTh13();
    th23True   = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
    dmsq32True = this_calc->GetDmsq32()*1000.;

    // Make the experiments
    std::cout<<"Building SingleSampleExperiments"<<std::endl;
    SingleSampleExperiment expNumuFhcReco(&predNumuFhcReco_wgt, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco(&predNueFhcReco_wgt, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco(&predNumuRhcReco_wgt, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco(&predNueRhcReco_wgt, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    expNumuFhcReco.SetMaskHist(0.5, 8.); 
    expNumuRhcReco.SetMaskHist(0.5, 8.);
    expNueFhcReco.SetMaskHist(0.5, 8.); 
    expNueRhcReco.SetMaskHist(0.5, 8.);

    std::cout<<"Built SingleSampleExperiments"<<std::endl;

    MultiExperiment multiExp({&expNumuFhcReco, &expNueFhcReco, &expNumuRhcReco, &expNueRhcReco, th13});
    SystShifts testSysts = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(1);
    std::vector<double> out = fitPoint(&multiExp, fitVars, fitsysts_norwt,
				       testOsc, testSysts, 
				       oscSeeds, this_calc);
    chi2      = out.at(out.size()-1);
    deltaFit  = out.at(3);
    deltaBias = out.at(out.size()-2);
    th13Fit   = out.at(1);
    dmsq32Fit = out.at(0);
    th23Fit   = out.at(2);

    // Now set up the fit itself                             
    SingleSampleExperiment expNumuFhcReco_norwt(&predNumuFhcReco, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_norwt(&predNueFhcReco, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_norwt(&predNumuRhcReco, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_norwt(&predNueRhcReco, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));     
    expNumuFhcReco_norwt.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_norwt.SetMaskHist(0.5, 8.);
    expNueFhcReco_norwt.SetMaskHist(0.5, 8.); 
    expNueRhcReco_norwt.SetMaskHist(0.5, 8.);

    osc::IOscCalculatorAdjustable* testOsc_norwt = NuFitOscCalc(1);
    SystShifts testSysts_norwt = kNoShift;
    // Prefits
    TH1* hpreNumuFhc_nowgt = makePrefit(Form("hpreNumuFhc_nowgt_%d", i), &predNumuFhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNueFhc_nowgt  = makePrefit(Form("hpreNueFhc_nowgt_%d", i), &predNueFhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNumuRhc_nowgt = makePrefit(Form("hpreNumuRhc_nowgt_%d", i), &predNumuRhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNueRhc_nowgt  = makePrefit(Form("hpreNueRhc_nowgt_%d", i), &predNueRhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    hpreNumuFhc_nowgt->Write();
    hpreNueFhc_nowgt->Write();
    hpreNumuRhc_nowgt->Write();
    hpreNueRhc_nowgt->Write();

    MultiExperiment multiExp_norwt({&expNumuFhcReco_norwt, &expNueFhcReco_norwt, &expNumuRhcReco_norwt, &expNueRhcReco_norwt, th13});
    std::vector<double> out_norwt = fitPoint(&multiExp_norwt, fitVars, fitsysts_norwt,
					     testOsc_norwt, testSysts_norwt, 
					     oscSeeds, this_calc);
    chi2_norwt      = out_norwt.at(out_norwt.size()-1);
    deltaFit_norwt  = out_norwt.at(3);
    deltaBias_norwt = out_norwt.at(out_norwt.size()-2);
    th13Fit_norwt   = out_norwt.at(1);
    dmsq32Fit_norwt = out_norwt.at(0);
    th23Fit_norwt   = out_norwt.at(2);

    // Data
    TH1* hdataNumuFhc_nowgt = makeData(Form("hdataNumuFhc_nowgt_%d", i), &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueFhc_nowgt  = makeData(Form("hdataNueFhc_nowgt_%d", i), &predNueFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNumuRhc_nowgt = makeData(Form("hdataNumuRhc_nowgt_%d", i), &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueRhc_nowgt  = makeData(Form("hdataNueRhc_nowgt_%d", i), &predNueRhcReco, fakedatashift, this_calc, pot_fd);
    hdataNumuFhc_nowgt->Write();
    hdataNueFhc_nowgt->Write();
    hdataNumuRhc_nowgt->Write();
    hdataNueRhc_nowgt->Write();

    // Postfits
    TH1* hpostNumuFhc_nowgt = makePostfit(Form("hpostNumuFhc_nowgt_%d", i), &predNumuFhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNueFhc_nowgt  = makePostfit(Form("hpostNueFhc_nowgt_%d", i), &predNueFhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNumuRhc_nowgt = makePostfit(Form("hpostNumuRhc_nowgt_%d", i), &predNumuRhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNueRhc_nowgt  = makePostfit(Form("hpostNueRhc_nowgt_%d", i), &predNueRhcReco, testSysts_norwt, testOsc_norwt, pot_fd);
    hpostNumuFhc_nowgt->Write();
    hpostNueFhc_nowgt->Write();
    hpostNumuRhc_nowgt->Write();
    hpostNueRhc_nowgt->Write();

    THStack *hsNumuFhc_nowgt = makeFitStack(Form("hsNumuFhc_nowgt_%d", i), "#nu_{#mu} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuFhc_nowgt, hdataNumuFhc_nowgt, hpostNumuFhc_nowgt);
    THStack *hsNueFhc_nowgt  = makeFitStack(Form("hsNueFhc_nowgt_%d", i), "#nu_{e} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueFhc_nowgt, hdataNueFhc_nowgt, hpostNueFhc_nowgt);
    THStack *hsNumuRhc_nowgt = makeFitStack(Form("hsNumuRhc_nowgt_%d", i), "#nu_{#mu} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuRhc_nowgt, hdataNumuRhc_nowgt, hpostNumuRhc_nowgt);
    THStack *hsNueRhc_nowgt  = makeFitStack(Form("hsNueRhc_nowgt_%d", i), "#nu_{e} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueRhc_nowgt, hdataNueRhc_nowgt, hpostNueRhc_nowgt);
    hsNumuFhc_nowgt->Write();
    hsNueFhc_nowgt->Write();
    hsNumuRhc_nowgt->Write();
    hsNueRhc_nowgt->Write();

    // theta13_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), this_calc->GetTh13()-fPostFitValues.at(1));
    // theta23_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23())-fPostFitValues.at(2));

    // Fit with no NuFit penalty
    MultiExperiment multiExp_noTh13({&expNumuFhcReco, &expNueFhcReco, &expNumuRhcReco, &expNueRhcReco});
    SystShifts testSysts_noTh13 = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_noTh13 = NuFitOscCalc(1);
    std::vector<double> out_noTh13 = fitPoint(&multiExp_noTh13, fitVars, fitsysts_norwt,
					      testOsc_noTh13, testSysts_noTh13, 
					      oscSeeds, this_calc);
    chi2_noTh13      = out_noTh13.at(out_noTh13.size()-1);
    deltaFit_noTh13  = out_noTh13.at(3);
    deltaBias_noTh13 = out_noTh13.at(out_noTh13.size()-2);
    th13Fit_noTh13   = out_noTh13.at(1);
    dmsq32Fit_noTh13 = out_noTh13.at(0);
    th23Fit_noTh13   = out_noTh13.at(2);

    // Fit with no NuFit penalty and no reweight
    MultiExperiment multiExp_norwt_noTh13({&expNumuFhcReco_norwt, &expNueFhcReco_norwt, &expNumuRhcReco_norwt, &expNueRhcReco_norwt});
    SystShifts testSysts_norwt_noTh13 = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_norwt_noTh13 = NuFitOscCalc(1);
    std::vector<double> out_norwt_noTh13 = fitPoint(&multiExp_norwt_noTh13, fitVars, fitsysts_norwt,
						    testOsc_norwt_noTh13, testSysts_norwt_noTh13, 
						    oscSeeds, this_calc);
    chi2_norwt_noTh13      = out_norwt_noTh13.at(out_norwt_noTh13.size()-1);
    deltaFit_norwt_noTh13  = out_norwt_noTh13.at(3);
    deltaBias_norwt_noTh13 = out_norwt_noTh13.at(out_norwt_noTh13.size()-2);
    th13Fit_norwt_noTh13   = out_norwt_noTh13.at(1);
    dmsq32Fit_norwt_noTh13 = out_norwt_noTh13.at(0);
    th23Fit_norwt_noTh13   = out_norwt_noTh13.at(2);

    SingleSampleExperiment expNumuFhcReco_lar(&predNumuFhcReco_wgt_lar, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_lar(&predNueFhcReco_wgt_lar, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_lar(&predNumuRhcReco_wgt_lar, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_lar(&predNueRhcReco_wgt_lar, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    expNumuFhcReco_lar.SetMaskHist(0.5, 8.); 
    expNumuRhcReco_lar.SetMaskHist(0.5, 8.);
    expNueFhcReco_lar.SetMaskHist(0.5, 8.); 
    expNueRhcReco_lar.SetMaskHist(0.5, 8.);
    MultiExperiment multiExp_lar({&expNumuFhcReco_lar, &expNueFhcReco_lar, &expNumuRhcReco_lar, &expNueRhcReco_lar});
    SystShifts testSysts_lar = kNoShift;
    osc::IOscCalculatorAdjustable* testOsc_lar = NuFitOscCalc(1);
    std::vector<double> out_lar = fitPoint(&multiExp_lar, fitVars, fitsysts_norwt,
					   testOsc_lar, testSysts_lar, 
					   oscSeeds, this_calc);
    chi2_lar      = out_lar.at(out_lar.size()-1);
    deltaFit_lar  = out_lar.at(3);
    deltaBias_lar = out_lar.at(out_lar.size()-2);
    th13Fit_lar   = out_lar.at(1);
    dmsq32Fit_lar = out_lar.at(0);
    th23Fit_lar   = out_lar.at(2);

    trFits->Fill();
  }
  dCPTrue->SetTitle("True #delta_{CP}; Point; #delta_{CP}");
  dCPBias->SetTitle("#delta_{CP} bias with NuWroRWT; True #delta_{CP}; #delta_{CP} bias");
  dCPBias_nowgt->SetTitle("#delta_{CP} bias without NuWroRWT; True #delta_{CP}; #delta_{CP} bias");
  dCPTrue->Write("dCPTrue");
  dCPBias->Write("dCPBias");
  dCPBias_nowgt->Write("dCPBias_nowgt");
  dCPBias_noCorr->Write("dCPBias_noCorr");

  theta13->Write("theta13");
  theta13_nowgt->Write("theta13_nowgt");
  theta23->Write("theta23");
  theta23_nowgt->Write("theta23_nowgt");

  trFits->Write();

  fout->Close();
} // fdFitWithRwgt
		  
