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
std::vector<double> dmsq32Seeds  = {kNuFitDmsq32CVNH};//, kNuFitDmsq32CVIH}; 
// Seed minimal and maximal delta
std::vector<double> deltaCPSeeds = {0., 0.5, 1., 1.5}; // Pi units

NuFitPenalizer penalty; // NuFit penalizer

void fdFitWithRwgt(const char* outfile, 
		   const int firstPoint=0, const int lastPoint=20, const int nPoints=20,
		   /*const int hie=1, const int oct=1,*/ bool makeFDInterps=false,
		   const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/",
		   const char* cafs="/dune/data/users/marshalc/CAFs/nuwroRW/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  assert(firstPoint>=0 && firstPoint<=lastPoint && lastPoint<=nPoints);

  oscSeeds.insert({&kFitSinSqTheta23, theta23Seeds});
  // oscSeeds.insert({&kFitDmSq32, dmsq32Seeds});
  // oscSeeds.insert({&kFitDeltaInPiUnits, deltaCPSeeds});

  osc::IOscCalculatorAdjustable* this_calc = /*DefaultOscCalc();*/NuFitOscCalc(1);
  // this_calc->SetTh13(kNuFitTh13CVNH);  
  // this_calc->SetTh23(kNuFitTh23CVNH);  
  // this_calc->SetDmsq21(kNuFitDmsq21CV);
  // this_calc->SetDmsq32(kNuFitDmsq32CVNH);
  // this_calc->SetL(kBaseline);
  // this_calc->SetRho(kEarthDensity);

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

  std::vector<const ISyst*> fitsysts_norwt = GetListOfSysts(true, true, 
							    true, false, true,
							    false, false,
							    false, 20, 
							    true);
  RemoveSysts(fitsysts_norwt, {"NuWroRWT", "NuWroReweightFakeData"});
  RemoveSysts(fitsysts, {"NuWroRWT", "NuWroReweightFakeData"});
  std::vector<const ISyst*> ndwgt = GetXSecSysts({"NuWroRWT"});
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  // Sets up the systematics correctly
  SystShifts ndshift(ndwgt.at(ndwgt.size()-1), 1); // ND derived weight
  // Little bit of printout to make sure everything is alright
  SystShifts fakedatashift(fakedata.at(fakedata.size()-1), 1);
  std::cout<<"fakedatashift is using "<<fakedata.at(fakedata.size()-1)->ShortName()<<std::endl;
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
  PredictionInterp& predNueRhcReco_wgt  = *ana::LoadFrom<PredictionInterp>(fFDrhc_wgt->GetDirectory("fd_interp_nue_rhc_wgt")).release();
  fFDrhc_wgt->Close();
  delete fFDrhc_wgt;

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

  TGraph *dCPTrue = new TGraph();
  TGraph *dCPBias = new TGraph();
  TGraph *dCPBias_nowgt = new TGraph();
  TGraph *dCPBias_noCorr = new TGraph();

  TGraph *theta13 = new TGraph();
  TGraph *theta23 = new TGraph();
  TGraph *rho     = new TGraph();
  TGraph *theta13_nowgt = new TGraph();
  TGraph *theta23_nowgt = new TGraph();
  dCPBias_nowgt->SetLineColor(kRed);
  dCPBias_nowgt->SetMarkerColor(kRed);
  dCPBias_nowgt->SetLineWidth(2);
  dCPBias->SetLineColor(kBlack);
  dCPBias->SetLineWidth(2);

  for (int i=firstPoint; i<=lastPoint; i++) {
    double dCP = ((double)i/(double)nPoints) * 2. * TMath::Pi();
    std::cout<<"\nRunning delta CP point "<<i<<" of "<<nPoints<<". dCP="<<dCP<<std::endl;
    this_calc->SetdCP(dCP);
    dCPTrue->SetPoint(i, i, dCP);
    // Make the experiments
    std::cout<<"Building SingleSampleExperiments"<<std::endl;
    SingleSampleExperiment expNumuFhcReco(&predNumuFhcReco_wgt, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco(&predNueFhcReco_wgt, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco(&predNumuRhcReco_wgt, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco(&predNueRhcReco_wgt, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
 
    std::cout<<"Built SingleSampleExperiments"<<std::endl;

    MultiExperiment multiExp({&expNumuFhcReco, &expNueFhcReco, &expNumuRhcReco, &expNueRhcReco});
    // multiExp.Add(&penalty);

    SystShifts fitsyst = kNoShift;
    SystShifts junk = kNoShift;

    auto start_fit = std::chrono::system_clock::now();
    // Now set up the fit itself                         
    osc::IOscCalculatorAdjustable* testOsc = /*DefaultOscCalc();*/NuFitOscCalc(1);
    // testOsc->SetTh13(kNuFitTh13CVNH);  
    // testOsc->SetTh23(kNuFitTh23CVNH);  
    // testOsc->SetDmsq21(kNuFitDmsq21CV);
    // testOsc->SetDmsq32(kNuFitDmsq32CVNH);
    // testOsc->SetdCP(dCP);
    // testOsc->SetL(kBaseline);
    // testOsc->SetRho(kEarthDensity);

    // Prefits
    TH1* hpreNumuFhc = makePrefit(Form("hpreNumuFhc_%d", i), &predNumuFhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpreNueFhc  = makePrefit(Form("hpreNueFhc_%d", i), &predNueFhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpreNumuRhc = makePrefit(Form("hpreNumuRhc_%d", i), &predNumuRhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpreNueRhc  = makePrefit(Form("hpreNueRhc_%d", i), &predNueRhcReco_wgt, junk, testOsc, pot_fd);
    hpreNumuFhc->Write();
    hpreNueFhc->Write();
    hpreNumuRhc->Write();
    hpreNueRhc->Write();
    // Data
    TH1* hdataNumuFhc = makeData(Form("hdataNumuFhc_%d", i), &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueFhc  = makeData(Form("hdataNueFhc_%d", i), &predNueFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNumuRhc = makeData(Form("hdataNumuRhc_%d", i), &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueRhc  = makeData(Form("hdataNueRhc_%d", i), &predNueRhcReco, fakedatashift, this_calc, pot_fd);
    hdataNumuFhc->Write();
    hdataNueFhc->Write();
    hdataNumuRhc->Write();
    hdataNueRhc->Write();

    std::cerr << "[INFO]: Beginning fit. " << BuildLogInfoString();
    Fitter fd_only_fit(&multiExp, fitVars, fitsysts_norwt);
    double fd_only_chisq = fd_only_fit.Fit(testOsc, junk, oscSeeds, {}, Fitter::kVerbose);
    auto end_fit = std::chrono::system_clock::now();
    std::time_t end_fit_time = std::chrono::system_clock::to_time_t(end_fit);
    std::cerr << "[FIT]: Finished fit in "
	      << std::chrono::duration_cast<std::chrono::seconds>(end_fit -
								  start_fit).count()
	      << " s after " << fd_only_fit.GetNFCN() << " iterations "
	      << BuildLogInfoString();

    // Postfits
    TH1* hpostNumuFhc = makePostfit(Form("hpostNumuFhc_%d", i), &predNumuFhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpostNueFhc  = makePostfit(Form("hpostNueFhc_%d", i), &predNueFhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpostNumuRhc = makePostfit(Form("hpostNumuRhc_%d", i), &predNumuRhcReco_wgt, junk, testOsc, pot_fd);
    TH1* hpostNueRhc  = makePostfit(Form("hpostNueRhc_%d", i), &predNueRhcReco_wgt, junk, testOsc, pot_fd);
    hpostNumuFhc->Write();
    hpostNueFhc->Write();
    hpostNumuRhc->Write();
    hpostNueRhc->Write();

    THStack *hsNumuFhc = makeFitStack(Form("hsNumuFhc_%d", i), "#nu_{#mu} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuFhc, hdataNumuFhc, hpostNumuFhc);
    THStack *hsNueFhc  = makeFitStack(Form("hsNueFhc_%d", i), "#nu_{e} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueFhc, hdataNueFhc, hpostNueFhc);
    THStack *hsNumuRhc = makeFitStack(Form("hsNumuRhc_%d", i), "#nu_{#mu} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuRhc, hdataNumuRhc, hpostNumuRhc);
    THStack *hsNueRhc  = makeFitStack(Form("hsNueRhc_%d", i), "#nu_{e} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueRhc, hdataNueRhc, hpostNueRhc);
    hsNumuFhc->Write();
    hsNueFhc->Write();
    hsNumuRhc->Write();
    hsNueRhc->Write();

    THStack* hsNumuFhcComp = dataMCWgtComp(Form("hsNumuFhcComp_%d", i), "#nu_{#mu} (FHC); E_{#nu, reco} / GeV; Events / GeV", &predNumuFhcReco_wgt, &predNumuFhcReco, &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
    THStack* hsNueFhcComp = dataMCWgtComp(Form("hsNueFhcComp_%d", i), "#nu_{e} (FHC); E_{#nu, reco} / GeV; Events / GeV", &predNueFhcReco_wgt, &predNueFhcReco, &predNueFhcReco, fakedatashift, this_calc, pot_fd);
    THStack* hsNumuRhcComp = dataMCWgtComp(Form("hsNumuRhcComp_%d", i), "#nu_{#mu} (RHC); E_{#nu, reco} / GeV; Events / GeV", &predNumuRhcReco_wgt, &predNumuRhcReco, &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
    THStack* hsNueRhcComp = dataMCWgtComp(Form("hsNueRhcComp_%d", i), "#nu_{e} (RHC); E_{#nu, reco} / GeV; Events / GeV", &predNueRhcReco_wgt, &predNueRhcReco, &predNueRhcReco, fakedatashift, this_calc, pot_fd);
    hsNumuFhcComp->Write();
    hsNueFhcComp->Write();
    hsNumuRhcComp->Write();
    hsNueRhcComp->Write();

    TTree *fd_only_tree = new TTree(Form("fd_only_%d",i), Form("fd_only_%d",i));
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

    double bias = (dCP - fPostFitValues.at(3)*TMath::Pi())*(180./TMath::Pi());
    if (bias<-180) bias+=360;
    else if (bias>360) bias-=360;
    dCPBias->SetPoint(i, dCP*(180./TMath::Pi()), bias);
    dCPBias_noCorr->SetPoint(i, dCP*(180./TMath::Pi()), (dCP - fPostFitValues.at(3)*TMath::Pi())*(180./TMath::Pi()));

    fTrueValues.resize(fParamNames.size(), 0.);
    fTrueValues.at(0) = this_calc->GetDmsq32()*1000.;
    fTrueValues.at(1) = this_calc->GetTh13();
    fTrueValues.at(2) = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
    fTrueValues.at(3) = this_calc->GetdCP();
    // fTrueValues.at(4) = this_calc->GetRho();

    fd_only_tree->Fill();
    fd_only_tree->Write();

    // rho->SetPoint(i, this_calc->GetRho()-fPostFitValues.at(4), dCP*(180./TMath::Pi()));
    theta13->SetPoint(i, this_calc->GetTh13()-fPostFitValues.at(1), dCP*(180./TMath::Pi()));
    theta23->SetPoint(i, TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23())-fPostFitValues.at(2), dCP*(180./TMath::Pi()));

    TMatrixDSym *covar = (TMatrixDSym *)fd_only_fit.GetCovariance();
    TH2D hist_covar = TH2D(*covar);
    hist_covar.SetName("covar");                                                                      
    TH2D hist_corr = *make_corr_from_covar(&hist_covar);
    hist_covar.Write();
    hist_corr.Write();
    covar->Write("covar_mat");

    // Do the same thing without the NuWroRWT CV weight  
    std::cout<<"\nNow making without CV weight"<<std::endl;
    SystShifts fitsyst_norwt = kNoShift;
    SystShifts junk_norwt    = kNoShift;

    // Make Experiments
    SingleSampleExperiment expNumuFhcReco_nowgt(&predNumuFhcReco, predNumuFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueFhcReco_nowgt(&predNueFhcReco, predNueFhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNumuRhcReco_nowgt(&predNumuRhcReco, predNumuRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));
    SingleSampleExperiment expNueRhcReco_nowgt(&predNueRhcReco, predNueRhcReco.PredictSyst(this_calc, fakedatashift).FakeData(pot_fd));

    std::cout<<"Built SingleSampleExperiments"<<std::endl;

    MultiExperiment multiExp_nowgt({&expNumuFhcReco_nowgt, &expNueFhcReco_nowgt, &expNumuRhcReco_nowgt, &expNueRhcReco_nowgt});
    // multiExp_nowgt.Add(&penalty);

    auto start_fit_norwt = std::chrono::system_clock::now();
    // Now set up the fit itself                                  
    osc::IOscCalculatorAdjustable* testOsc_norwt = /*DefaultOscCalc();*/NuFitOscCalc(1);
    // testOsc_norwt->SetTh13(kNuFitTh13CVNH);  
    // testOsc_norwt->SetTh23(kNuFitTh23CVNH);  
    // testOsc_norwt->SetDmsq21(kNuFitDmsq21CV);
    // testOsc_norwt->SetDmsq32(kNuFitDmsq32CVNH);
    // testOsc_norwt->SetdCP(dCP);
    // testOsc_norwt->SetL(kBaseline);
    // testOsc_norwt->SetRho(kEarthDensity);

    // Prefits
    TH1* hpreNumuFhc_nowgt = makePrefit(Form("hpreNumuFhc_nowgt_%d", i), &predNumuFhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNueFhc_nowgt  = makePrefit(Form("hpreNueFhc_nowgt_%d", i), &predNueFhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNumuRhc_nowgt = makePrefit(Form("hpreNumuRhc_nowgt_%d", i), &predNumuRhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpreNueRhc_nowgt  = makePrefit(Form("hpreNueRhc_nowgt_%d", i), &predNueRhcReco, junk_norwt, testOsc_norwt, pot_fd);
    hpreNumuFhc_nowgt->Write();
    hpreNueFhc_nowgt->Write();
    hpreNumuRhc_nowgt->Write();
    hpreNueRhc_nowgt->Write();
    // Data
    TH1* hdataNumuFhc_nowgt = makeData(Form("hdataNumuFhc_nowgt_%d", i), &predNumuFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueFhc_nowgt  = makeData(Form("hdataNueFhc_nowgt_%d", i), &predNueFhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNumuRhc_nowgt = makeData(Form("hdataNumuRhc_nowgt_%d", i), &predNumuRhcReco, fakedatashift, this_calc, pot_fd);
    TH1* hdataNueRhc_nowgt  = makeData(Form("hdataNueRhc_nowgt_%d", i), &predNueRhcReco, fakedatashift, this_calc, pot_fd);
    hdataNumuFhc_nowgt->Write();
    hdataNueFhc_nowgt->Write();
    hdataNumuRhc_nowgt->Write();
    hdataNueRhc_nowgt->Write();

    std::cerr << "[INFO]: Beginning fit. " << BuildLogInfoString();
    Fitter fd_only_fit_norwt(&multiExp_nowgt, fitVars, fitsysts_norwt);
    double fd_only_chisq_norwt = fd_only_fit_norwt.Fit(testOsc_norwt, junk_norwt, oscSeeds, {}, Fitter::kVerbose);
    auto end_fit_norwt = std::chrono::system_clock::now();
    std::time_t end_fit_time_norwt = std::chrono::system_clock::to_time_t(end_fit);
    std::cerr << "[FIT]: Finished fit in "
	      << std::chrono::duration_cast<std::chrono::seconds>(end_fit_norwt -
								  start_fit_norwt).count()
	      << " s after " << fd_only_fit_norwt.GetNFCN() << " iterations "
	      << BuildLogInfoString();

    // Postfits
    TH1* hpostNumuFhc_nowgt = makePostfit(Form("hpostNumuFhc_nowgt_%d", i), &predNumuFhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNueFhc_nowgt  = makePostfit(Form("hpostNueFhc_nowgt_%d", i), &predNueFhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNumuRhc_nowgt = makePostfit(Form("hpostNumuRhc_nowgt_%d", i), &predNumuRhcReco, junk_norwt, testOsc_norwt, pot_fd);
    TH1* hpostNueRhc_nowgt  = makePostfit(Form("hpostNueRhc_nowgt_%d", i), &predNueRhcReco, junk_norwt, testOsc_norwt, pot_fd);
    hpostNumuFhc->Write();
    hpostNueFhc->Write();
    hpostNumuRhc->Write();
    hpostNueRhc->Write();

    THStack *hsNumuFhc_nowgt = makeFitStack(Form("hsNumuFhc_nowgt_%d", i), "#nu_{#mu} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuFhc_nowgt, hdataNumuFhc_nowgt, hpostNumuFhc_nowgt);
    THStack *hsNueFhc_nowgt  = makeFitStack(Form("hsNueFhc_nowgt_%d", i), "#nu_{e} (FHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueFhc_nowgt, hdataNueFhc_nowgt, hpostNueFhc_nowgt);
    THStack *hsNumuRhc_nowgt = makeFitStack(Form("hsNumuRhc_nowgt_%d", i), "#nu_{#mu} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNumuRhc_nowgt, hdataNumuRhc_nowgt, hpostNumuRhc_nowgt);
    THStack *hsNueRhc_nowgt  = makeFitStack(Form("hsNueRhc_nowgt_%d", i), "#nu_{e} (RHC); E_{#nu, reco} / GeV; Events / GeV", hpreNueRhc_nowgt, hdataNueRhc_nowgt, hpostNueRhc_nowgt);
    hsNumuFhc_nowgt->Write();
    hsNueFhc_nowgt->Write();
    hsNumuRhc_nowgt->Write();
    hsNueRhc_nowgt->Write();

    TTree *fd_only_tree_norwt = new TTree(Form("fd_only_norwt_%d",i), Form("fd_only_norwt_%d",i));
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
double 
    bias = (dCP - fPostFitValues.at(3)*TMath::Pi())*(180./TMath::Pi());
    if (bias<-180) bias+=360;
    else if (bias>360) bias-=360;
    dCPBias_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), bias);

    fTrueValues.resize(fParamNames.size(), 0.);
    fTrueValues.at(0) = this_calc->GetDmsq32()*1000.;
    fTrueValues.at(1) = this_calc->GetTh13();
    fTrueValues.at(2) = TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23());
    fTrueValues.at(3) = this_calc->GetdCP();
    // fTrueValues.at(4) = this_calc->GetRho();

    fd_only_tree_norwt->Fill();
    fd_only_tree_norwt->Write();

    // rho_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), this_calc->GetRho()-fPostFitValues.at(4));
    theta13_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), this_calc->GetTh13()-fPostFitValues.at(1));
    theta23_nowgt->SetPoint(i, dCP*(180./TMath::Pi()), TMath::Sin(this_calc->GetTh23())*TMath::Sin(this_calc->GetTh23())-fPostFitValues.at(2));
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

  fout->Close();
} // fdFitWithRwgt
		  
