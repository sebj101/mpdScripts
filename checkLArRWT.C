// checkLArRWT.C
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

  hmc_wgt->Write(Form("%s_mc_wgt", name));
  hmc_lar_wgt->Write(Form("%s_mc_lar_wgt", name));
  hdata->Write(Form("%s_data", name));
  hmc_no_wgt->Write(Form("%s_mc_no_wgt", name));

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
const HistAxis axisnumu("E_{#nu, reco} (GeV)", binsEreco, Enu_reco_numu);
const HistAxis axisnue("E_{#nu, reco} (GeV)", binsEreco, Enu_reco_nue);
// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;

//const char* cafs="/dune/data/users/sbjones/CAFs/nuwroRW/";
const char* cafs="/dune/data/users/marshalc/CAFs/nuwroRW/";

void checkLArRWT(const char* outfile)
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);

  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true);

  std::vector<const ISyst*> systs = {};
  std::vector<const ISyst*> ndwgtQ2  = GetXSecSysts({"NuWroQ2"});
  std::vector<const ISyst*> ndwgtW   = GetXSecSysts({"NuWroW"});
  std::vector<const ISyst*> ndwgtlar = GetXSecSysts({"NuWroRWTLar"});
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});

  std::cout<<fakedata.at(fakedata.size()-1)->ShortName()<<" is fake data"<<std::endl;
  std::cout<<ndwgtQ2.at(ndwgtQ2.size()-4)->ShortName()<<" is Q2 rwt data"<<std::endl;
  std::cout<<ndwgtW.at(ndwgtW.size()-3)->ShortName()<<" is W rwt data"<<std::endl;
  std::cout<<ndwgtlar.at(ndwgtlar.size()-2)->ShortName()<<" is lar rwt data"<<std::endl;

  systs.push_back(fakedata.at(fakedata.size()-1));
  systs.push_back(ndwgtQ2.at(ndwgtQ2.size()-4));
  systs.push_back(ndwgtW.at(ndwgtW.size()-3));
  systs.push_back(ndwgtlar.at(ndwgtlar.size()-2));

  for (unsigned int i=0; i<systs.size(); i++) {
    std::cout<<systs.at(i)->ShortName()<<std::endl;
  }

  SystShifts fakeshift(fakedata.at(fakedata.size()-1), 1);
  SystShifts wgtshift(ndwgtQ2.at(ndwgtQ2.size()-3), 1);
  SystShifts wgtshiftW(ndwgtW.at(ndwgtW.size()-3), 1);
  SystShifts wgtlarshift(ndwgtlar.at(ndwgtlar.size()-2), 1);

  Loaders loadersFHC, loadersRHC;
  loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
  loadersFHC.SetLoaderPath(Form("%s/FD_FHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);

  loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
  loadersRHC.SetLoaderPath(Form("%s/FD_RHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);

  NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
  NoExtrapGenerator gennuereco (axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 

  PredictionInterp predNumuFHCreco(systs, this_calc, gennumureco, loadersFHC); 
  PredictionInterp predNumuFHCreco_wgt(systs, this_calc, gennumureco, loadersFHC, wgtshift); 
  PredictionInterp predNumuFHCreco_wgt_w(systs, this_calc, gennumureco, loadersFHC, wgtshiftW); 
  PredictionInterp predNumuFHCreco_lar_wgt(systs, this_calc, gennumureco, loadersFHC, wgtlarshift); 

  PredictionInterp predNueFHCreco(systs, this_calc, gennuereco, loadersFHC); 
  PredictionInterp predNueFHCreco_wgt(systs, this_calc, gennuereco, loadersFHC, wgtshift); 
  PredictionInterp predNueFHCreco_wgt_w(systs, this_calc, gennuereco, loadersFHC, wgtshiftW); 
  PredictionInterp predNueFHCreco_lar_wgt(systs, this_calc, gennuereco, loadersFHC, wgtlarshift); 

  PredictionInterp predNumuRHCreco(systs, this_calc, gennumureco, loadersRHC); 
  PredictionInterp predNumuRHCreco_wgt(systs, this_calc, gennumureco, loadersRHC, wgtshift); 
  PredictionInterp predNumuRHCreco_wgt_w(systs, this_calc, gennumureco, loadersRHC, wgtshiftW); 
  PredictionInterp predNumuRHCreco_lar_wgt(systs, this_calc, gennumureco, loadersRHC, wgtlarshift); 

  PredictionInterp predNueRHCreco(systs, this_calc, gennuereco, loadersRHC); 
  PredictionInterp predNueRHCreco_wgt(systs, this_calc, gennuereco, loadersRHC, wgtshift); 
  PredictionInterp predNueRHCreco_wgt_w(systs, this_calc, gennuereco, loadersRHC, wgtshiftW); 
  PredictionInterp predNueRHCreco_lar_wgt(systs, this_calc, gennuereco, loadersRHC, wgtlarshift); 

  loadersFHC.Go();
  loadersRHC.Go();

  TFile *fout = new TFile(outfile, "recreate");
  fout->cd();
  THStack *hsnumufhc = makeBigStacks("hsnumufhc", "#nu_{#mu} FHC; E_{#nu, reco}; Events / GeV", 
   				     &predNumuFHCreco_wgt, &predNumuFHCreco_lar_wgt, &predNumuFHCreco, 
   				     fakeshift, this_calc, pot_fd);
  THStack *hsnumurhc = makeBigStacks("hsnumurhc", "#nu_{#mu} RHC; E_{#nu, reco}; Events / GeV", 
   				     &predNumuRHCreco_wgt, &predNumuRHCreco_lar_wgt, &predNumuRHCreco, 
   				     fakeshift, this_calc, pot_fd);
  THStack *hsnuefhc = makeBigStacks("hsnuefhc", "#nu_{e} FHC; E_{#nu, reco}; Events / GeV", 
   				     &predNueFHCreco_wgt, &predNueFHCreco_lar_wgt, &predNueFHCreco, 
   				     fakeshift, this_calc, pot_fd);
  THStack *hsnuerhc = makeBigStacks("hsnuerhc", "#nu_{e} RHC; E_{#nu, reco}; Events / GeV", 
   				     &predNueRHCreco_wgt, &predNueRHCreco_lar_wgt, &predNueRHCreco, 
   				     fakeshift, this_calc, pot_fd);
  hsnumufhc->Write();

  TH1 *hnumufhc_w = predNumuFHCreco_wgt_w.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kGreen+2, kDashed, kPOT, kBinDensity);
  TH1 *hnumurhc_w = predNumuRHCreco_wgt_w.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kGreen+2, kDashed, kPOT, kBinDensity);
  TH1 *hnuefhc_w  = predNueFHCreco_wgt_w.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kGreen+2, kDashed, kPOT, kBinDensity);
  TH1 *hnuerhc_w  = predNueRHCreco_wgt_w.Predict(this_calc).FakeData(pot_fd).ToTH1(pot_fd, kGreen+2, kDashed, kPOT, kBinDensity);
  hnumufhc_w->Write("hnumufhc_w");
  hnumurhc_w->Write("hnumurhc_w");
  hnuefhc_w->Write("hnuefhc_w");
  hnuerhc_w->Write("hnuerhc_w");

  fout->Close();

} // checkLArRWT
