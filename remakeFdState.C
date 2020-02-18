// remakeFdState.C
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

const Var kFHC = SIMPLEVAR(dune.isFHC);
const Var Enu_reco_numu = SIMPLEVAR(dune.Ev_reco_numu);
const Var Enu_reco_nue  = SIMPLEVAR(dune.Ev_reco_nue);

std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
const Binning binsEreco  = Binning::Custom(binEEdges);

const HistAxis axisnumu("Reco #nu energy (GeV)", binsEreco, Enu_reco_numu);
const HistAxis axisnue("Reco #nu energy (GeV)", binsEreco, Enu_reco_nue);

// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;
// Oscillation variables to fit
std::vector<const IFitVar*> fitVars = {&kFitDmSq32NHScaled, &kFitTheta13, &kFitSinSqTheta23, &kFitDeltaInPiUnits, &kFitRho};

void remakeFdState(bool isFhc, bool isWeighted, bool isLar,
		   const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2/",
		   const char* cafs="/dune/data/users/sbjones/CAFs/nuwroRW/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  
  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true); 
  
  //  std::vector<const ISyst*> systlist = {};
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});//, NuWroRWT"});
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  std::vector<const ISyst*> nuwrorwt    = GetXSecSysts({"NuWroRWT"});
  std::vector<const ISyst*> nuwrorwtlar = GetXSecSysts({"NuWroRWTLar"});
  systlist.insert(systlist.end(), nuwrorwtlar.end()-3, nuwrorwtlar.end()-1);
  std::cout<<"Systs to make state files"<<std::endl;
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist.at(i)->ShortName()<<std::endl;
  }

  SystShifts ndWeight(nuwrorwt.at(nuwrorwt.size()-3), 1);
  SystShifts ndWeightlar(nuwrorwtlar.at(nuwrorwtlar.size()-2), 1);
  std::cout<<"ndweight syst is "<<nuwrorwt.at(nuwrorwt.size()-3)->ShortName()<<std::endl;
  std::cout<<"ndweightlar syst is "<<nuwrorwtlar.at(nuwrorwtlar.size()-2)->ShortName()<<std::endl;

  if (isFhc) {
    Loaders loadersFHC;
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersFHC.SetLoaderPath(Form("%s/FD_FHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    if (isWeighted) {
      if (isLar) {
 	PredictionInterp predNumuFHCreco_wgt(systlist, this_calc, gennumureco, loadersFHC, ndWeightlar); 
	PredictionInterp predNueFHCreco_wgt(systlist, this_calc, gennuereco, loadersFHC, ndWeightlar); 
	loadersFHC.Go();

	TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC_wgt_lar.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted FD numu FHC to "<<Form("%s/state_FD_FHC_wgt_lar.root", stateFileDir)<<std::endl;
	predNumuFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc_wgt"));
	std::cout<<"Saving weighted FD nue FHC to "<<Form("%s/state_FD_FHC_wgt_lar.root", stateFileDir)<<std::endl;
	predNueFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc_wgt"));
	fFDfhc->Close();
	delete fFDfhc;
	std::cout<<"Made & saved all PredictionInterps"<<std::endl;
      }
      else {
	PredictionInterp predNumuFHCreco_wgt(systlist, this_calc, gennumureco, loadersFHC, ndWeight); 
	PredictionInterp predNueFHCreco_wgt(systlist, this_calc, gennuereco, loadersFHC, ndWeight); 
	loadersFHC.Go();

	TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC_wgt.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted FD numu FHC to "<<Form("%s/state_FD_FHC_wgt.root", stateFileDir)<<std::endl;
	predNumuFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc_wgt"));
	std::cout<<"Saving weighted FD nue FHC to "<<Form("%s/state_FD_FHC_wgt.root", stateFileDir)<<std::endl;
	predNueFHCreco_wgt.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc_wgt"));
	fFDfhc->Close();
	delete fFDfhc;
	std::cout<<"Made & saved all PredictionInterps"<<std::endl;
      }      
    }
    else {
      PredictionInterp predNumuFHCreco(systlist, this_calc, gennumureco, loadersFHC); 
      PredictionInterp predNueFHCreco(systlist, this_calc, gennuereco, loadersFHC); 
      loadersFHC.Go();

      TFile *fFDfhc = new TFile(Form("%s/state_FD_FHC.root", stateFileDir), "recreate"); 
      std::cout<<"Saving FD numu FHC to "<<Form("%s/state_FD_FHC.root", stateFileDir)<<std::endl;
      predNumuFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_numu_fhc"));
      std::cout<<"Saving FD nue FHC to "<<Form("%s/state_FD_FHC.root", stateFileDir)<<std::endl;
      predNueFHCreco.SaveTo(fFDfhc->mkdir("fd_interp_nue_fhc"));
      fFDfhc->Close();
      delete fFDfhc;
      std::cout<<"Made & saved all PredictionInterps"<<std::endl;
    }
  }
  else {
    Loaders loadersRHC;
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nonswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_nueswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    loadersRHC.SetLoaderPath(Form("%s/FD_RHC_tauswap.root", cafs) , caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    NoExtrapGenerator gennumureco(axisnumu, kPassFD_CVN_NUMU && kIsTrueFV, kCVXSecWeights); 
    NoExtrapGenerator gennuereco(axisnue, kPassFD_CVN_NUE && kIsTrueFV, kCVXSecWeights); 
    if (isWeighted) {
      if (isLar) {
	PredictionInterp predNumuRHCreco_wgt(systlist, this_calc, gennumureco, loadersRHC, ndWeightlar); 
	PredictionInterp predNueRHCreco_wgt(systlist, this_calc, gennuereco, loadersRHC, ndWeightlar);
	loadersRHC.Go();

	TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC_wgt_lar.root", stateFileDir), "recreate");
	std::cout<<"Saving weighted FD numu RHC to "<<Form("%s/state_FD_RHC_wgt_lar.root", stateFileDir)<<std::endl;
	predNumuRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc_wgt"));
	std::cout<<"Saving weighted FD nue RHC to "<<Form("%s/state_FD_RHC_wgt_lar.root", stateFileDir)<<std::endl;
	predNueRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc_wgt"));
	fFDrhc->Close();
	delete fFDrhc;
	std::cout<<"Made & saved all PredictionInterps"<<std::endl;
      }
      else {
	PredictionInterp predNumuRHCreco_wgt(systlist, this_calc, gennumureco, loadersRHC, ndWeight); 
	PredictionInterp predNueRHCreco_wgt(systlist, this_calc, gennuereco, loadersRHC, ndWeight);
	loadersRHC.Go();

	TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC_wgt.root", stateFileDir), "recreate");
	std::cout<<"Saving weighted FD numu RHC to "<<Form("%s/state_FD_RHC_wgt.root", stateFileDir)<<std::endl;
	predNumuRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc_wgt"));
	std::cout<<"Saving weighted FD nue RHC to "<<Form("%s/state_FD_RHC_wgt.root", stateFileDir)<<std::endl;
	predNueRHCreco_wgt.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc_wgt"));
	fFDrhc->Close();
	delete fFDrhc;
	std::cout<<"Made & saved all PredictionInterps"<<std::endl;
      }
    }
    else {
      PredictionInterp predNumuRHCreco(systlist, this_calc, gennumureco, loadersRHC); 
      PredictionInterp predNueRHCreco(systlist, this_calc, gennuereco, loadersRHC);
      loadersRHC.Go();

      TFile *fFDrhc = new TFile(Form("%s/state_FD_RHC.root", stateFileDir), "recreate");
      std::cout<<"Saving FD numu RHC to "<<Form("%s/state_FD_RHC.root", stateFileDir)<<std::endl;
      predNumuRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_numu_rhc"));
      std::cout<<"Saving FD nue RHC to "<<Form("%s/state_FD_RHC.root", stateFileDir)<<std::endl;
      predNueRHCreco.SaveTo(fFDrhc->mkdir("fd_interp_nue_rhc"));
      fFDrhc->Close();
      delete fFDrhc;
      std::cout<<"Made & saved all PredictionInterps"<<std::endl;
    }
  }
} // remakeFdState
