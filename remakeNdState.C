// remakeNdState.C
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
const Var Enu_reco = SIMPLEVAR(dune.Ev_reco);

std::vector<double> binEEdges = {0., 0.5, 1., 1.25, 1.5, 1.75,
                                 2., 2.25, 2.5, 2.75, 
				 3., 3.25, 3.5, 3.75, 
				 4., 5., 6., 10.};
std::vector<double> binYEdges = {0, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0};

const Binning binsEreco  = Binning::Custom(binEEdges);
const Binning binsYreco  = Binning::Custom(binYEdges);
const HistAxis axND("E_{#nu, reco} (GeV)", binsEreco, Enu_reco,
		    "y_{reco}", binsYreco, kRecoYND);

void remakeNdState(bool isFhc, bool isWeighted, bool isLar,
		   const char* stateFileDir="/dune/data/users/sbjones/stateFiles/withNuWroRWT/highQ2/",
		   const char* cafs="/dune/data/users/sbjones/CAFs/nuwroRW/withEAvail/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);
  
  std::vector<const ISyst*> systlist = GetListOfSysts(true, true, 
						      true, false, true,
						      false, false,
						      false, 20, 
						      true); 
  
  std::vector<const ISyst*> fakedata = GetXSecSysts({"NuWroReweightFakeData"});
  systlist.insert(systlist.end(), fakedata.end()-1, fakedata.end());
  std::vector<const ISyst*> nuwrorwtq0q3 = GetXSecSysts({"NuWroQ0Q3"});
  std::cout<<"Systs to make state files"<<std::endl;
  for (unsigned int i=0; i<systlist.size(); i++) {
    std::cout<<systlist.at(i)->ShortName()<<std::endl;
  }

  SystShifts ndWeightq0q3(nuwrorwtq0q3.at(nuwrorwtq0q3.size()-3), 1);
  SystShifts ndWeightq0q3CC(nuwrorwtq0q3.at(nuwrorwtq0q3.size()-4), 1);
  std::cout<<"ndweight Q0Q3 CC syst is "<<nuwrorwtq0q3.at(nuwrorwtq0q3.size()-4)->ShortName()<<std::endl;
  std::cout<<"ndweight Q0Q3 syst is "<<nuwrorwtq0q3.at(nuwrorwtq0q3.size()-3)->ShortName()<<std::endl;

  Loaders loadersFHC;
  Loaders loadersRHC;
  SpectrumLoader loaderFHC(Form("%s/ND_FHC_CAF.root", cafs), kBeam);
  SpectrumLoader loaderRHC(Form("%s/ND_RHC_CAF.root", cafs), kBeam);
  loadersFHC.AddLoader(&loaderFHC, caf::kNEARDET, Loaders::kMC);
  loadersRHC.AddLoader(&loaderRHC, caf::kNEARDET, Loaders::kMC);

  if (isFhc) {
    if (isWeighted) {
      if (isLar) {
	NoOscPredictionGenerator genFhc(axND, kPassND_FHC_NUMU && kIsTrueFV, kCVXSecWeights); 
	PredictionInterp predFhc(systlist, 0, genFhc, loadersFHC, ndWeightq0q3CC);
	loadersFHC.Go();
	TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC_wgt_cc.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted ND numu FHC to "<<Form("%s/state_ND_FHC_wgt_cc.root", stateFileDir)<<std::endl;
	predFhc.SaveTo(fNDfhc->mkdir("nd_interp_numu_fhc_wgt"));
	fNDfhc->Close();
	delete fNDfhc;
	std::cout<<"Made & saved PredictionInterps"<<std::endl;
      }
      else {
	NoOscPredictionGenerator genFhc(axND, kPassND_FHC_NUMU && kIsTrueFV, kCVXSecWeights); 
	PredictionInterp predFhc(systlist, 0, genFhc, loadersFHC, ndWeightq0q3);
	loadersFHC.Go();
	TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC_wgt_q0q3.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted ND numu FHC to "<<Form("%s/state_ND_FHC_wgt_q0q3.root", stateFileDir)<<std::endl;
	predFhc.SaveTo(fNDfhc->mkdir("nd_interp_numu_fhc_wgt"));
	fNDfhc->Close();
	delete fNDfhc;
	std::cout<<"Made & saved PredictionInterps"<<std::endl;
      }
    }
    else {
      NoOscPredictionGenerator genFhc(axND, kPassND_FHC_NUMU && kIsTrueFV, kCVXSecWeights); 
      PredictionInterp predFhc(systlist, 0, genFhc, loadersFHC);
      loadersFHC.Go();
      TFile *fNDfhc = new TFile(Form("%s/state_ND_FHC.root", stateFileDir), "recreate"); 
      std::cout<<"Saving unweighted ND numu FHC to "<<Form("%s/state_ND_FHC.root", stateFileDir)<<std::endl;
      predFhc.SaveTo(fNDfhc->mkdir("nd_interp_numu_fhc"));
      fNDfhc->Close();
      delete fNDfhc;
      std::cout<<"Made & saved PredictionInterps"<<std::endl;
    }
  }
  else {
    if (isWeighted) {
      if (isLar) {
	NoOscPredictionGenerator genRhc(axND, kPassND_RHC_NUMU && kIsTrueFV, kCVXSecWeights); 
	PredictionInterp predRhc(systlist, 0, genRhc, loadersRHC, ndWeightq0q3CC); 
	loadersRHC.Go();
	TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC_wgt_cc.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted ND numu RHC to "<<Form("%s/state_ND_RHC_wgt_cc.root", stateFileDir)<<std::endl;
	predRhc.SaveTo(fNDrhc->mkdir("nd_interp_numu_rhc_wgt"));
	fNDrhc->Close();
	delete fNDrhc;
	std::cout<<"Made & saved PredictionInterps"<<std::endl;
      }
      else {
	NoOscPredictionGenerator genRhc(axND, kPassND_RHC_NUMU && kIsTrueFV, kCVXSecWeights); 
	PredictionInterp predRhc(systlist, 0, genRhc, loadersRHC, ndWeightq0q3); 
	loadersRHC.Go();
	TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC_wgt_q0q3.root", stateFileDir), "recreate"); 
	std::cout<<"Saving weighted ND numu RHC to "<<Form("%s/state_ND_RHC_wgt_q0q3.root", stateFileDir)<<std::endl;
	predRhc.SaveTo(fNDrhc->mkdir("nd_interp_numu_rhc_wgt"));
	fNDrhc->Close();
	delete fNDrhc;
	std::cout<<"Made & saved PredictionInterps"<<std::endl;
      }
    }
    else {
      NoOscPredictionGenerator genRhc(axND, kPassND_RHC_NUMU && kIsTrueFV, kCVXSecWeights); 
      PredictionInterp predRhc(systlist, 0, genRhc, loadersRHC); 
      loadersRHC.Go();
      TFile *fNDrhc = new TFile(Form("%s/state_ND_RHC.root", stateFileDir), "recreate"); 
      std::cout<<"Saving unweighted ND numu RHC to "<<Form("%s/state_ND_RHC.root", stateFileDir)<<std::endl;
      predRhc.SaveTo(fNDrhc->mkdir("nd_interp_numu_rhc"));
      fNDrhc->Close();
      delete fNDrhc;
      std::cout<<"Made & saved PredictionInterps"<<std::endl;
    }
  }

}
