// transverseMom.C
// Transverse momentum imbalance in CCQE interactions

// CAFAna includes
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
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
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/common_fit_definitions.h"
// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"

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

// POT for 3.5 years
const double pot_fd = 3.5 * POT120 * 40/1.13;
const double pot_nd = /*3.5*/ 0.176 * POT120;

const Var kNuMomY  = SIMPLEVAR(dune.NuMomY);
const Var kLepMomX = SIMPLEVAR(dune.LepMomX);
const Var kLepMomY = SIMPLEVAR(dune.LepMomY);
const Var kLepMomZ = SIMPLEVAR(dune.LepMomZ);
const Var kLepNuAngle = SIMPLEVAR(dune.LepNuAngle);
const Var kProMomX = SIMPLEVAR(dune.gastpc_ProMomX);
const Var kProMomY = SIMPLEVAR(dune.gastpc_ProMomY);
const Var kProMomZ = SIMPLEVAR(dune.gastpc_ProMomZ);
const Var kRecoProMomX = SIMPLEVAR(dune.gastpc_RecoProMomX);
const Var kRecoProMomY = SIMPLEVAR(dune.gastpc_RecoProMomY);
const Var kRecoProMomZ = SIMPLEVAR(dune.gastpc_RecoProMomZ);
const Var kNRecoFS = SIMPLEVAR(dune.gastpc_nRecoFS);
const Var kNFS = SIMPLEVAR(dune.nN)+SIMPLEVAR(dune.nP)+SIMPLEVAR(dune.nipip)+SIMPLEVAR(dune.nipim)+SIMPLEVAR(dune.nipi0)+SIMPLEVAR(dune.nikp)+SIMPLEVAR(dune.nikm)+SIMPLEVAR(dune.nik0)+SIMPLEVAR(dune.nNucleus);

const Var kTotalFSMomX({"dune.gastpc_ProMomX", "dune.LepMomX"},
		       [](const caf::StandardRecord* sr) {
			 double mom = 0;
			 mom =  sr->dune.gastpc_ProMomX + sr->dune.LepMomX;
			 return mom;
		       });
const Var kTotalFSMomY({"dune.gastpc_ProMomY", "dune.LepMomY"},
		       [](const caf::StandardRecord* sr) {
			 double mom = 0;
			 mom =  sr->dune.gastpc_ProMomY + sr->dune.LepMomY;
			 return mom;
		       });
const Var kTotalFSMomZ({"dune.gastpc_ProMomZ", "dune.LepMomZ"},
		       [](const caf::StandardRecord* sr) {
			 double mom = 0;
			 mom =  sr->dune.gastpc_ProMomZ + sr->dune.LepMomZ;
			 return mom;
		       });

const Binning simpleMomBins  = Binning::Simple(30, 0, 3);
const Binning simpleMomBins2 = Binning::Simple(60, -1.5, 1.5);
const Binning simpleProMomBins  = Binning::Simple(30, 0, 1.5);
const Binning simpleProMomBins2 = Binning::Simple(60, -1.5, 1.5);
const HistAxis axLepMomX("P_{x, #mu} / GeV", simpleMomBins2, kLepMomX);
const HistAxis axLepMomY("P_{y, #mu} / GeV", simpleMomBins2, kLepMomY);
const HistAxis axLepMomZ("P_{z, #mu} / GeV", simpleMomBins, kLepMomZ);
const HistAxis axLepNuAngle("#theta_{l,#nu}", simpleMomBins, kLepNuAngle);
const HistAxis axProMomX("P_{x, p} / GeV", simpleProMomBins2, kProMomX);
const HistAxis axProMomY("P_{y, p} / GeV", simpleProMomBins2, kProMomY);
const HistAxis axProMomZ("P_{z, p} / GeV", simpleProMomBins, kProMomZ);
const HistAxis axRecoProMomX("P_{x, p} / GeV", simpleProMomBins2, kRecoProMomX);
const HistAxis axRecoProMomY("P_{y, p} / GeV", simpleProMomBins2, kRecoProMomY);
const HistAxis axRecoProMomZ("P_{z, p} / GeV", simpleProMomBins, kRecoProMomZ);
const HistAxis axNuMomY("P_{y, #nu} / GeV", simpleMomBins2, kNuMomY);
const HistAxis axTotalFSMomX("P_{x, FS} / GeV", simpleMomBins2, kTotalFSMomX);
const HistAxis axTotalFSMomY("P_{y, FS} / GeV", simpleMomBins2, kTotalFSMomY);
const HistAxis axTotalFSMomZ("P_{z, FS} / GeV", simpleMomBins, kTotalFSMomZ);
const HistAxis axFSvsNuMomY("P_{y, FS} / GeV", simpleMomBins2, kTotalFSMomY, 
			    "P_{y, #nu} / GeV", simpleMomBins2, kNuMomY);

void transverseMom(const char* outFile, const char* garDir="/dune/data/users/sbjones/gasTpcCAF/v6/")
{
  gROOT->SetBatch(kTRUE);
  rootlogon();

  TFile *fout = new TFile(outFile, "recreate");

  Loaders loadersGArFHC;
  Loaders loadersGArRHC;
  SpectrumLoader loaderGArFHC(Form("%s/CAF_FHC.root", garDir), kBeam);
  SpectrumLoader loaderGArRHC(Form("%s/CAF_RHC.root", garDir), kBeam);
  loadersGArFHC.AddLoader(&loaderGArFHC, caf::kNEARDET, Loaders::kMC);
  loadersGArRHC.AddLoader(&loaderGArRHC, caf::kNEARDET, Loaders::kMC);

  NoOscPredictionGenerator genNuMomYFhc(axNuMomY, kIsTrueGasFV && kProMomX!=0);
  NoOscPredictionGenerator genLepMomXFhc(axLepMomX, kIsTrueGasFV && kProMomX!=0);
  NoOscPredictionGenerator genLepMomYFhc(axLepMomY, kIsTrueGasFV && kProMomY!=0);
  NoOscPredictionGenerator genLepMomZFhc(axLepMomZ, kIsTrueGasFV && kProMomZ!=0);
  NoOscPredictionGenerator genProMomXFhc(axProMomX, kIsTrueGasFV && kProMomX!=0);
  NoOscPredictionGenerator genProMomYFhc(axProMomY, kIsTrueGasFV && kProMomY!=0);
  NoOscPredictionGenerator genProMomZFhc(axProMomZ, kIsTrueGasFV && kProMomZ!=0);
  NoOscPredictionGenerator genTotalFSMomXFhc(axTotalFSMomX, kIsTrueGasFV && kProMomX!=0);
  NoOscPredictionGenerator genTotalFSMomYFhc(axTotalFSMomY, kIsTrueGasFV && kProMomY!=0);
  NoOscPredictionGenerator genTotalFSMomZFhc(axTotalFSMomZ, kIsTrueGasFV && kProMomZ!=0);
  NoOscPredictionGenerator genLepNuAngleFhc(axLepNuAngle, kIsTrueGasFV && kProMomZ!=0);
  NoOscPredictionGenerator genFSvsNuMomY(axFSvsNuMomY, kIsTrueGasFV && kProMomZ!=0);


  PredictionInterp predNuMomYFhc({}, 0, genNuMomYFhc, loadersGArFHC);
  PredictionInterp predLepMomXFhc({}, 0, genLepMomXFhc, loadersGArFHC);
  PredictionInterp predLepMomYFhc({}, 0, genLepMomYFhc, loadersGArFHC);
  PredictionInterp predLepMomZFhc({}, 0, genLepMomZFhc, loadersGArFHC);
  PredictionInterp predProMomXFhc({}, 0, genProMomXFhc, loadersGArFHC);
  PredictionInterp predProMomYFhc({}, 0, genProMomYFhc, loadersGArFHC);
  PredictionInterp predProMomZFhc({}, 0, genProMomZFhc, loadersGArFHC);
  PredictionInterp predTotalFSMomXFhc({}, 0, genTotalFSMomXFhc, loadersGArFHC);
  PredictionInterp predTotalFSMomYFhc({}, 0, genTotalFSMomYFhc, loadersGArFHC);
  PredictionInterp predTotalFSMomZFhc({}, 0, genTotalFSMomZFhc, loadersGArFHC);
  PredictionInterp predLepNuAngleFhc({}, 0, genLepNuAngleFhc, loadersGArFHC);
  PredictionInterp predFSvsNuMomY({}, 0, genFSvsNuMomY, loadersGArFHC);

  loadersGArFHC.Go();
  // loadersGArRHC.Go();

  TH1 *hNuMomYFhc = predNuMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hLepMomXFhc = predLepMomXFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hLepMomYFhc = predLepMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hLepMomZFhc = predLepMomZFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hProMomXFhc = predProMomXFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hProMomYFhc = predProMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hProMomZFhc = predProMomZFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hTotalFSMomXFhc = predTotalFSMomXFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hTotalFSMomYFhc = predTotalFSMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hTotalFSMomZFhc = predTotalFSMomZFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hLepNuAngleFhc = predLepNuAngleFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH2 *hFSvsNuMomY = predFSvsNuMomY.Predict(0).FakeData(pot_nd).ToTH2(pot_nd);
  setHistAttr(hNuMomYFhc);
  setHistAttr(hLepMomXFhc);
  setHistAttr(hLepMomYFhc);
  setHistAttr(hLepMomZFhc);
  setHistAttr(hProMomXFhc);
  setHistAttr(hProMomYFhc);
  setHistAttr(hProMomZFhc);
  setHistAttr(hTotalFSMomXFhc);
  setHistAttr(hTotalFSMomYFhc);
  setHistAttr(hTotalFSMomZFhc);
  setHistAttr(hLepNuAngleFhc);
  setHistAttr(hFSvsNuMomY);

  fout->cd();
  hNuMomYFhc->Write("hNuMomYFhc");
  hLepMomXFhc->Write("hLepMomXFhc");
  hLepMomYFhc->Write("hLepMomYFhc");
  hLepMomZFhc->Write("hLepMomZFhc");
  hProMomXFhc->Write("hProMomXFhc");
  hProMomYFhc->Write("hProMomYFhc");
  hProMomZFhc->Write("hProMomZFhc");
  hTotalFSMomXFhc->Write("hTotalFSMomXFhc");
  hTotalFSMomYFhc->Write("hTotalFSMomYFhc");
  hTotalFSMomZFhc->Write("hTotalFSMomZFhc");
  hLepNuAngleFhc->Write("hLepNuAngleFhc");
  hFSvsNuMomY->Write("hFSvsNuMomY");

  fout->Close();
  delete fout;
} // transverseMom
