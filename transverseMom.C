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

// POT for n years
const double years = 1.;
const double pot_fd = years * POT120 * 40/1.13;
const double pot_nd = years * POT120;

const Var kRecoP   = SIMPLEVAR(dune.gastpc_pro_mult);
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
const Var kRecoLepMomX = SIMPLEVAR(dune.gastpc_RecoLepMomX);
const Var kRecoLepMomY = SIMPLEVAR(dune.gastpc_RecoLepMomY);
const Var kRecoLepMomZ = SIMPLEVAR(dune.gastpc_RecoLepMomZ);
const Var kNRecoFS = SIMPLEVAR(dune.gastpc_nRecoFS);
const Var kNFS = SIMPLEVAR(dune.nN)+SIMPLEVAR(dune.nP)+SIMPLEVAR(dune.nipip)+SIMPLEVAR(dune.nipim)+SIMPLEVAR(dune.nipi0);

const Var kRecoTotalFSMomX = SIMPLEVAR(dune.gastpc_RecoLepMomX) + SIMPLEVAR(dune.gastpc_RecoProMomX);
const Var kRecoTotalFSMomY = SIMPLEVAR(dune.gastpc_RecoLepMomY) + SIMPLEVAR(dune.gastpc_RecoProMomY);
const Var kRecoTotalFSMomZ = SIMPLEVAR(dune.gastpc_RecoLepMomZ) + SIMPLEVAR(dune.gastpc_RecoProMomZ);

// Reco transverse mom
const Var kRecoTransMomAng({"dune.gastpc_RecoLepMomX", "dune.gastpc_RecoLepMomY", "dune.gastpc_RecoProMomX", "dune.gastpc_RecoProMomY"},
			   [](const caf::StandardRecord* sr) {
			     double momX = sr->dune.gastpc_RecoLepMomX + sr->dune.gastpc_RecoProMomX;
			     double momY = sr->dune.gastpc_RecoLepMomY + sr->dune.gastpc_RecoProMomY;
			     double magT = sqrt(momX*momX + momY*momY);
			     double magMu = sqrt(sr->dune.gastpc_RecoLepMomX*sr->dune.gastpc_RecoLepMomX + sr->dune.gastpc_RecoLepMomY*sr->dune.gastpc_RecoLepMomY);
			     double alpha = acos(-1.*(momX*sr->dune.gastpc_RecoLepMomX +
						      momY*sr->dune.gastpc_RecoLepMomY)/(magT*magMu));
			     return alpha;
			   });
const Var kRecoTransMomMag({"dune.gastpc_RecoLepMomX", "dune.gastpc_RecoLepMomY", "dune.gastpc_RecoProMomX", "dune.gastpc_RecoProMomY"},
			   [](const caf::StandardRecord* sr) {
			     double momX = sr->dune.gastpc_RecoLepMomX + sr->dune.gastpc_RecoProMomX;
			     double momY = sr->dune.gastpc_RecoLepMomY + sr->dune.gastpc_RecoProMomY;
			     double mag = sqrt(momX*momX + momY*momY);
			     return mag;
			   });
// True transverse mom
const Var kTransMomAng({"dune.LepMomX", "dune.LepMomY", "dune.gastpc_ProMomX", "dune.gastpc_ProMomY"},
			   [](const caf::StandardRecord* sr) {
			     double momX = sr->dune.LepMomX + sr->dune.gastpc_ProMomX;
			     double momY = sr->dune.LepMomY + sr->dune.gastpc_ProMomY;
			     double magT = sqrt(momX*momX + momY*momY);
			     double magMu = sqrt(sr->dune.LepMomX*sr->dune.LepMomX + sr->dune.LepMomY*sr->dune.LepMomY);
			     double alpha = acos(-1.*(momX*sr->dune.LepMomX +
						      momY*sr->dune.LepMomY)/(magT*magMu));
			     return alpha;
			   });
const Var kTransMomMag({"dune.LepMomX", "dune.LepMomY", "dune.gastpc_ProMomX", "dune.gastpc_ProMomY"},
			   [](const caf::StandardRecord* sr) {
			     double momX = sr->dune.LepMomX + sr->dune.gastpc_ProMomX;
			     double momY = sr->dune.LepMomY + sr->dune.gastpc_ProMomY;
			     double mag = sqrt(momX*momX + momY*momY);
			     return mag;
			   });


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
const Binning binsAng = Binning::Simple(40, 0., TMath::Pi());
const Binning binsMag = Binning::Simple(40, 0., 1.);
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
const HistAxis axRecoLepMomX("P_{x, #mu} / GeV", simpleProMomBins2, kRecoLepMomX);
const HistAxis axRecoLepMomY("P_{y, #mu} / GeV", simpleProMomBins2, kRecoLepMomY);
const HistAxis axRecoLepMomZ("P_{z, #mu} / GeV", simpleProMomBins, kRecoLepMomZ);
const HistAxis axNuMomY("P_{y, #nu} / GeV", simpleMomBins2, kNuMomY);
const HistAxis axTotalFSMomX("P_{x, FS} / GeV", simpleMomBins2, kTotalFSMomX);
const HistAxis axTotalFSMomY("P_{y, FS} / GeV", simpleMomBins2, kTotalFSMomY);
const HistAxis axTotalFSMomZ("P_{z, FS} / GeV", simpleMomBins, kTotalFSMomZ);
const HistAxis axFSvsNuMomY("P_{y, FS} / GeV", simpleMomBins2, kTotalFSMomY, 
			    "P_{y, #nu} / GeV", simpleMomBins2, kNuMomY);
const HistAxis axRecoTransMomAng("#delta#alpha_{T, reco} / radians", binsAng, kRecoTransMomAng);
const HistAxis axTransMomAng("#delta#alpha_{T} / radians", binsAng, kTransMomAng);
const HistAxis axRecoTransMomMag("|p_{T, reco}| / GeV", binsMag, kRecoTransMomMag);
const HistAxis axTransMomMag("|p_{T}| / GeV", binsMag, kTransMomMag);

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

  // True kinematics
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

  // Reconstructed proton and lepton kinematics
  NoOscPredictionGenerator genRecoProMomXFhc(axRecoProMomX, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genRecoProMomYFhc(axRecoProMomY, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genRecoProMomZFhc(axRecoProMomZ, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  PredictionInterp predRecoProMomXFhc({}, 0, genRecoProMomXFhc, loadersGArFHC);
  PredictionInterp predRecoProMomYFhc({}, 0, genRecoProMomYFhc, loadersGArFHC);
  PredictionInterp predRecoProMomZFhc({}, 0, genRecoProMomZFhc, loadersGArFHC);
  NoOscPredictionGenerator genRecoLepMomXFhc(axRecoLepMomX, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genRecoLepMomYFhc(axRecoLepMomY, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genRecoLepMomZFhc(axRecoLepMomZ, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  PredictionInterp predRecoLepMomXFhc({}, 0, genRecoLepMomXFhc, loadersGArFHC);
  PredictionInterp predRecoLepMomYFhc({}, 0, genRecoLepMomYFhc, loadersGArFHC);
  PredictionInterp predRecoLepMomZFhc({}, 0, genRecoLepMomZFhc, loadersGArFHC);

  NoOscPredictionGenerator genRecoTransMomMag(axRecoTransMomMag, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genRecoTransMomAng(axRecoTransMomAng, kIsTrueGasFV && kPassND_FHC_NUMU && kNRecoFS==2 && kRecoP==1);
  NoOscPredictionGenerator genTransMomMag(axTransMomMag, kIsTrueGasFV && kProMomX!=0);
  NoOscPredictionGenerator genTransMomAng(axTransMomAng, kIsTrueGasFV && kProMomX!=0);
  PredictionInterp predRecoTransMomMag({}, 0, genRecoTransMomMag, loadersGArFHC);
  PredictionInterp predRecoTransMomAng({}, 0, genRecoTransMomAng, loadersGArFHC);
  PredictionInterp predTransMomMag({}, 0, genTransMomMag, loadersGArFHC);
  PredictionInterp predTransMomAng({}, 0, genTransMomAng, loadersGArFHC);

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

  // Reco kinematics
  TH1 *hRecoLepXFhc = predRecoLepMomXFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRecoLepYFhc = predRecoLepMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRecoLepZFhc = predRecoLepMomZFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hRecoLepXFhc);
  setHistAttr(hRecoLepYFhc);
  setHistAttr(hRecoLepZFhc);
  hRecoLepXFhc->Write("hRecoLepXFhc");
  hRecoLepYFhc->Write("hRecoLepYFhc");
  hRecoLepZFhc->Write("hRecoLepZFhc");
  TH1 *hRecoProXFhc = predRecoProMomXFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRecoProYFhc = predRecoProMomYFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRecoProZFhc = predRecoProMomZFhc.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hRecoProXFhc);
  setHistAttr(hRecoProYFhc);
  setHistAttr(hRecoProZFhc);
  hRecoProXFhc->Write("hRecoProXFhc");
  hRecoProYFhc->Write("hRecoProYFhc");
  hRecoProZFhc->Write("hRecoProZFhc");

  // Momentum imbalance
  TH1 *hRecoTransMomMag = predRecoTransMomMag.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hRecoTransMomAng = predRecoTransMomAng.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hTransMomMag = predTransMomMag.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  TH1 *hTransMomAng = predTransMomAng.Predict(0).FakeData(pot_nd).ToTH1(pot_nd);
  setHistAttr(hRecoTransMomMag);
  setHistAttr(hRecoTransMomAng);
  setHistAttr(hTransMomMag);
  setHistAttr(hTransMomAng);
  hRecoTransMomMag->SetTitle(Form("Reconstructed transverse momentum magnitude in HPgTPC (FHC): %.2g years", years));
  hRecoTransMomAng->SetTitle(Form("Reconstructed transverse momentum angle in HPgTPC (FHC): %.2g years", years));
  hTransMomMag->SetTitle(Form("Transverse momentum magnitude in HPgTPC (FHC): %.2g years", years));
  hTransMomAng->SetTitle(Form("Transverse momentum angle in HPgTPC (FHC): %.2g years", years));
  hRecoTransMomMag->Write("hRecoTransMomMag");
  hRecoTransMomAng->Write("hRecoTransMomAng");
  hTransMomMag->Write("hTransMomMag");
  hTransMomAng->Write("hTransMomAng");

  fout->Close();
  delete fout;
} // transverseMom
