// plotNDGAr.C
void setHistAttr(TH2D *h2)
{
  h2->Sumw2();
  h2->GetXaxis()->SetTitleSize(.05);
  h2->GetYaxis()->SetTitleSize(.05);
  h2->GetZaxis()->SetTitleSize(.05);
  h2->GetXaxis()->SetLabelSize(.05);
  h2->GetYaxis()->SetLabelSize(.05);
  h2->GetZaxis()->SetLabelSize(.05);
  h2->SetOption("colz");
}

void setHistAttr(TH1D *h)
{
  h->Sumw2();
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitleSize(.05);
  h->GetYaxis()->SetTitleSize(.05);
  h->GetXaxis()->SetLabelSize(.05);
  h->GetYaxis()->SetLabelSize(.05);
  h->SetOption("hist");
}

bool IsInFV(const double x, const double y, const double z)
{
  double r = sqrt( (y+75.)*(y+75.) + (z-955.)*(z-955.) );
  return ( r < 200. && abs(x) < 200. );
}

void plotNDGAr(const char *outFile, const char* cafDir="/dune/data/users/sbjones/gasTpcCAF/v3/")
{
  TFile *fout = new TFile(outFile, "recreate");
  TChain *t = new TChain("caf");
  t->Add(Form("%s/CAF_FHC.root", cafDir));
  t->Add(Form("%s/CAF_RHC.root", cafDir));

  int isFHC, isCC, reco_q, reco_numu, muon_contained, muon_tracker, nuPDG;
  int nipip, nipim, nipi0;
  int gastpc_pi_pl_mult, gastpc_pi_min_mult, gastpc_pi_0_mult;
  double Ehad_veto;
  double Ev, LepE, Ev_reco, Elep_reco;
  double vtxx, vtxy, vtxz;
  int nFSP;
  int pdg[70];
  double partEvReco[70];
  double ptrue[70];
  t->SetBranchAddress("reco_numu", &reco_numu);
  t->SetBranchAddress("reco_q", &reco_q);
  t->SetBranchAddress("isCC", &isCC);
  t->SetBranchAddress("isFHC", &isFHC);
  t->SetBranchAddress("muon_contained", &muon_contained);
  t->SetBranchAddress("muon_tracker", &muon_tracker);
  t->SetBranchAddress("nuPDG", &nuPDG);
  t->SetBranchAddress("nipip", &nipip);
  t->SetBranchAddress("nipim", &nipim);
  t->SetBranchAddress("nipi0", &nipi0);
  t->SetBranchAddress("vtx_x", &vtxx);
  t->SetBranchAddress("vtx_y", &vtxy);
  t->SetBranchAddress("vtx_z", &vtxz);
  t->SetBranchAddress("Ev", &Ev);
  t->SetBranchAddress("Ev_reco", &Ev_reco);
  t->SetBranchAddress("Elep_reco", &Elep_reco);
  t->SetBranchAddress("LepE", &LepE);
  t->SetBranchAddress("nFSP", &nFSP);
  t->SetBranchAddress("pdg", pdg);
  t->SetBranchAddress("ptrue", ptrue);
  t->SetBranchAddress("partEvReco", partEvReco);
  t->SetBranchAddress("gastpc_pi_pl_mult", &gastpc_pi_pl_mult);
  t->SetBranchAddress("gastpc_pi_min_mult", &gastpc_pi_min_mult);
  t->SetBranchAddress("gastpc_pi_0_mult", &gastpc_pi_0_mult);

  // True vs reconstructed pion
  TH2D *h2PiTrueReco = new TH2D("h2PiTrueReco", "Reconstructed vs true #pi multiplicity; N_{#pi, true}; N_{#pi, reco}", 5, -0.5, 4.5, 5, -0.5, 4.5);
  // True vs reco
  TH2D *h2EReco = new TH2D("h2EReco", "Reconstructed vs. true neutrino energy; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco);
  TH2D *h2EReco0Pi = new TH2D("h2EReco0Pi", "Reconstructed vs. true neutrino energy for CC0#pi; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco0Pi);
  TH2D *h2EReco1Pi0 = new TH2D("h2EReco1Pi0", "Reconstructed vs. true neutrino energy for CC1#pi^{0}; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco1Pi0);
  TH2D *h2EReco1Pi = new TH2D("h2EReco1Pi", "Reconstructed vs. true neutrino energy for CC1#pi; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco1Pi);
  TH2D *h2EReco2Pi = new TH2D("h2EReco2Pi", "Reconstructed vs. true neutrino energy for CC2#pi; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco2Pi);
  TH2D *h2EReco3Pi = new TH2D("h2EReco3Pi", "Reconstructed vs. true neutrino energy for CC3#pi; E_{#nu, true} / GeV; E_{#nu, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco3Pi);
  // Reco - true
  TH1D *hTrueReco = new TH1D("hTrueReco", "E_{#nu, reco} - E_{#nu, true}; E_{#nu, reco} - E_{#nu, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueReco);
  TH1D *hTrueRecoCC0Pi = new TH1D("hTrueRecoCC0Pi", "E_{#nu, reco} - E_{#nu, true} for CC0#pi; E_{#nu, reco} - E_{#nu, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoCC0Pi);
  TH1D *hTrueRecoCC1Pi = new TH1D("hTrueRecoCC1Pi", "E_{#nu, reco} - E_{#nu, true} for CC1#pi; E_{#nu, reco} - E_{#nu, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoCC1Pi);
  TH1D *hTrueRecoCC2Pi = new TH1D("hTrueRecoCC2Pi", "E_{#nu, reco} - E_{#nu, true} for CC2#pi; E_{#nu, reco} - E_{#nu, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoCC2Pi);
  TH1D *hTrueRecoCC3Pi = new TH1D("hTrueRecoCC3Pi", "E_{#nu, reco} - E_{#nu, true} for CC3#pi; E_{#nu, reco} - E_{#nu, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoCC3Pi);
  // Hadron and lepton separately
  TH2D *h2EHadReco0Pi = new TH2D("h2EHadReco0Pi", "Reconstructed vs. true hadronic energy for CC0#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2EHadReco0Pi);
  TH2D *h2EHadReco1Pi0 = new TH2D("h2EHadReco1Pi0", "Reconstructed vs. true hadronic energy for CC1#pi^{0}; E_{had, true} / GeV; E_{had, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2EHadReco1Pi0);
  TH2D *h2EHadReco1Pi = new TH2D("h2EHadReco1Pi", "Reconstructed vs. true hadronic energy for CC1#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2EHadReco1Pi);
  TH2D *h2EHadReco2Pi = new TH2D("h2EHadReco2Pi", "Reconstructed vs. true hadronic energy for CC2#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2EHadReco2Pi);
  TH2D *h2EHadReco3Pi = new TH2D("h2EHadReco3Pi", "Reconstructed vs. true hadronic energy for CC3#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2EHadReco3Pi);
  // Lepton
  TH2D *h2ELepReco0Pi = new TH2D("h2ELepReco0Pi", "Reconstructed vs. true leptonic energy for CC0#pi; E_{lep, true} / GeV; E_{lep, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2ELepReco0Pi);
  TH2D *h2ELepReco1Pi0 = new TH2D("h2ELepReco1Pi0", "Reconstructed vs. true leptonic energy for CC1#pi^{0}; E_{lep, true} / GeV; E_{lep, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2ELepReco1Pi0);
  TH2D *h2ELepReco1Pi = new TH2D("h2ELepReco1Pi", "Reconstructed vs. true leptonic energy for CC1#pi; E_{lep, true} / GeV; E_{lep, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2ELepReco1Pi);
  TH2D *h2ELepReco2Pi = new TH2D("h2ELepReco2Pi", "Reconstructed vs. true leptonic energy for CC2#pi; E_{lep, true} / GeV; E_{lep, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2ELepReco2Pi);
  TH2D *h2ELepReco3Pi = new TH2D("h2ELepReco3Pi", "Reconstructed vs. true leptonic energy for CC3#pi; E_{lep, true} / GeV; E_{lep, reco}; Events", 20., 0., 4., 20., 0., 4.);
  setHistAttr(h2ELepReco3Pi);
  // True  - reco
  TH1D *hTrueRecoHadCC0Pi = new TH1D("hTrueRecoHadCC0Pi", "E_{had, reco} - E_{had, true} for CC0#pi; E_{had, reco} - E_{had, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoHadCC0Pi);
  TH1D *hTrueRecoHadCC1Pi = new TH1D("hTrueRecoHadCC1Pi", "E_{had, reco} - E_{had, true} for CC1#pi; E_{had, reco} - E_{had, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoHadCC1Pi);
  TH1D *hTrueRecoHadCC2Pi = new TH1D("hTrueRecoHadCC2Pi", "E_{had, reco} - E_{had, true} for CC2#pi; E_{had, reco} - E_{had, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoHadCC2Pi);
  TH1D *hTrueRecoHadCC3Pi = new TH1D("hTrueRecoHadCC3Pi", "E_{had, reco} - E_{had, true} for CC3#pi; E_{had, reco} - E_{had, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoHadCC3Pi);
  // Lepton
  TH1D *hTrueRecoLepCC0Pi = new TH1D("hTrueRecoLepCC0Pi", "E_{lep, reco} - E_{lep, true} for CC0#pi; E_{lep, reco} - E_{lep, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoLepCC0Pi);
  TH1D *hTrueRecoLepCC1Pi = new TH1D("hTrueRecoLepCC1Pi", "E_{lep, reco} - E_{lep, true} for CC1#pi; E_{lep, reco} - E_{lep, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoLepCC1Pi);
  TH1D *hTrueRecoLepCC2Pi = new TH1D("hTrueRecoLepCC2Pi", "E_{lep, reco} - E_{lep, true} for CC2#pi; E_{lep, reco} - E_{lep, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoLepCC2Pi);
  TH1D *hTrueRecoLepCC3Pi = new TH1D("hTrueRecoLepCC3Pi", "E_{lep, reco} - E_{lep, true} for CC3#pi; E_{lep, reco} - E_{lep, true} / GeV; Events", 100, -2, 2);
  setHistAttr(hTrueRecoLepCC3Pi);
  // Spatial distribution
  TH2D *h2TrueRecoHad0Pi = new TH2D("h2TrueRecoHad0Pi", "E_{had, reco} - E_{had, true} for CC0#pi; E_{had, reco} - E_{had, true} / GeV; r / cm; Events", 100, -2, 2, 100, 0, 200);
  setHistAttr(h2TrueRecoHad0Pi);
  TH2D *h2TrueRecoHad1Pi = new TH2D("h2TrueRecoHad1Pi", "E_{had, reco} - E_{had, true} for CC1#pi; E_{had, reco} - E_{had, true} / GeV; r / cm; Events", 100, -2, 2, 100, 0, 200);
  setHistAttr(h2TrueRecoHad1Pi);
  TH2D *h2TrueRecoHad2Pi = new TH2D("h2TrueRecoHad2Pi", "E_{had, reco} - E_{had, true} for CC2#pi; E_{had, reco} - E_{had, true} / GeV; r / cm; Events", 100, -2, 2, 100, 0, 200);
  setHistAttr(h2TrueRecoHad2Pi);
  TH2D *h2TrueRecoHad3Pi = new TH2D("h2TrueRecoHad3Pi", "E_{had, reco} - E_{had, true} for CC3#pi; E_{had, reco} - E_{had, true} / GeV; r / cm; Events", 100, -2, 2, 100, 0, 200);
  setHistAttr(h2TrueRecoHad3Pi);

  // True - reco for particle types
  TH1D *hTrueRecoPro = new TH1D("hTrueRecoPro", "Proton reco - true; E_{P, reco} - E_{P, true} / GeV; Events", 100, -2, 2);
  TH1D *hTrueRecoProLowE = new TH1D("hTrueRecoProLowE", "Proton reco - true (low energy); E_{P, reco} - E_{P, true} / GeV; Events", 100, -2, 2);
  TH1D *hTrueRecoPip = new TH1D("hTrueRecoPip", "Pion reco - true; E_{#pi^{+}, reco} - E_{#pi^{+}, true} / GeV; Events", 100, -2, 2);
  TH1D *hTrueRecoPim = new TH1D("hTrueRecoPim", "Pion reco - true; E_{#pi^{-}, reco} - E_{#pi^{-}, true} / GeV; Events", 100, -2, 2);
  TH1D *hTrueRecoPi0 = new TH1D("hTrueRecoPi0", "Pion reco - true; E_{#pi^{0}, reco} - E_{#pi^{0}, true} / GeV; Events", 100, -2, 2);

  for (size_t i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    if (i % 100000 == 0) cout<<i<<" of "<<t->GetEntries()<<endl;
    // TO DO: better selection cuts
    if (isCC && abs(nuPDG)==14 && ((isFHC==1 && reco_q == -1) || (isFHC==0 && reco_q == +1)) && IsInFV(vtxx, vtxy, vtxz) && reco_numu==1 && muon_tracker==1) {
      h2EReco->Fill(Ev, Ev_reco);
      hTrueReco->Fill(Ev_reco - Ev);
      int nPi = nipip + nipim + nipi0;
      double r = sqrt( (vtxy+75.)*(vtxy+75.) + (vtxz-955.)*(vtxz-955.) );

      for (int n=0; n<nFSP; n++) {
	if (pdg[n] == 211) 
	  hTrueRecoPip->Fill(partEvReco[n] - sqrt(ptrue[n]*ptrue[n] + 0.1395*0.1395));
	else if (pdg[n] == -211) 
	  hTrueRecoPim->Fill(partEvReco[n] - sqrt(ptrue[n]*ptrue[n] + 0.1395*0.1395));
	else if (pdg[n] == 111) 
	  hTrueRecoPi0->Fill(partEvReco[n] - sqrt(ptrue[n]*ptrue[n] + 0.135*0.135));
       	else if (pdg[n] == 2212) {
	  hTrueRecoPro->Fill(partEvReco[n] - sqrt(ptrue[n]*ptrue[n]));
	  if (ptrue[n] < 1.) {
	    hTrueRecoProLowE->Fill(partEvReco[n] - sqrt(ptrue[n]*ptrue[n]));
	  }
	}
      }

      h2PiTrueReco->Fill(nipip+nipim+nipi0, gastpc_pi_pl_mult+gastpc_pi_min_mult+gastpc_pi_0_mult);

      if (nPi == 0) {
	h2EReco0Pi->Fill(Ev, Ev_reco);
	h2EHadReco0Pi->Fill(Ev - LepE, Ev_reco - Elep_reco);
	hTrueRecoCC0Pi->Fill(Ev_reco - Ev);
	hTrueRecoHadCC0Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE));
	h2ELepReco0Pi->Fill(LepE, Elep_reco);
	hTrueRecoLepCC0Pi->Fill(Elep_reco - LepE);
	h2TrueRecoHad0Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE), r);
      }
      if (nPi == 1) {
	h2EReco1Pi->Fill(Ev, Ev_reco);
	h2EHadReco1Pi->Fill(Ev - LepE, Ev_reco - Elep_reco);
	hTrueRecoCC1Pi->Fill(Ev_reco - Ev);
	hTrueRecoHadCC1Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE));
	h2ELepReco1Pi->Fill(LepE, Elep_reco);
	hTrueRecoLepCC1Pi->Fill(Elep_reco - LepE);
	h2TrueRecoHad1Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE), r);
	if (nipi0 == 1) h2EReco1Pi0->Fill(Ev, Ev_reco);
      }
      if (nPi == 2) {
	h2EReco2Pi->Fill(Ev, Ev_reco);
	h2EHadReco2Pi->Fill(Ev - LepE, Ev_reco - Elep_reco);
	hTrueRecoCC2Pi->Fill(Ev_reco - Ev);
	hTrueRecoHadCC2Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE));
	h2ELepReco2Pi->Fill(LepE, Elep_reco);
	hTrueRecoLepCC2Pi->Fill(Elep_reco - LepE);
	h2TrueRecoHad2Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE), r);
      }
      if (nPi == 3) {
	h2EReco3Pi->Fill(Ev, Ev_reco);
	h2EHadReco3Pi->Fill(Ev - LepE, Ev_reco - Elep_reco);
	hTrueRecoCC3Pi->Fill(Ev_reco - Ev);
	hTrueRecoHadCC3Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE));
	h2ELepReco3Pi->Fill(LepE, Elep_reco);
	hTrueRecoLepCC3Pi->Fill(Elep_reco - LepE);
	h2TrueRecoHad3Pi->Fill((Ev_reco-Elep_reco) - (Ev-LepE), r);
      }
    }
  }  
  TLine *l = new TLine(0, 0, 5, 5);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  fout->cd();
  l->Write("line");

  h2PiTrueReco->Write();

  h2EReco->Write();
  h2EReco0Pi->Write();
  h2EReco1Pi0->Write();
  h2EReco1Pi->Write();
  h2EReco2Pi->Write();
  h2EReco3Pi->Write();

  h2EHadReco0Pi->Write();
  h2EHadReco1Pi->Write();
  h2EHadReco2Pi->Write();
  h2EHadReco3Pi->Write();

  hTrueReco->Write();
  hTrueRecoCC0Pi->Write();
  hTrueRecoCC1Pi->Write();
  hTrueRecoCC2Pi->Write();
  hTrueRecoCC3Pi->Write();

  hTrueRecoHadCC0Pi->Write();
  hTrueRecoHadCC1Pi->Write();
  hTrueRecoHadCC2Pi->Write();
  hTrueRecoHadCC3Pi->Write();

  h2TrueRecoHad0Pi->Write();
  h2TrueRecoHad1Pi->Write();
  h2TrueRecoHad2Pi->Write();
  h2TrueRecoHad3Pi->Write();

  hTrueRecoLepCC0Pi->Write();
  hTrueRecoLepCC1Pi->Write();
  hTrueRecoLepCC2Pi->Write();
  hTrueRecoLepCC3Pi->Write();

  h2ELepReco0Pi->Write();
  h2ELepReco1Pi->Write();
  h2ELepReco2Pi->Write();
  h2ELepReco3Pi->Write();

  hTrueRecoPip->Write();
  hTrueRecoPim->Write();
  hTrueRecoPi0->Write();
  hTrueRecoPro->Write();
  hTrueRecoProLowE->Write();

  fout->Close();
  delete fout;
} // plotNDGAr
