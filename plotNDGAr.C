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

void plotNDGAr(const char *outFile, const char* cafDir="/dune/data/users/sbjones/gasTpcCAF/")
{
  TFile *fout = new TFile(outFile, "recreate");
  TChain *t = new TChain("caf");
  t->Add(Form("%s/CAF_FHC.root", cafDir));
  t->Add(Form("%s/CAF_RHC.root", cafDir));

  int isFHC, isCC, reco_q, reco_numu, muon_contained, muon_tracker, nuPDG;
  int nipip, nipim, nipi0;
  double Ehad_veto;
  double Ev, LepE, Ev_reco, Elep_reco;
  double vtxx, vtxy, vtxz;
  t->SetBranchAddress("reco_numu", &reco_numu);
  t->SetBranchAddress("reco_q", &reco_q);
  t->SetBranchAddress("isCC", &isCC);
  t->SetBranchAddress("isFHC", &isFHC);
  t->SetBranchAddress("muon_contained", &muon_contained);
  t->SetBranchAddress("muon_tracker", &muon_tracker);
  t->SetBranchAddress("Ehad_veto", &Ehad_veto);
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

  for (size_t i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    if (i % 100000 == 0) cout<<i<<" of "<<t->GetEntries()<<endl;
    // TO DO: better selection cuts
    if (isCC && abs(nuPDG)==14 && ((isFHC==1 && reco_q == -1) || (isFHC==0 && reco_q == +1)) && IsInFV(vtxx, vtxy, vtxz) && reco_numu==1 && muon_tracker==1) {
      h2EReco->Fill(Ev, Ev_reco);
      hTrueReco->Fill(Ev_reco - Ev);
      int nPi = nipip + nipim + nipi0;
      if (nPi == 0) {
	h2EReco0Pi->Fill(Ev, Ev_reco);
	hTrueRecoCC0Pi->Fill(Ev_reco - Ev);
      }
      if (nPi == 1) {
	h2EReco1Pi->Fill(Ev, Ev_reco);
	hTrueRecoCC1Pi->Fill(Ev_reco - Ev);
	if (nipi0 == 1) h2EReco1Pi0->Fill(Ev, Ev_reco);
      }
      if (nPi == 2) {
	h2EReco2Pi->Fill(Ev, Ev_reco);
	hTrueRecoCC2Pi->Fill(Ev_reco - Ev);
      }
      if (nPi == 3) {
	h2EReco3Pi->Fill(Ev, Ev_reco);
	hTrueRecoCC3Pi->Fill(Ev_reco - Ev);
      }
    }
  }  
  TLine *l = new TLine(0, 0, 5, 5);
  l->SetLineWidth(3);
  l->SetLineStyle(2);
  fout->cd();
  l->Write("line");

  h2EReco->Write();
  h2EReco0Pi->Write();
  h2EReco1Pi0->Write();
  h2EReco1Pi->Write();
  h2EReco2Pi->Write();
  h2EReco3Pi->Write();

  hTrueReco->Write();
  hTrueRecoCC0Pi->Write();
  hTrueRecoCC1Pi->Write();
  hTrueRecoCC2Pi->Write();
  hTrueRecoCC3Pi->Write();

  fout->Close();
  delete fout;
} // plotNDGAr
