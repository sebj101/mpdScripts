// plotNDLAr.C
// Replicates a few basic plots from the ND CAFs to do with pion multiplicity in the LAr
bool IsInNDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) {
  bool inDeadRegion = false;
  for (int i = -3; i <= 3; ++i) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i * 102.1;
    if (pos_x_cm > cathode_center - 0.75 && pos_x_cm < cathode_center + 0.75)
      inDeadRegion = true;

    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel
    // plane, x2) don't worry about outer boundary because events are only
    // generated in active Ar + insides
    double module_boundary = i * 102.1 + 51.05;
    if (i <= 2 && pos_x_cm > module_boundary - 1.3 &&
	pos_x_cm < module_boundary + 1.3)
      inDeadRegion = true;
  }
  for (int i = 1; i <= 4; ++i) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module
    // wall, x2) module is 102.1cm wide, but only 101.8cm long due to cathode
    // (0.5cm) being absent in length but ArCLight is 0.1cm thicker than pixel
    // plane so it's 0.3cm difference positions are off-set by 0.6 because I
    // defined 0 to be the upstream edge based on the active volume by
    // inspecting a plot, and aparently missed by 3 mm, but whatever add 8mm =
    // 2 pad buffer due to worse position resolution in spatial dimension z
    // compared to timing direction x so total FV gap will be 1.8 + 2*0.8
    // = 3.4cm
    double module_boundary = i * 101.8 - 0.6;
    if (pos_z_cm > module_boundary - 1.7 && pos_z_cm < module_boundary + 1.7)
      inDeadRegion = true;
  }

  return (abs(pos_x_cm) < 300 && abs(pos_y_cm) < 100 && pos_z_cm > 50 &&
	  pos_z_cm < 350 && !inDeadRegion);
}

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

void plotNDLAr(const char *outFile, const char* cafDir="/pnfs/dune/persistent/users/picker24/CAFv4/") 
{
  TFile *fout = new TFile(outFile, "recreate");
  TChain *t = new TChain("caf");
  t->Add(Form("%s/ND_FHC_CAF.root", cafDir));
  t->Add(Form("%s/ND_RHC_CAF.root", cafDir));
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

  TH2D *h2EReco = new TH2D("h2EReco", "Reconstructed vs. true neutrino energy; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco);
  TH2D *h2EReco0Pi = new TH2D("h2EReco0Pi", "Reconstructed vs. true neutrino energy for CC0#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco0Pi);
  TH2D *h2EReco1Pi0 = new TH2D("h2EReco1Pi0", "Reconstructed vs. true neutrino energy for CC1#pi^{0}; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco1Pi0);
  TH2D *h2EReco1Pi = new TH2D("h2EReco1Pi", "Reconstructed vs. true neutrino energy for CC1#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco1Pi);
  TH2D *h2EReco2Pi = new TH2D("h2EReco2Pi", "Reconstructed vs. true neutrino energy for CC2#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco2Pi);
  TH2D *h2EReco3Pi = new TH2D("h2EReco3Pi", "Reconstructed vs. true neutrino energy for CC3#pi; E_{had, true} / GeV; E_{had, reco}; Events", 20, 0., 5., 20, 0., 5.);
  setHistAttr(h2EReco3Pi);

  for (size_t i=0; i<t->GetEntries(); i++) {
    t->GetEntry(i);
    if (i % 100000 == 0) cout<<i<<" of "<<t->GetEntries()<<endl;
    // Is a reco numu and is true numu CC
    if (reco_numu && (muon_contained || muon_tracker) && Ehad_veto<30 && 
	IsInNDFV(vtxx, vtxy, vtxz) && isCC && abs(nuPDG)==14 &&
	((isFHC && reco_q == -1) || (isFHC==0 && reco_q == +1))) {
      h2EReco->Fill(Ev, Ev_reco);
      int nPi = nipip + nipim + nipi0;
      if (nPi == 0) h2EReco0Pi->Fill(Ev, Ev_reco);
      if (nPi == 1) {
	h2EReco1Pi->Fill(Ev, Ev_reco);
	if (nipi0 == 1) h2EReco1Pi0->Fill(Ev, Ev_reco);
      }
      if (nPi == 2) h2EReco2Pi->Fill(Ev, Ev_reco);
      if (nPi == 3) h2EReco3Pi->Fill(Ev, Ev_reco);
    }
  } // Loop over entries
  
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
  fout->Close();
  delete fout;
  /*
  infhc->Close();
  inrhc->Close();
  delete infhc;
  delete inrhc;
  */
} // plotNDLAr
