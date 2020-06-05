// acceptanceGAr.C

bool IsInFV(const double x, const double y, const double z)
{
  double r = sqrt( (y+75.)*(y+75.) + (z-955.)*(z-955.) );
  return ( r < 200. && abs(x) < 200. );
}


void acceptanceGAr(const char* outfile)
{
  TFile *fout = new TFile(outfile, "recreate");

  TChain *tch = new TChain("caf");
  tch->Add("/dune/data/users/sbjones/gasTpcCAF/v10/CAF_FHC.root");
  tch->Add("/dune/data/users/sbjones/gasTpcCAF/v10/CAF_RHC.root");

  int mode;
  int reco_numu;
  int LepPDG;
  double vtx_x, vtx_y, vtx_z;
  double Ev, Q2;
  tch->SetBranchAddress("mode", &mode);
  tch->SetBranchAddress("reco_numu", &reco_numu);
  tch->SetBranchAddress("vtx_x", &vtx_x);
  tch->SetBranchAddress("vtx_y", &vtx_y);
  tch->SetBranchAddress("vtx_z", &vtx_z);
  tch->SetBranchAddress("Ev", &Ev);
  tch->SetBranchAddress("Q2", &Q2);
  tch->SetBranchAddress("LepPDG", &LepPDG);

  TH2D *h2TrueCCQE     = new TH2D("h2TrueCCQE", "CCQE, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2TrueCCQE_fhc = new TH2D("h2TrueCCQE_fhc", "CCQE, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2TrueCCQE_rhc = new TH2D("h2TrueCCQE_rhc", "CCQE, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2RecoCCQE     = new TH2D("h2RecoCCQE", "CCQE, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2RecoCCQE_fhc = new TH2D("h2RecoCCQE_fhc", "CCQE, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2RecoCCQE_rhc = new TH2D("h2RecoCCQE_rhc", "CCQE, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);

  TH2D *h2True     = new TH2D("h2True", "#nu_{#mu} CC, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2True_fhc = new TH2D("h2True_fhc", "#nu_{#mu} CC, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2True_rhc = new TH2D("h2True_rhc", "#nu_{#mu} CC, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Reco     = new TH2D("h2Reco", "#nu_{#mu} CC, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Reco_fhc = new TH2D("h2Reco_fhc", "#nu_{#mu} CC, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Reco_rhc = new TH2D("h2Reco_rhc", "#nu_{#mu} CC, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);

  for (int i=0; i<tch->GetEntries(); i++) {
    tch->GetEntry(i);
    if (i % 10000 == 0) cout<<i<<" of "<<tch->GetEntries()<<endl;
    if (abs(LepPDG)==13 && IsInFV(vtx_x, vtx_y, vtx_z)) {
      if (mode == 1) {
	h2TrueCCQE->Fill(Ev, Q2);
	(i < tch->GetEntries()/2.) ? h2TrueCCQE_fhc->Fill(Ev, Q2) : h2TrueCCQE_rhc->Fill(Ev, Q2);
	if (reco_numu == 1) {
	  h2RecoCCQE->Fill(Ev, Q2);
	  (i < tch->GetEntries()/2.) ? h2RecoCCQE_fhc->Fill(Ev, Q2) : h2RecoCCQE_rhc->Fill(Ev, Q2);
	}
      } // Is QE
      h2True->Fill(Ev, Q2);
      (i < tch->GetEntries()/2.) ? h2True_fhc->Fill(Ev, Q2) : h2True_rhc->Fill(Ev, Q2);
      if (reco_numu == 1) {
	h2Reco->Fill(Ev, Q2);
	(i < tch->GetEntries()/2.) ? h2Reco_fhc->Fill(Ev, Q2) : h2Reco_rhc->Fill(Ev, Q2);
      }
    } // True numu CC and in FV
  } // Loop over entries

  TH2D *h2AccCCQE = new TH2D("h2AccCCQE", "CCQE, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2AccCCQE_fhc = new TH2D("h2AccCCQE_fhc", "CCQE, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2AccCCQE_rhc = new TH2D("h2AccCCQE_rhc", "CCQE, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Acc = new TH2D("h2Acc", "#nu_{#mu} CC, HPgTPC Near Detector; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Acc_fhc = new TH2D("h2Acc_fhc", "#nu_{#mu} CC, HPgTPC Near Detector, FHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);
  TH2D *h2Acc_rhc = new TH2D("h2Acc_rhc", "#nu_{#mu} CC, HPgTPC Near Detector, RHC; True E_{#nu}; True Q^{2}", 25, 0, 5, 25, 0, 5);

  h2AccCCQE->Divide(h2RecoCCQE, h2TrueCCQE);
  h2AccCCQE_fhc->Divide(h2RecoCCQE_fhc, h2TrueCCQE_fhc);
  h2AccCCQE_rhc->Divide(h2RecoCCQE_rhc, h2TrueCCQE_rhc);
  h2Acc->Divide(h2Reco, h2True);
  h2Acc_fhc->Divide(h2Reco_fhc, h2True_fhc);
  h2Acc_rhc->Divide(h2Reco_rhc, h2True_rhc);

  fout->cd();
  h2TrueCCQE->Write();
  h2TrueCCQE_fhc->Write();
  h2TrueCCQE_rhc->Write();
  h2RecoCCQE->Write();
  h2RecoCCQE_fhc->Write();
  h2RecoCCQE_rhc->Write();

  h2True->Write();
  h2True_fhc->Write();
  h2True_rhc->Write();
  h2Reco->Write();
  h2Reco_fhc->Write();
  h2Reco_rhc->Write();

  h2AccCCQE->Write();
  h2AccCCQE_fhc->Write();
  h2AccCCQE_rhc->Write();

  h2Acc->Write();
  h2Acc_fhc->Write();
  h2Acc_rhc->Write();


  fout->Close();
  delete fout;
}
