//martinillo
#include <iostream>
using namespace std;

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
//
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TProfile.h"
//

int main() 
{
  TFile* post_histos = new TFile("post_histos.root","recreate");

  TFile* mydata      = TFile::Open("mergedHistos.root");

  TH2D* hist_pTdelta = (TH2D*) mydata->Get("hist_pTdelta"); 

  TH2D* hist_final   = new TH2D("hist_final", "Whatever this is", 10, 0, 1.8, 1000, -300, 300); 


  double delta_edges[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8 };

  for (int n = 0; n <= 9; ++n)
  {
    hist_pTdelta->GetXaxis()->SetRange(delta_edges[n],delta_edges[n+1]);//to do this we first define each ring
    hist_final->Fill( ( delta_edges[n]+delta_edges[n+1] )/2, hist_pTdelta->GetMean(2) );
  }

  hist_final->Sumw2();

  TH1D* hist_projX = hist_final->ProjectionX();
  TH1D* hist_projY = hist_final->ProjectionY();

  hist_projX->SetOption("bar");
  hist_projY->SetOption("bar");

  post_histos->cd();
  hist_projX->Write();
  hist_projY->Write();
  post_histos->Close();

  //clean up
  delete hist_projX;
  delete hist_projY;
  delete hist_final;
  delete hist_pTdelta;
  
  return 0;
}