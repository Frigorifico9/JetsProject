//pratchett
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

  int filas = 12;
  int columnas = 3;

  TH2D* all_histos[filas][columnas];
  TH1D* hist_final[filas][columnas];
  int mult_range[]   = {234, 204, 188, 168, 152, 140, 128, 116, 104, 92, 76, 0};
  const Char_t* Aj_range[] = {"Inc", "A", "B"};
  const Char_t* lol[] = {"Inclussive Aj", "Aj>0.22", "Aj<0.22"};

  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      all_histos[f][c] = (TH2D*) mydata->Get(Form("hist_mult>=%d_Aj%s",mult_range[f],Aj_range[c]));
      hist_final[f][c] = new TH1D( Form("hist_final_mult>=%d_Aj%s",mult_range[f],Aj_range[c]), Form("Projected pT vs Delta with %s",lol[c]), 9, 0, 1.8);
    }
  }
  cout<<1<<endl;

  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      for (int n = 0; n < 9; ++n)
      {
        all_histos[f][c]->GetXaxis()->SetRangeUser( 0.2*n, 0.2*(n+1) );//we first define each ring
        hist_final[f][c]->SetBinContent( n+1, all_histos[f][c]->GetMean(2) ); //then take the average in that ring
        hist_final[f][c]->SetBinError( n+1, all_histos[f][c]->GetMeanError(2) );
      }
    }
  }
  cout<<2<<endl;

  post_histos->cd();
  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      hist_final[f][c] ->Write();
    }
  }
  cout<<3<<endl;
  post_histos->Close();
  
  return 0;
}