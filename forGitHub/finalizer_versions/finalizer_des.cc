//sorrow
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
#include "TTree.h"
//

int main() 
{
  TFile* post_histos = new TFile("post_histos.root","recreate");

  TFile* mydata = TFile::Open("mergedHistos.root");

  //vector< vector<double> > dP(9, vector<double>(5));
  vector< vector<double> > TdP(9, vector<double>(5));

  TTree* tree = (TTree*) mydata->Get("tree");

  TBranch *bdP = tree->SetBranchAddress("dP",&dP);

  dP = bdP->Get("values");

  for (int d = 0; d <= 9; ++d)
  {
    for (int p = 0; p < 5; ++p)
    {
      cout<<dP[d][p];
    }
    cout<<endl;
  }

/*/
  TdP = (vector< vector<double> >) mydata->Get("TdP");

  double delta_edges[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8};

  TH1D* hist_final  = new TH1D("hist_final", "Whatever this is", 10, delta_edges ); 
  
  double final[9];
  for (int d = 0; d <= 9; ++d)
  {
    for (int p = 0; p < 5; ++p)
    {
      dP[d][p] = dP[d][p]/TdP[d][p];
      final[d] = final[d] + dP[d][p]/TdP[d][p];
    }
    hist_final->SetBinContent(bien,final[d]);
  }

  hist_final->Sumw2();
  //hist_final->SetOption("bar");

  post_histos->cd();
  hist_final->Write();
  post_histos->Close();

  //clean up
  delete hist_final;

  /*/
  
  return 0;
}