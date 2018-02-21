//halo 1
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

using namespace Pythia8;
using namespace fastjet;

int main() 
{
  // Generator. Shorthand for event.
  Pythia pythia;
  pythia.readFile("seetings.cmnd");
  int Nev = pythia.mode("Main:numberOfEvents");
  Event& event = pythia.event;
  pythia.init();

  TFile* histos = new TFile("histos.root","recreate");
  int N = 10000; //number of bins in hist_mult
  TH1D*  hist_mult = new TH1D("hist_mult", "The multiplicity", N, 0, 500);
//First we want the multiplicities, the mounf of particles with pT>0 and |eta|>2 in each event
  int mult = 0;
  for (int iEvent = 0; iEvent < Nev; ++iEvent) 
  {
    mult = 0;
    if (!pythia.next()) continue;

    for (int i = 0; i < event.size(); ++i) if (event[i].isCharged() && event[i].isFinal())
    {
      ++mult;
    }
    hist_mult->Fill(mult);
  }

  //Now we have to define the sets of events according to their multiplicity
  double sets[] = {0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; //this are multiplicity sets we want to find
  int nset = 1; 
  double I = hist_mult->Integral(); //this is the 100%
  double II = 0;
  double ii = 0; 
  cout<<"The integral is "<<I<<endl;

  for (int i = N-1; i > 0; --i) //we'll loop through the bins in hist_mult
  {
    hist_mult->GetXaxis()->SetRange(i,N); //first we declare the range to consider
    II = hist_mult->Integral();
    hist_mult->GetXaxis()->SetRange(i+1,N);
    ii = hist_mult->Integral();;
    if ( ii < (I*sets[nset+1]) && II > (I*sets[nset]) )
    //if ( hist_mult->Integral() > (I*sets[nset]) ) 
    {
      //hist_mult->GetXaxis()->SetRange(i-1,N);
      cout<<"La integral en este rango es "<<ii<<" y se supone que es menor o igual a "<<I*sets[nset]<<endl;
      cout<<"El "<<sets[nset]*100<<" porciento de los eventos tiene una multiplicidad mayor a: "<<endl;
      cout<<(i+1)*500/N<<endl;
      cout<<"Lo cual corresponde al bin "<<i+1<<endl;
      cout<<endl;
      ++nset;
    }
  }

  //
  pythia.stat();
  hist_mult->Sumw2();
  //

  histos->cd();
  hist_mult->Write();
  histos->Close();

  return 0;
}
