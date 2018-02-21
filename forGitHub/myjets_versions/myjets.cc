//dolor
#include <iostream>
using namespace std;

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/FastJet3.h"
//
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TAxis.h"
//

using namespace Pythia8;
using namespace fastjet;

int main() 
{
  TFile* histos             = new TFile("histos.root","recreate");
  TH1F*  hist_leading_pT    = new TH1F("hist_leading_pT", "Momentum of leading Jet", 1000, 0, 500);
  TH1F*  hist_subleading_pT = new TH1F("hist_subleading_pT", "Momentum of subleading Jet", 1000, 0, 500);
  TH1F*  hist_axis_pT       = new TH1F("hist_axis_pT", "Momentum of axis", 1000, 0, 500);
  //TH1F*  hist_jetTrack_pT   = new TH1F("hist_jetTrack_pT", "Momentum of tracks inside jets", 1000, 0, 300);
  //TH1F*  hist_track_pT      = new TH1F("hist_track_pT", "Momentum of tracks outside jets", 1000, 0, 500);
  TH1F*  hist_missing_pT    = new TH1F("hist_missing_pT", "Missing pT", 1000, 0, 500);
  TH1F*  hist_all_phi       = new TH1F("hist_all_phi", "Phi of all jets", 1000, 0, 7);
  TH1F*  hist_selected_phi  = new TH1F("hist_selected_phi", "Phi of selected jets", 1000, 0, 7);

  // Select common parameters for SlowJet and FastJet analyses.
  int    power          = -1; // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R              = 0.4;    
  double pTMin          = 0.0;
  double eta            = 2.0;
  double delta          = 0;
  double deltaMin       = 0.;
  double deltaMax       = 3.6;
  double pi             = 3.14159265359;
  double shadow         = 0;
  double pT_leading     = 120;
  double pT_subleading  = 50;
  bool   notdijet       = true;
  double subleading_r   = 0;
  double subleading_t   = 0;
  double subleading_phi = 0;
  double subleading_px  = 0;
  double subleading_py  = 0;
  double subleading_pz  = 0;
  double missingPT      = 0;

  // Generator. Shorthand for event.
  Pythia pythia;
  pythia.readFile("seetings.cmnd");
  int Nev = pythia.mode("Main:numberOfEvents");
  Event& event = pythia.event;
  pythia.init();

  JetDefinition jetDef(genkt_algorithm, R, power);
  vector <PseudoJet> fjInputs;

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < Nev; ++iEvent) 
  {
    if (!pythia.next()) continue;
    fjInputs.resize(0);

    shadow         = 0;
    notdijet       = true;
    subleading_r   = 0;
    subleading_t   = 0;
    subleading_phi = 0;
    delta          = 0;
    PseudoJet leading;
    PseudoJet subleading;
    vector <PseudoJet> inclusiveJets, sortedJets;
    PseudoJet axis;

    //First we need to create the jets
    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) 
    {
      if (!event[i].isVisible() || event[i].isNeutral() ) continue;
      PseudoJet particleTemp = event[i];
      fjInputs.push_back(particleTemp);
    }
    //cout<<"The particles' data was stored"<<endl;

    ClusterSequence clustSeq(fjInputs, jetDef); //it takes what we found, fjInuts, and transforms it into jet that fit our jetDef
    inclusiveJets = clustSeq.inclusive_jets(pTMin);
    sortedJets    = sorted_by_pt(inclusiveJets); //return a vector of jets sorted into decreasing pT2
    //cout<<"The Jets have been created and sorted"<<endl;

    if (sortedJets.size()<2) continue;

    //find most energetic jet
    for (int j = 0; j < sortedJets.size(); ++j)
    {
      leading  = sortedJets[j];
      if ( leading.perp()>pT_leading && abs(leading.pseudorapidity())<=eta )
      {
        //cout<<"The leading Jet was found"<<endl;
        //once we found the most energetic is time to infd the second most energetic
        for (int k = j+1; k < sortedJets.size(); ++k)
        {
          subleading = sortedJets[k];
          hist_all_phi->Fill(leading.phi()-sortedJets[k].phi());
          if ( subleading.perp()>pT_subleading && abs(subleading.pseudorapidity())<=eta && abs(leading.phi()-subleading.phi())>((5*pi)/6) )
          {
            //cout<<"The subleading Jet was found"<<endl;
            hist_leading_pT->Fill(leading.perp());
            hist_subleading_pT->Fill(subleading.perp());
            hist_selected_phi->Fill(leading.phi()-subleading.phi());
            notdijet = false; //we have found the two jets, hence this is a dijet event
            //we have what we've come for
            break;
          }
        }
        //we have what we've come for
        break;
      }
    }

    if (notdijet) continue;//lets avoid nondijet events

//####################################################################################################

    //Now, we have to reflect subleading around the origin, which requieres a bit of math
    //first the spherical components
    subleading_r   = sqrt(pow(subleading.px(),2)+pow(subleading.py(),2)+pow(subleading.pz(),2));
    subleading_t   = atan(subleading.py()/subleading.px());
    subleading_phi = subleading.phi();
    //then flip phi
    subleading_phi = subleading_phi + pi;
    //then recalculate the cartesian components using the new phi
    subleading_px  = subleading_r*cos(subleading_t)*sin(subleading_phi);
    subleading_py  = subleading_r*sin(subleading_t)*sin(subleading_phi);
    subleading_pz  = subleading_r*cos(subleading_phi);
    //we can't redefine subleading, but we can make it equal to some other PseudoJet
    PseudoJet temporal1(subleading_px,subleading_py,subleading_pz,subleading.E());
    subleading = temporal1;

    //find their bisector, this is the average direction of leading and subleading
    axis = subleading.m()*leading + leading.m()*subleading;
    hist_axis_pT->Fill(axis.perp());
    
    int tracks = 0;
    //now we have to analyze particles individually
    for (int u = 0; u < event.size(); ++u) if (event[u].isFinal()) 
    {
      delta = sqrt(pow(event[u].eta() - axis.pseudorapidity(),2) + pow(event[u].phi() - axis.phi(),2)); 
      //This are the particles we expect to form part of the Jet
      if (event[u].isVisible() && event[u].isCharged() && deltaMin<=delta && delta<=deltaMax)
      {
        //Now we have to ptoject the particle into the axis
        //it is important that at the end we divide by the magnitude of the axis vector
        //if (axis.m()!=0)
        //{
        shadow = shadow + event[u].pT() * cos( event[u].phi() - axis.phi() );
        cout<<"The projection is "<<shadow<<endl;
        //}
        ++tracks;
      }
    }
    missingPT = missingPT + shadow/tracks;
    cout<<"The missing pT is "<<missingPT<<endl;
    hist_missing_pT->Fill(missingPT);

  // End of event loop.
  }

  cout<<"We finished recolecting data"<<endl;

  // Statistics. Histograms.
  pythia.stat();


  hist_leading_pT->Sumw2();
  hist_subleading_pT->Sumw2();
  hist_axis_pT->Sumw2();
  hist_missing_pT->Sumw2();
  hist_all_phi->Sumw2();
  hist_selected_phi->Sumw2();

  histos->cd();
  hist_leading_pT->Write();
  hist_subleading_pT->Write();
  hist_axis_pT->Write();
  hist_missing_pT->Write();
  hist_all_phi->Write();
  hist_selected_phi->Write();
  histos->Close();

  //clean up: delete all pointers
  delete hist_leading_pT;
  delete hist_subleading_pT;
  delete hist_axis_pT;
  delete hist_missing_pT;
  delete hist_all_phi;
  delete hist_selected_phi;

  return 0;
}
