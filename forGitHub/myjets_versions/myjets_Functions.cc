  //row row
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

  //Functions definitions
  //FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  //This function receives the particle data and returns jets sorted in decreasing pT

  vector<PseudoJet> getSortedJets(Event & event, double pTMin, int power, double R)
  {
      vector <PseudoJet> inclusiveJets, fjInputs;

      fjInputs.resize(0);

      for (int i = 0; i <= event.size(); ++i) if (event[i].isFinal()) 
      {
        if (!event[i].isVisible() || event[i].isNeutral() ) continue; //basically, only use charged events
        PseudoJet particleTemp = event[i];
        fjInputs.push_back(particleTemp);
      }

      JetDefinition jetDef(genkt_algorithm, R, power);
      ClusterSequence clustSeq(fjInputs, jetDef); //it takes what we found, fjInuts, and transforms it into jet that fit our jetDef
      inclusiveJets = clustSeq.inclusive_jets(pTMin);

      return sorted_by_pt(inclusiveJets); //return a vector of jets sorted into decreasing pT2
  }

  //This function holds the critaria for the leading Jet
  bool leadingCondition(PseudoJet leading, double pT_leading, double eta)
  {
    return ( leading.perp()>=pT_leading && abs(leading.eta())<=eta ) ;
  }

  //This function holds the criteria for the subleading Jet
  bool subleadingCondition(PseudoJet leading, PseudoJet subleading, double pT_subleading, double eta, double d_phi )
  {
    return ( subleading.perp()>=pT_subleading && abs(subleading.eta())<=eta && abs(leading.phi()-subleading.phi())>d_phi ) ;
  }

  //This function returns the index of the leading Jet
  int findMax(vector<PseudoJet> & sortedJets, double pT_leading, double eta)
  {
    PseudoJet leading;
    for (int j = 0; j <= sortedJets.size(); ++j)
    {
      leading  = sortedJets[j];
      if ( leadingCondition( leading, pT_leading, eta) )
      {
        return j;
      }
    }
    int j=sortedJets.size();
    return j;
  }

  //This function returns the index of the subleading Jet
  int findNextMax( vector<PseudoJet> & sortedJets, PseudoJet & leading, int j, double pT_subleading, double eta, double d_phi )
  {
    PseudoJet subleading;
    for (int k = j+1; k <= sortedJets.size(); ++k)
    {
      subleading  = sortedJets[k];
      if ( subleadingCondition( leading, subleading, pT_subleading, eta, d_phi ) )
      {
        return k;
      }
    }
    int k=0;
    return k;
  }
  //FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF


  int main() 
  {
    TFile* histos             = new TFile("histos.root","recreate");
    TH1F*  hist_leading_pT    = new TH1F("hist_leading_pT", "Momentum of leading Jet", 1000, 0, 500);
    TH1F*  hist_subleading_pT = new TH1F("hist_subleading_pT", "Momentum of subleading Jet", 1000, 0, 500);
    TH1F*  hist_axis_pT       = new TH1F("hist_axis_pT", "Momentum of axis", 1000, 0, 500);

    double delta_edges[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 3.6 };
    double pT_edges[]    = { 0.5, 1, 2, 4, 8, 300 };
    TH2D*  hist_pTdelta  = new TH2D("hist_pTdelta", "The momentum of particles in each Delta Ring", 5, pT_edges, 10, delta_edges );

    
    // All the things you can change
    int    power          = -1; // -1 = anti-kT; 0 = C/A; 1 = kT.
    double R              = 0.4;    
    double pTMin          = 0.0;
    double eta            = 2.0;
    double delta          = 0;
    double deltaMin       = 0.2;
    double deltaMax       = 0.4;
    double pi             = 3.14159265359;
    double d_phi          = (5*pi)/6; 
    double jetshadow      = 0;
    double nonjetshadow   = 0;
    double pT_leading     = 120;
    double pT_subleading  = 50;
    double subleading_r   = 0;
    double subleading_t   = 0;
    double subleading_phi = 0;
    double subleading_px  = 0;
    double subleading_py  = 0;
    double subleading_pz  = 0;

    // Generator. Shorthand for event.
    Pythia pythia;
    pythia.readFile("seetings.cmnd");
    int Nev = pythia.mode("Main:numberOfEvents");
    Event& event = pythia.event;
    pythia.init();

    JetDefinition jetDef(genkt_algorithm, R, power);

    // Begin event loop. Generate event. Skip if error.
    for (int iEvent = 0; iEvent < Nev; ++iEvent) 
    {
      if (!pythia.next()) continue;

      //Reset variables
      jetshadow      = 0;
      nonjetshadow   = 0;
      delta          = 0;
      PseudoJet leading, subleading, axis, particle;
      vector <PseudoJet> sortedJets;

      sortedJets = getSortedJets(event,pTMin,power,R);

      if (sortedJets.size()<2) continue; //there were two few Jets in this event

      int j = findMax( sortedJets, pT_leading, eta);

      if (j==sortedJets.size()) continue; //there are no suitable Jets in this event

      leading = sortedJets[j];

      hist_leading_pT->Fill(leading.perp());

      int k = findNextMax( sortedJets, leading, j, pT_subleading, eta, d_phi ); 

      if (k==0) continue; // this was not a dijet event

      subleading = sortedJets[k]; 

      hist_subleading_pT->Fill(subleading.perp());

  //####################################################################################################

      //Now, we have to reflect subleading around the origin, which requieres a bit of math
      //first the spherical components
      
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

      //find their bisector, this is the average direction of leading and subleading
      axis = subleading.m()*leading + leading.m()*subleading;
      
      //now we have to analyze particles individually
      
      for (int n = 0; n <= event.size(); ++n)
      {
        PseudoJet particle = event[n];
        delta = sqrt(pow(particle.eta()-axis.eta(),2)+pow(particle.phi()-axis.phi(),2));
        if (delta>3.6 || particle.perp()>300) continue;
        pTll = -particle.perp()*cos( particle.phi()-axis.phi() );
        hist_pTdelta->Fill(pTll,delta); //we fill the histogram in the same order that we declared it, first pT and then delta
      }

      TProfile* hist_test = TH2::ProfileY(hist_pTdelta);

      for (int y = 0; y < 9; ++y)
      {
        for (int x = 0; x <4 ; ++x)
        {
          bin = hist_pTdelta->GetBin(x,y);
        }
      }
    // End of event loop.
    }

    // Statistics. Histograms.
    pythia.stat();


    hist_leading_pT->Sumw2();
    hist_subleading_pT->Sumw2();
    hist_axis_pT->Sumw2();
    hist_pTdelta->Sumw2();

    histos->cd();
    hist_leading_pT->Write();
    hist_subleading_pT->Write();
    hist_axis_pT->Write();
    hist_pTdelta->Write();
    histos->Close();

    return 0;
  }
