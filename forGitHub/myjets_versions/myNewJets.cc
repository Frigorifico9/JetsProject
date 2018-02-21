//the town
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
//CLASES
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//This class holds the cartesian coordinates of the PseudoJet
class Vector3
{
public:
  double x, y, z;
  PseudoJet & pseudoJet;

  Vector3( PseudoJet* pseudoJet)
  {
    this->pseudoJet = *pseudoJet;
    this->x = pseudoJet.px();
    this->y = pseudoJet.py();
    this->z = pseudoJet.pz();
  }
 
//This method lets you access the PesudoJet
  PseudoJet getPseudoJet() 
  {
    return this->pseudoJet;
  }
};

class VectorS
{
public:
  double r, theta, phi;
  PseudoJet & pseudoJet;

  VectorS( PseudoJet* pseudoJet )
  {
    this->pseudoJet = *pseudoJet;
    Vector3 temporal(pseudoJet)
    double x  = temporal.x();
    double y  = temporal.y();
    double z  = temporal.z();
    this->r = sqrt( pow( x, 2) + pow( y, 2) + pow( z, 2) );
    this->theta = atan(y/x);
    this->phi = this->pseudoJet.phi();
  }
//Method to shift phi by an angle a 
  void shiftPhi( double a = 3.14159265359 )
  {
    this->phi = this->pseudoJet.phi() + a; //we shift phi
    double x  = this->r*cos(this->theta)*sin(this->phi);
    double y  = this->r*sin(this->theta)*sin(this->phi);
    double z  = this->r*cos(this->phi);
    this->pseudoJet = (new PseudoJet( x, y, z, pseudoJet.E())); //we create the PseudoJet with the shifted phi
  }

//This method lets you access the PesudoJet
  PseudoJet getPseudoJet() 
  {
    return this->pseudoJet;
  }
};

class VectorX
{
public:
  double r, theta, phi, x, y, z;
  PseudoJet & pseudoJet;

//constructor
  VectorX( PseudoJet* pseudoJet )
  {
    this->pseudoJet = *pseudoJet;
    this->x = pseudoJet.px();
    this->y = pseudoJet.py();
    this->z = pseudoJet.pz();
  
    this->r = sqrt( pow( this->x, 2) + pow( this->y, 2) + pow( this->z, 2) );
    this->theta = atan(this->y/this->x);
    this->phi = this->pseudoJet.phi();
  }

//Method to go from spherical coordinates to cartesian coordinates
  Vector3 toVector3()
  {
    double x  = this->r*cos(this->theta)*sin(this->phi);
    double y  = this->r*sin(this->theta)*sin(this->phi);
    double z  = this->r*cos(this->phi);
    return Vector3(x,y,z,this->pseudoJet);
  }

//Method to go from cartesian coordinates to cartesian spherical
  VectorS toVectorS()
  {
    double r = sqrt( pow( this->x, 2) + pow( this->y, 2) + pow( this->z, 2) );
    double theta = atan(this->y/this->x);
    double phi = this->pseudoJet.phi();
    return VectorS(r,theta,phi,this->pseudoJet);
  }

};

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//FUNCTIONS
//FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
//This function receives the particle data and returns jets sorted in decreasing pT
vector<PseudoJet> getSortedJets(Event & event, double pTMin)
{
    vector <PseudoJet> inclusiveJets, fjInputs;

    fjInputs.resize(0);

    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) 
    {
      if (!event[i].isVisible() || event[i].isNeutral() ) continue; //basically, only use charged events
      PseudoJet particleTemp = event[i];
      fjInputs.push_back(particleTemp);
    }
    //cout<<"The particles' data was stored"<<endl;

    ClusterSequence clustSeq(fjInputs, jetDef); //it takes what we found, fjInuts, and transforms it into jet that fit our jetDef
    inclusiveJets = clustSeq.inclusive_jets(pTMin);

    return sorted_by_pt(inclusiveJets); //return a vector of jets sorted into decreasing pT2
}

//This function holds the critaria for the leading Jet
bool leadingCondition(PseudoJet leading, double pT_leading, double eta)
{
  return ( leading.perp()>=pT_leading && abs(leading.pseudorapidity())<=eta ) ;
}

//This function holds the criteria for the subleading Jet
bool subleadingCondition(PseudoJet leading, PseudoJet subleading, double pT_subleading, double eta, double d_phi )
{
  return ( subleading.perp()>=pT_subleading && abs(subleading.pseudorapidity())<=eta && abs(leading.phi()-subleading.phi())>d_phi ) ;
}

//This function returns the index of the leading Jet
int findMax(vector<PseudoJet> & sortedJets, double pT_leading, double eta)
{
  PseudoJet leading;
  for (int j = 0; j < sortedJets.size(); ++j)
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
  for (int k = j+1; k < sortedJets.size(); ++k)
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

//This function projects event onto axis using the 1 -1 -1 -1 metric 
double getShadow(int u, Event & event, PseudoJet axis)
{
  PseudoJet projected = event[u];
  return (axis.E()*projected.E()-axis.px()*projected.px()-axis.py()*projected.py()-axis.pz()*projected.pz())/axis.m();
}
//FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

int main() 
{
  TFile* histos             = new TFile("histos.root","recreate");
  TH1F*  hist_leading_pT    = new TH1F("hist_leading_pT", "Momentum of leading Jet", 1000, 0, 1000000);
  TH1F*  hist_subleading_pT = new TH1F("hist_subleading_pT", "Momentum of subleading Jet", 1000, 0, 1000000);
  TH1F*  hist_axis_pT       = new TH1F("hist_axis_pT", "Momentum of axis", 1000, 0, 1000000);
  TH1F*  hist_jetTrack_pT   = new TH1F("hist_jetTrack_pT", "Momentum of tracks inside jets", 1000, 0, 1000000);
  TH1F*  hist_track_pT      = new TH1F("hist_track_pT", "Momentum of tracks outside jets", 1000, 0, 1000000);
  
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
    PseudoJet leading, subleading, axis;
    vector <PseudoJet> sortedJets;

    sortedJets = getSortedJets(event,pTMin);

    if (sortedJets.size()<2) continue; //there were two few Jets in this event

    int j = findMax( sortedJets, pT_leading, eta);

    if (j==sortedJets.size()) continue; //there are no suitable Jets in this event

    leading = sortedJets[j];

    int k = findNextMax( sortedJets, leading, j, pT_subleading, eta, d_phi ); 

    if (k==0) continue; // this was not a dijet event

    subleading = sortedJets[k]; 

//####################################################################################################

    //Now, we have to reflect subleading around the origin, which requieres a bit of math
    //first the spherical components
    VectorX x_subleading(subleading);
    VectorS s_subleading = x_subleading.toSpherical();
    s_subleading.shiftPhi(pi);
    subleading = s_coordinates.getPseudoJet();

    //find their bisector, this is the average direction of leading and subleading
    axis = subleading.m()*leading + leading.m()*subleading;
    
    //now we have to analyze particles individually
    for (int u = 0; u < event.size(); ++u) if (event[u].isFinal()) 
    {
      delta = sqrt(pow(event[u].eta() - axis.pseudorapidity(),2) + pow(event[u].phi() - axis.phi(),2)); 
      //This are the particles we expect to form part of the Jet
      if (event[u].isVisible() && event[u].isCharged() && deltaMin<=delta && delta<=deltaMax)
      {
        //Now we have to project the particle into the axis
        //it is important that at the end we divide by the magnitude of the axis vector
        jetshadow = jetshadow + getShadow(u, event, axis);
      }
      //This are the particles we don't expect to form part of the Jet
      else
      {
        nonjetshadow = nonjetshadow + getShadow(u, event, axis);
      }
    }

  // End of event loop.
  }

  // Statistics. Histograms.
  pythia.stat();


  hist_leading_pT->Sumw2();
  hist_subleading_pT->Sumw2();
  hist_axis_pT->Sumw2();
  hist_jetTrack_pT->Sumw2();
  hist_track_pT->Sumw2();

  histos->cd();
  hist_leading_pT->Write();
  hist_subleading_pT->Write();
  hist_axis_pT->Write();
  hist_jetTrack_pT->Write();
  hist_track_pT->Write();
  histos->Close();

  return 0;
}
