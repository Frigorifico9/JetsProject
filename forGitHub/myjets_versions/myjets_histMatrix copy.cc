//philip
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
#include "TObjArray.h"
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

    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) 
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
  TH2D*  hist_pTdelta_AjInc = new TH2D("hist_pTdelta_AjInc", "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20); //Aj inclusive
  TH2D*  hist_pTdelta_AjA   = new TH2D("hist_pTdelta_AjA", "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20); //Aj above Aj_V 
  TH2D*  hist_pTdelta_AjB   = new TH2D("hist_pTdelta_AjB", "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20); //Aj bellow Aj_V

/*/
  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      all_histos[f][c] = new TH2D( Form("h%d%d",f,c), "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20);
    }
  }
/*/

  //All the things is best not to change
  double subleading_r       = 0;
  double subleading_t       = 0;
  double subleading_phi     = 0;
  double subleading_px      = 0;
  double subleading_py      = 0;
  double subleading_pz      = 0;
  double delta              = 0;
  double pi                 = 3.14159265359;
  double pTll               = 0.0;
  double Aj                 = 0;
  int    mult_range[]       = {234, 204, 188, 168, 152, 140, 128, 116, 104, 92, 76, 0}; //multiplicity ranges
  const  Char_t* Aj_range[] = {"Inc", "A", "B"};
  int    mult               = 0;
  int    filas              = 12;
  int    columnas           = 3;
  TH2D*  all_histos[filas][columnas];

  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      //all_histos[f][c] = new TH2D( Form("hist_mult>=%d_Aj%s",mult_range[f],Aj_range[c]), "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20);
      all_histos[f][c] = new TH2D( Form("hist_mult>=%d_Aj%s",mult_range[f],Aj_range[c]), "The momentum of particles in each Delta Ring", 9, 0, 1.8, 1000, -20, 20);
    }
  }
 
  // All the things you can change  
  int    power          = -1; // -1 = anti-kT; 0 = C/A; 1 = kT.
  double R              = 0.4; //radius of the jet cone   
  double pTMin          = 0.0; //minimum momentum of the jets
  double eta            = 2.0; //maximum eta of jets
  double d_phi          = (5*pi)/6; //minimum angular separation of leading and subleading jets 
  double pT_leading     = 120; //minimum pomentum of leadin jet
  double pT_subleading  = 50; //minimum momentum of subleading jet
  double Aj_V           = 0.22; //Aj value we are using to discriminate events

  // Generator. Shorthand for event.
  Pythia pythia;
  pythia.readFile("seetings.cmnd");
  int Nev = pythia.mode("Main:numberOfEvents");
  Event& event = pythia.event;
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < Nev; ++iEvent) 
  {
    if (!pythia.next()) continue;

    hist_pTdelta_AjInc->Reset();
    hist_pTdelta_AjA->Reset();
    hist_pTdelta_AjB->Reset();

    //hist_pTdelta_AjInc->Sumw2(false);
    //hist_pTdelta_AjA->Reset(false);
    //hist_pTdelta_AjB->Reset(false);    

    //Reset variables just because I'm a paranoid
    delta          = 0;
    pTll           = 0;
    Aj             = 0;
    subleading_r   = 0;
    subleading_t   = 0;
    subleading_phi = 0;
    subleading_px  = 0;
    subleading_py  = 0;
    subleading_pz  = 0;
    delta          = 0;
    pTll           = 0;
    PseudoJet leading, subleading, axis, particle;
    vector <PseudoJet> sortedJets;

    sortedJets = getSortedJets(event,pTMin,power,R);

    if (sortedJets.size()<2) continue; //there were too few Jets in this event

    int j = findMax( sortedJets, pT_leading, eta);

    if (j==sortedJets.size()) continue; //there are no suitable Jets in this event

    leading = sortedJets[j];

    int k = findNextMax( sortedJets, leading, j, pT_subleading, eta, d_phi ); 

    if (k==0) continue; // this was not a dijet event

    subleading = sortedJets[k];

    Aj = (leading.perp()-subleading.perp())/(leading.perp()+subleading.perp());

//####################################################################################################

    //Now, we have to reflect subleading around the origin, which requieres a bit of math
    //first the spherical components
    
    subleading_r   = sqrt(pow(subleading.px(),2)+pow(subleading.py(),2)+pow(subleading.pz(),2));
    subleading_t   = acos(subleading.pz()/subleading_r);
    subleading_phi = subleading.phi();
    //then flip phi
    subleading_phi = subleading_phi + pi;
    //then recalculate the cartesian components using the new phi
    subleading_px  = subleading_r*cos(subleading_phi)*sin(subleading_t);
    subleading_py  = subleading_r*sin(subleading_phi)*sin(subleading_t);
    subleading_pz  = subleading_r*cos(subleading_t);
    //we can't redefine subleading, but we can make it equal to some other PseudoJet
    PseudoJet temporal1(subleading_px,subleading_py,subleading_pz,subleading.E());
    subleading = temporal1;
    //find their bisector, this is the average direction of leading and subleading
    axis = subleading.m()*leading + leading.m()*subleading;

    //now we analize individual events using the axis

    mult = 0;
    for (int n = 0; n <= event.size(); ++n) if (event[n].isCharged() && event[n].isFinal())
    {
      ++mult;
      particle = event[n]; //we declare a particle
      //delta = sqrt(pow(particle.eta()-axis.eta(),2)+pow(particle.phi()-axis.phi(),2)); //we get the delta of that particle with the aixs
      delta = sqrt(pow(event[n].eta()-axis.eta(),2)+pow(event[n].phi()-axis.phi(),2));
      //if (delta>3.6) continue;
      if ( !( abs(particle.perp())<300 && (delta<1.8) ) && !( abs(particle.perp())>300 && (3.6>delta>1.8) ) ) continue; //if there are charged particles outside our pT range but have a delta between 3.6 and 1.8 they go to the overflow bin
      pTll = -particle.perp()*cos( particle.phi()-axis.phi() ); //we calculate the "projection"
      //pTll = -event[n].pT()*cos( event[n].phi()-axis.phi() );
      hist_pTdelta_AjInc->Fill(delta,pTll); //we fill the histogram in the same order that we declared it, first pT and then delta
      if (Aj>Aj_V)  hist_pTdelta_AjA->Fill(delta,pTll);
      if (Aj<Aj_V)  hist_pTdelta_AjB->Fill(delta,pTll);
    }

    //hist_pTdelta_AjInc->Sumw2();
    //hist_pTdelta_AjA->Sumw2();
    //hist_pTdelta_AjB->Sumw2();

    //Now that we have everything, lets put it in the right collection
    
    for (int n = 0; n < 12; ++n)
    {
      if (mult>=mult_range[n]) 
      {
        all_histos[n][0]->Add(hist_pTdelta_AjInc,1);
        all_histos[n][1]->Add(hist_pTdelta_AjA,1);
        all_histos[n][2]->Add(hist_pTdelta_AjB,1);
        break;
      }
    }
  // End of event loop.
  }

  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      all_histos[f][c]->Sumw2();
    }
  }

  // Statistics. Histograms.
  pythia.stat();

  histos->cd();
  for (int f = 0; f < filas; ++f)
  {
    for (int c = 0; c < columnas; ++c)
    {
      all_histos[f][c] ->Write();
    }
  }
  histos->Close();

  return 0;
}