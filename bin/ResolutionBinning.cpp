#include <iostream>
#include <cassert>
#include <TROOT.h>
#include <TSystem.h>
#include <Cintex/Cintex.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <vector>
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"


using namespace std;

//histograms
TH2F *lJetEtEtaBinning;
TH2F *bJetEtEtaBinning;
TH2F *muonEtEtaBinning;
TH2F *metEtEtaBinning;
TFile *outfile;

//variables
const static int nrEtbins = 10;
const static int nrEtabins = 10;
const static int maxEt = 200;
const static int maxEta = 200;
const static float truemaxEta = 2.5;
int totljets;
int totbjets;
int totm;
int totmet;

//
// create and define the histograms
//

void defineSolutionHistos() {
  lJetEtEtaBinning = new TH2F("lJetEtEtaBinning","to bin the Et-Eta of light jets",40,0,200,20,0,3);
  bJetEtEtaBinning = new TH2F("bJetEtEtaBinning","to bin the Et-Eta of b jets",40,0,200,20,0,3);
  muonEtEtaBinning = new TH2F("muonEtEtaBinning","to bin the Et-Eta of muons",40,0,200,20,0,3);
  metEtEtaBinning  = new TH2F("metEtEtaBinning","to bin the Et-Eta of the neutrino",40,0,200,20,0,3);
}

void defineVariables(){
  totljets = 0;
  totbjets = 0;
  totm = 0;
  totmet = 0;
}
	
//
// Analysis for 1 event
//
void AnalyseSolution(TtSemiEvtSolution &sol) {
 	 lJetEtEtaBinning  -> Fill(sol.getCalHadp().et(),fabs(sol.getCalHadp().eta()));
	 lJetEtEtaBinning  -> Fill(sol.getCalHadq().et(),fabs(sol.getCalHadq().eta()));
	 bJetEtEtaBinning  -> Fill(sol.getCalHadb().et(),fabs(sol.getCalHadb().eta()));
	 bJetEtEtaBinning  -> Fill(sol.getCalLepb().et(),fabs(sol.getCalLepb().eta()));
	 muonEtEtaBinning  -> Fill(sol.getRecLepm().et(),fabs(sol.getRecLepm().eta()));
	 metEtEtaBinning   -> Fill(sol.getRecLepn().et(),fabs(sol.getRecLepn().eta()));
}  


int main() {
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  outfile = new TFile("../data/EtEtaBinning.root", "RECREATE");;
  defineSolutionHistos();
  defineVariables();
   
  //lJet variables
  float EtlJetbinup[nrEtbins+1];
  EtlJetbinup[0]=0.;
  EtlJetbinup[nrEtbins]=(float)maxEt;
  int countEtlJet[maxEt+1];
  float EtalJetbinup[nrEtabins+1];
  EtalJetbinup[0]=0.;
  EtalJetbinup[nrEtabins]=truemaxEta;
  int countEtalJet[maxEta+1];
  float etljets[2],etaljets[2];
  
  //bJet variables
  float EtbJetbinup[nrEtbins+1];
  EtbJetbinup[0]=0.;
  EtbJetbinup[nrEtbins]=(float)maxEt;
  int countEtbJet[maxEt+1];
  float EtabJetbinup[nrEtabins+1];
  EtabJetbinup[0]=0.;
  EtabJetbinup[nrEtabins]=truemaxEta;
  int countEtabJet[maxEta+1];
  float etbjets[2],etabjets[2];
  
  //muon variables
  float Etmbinup[nrEtbins+1];
  Etmbinup[0]=0.;
  Etmbinup[nrEtbins]=(float)maxEt;
  int countEtm[maxEt+1];
  float Etambinup[nrEtabins+1];
  Etambinup[0]=0.;
  Etambinup[nrEtabins]=truemaxEta;
  int countEtam[maxEta+1];
  float etm,etam;
   
  //neutrino variables
  float Etmetbinup[nrEtbins+1];
  Etmetbinup[0]=0.;
  Etmetbinup[nrEtbins]=(float)maxEt;
  int countEtmet[maxEt+1];
  float etmet;
  
  //initialisation of the variables
  for(int i = 0; i <= maxEt; i++){
    countEtlJet[i] = 0;
    countEtbJet[i] = 0;
    countEtm[i] = 0;
    countEtmet[i] = 0;
  }
  for(int i = 0; i <= maxEta; i++){
    countEtalJet[i] = 0;
    countEtabJet[i] = 0;
    countEtam[i] = 0;
  }
  
  TString ft = "file:/beo5/pvmulder/CMSSW/src/myTopAnalyses/TopEtEtaBinning/test/LooserCutsEtEtaBinning.root"; 
  if (!gSystem->AccessPathName(ft)) {
    TFile *file = TFile::Open(ft);
    TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
    assert( events != 0 );
    TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_solutions__LooserCutsEtEtaBinning.obj" );
    assert( solsbranch != 0 );

    vector<TtSemiEvtSolution> sols;
    solsbranch->SetAddress( & sols );
      
    int nev = events->GetEntries();
    for( int ev = 0; ev < nev; ++ ev ) { 
      solsbranch -> GetEntry( ev );
      if(sols.size()>1){
        //get bestSol
        int bestSol = sols[0].getMCBestJetComb();
        if(bestSol>-1){ 
	  AnalyseSolution(sols[bestSol]);
	  
	  //lJet et and eta
	  etljets[0] = sols[bestSol].getCalHadp().et();
	  etljets[1] = sols[bestSol].getCalHadq().et();
	  etaljets[0] = fabs(sols[bestSol].getCalHadp().eta());
	  etaljets[1] = fabs(sols[bestSol].getCalHadq().eta());
	  //bJet et and eta
	  etbjets[0] = sols[bestSol].getCalHadb().et();
	  etbjets[1] = sols[bestSol].getCalLepb().et();
	  etabjets[0] = fabs(sols[bestSol].getCalHadb().eta());
	  etabjets[1] = fabs(sols[bestSol].getCalLepb().eta());
	  //muon et and eta
	  etm = sols[bestSol].getRecLepm().et();
	  etam = fabs(sols[bestSol].getRecLepm().eta());
	  //neutrino et
	  etmet = sols[bestSol].getRecLepn().et();
	 
	  for(int j = 0; j <= maxEt; j++){
	    for(int i=0; i<2; i++){
	      if((etljets[i] <= (float)j) && (etaljets[i] <= truemaxEta)) countEtlJet[j]++;
	      if((etbjets[i] <= (float)j) && (etabjets[i] <= truemaxEta)) countEtbJet[j]++;
	    }
	    if((etm <= (float)j) && (etam <= truemaxEta)) countEtm[j]++;
	    if((etmet <= (float)j)) countEtmet[j]++;
	  }     
 	  for(int j = 0; j <= maxEta; j++){
	    float trueEtaCut = (float)j*truemaxEta/maxEta;
	    for(int i=0; i<2; i++){
	      if((etaljets[i] <= trueEtaCut) && (etljets[i] <= (float)maxEt)) countEtalJet[j]++;
	      if((etabjets[i] <= trueEtaCut) && (etbjets[i] <= (float)maxEt)) countEtabJet[j]++;
	    }
	    if((etam <= trueEtaCut) && (etm <= (float)maxEt)) countEtam[j]++;
	  } 
	  //total number of ljets
	  if((etljets[0] <= (float)maxEt) && (etaljets[0] <= truemaxEta)) totljets++;
	  if((etljets[1] <= (float)maxEt) && (etaljets[1] <= truemaxEta)) totljets++;
	  //total number of bjets
	  if((etbjets[0] <= (float)maxEt) && (etabjets[0] <= truemaxEta)) totbjets++;
	  if((etbjets[1] <= (float)maxEt) && (etabjets[1] <= truemaxEta)) totbjets++;
	  //total number of muons
	  if((etm <= (float)maxEt) && (etam <= truemaxEta)) totm++;
	  //total number of neutrinos
	  if((etmet <= (float)maxEt)) totmet++;
	} 
      } 
    }
    file->Close();
  } 
  
  for(int i=1; i<nrEtbins; i++){
    EtlJetbinup[i]=0.;
    EtbJetbinup[i]=0.;
    Etmbinup[i]=0.;
    Etmetbinup[i]=0.;
    for(int Etcut = 0; Etcut <= maxEt; Etcut++){
      if(((float)(countEtlJet[Etcut])/((float)totljets) > ((float)(i)/(float)(nrEtbins))) && (EtlJetbinup[i]==0.)) EtlJetbinup[i]= (float)Etcut;
      if(((float)(countEtbJet[Etcut])/((float)totbjets) > ((float)(i)/(float)(nrEtbins))) && (EtbJetbinup[i]==0.)) EtbJetbinup[i]= (float)Etcut;
      if(((float)(countEtm[Etcut])/((float)totm) > ((float)(i)/(float)(nrEtbins))) && (Etmbinup[i]==0.)) Etmbinup[i]= (float)Etcut;
      if(((float)(countEtmet[Etcut])/((float)totmet) > ((float)(i)/(float)(nrEtbins))) && (Etmetbinup[i]==0.)) Etmetbinup[i]= (float)Etcut;
    }
  }
  for(int i=1; i<nrEtabins; i++){
    EtalJetbinup[i]=0.;
    EtabJetbinup[i]=0.;
    Etambinup[i]=0.;
    for(int j = 0; j <= maxEta; j++){
      if(((float)(countEtalJet[j])/((float)totljets) > ((float)(i)/(float)(nrEtabins))) && (EtalJetbinup[i]==0.)) EtalJetbinup[i]=(float)j*truemaxEta/maxEta;
      if(((float)(countEtabJet[j])/((float)totbjets) > ((float)(i)/(float)(nrEtabins))) && (EtabJetbinup[i]==0.)) EtabJetbinup[i]=(float)j*truemaxEta/maxEta;
      if(((float)(countEtam[j])/((float)totm) > ((float)(i)/(float)(nrEtabins))) && (Etambinup[i]==0.)) Etambinup[i]=(float)j*truemaxEta/maxEta;
    }
  } 
  cout << "in all the EtEtabins there is ~"<< 100/(float)(nrEtbins*nrEtabins) <<"% of the objects" << endl;
  cout << " " << endl;
  cout << "light jets | Etbins | Etabins"<< endl; 
  cout << "----------------------------" << endl; 
  for(int i=0; i<nrEtbins; i++){
    cout << "lJetsbin       " << EtlJetbinup[i] << "         " << EtalJetbinup[i] << endl;
  }
  cout << " " << endl;
  cout << "b jets     | Etbins | Etabins"<< endl; 
  cout << "----------------------------" << endl; 
  for(int i=0; i<nrEtbins; i++){
    cout << "bJetsbin       " << EtbJetbinup[i] << "         " << EtabJetbinup[i] << endl;
  } 
  cout << " " << endl;
  cout << "muons      | Etbins | Etabins"<< endl; 
  cout << "----------------------------" << endl; 
  for(int i=0; i<nrEtbins; i++){ 
    cout << "muonbin        " << Etmbinup[i]    << "         " << Etambinup[i]    << endl;
  }
  cout << " " << endl;
  cout << "neutrinos  | Etbins | Etabins"<< endl; 
  cout << "----------------------------" << endl; 
  for(int i=0; i<nrEtbins; i++){ 
    cout << "metbin         " << Etmetbinup[i]  << "    no eta-bin"    << endl;
  }
    
  gROOT->SetBatch();
  gROOT->SetStyle("Plain");  
  
  outfile->cd();
  outfile->Write();
  outfile->Close();

  return 0;
}
