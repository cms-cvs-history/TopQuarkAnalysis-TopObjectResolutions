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
#include "FWCore/FWLite/src/AutoLibraryLoader.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "AnalysisDataFormats/TopObjects/interface/TtSemiEvtSolution.h"


using namespace std;

//histograms
TH2F *lJetEtEtaBinning;
TH2F *bJetEtEtaBinning;
TH2F *muonEtEtaBinning;
TH2F *metEtEtaBinning;
TH1F *lJetsBinning;//to check if code is working
TH1F *bJetsBinning;//to check if code is working
TH1F *mBinning;//to check if code is working
TH1F *metBinning;//to check if code is working
TFile *outfile;

//variables
const static int nrEtbins = 10;
const static int nrEtabins = 10;
const static int maxEt = 200;
const static int maxEta = 200;
const static float truemaxEta = 2.5;
const static int nrsel = 100000; // fix number, otherwise array possibly to big !! 
int icount;
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
  lJetsBinning 	   = new TH1F("lJetsBinning","same number of lJets in each bin ",50,0,0.05);//to check if code is working
  bJetsBinning 	   = new TH1F("bJetsBinning","same number of bJets in each bin ",50,0,0.05);//to check if code is working
  mBinning 	   = new TH1F("mBinning","same number of muons in each bin ",50,0,0.05);//to check if code is working
  metBinning 	   = new TH1F("metBinning","same number of neutrinos in each bin ",50,0,0.5);//to check if code is working
}

void defineVariables(){
  icount = 0;
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
  
  cout << "I am defining the variables" << endl;
  
  //lJet variables
  float EtlJetbinup[nrEtbins+1];
  EtlJetbinup[0]=0.;
  EtlJetbinup[nrEtbins]=(float)maxEt;
  int countEtlJet[maxEt+1];
  float EtalJetbinup[nrEtabins+1];
  EtalJetbinup[0]=0.;
  EtalJetbinup[nrEtabins]=truemaxEta;
  int countEtalJet[maxEta+1];
  float etljets[2][nrsel],etaljets[2][nrsel];
  int countEtEtalJet[nrEtbins][nrEtabins]; //to check if code is working
  
  //bJet variables
  float EtbJetbinup[nrEtbins+1];
  EtbJetbinup[0]=0.;
  EtbJetbinup[nrEtbins]=(float)maxEt;
  int countEtbJet[maxEt+1];
  float EtabJetbinup[nrEtabins+1];
  EtabJetbinup[0]=0.;
  EtabJetbinup[nrEtabins]=truemaxEta;
  int countEtabJet[maxEta+1];
  float etbjets[2][nrsel],etabjets[2][nrsel];
  int countEtEtabJet[nrEtbins][nrEtabins];//to check if code is working
  
  //muon variables
  float Etmbinup[nrEtbins+1];
  Etmbinup[0]=0.;
  Etmbinup[nrEtbins]=(float)maxEt;
  int countEtm[maxEt+1];
  float Etambinup[nrEtabins+1];
  Etambinup[0]=0.;
  Etambinup[nrEtabins]=truemaxEta;
  int countEtam[maxEta+1];
  float etm[nrsel],etam[nrsel];
  int countEtEtam[nrEtbins][nrEtabins];//to check if code is working
  
  //neutrino variables
  float Etmetbinup[nrEtbins+1];
  Etmetbinup[0]=0.;
  Etmetbinup[nrEtbins]=(float)maxEt;
  int countEtmet[maxEt+1];
  float etmet[nrsel];
  int countEtmettest[nrEtbins];//to check if code is working
  
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
    cout << "I am opening the rootfile" << endl;
    TTree * events = dynamic_cast<TTree*>( file->Get( "Events" ) );
    assert( events != 0 );
    TBranch * solsbranch = events->GetBranch( "TtSemiEvtSolutions_solutions__LooserCutsEtEtaBinning.obj" );
    assert( solsbranch != 0 );

    vector<TtSemiEvtSolution> sols;
    solsbranch->SetAddress( & sols );
    
    cout << "now I will choose the best solution for each event and calculate the nr of objects for certain Et and Eta" << endl; 
    
    int nev = events->GetEntries();
    for( int ev = 0; ev < nev; ++ ev ) { 
      solsbranch -> GetEntry( ev );
      if(sols.size()>1){
        //get bestSol
        int bestSol = sols[0].getMCBestJetComb();
        if(bestSol>-1){ 
	  AnalyseSolution(sols[bestSol]);
	  icount++;
	  
	  if(icount < nrsel){
	    //lJet et and eta
	    etljets[0][icount] = sols[bestSol].getCalHadp().et();
	    etljets[1][icount] = sols[bestSol].getCalHadq().et();
	    etaljets[0][icount] = fabs(sols[bestSol].getCalHadp().eta());
	    etaljets[1][icount] = fabs(sols[bestSol].getCalHadq().eta());
	    //bJet et and eta
	    etbjets[0][icount] = sols[bestSol].getCalHadb().et();
	    etbjets[1][icount] = sols[bestSol].getCalLepb().et();
	    etabjets[0][icount] = fabs(sols[bestSol].getCalHadb().eta());
	    etabjets[1][icount] = fabs(sols[bestSol].getCalLepb().eta());
	    //muon et and eta
	    etm[icount] = sols[bestSol].getRecLepm().et();
	    etam[icount] = fabs(sols[bestSol].getRecLepm().eta());
	    //neutrino et
	    etmet[icount] = sols[bestSol].getRecLepn().et();
	  
	    for(int j = 0; j <= maxEt; j++){
	      for(int i=0; i<2; i++){
	        if((etljets[i][icount] <= (float)j) && (etaljets[i][icount] <= truemaxEta)) countEtlJet[j]++;
	        if((etbjets[i][icount] <= (float)j) && (etabjets[i][icount] <= truemaxEta)) countEtbJet[j]++;
	      }
	      if((etm[icount] <= (float)j) && (etam[icount] <= truemaxEta)) countEtm[j]++;
	      if((etmet[icount] <= (float)j)) countEtmet[j]++;
	    }     
 	    for(int j = 0; j <= maxEta; j++){
	      float trueEtaCut = (float)j*truemaxEta/maxEta;
	      for(int i=0; i<2; i++){
	        if((etaljets[i][icount] <= trueEtaCut) && (etljets[i][icount] <= (float)maxEt)) countEtalJet[j]++;
	        if((etabjets[i][icount] <= trueEtaCut) && (etbjets[i][icount] <= (float)maxEt)) countEtabJet[j]++;
	      }
	      if((etam[icount] <= trueEtaCut) && (etm[icount] <= (float)maxEt)) countEtam[j]++;
	    } 
	    //total number of ljets
	    if((etljets[0][icount] <= (float)maxEt) && (etaljets[0][icount] <= truemaxEta)) totljets++;
	    if((etljets[1][icount] <= (float)maxEt) && (etaljets[1][icount] <= truemaxEta)) totljets++;
	    //total number of bjets
	    if((etbjets[0][icount] <= (float)maxEt) && (etabjets[0][icount] <= truemaxEta)) totbjets++;
	    if((etbjets[1][icount] <= (float)maxEt) && (etabjets[1][icount] <= truemaxEta)) totbjets++;
	    //total number of muons
	    if((etm[icount] <= (float)maxEt) && (etam[icount] <= truemaxEta)) totm++;
	    //total number of neutrinos
	    if((etmet[icount] <= (float)maxEt)) totmet++;
	  }
	} 
      } 
    }
    file->Close();
  } 
  
  if(icount > nrsel) icount = nrsel;
  
  cout << "now I will calculate the EtEtabins" << endl; 
    
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
  
  cout << "now I will make the histograms which show that there is a peak around  ~"<< 1/(float)(nrEtbins*nrEtabins) << endl;
  
  //to check if code is working: lines 256 -> 287
  for(int iEtbin=1; iEtbin < nrEtbins; iEtbin++){
    for(int iEtabin=1; iEtabin < nrEtabins; iEtabin++){
      countEtEtalJet[iEtbin][iEtabin]=0;
      countEtEtabJet[iEtbin][iEtabin]=0;
      countEtEtam[iEtbin][iEtabin]=0;
    }
    countEtmettest[iEtbin]=0;
  }
	
  for(int iEtbin=1; iEtbin <= nrEtbins; iEtbin++){
    for(int iEtabin=1; iEtabin <= nrEtabins; iEtabin++){
      for(int i=1; i <= icount; i++){
	for(int ijet=0; ijet<2; ijet++){
	  if((etljets[ijet][i] > EtlJetbinup[iEtbin-1]) && (etljets[ijet][i] < EtlJetbinup[iEtbin]) && (etaljets[ijet][i] > EtalJetbinup[iEtabin-1]) && (etaljets[ijet][i] <
	    EtalJetbinup[iEtabin])) countEtEtalJet[iEtbin][iEtabin]++;  
	  if((etbjets[ijet][i] > EtbJetbinup[iEtbin-1]) && (etbjets[ijet][i] < EtbJetbinup[iEtbin]) && (etabjets[ijet][i] > EtabJetbinup[iEtabin-1]) && (etabjets[ijet][i] <
	    EtabJetbinup[iEtabin])) countEtEtabJet[iEtbin][iEtabin]++;
	}
	if((etm[i] > Etmbinup[iEtbin-1]) && (etm[i] < Etmbinup[iEtbin]) && (etam[i] > Etambinup[iEtabin-1]) && (etam[i] < Etambinup[iEtabin])) countEtEtam[iEtbin][iEtabin]++;
      }
      lJetsBinning->Fill((float)countEtEtalJet[iEtbin][iEtabin]/(float)totljets);
      bJetsBinning->Fill((float)countEtEtabJet[iEtbin][iEtabin]/(float)totbjets);   
      mBinning    ->Fill((float)countEtEtam[iEtbin][iEtabin]/(float)totm);   
    }
  }
  for(int iEtbin=1; iEtbin <= nrEtbins; iEtbin++){
    for(int i=1; i <= icount; i++){
      if((etmet[i] > Etmetbinup[iEtbin-1]) && (etmet[i] < Etmetbinup[iEtbin])) countEtmettest[iEtbin]++;
    }
    metBinning    ->Fill((float)countEtmettest[iEtbin]/(float)totmet); 
  } 
    
  gROOT->SetBatch();
  gROOT->SetStyle("Plain");  
  
  outfile->cd();
  outfile->Write();
  outfile->Close();

  return 0;
}
