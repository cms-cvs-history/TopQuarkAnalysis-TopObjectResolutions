#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMultiLayerPerceptron.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TKey.h"
#include "TObjString.h"
#include <iostream>
#include <vector>

/////////////////////////////////////////////////////////////////
// MACRO TO DRAW AND EVALUATE THE RESOLUTION PLOTS             //
/////////////////////////////////////////////////////////////////

void ResolutionEvaluation(TString fileName="Resolutions_tau.root", bool doNN=true, TString output="ResolutionEvaluation.ps") {
  
  // define all histo's
  int Nepoch = 10000;
  TH1F  *hResEtEtaBin[10][20][20], *hResEtaBin[10][20];
  TString  resObsName[8] = {"_ares","_bres","_cres","_dres","_thres","_phres","_etres","_etares"};
  TString objectType = fileName;
  objectType.Remove(0,objectType.Index("_")+1);
  objectType.Remove(objectType.Index("."), objectType.Length());
  if(objectType.Index("_")>0) objectType.Remove(objectType.Index("_"),objectType.Length());

  // read in histo's
  TFile *fIn = new TFile(fileName,"UPDATE");
  int nrETBins = 0; 
  std::vector<Double_t> etabin;
  TH1F *tmpEta = (TH1F*) (fIn->GetKey("hEtaBins")->ReadObj());
  for(Int_t b=1; b<=tmpEta->GetNbinsX(); b++) etabin.push_back(tmpEta->GetXaxis()->GetBinLowEdge(b));
  etabin.push_back(tmpEta->GetXaxis()->GetBinUpEdge(tmpEta->GetNbinsX()));
  int nrEtaBins = etabin.size()-1;
  std::cout<<"Found "<<nrEtaBins<< " eta-bins with edges: ( ";
  for(size_t u=0; u<etabin.size(); u++) cout<<etabin[u]<<", ";
  std::cout<<"\b\b )"<<std::endl;
  
  for(int ro=0; ro<8; ro++) {
    for(int etab=0; etab<nrEtaBins; etab++) {	
      TString resEtabinName = objectType; resEtabinName += resObsName[ro]; resEtabinName += "_etabin"; resEtabinName += etab;
      hResEtaBin[ro][etab] = (TH1F*) (fIn->GetKey(resEtabinName)->ReadObj());
      nrETBins = (Int_t) hResEtaBin[ro][etab]->GetEntries();
      for(int etb=0; etb<nrETBins; etb++) {
        TString eTEtabinName = resEtabinName; eTEtabinName += "_etbin";  eTEtabinName+= etb;
        hResEtEtaBin[ro][etab][etb] = (TH1F*) (fIn->GetKey(eTEtabinName)->ReadObj());
      }
    }
  }
  
    
  // Draw resolutions in (ET,Eta)-bins & resolutions vs ET for different Eta bins
  TCanvas *c = new TCanvas("dummy","",1);
  c -> Print(output+"[","landscape");
  TCanvas *cEtEta[10], *cEta[10];
  for(Int_t ro=0; ro<6; ro++) { //loop on variables
    TString cName   = "Res_vs_ET_for_different_Eta_bins_for_ResoObs_"; cName += resObsName[ro];
    TString cName2  = "EtEtaBinResolutions_for_ResoObs_"; cName2 += resObsName[ro];
    cEtEta[ro] =  new TCanvas(cName,cName,1);
    cEtEta[ro] -> Divide(nrEtaBins,nrETBins);
    cEta[ro]   =  new TCanvas(cName2,cName2,1);
    cEta[ro]   -> Divide(int(sqrt(nrEtaBins)),int(ceil(nrEtaBins/float(int(sqrt(nrEtaBins))))));
    // do the NN part
    TMultiLayerPerceptron* mlp = NULL;
    TTree* tResVar = NULL;
    if(doNN){
      tResVar = (TTree*)fIn->GetKey("tResVar")->ReadObj();
      TH1D tmp("tmpb","tmpb",1,-FLT_MAX,FLT_MAX);
      tResVar->Draw("value>>tmpb",Form("ro==%d",ro),"goff");
      if(nrEtaBins>1) {
        mlp = new TMultiLayerPerceptron("@Et,@Eta:6:3:value","1./error",tResVar,
        //mlp = new TMultiLayerPerceptron("@Et,@Eta:6:3:value","1.",tResVar,
  	                                Form("ro==%d",ro),
	                                Form("ro==%d",ro));
	Nepoch = 10000;
      }
      else {
        mlp = new TMultiLayerPerceptron("@Et:2:value","1./error",tResVar,
  	                                Form("ro==%d",ro),
	                                Form("ro==%d",ro));
	Nepoch = 1000;
      }
      mlp->Train(Nepoch,"graph update=10");
      gPad->Print(output,"landscape");
    }
    //find the axis limits
    tResVar->Draw("value:Et:Eta>>new3d","ro==0","goff");
    TProfile2D* prof2D = ((TH3F*)gDirectory->GetObjectUnchecked("new3d"))->Project3DProfile("xy");
    Double_t maxX = prof2D->GetXaxis()->GetXmax();
    Double_t maxY = prof2D->GetYaxis()->GetXmax();
    Double_t minX = prof2D->GetXaxis()->GetXmin();
    Double_t minY = prof2D->GetYaxis()->GetXmin();
    Double_t maxZ = 0.;
    Double_t minZ = FLT_MAX;
    // loop on eta bins
    for(Int_t etab=0; etab<nrEtaBins; etab++) {
      // plot individual histograms used to measure resolutions in that eta bin
      for(Int_t etb=0; etb<nrETBins; etb++) { 
        cEtEta[ro] -> cd(etb*nrEtaBins+etab+1);
        hResEtEtaBin[ro][etab][etb] -> Draw();
      }
      // plot resolution curves in eta bin
      cEta[ro] -> cd(etab+1);
      hResEtaBin[ro][etab] -> Draw();
      maxZ = maxZ < hResEtaBin[ro][etab]->GetMaximum() ? hResEtaBin[ro][etab]->GetMaximum() : maxZ;
      minZ = minZ > hResEtaBin[ro][etab]->GetMinimum() ? hResEtaBin[ro][etab]->GetMinimum() : minZ;
      // draw the NN 1D curve
      Double_t v[2];
      Double_t* vx = new Double_t[hResEtaBin[ro][etab]->GetNbinsX()+2];
      Double_t* vy = new Double_t[hResEtaBin[ro][etab]->GetNbinsX()+2];
      Double_t* out= new Double_t[hResEtaBin[ro][etab]->GetNbinsX()+2];
      v[0]=0.;
      v[1]=etabin[etab];
      vx[0]=v[0];
      vy[0]=v[1];
      out[0]=mlp->Evaluate(0, v);
      for (Int_t ix=1; ix<=hResEtaBin[ro][etab]->GetNbinsX(); ix++) {
        v[0]=hResEtaBin[ro][etab]->GetBinCenter(ix); //et
        v[1]=etabin[etab];               //eta
        vx[ix]=v[0];
        vy[ix]=v[1];
        out[ix]=mlp->Evaluate(0, v);
      }
      v[0]=hResEtaBin[ro][etab]->GetBinLowEdge(hResEtaBin[ro][etab]->GetNbinsX()+1);
      v[1]=etabin[etab];
      vx[hResEtaBin[ro][etab]->GetNbinsX()+1]=v[0];
      vy[hResEtaBin[ro][etab]->GetNbinsX()+1]=v[1];
      out[hResEtaBin[ro][etab]->GetNbinsX()+1]=mlp->Evaluate(0, v);
      TGraph* gExtrapolate=new TGraph(hResEtaBin[ro][etab]->GetNbinsX()+2, vx,out);
      gExtrapolate->SetLineWidth(3);
      gExtrapolate->SetLineColor(kRed);
      gExtrapolate->Draw("l");
      delete[] vx;
      delete[] vy;
      delete[] out;
    }
    for(Int_t etab=0; etab<nrEtaBins; etab++) {
      hResEtaBin[ro][etab] -> GetYaxis()->SetRangeUser(minZ*0.9,maxZ*1.1);
    }
    // save the result in ps file
    cEtEta[ro] -> Print(output,"landscape");
    cEta[ro]   -> Print(output,"landscape");
    // draw the 2D curves
    if(nrEtaBins>1) {
      new TCanvas;
      Double_t vvx[225];
      Double_t vvy[225];
      Double_t vout[225];
      Double_t v[2];
      for (Int_t ix=0; ix<15; ix++) {
          v[0]=minX+(maxX-minX)/15.*ix; //et
          for (Int_t iy=0; iy<15; iy++) {
             v[1]=minY+(maxY-minY)/15.*iy; //eta
             Int_t idx=ix*15+iy;
             vvx[idx]=v[0];
             vvy[idx]=v[1];
             vout[idx]=mlp->Evaluate(0, v);
          }
      }
      TGraph2D* g2Extrapolate=new TGraph2D("ANN fit","ANN fit",225, vvx, vvy, vout);
      g2Extrapolate->Draw("surf3");
      gPad->Print(output,"landscape");
    }
    // save the NN
    TString NNoutput = fileName;
    NNoutput.Remove(NNoutput.Index("."),NNoutput.Length());
    NNoutput += resObsName[ro]; NNoutput += "_NN";
    mlp->Export(NNoutput,"c++"); // C++ code
    NNoutput = objectType; NNoutput += resObsName[ro]; NNoutput += "_NN";
    mlp->Write(NNoutput,TObject::kOverwrite); // add it to the ROOT file
  } //loop on variables
  // close the ps file
  c -> Print(output+"]","landscape");

} 

