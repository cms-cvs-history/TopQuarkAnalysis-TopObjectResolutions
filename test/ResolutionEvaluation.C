#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMultiLayerPerceptron.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TKey.h"

/////////////////////////////////////////////////////////////////
// MACRO TO DRAW AND EVALUATE THE RESOLUTION PLOTS             //
/////////////////////////////////////////////////////////////////

void ResolutionEvaluation(TString fileName="Resolutions_tau.root", int absOrRel=0, bool doNN=true, TString output="ResolutionEvaluation.ps") {
  
  //TString fileName   = "Resolutions_tau.root";
  //int absOrRel       = 0;	// 0:abs, 1:rel, 2: both
  //bool doNN          = true;

  // define all histo's
  Double_t etabin[10] = {0.0875,0.26875,0.45,0.63125,0.825,1.0375,1.275,1.55,1.8875,2.2875};
  int Nepoch = 1000;
  TH1F  *hResEtEtaBin[10][20][20][2], *hResEtaBin[10][20][2],  *hResPar[10][2][3];
  TString objectType = fileName;
  objectType.Remove(0,objectType.Index("_")+1);
  if(objectType.Index("_")>0){
    objectType.Remove(objectType.Index("_"),objectType.Length());
  } else {
    objectType.Remove(objectType.Index("."),objectType.Length());    
  }
  int aorBegin = absOrRel;
  int aorEnd   = absOrRel+1;
  if(absOrRel == 2) { aorBegin = 0; aorEnd = 2; };
  
  // read in histo's
  TFile *fIn = new TFile(fileName,"UPDATE");
  int     nrEtaBins, nrETBins; 
  TString  resObsName[6] = {"pres","eres","thres","phres","etres","etares"};
  TString def[2] = {"_abs","_rel"};
  for(Int_t ro=0; ro<6; ro++) {
    for(Int_t aor=aorBegin; aor<aorEnd; aor++) {
      if(objectType != "met"){
        for(Int_t par=0; par<3; par++) {
          TString obsName3 = objectType; obsName3 += resObsName[ro]; obsName3 += "_par"; 
	  obsName3 += par; obsName3 += def[aor];
          hResPar[ro][aor][par] = (TH1F*) (fIn->GetKey(obsName3)->ReadObj());
          nrEtaBins = hResPar[ro][aor][par]->GetEntries();
        }
      } else {
        nrEtaBins = 1;
      }
      for(Int_t etab=0; etab<nrEtaBins; etab++) {	
        TString obsName2 = objectType; obsName2 += resObsName[ro]; obsName2 += "_etabin"; 
	        obsName2 += etab; obsName2 += def[aor];
        hResEtaBin[ro][etab][aor] = (TH1F*) (fIn->GetKey(obsName2)->ReadObj());
        nrETBins = hResEtaBin[ro][etab][aor]->GetEntries();
        for(Int_t etb=0; etb<nrETBins; etb++) {
          TString obsName = objectType; obsName += resObsName[ro]; obsName += "_etabin"; 
	          obsName += etab; obsName += "_etbin"; obsName += etb; obsName += def[aor];
	  hResEtEtaBin[ro][etab][etb][aor] = (TH1F*) (fIn->GetKey(obsName)->ReadObj());
        } 
      }
    }
  }
  
  // Draw resolutions in (ET,Eta)-bins & resolutions vs ET for different Eta bins
  TCanvas *c = new TCanvas("dummy","",1);
  c -> Print(output+"[","landscape");
  TCanvas *cEtEta[10][2], *cEta[10][2], *cPar[10][2];
  for(Int_t ro=0; ro<6; ro++) { //loop on variables
    for(Int_t aor=aorBegin; aor<aorEnd; aor++) { //loop on aor
      TString cName   = "Res_vs_ET_for_different_Eta_bins_for_ResoObs_"; cName += resObsName[ro]; cName += "_"; cName += def[aor];
      TString cName2  = "EtEtaBinResolutions_for_ResoObs_"; cName2 += resObsName[ro]; cName2 += "_"; cName2 += def[aor];
      cEtEta[ro][aor] =  new TCanvas(cName,cName,1);
      cEtEta[ro][aor] -> Divide(nrEtaBins,nrETBins);
      cEta[ro][aor]   =  new TCanvas(cName2,cName2,1);
      cEta[ro][aor]   -> Divide(int(sqrt(nrEtaBins)),int(ceil(nrEtaBins/float(int(sqrt(nrEtaBins))))));
      // do the NN part
      TMultiLayerPerceptron* mlp = NULL;
      TTree* tResVar = NULL;
      if(doNN) {
        tResVar = (TTree*)fIn->GetKey("tResVar")->ReadObj();
        mlp = new TMultiLayerPerceptron("Et,Eta:5:2:value","1/error",tResVar,
	                                Form("aor==%d && ro==%d",aor,ro),
	                                Form("aor==%d && ro==%d",aor,ro));
        mlp->Train(Nepoch,"graph update=10");
	gPad->Print(output,"landscape");
      }
      // loop on eta bins
      for(Int_t etab=0; etab<nrEtaBins; etab++) {
        // plot individual histograms used to measure resolutions in that eta bin
        for(Int_t etb=0; etb<nrETBins; etb++) { 
          cEtEta[ro][aor] -> cd(etb*nrEtaBins+etab+1);
	  hResEtEtaBin[ro][etab][etb][aor] -> Draw();
	}
	// plot resolution curves in eta bin
        cEta[ro][aor] -> cd(etab+1);
        hResEtaBin[ro][etab][aor] -> Draw();
	// draw the NN 1D curve
	Double_t v[2];
	Double_t* vx = new Double_t[hResEtaBin[ro][etab][aor]->GetNbinsX()+2];
	Double_t* vy = new Double_t[hResEtaBin[ro][etab][aor]->GetNbinsX()+2];
	Double_t* out= new Double_t[hResEtaBin[ro][etab][aor]->GetNbinsX()+2];
        v[0]=0.;
        v[1]=etabin[etab];
        vx[0]=v[0];
        vy[0]=v[1];
        out[0]=mlp->Evaluate(0, v);
        for (Int_t ix=1; ix<=hResEtaBin[ro][etab][aor]->GetNbinsX(); ix++) {
          v[0]=hResEtaBin[ro][etab][aor]->GetBinCenter(ix); //et
          v[1]=etabin[etab];               //eta
          vx[ix]=v[0];
          vy[ix]=v[1];
          out[ix]=mlp->Evaluate(0, v);
        }
        v[0]=hResEtaBin[ro][etab][aor]->GetBinLowEdge(hResEtaBin[ro][etab][aor]->GetNbinsX()+1);
        v[1]=etabin[etab];
        vx[hResEtaBin[ro][etab][aor]->GetNbinsX()+1]=v[0];
        vy[hResEtaBin[ro][etab][aor]->GetNbinsX()+1]=v[1];
        out[hResEtaBin[ro][etab][aor]->GetNbinsX()+1]=mlp->Evaluate(0, v);
        TGraph* gExtrapolate=new TGraph(hResEtaBin[ro][etab][aor]->GetNbinsX()+2, vx,out);
        gExtrapolate->SetLineWidth(3);
        gExtrapolate->SetLineColor(kRed);
        gExtrapolate->Draw("l");
	delete[] vx;
	delete[] vy;
	delete[] out;
      }
      // draw the additional canvas with variables resolutions
      if(nrEtaBins > 1){
        TString cName3  = "Par_vs_Eta_for_ResoObs_"; cName3 += resObsName[ro]; cName3 += "_"; cName3 += def[aor];
        cPar[ro][aor]   =  new TCanvas(cName3,cName3,1);
        cPar[ro][aor]   -> Divide(2,2);
        for(Int_t par=0; par<3; par++) {	
          cPar[ro][aor] -> cd(par+1);
          hResPar[ro][aor][par] -> Draw();
        }
      }
      // save the result in ps file
      cEtEta[ro][aor] -> Print(output,"landscape");
      cEta[ro][aor]   -> Print(output,"landscape");
      if(nrEtaBins > 1){
        cPar[ro][aor] -> Print(output,"landscape");
      }
      // draw the 2D curves
      new TCanvas;
      tResVar->Draw("value:Et:Eta>>new3d","ro==0","goff");
      TProfile2D* prof2D = ((TH3F*)gDirectory->GetObjectUnchecked("new3d"))->Project3DProfile("xy");
      Double_t maxX = prof2D->GetXaxis()->GetXmax();
      Double_t maxY = prof2D->GetYaxis()->GetXmax();
      Double_t minX = prof2D->GetXaxis()->GetXmin();
      Double_t minY = prof2D->GetYaxis()->GetXmin();
      Double_t vx[225];
      Double_t vy[225];
      Double_t out[225];
      Double_t v[2];
      for (Int_t ix=0; ix<15; ix++) {
          v[0]=minX+(maxX-minX)/15.*ix; //et
          for (Int_t iy=0; iy<15; iy++) {
             v[1]=minY+(maxY-minY)/15.*iy; //eta
             Int_t idx=ix*15+iy;
             vx[idx]=v[0];
             vy[idx]=v[1];
             out[idx]=mlp->Evaluate(0, v);
          }
      }
      TGraph2D* g2Extrapolate=new TGraph2D("ANN fit","ANN fit",225, vx, vy, out);
      g2Extrapolate->Draw("surf3");
      gPad->SetLogz();
      gPad->Print(output,"landscape");
      // save the NN
      TString NNoutput = fileName;
      NNoutput.Remove(NNoutput.Index("."),NNoutput.Length());
      NNoutput += resObsName[ro]; NNoutput += def[aor]; NNoutput += "_NN";
      mlp->Export(NNoutput,"c++"); // C++ code
      NNoutput = objectType; NNoutput += resObsName[ro]; NNoutput += "_NN"; NNoutput += def[aor];
      mlp->Write(NNoutput); // add it to the ROOT file
    } //loop on aor
  } //loop on variables
  // close the ps file
  c -> Print(output+"]","landscape");
} 

