{
  /////////////////////////////////////////////////////////////////
  // MACRO TO DRAW AND EVALUATE THE RESOLUTION PLOTS		 //
  /////////////////////////////////////////////////////////////////
  
  TString fileName   = "../data/Resolutions_lJets_MCJetCorJetIcone5.root";
  int absOrRel       = 0;	// 0:abs, 1:rel, 2: both

  /////////////////////////////////////////////////////////////////
  
  
  
  // define all histo's
  TH1F  *hResEtEtaBin[10][20][20][2], *hResEtaBin[10][20][2],  *hResPar[10][2][3];
  
  TString objectType = fileName;
  objectType.Remove(0,objectType.Index("_")+1);
  if(objectType.Index("_")>0){
    objectType.Remove(objectType.Index("_"),objectType.Length());
  }
  else
  {
    objectType.Remove(objectType.Index("."),objectType.Length());    
  }
  
  int aorBegin = absOrRel;
  int aorEnd   = absOrRel+1;
  if(absOrRel == 2) { aorBegin = 0; aorEnd = 2; };
  
  
  // read in histo's
  TFile *fIn = new TFile(fileName);
  int     nrEtaBins, nrETBins; 
  TString  resObsName[6] = {"pres","eres","thres","phres","etres","etares"};
  TString def[2] = {"_abs","_rel"};
  for(Int_t ro=0; ro<6; ro++) {
    for(Int_t aor=aorBegin; aor<aorEnd; aor++) {
      if(objectType != "met"){
        for(Int_t par=0; par<3; par++) {
          TString obsName3 = objectType; obsName3 += resObsName[ro]; obsName3 += "_par"; obsName3 += par; obsName3 += def[aor];
          hResPar[ro][aor][par] = (TH1F*) (fIn->GetKey(obsName3)->ReadObj());
          nrEtaBins = hResPar[ro][aor][par]->GetEntries();
        }
      }
      else
      {
        nrEtaBins = 1;
      }
      
      for(Int_t etab=0; etab<nrEtaBins; etab++) {	
        TString obsName2 = objectType; obsName2 += resObsName[ro]; obsName2 += "_etabin"; obsName2 += etab; obsName2 += def[aor];
        hResEtaBin[ro][etab][aor] = (TH1F*) (fIn->GetKey(obsName2)->ReadObj());
        nrETBins = hResEtaBin[ro][etab][aor]->GetEntries();
        for(Int_t etb=0; etb<nrETBins; etb++) {
          TString obsName = objectType; obsName += resObsName[ro]; obsName += "_etabin"; obsName += etab; obsName += "_etbin"; obsName += etb; obsName += def[aor];           
	  hResEtEtaBin[ro][etab][etb][aor] = (TH1F*) (fIn->GetKey(obsName)->ReadObj());
        } 
      }
    }
  }
  
  // Draw resolutions in (ET,Eta)-bins & resolutions vs ET for different Eta bins
  TCanvas *c = new TCanvas("dummy","",1);
  c -> Print("ResolutionEvaluation.ps[","landscape");
  TCanvas *cEtEta[10][2], *cEta[10][2], *cPar[10][2];
  for(Int_t ro=0; ro<6; ro++) {
    for(Int_t aor=aorBegin; aor<aorEnd; aor++) {
      TString cName   = "Res_vs_ET_for_different_Eta_bins_for_ResoObs_"; cName += resObsName[ro]; cName += "_"; cName += def[aor];
      TString cName2  = "EtEtaBinResolutions_for_ResoObs_"; cName2 += resObsName[ro]; cName2 += "_"; cName2 += def[aor];
      cEtEta[ro][aor] =  new TCanvas(cName,cName,1);
      cEtEta[ro][aor] -> Divide(nrEtaBins,nrETBins);
      cEta[ro][aor]   =  new TCanvas(cName2,cName2,1);
      cEta[ro][aor]   -> Divide((1+(int)sqrt(nrEtaBins)),(1+(int)sqrt(nrEtaBins)));
      for(Int_t etab=0; etab<nrEtaBins; etab++) {	
        for(Int_t etb=0; etb<nrETBins; etb++) {
          cEtEta[ro][aor] -> cd(etb*nrEtaBins+etab+1);
	  hResEtEtaBin[ro][etab][etb][aor] -> Draw();
	}
        cEta[ro][aor] -> cd(etab+1);
        hResEtaBin[ro][etab][aor] -> Draw();
      }
      if(nrEtaBins > 1){
        TString cName3  = "Par_vs_Eta_for_ResoObs_"; cName3 += resObsName[ro]; cName3 += "_"; cName3 += def[aor];
        cPar[ro][aor]   =  new TCanvas(cName3,cName3,1);
        cPar[ro][aor]   -> Divide(2,2);
        for(Int_t par=0; par<3; par++) {	
          cPar[ro][aor] -> cd(par+1);
          hResPar[ro][aor][par] -> Draw();
        }
      }
      cEtEta[ro][aor] -> Print("ResolutionEvaluation.ps","landscape");
      cEta[ro][aor]   -> Print("ResolutionEvaluation.ps","landscape");
      if(nrEtaBins > 1){
        cPar[ro][aor]   -> Print("ResolutionEvaluation.ps","landscape");
      }
    }
  }

  c -> Print("ResolutionEvaluation.ps]","landscape");
      	
} 
