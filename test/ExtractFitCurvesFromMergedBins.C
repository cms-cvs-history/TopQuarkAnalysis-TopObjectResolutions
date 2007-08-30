// simple script to extract the resolution versus ET graphs for the different eta bins 
// starting from the (Et,Eta)-resolution bins.
// Also the tree will be recalculated.
// By default, no (Et,Eta)-bins are saved to the final root-file to reduce its size
// Author: J.Heyninck
// run as :
// root -l -q 'ExtractFitCurvesFromMergedBins.C("../data/ResolutionBins.root","../data/Resolutions.root")'

void ExtractFitCurvesFromMergedBins(const char* fileIn, const char* fileOut) {

  bool      saveBinPlots 	= true;
  TString   resObsName[8] 	= {"_ares","_bres","_cres","_dres","_thres","_phres","_etres","_etares"};

  TFile input(fileIn,"READ");
  TFile output(fileOut,"RECREATE");
  input.cd();
  
  // find et and eta bins values
  vector<double> eTbinVals_, etabinVals_;
  TString objectName;
  //vector<TString> resVsEtNames;
  TList* keys = input.GetListOfKeys();
  TIter nextitem(keys);
  TKey* key = NULL;
  bool found = false;
  while((key = (TKey*)nextitem()) && !found) {
    TString name = key->GetName();
    if( name.Contains("_ares") && name.Contains("etabin0") && (!name.Contains("etbin"))) {
      objectName = name; objectName.Remove(objectName.Index("_"),objectName.Length());
      TH1F *h = (TH1F*) (input.GetKey(key->GetName())->ReadObj());
      for(size_t b=1; b<=h->GetNbinsX(); b++) eTbinVals_.push_back(h->GetXaxis()->GetBinLowEdge(b));
      eTbinVals_.push_back(h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));
      found = true;
    }
  }
  TH1F *tmpEta = (TH1F*) (input.GetKey("hEtaBins")->ReadObj());
  for(size_t b=1; b<=tmpEta->GetNbinsX(); b++) etabinVals_.push_back(tmpEta->GetXaxis()->GetBinLowEdge(b));
  etabinVals_.push_back(tmpEta->GetXaxis()->GetBinUpEdge(tmpEta->GetNbinsX()));
  int etanrbins = etabinVals_.size()-1;
  int etnrbins  = eTbinVals_.size()-1;
  
  cout<<"Found "<<etanrbins<< " eta-bins with edges: ( ";
  for(size_t u=0; u<etabinVals_.size(); u++) cout<<etabinVals_[u]<<", ";
  cout<<"\b\b )"<<endl;
  cout<<"Found "<<etnrbins<< " eT-bins with edges: ( ";
  for(size_t u=0; u<eTbinVals_.size(); u++) cout<<eTbinVals_[u]<<", ";
  cout<<"\b\b )"<<endl;
 
 
  
  // Remake the res_vs_ET histograms for all etabins and observables
  Int_t ro=0;
  Double_t et=0.;
  Double_t eta=0.;
  Double_t value,error;
  // CD: create the output tree
  TTree* tResVar = new TTree("tResVar","Resolution tree");
  tResVar->Branch("Et",&et,"Et/D");
  tResVar->Branch("Eta",&eta,"Eta/D");
  tResVar->Branch("ro",&ro,"ro/I");
  tResVar->Branch("value",&value,"value/D");
  tResVar->Branch("error",&error,"error/D");
  for(ro=0; ro<8; ro++) {
    for(int etab=0; etab<etanrbins; etab++) {	
      //CD set eta at the center of the bin
      eta = (etanrbins > 1) ? ((etabinVals_[etab]+etabinVals_[etab+1])/2.) : 2.5 ; 
      TString resEtabinName = objectName; resEtabinName += resObsName[ro]; resEtabinName += "_etabin"; resEtabinName += etab;
      TH1F *res = (TH1F*) (input.GetKey(resEtabinName)->ReadObj());
      res->Reset();
      for(int etb=0; etb<etnrbins; etb++) {
	//CD set et at the center of the bin
	et = (eTbinVals_[etb]+eTbinVals_[etb+1])/2.;
        TString eTEtabinName = resEtabinName; eTEtabinName += "_etbin";  eTEtabinName+= etb;
        TH1F *h = (TH1F*) (input.GetKey(eTEtabinName)->ReadObj());
        double maxcontent = 0.;
        int maxbin = 0;
        for(int nb=1; nb<h->GetNbinsX(); nb ++){
 	  if (h->GetBinContent(nb)>maxcontent) {
	    maxcontent = h->GetBinContent(nb);
	    maxbin = nb;
          }
        }
        int range = (int)(h->GetNbinsX()/6); //in order that ~1/3 of X-axis range is fitted
        TF1 f("f","gaus",h->GetBinCenter(maxbin-range),h->GetBinCenter(maxbin+range));
        f.SetParameters(h -> GetMaximum(), h -> GetMean(), h -> GetRMS());
        h -> Fit("f","RQ");
        res -> SetBinContent(etb+1,f.GetParameter(2));
        res -> SetBinError(etb+1,f.GetParError(2));
	//CD: Fill the tree
	value = f.GetParameter(2); //parameter value
	error = f.GetParError(2);  //parameter error
	tResVar->Fill();
        if(saveBinPlots) h->SetDirectory(&output); 
      }
      //CD: add a fake entry in et=0 for the NN training
      // for that, use a linear extrapolation.
      et = 0.;
      value = ((eTbinVals_[0]+eTbinVals_[1])/2.)*(res->GetBinContent(1)-res->GetBinContent(3))/((eTbinVals_[2]-eTbinVals_[0])/2.)+res->GetBinContent(1);
      error = res->GetBinContent(1)+res->GetBinContent(2);
      tResVar->Fill();
      // standard fit
      TString fname = "F_"; fname += resEtabinName;
      TF1 f2(fname,"[0]+[1]*exp(-[2]*x)",res->GetXaxis()->GetXmin(),res->GetXaxis()->GetXmax());
      res->Fit(f2.GetName(),"RQ");
      res->SetDirectory(&output);
    }
  }
  // write the output
  output.cd();
  tmpEta->Write();
  tResVar->Write();
  output.Write();
  output.Close();
  input.Close();

}
