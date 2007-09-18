// simple script to skim the resolution files.
// Author: C. Delaere.
// run as :
// root -l -q 'skimResolutionFile.C("../data/Resolutions.root","../data/skimmed_Resolutions.root")'
void skimResolutionFile(const char* fileIn, const char* fileOut) {
  TFile input(fileIn,"READ");
  TFile output(fileOut,"RECREATE");
  input.cd();
  // loop over keys in input
  TList* keys = input.GetListOfKeys();
  TIter nextitem(keys);
  TKey* key = NULL;
  while(key = (TKey*)nextitem()) {
    // select based on the name and on the cycle (select cycle 2 (last one))
    TString name = key->GetName();
    std::cout << name << std::endl;
    if(name.Contains("etabin") && !name.Contains("etbin")) {
      //copy to the output file
      ((TH1*)key->ReadObj())->SetDirectory(&output);
      std::cout << "Kept" << std::endl;
    }
    if(name.Contains("NN")) {
      output.cd();
      ((TMultiLayerPerceptron*)key->ReadObj())->Write(name);
      std::cout << "Kept" << std::endl;
    }
    if(name.Contains("hEtaBins")) {
      output.cd();
      ((TH1*)key->ReadObj())->SetDirectory(&output);
      std::cout << "Kept" << std::endl;
    }
  }
  // write the output
  std::cout << "Save output" << std::endl;
  output.Write();
  output.Close();
  input.Close();
}
