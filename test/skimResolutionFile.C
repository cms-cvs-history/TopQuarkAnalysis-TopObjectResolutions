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
    if(key->GetCycle()==2 && name.Contains("etabin") && !name.Contains("etbin")) {
      //copy to the output file
      ((TH1*)key->ReadObj())->SetDirectory(&output);
    }
  }
  // write the output
  output.Write();
  output.Close();
  input.Close();
}
