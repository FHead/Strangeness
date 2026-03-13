#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"

int main(int argc, char* argv[]) {
  const std::string nominalFileName =
      (argc > 1) ? argv[1] : "D0SignalOnlyHistograms.root";
  const std::string swappedFileName =
      (argc > 2) ? argv[2] : "D0WrongAssignmentHistograms.root";
  const std::string outputFileName =
      (argc > 3) ? argv[3] : "D0CombinedAssignmentHistograms.root";

  TFile nominalFile(nominalFileName.c_str(), "READ");
  TFile swappedFile(swappedFileName.c_str(), "READ");
  if (nominalFile.IsZombie() || swappedFile.IsZombie()) {
    std::cerr << "Error opening input files" << std::endl;
    return 1;
  }

  TH1D* hNominalKaonTag = nullptr;
  TH1D* hNominalKaonPionTag = nullptr;
  TH1D* hSwappedNoPionTag = nullptr;
  TH1D* hSwappedWithPionTag = nullptr;
  nominalFile.GetObject("hD0MassKaonTag", hNominalKaonTag);
  nominalFile.GetObject("hD0MassKaonPionTag", hNominalKaonPionTag);
  swappedFile.GetObject("hSwappedNoPionTag", hSwappedNoPionTag);
  swappedFile.GetObject("hSwappedWithPionTag", hSwappedWithPionTag);
  if (hNominalKaonTag == nullptr || hNominalKaonPionTag == nullptr ||
      hSwappedNoPionTag == nullptr || hSwappedWithPionTag == nullptr) {
    std::cerr << "Error: required histograms not found" << std::endl;
    return 1;
  }

  TH1D hCombinedKaonTag(*hNominalKaonTag);
  hCombinedKaonTag.SetName("hD0MassKaonTag");
  hCombinedKaonTag.SetTitle("D^{0} signal-only MC, combined kaon tag; m(K#pi) [GeV]; Candidates / bin");
  hCombinedKaonTag.Add(hSwappedNoPionTag);

  TH1D hCombinedKaonPionTag(*hNominalKaonPionTag);
  hCombinedKaonPionTag.SetName("hD0MassKaonPionTag");
  hCombinedKaonPionTag.SetTitle("D^{0} signal-only MC, combined kaon+pion tag; m(K#pi) [GeV]; Candidates / bin");
  hCombinedKaonPionTag.Add(hSwappedWithPionTag);

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hCombinedKaonTag.Write();
  hCombinedKaonPionTag.Write();
  TNamed selection("SelectionSummary",
                   "Combined D0 signal-only histograms from nominal and swapped assignments");
  selection.Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Combined kaon-tag entries:      " << hCombinedKaonTag.Integral() << std::endl;
  std::cout << "  Combined kaon+pion-tag entries: " << hCombinedKaonPionTag.Integral() << std::endl;
  return 0;
}
