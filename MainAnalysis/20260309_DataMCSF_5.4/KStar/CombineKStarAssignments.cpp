#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"

int main(int argc, char* argv[]) {
  const std::string inputFileName =
      (argc > 1) ? argv[1] : "KStarWrongAssignmentHistograms.root";
  const std::string outputFileName =
      (argc > 2) ? argv[2] : "KStarCombinedAssignmentHistograms.root";

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
    return 1;
  }

  TH1D* hNominalNoPionTag = nullptr;
  TH1D* hSwappedNoPionTag = nullptr;
  TH1D* hNominalWithPionTag = nullptr;
  TH1D* hSwappedWithPionTag = nullptr;
  inputFile.GetObject("hNominalNoPionTag", hNominalNoPionTag);
  inputFile.GetObject("hSwappedNoPionTag", hSwappedNoPionTag);
  inputFile.GetObject("hNominalWithPionTag", hNominalWithPionTag);
  inputFile.GetObject("hSwappedWithPionTag", hSwappedWithPionTag);
  if (hNominalNoPionTag == nullptr || hSwappedNoPionTag == nullptr ||
      hNominalWithPionTag == nullptr || hSwappedWithPionTag == nullptr) {
    std::cerr << "Error: required histograms not found in " << inputFileName << std::endl;
    return 1;
  }

  TH1D hCombinedKaonTag(*hNominalNoPionTag);
  hCombinedKaonTag.SetName("hKStarMassKaonTag");
  hCombinedKaonTag.SetTitle("K^{*} signal-only MC, combined kaon tag; m(K#pi) [GeV]; Candidates / bin");
  hCombinedKaonTag.Add(hSwappedNoPionTag);

  TH1D hCombinedKaonPionTag(*hNominalWithPionTag);
  hCombinedKaonPionTag.SetName("hKStarMassKaonPionTag");
  hCombinedKaonPionTag.SetTitle("K^{*} signal-only MC, combined kaon+pion tag; m(K#pi) [GeV]; Candidates / bin");
  hCombinedKaonPionTag.Add(hSwappedWithPionTag);

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hCombinedKaonTag.Write();
  hCombinedKaonPionTag.Write();
  TNamed selection("SelectionSummary",
                   "Combined KStar signal-only histograms from nominal and swapped assignments");
  selection.Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Combined kaon-tag entries:      " << hCombinedKaonTag.Integral() << std::endl;
  std::cout << "  Combined kaon+pion-tag entries: " << hCombinedKaonPionTag.Integral() << std::endl;
  return 0;
}
