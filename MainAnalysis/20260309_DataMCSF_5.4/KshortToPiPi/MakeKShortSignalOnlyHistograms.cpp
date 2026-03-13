#include <cmath>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TString.h"
#include "TTree.h"

namespace {
constexpr double kPionMass = 0.13957039;
constexpr double kMassWindowMin = 0.30;
constexpr double kMassWindowMax = 1.00;
constexpr int kMassBins = 280;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMatchAngleMax = 0.025;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;
constexpr int kMaxKShort = 4096;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
};

double buildMass(const TrackKinematics& t1, const TrackKinematics& t2) {
  const double p1sq = t1.px * t1.px + t1.py * t1.py + t1.pz * t1.pz;
  const double p2sq = t2.px * t2.px + t2.py * t2.py + t2.pz * t2.pz;
  const double e1 = std::sqrt(p1sq + kPionMass * kPionMass);
  const double e2 = std::sqrt(p2sq + kPionMass * kPionMass);
  const double px = t1.px + t2.px;
  const double py = t1.py + t2.py;
  const double pz = t1.pz + t2.pz;
  const double e = e1 + e2;
  const double m2 = e * e - (px * px + py * py + pz * pz);
  return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

bool passAcceptance(const TrackKinematics& t) {
  const double p = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz);
  if (p <= 0.0) return false;
  const double absCosTheta = std::fabs(t.pz / p);
  return (absCosTheta >= kAbsCosMin && absCosTheta <= kAbsCosMax);
}

std::string getArgument(int argc, char* argv[], const std::string& option,
                        const std::string& defaultValue) {
  for (int i = 1; i + 1 < argc; ++i)
    if (argv[i] == option) return argv[i + 1];
  return defaultValue;
}

double getDoubleArgument(int argc, char* argv[], const std::string& option, double defaultValue) {
  const std::string value = getArgument(argc, argv, option, "");
  return value.empty() ? defaultValue : std::stod(value);
}
}  // namespace

int main(int argc, char* argv[]) {
  const std::string inputFileName =
      getArgument(argc, argv, "--input", "../../../../Samples/merged_mc_v2.3.root");
  const std::string outputFileName =
      getArgument(argc, argv, "--output", "KShortSignalOnlyHistograms.root");
  const std::string treeName = getArgument(argc, argv, "--tree", "Tree");
  const double massMin = getDoubleArgument(argc, argv, "--mass-min", kMassWindowMin);
  const double massMax = getDoubleArgument(argc, argv, "--mass-max", kMassWindowMax);
  const double matchAngleMax = getDoubleArgument(argc, argv, "--match-angle-max", kMatchAngleMax);

  TFile inputFile(inputFileName.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Error: cannot open input file " << inputFileName << std::endl;
    return 1;
  }

  TTree* tree = nullptr;
  inputFile.GetObject(treeName.c_str(), tree);
  if (tree == nullptr) {
    std::cerr << "Error: cannot find tree '" << treeName << "' in " << inputFileName << std::endl;
    return 1;
  }

  long long nReco = 0;
  double recoPx[kMaxReco] = {0.0};
  double recoPy[kMaxReco] = {0.0};
  double recoPz[kMaxReco] = {0.0};
  double recoCharge[kMaxReco] = {0.0};
  long long recoPIDPion[kMaxReco] = {0};
  long long recoGoodTrack[kMaxReco] = {0};

  long long nKShort = 0;
  long long ksReco1ID[kMaxKShort] = {0};
  long long ksReco2ID[kMaxKShort] = {0};
  double ksReco1Angle[kMaxKShort] = {0.0};
  double ksReco2Angle[kMaxKShort] = {0.0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  tree->SetBranchAddress("NKShort", &nKShort);
  tree->SetBranchAddress("KShortReco1ID[NKShort]", ksReco1ID);
  tree->SetBranchAddress("KShortReco2ID[NKShort]", ksReco2ID);
  tree->SetBranchAddress("KShortReco1Angle[NKShort]", ksReco1Angle);
  tree->SetBranchAddress("KShortReco2Angle[NKShort]", ksReco2Angle);

  TH1D hMass1Tag("hKShortMass1Tag",
                 "K^{0}_{S} signal-only MC, 1-tag; m(#pi^{+}#pi^{-}) [GeV]; Candidates / bin",
                 kMassBins, massMin, massMax);
  TH1D hMass2Tag("hKShortMass2Tag",
                 "K^{0}_{S} signal-only MC, 2-tag; m(#pi^{+}#pi^{-}) [GeV]; Candidates / bin",
                 kMassBins, massMin, massMax);
  TH1D hMassAccepted("hKShortMassAccepted",
                     "K^{0}_{S} signal-only MC, accepted; m(#pi^{+}#pi^{-}) [GeV]; Candidates / bin",
                     kMassBins, massMin, massMax);

  long long totalCandidates = 0;
  long long passValidRecoID = 0;
  long long passMatchAngle = 0;
  long long passGoodTrack = 0;
  long long passAcceptanceBoth = 0;
  long long passOppositeCharge = 0;
  long long count1Tag = 0;
  long long count2Tag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    for (long long iKS = 0; iKS < nKShort; ++iKS) {
      totalCandidates++;

      const long long reco1 = ksReco1ID[iKS];
      const long long reco2 = ksReco2ID[iKS];
      if (reco1 < 0 || reco2 < 0 || reco1 >= nReco || reco2 >= nReco) continue;
      passValidRecoID++;

      if (ksReco1Angle[iKS] >= matchAngleMax || ksReco2Angle[iKS] >= matchAngleMax)
        continue;
      passMatchAngle++;

      if (recoGoodTrack[reco1] != 1 || recoGoodTrack[reco2] != 1) continue;
      passGoodTrack++;

      const TrackKinematics track1{recoPx[reco1], recoPy[reco1], recoPz[reco1]};
      const TrackKinematics track2{recoPx[reco2], recoPy[reco2], recoPz[reco2]};
      if (!passAcceptance(track1) || !passAcceptance(track2)) continue;
      passAcceptanceBoth++;

      if (recoCharge[reco1] * recoCharge[reco2] >= 0) continue;
      passOppositeCharge++;

      const double mass = buildMass(track1, track2);
      hMassAccepted.Fill(mass);

      int nTagged = 0;
      if (recoPIDPion[reco1] >= kPionTagThreshold) nTagged++;
      if (recoPIDPion[reco2] >= kPionTagThreshold) nTagged++;

      if (nTagged == 1) {
        hMass1Tag.Fill(mass);
        count1Tag++;
      }
      if (nTagged == 2) {
        hMass2Tag.Fill(mass);
        count2Tag++;
      }
    }
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hMass1Tag.Write();
  hMass2Tag.Write();
  hMassAccepted.Write();

  TNamed selection("SelectionSummary",
                   Form("Signal-only KShort->pipi from MC: valid KShortReco IDs, "
                        "KShortReco1Angle<%.4f, KShortReco2Angle<%.4f, both RecoGoodTrack==1, "
                        "both 0.15<=|cos(theta)|<=0.675, opposite charge, pion mass hypothesis "
                        "on both daughters, tag if RecoPIDPion>=2, hist range %.3f-%.3f GeV",
                        matchAngleMax, matchAngleMax, massMin, massMax));
  selection.Write();

  TParameter<long long>("TotalKShortCandidates", totalCandidates).Write();
  TParameter<long long>("PassValidRecoID", passValidRecoID).Write();
  TParameter<long long>("PassMatchAngle", passMatchAngle).Write();
  TParameter<long long>("PassGoodTrack", passGoodTrack).Write();
  TParameter<long long>("PassAcceptanceBoth", passAcceptanceBoth).Write();
  TParameter<long long>("PassOppositeCharge", passOppositeCharge).Write();
  TParameter<long long>("Count1Tag", count1Tag).Write();
  TParameter<long long>("Count2Tag", count2Tag).Write();
  TParameter<double>("MassMin", massMin).Write();
  TParameter<double>("MassMax", massMax).Write();
  TParameter<double>("MatchAngleMax", matchAngleMax).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Total KShort candidates: " << totalCandidates << std::endl;
  std::cout << "  Pass valid reco IDs:     " << passValidRecoID << std::endl;
  std::cout << "  Pass match angles:       " << passMatchAngle << std::endl;
  std::cout << "  Pass good tracks:        " << passGoodTrack << std::endl;
  std::cout << "  Pass acceptance:         " << passAcceptanceBoth << std::endl;
  std::cout << "  Pass opposite charge:    " << passOppositeCharge << std::endl;
  std::cout << "  1-tag entries:           " << count1Tag << std::endl;
  std::cout << "  2-tag entries:           " << count2Tag << std::endl;
  return 0;
}
