#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace {
constexpr double kPionMass = 0.13957039;
constexpr double kKShortMassWindowMin = 0.30;
constexpr double kKShortMassWindowMax = 1.00;
constexpr int kKShortMassBins = 280;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double charge = 0.0;
  long long pionTag = 0;
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
      getArgument(argc, argv, "--output", "KShortSBHistograms.root");
  const std::string treeName = getArgument(argc, argv, "--tree", "Tree");
  const double massMin = getDoubleArgument(argc, argv, "--mass-min", kKShortMassWindowMin);
  const double massMax = getDoubleArgument(argc, argv, "--mass-max", kKShortMassWindowMax);

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

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  TH1D hMass1Tag("hKShortSBMass1Tag",
                 "K_{S}^{0} same-event reco pairs, 1-tag; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin",
                 kKShortMassBins, massMin, massMax);
  TH1D hMass2Tag("hKShortSBMass2Tag",
                 "K_{S}^{0} same-event reco pairs, 2-tag; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin",
                 kKShortMassBins, massMin, massMax);
  TH1D hMassAccepted("hKShortSBMassAccepted",
                     "K_{S}^{0} same-event reco pairs, accepted; m(#pi^{+}#pi^{-}) [GeV]; Pairs / bin",
                     kKShortMassBins, massMin, massMax);

  long long acceptedTracks = 0;
  long long totalOppositeSignPairs = 0;
  long long count1Tag = 0;
  long long count2Tag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    std::vector<TrackKinematics> tracks;
    tracks.reserve(nReco);

    for (long long i = 0; i < nReco; ++i) {
      if (recoGoodTrack[i] != 1) continue;
      if (recoCharge[i] == 0) continue;
      TrackKinematics t{recoPx[i], recoPy[i], recoPz[i], recoCharge[i], recoPIDPion[i]};
      if (!passAcceptance(t)) continue;
      tracks.push_back(t);
      acceptedTracks++;
    }

    const int trackCount = static_cast<int>(tracks.size());
    for (int i = 0; i < trackCount; ++i) {
      for (int j = i + 1; j < trackCount; ++j) {
        if (tracks[i].charge * tracks[j].charge >= 0) continue;
        totalOppositeSignPairs++;

        const double mass = buildMass(tracks[i], tracks[j]);
        hMassAccepted.Fill(mass);

        int nTagged = 0;
        if (tracks[i].pionTag >= kPionTagThreshold) nTagged++;
        if (tracks[j].pionTag >= kPionTagThreshold) nTagged++;

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
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hMass1Tag.Write();
  hMass2Tag.Write();
  hMassAccepted.Write();

  TNamed selection("SelectionSummary",
                   Form("Reco-only Kshort S+B pairs from same event: RecoGoodTrack==1, nonzero "
                        "charge, 0.15<=|cos(theta)|<=0.675 on both tracks, opposite charge, pion "
                        "mass hypothesis on both tracks, tag if RecoPIDPion>=2, hist range %.3f-%.3f GeV",
                        massMin, massMax));
  selection.Write();
  TParameter<long long>("AcceptedTracks", acceptedTracks).Write();
  TParameter<long long>("TotalOppositeSignPairs", totalOppositeSignPairs).Write();
  TParameter<long long>("Count1Tag", count1Tag).Write();
  TParameter<long long>("Count2Tag", count2Tag).Write();
  TParameter<double>("MassMin", massMin).Write();
  TParameter<double>("MassMax", massMax).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Accepted tracks:      " << acceptedTracks << std::endl;
  std::cout << "  Opposite-sign pairs:  " << totalOppositeSignPairs << std::endl;
  std::cout << "  1-tag pairs:          " << count1Tag << std::endl;
  std::cout << "  2-tag pairs:          " << count2Tag << std::endl;
  return 0;
}
