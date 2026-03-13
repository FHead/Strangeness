#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

namespace {
constexpr double kKaonMass = 0.493677;
constexpr double kPionMass = 0.13957039;
constexpr double kMassWindowMin = 0.70;
constexpr double kMassWindowMax = 1.10;
constexpr int kMassBins = 320;
constexpr double kAbsCosMin = 0.15;
constexpr double kAbsCosMax = 0.675;
constexpr double kMatchAngleMax = 0.01;
constexpr double kTruthEnergyMatchMax = 0.025;
constexpr double kTruthMomentumAngleMax = 0.025;
constexpr long long kKaonTagThreshold = 2;
constexpr long long kPionTagThreshold = 2;
constexpr int kMaxReco = 10000;
constexpr int kMaxPhi = 4096;
constexpr int kMaxGen = 20000;

struct TrackKinematics {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
};

double angleBetween(double px1, double py1, double pz1,
                    double px2, double py2, double pz2) {
  const double p1 = std::sqrt(px1 * px1 + py1 * py1 + pz1 * pz1);
  const double p2 = std::sqrt(px2 * px2 + py2 * py2 + pz2 * pz2);
  if (p1 <= 0.0 || p2 <= 0.0) return std::numeric_limits<double>::infinity();
  double c = (px1 * px2 + py1 * py2 + pz1 * pz2) / (p1 * p2);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  return std::acos(c);
}

double buildMass(const TrackKinematics& kaon, const TrackKinematics& pion) {
  const double pK2 = kaon.px * kaon.px + kaon.py * kaon.py + kaon.pz * kaon.pz;
  const double pPi2 = pion.px * pion.px + pion.py * pion.py + pion.pz * pion.pz;
  const double eK = std::sqrt(pK2 + kKaonMass * kKaonMass);
  const double ePi = std::sqrt(pPi2 + kPionMass * kPionMass);
  const double px = kaon.px + pion.px;
  const double py = kaon.py + pion.py;
  const double pz = kaon.pz + pion.pz;
  const double e = eK + ePi;
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
      getArgument(argc, argv, "--output", "PhiWrongAsKStarHistograms.root");
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
  long long recoPIDKaon[kMaxReco] = {0};
  long long recoPIDPion[kMaxReco] = {0};
  long long recoGoodTrack[kMaxReco] = {0};

  long long nPhi = 0;
  long long phiReco1ID[kMaxPhi] = {0};
  long long phiReco2ID[kMaxPhi] = {0};
  double phiReco1Angle[kMaxPhi] = {0.0};
  double phiReco2Angle[kMaxPhi] = {0.0};

  long long nGen = 0;
  double genPx[kMaxGen] = {0.0};
  double genPy[kMaxGen] = {0.0};
  double genPz[kMaxGen] = {0.0};
  double genE[kMaxGen] = {0.0};
  long long genID[kMaxGen] = {0};
  long long genStatus[kMaxGen] = {0};
  long long genMatchIndex[kMaxGen] = {0};
  double genMatchAngle[kMaxGen] = {0.0};

  tree->SetBranchAddress("NReco", &nReco);
  tree->SetBranchAddress("RecoPx", recoPx);
  tree->SetBranchAddress("RecoPy", recoPy);
  tree->SetBranchAddress("RecoPz", recoPz);
  tree->SetBranchAddress("RecoCharge", recoCharge);
  tree->SetBranchAddress("RecoPIDKaon", recoPIDKaon);
  tree->SetBranchAddress("RecoPIDPion", recoPIDPion);
  tree->SetBranchAddress("RecoGoodTrack", recoGoodTrack);

  tree->SetBranchAddress("NPhi", &nPhi);
  tree->SetBranchAddress("PhiReco1ID[NPhi]", phiReco1ID);
  tree->SetBranchAddress("PhiReco2ID[NPhi]", phiReco2ID);
  tree->SetBranchAddress("PhiReco1Angle[NPhi]", phiReco1Angle);
  tree->SetBranchAddress("PhiReco2Angle[NPhi]", phiReco2Angle);
  tree->SetBranchAddress("NGen", &nGen);
  tree->SetBranchAddress("GenPx", genPx);
  tree->SetBranchAddress("GenPy", genPy);
  tree->SetBranchAddress("GenPz", genPz);
  tree->SetBranchAddress("GenE", genE);
  tree->SetBranchAddress("GenID", genID);
  tree->SetBranchAddress("GenStatus", genStatus);
  tree->SetBranchAddress("GenMatchIndex", genMatchIndex);
  tree->SetBranchAddress("GenMatchAngle", genMatchAngle);

  TH1D hMassKaonTag("hPhiWrongAsKStarMassKaonTag",
                    "#phi#rightarrow K^{+}K^{-} treated as K#pi, kaon tag; m(K#pi) [GeV]; Candidates / bin",
                    kMassBins, massMin, massMax);
  TH1D hMassKaonPionTag("hPhiWrongAsKStarMassKaonPionTag",
                        "#phi#rightarrow K^{+}K^{-} treated as K#pi, kaon+pion tag; m(K#pi) [GeV]; Candidates / bin",
                        kMassBins, massMin, massMax);
  TH1D hMassAccepted("hPhiWrongAsKStarMassAccepted",
                     "#phi#rightarrow K^{+}K^{-} treated as K#pi, accepted; m(K#pi) [GeV]; Candidates / bin",
                     kMassBins, massMin, massMax);

  long long totalCandidates = 0;
  long long passValidRecoID = 0;
  long long passMatchAngle = 0;
  long long passGoodTrack = 0;
  long long passAcceptanceBoth = 0;
  long long passOppositeCharge = 0;
  long long passKaonTag = 0;
  long long passKaonPionTag = 0;

  const long long entryCount = tree->GetEntries();
  for (long long entry = 0; entry < entryCount; ++entry) {
    tree->GetEntry(entry);

    if (nPhi > 0) {
      for (long long iPhi = 0; iPhi < nPhi; ++iPhi) {
        totalCandidates++;

        const long long reco1 = phiReco1ID[iPhi];
        const long long reco2 = phiReco2ID[iPhi];
        if (reco1 < 0 || reco2 < 0 || reco1 >= nReco || reco2 >= nReco) continue;
        passValidRecoID++;

        if (phiReco1Angle[iPhi] >= matchAngleMax || phiReco2Angle[iPhi] >= matchAngleMax) continue;
        passMatchAngle++;

        if (recoGoodTrack[reco1] != 1 || recoGoodTrack[reco2] != 1) continue;
        passGoodTrack++;

        const TrackKinematics assumedKaon{recoPx[reco1], recoPy[reco1], recoPz[reco1]};
        const TrackKinematics assumedPion{recoPx[reco2], recoPy[reco2], recoPz[reco2]};
        if (!passAcceptance(assumedKaon) || !passAcceptance(assumedPion)) continue;
        passAcceptanceBoth++;

        if (recoCharge[reco1] * recoCharge[reco2] >= 0) continue;
        passOppositeCharge++;

        const double mass = buildMass(assumedKaon, assumedPion);
        hMassAccepted.Fill(mass);

        if (recoPIDKaon[reco1] < kKaonTagThreshold) continue;
        passKaonTag++;
        hMassKaonTag.Fill(mass);

        if (recoPIDPion[reco2] < kPionTagThreshold) continue;
        passKaonPionTag++;
        hMassKaonPionTag.Fill(mass);
      }
    } else {
      std::vector<long long> kPlus;
      std::vector<std::pair<double, long long>> kMinus;
      kPlus.reserve(64);
      kMinus.reserve(64);
      for (long long i = 0; i < nGen; ++i) {
        if (genStatus[i] != 1) continue;

        const long long reco = genMatchIndex[i];
        if (reco < 0 || reco >= nReco) continue;
        if (genMatchAngle[i] >= matchAngleMax) continue;
        if (recoGoodTrack[reco] != 1) continue;
        const TrackKinematics t{recoPx[reco], recoPy[reco], recoPz[reco]};
        if (!passAcceptance(t)) continue;

        if (genID[i] == 321 && recoCharge[reco] > 0) kPlus.push_back(i);
        if (genID[i] == -321 && recoCharge[reco] < 0) kMinus.push_back({genE[i], i});
      }
      std::sort(kMinus.begin(), kMinus.end());

      for (long long iPhi = 0; iPhi < nGen; ++iPhi) {
        if (genID[iPhi] != 333) continue;
        totalCandidates++;

        long long bestKPlus = -1;
        long long bestKMinus = -1;
        double bestAbsDeltaE = std::numeric_limits<double>::infinity();
        double bestAngle = std::numeric_limits<double>::infinity();

        for (long long i : kPlus) {
          const double targetE = genE[iPhi] - genE[i];
          auto beginIt =
              std::lower_bound(kMinus.begin(), kMinus.end(), std::make_pair(targetE - kTruthEnergyMatchMax, -1LL));
          auto endIt =
              std::upper_bound(kMinus.begin(), kMinus.end(), std::make_pair(targetE + kTruthEnergyMatchMax, std::numeric_limits<long long>::max()));
          for (auto it = beginIt; it != endIt; ++it) {
            const long long j = it->second;
            const double absDeltaE = std::fabs(genE[i] + genE[j] - genE[iPhi]);
            if (absDeltaE >= bestAbsDeltaE) continue;
            const double sx = genPx[i] + genPx[j];
            const double sy = genPy[i] + genPy[j];
            const double sz = genPz[i] + genPz[j];
            const double angle = angleBetween(sx, sy, sz, genPx[iPhi], genPy[iPhi], genPz[iPhi]);
            bestAbsDeltaE = absDeltaE;
            bestAngle = angle;
            bestKPlus = i;
            bestKMinus = j;
          }
        }

        if (bestKPlus < 0 || bestKMinus < 0) continue;
        if (bestAbsDeltaE >= kTruthEnergyMatchMax) continue;
        if (bestAngle >= kTruthMomentumAngleMax) continue;

        const long long reco1 = genMatchIndex[bestKPlus];
        const long long reco2 = genMatchIndex[bestKMinus];
        if (reco1 < 0 || reco2 < 0 || reco1 >= nReco || reco2 >= nReco) continue;
        passValidRecoID++;

        if (genMatchAngle[bestKPlus] >= matchAngleMax || genMatchAngle[bestKMinus] >= matchAngleMax) continue;
        passMatchAngle++;

        if (recoGoodTrack[reco1] != 1 || recoGoodTrack[reco2] != 1) continue;
        passGoodTrack++;

        const TrackKinematics assumedKaon{recoPx[reco1], recoPy[reco1], recoPz[reco1]};
        const TrackKinematics assumedPion{recoPx[reco2], recoPy[reco2], recoPz[reco2]};
        if (!passAcceptance(assumedKaon) || !passAcceptance(assumedPion)) continue;
        passAcceptanceBoth++;

        if (recoCharge[reco1] * recoCharge[reco2] >= 0) continue;
        passOppositeCharge++;

        const double mass = buildMass(assumedKaon, assumedPion);
        hMassAccepted.Fill(mass);

        if (recoPIDKaon[reco1] < kKaonTagThreshold) continue;
        passKaonTag++;
        hMassKaonTag.Fill(mass);

        if (recoPIDPion[reco2] < kPionTagThreshold) continue;
        passKaonPionTag++;
        hMassKaonPionTag.Fill(mass);
      }
    }
  }

  TFile outputFile(outputFileName.c_str(), "RECREATE");
  hMassKaonTag.Write();
  hMassKaonPionTag.Write();
  hMassAccepted.Write();

  TNamed selection("SelectionSummary",
                   Form("Phi-as-KStar wrong-treatment study: use Phi* branches when available, otherwise for each "
                        "truth phi choose the final-state K+K- pair minimizing |E(KK)-E(phi)|, require "
                        "|E(KK)-E(phi)|<%.3f and angle(KK,phi)<%.3f, then require GenMatchAngle<%.4f on both legs; "
                        "both RecoGoodTrack==1, both 0.15<=|cos(theta)|<=0.675, opposite charge, treat daughter1 as "
                        "kaon and daughter2 as pion in mass build, require RecoPIDKaon>=2 on daughter1, and "
                        "optionally RecoPIDPion>=2 on daughter2, hist range %.3f-%.3f GeV",
                        kTruthEnergyMatchMax, kTruthMomentumAngleMax, matchAngleMax, massMin, massMax));
  selection.Write();

  TParameter<long long>("TotalPhiCandidates", totalCandidates).Write();
  TParameter<long long>("PassValidRecoID", passValidRecoID).Write();
  TParameter<long long>("PassMatchAngle", passMatchAngle).Write();
  TParameter<long long>("PassGoodTrack", passGoodTrack).Write();
  TParameter<long long>("PassAcceptanceBoth", passAcceptanceBoth).Write();
  TParameter<long long>("PassOppositeCharge", passOppositeCharge).Write();
  TParameter<long long>("PassKaonTag", passKaonTag).Write();
  TParameter<long long>("PassKaonPionTag", passKaonPionTag).Write();
  TParameter<double>("MassMin", massMin).Write();
  TParameter<double>("MassMax", massMax).Write();
  TParameter<double>("MatchAngleMax", matchAngleMax).Write();
  outputFile.Close();

  std::cout << "Wrote " << outputFileName << std::endl;
  std::cout << "  Total phi candidates:       " << totalCandidates << std::endl;
  std::cout << "  Pass valid reco IDs:        " << passValidRecoID << std::endl;
  std::cout << "  Pass match angles:          " << passMatchAngle << std::endl;
  std::cout << "  Pass good tracks:           " << passGoodTrack << std::endl;
  std::cout << "  Pass acceptance:            " << passAcceptanceBoth << std::endl;
  std::cout << "  Pass opposite charge:       " << passOppositeCharge << std::endl;
  std::cout << "  Pass kaon tag:              " << passKaonTag << std::endl;
  std::cout << "  Pass kaon+pion tag:         " << passKaonPionTag << std::endl;
  std::cout << "  Accepted entries in window: " << static_cast<long long>(hMassAccepted.GetEntries()) << std::endl;
  std::cout << "  Kaon-tag entries in window: " << static_cast<long long>(hMassKaonTag.GetEntries()) << std::endl;
  std::cout << "  K+pi-tag entries in window: " << static_cast<long long>(hMassKaonPionTag.GetEntries()) << std::endl;
  return 0;
}
