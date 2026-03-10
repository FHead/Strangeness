This directory holds systematic-variation outputs and notes for the `phi -> KK` efficiency study.

Current variation families:

- `SignalFunction/`: compare the nominal `GaussPlusRightTailCB` signal model against `TripleGaussian`.
- `FitRange/`: compare the nominal fit range `0.99-1.06` against the narrower `1.00-1.05`.
- `MatchingAngle/`: compare the nominal generator-reco matching requirement `0.01` against `0.025`.

Each concrete variation gets its own output directory with:

- `PhiSignalOnlyHistograms.root`
- `SignalOnlyFits/`
- `DataFits/`
- `manifest.txt`

Use `../RunPhiSystematics.sh` from the `PhiToKK` directory to populate these outputs.
