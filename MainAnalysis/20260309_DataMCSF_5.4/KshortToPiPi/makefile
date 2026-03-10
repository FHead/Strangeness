ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

default: ExecuteMakeKShortSignalOnlyHistograms ExecuteFitKShortSignalOnlyShapes ExecuteMakeKShortSBHistograms ExecuteFitKShortSB

ExecuteMakeKShortSignalOnlyHistograms: MakeKShortSignalOnlyHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKShortSignalOnlyHistograms.cpp \
		-o ExecuteMakeKShortSignalOnlyHistograms \
		$(ROOTLIBS)

ExecuteFitKShortSignalOnlyShapes: FitKShortSignalOnlyShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKShortSignalOnlyShapes.cpp \
		-o ExecuteFitKShortSignalOnlyShapes \
		$(ROOTLIBS)

ExecuteMakeKShortSBHistograms: MakeKShortSBHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKShortSBHistograms.cpp \
		-o ExecuteMakeKShortSBHistograms \
		$(ROOTLIBS)

ExecuteFitKShortSB: FitKShortSB.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKShortSB.cpp \
		-o ExecuteFitKShortSB \
		$(ROOTLIBS)
