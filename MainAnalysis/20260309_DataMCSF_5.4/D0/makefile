ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

default: ExecuteMakeD0SignalOnlyHistograms ExecuteMakeD0WrongAssignmentHistograms ExecuteCombineD0Assignments ExecuteFitD0SignalOnlyShapes ExecuteMakeD0SBHistograms ExecuteFitD0SB

ExecuteMakeD0SignalOnlyHistograms: MakeD0SignalOnlyHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeD0SignalOnlyHistograms.cpp \
		-o ExecuteMakeD0SignalOnlyHistograms \
		$(ROOTLIBS)

ExecuteMakeD0WrongAssignmentHistograms: MakeD0WrongAssignmentHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeD0WrongAssignmentHistograms.cpp \
		-o ExecuteMakeD0WrongAssignmentHistograms \
		$(ROOTLIBS)

ExecuteCombineD0Assignments: CombineD0Assignments.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		CombineD0Assignments.cpp \
		-o ExecuteCombineD0Assignments \
		$(ROOTLIBS)

ExecuteFitD0SignalOnlyShapes: FitD0SignalOnlyShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitD0SignalOnlyShapes.cpp \
		-o ExecuteFitD0SignalOnlyShapes \
		$(ROOTLIBS)

ExecuteMakeD0SBHistograms: MakeD0SBHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeD0SBHistograms.cpp \
		-o ExecuteMakeD0SBHistograms \
		$(ROOTLIBS)

ExecuteFitD0SB: FitD0SB.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitD0SB.cpp \
		-o ExecuteFitD0SB \
		$(ROOTLIBS)
