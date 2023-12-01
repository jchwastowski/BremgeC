#!/usr/bin/bash
mINPUT=$1".hepmc3"
mOUTPUT=$1".edm4hep.root"
mEvents=$2
#
echo " inputfile: $mINPUT outputfile: $mOUTPUT  number of events: $mEvents"
#
ddsim --compactFile $DETECTOR_PATH/epic_lumi_only.xml -I $mINPUT -O $mOUTPUT -N $mEvents --output.kernel 6 --output.inputStage 6
#
