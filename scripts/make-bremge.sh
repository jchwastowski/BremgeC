#!/bin/bash
#
echo "    Found: HepMC3  version: "`HepMC3-config --version`" at "`root-config --prefix`
echo "              GSL  version: "`gsl-config --version`   " at "`gsl-config --prefix`
echo "             ROOT  version: "`root-config --version`  " at "`root-config --prefix`
#################################
hepmcinc=`HepMC3-config --includedir`
rootinc=`root-config --incdir`
rootlib=`root-config --libdir`
rootlibs=`root-config --libs`
DSMEAR=""
DPHOTO=""
BIN=`pwd`"/bin/"
SRC=`pwd`"/src/"
#
mkdir -p ${BIN}
#
let case=$1+1
echo "   case: "$case
case $case in
    1) echo "     regular HepMC3 output ";;
    2) echo "     regular HepMC3 output with beam divergency"
	    DSMEAR="-DBEAM_SPREAD_THETA";;
    3) echo "     only radiated photons on HepMC3 output "
            DPHOTO="-DLUMI_DIRECT_TEST";;
    4) echo "     only radiated photons with beam divergency efect on HepMC3 output "
 	    DSMEAR="-DBEAM_SPREAD_THETA"
	    DPHOTO="-DLUMI_DIRECT_TEST";;
    4) echo "      ERROR: unknown option: $case "       
       exit;;
esac
#
echo "     compiler options: $DSMEAR  $DPHOTO "
#
c++  -o ${BIN}bremgenC ${SRC}Main.C ${DSMEAR} -Wno-cpp ${DPHOTO} -I./ -I${hepmcinc} -I${rootinc} -L/usr/lib64 -lHepMC3 -lgsl -lgslcblas -L${rootlib} ${rootlibs}		  
#
