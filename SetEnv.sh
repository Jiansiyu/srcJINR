#!/bin/bash

SIMMACHINE=$(uname -n)

# change SIMPATH if you installed FairSoft not to the '/opt/fairsoft/install' directory
export SIMPATH=/nobackup1/segarrae/software/fairsoft
# change FAIRROOTPATH if you installed FairRoot not to the '/opt/fairroot/install' directory
export FAIRROOTPATH=/nobackup1/segarrae/software/fairroot

export ROOTSYS=$SIMPATH

export PATH=$SIMPATH/bin:$PATH
export LD_LIBRARY_PATH=$SIMPATH/lib:$SIMPATH/lib/root:$LD_LIBRARY_PATH

source geant4.sh

platform=$(root-config --arch)
echo SIMPATH is pointing to $SIMPATH
echo FAIRROOTPATH is pointing to $FAIRROOTPATH
