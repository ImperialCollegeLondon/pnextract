#!/bin/bash

if [  -n "$PoreScaleDir" ]; then
	echo "Info: PoreScaleDir(=$PoreScaleDir) is NOT reset from "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
else

	export PoreScaleDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
	export myBinDir=$PoreScaleDir/bin


	export PATH=$PATH:$PoreScaleDir/scripts
	export PATH=$PATH:$PoreScaleDir/scripts/misc
	export PATH=$PATH:$PoreScaleDir/scripts/twoPhase
	export PATH=$PATH:$PoreScaleDir/scripts/networkModel
	export PATH=$PATH:$PoreScaleDir/scripts/singlePhase
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$myBinDir
	export PATH=$PATH:$myBinDir


	if ! [ -d $myBinDir ]; then
		 mkdir -p $myBinDir
	fi
fi

# to change based on the version of poreFoam, used by AllRun* scripts
# single-phase codes and preprocessing codes for two-phase flow solver
export poreFoamBashrc="$PoreScaleDir/poreFoam/bashrc"

# two-phase flow solver
export poreFoamBashrcExtend="$PoreScaleDir/interfoamccf/bashrc"

