# doAnalysis.sh
# Copyright (C) 2016, 2017 flossCoder
#
# This file is part of largeDevAnalyse.
#
# largeDevAnalyse is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# largeDevAnalyse is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

#!/bin/bash

# usage: bash doAnalysis.sh
#  "directory of the /large-deviation-properties-of-the-stochastic-block-model/analysis/main.R-file"
#  "filename of the root directory of the simulation results"

# for all subdirectories, run executeR, if settings.R exists
for dir in $2/*
do
	if [ -e "$dir/settings.R" ]
	then
		bash $1large-deviation-properties-of-the-stochastic-block-model/analysis/executeR.sh $1 $dir/settings.R
	fi
done
