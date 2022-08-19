# runParamList.sh
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

# usage: bash runParamList.sh
#  "root directory of the experiment"
#  "filename of the parameter list, which should be recalculated"
#  "directory of the /large-deviation-properties-of-the-stochastic-block-model/analysis/main.R-file"

# for all lines in the parameter list, run executeR, if settings.R exists

while read p; do
    echo $1$p/settings.R
    if [ -e "$1$p/settings.R" ]
	then
        bash $3large-deviation-properties-of-the-stochastic-block-model/analysis/executeR.sh $3 $1$p/settings.R
    fi
done <$1$2
