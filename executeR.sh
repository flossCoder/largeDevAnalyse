# executeR.sh
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

# usage: bash executeR.sh
#  "directory of the /large-deviation-properties-of-the-stochastic-block-model/analysis/main.R-file"
#  "directory and filename of the settings script"

# run the given arguments
R -e "arg = '$2'; source('$1large-deviation-properties-of-the-stochastic-block-model/analysis/main.R')"
