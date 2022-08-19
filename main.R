# main.R
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

# Load packages.
require("Rmpfr")

# Set up R-file sources.

# Function from the ?source documentation.
sourceDir <- function(path, thisFile, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if (trace) {
      cat(file.path(path, nm))
    }
    if (file.path(path, nm) != file.path(path, thisFile)) {
      source(file.path(path, nm))
    }
    if (trace) {
      cat("\n")
    }
  }
}

sourceDir(dirname(sys.frame(1)$ofile), "main.R")

# for none scientific file processing => needed for file names
options(scipen = 999)

# set default values for some arguments (can be changed in the settings file)
precision <- 1000 # number of bits precision for the calculation
outputWithLog <- FALSE # do extra small values
algorithm <- 1 # 1: Metropolis; 2: Wang-Landau
dataMetropolis <- NA
dataWangLandau <- NA
xLabel <- "S"
yLabel <- "P(S)"
noGlueing <- FALSE # TRUE: no histogram glueing is perfomed
noPlotting <- FALSE # TRUE: no plotting is performed

# source given settings file
source(arg)

# do the whole analysis
doAnalysis(directory, numberOfVertices, numberOfGraphs, dataMetropolis, dataWangLandau,
           precision, outputWithLog, noGlueing)

# plot the whole stuff
if (!noPlotting) {
  doPlotting(directory, numberOfVertices, numberOfGraphs, dataMetropolis, dataWangLandau,
             xLabel, yLabel, outputWithLog, noGlueing)
}
