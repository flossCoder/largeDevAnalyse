# doAnalysis.R
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

#' Do the whole analysis including histogram glueing, saving the results, ....
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param numberOfGraphs The number of guessed graphs during simple sampling.
#' @param dataMetropolis A matrix containing all positive and the number of sweeps for each temperature.
#' @param dataWangLandau A list containing all lower- and upper bounds of the given
#'                       histograms determined by the Wang-Landau algorithm.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' @param noGlueing Do not glue the histograms.
doAnalysis <- function(directory, numberOfVertices, numberOfGraphs,
                       dataMetropolis, dataWangLandau, precision, outputWithLog, noGlueing) {
  # process simple sampling
  simpleSampling <- processSimpleSamplingHistogram(directory, numberOfVertices, numberOfGraphs,
                                                   precision, outputWithLog)
  
  metropolisHistograms <- NA
  if (!is.null(dim(dataMetropolis))) {
    metropolisHistograms <- processMetropolisAlgorithmListOfHistogram(directory, numberOfVertices,
                                                                            dataMetropolis, precision,
                                                                            outputWithLog)
  }
  
  wangLandauHistograms <- NA
  # process the histograms determined by the Wang-Landau algorithm
  if (!is.null(dim(dataWangLandau))) {
    wangLandauHistograms <- processWLHistograms(directory, numberOfVertices, dataWangLandau,
                                                precision, outputWithLog)
  }
  
  if (!noGlueing) {
    # glue the histograms together
    resultingHistogram <- histogramGlueing(directory, numberOfVertices, simpleSampling,
                                           metropolisHistograms, dataMetropolis,
                                           wangLandauHistograms, dataWangLandau,
                                           precision, outputWithLog)
    # calculate the empirical rate function for the resulting histogram
    rateFunction <- calculateRateFunction(directory, numberOfVertices, resultingHistogram,
                                          outputWithLog)
  }
}
