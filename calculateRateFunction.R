# calculateRateFunction.R
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

#' Calculate the empirical rate function for the given histogram.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param histogram The given histogram.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return The empirical rate function.
calculateRateFunction <- function(directory, numberOfVertices, histogram, outputWithLog) {
  # calculate the rate function
  histogram[, 1] <- histogram[, 1] / numberOfVertices
  histogram[, 2] <- -1 / numberOfVertices * log10(histogram[, 2])
  histogram[, 3] <- -(histogram[, 2] != 0) / (numberOfVertices * histogram[, 2]) * histogram[, 3]
  
  # save the rate function
  saveHistogram(histogram, directory,
                paste("hist_", numberOfVertices, "_rate-function", sep = ""),
                FALSE) # never take the log of the rate function
  
  return(histogram)
}
