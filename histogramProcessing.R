# histogramProcessing.R
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

#' Process a histogram obtained by simple sampling:
#'  1.) open the data
#'  2.) calculate the PDF
#' 
#' @param directory Where the file is saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param numberOfGraphs The number of guessed graphs during simple sampling.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return The histogram representing the unbiased PDF of the given raw data.
processSimpleSamplingHistogram <- function(directory, numberOfVertices, numberOfGraphs, precision,
                                           outputWithLog) {
  # open importance sampling histogram
  histogram <- readHistogram(directory,
                             paste("hist_ss_", numberOfVertices, "_", numberOfGraphs, sep = ""),
                             precision)
  # calculate PDF of the histogram
  histogram <- calculatePDFForHistogram(histogram)
  # save PDF histogram
  saveHistogram(histogram, directory,
                paste("hist_ss_", numberOfVertices, "_", numberOfGraphs, "_PDF", sep = ""),
                outputWithLog)
  # determine the average value of the histogram and save it to file
  res <- averageHistogram(histogram)
  cat(sprintf("%s %s\n", toString(format(res[1], drop0 = TRUE, digits = 22)),
                         toString(format(res[2], drop0 = TRUE, digits = 22))),
      file = paste(directory, "/avg_ss_", numberOfVertices, "_", numberOfGraphs, ".dat", sep = ""),
      append = FALSE)
  return(histogram)
}

#' Process the histograms obtained by the Metropolis algorithm and return the result as
#' a big list for further processing.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param data A list containing all temperatures (correctly ordered) and the number of sweeps
#'             for each temperature.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return A list containing all histograms to the given temperatures.
processMetropolisAlgorithmListOfHistogram <- function(directory, numberOfVertices, 
                                                      data, precision, outputWithLog) {
  result <- list()
  out <- paste(directory, "/avg_metropolis_", numberOfVertices, ".dat", sep = "")
  if (file.exists(out)) {
    file.remove(out)
  }
  for (i in 1:dim(data)[1]) {
    if (!is.na(data[i, 1]) && !is.na(data[i, 2])) {
      result[[length(result) + 1]] <- processMetropolisAlgorithmHistogram(directory,
                                                                          numberOfVertices,
                                                                          data[i, 2],
                                                                          data[i, 1],
                                                                          precision,
                                                                          outputWithLog)
      # determine the average value of the histogram and save it to file
      res <- averageHistogram(result[[length(result)]])
      cat(sprintf("%s %s %s\n", toString(data[i, 1]),
                                toString(format(res[1], drop0 = TRUE, digits = 22)),
                                toString(format(res[2], drop0 = TRUE, digits = 22))),
          file = out, append = TRUE)
    }
  }
  return(result)
}

#' Process a histogram obtained by the Metropolis algorithm:
#'  1.) open the data
#'  2.) calculate the PDF
#'  3.) unbias the PDF
#' 
#' @param directory Where the file is saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param sweeps How many sweeps should be recorded.
#' @param temperature The artificial temperature.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return The histogram representing the unbiased PDF of the given raw data.
processMetropolisAlgorithmHistogram <- function(directory, numberOfVertices, sweeps,
                                                temperature, precision, outputWithLog) {
  # open importance sampling histogram
  histogram <- readHistogram(directory,
                             paste("hist_is_", numberOfVertices, "_", sweeps, "_",
                                   temperature, sep = ""), precision)
  
  # calculate PDF of the histogram
  histogram <- calculatePDFForHistogram(histogram)
  # save PDF histogram
  saveHistogram(histogram, directory,
                paste("hist_is_", numberOfVertices, "_", sweeps, "_", temperature,
                      "_PDF", sep = ""), outputWithLog)
  # unbias PDF histogram
  histogram <- unbiasHistogram(histogram, temperature)
  # save unbiased histogram
  saveHistogram(histogram, directory,
                paste("hist_is_", numberOfVertices, "_", sweeps, "_", temperature,
                      "_unbiased", sep = ""), outputWithLog)
  return(histogram)
}

#' Process the given histograms obtained by the Wang-Landau algorithm.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param data A list containing all lower- and upper histogram bounds.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return A list containing all histograms given in data.
processWLHistograms <- function(directory, numberOfVertices, data, precision, outputWithLog) {
  result <- list()
  out <- paste(directory, "/avg_wl_", numberOfVertices, ".dat", sep = "")
  if (file.exists(out)) {
    file.remove(out)
  }
  for (i in 1:dim(data)[1]) {
    if (!is.na(data[i, 1]) && !is.na(data[i, 2])) {
      result[[length(result) + 1]] <- processWLHistogram(directory, numberOfVertices,
                                                         data[i, 1], data[i, 2], precision,
                                                         outputWithLog)
      res <- averageHistogram(result[[length(result)]])
      cat(sprintf("%s %s %s %s\n", toString(data[i, 1]), toString(data[i, 2]),
                                   toString(format(res[1], drop0 = TRUE, digits = 22)),
                                   toString(format(res[2], drop0 = TRUE, digits = 22))),
          file = out, append = TRUE)
    }
  }
  return(result)
}

#' Process a histogram obtained by the Wang-Landau algorithm:
#'  1.) open the data
#'  2.) calculate the PDF
#' 
#' @param directory Where the file is saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param lowerBound The minimum interval value of the obtained histogram.
#' @param upperBound The maximum interval value of the obtained histogram.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of
#'        the output or not.
#' 
#' @return The histogram representing the unbiased PDF of the given raw data.
processWLHistogram <- function(directory, numberOfVertices, lowerBound, upperBound,
                               precision, outputWithLog) {
  # open the histogram
  histogram <- readHistogram(directory,
                             paste("density_", numberOfVertices, "_", lowerBound, "_",
                                   upperBound, sep = ""), precision)
  # calculate PDF of the histogram
  histogram <- calculatePDFForHistogram(histogram)
  # save PDF histogram
  saveHistogram(histogram, directory,
                paste("density_wl_", numberOfVertices, "_", lowerBound, "_", upperBound,
                      "_PDF", sep = ""), outputWithLog)
  return(histogram)
}

#' Remove the bias from the histogram.
#' 
#' @param histogram The pdf of the histogram with bins and errors.
#' @param temperature The artificial temperature of the system.
#' 
#' @return The unbiased histogram.
unbiasHistogram <- function(histogram, temperature) {
  histogram[, 2] <- histogram[, 2] * exp(histogram[, 1] / temperature)
  histogram[, 3] <- histogram[, 3] * exp(histogram[, 1] / temperature)
  return(histogram)
}

#' Calculate the pdf for an histogram.
#' 
#' @param histogram The given histogram with integer bins and errors.
#' 
#' @return The pdf of the histogram.
calculatePDFForHistogram <- function(histogram) {
  # divide histogram through the number of counts, if the number of counts != 0
  numberOfSamples <- sum(histogram[, 2], na.rm = TRUE)
  if (numberOfSamples != 0) {
    histogram[, 2] <- histogram[, 2] / numberOfSamples
    histogram[, 3] <- histogram[, 3] / numberOfSamples
  }
  return(histogram)
}
