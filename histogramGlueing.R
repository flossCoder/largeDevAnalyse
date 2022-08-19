# histogramGlueing.R
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

# some required constants
USEMETROPOLISALGORITHM <- 1
USEWANGLANDAUALGORITHM <- 2

#' Do the histogram glueing and determine the resulting histogram.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param simpleSampling The histogram obtained by simple sampling.
#' @param metropolisHistograms A list containing all histograms obtained by the Metropolis algorithm.
#' @param dataMetropolis A list containing meta informations fot the histograms obtained by the
#'                       Metropolis algorithm.
#' @param wangLandauHistograms A list containing all histograms obtained by the
#'                             Wang-Landau algorithm.
#' @param dataWangLandau A list containing all lower- and upper bounds of the given
#'                       histograms obtained by the Wang-Landau algorithm.
#' @param precision The number of bits for the Rmpfr data-types.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return The resulting histogram after adding the given histograms.
histogramGlueing <- function(directory, numberOfVertices, simpleSampling,
                             metropolisHistograms, dataMetropolis,
                             wangLandauHistograms, dataWangLandau,
                             precision, outputWithLog) {
  # remove all glueing files containing the z values, if they exist
  outMetropolis <- paste(directory, "/z_metropolis_", numberOfVertices, ".dat", sep = "")
  if (file.exists(outMetropolis)) {
    file.remove(outMetropolis)
  }
  outWL <- paste(directory, "/z_wl_", numberOfVertices, ".dat", sep = "")
  if (file.exists(outWL)) {
    file.remove(outWL)
  }
  
  # set up the resulting histogram
  resultingHistogram <- matrix(NA, nrow=numberOfVertices, ncol=3)
  for (i in 1:numberOfVertices) {
    resultingHistogram[i, 1] <- i
    resultingHistogram[i, 2] <- 0
    resultingHistogram[i, 3] <- Inf
  }
  
  # convert to mpfrMatrix
  resultingHistogram <- new("mpfrMatrix", mpfr(resultingHistogram, precision),
                            Dim = c(dim(resultingHistogram)[1], dim(resultingHistogram)[2]))
  
  # insert simple sampling into the resulting matrix
  for (i in 1:dim(simpleSampling)[1]) {
    resultingHistogram[as.numeric(simpleSampling[i, 1]), 2] <- simpleSampling[i, 2]
    resultingHistogram[as.numeric(simpleSampling[i, 1]), 3] <- simpleSampling[i, 3]
  }
  
  # as long as there are histograms that have not been glued jet, glue the best fitting histogram
  # to the resulting one
  glueCounter <- 0
  if (!is.null(dim(dataMetropolis))) {
    notGluedMetropolis <- matrix(TRUE, nrow=dim(dataMetropolis)[1])
  }
  if (!is.null(dim(dataWangLandau))) {
    notGluedWangLandau <- matrix(TRUE, nrow=dim(dataWangLandau)[1])
  }
  
  lengthWangLandau <- 0
  if (!is.null(dim(dataWangLandau))) {
    lengthWangLandau <- sum(!is.na(dataWangLandau[, 1]), na.rm = TRUE)
  }
  lengthMetropolis <- 0
  if (!is.null(dim(dataMetropolis))) {
    lengthMetropolis <- sum(!is.na(dataMetropolis[, 1]), na.rm = TRUE)
  }
  
  while ((glueCounter < sum(lengthWangLandau, lengthMetropolis, na.rm = TRUE))) {
    glueCounter <- glueCounter + 1
    # determine the best fitting histogram with the resulting one
    currentAlgorithm <- NA # 0: Metropolis, 1: Wang-Landau
    currentIndex <- NA
    currentOverlap <- c()
    
    # first check all histograms determined by the Metropolis algorithm
    if (!is.null(dim(dataMetropolis))) {
      for (i in which(notGluedMetropolis & !is.na(dataMetropolis[, 1]), arr.ind=TRUE)[, 1]) {
        # calculate the overlap for the current histogram
        auxOverlap <- calculateOverlap(resultingHistogram, metropolisHistograms[[i]])
        if ((!is.null(dim(auxOverlap)[1]) && is.null(dim(currentOverlap)[1])) ||
            (!is.null(dim(auxOverlap)[1]) && !is.null(dim(currentOverlap)[1]) && 
             (dim(auxOverlap)[1] >= dim(currentOverlap)[1]))) {
          # exchange the currentOverlap
          currentAlgorithm <- USEMETROPOLISALGORITHM
          currentIndex <- i
          currentOverlap <- auxOverlap
        }
      }
    }
    
    # second check all histograms determined by the Wang-Landau algorithm
    if (!is.null(dim(dataWangLandau))) {
      for (i in which(notGluedWangLandau & !is.na(dataWangLandau[, 1]), arr.ind=TRUE)[, 1]) {
        # calculate the overlap for the current histogram
        auxOverlap <- calculateOverlap(resultingHistogram, wangLandauHistograms[[i]])
        if ((!is.null(dim(auxOverlap)[1]) && is.null(dim(currentOverlap)[1])) ||
            (!is.null(dim(auxOverlap)[1]) && !is.null(dim(currentOverlap)[1]) && 
             (dim(auxOverlap)[1] >= dim(currentOverlap)[1]))) {
          # exchange the currentOverlap
          currentAlgorithm <- USEWANGLANDAUALGORITHM
          currentIndex <- i
          currentOverlap <- auxOverlap
        }
      }
    }
    
    if (is.na(currentAlgorithm)) {
      # there are histograms left, that can't be glued because of a missing overlap
      warning("failed to glue the histograms")
      glueCounter <- Inf
    } else {
      if (currentAlgorithm == USEMETROPOLISALGORITHM) {
        # glue the best histogram obtained by the Metropolis algorithm
        resultingHistogram <- glueHistogramToResult(directory, numberOfVertices, resultingHistogram,
                                                    metropolisHistograms[[currentIndex]],
                                                    dataMetropolis[currentIndex,],
                                                    currentOverlap, USEMETROPOLISALGORITHM,
                                                    outputWithLog)
        # mark the current histogram as glued
        notGluedMetropolis[currentIndex] <- FALSE
      } else if (currentAlgorithm == USEWANGLANDAUALGORITHM) {
        # glue the best histogram obtained by the Wang-Landau algorithm
        resultingHistogram <- glueHistogramToResult(directory, numberOfVertices, resultingHistogram,
                                                    wangLandauHistograms[[currentIndex]],
                                                    dataWangLandau[currentIndex,],
                                                    currentOverlap, USEWANGLANDAUALGORITHM,
                                                    outputWithLog)
        # mark the current histogram as glued
        notGluedWangLandau[currentIndex] <- FALSE
      } else {
        warning(paste("invallid current algorithm =", algorithm))
      }
    }
  }
  
  # save results
  saveHistogram(resultingHistogram, directory, paste("hist_", numberOfVertices, "_result", sep = ""),
                outputWithLog)
  # determine the average value of the histogram and save it to file
  res <- averageHistogram(resultingHistogram)
  cat(sprintf("%s %s\n", toString(format(res[1], drop0 = TRUE, digits = 22)),
              toString(format(res[2], drop0 = TRUE, digits = 22))),
      file = paste(directory, "/avg_", numberOfVertices, ".dat", sep = ""), append = FALSE)
  
  return(resultingHistogram)
}

#' Glue the given histogram to the resulting histogram, save the z factor plus the shifted histogram.
#' Finally the shifted histogram is added to the resulting histogram.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param resultingHistogram The resulting histogram.
#' @param histogram The histogram which will be glued to the resulting histogram.
#' @param data Providing information for saving the histogram under an correct filename.
#' @param overlap A matrix containing the overlap of the two histograms.
#' @param algorithm The number of the algorithm (e. g. 1 = Metropolis).
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' 
#' @return The resulting histogram after adding the given histogram.
glueHistogramToResult <- function(directory, numberOfVertices, resultingHistogram, histogram, data,
                                  overlap, algorithm, outputWithLog) {
  z <- fitHistograms(resultingHistogram, histogram, overlap)
  histogram[, 2] <- histogram[, 2] * z[1]
  histogram[, 3] <- abs(histogram[, 3] * z[1])# + abs(histogram[, 2] * z[2])
  # add the histogram to the resulting histogram
  resultingHistogram <- mergeHistogram(resultingHistogram, histogram)
  # save the histogram
  outname <- NA
  if (algorithm == USEMETROPOLISALGORITHM) {
    outname <- paste("hist_is_", numberOfVertices, "_", data[2], "_", data[1], "_shifted", sep = "")
  } else if (algorithm == USEWANGLANDAUALGORITHM) {
    outname <- paste("density_wl_", numberOfVertices, "_", data[1], "_", data[2], "_shifted", sep = "")
  }
  
  saveHistogram(histogram, directory, outname, outputWithLog)
  
  # save z
  if (algorithm == USEMETROPOLISALGORITHM) {
    # out format: temperature z Delta-z
    out <- paste(directory, "/z_metropolis_", numberOfVertices, ".dat", sep = "")
    if (outputWithLog) {
      cat(sprintf("%s %s %s\n", toString(format(data[1], drop0 = TRUE)),
                  toString(format(log10(z[1]), digits = 22)),
                  toString(format(log10(z[2]), digits = 22))), file = out, append = TRUE)
    } else {
      cat(sprintf("%s %s %s\n", toString(format(data[1], drop0 = TRUE)),
                  toString(format(z[1], digits = 22)),
                  toString(format(z[2], digits = 22))), file = out, append = TRUE)
    }
  } else if (algorithm == USEWANGLANDAUALGORITHM) {
    # out format: min-interval max-interval z Delta-z
    out <- paste(directory, "/z_wl_", numberOfVertices, ".dat", sep = "")
    if (outputWithLog) {
      cat(sprintf("%s %s %s %s\n", toString(format(data[1], drop0 = TRUE)),
                  toString(format(data[2], drop0 = TRUE)),
                  toString(format(log10(z[1]), digits = 22)),
                  toString(format(log10(z[2]), digits = 22))), file = out, append = TRUE)
    } else {
      cat(sprintf("%s %s %s %s\n", toString(format(data[1], drop0 = TRUE)),
                  toString(format(data[2], drop0 = TRUE)),
                  toString(format(z[1], digits = 22)),
                  toString(format(z[2], digits = 22))), file = out, append = TRUE)
    }
  }
  return(resultingHistogram)
}

#' Glue the two histograms together.
#' 
#' @param histogram1 The basis histogram.
#' @param histogram2 This histogram will be fitted to histogram1.
#' @param overlap A matrix containing the overlap of the two histograms.
#' 
#' @return The normalisation factor z.
fitHistograms <- function(histogram1, histogram2, overlap) {
  result <- 0
  
  # calculate the average of the z value
  for (i in 1:dim(overlap)[1]) {
    if (histogram2[overlap[i, 2], 2] != 0) {
      # security check to prevent from division through zero
      result <- result + histogram1[overlap[i, 1], 2] / histogram2[overlap[i, 2], 2]
    }
  }
  avg <- result / dim(overlap)[1]
  
  result2 <- 0
  # calculate the variance
  for (i in 1:dim(overlap)[1]) {
    if (histogram2[overlap[i, 2], 2] != 0) {
      # security check to prevent from division through zero
      result2 <- result2 + (histogram1[overlap[i, 1], 2] / histogram2[overlap[i, 2], 2] - avg)^2
    }
  }
  var <- result2 / dim(overlap)[1]
  
  # calculate the error
  err <- sqrt(var / (dim(overlap)[1] - 1))
  
  return(c(avg, err))
}

#' Calculate the indices of the overlapping region.
#' 
#' @param resultingHistogram The resulting histogram.
#' @param histogram The histogram which will be glued to the resulting histogram.
#' 
#' @return A matrix containing the indices of the overlapping region.
calculateOverlap <- function(resultingHistogram, histogram) {
  # calculate the intersection of the two lists
  intersection <- intersect(as.numeric(resultingHistogram[, 1] * (resultingHistogram[, 2] != 0)),
                            as.numeric(histogram[, 1]))
  
  # explanation of the following code:
  #   y = is.element(x, intersection) calculates a boolean list containing
  #     TRUE if the corresponding element of x is in intersection
  #     FALSE otherwise
  boolResultingHistogram <- is.element(as.numeric(resultingHistogram[, 1]), intersection)
  boolHistogram <- is.element(as.numeric(histogram[, 1]), intersection)
  
  # calculate the indices of the two histograms (if possible)
  if (any(boolResultingHistogram) && any(boolHistogram)) {
    # which(y, arr.ind = TRUE) calculates a list of all indices where y is TRUE
    indexResultingHistogram <- which(boolResultingHistogram, arr.ind = TRUE)
    indexHistogram <- which(boolHistogram, arr.ind = TRUE)
    # return the result
    return(cbind(indexResultingHistogram, indexHistogram))
  } else {
    # there is no intersection
    return(NA)
  }
}

#' Insert the given histogram into the one which should be created as a final result.
#' 
#' @param resultingHistogram The resulting histogram.
#' @param histogram The given histogram.
#' 
#' @return The resulting histogram after processing the given histogram.
mergeHistogram <- function(resultingHistogram, histogram) {
  for (i in 1:dim(histogram)[1]) {
    if ((asNumeric(resultingHistogram[asNumeric(histogram[i, 1]), 3]) != 0) &&
        (asNumeric(histogram[i, 3]) != 0)) {
      # use the error for determining the best fitting histogram
      if (asNumeric(resultingHistogram[asNumeric(histogram[i, 1]), 3]) > asNumeric(histogram[i, 3])) {
        resultingHistogram[asNumeric(histogram[i, 1]), 2] <- histogram[i, 2]
        resultingHistogram[asNumeric(histogram[i, 1]), 3] <- histogram[i, 3]
        
      }
    } else if (calculateNeighbourSlope(i, histogram) <
               calculateNeighbourSlope(asNumeric(histogram[i, 1]), resultingHistogram) &&
              (asNumeric(histogram[i, 2]) != 0)) {
      # can´t use error => use slope instead
      resultingHistogram[asNumeric(histogram[i, 1]), 2] <- histogram[i, 2]
      resultingHistogram[asNumeric(histogram[i, 1]), 3] <- histogram[i, 3]
    }
  }
  return(resultingHistogram)
}

#' Calculate the slope of the given histogram in the index point.
#' In case the slope (or the index point) is zero the function returns 1000.
#' 
#' @param index The index where the slope shall be determinded.
#' @param histogram The histogram to evaluate.
#' 
#' @return The slope or 1000.
calculateNeighbourSlope <- function(index, histogram) {
  dX <- 0
  dY <- 0
  if (asNumeric(histogram[index, 2]) == 0) {
    return(1000)
  }
  if (index > 1) {
    dX <- abs(asNumeric(histogram[index, 1]) - asNumeric(histogram[(index - 1), 1]))
    dY <- abs(asNumeric(histogram[index, 2]) - asNumeric(histogram[(index - 1), 2]))
  }
  if (index < dim(histogram)[1]) {
    dX <- dX + abs(asNumeric(histogram[(index + 1), 1]) - asNumeric(histogram[index, 1]))
    dY <- dY + abs(asNumeric(histogram[(index + 1), 2]) - asNumeric(histogram[index, 2]))
  }
  res <- 1000
  if (dX != 0 && dY != 0) {
    res <- (dY / dX)
  }
  return(res)
}
