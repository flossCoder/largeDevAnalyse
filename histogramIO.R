# histogramIO.R
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

#' Read the data from a file as matrix.
#' 
#' @param directory Where the file is saved.
#' @param filename Name of the file (without file extension).
#' @param precision The number of bits for the Rmpfr data-types.
#' 
#' @return A matrix containing the histogram data.
readHistogram <- function(directory, filename, precision) {
  # read data as data.frame
  frame <- read.table(paste(directory, "/", filename, ".dat", sep=""), sep = " ", header = FALSE)
  # convert data to matrix
  data <- as.matrix(frame)
  if (dim(data)[2] == 2) {
    # if there are just two columns in data, add a third one with zeros
    data <- cbind(data, matrix(0, ncol = 1, nrow = dim(data)[1]))
  }
  # return as mpfrMatrix
  return(new("mpfrMatrix", mpfr(data, precision), Dim = c(dim(data)[1], dim(data)[2])))
}

#' Save the given histogram to a dat file.
#' 
#' @param histogram The histogram, which has to be saved to file.
#' @param directory Where the file is saved.
#' @param filename Name of the file (without file extension).
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
saveHistogram <- function(histogram, directory, filename, outputWithLog) {
  out = paste(directory, "/", filename, ".dat", sep="")
  if (file.exists(out)){
    file.remove(out)
  }
  for (i in 1:dim(histogram)[1]) {
    if (outputWithLog) {
      value <- log10(histogram[i, 2])
      error <- histogram[i, 3]
      if (error == -Inf) {
        # set a very small number as the error
        error <- log10(mpfr(10^(-323), precision))
      }
      if (value != -Inf) {
        # Just print the result, when it is not -Inf.
        cat(sprintf("%s %s %s\n", toString(format(histogram[i, 1], digits = 3)),
                    toString(format(value, digits = 22)), toString(format(error, digits = 22))),
            file = out, append = TRUE)
      }
    } else {
      if (histogram[i, 2] != 0) {
        cat(sprintf("%s %s %s\n", toString(format(histogram[i, 1], digits = 3)),
                    toString(format(histogram[i, 2], digits = 22)),
                    toString(format(histogram[i, 3], digits = 22))), file = out, append = TRUE)
      }
    }
  }
}
