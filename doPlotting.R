# doPlotting.R
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

#' Do the whole plotting business in gnuplot.
#' 
#' @param directory Where the files are saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param numberOfGraphs The number of guessed graphs during simple sampling.
#' @param dataMetropolis A matrix containing all positive and the number of sweeps for each temperature.
#' @param dataWangLandau A list containing all lower- and upper bounds of the given
#'                       histograms determined by the Wang-Landau algorithm.
#' @param xLabel Label for the abscissa.
#' @param yLabel Label for the ordinate.
#' @param outputWithLog A boolean variable indicating whether to save the log of the output or not.
#' @param noGlueing Do not glue the histograms.
doPlotting <- function(directory, numberOfVertices, numberOfGraphs, dataMetropolis, dataWangLandau,
                       xLabel, yLabel, outputWithLog, noGlueing) {
  # basic set up for writing out the gnuplot file
  simpleSampling <- paste('"hist_ss_', numberOfVertices, '_', numberOfGraphs,
                          '_PDF.dat" using 1:2 title "simple sampling"', sep = "")
  listOfNamesMetropolis <- NA
  if (!is.null(dim(dataMetropolis))) {
    listOfNamesMetropolis <- prepareListOfNamesMetropolis(numberOfVertices, dataMetropolis)
  }
  listOfNamesWangLandau <- NA
  if (!is.null(dim(dataWangLandau))) {
    listOfNamesWangLandau <- prepareListOfNamesWangLandau(numberOfVertices, dataWangLandau)
  }
  
  # prepare gnuplot file
  filename <- paste(directory, "/plot.gpl", sep = "")
  # write first line of gnuplot script
  cat("#! /usr/bin/gnuplot\n\n", file = filename)
  cat('set pointsize 2\n\n', file = filename, append = TRUE)
  # set the output terminal to postscript enhanced color
  cat('set terminal postscript enhanced color font "Helvetica, 18" lw "2pt"\n\n', file = filename,
      append = TRUE)
  
  # set format of ordinate
  if (outputWithLog) {
    cat('set format y "10^{%.s}"\n\n', file = filename, append = TRUE)
  } else {
    cat('set format y "10^{%T}"\n\n', file = filename, append = TRUE)
  }
  
  # set ordinate axis logarithmic
  if (!outputWithLog) {
    cat("set logscale y\n\n", file = filename, append = TRUE)
  }
  
  # set x label
  cat("set xlabel '{/Helvetica-Italic ", xLabel, "}'\n\n", file = filename, append = TRUE)
  # set y label
  cat("set ylabel '{/Helvetica-Italic ", yLabel, "}'\n\n", file = filename, append = TRUE)
  # for compatibility reasons is the legend outside of the plot
  if (sum(dim(dataWangLandau)[1], dim(dataMetropolis)[1], na.rm = TRUE) > 20) {
    cat("unset key\n\n", file = filename, append = TRUE)
  } else {
    cat("set key outside\n\n", file = filename, append = TRUE)
  }
  
  # concatenate list of names
  listOfNames <- NA
  if (!is.null(dim(listOfNamesMetropolis)) && !is.null(dim(listOfNamesWangLandau))) {
    listOfNames <- rbind(listOfNamesMetropolis, listOfNamesWangLandau)
  } else if (!is.null(dim(listOfNamesMetropolis))) {
    listOfNames <- listOfNamesMetropolis
  } else if (!is.null(dim(listOfNamesWangLandau))) {
    listOfNames <- listOfNamesWangLandau
  } else {
    return()
  }
  
  # save instructions for plotting
  savePlottingInstruction(filename, numberOfVertices, simpleSampling, listOfNames, "PDF")
  cat("\n\n", file = filename, append = TRUE)
  if (!is.null(dim(listOfNamesMetropolis))) {
    savePlottingInstruction(filename, numberOfVertices, simpleSampling, listOfNamesMetropolis,
                            "unbiased")
    cat("\n\n", file = filename, append = TRUE)
  }
  if (! noGlueing) {
    # plot shifted histograms
    savePlottingInstruction(filename, numberOfVertices, simpleSampling, listOfNames, "shifted")
    cat("\n\n", file = filename, append = TRUE)
    
    # plot the result with logarithmic y-axis
    cat(paste('set output "', numberOfVertices, '_result.eps"\n', sep = ""),
        file = filename, append = TRUE)
    cat(paste('plot [:][:] "hist_', numberOfVertices, '_result.dat" using 1:2 title "result"\n\n', sep = ""),
        file = filename, append = TRUE)
    
    # plot rate function
    cat(paste('set output "', numberOfVertices, '_rate-function.eps"\n', sep = ""),
        file = filename, append = TRUE)
    cat(paste('plot [:][:] "hist_', numberOfVertices, '_rate-function.dat" using 1:2 title "rate-function"\n\n', sep = ""),
        file = filename, append = TRUE)
    
    # plot the result with non logarithmic y-axis
    cat("unset logscale y\n\n", file = filename, append = TRUE)
    cat('unset format y\n\n', file = filename, append = TRUE)
    cat(paste('set output "', numberOfVertices, '_result_no-log.eps"\n\n', sep = ""),
        file = filename, append = TRUE)
    cat(paste('plot [:][:] "hist_', numberOfVertices, '_result.dat" using 1:2 title "result"\n\n', sep = ""),
        file = filename, append = TRUE)
  }
  
  # call gnuplot
  system(paste("cd ", directory, "; gnuplot plot.gpl", sep = ""), intern = FALSE)
}

#' Generate a plot command for the given parameter and save it to the file.
#' 
#' @param filename The filename where the result should be saved.
#' @param numberOfVertices The number of vertices of the graph.
#' @param simpleSampling Call for simple sampling.
#' @param listOfNames A list containing entries for plotting importance sampling results.
#' @param kind The kind of plot (e. g. PDF, unbiased or shifted).
savePlottingInstruction <- function(filename, numberOfVertices, simpleSampling, listOfNames, kind) {
  cat(paste('set output "', numberOfVertices, '_', kind, '.eps"\n', sep = ""),
      file = filename, append = TRUE)
  cat(preparePlotCall(simpleSampling, listOfNames, kind), file = filename, append = TRUE)
}

#' Prepare a plot call in gnuplot.
#' 
#' @param simpleSampling Call for simple sampling.
#' @param listOfNames A list containing entries for plotting importance sampling results.
#' @param kind The kind of plot (e. g. PDF, unbiased or shifted).
#' 
#' @return The plot instruction in gnuplot.
preparePlotCall <- function(simpleSampling, listOfNames, kind) {
  res <- paste('plot [:][:] ', simpleSampling, sep = "")
  for (i in 1:dim(listOfNames)[1]) {
    res <- paste(res, ', "', listOfNames[i, 2], '_', kind,
                 '.dat" using 1:2 title "', listOfNames[i, 1], '"', sep = "")
  }
  return(res)
}

#' Prepare a list with entries for plotting for all data obtained by the Metropolis algorithm.
#' 
#' @param numberOfVertices The number of vertices of the graph.
#' @param data A matrix containing the data to process.
#' 
#' @return The calculated list containing entries for plotting.
prepareListOfNamesMetropolis <- function(numberOfVertices, data) {
  listOfNames <- c()
  for (i in 1:dim(data)[1]) {
    if (!is.na(data[i, 1]) && !is.na(data[i, 2])) {
      listOfNames <- rbind(listOfNames, c(paste("{/Helvetica-Italic T} = ", data[i, 1], sep = ""),
                           paste('hist_is_', numberOfVertices, '_',
                           data[i, 2], '_', data[i, 1], sep = "")))
    }
  }
  return(listOfNames)
}

#' Prepare a list with entries for plotting for all data obtained by the Wang-Landau algorithm.
#' 
#' @param numberOfVertices The number of vertices of the graph.
#' @param data A matrix containing the data to process.
#' 
#' @return The calculated list containing entries for plotting.
prepareListOfNamesWangLandau <- function(numberOfVertices, data) {
  listOfNames <- c()
  for (i in 1:dim(data)[1]) {
    if (!is.na(data[i, 1]) && !is.na(data[i, 2])) {
      listOfNames <- rbind(listOfNames, c(paste("{/Helvetica-Italic min} = ", data[i, 1],
                                                " {/Helvetica-Italic max} = ", data[i, 2], sep = ""),
                                          paste('density_wl_', numberOfVertices, '_',
                                          data[i, 1], '_', data[i, 2], sep = "")))
    }
  }
  return(listOfNames)
}
