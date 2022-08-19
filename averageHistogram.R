# averageHistogram.R
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

#' Calculate average value plus standard error for the given histogram.
#' 
#' @param histogram The pdf of the given histogram.
#' 
#' @return Average value and standard error of the given histogram.
averageHistogram <- function(histogram) {
  avg <- sum(histogram[, 1] * histogram[, 2]) / sum(histogram[, 2])
  var <- sum(histogram[, 1]^2 * histogram[, 2]) / sum(histogram[, 2]) - avg^2
  err <- sqrt(var / (length(histogram[, 2]) - 1))
  return(c(avg, err))
}
