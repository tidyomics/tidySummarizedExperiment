#' tidySummarizedExperiment: Bringing SummarizedExperiment to the Tidyverse
#'
#' @description
#' The tidySummarizedExperiment package provides a bridge between Bioconductor
#' SummarizedExperiment objects and the tidyverse. It creates an invisible layer
#' that enables viewing the Bioconductor SummarizedExperiment object as a tidyverse
#' tibble, and provides SummarizedExperiment-compatible dplyr, tidyr, ggplot2 and
#' plotly functions. This allows users to get the best of both Bioconductor and
#' tidyverse worlds.
#'
#' @details
#' The main functions and methods are:
#'
#' \itemize{
#' \item \code{\link{tidy}} - Convert SummarizedExperiment to tidy format
#' \item \code{\link{as_tibble}} - Convert to tibble representation
#' \item \code{\link{filter}}, \code{\link{select}}, \code{\link{mutate}}, \code{\link{arrange}} - dplyr verbs
#' \item \code{\link{pivot_longer}}, \code{\link{pivot_wider}}, \code{\link{nest}}, \code{\link{unnest}} - tidyr functions
#' \item \code{\link{ggplot}} - ggplot2 visualization
#' \item \code{\link{left_join}}, \code{\link{inner_join}}, \code{\link{right_join}}, \code{\link{full_join}} - dplyr joins
#' }
#'
#' For detailed information on usage, see the package vignette, by typing
#' \code{vignette("introduction", package = "tidySummarizedExperiment")}.
#'
#' All software-related questions should be posted to the Bioconductor Support Site:
#'
#' \url{https://support.bioconductor.org}
#'
#' The code can be viewed at the GitHub repository:
#'
#' \url{https://github.com/stemangiola/tidySummarizedExperiment}
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics
#' ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024).
#' \doi{10.1038/s41592-024-02299-2}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/stemangiola/tidySummarizedExperiment}
#'   \item Report bugs at \url{https://github.com/stemangiola/tidySummarizedExperiment/issues}
#' }
#'
#' @author Stefano Mangiola
#'
#' @docType package
#' @name tidySummarizedExperiment-package
#' @aliases tidySummarizedExperiment
#' @keywords internal
"_PACKAGE"

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    attached <- tidyverse_attach()
    
    # Print loading message about printing
    packageStartupMessage(
        "tidySummarizedExperiment says: Printing is now handled externally. ",
        "If you want to visualize the data in a tidy way, do library(tidyprint). ",
        "See https://github.com/tidyomics/tidyprint for more information."
    )
}