#' @importFrom methods getMethod
setClass("tidySummarizedExperiment",
    contains=c("SummarizedExperiment", "RangedSummarizedExperiment"))

#' @name tidy
#' @rdname tidy
#' @title tidy for `Seurat`
#'
#' @param object A `Seurat` object.
#' @return A `tidyseurat` object.
#'
#' @examples
#' data(pasilla)
#' pasilla %>% tidy()
#'
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#'
#' @importFrom generics tidy

#' @importFrom lifecycle deprecate_warn
tidy_ <- function(object) {
    
    # DEPRECATE
    deprecate_warn(
        when = "1.1.1",
        what = "tidy()",
        details = "tidySummarizedExperiment says: tidy() is not needed anymore."
    )
    
    object
}

#' @importFrom methods as
#' @rdname tidy
#' @param object A SummarizedExperiment object
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
tidy.SummarizedExperiment <- tidy_

#' @importFrom methods as
#' @rdname tidy
#' @param object A SummarizedExperiment object
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
tidy.RangedSummarizedExperiment <- tidy_