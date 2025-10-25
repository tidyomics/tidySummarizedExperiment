#' @name as_tibble
#' @rdname as_tibble
#' @inherit tibble::as_tibble
#' @return `tibble`
#' 
#' @examples
#' tidySummarizedExperiment::pasilla %>%
#'     as_tibble()
#'     
#' tidySummarizedExperiment::pasilla %>%
#'     as_tibble(.subset=-c(condition, type))
#' 
#' @importFrom purrr reduce
#' @importFrom purrr map
#' @importFrom purrr when
#' @importFrom tidyr spread
#' @importFrom stats setNames
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom pkgconfig get_config
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' 
#' Müller, K., Wickham, H. (2023). tibble: Simple Data Frames. R package version 3.2.1, https://CRAN.R-project.org/package=tibble
#' @export
as_tibble.SummarizedExperiment <- function(x, ...,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {

    .as_tibble_optimised(x = x, ...,
        .name_repair=.name_repair, rownames=rownames)

}

.as_tibble_optimised <- function(x, skip_GRanges=FALSE, .subset=NULL,
    .name_repair=c("check_unique", "unique", "universal", "minimal"),
    rownames=pkgconfig::get_config("tibble::rownames", NULL)) {
  
    .subset <- enquo(.subset)
  
    sample_info <-
        colData(x) %>%
    
        # If reserved column names are present add .x
        change_reserved_column_names(x) %>% 

        # Convert to tibble
        tibble::as_tibble(rownames=s_(x)$name) %>% 
        setNames(c(s_(x)$name, colnames(colData(x))))
  
    range_info <-
        skip_GRanges %>%
        when(
            (.) ~ tibble() %>% list,
            ~  get_special_datasets(x) 
        ) %>%
        reduce(left_join, by="coordinate") 
    
    gene_info <-
        rowData(x) %>% 
    
        # If reserved column names are present add .x
        change_reserved_column_names(x)%>% 
  
        # Convert to tibble
        tibble::as_tibble(rownames=f_(x)$name) %>% 
        setNames(c(f_(x)$name, colnames(rowData(x))))
  
    count_info <- get_count_datasets(x)
  
    # Return 
    if (quo_is_null(.subset))
    
        # If I want to return all columns
        count_info %>%
            full_join(sample_info, by=s_(x)$name) %>%
            full_join(gene_info, by=f_(x)$name) %>%
            when(nrow(range_info) > 0 ~ 
                (.) %>% left_join(range_info) %>% suppressMessages(),
                ~ (.)) 
    
    # This function outputs a tibble after subsetting the columns
    else subset_tibble_output(x, count_info, sample_info,
        gene_info, range_info, !!.subset)
}
