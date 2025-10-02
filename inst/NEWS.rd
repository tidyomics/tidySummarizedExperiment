\name{NEWS}
\title{News for Package \pkg{tidySummarizedExperiment}}

\section{Changes in version 1.19.7}{
\itemize{
    \item Refactored \code{mutate()}, \code{filter()}, \code{select()} and join methods for \code{SummarizedExperiment} with query-scope analysis and operation decomposition.
    \item Added query-scope analyzers (e.g., \code{analyze_query_scope_mutate()}, \code{analyze_query_scope_filter()}) with mixed-scope handling and recording of the latest mutate scope.
    \item Introduced utilities to link query results to \code{SummarizedExperiment} slots and mapping helpers.
    \item Added helpers to extract column and row data from mutated tibbles.
    \item Split and modernized join routines; updated documentation including \code{left_join()} examples.
    \item Replaced deprecated \code{mutate_features()} and \code{mutate_samples()} with internal helpers \code{modify_features()} and \code{modify_samples()}, and added internal helper \code{mutate_assay()} for assay-only mutations; user-facing \code{mutate()} behavior is unchanged. Removed deprecated documentation for the old helpers.
    \item Added an experimental plyxp-backed execution path for mixed-scope operations (across assays, \code{colData}, and \code{rowData}) to improve correctness and performance.
    \item Performance: mutate benchmarks show substantial speed-ups versus \code{master}. Typical median elapsed times drop from approximately 280-350 ms to 10-35 ms for assay- and colData-only scenarios and for chained mutates (about 8x-30x faster). Mixed assay+colData scenarios improve from approximately 290 ms to 25-40 ms. Grouped mean shows parity around 80-100 ms. Benchmarks are in \code{dev/benchmark-scoping.Rmd} with cached results in \code{vignettes/benchmark_results.rda}.
    \item Improved decomposition in \code{decompose_tidy_operations()}, mapping, and related utilities for robustness and performance.
    \item Query decomposition utilities: added internal \code{is_composed()} and \code{decompose_tidy_operation()} to break down complex dplyr operations into sequential steps for granular execution and introspection.
    \item Filter analysis: added internal \code{analyze_query_scope_filter()} and updated \code{filter.SummarizedExperiment} to use scope analysis and to record latest filter scope metadata.
    \item Join improvements: added internal \code{analyze_query_scope_join()} and refactored \code{left_join()}, \code{inner_join()}, \code{right_join()}, and \code{full_join()} for smarter dispatch with fallbacks.
    \item Namespace and metadata: updated \code{NAMESPACE} imports (including tidyselect helpers), removed deprecated exports, and bumped \pkg{plyxp} (Imports) and \pkg{airway} (Suggests).
    \item Vignette updates: fixed \code{count()} example, clarified \code{rename()} usage, updated website links, cached benchmark results, and moved benchmarking vignette to \code{dev/}.
    \item Metadata updates: added \pkg{airway} to \code{Suggests} and \pkg{plyxp} to \code{Imports}.
    \item Expanded unit tests for mutate, filter, select, joins, and query-to-slot linking.
}}

\section{Changes in version 1.19.5}{
\itemize{
    \item Soft deprecated \code{bind_rows()} in favor of \code{append_samples()} from ttservice.
    \item Added \code{append_samples()} method for SummarizedExperiment objects.
    \item \code{bind_rows()} is not a generic method in dplyr and may cause conflicts.
    \item Users are encouraged to use \code{append_samples()} instead.
}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item Updated documentation to properly reflect S3 method structure.
    \item Simplified ggplot2 compatibility - S3 method continues to work with ggplot2 4.0.0.
    \item Users are now directed to library(tidyprint) for tidy visualization and https://github.com/tidyomics/tidyprint for more information.
}}

\section{Changes in version 1.4.0, Bioconductor 3.14 Release}{
\itemize{
    \item Improved join_*() functions.
    \item Changed special column names with a starting "." to avoid conflicts with pre-existing column names.
    \item Improved all method for large-scale datasets.
}}

\section{Changes in version 1.5.3, Bioconductor 3.15 Release}{
\itemize{
    \item Speed-up nest.
    \item Adaptation to Ranged-SummarizedExperiment.
}}

\section{Changes in version 1.7.3, Bioconductor 3.16 Release}{
\itemize{
    \item Fixed as_tibble edge case
    \item Fixed print for DelayedArray
    \item Improve performance for large-scale datasets
    \item Fixed filter is the result is a no-gene dataset, and improve performance of filtering
}}

