#' Blueprint+Encode reference dataset for human
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"Blue.comb"

#' Human Primary Cell Atlas (HPCA) reference dataset for human
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"hpca.comb"

#' Gene order
#'
#' @format A dataframe with the following items:
#' \describe{
#'   \item{seqnames}{chromosome names}
#'   \item{start}{gene start}
#'   \item{end}{gene end}
#' }
"geneOrder"

#' MMR gene expression matrix
#'
#' @format A matrix with the following items:
#' \describe{
#'   \item{colnames}{sample id}
#'   \item{rownames}{MMR genes}
#'   \item{value}{mean expression in malignant cells}
#' }
"mmr_matrix"

