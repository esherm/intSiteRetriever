library(BSgenome)

# generate uniform random position on the reference genome
# for all chromosomes(but  gender specific)

#' return BS genome object for human readable UCSC format
#' 
#' format examples are: hg18, hg19, hg38 for human
#'                      mm8, mm9 for mouse
#' @return BS object
get_reference_genome <- function(reference_genome) {
    pattern <- paste0("\\.", reference_genome, "$")
    match_index <- which(grepl(pattern, installed.genomes()))
    if (length(match_index) != 1) {
        write("Installed genomes are:", stderr())
        write(installed.genomes(), stderr())
        stop(paste("Cannot find unique genome for", reference_genome))
    }
    BS_genome_full_name <- installed.genomes()[match_index]
    get(BS_genome_full_name)
}

#' for a given reference genome and gender generate random positions 
#'
#' @param number_of_positions total number of random positions to generate
#' @param reference_genome BS object reference genome(@see get_reference_genome)
#' @param gender 'm' or 'f'
#' @param male_chr list of male-specific chromosomes
#' @return dataframe with columns: chr, strand, position
get_random_positions <- function(
    number_of_positions, reference_genome, gender, male_chr=c("chrY")
) {
    
    # check gender is only: 'm' or 'f'

    # get chromosome lengths

    # get random sites
}

#' generate random controls for samples and sites
#'
#' @param sites_meta dataframe with columns: siteID, gender, sampleName 
#' @param reference_genome BS genome object. All sites have the same genome.
#' @param number_mrcs_per_site number of MRCs to generate for each site 
#' @return dataframe with columns: siteID, chr, strand, position
#'
#' @note siteID are the same as given by sites df
get_N_MRCs <- function(sites_meta, reference_genome, number_mrcs_per_site=3) {

    # check names are as expected

    # number_mrcs_per_site > 0

    # check gender is only: 'm' or 'f'

    # for each site generate number_mrcs_per_site MRCS

}


#' gender can only be male('m') or female('f')
check_gender <- function(gender) {

}
