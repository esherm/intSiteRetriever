library(BSgenome)
library(dplyr)

# generate uniform random position on the reference genome
# for all chromosomes(but  gender specific)
# gender is encoded as: 'm' or 'f'
# @note length of chromosome is used as a weight

#' return BS genome object for human readable UCSC format
#' 
#' format examples are: hg18, hg19, hg38 for human
#'                      mm8, mm9 for mouse
#' @stop if cannot find unique genome from installed BSgenome
get_reference_genome <- function(reference_genome) {
    pattern <- paste0("\\.", reference_genome, "$")
    match_index <- which(grepl(pattern, installed.genomes()))
    if (length(match_index) != 1) {
        write("Installed genomes are:", stderr())
        write(installed.genomes(), stderr())
        stop(paste("Cannot find unique genome for", reference_genome))
    }
    BS_genome_full_name <- installed.genomes()[match_index]
    library(BS_genome_full_name, character.only=T)
    get(BS_genome_full_name)
}

#' for a given reference genome and gender generate random positions 
#'
#' @param siteIDs vector of unique siteIDs for use as random seed
#' @param reference_genome BS object reference genome(@seealso get_reference_genome)
#' @param gender 'm' or 'f
#' @param number_of_positions total number of random positions to generate
#' @param male_chr list of male-specific chromosomes prefixes(only 1 prefix is allowed at present)
#' @return dataframe with columns: siteID, chr, strand, position
get_random_positions <- function(siteIDs, reference_genome, gender,
                                 number_of_positions=3, male_chr=c("chrY")){
  stopifnot(length(male_chr) == 1)
  stopifnot(length(gender) == 1)
  stopifnot(check_gender(gender))

  chr_len <- seqlengths(reference_genome)
  stopifnot(any(grepl(male_chr, names(chr_len)))) # male chomosome is in genome
  chr_len <- get_gender_specific_chr(chr_len, gender, male_chr)
  chr_len <- chr_len[names(chr_len) != "chrM"] #remove mitochondria

  cs <- c(0,cumsum(as.numeric(chr_len)))
  genomeLength <- max(cs)

  seed <- .Random.seed #don't want to screw up global randomness

  mrcs <- sapply(siteIDs, function(x){
    set.seed(x)
    rands <- round(runif(number_of_positions, 1, genomeLength*2)-genomeLength)
    cuts <- cut(abs(rands), breaks=cs, labels=names(chr_len))

    #outputs in format of "siteID", "chr", "strand", "position"
    list(rep(x,number_of_positions),
         as.character(cuts),
         as.character(cut(sign(rands), breaks=c(-1,0,1), labels=c("-", "+"), include.lowest=T)),
         abs(rands) - cs[match(cuts, names(chr_len))])
  })

  .Random.seed <- seed #resetting the seed

  #this is ugly and requires as.character above but ~5X faster than outputting a list of
  #data.frames and then rbinding them all together
  mrcs <- data.frame(matrix(unlist(t(mrcs)), nrow=length(siteIDs)*number_of_positions))
  names(mrcs) <- c("siteID", "chr", "strand", "position")

  mrcs
}

#' generate random controls for sites
#'
#' @param sites_meta dataframe with columns: siteID, gender
#' @param reference_genome BS genome object. All sites have the same genome.
#' @param number_mrcs_per_site number of MRCs to generate for each site 
#' @return dataframe with columns: siteID, chr, strand, position
#'
#' @note siteID are the same as given by sites_meta df
get_N_MRCs <- function(sites_meta, reference_genome, number_mrcs_per_site=3) {
    stopifnot(setequal(names(sites_meta), c("siteID", "gender")))
    stopifnot(number_mrcs_per_site > 0)
    stopifnot(check_gender(sites_meta$gender))

    num_sites <- nrow(sites_meta)
    tot_num_mrcs <- num_sites * number_mrcs_per_site

    mrcs <- lapply(split(sites_meta, sites_meta$gender), function(sites){
      get_random_positions(sites$siteID, reference_genome, sites[1,"gender"],
                           number_mrcs_per_site)
    })

    unrowname(do.call(rbind, mrcs))
}

#' from vector of chromosome lengths with names creates vector for male or female
#' @param all_chromosomes vector with length, names(all_chromosomes) are actual names of chromosome
get_gender_specific_chr <- function(all_chromosomes, gender, male_chr) {
    stopifnot( ! is.null(names(all_chromosomes)))
    stopifnot(length(male_chr) == 1)
    if (gender == 'm') {
        return(all_chromosomes)
    }
    stopifnot(gender =='f')
    chromosome_names <- names(all_chromosomes)
    female_chromosomes <- all_chromosomes[ ! grepl(male_chr, chromosome_names)]
    female_chromosomes
}

#' gender can only be male('m') or female('f')
check_gender <- function(gender) {
    valid <- c('m', 'f')
    values <- unique(gender)
    all(values %in% valid)
}

