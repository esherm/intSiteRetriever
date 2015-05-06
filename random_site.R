library(BSgenome)

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
    get(BS_genome_full_name)
}

#' for a given reference genome and gender generate random positions 
#'
#' @param number_of_positions total number of random positions to generate
#' @param reference_genome BS object reference genome(@seealso get_reference_genome)
#' @param gender 'm' or 'f'
#' @param male_chr list of male-specific chromosomes prefixes(only 1 prefix is allowed at present)
#' @return dataframe with columns: chr, strand, position
get_random_positions <- function(
    number_of_positions, reference_genome, gender, male_chr=c("chrY")
) {
    stopifnot(length(male_chr) == 1)
    stopifnot(length(gender) == 1)
    stopifnot(check_gender(gender))

    chr_len <- seqlengths(reference_genome)
    stopifnot(any(grepl(male_chr, names(chr_len)))) # male chomosome is in genome
    chr_len <- get_gender_specific_chr(chr_len, gender, male_chr)

    random_sites <- data.frame(
        chr=character(number_of_positions), 
        position=numeric(number_of_positions),
        strand=character(number_of_positions),
        stringsAsFactors=FALSE
    )

    random_sites$chr <- get_random_chromosome(number_of_positions, chr_len)
    random_sites$position <- get_random_positions_on_chromosome(
        random_sites$chr, chr_len)
    random_sites$strand <- get_random_strand(number_of_positions)
    random_sites
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

    mrcs <- data.frame(
        siteID=sites_meta$siteID,
        chr=character(tot_num_mrcs), 
        position=numeric(tot_num_mrcs),
        strand=character(tot_num_mrcs),
        stringsAsFactors=FALSE
    )

    # for each site generate number_mrcs_per_site MRCS
    sapply(1:nrow(sites_meta), function(i) {
        start_index <- (i - 1) * number_mrcs_per_site + 1
        end_index <- start_index + number_mrcs_per_site - 1
        mrcs_for_site <- get_random_positions(
            number_mrcs_per_site, reference_genome, sites_meta$gender[i])
        mrcs_for_site <- cbind(siteID=sites_meta$siteID[i], mrcs_for_site)
        mrcs[start_index:end_index, ] <<- mrcs_for_site
    })
    mrcs
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

#' @param chr_len is a vector with names of chomosomes and values equal to length of the chr
#' @param size number of random chromosomes to generate
#' @return vector of chromosome names with weigths proportional to length of the chromosome
get_random_chromosome <- function(size, chr_len) {
    stopifnot( ! is.null(names(chr_len)))
    sample(names(chr_len), size, replace=TRUE, prob=chr_len)
}

#' @param chr_length is a vector with names of chomosomes and values equal to length of the chr
#' @param chromosomes vector of chromosome names for which random positions are generated
#' @return vector of positions
get_random_positions_on_chromosome <- function(chromosomes, chr_len) {
    stopifnot( ! is.null(names(chr_len)))
    stopifnot(chromosomes %in% names(chr_len))
    random_pos <- sapply(chromosomes , function(chr) {
        floor(runif(n=1, min=1, max=chr_len[[chr]] + 1))
    })
    random_pos
}

#' @return vector of "+" an "-" with the same probability
get_random_strand <- function(size) {
    valid_strand <- c("+", "-")
    sample(valid_strand, size, replace=TRUE)
}

#' gender can only be male('m') or female('f')
check_gender <- function(gender) {
    valid <- c('m', 'f')
    values <- unique(gender)
    all(values %in% valid)
}

