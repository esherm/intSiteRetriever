source("random_site.R")

context("Testing random positions generation.")

reference <- get_reference_genome("hg18")

test_that("generate_random_positions returns df", {
    randoms <- get_random_positions(1, reference, 'm')
    expect_named(randoms, 
        c("siteID", "chr", "position", "strand"), ignore.order=TRUE)
})

test_that("generate_random_positions only accepts male or female", {
    expect_error(get_random_positions(1, reference, 'xxx'))
    expect_error(get_random_positions(10, reference, 'A'))
})

test_that("size of the df is number_of_positions", {
    sapply(42, function(x)
        expect_equal(nrow(get_random_positions(1, reference, 'f', x)), x)
    )
})

test_that("male-specific chromosome in the genome", {
    expect_error(get_random_positions(
        3, reference, 'f', male_chr=c("chr_EVER_UNKNOWN")))
})

sites_meta <- data.frame( 
    siteID=c(1, 2, 3),
    gender=c('f', 'm', 'm')
)

test_that("get_N_MRCs have all required columns", {
    mrcs <- get_N_MRCs(sites_meta, reference, 3)
    expect_named(mrcs, c("siteID", "chr", "position", "strand"), ignore.order=T)
})

test_that("get_N_MRCs only uses given siteIDs and do not introduce new one", {
    mrcs <- get_N_MRCs(sites_meta, reference, 3)
    expect_true(setequal(mrcs$siteID, sites_meta$siteID))
})

test_that("get_N_MRCs return df", {
    N <- 7
    mrcs <- get_N_MRCs(sites_meta, reference, N)
    expect_equal(nrow(mrcs), c(nrow(sites_meta)*N))
    expect_equal(ncol(mrcs), 4)
})

test_that("get_N_MRCs return the same result for the same siteIDs", {
    mrcs1 <- get_N_MRCs(sites_meta, reference, 9)
    mrcs2 <- get_N_MRCs(sites_meta, reference, 9)
    expect_equal(mrcs1, mrcs2)
})
