source("random_site.R")

context("Testing random positions generation.")

reference <- get_reference_genome("hg18")

test_that("generate_random_positions returns df", {
    randoms <- get_random_positions(1, reference, 'm')
    has_names(randoms, c("chr", "position", "strand"))
})

test_that("generate_random_positions only accepts male or female", {
    expect_error(get_random_positions(1, reference, 'xxx'))
    expect_error(get_random_positions(10, reference, 'A'))
})

test_that("size of the df is number_of_positions", {
    number_of_positions <- seq(1:10)
    sapply(42, function(x)
        expect_equal(nrow(get_random_positions(x, reference, 'f')), x)
    )
})

test_that("male-specific chromosome in the genome", {
    expect_error(get_random_positions(
        3, reference, 'f', c("chr_EVER_UNKNOWN")))
})

sites_meta <- data.frame( 
    siteID=c(1, 2, 3),
    gender=c('f', 'm', 'm'),
    sampleName=c("subject1", "subject1", "subject2")
)

test_that("get_N_MRCs have all required columns", {
    mrcs <- get_N_MRCs(sites_meta, reference, 3)
    has_names(mrcs, c("siteID", "chr", "position", "strand"))
})

test_that("get_N_MRCs only uses given siteIDs and do not introduce new one", {
    mrcs <- get_N_MRCs(sites_meta, reference, 3)
    expect_true(setequal(mrcs$siteID, sites_meta$siteID))
})

test_that("get_N_MRCs return df", {
    N <- 7
    mrcs <- get_N_MRCs(sites_meta, reference, N)
    expect_equal(dim(mrcs), c(nrow(sites_meta)*N), 4)
})

test_that("get_N_MRCs return the same result for the same siteIDs", {
    mrcs1 <- get_N_MRCs(sites_meta, reference, 9)
    mrcs2 <- get_N_MRCs(sites_meta, reference, 9)
    expect_equal(mrcs1, mrcs2)
})
