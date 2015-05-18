.connectToDB <- function(dbConn){
  if(is.null(dbConn)){
    library("RMySQL") #also loads DBI
    all_cons <- dbListConnections(MySQL())
    for (con in all_cons) {
      discCon <- dbDisconnect(con)
    }
    dbConn <- dbConnect(MySQL(), group="intSitesDEV-dev") #~/.my.cnf must be present 
    dbConn
  }else{
    dbConn
  }
}

#this is a hack-around version of dbQuoteString that does the exact same thing
#but does not require an instance of DBIConnection to be passed along.  This
#gets around the issue where we have to SQL quote before a DB connection is
#established
.quoteForMySQL <- function(x){
  x <- gsub("'", "''", x, fixed=TRUE)
  str <- paste("'", encodeString(x), "'", sep = "")
  str[is.na(x)] <- "NULL"
  SQL(str)
}

.disconnectFromDB <- function(dbConn, conn){ #conn is passed in from user, dbConn is the actual connection
  if(is.null(conn)){ #we made a temp connection which needs to be closed
    dbDisconnect(dbConn)
  }
}

#parse a vector of strings into a piped together list ready for SQL REGEXP
#is tolerant of MySQL's '%' wildcard which is unfortunately used pretty extensively in legacy code
#allows single distinct queries
.parseSetNames <- function(setName){
  .quoteForMySQL(paste0(gsub("%", "(.*)", paste0("^", setName, "$")), collapse="|"))
}

.intSiteRetrieverQuery <- function(command, conn){
  dbConn <- .connectToDB(conn)
  res <- suppressWarnings(dbGetQuery(dbConn, command))
  .disconnectFromDB(dbConn, conn)
  res
}

getUniqueSites <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        sites.position,
                                        samples.sampleName
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                 " AND sites.multihitID IS NULL;"), conn)
}

getMRCs <- function(setName, conn=NULL){
  sites.metadata <- .intSiteRetrieverQuery(paste0("SELECT sites.siteID,
                                                         samples.refGenome,
                                                         samples.gender,
                                                         samples.sampleName
                                                  FROM sites, samples
                                                  WHERE sites.sampleID = samples.sampleID
                                                  AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                                  " AND sites.multihitID IS NULL;"), conn)
  
  sites.metadata <- split(sites.metadata, with(sites.metadata, paste0(refGenome, ".", gender)))
  
  #hardcoding reference genome data/info for MRC creation might not be the best idea
  #but it will work for now.  Would it be appropriate to put into the database?
  #it would be a tiny, infrequently updated table
  
  refGenomeInfo <- list("hg18.M"=Seqinfo(seqnames=c("chr10", "chr11", "chr12", "chr13", "chr14",
                                                   "chr15", "chr16", "chr17", "chr18", "chr19",
                                                   "chr1", "chr20", "chr21", "chr22", "chr2",
                                                   "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                                   "chr9", "chrM", "chrX", "chrY"),
                                        seqlengths=c(135374737, 134452384, 132349534, 114142980,
                                                     106368585, 100338915, 88827254, 78774742,
                                                     76117153, 63811651, 247249719, 62435964,
                                                     46944323, 49691432, 242951149, 199501827,
                                                     191273063, 180857866, 170899992, 158821424,
                                                     146274826, 140273252, 16571, 154913754, 57772954),
                                        isCircular=rep(F,25),
                                        genome=rep("hg18.M",25)),
                       "hg18.F"=Seqinfo(seqnames=c("chr10", "chr11", "chr12", "chr13", "chr14",
                                                   "chr15", "chr16", "chr17", "chr18", "chr19",
                                                   "chr1", "chr20", "chr21", "chr22", "chr2",
                                                   "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                                                   "chr9", "chrM", "chrX"),
                                        seqlengths=c(135374737, 134452384, 132349534, 114142980,
                                                     106368585, 100338915, 88827254, 78774742,
                                                     76117153, 63811651, 247249719, 62435964,
                                                     46944323, 49691432, 242951149, 199501827,
                                                     191273063, 180857866, 170899992, 158821424,
                                                     146274826, 140273252, 16571, 154913754),
                                        isCircular=rep(F,24),
                                        genome=rep("hg18.F",24)))
  
  seed <- .Random.seed #don't want to screw up global randomness
  
  #could possibly be made faster by sticking refGenomeInfo in a DB and not having to worry
  #about all this stuff... currently benchmarking at .372sec/195 sites worth of MRCs
  allMRCs <- do.call(rbind, lapply(sites.metadata, function(sites){
    refGenome <- sites[1,"refGenome"]
    gender <- sites[1,"gender"]
    indexSeqInfo <- refGenomeInfo[[paste0(refGenome, ".", gender)]]
    
    cs <- c(0,cumsum(as.numeric(seqlengths(indexSeqInfo))))
    genomeLength <- max(cs)
    
    mrcs <- sapply(sites$siteID, function(x){
      set.seed(x)
      rands <- round(runif(3, 1, genomeLength*2)-genomeLength)
      cuts <- cut(abs(rands), breaks=cs, labels=seqnames(indexSeqInfo))
      
      #outputs in format of "siteID", "chr", "strand", "position"
      list(rep(x,3),
           as.character(cuts),
           as.character(cut(sign(rands), breaks=c(-1,0,1), labels=c("-", "+"), include.lowest=T)),
           abs(rands) - cs[match(cuts, seqnames(indexSeqInfo))])
      
    })
    #this is ugly and requires as.character above but ~5X faster than outputting a list of
    #data.frames and then rbinding them all together
    mrcs <- data.frame(matrix(unlist(t(mrcs)), nrow=nrow(sites)*3))
    names(mrcs) <- c("siteID", "chr", "strand", "position")
    mrcs
  }))
  
  #keep output consistant across functions
  allMRCs$siteID <- as.numeric(levels(allMRCs$siteID))[allMRCs$siteID]
  allMRCs$position <- as.numeric(levels(allMRCs$position))[allMRCs$position]
  allMRCs$strand <- as.character(allMRCs$strand)
  allMRCs$chr <- as.character(allMRCs$chr)
  
  .Random.seed <- seed #resetting the seed
  merge(allMRCs, do.call(rbind, sites.metadata)[c("siteID", "sampleName")])
}

getMultihits <- function(setName, conn=NULL){
  stop("FUNCTION NOT IMPLEMENTED") 
}

getUniquePCRbreaks <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT pcrbreakpoints.breakpoint,
                                        pcrbreakpoints.count,
                                        sites.position AS integration,
                                        sites.siteID,
                                        sites.chr,
                                        sites.strand,
                                        samples.sampleName
                                 FROM sites, samples, pcrbreakpoints
                                 WHERE (sites.sampleID = samples.sampleID AND
                                        pcrbreakpoints.siteID = sites.siteID)
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                 " AND sites.multihitID IS NULL;"), conn)
}

#setNameExists shouldn't support MySQL % syntax, thus it needs its own dbQuoteString
setNameExists <- function(setName, conn=NULL){
  stopifnot(!any(grepl("%", setName)))
  
  res <- .intSiteRetrieverQuery(paste0("SELECT DISTINCT sampleName
                                                    FROM samples
                                                    WHERE sampleName REGEXP ", .quoteForMySQL(paste0("^", setName, "$", collapse="|")), ";"), conn)
  setName %in% res$sampleName
}

getRefGenome <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        samples.refGenome
                                 FROM samples
                                 WHERE samples.sampleName REGEXP ", .parseSetNames(setName), ";"), conn)
}

getReadCounts <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        SUM(pcrbreakpoints.count) AS readCount
                                 FROM sites, samples, pcrbreakpoints
                                 WHERE (sites.sampleID = samples.sampleID AND
                                        pcrbreakpoints.siteID = sites.siteID)
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                 " AND sites.multihitID IS NULL
                                 GROUP BY sites.sampleID;"), conn)
}

getUniqueSiteCounts <- function(setName, conn=NULL){
  .intSiteRetrieverQuery(paste0("SELECT samples.sampleName,
                                        COUNT(*) AS uniqueSites
                                 FROM sites, samples
                                 WHERE sites.sampleID = samples.sampleID
                                 AND samples.sampleName REGEXP ", .parseSetNames(setName),
                                 " AND sites.multihitID IS NULL
                                 GROUP BY sites.sampleID;"), conn)
}