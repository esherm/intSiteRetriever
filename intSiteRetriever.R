.connectToDB = function(dbConn){
  if(is.null(dbConn)){
    library("RMySQL") #also loads DBI
    all_cons <- dbListConnections(MySQL())
    for (con in all_cons) {
      discCon <- dbDisconnect(con)
    }
    dbConn = dbConnect(MySQL(), group="intSitesDEV-dev") #~/.my.cnf must be present 
    dbConn
  }else{
    dbConn
  }
}

.disconnectFromDB = function(dbConn, conn){ #conn is passed in from user, dbConn is the actual connection
  if(is.null(conn)){ #we made a temp connection which needs to be closed
    dbDisconnect(dbConn)
  }
}

#parse a vector of strings into a piped together list ready for SQL REGEXP
#is tolerant of MySQL's '%' wildcard which is unfortunately used pretty extensively in legacy code
#allows single distinct queries
#'private' method so we can trust that dbConn is valid when passed in
.parseSetNames = function(setName, dbConn){
  dbQuoteString(dbConn, paste0(gsub("%", "(.*)", paste0("^", setName, "$")), collapse="|"))
}


getUniqueSites = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT sites.siteID,
                                                           sites.chr,
                                                           sites.strand,
                                                           sites.position,
                                                           samples.sampleName
                                                    FROM sites, samples
                                                    WHERE sites.sampleID = samples.sampleID
                                                    AND samples.sampleName REGEXP ", .parseSetNames(setName, dbConn),
                                                  " AND sites.multihitID IS NULL;")))
  .disconnectFromDB(dbConn, conn)
  res
}

getMRCs = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT sites.siteID,
                                                           samples.refGenome,
                                                           samples.gender,
                                                           samples.sampleName
                                                    FROM sites, samples
                                                    WHERE sites.sampleID = samples.sampleID
                                                    AND samples.sampleName REGEXP ", .parseSetNames(setName, dbConn),
                                                  " AND sites.multihitID IS NULL;")))

  res = split(res, with(res, paste0(refGenome, ".", gender)))
  
  #hardcoding reference genome data/info for MRC creation might not be the best idea
  #but it will work for now.  Would it be appropriate to put into the database?
  #it would be a tiny, infrequently updated table
  
  refGenomeInfo = list("hg18.M"=Seqinfo(seqnames=c("chr10", "chr11", "chr12", "chr13", "chr14",
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
  
  seed = .Random.seed #don't want to screw up global randomness
  
  #could possibly be made faster by sticking refGenomeInfo in a DB and not having to worry
  #about all this stuff... currently benchmarking at .372sec/195 sites worth of MRCs
  MRCs = do.call(rbind, lapply(res, function(sites){
    refGenome = sites[1,"refGenome"]
    gender = sites[1,"gender"]
    indexSeqInfo = refGenomeInfo[[paste0(refGenome, ".", gender)]]
    
    cs = c(0,cumsum(as.numeric(seqlengths(indexSeqInfo))))
    genomeLength = max(cs)
        
    do.call(rbind, lapply(sites$siteID, function(x){
      set.seed(x)
      rands = round(runif(3, 1, genomeLength*2)-genomeLength)
      cuts = cut(abs(rands), breaks=cs, labels=seqnames(indexSeqInfo))
      
      data.frame("siteID"=x,
                 "chr"=cuts,
                 "strand"=cut(sign(rands), breaks=c(-1,0,1), labels=c("-", "+"), include.lowest=T),
                 "position"=abs(rands) - cs[match(cuts, seqnames(indexSeqInfo))])
      
    }))
  }))
    
  .Random.seed = seed #resetting the seed
  .disconnectFromDB(dbConn, conn)
  merge(MRCs, do.call(rbind, res)[c("siteID", "sampleName")]) #not quite right
}

getMultihits = function(setName, conn=NULL){
  stop("FUNCTION NOT IMPLEMENTED") 
  dbConn = .connectToDB(conn)
 
  .disconnectFromDB(dbConn, conn)
}

getUniquePCRbreaks = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT pcrbreakpoints.breakpoint,
                                                           pcrbreakpoints.count,
                                                           sites.position AS integration,
                                                           sites.siteID,
                                                           sites.chr,
                                                           sites.strand,
                                                           samples.sampleName
                                                    FROM sites, samples, pcrbreakpoints
                                                    WHERE (sites.sampleID = samples.sampleID AND
                                                          pcrbreakpoints.siteID = sites.siteID)
                                                    AND samples.sampleName REGEXP ", .parseSetNames(setName, dbConn),
                                                  " AND sites.multihitID IS NULL;")))
  .disconnectFromDB(dbConn, conn)
  res
}

#setNameExists shouldn't support MySQL % syntax, thus it needs its own dbQuoteString
setNameExists = function(setName, conn=NULL){
  stopifnot(!any(grepl("%", setName)))
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT DISTINCT sampleName
                                                    FROM samples
                                                    WHERE sampleName REGEXP ", dbQuoteString(dbConn, paste0("^", setName, "$", collapse="|")), ";")))
  .disconnectFromDB(dbConn, conn)
  setName %in% res$sampleName
}

getRefGenome = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT samples.sampleName,
                                                           samples.refGenome                                                    FROM samples
                                                    WHERE samples.sampleName REGEXP ", .parseSetNames(setName, dbConn), ";")))
  .disconnectFromDB(dbConn, conn)
  res
}

