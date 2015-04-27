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
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT sites.position,
                                                           sites.chr,
                                                           sites.strand,
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
  
  .disconnectFromDB(dbConn, conn)
}

getMultihits = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
 
  .disconnectFromDB(dbConn, conn)
}

getUniquePCRbreaks = function(setName, conn=NULL){
  dbConn = .connectToDB(conn)
  res = suppressWarnings(dbGetQuery(dbConn, paste0("SELECT pcrbreakpoints.breakpoint,
                                                           pcrbreakpoints.count,
                                                           sites.position AS integration,
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

