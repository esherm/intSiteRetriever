connectToDB = function(){
  library("RMySQL") #also loads DBI
  all_cons <- dbListConnections(MySQL())
  for (con in all_cons) {
    discCon <- dbDisconnect(con)
  }
  dbConn = dbConnect(MySQL(), group="intSitesDEV-dev") #~/.my.cnf must be present 
  dbConn
}


#parse a vector of strings into a piped together list ready for SQL REGEXP
#is tolerant of MySQL's '%' wildcard which is unfortunately used pretty extensively in legacy code
#allows single distinct queries
parseSetNames = function(setName){
  paste0(gsub("%", "(.*)", paste0("^", setName, "$")), collapse="|")
}


getUniqueSites = function(setName){
  dbConn = connectToDB() 
  suppressWarnings(dbGetQuery(dbConn, paste0("SELECT sites.position, sites.chr, sites.strand, samples.sampleName
                                              FROM sites, samples
                                              WHERE sites.sampleID = samples.sampleID
                                              AND samples.sampleName REGEXP ", dbQuoteString(dbConn, parseSetNames(setName)), "AND sites.multihitID IS NULL;"))) 
}

getMRCs = function(){
  dbConn = connectToDB() 
  
}

getMultihits = function(){
  dbConn = connectToDB() 
  
}

getPCRbreaks = function(){
  dbConn = connectToDB() 
  
}

setNameExists = function(setName){
  dbConn = connectToDB() 
  
}