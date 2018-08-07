downGSE <- function(studyID = "GSE1009", destdir = ".") {
  
  library(GEOquery)
  eSet <- getGEO(studyID, destdir = destdir, getGPL = F)
  
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  
  write.csv(exprSet, paste0(studyID, "_exprSet.csv"))
  write.csv(pdata, paste0(studyID, "_metadata.csv"))
  return(eSet)
  
}

downGSE(studyID = "GSE75478")


get_GSE_links <- 
  function(studyID = "GSE1009", down = F, destdir = "./") {
    ## studyID destdir ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/matrix/GSE1009_series_matrix.txt.gz
    ## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/suppl/GSE1009_RAW.tar
    ## http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&mode=csv&series=1009
    
    supp_link = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", substr(studyID, 1, nchar(studyID) - 3), "nnn/", 
                       studyID, "/suppl/", studyID, "_RAW.tar")
    meta_link = paste0("http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&mode=csv&series=", substr(studyID, 
                                                                                                       4, nchar(studyID)))
    matrix_link = paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/", substr(studyID, 1, nchar(studyID) - 3), "nnn/", 
                         studyID, "/matrix/", studyID, "_series_matrix.txt.gz")
    
    print(paste0("The URL for raw data is : ", supp_link))
    print(paste0("The URL for metadata is : ", meta_link))
    print(paste0("The URL for matrix   is : ", matrix_link))
    
  }
get_GSE_links("GSE75478")


GSE75478 = GEOquery::getGEO("GSE75478")
length(GSE75478)
GSE75478_Data = GSE75478[[1]]
exprMatrix=exprs(GSE75478_Data)
dim(exprMatrix)


eList <- GEOquery::getGEO("GSE11675")
length(eList)
eData = eList[[1]]
eData
names(pData(eData))
exprMatrix=exprs(eData)
dim(exprMatrix)
