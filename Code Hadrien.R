## Load dataset
## ------------------
cat("Load dataset Module analysis \n")
## ------------------
sort.by <- "none"
daysIn <- c(0,1,3,7)
thCount <- 5
removeNoCompTemp <- T
group <- c(1,2)
formula <- ~id_num+day
dayToFactor <- F
removePat2and7 <- F
removeNoCounts <- T
removeOthers <- c("X101_d1",
                  "X102_d0","X102_d1",
                  "X103_d0","X103_d1","X103_d3",
                  "X104_d0",
                  "X107_d0",
                  paste("X109_d",daysIn,sep=""),
                  "X111_d3","X111_d7",
                  paste("X114_d",daysIn,sep=""),
                  "X115_d1",
                  "X116_d3","X116_d7",
                  "X117_d0","X117_d1","X117_d7")

# ---------------------New object DESeq function
newObjectDESeq <-
  function(daysIn=c(0,1),formula=~day,dayToFactor=F,removePat2and7=T,removeNoCounts=T,group=c(1,2),
           removeOthers=NULL,thCount=0){
    return(list(daysIn=daysIn,formula=formula,removeOthers=removeOthers,thCount=thCount,rVersion=getRversion(),
                dayToFactor=dayToFactor,removePat2and7=removePat2and7,removeNoCounts =
                  removeNoCounts,id=sample(x=1:10^6,size = 1)))
  }
# ---------------------createNewMetas function

createNewMetas <- function(newMetasDESeq){
  test <- isNewMetasInAllMetas(newMetasDESeq)
  nameData <- paste("../../data/",test,".rds",sep="")
  out <- tryCatch(readRDS(nameData),
                  warning=function(w) {
                    print("Create new dataset")
                    out <- 
                      runDESeq(daysIn=newMetasDESeq$daysIn,formula=newMetasDESeq$formula,
                               removeOthers=newMetasDESeq$removeOthers,thCount=newMetasDESeq$thCount,
                               dayToFactor=newMetasDESeq$dayToFactor,removePat2and7=newMetasDESeq$removePat2and7,
                               removeNoCounts = 
                                 newMetasDESeq$removeNoCounts)
                    saveRDS(object=out,file=nameData)
                    allMetas <- openAllMetas()
                    allMetas[[length(allMetas)+1]] <- newMetasDESeq
                    saveRDS(object = allMetas,file = 
                              "../../data/allMetas.rds")
                    return(out)
                  }
  )
  return(out)
}


# --------------isNewMetasInAllMetas function
isNewMetasInAllMetas <- function(newMetasDESeq){
  allMetas <- openAllMetas()
  tests <- 0
  m <- 0
  out <- newMetasDESeq$id
  tests <- NULL
  if(length(allMetas)>0){
    for(m in 1:length(allMetas)){
      tests <- c(tests,compareTwoMetas(newMetasDESeq,allMetas[[m]]))
    }
    if(sum(tests)<length(allMetas)){
      out <- allMetas[[which(tests==0)]]$id
    }
  }
  return(out)
}




#------------ Continue analysis
day0 <- 0
day1 <- 1
N <- 10
type <- "DC"
newMetasDESeq.In <-
  newObjectDESeq(daysIn,formula,dayToFactor,removePat2and7,
                 removeNoCounts,group,removeOthers,thCount)
DESeqOut <- createNewMetas(newMetasDESeq.In)



## Get meta information
dds <- DESeqOut$dds
day_ok <- DESeqOut$day_ok
id_num_ok <- DESeqOut$id_num_ok
pos_rows_ok <- DESeqOut$pos_rows_ok
se <- DESeqOut$se
counts <- DESeq::counts(dds)
toKeep <- which(SummarizedExperiment::colData(se)$day %in% c(0,1,3,7))
coldataSE <- SummarizedExperiment::colData(se)[toKeep,]
day <- factor(coldataSE$day)#,labels = c("_d0","_d1","_d3","_d7"))
id_num <- factor(coldataSE$id_num)
group <- rep(1,length(id_num))
group[which(as.numeric(levels(id_num)[id_num])>10)] <-2
design <- model.matri~0+day);colnames(design);rownames(design) <-
  colnames(counts)