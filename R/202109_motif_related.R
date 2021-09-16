
#' de novo motif and find uniquemotif from fasta file
#'
#'
#' @param indir  the directory containing the reads for motif discovery, can be either a fastq, fasta, or plain text file,must .gz
#' @param ic_threshold use ic to chose motif
#' @param icscore_threshold use icscore to chose motif
#' @param iclength_treshold use icscore / consensus_length to chose motif
#' @param threshold use Pearson to compare two motif
#' @param thread the number of threads
#'
#'
#' @return
#'the results of de novo motifs and unique motifs,and a html table
#' @export
#'
#' @examples
#' uniquemotif("/wrk/wenchenjin/work/Chenjin1_3__ATI_30N101N/pre_processed/cleaned","/wrk/yuanzhen")
#'
#'
yz_uniquemotif <- function(files,ic_threshold=1,icscore_threshold = 3,iclength_treshold=0.5,threshold=0.95,thread = 30){
  library(TFBSTools)
  library(universalmotif)

  cat("we will deal with your data by four steps","\n")
  cat("1. de novo motif discovry","\n")
  cat("2. find the unique motif","\n")
  cat("3. detect the TF name of the unique motif from JASPAR2020","\n")
  cat("4. visualize the results","\n")

  ###Now we will the first steps to de novo motif discovry
  cat("Now we will the first steps to de novo motif discovry","/n")

  filenamedisc <- Sys.glob(files)


  outdir <- tempdir()
  motifdicr <- paste0(outdir,"/","motifdicr")
  dir.create(motifdicr)


  f <- function(x){
    hbmcbioR::autoseed_motif_disc(x)
  }

  dtdicr <- BiocParallel::bplapply(filenamedisc, f, BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))

  write_rds(dtdicr,paste0(outdir,"/","total_matrix",".Rds"))

  for(c in 1:length(dtdicr)){
    for(s in 1:length(dtdicr[[c]])){
      write_rds(dtdicr[[c]][[s]], paste0(motifdicr,"/",c,"_",s,".Rds"))
    }
  }

  motifrawdir <- paste0(outdir,"/","motifrawdir")
  txttempdir <- paste0(outdir,"/", "txttempdir")
  dir.create(motifrawdir)
  dir.create(txttempdir)


  sum_comse <- function(con){
    nub <- con %>% str_split("") %>%.[[1]] %>% as.data.frame()
    colnames(nub) <- "base"
    nub$nuber <- 1:length(nub$base)
    for(i in 1:length(nub$base)){
      if(nub$base[i] == "A"){
        nub[i,3] <- 1
      }else if(nub$base[i] == "C"){
        nub[i,3] <- 2
      }else if(nub$base[i] == "G"){
        nub[i,3] <- 3
      }else if(nub$base[i] == "T"){
        nub[i,3] <- 4
      }else if(nub$base[i] == "N"){
        nub[i,3] <- 5
      }else if(nub$base[i] == "Y"){
        nub[i,3] <- 6
      }else if(nub$base[i] == "W"){
        nub[i,3] <- 7
      }else if(nub$base[i] == "S"){
        nub[i,3] <- 8
      }else if(nub$base[i] == "M"){
        nub[i,3] <- 9
      }else if(nub$base[i] == "K"){
        nub[i,3] <- 10
      }else if(nub$base[i] == "H"){
        nub[i,3] <- 11
      }else if(nub$base[i] == "B"){
        nub[i,3] <- 12
      }else if(nub$base[i] == "V"){
        nub[i,3] <- 13
      }else if(nub$base[i] == "D"){
        nub[i,3] <- 14
      }else if(nub$base[i] == "R"){
        nub[i,3] <- 15
      }else{
        next()
      }
    }
    colnames(nub) <- c("base","nuber","value")
    nub$multi <- nub$nuber * nub$value
    sum(nub$multi)
  }


  matrixlist <- list()
  filedressmatrix <- list.files(motifdicr)
  filenamematrix <- paste0(motifdicr, "/", filedressmatrix)

  for(v in 1:length(filedressmatrix)){
    n_name <- gsub('?.Rds','',filedressmatrix[v])
    a <- read_rds(filenamematrix[v])
    if(is.null(a)) next()
    if(is.null(nrow(a))) next()
    b <- as.matrix(a)
    colnames(b) <- NULL
    rownames(b) <- c("A", "C", "G", "T")
    motifyz <- universalmotif::convert_motifs(b)
    motifyzrc <- universalmotif::motif_rc(motifyz)
    nub <- sum_comse(motifyz["consensus"])
    nubrc <- sum_comse(motifyzrc["consensus"])
    if(nub > nubrc){
      motifyzrc["name"] <- n_name
      universalmotif::write_motifs(motifyzrc, paste0(motifrawdir, "/", n_name),overwrite = TRUE)
    }else{
      motifyz["name"] <- n_name
      universalmotif::write_motifs(motifyz, paste0(motifrawdir, "/", n_name),overwrite = TRUE)
    }
  }

  cat("2. find the unique motif","\n")
  IC_calc<-function(consensus,pseudo_freq=0.0000000000001)
  {
    letters=consensus %>% str_split("") %>% .[[1]]
    letter_freqs= letters %>% table() %>% prop.table()
    for (letter in c("A","C","G","T"))
    {
      if (is.na(letter_freqs[letter])) letter_freqs[letter]=pseudo_freq
    }
    ic=2-sum(-letter_freqs*(log2(letter_freqs)))
    ic
  }

  filedress <- list.files(motifrawdir)
  filename <- paste0(motifrawdir, "/",filedress)

  for(y in 1:length(filename)){
    motiftrim <- read_motifs(filename[y])
    if(str_detect(motiftrim["consensus"],"^N+[^ATCGRYWSMKHBVD]$")){
      cat("we remove motif", filename[y],"\n","because of repeat N consensus:","",motiftrim["consensus"],"\n")
      unlink(filename[y])
    }else if(str_detect(motiftrim["consensus"],"^N*[[RYWSMKHBVD]+N*[RYWSMKHBVD]+]+N*$")){
      cat("we remove motif", filename[y],"\n","because of unreliable consensus:","",motiftrim["consensus"],"\n")
      unlink(filename[y])
    }else if(IC_calc(motiftrim["consensus"]) > ic_threshold){
      cat("we remove motif", filename[y],"\n","because of IC too hight","",motiftrim["consensus"],"\n")
    }else{
      y=y+1
    }
  }


  for (n in 1:length(filename)){
    if(file.exists(filename[n])){
      t1 <- read_motifs(filename[n])
      t1rc <- motif_rc(t1)
      pwm_1 <- convert_motifs(t1, "TFBSTools-PWMatrix")
      pwm_1rc <- convert_motifs(t1rc, "TFBSTools-PWMatrix")
      cat("we deal with motif",filename[n],"just step one","\n")
      p <- n+1
      leqm <- length(filename)
      for (m in p:leqm){
        if(file.exists(filename[m])){
          t2 <- read_motifs(filename[m])
          pwm_2 <- convert_motifs(t2, "TFBSTools-PWMatrix")
          if(((PWMSimilarity(pwm_1, pwm_2, method="Pearson")) > threshold)){
            unlink(filename[m])
            cat("we remove motif", filename[m],"\n")
          } else if((PWMSimilarity(pwm_1rc, pwm_2, method="Pearson")) > threshold) {
            unlink(filename[m])
            cat("we remove motif", filename[m],"\n")
          }else{
            m=m+1
          }
        }
      }
    }else{
      n=n+1
    }
  }

  motifunique <- list()
  fileunique <- list.files(motifrawdir)
  filenameunique <- paste0(motifrawdir,"/",fileunique)
  for(h in 1:length(fileunique)){
    motifunique <- c(motifunique,list(read_motifs(filenameunique[h])))
  }
  motif_dt_unique <- to_df(motifunique)

  motif_dt_unique_1 <- mutate(motif_dt_unique,iclength = icscore / str_length(consensus))


  motiffilterdir <- paste0(outdir,"/","motiffilterdir")
  dir.create(motiffilterdir)
  motiffilter_dt <- filter(motif_dt_unique_1,icscore > icscore_threshold)
  motiffilter_dt_1 <- filter(motiffilter_dt,iclength > iclength_treshold)
  name_filter <- motiffilter_dt_1$name
  sapply(name_filter,function(x){file.copy(paste(motifdicr,paste0(x,".Rds"),sep="/"),motiffilterdir)})

  matrixlistlast <- list()
  filedressmatrixlast <- list.files(motiffilterdir)
  filenamematrixlast <- paste0(motiffilterdir, "/", filedressmatrixlast)
  for(j in 1:length(filedressmatrixlast)){
    matrixlistlast[[j]] <- read_rds(filenamematrixlast[j])
  }
  unlink(motiffilterdir)
  unlink(motifrawdir)
  unlink(motifdicr)
  unlink(outdir,recursive = T)
  matrixlistlast

}

library(sysbiofj)
motifs=yz_uniquemotif("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B6_18_Mock__root__conc_0nM__time_30min__30n.gz")
