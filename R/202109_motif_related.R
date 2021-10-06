# convert a motif list to PFMatrix
motif_list_to_PFMatrix<- function(motif_list){map(motif_list, function(x){PFMatrix(profileMatrix = x %>% set_rownames(c("A","C","G","T")))} )}


# aux function to help orienting the motifs
.sum_comse <- function(con){
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

.IC_calc<-function(consensus,pseudo_freq=0.0000000000001)
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

#' For input file(s),
#' 1. De novo motif discovery
#' 2. filtering for the "good" motifs
#' 3. retaining only 1 copy of similar motifs
#'
#'
#' @param files  the target sequence file(s) for motif discovery, can be fastq, fasta, or plain text, either compressed as .gz or not
#' @param outdir the output folder, if specified, to store all motifs
#' @param consensus_ic_threshold to filter out motifs whose consensus has too high ic, i.e., lack complexity
#' @param total_ic_threshold_upper to filter out motifs whose total ic is too high, i.e., motif from fixed sequences
#' @param total_ic_threshold to filter out motifs whose total ic is too low, i.e., signal strength is weak
#' @param iclength_treshold to filter out motifs whose ic per base is too low, i.e., signal strength is weak
#' @param similarity_corr_threshold to retain only 1 copy of duplicating motifs
#' @param thread the number of threads
#'
#'
#' @return
#' a list of
#' @export
#'
#' @examples
#' uniquemotif("/wrk/wenchenjin/work/Chenjin1_3__ATI_30N101N/pre_processed/cleaned","/wrk/yuanzhen")
#'
#'
uniquemotif_files <- function(files,outdir=NULL,consensus_ic_threshold=1,total_ic_threshold = 3,total_ic_threshold_upper=30, iclength_treshold=0.3,similarity_corr_threshold=0.90,thread = 30){
  library(TFBSTools)
  library(universalmotif)

  ###Now we will the first steps to de novo motif discovry
  filenamedisc <- Sys.glob(files)
    f <- function(x){
      hbmcbioR::autoseed_motif_disc(x)
    }

    rdsfile=paste0(outdir,"/","total_matrix",".Rds")
  if (!file.exists(rdsfile)){
    dtdicr <- BiocParallel::bplapply(filenamedisc, f, BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))
    if(!is.null(outdir)) {dir.create(outdir,recursive = T); write_rds(dtdicr,rdsfile)}
  }else{
    dtdicr <-readRDS(rdsfile)
  }


  motifs=do.call(c,dtdicr)
  PCMs=map(motifs, function(x){PFMatrix(profileMatrix = x %>% set_rownames(c("A","C","G","T")))} )
  uni_PCMs=map(PCMs,function(x){convert_motifs(x)})
  uni_PCMs=map(uni_PCMs,function(x){ x["name"]=x["consensus"]; x})

    uni_rcPCMs=map(uni_PCMs,function(x){motif_rc(x)})

  # orient the motifs
  uni_PCMs_oriented=map2(uni_PCMs,uni_rcPCMs, function(pcm,rcpcm){
      if(.sum_comse(pcm["consensus"])> .sum_comse(rcpcm["consensus"])){
        return(rcpcm)}else{return(pcm)}
    })


  # filter for low quality motifs by consensus
  PCMs_o_f=map(uni_PCMs_oriented,function(x){
    if(str_detect(x["consensus"],"^N+[^ATCGRYWSMKHBVD]$"))return(NULL)  #repeat N consensus
    if(str_detect(x["consensus"],"^N*[[RYWSMKHBVD]+N*[RYWSMKHBVD]+]+N*$"))return(NULL) #unreliable consensus
    if(.IC_calc(x["consensus"]) > consensus_ic_threshold)return(NULL) #consensus lacks variety
    x
  })
      null_filter=map(PCMs_o_f,is.null) %>% unlist()
      PCMs_o_f[null_filter %>% which]=NULL

  # retain only uniq motifs
  for (n in 1:length(PCMs_o_f)){
      if(is.null(PCMs_o_f[n] %>% unlist()))next
      t1 <- PCMs_o_f[[n]]
      t1rc <- motif_rc(t1)
      pwm_1 <- convert_motifs(t1, "TFBSTools-PWMatrix")
      pwm_1rc <- convert_motifs(t1rc, "TFBSTools-PWMatrix")
      # print("n=" %>% paste0(n))
      p <- n+1
      leqm <- length(PCMs_o_f)
      for (m in p:leqm){
          if(is.null(PCMs_o_f[m] %>% unlist()))next
          t2 <- PCMs_o_f[[m]]
          pwm_2 <- convert_motifs(t2, "TFBSTools-PWMatrix")
          if(((PWMSimilarity(pwm_1, pwm_2, method="Pearson")) > similarity_corr_threshold)){
            if (PCMs_o_f[[m]]["icscore"] < PCMs_o_f[[n]]["icscore"]){ PCMs_o_f[[m]]=NULL } else { PCMs_o_f[[n]]=NULL }
          } else if((PWMSimilarity(pwm_1rc, pwm_2, method="Pearson")) > similarity_corr_threshold) {
            if (PCMs_o_f[[m]]["icscore"] < PCMs_o_f[[n]]["icscore"]){ PCMs_o_f[[m]]=NULL } else { PCMs_o_f[[n]]=NULL }
          }
      }
  }

  # filter motifs by IC content
  motif_dt_unique <- to_df(PCMs_o_f) %>%  mutate(iclength = icscore / str_length(consensus))

  motiffilter_dt <- filter(motif_dt_unique,icscore > total_ic_threshold) %>% filter(icscore < total_ic_threshold_upper)
  motiffilter_dt <- filter(motiffilter_dt,iclength > iclength_treshold)

  motiffilter_dt[[1]]
}

# library(fjComm)
# motifs=uniquemotif_yz("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B6_18_Mock__root__conc_0nM__time_30min__30n.gz")
# # motifs %>% ggseqlogo_lab_list() %>% print()



#' Still buggy, have to run 2-3 times to get real uniq motifs
#'
#' @param uni_motif_list
#' @param similarity_corr_threshold
#'
#' @return
#' @export
#'
#' @examples
motifs_keep_uniq<-function(uni_motif_list, similarity_corr_threshold=0.9)
{
  # retain only uniq motifs
  for (n in 1:length(uni_motif_list)){
    # if(is.null(uni_motif_list[n] %>% unlist()))next
    if(is.na(uni_motif_list[n]))next
    t1 <- uni_motif_list[[n]]
    t1rc <- motif_rc(t1)
    pwm_1 <- convert_motifs(t1, "TFBSTools-PWMatrix")
    pwm_1rc <- convert_motifs(t1rc, "TFBSTools-PWMatrix")
    p <- n+1
    for (m in p:length(uni_motif_list)){
      if(is.null(uni_motif_list[m] %>% unlist()))next
      if(is.na(uni_motif_list[m]))next
      if(m==n)next
      t2 <- uni_motif_list[[m]]
      pwm_2 <- convert_motifs(t2, "TFBSTools-PWMatrix")
      if(((PWMSimilarity(pwm_1, pwm_2, method="Pearson")) > similarity_corr_threshold)){
        if (uni_motif_list[[m]]["icscore"] < uni_motif_list[[n]]["icscore"]){ uni_motif_list[[m]]=NA} else { uni_motif_list[[n]]=uni_motif_list[[m]]; uni_motif_list[[m]]=NA}
      } else if((PWMSimilarity(pwm_1rc, pwm_2, method="Pearson")) > similarity_corr_threshold) {
        if (uni_motif_list[[m]]["icscore"] < uni_motif_list[[n]]["icscore"]){ uni_motif_list[[m]]=NA} else { uni_motif_list[[n]]=uni_motif_list[[m]]; uni_motif_list[[m]]=NA}
      }
    }
  }
  uni_motif_list[is.na(uni_motif_list)]=NULL
  uni_motif_list
}


#' orient the motifs
#'
#' @param uni_motif_list
#'
#' @return
#' @export
#'
#' @examples
motifs_orient<-function(uni_motif_list){
  uni_PCMs=uni_motif_list
  uni_rcPCMs=map(uni_PCMs,function(x){motif_rc(x)})
  map2(uni_PCMs,uni_rcPCMs, function(pcm,rcpcm){
    if(.sum_comse(pcm["consensus"])> .sum_comse(rcpcm["consensus"])){
      return(rcpcm)}else{return(pcm)}
  })
}




#' Filter for low quality motifs by consensus
#'
#' @param uni_motif_list
#' @param consensus_ic_threshold
#'
#' @return
#' @export
#'
#' @examples
motifs_filter_by_quality<-function(uni_motif_list,consensus_ic_threshold=1){
  uni_motif_list=map(uni_motif_list,function(x){
    if(str_detect(x["consensus"],"^N+[^ATCGRYWSMKHBVD]$"))return(NULL)  #repeat N consensus
    if(str_detect(x["consensus"],"^N*[[RYWSMKHBVD]+N*[RYWSMKHBVD]+]+N*$"))return(NULL) #unreliable consensus
    if(.IC_calc(x["consensus"]) > consensus_ic_threshold)return(NULL) #consensus lacks variety
    x
  })
  null_filter=map(uni_motif_list,is.null) %>% unlist()
  uni_motif_list[null_filter %>% which]=NULL
  uni_motif_list
}




#' filter motifs by IC content
#'
#' @param uni_motif_list
#' @param total_ic_threshold
#' @param total_ic_threshold_upper
#' @param iclength_treshold
#'
#' @return
#' @export
#'
#' @examples
motifs_filter_by_IC<-function(uni_motif_list, total_ic_threshold = 3, total_ic_threshold_upper=30, iclength_treshold=0.3){
  motif_dt_unique <- to_df(uni_motif_list) %>%  mutate(iclength = icscore / str_length(consensus))
  motiffilter_dt <- filter(motif_dt_unique,icscore > total_ic_threshold) %>% filter(icscore < total_ic_threshold_upper)
  motiffilter_dt <- filter(motiffilter_dt,iclength > iclength_treshold)
  motiffilter_dt %>% to_list()
}


