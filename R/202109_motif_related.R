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
  uniquemotif_motiflist(consensus_ic_threshold=consensus_ic_threshold,total_ic_threshold = consensus_ic_threshold,total_ic_threshold_upper=total_ic_threshold_upper, iclength_treshold=iclength_treshold, similarity_corr_threshold=similarity_corr_threshold)

}

uniquemotif_motiflist<-function(motif_list,consensus_ic_threshold=1,total_ic_threshold = 3,total_ic_threshold_upper=30, iclength_treshold=0.3,similarity_corr_threshold=0.90)
{
  # convert to PFM list
  uni_PCMs=.list_to_PFMatrixL(motif_list)

  # orient the motifs
  uni_PCMs_oriented=.motifs_orient(uni_PCMs)

  # filter for low quality motifs by consensus
  PCMs_o_f=.motifs_filter_by_quality(uni_PCMs_oriented,consensus_ic_threshold=consensus_ic_threshold)

  # retain only uniq motifs
  PCMs_o_f=.motifs_keep_uniq(PCMs_o_f, similarity_corr_threshold=similarity_corr_threshold)

  # filter motifs by IC content
  motif_dt_unique <- to_df(PCMs_o_f) %>%  mutate(iclength = icscore / str_length(consensus))

  motiffilter_dt <- filter(motif_dt_unique,icscore > total_ic_threshold) %>% filter(icscore < total_ic_threshold_upper)
  motiffilter_dt <- filter(motiffilter_dt,iclength > iclength_treshold)

  # order by similarity
  motiffilter_dt[[1]] %>% .motifs_order_by_similarity()
}










# library(fjComm)
# motifs=uniquemotif_yz("/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/B6_18_Mock__root__conc_0nM__time_30min__30n.gz")
# # motifs %>% ggseqlogo_lab_list() %>% print()


#' Title
#'
#' @param motif_list
#'
#' @return
#' @export
#'
#' @examples
.list_to_PFMatrixL<-function(motif_list)
{
  library(TFBSTools)
  library(universalmotif)
  PCMs=map(motif_list, function(x){PFMatrix(profileMatrix = x %>% set_rownames(c("A","C","G","T")))} )
  uni_PCMs=map(PCMs,function(x){convert_motifs(x)})
  uni_PCMs=map(uni_PCMs,function(x){ x["name"]=x["consensus"]; x})
  uni_PCMs
}


#' Still buggy, have to run 2-3 times to get real uniq motifs
#'
#' @param uni_motif_list
#' @param similarity_corr_threshold
#'
#' @return
#' @export
#'
#' @examples
.motifs_keep_uniq<-function(uni_motif_list, similarity_corr_threshold=0.9)
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
.motifs_orient<-function(uni_motif_list){
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
.motifs_filter_by_quality<-function(uni_motif_list,consensus_ic_threshold=1){
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
.motifs_filter_by_IC<-function(uni_motif_list, total_ic_threshold = 3, total_ic_threshold_upper=30, iclength_treshold=0.3){
  motif_dt_unique <- to_df(uni_motif_list) %>%  mutate(iclength = icscore / str_length(consensus))
  motiffilter_dt <- filter(motif_dt_unique,icscore > total_ic_threshold) %>% filter(icscore < total_ic_threshold_upper)
  motiffilter_dt <- filter(motiffilter_dt,iclength > iclength_treshold)
  motiffilter_dt %>% to_list()
}


#' make a pcm to the most similar known TFs motif name
#' @param pcm  include your motif pcm
#' @param species your interested species(the JASPAR included),default is Arabidopsis thaliana
#' @param n the number of similar motif name,default is 1
#'
#' @export  list of motif_name
#'
#' @example
#' my_know <- motif_from_JP(root_motifs[[1]])
motif_to_TF <- function(PFMatrix,species='Arabidopsis thaliana',n=1){
  #to get all of the TFs PWM from JASPAR2020
  library(JASPAR2020)
  library(TFBSTools)
  pwm_library <- getMatrixSet(
    JASPAR2020,
    opts=list(
      species =  species,
      matrixtype = 'PWM'
    ))

  motiflogo <- universalmotif::convert_motifs(PFMatrix, "TFBSTools-PWMatrix")
  # find the most similar motif to our motif
  pwm_sim <- PWMSimilarity(pwm_library, motiflogo, method = 'Pearson')
  pwm_sim_rc <- PWMSimilarity(pwm_library, universalmotif::motif_rc(motiflogo), method = 'Pearson')
  pwm_sim <- pmax(pwm_sim,pwm_sim_rc)

  # extract the motif names from the pwm library
  pwm_library_list = lapply(pwm_library, function(x){
    data.frame(ID = ID(x))})

  # combine the list into one data frame
  pwm_library_dt = dplyr::bind_rows(pwm_library_list)

  # fetch the similarity of each motif to our unknown motif
  pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

  # find the most similar motif in the library
  pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
  x_n <- head(pwm_library_dt,n)
  pfm_all <- sapply(x_n$ID,function(x){TFBSTools::getMatrixByID(JASPAR2020, ID = x)})
  lapply(pfm_all, function(x)x@name)
}


#' Title
#'
#' @param uni_motif_list
#'
#' @return
#' @export
#'
#' @examples
.motifs_order_by_similarity<-function(uni_motif_list)
{
  library(universalmotif)
  dist_mat = compare_motifs(uni_motif_list,min.mean.ic = 0,method = "PCC")
  tree = pheatmap:::cluster_mat(dist_mat, distance = "euclidean", method = "complete")
  tree = pheatmap:::identity2(tree, dist_mat)

  # dist_mat=dist_mat[tree$order,tree$order]
  # pheatmap::pheatmap(dist_mat,cluster_rows = F,cluster_cols = F)

  uni_motif_list[tree$order]
}




ggseqlogo_add_motif_names<-function(motif_list_plot,TFnames)
{
  names(TFnames)=1:length(TFnames)
  motif_list_plot+facet_wrap(~seq_group,scales = "free_x",labeller = labeller(seq_group=TFnames))
}

