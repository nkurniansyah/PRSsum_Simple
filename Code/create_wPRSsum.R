
#' Create Weighted PRSsum
#'
#' @param Summary_stat1 : Summary statitsic 1, this can be based on ancestry or GWAS name 
#' @param Summary_stat2 : Summary statitsic 2, this can be based on ancestry or GWAS name
#' @param Summary_stat3 : Summary statitsic 3, this can be based on ancestry or GWAS name 
#' @param Summary_stat4 : Summary statitsic 4, this can be based on ancestry or GWAS name 
#' @param Summary_stat5 : Summary statitsic 5, this can be based on ancestry or GWAS name
#' @param n_sample_Summary_stat1 : Number of SNP to include in summary stat1. (n_sample_Summary_stat1=nrow(Summary_stat1))
#' @param n_sample_Summary_stat2 : Number of SNP to include in summary stat2. (n_sample_Summary_stat1=nrow(Summary_stat2))
#' @param n_sample_Summary_stat3 : Number of SNP to include in summary stat3. (n_sample_Summary_stat1=nrow(Summary_stat3))
#' @param n_sample_Summary_stat4 : Number of SNP to include in summary stat4. (n_sample_Summary_stat1=nrow(Summary_stat4))
#' @param n_sample_Summary_stat5 : Number of SNP to include in summary stat5. (n_sample_Summary_stat1=nrow(Summary_stat5))
#' @param scaling_file : This is based on TOPmed Scaling, make sure the study name match with summary statistics.
#                         e.g: if your summary Summary_stat1= FINNGEN, in scaling file.csv, 
#                              in first row of scaling in Study column should Summary_stat1 mean and SD of FINNGEN.
#                              see example file
#'
#' @param weight_file : This is based on MGB Weight, make sure the study name match with summary statistics.
#                         e.g: if your summary Summary_stat1= FINNGEN, in weight file.csv, 
#                              in first row of weight in Study column should Summary_stat1 prs_effect FINNGEN in MGB
#                              see example file formating
#'
#' @return  (data frmae of PRssum)
#' @export
#'
#' @examples
#' 
create_wprsum<-function(Summary_stat1, Summary_stat2=NA,Summary_stat3=NA,
                       Summary_stat4=NA, Summary_stat5=NA, n_sample_Summary_stat1, 
                       n_sample_Summary_stat2=NA, n_sample_Summary_stat3=NA,
                       n_sample_Summary_stat4=NA, n_sample_Summary_stat5=NA,
                       scaling_file, weight_file=NA){
  
  
  list_n_samples<-list(n_sample_Summary_stat1=n_sample_Summary_stat1, n_sample_Summary_stat2=n_sample_Summary_stat2,
                       n_sample_Summary_stat3=n_sample_Summary_stat3, n_sample_Summary_stat4=n_sample_Summary_stat4,
                       n_sample_Summary_stat5=n_sample_Summary_stat5)
  
  
  
  list_summary_stat<-list(Summary_stat1=Summary_stat1,Summary_stat2=Summary_stat2, Summary_stat3=Summary_stat3, 
                          Summary_stat4=Summary_stat4,
                          Summary_stat5=Summary_stat5)
  
  
  
  scaling<- read.csv(scaling_file)
  weight<-read.csv(weight_file)
  list_summary_stat<-list_summary_stat[!is.na(list_summary_stat)]

  #makes sure your study name is Summary_stat1,Summary_stat2, ...
  # make sure the name weight is prs_effect (Follow the example)
  study_names<- names(list_summary_stat)
  
  out_beta<-list()
  out_snp_include<-list()
  for(study in study_names){
    summary_stat_df<- fread(list_summary_stat[[study]], data.table = F)
    
    
    list_snp<-summary_stat_df %>% mutate(SNP=paste0(CHR,":",POS_hg38,":",A1,":",A2)) %>% dplyr::select(rsID, SNP)
    head(summary_stat_df)
    
    scaling_study<-scaling[which(scaling$Study==study),]
    mean_study<- as.numeric(scaling_study$Mean)
    sd_study<- as.numeric(scaling_study$SD)
    
    weight_df<-weight[which(weight$Study==study),]
    weight_study<-as.numeric(weight_df$prs_effect)
    n_study<-list_n_samples[[paste0("n_sample_",study)]]
    beta<- data.frame(rsID=summary_stat_df$rsID,
                      beta_to_add=((weight_study/sd_study)*summary_stat_df$BETA)/(2*n_study),
                      beta_to_substract=((weight_study * mean_study) /sd_study) /(2*n_study))
    head(beta)
    
    colnames(beta)<-c("rsID",paste0(study,"_add"), paste0(study,"_sub"))
    out_snp_include[[study]]<- list_snp
    out_beta[[study]]<- beta
    
  }
 
  
  snp_to_include<- do.call(rbind, out_snp_include)
  beta_clean<- purrr::reduce(out_beta, full_join, "rsID")
  head(beta_clean)
  # repalce NA with 0
  
  beta_clean[is.na(beta_clean)]<- 0
  

  beta_scaling_add<- rowSums(beta_clean[, grepl("_add", colnames(beta_clean))])
  beta_scaling_substract<- rowSums(beta_clean[, grepl("_sub", colnames(beta_clean))])
  
  final_beta<- data.frame(rsID=beta_clean$rsID, BETA=beta_scaling_add-beta_scaling_substract)
  head(final_beta)
  dim(final_beta)
  head(snp_to_include)
  
  snp_to_include<-snp_to_include[!duplicated(snp_to_include$SNP),]
  dim(snp_to_include)
  
  final_beta_annot<- left_join(snp_to_include, final_beta, by="rsID")
  head(final_beta_annot)
  
  final_beta_annot<- final_beta_annot %>% separate(SNP, c("CHR","POS","A1","A2"),":")%>% 
                                          dplyr::select(CHR,rsID,POS,A1,A2,BETA)
  
  return(final_beta_annot)
  
}


