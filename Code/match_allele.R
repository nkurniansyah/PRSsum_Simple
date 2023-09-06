
#' Title : Matching allele accross study
#'
#' @param refrence_snp : A data frame for SNP reference data, which can be derived from a full summary statistics dataset. This file should include the following columns: 
#'                       Chromosome, Position, Allele1_ref, and Allele2_ref. Ensure that the naming conventions are consistently followed.
#' @param list_variants_weight : list data frame of variants weight, and each study must have following columns:
#'                       Chromosome, Position, Allele1, and Allele2, BETA
#' @param match_by_position : If match_by_position is set to TRUE, the reference dataset will be merged with the variant weights dataset using the Chromosome and Position columns.
#' @return a data frame of all varaints weight with annotation from the reference SNP
#' 
#' @export
#'
#' @examples
#' 
match_allele<- function(refrence_snp, 
                        list_variants_weight,
                        match_by_position=TRUE){
  
  if(!is.list(list_variants_weight)){
    message("variants weight have to be in list")
  }
  
  study_names<-names(study_list)
  out<-list()
  for(study in study_names){
    selected_variant_weight<-list_variants_weight[[study]] 

    if(match_by_position){
      message("variant weight will match using chromosome and position")
      
      variant_weight_annot<-left_join(selected_variant_weight,refrence_snp,by=c("Chromosome","Position"))
      ind_flip_beta<- which(variant_weight_annot$Allele1_ref!=variant_weight_annot$Allele1)
      
      if(length(ind_flip_beta)> 0){
        message(paste0("fliping allele with total ", length(ind_flip_beta)))
        variant_weight_annot$BETA[ind_flip_beta]<- variant_weight_annot$BETA[ind_flip_beta]*-1
      }
      
      ind_allele_na<-which(is.na(variant_weight_annot$Allele1_ref))
      if(length(ind_allele_na)> 0){
        variant_weight_annot$Allele1_ref<-variant_weight_annot$Allele1
        variant_weight_annot$Allele2_ref<-variant_weight_annot$Allele2
      }
      
      variant_weight_annot<-variant_weight_annot %>% dplyr::select(Chromosome, Position, Allele1_ref, Allele2_ref )
      head(variant_weight_annot)
      
      variant_weight_df<- variant_weight_annot %>% mutate(variant_id= paste0(Chromosome,":",Position,":",Allele1_ref,":",Allele2_ref))
      head(variant_weight_df)
      variant_weight_df<-variant_weight_df %>% dplyr::select(variant_id, BETA)
      colnames(variant_weight_df)<-c("variant_id",study )
      
    }else{
      message("variant weight will match using rsID")
      
      variant_weight_annot<-left_join(selected_variant_weight,refrence_snp,by=c("rsID"))
      ind_flip_beta<- which(variant_weight_annot$Allele1_ref!=variant_weight_annot$Allele1)
      
      if(length(ind_flip_beta)> 0){
        message(paste0("fliping allele with total ", length(ind_flip_beta)))
        variant_weight_annot$BETA[ind_flip_beta]<- variant_weight_annot$BETA[ind_flip_beta]*-1
      }
      
      ind_allele_na<-which(is.na(variant_weight_annot$Allele1_ref))
      if(length(ind_allele_na)> 0){
        variant_weight_annot$Allele1_ref<-variant_weight_annot$Allele1
        variant_weight_annot$Allele2_ref<-variant_weight_annot$Allele2
      }
      variant_weight_annot<-variant_weight_annot %>% dplyr::select(rsID,Chromosome_ref, Position_ref, Allele1_ref, Allele2_ref,BETA )
      head(variant_weight_annot)
      
      variant_weight_df<- variant_weight_annot %>% mutate(variant_id= paste0(Chromosome_ref,":",Position_ref,":",Allele1_ref,":",Allele2_ref))
      head(variant_weight_df)
      variant_weight_df<-variant_weight_df %>% dplyr::select(variant_id,rsID, BETA)
      colnames(variant_weight_df)<-c("variant_id","rsID",study )
    }
    
    out[[study]]<- variant_weight_df
  }
  
  all_variant_weight<- purrr::reduce(out, full_join, by="variant_id")
  head(all_variant_weight)
  
  all_variant_weight<- all_variant_weight %>%separate(variant_id, c("Chromosome","Position","Allele1","Allele2"),":")
  
  return(all_variant_weight)
  
}