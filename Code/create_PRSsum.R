
#' Title: Create PRSsum
#'
#' @param variant_weights: A data frame that contains information on chromosome (chr), position (pos), allele, 
#'                         and the corresponding weights from various studies.
#'                         NOTE: The names of the variants' weights must match the study names in the PRSsum scaling file. 
#'                         For example, if you obtain variant weights from both MVP and FINNGEN, their names must correspond accordingly.
#' @param PRSsum_scaling: A data frame that includes the study names matching the column names in "variants_weights," along with the Mean and standard deviation (SD) 
#'                        from the training dataset, as well as the number of variants selected to generate PRS for each study.
#' @param weight_file : A path file containing information on weights will be implemented in 'variants_weight' for each study. 
#'                      The column names must match with 'prs_effect' and the study name in the same column as 'variants_weight' and 'PRSsum_scaling'
#' @param chr_col_name : Chromosome column name in "variant_weights"
#' @param pos_col_name : Position column name in "variant_weights"
#' @param effect_allele_col_name : Effect allele column name in "variant_weights"
#' @param other_allele_col_name : Other allele column name in "variant_weights" 
#' @param rsID_col_name : rsID column name in "variant_weights"  if available. leave it NULL if its not available
#'
#' @return A data frame of annotation (rsID, chormosome, position, effect allele and othe allele) and variant weight of PRSsum.
#' @export
#'
#' @examples
create_prsum<-function(variant_weights, 
                       PRSsum_scaling, 
                       weight_file=NULL,
                       chr_col_name, 
                       pos_col_name, 
                       effect_allele_col_name,
                       other_allele_col_name,
                       rsID_col_name=NULL){
  
  variant_weights[is.na(variant_weights)] <- 0
  
  if(!is.null(weight_file)){
    message("Run weighted PRSsum")
    
    weight<- fread(weight_file, data.table=F)
    PRSsum_scaling<-left_join(PRSsum_scaling, weight, by="Study")
    PRSsum_scaling$Scaling <- (PRSsum_scaling$prs_effect/PRSsum_scaling$SD)/(2*PRSsum_scaling$N_variants)
  }else{
    message("Run unweighted PRSsum")
    PRSsum_scaling$Scaling <- (1/PRSsum_scaling$SD)/(2*PRSsum_scaling$N_variants)
    
  }
  
  
  variant_weights[,PRSsum_scaling$Study]<- (mapply("*",as.data.frame(variant_weights[,PRSsum_scaling$Study]),as.numeric(PRSsum_scaling$Scaling)))
  
  head(variant_weights)
  
  variant_weights$PRSsum<- rowSums(variant_weights[,PRSsum_scaling$Study])
  
  head(variant_weights)
  
  if(!is.null(rsID_col_name)){
    variant_weights_clean<- variant_weights[,c(rsID_col_name,chr_col_name,pos_col_name,effect_allele_col_name, other_allele_col_name,"PRSsum")]
    colnames(variant_weights_clean)<-c("rsID","chr_name","chr_position","effect_allele","other_allele","effect_weight")
  }else{
    variant_weights_clean<- variant_weights[,c(chr_col_name,pos_col_name,effect_allele_col_name, other_allele_col_name,"PRSsum")]
    colnames(variant_weights_clean)<-c("chr_name","chr_position","effect_allele","other_allele","effect_weight")
    
  }
  

  return(variant_weights_clean)
  
}

