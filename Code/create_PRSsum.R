
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
    message("Runweighted PRSsum")
    
    weight<- fread(weight_file, data.table=F)
    PRSsum_scaling<-left_join(PRSsum_scaling, variant_weights, by="Study")
    PRSsum_scaling$Scaling <- (PRSsum_scaling$prs_effect/PRSsum_scaling$SD)/(2*PRSsum_scaling$N)
  }else{
    message("Run unweighted PRSsum")
    PRSsum_scaling$Scaling <- (1/PRSsum_scaling$SD)/(2*PRSsum_scaling$N)
    
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
  
  head(variant_weights_clean)
  
  
}

