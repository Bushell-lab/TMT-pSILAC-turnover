library(dplyr)
library(tidyverse)
library(limma)
library(minpack.lm)
library(BSDA)


# loading variables
home_folder = "~/data/TSCHMIDT/backup/temp data/pSILAC data/R analysis/"
treatment = "a2ko"


#load SILAC file with header and row names, link to path of file
master = read.table(moore_pardo.txt",sep=""),header=TRUE,sep="\t")
row.names(master)=make.names(master$gene_name, unique=TRUE)
master = master %>% select(-contains("POOL"),-contains("gene_name")) %>% drop_na() # remove empty lanes
master = 2^master #unlog data

#load sample sheet, link to path of file
ss = read.table(samplesheet.txt",sep=""),sep="\t") 
#### assuming that the order of columns in the result file is t0|2h control|2h treatment|4h...|

# generate vector with new order columns to replicates and groups next to each other ####
exp_groups = 7 #number of sample groups, i.e. 0h, and 2 4 8 for each dmso and hipp
replicates = 4 # number of replicates
new_order = 1
index=0
for (replicate_order in 1:exp_groups) {
  for (replicate_column in 1:replicates)   {
    index = index + 1
    new_order[index] = replicate_order + (replicate_column-1)*7
  }
}

# split data into tables for Fraction (turnover), heavy intensities (synthesis) and light (degradation) intensities only, 
# pre-filters degradation has to shown downwards trend over time,
# and synthesis has to shown upwards trend over time

for (int in c("heavy","light","fraction")){
  if (int == "fraction")
  {
    
    master_heavy = subset(master_heavy, row.names(master_heavy)%in%row.names(master_light))
    master_light = subset(master_light, row.names(master_light)%in%row.names(master_heavy))
    master_heavy =na.omit(master_heavy)
    master_light =na.omit(master_light)
    master_fraction = (na.omit(master_heavy / (master_light + master_heavy)))
  }
  else
  {
    df = (select(master, contains(int), -contains("pool")))
    df = df[,new_order]
    names(df)=ss[1:28,"V1"]
    df= df %>% mutate(mean_t0 = rowMeans(across(contains("t0")&contains(c("r1","r3","r4")))),
                      mean_t24 = rowMeans(across(contains("t24")&contains("c"))),
                      mean_t24_treat = rowMeans(across(contains("t24")&contains("a2ko")&contains(c("r2","r3","r4")))))
    if (int == "light") {df = df %>% filter((mean_t24-mean_t0)<0 & (mean_t24_treat-mean_t0)<0) %>% select(-contains("mean"))}
    else {df = df %>% filter((mean_t24-mean_t0)>0 & (mean_t24_treat-mean_t0)>0 ) %>% select(-contains("mean"))}
    assign(paste0("master_",int), df)
  }
}



#normalises degradation and synthesis tables to first and last time point respectively
for (i in 1:4){
  columns = c(0,4,8,12,16,20,24)+i
  columns_hNT = c(0,4,12,20)+i
  columns_hA2 =c(0,8,16,24)+i
  master_light[,(columns)] = master_light[,(columns)]/ master_light[,i]
  master_heavy[,(columns_hNT)] = master_heavy[,(columns_hNT)]/ master_heavy[,(20+i)]
  master_heavy[,(columns_hA2)] = master_heavy[,(columns_hA2)]/ master_heavy[,(24+i)]
}

# batch correction using LIMMA
batch <- rep(c("A","B","C","D"),7)
master_light <- removeBatchEffect(master_light, batch) %>% as.data.frame()
master_heavy <- removeBatchEffect(master_heavy, batch) %>% as.data.frame()
master_fraction <- removeBatchEffect(master_fraction, batch) %>% as.data.frame()


#split data into DMSO and TREATMENT for TURNOVER table
master_list = list(master_fraction=master_fraction)
master_list_split = list()
counter=1
for (int in 1){
  df=data.frame(t(data.frame(master_list[[int]])))
  df$timepoint = ss[c(1:28),"V5"]
  df=data.frame(t(df))
  assign(paste0(names(master_list[int]),"_",treatment),data.frame(t(select(df, contains("a2ko"),contains("t0")))))
  assign(paste0(names(master_list[int]),"_dmso"),data.frame(t(select(df, contains("c"),contains("t0")))))
  master_list_split[[counter]] = data.frame(t(select(df, contains("c"),contains("t0"))))
  names(master_list_split)[counter]=paste0(names(master_list[int]),"_dmso")
  counter=counter+1
  master_list_split[[counter]] = data.frame(t(select(df,contains("a2ko"),contains("t0"))))
  names(master_list_split)[counter]=paste0(names(master_list[int]),"_",treatment)
  counter=counter+1
  
}

#filters for steadiness in TURNOVER, ie exponential growth
#split data (above) before use this part
for (int in 1:length(master_list_split)){
  df=master_list_split[[int]]
  df2=df
  for (col_number in 1:(ncol(df)-1))
  {
    stats_temp = summary(lm(df[1:(nrow(df)-4),col_number] ~ df[1:(nrow(df)-4),"timepoint"],df[1:(nrow(df)-4),]))
    stats_temp_2 = summary(lm(df[1:(nrow(df)),col_number] ~ df[1:(nrow(df)),"timepoint"],df[1:(nrow(df)),]))
    
    if (  stats_temp$coefficients[2,4]>0.1 |stats_temp_2$coefficients[2,4]>0.3  | stats_temp$coefficients[2,1]<0 | stats_temp_2$coefficients[2,1]<0   )
    {
      df2=df2 %>% select(-contains(names(df)[col_number]))
    }
  }
  master_list_split[[int]] = df2
}

# calculates rate of TURNOVER,  p value, log2 and fdr #####
# filter data as above before running this

for (int in 1:length(master_list_split)){
  fraction_slopes_stats=data.frame()
  df=master_list_split[[int]]
  for (col_number in 1:(ncol(df)-1))
  {
    gene_lm=names(df[col_number]) 
    DMSO_lm =  nlsLM(data=df, df[,col_number] ~ a + (b-a)*(1-exp(-c*timepoint)) ,start=list(a=mean(df[13:16,col_number])*0.9,b=mean(df[1:12,col_number])*0.9,c=(log(mean(df[1:4,col_number]))-log(mean(df[13:16,col_number])))/8*0.9  ), na.action = "na.exclude", control=nls.lm.control(maxiter = 1000))
    fraction_slopes_stats_temp = as.data.frame(summary(DMSO_lm)$coefficients)
    fraction_slopes_stats_temp$gene = gene_lm
    fraction_slopes_stats= rbind(fraction_slopes_stats,fraction_slopes_stats_temp["c",c(1,2,4,5)])
    row.names(fraction_slopes_stats)[nrow(fraction_slopes_stats)]=gene_lm
   
 
    if (col_number == ncol(df)-1)
    {names(fraction_slopes_stats)=c("rate","std_error","pvalue","gene_name")}
  }
  
  assign(paste0(names(master_list_split[int]),"_slope"),fraction_slopes_stats)
  if (int%%2==0)
  {
    df = merge(get(paste0(names(master_list_split[int-1]),"_slope")),get(paste0(names(master_list_split[int]),"_slope")), by.x="gene_name",by.y="gene_name",suffixes = c(".dmso",paste0(".",treatment)))
    df$log2FC = log(df[,paste0("rate.",treatment)]/df$rate.dmso,2)
    rate_test = tsum.test(mean.x=df[,"rate.dmso"],   s.x=df[,"std_error.dmso"], n.x=replicates,
                           mean.y=df[,paste0("rate.",treatment)], s.y=df[,paste0("std_error.",treatment)], n.y=replicates)
    df[,"pvalue_delta"] = rate_test$p.value
    df=df[order(df$pvalue_delta, decreasing = F),] 
    
    df$fdr_delta=p.adjust(as.numeric(df$pvalue_delta), method="fdr")
    for (row_number in 1:nrow(df))
    {
      if (df[row_number,"pvalue_delta"]<0.05)
      {df[row_number,"sig"] = "eIF4A2-dependent"}
      else if (df[row_number,"pvalue_delta"]>0.7)
      {df[row_number,"sig"] = "eIF4A2-independent"}
      else {df[row_number,"sig"] = "non-sig"}
    }
    df$sig = factor(df$sig, levels=c("eIF4A2-independent","non-sig","eIF4A2-dependent"))
    
    
    assign(paste0(names(master_list[int/2]),"_complete"), df)
    df2=subset(df,df$pvalue.dmso<0.1&df[,paste0("pvalue.",treatment)]<0.1)
    assign(paste0(names(master_list[int/2]),"_complete_sig"), df2)
    assign(paste0(treatment),df2)
  }
}

#adds half life to the table
a2ko$halflife_nt = log(2)/a2ko$rate.dmso
a2ko$halflife_a2ko = log(2)/a2ko$rate.a2ko
a2ko$log2FC_halflife= log(a2ko$halflife_a2ko/a2ko$halflife_nt,2)
master_fraction=master_fraction %>% filter( row.names(master_fraction)%in%names(master_list_split[[1]]) )





# split data into DMSO and TREATMENT for DEGRADATION and SYNTHESIS 


master_heavy= master_heavy %>%  
  mutate(diff_t8_t16_c = (rowMeans(across(contains("t16")&contains("c"))) - rowMeans(across(contains("t8")&contains("c"))))/rowMeans(across(contains("t8")&contains("c"))),
         diff_t8_t16_treat = (rowMeans(across(contains("t16")&contains("a2ko"))) - rowMeans(across(contains("t8")&contains("a2ko"))))/rowMeans(across(contains("t8")&contains("a2ko")))    ) %>%
  filter((diff_t8_t16_c >(0.1) &(diff_t8_t16_treat)>(0.1) )) %>% 
  select(-contains("diff"))

master_heavy = master_heavy %>%
  mutate(diff_t16_t24_c = (rowMeans(across(contains("t24")&contains("c"))) - rowMeans(across(contains("t16")&contains("c")))),
         diff_t16_t24_treat = (rowMeans(across(contains("t24")&contains("a2ko"))) - rowMeans(across(contains("t16")&contains("a2ko"))))    ) %>%
  filter((diff_t16_t24_c >(0) &(diff_t16_t24_treat)>(0) )) %>% 
  select(-contains("diff"))


master_list = list(master_light=master_light,
                   master_heavy=master_heavy)
master_list_split = list()
counter=1
for (int in 1:length(master_list)){
  df=data.frame(t(data.frame(master_list[[int]])))
  df$timepoint = ss[c(1:28),"V5"]
  df=data.frame(t(df))
  assign(paste0(names(master_list[int]),"_",treatment),data.frame(t(select(df, contains("a2ko"),contains("t0")))))
  assign(paste0(names(master_list[int]),"_dmso"),data.frame(t(select(df, contains("c"),contains("t0")))))
  master_list_split[[counter]] = data.frame(t(select(df, contains("c"),contains("t0"))))
  names(master_list_split)[counter]=paste0(names(master_list[int]),"_dmso")
  counter=counter+1
  master_list_split[[counter]] = data.frame(t(select(df,contains("a2ko"),contains("t0"))))
  names(master_list_split)[counter]=paste0(names(master_list[int]),"_",treatment)
  counter=counter+1
  
}

# filters for steadiness in DEGRADATION, ie exponential decay
# filters for steadiness in SYNTHESIS, ie exponential growth

for (int in 1:4){
  df=master_list_split[[int]]
  df2=df
  for (col_number in 1:(ncol(df)-1))
  {
    stats_temp = summary(lm(df[1:(nrow(df)-4),col_number] ~ df[1:(nrow(df)-4),"timepoint"],df[1:(nrow(df)-4),]))
    stats_temp_2 = summary(lm(df[1:(nrow(df)),col_number] ~ df[1:(nrow(df)),"timepoint"],df[1:(nrow(df)),]))
    
  if (  (stats_temp$coefficients[2,4]>0.4 |stats_temp_2$coefficients[2,4]>0.4)  | ((stats_temp$coefficients[2,1]>0 | stats_temp_2$coefficients[2,1]>0)  & (int==1 | int==2))  |  ((stats_temp$coefficients[2,1]<0 | stats_temp_2$coefficients[2,1]<0)  & (int==3 | int==4)) )
    {
    df2=df2 %>% select(-contains(names(df)[col_number]))
    }
  }
  master_list_split[[int]] = df2
}
  


# calulates  rates for SYNTHESIS and DEGRADATION,  p value, log2 and fdr #####

for (int in 1:4) {
  fraction_slopes_stats=data.frame()
  df=master_list_split[[int]]
  for (col_number in 1:(ncol(df)-1))
  {
    if   (grepl("light",names(master_list_split)[int],fixed=T) == T) 
    {
      start_list=list(a=(mean(df[13:16,col_number]))*0.9,b=(mean(df[9:12,col_number]))*0.8,c=abs((log(mean(df[1:4,col_number]))-log(mean(df[13:16,col_number])))/8*0.8)  )
      my_formula =  formula("df[,col_number] ~ a + (b-a)*(1-exp(-c*timepoint))")
    }
    else
    {
      start_list=list(a=mean(df[9:12,col_number])*0.9,b=(mean(df[13:16,col_number]))*0.8,c=abs((log(mean(df[1:4,col_number]))-log(mean(df[13:16,col_number])))/8*0.71)  )
      my_formula =  formula("df[,col_number] ~ b + (a-b)*(1-exp(-c*timepoint))")
    }
    gene_lm=names(df[col_number]) 
    DMSO_lm =  nlsLM(data=df, my_formula  ,start=start_list, na.action = "na.exclude", control=nls.lm.control(maxiter = 1000))
    
    fraction_slopes_stats_temp = as.data.frame(summary(DMSO_lm)$coefficients)
    fraction_slopes_stats_temp$gene = gene_lm
    fraction_slopes_stats= rbind(fraction_slopes_stats,fraction_slopes_stats_temp["c",c(1,2,4,5)])
    row.names(fraction_slopes_stats)[nrow(fraction_slopes_stats)]=gene_lm
    
    
    if (col_number == ncol(df)-1)
    {names(fraction_slopes_stats)=c("rate","std_error","pvalue","gene_name")}
  }
  
  assign(paste0(names(master_list_split[int]),"_slope"),fraction_slopes_stats)
  if (int%%2==0)
  {
    df = merge(get(paste0(names(master_list_split[int-1]),"_slope")),get(paste0(names(master_list_split[int]),"_slope")), by.x="gene_name",by.y="gene_name",suffixes = c(".dmso",paste0(".",treatment)))
    df$log2FC = log(df[,paste0("rate.",treatment)]/df$rate.dmso,2)
    rate_test = tsum.test(mean.x=df[,"rate.dmso"],   s.x=df[,"std_error.dmso"], n.x=replicates,
                          mean.y=df[,paste0("rate.",treatment)], s.y=df[,paste0("std_error.",treatment)], n.y=replicates)
    df[,"pvalue_delta"] = rate_test$p.value
    df=df[order(df$pvalue_delta, decreasing = F),] 
    
    df$fdr_delta=p.adjust(as.numeric(df$pvalue_delta), method="fdr")
    for (row_number in 1:nrow(df))
    {
      if (df[row_number,"pvalue_delta"]<0.05)
      {df[row_number,"sig"] = "eIF4A2-dependent"}
      else if (df[row_number,"pvalue_delta"]>0.7)
      {df[row_number,"sig"] = "eIF4A2-independent"}
      else {df[row_number,"sig"] = "non-sig"}
    }
    df$sig = factor(df$sig, levels=c("non-sig","eIF4A2-independent","eIF4A2-dependent"))
    
    
    assign(paste0(names(master_list[int/2]),"_complete"), df)
    df2=subset(df,df$pvalue.dmso<0.1&df[,paste0("pvalue.",treatment)]<0.1)
    assign(paste0(names(master_list[int/2]),"_complete_sig"), df2)

  }
}

#adds half life
master_light_complete_sig$halflife_nt = log(2)/master_light_complete_sig$rate.dmso
master_light_complete_sig$halflife_a2ko = log(2)/master_light_complete_sig$rate.a2ko
master_light_complete_sig$log2FC_halflife= log(master_light_complete_sig$halflife_a2ko/master_light_complete_sig$halflife_nt,2)
master_heavy_complete_sig$halflife_nt = log(2)/master_heavy_complete_sig$rate.dmso
master_heavy_complete_sig$halflife_a2ko = log(2)/master_heavy_complete_sig$rate.a2ko
master_heavy_complete_sig$log2FC_halflife= log(master_heavy_complete_sig$halflife_a2ko/master_heavy_complete_sig$halflife_nt,2)

