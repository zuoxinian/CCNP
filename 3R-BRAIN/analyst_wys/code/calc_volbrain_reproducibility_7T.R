rm(list=ls())
workdir <- '~/R3BRAIN_for_OHBM/'  #change it to your data path accordingly
setwd(workdir)
library("dplyr")
library("effsize")
library("pwr")
#load 3T reliability results
data_csv <- paste0(workdir,'data/','volbrain_replicability_7T.csv')
data_df <- read.csv(data_csv)
volbrain_ipcas <- data_df %>% filter(scanner=='IPCAS_MR750')
#volbrain_ipcas <- data_df[which(data_df$scanner=='IPCAS_MR750'),]


# orginaze dataframe
## right hemi
df_names <- colnames(volbrain_ipcas)
R_volume_names <- regexpr('.R.[A-Z][A-Z].cm3|.Right.cm3',df_names) #find metric loc have Right hemi volumes name 
R_volume_names_index <- ifelse(R_volume_names==-1,FALSE,TRUE) #
R_df_names <- df_names[R_volume_names_index]
R_volume_str_end <- R_volume_names[R_volume_names >0] - 1
struct_names <- substr(R_df_names,1,R_volume_str_end) #get structure names
struct_names[c(4,5,7,8)] <- c('Cerebelum.GM','Cerebelum.WM','Cerebrum.GM','Cerebrum.WM') #manually edit
## left hemi 
L_volume_names <- regexpr('.L.[A-Z][A-Z].cm3|.Left.cm3',df_names) #find metric loc have left hemi volumes name 
L_volume_names_index <- ifelse(L_volume_names==-1,FALSE,TRUE) #
L_df_names <- df_names[L_volume_names_index]

nmeasures <- length(L_df_names)

#centralize data
center_cols_names <- c("age","height","weight","education")
chosen_cols <- volbrain_ipcas[,center_cols_names]
centered_values <- scale(chosen_cols,center = TRUE)
centered_cols_names <- c("cage","cheight","cweight","ceducation")
volbrain_ipcas[,centered_cols_names] <- centered_values

#data_frame_values
rownames <- c('t','pvalue','mean_diff','lower_bound','upper_bound','effsize','eff_low','eff_up','pwr')
df_lr <- data.frame(matrix(0,nrow = length(rownames),ncol = nmeasures ),row.names = rownames)
colnames(df_lr) <- struct_names
for (i in 1:nmeasures){
  print(struct_names[[i]])
  L_metric <- volbrain_ipcas[,L_df_names[i]]
  R_metric <- volbrain_ipcas[,R_df_names[i]]
  #paired t-test
  result <- t.test(L_metric,R_metric,paired = T)
  effsize <- cohen.d(L_metric,R_metric,paired=T)
  #cohen's d and conf.int
  cohend <- effsize$estimate
  cohend_conf.int <- effsize$conf.int
  stpwoer <- pwr.t.test(n=length(L_metric),d=cohend,sig.level=0.05,power = NULL,type = "paired",alternative = "two.sided")
  df_lr[,i] <- c(result[['statistic']],result[['p.value']],result[['estimate']],result[["conf.int"]][1],result[["conf.int"]][2],cohend,cohend_conf.int[1],cohend_conf.int[2],stpwoer$power)
}
fout <- paste0(workdir,'results/','reproducibility/','volbrain_ipcas_for-7T_LR_diff_by_wys.csv')
write.csv(df_lr,fout)


####################### gender differences########################################
#calculat AI
cal_AI <- function(L,R){((L-R)/(L+R))*2} #define lr asymmetry function

volbrain_ipcas$gender <- as.factor(volbrain_ipcas$gender)
ICV <- volbrain_ipcas$Tissue.IC.cm3
volbrain_ipcas$cICV <- scale(ICV,center=T)

#cal gender and age diff
gender_distribute_diff <- chisq.test(summary(volbrain_ipcas$gender))
volbrain_ipcas_male <- volbrain_ipcas %>% filter(gender=='1')
volbrain_ipcas_female <- volbrain_ipcas %>% filter(gender=='0')
age_diff <- t.test(volbrain_ipcas_male$age,volbrain_ipcas_female$age)
education_diff <- t.test(volbrain_ipcas_male$education,volbrain_ipcas_female$education)
#save gender age education p
matched_test <- data.frame(gender_chisq_p=gender_distribute_diff$p.value,age_diff_p = age_diff$p.value,education_diff_p=education_diff$p.value)
fout <-  paste0(workdir,'results/','reproducibility/','gender_match_test-7T.csv')
write.csv(matched_test,fout)

#loop each metric calculate the residual p value after regressed age and icv
idx_male <- which(volbrain_ipcas['gender']=='1')
idx_female <- which(volbrain_ipcas['gender']=='0')
rownames <- c('t','pvalue','mean_diff_male','mean_diff_female','lower_bound','upper_bound','effsize','eff_low','eff_up','pwr')
df_gender_diff <- data.frame(matrix(0,nrow = length(rownames),ncol = nmeasures ),row.names = rownames)
colnames(df_gender_diff) <- struct_names

for (k in 1:nmeasures){
  print(struct_names[k])
  metric_AI <- cal_AI(volbrain_ipcas[,L_df_names[k]],volbrain_ipcas[,R_df_names[k]])
  tmp_df <- data.frame(metric_AI=metric_AI,cICV=volbrain_ipcas$cICV)
  residuals <- lm(metric_AI ~ cICV,data = tmp_df) %>% residuals  #regress age & icv
  result <- t.test(residuals[idx_male],residuals[idx_female])
  effsize <- cohen.d(residuals[idx_male],residuals[idx_female])
  #cohen's d and conf.int
  cohend <- effsize$estimate
  cohend_conf.int <- effsize$conf.int
  nsubj <- result$n
  stpwoer <- pwr.t.test(n=length(residuals),d=cohend,sig.level=0.05,power = NULL,type = "two.sample",alternative = "two.sided")
  df_gender_diff[,k] <- c(result[['statistic']],result[['p.value']],result[['estimate']][1],result[['estimate']][2],result[["conf.int"]][1],result[["conf.int"]][2],cohend,cohend_conf.int[1],cohend_conf.int[2],stpwoer$power)
}
fout <- paste0(workdir,'results/','reproducibility/','volbrain_ipcas_for-7T_gender-diff_by_wys.csv')
write.csv(df_gender_diff,fout)
