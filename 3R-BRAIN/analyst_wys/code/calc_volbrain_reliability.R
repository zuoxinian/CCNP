################################################################################################################################
# The between session, between scanner and total reliability analysis
# Between session reliability is defined as the ratio of between participants variance to variance caused by different sessions
# Between scanner reliability is defined as the ratio of between participants variance to variance caused by different scanners
# Total reliability is defined as the ratio of between participants variance to variance we observed
# sigma_subid square is the variance between participants
# sigma_scanner square is the variance between scanners
# sigma_session square is the variance between sessions
# sigma_scanner_session square is the variance caused by different sessions and different scanner
# sigma_subid_session square is the variance by inter individual differences differ by sessions
# sigma_subid_scanner square is the variance by inter individual differences differ by scanners
# Asymmetry Index(AI) is defined as ((LeftHemi - RightHemi)/(LeftHemi + RightHemi))*2
###############################################################################################################################
rm(list=ls())
#load packages
library("lmerTest")           
library("boot")

#set path
workdir <- '/Home/3R-BRAIN/'      #Should change it to your path accordingly
setwd(workdir)

#load data file
fname <- paste0(workdir,'/data/','volbrain_reliability_3R-BRAIN.csv')
volbrain_df <- read.csv(fname,header = TRUE)

#define bootstrap function to compute 3 different iccs
##########################################functions###########################################################################
calc.icc_session <- function(fm){
  results <- summary(fm)
  varcor <- as.data.frame(results$varcor)                   
  #variance of each factor
  sigma_subid <- varcor[3,4]
  sigma_scanner <- varcor[5,4]
  sigma_session <- varcor[6,4]
  sigma_scanner_session <-varcor[4,4]
  sigma_subid_scanner <- varcor[1,4]
  sigma_subid_session <- varcor[2,4]
  sigma_error <- varcor[7,4]
  sigma_all <- sum(varcor$vcov)
  vary <- var(tmpdf$y)
  sigma_subid/(sigma_subid+sigma_session+sigma_subid_session+sigma_error)
}
calc.icc_scanner <- function(fm){
  results <- summary(fm)
  varcor <- as.data.frame(results$varcor)                   
  #variance of each factor
  sigma_subid <- varcor[3,4]
  sigma_scanner <- varcor[5,4]
  sigma_session <- varcor[6,4]
  sigma_scanner_session <-varcor[4,4]
  sigma_subid_scanner <- varcor[1,4]
  sigma_subid_session <- varcor[2,4]
  sigma_error <- varcor[7,4]
  sigma_all <- sum(varcor$vcov)
  vary <- var(tmpdf$y)
  sigma_subid/(sigma_subid+sigma_scanner+sigma_subid_scanner+sigma_error)
}

calc.icc <- function(fm){
  results <- summary(fm)
  varcor <- as.data.frame(results$varcor) 
  sigma_subid <- varcor[3,4]
  sigma_scanner <- varcor[5,4]
  sigma_session <- varcor[6,4]
  sigma_scanner_session <-varcor[4,4]
  sigma_subid_scanner <- varcor[1,4]
  sigma_subid_session <- varcor[2,4]
  sigma_error <- varcor[7,4]
  sigma_all <- sum(varcor$vcov)
  vary <- var(tmpdf$y)
  icc <- sigma_subid/sigma_all
}


############################################ calculate icc that gained by volbrain ############################################

#centralize data
center_cols_names <- c("interval_baseline","age","height","weight","education")
chosen_cols <- volbrain_df[,center_cols_names]
centered_values <- scale(chosen_cols,center = TRUE)
centered_cols_names <- c("cinterval_baseline","cage","cheight","cweight","ceducation")
volbrain_df[,centered_cols_names] <- centered_values
volbrain_df$sex <- as.factor(volbrain_df$sex)

#create results dataframe
metric_names <- names(volbrain_df)[18:77]
nmetrics <- length(metric_names)
rownames <- c('icc_session','icc_scanner','icc','vary','var_subid','var_scanner','var_session','var_scanner_session','var_subid_scanner','var_subid_session','varw')
nrows <- length(rownames)
iccs <- data.frame(matrix(0,nrows,nmetrics))
names(iccs) <- metric_names
row.names(iccs) <- rownames
p_values <- data.frame(matrix(0,2,nmetrics))
names(p_values) <- metric_names
row.names(p_values) <- c('Intercept','pvalue')
ci_2.5 <- data.frame(matrix(0,3,nmetrics))
names(ci_2.5) <- metric_names
row.names(ci_2.5) <- c('icc_session','icc_scanner','icc')
ci_97.5 <- data.frame(matrix(0,3,nmetrics))
names(ci_97.5) <- metric_names
row.names(ci_97.5) <- c('icc_session','icc_scanner','icc')

#fit the 3-level lme model
for (i in 1:nmetrics){
  tmpdf <- volbrain_df
  print(metric_names[i])
  metric <- metric_names[i]
  tmpdf$y <- tmpdf[[metric]]
  index_outliers <- !tmpdf$y %in% boxplot.stats(tmpdf$y)$out  #exclude outliers
  tmpdf <- tmpdf[index_outliers,]
  try({
  #use lmer function fit linear mixed effect model
  fm <- lmer(y ~cage + sex + cheight + cweight + ceducation + cinterval_baseline + (1|subid)+(1|session)+(1|scanner)+(1|subid:session)+(1|subid:scanner)+(1|session:scanner),data=tmpdf)
  results <- summary(fm)
  varcor <- as.data.frame(results$varcor)
  p_value <- results$coefficients[1,c(1,5)]
  #variance of each factor
  sigma_subid <- varcor[3,4]
  sigma_scanner <- varcor[5,4]
  sigma_session <- varcor[6,4]
  sigma_scanner_session <-varcor[4,4]
  sigma_subid_scanner <- varcor[1,4]
  sigma_subid_session <- varcor[2,4]
  sigma_error <- varcor[7,4]
  sigma_all <- sum(varcor$vcov)
  vary <- var(tmpdf$y)
  
  var_subid <- sigma_subid/vary  
  var_scanner <- sigma_scanner/vary
  var_session <- sigma_session/vary
  var_scanner_session <- sigma_scanner_session/vary
  var_subid_scanner <- sigma_subid_scanner/vary
  var_subid_session <- sigma_subid_session/vary
  varw <- sigma_error/vary
  #calculate icc
  icc <- sigma_subid/sigma_all
  icc_session <- sigma_subid/(sigma_subid+sigma_session+sigma_subid_session+sigma_error)
  icc_scanner <- sigma_subid/(sigma_subid+sigma_scanner+sigma_subid_scanner+sigma_error)
  
  
  #assign values
  iccs[,metric] <- c(icc_session,icc_scanner,icc,vary,var_subid,var_scanner,var_session,
                     var_scanner_session,var_subid_scanner,var_subid_session,varw)
  
  p_values[,metric] <- p_value
  #boot strap
  boot.icc_sessions = bootMer(fm,calc.icc_session,1000)
  ci <- quantile(boot.icc_sessions$t,c(0.025,0.975))
  ci_2.5['icc_session',metric] <- ci[["2.5%"]]
  ci_97.5['icc_session',metric] <- ci[["97.5%"]]
  
  boot.icc_scanner = bootMer(fm,calc.icc_scanner,1000)
  ci <- quantile(boot.icc_scanner$t,c(0.025,0.975))
  ci_2.5['icc_scanner',metric] <- ci[["2.5%"]]
  ci_97.5['icc_scanner',metric] <- ci[["97.5%"]]
  
  boot.icc = bootMer(fm,calc.icc,1000)
  ci <- quantile(boot.icc$t,c(0.025,0.975))
  ci_2.5['icc',metric] <- ci[["2.5%"]]
  ci_97.5['icc',metric] <- ci[["97.5%"]]
  },silent=FALSE)
}

#check output dir
output_dir <- paste0(workdir,'/results/','reliability/')
ifelse(dir.exists(output_dir),print("Saving results into reliability folder"),dir.create(output_dir))

#save results
fout <- paste0(output_dir,'volbrain_icc.csv')
write.csv(iccs,fout)

fout <- paste0(output_dir,'volbrain_pvalue.csv')
write.csv(p_values,fout)

fout <- paste0(output_dir,'volbrain_icc_ci25.csv')
write.csv(ci_2.5,fout)

fout <- paste0(output_dir,'volbrain_icc_ci975.csv')
write.csv(ci_97.5,fout)


#######################################cal AI icc############################################################################
#calculat AI
cal_AI <- function(L,R){((L-R)/(L+R))*2} #define lr asymmetry function
## right hemi
df_names <- colnames(volbrain_df)
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
#setup dataframe for AI analysis
rownames <- c('icc_session','icc_scanner','icc','vary','var_subid','var_scanner','var_session','var_scanner_session','var_subid_scanner','var_subid_session','varw')
nrows <- length(rownames)
iccs <- data.frame(matrix(0,nrows,nmeasures))

names(iccs) <- struct_names
row.names(iccs) <- rownames
p_values <- data.frame(matrix(0,2,nmeasures))
names(p_values) <- struct_names
row.names(p_values) <- c('Intercept','pvalue')
ci_2.5 <- data.frame(matrix(0,3,nmeasures))
names(ci_2.5) <- struct_names
row.names(ci_2.5) <- c('icc_session','icc_scanner','icc')
ci_97.5 <- data.frame(matrix(0,3,nmeasures))
names(ci_97.5) <- struct_names
row.names(ci_97.5) <- c('icc_session','icc_scanner','icc')

for (i in 1:nmeasures){
  tmpdf <- volbrain_df
  L <- volbrain_df[[L_df_names[i]]]
  R <- volbrain_df[[R_df_names[i]]]
  y <- cal_AI(L,R)
  print(struct_names[i])
  metric <- struct_names[i]
  index_outliers <- !y %in% boxplot.stats(y)$out  #exclude outliers
  tmpdf$y <- y
  tmpdf <- tmpdf[index_outliers,]
  try({
    fm <- lmer(y~cage + sex + cheight + cweight + ceducation + cinterval_baseline + (1|subid)+(1|session)+(1|scanner)+(1|subid:session)+(1|subid:scanner)+(1|session:scanner),data=tmpdf)
    results <- summary(fm)
    varcor <- as.data.frame(results$varcor)
    p_value <- results$coefficients[1,c(1,5)]
    #variance of each factor
    sigma_subid <- varcor[3,4]
    sigma_scanner <- varcor[5,4]
    sigma_session <- varcor[6,4]
    sigma_scanner_session <-varcor[4,4]
    sigma_subid_scanner <- varcor[1,4]
    sigma_subid_session <- varcor[2,4]
    sigma_error <- varcor[7,4]
    sigma_all <- sum(varcor$vcov)
    vary <- var(tmpdf$y)
    
    var_subid <- sigma_subid/vary
    var_scanner <- sigma_scanner/vary
    var_session <- sigma_session/vary
    var_scanner_session <- sigma_scanner_session/vary
    var_subid_scanner <- sigma_subid_scanner/vary
    var_subid_session <- sigma_subid_session/vary
    varw <- sigma_error/vary
    #calculate icc
    icc <- sigma_subid/sigma_all
    icc_session <- sigma_subid/(sigma_subid+sigma_session+sigma_subid_session+sigma_error)
    icc_scanner <- sigma_subid/(sigma_subid+sigma_scanner+sigma_subid_scanner+sigma_error)
    
    
    #assign values
    iccs[,metric] <- c(icc_session,icc_scanner,icc,vary,var_subid,var_scanner,var_session,
                       var_scanner_session,var_subid_scanner,var_subid_session,varw)
    
    p_values[,metric] <- p_value
    #boot strap
    boot.icc_sessions = bootMer(fm,calc.icc_session,1000)
    ci <- quantile(boot.icc_sessions$t,c(0.025,0.975))
    ci_2.5['icc_session',metric] <- ci[["2.5%"]]
    ci_97.5['icc_session',metric] <- ci[["97.5%"]]
    
    boot.icc_scanner = bootMer(fm,calc.icc_scanner,1000)
    ci <- quantile(boot.icc_scanner$t,c(0.025,0.975))
    ci_2.5['icc_scanner',metric] <- ci[["2.5%"]]
    ci_97.5['icc_scanner',metric] <- ci[["97.5%"]]
    
    boot.icc = bootMer(fm,calc.icc,1000)
    ci <- quantile(boot.icc$t,c(0.025,0.975))
    ci_2.5['icc',metric] <- ci[["2.5%"]]
    ci_97.5['icc',metric] <- ci[["97.5%"]]
    
  },silent=FALSE)
}

#save results
fout <- paste0(output_dir,'volbrain_icc_subcortical_AI.csv')
write.csv(iccs,fout)

fout <- paste0(output_dir,'volbrain_pvalue_subcortical_AI.csv')
write.csv(p_values,fout)

fout <- paste0(output_dir,'volbrain_icc_ci25_subcortical_AI.csv')
write.csv(ci_2.5,fout)

fout <- paste0(output_dir,'volbrain_icc_ci975_subcortical_AI.csv')
write.csv(ci_97.5,fout)

print("FINISH ALL")
