#package: stringr

library(stringr)
library(enviPat)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("LICAR_functions.R")
# Warning:need to add 2 class in one line
# DLPC is changed manually

### raw data file input
#aimpath <- readline("Please enter the aim file :")
content <- c("TL_pos_median_concentration_by_biorepl.csv","TL_neg_median_concentration_by_biorepl.csv")
for (i in 1:length(content)){
  aimpath <- content[i]

# determine the mode
if (grepl("pos",aimpath,ignore.case=TRUE)==TRUE) {
  mode <- "POS" 
} else if (grepl("neg",aimpath,ignore.case=TRUE)==TRUE) {
  mode <- "NEG"
}
### read file and standardise lipidclass
raw <- read.csv(aimpath,header=TRUE)


num_index <- c(str_locate(raw$ion,pattern = '[0-9]')[,'start']) #find the number position

raw$ion2 <- paste(substr(raw$ion,start=1,stop = num_index-1),substr(raw$ion,start=num_index,stop=length(raw$ion)-num_index+1)) #insert space between lipid class and carbon number

#standardize the lipid class name
std_class <- function(x){
for (i in 1:length(x)){
  if (grepl('lyso',x[i],ignore.case = TRUE)==TRUE) {
    x[i]<- gsub('lyso','L',x[i],ignore.case = TRUE)
  } else if (grepl('(o-)',x[i],ignore.case = TRUE)==TRUE){
    x[i] <- gsub('\\(O-\\)','O',x[i],ignore.case = TRUE)
    } else if (grepl('_C',x[i])==TRUE){
      x[i] <- gsub('_C','_',x[i])
} else next}
return(x)
  }
raw$ion2 <- std_class(unlist(raw$ion2))
raw$LipidClass <- lapply(strsplit(raw$ion2," "),function(x){x[1]})


### determine the lipid group
pos_HG <- c('PC','PE','LPC','LPE','PCO','PEO')
neg_HG <- c('PI','PS','LPI','LPS','PIO')
neg_FA <- c('CL')
raw$LipidClass <- factor(unlist(raw$LipidClass))
l_num <- levels(raw$LipidClass)

#format raw to pre_df
rownames(raw) <- raw[,'ion2']
for (i in 1:length(l_num)) {
df <- raw[raw$LipidClass==l_num[i],]
drops <- c("ion","standard_ion","LipidClass","LipidClass2","ion2")
df <- df[ ,!(names(df) %in% drops)]
colnames(df)[1] <- "Precursor"
colnames(df)[2] <- "Product"
rawdata <- df[order(df$Precursor, df$Product), ]
if (i==1) {
  tot_res <- rawdata[0,]
  tot_ratio <- rawdata[0,]
  }
get_LG <- function(lipidclass, mode){
  if (mode == "POS") {
    if (lipidclass %in% pos_HG){
      lipidgroup <- "Head Group"
    }
  }
  if (mode=="NEG") {
    if (lipidclass %in% neg_HG){
      lipidgroup <- "Head Group"
    } else if (lipidclass %in% neg_FA){
      lipidgroup <- "FA"
    }
  }
  return(lipidgroup)
}
lipidgroup <- get_LG(l_num[i],mode)
if (mode=='NEG') {
if (lipidgroup=='FA') {
  l_num[i] <- paste(l_num[i],'N',lipidgroup,sep='')
} else if (lipidgroup=='Head Group'){
  if (l_num[i]=='PI'|l_num[i]=='PS'){
    l_num[i] <- paste(l_num[i],'NHG',sep='')
  }
}
}
# execute isotope correction


# check lipid class by using productIon-param3
class_name <- sub( " .*", "", rownames(rawdata))
if((length(table(class_name)) > 1) && (lipidgroup %in% c("Head Group", "FA", "LCB", "Neutral"))) {
  stop(paste("Lipid class is not unique, please check the data! Continue if it belongs to group RPLC."))}

# isotope correction--lipidclass-param2
correctedResults <- isoCorrect(rawdata, lipidClass=l_num[i], lipidGroup=lipidgroup)
tot_res <- rbind(tot_res,correctedResults)


# calculate the ratio
ratio <- correctedResults/rawdata
ratio[, c("Precursor", "Product")] <- rawdata[,c("Precursor", "Product")]
tot_ratio <- rbind(tot_ratio,ratio)
}

# add names
tot <- list("Correction Res"=tot_res,"Correction Ratio"=tot_ratio)
tot <- lapply(tot, function(x){
  names <- rownames(x)
  rownames(x) <- NULL
  x <- cbind(names,x)
})

# add class to tot_res
tot <- lapply(tot,function(x){
res_class <- unlist(lapply(strsplit(x$names," "),function(x){x[1]}))
x <- cbind(x,res_class)
})


# write output
outpath <- paste(str_sub(aimpath,start=1,end=-5),"_after_corr.xlsx",sep='')


writexl::write_xlsx(tot,
outpath)
}