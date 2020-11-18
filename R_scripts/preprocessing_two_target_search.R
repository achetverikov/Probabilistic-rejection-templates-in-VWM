library(apastats)
load.libs(c('data.table','stringr','plyr','dplyr','Hmisc'))
angle_dist<-function(a,b){
  c = a - b
  (c+90)%%180 - 90
}

data_e2 <- fread('../bimodal_two_target_search.csv')

# exclude subjects with very low accuracy or very high RT
data_e2 <- data_e2[subjectId %nin% data_e2[blockType=='probe',mymean(rt_2),by=.(subjectId)][V1>1400, subjectId]]

# reshape data to long format
datam<-reshape(data_e2, varying=c('answer_1','correct_1','rt_1','answer_2','correct_2','rt_2'), dir='long', sep='_')

# determine which target was responded to
datam[answer==corrResp_1, targetOri:=targetOri_1]
datam[answer==corrResp_2, targetOri:=targetOri_2]

datam[,targetDist:=angle_dist(distrMean, targetOri)]
datam[,ctpd:=angle_dist(prevDistrMean, targetOri)]
datam[,actpd:=abs(round(ctpd))]

# determine target type for test trials
datam[actpd%in%c(25,30), stimType:='peak']
datam[actpd%in%c(0), stimType:='inbetween']
datam[actpd%in%c(50,60), stimType:='out']

datam[,lrt:=log(rt)]

# reshape long to wide
data_r<-dcast(datam[blockType=='probe'],subjectId+totCorrect+session+block+blockType+trial+targets~time, value.var= c('answer','correct','rt','lrt','ctpd','actpd','targetOri','stimType','targetDist'))

# determine target type for trials with only one correct response 
data_r[correct_2==0&correct_1==1, stimType_2:=setdiff(stimType_1,str_split(targets,'_')),by=.(subjectId, block, blockType, trial)]
data_r[correct_1==0&correct_2==1, stimType_1:=setdiff(stimType_2,str_split(targets,'_')),by=.(subjectId, block, blockType, trial)]

# reshape back to long
datam2<-reshape(data_r, varying=c('actpd_1','actpd_2','answer_1','correct_1','rt_1','answer_2','correct_2','rt_2','stimType_1','stimType_2','lrt_1','lrt_2'), dir='long', sep='_')
datam2[,stimType:=factor(stimType, levels=c('peak','inbetween','out'), labels=c('Peak','In-between','Outside'))]

datam2[,actpdf:=factor(actpd)]
datam2[,timef:=factor(time)]

# is it the first or the second response?
datam2[,p_first:=1-(time-1)]

data_r[,atarget_diff:=abs(angle_dist(targetOri_1, targetOri_2))]
data_r[,targets2:=relevel(factor(targets), ref='Peak + Peak')]

datam[,stimTypef:=factor(stimType, levels=c('out','inbetween','peak'),labels=c('Outside','Between','Peak'))]
datam[, expName:='Main Experiment']
datam_e2 <- datam
data_r_e2 <- data_r
rm(datam)
rm(data_r)
# preparations for simulations and plots

data_r_e2[totCorrect==2&trial==0,grandmean:=mymean(lrt_2)]
data_r_e2[totCorrect==2,lrt_2c:=scale(lrt_2, scale=F), by=subjectId]
data_r_e2[,targets_fp:=factor(targets, levels=c('Peak + Outside', 'Peak + In-between','Outside + In-between','Peak + Peak'), labels = c('P+O','P+B','O+B','P+P'))]
data_r_e2[,short_exp:='e2']

datam2<-reshape(data_r_e2, varying=c('actpd_1','actpd_2','answer_1','correct_1','rt_1','answer_2','correct_2','rt_2','stimType_1','stimType_2','lrt_1','lrt_2'), dir='long', sep='_')
datam2[,stimType:=factor(stimType, levels=c('peak','inbetween','out'), labels=c('Peak','Between','Outside'))]

datam2[,actpdf:=factor(actpd)]
datam2[,timef:=factor(time)]

datam2[,p_first:=2-(time)]
datam2[,ctpd:=NULL]
datam2[,ctpd:=ifelse(stimType=='Outside', ifelse(abs(ctpd_1)==50,ctpd_1, ctpd_2),NA)]
datam2[,ctpd_other:=ifelse(stimType=='Outside', ifelse(abs(ctpd_1)!=50,ctpd_1, ctpd_2),NA)]
datam2[is.na(ctpd), ctpd:=ifelse(stimType=='Between', ifelse(ctpd_1==0,ctpd_1, ctpd_2),NA)]
datam2[is.na(ctpd_other), ctpd_other:=ifelse(stimType=='Between', ifelse(ctpd_1!=0,ctpd_1, ctpd_2),NA)]

datam2[,ctpd_comb:=paste(pmin(ctpd_1, ctpd_2), pmax(ctpd_1, ctpd_2), sep='.')]
datam2_fp<-datam2[totCorrect==2&trial==0&targets!='Peak + Peak'&(stimType=='Outside'|targets=='Peak + In-between')&stimType!='Peak']

datam2_fp[,targets_stimType:=factor(interaction(stimType,targets))]
datam2_fp[,targets_stimType:=factor(targets_stimType, levels=c('Outside.Peak + Outside','Between.Peak + In-between','Outside.Outside + In-between'), labels=c('O|P+O','B|P+B','O|O+B'))]

