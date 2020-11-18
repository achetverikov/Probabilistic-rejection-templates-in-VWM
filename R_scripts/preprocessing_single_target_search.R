library(apastats)
load.libs(c('data.table','stringr','plyr','dplyr','Hmisc'))
angle_dist<-function(a,b){
  c = a - b
  (c+90)%%180 - 90
}

data_e1 <- fread(file = '../bimodal_single_target_search.csv')

# exclude subjects with very low accuracy or very high RT
data_e1 <- data_e1[subjectId %nin% data_e1[,.(acc=mymean(correct), mean_rt=mymean(rt)), by=.(subjectId)][acc < 0.75 | mean_rt > 1400, subjectId]
]

data_e1[,lrt:=log(rt)]
ds_data_e1<-data_e1[blockType=='probe']
ds_data_e1[,actpd:=abs(t_dist_to_prev_d)]
ds_data_e1[,peak_type:=case_when(actpd>50 ~ 'outside', actpd>20&actpd<40 ~ 'peak', actpd<15 ~ 'in-between')]

ds_data_e1[,stimTypef:=factor(peak_type, levels=c('outside','in-between','peak'),labels=c('Outside','Between','Peak'))]
ds_data_e1[, expName:='Suppl. Experiment']

