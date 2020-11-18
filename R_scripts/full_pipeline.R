#' ---
#' title: "Analyses for the manuscript \"Probabilistic working memory templates guide visual attention\" by Chetverikov, Campana, & Kristjansson" 
#' author: "Andrey Chetverikov"
#' date: '`r format(Sys.time(), "%d %B, %Y")`'
#' output: html_document
#' ---
#' 

#' Set the working directory to be the same as the script directory. 
#' The scripts assume that the data is in the folder one level above

set.seed(657645644) # for bootstrapping reproducibility

# load libraries; apastats can be install from here: https://github.com/achetverikov/APAstats (it's mostly used for outputting statistics)

library(apastats)
load.libs(c('circular','doSNOW','foreach','doParallel','mgcv','ggplot2',
            'data.table', 'bbmle','stringr','plyr','dplyr','scales','Hmisc', 
            'ggthemes', 'cowplot','lmSupport','lme4', 'lmerTest','ez'))


# temporarily suppress warnings
oldw <- getOption("warn")
options(warn = -1)

source('preprocessing_single_target_search.R')
source('preprocessing_two_target_search.R')

# combine single-target results for two experiments
single_target_data<-rbind(
  datam_e2[totCorrect==2&trial==0&blockType=='probe', 
           .(lrt, stimTypef, ctpd=ctpd, subjectId, expName)],
  ds_data_e1[correct==1&trial==0, 
             .(lrt, stimTypef, ctpd=t_dist_to_prev_d, subjectId, expName)])
single_target_data[,short_exp:=factor(expName, levels = c('Main Experiment','Suppl. Experiment'), labels=c('e2','e1'))]
single_target_data[,grandmean:=mymean(lrt),by=.(short_exp)]
single_target_data[,lrtc:=scale(lrt, scale=F), by=.(subjectId, expName)]

single_target_data[short_exp=='e2',mean_ctpd:=ctpd]
single_target_data[stimTypef=='Between', mean_ctpd:=0]
single_target_data[stimTypef=='Outside'&short_exp=='e1', mean_ctpd:=90]
single_target_data[stimTypef=='Peak'&short_exp=='e1', mean_ctpd:=sign(ctpd)*30]

# simulations take time, so we load the pre-computed results instead
# NB: this script uses multicore parallel processing set with SNOW that might 
# not work everywhere; in this is the case, simply replace %dopar% with %do% in
# 'foreach' call
#
# source('simulations.R') 
load('simulations_results.RData')

template_sds[,1:3] # estimated template SDs 
bic_tbl # models BICs
bic_tbl[,.(peaks, 
           bic_rt_st=bic_rt_st-bic_rt_st[which(peaks=='two')], 
           bic_pf=bic_pf-bic_pf[which(peaks=='two')], 
           bic_rt_2t=bic_rt_2t-bic_rt_2t[which(peaks=='two')]), 
        by=exp] # differences in BICs

#' # Statistics for Supplementary Experiment

#' ## Overall performance

#' ### Block type (prime vs test)
data_e1[correct==1,describe.mean.and.t(lrt, blockType, 3, paired=T, eff.size=T, aggregate_by = subjectId, transform_means=exp)]
data_e1[,describe.mean.and.t(correct, blockType, 3, paired=T, eff.size=T, aggregate_by = subjectId)]

#' ### Learning within prime streaks
data_prime_e1<-data_e1[blockType=='prime']
data_prime_e1[,trial:=factor(trial)]
contrasts(data_prime_e1$trial) <- varContrasts(data_prime_e1$trial, Type='HELMERT', RefLevel = 'Highest')

fit_rep_rt<-lmer(lrt~trial+(1+trial|subjectId),data_prime_e1[correct==1])
fit_rep_acc<-lmer(correct~trial+(1+trial|subjectId),data_prime_e1) #note: binomial regression gives the same results, but throws convergence warnings, so we use Gaussian model 
describe.glm(fit_rep_rt,'trial0_v_Later', dtype = 3, test.df=T)
describe.glm(fit_rep_acc,'trial0_v_Later', dtype = 3, test.df=T)

#' ## Test trials performance
ds_data_1<-ds_data_e1[correct==1&trial==0&prev_correct==1, .(lrt=mymean(lrt)), keyby=.(subjectId, peak_type)]
ds_data_1[peak_type!='outside',describe.mean.and.t(lrt, factor(peak_type), 0, eff.size = T, paired=T)]
ds_data_1[peak_type!='peak',describe.mean.and.t(lrt, factor(peak_type), 0, eff.size = T, paired=T)]

#' # Statistics for Main Experiment

#' ## Overal performance

#' ### Block type (prime vs test)

acc_table<-dcast(data_e2[,data.table(prop.table(table(totCorrect))), by=.(subjectId, blockType)][,totCorrect:=paste0('t',totCorrect)], subjectId+blockType~totCorrect)
acc_table[is.na(t0),t0:=0]
perf_table<-merge(acc_table, data_e2[totCorrect==2, .(rt_2=exp(mymean(log(rt_2))), rt_1=exp(mymean(log(rt_1)))), by=.(subjectId, blockType)])
perf_table[,blockType:=factor(blockType, labels=c('Learning','Test'))]

perf_table[blockType=='Learning',describe.mean.conf(t2)]
perf_table[blockType=='Learning',describe.mean.conf(t1)]
perf_table[blockType=='Learning',describe.mean.conf(t0)]

perf_table[blockType=='Test',describe.mean.conf(t2)]

data_e2[totCorrect==2,mymean(((rt_2-rt_1))),by=.(subjectId,blockType)][,describe.mean.and.t(V1,blockType, 3, paired = T, digits=0)]
data_e2[totCorrect==2,mymean(log(rt_1)),keyby=.(subjectId,blockType)][,describe.mean.and.t(V1,blockType, 3, paired = T, transform=T, digits=0)]

#' ### Learning effects
data_prime_2nd<-datam_e2[time==2&blockType=='prime']
data_prime_2nd[,trial:=factor(trial)]
contrasts(data_prime_2nd$trial)<-varContrasts(data_prime_2nd$trial, Type='HELMERT', RefLevel = 'Highest')

fit_rep_rt<-lmer(lrt~trial+(1+trial|subjectId),data_prime_2nd[totCorrect==2])
fit_rep_acc<-lmer(totCorrect~trial+(1+trial|subjectId),data_prime_2nd)

ins.lmer(fit_rep_rt,'trial0_v_Later')
ins.lmer(fit_rep_acc,'trial0_v_Later')

#' ## Test trials

#' ### Search time for a single target
st_ez<-ezANOVA(datam_e2[totCorrect==2&blockType=='probe'&trial==0], lrt, wid = subjectId, within = stimTypef)
describe.ezanova(st_ez, 'stimTypef')

#' ### Total search times for two targets
t2ez<-ezANOVA(data_r_e2[totCorrect==2], lrt_2, wid = subjectId, within = targets)
describe.ezanova(t2ez, 'targets')

data_ra<-data_r_e2[totCorrect==2&trial==0,.(lrt_2=mymean(lrt_2)), by=.(subjectId, targets)]

data_ra[targets%in%c('Peak + Peak','Outside + In-between'), describe.ttest(t.test(lrt_2~targets, .SD,paired=T))]
data_ra[targets%in%c('Peak + Outside','Peak + In-between'), describe.ttest(t.test(lrt_2~targets, .SD,paired=T))]
data_ra[targets%in%c('Peak + Peak','Peak + In-between'), describe.ttest(t.test(lrt_2~targets, .SD,paired=T))]

#' ### Probabilities of reporting first
datam2a<- datam2[totCorrect==2&trial==0&!is.na(stimType)&targets!='Peak + Peak',]

pfirst_glmer<-datam2a[,describe.glm(glmer(p_first~(1|subjectId),.SD, family = 'binomial')), keyby=.(targets, stimType)]

pfirst_glmer[stimType=='Peak'&targets=='Peak + In-between', str]
pfirst_glmer[stimType=='Peak'&targets=='Peak + Outside', str]
pfirst_glmer[stimType=='Outside'&targets=='Outside + In-between', str]

#' # Plot Figure 2 from the paper
source('final_plots.R')

#+ fig.width=8.225, fig.height=7, out.width='705px', out.height='600px', dpi=320
plot_grid(plot_grid(ggdraw() + draw_label("Results", fontface='bold', size = 11),
                    plot_grid(p1r, p2r,p3r,nrow=1, labels=c('A','B','C'), align = 'h', axis="tblr",label_size = 10), ncol=1, rel_heights = c(0.1,0.9)),
          plot_grid(ggdraw() + draw_label("Models' best-fitting predictions", fontface='bold', size = 11),
                    plot_grid(p2+theme(legend.position = 'none') + coord_cartesian(ylim=(ggplot_build(p1r)$layout$panel_ranges[[1]]$y.range)),
                              p1+theme(legend.position = 'none') + coord_cartesian(ylim=(ggplot_build(p2r)$layout$panel_ranges[[1]]$y.range)),
                              p3+theme(legend.position = c(1,1), legend.justification = c(1,1)) + coord_cartesian(ylim=(ggplot_build(p3r)$layout$panel_ranges[[1]]$y.range)), 
                              nrow=1, labels=c('D','E','F'), align = 'h',
                              axis="tblr", label_size = 10),
                    ncol=1, rel_heights = c(0.1,0.9)),
          nrow=2)
# turn warnings back on
options(warn = oldw)

