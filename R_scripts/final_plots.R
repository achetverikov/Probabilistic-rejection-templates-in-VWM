default_colors <-c('#3498db','#009f06','#AA2255','#FF7F00')
scale_colour_discrete <- function(...) scale_color_manual(values=default_colors, ...)
scale_fill_discrete <- function(...) scale_fill_manual(values=default_colors, ...)
scale_shape_discrete <- function(...) scale_shape_ac(...)
update_geom_defaults("point", list(shape=I(21), size=I(4),color='white'))

theme_set(plots_theme<-theme_light(base_family = 'sans')+theme(panel.grid.major = element_line(color='gray', linetype=2), legend.position = "right", strip.background = element_blank(),strip.text = element_text(color='black'), panel.border = element_blank(), plot.title = element_text(hjust=0, size = rel(1)), axis.title = element_text(size=rel(0.9), color=I('gray30')), legend.margin = margin(0)))

rt_scaling[,peaks_fp:=factor(peaks, levels=c('single', 'two', 'avgd'),labels=c('Single-peak','Probabilistic','Averaged'))]
p_first_scaling[,peaks_fp:=factor(peaks, levels=c('single', 'two', 'avgd'),labels=c('Single-peak','Probabilistic','Averaged'))]

rt_two_t_scaling[,peaks_fp:=factor(peaks, levels=c('single', 'two', 'avgd'),labels=c('Single-peak','Probabilistic','Averaged'))]

p1<-ggplot(rt_two_t_scaling[,.(value=mean(pred_rt+grandmean)),by=.( peaks_fp, targets_fp)], aes(group=peaks_fp, x=targets_fp, label=peaks_fp, color=peaks_fp, fill=peaks_fp, y=value))+geom_line()+geom_point(stat='summary',fun.y=mean, color=I('white'))+labs(x=NULL, y='RT for two targets (ms)', color=NULL, fill=NULL)+guides(fill=guide_legend(ncol=1,title = NULL),color=guide_legend(ncol=1,title = NULL))+theme(legend.position = c(1,1), legend.justification = c(1,1))+scale_y_exp()+scale_y_exp()+scale_x_discrete(breaks = c('O+B','P+O','P+B','P+P'))


p2<-ggplot(rt_scaling[exp=='e2',.(value=mean(y+grandmean)), by=.(peaks_fp,stimTypef)], aes(group=peaks_fp, x=stimTypef, fill=peaks_fp, color=peaks_fp, y=value))+geom_line()+geom_point(stat='summary',fun.y=mean,color=I('white'))+labs(x=NULL, y='RT for a single target (ms)', color=NULL, fill=NULL)+theme(legend.position = 'top')+guides(fill=guide_legend(ncol=1),color=guide_legend(ncol=1))+theme(legend.position = c(1,1), legend.justification = c(1,1))+scale_y_exp()


p3<-ggplot(p_first_scaling[,.(value=mean(as.numeric(prob))),by=.(peaks_fp, typefp)], aes(group=peaks_fp, x=typefp,fill=peaks_fp, color=peaks_fp, y=as.numeric(value)))+geom_line()+geom_point(color=I('white'))+labs(y='P(find first)', x=NULL, color=NULL, fill=NULL)+theme(legend.position = 'top')+guides(fill=guide_legend(ncol=1),color=guide_legend(ncol=1))+theme(legend.position = c(1,0), legend.justification = c(1,0))+geom_hline(yintercept = 0.5, linetype=2)

#predictions_plots<-ggarrange(p2,p1,p3, ncol=1)

single_target_data_by_subj <- single_target_data[,.(lrt = (mymean(lrtc+grandmean))), by=.(stimTypef,subjectId,expName)]
single_target_data[,lrtc_with_grandmean:=lrtc+grandmean]
p1r<-plot.pointrange(single_target_data, aes(x=stimTypef, y=lrtc_with_grandmean, group=1), wid = 'subjectId', 
                     within_subj = T, do_aggregate = T,exp_y = F, pointshape = NULL, 
                     print_aggregated_data = F, pointsize = NULL,pretty_y_axis = T, betweenvars = 'expName',
                     connecting_line = T, 
                     custom_geom = list(geom_violin(data=single_target_data_by_subj, aes(y=lrt, x=stimTypef), fill=I('black'),inherit.aes = F, alpha=I(0.1), color=I('white')),
                                        geom_jitter(data=single_target_data_by_subj, aes(y=lrt, x=stimTypef), color=I('black'), width=0.15, alpha=I(0.9), size=I(2),inherit.aes = F, fill=I('white'))))+
  labs(y='RT for a single target (ms)', x='Target Type')+
  scale_shape_ac()+
  guides(color='none', fill='none')+
  facet_grid(expName~., scale='free')+
  scale_y_exp()

data_r_e2_by_subj <- data_r_e2[totCorrect==2&trial==0,.(lrt = (mymean(lrt_2c+grandmean))), by=.(targets_fp,subjectId)]
data_r_e2[,lrt_2c_with_grandmean:=lrt_2c+grandmean]

p2r<-plot.pointrange(data_r_e2[totCorrect==2&trial==0], aes(x=targets_fp, y=lrt_2c_with_grandmean, group=1), 
                     wid = 'subjectId', within_subj = T, do_aggregate = T, connecting_line = T,
                     pointshape = NULL, pointsize = NULL, pretty_y_axis = T, 
                     custom_geom = list(
                       geom_violin(data=data_r_e2_by_subj, aes(y=lrt, x=targets_fp), fill=I('black'),inherit.aes = F, alpha=I(0.1), color=I('white')),
                       geom_jitter(data=data_r_e2_by_subj, aes(y=lrt, x=targets_fp), color=I('black'), width=0.15, alpha=I(0.9), size=I(2),inherit.aes = F, fill=I('white'))))+
  labs(y='RT for two targets (ms)', x='Condition')+
  scale_shape_ac()+
  guides(color='none', fill='none')+
  scale_x_discrete(limits = c('O+B','P+O','P+B','P+P'))+
  scale_y_exp()

scaleFUN <- function(x) sprintf("%.2f", x)
datam2_fp_bysubj<-datam2_fp[,.(p_first=mymean(p_first)),by=.(targets_stimType, subjectId)]

p3r<-ggplot(datam2_fp_bysubj[,data.frame(mean_cl_boot(p_first)), by=.(targets_stimType)], 
            aes(x=targets_stimType, y=y, ymin=ymin, ymax=ymax, group = 1))+
  list(geom_hline(yintercept = 0.5, linetype=2), 
       geom_violin(data=datam2_fp_bysubj, aes(y=p_first, x=targets_stimType), fill=I('black'),inherit.aes = F, alpha=I(0.1), color=I('white')), 
       geom_jitter(data=datam2_fp_bysubj, aes(y=p_first, x=targets_stimType), color=I('black'), width=0.15, alpha=I(0.9), size=I(2),inherit.aes = F, fill=I('white')))+
  labs(y='P(find first)', x='Target Type | Condition')+
  scale_shape_ac()+
  scale_x_discrete()+
  guides(color='none', fill='none')+
  scale_y_continuous()+
  expand_limits(y=c(0.45, 0.7))+
  geom_linerange(size=I(1))+
  geom_line(size=I(1))+
  geom_point(fill=I('white'))+
  coord_cartesian(ylim=c(0.4,0.8))

png(filename = 'Figure2.png', type='cairo', width = 9.4*280, height = 7*320, res = 320)
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
dev.off()

