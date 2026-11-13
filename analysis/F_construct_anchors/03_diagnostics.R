# Author / creator / last modified by: Chunhao Guo
# post hoc follow-up (specified after modules B2 and D)
source('scripts/00_common.R')
suppressPackageStartupMessages(library(ggplot2))
# Original-design, single-panel diagnostics; no publication-reference reproduction claim.
theme_set(theme_classic(base_size=11,base_family='Helvetica')+theme(text=element_text(colour='black'),plot.margin=margin(10,12,10,10)))
d<-rd('results/F1_GENOTYPE_ALL.csv')[signature=='Angio36'&operator=='O2']
d[,comparison:=paste(cohort,anchor,sep=' / ')];d[,comparison:=factor(comparison,levels=rev(comparison))]
p<-ggplot(d,aes(difference_vs_O1,comparison))+geom_vline(xintercept=0,linetype=2,colour='grey55')+
 geom_errorbar(aes(xmin=difference_lower,xmax=difference_upper),width=.18,orientation='y',linewidth=.5)+
 geom_point(size=2)+labs(x='Angio36 O2 minus O1: difference in standardised genotype contrast',y=NULL)
s<-rd('results/F3_GROUP_ANCHORS_ALL.csv')[cohort=='CPTAC']
s[,gene_group:=factor(gene_group,levels=rev(c('Shared 6','Angio16-only','Angio36-only FGF/FGFR','Angio36-only other')))]
q<-ggplot(s,aes(estimate,gene_group,shape=operator,colour=operator))+
 geom_errorbar(aes(xmin=lower,xmax=upper),width=.18,position=position_dodge(width=.5),orientation='y',linewidth=.5)+
 geom_point(size=2,position=position_dodge(width=.5))+scale_colour_manual(values=c(O1='#267A9D',O2='#BD6719'))+
 labs(x='Spearman correlation with the four-protein vascular proxy (78 cases)',y=NULL,colour=NULL,shape=NULL)+theme(legend.position='top')
for(nm in c('F1_Angio36_O2_minus_O1','F3_CPTAC_group_anchors')){
 fig<-(if(grepl('^F1',nm))p else q)+labs(subtitle=LABEL)
 grDevices::pdf(file.path('plots',paste0(nm,'.pdf')),width=9,height=4.7,family='Helvetica',useDingbats=FALSE);print(fig);dev.off()
 ragg::agg_png(file.path('plots',paste0(nm,'.png')),width=9,height=4.7,units='in',res=180);print(fig);dev.off()
}
wr(d,'PLOT_SOURCE_F1.csv');wr(s,'PLOT_SOURCE_F3.csv')
writeLines(capture.output(sessionInfo()),'logs/03_plot_sessionInfo.txt')
