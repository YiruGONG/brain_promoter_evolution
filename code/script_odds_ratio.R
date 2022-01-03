library(ggplot2)

# filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
# filenames <- list.files(filepath,pattern='*.divergent.counts.window.txt$')
# filenames <- c(filenames,list.files(filepath,pattern='*.matched.counts.window.txt$'))
# names = gsub('.counts.window.txt','',filenames)
# names = gsub('hg19.','',names)



### counts.txt
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/'
filenames <- list.files(filepath,pattern='hg19.*.counts.txt$')
names = gsub('.counts.txt','',filenames)
names = gsub('hg19.','',names)

example = read.table('hg19.amygdala.expression.txt',stringsAsFactors = F,header=T)
exp.col = colnames(example)
count.col = c('chr','start','end','name','PASS','Reference','Alternate','Derived','Derived_freq','Alternate_freq','Type')

rare_genome    = 12205724
common_genome  = 5263405




filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/'
filenames <- list.files(filepath,pattern='hg19.*.counts.window.txt$')
filenames <- filenames[-grep('divergent',filenames)]
filenames <- filenames[-grep('matched',filenames)]
names = gsub('.counts.window.txt','',filenames)
names = gsub('hg19.','',names)

example = read.table('hg19.amygdala.expression.txt',stringsAsFactors = F,header=T)
exp.col = colnames(example)
count.col = c('chr','start','end','name','PASS','Reference','Alternate','Derived','Derived_freq','Alternate_freq','Type')

rare_genome    = 12205724
common_genome  = 5263405

result=NULL
for (idx in 1:length(filenames)){
  data = read.table(filenames[idx],stringsAsFactors = F,header=F)
  colnames(data) = c(count.col,exp.col)
  rare_region = sum(data$Derived_freq<0.015)
  common_region  = sum(data$Derived_freq>0.05)
  odds.ratio = (rare_region/common_region)/(rare_genome/common_genome)
  fisher = fisher.test(matrix(c(rare_region, common_region, rare_genome, common_genome), ncol = 2))
  reg.result = cbind(rare_region,common_region,rare_genome,common_genome,odds_ratio=fisher$estimate[["odds ratio"]],lower_conf=fisher$conf.int[1],upper_conf=fisher$conf.int[2],p.value=fisher$p.value)
  result = rbind(result,reg.result)
}
result = as.data.frame(result)
result = cbind(brain_region=names,result)
result$brain_region = as.character(result$brain_region)
# result = result[order(result$brain_region),]
result$p.adjust = p.adjust(result$p.value,method='bonferroni')
result$region = 'region'
result$region[grep('brain',result$brain_region)]= 'brain'
result= result[order(result$region,result$brain_region),]
result$brain_region = factor(result$brain_region,levels=result$brain_region)
save(result,file='hg19.OR.window.Rdata')

# plot(result$odds_ratio,result$brain_region,las=1,pch='|')

# ggplot(result,aes(odds_ratio,brain_region,xmin=lower_conf,xmax=upper_conf))+
#   geom_pointrange(shape='|',color='black')+
#   geom_vline(xintercept = 1,linetype=2,color='lightblue')+
#   theme_classic()+
#   labs(y=NULL,x='Derived allele frequency\n odds ratio')
#   # coord_flip()

# pdf(file='/exports/cmvm/datastore/sbms/groups/young-lab/yiru/plots/hg19.regions.odds_ratio.pdf',width=6,height=6)
# ggplot(result,aes(brain_region,odds_ratio,ymin=lower_conf,ymax=upper_conf))+
#   geom_pointrange(shape='|')+  #,color='black')+
#   geom_hline(yintercept = 1,linetype=2,color='lightblue')+
#   theme_classic()+
#   scale_y_continuous(breaks=seq(0,4,0.5))+
#   labs(x=NULL,y='Derived allele frequency\n odds ratio')+
#   coord_flip()
# dev.off()


#### own computer
load('C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.window.Rdata')
result
result$p.adjust = p.adjust(result$p.value,method='bonferroni')
result$region = 'region'
result$region[grep('brain',result$brain_region)]= 'brain'
result= result[order(result$region,result$brain_region),]
result$brain_region = factor(result$brain_region,levels=result$brain_region)
save(result,file='C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.window.Rdata')
idx = which(result$p.value<0.05)

pdf(file='C:/Users/lenovo/Desktop/Young lab/output/hg19_region_OR.pdf',width=7,height=6)
ggplot(result,aes(brain_region,odds_ratio,ymin=lower_conf,ymax=upper_conf,color=region))+
  geom_pointrange(shape='|',size=0.8)+  #,color='black')+
  geom_hline(yintercept = 1,linetype=2,color='lightblue')+
  geom_vline(xintercept = 2.5,linetype=1,color='grey')+
  annotate('text',x=idx-0.05,y=result$upper_conf[idx]+0.25,label='*')+
  theme_classic()+
  scale_y_continuous(breaks=seq(0,4,0.5))+
  scale_color_manual(values=c('blue','black'))+
  theme(legend.position = 'none')+
  labs(x=NULL,y='Derived allele frequency\n odds ratio')+
  coord_flip()
dev.off()


##### manual adding window #####
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data'
filenames <- list.files(filepath,pattern='hg19.*.expression.txt$')
# names = gsub('.expression.txt','',filenames)
# names = gsub('hg19.','',names)

for (idx in 1:length(filenames)){
  data = read.table(filenames[idx],stringsAsFactors = F,header=T)
  data[which(data$Strand=='+'),2] = data[which(data$Strand=='+'),2]-150
  data[which(data$Strand=='+'),3] = data[which(data$Strand=='+'),3]+50
  data[which(data$Strand=='-'),2] = data[which(data$Strand=='-'),2]-50
  data[which(data$Strand=='-'),3] = data[which(data$Strand=='-'),3]+150
  
  data = data[order(data[,1],data[,2],data[,3]),]
  write.table(data,gsub('expression','expression.window',filenames[idx]),row.names = F,quote=F, sep='\t')
}

# ## function mode
# window = function(data,l=150,r=50,strand=T){
#   if (strand){
#     data[which(data$Strand=='+'),2] = data[which(data$Strand=='+'),2]-l
#     data[which(data$Strand=='+'),3] = data[which(data$Strand=='+'),3]+r
#     data[which(data$Strand=='-'),2] = data[which(data$Strand=='-'),2]-r
#     data[which(data$Strand=='-'),3] = data[which(data$Strand=='-'),3]+l
#   } else {
#     data[,2] = data[,2]-l
#     data[,3] = data[,3]+r
#   }
#   data[order(data[,1],data[,2],data[,3]),]
# }
# 
# filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
# filenames <- list.files(filepath,pattern='hg19.*.expression.txt$')
# names = gsub('.expression.txt','',filenames)
# names = gsub('hg19.','',names)
# 
# for (idx in 1:length(filenames)){
#   data = read.table(filenames[idx],stringsAsFactors = F,header=T)
#   data = window(data,l=150,r=50,strand=T)
#   write.table(data,paste0('hg19.',names[idx],'.expression.window.txt'),row.names = F,quote=F, sep='\t')
# }



###### conserved vs others ########
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
filenames <- list.files(filepath,pattern='hg19.*.counts.window.txt$')
filenames <- filenames[-grep('divergent',filenames)]
filenames <- filenames[-grep('matched',filenames)]
names = gsub('.counts.window.txt','',filenames)
# names = gsub('hg19.','',names)

for (idx in 1:length(filenames)){
  data = read.table(filenames[idx],stringsAsFactors = F,header=F)
  matched = data[which(data[,19]=='conserved'),]
  divergent = data[which(data[,19]=='aligned'),]
  write.table(matched,paste0(names[idx],'.matched.counts.window.txt'),row.names = F,col.names=F,quote=F, sep='\t')
  write.table(divergent,paste0(names[idx],'.divergent.counts.window.txt'),row.names = F,col.names = F,quote=F, sep='\t')
}


### own computer plot
load('C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.matched.window.Rdata')
result
result[which(result$p.value<0.05),]
# result$region = 'region'
# result$region[grep('brain',result$brain_region)]= 'brain'
# result= result[order(result$region,result$brain_region),]
# result$brain_region = factor(result$brain_region,levels=result$brain_region)
# save(result,file='C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.matched.window.Rdata')

idxs = strsplit(as.character(result$brain_region),"\\.")
# region = sapply(idxs,function(x) x[1])
promoter = sapply(idxs,function(x) x[2])
result = cbind(result,promoter)
result$Promoter = paste(result$region,result$promoter)
# result$promoter = 'divergent'
# result[grep('matched',result$brain_region),'promoter']='matched'
# result$promoter = factor(result$promoter,levels=c('divergent','matched'))

pdf(file='C:/Users/lenovo/Desktop/Young lab/output/hg19_matched_OR_window.pdf',width=10,height=8)
ggplot(result,aes(brain_region,odds_ratio,ymin=lower_conf,ymax=upper_conf,color=Promoter))+ # ,color=region,type=promoter))+
  geom_pointrange(shape='|',size=1)+  #,color='black')+
  geom_hline(yintercept = 1,linetype=2,color='lightblue')+
  geom_vline(xintercept = 4.5,linetype=1,color='grey')+
  theme_classic()+
  scale_y_continuous(breaks=c(seq(0,5,1),seq(10,25,5)))+
  labs(x=NULL,y='Derived allele frequency\n odds ratio')+
  coord_flip()
  # scale_linetype_manual(values=c(3,1))
  # scale_size_manual(values=c(0.5,1))
dev.off()



### change the background
##on eddie
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
filenames <- list.files(filepath,pattern='*.divergent.counts.window.txt$')
filenames <- c(filenames,list.files(filepath,pattern='*.matched.counts.window.txt$'))
# names = gsub('.counts.window.txt','',filenames)
# names = gsub('hg19.','',names)
idxs = strsplit(filenames,"\\.")
names = unique(sapply(idxs,function(x) x[2]))

example = read.table('hg19.amygdala.expression.txt',stringsAsFactors = F,header=T)
exp.col = colnames(example)
count.col = c('chr','start','end','name','PASS','Reference','Alternate','Derived','Derived_freq','Alternate_freq','Type')

result=NULL
for (idx in 1:length(names)){
  matched = read.table(grep(paste0(names[idx],'.matched'),filenames,value=T),stringsAsFactors = F,header=F)
  divergent = read.table(grep(paste0(names[idx],'.divergent'),filenames,value=T),stringsAsFactors = F,header=F)
  colnames(matched) <- colnames(divergent) <- c(count.col,exp.col)
  rare_matched = sum(matched$Derived_freq<0.015)
  common_matched  = sum(matched$Derived_freq>0.05)
  rare_divergent = sum(divergent$Derived_freq<0.015)
  common_divergent  = sum(divergent$Derived_freq>0.05)
  
  odds.ratio = (rare_matched/common_matched)/(rare_divergent/common_divergent)
  fisher = fisher.test(matrix(c(rare_matched, common_matched, rare_divergent, common_divergent), ncol = 2))
  reg.result = cbind(rare_matched,common_matched,rare_divergent,common_divergent,odds_ratio=fisher$estimate[["odds ratio"]],lower_conf=fisher$conf.int[1],upper_conf=fisher$conf.int[2],p.value=fisher$p.value)
  result = rbind(result,reg.result)
}
result = as.data.frame(result)
result = cbind(brain_region=names,result)
result$brain_region = as.character(result$brain_region)
# result$p.adjust = p.adjust(result$p.value,method = 'bonferroni')
save(result,file='hg19.OR.matched_divergent.window.Rdata')

##own computer
load('C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.matched_divergent.window.Rdata')
# result$region = 'region'
# result$region[grep('brain',result$brain_region)]= 'brain'
# result= result[order(result$region,result$brain_region),]
# result$brain_region = factor(result$brain_region,levels=result$brain_region)
# save(result,file='C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.matched_divergent.window.Rdata')

idx = which(result$p.value<0.05)
pdf(file='C:/Users/lenovo/Desktop/Young lab/output/hg19_matched_divergent_OR_window.pdf',width=7,height=5)
ggplot(result,aes(brain_region,odds_ratio,ymin=lower_conf,ymax=upper_conf,color=region))+
  geom_pointrange(shape='|')+  #,color='black')+
  geom_hline(yintercept = 1,linetype=2,color='lightblue')+
  geom_vline(xintercept = 4.5,linetype=1,color='grey')+
  annotate('text',x=idx-0.05,y=result$upper_conf[idx]+1,label='*')+
  theme_classic()+
  scale_color_manual(values=c('blue','black'))+
  scale_y_continuous(breaks=c(seq(0,5,1),seq(10,25,5)))+
  labs(x=NULL,y='Derived allele frequency\n Matched/Divergent odds ratio')+
  coord_flip()
dev.off()





####### all regions compared ########
## method 1 - all_brain
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
filenames <- list.files(filepath,pattern='hg19.*.expression.window.txt$')
filenames <- filenames[-grep('brain_general',filenames)]
filenames <- filenames[-grep('all_brain',filenames)]
names = gsub('.expression.window.txt','',filenames)
names = gsub('hg19.','',names)

# all=NULL
# for (idx in 1:length(names)){
#   data = read.table(filenames[idx],stringsAsFactors = F,header=T)
#   data = cbind(data,region=names[idx])
#   all = rbind(all,data)
# }
all = lapply(filenames,function(x) read.table(x,stringsAsFactors = F,header=T))
all = do.call(rbind,all)
# both = all[!duplicated(all$Id),]
both=NULL
for (Id in unique(all$Id)){
  data = all[which(all$Id ==Id),]
  if (nrow(data)>1){
    if ( sum(data$Macaque=='conserved')>0 ) {data = data[which(data$Macaque=='conserved'),][1,]
    }else if (sum(data$Macaque=='unconserved')>0) {data = data[which(data$Macaque=='unconserved'),][1,]
    } else data = unique(data[which(data$Macaque=='aligned'),][1,])
  }
  both = rbind(both,data)
}
both = both[order(both[,1],both[,2],both[,3]),]
write.table(both,'hg19.all_brain.expression.txt',row.names = F,quote=F, sep='\t')

matched = all[which(all$Macaque=='conserved'),]
matched = matched[!duplicated(matched$Id),]
divergent = all[which(all$Macaque == 'aligned'),]
divergent = divergent[!duplicated(divergent$Id),]
divergent = divergent[-which(divergent$Id %in% matched$Id),]
matched = matched[order(matched[,1],matched[,2],matched[,3]),]
divergent = divergent[order(divergent[,1],divergent[,2],divergent[,3]),]
write.table(matched,'hg19.all_brain.matched.expression.window.txt',row.names = F,quote=F, sep='\t')
write.table(divergent,'hg19.all_brain.divergent.expression.window.txt',row.names = F,quote=F, sep='\t')


## method 2 - brain_general


### downsampling of brain region
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/'
filenames <- list.files(filepath,pattern='hg19.*.expression.txt$')
filenames1 <- filenames[-grep('brain',filenames)]
filenames2 <- filenames[grep('brain',filenames)]

num = sapply(filenames1,function(file) {
  data = read.table(file,stringsAsFactors = F,header=T)
  nrow(data)
})
avenum = as.integer(mean(num))

idxs = strsplit(filenames2,"\\.")
names = unique(sapply(idxs,function(x) x[2]))
for (x in 1:length(filenames2)){
  for (n in 1:5){
    data = read.table(filenames2[x],stringsAsFactors = F,header=T)
    Ids = sample(unique(data$Id),avenum)
    down = data[which(data$Id %in% Ids),]
    write.table(down,paste0('./downsampled/hg19.',names[x],'.downsampled',n,'.expression.txt'),row.names = F,quote=F, sep='\t')
  }
}


### own computer
load('C:/Users/lenovo/Desktop/Young lab/data/hg19.OR.brain.downsampled.window.Rdata')

