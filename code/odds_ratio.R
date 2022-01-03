
############### function ##################

odds_ratio = function(filenames,names,background=c('matched','divergent','genome'),rare_genome=12205724,common_genome=5263405){
  result=NULL
  if (background=='genome'){
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
  } else if (background=='divergent'){
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
  } else if (background=='matched'){
    for (idx in 1:length(names)){
      matched = read.table(grep(paste0(names[idx],'.matched'),filenames,value=T),stringsAsFactors = F,header=F)
      divergent = read.table(grep(paste0(names[idx],'.divergent'),filenames,value=T),stringsAsFactors = F,header=F)
      colnames(matched) <- colnames(divergent) <- c(count.col,exp.col)
      rare_matched = sum(matched$Derived_freq<0.015)
      common_matched  = sum(matched$Derived_freq>0.05)
      rare_divergent = sum(divergent$Derived_freq<0.015)
      common_divergent  = sum(divergent$Derived_freq>0.05)
      
      odds.ratio = (rare_divergent/common_divergent)/(rare_matched/common_matched)
      fisher = fisher.test(matrix(c(rare_divergent, common_divergent,rare_matched, common_matched), ncol = 2))
      reg.result = cbind(rare_matched,common_matched,rare_divergent,common_divergent,odds_ratio=fisher$estimate[["odds ratio"]],lower_conf=fisher$conf.int[1],upper_conf=fisher$conf.int[2],p.value=fisher$p.value)
      result = rbind(result,reg.result)
    }
  } else warning('Uncorrect Background selection')
  
  result = as.data.frame(result)
  result = cbind(brain_region=names,result)
  result$brain_region = as.character(result$brain_region)
  # result = result[order(result$brain_region),]
  result$p.adjust = p.adjust(result$p.value,method='bonferroni')
  result$region = 'region'
  result$region[grep('brain',result$brain_region)]= 'brain'
  result= result[order(result$region,result$brain_region),]
  result$brain_region = factor(result$brain_region,levels=result$brain_region)
  result
}


################## application #######################
setwd('/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/')
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/'

example = read.table('hg19.amygdala.expression.txt',stringsAsFactors = F,header=T)
exp.col = colnames(example)
count.col = c('chr','start','end','name','PASS','Reference','Alternate','Derived','Derived_freq','Alternate_freq','Type')

### .count.txt
filenames <- list.files(filepath,pattern='hg19.*.counts.txt$')
names = gsub('.counts.txt','',filenames)
names = gsub('hg19.','',names)

rare_genome    = 12205724
common_genome  = 5263405
result = odds_ratio(filenames,names,background = 'genome',rare_genome=rare_genome,common_genome=common_genome)
save(result,file='hg19.OR.Rdata')


### .count.window.txt
filenames <- list.files(filepath,pattern='hg19.*.counts.window.txt$')
filenames <- filenames[-grep('divergent',filenames)]
filenames <- filenames[-grep('matched',filenames)]
names = gsub('.counts.window.txt','',filenames)
names = gsub('hg19.','',names)

rare_genome    = 12205724
common_genome  = 5263405
result = odds_ratio(filenames,names,background = 'genome',rare_genome,common_genome)
save(result,file='hg19.OR.window.Rdata')


### subsets
filenames <- list.files(filepath,pattern='*.divergent.counts.window.txt$')
filenames <- c(filenames,list.files(filepath,pattern='*.matched.counts.window.txt$'))
names = gsub('.counts.window.txt','',filenames)
names = gsub('hg19.','',names)

rare_genome    = 12205724
common_genome  = 5263405
result = odds_ratio(filenames,names,background = 'genome',rare_genome,common_genome)
save(result,file='hg19.OR.subsets.window.Rdata')


### subset - change background to divergent
filenames <- list.files(filepath,pattern='*.divergent.counts.window.txt$')
filenames <- c(filenames,list.files(filepath,pattern='*.matched.counts.window.txt$'))
names = gsub('.divergent.counts.window.txt','',filenames)
names = gsub('.matched.counts.window.txt','',names)
names = gsub('hg19.','',names)
names = unique(names)
# idxs = strsplit(filenames,"\\.")
# names = unique(sapply(idxs,function(x) x[2]))

result = odds_ratio(filenames,names,background = 'divergent')
save(result,file='hg19.OR.matched_divergent.window.Rdata')


