

#################################
## output data: human,mac,h.both,m.both,pairs

# rheMac3
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/macaque/data'
filenames <- list.files(filepath,pattern='rheMac3.*.outcomes.bed.gz$',full.names = T)
filenames <- filenames[-grep('rheMac3_',filenames)]
name = gsub(paste0(filepath,'/','rheMac3.'),'',filenames)
name = gsub('.outcomes.bed.gz','',name)

mac = lapply(filenames,function(x) read.table(x,stringsAsFactors = F))
# mac = lapply(mac,function(x) cbind(x,region=name))
for (i in 1:length(filenames)){
  mac[[i]]= cbind(mac[[i]],region=name[i])
}
mac = do.call(rbind,mac)
colnames(mac) = c('chr','start','end','Id','TPM','Strand','Brain_expression','Human','Mouse','region')
mac = mac[-which(mac$region %in% c('housekeeping','broadly_expressed','tissue_restricted')),]


##hg19
filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/macaque/data'
filenames <- list.files(filepath,pattern='hg19.*.outcomes.bed.gz$',full.names = T)
filenames <- filenames[-grep('hg19_',filenames)]
# filenames <- filenames[-grep('amygdala',filenames)]
name = gsub(paste0(filepath,'/','hg19.'),'',filenames)
name = gsub('.outcomes.bed.gz','',name)

human = lapply(filenames,function(x) read.table(x,stringsAsFactors = F))
# human = lapply(human,function(x) cbind(x,region=name))
for (i in 1:length(filenames)){
  human[[i]]= cbind(human[[i]],region=name[i])
}
human = do.call(rbind,human)
colnames(human) = c('chr','start','end','Id','TPM','Strand','Brain_expression','Body_expression','Macaque','Mouse','region')
human = human[-which(human$region %in% c('housekeeping','broadly_expressed','tissue_restricted')),]


# ######specificity
# #human
# hp = table(human$Id)
# special = names(hp[which(hp==1)])
# all = unique(human$Id)
# h.spe = cbind(all,expression='non-specific')
# h.spe[which(h.spe[,1] %in% special),2]='specific'
# try = human[!duplicated(human[,c(1:4,6)]),c(1:4,6)]
# rownames(try) = try[,4]
# h.spe = cbind(try[all,],expression=h.spe[,2])
# 
# region = sapply(all,function(idx){
#   each = human[which(human$Id==idx),]
#   region = Reduce(paste,each$region)
#   region = gsub(' ',';',region)
#   region
# })
# h.spe = cbind(h.spe,region)
# rownames(h.spe) = NULL

# #macaque
# mp = table(mac$Id)
# special = names(mp[which(mp==1)])
# all = unique(mac$Id)
# m.spe = cbind(all,expression='non-specific')
# m.spe[which(m.spe[,1] %in% special),2]='specific'
# try = mac[!duplicated(mac[,c(1:4,6)]),c(1:4,6)]
# rownames(try) = try[,4]
# m.spe = cbind(try[all,],expression=m.spe[,2])
# 
# region = sapply(all,function(idx){
#   each = mac[which(mac$Id==idx),]
#   region = Reduce(paste,each$region)
#   region = gsub(' ',';',region)
#   region
# })
# m.spe = cbind(m.spe,region)
# rownames(m.spe) = NULL


# ################################ all region version
# #####conservation
# name = name[-which(name %in% c('housekeeping','broadly_expressed','tissue_restricted'))]
# total = read.table('/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/macaque/data/hg19_expression_outcomes.bed.gz',stringsAsFactors=F)
# colnames(total) = c('chr','start','end','Id','TPM','Strand','Human','Expression','Macaque','Macaque_Id','Mouse','Mouse_Id')
# 
# pairs = NULL
# for (idx in name){
#   h.reg = human[which(human$region==idx),]
#   m.reg = mac[which(mac$region==idx),]
#   
#   m1 = lapply(m.reg$Id,function(x) total[grep(x,total$Macaque_Id),])
#   for (n in 1:nrow(m.reg)){
#     if (nrow(m1[[n]]) != 0) {m1[[n]]=cbind(m1[[n]],Macaque_region=m.reg$Id[n])}
#   }
#   m1 = do.call(rbind, m1)
#   
#   h1 = total[which(total$Id %in% h.reg$Id),]
#   both = m1[which(m1$Id %in% h1$Id),]
#   both = both[order(both$Id),]
#   
#   r.conserved = cbind(Human_Id=both$Id,Macaque_Id=as.character(both$Macaque_region),region=idx)
#   pairs = rbind(pairs,r.conserved)
# }
# 
# ##human
# conserved = unique(pairs[,1])
# Human_Id=unique(human$Id)
# h.con = cbind(Human_Id,alignment='unconserved')
# h.con[which(h.con[,1] %in% conserved),2]='conserved'
# 
# try = human[!duplicated(human[,c(1:4,6)]),c(1:4,6)]
# rownames(try) = try[,4]
# h.con = cbind(try[Human_Id,],alignment=h.con[,2])
# 
# Human_region = sapply(Human_Id,function(idx){
#   each = human[which(human$Id==idx),]
#   region = Reduce(paste,each$region)
#   region = gsub(' ',';',region)
#   region
# })
# h.con = cbind(h.con,Human_region)
# rownames(h.con) = NULL
# 
# h.con = cbind(h.con,Macaque_Id=NA,conserved_region=NA)
# # h.con2 = NULL
# for (idx in conserved){
#   each = pairs[which(pairs[,1]==idx),]
#   if (is.null(dim(each))) each=t(as.matrix(each))
#   
#   Macaque_Id = Reduce(paste,unique(each[,2]))
#   Macaque_Id = gsub(' ',';',Macaque_Id)
#   region = Reduce(paste,unique(each[,3]))
#   region = gsub(' ',';',region)
#   h.con[which(h.con$Id==idx),'Macaque_Id']=Macaque_Id
#   h.con[which(h.con$Id==idx),'conserved_region']=region
#   # this = cbind(Human_Id=idx,Macaque_Id,region)
#   # h.con2 = rbind(h.con2,this)
# }
# 
# ##macaque
# conserved = unique(pairs[,2])
# Macaque_Id=unique(mac$Id)
# m.con = cbind(Macaque_Id,alignment='unconserved')
# m.con[which(m.con[,1] %in% conserved),2]='conserved'
# 
# try = mac[!duplicated(mac[,c(1:4,6)]),c(1:4,6)]
# rownames(try) = try[,4]
# m.con = cbind(try[Macaque_Id,],alignment=m.con[,2])
# 
# Macaque_region = sapply(Macaque_Id,function(idx){
#   each = mac[which(mac$Id==idx),]
#   region = Reduce(paste,each$region)
#   region = gsub(' ',';',region)
#   region
# })
# m.con = cbind(m.con,Macaque_region)
# rownames(m.con) = NULL
# 
# m.con = cbind(m.con,Human_Id=NA,conserved_region=NA)
# # h.con2 = NULL
# for (idx in conserved){
#   each = pairs[which(pairs[,2]==idx),]
#   if (is.null(dim(each))) each=t(as.matrix(each))
#   
#   Human_Id = Reduce(paste,unique(each[,1]))
#   Human_Id = gsub(' ',';',Human_Id)
#   region = Reduce(paste,unique(each[,3]))
#   region = gsub(' ',';',region)
#   m.con[which(m.con$Id==idx),'Human_Id']=Human_Id
#   m.con[which(m.con$Id==idx),'conserved_region']=region
#   # this = cbind(Human_Id=idx,Macaque_Id,region)
#   # h.con2 = rbind(h.con2,this)
# }
# 
# ###synthesize
# h.both = cbind(h.spe[,1:6],h.con[,6:9])
# m.both = cbind(m.spe[,1:6],m.con[,6:9])
# 
# write.csv(h.both,'human_region_output.csv',row.names = F)
# write.csv(m.both,'macaque_region_output.csv',row.names = F)
# 
# save(human,mac,h.both,m.both,pairs,file='output.csv')




############################## single region version
name = name[-which(name %in% c('housekeeping','broadly_expressed','tissue_restricted'))]
name = c(name,'brain_general')
total = read.table('/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/macaque/data/hg19_expression_outcomes.bed.gz',stringsAsFactors=F)
colnames(total) = c('chr','start','end','Id','TPM','Strand','Human','Expression','Macaque','Macaque_Id','Mouse','Mouse_Id')
h.aligned = total[which(total$Macaque!='unaligned'),]
# m.id = strsplit(total$Macaque_Id,';')
# m.id = do.call(c,m.id)
# m.id = unique(na.omit(m.id))
# m.aligned = sapply(m.id,function(x){
#   
# })

#human
hp = table(human$Id)
h.special = names(hp[which(hp==1)])
mp = table(mac$Id)
m.special = names(mp[which(mp==1)])

all = NULL
for (region in name){
  if (region != 'brain_general'){
    h.reg = human[which(human$region==region),]
    m.reg = mac[which(mac$region==region),]
  } else if (region == 'brain_general'){
    h.reg = human[!duplicated(human$Id),]
    m.reg = mac[!duplicated(mac$Id),]
  }
  
  h.both = cbind(h.reg[,1:6],Brain_expression='non_specific')
  h.both[,7] = as.character(h.both[,7])
  h.both[which(h.both[,4] %in% h.special),7]='region_specific'
  m.both = cbind(m.reg[,1:6],Brain_expression='non_specific')
  m.both[,7] = as.character(m.both[,7])
  m.both[which(m.both[,4] %in% m.special),7]='region_specific'
  
  m1 = lapply(m.reg$Id,function(x) total[grep(x,total$Macaque_Id),])
  for (n in 1:nrow(m.reg)){
    if (nrow(m1[[n]]) != 0) {m1[[n]]=cbind(m1[[n]],Macaque_region=m.reg$Id[n])}
  }
  m1 = do.call(rbind, m1)
  h1 = total[which(total$Id %in% h.reg$Id),]
  both = m1[which(m1$Id %in% h1$Id),]
  both = both[order(both$Id),]
  r.conserved = cbind(Human_Id=both$Id,Macaque_Id=as.character(both$Macaque_region))
  
  h.both = cbind(h.both,Macaque='unconserved')
  h.both[,8] = as.character(h.both[,8])
  h.both[which(h.both[,4] %in% h.aligned$Id),8] = 'aligned'
  h.both[which(h.both[,4] %in% r.conserved[,1]),8]='conserved'
  m.both = cbind(m.both,Human='unconserved')
  m.both[,8] = as.character(m.both[,8])
  m.both[which(m.both[,4] %in% r.conserved[,2]),8]='conserved'
  
  h.both = cbind(h.both,Macaque_Id=NA)
  # Macaque_Id = lapply(unique(r.conserved[,1]),function(idx){
  for (idx in unique(r.conserved[,1])){
    each = r.conserved[which(r.conserved[,1]==idx),2]
    Macaque_Id = Reduce(paste,each)
    Macaque_Id = gsub(' ',';',Macaque_Id)
    h.both[which(h.both[,4]==idx),9]=Macaque_Id
  }
  m.both = cbind(m.both,Human_Id=NA)
  for (idx in unique(r.conserved[,2])){
    each = r.conserved[which(r.conserved[,2]==idx),1]
    Human_Id = Reduce(paste,each)
    Human_Id = gsub(' ',';',Human_Id)
    m.both[which(m.both[,4]==idx),9]=Human_Id
  }
  
  h.both = h.both[order(h.both[,1],h.both[,2],h.both[,3]),]
  h.both$chr = gsub('chr','',h.both$chr)
  m.both = m.both[order(m.both[,1],m.both[,2],m.both[,3]),]
  m.both$chr = gsub('chr','',m.both$chr)
  
  all = rbind(all,h.both)
  write.table(h.both,paste('hg19',region,'expression.txt',sep='.'),row.names = F,quote=F, sep='\t')
  write.table(m.both,paste('rheMac3',region,'expression.txt',sep='.'),row.names = F,quote=F, sep='\t')
}

## method 1 - all_brain
# filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/window'
# filenames <- list.files(filepath,pattern='hg19.*.expression.window.txt$')
# filenames <- filenames[-grep('brain_general',filenames)]
# # filenames <- filenames[-grep('all_brain',filenames)]
# names = gsub('.expression.window.txt','',filenames)
# names = gsub('hg19.','',names)

# # all=NULL
# # for (idx in 1:length(names)){
# #   data = read.table(filenames[idx],stringsAsFactors = F,header=T)
# #   data = cbind(data,region=names[idx])
# #   all = rbind(all,data)
# # }
# 
# all = lapply(filenames,function(x) read.table(x,stringsAsFactors = F,header=T))
# all = do.call(rbind,all)
# # both = all[!duplicated(all$Id),]
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

# matched = all[which(all$Macaque=='conserved'),]
# matched = matched[!duplicated(matched$Id),]
# divergent = all[which(all$Macaque == 'aligned'),]
# divergent = divergent[!duplicated(divergent$Id),]
# divergent = divergent[-which(divergent$Id %in% matched$Id),]
# matched = matched[order(matched[,1],matched[,2],matched[,3]),]
# divergent = divergent[order(divergent[,1],divergent[,2],divergent[,3]),]
# write.table(matched,'hg19.all_brain.matched.expression.window.txt',row.names = F,quote=F, sep='\t')
# write.table(divergent,'hg19.all_brain.divergent.expression.window.txt',row.names = F,quote=F, sep='\t')


###
# filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data/output1'
# filenames <- list.files(filepath,pattern='*.expression.txt$')
# # name = gsub('.expression.txt','',filenames)
# 
# for (idx in 1:length(filenames)){
#   data = read.table(filenames[idx],stringsAsFactors = F,header=T)
#   data = data[order(data[,1],data[,2],data[,3]),]
#   data$chr = gsub('chr','',data$chr)
#   write.table(data,filenames[idx],row.names = F,quote=F, sep='\t')
# }

# exp.col = colnames(h.both)
# count.col = c('chr','start','end','name','PASS','Reference','Alternate','Derived','Derived_freq','Alternate_freq','Type')
# count.col = matrix(count.col,nrow=1)
# write.table(count.col,'hg19.counts.colnames.txt',row.names = F,col.names = F,quote=F, sep='\t')

