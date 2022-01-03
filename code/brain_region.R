
#################################
## output data: human,mac,h.both,m.both,pairs,both

setwd('/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/data')

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



############################## single region version
name = name[-which(name %in% c('housekeeping','broadly_expressed','tissue_restricted'))]
# name = c(name,'brain_general')
total = read.table('/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/macaque/data/hg19_expression_outcomes.bed.gz',stringsAsFactors=F)
colnames(total) = c('chr','start','end','Id','TPM','Strand','Human','Expression','Macaque','Macaque_Id','Mouse','Mouse_Id')
# h.aligned = total[which(total$Macaque!='unaligned'),]
brain = total[which(total$Macaque=='brain'),'Id']

hp = table(human$Id)
h.special = names(hp[which(hp==1)])
mp = table(mac$Id)
m.special = names(mp[which(mp==1)])

all <- num <- NULL
for (region in name){
  # if (region != 'brain_general'){
  #   h.reg = human[which(human$region==region),]
  #   m.reg = mac[which(mac$region==region),]
  # } else if (region == 'brain_general'){
  #   h.reg = human[!duplicated(human$Id),]
  #   m.reg = mac[!duplicated(mac$Id),]
  # }
  h.reg = human[which(human$region==region),]
  m.reg = mac[which(mac$region==region),]
  
  h.both = cbind(h.reg[,1:6],Brain_expression='non_specific')
  h.both[,7] = as.character(h.both[,7])
  h.both[which(h.both[,4] %in% h.special),7]='region_specific'
  m.both = cbind(m.reg[,1:6],Brain_expression='non_specific')
  m.both[,7] = as.character(m.both[,7])
  m.both[which(m.both[,4] %in% m.special),7]='region_specific'
  
  ##human promoter conservation
  rownames(total) = total$Id
  h.both$Macaque = 'unconserved'
  h.both$Macaque_Id = total[h.both$Id,]$Macaque_Id
  h.both$matched_Id = NA
  h.both[which(h.both$Id %in% brain),'Macaque'] = 'divergent'
  for (n in which(h.both$Macaque=='divergent')){
    mac_Id = total[h.both$Id[n],]$Macaque_Id
    mac_Id = strsplit(mac_Id,';')[[1]]
    idx = sapply(mac_Id,function(Id) Id %in% m.both$Id)
    # idx = ifelse(length(mac_Id)>1,idx[,1],idx)
    if ( sum(idx)>0 ) {
      h.both$Macaque[n] = 'matched'
      h.both$matched_Id[n] = Reduce(paste,mac_Id[which(idx)])
    } 
  }
  h.both$matched_Id = gsub(' ',';',h.both$matched_Id)
  
  ## macaque promoter conservation
  rownames(total) = NULL
  m.both$Human = 'unconserved'
  m.both$matched_Id <-m.both$Human_Id <-  NA
  
  m1 = lapply(m.reg$Id,function(x) total[grep(x,total$Macaque_Id),])
  for (n in 1:nrow(m.reg)){
    if (nrow(m1[[n]]) != 0) {
      m.both$Human[n] = 'divergent'
      Human_Id = Reduce(paste,m1[[n]]$Id)
      m.both$Human_Id[n] = gsub(' ',';',Human_Id)
      
      idx = sapply(m1[[n]]$Id,function(Id) Id %in% h.both$Id)
      if ( sum(idx)>0 ) {
        m.both$Human[n] = 'matched'
        m.both$matched_Id[n] = Reduce(paste,m1[[n]][which(idx),]$Id)
      }
    }
  }
  m.both$matched_Id = gsub(' ',';',m.both$matched_Id)
  
  ##order
  h.both = h.both[order(h.both[,1],h.both[,2],h.both[,3]),]
  h.both$chr = gsub('chr','',h.both$chr)
  m.both = m.both[order(m.both[,1],m.both[,2],m.both[,3]),]
  m.both$chr = gsub('chr','',m.both$chr)
  
  # if (region != 'brain_general'){
    all = rbind(all,h.both)
    num = c(num,nrow(h.both))
  # }
  
  write.table(h.both,paste('hg19',region,'expression.txt',sep='.'),row.names = F,quote=F, sep='\t')
  write.table(m.both,paste('rheMac3',region,'expression.txt',sep='.'),row.names = F,quote=F, sep='\t')
}


########################## method 1 - all_brain
both=NULL
for (Id in unique(all$Id)){
  data = all[which(all$Id ==Id),]
  if (nrow(data)>1){
    if ( sum(data$Macaque=='matched')>0 ) {data = data[which(data$Macaque=='matched'),][1,]
    }else if (sum(data$Macaque=='unconserved')>0) {data = data[which(data$Macaque=='unconserved'),][1,]
    } else data = unique(data[which(data$Macaque=='divergent'),][1,])
  }
  both = rbind(both,data)
}
both = both[order(both[,1],both[,2],both[,3]),]
write.table(both,'hg19.all_brain.expression.txt',row.names = F,quote=F, sep='\t')


######## downsampled all_brain
avenum = as.integer(mean(num))
Ids = sample(unique(both$Id),avenum)
down = data[which(data$Id %in% Ids),]
write.table(down,'hg19.all_brain_downsampled.expression.txt',row.names = F,quote=F, sep='\t')
for (n in 1:4){
  Ids = sample(unique(both$Id),avenum)
  down = data[which(data$Id %in% Ids),]
  write.table(down,paste0('./downsampled/hg19.all_brain_downsampled',n,'.expression.txt'),row.names = F,quote=F, sep='\t')
}

# # brain_general
# data = read.table('hg19.brain_general.expression.txt',stringsAsFactors = F,header=T)
# Ids = sample(unique(data$Id),avenum)
# down = data[which(data$Id %in% Ids),]
# write.table(down,'hg19.brain_general_downsampled.expression.txt',row.names = F,quote=F, sep='\t')
# for (n in 1:4){
#   Ids = sample(unique(data$Id),avenum)
#   down = data[which(data$Id %in% Ids),]
#   write.table(down,paste0('./downsampled/hg19.brain_general_downsampled',n,'.expression.txt'),row.names = F,quote=F, sep='\t')
# }

