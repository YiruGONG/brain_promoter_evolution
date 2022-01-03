filepath <- '/mnt/nfsstaging/cmvm/datastore/sbms/groups/young-lab/yiru/test'
filenames <- list.files(filepath,pattern='hg19.*.expression.txt$')

for (idx in 1:length(filenames)){
  data = read.table(filenames[idx],stringsAsFactors = F,header=T)
  data[which(data$Strand=='+'),2] = data[which(data$Strand=='+'),2]-150
  data[which(data$Strand=='+'),3] = data[which(data$Strand=='+'),3]+50
  data[which(data$Strand=='-'),2] = data[which(data$Strand=='-'),2]-50
  data[which(data$Strand=='-'),3] = data[which(data$Strand=='-'),3]+150
  
  data = data[order(data[,1],data[,2],data[,3]),]
  write.table(data,gsub('expression','expression.window',filenames[idx]),row.names = F,quote=F, sep='\t')
}