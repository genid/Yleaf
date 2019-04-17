require(mice, quietly = T)
require(plyr, quietly = T)

options(stringsAsFactors = FALSE)
df_no_informative <- data.frame()

filelog = file(log_output, "w")
allow_base = c('a', 't', 'c', 'g', 'A', 'T', 'C', 'G', '.')
# header_names = c("chr", "name", "hg", "hierach", "pos", "mut", "anc", "der")
# Create a table/dataframe with the positions.txt
header_names = c("chr", "name", "hg", "pos", "mut", "anc", "der")
marker_table = read.table(Markerfile, col.names=header_names, fill=TRUE, na.strings="", comment.char="")

columnames <- c("chr","pos","marker_name","haplogroup","mutation","anc","der","reads","call_perc",
                "called_base","state","description")

marker_table = cc(marker_table) #What is the meaning of this 'cc'

tmp = (sprintf("Total of reads %s.", (total_reads)))
message(tmp)
write(tmp, filelog, append= T)
tmp = (sprintf("%d valid markers provided.", nrow(marker_table)))
message(tmp)
write(tmp, filelog, append= T)

# Create a table/dataframe with the PileUp
header_names = c('chr', 'pos', 'refbase', 'reads', 'align', 'quality')

pileup_table = read.table(Pileupfile, col.names=header_names, fill=TRUE, na.strings='', comment.char="", quote="")


comb_indel_begin = function(align)
{
  while(any(align == '^'))
  {
    p0 = which(align == '^')[1]
    p1 = p0 + 1
    # begin = paste(align[p0:p1], collapse='')
    align = align[-(p0:p1)]
  }
  while(any(align == '+'))
  {
    p0 = which(align == '+')[1]
    nDigits = sum(align[(p0+1):(p0+3)] %in% 0:9)
    nIns = as.integer(paste(align[(p0+1):(p0+nDigits)], collapse=''))
    p1 = p0 + nDigits + nIns
    align = align[-(p0:p1)]
  }
  while(any(align == '-'))
  {
    p0 = which(align == '-')[1]
    nDigits = sum(align[(p0+1):(p0+3)] %in% 0:9)
    nIns = as.integer(paste(align[(p0+1):(p0+nDigits)], collapse=''))
    p1 = p0 + nDigits + nIns
    align = align[-(p0:p1)]
  }
  align = align[which(align != '$')]
  return(align)
}

getbase = function(mydat, Base_majority)
{
  
  myalign = strsplit(mydat$ALIGN, '')[[1]]
  myalign[myalign == ','] = '.'
  myalign = comb_indel_begin(myalign)
  qualraw = charToRaw(mydat$quality)
  if(length(myalign) != length(qualraw)) {
    stop('align and quality not of the same length!\n')
  }
  
  align_tab = prop.table(table(myalign))
  max_perc = max(align_tab) * 100
  max_base = names(align_tab)[which.max(align_tab)]
  if(max_perc >= Base_majority & max_base %in% allow_base)
  {
    if(max_base == '.')
    {
      base = mydat$REFBASE
    }
    else
    {
      base = max_base
    }
  }else{
    base = NA
  }
  myalign = paste(myalign, collapse='')
  f_qual = paste('Q', rawToChar(qualraw))
  # print(align_tab)
  max_perc <- as.integer(max_perc)
  res = data.frame(f_align=myalign,
                   f_qual=f_qual,
                   f_len=nchar(myalign),
                   max_base=paste('M', max_base),
                   max_perc=max_perc, base=base)
  return(res)
}

pi <-  pileup_table[1:2] 
po <-  marker_table[c(1,4)] 
pi <- as.array(pi$pos)
po <- as.array(po$pos)

diff_position_pileup <- setdiff(po,pi)
diff_position_pileup <- as.data.frame(diff_position_pileup)
diff_position_pileup <- cbind(rep("chrY",nrow(diff_position_pileup)), diff_position_pileup$diff_position_pileup)
colnames(diff_position_pileup) <- c("chr", "pos")

diff_position_pileup <- merge(marker_table,diff_position_pileup )

if(nrow(diff_position_pileup) > 0){
  tmp = (sprintf("Position with zero reads"))
  diff_position_pileup$reads <- 0
  diff_position_pileup$max_perc <- NA
  diff_position_pileup$base <- NA
  diff_position_pileup$branch <- NA
  diff_position_pileup$description <- tmp
}

positions_zero <- pileup_table[which(pileup_table$reads == 0),]
df_zero_reads = merge(marker_table, positions_zero)
if(nrow(df_zero_reads) > 0){
  tmp = (sprintf("Position with zero reads"))
  df_zero_reads$max_perc <- NA
  df_zero_reads$base <- NA
  df_zero_reads$branch <- NA
  df_zero_reads$description <- tmp
  df_zero_reads = cbind(df_zero_reads[c(1:7,9,12:15)])
  pileup_table <- pileup_table[which(pileup_table$reads > 0),] 
}

if(nrow(diff_position_pileup) > 0 & nrow(df_zero_reads) > 0){
  df_zero_reads = rbind(diff_position_pileup,df_zero_reads)
  tmp = (sprintf("%d Position with zero reads", nrow(df_zero_reads)))
  message(tmp)
  write(tmp, filelog, append= T)
  df_no_informative <- rbind(df_no_informative, df_zero_reads)
  
}else if(nrow(diff_position_pileup)== 0 & nrow(df_zero_reads) > 0){
  tmp = (sprintf("%d Position with zero reads", nrow(df_zero_reads)))
  message(tmp)
  write(tmp, filelog, append= T)
  df_no_informative <- rbind(df_no_informative, df_zero_reads)
  
}else if(nrow(diff_position_pileup) > 0 & nrow(df_zero_reads) == 0){
  df_zero_reads <- diff_position_pileup
  tmp = (sprintf("%d Position with zero reads", nrow(df_zero_reads)))
  message(tmp)
  write(tmp, filelog, append= T)
  df_no_informative <- rbind(df_no_informative, df_zero_reads)
}

pileup_table_below_reads <- pileup_table[which(pileup_table$reads < Reads_thresh),]

df_reads_below = merge(marker_table, pileup_table_below_reads)

if(nrow(df_reads_below) > 0){
  tmp = (sprintf("%d positions below the r (reads) %d threshold", nrow(df_reads_below), Reads_thresh))
  message(tmp)
  write(tmp, filelog, append= T)
  
  tmp = (sprintf("Position below the r (reads) %d threshold", Reads_thresh))
  
  df_reads_below$ANC = toupper(df_reads_below$anc)
  df_reads_below$DER = toupper(df_reads_below$der)
  df_reads_below$REFBASE = toupper(df_reads_below$refbase)
  df_reads_below$ALIGN = toupper(df_reads_below$align)
  
  nmarkers = nrow(df_reads_below)
  baseinfo = NULL
  
  for(i in 1:nmarkers) {
    baseinfo_i = getbase(df_reads_below[i, ], 51)
    baseinfo = rbind(baseinfo, baseinfo_i)
  }
  df_reads_below = cbind(df_reads_below, baseinfo)
  #d <- df_reads_below
  #df_reads_below <- d
  eqANC = (df_reads_below$base == df_reads_below$ANC)
  eqDER = (df_reads_below$base == df_reads_below$DER)
  df_reads_below$branch = ifelse(eqANC & eqDER, 'E', ifelse(eqANC, 'A', ifelse(eqDER, 'D', NA)))
  df_reads_below$quality = paste('Q', df_reads_below$quality)
  df_reads_below$description <- tmp
  df_reads_below <- cbind(df_reads_below[c(1:7,9,20:23)])
  df_no_informative <- rbind(df_no_informative, df_reads_below)
  
}

marker_pileup = merge(marker_table, pileup_table)

if(nrow(marker_pileup) == 0) {
  tmp = "None of the positions is in your markers list. Exit..."
  write(tmp, filelog, append=T)
  stop(tmp)
}else{
  tmp = (sprintf("%d position(s) is in your markers list.", nrow(marker_pileup)))
  message(tmp)
  write(tmp, filelog, append= T)
}

pileup_table <- pileup_table[(pileup_table$chr %in% unique(marker_table$chr)) & 
                               (pileup_table$reads >= Reads_thresh), ]


if(nrow(pileup_table) == 0) {
  tmp = "No positions gave haplogroup information..."
  write(tmp, filelog, append=T)
  df = rbind(df_zero_reads, df_reads_below)  
  write.table(df,file=emf_output, sep='\t', quote=FALSE, row.names=FALSE)
  stop(tmp)
}

marker_pileup = merge(marker_table, pileup_table)

marker_pileup$ANC = toupper(marker_pileup$anc)
marker_pileup$DER = toupper(marker_pileup$der)
marker_pileup$REFBASE = toupper(marker_pileup$refbase)
marker_pileup$ALIGN = toupper(marker_pileup$align)

nmarkers = nrow(marker_pileup)
baseinfo = NULL
base_info_below = 0
discordant_genotype = 0

if(nmarkers > 0){
  for(i in 1:nmarkers) {
    baseinfo_i = getbase(marker_pileup[i, ], Base_majority)
    #print(baseinfo_i)
    baseinfo = rbind(baseinfo, baseinfo_i)
  }
  marker_pileup = cbind(marker_pileup, baseinfo)
  base_info_below = marker_pileup[which((marker_pileup$max_perc < Base_majority)),]
  
}

if(typeof(base_info_below) == "list"){
  
  if(nrow( base_info_below) > 0){
    base_info_below <- cbind(base_info_below[c(1:15)])
    nmarkers = nrow(base_info_below)
    baseinfo = NULL
    
    for(i in 1:nmarkers) {
      baseinfo_i = getbase(base_info_below[i, ], 51)
      #print(baseinfo_i)
      baseinfo = rbind(baseinfo, baseinfo_i)
    }
    
    base_info_below = cbind(base_info_below, baseinfo)
    eqANC = (base_info_below$base == base_info_below$ANC)
    eqDER = (base_info_below$base == base_info_below$DER)
    base_info_below$branch = ifelse(eqANC & eqDER, 'E', ifelse(eqANC, 'A', ifelse(eqDER, 'D', NA)))
    base_info_below$quality = paste('Q', base_info_below$quality)
    tmp = (sprintf("Below the b threshold %d", Base_majority))
    
    base_info_below$description <- tmp
    base_info_below <- cbind(base_info_below[c(1:7,9,20:23)])
    df_no_informative <- rbind(df_no_informative, base_info_below)
  }
}else{
  tmp = (sprintf("%d Below the b threshold %d", nrow(base_info_below), Base_majority))
  message(tmp)
  write(tmp, filelog, append= T)
}

if(nrow(marker_pileup) > 0){
  
  eqANC = (marker_pileup$base == marker_pileup$ANC)
  eqDER = (marker_pileup$base == marker_pileup$DER)
  marker_pileup$branch = ifelse(eqANC & eqDER, 'E', ifelse(eqANC, 'A', ifelse(eqDER, 'D', NA)))
  marker_pileup$quality = paste('Q', marker_pileup$quality)
  
  marker_pileup = marker_pileup[which((marker_pileup$max_perc >= Base_majority)),]
  marker_pileup = marker_pileup[, ! names(marker_pileup) %in% c('ANC', 'DER', 'REFBASE', 'ALIGN')]
  discordant_genotype = marker_pileup[which( is.na(marker_pileup$branch)), ]
  marker_pileup = marker_pileup[which(! is.na(marker_pileup$branch)), ]
  tmp = (sprintf("%d markers gives you haplogroup information", nrow(marker_pileup)))
  message(tmp)
  write(tmp, filelog, append = T)
  
  colnames(marker_pileup)[3]  <- "marker_name"
  colnames(marker_pileup)[4]  <- "haplogroup"
  colnames(marker_pileup)[5]  <- "mutation"
  colnames(marker_pileup)[14] <- "reads"
  colnames(marker_pileup)[16] <- "call_perc"
  colnames(marker_pileup)[17] <- "called_base"
  colnames(marker_pileup)[18] <- "state"
  new_marker_pileup <- marker_pileup[c(1:7,14,16:18)]
  new_marker_pileup <- new_marker_pileup[order(new_marker_pileup$haplogroup),]
  write.table(new_marker_pileup, file=Outputfile, sep='\t', quote=FALSE, row.names=FALSE)  
  
}else{
  tmp = (sprintf("0 markers gives you haplogroup information"))
  message(tmp)
  write(tmp, filelog, append = T)
}

if(typeof(discordant_genotype) == "list"){
  if(nrow(discordant_genotype)){
    tmp = (sprintf("%d Discordant genotype(s)" , nrow(discordant_genotype)))
    message(tmp)
    write(tmp, filelog, append= T)
    discordant_genotype <- discordant_genotype[c(1:7,14,16:18)]
    colnames(discordant_genotype)[8] <- "reads"
    discordant_genotype$description <- "Discordant genotype"
    df_no_informative <- rbind(df_no_informative, discordant_genotype)  
  }
}

if(nrow(df_no_informative) > 0){
  colnames(df_no_informative) <- columnames
  write.table(df_no_informative,file=emf_output, sep='\t', quote=FALSE, row.names=FALSE)
  tmp = (sprintf("%d markers did not give information haplogroup", nrow(df_no_informative)))
  write(tmp, filelog, append = T)
  message(tmp)
  close(filelog)    
  
}else{
  tmp = (sprintf("No info. available for this sample." ))
  message(tmp)
  write(tmp, emf_output, append = T)
  close(emf_output)
}