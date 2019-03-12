

require(mice)
require(plyr)
options(stringsAsFactors = FALSE)

allow_base = c('a', 't', 'c', 'g', 'A', 'T', 'C', 'G', '.')
# header_names = c("chr", "name", "hg", "hierach", "pos", "mut", "anc", "der")
# Create a table/dataframe with the positions.txt
header_names = c("chr", "name", "hg", "pos", "mut", "anc", "der")
marker_table = read.table(Markerfile, col.names=header_names, fill=TRUE, na.strings="", comment.char="")
marker_table = cc(marker_table) #What is the meaning of this 'cc'
message(sprintf("%d valid markers provided.", nrow(marker_table)))
# Create a table/dataframe with the PileUp
header_names = c('chr', 'pos', 'refbase', 'reads', 'align', 'quality')

pileup_table = read.table(Pileupfile, col.names=header_names, fill=TRUE, na.strings='', comment.char="", quote="")
message(sprintf("%d positions after pileup Y.", nrow(pileup_table)))

pileup_table = pileup_table[(pileup_table$chr %in% unique(marker_table$chr)) & 
(pileup_table$reads >= Reads_thresh), ]


message(sprintf("%d positions passed reads number filter.", nrow(pileup_table)))
marker_pileup = merge(marker_table, pileup_table)
message(sprintf("%d positions are in your markers list.", nrow(marker_pileup)))

if(nrow(marker_pileup) == 0) {
    stop("None of the positions is in your markers list. Exit...")
}

marker_pileup$ANC = toupper(marker_pileup$anc)
marker_pileup$DER = toupper(marker_pileup$der)
marker_pileup$REFBASE = toupper(marker_pileup$refbase)
marker_pileup$ALIGN = toupper(marker_pileup$align)


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

getbase = function(mydat)
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

    if(max_perc > Base_majority & max_base %in% allow_base)
    {
        if(max_base == '.')
        {
            base = mydat$REFBASE
        }
        else
        {
            base = max_base
        }
    } else {
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

nmarkers = nrow(marker_pileup)
#-----------here still have 1 reads---------
baseinfo = NULL
for(i in 1:nmarkers) {
    baseinfo_i = getbase(marker_pileup[i, ])
    baseinfo = rbind(baseinfo, baseinfo_i)
}

marker_pileup = cbind(marker_pileup, baseinfo)

# ---------------------------------
eqANC = (marker_pileup$base == marker_pileup$ANC)
eqDER = (marker_pileup$base == marker_pileup$DER)
marker_pileup$branch = ifelse(eqANC & eqDER, 'E', ifelse(eqANC, 'A', ifelse(eqDER, 'D', NA)))
marker_pileup$quality = paste('Q', marker_pileup$quality)
marker_pileup = marker_pileup[, ! names(marker_pileup) %in% c('ANC', 'DER', 'REFBASE', 'ALIGN')]
marker_pileup = marker_pileup[which(! is.na(marker_pileup$branch)), ]
message(sprintf("%d markers gives you haplogroup information.", nrow(marker_pileup)))

colnames(marker_pileup)[3]  <- "marker_name"
colnames(marker_pileup)[4]  <- "haplogroup"
colnames(marker_pileup)[5]  <- "mutation"
colnames(marker_pileup)[14] <- "reads"
colnames(marker_pileup)[16] <- "call_perc"
colnames(marker_pileup)[17] <- "called_base"
colnames(marker_pileup)[18] <- "state"
new_marker_pileup <- marker_pileup[c(1:7,14,16:18)]

write.table(new_marker_pileup, file=Outputfile, sep='\t', quote=FALSE, row.names=FALSE)
