library(data.table)

online.monitor.dir <- '/home/sasha/EcolePolytechique/Code_R/ModLanGauss'
pedestal.suppression <- 0.05

map <- fread(paste0(online.monitor.dir,'/fev10_chip_channel_x_y_mapping.txt')) # data.table with chip:channel <-> X-Y mapping
## load(file=paste0(online.monitor.dir,'/map.RData'))
map[,i:=channel]

load.raw <- function(file) {
    cat('Loading file ',file,'\n')
    df.names <- c('acq','bx','sca','chip','i','adc','trig')
    df.classes <- c(rep('integer',6), 'logical')

    ## hits.local <- fread(paste0(online.monitor.dir,'/raw/raw ',file),
    ##                     col.names=c('acq','bx','chip','sca',paste0(c('h','hadc','l','ladc'),'.',rep(1:64,each=4))))
    ## hits.local <- melt(hits.local,id=1:4,measure.vars=patterns('^h\\.','^hadc','^l\\.','^ladc'),variable.name='i')
    ## setnames(hits.local, c(names(hits.local)[1:5],'trig','adc','l.trig','l.adc'))
    ## setkey(hits.local, acq,chip,sca,bx)
    ## hits.local[,i:=as.integer(i)]

    hits.local <- fread(paste0(online.monitor.dir,'/raw/raw ',file,' 0 high_gain_triggers ', pedestal.suppression),
                        col.names=df.names, colClasses=df.classes, logical01=TRUE,
                        key='acq,chip,sca,bx')

    ## hits.local <- data.table(read.table(pipe(paste0(online.monitor.dir,'/raw/raw ',file,' 0 high_gain_triggers ',
    ##                                                 pedestal.suppression)),
    ##                                     col.names=df.names, colClasses=df.classes, comment.char=''),
    ##                          key='acq,chip,sca,bx')
    if (nrow(hits.local) > 0) {
	hits.local[,bx.cor:=cumsum( {
	    dbx = c(0,diff(bx))
	    dbx < 0 | ( dbx == 0 & c(0,diff(sca))>0 )
	       })*4096+bx, by=list(acq,chip)] # For each spill,chip: sort in SCA: bx.cor = bx + N*4096, where N - number of detected oveflows,
					# ie. number of times when BX goes down, eg. maximally: from 4095 to 0.
					# Note, if bx jumps by >4096 crossings, this is not correct (but this is the best one can do with only 12 bits for BX).
					# dbx < 0 | ( dbx == 0 & c(0,diff(sca))>0 ): treat a jump by 4096 BXs as a special case:
					# in this case BX stays the same (dbx==0) but SCA changes

	## bx.cor correction is done per chip; BX in different chips may still be compared only modulo(4096): eg. in the same event the
	## jump may be detected in one chip, but not in the other. Ie. across chips one uses bx (== bx.cor modulo(4096)), within the chip - bx.cor.

	ev.local <- hits.local[,list(n.trig=sum(trig)), keyby=list(acq,bx)]  # Find retriggers across all chips in dif
	ev.local[,bx.group:=cumsum( c(0,diff(bx)) > 4 ), by=list(acq)]  # c(0,diff(bx)) = distances to previous BX;
	## nice trick to group bxid's differing by at most 3 ("successive")
	ev.local[,`:=`(nbx=.N, ibx=1:.N), by=list(acq,bx.group)]   # finally, number of "successive" BXs in each group, ibx>1 means retriggerings

	## In FEV10, sometimes there are events with negative signals (eg. ADC==4). Pedestals in such events
	## are shifted. To avoid errors, pedestals are calculated only if there is no any channel with "negative" signal in the same event.
	##
	## Algorithm: first, find approximate pedestal values (biased by "negative" signals) from untriggered and not retriggered data per chip, channel, SCA.
        ## Then, take a median pedestal inside every chip and, finally, find the minimum across the chips.
	## In the following "unbiased" pedestal calculation, only entries above (0.75 * this value) will be considered.
	## Note, for the chip with the minimal pedestals, this assumes that all its pedestals > 0.75 * (chip average).
	negative.adc.threshold <- 0.75 * min(merge(hits.local,
						   ev.local[,list(acq,bx,nbx)], # add nbx to hits.local to require nbx==1
						   by=c('acq','bx'))[trig==FALSE & nbx==1
						       ][, list(ped=as.double(median(adc))), by=list(chip,i,sca)
							 ][,list(ped=median(ped)),by=list(chip)]$ped)
	## cat('Negitive threshold:',negative.adc.threshold, '\n')

	## same per chip (name with .chip), events for one chip are sorted first in SCA
	ev.local.chip <- hits.local[,list(n.trig.chip=sum(trig),
					  n.neg.trig.chip=sum(adc<negative.adc.threshold)), keyby=list(acq,chip,sca,bx.cor,bx)]
	ev.local.chip[,bx.group.chip:=cumsum( c(0,diff(bx.cor)) > 4 ), by=list(acq,chip)]
	ev.local.chip[,`:=`(nbx.chip=.N, ibx.chip=1:.N), by=list(acq,chip,bx.group.chip)]

	## do the same, but taking into account only events with at least one trigger
	## (in SKIROC, if trigger arrives at the clock edge, both BX and BX+1 are written, but the latter without any triggers, if there is no retriggering)
	any.event.wo.trig <- any(ev.local.chip$n.trig.chip == 0)
	if (any.event.wo.trig) {
	    ev.local.trig <- ev.local.chip[n.trig.chip>0, list(acq,chip,sca,bx.cor,bx)] # remove events wo triggers, name with suffix .trig;
	    ## note: this automatically means per chip
	    ev.local.trig[,bx.group.trig:=cumsum( c(0,diff(bx.cor)) > 4 ), by=list(acq,chip)]  # find retriggers per chip AND only for events with triggers
	    ev.local.trig[,`:=`(nbx.trig=.N, ibx.trig=1:.N), by=list(acq,chip,bx.group.trig)]
	    ev.local.chip <- ev.local.trig[ev.local.chip]
	} else ev.local.chip[,`:=`(bx.group.trig=bx.group.chip, nbx.trig=nbx.chip, ibx.trig=ibx.chip)]

	ev.local.chip <- merge(ev.local, ev.local.chip, by=c('acq','bx'), all=TRUE)
	hits.local <- merge(hits.local, ev.local.chip, by=c('acq','chip','sca','bx.cor','bx'))

#        hits.local <- hits.local[ibx==1] # remove retriggers

	hits.local <- merge(hits.local, map[,list(chip,i,x,y)], by=c('chip','i'), all.x=TRUE)

	## subtract pedestals, find them from trig==FALSE & ibx==1 & (no channel in the same chip with ADC < negative.adc.threshold), otherwise set to NA
        ## change in Feb.2016: require nbx==1 as well to remove second pedestal peak due to retrigger, then ibx==1 automatically
	hits.local[,a := {
		       unbiased.pedestal <- trig==FALSE & n.neg.trig.chip == 0 & nbx==1
		       as.double(adc) - if(any(unbiased.pedestal)) median(adc[unbiased.pedestal]) else NA
		     }, by=list(chip,i,sca)]
	## hits.local[,l.a := {
	## 	       unbiased.pedestal <- trig==FALSE & n.neg.trig.chip == 0 & nbx==1
	## 	       as.double(l.adc) - if(any(unbiased.pedestal)) median(l.adc[unbiased.pedestal]) else NA
	## 	     }, by=list(chip,i,sca)]
	setkey(hits.local, acq, chip, bx.cor, i)
	setkey(ev.local.chip, acq, chip, bx.cor)
	setkey(ev.local, acq, bx)

	hits.local
##	ev <<- ev.local
##	ev.chip <<- ev.local.chip
    } else NULL
}

data.dir <- 'data'
files <- list.files(data.dir, pattern='\\.raw$')

mip <- data.table(file=files)[ ,{
    load.raw(paste0(data.dir, '/', file))[trig==TRUE & ibx==1 & !is.na(a)]
}, by=file]

write.table(mip, "MIP_all.dat", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)
