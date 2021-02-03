library(VSURF) # version 1.0.3 was used in publication but results remain stable with newer versions (up to version 1.1. tested)
library(data.table)
installed.packages()[which(installed.packages()[,"Package"]=="VSURF"),"Version"]

allfeat <- readRDS("dat/features_prediction_b150.RDS")
data <- fread("dat/bacterial_load_persistence_pop-growth.csv")
meta <- fread("dat/meta_pwy.tbl")

idx <- match(data$strain, rownames(allfeat))
selfeat <- allfeat[idx,]
lmdat <- selfeat[,apply(selfeat, 2, function(x){(length(unique(x))>1)})] # remove constant cols

# mean population growth
target1 <- data$mean_pop
rf1 <- VSURF(x=lmdat, y=target1, parallel=T)
rf1.imp <- data.table(nr=rf1$imp.mean.dec.ind, id=colnames(lmdat)[rf1$imp.mean.dec.ind], imp.mean=rf1$imp.mean.dec, imp.sd=rf1$imp.sd.dec)
rf1.imp[,name:=meta$name[match(id, meta$id)]]
rf1.imp[nr %in% rf1$varselect.interp] # highest importance: PWY-7515 (trans-3-hydroxy-L-proline degradation)

# bacterial load
target2 <- log10(data$mean_cfu_load)
rf2  <- VSURF(x=lmdat, y=target2, parallel=T)
rf2.imp <- data.table(nr=rf2$imp.mean.dec.ind, id=colnames(lmdat)[rf2$imp.mean.dec.ind], imp.mean=rf2$imp.mean.dec, imp.sd=rf2$imp.sd.dec)
rf2.imp[,name:=meta$name[match(id, meta$id)]]
rf2.imp[nr %in% rf2$varselect.interp] # highest importance: PWY-6389 (pyruvate fermentation to (S)-acetoin)
