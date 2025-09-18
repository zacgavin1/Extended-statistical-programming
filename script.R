x <- c("10","2","7","89","43","1") ## example vector
ii <- which(nchar(x)>1) # where the double digits are
xs <- rep("", length(ii)+length(x))
iis <- ii+ 1:length(ii) # each element gets shifted one more over than the previous one
xs[iis]<- substr(x[ii],2,2)
xs[-iis]<- substr(x,1,1)
xs