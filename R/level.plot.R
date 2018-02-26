`level.plot` <-
function(data.in, XY, color.gradient='red', cex=1, level.range=c(min(data.in),max(data.in)), show.scale=TRUE, title="level plot", SRC=FALSE, save.file="no", ImageSize="small", AddPresAbs=NULL, PresAbsSymbol=c(cex*0.8,16,4),...){
  
  extra.args <- list(...)
  if(!is.null(extra.args$multiple.plot)){
    multiple.plot <- extra.args$multiple.plot
  } else {
    multiple.plot <- FALSE
  }
  
    
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n") 
    if(ncol(XY)!=2) stop("\n wrong coordinates given in 'XY' : there should be two columns \n")
    if(nrow(XY)!=length(data.in)) stop("\n data and coordinates should be of the same length \n")

#     if() multiple.plot <- TRUE  else multiple.plot <- FALSE


    if(SRC){
        if(length(unique(data.in))>4){
            cat("\n not possible to render SRC plot -> more than four different values in data ")
            SRC <- F
        } else {    
            SRCvalues <- sort(unique(data.in))
			      col_id <- data.in + 3
            color.system <- c("red", "lightgreen", "grey", "darkgreen")
            title <- paste("SRC plot ", title, sep="")
        }
    } else
	{    if(color.gradient=='grey') {
#           color.system <- c()
#           for(i in seq(93,10,length.out=100)) color.system <- c(color.system, gray(i/100))
#           color.system <- c(gray(0.93), color.system, gray(0))
	        color.system <- gray(seq(0.95,0, length.out=102))
            
        }
        if(color.gradient=='blue') {
          color.system <- c('grey88',
          rainbow(45, start=0.5, end=0.65),                       
          rainbow(10, start=0.65, end=0.7),
          rainbow(45, start=0.7, end=0.85),
          'red')
        }
        if(color.gradient=='red') {    
          color.system <- c(
          'grey88',
          c(rep(c(colors()[c(417,417,515)]), each=5),
          rev(rainbow(55, start=0.13, end=0.23 )),
          rev(rainbow(50, start=0.08, end=0.13 )[seq(1,50,length.out=15)]),
          rev(rainbow(50, end=0.08)[seq(1,50,length.out=15)])), 
          'brown2')
        } 
    
	
      #if range wanted is broader than possible, set to actual range limits
      #if(level.range[1]<min(data.in)) level.range[1] <- min(data.in)  
      #if(level.range[2]>max(data.in)) level.range[2] <- max(data.in)  
    
      #determine the color code to assess to each value
      
      # define a vector for make a correspundance between values and colors
      val_seq <- c(seq(level.range[1], level.range[2], length.out=101),Inf)  
      col_id <- sapply(data.in, function(x){return(which(x<=val_seq)[1])})
       
       
#       g <- gg <- data.in
#       gg[gg <= level.range[1]]  <- level.range[1]
#       gg[gg >= level.range[2]] <- level.range[2]
#       gg <- gg-min(g) 
#       gg <- gg/max(gg)*100 + 1
#    
#       #over and under-ranged values set to limits of color range
#       gg[g < level.range[1]] <- 1
#       gg[g > level.range[2]] <- 102
    }
    
    # define plotting symbols for presence and absences if required by user
    if(!is.null(AddPresAbs)){
        cex2 <- PresAbsSymbol[1]
        pchPres <- PresAbsSymbol[2]
        pchAbs <- PresAbsSymbol[3]
    }
    
    #define image size for JPEG and TIFF
    if(ImageSize=="small") {SizeInPix <- 480; FontSize=12} else if(ImageSize=="standard") {SizeInPix <- 1000; FontSize=22} else if(ImageSize=="large") {SizeInPix <- 2000; FontSize=44}
    
    if(save.file == "pdf") pdf(paste(title, ".pdf", sep=""))
    if(save.file == "jpeg") jpeg(paste(title, ".jpeg", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize, quality=85)
    if(save.file == "tiff") tiff(paste(title, ".tiff", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize)
    if(save.file == "postscript") postscript(paste(title, ".eps", sep=""))
  
    if(show.scale){
        
        if(!multiple.plot) layout(matrix(c(1,2),nrow=1), widths=c(5,1), heights=c(1,1))
        plot(XY[,2]~XY[,1], col=color.system[col_id], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
        #Add Presence and Absence locations if requested by user:
        if(!is.null(AddPresAbs)){points(AddPresAbs[AddPresAbs[,3]==1,1:2], col="black", pch=pchPres, cex=cex2); points(AddPresAbs[AddPresAbs[,3]==0,1:2], col="black", pch=pchAbs, cex=cex2)}
     
        par(mar=c(0.1,0.1,0.1,0.1))
        plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE) 
        polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
        
        if(SRC){ legend(0,0.8,legend=list(' (1) new', ' (0) stable', '(-1) kept', '(-2) lost'),cex=1, fill=rev(color.system),bty='n') 
         } else {
          if(level.range[1] == min(data.in)) lmin <- round(level.range[1], digits=2) else lmin <- paste(round(level.range[1], digits=2), " or lower", sep="")
          if(level.range[2] == max(data.in)) lmax <- round(level.range[2], digits=2) else {lmax <- paste(round(level.range[2], digits=2), " or over", sep="") }

          if(!multiple.plot){
              legend(0.2,0.92,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
              '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin),cex=1, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')
          
          } else legend(0.2,1.05,legend=list(lmax,'','','','',round((3*level.range[2]+level.range[1])/4, digits=2),'','','','',round(sum(level.range)/2, digits=2),
          '','','','',round((level.range[2]+3*level.range[1])/4, digits=2),'','','','',lmin), cex=cex, fill=rev(color.system[c(1,seq(2,101,length.out=19),102)]),bty='n')
        
        }
    }
     else{
          plot(XY[,2]~XY[,1], col=color.system[col_id], cex=cex, pch=19, xlab='', ylab='', xaxt='n', yaxt='n', main=title)
          #Add Presence and Absence locations if requested by user:
          if(!is.null(AddPresAbs)){points(AddPresAbs[AddPresAbs[,3]==1,1:2], col="black", pch=pchPres, cex=cex2); points(AddPresAbs[AddPresAbs[,3]==0,1:2], col="black", pch=pchAbs, cex=cex2)}
     }
    
    if(!is.null(AddPresAbs)) rm(cex2,pchPres,pchAbs)
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
     
}
