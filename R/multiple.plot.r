`multiple.plot` <-
function(Data, coor, color.gradient='red', plots.per.window=9, cex=1, save.file="no", name="multiple plot", ImageSize="small", AddPresAbs=NULL, PresAbsSymbol=c(cex*0.8,16,4)){

    if(nrow(coor) != nrow(Data)) stop("Uncorrect mapping coordinates : coor and Data are not of the same length")
    if(color.gradient!='grey' && color.gradient!='red' && color.gradient!='blue') stop("\n color.gradient should be one of 'grey', 'red' or 'blue' \n")

#     assign("multiple", 564, pos=1) #to communicate to level.plot that a multiple plot is wanted (564 is just random)
    multiple <- 564
    
    #function plotting color boxes
    pbox <- function(co){ 
        plot(x=c(-1,1), y=c(0,1), xlim=c(0,1), ylim=c(0,1), type="n", axes=FALSE) 
        polygon(x=c(-2,-2,2,2), y=c(-2,2,2,-2), col=co, border=NA) 
    }

    #Take off NA data
    Data <- Data[ ,!is.na(Data[1,])]
    

    #calculating the number of windows to open    
    NbPlots <- ncol(Data)
    NbWindows <- ceiling(NbPlots/plots.per.window)
    if(NbWindows==1) plots.per.window <- NbPlots
    
    #define image size for JPEG and TIFF
    if(ImageSize=="small") {SizeInPix <- 480; FontSize=12} else if(ImageSize=="standard") {SizeInPix <- 1000; FontSize=22} else if(ImageSize=="large") {SizeInPix <- 2000; FontSize=44}

    if(save.file=="pdf") pdf(paste(name, ".pdf", sep=""))
    if(save.file=="jpeg") jpeg(paste(name, ".jpeg", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize, quality=85)
    if(save.file=="tiff") tiff(paste(name, ".tiff", sep=""), width=SizeInPix, height=SizeInPix, pointsize=FontSize)
    if(save.file=="postscript") postscript(paste(name, ".eps", sep=""))    
    
    for(W in 1:NbWindows){
        if(save.file=="no") dev.new()
        
        Wstart <- (W-1)*plots.per.window + 1
        if(W*plots.per.window > NbPlots) Wfinal <- NbPlots  else Wfinal <- W*plots.per.window
        DataW <- as.data.frame(Data[ ,Wstart:Wfinal])
        colnames(DataW) <- colnames(Data)[Wstart:Wfinal]
        
        #determine the organisation of the plots on the window
        W.width <- ceiling(sqrt(plots.per.window))
        W.height <- ceiling(plots.per.window/W.width)

        #create object for scaling the legend
        legendcex <- 0.64+1/exp(W.height)
#         assign("legendcex", 0.64+1/exp(W.height), pos=1)
        
        #matrix of indexes for ordering the layout
        mat <- c(1,2)
        for(i in 1:(W.width-1))  mat <- c(mat, mat[1:2] + 4*i)
        mat <- rbind(mat, mat+2)
        for(i in 1:(W.height-1))  mat <- rbind(mat, mat[1:2,] + W.width*4*i)  
        
        layout(mat, widths=rep(c(1,0.3), W.width), heights=rep(c(0.2,1), W.height))
        
        par(mar = c(0.1,0.1,0.1,0.1))
        for(i in 1:(Wfinal-Wstart+1)){
             pbox("grey98")
             text(x=0.5, y=0.8, pos=1, cex=1.6, labels=colnames(DataW)[i], col="#4c57eb")
             pbox("grey98")
             level.plot(DataW[,i], XY=coor, color.gradient=color.gradient, cex=cex, title="", AddPresAbs=AddPresAbs, PresAbsSymbol=PresAbsSymbol, multiple.plot=TRUE) 
        }
        
        #fill gaps by grey boxes
        if(W.width*W.height-plots.per.window != 0) for(i in 1:((W.width*W.height-plots.per.window)*4)) pbox("grey98")                    
            
    } #W loop   
    
#     rm('legendcex', 'multiple', pos=1)
    if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
}
