"trad.scatter.plot" <-
function(x,y,add=F,fc.lines=c(1,2,3,4),draw.fc.lines=T,draw.fc.line.labels=T,fc.line.col="lightgrey",pch=20,...) {
   if(!add) {
     mx <- max(x,y);     
     mn <- min(x,y);
     plot(x,y,xlim=range(mn,mx),ylim=range(mn,mx),pch=pch,...);
     if(draw.fc.lines) {
       for( lne in fc.lines) {
         xc <- c(mn,mx-lne);
         yc <- c(lne + mn, mx);
         lines(xc,yc,col=fc.line.col); 
         if(draw.fc.line.labels) {text(mn+0.25,lne+mn,2^lne,col=fc.line.col); }
         xc <- c(lne + mn, mx);
         yc <- c(mn,mx-lne);
         lines(xc,yc,col=fc.line.col); 
          if(draw.fc.line.labels) {text(lne+mn,mn+0.25,2^lne,col=fc.line.col); }
       }
     }
   }
   else {
     points(x,y,pch=pch,...); 
   }
}
