"trad.scatter.plot" <-
function(x,y,add=FALSE,fc.lines=log2(c(2,4,6,8)),draw.fc.lines=TRUE,draw.fc.line.labels=TRUE,fc.line.col="lightgrey",pch=20,...) {
   if(!add) {
     mx <- max(x,y);     
     mn <- min(x,y);
     plot(x,y,xlim=range(mn,mx),ylim=range(mn,mx),pch=pch,...);
     if(draw.fc.lines) {
       for( lne in fc.lines) {
         xc <- c(mn,mx);
         yc <- c(mn+lne,mx+lne);
         lines(xc,yc,col=fc.line.col); 
         if(draw.fc.line.labels) {text(mn+0.25,lne+mn,2^lne,col=fc.line.col); }
         xc <- c(mn, mx);
         yc <- c(mn-lne,mx-lne);
         lines(xc,yc,col=fc.line.col); 
          if(draw.fc.line.labels) {text(lne+mn,mn+0.25,2^lne,col=fc.line.col); }
       }
     }
   }
   else {
     points(x,y,pch=pch,...); 
   }
}
