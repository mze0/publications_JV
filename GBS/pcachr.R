pcachr <- function(dataPCAChr,Conf,Alpha,ind){

Lab <- 2.4
Main <- 4
Ax <- 2.5
	
	
ResPCA <- PCA(dataPCAChr, scale = T,graph=F)

vPC1 <- ResPCA$var$coord[,1]
vPC2 <- ResPCA$var$coord[,2]
vlabs <- rownames(ResPCA$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))

Names<-rownames(ResPCA$var$coord)
Names[grep("MDA",Names)]<-"MDA"
Names[grep("Pic",Names)]<-"Pic"
Names[grep("MC",Names)]<-"MC"

LabPCA <- cbind.data.frame(Names,ResPCA$var$coord[,1],ResPCA$var$coord[,2])
colnames(LabPCA)<-c("WGA","P1","P2")
LabPCA[nrow(LabPCA),"WGA"] <- "MDA"

plot(vPCs$vPC1, vPCs$vPC2, pch = 18, col = "black",main="PCA logR-values", xlab="PC1", ylab="PC2",cex=1,cex.lab=Lab,cex.main=3,cex.axis=Ax,frame=F,yaxt="n",xaxt="n")
axis(1, tck = 1, col = "#BABABA70", lty = 2,cex.axis=Ax)
axis(2, tck = 1, col = "#BABABA70", lty = 2,cex.axis=Ax)
abline(h=0, v=0, col = "grey")


ElipMDA<-dataEllipse(LabPCA[LabPCA[,"WGA"]=="MDA","P1"],LabPCA[LabPCA[,"WGA"]=="MDA","P2"],levels=Conf,col="blue",robust=T,fill = T, fill.alpha = Alpha,add=T,center.pch=F,lwd=0)
ElipPic<-dataEllipse(LabPCA[LabPCA[,"WGA"]=="Pic","P1"],LabPCA[LabPCA[,"WGA"]=="Pic","P2"],levels=Conf,col="red",robust=T,fill = T, fill.alpha = Alpha,add=T,center.pch=F,lwd=0)
ElipMC<-dataEllipse(LabPCA[LabPCA[,"WGA"]=="MC","P1"],LabPCA[LabPCA[,"WGA"]=="MC","P2"],levels=Conf,col="green",robust=T,fill = T, fill.alpha = Alpha,add=T,center.pch=F,lwd=0)

points(LabPCA[nrow(LabPCA),"P1"],LabPCA[nrow(LabPCA),"P2"],col="black",pch=10, cex=5)

Points<-expand.grid(LabPCA[ind,"P1"],LabPCA[ind,"P2"])
PolyPoints<-round(ElipMDA)
out = pnt.in.poly(Points,PolyPoints)
if(out$pip==1){
str<-paste(blast[j],"chr_",chrnr,"inside_blue",out$pip,sep=" ")

write(str,paste(outPath,Family,"_ChrSpecOutPCA.txt"),append=T,sep="\t")
print(str)
}

}