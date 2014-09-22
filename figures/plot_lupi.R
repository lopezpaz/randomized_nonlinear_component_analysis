library(tikzDevice)

lupi <- rbind(
c(55,62.45),
c(56.85,64.55),
c(56.25,64.55),
c(53.9,61.65),
c(53.05,62.15),
c(54,62.45),
c(55.8,60.6),
c(56.65,65.6),
c(54.85,61.55),
c(55.7,61.6),
c(54.5,61.75),
c(56.6,62.9),
c(54.2,61.8),
c(54.8,64.7))

tikz("lupi.tex", height=2, width=5, standAlone=TRUE)
par(mar=c(2,4,0.5,0.5),cex=1.3)

plot(lupi[,1], t="o", pch=0, ylim=c(50,70), ylab="Classification Acc.",lwd=2)
points(lupi[,2], t="o", pch=1,lwd=2)

legend("topright", c("SURF", "RCCA"), lty=rep(1,2), pch=0:1,lwd=2,bg="white")
dev.off()

system("pdflatex lupi.tex; evince lupi.pdf")
