library(tikzDevice)

options( tikzLatexPackages = c(
  getOption( "tikzLatexPackages" ),
  "\\usepackage{bm,times}"
))

my_mar <- c(3.3,3.3,1,1)
my_mgp <- c(2.1,1,0)
my_cex <- 1.5

r <- rbind(c(1e-04,5.792990e+04),
           c(1e-03,6.431905e+03),
           c(1e-02,6.014578e+02),
           c(1e-01,6.245001e+01),
           c(1e+00,6.651322e+00),
           c(1e+01,7.249647e-01),
           c(1e+02,8.783448e-02),
           c(1e+03,1.341658e-02))

tikz("norm_cca_g.tex",standAlone=TRUE,width=3,height=3)
par(mar=my_mar, mgp=my_mgp,cex=my_cex) 
plot(log(r),t="o",xlab="$\\log(\\gamma)$",ylab="$\\log(\\|\\hat{\\bm R}^{-1}\\hat{\\bm L}-\\bm R^{-1}\\bm L\\|$)",lwd=3) 
dev.off()

system("pdflatex norm_cca_g.tex")

r <- rbind(c( 100, 60.574517, 19657.962),
           c( 250, 38.365135, 11739.667),
           c( 500, 27.343299,  8340.034),
           c(1000, 18.891724,  5904.768),
           c(1500, 14.984350,  4727.440),
           c(2000, 14.091385,  4141.354),
           c(2500, 12.474679,  3709.368),
           c(3000, 11.025231,  3372.574),
           c(4000,  9.240461,  2910.546),
           c(5000,  8.677846,  2603.916))

tikz("norm_cca_m.tex",standAlone=TRUE,width=3,height=3);
par(mar=my_mar, mgp=my_mgp,cex=my_cex) 
plot(r[,1],r[,3],t="o",xlab="$m$",ylab="$\\|\\hat{\\bm R}^{-1}\\hat{\\bm L}-\\bm R^{-1}\\bm L\\|$",lwd=3) 
dev.off()

system("pdflatex norm_cca_m.tex")

r <- rbind(c( 100,  1.961450,   393.3971),
           c( 250,  4.850861,  1295.4472),
           c( 500,  9.189258,  2782.4135),
           c(1000, 21.171350,  5834.5632),
           c(1500, 32.406314,  9154.8163),
           c(2000, 34.329282, 12056.4939),
           c(2500, 47.685869, 15105.3356),
           c(3000, 63.446192, 17977.8319),
           c(4000, 71.370321, 24444.7370),
           c(5000, 92.434716, 30033.6785))
tikz("norm_cca_n.tex",standAlone=TRUE,width=3,height=3);
par(mar=my_mar, mgp=my_mgp,cex=my_cex) 
plot(r[,1],r[,3],t="o",xlab="$n$",ylab="$\\|\\hat{\\bm R}^{-1}\\hat{\\bm L}-\\bm R^{-1}\\bm L\\|$",lwd=3)
dev.off()

system("pdflatex norm_cca_n.tex")

r <- rbind(c( 100, 60.574517, 19657.962),
           c( 250, 38.365135, 11739.667),
           c( 500, 27.343299,  8340.034),
           c(1000, 18.891724,  5904.768),
           c(1500, 14.984350,  4727.440),
           c(2000, 13.591385,  4141.354),
           c(2500, 12.474679,  3709.368),
           c(3000, 11.025231,  3372.574),
           c(4000,  9.240461,  2910.546),
           c(5000,  8.677846,  2603.916))
fff <- 
c(60.623741,
  38.305267,
  27.056785,
  19.102907,
  15.579209,
  13.478666,
  12.045186,
  10.987035,
   9.501727,
   8.488104)

tikz("norm_pca_m.tex",standAlone=TRUE,width=3,height=3);
par(mar=my_mar, mgp=my_mgp,cex=my_cex) 
plot(r[,1],r[,2],t="o",xlab="$m$",ylab="$\\|\\hat{\\bm K}-\\bm K\\|$",lwd=3) 
lines(r[,1],fff,col="red",lwd=3)
dev.off()

system("pdflatex norm_pca_m.tex")

r <- rbind(c( 100, 1.874253),
           c( 250, 4.923609),
           c( 500, 9.114774),
           c(1000,18.907895),
           c(1500,27.906247),
           c(2000,39.030703),
           c(2500,45.969439),
           c(3000,58.737486),
           c(4000,75.526354),
           c(5000,95.037361))

tikz("norm_pca_n.tex",standAlone=TRUE,width=3,height=3);
par(mar=my_mar, mgp=my_mgp,cex=my_cex) 
plot(r[,1],r[,2],t="o",xlab="$n$",ylab="$\\|\\hat{\\bm K}-\\bm K\\|$",lwd=3)
dev.off()

system("pdflatex norm_pca_n.tex")
