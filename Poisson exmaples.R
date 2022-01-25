require(ggplot2)

zf <- read.table("C:/Users/u0111580/Desktop/LOSdata.csv", header = TRUE, sep = ";")
# modify the above line as needed
Pop  = data.frame(LOS=zf$LOS,Age=zf$Age,Sex=zf$Sexe,DRG=zf$DRG,MDC=zf$MDC,Nbdg=zf$NbDiag,Nbtt=zf$NbAct,CW=zf$CWeff)
table(Pop$MDC) # These are the DRG groups

