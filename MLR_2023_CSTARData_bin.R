
library(dplyr)
library(MBC) #quantile mapping
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(grid)

print.stat.lm=function(x){
    mab=mean(abs(x$model[,1]-x$model[,2]))
    mapb=100*mean(abs(x$model[,1]-x$model[,2]))/mean(x$model[,2])
    r.squared=cor(x$model[,1],x$model[,2])^2
    if (sign(x$coefficients[2])>0){
    c(paste("y=", format(round(x$coefficients[1],0), nsmall=0), "+", format(round(x$coefficients[2],2), nsmall=2), "x; R-sq=", format(round(r.squared,2), nsmall=2), sep=""),
    paste("MAB=", format(round(mab,1), nsmall=1), ", MAPB=", format(round(mapb,1), nsmall=1), "%", sep=""))
    } else {
    c(paste("y=", format(round(x$coefficients[1],0), nsmall=0), "-", format(abs(round(x$coefficients[2],2)), nsmall=2), "x; R-sq=", format(round(r.squared,2), nsmall=2), sep=""),
    paste("MAB=", format(round(mab,1), nsmall=1), ", MAPB=", format(round(mapb,1), nsmall=1), "%", sep=""))
    }
}

project_home="C:\Users\maris\OneDrive - UCB-O365\Hannigan Research Lab\CSTAR\Shed Test 2025\calibrated fire activities"
setwd(project_home)
ft=readxl::read_excel("2023_fire_activities_calibrated_master.xlsx", sheet="Sheet1")


output=list()
id=c(1150,1152,1154,1155,1156,1157,1159)
tab=matrix(NA,ncol=4, nrow=24)
colnames(tab)=c("intercept","PM_raw","Temp","RH")
rownames(tab)=c(paste0("sensor",id,"_0_100"), paste0("sensor",id,"_100_1000"), paste0("sensor",id,"_over_1000"),"pooled_0_100","pooled_100_1000","pooled_over_1000")
tab2=matrix(NA,ncol=5, nrow=24)
colnames(tab2)=c("intercept","PM_raw","PMM_raw_sq","Temp","RH")
rownames(tab2)=c(paste0("sensor",id,"_0_100"), paste0("sensor",id,"_100_1000"), paste0("sensor",id,"_over_1000"),"pooled_0_100","pooled_100_1000","pooled_over_1000")


for (s in id){
    pw1=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet="0-100")
    pw2=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet="100-1000")
    pw3=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet=">1000")

    k=which(id==s)
    eq1=as.formula(paste0("POPS_PM25~PAM_PM25+TempC+RH"))
    eq2=as.formula(paste0("POPS_PM25~PAM_PM25+I(PAM_PM25^2)+TempC+RH"))
    
    #MLR (m1), log-MLR (m2) & RF (m3) by segments
    m1.pw1=lm(eq1, pw1)
    m1.pw2=lm(eq1, pw2)
    m2.pw2=lm(eq2, pw2)
    m1.pw3=lm(eq1, pw3)
    m2.pw3=lm(eq2, pw3)
    
    tab[k,]=coefficients(m1.pw1)
    tab[k+length(id),]=coefficients(m1.pw2)
    tab[k+2*length(id),]=coefficients(m1.pw3)
    tab2[k+length(id),]=coefficients(m2.pw2)
    tab2[k+2*length(id),]=coefficients(m2.pw3)
    
    
    png(paste("pm_sensor",s,"_pw1.png", sep=""), height=3, width=7.5, units="in", res=400)
    par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
    x=m1.pw1$model[,1]
    y=m1.pw1$model[,2]
    plot(x,y, xlim=c(-10,100),ylim=c(-10,100), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)
    
    x=m1.pw1$model[,1]
    y=m1.pw1$fitted.values
    plot(x,y, xlim=c(-10,100),ylim=c(-10,100), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)
    
    mtext("(a) raw data", side=3,line=0,outer=T, at=0.175)
    mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
    mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
    par(las=0)
    mtext(paste("PAM Sensor ", s, sep=""),side=2,line=1.3,outer=T,cex=0.9)
    dev.off()



    png(paste("pm_sensor",s,"_pw2.png", sep=""), height=3, width=7.5, units="in", res=400)
    par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
    x=m1.pw2$model[,1]
    y=m1.pw2$model[,2]
    plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)
    
    x=m1.pw2$model[,1]
    y=m1.pw2$fitted.values
    plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

    x=m1.pw2$model[,1]
    y=m2.pw2$fitted.values
    plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

    mtext("(a) raw data", side=3,line=0,outer=T, at=0.165)
    mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
    mtext("(c) quadratic regression", side=3,line=0,outer=T, at=0.835)
    mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
    par(las=0)
    mtext(paste("PAM Sensor ", s, sep=""),side=2,line=1.3,outer=T,cex=0.9)
    dev.off()



    png(paste("pm_sensor",s,"_pw3.png", sep=""), height=3, width=7.5, units="in", res=400)
    par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
    x=m1.pw3$model[,1]
    y=m1.pw3$model[,2]
    plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)
    
    x=m1.pw3$model[,1]
    y=m1.pw3$fitted.values
    plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

    x=m1.pw3$model[,1]
    y=m2.pw3$fitted.values
    plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
    mod=lm(y~x)
    abline(a=0, b=1, lty=3, lwd=2)
    abline(mod, col=2, lty=2, lwd=2)
    legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

    mtext("(a) raw data", side=3,line=0,outer=T, at=0.165)
    mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
    mtext("(c) quadratic regression", side=3,line=0,outer=T, at=0.835)
    mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
    par(las=0)
    mtext(paste("PAM Sensor ", s, sep=""),side=2,line=1.3,outer=T,cex=0.9)
    dev.off()

   
    
    if (s %in% unique(ft$sensor)){
        sensor=ft[ft$sensor==s,]
        sensor$PM_lin_QMcalibrated=NA
        sensor$PM_qua_QMcalibrated=NA
        
        #do quantile mapping based on m1
        qm=QDM(m.c=m1.pw1$fitted.values, o.c=m1.pw1$model[,2], m.p=c(sensor[sensor[,7] < 100,7])[[1]])
        sensor[sensor[,7] < 100,"PM_lin_QMcalibrated"]=qm$mhat.p
        
        qm=QDM(m.c=m1.pw2$fitted.values, o.c=m1.pw2$model[,2], m.p=c(sensor[sensor[,7] >= 100 & sensor[,7] < 1000,7])[[1]])
        sensor[sensor[,7] >= 100 & sensor[,7] < 1000,"PM_lin_QMcalibrated"]=qm$mhat.p
        
        qm=QDM(m.c=m1.pw3$fitted.values, o.c=m1.pw3$model[,2], m.p=c(sensor[sensor[,7] >= 1000,7])[[1]])
        sensor[sensor[,7] >= 1000,"PM_lin_QMcalibrated"]=qm$mhat.p
                
        #do quantile mapping based on m2
        qm=QDM(m.c=m1.pw1$fitted.values, o.c=m1.pw1$model[,2], m.p=c(sensor[sensor[,7] < 100,7])[[1]])
        sensor[sensor[,7] < 100,"PM_qua_QMcalibrated"]=qm$mhat.p
        
        qm=QDM(m.c=m2.pw2$fitted.values, o.c=m2.pw2$model[,2], m.p=c(sensor[sensor[,7] >= 100 & sensor[,7] < 1000,7])[[1]])
        sensor[sensor[,7] >= 100 & sensor[,7] < 1000,"PM_qua_QMcalibrated"]=qm$mhat.p
        
        qm=QDM(m.c=m2.pw3$fitted.values, o.c=m2.pw3$model[,2], m.p=c(sensor[sensor[,7] >= 1000,7])[[1]])
        sensor[sensor[,7] >= 1000,"PM_qua_QMcalibrated"]=qm$mhat.p
        
        #write.csv(sensor,paste("2023_fire_activities_calibrated_master_sensor", s, ".csv", sep=""), row.names=F)
        output[[k]]=sensor
    }
    print(paste("Sensor =", s))
}
output=bind_rows(output)




#sensor 1151, 1158
pw1=list()
pw2=list()
pw3=list()
for (s in id){
    k=which(id==s)
    pw1[[k]]=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet="0-100")
    pw2[[k]]=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet="100-1000")
    pw3[[k]]=readxl::read_excel(paste0("PAMs_Regime_Bins_Updated/PAM",s,"_Regime_Bins.xlsx"), sheet=">1000")
}
pw1=bind_rows(pw1)
pw2=bind_rows(pw2)
pw3=bind_rows(pw3)

#MLR (m1), log-MLR (m2) & RF (m3) by segments
m1.pw1=lm(eq1, pw1)
m1.pw2=lm(eq1, pw2)
m2.pw2=lm(eq2, pw2)
m1.pw3=lm(eq1, pw3)
m2.pw3=lm(eq2, pw3)

tab[22,]=coefficients(m1.pw1)
tab[23,]=coefficients(m1.pw2)
tab[24,]=coefficients(m1.pw3)
tab2[23,]=coefficients(m2.pw2)
tab2[24,]=coefficients(m2.pw3)






s="_pooled"
png(paste("pm_sensor",s,"_pw1.png", sep=""), height=3, width=7.5, units="in", res=400)
par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
x=m1.pw1$model[,1]
y=m1.pw1$model[,2]
plot(x,y, xlim=c(-10,100),ylim=c(-10,100), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

x=m1.pw1$model[,1]
y=m1.pw1$fitted.values
plot(x,y, xlim=c(-10,100),ylim=c(-10,100), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

mtext("(a) raw data", side=3,line=0,outer=T, at=0.175)
mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
par(las=0)
mtext(paste("Pooled PAM Sensors"),side=2,line=1.3,outer=T,cex=0.9)
dev.off()



png(paste("pm_sensor",s,"_pw2.png", sep=""), height=3, width=7.5, units="in", res=400)
par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
x=m1.pw2$model[,1]
y=m1.pw2$model[,2]
plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

x=m1.pw2$model[,1]
y=m1.pw2$fitted.values
plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

x=m1.pw2$model[,1]
y=m2.pw2$fitted.values
plot(x,y, xlim=c(100,1000),ylim=c(100,1000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

mtext("(a) raw data", side=3,line=0,outer=T, at=0.165)
mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
mtext("(c) quadratic regression", side=3,line=0,outer=T, at=0.835)
mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
par(las=0)
mtext(paste("Pooled PAM Sensors"),side=2,line=1.3,outer=T,cex=0.9)
dev.off()



png(paste("pm_sensor",s,"_pw3.png", sep=""), height=3, width=7.5, units="in", res=400)
par(mar=c(2, 1.8, 0.5, 1.5), las=1, oma=c(2,2.5,1.5,0), mfrow=c(1,3))
x=m1.pw3$model[,1]
y=m1.pw3$model[,2]
plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

x=m1.pw3$model[,1]
y=m1.pw3$fitted.values
plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

x=m1.pw3$model[,1]
y=m2.pw3$fitted.values
plot(x,y, xlim=c(1000,5000),ylim=c(1000,5000), xlab="", ylab="", main="", col=rgb(0.5,0.5,0.5,0.5), pch=20)
mod=lm(y~x)
abline(a=0, b=1, lty=3, lwd=2)
abline(mod, col=2, lty=2, lwd=2)
legend("bottomright", print.stat.lm(mod), text.col=2, lty=0, bty='n', ncol=1)

mtext("(a) raw data", side=3,line=0,outer=T, at=0.165)
mtext("(b) linear regression", side=3,line=0,outer=T, at=0.5)
mtext("(c) quadratic regression", side=3,line=0,outer=T, at=0.835)
mtext(expression(paste("POPS Referenced ", PM[2.5], " [", mu, "g/", m^3, "]", sep="")), side=1,line=1,outer=T,cex=0.9)
par(las=0)
mtext(paste("Pooled PAM Sensors"),side=2,line=1.3,outer=T,cex=0.9)
dev.off()



additions=list()
id=c(1151,1158)
for (s in id){
    k=which(id==s)
    sensor=ft[ft$sensor==s,]
    sensor$PM_lin_QMcalibrated=NA
    sensor$PM_qua_QMcalibrated=NA
    
    #do quantile mapping based on m1
    qm=QDM(m.c=m1.pw1$fitted.values, o.c=m1.pw1$model[,2], m.p=c(sensor[sensor[,7] < 100,7])[[1]])
    sensor[sensor[,7] < 100,"PM_lin_QMcalibrated"]=qm$mhat.p
    
    qm=QDM(m.c=m1.pw2$fitted.values, o.c=m1.pw2$model[,2], m.p=c(sensor[sensor[,7] >= 100 & sensor[,7] < 1000,7])[[1]])
    sensor[sensor[,7] >= 100 & sensor[,7] < 1000,"PM_lin_QMcalibrated"]=qm$mhat.p
    
    qm=QDM(m.c=m1.pw3$fitted.values, o.c=m1.pw3$model[,2], m.p=c(sensor[sensor[,7] >= 1000,7])[[1]])
    sensor[sensor[,7] >= 1000,"PM_lin_QMcalibrated"]=qm$mhat.p
            
    #do quantile mapping based on m2
    qm=QDM(m.c=m1.pw1$fitted.values, o.c=m1.pw1$model[,2], m.p=c(sensor[sensor[,7] < 100,7])[[1]])
    sensor[sensor[,7] < 100,"PM_qua_QMcalibrated"]=qm$mhat.p
    
    qm=QDM(m.c=m2.pw2$fitted.values, o.c=m2.pw2$model[,2], m.p=c(sensor[sensor[,7] >= 100 & sensor[,7] < 1000,7])[[1]])
    sensor[sensor[,7] >= 100 & sensor[,7] < 1000,"PM_qua_QMcalibrated"]=qm$mhat.p
    
    qm=QDM(m.c=m2.pw3$fitted.values, o.c=m2.pw3$model[,2], m.p=c(sensor[sensor[,7] >= 1000,7])[[1]])
    sensor[sensor[,7] >= 1000,"PM_qua_QMcalibrated"]=qm$mhat.p
    additions[[k]]=sensor
}

additions=bind_rows(additions)
output=rbind(output,additions)
write.csv(output,"2023_fire_activities_calibrated_master_QM.csv", row.names=F)

tp=tableGrob(round(tab,2))
my_title = textGrob("linear regression coefficients by sensors and segments", gp=gpar(fontsize=14, fontface="bold", color="darkblue"))
png("reg_table.png", height=7.5, width=6, units="in", res=350)
grid.arrange(tp, top= my_title)
dev.off()

tp=tableGrob(round(tab2,6))
my_title = textGrob("quadratic regression coefficients by sensors and segments", gp=gpar(fontsize=14, fontface="bold", color="darkblue"))
png("quareg_table.png", height=7.5, width=8, units="in", res=350)
grid.arrange(tp, top= my_title)
dev.off()


