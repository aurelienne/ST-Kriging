library(gridExtra)
library(automap)
library(geoR)
library(sp)
library(spacetime)
library(rlist)
library(gstat)
library(Metrics)
library(ggplot2)

# Area
deltaY = 1
deltaX = 1
lonIni = 1
latIni = 1
lonMax = 27
latMax = 29
nX = 27
nY = 29

lats <- seq(latIni, latIni+nY*deltaY, deltaY)[1:nY]
lons <- seq(lonIni, lonIni+nX*deltaX, deltaX)[1:nX]
coords <- merge(lons, lats, all.x = TRUE, all.y = TRUE)
sp <- SpatialPoints(coords = coords, proj4string = CRS(as.character(NA)))

# Reading data
dt_list <- list()
bin_list <- list()
coords_list <- list()
prec_list <- list()
df_series <- data.frame(Obs=double())
dirname <- "/dados/radar/saoroque/cappi/cappi3km_tamanduatei/2019/caso_1003/"
for (filename in list.files(dirname)){
  print(filename)
  read.filename <- file(paste(dirname, filename, sep=""), "rb")
  bindata <- readBin(read.filename, "double",  n=nX*nY, size=4)
  bindata[bindata == -99] <- 0
  df <- data.frame(bindata)
  df_series <- rbind(df_series, df)
  coords_list <- append(coords_list, coords)
  prec_list <- append(prec_list, sum(bindata))
  dt <- substr(filename, 11, 22)
  print(dt)
  z <- as.POSIXct(dt,format="%Y%m%d%H%M")
  dt_list <- append(dt_list, z)
}
stf <- STFDF(sp, dt_list, data = df_series)
stplot(stf)
plot(dt_list, prec_list, type="l", col="blue", xlab="Horario", ylab="Prec. Tot. (mm/h)", 
     main="Precipitação Total na bacia do Tamanduateí", sub="Datas: 10/03/2019 - 11/03/2019")

# Sample Variogram
breaks <- seq(0,40,1)
vst <- variogramST(formula = bindata~1, locations = stf, data = stf, option = "bin", 
                   tlags = tlags, na.ommit=F, assumeRegular=T, cutoff = 15
#                   alpha=330, tol.hor=30, tol.ver=90
                   )
plot(vst)
plot(vst, map=F)
plot(vst[vst$timelag==0,], map=F)
boxplot(vst$gamma ~ vst$timelag, xlab="timelag", ylab="gamma")
plot(vst,wireframe=T)

# Variogram Model
sumMetric <- vgmST("sumMetric", space = vgm(psill=120,"Exp", range=14, nugget=0),
                   time = vgm(psill=50, "Exp", range=40, nugget=0), 
                   joint = vgm(psill=90, "Exp", range=40, nugget=0), 
                   stAni=14)
plot(vst, sumMetric, map=F, all=T)
sumMetric_Vgm <- fit.StVariogram(vst, sumMetric, stAni = 14)
plot(vst, sumMetric_Vgm, map=F, all=T)
attr(sumMetric_Vgm, "MSE")

# Cycle Prediction
breaks <- seq(0,40,1)
tlags <- seq(0,20,1)
N <- length(stf@time)-6
mae_pred_list <- vector("list", N)
mae_pred2_list <- vector("list", N)
rmse_pred_list <- vector("list", N)
rmse_pred2_list <- vector("list", N)
mae_persist_list <- vector("list", N)
mae_persist2_list <- vector("list", N)
rmse_persist_list <- vector("list", N)
rmse_persist2_list <- vector("list", N)
errormap <- data.frame()
mse <- 0
for (i in seq(1, N)) {
  print(i)
  f = i+4
  print(f)
  stf_sample <- stf[,i:f]
  stplot(stf_sample)
  vst <- variogramST(formula = bindata~1, locations = stf_sample, data = stf_sample, option = "bin", breaks = breaks, tlags = tlags, na.ommit=T, assumeRegular=T, cutoff=15)
  plot(vst)
  plot(vst, map=F)
  boxplot(vst$gamma ~ vst$timelag, xlab="timelag", ylab="gamma")
  plot(vst,wireframe=T)
  
  plot(vst, sumMetric, map=F, all=T)
  plot(vst, sumMetric, wireframe=T, all=T)

  # Fitting
  sumMetric_Vgm <- fit.StVariogram(vst, sumMetric, stAni = 14)
  plot(vst, sumMetric_Vgm, map=F, all=T)
  plot(vst, sumMetric_Vgm, wireframe=T, all=T)
  mse <- mse + attr(sumMetric_Vgm, "MSE")
  print(mse)

  # Predicting
  f1 = f+1
  f2 = f+2
  new_stf <- STF(sp, stf@time$timeIndex[f1:f2])
  pred <- krigeST(bindata~1, data=stf_sample, modelList=sumMetric_Vgm, newdata=new_stf)
  stplot(stf_sample)
  stplot(pred)
  
  # Calculating errors
  mae_pred <- mae(stf[,f1]@data$bindata, pred[,1]@data$var1.pred) 
  mae_pred2 <- mae(stf[,f2]@data$bindata, pred[,2]@data$var1.pred) 
  mae_persist <- mae(stf[,f1]@data$bindata, stf[,f]@data$bindata)  
  mae_persist2 <- mae(stf[,f2]@data$bindata, stf[,f]@data$bindata)  
  rmse_pred <- rmse(stf[,f1]@data$bindata, pred[,1]@data$var1.pred) 
  rmse_pred2 <- rmse(stf[,f2]@data$bindata, pred[,2]@data$var1.pred) 
  rmse_persist <- rmse(stf[,f1]@data$bindata, stf[,f]@data$bindata)  
  rmse_persist2 <- rmse(stf[,f2]@data$bindata, stf[,f]@data$bindata)  
  mae_pred_list[i] <- mae_pred
  mae_pred2_list[i] <- mae_pred2
  mae_persist_list[i] <- mae_persist
  mae_persist2_list[i] <- mae_persist2
  rmse_pred_list[i] <- rmse_pred
  rmse_pred2_list[i] <- rmse_pred2
  rmse_persist_list[i] <- rmse_persist
  rmse_persist2_list[i] <- rmse_persist2
  
  print(paste("MAE pred: ", mae_pred))
  print(paste("MAE PERSIST: ", mae_persist))
  print(paste("RMSE pred: ", rmse_pred))
  print(paste("RMSE PERSIST: ", rmse_persist))
  if (i==1){
    errormap_pred <- data.frame(rep(0, nX*nY))
    errormap_persist <- data.frame(rep(0, nX*nY))
  }
  errormap_pred <- errormap_pred + abs(pred[,1]@data-stf[,f1]@data)
  errormap_persist <- errormap_persist + abs(stf[,f]@data-stf[,f1]@data)
}

# Ploting Errors
mse <- mse/i
print(paste("MSE = ",mse))
  
errormap <- rbind(errormap_pred, errormap_persist)
errormap_pred_stf <- STFDF(sp, stf@time$timeIndex[f1:f2], data = errormap)
stplot(errormap_pred_stf, cuts=6, col.regions=c("light blue", "yellow","orange","red","brown"))

plot(seq_along(mae_pred_list), mae_pred_list, type="l", col="blue", ylim=c(4,12), 
     xlab="Previsões", ylab="MAE", main="Previsão 10 minutos")
lines(seq_along(mae_pred_list), mae_persist_list, col="red")
legend(0.5, 12, legend=c("Persistido", "Previsto"),
       col=c("red", "blue"), lty=1, cex=0.8)
plot(seq_along(mae_pred_list), mae_pred2_list, type="l", col="blue", ylim=c(4,12), 
     xlab="Previsões", ylab="MAE", main="Previsão 20 minutos")
lines(seq_along(mae_pred_list), mae_persist2_list, col="red")
legend(0.5, 12, legend=c("Persistido", "Previsto"),
       col=c("red", "blue"), lty=1, cex=0.8)

plot(seq_along(rmse_pred_list), rmse_pred_list, type="l", col="blue", ylim=c(6,14), 
     xlab="Previsões", ylab="RMSE", main="Previsão 10 minutos")
lines(seq_along(rmse_pred_list), rmse_persist_list, col="red")
legend(0.5, 14, legend=c("Persistido", "Previsto"),
       col=c("red", "blue"), lty=1, cex=0.8)
plot(seq_along(rmse_pred_list), rmse_pred2_list, type="l", col="blue", ylim=c(6,16), 
     xlab="Previsões", ylab="RMSE", main="Previsão 20 minutos")
lines(seq_along(rmse_pred_list), rmse_persist2_list, col="red")
legend(0.5, 16, legend=c("Persistido", "Previsto"),
       col=c("red", "blue"), lty=1, cex=0.8)
