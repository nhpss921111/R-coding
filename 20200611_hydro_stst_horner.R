# 水文統計
# 雨量資料處理
# 撰寫時間：2020.05.29
# BY連育成

# 待改進地方
# 1. 起始年和結束年用list表示
# 2. 1~24延時改成用list表示
# 3. 重現期用list表示
# 4. 最適模型判斷：ks檢定的p-value，要再加入AIC或BIC指標

rm(list = ls()) # 刪除資料

library(fitdistrplus) #
library(FAdist) 
library(openxlsx) # 可以讀取資料量大的檔案
library(dplyr)
library(grDevices)
library(stats) # 機率分布套件
library(actuar) # 機率分布套件
library(gumbel)
library(EnvStats) # 估計gumble的參數
library(goft)
library(ggplot2) #畫圖用

# ---------------------- 參數設定區 --------------




# ----------------------------------------


setwd("E:/R_reading")  # data資料夾路徑
# 全部讀取
rowdata <- read.xlsx(file.path(getwd(), "daton_1987-2019.xlsx"), 
                     sheet= "rainfall" ,skipEmptyRows = TRUE,skipEmptyCols = TRUE,
                     na.strings = "NA") 
# 截取需要的資料
data <- rowdata[ ,7:30]

# 檢查資料
sum(is.na(data))
sum(data=="null")
sum(data==-999997)
data[is.na(data)] <- 0 #把NA值換成0
data[data<0] <- 0  #把負數換成0

# 閏年的判斷函數
leap_year <- function(date) {
  if (is.numeric(date)) {
    year <- date
  } else {
    year <- year(date)
  }
  (year %% 4 == 0) & ((year %% 100 != 0) | (year %% 400 == 0))
}


#  資料分割

yi <- 1987 #請輸入:起始年分
yend <- 2019 #請輸入:結束年分
y <- c(yi:yend)
daycut <- c()
d <- 1 # 第一筆雨量資料位置
for(i in c(yi:yend)){
  if(leap_year(i)){
    delta_d <- 365 # 閏年+365天
    d <- d + delta_d
  }else{
    delta_d <- 364 # 平年+364天
    d <- d + delta_d}
  
  #print(delta_d)
  daycut <- append(daycut,d) #每一年的最後一天
  d <- d + 1
  
}

# 分成每年一組數據，計算各延時下最大降雨深度
# 儲存在m矩陣裡

m <- matrix(nrow=33, ncol=24)
a <- 1
j <- 0
for (i in 1:length(daycut)){  #以每年的最後一天當作切點，把每年的資料分割出來
  y <- slice(data, a:daycut[i])
  x <- c()
  for (j in 1:length(y[,1])){ # j 為當年天數:365 or 366
    x <- append(x, y[j,]) # 已加入第 j 天資料到 x，x為每年的數據
    j <- j+1
  }
  yi <- c()
  yi <- as.matrix(x)
  yi <- as.numeric(yi) # 每年有 (一年天數*24小時)筆資料
  print(paste0("已分割出第",i,"年資料"))
  print(paste0("當年資料數量",length(yi),"筆")) 
  for(h in c(0:23)){  # 延時1小時到延時24小時
    print(paste0("延時",h+1,"小時"))
    f <- c()
    for(n in 1:length(yi)){  # 位置編號 1:(365*24)
      if (h==0){
        f <- append(f,yi[n]) # 連續一小時
      }else{
        yi[n][is.na(yi[n])] <- 0 # 如果是NA，給值0
        yi[n+h][is.na(yi[n+h])] <- 0 # 如果是NA，給值0
        f <- append(f, sum(yi[n:(n+h)]))# 第n個位置往後h延時，加總
        ##可惡的 "(" n+h ")" ，害我搞了三天!!!!!!!!!!!!!!
      }
    }
    m[i,h+1] <- max(f,na.rm=TRUE)  #取每個延時的最大值，且捨棄NA值
  }
  print(paste0("第",i,"年延時",h+1,"小時完成"))
  a <- daycut[i] + 1 # 為下一年的起始日 +1
}


m <- m[rowSums(m==0)==0,] # 把有0的那年資料刪除


# ----------------------------------------------------------------------

print("參數估計")
candidate <- c("norm","lnorm","gumbel","gamma3","lgamma3")

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))


dist.table <- matrix(nrow=length(candidate)+1,ncol=24)
rownames(dist.table) <- c("norm","lnorm","gumbel","gamma3","lgamma3","good dist")
#dist.table[(length(candidate)+1),] <- candidate[1] #假設最適分布為第一個候選分布
for(i in 1:24){
  var <- m[,i]
  print(paste0("延時",i,"小時"))
  # By Maximun Likelihood Estimate Method
  for(dist in c(1:length(candidate))){
    # -------------------------- parameter estimate -----------------------
    print(candidate[dist])
    dist.char <- c(candidate[dist],
                   paste0("d", candidate[dist]), 
                   paste0("p", candidate[dist]),
                   paste0("q", candidate[dist]))
    
    #md <- fitdistr(x, distribution, start = list(parameter1 = 1, parameter2 = 1))
    if(candidate[dist] == "norm"){
      md <- fitdist(var, dist = dist.char[1])
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "lnorm"){
      md <- fitdist(var, dist = dist.char[1], start = list(meanlog=1, sdlog=1))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "gumbel"){
      md <- eevd(var,method = "mle") 
      par1 <- md$parameters[1] #fitting參數1
      par2 <- md$parameters[2] #參數2
      print(c(par1, par2))
    }else if(candidate[dist] == "gamma3"){
      md <- fitdist(var, dist = dist.char[1], start=list(shape=1,scale=1,thres=1),lower=c(0,0))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      par3 <- md$estimate[3] #參數3
      print(c(par1, par2, par3))
      plot(md)
    }else if(candidate[dist] == "lgamma3"){
      md <- fitdist(var, dist = dist.char[1],start=list(shape=1,scale=1,thres=1),lower=c(0,0))
      par1 <- md$estimate[1] #fitting參數1
      par2 <- md$estimate[2] #參數2
      par3 <- md$estimate[3] #參數3
      print(c(par1, par2, par3))
      plot(md)
    }else{break}
    
    
    
    
    # ------------------------------- K-S test ----------------------
    print("KS test")
    if(candidate[dist] == "norm"){xi.cdf <- get(dist.char[3])(var, par1,par2)}
    if(candidate[dist] == "lnorm"){xi.cdf <- get(dist.char[3])(var, meanlog = par1, sdlog = par2)}
    if(candidate[dist] == "gumbel"){xi.cdf <- get(dist.char[3])(var, par1, par2)}
    if(candidate[dist] == "gamma3"){xi.cdf <- get(dist.char[3])(var, shape=par1, scale=par2, thres=par3)}
    if(candidate[dist] == "lgamma3"){xi.cdf <- get(dist.char[3])(var, shape=par1, scale=par2, thres=par3)}
    
    result <- ks.test(var, dist.char[3], par1,par2)
    print(paste0(candidate[dist], "的KS檢定P-value: ", result$p.value))
    
    # ------------- 將P-value整理成表格-------------------
    
    dist.table[dist,i] <- result$p.value
    
  }
  # 每個延時P-value排序，變成數值再排序
  dist.choice <- rank(as.numeric(dist.table[(1:dist),i])) 
  dist.table[length(candidate)+1,i] <- candidate[which.max(dist.choice)] #最大的P-value對應的機率分布
  
}


# ----------------------- 輸出資料 ----------------

file <- paste("E:/R_output/", "rainfall_return.xlsx", sep="")
write.xlsx(m, file) # 輸出檔：每年各延時的降雨深度(mm)
file <- paste("E:/R_output/", "dist_choice.xlsx", sep="")
write.xlsx(dist.table,file) #輸出檔：K-S test 的p-value 及各延時最適合機率分布


# 決定適合分布：gumbel
## 各重現期的水文量

#輸入分布名稱："q"+gunbel

Ti <- c()
Ty <- c() # 各延時不同重現期對應的水文量
for(i in 1:24){
  var <- m[,i]
  dist.char <- c("qgumbel")
  md <- eevd(var,method = "mle") 
  par1 <- md$parameters[1] #參數1
  par2 <- md$parameters[2] #參數2
  
  Q2 <- get(dist.char)(0.5,par1,par2)
  print(paste0("2年重現期洪峰流量",Q2))
  Q10 <- get(dist.char)(0.9,par1,par2)
  print(paste0("10年重現期洪峰流量",Q10))
  Q50 <- get(dist.char)(0.98,par1,par2)
  print(paste0("50年重現期洪峰流量",Q50))
  Q100 <- get(dist.char)(0.99,par1,par2)
  print(paste0("100年重現期洪峰流量",Q100))
  
  Ti <- rbind(Q2,Q10,Q50,Q100)
  Ty <- append(Ty,Ti)
}

# 數據整理
hydro_value <- matrix(Ty,ncol=24)


# 累積雨量深度/各延時 = 降雨強度(mm/hr)
for (i in c(1:24)){
  hydro_value[,i] <- hydro_value[,i]/i
}

colnames(hydro_value) <- c(paste0("h",1:24))
rownames(hydro_value) <- c("2-","10-","25-","100-")

file <- paste("E:/R_output/", "rainfall_intensity_return_period.xlsx", sep="")
write.xlsx(hydro_value,file)

#
horner_para <- c() #長度待修改
dura <- c((1:24)*60)  #延時：小時->分鐘
for (i in c(1:4)){ #重現期數量
  hydro_value_i <- hydro_value[i,]
  horner <- nls(log((hydro_value_i)) ~ log(a)-c*log(dura+b), #Hornor empirical formula 
                start = list(a = 1, b = 1, c = 1), control = list(maxiter = 50), trace = T) 
  para <- coef(horner)
  horner_para <- append(horner_para,para)
  #c <- cbind(horner_para)
}
#c <- cbind(horner_para)
horner_para <- rbind(horner_para[1:3],horner_para[4:6],horner_para[7:9],horner_para[10:12])
rownames(horner_para) <- c(paste0("重現期",c(2,10,50,100),"年"))
horner_para <- round(horner_para, digits = 2)
print(horner_para)

file <- paste("E:/R_output/", "dist_choice.xlsx", sep="")
write.xlsx(dist.table,file)

#----繪 IDF 曲線 plot IDF Curve---- 
color <- c("Red", "Black", "Blue", "aquamarine3")
Tyr <- c(2,10,50,100) # 重現期
for(i in c(1:length(Tyr))){ 
  a <- horner_para[i, 1]
  b <- horner_para[i, 2]
  c <- horner_para[i, 3]
  if(i == 1 ){ 
    plot(x = dura, y = a/(dura+b)^c, type = "l", lwd = 2, col = color[i], xlim = c(0, 1440), 
         ylim = c(0, 180), 
         xlab = "duration(min)", ylab = "intensity(mm/hr)", main = "Horner IDF Curve", 
         cex.lab = 1.5,cex.axis = 1.5, cex.main = 2) 
    par(new = T)
    plot(x = dura, y = hydro_value[i,], type = "p", lwd = 2, col = color[i], xlim = c(0, 1440),
         ylim = c(0, 180),
         xlab = "duration(min)", ylab = "intensity(mm/hr)", main = "Horner IDF Curve",
         cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
  }else{ 
    plot(x = dura, y = a/(dura+b)^c, type = "l", lwd = 2, col = color[i], xlim = c(0, 1440), 
         ylim = c(0, 180), xlab = "", ylab = "", axes = F)
    par(new = T)
    plot(x = dura, y = hydro_value[i,], type = "p", lwd = 2, col = color[i], xlim = c(0, 1440),
         ylim = c(0, 180), xlab = "", ylab = "", axes = F)
  } 
  par(new = T) 
} 
legend(legend = (paste0(Tyr, "yr")), lty = 1, lwd = 5, col = color,"topright", cex = 1) 
legend(legend = c("empirical"), pch=1, x=200,y=140, cex = 1)
legend(legend = c("expectation"), lty = 1 , x=200,y=180, cex = 1)
par(new = F) 

file <- paste("E:/R_output/", "horner_para.xlsx", sep="")
write.xlsx(horner_para,file)
