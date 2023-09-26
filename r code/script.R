library(rEDM)
library(tidyverse)
library(latex2exp)
library(ggplot2)
library(plotly)
library(tayloRswift)

plotBM <- function(BM, T, plot=1, type="BM"){
  df <- data.frame(
    Time = seq(0,T, T/(length(BM[1,])-1)),
    bm = t(BM)
  )
  theme_set(theme_bw())
  if(plot==1){
    data <- df %>%
      pivot_longer(!Time)
    return(ggplot(data=data, aes(x=Time,y=value, color=name)) + 
             geom_line(linewidth=0.3, show.legend=FALSE) +
             xlab("t") + 
             ylab(ifelse(type=="BM",unname(latex2exp::TeX("$B_t")),
                         unname(latex2exp::TeX("$X_t")))) + 
             scale_color_taylor(palette="lover")
    )
  }
  if(plot==2){
    return(ggplot(data=df, aes(x=bm.1, y=bm.2,color=Time))+geom_path(linewidth=0.3) + 
             xlab(unname(latex2exp::TeX("$B^{(1)}"))) + 
             ylab(unname(latex2exp::TeX("$B^{(2)}"))))+ 
      scale_color_taylor(palette="lover")
  }
  if(plot==3){
    return(plot_ly(data=df,x=~bm.1,y=~bm.2,z=~bm.3,
                   type='scatter3d',mode='lines',color=~Time))
  }
  else{
    return("Error")
  }
}


simOU <- function(dim, rho, driftMat, time, length, X0){
  df <- t(matrix(rep(X0,length),ncol=dim,nrow=length))
  for(i in 1:(length-1)){
    df[,i+1] <- as.vector( df[,i] - driftMat %*% df[,i] * time/(length-1) + 
                             rho*rnorm(dim, mean = 0, sd = sqrt(time/(length-1)))
    )
  }
  return(df)
}
simLog2 <- function(xtoy, ytox, rx, ry, length, X0){
  df <- t(matrix(rep(X0,length), ncol=2,nrow=length))
  for(i in 1:(length-1)){
    df[,i+1] <- (df[,i])*(c(rx,ry)*(1-df[,i])-c(ytox,xtoy)*rev(df[,i]))
  }
  return(df)
}

M <- matrix(c(1,0,0,1,1,0,1,0,1),ncol=3,byrow=T)
X <- simOU(dim=3, rho=0.05, driftMat = M,
           time = 100, length = 1e5, X0 = rep(0,3))
J <- simLog2(xtoy=0,ytox=0.05,rx=3.65,ry=3.77,length=1e5, X0=rnorm(2))

plotBM(X, T = 100, type="OU")

df <- data.frame(Var1 = X[1,3000:4000],
                 Var2 = X[2,3000:4000],
                 Var3 = X[3,3000:4000])

simplex_output <- simplex(df, c(1,500),c(701,900))
par(mar = c(4, 4, 1, 1), mgp = c(2.5, 1, 0))  # set margins for plotting
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", 
     ylab = "Forecast Skill (rho)")


var1_xmap <- ccm(df, E=1, lib_column = "Var1",
                 target_column = "Var2", lib_sizes = seq(10, 140, by = 10), num_samples = 100, 
                 random_libs = TRUE, replace = TRUE, silent = TRUE)
var2_xmap <- ccm(df, E=1, lib_column = "Var2",
                 target_column = "Var1", lib_sizes = seq(10, 140, by = 10), num_samples = 100, 
                 random_libs = TRUE, replace = TRUE, silent = TRUE)

plot(var1_xmap$LibSize, var1_xmap$`Var1:Var2`, type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0,0.6))
lines(var1_xmap$LibSize, var1_xmap$`Var2:Var1`, col = "blue")
legend(x = "topleft", legend = c("var1 xmap", "var2 xmap"), col = c("red", 
                                                                                  "blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)




## Sardine data
library(readr)
NewPortTemp <- read_csv("data/Newport Pier/NewportBeach_TEMP_1924-202303.csv", 
                        col_types = cols(...8 = col_skip(), ...9 = col_skip(), 
                                         ...10 = col_skip(), ...11 = col_skip(), 
                                         TIME_PST = col_skip(), TIME_FLAG = col_skip(), 
                                         TEMP_FLAG = col_skip()), skip = 44)
View(NewPortTemp)

ScrippsTemp <- read_csv("data/Scripps Pier/LaJolla_TEMP_1916-202303.csv", 
                        col_types = cols(...10 = col_skip(), 
                                         ...11 = col_skip(), ...12 = col_skip(), 
                                         ...13 = col_skip(), ...14 = col_skip(), 
                                         TIME_PST = col_skip(), TIME_FLAG = col_skip(), 
                                         BOT_TEMP_C = col_skip(), BOT_FLAG = col_skip()), 
                        skip = 45)
View(ScrippsTemp)

anchovies <- read_csv("data/anchovies.csv")
anchovies <- anchovies[2:nrow(anchovies),]
View(anchovies)

sardines <- read_csv("data/sardines.csv")
sardines <- sardines[2:nrow(sardines),]
View(sardines)

library(readxl)
NewFish <- read_excel("data/Marine Fisheries Data Explorer Extract - 2023-09-25 11.10.49.xlsx")
View(NewFish)
