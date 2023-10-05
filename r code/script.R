library(rEDM)
library(tidyverse)
library(latex2exp)
library(ggplot2)
library(plotly)
library(tayloRswift)
library(extrafont)

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


sardines <- read_csv("data/sardines.csv")
sardines <- sardines[2:nrow(sardines),]


library(readxl)
NewFish <- read_excel("data/Marine Fisheries Data Explorer Extract - 2023-09-25 11.10.49.xlsx")
new_sardines <- NewFish[NewFish$`Species Name`=="Sardine, Pacific",1:4]
new_anchovies <- NewFish[NewFish$`Species Name`=="Anchovy, northern",1:4]


tail(anchovies)
anchovies$month <- sardines$month <- rep(1:12,75)
comb <- merge(new_sardines, new_anchovies, by=c("Year", "Month"))

df <- data.frame(
  year <- c(anchovies$year, comb$Year[comb$Year>2002]),
  month <- c(anchovies$month, comb$Month[comb$Year>2002]),
  anchovies <- c(anchovies$landings, comb$Pounds.y[comb$Year > 2002]),
  sardines <- c(sardines$landings, comb$Pounds.x[comb$Year > 2002])
)
colnames(df) <- c("year", "month", "anchovies", "sardines")
df <- as_tibble(df)
df

library(zoo)
roll_anchovies <- zoo::rollmean(na.omit(as.numeric(df$anchovies)), k=48, fill=NA)
roll_sardines <- zoo::rollmean(na.omit(as.numeric(df$sardines)), k=48, fill=NA)

NewPortTemp <- group_by(NewPortTemp,  YEAR,MONTH)
NewPortTemp <- summarise(NewPortTemp, result=mean(SURF_TEMP_C))
roll_temp <- zoo::rollmean(na.omit(as.numeric(NewPortTemp$result)), k=36, fill=NA)
ScrippsTemp <- group_by(ScrippsTemp,  YEAR,MONTH)
ScrippsTemp <- summarise(ScrippsTemp, result=mean(SURF_TEMP_C))
roll_temp2 <- zoo::rollmean(na.omit(as.numeric(ScrippsTemp$result)), k=36, fill=NA)


yearFish <- group_by(df[,1:3], year)
yearFish$anchovies <- as.numeric(yearFish$anchovies)
yearFish <- na.omit(yearFish)
yearAnchovies <- summarise(yearFish, result=mean(as.numeric(anchovies)))


yearFish <- group_by(df[,c(1,2,4)], year)
yearFish$sardines <- as.numeric(yearFish$sardines)
yearFish <- na.omit(yearFish)
yearSardines <- summarise(yearFish, result=mean(as.numeric(sardines)))


tempNA <- na.omit(NewPortTemp)
temp2NA <- na.omit(ScrippsTemp)
tempDF <- merge(tempNA,temp2NA, by=c("YEAR", "MONTH"))
roll_temp <- zoo::rollmean(tempDF$result.x, k=36, fill=NA)
roll_temp2 <- zoo::rollmean(tempDF$result.y, k=36, fill=NA)
cor(na.omit(roll_temp),na.omit(roll_temp2))

library(ggtext)

p <- qplot((1928+1:1090*(94/1090)),roll_anchovies[1:1090]/1e6,geom="line", color="#ed713a") + 
  theme(
    text = element_text(family = "IBM Plex Sans Condensed", size = 10)
  ) + 
  theme(axis.title = element_text(), legend.position = "none",
        plot.subtitle = element_markdown(),
        axis.title.y = element_markdown(),
        axis.title.y.right = element_markdown()) +
  geom_line(aes(y=roll_sardines[1:1090]/5e6),color="#3fc1c9") +
  geom_line(aes(x=1928+1:740*(94/740),y=(roll_temp[1:740]+roll_temp2)*5-1.32e2), color="#77AB43") +
  scale_y_continuous(
    # Features of the first axis
    name = "Landings of <span style = 'color: #ed713a;'>anchovy</span> (million pounds)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*5, name="Landings of <span style = 'color: #3fc1c9;'>sardine</span> (million pounds)")
  ) +
  geom_bar(aes(x=1928:2022, y=yearAnchovies$result/1e6), stat="identity", fill="#ed713a", alpha=0.2, linetype=0) +
  geom_bar(aes(x=1928:2022, y=yearSardines$result/5e6), stat="identity", fill="#3fc1c9", alpha=0.2, linetype=0)  +
  geom_hline(yintercept = 27.06474, linetype=2, color ="#77AB43", alpha=0.5) +
  geom_hline(yintercept = 45.46639, linetype=2, color ="#77AB43", alpha=0.5) +
  annotate("richtext", x=2015, y=25, 
           label="<span style = 'color: #77AB43;'>15.91 °C</span>",  
           fill = NA, label.color = NA,size=3.4, 
           family="IBM Plex Sans Condensed") +
  annotate("richtext", x=1935, y=43.4, 
           label="<span style = 'color: #77AB43;'>17.75 °C</span>",  
           fill = NA, label.color = NA,size=3.4, 
           family="IBM Plex Sans Condensed") +
  labs(
    title="Sardine and Anchovy Landings",
    subtitle="4-year and yearly averages landings of 
      <span style = 'color: #3fc1c9;'>sardine</span> and
      <span style = 'color: #ed713a;'>anchovy</span> compared with 4-year averages of
    <br> <span style = 'color: #77AB43;'>sea surface temperature</span> from 1928-2022 in the California Current System.",
    caption="Source: Calif. Dept. of Wildlife and Game",
    x="Year"
  )

p

ggsave(
  "sardine.png",
  p,
  width = 6.5,
  height = 4.5
)
  

max(na.omit((roll_temp[1:740]+roll_temp2)/2))


qplot(1:883,diff(roll_anchovies[1:884]),geom="line", color="red", alpha=0.8) + geom_line(aes(y=diff(roll_sardines[1:884])/5),color="blue",alpha=0.8) +
  geom_line(aes(y=(diff(roll_temp[1:884]+roll_temp2))*1e7/2), color="green", alpha=0.3)

rollCor <- zoo::rollapply(data.frame(df1 <- roll_temp,
                                     df2 <- roll_temp2)
                          , fill=NA,width=10, 
                          FUN=function(df) cor(df[,1],df[,2]), 
                          by.column=FALSE)
rollCorAnchSar <- zoo::rollapply(data.frame(df1 <- roll_anchovies,
                                            df2 <- roll_sardines[1:1097])
                                 , fill=NA,width=360, 
                                 FUN=function(df) cor(df[,1],df[,2]), 
                                 by.column=FALSE)
rollCorTempSar <- zoo::rollapply(data.frame(df1 <- roll_temp,
                                            df2 <- roll_sardines[sort(s)])
                                 , fill=NA,width=360, 
                                 FUN=function(df) cor(df[,1],df[,2]), 
                                 by.column=FALSE)
rollCorTempSar
s <- sample(1:1119,740)
sort(s)

cor(diff(na.omit(roll_anchovies)[800:1000]), diff(na.omit(roll_sardines)[800:1000]))
roll_sardines[1090]
roll_anchovies[800]
summary(tempDF)
summary(df)


# Sample data
data <- data.frame(
  Category = c("A", "B", "C", "D", "E"),
  Value = c(25, 40, 30, 15, 20)
)

# Create the ggplot
p <- ggplot(data, aes(x = Category, y = Value)) +
  geom_bar(stat = "identity", fill = "#008FD5") +
  
  # Adjust the theme to mimic FiveThirtyEight style
  theme_minimal() +
  theme(
    text = element_text(family = "IBM Plex Sans Condensed", size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(vjust = 0.5, hjust = 0.5),
    panel.grid.major = element_line(size = 0.2, color = "gray"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "#F0F0F0", color = NA),
    plot.background = element_rect(fill = "#F0F0F0", color = NA),
    legend.position = "none"
  ) +
  
  # Add a title and subtitle
  labs(
    title = "FiveThirtyEight-Inspired Bar Plot",
    subtitle = "A demonstration of a ggplot2 plot in FiveThirtyEight style",
    caption = "Source: Your Data Source"
  )

# Display the plot
print(p)

library(ggthemes)
theme_set(theme_fivethirtyeight())
