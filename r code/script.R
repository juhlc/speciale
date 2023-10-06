## Base packages
{
library(readr)  ## Loading data
library(readxl) ## Loading Excel data

library(rEDM)  ## CCM and other related methods

library(tidyverse)  ## Data manipulation
library(zoo)  ## Computing functions rolling

library(ggplot2) ## Plotting
library(latex2exp)  ## Adding latex in plots
library(ggthemes) ## Using themes in ggplot
library(extrafont)  ## Adding fonts in R
library(ggtext)  ## Using fonts with ggplot
library(tayloRswift) ## Coloring
library(plotly) ## 3d plotting
}

## Functions for simulation of OU- and coupled deterministic logistic processes
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

## Simulation and test of rEDM package
{
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
                 target_column = "Var2", lib_sizes = seq(10, 140, by = 10), 
                 num_samples = 100, random_libs = TRUE, replace = TRUE, silent = TRUE)
var2_xmap <- ccm(df, E=1, lib_column = "Var2",
                 target_column = "Var1", lib_sizes = seq(10, 140, by = 10), 
                 num_samples = 100, random_libs = TRUE, replace = TRUE, silent = TRUE)

plot(var1_xmap$LibSize, var1_xmap$`Var1:Var2`, type = "l", col = "red", 
     xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0,0.6))
lines(var1_xmap$LibSize, var1_xmap$`Var2:Var1`, col = "blue")
legend(x = "topleft", legend = c("var1 xmap", "var2 xmap"), 
       col = c("red","blue"), lwd = 1, bty = "n", 
       inset = 0.02, cex = 0.8)
}


                  ############################
                  ##      Sardine data      ##
                  ############################

lag <- 3


## Reading and cleaning data
{
NewPortTemp <- read_csv("data/Newport Pier/NewportBeach_TEMP_1924-202303.csv", 
                        col_types = cols(...8 = col_skip(), ...9 = col_skip(), 
                                         ...10 = col_skip(), ...11 = col_skip(), 
                                         TIME_PST = col_skip(), TIME_FLAG = col_skip(), 
                                         TEMP_FLAG = col_skip()), skip = 44)
ScrippsTemp <- read_csv("data/Scripps Pier/LaJolla_TEMP_1916-202303.csv", 
                        col_types = cols(...10 = col_skip(), 
                                         ...11 = col_skip(), ...12 = col_skip(), 
                                         ...13 = col_skip(), ...14 = col_skip(), 
                                         TIME_PST = col_skip(), TIME_FLAG = col_skip(), 
                                         BOT_TEMP_C = col_skip(), BOT_FLAG = col_skip()), 
                        skip = 45)


anchovies <- read_csv("data/anchovies.csv")
anchovies <- anchovies[2:nrow(anchovies),]

sardines <- read_csv("data/sardines.csv")
sardines <- sardines[2:nrow(sardines),]

NewFish <- read_excel("data/Marine Fisheries Data Explorer Extract - 2023-09-25 11.10.49.xlsx")
new_sardines <- NewFish[NewFish$`Species Name`=="Sardine, Pacific",1:4]
new_anchovies <- NewFish[NewFish$`Species Name`=="Anchovy, northern",1:4]

anchovies$month <- sardines$month <- rep(1:12,75)
comb <- merge(new_sardines, new_anchovies, by=c("Year", "Month"))

df <- data.frame(
  year = c(anchovies$year, comb$Year[comb$Year>2002]),
  month = c(anchovies$month, comb$Month[comb$Year>2002]),
  anchovies = c(anchovies$landings, comb$Pounds.y[comb$Year > 2002]),
  sardines = c(sardines$landings, comb$Pounds.x[comb$Year > 2002])
)
df <- as_tibble(df)

## Missing values imputed by 0 - hopefully only producing a small bias
anchovies_by_year <- summarise(df,.by=year, result=sum(na.omit(as.numeric(anchovies))))
sardines_by_year <- summarise(df,.by=year, result=sum(na.omit(as.numeric(sardines))))

roll_anchovies <- rollmean(anchovies_by_year$result, k=lag, fill=NA)
roll_sardines <- rollmean(sardines_by_year$result, k=lag, fill=NA)

## Temperature data
NewPortTemp <- NewPortTemp[NewPortTemp$YEAR >= 1928 & NewPortTemp$YEAR <=2022,]
ScrippsTemp <- ScrippsTemp[ScrippsTemp$YEAR >= 1928 & ScrippsTemp$YEAR <=2022,]
mean(ScrippsTemp[,1:3] == NewPortTemp[,1:3])

## Impute NA values by measured values from other variable
AvTemp <- cbind(NewPortTemp[,1:3], 
                rowMeans(
                  cbind(
                    NewPortTemp$SURF_TEMP_C,
                    ScrippsTemp$SURF_TEMP_C),
                  na.rm = TRUE)
                )
colnames(AvTemp) <- c("year", "month", "day", "avg_temp")
# AvTemp[is.na(AvTemp$avg_temp),]
AvTemp <- as_tibble(AvTemp)
AvTemp <- group_by(AvTemp,  year, month)

## Remaining missing values imputed by the average temperature of that month
AvTemp <- summarise(AvTemp, result=mean(na.omit(avg_temp)))
temp_by_year <- summarise(AvTemp, result=mean(result))
roll_temp <- rollmean(AvTemp$result, k=lag*12, fill=NA)
roll_temp_year <- rollmean(temp_by_year$result, k=lag, fill=NA)


## More sophisticated imputation
AvTemp2<- cbind(NewPortTemp[,1:3], 
                rowMeans(
                  cbind(
                    NewPortTemp$SURF_TEMP_C,
                    ScrippsTemp$SURF_TEMP_C),
                  na.rm = FALSE)
)
colnames(AvTemp2) <- c("year", "month", "day", "temp")
AvTemp2[1,4] <- AvTemp2[2,4] + ScrippsTemp[1,4] - ScrippsTemp[2,4]
for(i in 2:nrow(AvTemp2)){
  if(AvTemp2$temp[i] == "NaN"){
    list <- rev(1:i)
    counter <- 0
    check <- FALSE
    while(check == FALSE){
      counter <- counter + 1
      while(AvTemp2[list[counter],4] == "NaN"){
        counter = counter + 1
      }
      check <- (NewPortTemp[i,4] != "NaN" & NewPortTemp[list[counter],4] != "NaN") |
        (ScrippsTemp[i,4] != "NaN" & ScrippsTemp[list[counter],4] != "NaN") |
        (ScrippsTemp[i,4] == "NaN" & NewPortTemp[i,4] == "NaN")
    }
    AvTemp2$temp[i] <- as.numeric(AvTemp2[list[counter], 4]) + 
      as.numeric(ifelse(ScrippsTemp[i,4] == "NaN", 0, - ScrippsTemp[list[counter], 4] + ScrippsTemp[i,4])) +
      as.numeric(ifelse(NewPortTemp[i,4] == "NaN", 0, - NewPortTemp[list[counter], 4] + NewPortTemp[i,4]))
  }
}
AvTemp2[is.na(AvTemp2$avg_temp),]
AvTemp2 <- as_tibble(AvTemp2)
AvTemp2 <- group_by(AvTemp2,  year, month)

## Remaining missing values imputed by the average temperature of that month
AvTemp2 <- summarise(AvTemp2, result=mean(temp))
temp_by_year2 <- summarise(AvTemp2, result=mean(result))
roll_temp2 <- rollmean(AvTemp2$result, k=lag*12, fill=NA)
roll_temp_year2 <- rollmean(temp_by_year2$result, k=lag, fill=NA)

# cor(AvTemp$result,AvTemp2$result)
# cor(temp_by_year$result, temp_by_year2$result)
# cor(na.omit(roll_temp), na.omit(roll_temp2))
# cor(na.omit(roll_temp_year), na.omit(roll_temp_year2))
}

## Structuring data for plotting and analysis
{
plot_data <- data.frame(
  year = 1928:2022,
  roll_anchovies = roll_anchovies,
  roll_sardines = roll_sardines,
  roll_temp = roll_temp_year2,
  roll_temp_alt = roll_temp_year,
  anchovies = anchovies_by_year$result,
  sardines = sardines_by_year$result,
  temp = temp_by_year2$result,
  temp_alt = temp_by_year$result
)
plot_data <- as_tibble(plot_data)

temp_data <- data.frame(
  year = 1:length(roll_temp2)/12+1928,
  temp = roll_temp2
)
temp_data <- as_tibble(temp_data)

rm(list=setdiff(ls(), c("temp_data", "plot_data", "lag")))
}

## Creating basic plot of sardines and anchovies
{
loadfonts(quiet=T)
theme_set(theme_fivethirtyeight())

p <- ggplot(data=plot_data, aes(x=year)) +
  theme(
    text = element_text(family = "IBM Plex Sans Condensed", size = 10),
    axis.title = element_text(),
    plot.subtitle = element_markdown(),
    axis.title.y = element_markdown(),
    axis.title.y.right = element_markdown(),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.position = "none") +
  scale_y_continuous(
    # Features of the first axis
    name = "Landings of <span style = 'color: #ed713a;'>anchovy</span> (million pounds)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*5, name="Landings of <span style = 'color: #3fc1c9;'>sardine</span> (million pounds)")
  ) +
  geom_line(aes(y=roll_sardines/5e6), color = "#3fc1c9") +
  geom_bar(aes(y=sardines/5e6), stat="identity", fill="#3fc1c9", alpha=0.2, linetype=0)  +
  geom_line(aes(y=roll_anchovies/1e6),color="#ed713a") +
  geom_bar(aes(y=anchovies/1e6), stat="identity", fill="#ed713a", alpha=0.2, linetype=0) +
  geom_line(aes(y=plot_data$roll_temp*60-600), color="#77AB43") +
  geom_hline(yintercept = min(na.omit(plot_data$roll_temp))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  geom_hline(yintercept = max(na.omit(plot_data$roll_temp))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  # geom_line(data=temp_data, aes(x=year, y=temp*60-600), color="#77AB43") +
  # geom_hline(yintercept = min(na.omit(roll_temp2))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  # geom_hline(yintercept = max(na.omit(roll_temp2))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  annotate("richtext", x=2015, y=325, 
           label="<span style = 'color: #77AB43;'>15.88 °C</span>",  
           fill = NA, label.color = NA,size=3.4, 
           family="IBM Plex Sans Condensed") +
  annotate("richtext", x=1935, y=510, 
           label="<span style = 'color: #77AB43;'>18.85 °C</span>",  
           fill = NA, label.color = NA,size=3.4, 
           family="IBM Plex Sans Condensed") +
  labs(
    title="Sardine and Anchovy Landings",
    subtitle="The yearly and 3-year averages of landings of 
      <span style = 'color: #3fc1c9;'>sardine</span> and
      <span style = 'color: #ed713a;'>anchovy</span> compared with 3-year
    <br> averages of <span style = 'color: #77AB43;'>sea surface temperature</span> from 1928-2022 in the California Current System.",
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
}

p

## First difference approach
{
corr_lag <- 30
diff_plot_data <- data.frame(
  year = 1929:2022,
  roll_anchovies = diff(plot_data$roll_anchovies),
  roll_sardines = diff(plot_data$roll_sardines),
  roll_temp = diff(plot_data$roll_temp),
  roll_temp_undiff = plot_data$roll_temp[2:nrow(plot_data)],
  anchovies = diff(plot_data$anchovies),
  sardines = diff(plot_data$sardines),
  roll_corr_a_s = rollapply(plot_data,
                             fill=NA, width=corr_lag,
                             FUN=function(df) cor(df[,"roll_anchovies"],
                                                  df[,"roll_sardines"]), 
                                 by.column=FALSE)[2:nrow(plot_data)],
  roll_corr_a_t = rollapply(plot_data,
                             fill=NA, width=corr_lag,
                             FUN=function(df) cor(df[,"roll_anchovies"],
                                                  df[,"roll_temp"]), 
                             by.column=FALSE)[2:nrow(plot_data)],
  roll_corr_s_t = rollapply(plot_data,
                             fill=NA, width=corr_lag,
                             FUN=function(df) cor(df[,"roll_temp"],
                                                  df[,"roll_sardines"]), 
                             by.column=FALSE)[2:nrow(plot_data)]
)
}

## Expanded plot
{
q <- ggplot(data=diff_plot_data, aes(x=year)) +
  theme(
    text = element_text(family = "IBM Plex Sans Condensed", size = 10),
    axis.title = element_text(),
    plot.subtitle = element_markdown(),
    axis.title.y = element_markdown(),
    axis.title.y.right = element_markdown(),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    legend.position = "none") +
  scale_y_continuous(
    # Features of the first axis
    name = "Net difference in landings of <span style = 'color: #ed713a;'>anchovy</span> (million pounds)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis( trans=~.*5, name="Net difference in landings of <span style = 'color: #3fc1c9;'>sardine</span> (million pounds)")
  ) +
  geom_line(aes(y=roll_sardines/5e6), color = "#3fc1c9") +
  geom_bar(aes(y=sardines/5e6), stat="identity", fill="#3fc1c9", alpha=0.2, linetype=0)  +
  geom_line(aes(y=roll_anchovies/1e6),color="#ed713a") +
  geom_bar(aes(y=anchovies/1e6), stat="identity", fill="#ed713a", alpha=0.2, linetype=0) +
  geom_line(aes(y=roll_corr_s_t*100),color="#77AB43") +
  #geom_line(aes(y=plot_data$roll_temp*60-600), color="#77AB43") +
  #geom_hline(yintercept = min(na.omit(plot_data$roll_temp))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  #geom_hline(yintercept = max(na.omit(plot_data$roll_temp))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  # geom_line(data=temp_data, aes(x=year, y=temp*60-600), color="#77AB43") +
  # geom_hline(yintercept = min(na.omit(roll_temp2))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  # geom_hline(yintercept = max(na.omit(roll_temp2))*60-600, linetype=2, color ="#77AB43", alpha=0.5) +
  #annotate("richtext", x=2015, y=325, 
   #        label="<span style = 'color: #77AB43;'>15.88 °C</span>",  
    #       fill = NA, label.color = NA,size=3.4, 
     #      family="IBM Plex Sans Condensed") +
  #annotate("richtext", x=1935, y=510, 
   #        label="<span style = 'color: #77AB43;'>18.85 °C</span>",  
    #       fill = NA, label.color = NA,size=3.4, 
     #      family="IBM Plex Sans Condensed") +
  labs(
    title="Sardine and Anchovy Landings",
    subtitle="The yearly and 3-year averages of landings of 
      <span style = 'color: #3fc1c9;'>sardine</span> and
      <span style = 'color: #ed713a;'>anchovy</span> compared with 3-year
    <br> averages of <span style = 'color: #77AB43;'>sea surface temperature</span> from 1928-2022 in the California Current System.",
    caption="Source: Calif. Dept. of Wildlife and Game",
    x="Year"
  )

q


ggsave(
  "diff.png",
  q,
  width = 6.5,
  height = 4.5
)
}


#### ChatGPT Style Guide
{
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
}
