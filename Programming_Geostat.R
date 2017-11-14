########Uebung 1 - 25.10.2017################

##use as a calculator
5+6
5*5^4-88*15
0/0
1/0

result1 <- 5+6
result1

result2 <- 10*5.2
result2



##operators
seq(1,100,by=2)         #generates a user defined sequence of data - 1 bis 100 in 2er Schritten

plot(seq(1,100,by=2)) #gibt Grafik aus

c("A", 1:100) #c() Umwandlung in Vektoren, kombiniert A und 1-100

getwd()




##Precipitation Data

prec_avg <- c(56,46,50,53,69,83,83,80,62,55,60,63)
plot(prec_avg)
plot(prec_avg, pch=19, cex=2, col="#00ff0060")
lines(lowess(prec_avg,f=.2))


random <-c (1,5,32,42,32,21,2,32,34)
plot(random, pch=15, cex=1.4, col="#00ff0060")
plot(log(random), pch=15, cex=1.4, col="#00ff0060")



##Niederschlag Deutschland 
install.packages("raster")
library(raster)

germany <- getData("GADM",country="DEU",level=2) #get country borders
plot(germany)

prec_ger <- getData("worldclim",var="prec",res=.5,lon=10,lat=51)
plot(prec_ger)

prec_ger1 <- crop(prec_ger,germany)     #crop to extent of Germany
spplot(prec_ger1)

prec_ger2 <- mask(prec_ger1,germany)    #mask to shape of Germany
spplot(prec_ger2)

prec_avg <- cellStats(prec_ger2,stat="mean")
plot(prec_avg)


##Rmarkdown & Knitr installieren
install.packages("rmarkdown")
install.packages("knitr")
install.packages("officer")
library("rmarkdown")
library("knitr")
library("officer")

eval=false #code not executed
echo=false #code not shown

##ggplot2
install.packages("ggplot2")
library("ggplot2")

x11()
x <- data.frame(x=1,y=2,label="ggplot2 Test /n@ EAGLE")
ggplot(data=x,aes(x=x,y=y))+geom_text(aes(label=label),size=15)

head(mpg)
ggplot(mpg, aes(x=displ, y=hwy))+geom_point()
ggplot(mpg,aes(displ,cty,colour=class)) + geom_point() + geom_smooth()
ggplot(mpg,aes(displ,hwy))+geom_point()+facet_wrap("class")+geom_smooth()

myPlot<-ggplot(mpg,aes(x=displ,y=hwy))+geom_point()
myPlot + geom_smooth()

ggplot()+geom_point(data=mpg,aes(x=displ,y=hwy))

ggplot(mpg,aes(drv,hwy))+geom_jitter()
ggplot(mpg,aes(drv,hwy))+geom_boxplot()
ggplot(mpg,aes(drv,hwy))+geom_violin()

ggplot(mpg,aes(drv,hwy))+geom_violin()+geom_jitter()+geom_boxplot()
ggplot(mpg,aes(drv,hwy))+geom_violin()+geom_jitter(aes(alpha=.7,size=2),colour="blue")

a <- ggplot()+geom_point(data=mpg,aes(x=displ, y=hwy, colour=class))
a+theme_bw()
theme_set(theme_bw())

ggplot()+geom_point(data=mpg,aes(x=displ,y=hwy,colour=class))+
  facet_grit(manufacturer"class")+
  ggtitle("TEST chart")+
  theme(plot.title=element_text(angle=0,size=22,colour="hotpink"))+
  scale_colour_discrete(name="type")


#########Uebung 2 - 07.11.2017##############
getwd()
setwd("D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day2")

read.table("bio_data_forest.csv")    #wieso klappt das hier nicht?
read.csv("bio_data_forest.csv")
read.table('tabelle.txt')
my.df <- read.table('tabelle.txt'),header=TRUE,sep='/')


##Indexing - Matrix
x <- matrix(c(4,7,3,8,9,2),nrow=2)
x
x[2,2]
x[,2]

'create a data frame or matrix'
numbers_1 <- rnorm(80,mean=0,sd=1)
mat_1 <- matrix(numbers_1,nrow=20,ncol=4)
mat_1

df_1 <- data.frame(mat_1)
names(df_1) <- c('var1','var2','var3','var4')
head(df_1)

##Indexing - Vector
x <- seq(1,100,by=2.5)
x
x[5]
x[4:10]
x[length(x)]
x[length(x)-1]

####Data Frame vs. Matrix: Matrix nur numerische Werte, keine Variablennamen, Data Frame kann auch verschiedene Werte und Namen

##Indexing - Data Frame
test <- data.frame(A=c(1,2,3),B=c('aB1','aB2','aB3'))
test[,1]
test[,'A']                                   #same thing

df <- data.frame(plot='location_name_1',measure1=runif(100)*1000,measure2=round(runif(100)*100),value=rnorm(100,2,1),ID=rep(LETTERS,100))
df_2 <- data.frame(plot='location_name_2',measure1=runif(50)*100,measure2=round(runif(50)*10),value=rnorm(50),ID=rep(LETTERS,50))
df <- rbind(df,df_2)                         #data merging
summary(df)
str(df)
head(df)
length(df$measure1)
df[,c('plot','measure1','measure2')]
a <- df[,c('plot','measure1','measure2')]
a
df[100,c('plot','measure1','measure2')]

##quering single values
prec_avg[7]       #prec July
prec_avg[4:9]     #prec April-June

##task for next week
getwd()
#setwd("D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day2")
setwd('E:/Introduction to R and Geostatistics/day2')

x <- read.table('pH.csv',sep=';',dec=',')         #seperator definiert, genauso wie Komma für Dezimalstellen
names(x) <- c('H2O','KCL','CaCL2')                #gibt es auch was, um die Zeilen zu benennen?
x
summary(x)
str(x)
mode(x)
cut(x)
sort('pH.csv',x,decreasing=FALSE)                 #was muss hier verändert werden?



########Uebung3 - 14.11.17#########

#erst einmal Wiederholung der letzten Stunde

library(car)
x <- seq(1,100,by=1)
x
recode(x,'0:30=1;31:70=2;else=3')                 #einfacher als manuell: recode aus car package

##Indexing again
m1 <- matrix(c(4,7,3,8,9,2),nrow=2)
m1
m2 <- matrix(c(2,4,3,1,5,7),nrow=2,ncol=3,byrow=TRUE)
m2

'Fortsetzung der letzten Woche'
df <- data.frame(plot='location_name_1',measure1=runif(100)*1000,measure2=round(runif(100)*100),value=rnorm(100,2,1),ID=rep(LETTERS,100))
df_2 <- data.frame(plot='location_name_2',measure1=runif(50)*100,measure2=round(runif(50)*10),value=rnorm(50),ID=rep(LETTERS,50))
df <- rbind(df,df_2)                         #data merging
summary(df)
str(df)
head(df)
length(df$measure1)
df[,c('plot','measure1','measure2')]
a <- df[,c('plot','measure1','measure2')]
b <- df[,c('measure1')]
plot(b)
c <- b <- df[,c('measure2')]
plot(c)
d <- df[,c('plot','measure1')]
plot(d)
