5+6

5*5^4-88*15

0/0

1/0

result1<-5+6
result1

result2<-10*5.2
result2

#the same, but confusing:

5+6->result3
result3

sum(5,6)

max(5,6)

log(5,6)

seq(1,100,by=2)

plot(seq(1,100, by=20))

result5<-c("A", 1:100) #c() umwandeln in Vektoren
result5
summary(result5)
str(result5)
head(result5)
ls()

rm(result2)
ls()

Niederschlag<-c(56,46,50,53,69,83,83,80,62,55,60,63)
Niederschlag
plot(Niederschlag)
plot(Niederschlag, pch=19, cex=2, col="#00ff0060")

lines(lowess(Niederschlag,f=.2))

############################################################

random<-c(1,5,32,42,32,21,2,32,34)
random
plot(random, pch=15, cex=1.4, col="#00ff0060")
plot(log(random), pch=15, cex=1.4, col="#00ff0060")

#####################Niederschlag D. Beispiel##############

install.packages("raster")
library(raster)

ger<-getData("GADM", country="DEU", level=2) #get country borders
plot(ger)

prec<-getData("worldclim",var="prec",res=.5,lon=10,lat=51)
plot(prec)

prec_ger1<-crop(prec,germany)
spplot(prec_ger1)

prec_ger2<-mask(prec_ger1, germany)
spplot(prec_ger2)

prec_avg<-cellStats(prec_ger2,stat="mean")
plot(prec_avg)

#################Dokument erstellen###########################

install.packages("rmarkdown")
install.packages("knitr")

install.packages("officer")
library("rmarkdown")
library("knitr")
library("officer")

#######################ggplot2##############################

ggplot(housing, aes(x = Home.Value)) + geom_histogram()

install.packages("ggplot2")
library("ggplot2")

x11()
x<-data.frame(x=1,y=2,label="ggplot2 Test /n@ EAGLE")
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

a<-ggplot()+geom_point(data=mpg,aes(x=displ, y=hwy, colour=class))
a+theme_bw()
theme_set(theme_bw())

ggplot()+geom_point(data=mpg,aes(x=displ,y=hwy,colour=class))+
  facet_grit(manufacturer"class")+
  ggtitle("TEST chart")+
  theme(plot.title=element_text(angle=0,size=22,colour="hotpink"))+
  scale_colour_discrete(name="type")
