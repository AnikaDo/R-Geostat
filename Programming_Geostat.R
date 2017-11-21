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

##Raster-Einführung

'create own raster'
install.packages("raster")
library(raster)

r1 <- raster(nrows=10,ncols=10)
r1[] <- rnorm(100)               #[]heißt: r1 wird nicht überschreiben, sondern die Werte werden eingefügt, werden 100 Zahlen eingefügt
plot(r1)

library(sp)
poi1 <- cbind(c(rnorm(10)),c(rnorm(10)))    #random set of values created as step-in coordinates
poi1
poi1.sp <- SpatialPoints(poi1)
plot(poi1.sp)                               #plot spatial point data set
df <- data.frame(attr1=c('a','b','z','d','e','q','w','r','z','y'),attr2=c(101:110))
df
poi1.spdf <- SpatialPointsDataFrame(poi1.sp,df)
plot(poi1.spdf)#hier je nach Intervall verschiedene Farben einfügen)      #wie geht das?

'Herumspielen'
rbPal <- colorRampPalette(c('red','blue'))
test <- rbPal(10)[as.numeric(cut(df$attr2,breaks=3))]
plot(test)

##Weiter arbeiten mit Rasterdaten
install.packages("RStoolbox")
library(RStoolbox)

lsat
plot(lsat$B1_dn)   #oder plot(lsat[[1]])
B2_B3 <- c(lsat$B2_dn,lsat$B3_dn)   #lsat[[2:3]]
B2_B3
B2_30_40 <- lsat$B2_dn[30:40]     #DNs zwischen 30 und 40 extrahiert
B2_30_40
plot(B2_30_40)
extract(r1,poi1.spdf)             #Extrahieren von Raster Values 


install.packages("move")
library(move)
data(leroy)
plot(leroy)
env <- raster(leroy,vals=rnorm(100))
x <- lsat[1:10,]
x <- lsat[]
x <- getValues(lsat)
x <- lsat[lsat=10]
x <- extract(env,leroy)
plot(x)

lsat[] <- rnorm(ncell(lsat))
lsat[lsat<0] <- NA
env[] <- 0
env[leroy] <- 1
plot(env)
plot(lsat)


##########Uebung 4 - 21.11.2017##############

a <- sqrt(2)
if(a*a != 2)               #!= ungleich
    {                      #if statement in spatial context, i.e. together with resolution - resample bla  
  print('R is great!')
    }

j <- 0
while (j<1)
    {
  j <- j+0.1;print(j)
    }

                            #macht irgendeinen Unterschied wie rum man es schreibt. noch einmal Felix fragen 
j <- 0
while (j<1)
{
print(j);j <- j+0.1
}

myfunction <- function(x,y){
  z <- x+y
  return(z)}              #return is just needed if you have many values in your function then just z is returned

myfunction(4,3)

myfunction <- function(x,y){
  z <- x+y
  }                       #so wird bloß das Ergebnis ausgegeben, unabhängig von anderen Variablen, da diese nicht weiter definiert

myfunction(4,3)           #why are functions important in RS? Calculation of different indices for example

fun_ndvi <- function(nir,red){(nir-red)/(nir+red)}

raster()  #single layer raster
brick()   #multi-layer raster from one file
stack()   #multi-layer raster from seperate files (same extent, resolution)

#raster imports just one band at a time
band_1 <- raster('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif', band=1)
band_2 <- raster('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif', band=2)
band_3 <- raster('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif', band=3)

#combine rasters of identical dimensions from raster objects
allbands <- stack(band_1,band_2,band_3)

#combine rasters of identical dimensions from raster objects
stacked <- stack(c(''))

#brick imports all bands of a single file
allbands <- brick('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif')
allbands

img <- brick('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif')
img

cellStats()       #wie geht das?
summary()
zonal()
quantile()
freq()

#displaying your multi-spectral datra on RGB
plotRGB(allbands,3,2,1)

#maybe colour stretch is needed
plotRGB(allbands,3,2,1,stretch='lin')

#a ggplot2 option using the commands provided by package 'RStoolbox'
library(ggplot2)
ggRGB(allbands,3,2,1,stretch='lin')

#single layer greyscale
ggR(allbands,layer=4,maxpixels=1e6,stretch='hist')

#single layer map to user defined legend
ggR(allbands,layer=1,stretch='lin',geom_raster=TRUE)+scale_fill_gradient(low='blue',high='green')

#export raster - overwrite if already existent
writeRaster(img,datatype='FLT4S',filename='new_data.tif',format='GTiff',overwrite=TRUE)

#export a picture to GoogleEarth, only works for geographic coordinates
KML(img,new_data.tif,col=rainbow(255),maxpixels=100000)       #irgendwas funktioniert hier nicht

plot(band_3)

#draw an extent on the monitor (NW corner and SE corner)
ext <- drawExtent()

#ext is an object of all class extent
band_3_crop <- crop(band_3, ext)

#grow and shrink extents by multiplying
ext*2  #grows extent in all four directions

#cropping funktioniert auch so

'Raster Algebra'

#band calculation
raster_sd <- calc(img,fun=sd)

#adding a calculation into a function
fun <- function(x){x/10}
raster_output <- calc(img,fun)

#set NA values to -999
fun <- function(x){x[is.na(x)] <- -999;return(x)}
raster_output <- calc(img, fun)
plot(raster_output)

#refer to single layers
raster_output <- calc(img,fun=function(x){x[1]+x[2]*x[3]})

#regression analysis
raster12 <- stack(raster_1,raster_2)
fun <- function(x){lm(x[1:5]~x[6:10])$coefficients[2]}
raster_output <- calc(raster12,fun)

#write your permanent results directly to disk when calculating them
calc(img,fun=sd,filename='img_sd.grd')


##NDVI calculation
#import all layers all point to single bands with [[]]
lsat <- brick('E:/R_und_Geostatistics/crop_p224r63_all_bands.tif')
ndvi <- (lsat[[4]]-lsat[[3]])/(lsat[[4]]+lsat[[3]])
plot(ndvi)

#or import single band rasters and point to individual objects
band_3 <- raster('E:/R_und_Geostatistics/raster_band3.tif')
band_4

'overlay'
ndvi <- overlay(lsat[[4]],lsat[[3]],fun=function(nir,red){(nir-red)/(nir+red)})
plot(ndvi)
'calc'
ndvi <- calc(lsat,fun=function(x){(x[,4]-x[,3])/(x[,4]+x[,3])})
plot(ndvi)
'exportimg resulting VI image'
savi <- overlay(lsat[[4]],lsat[[3]],fun=function(nir,red){(nir-red)/(nir+red+0.5)*(1+0.5)},filename='savi.tif',format='GTiff')

#we use defined NDVI function, so we dont have to write it several times
fun_ndvi <- function(nir,red){(nir-red)/(nir+red)}

#MSAVI Funktion
fun_msavi <- function(nir,red){(2*nir+1-sqrt((2*nir+1)^2-8*(nir-red)))/2}
msavi <- overlay(lsat[[4]],lsat[[3]],fun=fun_msavi)         'IRGENDWO IST HIER NOCH EIN FEHLER!'
plot(msavi)

#An alternative command aus der RStoolbox
spectralIndices()
