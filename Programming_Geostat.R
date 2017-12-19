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


#######Uebung 5 - 28.11.17#################################################################

'a list'
a <- runif(199)
b <- c("aa","bb","cc","dd","ee")
c <- list(a,b)
c

'indexing a list of two vectors'
c[2]            #ist das gleiche wie: c[[2]] <- indexing the second object
c[[2]][1]       #first entry of second object

'...of different size'
a <- list(obj_1=runif(100),obj_2=c("aa","bb"),obj_3=c(1,2,4))
a$obj_1         #call the object name
'oder' 
a[["obj_1"]]
'oder'    
a[[2]]

'a list with a matrix, vector and data frame of different sizes'
a <- list(m1=matrix(runif(50),nrow = 5),v1=c(1,6,10),df1=data.frame(a=runif(100),b=rnorm(100)))
a$df1[,1]    #index a data frame or matrix as known

#if-else statement
a <- 5
if(a>0)
{
  print("it is a positive number")
}

a <- 5
if(a !=5)
{
  print("number is not equal 5")
} else {
  print("number is equal 5")
}

'while statement'
j <- 0
while(j<1)
{
  j <- j+0.1; print(j)
}

'add commonly used analyses to a function'
#myfunction <- function(arg1,arg2,...){
#  statements
#  return(object)
#  }

myfunction <- function(x,y){
  z <- x+y
  return(z)         #ohne return(): last created object is returned
  }

myfunction(4,3)

#Vegetation Indices
install.packages("raster")
library(raster)
install.packages("rgdal")
library(rgdal)

lsat <- brick("D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/UG_Landsat.tif")
ndvi <- (lsat[[4]]-lsat[[3]])/(lsat[[4]]+lsat[[3]])
plot(ndvi)
'oder'
band_3 <- raster("D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/UG_Landsat_B3.tif")
band_4 <- raster("D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/UG_Landsat_B4.tif")
ndvi <- (band_4-band_3)/(band_4+band_3)
plot(ndvi)

#functions in band calculations using commands
ndvi <- overlay(band_4, band_3, fun=function(nir,red){(nir-red)/(nir+red)})

ndvi <- calc(lsat, fun=function(x){(x[,4]-x[,3])/(x[,4]+x[,3])})

savi <- overlay(band_4, band_3, fun=function(nir,red0){(nir-red)/(nir+red+0.5)*1+0.5}, filename="savi.tif", format="GTiff")
     #FEHLERMELDUNG: cannot use this formula, because it is not vectorized

fun_ndvi <- function(nir,red){(nir-red)/(nir+red)}
ndvi <- overlay(band_4, band_3, fun=fun_ndvi)

install.packages("RStoolbox")
library(RStoolbox)
ndvi <- spectralIndices(lsat, red="band_3", nir="band_4", indices="NDVI")  #wieso klappt das nicht?
VIs <- spectralIndices(lsat_ref, red="B3_tre", nir="B4_tre")               #hier auch nicht


##Classification
'unsupervised'
install.packages("raster")
library(raster)
install.packages("rgdal")
library(rgdal)
install.packages("RStoolbox")
library(RStoolbox)
#band_1 <- raster('D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/crop_p224r63_all_bands.tif', band=1)
#band_2 <- raster('D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/crop_p224r63_all_bands.tif', band=2)
#band_3 <- raster('D:/Studium/Master/3. WiSe 17-18/Introduction to R and Geostatistics/day5/crop_p224r63_all_bands.tif', band=3)
band_1 <- raster('E:/Introduction to R and Geostatistics/crop_p224r63_all_bands.tif', band=1)
band_2 <- raster('E:/Introduction to R and Geostatistics/crop_p224r63_all_bands.tif', band=2)
band_3 <- raster('E:/Introduction to R and Geostatistics/crop_p224r63_all_bands.tif', band=3)
allbands <- stack(band_1,band_2,band_3)
uc <- unsuperClass(allbands,nClasses=5)         #uc is a list containing auxilary information
plot(uc$map)

landsat_allbands.kmeans <- kmeans(allbands[],5) #is a clustering algorithm, dont use the spatial object, just use the data of the spatial object
kmeansraster <- raster(allbands)
kmeansraster[] <- landsat_allbands.kmeans$cluster #insert it into the slot of kmeansraster
plot(kmeansraster)

'create a copy, filled with NAs'
kmeansraster <- allbands[[1]]     #create a copy, filled with NAs
kmeansraster[] <- NA              
values <- getValues(allbands)     #extract values from raster layers
valid <- complete.cases(values)   #just use complete cases
allbands.kmeans <- kmeans(values[valid,],5,iter.max=100,nstart=3)     #run the kmeans clustering
kmeansraster[valid] <- allbands.kmeans$cluster   #populate empty vector with cluster values derived from kmeans
plot(kmeansraster)
click(kmeansraster, n=3)
arg <- list(at=seq(1,5,1),labels=c("none","none","water","forest","defo"))     #argument scheme
colour <- c("white","white","blue","green","brown")                            #colour schmeme
plot(kmeansraster,col=colour, axis.arg=arg)


##ggplot2 intro and examples
install.packages("ggplot2")
library(ggplot2)
x11()
x <- data.frame(x=1,y=1,label="ggplot 2 introduction \n@ EAGLE")
ggplot(data=x, aes(x=x,y=y))+geom_text(aes(label=label),size=15)

install.packages("devtools")
library(devtools)

install_bitbucket("EAGLE_MSc/steigerwald",username=AnikaDo)
library(steigerwald)
data("bio_data")
head(bio_data)
ggplot(bio_data$forest_short,aes(x=beech,y=ndvi))+geom_point()

ggplot(bio_data$forest_short,aes(beech,ndvi,colour=height))+geom_point()+geom_smooth()
ggplot(bio_data$forest_short,aes(beech,ndvi))+geom_point()+facet_wrap(~sub_basin)+geom_smooth()
ggplot(bio_data$forest_short,aes(sub_basin,ndvi))+geom_boxplot(alpha=.5)+geom_point(aes(color=height),alpha=.7,size=1.5,position=position_jitter(width=.25,height=0))
ggplot()+geom_point(data=bio_data$forest_short,aes(sub_basin,ndvi))
ggplot()+geom_point(data=bio_data$forest_short,aes(sub_basin,ndvi,color=height))

'themes and storing the plot'
a <- ggplot()+geom_point(data=bio_data$forest_short,aes(sub_basin,ndvi,color=height))
a+theme_bw()    #call a stored plot and add new options
theme_set(theme_bw())  #for all further plots
a         #global settings
a <- ggplot()+geom_point(data=mpg,aes(x=displ,y=hwy,colour=class)) #easily change the theme
a+theme_bw() #change the theme for one plot
theme_set(theme_bw()) #change the theme globally

'define your theme'
ggplot()+geom_point(data=mpg,aes(x=displ,y=hwy,colour=class))+facet_grid(manufacturer~class)+ggtitle('EAGLE chart')+theme(plot.title=element_text(angle=0,size=22,colour="hotpink"))+scale_color_discrete(name="type")


######Task for December 5th 2017#########
#plot some other data

install.packages("RCurl")
library(RCurl)
#getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vTbXxJqjfY-voU-9UWgWsLW09z4dzWsv9c549qxvVYxYkwbZ9RhGE4wnEY89j4jzR_dZNeiWECW9LyW/pub?gid=0&single=true&output=csv")
x <- read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vTbXxJqjfY-voU-9UWgWsLW09z4dzWsv9c549qxvVYxYkwbZ9RhGE4wnEY89j4jzR_dZNeiWECW9LyW/pub?gid=0&single=true&output=csv")))
x
summary(x)

install.packages("ggplot2")
library(ggplot2)
x11()
ggplot()


#########Uebung 6 - 05.12.2017#################################################################################################
#ggplot2 playground
install.packages("ggplot2")
library(ggplot2)
data=data.frame(value=rnorm(10000))
ggplot(data,aes(x=value))+geom_histogram()      #basic histogram
ggplot(data,aes(x=value))+geom_histogram(binwidth = 0.05)  #custom binning = giving the size of thge bin
ggplot(data,aes(x=value))+geom_histogram(binwidth = 0.2,color="white",fill=rgb(0.2,0.7,0.1,0.4))    #uniform colour
#ggplot(data, aes(x=value)) + geom_histogram(binwidth = 0.2, aes(fill = ..count..) )  #proportional colour


##Daten aus CSV File plotten
install.packages("RCurl")
library(RCurl)
#getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vTbXxJqjfY-voU-9UWgWsLW09z4dzWsv9c549qxvVYxYkwbZ9RhGE4wnEY89j4jzR_dZNeiWECW9LyW/pub?gid=0&single=true&output=csv")
x <- read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vTbXxJqjfY-voU-9UWgWsLW09z4dzWsv9c549qxvVYxYkwbZ9RhGE4wnEY89j4jzR_dZNeiWECW9LyW/pub?gid=0&single=true&output=csv")))
x
summary(x)
install.packages("reshape2")
library(reshape2)
x2 <- melt(data=x)
library(ggplot2)
ggplot(x2,aes(x=variable,y=value))+geom_boxplot()

x.cs <- data.frame(variable=names(x),cs=t(cumsum(x)[nrow(x),])) #cumulative sum of missed minutes is in last row and we index that
names(x.cs) <- c("variable","cumsum")    #change the names to fit the melt output and to be able to merge it later on
x2 <- melt(data=x)     #reshaping the data to see the differences
x3 <- merge(x.cs,x2,by.x="variable",all=T)   #merge the two data frames based on "variable"
ggplot(x3,aes(x=variable,y=value,colour=cumsum))+geom_point()    #plot the sum as colour
ggplot(x3,aes(x=variable,y=value,colour=cumsum))+geom_boxplot(alpha=.5)+geom_point(alpha=.7,size=1.5, position=position_jitter(width=0.25,height=.5))

install.packages("gender")
library(gender)
x.g <- gender(names(x))  #run the gender detection on names
colnames(x.g)[1] <- "variable"  #change the column name again for later merging
x4 <- merge(x3,x.g,by.x="variable",all=T) #merging w/ previously created data
ggplot(x4,aes(x=variable,y=value,colour=cumsum))+geom_boxplot()+facet_wrap(~gender)


#your theoretical field data format
fielddata_wide <- read.table(header=TRUE,text='
plot_id name Cover LAI DBH
1 Sophie 7.9 12.3 10.7
2 Achmed 6.3 10.6 11.1
3 Achmed 9.5 13.1 13.8
4 Sophie 11.5 13.4 12.9
')
library(reshape2)
fielddata_long <- melt(fielddata_wide, id.vars=c("plot_id","name"),measure.vars=c("Cover","LAI","DBH"),variable.name="method",value.name="measurement")
fielddata_long

data_wide <- dcast(fielddata_long,plot_id+name~sample,value.var="measurement")  #doesnt work :(

##ggplot2 and spatial data
install.packages("ggmap")
library(ggmap)
install.packages("mapproj")
library(mapproj)
map.do <- get_map("Dortmund",zoom=15)
map.wue <- get_map("Wurzburg")
map <- get_map("Bavaria",zoom=6)
ggmap(map)
ggmap(map.do)


lsat.df <- data.frame(coordinates(lsat),getValues(lsat))
lsat.df <- lsat.df[lsat.df$B1_dn!=0,]
ggplot(lsat.df)+geom_raster(aes(x=x,y=y,fill=B3_dn))+scale_fill_gradient(na.value = NA)+coord_equal()
ggplot(lsat.df)+geom_raster(aes(x=x,y=y,fill=B3_dn))+scale_fill_gradient(low="black",high="white",na.value = NA)+coord_equal()

a <- ggplot(lsat.df)+geom_raster(aes(x=x,y=y,fill=B3_dn))+scale_fill_gradient(low="black",high="white",na.value = NA)+coord_equal()
a
poly <- readRDS(system.file("external/trainingPolygons.rds",package="RStoolbox"))
plots <- as.data.frame(coordinates(poly))
a+guides(fill=guide_colourbar())+geom_point(data=plots,aes(x=V1,y=V2),shape=2,colour="orange",size=26)+theme(axis.title.x = element_blank())

#plotting raster objects using RStoolbox
install.packages("devtools")
library(devtools)

install_bitbucket("EAGLE_MSc/steigerwald",username=AnikaDo)
library(steigerwald)
data("bio_data")
head(bio_data)
ggplot(bio_data$forest_short,aes(x=beech,y=ndvi))+geom_point()


################Uebung 8 - 19.12.17######################################################################################

#need a break? play a game
install.packages("fun")
library(fun)
if(.Platform$OS.type=='windows')x11()else x11(type="Xlib")
mine_sweeper()

install.packages('sudoku')
if(.Platform$OS.type=='windows')x11()else x11(type="Xlib")
library(sudoku)
playSudoku()

if(!require(devtools)){install.packages(devtools)}
devtools::install_github("brooke-watson/BRRR")
library(BRRR)  
skrrrahh('drummaboy')
skrrrahh('snoop')
skrrrahh(41)
