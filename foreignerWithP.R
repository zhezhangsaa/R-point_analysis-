nsim= 1000
#create matrix to save simulated L function. Each row contain 1 foreigner building L function simulation. 
# What is the maximum distance for which you want to inspect the functions, typically not larger than 1/4 of the shorter side of the window
 #window size window: rectangle = [16760, 17867] x [49037, 50346] units shortest window side is 1107 meter. 1107*0.25=276

maxr<- 270
r=seq(0,maxr,length=1000)# length is number of sequence, 1000 points. points x value vary from 0 to 270. 
Fs <- matrix(nrow=nsim, ncol=length(r))


#read the DBf data from work space
myData<-read.dbf("D:/zhe/pointProcess/foreigner/helsinkiCenter.dbf", as.is=FALSE)

#process the data so give mark 1 to foreigner building and 0 to non-foreigner building
myData$ULKOKANS01 <- 1*(myData$ULKOKANS!=0)

#create ppp file from dataset
my_ppp <- ppp(myData$XKOO, myData$YKOO, window=owin(c(16760,17867),c(49037,50346)),marks=as.factor(myData$ULKOKANS01));

#select all the foreiner buildings and saved it to myData_foreigners
myData_foreigners <- myData[myData$ULKOKANS01==1,]

#create ppp object for foreigner buildings without marks
foreign_ppp <- ppp(myData_foreigners$XKOO, myData_foreigners$YKOO, window=owin(c(16760,17867),c(49037,50346))) # c(15870,21618),c(46716,51976))

#plot L function for foreinger buildings, this will be the observed function
LfuncOb<- Lest(foreign_ppp,correction="isotropic",r=r)
plot(LfuncOb)
#Creat 500 L function simulations after random labeling
for(i in 1:nsim){
          sim<- rlabel(my_ppp, labels =marks(my_ppp),permute=TRUE) 
          #select foreigner building from sim
          foreigner_sim<-split(sim)$"1"
          #create L function 
          Lfunc= Lest(foreigner_sim,correction="isotropic", r=r)
          Fs[i,]<- Lfunc$iso#each row in the matrix contain one simulation
          #plot each simulation of L function 
          
}
# min and maximum of Fs

lower <- apply(Fs, MARGIN = 2, FUN=min)
upper <- apply(Fs, MARGIN = 2, FUN=max)
theo <- apply(Fs, MARGIN = 2, FUN=mean)
#plot boundary of envelope
x<-r
plot(LfuncOb)
lines(x,lower,lty=2)
lines(x,upper,lty=2)
lines(x,theo,lty=3,col=3)

# Deviation test: Here we calculate the difference

dif <- matrix(nrow=nsim, ncol=1)
for(i in 1:nsim){
	sum <- sum(abs(Fs[i,]-theo))
	dif[i,] <- sum
}

# difference for observed data
L <- Lest(foreign_ppp, correction="isotropic", r=r)
data_dif <- sum(abs(L$iso-L$r))

#calculating p
a = 1 # this is the upstairs in page 458 p value
for(i in 1:nsim){
	if(dif[i,1]>data_dif[1]){
		a = a +1
	}
}
p <- a/(nsim+1)


