##---------------------------------------------------------------------------------------------
## CMSY and BSM analysis 
## Original Code Written by Rainer Froese, Gianpaolo Coro and Henning Winker
# ><> CMSY_ICCAT.R Version modified for SMA ICCAT Assessment by Henning Winker, Cape Town, May 2017 <><
##---------------------------------------------------------------------------------------------
library(R2jags)  # Interface with JAGS
library(coda) 
library("parallel")
library("foreach")
library("doParallel")
library("gplots")
library(fitdistrplus) # Additional package for prior plots #><>
  #-----------------------------------------
  # Some general settings
  #-----------------------------------------
  # set.seed(999) # use for comparing results between runs
  rm(list=ls(all=TRUE)) # clear previous variables etc
  options(digits=3) # displays all numbers with three significant digits as default
  graphics.off() # close graphics windows from previous sessions
  FullSchaefer <- F    # initialize variable; automatically set to TRUE if enough abundance data are available
  n.chains     <- ifelse(detectCores() > 2,3,2) # set 3 chains in JAGS if more than 2 cores are available
  ncores_for_computation=detectCores() # cores to be used for parallel processing of CMSY
  cl           <- makeCluster(ncores_for_computation)
  registerDoParallel(cl, cores = ncores_for_computation)
  
  
  setwd("C:/Work/Research/LargePelagics/ASSESSMENTS/SMA_ICCAT/CMSY_report_final")
  #-----------------------------------------
  # Required settings, File names
  #-----------------------------------------
  catch_file  <-  "O_SMA_Catch.csv" #"O_Stocks_Catch_8_Med.csv"    # name of file containing "stock", "yr", "ct", and optional "bt"
  id_file     <-  "O_SMA_ID.csv"    #"O_Stocks_ID_8_Med.csv"  # name of file containing stock-specific info and settings for the analysis
  outfile     <- paste("Out_",format(Sys.Date(),format="%B%d%Y_"),id_file,sep="") # default name for output file
  outfile.txt <- paste(outfile,".txt", sep="") 

  #----------------------------------------
  # Select stock to be analyzed
  #----------------------------------------
  stocks      <-NA
  # If the input files contain more than one stock, specify below the stock to be analyzed
  # If the line below is commented out (#), all stocks in the input file will be analyzed
  stocks <- c("SMA.NA.C1.proof","SMA.SA.C1.proof","SMA.SA.C1","SMA.SA.C2") [1:4]     #c("nep-6")     #c("cod-2224") #c("tur-kask")   #c("cod-scow") # c("MULLBAR_BA") #c("CORYHIP_SA")  #   c("Tur_BS")  #  #   "sal-2231"    #c("sol-iris") # c("cod-2224")   
  
  
  #-----------------------------------------
  # General settings for the analysis
  #-----------------------------------------
  dataUncert   <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
  sigmaR       <- 0.1 # overall process error for CMSY; SD=0.1 is the default
  n            <- 10000 # initial number of r-k pairs
  n.new        <- n # initialize n.new
  ni           <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
  nab          <- 2 # default=5; minimum number of years with abundance data to run BSM
  mgraphs      <- T # set to TRUE to produce additional graphs for management
  save.plots   <- T # set to TRUE to save graphs to JPEG files
  close.plots  <- T # set to TRUE to close on-screen plots after they are saved, to avoid "too many open devices" error in batch-processing
  write.output <- T # set to TRUE if table with results in output file is wanted; expects years 2004-2010 to be available
  force.cmsy   <- F # set to TRUE if CMSY results are to be preferred over BSM results
  select.yr    <- 2015 # option to display F, B, F/Fmsy and B/Bmsy for a certain year; default NA
  
  #----------------------------------------------
  #  FUNCTIONS
  #----------------------------------------------
  # Monte Carlo filtering with Schaefer Function
  #----------------------------------------------
  SchaeferParallelSearch<-function(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki,i, ri,int.yr.i, nstartbt, yr, end.yr, endbio, npoints, pt){
    ptm<-proc.time()
    # create vectors for viable r, k and bt
    inmemorytable <- vector()
    # parallelised for the points in the r-k space
    inmemorytable <- foreach (i = 1 : npoints, .combine='rbind', .packages='foreach', .inorder=TRUE) %dopar%{
    nsbt = length(startbt)
    VP   <- FALSE
     for(nj in 1:nsbt) {  
      # create empty vector for annual biomasses
      bt <- vector()
      j<-startbt[nj]
      # set initial biomass, including 0.1 process error to stay within bounds
      bt[1]=j*ki[i]*exp(rnorm(1,0, 0.1*sigR))  ## set biomass in first year
      # repeat test of r-k-startbt combination to allow for different random error
      for(re in 1:ni)   {
        #loop through years in catch time series
        for (t in 1:nyr)  {  # for all years in the time series
          xt=rnorm(1,0, sigR) # set new process error for every year  
          zlog.sd = sqrt(log(1+(duncert)^2))
          zt=rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
          # calculate biomass as function of previous year's biomass plus surplus production minus catch
          bt[t+1] <- ifelse(bt[t]/ki[i] >= 0.01,
                            bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt,
                              bt[t]+(4*bt[t]/ki[i])*ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt) # assuming reduced r at B/k < 0.01
  
          # if biomass < 0.01 k, discard r-k-startbt combination
          if(bt[t+1] < 0.01*ki[i]) { 
            break
          } # stop looping through years, go to next upper level
          # intermediate year check
          if ((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))) { 
            break 
          }  
        } # end of loop of years
          # if loop was broken or last biomass falls outside of expected ranges 
          # do not store results, go directly to next startbt
          if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i]) ) { 
              next 
            } else {
              #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
              inmemorytablerow<-c(i,j,ri[i],ki[i],bt[1:(nyr+1)]/ki[i])
              if (length(inmemorytablerow)==(4+nyr+1)){
               if (VP==FALSE)
                {
                 inmemorytable<-inmemorytablerow
                }
               else
                {
                    inmemorytable<-rbind(inmemorytable,inmemorytablerow)
                 }
               VP<-TRUE
              }
            }
        } # end of repetition for random error
       } # end of j-loop of initial biomasses 
    # instruction necessary to make the foreach loop see the variable:
    if (length(inmemorytable)==0)
        {inmemorytable<-vector(length=4+nyr+1)*NA}
      else
        {inmemorytable}
    }#end loop on points
  
    #create the output matrix
    mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1) 
    npointsmem = dim(inmemorytable)[1]
    npointscols = dim(inmemorytable)[2]
    #reconstruction of the processing matrix after the parallel search
    if (npointsmem>0 && npointscols>0){
     for (idxr in 1:npointsmem){
          i = inmemorytable[idxr,1]
          if (!is.na(i)){
            j = inmemorytable[idxr,2]
            mdatindex<-((i-1)*nstartbt)+which(startbt==j)  
            mdat[mdatindex,1]           <- inmemorytable[idxr,3]
            mdat[mdatindex,2]           <- inmemorytable[idxr,4]
            mdat[mdatindex,3:(2+nyr+1)] <- inmemorytable[idxr,5:(4+nyr+1)]
            if(pt==T) points(x=ri[i], y=ki[i], pch=".", cex=4, col="gray")
          }
      }
    }
    ptm<-proc.time()-ptm
    mdat <- na.omit(mdat)
    return(mdat)
  }
  
  SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni) {
    
    # create vector for initial biomasses
    startbt     <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
    nstartbt    <- length(startbt)
    npoints     <- length(ri)
    # get index of intermediate year
    int.yr.i     <- which(yr==int.yr) 
    
    #loop through r-k pairs with parallel search
    mdat<-SchaeferParallelSearch(ni, nyr,sigR,duncert,ct,int.yr,intbio, startbt, ki, i, ri, int.yr.i, nstartbt, yr, end.yr, endbio, npoints,pt)
  
    cat("\n")
    return(list(mdat))
  } # end of SchaeferMC function
  
  #-----------------------------------------------
  # Function for moving average
  #-----------------------------------------------
  ma    <- function(x){
    x.1    <-   filter(x,rep(1/3,3),sides=1)
    x.1[1] <- x[1]
    x.1[2] <- (x[1]+x[2])/2
    return(x.1)
  }
  #---------------------------------------------
  # END OF FUNCTIONS
  #---------------------------------------------
  
  #-----------------------------------------
  # Start output to screen
  #-----------------------------------------
  cat("-------------------------------------------\n")
  cat("CMSY Analysis,", date(),"\n")
  cat("-------------------------------------------\n")
  
  #------------------------------------------
  # Read data and assign to vectors
  #------------------------------------------
  # create headers for data table file
  if(write.output==T){
    outheaders = data.frame("Region", "Subregion","Name","SciName","Stock","start.yr","end.yr","btype",
                            "MaxCatch","LastCatch","MSY_BSM","lcl","ucl","r_BSM","lcl","ucl",
                            "k_BSM","lcl","ucl","q_BSM","lcl","ucl","rel_B_BSM","lcl","ucl","rel_F_BSM",
                            "r_CMSY","lcl","ucl","k_CMSY","lcl","ucl","MSY_CMSY","lcl","ucl",
                            "rel_B_CMSY","2.5th","97.5th","rel_F_CMSY",
                            "F_msy","lcl","ucl","curF_msy","lcl","ucl",
                            "MSY","lcl","ucl","Bmsy","lcl","ucl",
                            "B","lcl","ucl","B_Bmsy","lcl","ucl",
                            "F","lcl","ucl","F_Fmsy","lcl","ucl",
                            "sel_B","sel_B_Bmsy","sel_F","sel_F_Fmsy",
                            "c00","c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13","c14","c15",
                            "F.Fmsy00","F.Fmsy01","F.Fmsy02","F.Fmsy03","F.Fmsy04","F.Fmsy05","F.Fmsy06","F.Fmsy07","F.Fmsy08","F.Fmsy09","F.Fmsy10","F.Fmsy11","F.Fmsy12","F.Fmsy13","F.Fmsy14","F.Fmsy15",
                            "B00","B01","B02","B03","B04","B05","B06","B07","B08","B09","B10","B11","B12","B13","B14","B15")
   
     write.table(outheaders,file=outfile, append = T, sep=",",row.names=F,col.names=F)
  }
  
  cat("Parallel processing will use",ncores_for_computation,"cores\n")
  # Read data
  cdat         <- read.csv(catch_file, header=T, dec=".", stringsAsFactors = FALSE)
  cinfo        <- read.csv(id_file, header=T, dec=".", stringsAsFactors = FALSE)
  cat("Files", catch_file, ",", id_file, "read successfully","\n")
  
  #---------------------------------
  # Analyze stock(s)
  #---------------------------------
  if(is.na(stocks[1])==TRUE){
     #stocks         <- as.character(cinfo$Stock) # Analyze stocks in sequence of ID file
    # stocks         <- sort(as.character(cinfo$Stock)) # Analyze stocks in alphabetic order
    # stocks         <- as.character(cinfo$Stock[cinfo$Subregion=="Baltic Sea"]) # Analyze stocks in Region
  }
  
  # analyze one stock after the other
  for(stock in stocks) {
    cat("Processing",stock,",", as.character(cinfo$ScientificName[cinfo$Stock==stock]),"\n")
    # assign data from cinfo to vectors
    res          <- as.character(cinfo$Resilience[cinfo$Stock==stock])
    start.yr     <- as.numeric(cinfo$StartYear[cinfo$Stock==stock])
    end.yr       <- as.numeric(cinfo$EndYear[cinfo$Stock==stock])
    r.low        <- as.numeric(cinfo$r.low[cinfo$Stock==stock])
    r.hi         <- as.numeric(cinfo$r.hi[cinfo$Stock==stock])
    user.log.r   <- ifelse(is.na(r.low)==F & is.na(r.hi)==F,TRUE,FALSE)     
    stb.low      <- as.numeric(cinfo$stb.low[cinfo$Stock==stock])
    stb.hi       <- as.numeric(cinfo$stb.hi[cinfo$Stock==stock])
    int.yr       <- as.numeric(cinfo$int.yr[cinfo$Stock==stock])
    intb.low     <- as.numeric(cinfo$intb.low[cinfo$Stock==stock])
    intb.hi      <- as.numeric(cinfo$intb.hi[cinfo$Stock==stock])
    endb.low     <- as.numeric(cinfo$endb.low[cinfo$Stock==stock])
    endb.hi      <- as.numeric(cinfo$endb.hi[cinfo$Stock==stock])
    btype        <- as.character(cinfo$btype[cinfo$Stock==stock])
    force.cmsy   <- ifelse(force.cmsy==T,T,cinfo$force.cmsy[cinfo$Stock==stock])
    comment      <- as.character(cinfo$Comment[cinfo$Stock==stock])
    # set global defaults for uncertainty
    duncert      <- dataUncert
    sigR         <- sigmaR
  
    # check for common errors
    if (length(btype)==0){
      cat("ERROR: Could not find the stock in the ID input file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)")
      return (NA) }
    if(start.yr < cdat$yr[cdat$Stock==stock][1]){
      cat("ERROR: start year in ID file before first year in catch file\n")
      return (NA)}
    
    # extract data on stock
    yr           <- as.numeric(cdat$yr[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])
    
    if (length(yr)==0){
      cat("ERROR: Could not find the stock in the Catch input files - Please check that the code is written correctly")
      return (NA)
    }
    
    ct.raw           <- as.numeric(cdat$ct[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that catch is given in tonnes, transforms to '000 tonnes
    if(btype=="biomass" | btype=="CPUE" ) {
      bt <- as.numeric(cdat$bt[cdat$Stock==stock & cdat$yr >= start.yr & cdat$yr <= end.yr])/1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
    } else {bt <- NA}
    
    if(is.na(mean(ct.raw))){
    cat("ERROR: Missing value in Catch data; fill or interpolate\n")  
    }
    nyr          <- length(yr) # number of years in the time series
    
    # change catch to 3 years moving average where value is average of past 3 years 
    ct              <- ma(ct.raw)
  
    # initialize vectors for viable r, k, bt, and all in a matrix
    mdat.all    <- matrix(data=vector(),ncol=2+nyr+1)
    
    # initialize other vectors anew for each stock
    current.attempts <- NA
    
    #----------------------------------------------------
    # Determine initial ranges for parameters and biomass
    #----------------------------------------------------
    # initial range of r from input file
    if(is.na(r.low)==F & is.na(r.hi)==F) {
      start.r <- c(r.low,r.hi)
    } else {
      # initial range of r based on resilience
      if(res == "High") {
        start.r <- c(0.6,1.5)} else if(res == "Medium") {
          start.r <- c(0.2,0.8)}    else if(res == "Low") {
            start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
              start.r <- c(0.015,0.1)} 
    }
    
    # get index of years with lowest and highest catch between start+3 and end-3 years
    min.yr.i     <- which.min(ct[4:(length(ct)-3)])+3
    max.yr.i     <- which.max(ct[4:(length(ct)-3)])+3
    min.ct       <- ct[min.yr.i]
    max.ct       <- ct[max.yr.i]
    
    # use initial biomass range from input file if stated
    if(is.na(stb.low)==F & is.na(stb.hi)==F) {
      startbio <- c(stb.low,stb.hi)
    } else {
      # if start year < 1960 assume high biomass
      if(start.yr < 1960) {startbio <- c(0.5,0.9)} else {
        # else use medium prior biomass range
        startbio <- c(0.2,0.6)} }
    
    # use year and biomass range for intermediate biomass from input file
    if(is.na(intb.low)==F & is.na(intb.hi)==F) {
      int.yr   <- int.yr
      intbio   <- c(intb.low,intb.hi)
      
      # if contrast in catch is low, use initial range again in mid-year
    } else if(min(ct)/max(ct) > 0.6) {
      int.yr    <- as.integer(mean(c(start.yr, end.yr)))
      intbio    <- startbio 
      
      # else if year of minimum catch is after max catch then use min catch
    } else if(min.yr.i > max.yr.i) {
      int.yr    <- yr[min.yr.i-1]
      if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) & 
           (min.ct/max.ct) > 0.3) intbio <- c(0.2,0.6) else intbio <- c(0.01,0.4)
      
      # else use max catch  
    } else {
      # assume that biomass range in year before maximum catch was high or medium
      int.yr    <- yr[max.yr.i-1]
      intbio    <- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                        # ((min.ct/max.ct < 0.3 & (max.yr.i - min.yr.i) < 25))) c(0.5,0.9) else c(0.2,0.6) }
                        (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else c(0.2,0.6) } # if incease is steep, assume high, else medium
    # end of intbio setting
    
    # final biomass range from input file
    if(is.na(endb.low)==F & is.na(endb.hi)==F) {
      endbio   <- c(endb.low,endb.hi)
    } else {
      # else use mean final catch/max catch to estimate final biomass
      endbio  <- if(ct[nyr]/max(ct) > 0.7) {c(0.5,0.9)} else if(ct[nyr]/max(ct) < 0.3) {c(0.01,0.4)} else {c(0.2,0.6)}
    } # end of final biomass setting
    
    # initial prior range of k values, assuming min k will be larger than max catch / prior for r 
    if(mean(endbio) <= 0.5) {
      start.k <- c(max(ct)/start.r[2],4*max(ct)/start.r[1])} else {
        start.k <- c(2*max(ct)/start.r[2],12*max(ct)/start.r[1])} 
    # start.k <- c(start.k[1],3000)   
    cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
        ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
        ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
    
    #------------------------------------------------------------------
    # Uniform sampling of the r-k space
    #------------------------------------------------------------------
    # get random set of r and k from log space distribution 
    ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))  
    ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))  
    
    #-----------------------------------------------------------------
    # Plot data and progress
    #-----------------------------------------------------------------
    # check for operating system, open separate window for graphs if Windows
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(10,5)}
    par(mfrow=c(1,2))
    # plot catch
    plot(x=yr, y=ct.raw, 
         ylim=c(0,max(ifelse(substr(id_file,1,3)=="Sim",
                             1.1*true.MSY,0),1.2*max(ct.raw))),
         type ="l", bty="l", main=paste("A: ",stock,"catch"), xlab="Year", ylab="Catch", lwd=2)
    lines(x=yr,y=ct,col="blue", lwd=1)
    points(x=yr[max.yr.i], y=max.ct, col="red", lwd=2)
    points(x=yr[min.yr.i], y=min.ct, col="red", lwd=2)
  
    # plot r-k graph
    plot(x=ri1, y=ki1, xlim = start.r, ylim = start.k, log="xy", xlab="r", ylab="k", 
         main="B: Finding viable r-k", pch=".", cex=3, bty="l", col="gray95")
  
    #---------------------------------------------------------------------
    # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
    #---------------------------------------------------------------------
    cat("First Monte Carlo filtering of r-k space with ",n," points...\n")
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                       pt=T, duncert=dataUncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    # count viable trajectories and r-k pairs 
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
    cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
    
    #----------------------------------------------------------------------- 
    # 2 - if the lower bound of k is too high, reduce it by half and rerun
    #-----------------------------------------------------------------------
    if(length(kv.all[kv.all < 1.1*start.k[1] & rv.all < mean(start.r)]) > 10) {
      cat("Reducing lower bound of k, resampling area with",n,"additional points...\n")
      start.k <- c(0.5*start.k[1],start.k[2])
      ri1 = exp(runif(n, log(start.r[1]), log(start.r[2])))  
      ki1 = exp(runif(n, log(start.k[1]), log(start.k[2])))  
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=T, duncert=dataUncert, startbins=10, ni=ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
    }
    
    #-------------------------------------------------------------------
    # 3 - if few points were found then resample and shrink the log k space
    #-------------------------------------------------------------------
    if (n.viable.b <= 1000){
      log.start.k.new  <- log(start.k) 
      max.attempts     <- 3
      current.attempts <- 1
      startbins        <- 10  
      while (n.viable.b <= 1000 && current.attempts <= max.attempts){
        if(n.viable.pt > 0) {
          log.start.k.new[1] <- mean(c(log(start.k[1]), min(log(kv.all))))
          log.start.k.new[2] <- mean(c(log.start.k.new[2], max(log(kv.all)))) }
        n.new <- n*current.attempts #add more points
        ri1 = exp(runif(n.new, log(start.r[1]), log(start.r[2])))  
        ki1 = exp(runif(n.new, log.start.k.new[1], log.start.k.new[2]))
        cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start.k.new[1]),",",exp(log.start.k.new[2]),"]\n")
        cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points...","\n")
        if(current.attempts==2 & n.viable.b < 50){
          duncert   <- 2*dataUncert
          sigR      <- 2*sigmaR
          startbins <- 20
          cat("Doubling startbins, catch and process error, and number of variability patterns \n")   
        }
        MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                           pt=T, duncert=duncert, startbins=startbins, ni=2*ni)
        mdat.all <- rbind(mdat.all,MCA[[1]])
        rv.all   <- mdat.all[,1]
        kv.all   <- mdat.all[,2]
        btv.all  <- mdat.all[,3:(2+nyr+1)]
        n.viable.b   <- length(mdat.all[,1])
        n.viable.pt <- length(unique(mdat.all[,1]))
        cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
        current.attempts=current.attempts+1 #increment the number of attempts
      }
      if(n.viable.b < 5) {
        cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")
        next
      }  
    }
    
    #------------------------------------------------------------------
    # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
    #------------------------------------------------------------------
    if(length(rv.all[rv.all > 0.9*start.r[2]]) < 5) { 
      l.sample.r        <- quantile(rv.all,0.6)
      add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
      cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points...\n")
      log.start.k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))
      
      ri1 = exp(runif(add.points, log(l.sample.r), log(start.r[2])))  
      ki1 = exp(runif(add.points, log.start.k.new[1], log.start.k.new[2]))
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=T, duncert=duncert, startbins=10, ni=ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
    }
  
    # ------------------------------------------------------------------
    # Bayesian analysis of catch & biomass (or CPUE) with Schaefer model
    # ------------------------------------------------------------------
    FullSchaefer <- F
    if(btype != "None" & length(bt[is.na(bt)==F])>=nab) {
      FullSchaefer <- T
   
      # set inits for r-k in lower right corner of log r-k space
      init.r      <- start.r[1]+0.8*(start.r[2]-start.r[1])
      init.k      <- start.k[1]+0.1*(start.k[2]-start.k[1])
      
      # vector with no penalty (=0) if predicted biomass is within viable range, else a penalty of 10 is set 
      pen.bk = pen.F = rep(0,length(ct))
      
      #><>><>
      # Add biomass priors
      b.yrs = c(1,length(start.yr:int.yr),length(start.yr:end.yr))
      b.prior = rbind(matrix(c(startbio[1],startbio[2],intbio[1],intbio[2],endbio[1],endbio[2]),2,3),rep(0,3)) # last row includes the 0 pen
      #><>><>
  
      cat("Running MCMC analysis....\n")
      if(btype == "biomass" ) {
        # Data to be passed on to JAGS
        jags.data        <- c('ct','bt','nyr', 'start.r','startbio','start.k',
                              'init.r','init.k', 'pen.bk','pen.F','b.yrs','b.prior')
        # Parameters to be returned by JAGS
        jags.save.params <- c('r','k','P','AI') # 
        
        # JAGS model
        Model = "model{
  	# to avoid crash due to 0 values
  	eps<- 0.01
  	penm[1] <- 0 # no penalty for first biomass
  	Pmean[1] <- log(alpha)
  	P[1] ~ dlnorm(Pmean[1],itau2)
  
  	for (t in 2:nyr) {
        Pmean[t] <- ifelse(P[t-1] > 0.01,
        log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps)),  # Process equation
        log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps))) # assuming reduced r at B/k < 0.01
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(k*P[t])-log(k*(eps+0.001)),ifelse(P[t]>1,log(k*P[t])-log(k*(0.99)),0)) # penalty if Pmean is outside viable biomass
        
      }
        
        # ><> Biomass priors/penalties are enforced as follows
        for (i in 1:3) {
        penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(k*P[b.yrs[i]])-log(k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(k*P[b.yrs[i]])-log(k*b.prior[2,i]),0)) 
        b.prior[3,i] ~ dnorm(penb[i],100)
        }
        
        for (t in 1:nyr){
        Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) #><> Penalty term on F > 1, i.e. ct>B
        pen.F[t]  ~ dnorm(Fpen[t],1000)
        pen.bk[t] ~ dnorm(penm[t],10000) 
        Bm[t] <- log(P[t]*k);
        bt[t]    ~ dlnorm(Bm[t],(isigma2+100));
        AI[t] <- P[t]*k #><> abundance index
        }
  
  	# priors
  	# search in the alpha space from the center of the range. Allow high variability
  	log.alpha               <- log((startbio[1]+startbio[2])/2)
  	sd.log.alpha            <- (log.alpha-log(startbio[1]))/5
  	tau.log.alpha           <- pow(sd.log.alpha,-2)
  	alpha                   ~  dlnorm(log.alpha,tau.log.alpha)
  
  	# search in the k space from 20% of the range
  	log.km              <- log(start.k[1]+0.2*(start.k[2]-start.k[1]))
  	sd.log.k            <- (log.km-log(start.k[1]))/4
  	tau.log.k           <- pow(sd.log.k,-2)
  	k                   ~  dlnorm(log.km,tau.log.k)
  
  	# define process (tau) and observation (sigma) variances as inversegamma priors
  	itau2 ~ dgamma(2,0.01)
  	tau2  <- 1/itau2
  	tau   <- pow(tau2,0.5)
  
  	isigma2 ~ dgamma(2,0.01)
  	sigma2 <- 1/isigma2
  	sigma <- pow(sigma2,0.5) 
  
    log.rm              <- mean(log(start.r))
    sigma.log.r         <- abs(log.rm - log(start.r[1]))/2
    tau.log.r           <- pow(sigma.log.r,-2)
    r                   ~  dlnorm(log.rm,tau.log.r)   
    } " # end of JAGS model for btype=="biomass"     
    
    # ---------------------------------------------------------------------
    # Schaefer model for Catch & CPUE
    # ---------------------------------------------------------------------
    } else {  
      # get prior for q from stable catch/biomass period, min 5 years; get range of years from input file
      q.start    <- cinfo$q.start[cinfo$Stock==stock]
      q.end      <- cinfo$q.end[cinfo$Stock==stock]
      if(is.na(q.start)==F & is.na(q.end)==F) {
        mean.last.ct      <-mean(ct[yr >= q.start & yr <= q.end], na.rm=T) # get mean catch of indicated years
        mean.last.cpue    <-mean(bt[yr >= q.start & yr <= q.end], na.rm=T) # get mean of CPUE of indicated years
      } else { 
        # get prior range for q from mean catch and mean CPUE in recent years 
        lyr               <- ifelse(mean(start.r)>=0.5,5,10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish    
        mean.last.ct      <-mean(ct[(nyr-lyr):nyr],na.rm=T) # get mean catch of last years
        mean.last.cpue    <-mean(bt[(nyr-lyr):nyr],na.rm=T) # get mean of CPUE of last years
        
      }
      gm.start.r      <- exp(mean(log(start.r))) # get geometric mean of prior r range
      if(mean(endbio) >= 0.5) {  # if biomass is high  
        q.1           <- mean.last.cpue*0.25*gm.start.r/mean.last.ct
        q.2           <- mean.last.cpue*0.5*start.r[2]/mean.last.ct
      } else {
        q.1           <- mean.last.cpue*0.5*gm.start.r/mean.last.ct
        q.2           <- mean.last.cpue*start.r[2]/mean.last.ct
      } 
      q.prior         <- c(q.1,q.2)
      init.q          <- mean(q.prior)
      
      # note priors
      Prior.range <- data.frame(Min=c(start.k[1],start.r[1],q.prior[1],b.prior[1,]),Max=round(c(start.k[2],start.r[2],q.prior[2],b.prior[2,]),4))
      row.names(Prior.range) <- c("k","r","q","B.init","B.int","B.end")
      write.csv(Prior.range,paste0("Prior_",stock,".csv"))
      
      # Data to be passed on to JAGS
      jags.data        <- c('ct','bt','nyr', 'start.r', 'start.k', 'startbio', 'q.prior',
                            'init.q','init.r','init.k','pen.bk','pen.F','b.yrs','b.prior')
      # Parameters to be returned by JAGS
      jags.save.params <- c('r','k','q', 'P', "AI","beta") 
      
      # JAGS model
      Model = "model{
      # to reduce chance of non-convergence, Pmean[t] values are forced >= eps
      eps<-0.01
      penm[1] <- 0 # no penalty for first biomass
      Pmean[1] <- log(alpha)
      P[1] ~ dlnorm(Pmean[1],itau2)
  
      for (t in 2:nyr) {
        Pmean[t] <- ifelse(P[t-1] > 0.01,
        log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps)),  # Process equation
        log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
        P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
        penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*k*P[t])-log(q*k*(eps+0.001)),ifelse(P[t]>1,log(q*k*P[t])-log(q*k*(0.99)),0)) # penalty if Pmean is outside viable biomass
      }
      
      # ><> Biomass priors/penalties are enforced as follows
      for (i in 1:3) {
        penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[2,i]),0)) 
        b.prior[3,i] ~ dnorm(penb[i],100)
      }
      
      for (t in 1:nyr){
        Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) #><> Penalty term on F > 1, i.e. ct>B
        pen.F[t]  ~ dnorm(Fpen[t],1000)
        pen.bk[t] ~ dnorm(penm[t],10000) 
        cpuem[t]  <- log(q*P[t]*k);
        bt[t]     ~ dlnorm(cpuem[t],(isigma2));
        AI[t] <- q*P[t]*k #><> abundance index
     }
      
    # priors
    log.alpha               <- log((startbio[1]+startbio[2])/2) # needed for fit of first biomass
    sd.log.alpha            <- (log.alpha-log(startbio[1]))/4
    tau.log.alpha           <- pow(sd.log.alpha,-2)
    alpha                   ~  dlnorm(log.alpha,tau.log.alpha)
  
    # search in the k space starting from 20% of the range
    log.km              <- log(start.k[1]+0.2*(start.k[2]-start.k[1]))
    sd.log.k            <- (log.km-log(start.k[1]))/4
    tau.log.k           <- pow(sd.log.k,-2)
    k                   ~  dlnorm(log.km,tau.log.k)
      
    # set realistic prior for q
    log.qm              <- mean(log(q.prior))
    sd.log.q            <- (log.qm-log(q.prior[1]))/4
    tau.log.q           <- pow(sd.log.q,-2)
    q                   ~  dlnorm(log.qm,tau.log.q)
      
    # define process (tau) and observation (sigma) variances as inversegamma prios
    itau2 ~ dgamma(4,0.01)
    tau2  <- 1/itau2
    tau   <- pow(tau2,0.5)
      
    isigma2 ~ dgamma(2,0.01)
    sigma2 <- 1/isigma2+0.225
    sigma <- pow(sigma2,0.5)
  
    log.rm              <- mean(log(start.r))
    sigma.log.r         <- abs(log.rm - log(start.r[1]))/2
    tau.log.r           <- pow(sigma.log.r,-2)
    r                   ~  dlnorm(log.rm,tau.log.r)   
   } "    # end of JAGS model for CPUE                  
      
  } # end of else loop for Schaefer with CPUE
  
  # Write JAGS model to file
  cat(Model, file="r2jags.bug")  
  
  if(btype=="biomass") {
    j.inits     <- function(){list("r"=rnorm(1,mean=init.r,sd=0.2*init.r),
                                   "k"=rnorm(1,mean=init.k,sd=0.1*init.k),
                                   "itau2"=1000,
                                   "isigma2"=1000)}} else {
    j.inits <- function(){list("r"=rnorm(1,mean=init.r,sd=0.2*init.r),
                                    "k"=rnorm(1,mean=init.k,sd=0.1*init.k),
                                    "q"=rnorm(1,mean=init.q,sd=0.2*init.q),
                                    "itau2"=1000,
                                    "isigma2"=1000)}}
  # run model
  jags_outputs <- jags.parallel(data=jags.data, 
                                working.directory=NULL, inits=j.inits, 
                                parameters.to.save=jags.save.params, 
                                model.file="r2jags.bug", n.chains = n.chains, 
                                n.burnin = 30000, n.thin = 10, 
                                n.iter = 60000)
  
  # ------------------------------------------------------
  # Results from JAGS Schaefer
  # ------------------------------------------------------
  r_raw            <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$r))
  k_raw            <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$k))
  # Importance sampling: only accept r-k pairs where r is near the prior range
  r_out            <- r_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]   
  k_out            <- k_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]   
  
  mean.log.r.jags  <- mean(log(r_out))
  sd.log.r.jags    <- sd(log(r_out))
  r.jags           <- exp(mean.log.r.jags)
  lcl.r.jags       <- exp(mean.log.r.jags - 1.96*sd.log.r.jags)
  ucl.r.jags       <- exp(mean.log.r.jags + 1.96*sd.log.r.jags)
  mean.log.k.jags  <- mean(log(k_out))
  sd.log.k.jags    <- sd(log(k_out))
  k.jags           <- exp(mean.log.k.jags)
  lcl.k.jags       <- exp(mean.log.k.jags - 1.96*sd.log.k.jags)
  ucl.k.jags       <- exp(mean.log.k.jags + 1.96*sd.log.k.jags)
  MSY.posterior     <- r_out*k_out/4 # simpler
  mean.log.MSY.jags <- mean(log(MSY.posterior))
  sd.log.MSY.jags   <- sd(log(MSY.posterior))
  MSY.jags          <- exp(mean.log.MSY.jags)
  lcl.MSY.jags      <- exp(mean.log.MSY.jags - 1.96*sd.log.MSY.jags)
  ucl.MSY.jags      <- exp(mean.log.MSY.jags + 1.96*sd.log.MSY.jags)
  
  if(btype=="CPUE") {
    q_out           <- as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$q))
    mean.log.q      <- mean(log(q_out))
    sd.log.q        <- sd(log(q_out))
    mean.q          <- exp(mean.log.q)
    lcl.q           <- exp(mean.log.q-1.96*sd.log.q)
    ucl.q           <- exp(mean.log.q+1.96*sd.log.q)
    F.bt.cpue       <- mean.q*ct.raw/bt
    Fmsy.cpue       <- r.jags/2
  }
  
  # get F from observed biomass
  if(btype == "biomass") {
    F.bt       <- ct.raw/bt
    Fmsy.bt    <- r.jags/2
   }
  
  # get relative biomass P=B/k as predicted by BSM, including predictions for years with NA abundance
  all.P    <- jags_outputs$BUGSoutput$sims.list$P # matrix with P distribution by year
  quant.P  <- apply(all.P,2,quantile,c(0.025,0.5,0.975),na.rm=T)
  
  # get k, r posterior ><>
  all.k  <- jags_outputs$BUGSoutput$sims.list$k # matrix with P distribution by year
  all.r  <- jags_outputs$BUGSoutput$sims.list$r # matrix with P distribution by year
  quant.k <- apply(all.k,2,quantile,c(0.025,0.5,0.975),na.rm=T) 
  quant.r <- apply(all.r,2,quantile,c(0.025,0.5,0.975),na.rm=T) 
   
  # get B/Bmys posterior 
  all.b_bmsy=NULL
  for(t in 1:ncol(all.P)){
    all.b_bmsy  <- cbind(all.b_bmsy,all.P[,t]*2)}
  
  # get F/Fmys posterior ><>
  all.F_Fmsy=NULL
  for(t in 1:ncol(all.P)){
    all.F_Fmsy<- cbind(all.F_Fmsy,(ct.raw[t]/(all.P[,t]*all.k))/ifelse(all.P[,t]>0.01,all.r/2,all.r/2*4*all.P[,t]))}
  
  quant.F_Fmsy <- apply(all.F_Fmsy,2,quantile,c(0.025,0.5,0.975),na.rm=T) 
  quant.b_bmsy <- apply(all.b_bmsy,2,quantile,c(0.025,0.5,0.975),na.rm=T) 
  
  all.AI <-  jags_outputs$BUGSoutput$sims.list$AI
  quant.AI <- apply(all.AI,2,quantile,c(0.025,0.5,0.975))
  
  } # end of MCMC Schaefer loop 
  
    
   
  #------------------------------------
  # get results from CMSY
  #------------------------------------
  # get estimate of most probable r as 75th percentile of mid log.r-classes
  # get unique combinations of r-k
  unique.rk         <- unique(mdat.all[,1:2])
  # get remaining viable log.r and log.k 
  log.rs           <- log(unique.rk[,1])
  log.ks           <- log(unique.rk[,2])
  # get vectors with numbers of r and mid values in classes
  # determine number of classes as a function of r-width
  r.width         <- (max(unique.rk[,1])-start.r[1])/(start.r[2]-start.r[1])
  classes         <- ifelse(r.width>0.8,100,ifelse(r.width>0.5,50,ifelse(r.width>0.3,25,12)))
  hist.log.r      <- hist(x=log.rs, breaks=classes, plot=F)
  log.r.counts    <- hist.log.r$counts
  log.r.mids      <- hist.log.r$mids
  # get most probable log.r as 75th percentile of mids with counts > 0
  log.r.est       <- as.numeric(quantile(log.r.mids[which(log.r.counts > 0)],0.75))
  median.log.r    <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.50))
  lcl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.5125))
  ucl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.9875))
  sd.log.r.est    <- (ucl.log.r - log.r.est) / 1.96 
  r.est           <- exp(log.r.est)
  lcl.r.est       <- exp(log.r.est-1.96*sd.log.r.est)
  ucl.r.est       <- exp(log.r.est+1.96*sd.log.r.est)
  
  # get r-k pairs above median of mids
  rem            <- which(unique.rk[,1] > exp(median.log.r))
  rem.log.r      <- log(unique.rk[,1][rem])
  rem.log.k      <- log(unique.rk[,2][rem])
  # do linear regression of log k ~ log r with slope fixed to -1 (from Schaefer)
  reg            <- lm(rem.log.k ~ 1 + offset(-1*rem.log.r))
  int.reg        <- as.numeric(reg[1])
  sd.reg      <- sd(resid(reg))
  # get estimate of log(k) from y where x = log.r.est
  log.k.est      <- int.reg + (-1) * log.r.est
  # get estimates of ucl of log.k.est from y + SD where x = ucl.log.r  
  ucl.log.k     <- int.reg + (-1) * lcl.log.r + sd.reg
  # get estimates of sd.log.k.est from upper confidence limit of log.k.est  
  sd.log.k.est   <- (ucl.log.k - log.k.est) / 1.96 
  lcl.log.k      <- log.k.est - 1.96*sd.log.k.est
  ucl.log.k      <- log.k.est + 1.96*sd.log.k.est
  
  k.est       <- exp(log.k.est)
  lcl.k.est   <- exp(lcl.log.k)
  ucl.k.est   <- exp(ucl.log.k)
  
  # get MSY from remaining log r-k pairs
  log.MSY.est     <- mean(rem.log.r + rem.log.k - log(4))
  sd.log.MSY.est  <- sd(rem.log.r + rem.log.k - log(4))
  lcl.log.MSY.est <- log.MSY.est - 1.96*sd.log.MSY.est 
  ucl.log.MSY.est <- log.MSY.est + 1.96*sd.log.MSY.est
  MSY.est         <- exp(log.MSY.est)
  lcl.MSY.est     <- exp(lcl.log.MSY.est)
  ucl.MSY.est     <- exp(ucl.log.MSY.est)
  
  # get predicted biomass vectors as median and quantiles 
  # only use biomass trajectories from r-k pairs within the confidence limits
  rem.btv.all <- mdat.all[which(mdat.all[,1] > lcl.r.est & mdat.all[,1] < ucl.r.est 
                                & mdat.all[,2] > lcl.k.est & mdat.all[,2] < ucl.k.est),3:(2+nyr+1)]
  median.btv <- apply(rem.btv.all,2, median)
  median.btv.lastyr  <- median.btv[length(median.btv)-1]
  nextyr.bt  <- median.btv[length(median.btv)]
  lcl.btv    <- apply(rem.btv.all,2, quantile, probs=0.025)
  q.btv      <- apply(rem.btv.all,2, quantile, probs=0.25)
  ucl.btv    <- apply(rem.btv.all,2, quantile, probs=0.975)
  lcl.median.btv.lastyr <- lcl.btv[length(lcl.btv)-1]
  ucl.median.btv.lastyr <- ucl.btv[length(lcl.btv)-1]  
  lcl.nextyr.bt <- lcl.btv[length(lcl.btv)]
  ucl.nextyr.bt <- ucl.btv[length(lcl.btv)]
  
  # get biomass for relative abundance ><>
  median.bt = median.btv*k.est
  lcl.bt = lcl.btv*k.est
  ucl.bt = ucl.btv*k.est
  
  # Normalize for plotting ><>
  norm.mbt = median.bt/mean(median.bt)
  norm.lbt = lcl.bt/mean(median.bt)
  norm.ubt = ucl.bt/mean(median.bt)
  
  # get F derived from predicted CMSY biomass
  F.CMSY      <- ct/(median.btv[1:nyr]*k.est)
  Fmsy.CMSY  <- r.est/2 # Fmsy from CMSY 
  
  # --------------------------------------------
  # Get results for management
  # --------------------------------------------
    # Original code to keep intengrity
  if(FullSchaefer==F | force.cmsy==T) { # if only CMSY is available or shall be used
    MSY   <-MSY.est; lcl.MSY<-lcl.MSY.est; ucl.MSY<-ucl.MSY.est 
    Bmsy  <-k.est/2; lcl.Bmsy<-lcl.k.est/2; ucl.Bmsy<-ucl.k.est/2
    Fmsy  <-r.est/2; lcl.Fmsy<-lcl.r.est/2; ucl.Fmsy<-ucl.r.est/2
    B.Bmsy<-2*median.btv[1:nyr];lcl.B.Bmsy<-2*lcl.btv[1:nyr];ucl.B.Bmsy<-2*ucl.btv[1:nyr]
    if(is.na(select.yr)==F){B.Bmsy.sel<-2*median.btv[yr==select.yr]}
    
  } else {
    MSY   <-MSY.jags; lcl.MSY<-lcl.MSY.jags; ucl.MSY<-ucl.MSY.jags 
    Bmsy  <-k.jags/2; lcl.Bmsy<-lcl.k.jags/2; ucl.Bmsy<-ucl.k.jags/2
    Fmsy  <-r.jags/2; lcl.Fmsy<-lcl.r.jags/2; ucl.Fmsy<-ucl.r.jags/2
    B.Bmsy<-2*quant.P[2,];lcl.B.Bmsy<-2*quant.P[1,];ucl.B.Bmsy<-2*quant.P[3,]
    if(is.na(select.yr)==F) {B.Bmsy.sel<-2*quant.P[2,][yr==select.yr]}
  }
  B          <-B.Bmsy*Bmsy;lcl.B<-lcl.B.Bmsy*Bmsy;ucl.B<-ucl.B.Bmsy*Bmsy
  B.last     <-B[nyr];lcl.B.last<-lcl.B[nyr];ucl.B.last<-ucl.B[nyr]
  B.Bmsy.last<-B.Bmsy[nyr];lcl.B.Bmsy.last<-lcl.B.Bmsy[nyr];ucl.B.Bmsy.last<-ucl.B.Bmsy[nyr]
  
  Fm           <- ct.raw/B;lcl.F<-ct.raw/ucl.B;ucl.F<-ct.raw/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  F.Fmsy       <- Fm/Fmsy.vec; lcl.F.Fmsy<-lcl.F/Fmsy.vec; ucl.F.Fmsy<-ucl.F/Fmsy.vec
  
  F.last     <-Fm[nyr];lcl.F.last<-lcl.F[nyr];ucl.F.last<-ucl.F[nyr]
  Fmsy.last  <-Fmsy.vec[nyr];lcl.Fmsy.last<-lcl.Fmsy.vec[nyr];ucl.Fmsy.last<-ucl.Fmsy.vec[nyr]
  F.Fmsy.last<-F.Fmsy[nyr];lcl.F.Fmsy.last<-lcl.F.Fmsy[nyr];ucl.F.Fmsy.last<-ucl.F.Fmsy[nyr]
  
  if(is.na(select.yr)==F){
    B.sel<-B.Bmsy.sel*Bmsy
    F.sel<-ct.raw[yr==select.yr]/B.sel
    F.Fmsy.sel<-F.sel/Fmsy.vec[yr==select.yr]
  }
  
  
    #---------------------------------------------
    # NEW SEPERATION between CMSY and JAGS # ><>
    #---------------------------------------------
    MSY.cmsy   <-MSY.est; lcl.MSY.cmsy<-lcl.MSY.est; ucl.MSY.cmsy<-ucl.MSY.est 
    Bmsy.cmsy  <-k.est/2; lcl.Bmsy.cmsy<-lcl.k.est/2; ucl.Bmsy.cmsy<-ucl.k.est/2
    Fmsy.cmsy  <-r.est/2; lcl.Fmsy.cmsy<-lcl.r.est/2; ucl.Fmsy.cmsy<-ucl.r.est/2
    B.Bmsy.cmsy <-2*median.btv[1:nyr];lcl.B.Bmsy.cmsy<-2*lcl.btv[1:nyr];ucl.B.Bmsy.cmsy<-2*ucl.btv[1:nyr]
    if(is.na(select.yr)==F){B.Bmsy.sel<-2*median.btv[yr==select.yr]}
    B.cmsy          <-B.Bmsy.cmsy*Bmsy.cmsy;lcl.B.cmsy<-lcl.B.Bmsy.cmsy*Bmsy.cmsy;ucl.B.cmsy<-ucl.B.Bmsy.cmsy*Bmsy.cmsy
    Fm.cmsy           <- ct.raw/B.cmsy;lcl.F.cmsy<-ct.raw/ucl.B.cmsy;ucl.F.cmsy<-ct.raw/lcl.B.cmsy
    Fmsy.vec.cmsy     <- ifelse(B.Bmsy.cmsy>0.5,Fmsy.cmsy,Fmsy.cmsy*2*B.Bmsy.cmsy)
    lcl.Fmsy.vec.cmsy <- ifelse(B.Bmsy.cmsy>0.5,lcl.Fmsy.cmsy,lcl.Fmsy.cmsy*2*B.Bmsy.cmsy)
    ucl.Fmsy.vec.cmsy <- ifelse(B.Bmsy.cmsy>0.5,ucl.Fmsy.cmsy,ucl.Fmsy.cmsy*2*B.Bmsy.cmsy)
    F.Fmsy.cmsy       <- Fm.cmsy/Fmsy.vec.cmsy; lcl.F.Fmsy.cmsy<-lcl.F.cmsy/Fmsy.vec.cmsy; ucl.F.Fmsy.cmsy<-ucl.F.cmsy/Fmsy.vec.cmsy
    
    if(FullSchaefer==T){
    # NEW SEPERATION # ><>
    Bmsy.jags  <-k.jags/2; lcl.Bmsy.jags<-lcl.k.jags/2; ucl.Bmsy.jags<-ucl.k.jags/2
    Fmsy.jags  <-r.jags/2; lcl.Fmsy.jags<-lcl.r.jags/2; ucl.Fmsy.jags<-ucl.r.jags/2
    B.Bmsy.jags<-2*quant.P[2,];lcl.B.Bmsy.jags<-2*quant.P[1,];ucl.B.Bmsy.jags<-2*quant.P[3,]
    if(is.na(select.yr)==F) {B.Bmsy.sel<-2*quant.P[2,][yr==select.yr]}
    B.jags          <-B.Bmsy.jags*Bmsy.jags;lcl.B.jags<-lcl.B.Bmsy.jags*Bmsy.jags;ucl.B.jags<-ucl.B.Bmsy.jags*Bmsy.jags
    Fm.jags           <- ct.raw/B.jags;lcl.F.jags<-ct.raw/ucl.B.jags;ucl.F.jags<-ct.raw/lcl.B.jags
    Fmsy.vec.jags     <- ifelse(B.Bmsy.jags>0.5,Fmsy.jags,Fmsy.jags*2*B.Bmsy.jags)
    lcl.Fmsy.vec.jags <- ifelse(B.Bmsy.jags>0.5,lcl.Fmsy.jags,lcl.Fmsy.jags*2*B.Bmsy.jags)
    ucl.Fmsy.vec.jags <- ifelse(B.Bmsy.jags>0.5,ucl.Fmsy.jags,ucl.Fmsy.jags*2*B.Bmsy.jags)
    F.Fmsy.jags       <- Fm.jags/Fmsy.vec.jags; lcl.F.Fmsy.jags<-lcl.F.jags/Fmsy.vec.jags; ucl.F.Fmsy.jags<-ucl.F.jags/Fmsy.vec.jags
    } # End 
  
    
  # ------------------------------------------
  # print input and results to screen
  #-------------------------------------------
  cat("---------------------------------------\n")
  cat("Species:", cinfo$ScientificName[cinfo$Stock==stock], ", stock:",stock,"\n")
  cat(cinfo$Name[cinfo$Stock==stock], "\n")
  cat("Region:",cinfo$Region[cinfo$Stock==stock],",",cinfo$Subregion[cinfo$Stock==stock],"\n")
  cat("Catch data used from years", min(yr),"-", max(yr),", abundance =", btype, "\n")
  cat("Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n")
  cat("Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n")
  cat("Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,"default","expert"), "\n")
  cat("Prior range for r =", format(start.r[1],digits=2), "-", format(start.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),  
      ", prior range for k =", start.k[1], "-", start.k[2],"\n")
  # if Schaefer and CPUE, print prior range of q
  if(FullSchaefer==T & btype=="CPUE") {
    cat("Prior range of q =",q.prior[1],"-",q.prior[2],"\n")
  }
  # results of CMSY analysis
  cat("\nResults of CMSY analysis \n")
  cat("-------------------------\n")
  cat("Altogether", n.viable.b, "viable trajectories for", n.viable.pt," r-k pairs were found \n")
  cat("r   =", r.est,", 95% CL =", lcl.r.est, "-", ucl.r.est,", k =", k.est,", 95% CL =", lcl.k.est, "-", ucl.k.est,"\n")
  cat("MSY =", MSY.est,", 95% CL =", lcl.MSY.est, "-", ucl.MSY.est,"\n")
  cat("Relative biomass in last year =", median.btv.lastyr, "k, 2.5th perc =", lcl.median.btv.lastyr, 
      ", 97.5th perc =", ucl.median.btv.lastyr,"\n")
  cat("Exploitation F/(r/2) in last year =", (F.CMSY/Fmsy.CMSY)[nyr],"\n\n")
  
  # print results from full Schaefer if available
  if(FullSchaefer==T) {
    cat("Results from Bayesian Schaefer model (BSM) using catch &",btype,"\n")
    cat("------------------------------------------------------------\n")
    if(btype == "CPUE") cat("q   =", mean.q,", lcl =", lcl.q, ", ucl =", ucl.q,"\n")
    cat("r   =", r.jags,", 95% CL =", lcl.r.jags, "-", ucl.r.jags,", k =", k.jags,", 95% CL =", lcl.k.jags, "-", ucl.k.jags,"\n")
    cat("MSY =", MSY.jags,", 95% CL =", lcl.MSY.jags, "-", ucl.MSY.jags,"\n")
    cat("Relative biomass in last year =", quant.P[2,][nyr], "k, 2.5th perc =",quant.P[1,][nyr], 
        ", 97.5th perc =", quant.P[3,][nyr],"\n")
    cat("Exploitation F/(r/2) in last year =", (ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2) ,"\n\n")
  }
  
  # print results to be used in management
  cat("Results for Management (based on",ifelse(FullSchaefer==F | force.cmsy==T,"CMSY","BSM"),"analysis) \n")
  cat("-------------------------------------------------------------\n")
  if(force.cmsy==T) cat("Mangement results based on CMSY because abundance data seem unrealistic\n")
  cat("Fmsy =",Fmsy,", 95% CL =",lcl.Fmsy,"-",ucl.Fmsy,"(if B > 1/2 Bmsy then Fmsy = 0.5 r)\n") 
  cat("Fmsy =",Fmsy.last,", 95% CL =",lcl.Fmsy.last,"-",ucl.Fmsy.last,"(r and Fmsy are linearly reduced if B < 1/2 Bmsy)\n")
  cat("MSY  =",MSY,", 95% CL =",lcl.MSY,"-",ucl.MSY,"\n") 
  cat("Bmsy =",Bmsy,", 95% CL =",lcl.Bmsy,"-",ucl.Bmsy,"\n") 
  cat("Biomass in last year =",B.last,", 2.5th perc =", lcl.B.last, ", 97.5 perc =",ucl.B.last,"\n")
  cat("B/Bmsy in last year  =",B.Bmsy.last,", 2.5th perc =", lcl.B.Bmsy.last, ", 97.5 perc =",ucl.B.Bmsy.last,"\n")
  cat("Fishing mortality in last year =",F.last,", 2.5th perc =", lcl.F.last, ", 97.5 perc =",ucl.F.last,"\n")
  cat("Exploitation F/Fmsy  =",F.Fmsy.last,", 2.5th perc =", lcl.F.Fmsy.last, ", 97.5 perc =",ucl.F.Fmsy.last,"\n")
  
  # show stock status and exploitation for optional selected year
  if(is.na(select.yr)==F) {
    cat("\nStock status and exploitation in",select.yr,"\n")
    cat("Biomass =",B.sel, ", B/Bmsy =",B.Bmsy.sel,", F =",F.sel,", F/Fmsy =",F.Fmsy.sel,"\n") }
  
  # indicate if less than 5 years of biomass or CPUE are available
  if(btype !="None" & length(bt[is.na(bt)==F])<nab) {
    cat("Less than",nab,"years with abundance data available, shown on second axis\n")
  }
  cat("Comment:", comment,"\n")
  cat("----------------------------------------------------------\n")
  
  
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(10,5)}
  par(mfrow=c(1,2),mex=0.75) 
  # -----------------------------------------
  # Plot results 
  # -----------------------------------------
  # Analysis of viable r-k plot
  # ----------------------------
  max.y    <- max(c(ifelse(FullSchaefer==T,ucl.k.jags,NA), max(kv.all)), 
                  ifelse(substr(id_file,1,3)=="Sim",1.2*true.k,max(kv.all)),    
                  na.rm=T)
  min.y    <- min(c(ifelse(FullSchaefer==T,lcl.k.jags,NA), 0.9*min(kv.all)), 
                  ifelse(substr(id_file,1,3)=="Sim",0.8*true.k,0.9*min(kv.all)),    
                  na.rm=T)
  
  plot(x=rv.all, y=kv.all, xlim=start.r, 
       ylim=c(min.y,max.y), 
       pch=16, col="gray",log="xy", bty="l",
       xlab="r", ylab="k", main="A: Analysis of viable r-k",cex.main=.8)
  
  # plot r-k pairs from MCMC
  if(FullSchaefer==T) {points(x=r_out, y=k_out, pch=16,cex=0.5)}
  
  # plot best r-k from full Schaefer analysis
  if(FullSchaefer==T) {
    points(x=r.jags, y=k.jags, pch=19, col="red")  
    lines(x=c(lcl.r.jags, ucl.r.jags),y=c(k.jags,k.jags), col="red")
    lines(x=c(r.jags,r.jags),y=c(lcl.k.jags, ucl.k.jags), col="red")
  }
  
  # plot blue dot for CMSY r-k, with 95% CL lines 
  points(x=r.est, y=k.est, pch=19, col="blue")
  lines(x=c(lcl.r.est, ucl.r.est),y=c(k.est,k.est), col="blue")
  lines(x=c(r.est,r.est),y=c(lcl.k.est, ucl.k.est), col="blue")
  
  
  # Pred. abundance plot #><>
  #--------------------
  # determine k to use for red line in b/k plot 
  if(FullSchaefer==T)  {k2use <- k.jags} else {k2use <- k.est}
  # determine hight of y-axis in plot
  max.y  <- max(c(ifelse(btype=="biomass",max(bt/k2use,na.rm=T),NA),
                  ifelse(btype=="CPUE",max(bt/(mean.q*k2use),na.rm=T),NA),
                  max(ucl.btv),0.6,startbio[2], endbio[2]),
                ifelse(FullSchaefer==T & btype=="biomass",max(bt[is.na(bt)==F]/lcl.k.jags,na.rm=T),NA),
                ifelse(FullSchaefer==T & btype=="CPUE",1.1*max(bt/(mean.q*lcl.k.jags),na.rm=T),NA), na.rm=T)
  
  # Main plot of relative Adandunce
  plot(x=yr,y=norm.mbt[1:nyr], lwd=1.5, xlab="Year", ylab="Relative Abundance", type="n",
       ylim=ifelse(FullSchaefer==rep(T,2),range(c(norm.lbt,norm.ubt,quant.AI/mean(quant.AI[2,]))),range(c(norm.lbt,norm.ubt))), bty="l",  main="B: Relative Abundance ",cex.main=.8,col="blue")
  
  
  
  # if observed biomass is available, plot red biomass line (use non-smoothed bt)
  if(FullSchaefer==T) {
    polygon(c(yr,rev(yr)),c(quant.AI[1,]/mean(quant.AI[2,]),rev(quant.AI[3,]/mean(quant.AI[2,]))),border=0,col="grey")
    }
  polygon(c(yr,rev(yr)),c(norm.lbt[1:nyr],rev(norm.ubt[1:nyr])),border=0,col=grey(0.4,0.6))
  lines(x=yr, y=norm.mbt[1:nyr],type="l",col=4,lwd=2,lty=2)
  if(FullSchaefer==T){
    norm.bt <- ifelse(is.na(bt),NA,bt/mean(quant.AI[2,]))
    points(yr,norm.bt,pch=21,bg="white")
    lines(yr,quant.AI[2,]/mean(quant.AI[2,]),lwd=2,col=2)
  }
  
  
  
  #save analytic chart to JPEG file
  if (save.plots==TRUE) 
  {
    jpgfile<-paste(stock,"_AN.jpeg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 500, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  
  
  #---------------------------------------------
  # Plot Management-Graphs if desired
  #---------------------------------------------
  if(mgraphs==T) {
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,10)}
    par(mfrow=c(2,2),mex=0.85) 
    # make margins narrower
    par(mar=c(3.1,4.1,2.1,2.1))
    
    
    #----------------------------------------
    # plot estimated biomass relative to Bmsy
    #----------------------------------------
    # plot empty frame
    plot(yr, rep(0,nyr),type="n", ylim=c(0,max(c(2, max(ifelse(FullSchaefer==T,c(ucl.B.Bmsy.jags,ucl.B.Bmsy.cmsy),ucl.B.Bmsy.cmsy))))), ylab="",xlab="", main="A: Biomass", bty="l",cex.main=1.2) 
    mtext("B / Bmsy",side=2, line=2,cex=1.1)
    mtext("Year",side=1, line=2,cex=1.1)
    
    # plot gray area of uncertainty in predicted biomass
    polygon(c(yr,rev(yr)), c(lcl.B.Bmsy.cmsy,rev(ucl.B.Bmsy.cmsy)),col=rgb(0,0,1,0.5),border=NA)
    
    if(FullSchaefer==T){
    polygon(c(yr,rev(yr)), c(lcl.B.Bmsy.jags,rev(ucl.B.Bmsy.jags)),col=rgb(1,0,0,0.5),border=NA)
    lines(yr,B.Bmsy.jags,lwd=2,lty=1,col="red")
    }  
    # plot median biomass
    lines(yr,B.Bmsy.cmsy,lwd=2,lty=2,col="darkblue")
    
    # plot lines for Bmsy and 0.5 Bmsy
    lines(x=c(yr[1],yr[nyr]),y=c(1,1), lty="dashed", lwd=1.5)
    lines(x=c(yr[1],yr[nyr]),y=c(0.5,0.5), lty="dotted", lwd=1.5)
    
    
    #---------------------------
    # plot exploitation rate
    # -------------------------
    # plot empty frame
    plot(yr, rep(0,nyr),type="n", ylim=c(0,max(c(2, min(5,max(ifelse(FullSchaefer==rep(T,nyr*2),c(ucl.F.Fmsy.jags,ucl.F.Fmsy.cmsy),ucl.F.Fmsy.cmsy)))))), 
         ylab="",xlab="", main="B: Exploitation", bty="l",cex.main=1.2)  
    # plot gray area of uncertainty in predicted exploitation
    polygon(c(yr,rev(yr)), c(lcl.F.Fmsy.cmsy,rev(ucl.F.Fmsy.cmsy)),col=rgb(0,0,1,0.5),border=NA)
    
    if(FullSchaefer==T){
    polygon(c(yr,rev(yr)), c(lcl.F.Fmsy.jags,rev(ucl.F.Fmsy.jags)),col=rgb(1,0,0,0.5),border=NA)
    lines(yr,F.Fmsy.jags,lwd=2,lty=1,col="red")
    }
    mtext("H / Hmsy",side=2, line=2,cex=1.1)
    mtext("Year",side=1, line=2,cex=1.1)
    
    # plot median biomass
    lines(yr,F.Fmsy.cmsy,lwd=2,lty=2,col="darkblue")
    
    # plot line for u.msy
    lines(x=c(yr[1],yr[nyr]),y=c(1,1), lty="dashed", lwd=1.5)
    
    
    #---------------------
    # plot catch with MSY
    #---------------------
    max.y <- max(c(1.1*max(ct.raw),ifelse(FullSchaefer==T,ucl.MSY.jags,0),ucl.MSY.cmsy),na.rm=T)
    plot(x=yr,rep(0,nyr),type="n",ylim=c(0,max.y), bty="l", main=paste("C: Catch and MSY"),cex.main=1.2, 
         ylab="",xlab="")
    mtext("Catch in 1000 t",side=2, line=2,cex=1.1)
    mtext("Year",side=1, line=2,cex=1.1)
    
    rect(yr[1],lcl.MSY.cmsy,yr[nyr],ucl.MSY.cmsy,col=rgb(0,0,1,0.5), border=NA)
    
    if(FullSchaefer==T){
    rect(yr[1],lcl.MSY.jags,yr[nyr],ucl.MSY.jags,col=rgb(1,0,0,0.5), border=NA)
    lines(x=c(yr[1],yr[nyr]),y=c(MSY.jags,MSY.jags),lty="dashed", col="red", lwd=2)
    }
    lines(x=c(yr[1],yr[nyr]),y=c(MSY.cmsy,MSY.cmsy),lty="dashed", col="darkblue", lwd=2)
    lines(x=yr, y=ct.raw, lwd=2)
    
    
    #><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
    # plot stock-status graph
    #---------------------------
    if(FullSchaefer==T & force.cmsy==F) {x.F_Fmsy.jags = all.F_Fmsy[,nyr]
    y.b_bmsy.jags = all.b_bmsy[,nyr]} 
      
    log.sd.B.Bmsy = (log(ucl.B.Bmsy.cmsy[nyr]+0.0011)-log(lcl.B.Bmsy.cmsy[nyr]+0.001))/(2*1.96)   
    log.sd.F.Fmsy = (log(ucl.F.Fmsy.cmsy[nyr]+0.005)-log(lcl.F.Fmsy.cmsy[nyr]+0.001))/(2*1.96)   
    x.F_Fmsy.cmsy= rlnorm(20000,log(F.Fmsy.cmsy[nyr]+0.001),log.sd.F.Fmsy)
    y.b_bmsy.cmsy =rlnorm(20000,log(B.Bmsy.cmsy[nyr]+0.001),log.sd.B.Bmsy)
  
    # ICCAT Kobe plot
    # Kernel CMSY 
    kernelF.cmsy <- ci2d(y.b_bmsy.cmsy,x.F_Fmsy.cmsy,nbins=101,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none")
    if(FullSchaefer==T){
    kernelF.jags <- ci2d(y.b_bmsy.jags,x.F_Fmsy.jags,nbins=101,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none")}
    c1 <- c(-1,100)
    c2 <- c(1,1)
    
    if(FullSchaefer==T){
    max.y1   <- max(c(2, max(kernelF.cmsy$contours$"0.95"$y*1.05,kernelF.jags$contours$"0.95"$y*1.05),na.rm =T))
    max.y    <- ifelse(max.y1 > 5,min(max(5,F.Fmsy.cmsy*2,F.Fmsy.jags*2),12),max.y1)
    } else {
      max.y1   <- max(c(2, max(kernelF.cmsy$contours$"0.95"$y*1.05),na.rm =T))
      max.y    <- ifelse(max.y1 > 5,min(max(5,F.Fmsy.cmsy*2),8),max.y1)
      }
    if(FullSchaefer==T){max.x    <- max(2,quantile(c(y.b_bmsy.cmsy,y.b_bmsy.jags),0.96))} else {
    max.x    <- max(2,quantile(c(y.b_bmsy.cmsy),0.96))}
    
    
    #Create plot
    plot(1000,1000,type="b", xlim=c(0,max.x), ylim=c(0,max.y),lty=3,xlab="",ylab="", bty="l",main="D: Kobe Plot",cex.main=1.2)
    mtext("H / Hmsy",side=2, line=2,cex=1)
    mtext("B / Bmsy",side=1, line=2,cex=1)
    
    # and fill areas using the polygon function
    zb2 = c(0,1)
    zf2  = c(1,100)
    zb1 = c(1,100)
    zf1  = c(0,1)
    polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
    polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
    polygon(c(1,100,100,1),c(1,1,100,100),col="yellow",border=0)
    polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
    
    
    # extract interval information from ci2d object
    # and fill areas using the polygon function
    polygon(kernelF.cmsy$contours$"0.95",lty=2,border=NA,col="cornsilk4")
    polygon(kernelF.cmsy$contours$"0.8",border=NA,lty=2,col="grey")
    polygon(kernelF.cmsy$contours$"0.5",border=NA,lty=2,col="cornsilk2")
    
    if(FullSchaefer==TRUE){ #><> 
    # and fill areas using the polygon function
    polygon(kernelF.jags$contours$"0.95",lty=2,border=NA,col=grey(0.3,0.6))
    polygon(kernelF.jags$contours$"0.8",border=NA,lty=2,col="grey")
    polygon(kernelF.jags$contours$"0.5",border=NA,lty=2,col="cornsilk2")
    }
    
      ## Add points and trajectory lines
      lines(c1,c2,lty=3,lwd=0.7)
      lines(c2,c1,lty=3,lwd=0.7)
      lines(B.Bmsy.cmsy,F.Fmsy.cmsy, lty=1,lwd=3,col="white")
      lines(B.Bmsy.cmsy,F.Fmsy.cmsy, lty=1,lwd=1,col=4)
      points(B.Bmsy.cmsy,F.Fmsy.cmsy,cex=0.8,pch=4,col=4)
      
      if(FullSchaefer==TRUE){ #><> plot both with cross-hair here
        lines(B.Bmsy.jags,F.Fmsy.jags, lty=1,lwd=3.,col="white")
        lines(B.Bmsy.jags,F.Fmsy.jags, lty=1,lwd=1.,col=2)
        points(B.Bmsy.jags,F.Fmsy.jags,cex=0.8,pch=4,col=2)
      }  
      
      points(B.Bmsy.cmsy[1],F.Fmsy.cmsy[1],col="white",pch=22,bg=4,cex=1.9)
      points(B.Bmsy.cmsy[which(yr==int.yr)],F.Fmsy.cmsy[which(yr==int.yr)],col="white",pch=21,bg="blue",cex=1.9)
      points(B.Bmsy.cmsy[nyr],F.Fmsy.cmsy[nyr],col="white",pch=24,bg="blue",cex=1.9)
      
      if(FullSchaefer==TRUE){ #><> plot both with cross-hair here
      points(B.Bmsy.jags[1],F.Fmsy.jags[1],col="white",pch=22,bg="red",cex=1.9)
      points(B.Bmsy.jags[which(yr==int.yr)],F.Fmsy.jags[which(yr==int.yr)],col="white",pch=21,bg="red",cex=1.9)
      points(B.Bmsy.jags[nyr],F.Fmsy.jags[nyr],col="white",pch=24,bg="red",cex=1.9)
      }
      
      # Get Propability
      Pr.green = sum(ifelse(y.b_bmsy.cmsy>1 & x.F_Fmsy.cmsy<1,1,0))/length(y.b_bmsy.cmsy)*100
      Pr.red = sum(ifelse(y.b_bmsy.cmsy<1 & x.F_Fmsy.cmsy>1,1,0))/length(y.b_bmsy.cmsy)*100
      Pr.yellow = (sum(ifelse(y.b_bmsy.cmsy<1 & x.F_Fmsy.cmsy<1,1,0))+sum(ifelse(y.b_bmsy.cmsy>1 & x.F_Fmsy.cmsy>1,1,0)))/length(y.b_bmsy.cmsy)*100
      
      if(FullSchaefer==TRUE){
      ## Add legend
      legend('topright', c(paste(start.yr),paste(int.yr),paste(end.yr),"50% C.I.","80% C.I.","95% C.I."), 
               lty=c(1,1,1,-1,-1,-1),pch=c(22,21,24,22,22,22),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4"), 
               col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),1.7,1.7,1.7),bty="n",y.intersp = 0.9)
        
      } else {  
      ## Add legend
      legend('topright', 
             c(paste(start.yr),paste(int.yr),paste(end.yr),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.green),1),"%")), 
             lty=c(1,1,1,rep(-1,7)),pch=c(22,21,24,rep(22,7)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","green"), 
             col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.1,3)),bty="n")
      }
      
      
      
      
    #><> End of Biplot
    
    
  
  # save management chart to JPEG file
  if (save.plots==TRUE & mgraphs==TRUE) 
  {
    jpgfile<-paste(stock,"_MAN.jpeg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 850, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  } # end of management graphs
  
  #---------------------------------
  # Plot prior postrior
  #---------------------------------
  if(mgraphs==T) {
  # Priors ><>
  Prior.k <- log(start.k[1]+0.2*(start.k[2]-start.k[1])) 
  Prior.k[2] <- (Prior.k[1]-log(start.k[1]))/4
  Prior.r  <- mean(log(start.r))
  Prior.r[2] <- abs(Prior.r[1]- log(start.r[1]))/2
  
  # open window for plot of four panels
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,10)}
  par(mfrow=c(2,2),mex=0.95) 
  # make margins narrower
  par(mar=c(4.4,4.1,.1,.1),cex.lab=1.1)
  
  # Show Prior, Posterior ><>
  # k
  if(FullSchaefer==T){
  out.jags = as.numeric(all.k)
  parm.jags = fitdist(out.jags, "gamma")$estimate
  pdf.jags = dgamma(sort(out.jags),parm.jags[1],parm.jags[2])  
  }
  #CMSY
  out.cmsy = rlnorm(10000,log.k.est,sd.log.k.est) 
  parm.cmsy = fitdist(out.cmsy, "gamma")$estimate
  pdf.cmsy = dgamma(sort(out.cmsy),parm.cmsy[1],parm.cmsy[2])  
  # Prior
  rpr = rlnorm(10000,Prior.k[1],Prior.k[2]) 
  prior = dlnorm((sort(rpr)),Prior.k[1],Prior.k[2])   
  
  if(FullSchaefer==T){
    xlim = range(quantile(c(out.jags,out.cmsy,rpr),c(0.,0.99)))
    ylim = c(0,1.1*max(pdf.cmsy,pdf.jags,prior))} else {
    xlim = range(quantile(c(out.cmsy,rpr),c(0.,0.99)))
    ylim = c(0,1.1*max(pdf.cmsy,prior))} 
  plot(sort(out.cmsy),pdf.cmsy,type="n",xlim=xlim,ylim=ylim,yaxt="n",xlab="k",ylab="",xaxs="i",yaxs="i")
  polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.5,1))
  polygon(c(sort(out.cmsy),rev(sort(out.cmsy))),c(pdf.cmsy,rep(0,length(sort(out.cmsy)))),col=rgb(0,0,1,0.4),border="white")
  if(FullSchaefer==T){
  polygon(c(sort(out.jags),rev(sort(out.jags))),c(pdf.jags,rep(0,length(sort(out.jags)))),col=rgb(1,0,0,0.4),border="white")
  legend('topright',c("Prior","CMSY","BSM"),lwd=c(rep(-1,3)),pch=c(rep(22,3)),lty=c(rep(-1,3)),pt.cex=1.5,pt.bg = c(grey(0.5,1),rgb(0,0,1,0.5),rgb(1,0,0,0.5)),bty="n")
  abline(v=mean(out.jags),col=2,lwd=2,lty=2)} else {
    legend('topright',c("Prior","CMSY"),lwd=c(rep(-1,2)),pch=c(rep(22,2)),lty=c(rep(-1,2)),pt.cex=1.5,pt.bg = c(grey(0.5,1),rgb(1,0,0,0.5)),bty="n")
  }
  abline(v=mean(out.cmsy),col=4,lwd=2,lty=2)
  #mtext("k",side=1, line=2,cex=1.1)
  # r
  if(FullSchaefer==T){
  out.jags = as.numeric(all.r)
  parm.jags = fitdist(out.jags, "gamma")$estimate
  pdf.jags = dgamma(sort(out.jags),parm.jags[1],parm.jags[2])  
  }
  #CMSY
  out.cmsy = rlnorm(10000,log.r.est,sd.log.r.est) 
  parm.cmsy = fitdist(out.cmsy, "gamma")$estimate
  pdf.cmsy = dgamma(sort(out.cmsy),parm.cmsy[1],parm.cmsy[2])  
  # Prior
  rpr = rlnorm(10000,Prior.r[1],Prior.r[2]) 
  prior = dlnorm((sort(rpr)),Prior.r[1],Prior.r[2])   
  
  if(FullSchaefer==T){
  xlim = range(quantile(c(out.jags,out.cmsy,rpr),c(0.,0.99)))
  ylim = c(0,1.1*max(pdf.cmsy,pdf.jags,prior))} else {
  xlim = range(quantile(c(out.cmsy,rpr),c(0.,0.99)))
  ylim = c(0,1.1*max(pdf.cmsy,prior))} 
  plot(sort(out.cmsy),pdf.cmsy,type="n",xlim=xlim,ylim=ylim,yaxt="n",xlab="r",ylab="",xaxs="i",yaxs="i")
  polygon(c(sort(rpr),rev(sort(rpr))),c(prior,rep(0,length(sort(rpr)))),col=gray(0.5,1))
  polygon(c(sort(out.cmsy),rev(sort(out.cmsy))),c(pdf.cmsy,rep(0,length(sort(out.cmsy)))),col=rgb(0,0,1,0.4),border="white")
  
  if(FullSchaefer==T){
  polygon(c(sort(out.jags),rev(sort(out.jags))),c(pdf.jags,rep(0,length(sort(out.jags)))),col=rgb(1,0,0,0.4),border="white")
  abline(v=mean(out.jags),col=2,lwd=2,lty=2)
  }
  abline(v=mean(out.cmsy),col=4,lwd=2,lty=2)
  
  # B/Bmsy jags
  if(FullSchaefer==T){
  out.jags = y.b_bmsy.jags
  parm.jags = fitdist(out.jags, "gamma")$estimate
  pdf.jags = dgamma(sort(out.jags),parm.jags[1],parm.jags[2])  
  }
  #CMSY
  out.cmsy =  y.b_bmsy.cmsy
  parm.cmsy = fitdist(out.cmsy, "gamma")$estimate
  pdf.cmsy = dgamma(sort(out.cmsy),parm.cmsy[1],parm.cmsy[2])  
  
  if(FullSchaefer==T){
    xlim = range(quantile(c(out.jags,out.cmsy),c(0.,0.99)))
    ylim = c(0,1.1*max(pdf.cmsy,pdf.jags))} else {
      xlim = range(quantile(c(out.cmsy),c(0.,0.99)))
      ylim = c(0,1.1*max(pdf.cmsy))} 
  plot(sort(out.cmsy),pdf.cmsy,type="n",xlim=xlim,ylim=ylim,yaxt="n",xlab=paste0("B",end.yr," / Bmsy"),ylab="",xaxs="i",yaxs="i")
  polygon(c(sort(out.cmsy),rev(sort(out.cmsy))),c(pdf.cmsy,rep(0,length(sort(out.cmsy)))),col=rgb(0,0,1,0.4),border="white")
  if(FullSchaefer==T){
  polygon(c(sort(out.jags),rev(sort(out.jags))),c(pdf.jags,rep(0,length(sort(out.jags)))),col=rgb(1,0,0,0.4),border="white")
  abline(v=mean(out.jags),col=2,lwd=2,lty=2)
  }
  abline(v=mean(out.cmsy),col=4,lwd=2,lty=2)
  abline(v=1,col=1,lwd=2,lty=5)
  
  # F/Fmsy 
  if(FullSchaefer==T){
  out.jags = x.F_Fmsy.jags
  parm.jags = fitdist(out.jags, "gamma")$estimate
  pdf.jags = dgamma(sort(out.jags),parm.jags[1],parm.jags[2])  
  }
  #CMSY
  out.cmsy =  x.F_Fmsy.cmsy
  parm.cmsy = fitdist(out.cmsy, "gamma")$estimate
  pdf.cmsy = dgamma(sort(out.cmsy),parm.cmsy[1],parm.cmsy[2])  
  
  if(FullSchaefer==T){
    xlim = range(quantile(c(out.jags,out.cmsy),c(0.,0.99)))
    ylim = c(0,1.1*max(pdf.cmsy,pdf.jags))} else {
      xlim = range(quantile(c(out.cmsy),c(0.,0.99)))
      ylim = c(0,1.1*max(pdf.cmsy))} 
  plot(sort(out.cmsy),pdf.cmsy,type="n",xlim=xlim,ylim=ylim,yaxt="n",xlab=paste0("H",end.yr," / Hmsy"),ylab="",xaxs="i",yaxs="i")
  polygon(c(sort(out.cmsy),rev(sort(out.cmsy))),c(pdf.cmsy,rep(0,length(sort(out.cmsy)))),col=rgb(0,0,1,0.4),border="white")
  if(FullSchaefer==T){
  polygon(c(sort(out.jags),rev(sort(out.jags))),c(pdf.jags,rep(0,length(sort(out.jags)))),col=rgb(1,0,0,0.4),border="white")
  abline(v=mean(out.jags),col=2,lwd=2,lty=2)}
  abline(v=mean(out.cmsy),col=4,lwd=2,lty=2)
  abline(v=1,col=1,lwd=2,lty=5)
  
  # save management chart to JPEG file
  if (save.plots==TRUE & mgraphs==TRUE) 
  {
    jpgfile<-paste(stock,"_Prior.jpeg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 850, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  } # end of PRIOR/Posterior graphs
  
  #----------------------------------
  # Make single ICCAT CMSY KOBE Plot
  #---------------------------------
  
  # ICCAT Kobe plot
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,10)}
  par(mfrow=c(1,1),mex=0.95) 
  # make margins narrower
  par(mar=c(4.4,4.1,.1,.1),cex.lab=1.1)
  
  # Kernel CMSY 
  c1 <- c(-1,100)
  c2 <- c(1,1)
  
  max.y1   <- max(c(2, max(kernelF.cmsy$contours$"0.95"$y*1.05),na.rm =T))
  max.y    <- ifelse(max.y1 > 5,min(max(5,F.Fmsy.cmsy*2),8),max.y1)
  
  #Create plot
  plot(1000,1000,type="b", xlim=c(0,max.x), ylim=c(0,max.y),lty=3,xlab="",ylab="", bty="l",cex.main=1.2)
  mtext("H / Hmsy",side=2, line=2,cex=1.2)
  mtext("B / Bmsy",side=1, line=2,cex=1.2)
  
  # and fill areas using the polygon function
  zb2 = c(0,1)
  zf2  = c(1,100)
  zb1 = c(1,100)
  zf1  = c(0,1)
  polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
  polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
  polygon(c(1,100,100,1),c(1,1,100,100),col="yellow",border=0)
  polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
  
  
  # extract interval information from ci2d object
  # and fill areas using the polygon function
  polygon(kernelF.cmsy$contours$"0.95",lty=2,border=NA,col="cornsilk4")
  polygon(kernelF.cmsy$contours$"0.8",border=NA,lty=2,col="grey")
  polygon(kernelF.cmsy$contours$"0.5",border=NA,lty=2,col="cornsilk2")
  
  ## Add points and trajectory lines
  lines(c1,c2,lty=3,lwd=0.7)
  lines(c2,c1,lty=3,lwd=0.7)
  lines(B.Bmsy.cmsy,F.Fmsy.cmsy, lty=1,lwd=3,col="white")
  lines(B.Bmsy.cmsy,F.Fmsy.cmsy, lty=1,lwd=1,col=4)
  points(B.Bmsy.cmsy,F.Fmsy.cmsy,cex=0.8,pch=4,col=4)
  
  points(B.Bmsy.cmsy[1],F.Fmsy.cmsy[1],col="white",pch=22,bg=4,cex=1.9)
  points(B.Bmsy.cmsy[which(yr==int.yr)],F.Fmsy.cmsy[which(yr==int.yr)],col="white",pch=21,bg="blue",cex=1.9)
  points(B.Bmsy.cmsy[nyr],F.Fmsy.cmsy[nyr],col="white",pch=24,bg="blue",cex=1.9)
  
  # Get Propability
  Pr.green = sum(ifelse(y.b_bmsy.cmsy>1 & x.F_Fmsy.cmsy<1,1,0))/length(y.b_bmsy.cmsy)*100
  Pr.red = sum(ifelse(y.b_bmsy.cmsy<1 & x.F_Fmsy.cmsy>1,1,0))/length(y.b_bmsy.cmsy)*100
  Pr.yellow = (sum(ifelse(y.b_bmsy.cmsy<1 & x.F_Fmsy.cmsy<1,1,0))+sum(ifelse(y.b_bmsy.cmsy>1 & x.F_Fmsy.cmsy>1,1,0)))/length(y.b_bmsy.cmsy)*100
  
  
  ## Add legend
  legend('topright', 
         c(paste(start.yr),paste(int.yr),paste(end.yr),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.green),1),"%")), 
         lty=c(1,1,1,rep(-1,7)),pch=c(22,21,24,rep(22,7)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","green"), 
         col=1,lwd=1.1,cex=0.9,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.1,3)),bty="n")
  
  if (save.plots==TRUE & mgraphs==TRUE) 
  {
    jpgfile<-paste(stock,"_KOBE.jpeg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 850, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    }
  dev.off()
  
  # -------------------------------------
  ## Write some results into csv outfile
  # -------------------------------------
  if(write.output == TRUE) {
    
    # write data into csv file
    output = data.frame(as.character(cinfo$Region[cinfo$Stock==stock]),
                        as.character(cinfo$Subregion[cinfo$Stock==stock]),
                        as.character(cinfo$Name[cinfo$Stock==stock]),
                        cinfo$ScientificName[cinfo$Stock==stock], 
                        stock, start.yr, end.yr, btype,
                        max(ct.raw),ct.raw[nyr],
                        ifelse(FullSchaefer==T,MSY.jags,NA), # full Schaefer
                        ifelse(FullSchaefer==T,lcl.MSY.jags,NA),
                        ifelse(FullSchaefer==T,ucl.MSY.jags,NA),
                        ifelse(FullSchaefer==T,r.jags,NA), 
                        ifelse(FullSchaefer==T,lcl.r.jags,NA),
                        ifelse(FullSchaefer==T,ucl.r.jags,NA),
                        ifelse(FullSchaefer==T,k.jags,NA),                                          
                        ifelse(FullSchaefer==T,lcl.k.jags,NA),                    
                        ifelse(FullSchaefer==T,ucl.k.jags,NA),
                        ifelse(FullSchaefer==T & btype=="CPUE",mean.q,NA),
                        ifelse(FullSchaefer==T & btype=="CPUE",lcl.q,NA),
                        ifelse(FullSchaefer==T & btype=="CPUE",ucl.q,NA),
                        ifelse(FullSchaefer==T,quant.P[2,][nyr],NA), # last B/k JAGS
                        ifelse(FullSchaefer==T,quant.P[1,][nyr],NA),
                        ifelse(FullSchaefer==T,quant.P[3,][nyr],NA), 
                        ifelse(FullSchaefer==T,(ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2),NA), # last F/Fmsy JAGS
                        r.est, lcl.r.est, ucl.r.est, # CMSY r
                        k.est, lcl.k.est, ucl.k.est, # CMSY k     
                        MSY.est, lcl.MSY.est, ucl.MSY.est, # CMSY MSY
                        median.btv.lastyr, lcl.median.btv.lastyr,ucl.median.btv.lastyr, # CMSY B/k in last year with catch data
                        (F.CMSY/Fmsy.CMSY)[nyr],
                        Fmsy,lcl.Fmsy,ucl.Fmsy,Fmsy.last,lcl.Fmsy.last,ucl.Fmsy.last,
                        MSY,lcl.MSY,ucl.MSY,Bmsy,lcl.Bmsy,ucl.Bmsy,
                        B.last, lcl.B.last, ucl.B.last, B.Bmsy.last, lcl.B.Bmsy.last, ucl.B.Bmsy.last,
                        F.last, lcl.F.last, ucl.F.last, F.Fmsy.last, lcl.F.Fmsy.last, ucl.F.Fmsy.last,
                        ifelse(is.na(select.yr)==F,B.sel,NA),
                        ifelse(is.na(select.yr)==F,B.Bmsy.sel,NA),
                        ifelse(is.na(select.yr)==F,F.sel,NA),
                        ifelse(is.na(select.yr)==F,F.Fmsy.sel,NA),
                        ifelse(yr[1]>2000,NA,ct.raw[yr==2000]),ifelse(yr[1]>2001,NA,ct.raw[yr==2001]),ifelse(yr[1]>2002,NA,ct.raw[yr==2002]), ifelse(yr[1]>2003,NA,ct.raw[yr==2003]),# allow missing 2000-2002
                        ct.raw[yr==2004],ct.raw[yr==2005],ct.raw[yr==2006],ct.raw[yr==2007],ct.raw[yr==2008],ct.raw[yr==2009],ct.raw[yr==2010],
                        ifelse(yr[nyr]<2011,NA,ct.raw[yr==2011]),ifelse(yr[nyr]<2012,NA,ct.raw[yr==2012]),ifelse(yr[nyr]<2013,NA,ct.raw[yr==2013]),ifelse(yr[nyr]<2014,NA,ct.raw[yr==2014]),ifelse(yr[nyr]<2015,NA,ct.raw[yr==2015]), # allow missing 2011-2015
                        ifelse(yr[1]>2000,NA,F.Fmsy[yr==2000]),ifelse(yr[1]>2001,NA,F.Fmsy[yr==2001]),ifelse(yr[1]>2002,NA,F.Fmsy[yr==2002]),ifelse(yr[1]>2003,NA,F.Fmsy[yr==2003]), # allow missing 2000-2002
                        F.Fmsy[yr==2004],F.Fmsy[yr==2005],F.Fmsy[yr==2006],F.Fmsy[yr==2007],F.Fmsy[yr==2008],F.Fmsy[yr==2009],F.Fmsy[yr==2010],
                        ifelse(yr[nyr]<2011,NA,F.Fmsy[yr==2011]),ifelse(yr[nyr]<2012,NA,F.Fmsy[yr==2012]),ifelse(yr[nyr]<2013,NA,F.Fmsy[yr==2013]),ifelse(yr[nyr]<2014,NA,F.Fmsy[yr==2014]),ifelse(yr[nyr]<2015,NA,F.Fmsy[yr==2015]),# allow missing 2011-2015
                        ifelse(yr[1]>2000,NA,B[yr==2000]),ifelse(yr[1]>2001,NA,B[yr==2001]),ifelse(yr[1]>2002,NA,B[yr==2002]),ifelse(yr[1]>2003,NA,B[yr==2003]), # allow missing 2000-2002
                        B[yr==2004],B[yr==2005],B[yr==2006],B[yr==2007],B[yr==2008],B[yr==2009],B[yr==2010],
                        ifelse(yr[nyr]<2011,NA,B[yr==2011]),ifelse(yr[nyr]<2012,NA,B[yr==2012]),ifelse(yr[nyr]<2013,NA,B[yr==2013]),ifelse(yr[nyr]<2014,NA,B[yr==2014]),ifelse(yr[nyr]<2015,NA,B[yr==2015])) # allow missing 2011-2015
                        
    write.table(output, file=outfile, append = T, sep = ",", 
                dec = ".", row.names = FALSE, col.names = FALSE)
    
    
    # write screen text into text outfile.txt
    
    cat("Species:", cinfo$ScientificName[cinfo$Stock==stock], ", stock:",stock,"\n",
        cinfo$Name[cinfo$Stock==stock], "\n",
        "Source:",cinfo$Source[cinfo$Stock==stock],"\n",
        "Region:",cinfo$Region[cinfo$Stock==stock],",",cinfo$Subregion[cinfo$Stock==stock],"\n",
        "Catch data used from years", min(yr),"-", max(yr),", abundance =", btype, "\n",
        "Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n",
        "Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n",
        "Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,", default","expert"), "\n",
        "Prior range for r =", format(start.r[1],digits=2), "-", format(start.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),  
        ", prior range for k  =", start.k[1], "-", start.k[2],
        file=outfile.txt,append=T)
    
    if(FullSchaefer==T & btype=="CPUE") {
      cat("\n Prior range of q =",q.prior[1],"-",q.prior[2],file=outfile.txt,append=T)
    }
    cat("\n\n Results of CMSY analysis with altogether",n.viable.b, "viable trajectories for", n.viable.pt,"r-k pairs \n",
        "r =", r.est,", 95% CL =", lcl.r.est, "-", ucl.r.est,
        ", k =", k.est,", 95% CL =", lcl.k.est, "-", ucl.k.est,"\n",
        "MSY =", MSY.est,", 95% CL =", lcl.MSY.est, "-", ucl.MSY.est,"\n",
        "Relative biomass last year =", median.btv.lastyr, "k, 2.5th =", lcl.median.btv.lastyr, 
        ", 97.5th =", ucl.median.btv.lastyr,"\n",
        "Exploitation F/(r/2) in last year =", (F.CMSY/Fmsy.CMSY)[length(median.btv)-1],"\n",
        file=outfile.txt,append=T)  
  
      if(FullSchaefer==T) {
      cat("\n Results from Bayesian Schaefer model using catch &",btype,"\n",
          "r =", r.jags,", 95% CL =", lcl.r.jags, "-", ucl.r.jags,
          ", k =", k.jags,", 95% CL =", lcl.k.jags, "-", ucl.k.jags,"\n",
          "MSY =", MSY.jags,", 95% CL =", lcl.MSY.jags, "-", ucl.MSY.jags,"\n",
          "Relative biomass in last year =", quant.P[2,][nyr], "k, 2.5th perc =",quant.P[1,][nyr], 
              ", 97.5th perc =", quant.P[3,][nyr],"\n",
          "Exploitation F/(r/2) in last year =", (ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2) ,
          file=outfile.txt,append=T)
       if(btype == "CPUE") {cat("\n q =", mean.q,", lcl =", lcl.q, ", ucl =", ucl.q,
          file=outfile.txt,append=T)}
    }
    
    cat("\n\n Results for Management (based on",ifelse(FullSchaefer==F | force.cmsy==T,"CMSY","BSM"),"analysis) \n",
        "Fmsy =",Fmsy,", 95% CL =",lcl.Fmsy,"-",ucl.Fmsy,"(if B > 1/2 Bmsy then Fmsy = 0.5 r)\n",
        "Fmsy =",Fmsy.last,", 95% CL =",lcl.Fmsy.last,"-",ucl.Fmsy.last,"(r and Fmsy are linearly reduced if B < 1/2 Bmsy)\n",
        "MSY  =",MSY,",  95% CL =",lcl.MSY,"-",ucl.MSY,"\n", 
        "Bmsy =",Bmsy,",  95% CL =",lcl.Bmsy,"-",ucl.Bmsy,"\n", 
        "Biomass in last year  =",B.last,", 2.5th perc =", lcl.B.last, ", 97.5 perc =",ucl.B.last,"\n",
        "B/Bmsy in last year   =",B.Bmsy.last,", 2.5th perc =", lcl.B.Bmsy.last, ", 97.5 perc =",ucl.B.Bmsy.last,"\n",
        "Fishing mortality in last year =",F.last,", 2.5th perc =", lcl.F.last, ", 97.5 perc =",ucl.F.last,"\n",
        "F/Fmsy  =",F.Fmsy.last,", 2.5th perc =", lcl.F.Fmsy.last, ", 97.5 perc =",ucl.F.Fmsy.last,"\n",
        file=outfile.txt,append=T)
    
    # show stock status and exploitation for optional selected year
    if(is.na(select.yr)==F) {
      cat("\n Stock status and exploitation in",select.yr,"\n",
      "Biomass =",B.sel, ", B/Bmsy =",B.Bmsy.sel,", fishing mortality F =",F.sel,", F/Fmsy =",F.Fmsy.sel,"\n",
      file=outfile.txt,append=T) }
    
    if(btype !="None" & length(bt[is.na(bt)==F])<nab) {
      cat(" Less than",nab,"years with abundance data available, shown on second axis\n",file=outfile.txt,append=T) }
    
    cat(" Comment:", comment,"\n","----------------------------------------------------------\n\n",
        file=outfile.txt,append=T)
   } # end of loop to write text to file
  
   if(close.plots==T) graphics.off() # close on-screen graphics windows after files are saved
  

  
  } # end of stocks loop
  
  #stop parallel processing clusters
  stopCluster(cl)
