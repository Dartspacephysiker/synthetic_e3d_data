source("ISgeometryMODIFIED.R")

cartesianToSpherical <- function(loc,degrees=T){
  #
  # conversion between Cartesian and spherical coordinate systems
  #
  # INPUT:
  #
  # loc      position in Cartesian coordinates c(x,y,z), 
  # degrees  logical, TRUE if want lat,lon converted to degrees
  #
  # OUPUT:
  #
  # a vector c(r,lat,lon) of the corresponding spherical coordinates
  #
  
  r <- sqrt(loc[1]**2+loc[2]**2+loc[3]**2)
  theta <- acos(loc[3]/r)
  lon <- atan2(loc[2],loc[1])
  
  lat = pi/2-theta
  locsph <- c(r,lat,lon)
  if(degrees) locsph[2:3] <- locsph[2:3]*180/pi
  
  return(locsph)
  
} # cartesianToSpherical

#Read in data outputted from 1_select_beam_geometry.py
indir <- './output/'

outdir <- indir

txfile <- 'TX_ECEF.txt'
rxfile <- 'RX_ECEF.txt'
ptsfile <- 'scatterpoints_ECEF.txt'
txd <- read.csv(paste(indir,txfile,sep=''))
rxd <- read.csv(paste(indir,rxfile,sep=''))
ptd <- read.csv(paste(indir,ptsfile,sep=''))
npts <- dim(ptd)[1]

#Select 1st row, 3 columns
ecefcols <-c('xECEF','yECEF','zECEF')
ptd[1,ecefcols]

# Load location of transmitter (assuming there is one transmitter)
locTrans <- as.numeric(c(txd[1,ecefcols],use.names=FALSE))

# Load locations of receivers (assuming there are three receivers)
locRec <- list(as.numeric(c(rxd[1,ecefcols],use.names=FALSE)),
              as.numeric(c(rxd[2,ecefcols],use.names=FALSE)),
              as.numeric(c(rxd[3,ecefcols],use.names=FALSE)))

locTrans2 <- cartesianToSpherical(locTrans)[2:3]
locRec2 <- lapply(locRec,function(x) {cartesianToSpherical(x)[c(2,3)]})

RE = EarthRadius()
pts <- apply(ptd[,ecefcols],1,function(x) cartesianToSpherical(as.numeric(c(x,use.names=F))))
alts <- pts[1,]-RE

ptsxy <- apply(pts[2:3,],2,function(x) sphericalToPlanar.geographic(x[1],x[2],zeroLatitude=locTrans2[1],zeroMeridian=locTrans2[2]))
x <- unlist(lapply(ptsxy,function(x) as.numeric(x$x)),recursive=F)
y <- unlist(lapply(ptsxy,function(x) as.numeric(x$y)),recursive=F)

#Ionosphere properties, used for self-noise calculation
fwhmIonSlab <- 100 # ionospheric slab thickness [km]
Ne <- 1e12         # electron density [m^-3]

#Additional radar system properties
fradar <- 233e6  # Radar frequency [Hz]
tau0 <- 100      # ACF time-scale [us] (IS THIS REASONABLE?)
dutyCycle <- .25  # Transmitter duty cycle
RXduty <- 0.8      # Receiver duty cycle
Tnoise <- c(200) # Noise temperature for receiver sites
Pt <- 5.e6      # Devin says " The goal for the first stage implementation of the system is for 5 MW TX power (https://eiscat.se/eiscat3d-information/eiscat_3d-faq/)"

mineleTrans <- 30
mineleRec <- 30
fwhmTrans <- 1
fwhmRec <- 1
fwhmRange <- 3

#Ilkka says he used gamma_0 = 0.01 in creating the error tables, so...
gamma0 <- 0.01

dat = multistaticNoiseLevels(locTrans2,locTrans2,locRec2,locxy=F,
                       fwhmTrans=fwhmTrans,fwhmRec=fwhmRec,fwhmRange=fwhmRange,
                       x=x,y=y,heights=alts,
                       infinity=defaultInfinity(),
                       Pt=Pt,
                       Ne=Ne,
                       fwhmIonSlab=fwhmIonSlab,
                       Tnoise=Tnoise,
                       fradar=fradar,
                       tau0=tau0,
                       phArrTrans=T,
                       phArrRec=T,
                       RXduty=RXduty,
                       verbose=FALSE,
                       mineleTrans=mineleTrans,
                       mineleRec=mineleRec,
                       ptlist=T)

write.table(dat$noiseLevel.isotropic$noiseLevel/gamma0,
            paste(outdir,"isotropicnoise.txt",sep=''),
          row.names=FALSE,
          col.names=FALSE)
write.table(dat$noiseLevel.velocity$noiseLevel/gamma0,
          paste(outdir,"velocitynoise.txt",sep=''),
          row.names=FALSE,
          col.names=FALSE)
