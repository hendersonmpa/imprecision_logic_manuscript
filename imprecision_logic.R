library("tidyverse")
     library("lubridate")
   #  library("readxl")
     library("RODBC")
     library("xtable")
library("pander")
library("scales")
     options(warn=-1) ## options(warn=0) to turn back on
     ## Suppress summarise info
     today <- as.Date(now())
     source("credentials.r")

     ## rescale a vector from 0 to 1
     rescale <- function(x){
       (x-min(x))/(max(x)-min(x))
     }

     '%!in%' <- function(x,y)!('%in%'(x,y))

     ### accept data, initial and confirm thresholds
     ### return the area of the probability density polygon 
     densprob <- function(dens, lower, upper) {
       x <- dens$x
       y <- dens$y
       dx <- x[2] - x[1] ## determine the increment
       C <- sum(y) * dx ## total area should be very close to 1
       p.unscaled <- sum(y[x >= lower & x <= upper]) * dx 
       round(p.unscaled/C, digits = 5) ## scaled probablity
     }


   ## Calculate the joint probability of the sample distribution and the imprecsion distribution for each y from the initial threshold to 6 SD
   jointprob  <- function(dens, confirm, sd, lower, upper) {
     #dens <- density(pop_data)
     x <- dens$x
     y <- dens$y
     dx <- x[2] - x[1] ## determine the increment
     pop_dens_region <- y[x >= lower & x <= upper] ##trim the pop dens to the region of interest

     ## create the imprecision region
     start <- confirm - (6*sd)
     stop <- confirm + (6*sd)
     x2 <- seq(start,stop,dx)
     y2 <- dnorm(x2,confirm,sd)
     imp_dens_region <- y2[x2 >= lower & x2 <= upper] ##trim the imprecision dens to the region of interest
     #Create a dataframe with the Ys from both densities side by side
     sum(pop_dens_region * imp_dens_region) * dx
     }

  ## Testing
## jointprob(galtfilter$result, 1.5, 2.7, 3.2)
## jointprob(galtfilter$result, 1.5, 1.5, 3.2)

  ## accept data, confirmation threshold, sd at the threshold, factor expansion factor
     ## return factor, lower, upper, grey area samples, uncertain area samples

  denssamples <- function(dens, confirm, sd, factor , direction = "left", samples = 145000) {
       #dens <- density(data)
       umsd  <- factor * sd
       sevensd  <- 7 * sd
       if (direction == "left") {
	 ## initial threshold based on the sd and factor
	 initial <- confirm + umsd
	 end <- confirm + sevensd
	 ## grey area between the confirm and initial thresholds
	 grey_area <- densprob(dens, confirm, initial)
	 ## Calculate the joint probability of the sample distribution and the imprecsion distribution for each value from the initial threshold to 6 SD
	 imprecision_area <- jointprob(dens, confirm, sd, initial, end)

       } else {
	 ## right sided threshold
	 print("Right sided thresholds not implemented")
       }
       ## area of the probability density polygon between the initial and 6 sd above
       grey_samples <- grey_area * samples
       imprecision_samples <- imprecision_area * samples
       list(factor, initial, grey_samples, imprecision_samples)
   }

  ## Testing
  ##  denssamples(galtfilter$result, 1.5, 0.2, 1, direction = "left")

galtquery <- "select s.spcextcode1 as accession,
	 a.ansTimeMeasured as measured_time,
	 s.spcExtcode2 as form,
	 sd.sd2GestationAge as ga,
	 sd.sd2Weight as bw,
	 sd.sd2AgeAtCollection as aoc,
	 a.ansvalueplain as result,
	 va.ResultCode as result_code
	 from (select s.specimenid, a.testid, max(answerix) as answerindex
	 from Answer a inner join specimen s on s.SpecimenID = a.SpecimenID
	 where a.TestId = 13 
	 and a.ansStatus = 110
	 and s.spcextcode1 like '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'
	 and substring(s.spcextcode1,1,8) between '20170000' and '20190000'
	 and substring(s.spcextcode1,9,1) not in ('4', '7', '8')
	 group by s.specimenid, a.TestId) a1
	 inner join answer a on a1.SpecimenID = a.SpecimenID and a1.AnswerIndex = a.AnswerIX and a1.TestId = a.TestId
	 inner join specimen s on a1.specimenid = s.specimenid
	 inner join vw_Answers va on s.spcExtcode1 = va.AccessionNumber and a.TestId = va.TestID
	 inner join specimendetail2 sd on sd.SpecimenId = va.SpecimenID
	 order by s.spcextcode1"
## galtdata <- with_con(galtquery)
## write.csv(galtdata, file= paste0("./data/galt_data_", today, ".csv"))
galtdata <- read.csv("./data/galt_data_2022-04-26.csv", stringsAsFactors = FALSE)
galtdata$measured_time  <- ymd_hms(galtdata$measured_time)
galtdata <- na.omit(galtdata)
galtfilter <-  galtdata %>%
  filter( !result_code %in% c("GALT-C-01-100", "GALT-C-01-001", "GALT-C-01-012")) ## initial results only
  #filter(measured_time >= ymd_hms("2018-06-11 00:00:00") & !result_code %in% c("GALT-C-01-100", "GALT-C-01-001", "GALT-C-01-012")) ## initial result only

annual_volume <- 145000
start <- 0
end <- 6
mean <- 1.5
sd <- 0.2

## initialize the dataframe
galtarea <- data.frame(factor = double(), initial = double(),
		       grey = double(), imprecision = double(),
		       stringsAsFactors = FALSE)

## stochastic simulation
c <- 0
for (s in 1:1000) {
  data_sample <- sample(galtfilter$result, size = annual_volume, replace = TRUE)
  d <- density(data_sample)
  for (i in seq(from = start, to = end, by = 1)) {
    c <- c + 1
    galtarea[c,] <- denssamples(d, confirm = 1.5, sd = 0.2, factor = i)
    # denssamples(d, confirm = 1.5, sd = 0.2, factor = i)
  }
}

galtarea %>%
  group_by(factor, initial) %>%
  summarise(grey_p025 = quantile(grey,probs = c(0.025), type = 8, na.rm = TRUE),
	    grey_median = median(grey, na.rm = TRUE),
	    grey_p975 = quantile(grey,probs = c(0.975), type = 8, na.rm = TRUE),
	    imp_p025 = quantile(imprecision,probs = c(0.025), type = 8, na.rm = TRUE),
	    imp_median = median(imprecision, na.rm = TRUE),
	    imp_p975 = quantile(imprecision,probs = c(0.975), type = 8, na.rm = TRUE)) %>%
  ## rename("Standard Deviations" = factor,
  ##       	    "Initial Threshold" = initial,
  ##       	    "Grey Zone" = grey,
  ##       	    "Imprecision Zone" =  imprecision) %>%
  xtable(caption = "Initial Threshold Simulation Results. In each simulation the confirmation threshold is set to 1.5 U/g Hb and the initial thresholds is increased by the corresponding number of standard deviations",
	 label = "tab:imprecision", display = c("d", "d", "f","f", "f", "f", "g", "g","g")) %>%
  print(include.rownames = FALSE)

##pandoc.table(style = "rmarkdown", caption = "Initial Threshold Simulation Results. In each simulation the confirmation threshold is set to 1.5 U/g Hb and the initial thresholds is increased by the corresponding number of standard deviations")

dens <- density(galtfilter$result)
  sd <- 0.2 ##SD at postive confirm
  confirm  <- 1.5
  initial <- confirm + (1.1*sd)
  theight  <- max(dens$y[which(dens$x <= confirm)])
  bheight  <- max(dens$y[which(dens$x <= initial)])
  ## defining the region of FN uncertainty
  start  <- confirm - (6*sd)
  stop <- confirm + (6*sd)
  x2 <- seq(start,stop,0.01)
  y2 <- theight*rescale(dnorm(x2,confirm,sd))
  ## create indices for half of the UM distribution
  halfx2 <- seq(confirm,stop,0.01)
  halfy2 <- y2[length(halfx2):length(x2)]
  fnx2 <- seq(initial,stop,0.01)
  fny2 <- y2[(length(x2) - length(fnx2)):(length(x2) -1)]

plot(x= 0:2*confirm, y = 0:2*bheight, type = "n",
       xlab = "GALT Activity U/g Hb",
       ylab = "Probability Density")
### polygons
polygon(dens,col = "steelblue", border = "steelblue")
## imprecision zome
# purple50 <- adjustcolor("purple", alpha.f = 0.25)
with(dens, polygon(x=c(stop, stop, x[x < stop]), y=c(0, y[x=stop], y[x < stop]), col="goldenrod", border = "goldenrod"))
## grey zone
with(dens, polygon(x=c(initial, initial, x[x < initial]), y=c(0, y[x=initial], y[x < initial]), col="grey75", border = "grey75"))
## positive
  with(dens, polygon(x=c(confirm, confirm, x[x < confirm]), y=c(0, y[x=confirm], y[x < confirm]), col="black", border = "black"))

  ## measurement error distribution
  points(x2,y2,type="l",col="red", lwd = 4) ## region of uncertainty of measurment
  zeros <- rep(0,length(x2)) # create a vector of zeros
  #polygon(c(x2,rev(x2)),c(y2,zeros), border = NA, col="red")

  #polygon(c(halfx2,rev(halfx2)),c(halfy2,zeros), border = NA, col="red")
  fnzeros <- rep(0,length(fnx2)) # create a vector of zeros
  polygon(c(fnx2,rev(fnx2)),c(fny2,fnzeros), border = NA, col="red")
  #area <- 0.01 * sum(halfy2)
  #samples <- round(area *145000, digits = 0)
  #text(x = 0.55, y = 0.004, label= paste("Annual results in red area:",samples), side = 3)

  abline(v = confirm, col = "black" , lty = 1, lwd = 2)
  abline(v = initial, col = "black", lty = 2, lwd = 2)
  #abline(v = confirm + (1*sd), col = "black", lty = 2, lwd = 2) 

  legend("topleft",
	 legend = c("positive", "grey zone", "imprecision zone", "negative", 
		     "analytical imprecision","confirmation threshold",
		    "initial threshold"),
	 col = c("black", "grey75", "goldenrod", "steelblue" , "red", "black", "black"),
	 lty = c(NA, NA, NA, NA, "solid", "solid", "dashed"),
	 lwd = c(NA, NA, NA, NA, 2, 2, 2),
	 pch = c(15, 15, 15, 15, NA , NA, NA))

dens <- density(galtfilter$result)
sd <- 0.2 ##SD at postive confirm
confirm  <- 1.5
initial <- confirm + (1.1*sd)
theight  <- max(dens$y[which(dens$x <= confirm)])
bheight  <- max(dens$y[which(dens$x <= initial)])
## defining the region of FN uncertainty
start <- confirm - (6*sd)
stop <- confirm + (6*sd)
x2 <- seq(start,stop,0.01)
					#y2 <- theight*rescale(dnorm(x2,confirm,sd))
y2 <- dnorm(x2,confirm,sd)

## create indices for half of the UM distribution
halfx2 <- seq(confirm,stop,0.01)
halfy2 <- y2[length(halfx2):length(x2)]
fnx2 <- seq(initial,stop,0.01)
fny2 <- y2[(length(x2) - length(fnx2)):(length(x2) -1)]

plot(dens,type = "n", main = "",
     xlim = c(0,5),
     ylim = c(0,0.04),
					#x= 0:10*confirm, y = 0:10*bheight, 
     xlab = "GALT Activity U/g Hb",
     ylab = "Probability Density")
red50 <- alpha("red", 0.5)
### polygons
## measurement error distribution
polygon(x2,y2,type="l",col=red50, border = "red")
polygon(dens,col = "steelblue", border = "steelblue")
## imprecision zome
					# purple50 <- adjustcolor("purple", alpha.f = 0.25)
with(dens, polygon(x=c(stop, stop, x[x < stop]), y=c(0, y[x=stop], y[x < stop]), col="goldenrod", border = "goldenrod"))
## grey zone
with(dens, polygon(x=c(initial, initial, x[x < initial]), y=c(0, y[x=initial], y[x < initial]), col="grey75", border = "grey75"))
## positive
with(dens, polygon(x=c(confirm, confirm, x[x < confirm]), y=c(0, y[x=confirm], y[x < confirm]), col="black", border = "black"))
points(x2,y2,type="l",col="red", lwd = 2) ## region of uncertainty of measurment



##polygon(c(fnx2,rev(fnx2)),c(fny2,fnzeros), border = NA, col= alpha("red", 0.5))
#fnzeros <- rep(0,length(fnx2)) # create a vector of zeros
#polygon(c(fnx2,rev(fnx2)),c(fny2,fnzeros), border = NA, col= alpha("red", 0.5))
					# thresholds
abline(v = confirm, col = "black" , lty = 1, lwd = 2)
abline(v = initial, col = "black", lty = 2, lwd = 2)
legend("topright", inset = 0.02,
       legend = c("positive", "grey zone", "imprecision zone", "negative", 
		  "analytical imprecision","confirmation threshold",
		  "initial threshold"),
       col = c("black", "grey75", "goldenrod", "steelblue" , red50, "black", "black"),
       lty = c(NA, NA, NA, NA, NA, "solid", "dashed"),
       lwd = c(NA, NA, NA, NA, NA, 2, 2),
       pch = c(15, 15, 15, 15, 15 , NA, NA),
       bg = "white", 
       box.lty = 0)
