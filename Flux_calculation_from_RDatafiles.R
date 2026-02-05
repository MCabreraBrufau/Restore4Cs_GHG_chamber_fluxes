# ---
# Authors: Miguel Cabrera-Brufau & Camille Minaudo
# Project: "RESTORE4Cs"
# ---

#TO-DO------
#Functionally works, check naming and comments (format pretty).
#Sanity check: compare fluxes produced here with those published in LifeWatch

# --- Description----
#This script uses the goFlux pacakge to produce (from raw GHG concentration timeseries and CO2 and CH4 auxfiles) flux estimates (total.flux, LM.flux, HM.flux). Full calculation of RESTORE4Cs dataset takes ~35 minutes.

#Input files: 
  #Rdata files (named after each sampling day), containing CO2 and CH4 concentration timeseries
  #CH4 auxfile, containing UniqueID, start-time, duration and initial temperature and pressure of static chamber incubations
  #CO2 auxfile, containing UniqueID, start-time, duration and initial temperature and pressure of static chamber incubations

#Output files:
  #PDF files with incubation plots (for CO2 and CH4 separately) for each sampling day
  #co2_fluxes.csv: csv file with all CO2 flux estimates (3 estimates + ancillary data per incubation)
  #ch4_fluxes.csv: csv file with all CH4 flux estimates (3 estimates + ancillary data per incubation)


#3 flux estimate methods for each incubation: LM, HM (goFlux package) and total.flux (custom). 
#total.flux is calculated using the mean concentration difference (first vs last 10s, and the elapsed time between these means). This ensures an more appropriate estimate of total GHG flux when non-linear patterns caused by ebullition exist (Which make LM and HM flux estimates highly unreliable). total.flux estimate is accompanied by SE (standard error) for consistency with goFLux LM and HM estimates.


#Clear Global Environment
rm(list = ls())


# ---- Directories ----

#Root path
dropbox_root <- "C:/Users/Miguel/Dropbox/GHG_bestflux_calculation" # You have to make sure this is pointing to the right folder on your local machine

#Path to Rdata files, containing CO2 and CH4 timeseries (UTC time) named per sampling event (subsite) and gas analyzer (one duplicated sampling event where two instruments were used: "S2-RI-P2_LI-7810" and "S2-RI-P2_Picarro", total of 145 files, but LI-8710 only has data from extra incubations for inter-instrument calibration, not used)
RData_path <- paste0(dropbox_root, "/RData_timeseries") 

#Path to CO2 and CH4 auxfiles with corrected start.time and duration (plus initial pressure and temperature)
auxfile_path<- paste0(dropbox_root,"/Auxfiles/") 

#ORIGINAL auxfiles used for aquaGHG calculation
auxfile_path<- paste0("C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/GHG/Working data/Auxfiles_4aquaGHG/") 

#Set results_path for computed fluxes and incubation plots
results_path <- paste0(dropbox_root,"/Results/")
#Create subdirectory for plots
plots_path<- paste0(results_path,"Incubation_plots/")
if (!dir.exists(plots_path)) {
  dir.create(plots_path, recursive = TRUE)
}


#goFlux Update----
#Check goFlux version: and warn if not up-to-date 
{desc <- packageDescription("goFlux")

if (is.null(desc$RemoteSha)) {
  warning("goFlux was not installed from GitHub; cannot check latest version.")
} else {
  installed_sha <- desc$RemoteSha
  latest_sha <- gh::gh(
    "/repos/{owner}/{repo}/commits/{ref}",
    owner = "Qepanna",
    repo  = "goFlux",
    ref   = "HEAD"
  )$sha
  
  if (installed_sha != latest_sha) {
    stop(
      "goFlux is out of date.\n",
      "Installed: ", installed_sha, "\n",
      "Latest:    ", latest_sha, "\n",
      "Please update intentionally."
    )
  }else{
  message("OK: Latest goFlux version is already installed with sha: ", installed_sha)
    rm(installed_sha,latest_sha,desc)}
}
}
#For reproducibility: this script was developed using the goFlux sha: "67c276d87d984d55b70d7b756ad702d63dbf2038"


#Decide wheter to update goFlux explicitly: 
update_goFlux<- F

# if(update_goFlux==T){
# try(detach("package:goFlux", unload = TRUE), silent = TRUE)
# devtools::install_github("Qepanna/goFlux")}



#Load PKGs and functions-----
library(goFlux)
library(dplyr)
library(tidyr)
library(ggplot2)
# library(lubridate)
# library(pbapply)
# library(egg)
# library(purrr)
# library(ggnewscale)
# library(stringr)
# library(data.table)


#Function to load subsite data ()
load_subsite<- function(auxfile, RData_path){

gas <- unique(auxfile$gas_analiser)
subsite<- unique((auxfile$subsite))
setwd(RData_path)
if(gas== "LI-COR"){
  gs_suffix <- "LI-7810"
} else {
  gs_suffix <- gas
}
load(file = paste0(subsite,"_",gs_suffix,".RData"))
mydata <- mydata[,c("POSIX.time", 
                    "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                    "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
return(mydata)
}




#Create function to calculate average total flux (tf-t0 divided by elapsed time between them, using average initial and final concentrations with 10s window), using same dataframe input as goFlux function plus the output of the goFlux call itself (to get flux.term)
total.flux <- function(dataframe, #dataframe passed to goFlux call
                       goflux_result, #output from goFlux call (to get flux.term)
                       gastype="CO2_ppm" #name of gas column for which to calculate
                       ) {
  
  t.win<-10 #time window for averaging initial and final concentration 
  
  result.total.flux <-dataframe %>%
    # 1. Group calculations by UniqueID
    group_by(UniqueID) %>%
    # 2. Compute concentrations and elapsed time for each UniqueID using the gastype column in the arguments
    summarise(
      # Initial concentration (first t.win seconds)
      C0.10s = mean(.data[[gastype]][Etime < t.win], na.rm = TRUE),
      C0.10s.sd = sd(.data[[gastype]][Etime < t.win], na.rm = TRUE),
      # Final concentration (last t.win seconds)
      Cf.10s = mean(.data[[gastype]][Etime > max(Etime) - t.win], na.rm = TRUE),
      Cf.10s.sd = sd(.data[[gastype]][Etime > max(Etime) - t.win], na.rm = TRUE),
      # Elapsed time for total flux calculation
      elapsed_time = max(Etime) - t.win,
      .groups = "drop") %>%
    # 3. Get flux.term from goFlux output
    left_join(goflux_result %>% select(UniqueID, flux.term), by = "UniqueID") %>%
    # 4. Calculate delta concentration and errors
    mutate(
      # Change in concentration
      delta.conc = Cf.10s - C0.10s,
      # SD of delta concentration (error propagation)
      # delta.conc.sd = sqrt(C0.10s.sd^2 + Cf.10s.sd^2),
      # SE of delta concentration
      delta.conc.se = sqrt((C0.10s.sd^2 / t.win) +
                             (Cf.10s.sd^2 / t.win)),
      # Total flux calculation
      total.flux = (delta.conc / elapsed_time) * flux.term,
      # SD of total flux
      # total.flux.sd = abs(total.flux) * delta.conc.sd / abs(delta.conc),
      # SE of total flux
      total.flux.se = (flux.term / elapsed_time) * delta.conc.se
    ) %>%
    # 5. Select final output columns
    select(
      UniqueID,
      C0.10s, C0.10s.sd,
      Cf.10s, Cf.10s.sd,
      total.flux, total.flux.se
    )
  return(result.total.flux)
}


#Alternative function to load each incubation independently (adding required auxfile parameters).
#Load function
#This loads exactly the same data as the funtion in aquaGHG script

load_incubation <- function(auxfile_i, RData_path){
  
  # message("Loading data for ",auxfile_i$UniqueID)
  gas <- unique(auxfile_i$gas_analiser)
  
  setwd(RData_path)
  if(gas== "LI-COR"){
    gs_suffix <- "LI-7810"
  } else {
    gs_suffix <- gas
  }
  load(file = paste0(auxfile_i$subsite,"_",gs_suffix,".RData"))
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  
  mydata <- mydata[which(mydata$POSIX.time>=auxfile_i$start.time & 
                           mydata$POSIX.time<=auxfile_i$start.time+auxfile_i$duration),]
  if(dim(mydata)[1]>0){
    #Add auxfile data to rawdata of incubation:
    mydata$UniqueID <- auxfile_i$UniqueID
    mydata$flag=1 #Indicates goflux that all data is part of incubation
    mydata$Etime <- as.numeric(mydata$POSIX.time) - min(as.numeric(mydata$POSIX.time))
    mydata$Area<-auxfile_i$Area
    mydata$Vtot<-auxfile_i$Vtot
    mydata$Pcham<-auxfile_i$Pcham
    mydata$Tcham<-auxfile_i$Tcham   
  } else {
    warning(paste0("For ", auxfile_i$UniqueID, ", no measurements were found!"))
  }
  return(mydata)
}



  
  
# ---- Loading auxfiles ----


#TO-delete------


#Exact copy from aquaGHG script: 
#CO2 auxfile:
co2_auxfile <- read.csv(file = paste0(auxfile_path,"co2_auxfile.csv"))

co2_auxfile <- co2_auxfile %>%
  mutate(start.time = as.POSIXct(start.time,tz = "UTC"))#original start.time (not rounded)

# co2_auxfile <- co2_auxfile %>%
  # mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"secs")),#rounding transforms to class different to POSIXct, causing error, forze as.POSIXct to fix
#          obs.length=floor(duration))

#CH4 auxfile:
# ch4_auxfile <- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))
# 
# ch4_auxfile <- ch4_auxfile %>% 
#   mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"secs")),#rounding transforms to class different to POSIXct, causing error, forze as.POSIXct to fix
#          obs.length=floor(duration))



#Using corrected (rounded) auxfiles
# 
# #CO2 auxfile:
# co2_auxfile <- read.csv(file = paste0(auxfile_path,"correct_co2_auxfile.csv"))
# 
# co2_auxfile <- co2_auxfile %>% 
#   mutate(start.time=as.POSIXct(start.time,tz = "UTC"),
#          obs.length=floor(duration),
#          #Add sampling (subsite_analiser Identifier)
#          analiser_sufix=if_else(gas_analiser=="LI-COR","LI-7810",gas_analiser),
#          sampling=paste(subsite, analiser_sufix,sep = "_"))%>% 
#   select(-analiser_sufix)
# 
# 
# #CH4 auxfile:
# ch4_auxfile <- read.csv(file = paste0(auxfile_path,"correct_ch4_auxfile.csv"))
# 
# ch4_auxfile <- ch4_auxfile %>% 
#   mutate(start.time=as.POSIXct(start.time,tz = "UTC"),
#          obs.length=floor(duration),
#          #Add sampling (subsite_analiser Identifier)
#          analiser_sufix=if_else(gas_analiser=="LI-COR","LI-7810",gas_analiser),
#          sampling=paste(subsite, analiser_sufix,sep = "_")) %>% 
#   select(-analiser_sufix)






#-----subset to test (Optional)------

#select only a couple of subsites to test the script
pattern<- ""

# ch4_auxfile<- ch4_auxfile %>% 
#   filter(grepl(pattern, UniqueID))

co2_auxfile<- co2_auxfile %>% 
  filter(grepl(pattern, UniqueID))

#Calculate fluxes------


##CO2 loop per incubation----
#Adapt code in "recalculate fluxes with aquaGHG_miguel_edits", loading each incubation in the same way that is done there, but using only goflux functions (then re-check flux results). 
#Use load_incubation function to create data frame to feed into goflux (including auxfile details Pcham, Tcham, Area,Vtot,...), goflux function only takes a single data-frame full results should be included.
#Loop over subsites, then over incubations for flux calculation, print per-subsite plot (after creating plotlist and appending each incubation including its UniqueID as ggtitle)

for (s in unique(co2_auxfile$subsite)){
  
  #Print subsite progress:
  message("Now processing CO2 from subsite ",s," (",which(unique(co2_auxfile$subsite)==s)," of ",length(unique(co2_auxfile$subsite)),")")
  
  #subset subsite auxfile
  auxfile_co2_s<- co2_auxfile %>% filter(subsite==s)
  
  #Initiate pdf file to save plots from subsite s
  
  # plot_filename <- paste("CO2",s, as.character(as.Date(last(auxfile_co2_s$start.time))),sep="_")
  # pdf(file = paste0(plots_path,plot_filename,".pdf"),
  #     width = 8,height = 8)
      
  #Loop over incubations within subsite s, calculating fluxes and producing plots 
  for (i in auxfile_co2_s$UniqueID){
    
    #get Incubation auxfile
    auxfile_co2_i<- auxfile_co2_s %>% filter(UniqueID==i)
    
    #Print Overall progress (%):
    incubnum<- which(co2_auxfile$UniqueID==i)
    message(paste0("Processing CO2 incubation ",incubnum, " of ",length(co2_auxfile$UniqueID)," (",round(100*incubnum/length(co2_auxfile$UniqueID),0), " %)"))
    
    #Load Incubation rawdata (combining it with auxfile parameters)
    mydata_i<- load_incubation(auxfile=auxfile_co2_i, RData_path = RData_path)
    
    #Use goFlux function to obtain the calculated flux for LM and HM.
    #Call structure to avoid console printouts (Warnings and progress bars)
    co2_goflux_i <- local({
      invisible(capture.output(
        res <- suppressWarnings(
          suppressMessages(
            goFlux(mydata_i, "CO2dry_ppm")
          )
        )
      ))
      res
    })
    
    #Calculate total.flux using function defined above
    co2_total.flux_i<- total.flux(dataframe = mydata_i,
                                  goflux_result = co2_goflux_i,
                                  gastype = "CO2dry_ppm")
    
    #Join goflux output (co2_goflux_i) and total.flux output (co2_total.flux_i)
    all_co2flux_i<- auxfile_co2_i %>% select(UniqueID) %>%
      left_join(co2_goflux_i, by="UniqueID") %>% 
      left_join(co2_total.flux_i, by="UniqueID")
    
    #Make sure all HM. parameters are numeric (even when NA), otherwise flux.plot throws error
    co2_goflux_i<- co2_goflux_i %>% 
      mutate(across(starts_with("HM."), as.numeric))
    
    #Create UniqueID plot (adding ggtitle with UniqueID)
    
    # #Create plot (silent)
    # testplot <- local({
    #   invisible(capture.output(
    #     res<- suppressWarnings(flux.plot(
    #       flux.results = co2_goflux_i,
    #       dataframe = mydata_i,
    #       gastype = "CO2dry_ppm",
    #       shoulder = 10,
    #       plot.display = c("MDF", "prec", "nb.obs","flux.term"),
    #       best.model = FALSE,
    #       quality.check = FALSE
    #     ))
    #   ))
    #   res
    # })
    # 
    # #print plot (avoiding printing of "[[1]]" in console)
    # invisible(capture.output(print(testplot))) 
    # #ADD the uniqueID to the displayed plot (to save it in pdf): 
    # grid::grid.text(as.character(i),
    #                 x = 0.5, y = 0.975,
    #                 gp = grid::gpar(fontsize = 14, fontface = "bold"))
    # 
    #Append results from each UniqueID, 
    if(i==co2_auxfile$UniqueID[1]){
      all_co2flux<- all_co2flux_i #Create general data-frame in first loop execution
    }else{
      #For subsequent runs, append subsite results to general results
      all_co2flux<- rbind(all_co2flux, all_co2flux_i)
    }
    
  } #end of per-uniqueID loop
  
  # dev.off() #Close and save pdf file with plots
  
}#end of per-subsite loop

#Save all_co2_flux:
write.csv(all_co2flux, file = paste0(results_path,"co2_fluxes.csv"),row.names = F)


#Check flux differences-----

#Load objects created by aquaGHG run
load("C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/Data/GHG/Processed data/computed_flux/CO2_aquaGHG.Rdata")


old_long<-CO2_flux.auto %>% 
  select(UniqueID,LM.C0,LM.Ct,nb.obs,C0,Ct,flux.term,LM.flux,HM.flux,mean.total.flux) %>%
  rename(total.flux=mean.total.flux) %>% 
  pivot_longer(cols = -UniqueID, names_to = "name",values_to = "old")


new_co2<- read.csv(file=paste0(results_path,"co2_fluxes.csv"))

new_long<- new_co2 %>%
  select(UniqueID,LM.C0,LM.Ct,nb.obs,C0,Ct,flux.term,LM.flux,HM.flux,total.flux) %>% 
  pivot_longer(cols = -UniqueID, names_to = "name",values_to = "new")


comparison_long<- new_long %>% 
  left_join(old_long) %>% 
  left_join(co2_auxfile %>% select(UniqueID,gas_analiser,duration)) %>% 
  mutate(equal=old==new,
         dif=new-old,
         percent_dif=dif/old*100,
         duration=duration+1)

#check nb.obs new
comparison_long %>% 
  filter(name=="nb.obs") %>% 
ggplot(aes(y=duration-new,fill=gas_analiser))+geom_histogram()
#For Li-COR (resolution==1.0s), same duration and new nb.obs, 
#ISSUE----
#a few Li-COR incubations with strong discrepancy duration vs new.... WHY???

#For Picarro (resolution >1.0s), duration is always higher than nb.obs new, OK, makes sense
#For los Gaots(resolution <1.0s), duration is either the same or higher lower than nb.obs (rouding error, ok)

#Old script interpolated concentration data to produce same nb.obs than duration 
comparison_long %>% 
  filter(name=="nb.obs") %>% 
  ggplot(aes(y=duration-old,fill=gas_analiser))+geom_histogram()

comparison_long %>% 
  filter(name=="nb.obs") %>%filter(duration-old<0) #One weird incubation with 1 extra nb.obs old 

#Old vs new initial and final concentrations (C0,Ct): not the same concentrations (due to old script interpolation)
comparison_long %>% 
  filter(name=="C0") %>% 
  ggplot(aes(y=old-new,fill=gas_analiser))+geom_histogram()

comparison_long %>% 
  filter(name=="Ct") %>% 
  ggplot(aes(y=old-new,fill=gas_analiser))+geom_histogram()

#Are the differences in LM.flux correlated to those of nb.obs? (i.e. does the "inflation" of nb.obs lead to unreliable fluxes?
comparison_long %>% 
  filter(name%in%c("nb.obs","LM.flux")) %>%
  select(UniqueID,name,dif,percent_dif,new) %>% 
  pivot_wider(names_from = name, values_from = c(dif,percent_dif,new)) %>% 
  filter(new_LM.flux >0.05) %>% #exclude very low fluxes that would inflate percent difference
  ggplot(aes(x=percent_dif_nb.obs, y=percent_dif_LM.flux,col=new_LM.flux<0.1))+
  geom_point()+ggtitle("percent differneces")
#there is no correlation, this means that differences in fluxes are not dependent on nb.obs differences. Fluxes are different due to interpolation changing the GHG concentrations (slightly), nb.obs is not used as the "time" in regressing concentration vs time.

#New fluxes differ by  ~7% (maximum) from old ones BUT the discrepancy does not scale with difference in nb.obs 

#Check flux term differences: 
comparison_long %>% 
  filter(name=="flux.term") %>% 
  ggplot(aes(y=old-new,fill=gas_analiser))+geom_histogram()

#Maximum relative difference in flux.term is <0.03%....OK, potentially arising from  
comparison_long %>% 
  filter(name=="flux.term") %>%
  summarise(max(percent_dif))

#Check the maximum percent difference (old-new/new) for the key parameters:
comparison_long %>% 
  group_by(name) %>% 
  summarise(maxpercent_dif=max(percent_dif,na.rm = T),
            meanpercent_dif=mean(percent_dif,na.rm=T))

comparison_long %>% 
  filter(name%in%c("flux.term","HM.flux","LM.flux","total.flux")) %>% 
  ggplot(aes(x=percent_dif))+
  geom_histogram()+
  facet_wrap(.~name,scales="free")


#incubations with large >3% LM.flux discrepancy
incub_LM_different<- comparison_long %>% 
  filter(name=="LM.flux") %>% filter(abs(percent_dif)>3) %>% pull(UniqueID)
#These incubations have a very high relative error for LM flux (discrepancy is within )
new_co2 %>% filter(UniqueID%in%incub_LM_different) %>% 
  arrange(abs(LM.se.rel)) %>% 
  select(LM.se.rel,LM.r2)

#incubations with large >3% total.flux discrepancy
incub_total.fluxdifferent<- comparison_long %>% 
  filter(name=="total.flux") %>% filter(abs(percent_dif)>3) %>% pull(UniqueID)
#These incubations have a very high relative error for LM flux (discrepancy is within )
new_co2 %>% filter(UniqueID%in%incub_total.fluxdifferent) %>% 
  arrange(abs(LM.se.rel)) %>% 
  select(LM.se.rel,LM.r2)

#incubations with large >3% HM.flux discrepancy
incub_HM.fluxdifferent<- comparison_long %>% 
  filter(name=="HM.flux") %>% filter(abs(percent_dif)>3) %>% pull(UniqueID)
#These incubations have a very high relative error for LM flux (discrepancy is within )
new_co2 %>% filter(UniqueID%in%incub_HM.fluxdifferent) %>% 
  arrange(abs(LM.se.rel)) %>% 
  select(LM.se.rel,LM.r2,g.fact)


#LM.flux discrepancies:
comparison_long %>% 
  filter(name=="LM.flux") %>% 
  filter(abs(new)>0.1) %>% 
  ggplot(aes(x=abs(new), y=percent_dif,col=gas_analiser))+
  geom_point()+
  geom_hline(yintercept = 2.37)+
  facet_wrap(.~gas_analiser,scales="free")
#LI-COR has very minor deviations, no bias. CAN BE USED AS-IS.
#Los Gatos have variable deviations with a bias of ca. 2.4% (new fluxes are 2.4% higher than old ones). mean time difference between data points: 0.9768881 s. If we assumed separation was exactly 1s, this biass would represent 2.31%, this might be the culprit. 

#Picarro has variable deviations, some very large, not limited to small fluxes. 


#LI-COR and Picarro are centered at cero. 
#Los-Gatos have a consistent biass of ~2.5 % (not terrible)
#OK incubations with large discrepancies old vs new are always associated to low fluxes

#HM.flux discrepancies:
#OK incubations with large discrepancies old vs new are always associated to low fluxes
comparison_long %>% 
  filter(name=="HM.flux") %>% 
  filter(abs(new)>0.1) %>% 
  ggplot(aes(x=log(abs(new)), y=percent_dif,col=gas_analiser))+
  geom_point()


comparison_long %>% 
  filter(name=="total.flux") %>% 
  filter(abs(new)>0.01) %>% 
  ggplot(aes(x=log(abs(new)), y=percent_dif,col=gas_analiser))+
  geom_point()




#Check instrument data frequency (in RData files)
#LI-COR:
licordata<-load_subsite(auxfile = co2_auxfile %>%
                          filter(subsite=="S1-CA-A1"),RData_path=RData_path)
first(licordata$POSIX.time)
unique(diff(licordata$POSIX.time)) #Always integer time difference (1s, or instrument down-time)


#Picarro: Non-constant time-frequency of data, non-integer frequency, non-integer POSIX.time, rounding to nearest second is not particularly elegant, but shoudl not cause too much problems, although it will create (i) duplicated POSIX.time values and (ii) gaps in POSIX.time. 

picarrodata<-load_subsite(auxfile = co2_auxfile %>%
                          filter(subsite=="S1-CA-P1"),RData_path=RData_path)
first(picarrodata$POSIX.time)
unique(diff(picarrodata$POSIX.time)) 
mean(diff(picarrodata$POSIX.time))
median(diff(picarrodata$POSIX.time))
picarrodata %>% 
  mutate(dif_time=c(NA,diff(.$POSIX.time))) %>% 
ggplot(aes(x=POSIX.time,y=dif_time))+geom_line()

picarrodata %>% 
  mutate(dif_time=c(NA,diff(.$POSIX.time))) %>% 
  ggplot(aes(y=dif_time))+geom_boxplot(outliers = F)

#Los Gatos: higher than 1s frequency, rounding will result in 1 duplicated timestamp from time to time. 
lgsdata<-load_subsite(auxfile = co2_auxfile %>%
                            filter(subsite=="S3-DU-P1"),RData_path=RData_path)

first(lgsdata$POSIX.time)
unique(diff(lgsdata$POSIX.time)) 
mean(diff(lgsdata$POSIX.time))
median(diff(lgsdata$POSIX.time))
lgsdata %>% 
  mutate(dif_time=c(NA,diff(.$POSIX.time))) %>% 
  filter(dif_time<9838) %>% #remove shutdown gap
  ggplot(aes(x=POSIX.time,y=dif_time))+geom_line()

lgsdata %>% 
  mutate(dif_time=c(NA,diff(.$POSIX.time))) %>% 
  filter(dif_time<9838) %>% 
  ggplot(aes(y=dif_time))+geom_boxplot(outliers = F)


exampledata<- data.frame(exact_time=seq(from=0, to=1500, by=0.977)) %>% 
  mutate(rounded_time=round(exact_time,digits = 0),
         rownumber=row_number()-1,
         ghg_slope1=2*exact_time)

ggplot(exampledata)+
  geom_line(aes(x=exact_time, y=ghg_slope1, col="exact_time"))+
  geom_line(aes(x=rounded_time, y=ghg_slope1, col="rounded_time"))+
  geom_line(aes(x=rownumber, y=ghg_slope1, col="rownum"))


#Slope for exact time: 
exact_lm<- lm(formula = ghg_slope1~exact_time, data = exampledata)


#Slope with rounded time
rounded_lm<- lm(formula = ghg_slope1~rounded_time, data = exampledata)

#Slope for rownumber
rownum_lm<- lm(formula = ghg_slope1~rownumber, data=exampledata)

#slope relative difference (exact vs rounded)
(exact_lm$coefficients[2]-rounded_lm$coefficients[2])/rounded_lm$coefficients[2]*100

#slope relative difference (runded vs rownum)
(rounded_lm$coefficients[2]-rownum_lm$coefficients[2])/rownum_lm$coefficients[2]*100
#2.354216% "overestimation" using rounded time instead of rownumber (assuming 1s per row)
#this matches the biass observed for Los Gatos. 


#_____-------







