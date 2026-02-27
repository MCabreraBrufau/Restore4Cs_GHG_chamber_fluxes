#Fluxes_from_RData

# ---
# Authors: Miguel Cabrera-Brufau & Camille Minaudo
# Project: "RESTORE4Cs"
# ---


# --- Description----
#This script uses the goFlux package to produce (from raw GHG concentration time-series and CO2 and CH4 auxfiles) three flux estimates (total.flux, LM.flux, HM.flux). 
#NO modification is made to the raw GHG concentration time-series (RData load -> filter for incubation timespan -> flux calculation)

#Input files: 
  #Rdata files (named after each subsite + gas_analizer), containing CO2, CH4 and H20 concentration time-series
  #CH4 auxfile
  #CO2 auxfile

#Output files:
  #PDF files with incubation plots (for CO2 and CH4 separately) for each subsite
  #all_co2flux.csv: csv file with all CO2 flux estimates (3 estimates + ancillary data per incubation)
  #all_ch4flux.csv: csv file with all CH4 flux estimates (3 estimates + ancillary data per incubation)

#Ancillary data per incubation in all_'ghg'flux.csv files is composed of parameters describing the structure of raw-data imported and of data used for flux calculation (incubation duration, total datapoints, and summary of delta-time between datapoints: average, sd, min and max Dt, in seconds).

#3 flux estimate methods for each incubation: LM, HM (goFlux package) and total.flux (custom).

#total.flux is calculated using the mean concentration difference (first vs last 10s, and the elapsed time between these means). This ensures an more appropriate estimate of GHG flux when non-linear patterns caused by ebullition exist (Which make LM and HM flux estimates highly unreliable). total.flux estimate is accompanied by SE (standard error) for consistency with goFLux LM and HM estimates.


#Clear Global Environment
rm(list = ls())


#Options------

#1st. subset auxfiles to run only the script for a subset of the incubations
#Add in pattern something to match in auxfiles UniqueID (via grepl function), leave as empty pattern<- "" to calculate full dataset
pattern<- ""

#2nd. Produce plots for fluxes?
#set as FALSE to skip plotting of fluxes into pdf files (much quicker calculation)
produce_plots<- T


#3rd Should the script be run for both co2 and ch4? 
#You can decide to run only co2 fluxes, only ch4 fluxes or both. 
calculate_gas<- c("co2","ch4") #can be c("co2"), c("ch4") or c("co2","ch4")



# ---- Directories ----
#Root path: You have to make sure this is pointing to the right folder on your local machine, by default, the repository folder is selected. If you want a different folder, make sure it contains the appropriate subfolders "Auxfiles" and "RData_timeseries" with the corresponding required data files.
root_path <-dirname(rstudioapi::getSourceEditorContext()$path)

#Path to Rdata files, containing CO2 and CH4 timeseries (UTC time) named per sampling event (subsite) and gas analyzer.
RData_path <- paste0(root_path, "/RData_timeseries") 

#Path to CO2 and CH4 auxfiles
auxfile_path<- paste0(root_path,"/Auxfiles/") 

#Set results_path for computed fluxes
results_path <- paste0(root_path,"/Results/")
#Create folder if it does not exist
if (!dir.exists(results_path)) {
  dir.create(results_path, recursive = TRUE)
}
#set plots_path for plots
plots_path<- paste0(results_path,"Incubation_plots/")
#Create folder if it does not exist
if (!dir.exists(plots_path)) {
  dir.create(plots_path, recursive = TRUE)
}


#Is goFlux up-to-date?----
#Check goFlux version: and warn if not up-to-date 
#For reproducibility: this script was developed using the goFlux sha: "67c276d87d984d55b70d7b756ad702d63dbf2038"
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
    warning(
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




#Packages and functions-----
renv::restore()
library(goFlux)
library(dplyr)
library(tidyr)
library(ggplot2)


#Load function
#Function to load data from each RData sampling event (subsite+gas_analizer).
#Maintains original data frequency and exact concentration (no reframing/resampling/interpolation/rouding)
#Requires auxfiles to have a column named "sampling" that is built as paste0(subsite,"_",gs_suffix)

load_sampling <- function(auxfile_s, RData_path){
  
  #Set RData folder as working directory
  setwd(RData_path)
  #Get sampling (name of RData file)
  s <- unique(auxfile_s$sampling)
  #Load RData file
  load(file = paste0(s,".RData"))
  #Select only relevant columns
  mydata <- mydata[,c("POSIX.time", 
                      "CO2dry_ppm", "CH4dry_ppb", "H2O_ppm",
                      "CO2_prec",   "CH4_prec",   "H2O_prec"  )]
  #Return full sampling data.frame if it contains observations
  if(dim(mydata)[1]>0){
    return(mydata)  
  } else {
    stop(paste0("No measurements for sampling ",s))
  }
}


# ---- Load auxfiles ----

#CO2 auxfile:
co2_auxfile <- read.csv(file = paste0(auxfile_path,"co2_auxfile.csv"))
#CH4 auxfile:
ch4_auxfile <- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))

#Necessary modifications: 
#transform start.time (character) to POSIXct
#add "sampling" column to match RData file names
co2_auxfile <- co2_auxfile %>%
  mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"sec")),#round start.time to nearest second
         sampling=if_else(gas_analiser=="LI-COR", paste0(subsite,"_LI-7810"), #Modify for LI-COR
                          paste0(subsite,"_",gas_analiser))) %>% 
  rename(flux_decision=co2_decission)

ch4_auxfile <- ch4_auxfile %>%
  mutate(start.time = as.POSIXct(round(as.POSIXct(start.time,tz = "UTC"),"sec")),#round start.time to nearest second
         sampling=if_else(gas_analiser=="LI-COR", paste0(subsite,"_LI-7810"), #Modify for LI-COR
                          paste0(subsite,"_",gas_analiser)))%>% 
  rename(flux_decision=ch4_decission)


#subset auxfiles (Optional) ------

#Apply filter defined in section "OPTIONS" 
co2_auxfile<- co2_auxfile %>% 
  filter(grepl(pattern, UniqueID))

ch4_auxfile<- ch4_auxfile %>%
  filter(grepl(pattern, UniqueID))



#Prepare loop ------

#Obtain total number of incubations to calculate (sum of incubations in auxfiles, depending on which gas should be calculated), to use in loop progress messages
incub_tocalc <- 0

if ("co2" %in% calculate_gas) {
  incub_tocalc <- incub_tocalc + length(co2_auxfile$UniqueID)
}

if ("ch4" %in% calculate_gas) {
  incub_tocalc <- incub_tocalc + length(ch4_auxfile$UniqueID)
}

#Initialize number of incubations calculated (updates within loop)
incub_calculated<- 0

#GHG Flux loop------
for(gs in calculate_gas){
  
  #Load appropriate auxfile and set string for gastype
  if(gs=="co2"){
    gs_auxfile<- co2_auxfile
    goflux_gs_string <- "CO2dry_ppm"
  }
  else if (gs=="ch4"){
    gs_auxfile<- ch4_auxfile
    goflux_gs_string <- "CH4dry_ppb"}
  else{stop("Incorrect calculate_gas option")}

  ##Subsite loop ----
  for(subs in unique(gs_auxfile$subsite)){
    
    #Message
    message(paste0("Processing ",toupper(gs), " of subsite ",subs))
    
    #Subset auxfile for subsite:
    gs_auxfile_subs<- gs_auxfile %>% filter(subsite==subs)
    
    #Initialize subsite pdf for plots (if plots are to be produced)
    if(produce_plots){
      plot_filename <- paste(toupper(gs) ,subs, as.character(as.Date(last(gs_auxfile_subs$start.time))),sep="_")
      pdf(file = paste0(plots_path,plot_filename,".pdf"),
          width = 8,height = 8)
    }
    
    #Sampling loop ----
    #Sampling loop (some subsites contain measures performed with different analyzers)
    for(s in unique(gs_auxfile_subs$sampling)){
      
      #Subset auxfile for sampling s
      gs_auxfile_s<- gs_auxfile_subs %>% filter(sampling==s)
      
      #Load rawdata of sampling s
      mydata_s <-load_sampling(auxfile_s = gs_auxfile_s, RData_path = RData_path)
      
      # UniqueID loop ----
      for (i in gs_auxfile_s$UniqueID) {
        
        # get UniqueID i auxfile
        gs_auxfile_i <- gs_auxfile_s %>%
          filter(UniqueID == i)
        
        # get UniqueID i rawdata
        mydata_i <- mydata_s[
          mydata_s$POSIX.time >= gs_auxfile_i$start.time &
            mydata_s$POSIX.time <= gs_auxfile_i$start.time + gs_auxfile_i$duration, ]
        
        # If rawdata does not contain data within UniqueID limits, message and skip this UniqueID
        if (nrow(mydata_i) == 0) {
          message(paste0("CAUTION! No data for UniqueID ", i))
          next
        }
        
        #Prepare mydata_i for flux calculation:
        #Add auxfile parameters (Identification and chamber details for flux.term):
        mydata_i$UniqueID <- gs_auxfile_i$UniqueID
        mydata_i$Area<-gs_auxfile_i$Area
        mydata_i$Vtot<-gs_auxfile_i$Vtot
        mydata_i$Pcham<-gs_auxfile_i$Pcham
        mydata_i$Tcham<-gs_auxfile_i$Tcham
        #Add goFlux required parameters:
        #flag: Indicates that all data is part of incubation
        mydata_i$flag=1
        #Etime: elapsed time (in seconds), used as x-axis for regressions
        mydata_i$Etime <- as.numeric(mydata_i$POSIX.time) - min(as.numeric(mydata_i$POSIX.time))
        
        ##Timegrid summary ----
        #Calculate basic rawdata parameters for quality assurance
        og_timegrid_i<- mydata_i %>% 
          summarise(UniqueID=unique(UniqueID), #Identifier
                    duration_s=gs_auxfile_i$duration, #Timespan of incubation in seconds
                    mydata_N_obs= sum(!is.na(POSIX.time)), #number of observations
                    mydata_avg_Dt = mean(diff(POSIX.time), na.rm=T), #AVG time-step
                    mydata_sd_Dt = sd(diff(POSIX.time), na.rm=T), #SD time-step
                    mydata_min_Dt = min(diff(POSIX.time), na.rm=T), #MIN time-step
                    mydata_max_Dt = max(diff(POSIX.time), na.rm=T) #MAX time-step
                    )
        
        ##GoFlux call------
        #Use goFlux function to obtain the calculated flux for LM and HM.
        #Call structure to avoid console printouts (Warnings and progress bars)
        goflux_i <- local({
          invisible(capture.output(
            res <- suppressWarnings(
              suppressMessages(
                goFlux(dataframe = mydata_i,gastype = goflux_gs_string)
              )
            )
          ))
          res
        })
        
        #Make sure all HM. parameters are numeric (even when NA), otherwise flux.plot throws error
        goflux_i<- goflux_i %>% 
          mutate(across(starts_with("HM."), as.numeric))
        
        ##Plot flux----
        if(produce_plots){
          #Create plot object (sintax structure for silent execution)
          plot_i <- local({
            invisible(capture.output(
              res<- suppressWarnings(flux.plot(
                flux.results = goflux_i,
                dataframe = mydata_i,
                gastype = goflux_gs_string,
                shoulder = 10,#just to expand plot limits
                plot.display = c("MDF", "prec", "nb.obs","flux.term"),
                best.model = FALSE,
                quality.check = FALSE
              ))
            ))
            res
          })
          
          
          #Create plot title
          titleplot_i<- paste(i, gs_auxfile_i$flux_decision, gs_auxfile_i$gas_analiser, sep=", ")
          
          #print plot (sintax to avoid printing of "[[1]]" in console)
          invisible(capture.output(print(plot_i)))
          #ADD plot title
          grid::grid.text(label = titleplot_i,
                          x = 0.5, y = 0.975,
                          gp = grid::gpar(fontsize = 14, fontface = "bold"))
        }
        
        ##calculate total.flux-----
        #total.flux is calculated based on concentration difference using the averages during the first and last 10s of incubation, divided by the elapsed time between the center of those windows, and multiplied by goFlux flux term.
        
        #define time variables and flux.term
        t.win<- 10
        etime_max<- max(mydata_i$Etime, na.rm=T)
        elapsed_time<- etime_max-t.win
        flux.term_i<- goflux_i$flux.term
        
        #Extract values from initial and final timewindow
        vals_init  <- mydata_i[[goflux_gs_string]][mydata_i$Etime < t.win]
        vals_final <- mydata_i[[goflux_gs_string]][mydata_i$Etime > (etime_max - t.win)]
        
        #Avg Initial concentration (first t.win seconds)
        C0.10s    <- mean(vals_init, na.rm = TRUE)
        C0.10s.sd <- sd(vals_init, na.rm = TRUE)
        C0.10s.n  <- sum(!is.na(vals_init))
        
        #Avg Final concentration (last t.win seconds)
        Cf.10s    <- mean(vals_final, na.rm = TRUE)
        Cf.10s.sd <- sd(vals_final, na.rm = TRUE)
        Cf.10s.n  <- sum(!is.na(vals_final))
        
        #Delta concentration and SE
        delta.conc <- Cf.10s - C0.10s
        delta.conc.se <- sqrt( (C0.10s.sd^2 / C0.10s.n) + (Cf.10s.sd^2 / Cf.10s.n))
        
        #Total flux and SE (propagated)
        total.flux <- (delta.conc / elapsed_time) * flux.term_i
        total.flux.se <- (flux.term_i / elapsed_time) * delta.conc.se
        
        #Build total.flux_i data.frame
        total.flux_i<- data.frame(
          UniqueID = gs_auxfile_i$UniqueID,
          C0.10s = C0.10s,
          C0.10s.sd = C0.10s.sd,
          C0.10s.n = C0.10s.n,
          Cf.10s = Cf.10s,
          Cf.10s.sd = Cf.10s.sd,
          Cf.10s.n = Cf.10s.n,
          total.flux = total.flux,
          total.flux.se = total.flux.se
        )
          
        # Update incub progress count
        incub_calculated <- incub_calculated + 1
        # Display progress
        message(
          paste0("Overall progress: ", round((incub_calculated / incub_tocalc) * 100, 2), "%")
        )
          
        ##Join Outputs----
        #Combine fluxes for UniqueID i
        fluxes_i<- gs_auxfile_i %>% 
          select(UniqueID) %>%
          left_join(og_timegrid_i, by="UniqueID") %>% 
          left_join(total.flux_i, by="UniqueID") %>% 
          left_join(goflux_i, by="UniqueID")
        
        #Append results from UniqueID i to rest (gas-specific) 
        #Add fluxes_i to target combined data.frame (co2 or ch4), creating these objects in first execution of each gas.
        if (gs == "co2") {
          if (!exists("all_co2flux")) {
            all_co2flux <- fluxes_i
          } else {
            all_co2flux <- rbind(all_co2flux, fluxes_i)
          }
        } else if (gs == "ch4") {
          if (!exists("all_ch4flux")) {
            all_ch4flux <- fluxes_i
          } else {
            all_ch4flux <- rbind(all_ch4flux, fluxes_i)
          }
        }
        
      }#end of UniqueID loop

    }#end of sampling loop
    
    #Close and save subsite pdf for plots (if plots are to be produced)
    if(produce_plots){
      message(paste0("Saving plot for subsite ", subs))
      dev.off()
      }
    
  }#end of subsite loop
  
  
  ##Save GHG fluxes----
  #After all fluxes are calculated, save them as csv
  if (gs == "co2" && exists("all_co2flux")) {
    write.csv(all_co2flux,
              file.path(results_path, "all_co2flux.csv"),
              row.names = FALSE)
  }
  
  if (gs == "ch4" && exists("all_ch4flux")) {
    write.csv(all_ch4flux,
              file.path(results_path, "all_ch4flux.csv"),
              row.names = FALSE)
  }
  
  # Final message, print how many incubations for each gas have been saved. 
  if (incub_calculated == incub_tocalc) {
    
    msg_parts <- character()
    
    if (exists("all_co2flux")) {
      msg_parts <- c(
        msg_parts,
        paste0("CO2 saved (", nrow(all_co2flux), " rows)")
      )
    }
    
    if (exists("all_ch4flux")) {
      msg_parts <- c(
        msg_parts,
        paste0("CH4 saved (", nrow(all_ch4flux), " rows)")
      )
    }
    
    message(
      "Flux calculation completed: ",
      paste(msg_parts, collapse = " | ")
    )
  }
  
  
  
}#End of gas loop


