#Selection criteria for flux best.model


# ---
# Authors: Miguel Cabrera
# Project: "RESTORE4Cs"
# ---

# --- Description----
#Description: this script is used to set the selection criteria to chose the best.model for CH4 and CO2 instantaneous flux (criteria below). And to produce the 'ghg'_bestflux.csv files for CH4 and CO2 static chamber incubations from project Restore4Cs. Visual inspection of incubation timeseries (pdf plots created by Fluxes_from_RData.R script) is used to confirm that best-model selection criteria is appropriate. 

#Choosing a best.flux follows this logic (Hard-thresholds set after inspection performed in this script): 

#0. If incubation is wrong (discarded due to artefacts or wrong manipulation), no model is selected as as best.model and a note is provided in best.model.flag. 
#1. IF NO bubbles are found (Visual inspection) --> GOflux criteria 
#2. IF bubbles are found (visual inspection) and LM.r2 < 0.99 --> use total.flux as best.model
#3. IF bubbles are found (visual inspection) and LM.r2 >= 0.99--> use LM as best.model

#GOflux criteria: For HM to be chosen as best flux, all criteria below must be met.
#HM.flux must exist to be chosen (is.na(HM.flux)--> LM)
#MDF fluxes are always LM (LM.flux<=MDF.lim --> LM)
#Kappa must be below maximum (HM.k>=k.max   --> LM)
#g.fact must be below threshold (CO2.g.fact<4, CH4.g.fact<3, otherwise ---> LM)
#HM.MAE must be at least 5% better than LM.MAE (otherwise default to LM)
#HM.AICc of must be below LM.AICc (otherwise default to LM)


#______________----
#CH4 Selection-----

#Clear Global Environment
rm(list = ls())

#Directories---------
#root path: You have to make sure this is pointing to the write folder on your local machine:
root_path <- "C:/Users/Miguel/Dropbox/GHG_bestflux_calculation" 

#Path to Co2 and Ch4 auxfiles with corrected start.time and duration and with discard decisions:
auxfile_path<- paste0(root_path,"/Auxfiles/") 

#Path to calculated fluxes
results_path <- paste0(root_path,"/Results/")




#Packages and functions----------
library(tidyverse)
library(ggExtra)


#goFlux::best.flux() calculates a rounded version of the MDF column, the limit under which a flux cannot be said to be significantly different from cero (based on instrument precission, duration of incubation and chamber geometry). We reproduce this MDF.lim here using exactly the same approach used in the goFlux best.flux() function. 
nb.decimal = function(x) {
  if (length(x) == 0) 
    return(numeric())
  x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% 
    as.numeric()
  x_int = floor(x) %>% abs() %>% nchar()
  x_nchr = x_nchr - 1 - x_int
  x_nchr[x_nchr < 0] = 0
  x_nchr
}


#1. Load results------

#Load ch4 auxfile with incubation discard decissions (ok/discard)
full_auxfilech4<- read.csv(file = paste0(auxfile_path,"ch4_auxfile.csv"))

#Load results for CH4 incubations:
ch4_flux<- read.csv(paste0(results_path,"all_ch4flux.csv"))

#Load ebullition inspection table: UniqueID, discard decission (ok/discard), visual ebullition (T/F)
ch4_ebul_inspection<- read.csv(file = paste0(auxfile_path,"ebullition_inspection.csv"))

#Reproduce best.flux MDF.lim, using the nb.decimal function defined above
ch4_flux<- ch4_flux %>% 
  mutate(MDF.lim = ifelse(nb.decimal(signif(MDF, 1)) != 0, signif(MDF, 2), round(MDF, 1)))



#2. Criteria for goflux -----

#Filter results to decide goflux criteria:
#To decide the goflux hard thresholds (i.e. g.fact) we remove incubations that:
  #Are Wrong incubations (those marked due to containing strong artefact, containing only backflush after ebullition or had wrong chamber manipulation), marked as "discard" in auxfile
  #Incubations with detected ebullition (by visual inspection), marked as "visual_ebullition"==T

incub_discard<- ch4_ebul_inspection %>% filter(ch4_decission=="discard") %>% pull(UniqueID)
incub_ebu_visual<- ch4_ebul_inspection %>% filter(visual_ebullition) %>% pull(UniqueID)


ch4_4gofluxcriteria<- ch4_flux %>% 
  filter(!UniqueID%in%c(incub_discard,incub_ebu_visual))


##2.1. Inspect mdf----
#We want to make sure MDF fluxes do not arise from artefacts "flattening" a real flux, all MDF incubations will default to LM 

ch4_mdf_LM<- ch4_4gofluxcriteria %>%  
  filter(abs(LM.flux)<=MDF.lim)


#MDFs CH4 visually inspected: true "below detection", flat-ish white noise incubations:
mdf_inspected<- c("s1-ca-p2-15-v-t-14:54","s1-ca-p2-4-v-d-11:24","s1-ca-p2-4-v-t-11:17","s1-ca-r1-9-v-t-10:57","s1-ca-r1-17-b-d-12:49","s1-ca-r2-1-b-d-08:47","s1-ca-r2-6-v-t-09:50","s1-da-a1-1-b-d-08:22","s1-da-a1-2-v-d-08:46","s1-da-a1-2-v-t-08:38","s1-da-a1-5-v-d-10:01","s1-da-a1-9-b-d-12:23","s1-du-a1-11-b-t-13:12","s1-du-a1-12-v-t-13:27","s1-du-a1-13-v-t-14:07","s1-du-a1-14-v-d-14:32","s1-du-a1-15-v-d-15:01","s1-du-a1-15-v-t-14:52","s1-du-a2-10-b-t-11:21","s1-du-a2-11-v-d-11:49","s1-du-a2-12-v-t-12:06","s1-du-a2-14-b-t-12:50","s1-du-a2-6-v-d-09:18","s1-du-a2-6-v-t-09:09","s1-du-p1-14-v-d-13:13","s1-du-p1-9-v-d-11:37","s1-du-p2-2-v-t-07:39","s1-du-p2-3-o-d-07:56","s1-du-p2-5-b-t-08:39","s1-du-p2-6-v-d-09:00","s1-du-p2-7-o-d-09:16","s1-du-p2-14-v-d-11:51","s1-du-p2-11-v-d-10:43","s1-du-p2-11-v-t-10:36","s1-du-p2-13-v-d-11:31","s1-du-p2-13-v-t-11:24","s1-du-p2-14-v-t-11:43","s1-du-p2-15-v-d-12:17","s1-du-p2-15-v-t-12:09","s1-du-p2-2-v-d-07:48","s1-du-p2-4-v-d-08:23","s1-du-p2-4-v-t-08:12","s1-du-p2-6-v-t-08:52","s1-du-p2-8-v-d-09:36","s1-du-p2-9-b-t-09:44","s1-du-r1-3-v-t-08:24","s1-du-r2-1-v-d-07:29","s1-du-r2-10-v-d-11:03","s1-du-r2-12-v-t-11:34","s1-du-r2-13-v-d-12:21","s1-du-r2-14-v-t-12:33","s1-du-r2-15-v-d-12:56","s1-du-r2-4-b-d-08:20","s1-va-p1-11-v-d-11:52","s1-va-p1-11-v-t-11:45","s1-va-p1-13-v-d-12:48","s1-va-p1-13-v-t-12:41","s1-va-p1-3-b-d-09:20","s1-va-p1-4-b-d-09:30","s1-va-p1-9-v-d-11:22","s1-va-r2-2-b-d-09:33","s1-va-r2-3-v-t-09:41","s1-va-r2-5-v-d-10:22",
                  "s2-ca-a1-4-v-t-09:38","s2-ca-a1-5-v-d-10:15","s2-ca-a1-5-v-t-10:02","s2-ca-a1-6-v-d-10:37","s2-ca-a1-6-v-t-10:30","s2-ca-a1-7-v-d-10:56","s2-ca-a1-7-v-t-10:47","s2-ca-r1-14-v-d-12:21","s2-cu-r1-1-b-d-08:05","s2-cu-r1-2-b-d-08:15","s2-cu-r1-3-b-d-08:28","s2-da-a1-7-b-d-10:36","s2-da-a1-7-b-t-10:41","s2-da-a1-7-v-d-10:28","s2-da-a1-7-v-t-10:20","s2-da-a1-9-b-d-11:27","s2-da-a1-9-b-t-11:21","s2-du-a1-10-v-d-12:11","s2-da-r1-2-v-d-09:22","s2-du-a1-10-v-t-12:04","s2-du-a1-11-v-d-12:28","s2-du-a1-11-v-t-12:22","s2-du-a1-13-v-d-12:58","s2-du-a1-13-v-t-12:53","s2-du-a1-14-v-d-13:16","s2-du-a1-14-v-t-13:11","s2-du-a1-15-v-d-13:35","s2-du-a1-15-v-t-13:29","s2-du-a1-16-v-d-14:05","s2-du-a1-16-v-t-13:57","s2-du-a1-4-v-d-10:13","s2-du-a1-5-b-d-10:30","s2-du-a2-3-v-t-09:34","s2-du-a2-10-b-d-12:03","s2-du-a2-10-b-t-11:50","s2-du-a2-4-v-d-09:57","s2-du-a2-8-v-d-11:19","s2-du-a2-8-v-t-11:13","s2-du-p1-3-o-d-09:31","s2-du-p1-10-b-d-11:29","s2-du-p1-11-v-t-11:44","s2-du-p1-14-v-d-12:42","s2-du-p1-2-v-t-09:15","s2-du-p1-4-v-t-09:49","s2-du-p1-8-b-d-10:56","s2-du-p2-10-v-d-12:14","s2-du-p2-10-v-t-12:09","s2-du-p2-11-v-d-12:34","s2-du-p2-11-v-t-12:27","s2-du-p2-12-v-d-12:46","s2-du-p2-12-v-t-12:41","s2-du-p2-13-b-d-13:01","s2-du-p2-13-b-t-12:54","s2-du-p2-14-v-d-13:16","s2-du-p2-2-v-d-09:41","s2-du-p2-2-v-t-09:33","s2-du-p2-3-o-d-10:11","s2-du-p2-4-b-d-10:32","s2-du-p2-4-b-t-10:25","s2-du-p2-5-o-d-10:43","s2-du-p2-7-v-d-11:31","s2-du-p2-7-v-t-11:21","s2-du-p2-8-v-d-11:47","s2-du-p2-8-v-t-11:41","s2-du-p2-9-b-d-11:59","s2-du-p2-9-b-t-11:55","s2-du-r1-15-v-t-12:55","s2-du-r1-2-v-t-09:35","s2-du-r1-7-v-d-10:54","s2-du-r1-7-v-t-10:41","s2-du-r2-8-o-d-11:11","s2-du-r2-1-v-d-09:29","s2-du-r2-10-v-d-11:48","s2-du-r2-13-v-d-12:59","s2-du-r2-14-v-d-13:16","s2-du-r2-14-v-t-13:09","s2-du-r2-15-v-t-13:25","s2-du-r2-6-v-d-10:40","s2-ri-a2-3-b-d-08:59","s2-ri-a2-9-v-t-10:19","s2-ri-r1-9-v-t-12:33",
                  "s3-ca-a1-9-v-t-10:01","s3-ca-r1-10-v-d-09:53","s3-ca-r1-10-v-t-09:44","s3-ca-r1-15-b-d-11:00","s3-ca-r1-5-v-d-08:42","s3-ca-r1-5-v-t-08:35","s3-cu-r2-2-b-d-07:07","s3-cu-r2-3-b-d-07:18","s3-da-a1-1-b-d-07:07","s3-da-a1-3-b-d-07:51","s3-da-a1-4-v-t-08:05","s3-da-a1-6-v-t-09:10","s3-da-a1-6-b-d-09:23","s3-da-a1-11-b-d-11:29","s3-da-a1-11-v-d-11:21","s3-da-a1-11-v-t-11:12","s3-da-a1-12-b-d-12:09","s3-da-a1-12-v-d-12:00","s3-da-a1-12-v-t-11:51","s3-da-a1-2-v-t-07:17","s3-da-a1-4-b-d-08:22","s3-da-a1-5-b-d-08:55","s3-da-a1-6-v-d-09:17","s3-da-a1-7-b-d-09:40","s3-da-a1-8-b-d-09:59","s3-da-a1-9-b-d-10:19","s3-du-a1-10-v-d-10:42","s3-du-a1-10-v-t-10:33","s3-du-a1-12-v-d-11:30","s3-du-a1-13-v-d-11:47","s3-du-a1-13-v-t-11:41","s3-du-a1-14-v-d-12:05","s3-du-a1-14-v-t-11:58","s3-du-a1-7-b-t-09:28","s3-du-a1-9-v-d-10:19","s3-du-a1-9-v-t-10:09","s3-du-a1-12-v-t-11:23","s3-du-a2-11-b-d-10:55","s3-du-a2-11-b-t-10:49","s3-du-a2-12-v-t-11:09","s3-du-a2-7-v-d-09:23","s3-du-a2-7-v-t-09:14","s3-du-a2-8-v-t-09:37","s3-du-a2-6-v-t-08:53","s3-du-a2-6-v-d-09:01","s3-du-a2-9-b-d-10:09","s3-du-a2-9-b-t-10:02","s3-du-p1-2-o-d-07:47","s3-du-p1-15-v-d-11:31",  "s3-du-p2-11-v-d-11:59","s3-du-p2-11-v-t-11:51","s3-du-p2-12-v-d-12:24","s3-du-p2-12-v-t-12:16","s3-du-p2-15-v-d-13:16","s3-du-p2-15-v-t-13:09","s3-du-p2-2-v-d-08:33","s3-du-p2-2-v-t-08:22","s3-du-p2-3-v-d-08:57","s3-du-p2-4-v-d-09:18","s3-du-p2-4-v-t-09:10","s3-du-p2-8-v-d-11:01","s3-du-p2-8-v-t-10:52","s3-du-r1-2-v-d-07:36","s3-du-r1-2-v-t-07:30","s3-du-r2-5-b-d-09:44","s3-va-p1-13-v-d-11:58","s3-va-p1-5-v-d-09:24","s3-va-p1-5-v-t-09:19",
                  "s4-ca-r1-2-v-t-07:55","s4-ca-r1-3-b-d-08:12","s4-ca-r1-4-v-t-08:19","s4-cu-r1-4-b-d-07:34","s4-cu-r2-2-b-d-06:55","s4-da-a1-3-v-t-08:02","s4-da-a1-5-b-d-08:32","s4-da-a1-7-v-t-09:22","s4-da-a1-10-b-d-10:04","s4-da-a1-2-b-d-07:48","s4-da-a1-8-v-t-09:40","s4-du-p1-12-v-t-10:40","s4-du-p1-14-v-d-11:26","s4-du-p1-4-v-t-08:37","s4-du-p1-6-v-d-09:17","s4-du-p1-2-v-t-08:05","s4-du-p1-13-v-t-10:50","s4-du-r2-1-v-d-08:20","s4-du-r2-12-v-t-11:06","s4-du-r2-2-v-t-08:29","s4-du-r2-3-b-d-08:47","s4-du-r2-9-v-d-10:17","s4-ri-a2-1-b-d-06:52","s4-ri-a2-6-b-d-07:41","s4-ri-a2-9-b-d-08:27","s4-va-p1-12-v-t-10:11","s4-va-p1-7-v-t-08:56","s4-va-p1-9-v-t-09:20","s4-va-r2-1-b-d-08:01","s4-va-r2-11-v-t-10:08","s4-va-r2-14-v-d-10:46","s4-va-r2-15-v-d-11:01","s4-va-r2-15-v-t-10:54","s4-va-r2-4-v-t-08:35"
)

#any MDF flux not inspected?
ch4_mdf_LM %>% filter(!UniqueID%in%mdf_inspected) %>% pull(UniqueID)

#any mdf_inspected that shouldn't have been?
mdf_inspected[!mdf_inspected%in%ch4_mdf_LM$UniqueID]



##2.2. Inspect Kmax-----
#Inspect HM.k>=k.max () to ensure that the LM flux is not influenced by artefacts (as these incubations will default to LM estimate)
ch4_kmax<- ch4_4gofluxcriteria %>%
  filter(!UniqueID%in%mdf_inspected) %>% #Already inspected (due to MDF)
  filter(HM.k>=k.max)

#In kmax_inspected, incubations with k=kmax with good LM estimates:
kmax_inspected<-c("s1-ca-p1-3-v-d-10:31","s1-ca-p2-3-v-d-10:59","s1-ca-p1-5-v-t-11:02","s1-da-a1-4-v-t-09:15","s1-da-a1-4-v-d-09:22","s1-da-a1-7-v-d-11:35","s1-du-a1-10-v-t-12:29","s1-du-a2-4-v-t-07:53","s1-du-a2-4-v-d-08:01","s1-du-p1-1-o-d-08:11","s1-du-p1-6-v-d-10:30","s1-du-p1-9-v-t-11:30","s1-du-p1-13-v-t-12:47","s1-du-p1-10-v-d-12:04","s1-du-p1-15-v-t-13:27","s1-du-p1-15-v-d-13:33","s1-du-r1-9-v-d-10:00","s1-du-r2-1-v-t-07:23","s1-du-r2-6-v-t-08:51","s1-du-r2-10-v-t-10:57","s1-du-r2-11-v-d-11:22","s1-du-r2-12-v-d-11:40","s1-du-r2-15-v-t-12:50","s1-va-a1-8-b-d-12:46",
                  "s2-ca-a1-4-v-d-09:50","s2-ca-a2-7-v-t-09:40","s2-ca-p1-2-v-d-08:54","s2-ca-r1-8-v-d-10:45","s2-da-r1-2-v-t-09:15","s2-du-p1-2-v-d-09:22","s2-du-p1-9-v-t-11:10","s2-du-p1-6-o-d-10:23","s2-du-p1-7-o-d-10:40","s2-du-p1-12-v-t-12:00","s2-du-p1-14-v-t-12:34","s2-du-r2-2-v-t-09:37","s2-du-r2-2-v-d-09:45","s2-du-r2-7-v-d-10:55","s2-du-r2-11-v-t-11:58","s2-du-r2-13-v-t-12:52","s2-ri-a2-7-b-d-10:00","s2-ri-p1-6-v-t-10:17","s2-ri-r2-12-v-t-12:35","s2-va-p1-2-v-d-09:21","s2-va-r2-1-b-t-09:44",
                  "s3-ca-a1-7-v-d-09:35","s3-da-a1-2-v-d-07:25","s3-da-a1-4-v-d-08:12","s3-du-a2-3-v-t-07:43","s3-du-p1-6-v-d-08:51","s3-du-p1-15-v-t-11:24","s3-du-r1-1-v-t-07:15","s3-du-r1-9-v-t-09:18","s3-ri-r1-9-b-d-12:22","s3-va-p1-12-b-d-11:38",
                  "s4-ca-p1-14-v-d-11:52","s4-ca-r2-2-v-d-07:27","s4-cu-r1-2-b-d-07:13","s4-da-a1-4-v-t-08:16","s4-du-a2-15-v-d-13:04","s4-du-p1-4-v-d-08:45","s4-du-p1-11-b-d-10:30","s4-du-p2-14-v-d-12:08","s4-du-p2-15-v-d-12:25","s4-du-r2-2-v-d-08:36","s4-du-r2-10-v-t-10:26","s4-du-r2-11-v-d-10:56","s4-ri-a1-5-b-d-07:53","s4-ri-a2-7-v-d-07:58","s4-ri-p1-1-v-d-06:59","s4-ri-p1-4-v-t-07:38","s4-ri-p1-5-o-d-07:54","s4-ri-p1-11-v-t-09:30","s4-va-p1-12-v-d-10:18","s4-va-p1-14-v-d-10:47","s4-va-p2-10-v-t-08:42","s4-va-p2-14-v-t-09:27","s4-va-r2-3-v-d-08:27","s4-va-r2-11-v-d-10:15")


tocropand_reinspect<- c(
                        )

#Any kmax not inspected?
ch4_kmax %>% filter(!UniqueID%in%c(kmax_inspected,tocropand_reinspect)) %>% pull(UniqueID)

#Any kmax_inspected that shouldnt?
kmax_inspected[!kmax_inspected%in%ch4_kmax$UniqueID]
 


##2.3. Inspect g.fact ------

#Filter the dataset to remove already inspected incubations (mdf, kappamax, tocrop)
ch4_gfact<-ch4_4gofluxcriteria %>%
  filter(!UniqueID%in%c(mdf_inspected,kmax_inspected,tocropand_reinspect)) 

#lets see what is the highest g.fact of the best-performing HM models (best 50%)
ch4_gfact %>% 
  filter(HM.r2>quantile(HM.r2, 0.5, na.rm=T)) %>%
  ggplot( aes(x=HM.r2, y=g.fact))+
  geom_point()+
  geom_hline(yintercept = 3)


#Here we check that the hard-threshold chosen for g.fact (i.e. if g.fact > threshold --> LM as best flux) does not cause any truly non-linear pattern to be lost.


#Based on the g.fact distribution above, we will start the inspection with g.fact=3 as threshold and progressively lower it until clear and clean non-linear patterns appear

#These are All the incubations with unreasonably large g-fact. They have all been visually inspected and LM is more representative for them. 
gfact_inspected<- c("s1-ca-p1-5-v-d-11:08","s1-ca-p1-8-v-t-12:12","s1-ca-p1-11-v-d-13:35","s1-ca-p1-14-v-d-14:49","s1-ca-p1-15-v-t-15:02","s1-ca-r1-3-v-d-09:33","s1-ca-r2-6-v-d-09:56","s1-ca-r2-9-v-d-10:56","s1-da-r1-5-o-d-09:05","s1-du-p1-8-v-d-11:13","s1-du-r1-9-v-t-09:54","s1-ri-p2-1-v-d-09:37","s1-ri-r1-9-b-t-11:39","s1-ri-r1-10-v-t-11:56","s1-ri-r1-13-v-t-12:48","s1-va-a1-7-b-d-12:37","s1-va-p2-8-b-t-11:12","s2-ca-p1-6-v-t-10:06","s2-ca-p2-8-b-t-10:50","s2-ca-p2-9-v-d-11:09","s2-ri-a1-1-b-t-08:50","s2-ri-a1-13-b-t-10:57","s2-ri-p2-3-v-t-09:09","s2-ri-p2-3-v-d-09:15","s2-ri-p2-8-v-d-10:58","s3-ca-a1-10-v-d-10:23","s3-cu-p1-2-v-d-06:50","s3-du-p1-9-v-d-09:38","s3-ri-r1-8-v-d-11:36","s3-ri-r1-14-v-t-14:12","s4-ca-a1-12-v-d-12:55","s4-ca-p1-8-v-d-09:50","s4-ca-p2-10-v-t-12:50","s4-ca-r2-7-v-d-08:46","s4-du-a1-7-v-d-10:09","s4-du-a1-11-v-t-11:40","s4-du-a1-12-v-d-12:04","s4-du-p2-2-v-d-08:16","s4-ri-a1-10-b-d-08:40","s4-ri-r2-4-v-t-07:32","s4-va-p1-6-b-d-08:47","s4-va-p2-3-v-t-07:19","s4-va-p2-13-v-t-09:11","s4-du-a1-15-v-t-12:53","s4-va-a1-3-b-d-09:10","s1-ca-r2-15-v-d-13:41","s1-du-p1-5-v-t-09:56","s1-du-p1-14-v-t-13:07","s1-ri-p2-11-v-d-12:44","s1-ri-r1-2-v-d-09:58","s1-ri-r1-3-v-t-10:06","s2-ca-p1-1-v-t-08:27","s2-ri-a1-2-b-d-09:15","s2-ri-r1-7-v-t-11:59","s3-du-r2-11-v-d-11:30","s4-ca-p1-15-v-t-12:13","s4-ca-p2-11-v-t-13:08","s4-du-p1-6-v-t-09:11","s4-ri-a1-2-b-d-07:22","s4-ri-r1-5-v-t-08:11","s4-ri-r1-6-v-d-08:35","s4-va-p1-10-v-t-09:34","s1-ca-r2-15-v-t-13:34","s1-ca-p1-6-v-d-11:24","s1-cu-a1-13-v-d-10:47","s1-du-p1-7-b-d-10:45","s1-du-p1-10-v-t-11:57","s1-ri-a1-8-b-t-12:05","s1-ri-r1-15-v-d-13:35","s1-va-a1-1-o-d-09:49","s2-cu-p1-11-v-d-09:57","s2-du-a1-1-o-d-09:21","s2-va-p1-6-o-d-10:34","s2-va-r2-1-b-d-09:52","s3-da-p2-9-v-t-10:51","s3-du-r1-11-v-t-09:53","s3-du-r2-7-o-d-10:12","s3-va-a1-2-v-d-09:13","s4-du-a2-7-v-d-10:19","s4-du-a2-11-v-t-11:34","s4-du-p1-13-v-d-11:03","s4-du-r2-6-o-d-09:19","s1-ca-p1-1-o-d-09:25","s3-da-p1-14-v-d-11:58","s4-ca-a1-12-v-t-12:48","s4-du-a2-10-v-d-11:19","s4-du-a2-12-v-d-12:00","s4-du-a2-14-v-d-12:47","s4-ri-a1-3-b-d-07:32","s4-ri-a2-11-b-d-08:54")#gfact>=3


#To fix and re-inspect: too long incubations, aproaching asyntotic concentration
toreinspect_gfact<- c(
                      )


#Any g.fact to inspect?
ch4_gfact %>% 
  filter(!UniqueID%in%c(gfact_inspected,toreinspect_gfact)) %>% #already inspected
  filter(g.fact>=3) %>% 
  select(UniqueID, g.fact)

#any gfact_inspected that shouldnt?
gfact_inspected[!gfact_inspected%in%ch4_gfact$UniqueID]


#g.fact threshold of 3 is the optimal solution, many true non-linear with g.fact >2.5


##2.4. AICc & MAE improvement-----

#To decide for the remaining incubations (those non-ebullitive, non-mdf, non-kmax, non-gfact>3), we will use a combination of AICc weight (probability 0-1 of HM model being best) and relative improvement on MAE with HM.
#Compare distributions for each criteria, inspect critical decissions. 


#After checking that g.fact and mdf criteria are appropriate, subset dataset that needs a decision between LM and HM:
ch4_todecide<- ch4_4gofluxcriteria %>% #Non-wrong, non-ebullitive (auto or visual)
  filter(!is.na(HM.flux)) %>% #Non-NA HM (nothing to decide in those cases, default LM)
  filter(!UniqueID%in%c(mdf_inspected,kmax_inspected,tocropand_reinspect, gfact_inspected, toreinspect_gfact))#Non-mdf, non-kmax, non-gfact>3


#Inspect AICc for selection criteria (using AICc weight)
#Test criteria: using AICc weight (based on the relative difference in AICc, calculate the AICc weight metric, a probabilistic estimate of how likely it is that one model is better than the other). 
#ADD relative improvement in MAE as additional threshold criteria (for HM to be selected, it has to reduce model error by more than 5%).

ch4_todecide <- ch4_todecide %>%
  mutate(
    delta_AICc_LM = LM.AICc - pmin(LM.AICc, HM.AICc),  # Compare to the best AICc (min AICc)
    delta_AICc_HM = HM.AICc - pmin(LM.AICc, HM.AICc),
    # Calculate the Akaike weight for HM: AICc_weight_HM, how likely is that HM is better (0-1)?
    AICc_weight_HM = exp(-0.5 * delta_AICc_HM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM)),
    # Calculate relative improvement in MAE using HM
    rel_reduct_MAE_HM  = (LM.MAE  - HM.MAE)  / LM.MAE,
    k.ratio=HM.k/k.max
  )

#Inspect AICc_weight_HM
ch4_todecide %>%
  ggplot(aes(x=AICc_weight_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Most cases have very clear separation. 

#Inspect against LM r2 (rel_reduction in MAE as a guide)
ch4_todecide %>%
  ggplot(aes(x=AICc_weight_HM, y=LM.r2,col=rel_reduct_MAE_HM>0.05))+
  geom_point()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Many incubations that win with AICc weight, have poor LM r2 (indicating noise) and do not improve MAE by more than 5%


#Inspect rel_reduct_MAE_HM
ch4_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="Relative reduction in MAE with HM")

ch4_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, y=LM.r2,col=AICc_weight_HM>0.5))+
  geom_point()+
  geom_vline(xintercept=0.025)+
  scale_x_continuous(name="Relative reduction in MAE with HM")


#Inspect against k.ratio 
ggMarginal(
  ch4_todecide %>%
    ggplot(aes(x=AICc_weight_HM,y=k.ratio, col=rel_reduct_MAE_HM>0.05))+
    geom_point()+
    # geom_hline(yintercept=1.25)+
    labs(col=NULL)+
    scale_x_continuous(name="AICc weight for HM, probability of HM being best")
)



#Inspect cases with better AICc for HM, but low reduction in MAE (<5%) and with relatively big differences in flux (gfact>1.5)

#Inspect to decide 
ch4_lowMAEimprovement<-ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(g.fact>1.5)

#This cases have been inspected: typical behavior is linear with a few artefacts forcing a curve in the model. 
bestLM<- c("s1-ca-p1-4-v-t-10:43","s1-ca-p1-9-v-t-12:32","s1-ca-p1-15-v-d-15:09","s1-ca-r2-8-v-t-10:25","s1-da-p2-3-v-d-11:07","s1-da-r2-14-b-t-12:26","s1-du-p1-11-v-t-12:17","s1-va-p2-14-b-d-13:13","s2-ca-p1-4-v-t-09:26","s2-ca-p1-10-o-d-12:22","s2-ca-p2-2-b-t-09:18","s2-ca-p2-12-v-t-12:36","s2-ca-r1-8-v-t-10:39","s2-ca-r2-10-v-t-11:07","s2-du-r1-3-v-t-09:52","s2-du-r1-15-v-d-13:01","s2-du-r2-7-v-t-10:49","s2-du-r2-9-o-d-11:25","s2-ri-a1-10-b-t-10:30","s2-ri-p2-4-v-t-09:27","s2-ri-p2-11-v-t-12:17","s2-ri-r1-2-v-t-10:21","s2-ri-r1-3-v-d-10:51","s2-ri-r1-4-v-d-11:12","s2-ri-r1-5-v-t-11:20","s2-ri-r1-5-v-d-11:30","s2-ri-r1-6-v-d-11:47","s2-va-r2-4-o-d-10:34","s3-ca-p1-3-v-t-08:00","s3-da-a1-5-v-t-08:40","s3-du-r1-12-v-d-10:16","s3-du-r2-2-v-t-09:04","s3-du-r2-9-v-t-10:46","s3-ri-p2-2-v-t-09:18","s3-ri-r1-6-v-t-10:41","s4-ca-a1-13-v-d-13:31","s4-ca-p1-15-v-d-12:21","s4-ca-p2-9-v-t-12:29","s4-du-p1-1-b-d-07:55","s4-du-p1-8-b-d-09:43","s4-du-p2-13-v-t-11:40","s4-du-r1-3-v-t-08:52","s4-ri-a2-4-v-t-07:19","s4-ri-a2-4-v-d-07:25","s4-ri-a2-10-v-t-08:40","s4-ri-p1-7-o-d-08:26","s4-ri-p1-12-v-t-09:44","s4-ri-p1-12-v-d-09:49","s4-ri-p2-1-v-d-07:09","s4-ri-p2-8-v-t-09:04","s4-ri-p2-10-v-t-09:38","s4-ri-r1-1-v-d-07:11","s4-ri-r1-9-v-d-10:13","s4-ri-r2-12-v-t-10:12","s4-va-a2-11-b-d-10:18","s4-va-p2-5-v-t-07:42","s4-va-p2-14-v-d-09:34","s4-va-r2-3-v-t-08:19","s4-va-r2-4-v-d-08:42","s1-ri-p2-4-v-t-10:24","s1-ri-p2-12-o-d-12:51","s1-ri-r1-8-v-t-11:24","s2-du-p1-4-v-d-09:57","s2-ri-r1-8-v-t-12:16","s2-ri-r1-10-v-t-12:51","s3-ca-p1-13-v-d-12:07","s3-da-a1-5-v-d-08:47","s3-du-p1-8-b-d-09:20","s3-du-r2-9-v-d-10:53","s3-ri-p2-3-v-d-09:48","s3-va-r2-3-b-d-09:10","s4-ca-p1-10-v-t-10:44","s4-cu-r1-1-b-d-07:06","s4-du-p2-4-o-d-08:47","s4-ri-a1-4-b-d-07:42","s4-ri-p1-4-v-d-07:43","s4-ri-p1-6-v-d-08:15","s4-ri-p2-4-o-d-07:55","s4-ri-p2-11-v-t-10:06","s4-ri-r1-5-v-d-08:19","s4-ri-r1-8-v-t-09:22","s4-ri-r2-10-v-t-09:40","s4-ca-a1-8-v-t-10:52")

noclear<- c("s1-du-a1-6-b-t-10:30","s1-ri-r1-14-v-d-13:20","s2-ca-p1-9-o-d-12:12","s2-cu-p1-11-v-t-09:50","s2-va-p1-8-o-d-11:07","s3-ri-p2-7-v-t-10:58","s3-ri-r1-14-v-d-14:19","s3-va-r2-4-b-d-09:18","s3-va-r2-5-v-t-09:27","s4-ri-p1-10-v-t-09:10","s4-ri-p2-2-o-d-07:17","s4-ri-r2-11-v-t-09:53","s4-ri-r2-11-v-d-10:01","s1-ca-p1-13-o-d-14:13","s1-cu-r1-13-v-d-09:02","s1-ri-r1-14-v-t-13:12","s2-ca-r2-7-v-d-10:19","s2-va-r2-9-v-d-11:47","s3-ca-p1-2-v-d-07:47","s4-du-a1-1-o-d-08:13")

 
bestHM<-c()

discard_artefacts<-c()

crop_reinspect<- c()


#Any to inspect?
ch4_lowMAEimprovement %>% 
  filter(!UniqueID%in%c(bestLM,bestHM,noclear,discard_artefacts,crop_reinspect)) %>%
  select(UniqueID, g.fact,k.ratio, rel_reduct_MAE_HM) %>% head()

#Any inspected that shouldnt?
c(bestLM,noclear)[!c(bestLM,noclear)%in%ch4_lowMAEimprovement$UniqueID]


#These criteria work well for all, a few-low g.fact with unclear HM or LM (not critical), but overall very good. 

#quality check: inspect high k.ratio (strong curvature) of selected HM and inspect visually.
k_final_good<-c("s1-da-a1-3-b-d-08:56","s3-va-p1-11-b-d-11:28","s2-ca-a1-8-b-d-11:08")

ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(k.ratio>0.5) %>% 
  filter(!UniqueID%in%c(crop_reinspect, k_final_good)) %>% 
  select(UniqueID,g.fact,k.ratio)


#quality check: inspect MAE distribution of selected HM and LM. Visually review worse incubations. 
hm_apropriate<- c("s2-ca-p2-15-b-t-13:23","s4-cu-a2-3-o-d-08:08","s2-ri-r2-2-o-d-12:47","s4-cu-a2-2-v-d-07:54","s3-va-r2-7-b-d-10:02","s4-va-a2-1-b-d-07:13","s4-va-a2-2-o-d-07:21","s2-cu-r1-9-o-d-10:03","s4-cu-a2-2-v-t-07:45","s2-da-r1-1-o-d-09:00")
non_detected_ebullition<-c()


lm_apropriate<-c("s1-ri-r2-3-o-d-14:05","s3-ca-a1-2-v-d-07:56","s1-ca-a2-12-v-d-13:07")
#HM chosen: check high MAE, worse-scenario is ebullitive dynamics and HM chosen.
ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(HM.MAE>15) %>% 
  arrange(desc(HM.MAE)) %>% 
  filter(g.fact>1.1) %>% 
  filter(LM.flux>5) %>%
  select(UniqueID, HM.MAE, g.fact, k.ratio)%>% 
  filter(!UniqueID%in%c(hm_apropriate,lm_apropriate,non_detected_ebullition)) %>% head()




ch4_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(!UniqueID%in%c(hm_apropriate,lm_apropriate,non_detected_ebullition)) %>% 
  ggplot(aes(x=HM.MAE))+
  geom_histogram()


#quality check: inspect MAE of selected LM (not critical, it is what it is)
ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate)) %>% 
  ggplot(aes(x=LM.MAE))+
  geom_histogram()


#Quality check LM. inspect ratio between total flux and LM to detect ebullitive patterns that are concave (exponential curve) for which HM is never selected, but LM is also not appropriate. Re-classify to ebullition after visual inspection. 
non_detected_ebullition<-c(
                           )

lm_apropriate<- c("s1-cu-p1-3-o-d-08:33","s1-cu-r1-15-v-d-09:44","s1-ri-a1-16-o-d-14:00","s1-va-r1-7-o-d-10:38","s1-ca-r1-10-v-t-11:15","s1-cu-a1-15-v-d-11:47","s1-cu-r1-13-v-d-09:02","s1-da-p2-1-v-d-10:30","s1-da-p2-3-v-d-11:07","s1-du-r1-10-v-d-10:20","s1-du-r2-7-v-t-09:16","s1-va-a1-10-v-t-13:03","s1-va-r1-13-v-d-12:55","s2-da-p2-3-v-d-10:19","s2-da-p2-4-o-d-10:26","s2-da-r1-9-v-d-10:56","s2-va-a2-9-v-d-11:33","s3-cu-r1-11-o-d-09:28","s3-va-a1-8-v-d-11:44","s4-ca-a1-11-v-t-12:14","s4-cu-r1-9-v-d-09:10","s4-da-r1-3-v-t-09:25","s4-du-a1-4-v-d-09:05","s4-du-a2-4-v-t-08:34","s4-du-a2-4-v-d-08:43","s4-va-a1-12-v-d-12:25","s4-va-a2-7-v-d-09:05","s4-va-p2-3-v-d-07:25","s4-va-p2-14-v-d-09:34","s4-va-r1-12-v-t-11:28","s4-va-r1-12-v-d-11:39","s4-va-r2-9-v-t-09:44","s1-va-r1-15-v-d-13:44","s2-da-p2-12-v-d-11:57","s3-du-r1-8-v-t-09:00","s4-ca-p1-9-v-d-10:33","s2-cu-p1-6-v-t-08:39","s3-da-p2-9-v-d-10:58")

review_cropORdiscard<-c()


ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(LM.flux>2) %>% 
  mutate(lmtotalratio=LM.flux/total.flux) %>% 
  filter(!between(lmtotalratio,0.925,1.075)) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate,review_cropORdiscard)) %>% 
  select(UniqueID,LM.MAE, g.fact,LM.r2,lmtotalratio) %>% head()

ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(LM.flux>2) %>% 
  filter(between(LM.flux/total.flux,0.5,1.5)) %>% 
  filter(!UniqueID%in%c(non_detected_ebullition,lm_apropriate)) %>%
  ggplot(aes(x=LM.MAE,y=LM.flux/total.flux))+
  geom_point()+
  geom_hline(yintercept = 1.1)+
  geom_hline(yintercept = 0.9)




#Extra: check distribution of LM vs total.flux for very good LM, inspect biass in mean.total.flux
ch4_todecide %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(LM.flux>2) %>% 
  ggplot(aes(x=LM.flux, y=total.flux, col=between(LM.flux/total.flux, 0.90,1.1)))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

#Similar SE for both approaches. 
ch4_todecide %>% 
  filter(LM.r2>0.85) %>% 
  ggplot(aes(x=abs(total.flux.se), y=abs(LM.SE)))+
  geom_point()
  
ch4_todecide %>% 
  filter(LM.r2>0.99) %>% 
  ggplot(aes(x=LM.flux/total.flux))+
  geom_histogram()




#3. Criteria for ebullition-----


##3.1. LM.r2 threshold-----
#If ebullition is detected by visual inspection: 

  #IF LM.r2 > threshold ---> select LM
  #IF LM.r2 < threshold ---> select totalflux 

#Check threshold to chose LM over totalflux, when ebullition is present, total.flux should be more reliable than linear depending on the shape. check difference Linear/total to decide LM.r2 threshold 

ch4_ebullitive<- ch4_flux %>% 
  filter(!UniqueID%in%incub_discard) %>% #Not wrong incubations
  filter(UniqueID%in%c(incub_ebu_visual)) %>% #Contains ebullition (visual)
  mutate(linear_total_ratio=LM.flux/total.flux)

ggMarginal(
  ch4_ebullitive %>% 
    filter(LM.r2>0.95) %>% 
    ggplot(aes(x=LM.r2, y=linear_total_ratio))+
    geom_point()+
    geom_vline(xintercept = 0.99)
)

#Threshold to chose LM over total.flux will be set at LM.r2 0.99, as below this threshold there are a lot of incubations whose linear flux deviates more than 10% from the total.flux


#4. CH4 bestflux selection-------
#We apply the criteria above to select the best_model in every case (ebullition or not) and add best_model_flags with flags detailing the reason for best model.

##4.1. Selection & flags----
ch4_bestflux<- ch4_flux %>% 
  #Create logical vectors for all cases, to be able to inspect (if any) cases where the best_model is overriden, due to errors in the order of the criteria. 
  mutate(is_discard=UniqueID%in%incub_discard,
         is_ebullition=UniqueID%in%c(incub_ebu_visual),
         is_LMr2_99=LM.r2>=0.99,
         is_NA_HM=is.na(HM.flux),
         is_mdflinear=abs(LM.flux)<=MDF.lim,
         is_kmax=HM.k>=k.max,
         is_gfactexceed=g.fact>=3,
         is_lowmaeimprove=((LM.MAE  - HM.MAE)  / LM.MAE)<0.05,
         is_LMaiccbetter=LM.AICc<HM.AICc) %>% 
  #The conditions are evaluated sequentially (if multiple conditions are met, only the first is used), if no condition is met, default to "HM"
  mutate(best_model=case_when(is_discard~"None appropriate",
                              is_ebullition&is_LMr2_99~"LM",
                              is_ebullition&!is_LMr2_99~"total.flux",
                              is_NA_HM~"LM",
                              is_mdflinear~"LM",
                              is_kmax~"LM",
                              is_gfactexceed~"LM",
                              is_lowmaeimprove~"LM",
                              is_LMaiccbetter~"LM",
                              TRUE~"HM"),
         # Add quality flags (some override others)
         best_model_flags = pmap_chr(
           list(is_discard, is_ebullition, is_LMr2_99, is_NA_HM,
                is_mdflinear, is_kmax, is_gfactexceed, is_lowmaeimprove, is_LMaiccbetter),
           function(discard, ebu, lmr2, na_hm, mdf, kmax, gfact, lowmae, aicc) {
             # Priority override
             if (discard) return("Discard incubation, artefacts preclude flux estimate")
             
             # EBULLITION branch
             if (ebu) {
               ebullition_flags <- c("Ebullitive dynamics")
               if (lmr2) ebullition_flags <- c(ebullition_flags, "LM.r2 >= 0.99")
               if (!lmr2) ebullition_flags <- c(ebullition_flags, "LM.r2 < 0.99")
               return(str_c(ebullition_flags, collapse = ", "))
             }
             
             # NON-EBULLITION branch
             other_flags <- c()
             if (na_hm) other_flags <- c(other_flags, "HM.flux is NA")
             if (mdf) other_flags <- c(other_flags, "MDF")
             if (isTRUE(kmax)) other_flags <- c(other_flags, "HM.k exceeds maximum")
             if (isTRUE(gfact)) other_flags <- c(other_flags, "g.fact >= 3")
             if (isTRUE(lowmae)) other_flags <- c(other_flags, "HM does not improve MAE by > 5%")
             if (isTRUE(aicc)) other_flags <- c(other_flags, "AICc of LM is best")
             
             if (length(other_flags) == 0) return("HM meets all criteria")
             str_c(other_flags, collapse = " | ")
           }
         )
  ) %>% 
  mutate(best.flux=case_when(best_model=="LM"~LM.flux,
                             best_model=="HM"~HM.flux,
                             best_model=="total.flux"~total.flux,
                             best_model=="None appropriate"~NA_real_))



ch4_bestflux %>%
  separate_rows(best_model_flags, sep = " \\| ") %>%
  count(best_model, best_model_flags)

ch4_bestflux %>% 
  count(best_model)

  


#QUality check: compare which method is assigned to biggest fluxes (HM should not be present among the highest fluxes, as these most likely arise from ebullitive patterns)

#Check CH4 flux distribution according to models
ch4_bestflux %>% 
  filter(!is.na(best.flux)) %>% 
  ggplot(aes( x=abs(best.flux), fill=best_model))+
  geom_histogram(bins = 100)+
  scale_x_log10()+
  geom_vline(xintercept=25)+
  facet_wrap(~best_model, nrow=3)

#Check HM models highest relative uncertainty (make sure we are not assigning HM when there is ebullition). 
checked_HM_correct<- c("s1-cu-r1-15-v-d-09:44","s2-ca-p2-15-b-t-13:23","s1-ca-a2-12-v-t-13:01","s4-va-a2-2-o-d-07:21","s1-va-r2-12-v-t-12:12","s2-cu-p1-6-v-d-08:45","s2-da-r2-5-o-d-08:25","s2-va-a2-13-v-t-12:38","s4-va-a2-1-b-d-07:13","s3-cu-p1-2-v-t-06:41","s4-ca-r2-11-v-t-09:36","s4-du-r2-11-v-t-10:49","s1-cu-p1-18-v-d-14:35","s1-da-p2-3-v-t-10:58","s1-da-p2-1-v-t-10:22","s3-du-p1-13-v-t-10:46","s2-va-a2-5-v-t-10:11","s2-va-a2-4-v-t-09:48"
                       )

ebullition<- c(
               )

ch4_bestflux %>%
  filter(best_model=="HM") %>% 
  filter(abs(best.flux)>5) %>% 
  filter(HM.se.rel>5) %>% 
  filter(g.fact>1.1) %>% 
  filter(!UniqueID%in%c(checked_HM_correct,ebullition)) %>%
  arrange(desc(HM.se.rel)) %>% 
  select(UniqueID, best.flux, g.fact, best_model, HM.MAE,HM.se.rel)
  


##4.2. Final format ----
names(ch4_bestflux)
#provisionally, only best.flux to be able to work with the data. 
ch4_bestflux_formated<- ch4_bestflux %>% 
  #general
  select(UniqueID, duration_s, flux.term,MDF.lim,
         #totalflux results
         C0.10s,C0.10s.sd,Cf.10s,Cf.10s.sd,total.flux, total.flux.se, 
         #LM results
         LM.flux, LM.SE, LM.MAE, LM.RMSE, LM.AICc, LM.r2, LM.p.val,
         #HM.results
         HM.flux, HM.SE, HM.MAE, HM.RMSE, HM.AICc, HM.r2, HM.k, k.max,
         #Fluxchoice results
         g.fact, best_model, best_model_flags,
         ) %>% 
  rename(duration=duration_s,
         LM.flux.se=LM.SE,
         HM.flux.se=HM.SE,
         best.model=best_model, best.model.flags=best_model_flags)

names(ch4_bestflux_formated)

#save ch4_bestflux
write.csv(ch4_bestflux_formated, file = paste0(results_path, "ch4_bestflux.csv"),row.names = F)





#______________----
#CO2 Selection-----

#IMPLEMENT same goflux criteria that for CH4 re-defining hard thresholds for g.fact based on co2 incubation inspections.  

#Clean global environment, leaving the data-paths and "nb.decimal" function only: 
rm(list = setdiff(ls(), c("dropbox_root", "auxfile_path", "results_path","nb.decimal")))


#1. Load results------

#Load discard decissions from auxfile co2
full_auxfileco2<- read.csv(file = paste0(auxfile_path,"co2_auxfile.csv"))

#Load results for CO2 incubations:
co2_flux<- read.csv(paste0(results_path,"all_co2flux.csv"))


#Reproduce best.flux MDF.lim
#goFlux::best.flux() calculates a rounded version of the MDF column, the limit under which a flux cannot be said to be significantly different from cero (based on instrument precission, duration of incubation and chamber geometry). We reproduce this MDF.lim here using exactly the same approach used in the goFlux best.flux() function. 

nb.decimal = function(x) {
  if (length(x) == 0) 
    return(numeric())
  x_nchr = x %>% abs() %>% as.character() %>% nchar() %>% 
    as.numeric()
  x_int = floor(x) %>% abs() %>% nchar()
  x_nchr = x_nchr - 1 - x_int
  x_nchr[x_nchr < 0] = 0
  x_nchr
}

co2_flux<- co2_flux %>% 
  mutate(MDF.lim = ifelse(nb.decimal(signif(MDF, 1)) != 0, signif(MDF, 2), round(MDF, 1)))



#2. Criteria for goflux -----

#Filter results to decide goflux criteria:
#To decide the goflux hard thresholds (i.e. g.fact) we remove incubations marked as "discard" due to strong artefact, ebullition-caused artefact or wrong manipulation)

incub_discard<- full_auxfileco2 %>% filter(co2_decission=="discard") %>% pull(UniqueID)

co2_4gofluxcriteria<- co2_flux %>% 
  filter(!UniqueID%in%c(incub_discard))



##2.1. Inspect mdf----
#We want to make sure MDF fluxes do not arise from artefacts, all MDF incubations will default to LM 
co2_mdf_LM<- co2_4gofluxcriteria %>%  
  filter(abs(LM.flux)<=MDF.lim)

#MDFs CO2 inspected: true below detection flat-ish white noise incubations
mdf_inspected<- c("s1-ca-a2-1-v-d-10:44","s1-du-p2-12-b-d-11:01","s1-ri-p1-1-v-d-10:19","s1-ri-p1-3-v-d-10:51","s1-ri-r2-4-b-d-10:01","s2-du-a1-15-v-t-13:29","s2-du-a2-1-o-d-08:54","s3-du-p2-7-b-d-10:26","s3-va-p1-13-v-t-11:48","s1-ca-a2-13-o-d-13:17","s1-ca-a2-4-v-d-11:25","s1-ca-a2-5-o-d-11:33","s1-ca-a2-6-v-d-11:49","s1-ca-a2-8-o-d-12:12","s1-ca-r1-4-o-d-09:45","s1-ca-r1-6-o-d-10:17","s1-cu-a1-14-v-d-11:14","s1-cu-a2-14-o-d-11:32","s1-cu-a2-16-o-d-11:58","s1-cu-a2-4-o-d-07:53","s1-cu-a2-6-o-d-08:12","s1-cu-p2-9-o-d-10:56","s1-cu-r2-18-v-d-10:36","s1-cu-r2-21-b-d-11:20","s1-da-p1-3-o-d-09:11","s1-da-p1-5-o-d-09:33","s1-da-p1-8-o-d-10:05","s1-da-p1-9-o-d-10:17","s1-du-a1-2-o-d-09:26","s1-du-a1-3-o-d-09:38","s1-du-a1-4-o-d-09:50","s1-du-p2-5-b-t-08:39","s1-du-p2-6-v-t-08:52","s1-du-r2-8-v-t-09:44","s1-ri-a2-10-b-d-11:50","s1-ri-a2-12-b-d-12:05","s1-ri-a2-5-b-d-10:56","s1-ri-a2-6-b-d-11:01","s1-ri-p1-10-v-t-12:59","s1-ri-p1-4-v-d-11:11","s1-ri-r2-6-v-d-10:39","s1-ri-r2-7-b-d-10:50","s1-va-a2-1-o-d-09:02","s1-va-a2-4-v-d-10:01","s1-va-a2-6-o-d-10:35","s1-va-p1-13-v-t-12:41","s1-va-r2-13-o-d-12:34","s1-va-r2-14-o-d-12:53","s1-va-r2-15-o-d-13:13","s1-va-r2-8-b-d-11:08",
                  "s2-ca-a1-3-v-t-09:12","s2-ca-a1-8-b-d-11:08","s2-ca-a2-1-v-d-08:15","s2-ca-a2-10-v-d-10:58","s2-ca-a2-7-v-d-09:47","s2-ca-a2-5-v-t-09:15","s2-ca-a2-2-v-d-08:33","s2-ca-a2-9-v-d-10:31","s2-ca-r1-14-v-d-12:21","s2-cu-a1-4-v-t-08:59","s2-cu-a2-2-v-t-08:25","s2-cu-a2-11-o-d-10:30","s2-cu-a2-14-o-d-11:03","s2-cu-a2-15-o-d-11:11","s2-cu-a2-3-v-d-08:54","s2-cu-a2-6-o-d-09:42","s2-cu-a2-7-o-d-09:50","s2-cu-a2-8-o-d-10:04","s2-cu-a2-9-o-d-10:13","s2-cu-p1-2-v-t-07:51","s2-cu-p2-1-v-d-08:26","s2-cu-p2-1-v-t-08:18","s2-cu-p2-10-o-d-10:52","s2-cu-p2-11-o-d-11:08","s2-cu-p2-13-o-d-11:23","s2-cu-p2-14-o-d-11:31","s2-cu-p2-15-o-d-11:39","s2-cu-p2-2-v-d-08:46","s2-cu-p2-2-v-t-08:38","s2-cu-p2-3-v-d-09:08","s2-cu-p2-3-v-t-09:00","s2-cu-p2-9-o-d-10:45","s2-cu-r1-4-v-t-08:40","s2-cu-r2-1-b-d-08:08","s2-cu-r2-1-b-t-08:01","s2-cu-r2-10-o-d-11:10","s2-cu-r2-11-o-d-11:23","s2-cu-r2-12-o-d-11:37","s2-cu-r2-3-b-d-08:41","s2-cu-r2-4-v-d-09:02","s2-cu-r2-4-v-t-08:55","s2-cu-r2-5-v-d-09:23","s2-cu-r2-5-v-t-09:14","s2-cu-r2-6-v-t-09:40","s2-cu-r2-9-o-d-11:00","s2-da-a1-1-o-d-08:17","s2-da-a1-3-o-d-08:43","s2-da-a2-1-o-d-07:32","s2-da-a2-10-o-d-09:31","s2-da-a2-11-o-d-09:57","s2-da-a2-12-o-d-10:06","s2-da-a2-13-o-d-10:17","s2-da-a2-15-o-d-10:46","s2-da-a2-2-o-d-07:52","s2-da-a2-3-o-d-08:04","s2-da-a2-4-o-d-08:16","s2-da-a2-5-o-d-08:26","s2-da-a2-6-o-d-08:36","s2-da-a2-7-o-d-08:47","s2-da-p1-1-o-d-08:10","s2-da-p2-2-v-t-09:58","s2-da-r2-2-v-d-07:46","s2-du-a1-1-o-d-09:21","s2-du-a1-2-o-d-09:36","s2-du-a1-3-o-d-09:48","s2-du-a2-12-v-d-12:36","s2-du-a2-14-b-d-13:04","s2-du-p2-1-o-d-09:17","s2-du-p2-3-o-d-10:11","s2-du-r1-9-v-d-11:22","s2-ri-a2-1-b-d-08:35","s2-ri-a2-16-o-d-12:08","s2-ri-a2-2-b-d-08:51","s2-ri-a2-7-b-d-10:00","s2-ri-a2-8-b-d-10:10","s2-ri-p1-13-o-d-12:52","s2-ri-p1-14-o-d-13:19","s2-ri-p1-15-o-d-13:41","s2-ri-p1-2-v-d-09:11","s2-ri-p1-3-v-d-09:31","s2-ri-p1-4-v-d-09:54","s2-ri-p1-5-v-d-10:08","s2-ri-p1-7-o-d-10:36","s2-ri-r2-12-v-d-12:40","s2-ri-r2-2-o-d-12:47","s2-ri-r2-3-o-d-13:26","s2-va-a1-2-o-d-10:21","s2-va-a1-4-v-d-11:06","s2-va-a1-4-v-t-10:57","s2-va-a1-6-b-d-11:34","s2-va-a1-7-v-t-11:57","s2-va-p2-12-v-t-11:44","s2-va-p2-13-v-d-12:14","s2-va-p2-15-v-t-12:49","s2-va-p2-5-o-d-10:05","s2-va-p2-6-o-d-10:18","s2-va-p2-7-o-d-10:35","s2-va-p2-9-b-d-11:07","s2-va-r1-1-o-d-09:18","s2-va-r1-11-v-d-12:46","s2-va-r1-11-v-t-12:36","s2-va-r1-9-v-d-11:42","s2-va-r1-12-v-d-13:13","s2-va-r1-12-v-t-13:03","s2-va-r1-5-o-d-10:22","s2-va-r1-8-o-d-11:15",
                  "s3-ca-a2-10-o-d-09:29","s3-ca-a2-12-o-d-09:56","s3-ca-a2-14-o-d-10:25","s3-ca-a2-15-v-d-10:39","s3-ca-a2-16-v-d-10:57","s3-ca-a2-2-o-d-07:43","s3-ca-a2-3-v-d-08:05","s3-ca-a2-4-o-d-08:13","s3-ca-a2-5-v-d-08:27","s3-ca-a2-5-v-t-08:20","s3-ca-a2-6-o-d-08:36","s3-ca-a2-7-v-d-08:50","s3-da-a2-1-o-d-06:06","s3-da-a2-4-o-d-07:00","s3-du-a1-6-b-d-09:19","s3-du-a2-11-b-d-10:55","s3-du-a2-3-v-t-07:43","s3-du-p2-5-o-d-09:38","s3-du-r1-9-v-t-09:18","s3-ri-a2-11-b-d-10:58","s3-ri-a2-13-b-d-11:20","s3-ri-a2-14-b-d-11:30","s3-ri-a2-16-o-d-12:03","s3-ri-a2-17-o-d-12:23","s3-ri-a2-6-b-d-09:51","s3-ri-a2-8-b-d-10:21","s3-ri-a2-9-b-d-10:29","s3-ri-r1-13-b-d-14:04","s3-va-r2-1-o-d-08:29","s3-va-r2-3-b-t-09:04","s3-va-r2-3-b-d-09:10",
                  "s4-cu-a2-11-o-d-10:31","s4-cu-r2-12-o-d-10:10","s4-da-a2-4-o-d-07:21","s4-da-a2-12-o-d-08:51","s4-da-p1-12-o-d-09:55","s4-ri-a2-6-b-d-07:41","s4-ri-p1-5-o-d-07:54")


#any MDF flux not inspected?
co2_mdf_LM %>% filter(!UniqueID%in%mdf_inspected) %>% pull(UniqueID)

#any mdf_inspected that shouldnt?
mdf_inspected[!mdf_inspected%in%co2_mdf_LM$UniqueID]


##2.2. Inspect Kmax-----
#Inspect HM.k>=k.max () to ensure that the LM flux is not influenced by artefacts (as these incubations will default to LM estimate)

co2_kmax<- co2_4gofluxcriteria %>%
  filter(!UniqueID%in%mdf_inspected) %>% #Already inspected (due to MDF)
  filter(HM.k>=k.max)

#In kmax_inspected, incubations with k=kmax with good LM estimates:
kmax_inspected<-c("s1-cu-a2-8-o-d-08:33","s1-da-r2-15-b-t-12:35","s1-ri-p1-2-v-d-10:33","s1-ri-p1-6-v-d-11:41","s1-va-a2-4-v-t-09:54","s1-va-p1-15-v-t-13:05","s1-va-p2-4-b-t-09:56","s2-ca-a2-2-v-t-08:25","s2-cu-a1-9-o-d-10:00","s2-cu-a1-10-o-d-10:09","s2-cu-a1-15-v-d-11:02","s2-cu-a2-12-o-d-10:40","s2-cu-a2-13-o-d-10:54","s2-cu-r2-3-b-t-08:34","s2-da-p2-9-o-d-11:15","s2-du-r2-4-b-d-10:07","s2-ri-a2-3-b-d-08:59","s2-ri-a2-11-b-d-10:55","s2-ri-r2-1-o-d-12:32","s2-va-a2-14-v-d-13:06","s3-ca-a1-7-v-t-09:28","s3-cu-p2-7-o-d-10:36","s3-da-a1-1-b-d-07:07","s3-du-a1-7-b-d-09:35","s3-ri-p1-17-o-d-13:28","s4-ca-a1-8-v-t-10:52","s4-cu-a1-4-v-t-08:22","s4-cu-a1-10-o-d-09:59","s4-cu-r2-13-o-d-10:17","s4-cu-r2-14-o-d-10:25")


tocropand_reinspect<- c() 


#Any kmax not inspected?
co2_kmax %>% filter(!UniqueID%in%c(kmax_inspected,tocropand_reinspect)) %>% pull(UniqueID)

#Any kmax_inspected that shouldnt?
kmax_inspected[!kmax_inspected%in%co2_kmax$UniqueID]



##2.3. Inspect g.fact ------

#Filter the dataset to remove already inspected incubations (mdf, kappamax, tocrop)
co2_gfact<-co2_4gofluxcriteria %>%
  filter(!UniqueID%in%c(mdf_inspected,kmax_inspected,tocropand_reinspect)) 

#lets see what is the highest g.fact of the best-performing HM models (best 90%)
co2_gfact %>% 
  filter(HM.r2>quantile(HM.r2, 0.1, na.rm=T)) %>%
  ggplot( aes(x=HM.r2, y=g.fact))+
  geom_point()+
  geom_hline(yintercept = 4)

#Here we check that the hard-threshold chosen for g.fact (i.e. if g.fact > threshold --> LM as best flux) does not cause any truly non-linear pattern to be lost.

#Based on the g.fact distribution above, we will start the inspection with g.fact=4 as threshold and progressively lower it until clear and clean non-linear patterns appear.

#These are All the incubations with unreasonably large g-fact. They have all been visually inspected and LM is more representative for them. 
co2_gfact_inspected<- c("s1-va-r1-5-o-d-10:04","s3-va-r1-1-v-d-08:47","s4-ri-p2-6-o-d-08:28","s4-da-p2-5-o-d-09:34","s1-ca-p2-3-v-t-10:53","s4-da-r2-7-o-d-08:01","s4-du-a1-14-b-d-12:34","s1-va-r1-13-v-d-12:55","s1-ca-p1-14-v-d-14:49","s2-cu-p1-13-o-d-10:11","s2-cu-p1-11-v-t-09:50","s4-ri-p1-14-o-d-11:04","s1-cu-p1-17-v-t-14:01","s1-da-p2-11-o-d-12:38","s3-da-p2-3-v-d-09:45","s4-da-p2-1-v-t-08:45","s2-da-p2-1-o-d-09:44","s4-ca-p2-1-v-d-08:01","s2-du-p1-15-v-d-13:01","s3-da-r1-10-o-d-10:27","s4-ri-p1-7-o-d-08:26","s4-du-r2-10-v-t-10:26","s1-da-p2-8-o-d-12:04","s2-da-p2-12-v-d-11:57","s4-ri-a1-11-b-d-08:50")#gfact>4
                        
                        
#To fix and re-inspect: too long incubations, approaching asymptotic concentration (which causes flattening of LM and thus unreasonably large g-fact)
toreinspect_gfact<- c()
                      

#Any g.fact to inspect?
co2_gfact %>% 
  filter(!UniqueID%in%c(co2_gfact_inspected,toreinspect_gfact)) %>% #already inspected
  filter(g.fact>=4) %>% 
  select(UniqueID, g.fact) %>% arrange(desc(g.fact))

#any gfact_inspected that shouldnt?
co2_gfact_inspected[!co2_gfact_inspected%in%(co2_gfact %>% filter(g.fact>=4))$UniqueID]


#g.fact threshold of 4 is the optimal solution, many true non-linear with g.fact >3.5


##2.4. AICc & MAE improvement-----

#To decide for the remaining incubations (those non-mdf, non-kmax, non-gfact>3), we will use a combination of AICc weight (probability 0-1 of HM model being best) and relative improvement on MAE with HM.
#Compare distributions for each criteria, inspect critical decissions. 


#After checking that g.fact and mdf criteria are appropriate, subset dataset that needs a decision between LM and HM:
co2_todecide<- co2_4gofluxcriteria %>% #Non-wrong
  filter(!is.na(HM.flux)) %>% #Non-NA HM (nothing to decide in those cases, default LM)
  filter(abs(LM.flux)>MDF.lim) %>% #Non-mdf
  filter(HM.k<k.max) %>% #Non-kmax
  filter(g.fact<4) #Non-gfact>=4


#Inspect AICc for selection criteria (using AICc weight)
#Test criteria: using AICc weight (based on the relative difference in AICc, calculate the AICc weight metric, a probabilistic estimate of how likely it is that one model is better than the other). 
#ADD relative improvement in MAE as additional threshold criteria (for HM to be selected, it has to reduce model error by more than 5%).

co2_todecide <- co2_todecide %>%
  mutate(
    delta_AICc_LM = LM.AICc - pmin(LM.AICc, HM.AICc),  # Compare to the best AICc (min AICc)
    delta_AICc_HM = HM.AICc - pmin(LM.AICc, HM.AICc),
    # Calculate the Akaike weight for HM: AICc_weight_HM, how likely is that HM is better (0-1)?
    AICc_weight_HM = exp(-0.5 * delta_AICc_HM) / (exp(-0.5 * delta_AICc_LM) + exp(-0.5 * delta_AICc_HM)),
    # Calculate relative improvement in MAE using HM
    rel_reduct_MAE_HM  = (LM.MAE  - HM.MAE)  / LM.MAE,
    k.ratio=HM.k/k.max
  )

#Inspect AICc_weight_HM
co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Most cases have very clear separation. 

#Inspect against LM r2 (rel_reduction in MAE as a guide)
co2_todecide %>%
  ggplot(aes(x=AICc_weight_HM, y=LM.r2,col=rel_reduct_MAE_HM>0.05))+
  geom_point()+
  scale_x_continuous(name="AICc weight for HM, probability of HM being best")
#Many incubations that win with AICc weight, have poor LM r2 (indicating noise) and do not improve MAE by more than 5%


#Inspect rel_reduct_MAE_HM
co2_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, fill=AICc_weight_HM>0.5))+
  geom_histogram(bins=100)+
  scale_x_continuous(name="Relative reduction in MAE with HM")

co2_todecide %>%
  ggplot(aes(x=rel_reduct_MAE_HM, y=LM.r2,col=AICc_weight_HM>0.5))+
  geom_point()+
  geom_vline(xintercept=0.025)+
  scale_x_continuous(name="Relative reduction in MAE with HM")


#Inspect against k.ratio 
ggMarginal(
  co2_todecide %>%
    ggplot(aes(x=AICc_weight_HM,y=k.ratio, col=rel_reduct_MAE_HM>0.05))+
    geom_point()+
    # geom_hline(yintercept=1.25)+
    labs(col=NULL)+
    scale_x_continuous(name="AICc weight for HM, probability of HM being best")
)



#Inpsect cases with better AICc for HM, but low reduction in MAE (<5%) and with relatively big differences in flux (gfact>1.5)

#Inspect to decide 
co2_lowMAEimprovement<-co2_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM<0.05) %>% 
  filter(g.fact>1.5)

#This cases have been inspected: tipical behaviour is linear with a few artefacts forcing a curve in the model. 
bestLM<- c("s1-cu-a1-14-v-t-11:07","s1-cu-p2-8-o-d-10:42","s1-cu-r2-6-o-d-07:07","s1-da-p1-10-o-d-10:34","s1-da-p2-9-o-d-12:13","s1-da-p2-13-o-d-12:59","s1-du-r1-9-v-d-10:00","s1-du-r2-4-b-d-08:20","s1-ri-a1-8-b-d-12:16","s1-ri-p2-15-o-d-13:37","s1-va-a2-7-v-d-10:53","s1-va-p2-16-b-t-13:56","s1-va-p2-16-b-d-14:02","s2-ca-p1-4-v-t-09:26","s2-cu-a1-4-v-d-09:05","s2-cu-a1-5-o-d-09:12","s2-cu-a1-8-o-d-09:52","s2-cu-a1-13-o-d-10:34","s2-cu-p2-5-o-d-09:36","s2-da-a2-8-o-d-09:02","s2-da-p1-8-o-d-09:36","s2-da-p1-14-v-t-11:04","s2-da-r2-2-v-t-07:39","s2-da-r2-13-o-d-09:46","s2-da-r2-14-v-t-10:13","s2-du-p1-13-v-d-12:23","s2-du-p2-7-v-d-11:31","s2-du-r2-10-v-t-11:41","s2-du-r2-14-v-t-13:09","s2-ri-a1-2-b-t-09:08","s2-ri-a2-10-b-d-10:48","s2-va-a2-13-v-d-12:45","s2-va-p1-14-v-d-12:51","s2-va-r1-10-v-d-12:14","s2-va-r2-7-v-t-11:08","s3-ca-a2-11-v-d-09:44","s3-cu-r1-5-b-d-07:52","s3-du-r1-6-o-d-08:30","s3-ri-a2-2-b-d-09:06","s3-ri-a2-4-v-d-09:35","s4-ca-p1-5-v-t-08:24","s4-ca-p2-10-v-t-12:50","s4-cu-a1-3-o-d-08:09","s4-cu-r1-12-o-d-09:49","s4-da-a2-13-o-d-09:09","s4-da-p2-11-v-t-10:40","s4-da-r1-3-v-d-09:32","s4-da-r2-6-v-d-07:53","s4-da-r2-6-v-t-07:46","s4-du-a2-15-v-t-12:59","s4-du-r1-4-o-d-09:05","s4-ri-a2-12-b-d-09:02","s4-ri-p1-6-v-d-08:15","s4-va-r2-14-v-d-10:46")

noclear<- c("s1-ca-r1-11-o-d-11:34","s1-da-p1-13-v-t-11:16","s2-du-r1-2-v-d-09:42","s3-da-p1-8-v-d-09:49")

bestHM<-c()

discard_artefacts<-c()

crop_reinspect<- c()

#Any to inspect?
co2_lowMAEimprovement %>% 
  filter(!UniqueID%in%c(bestLM,bestHM,noclear,discard_artefacts,crop_reinspect)) %>%
  select(UniqueID, g.fact,k.ratio, rel_reduct_MAE_HM) %>% head()

#Any inspected that shouldnt?
c(bestLM,noclear)[!c(bestLM,noclear)%in%co2_lowMAEimprovement$UniqueID]


#These criteria work well for all, a few low-g.fact HM unclear, but overall very good. 

#quality check: inspect k.ratio of selected HM and inspect visually.
k_final_good<-c("s1-ca-a1-3-v-t-09:10","s1-ca-a1-12-v-t-11:41","s1-du-a1-13-v-t-14:07","s2-ca-a1-1-o-d-08:47","s2-du-p2-9-b-d-11:59","s2-va-a1-10-v-t-13:21","s2-va-p2-1-o-d-09:04","s3-cu-r2-2-b-d-07:07","s3-da-a2-14-v-t-10:08","s3-cu-r2-14-o-d-12:19")

co2_todecide %>% 
  filter(AICc_weight_HM>0.5) %>% 
  filter(rel_reduct_MAE_HM>=0.05) %>% 
  filter(k.ratio>0.5) %>% 
  filter(!UniqueID%in%c(crop_reinspect, k_final_good)) %>% 
  select(UniqueID,g.fact,k.ratio)



#3. CO2 bestflux selection  -------

##3.1. Selection and flags------

co2_bestflux<- co2_flux %>% 
  #Create logical vectors for all cases, to be able to inspect (if any) cases where the best_model is overriden, due to errors in the order of the criteria. 
  mutate(is_discard=UniqueID%in%incub_discard,
         is_NA_HM=is.na(HM.flux),
         is_mdflinear=abs(LM.flux)<=MDF.lim,
         is_kmax=HM.k>=k.max,
         is_gfactexceed=g.fact>=4,
         is_lowmaeimprove=((LM.MAE  - HM.MAE)  / LM.MAE)<0.05,
         is_LMaiccbetter=LM.AICc<HM.AICc) %>% 
  #The conditions are evaluated sequentially (if multiple conditions are met, only the first is used), if no condition is met, default to "HM"
  mutate(best_model=case_when(is_discard~"None appropriate",
                              is_NA_HM~"LM",
                              is_mdflinear~"LM",
                              is_kmax~"LM",
                              is_gfactexceed~"LM",
                              is_lowmaeimprove~"LM",
                              is_LMaiccbetter~"LM",
                              TRUE~"HM"),
         # Add quality flags (some override others)
         best_model_flags = pmap_chr(
           list(is_discard, is_NA_HM,
                is_mdflinear, is_kmax, is_gfactexceed, is_lowmaeimprove, is_LMaiccbetter),
           function(discard, na_hm, mdf, kmax, gfact, lowmae, aicc) {
             # Priority override
             if (discard) return("Discard incubation, artefacts preclude flux estimate")
             
             # Goflux branch
             other_flags <- c()
             if (na_hm) other_flags <- c(other_flags, "HM.flux is NA")
             if (mdf) other_flags <- c(other_flags, "MDF")
             if (isTRUE(kmax)) other_flags <- c(other_flags, "HM.k exceeds maximum")
             if (isTRUE(gfact)) other_flags <- c(other_flags, "g.fact >= 4")
             if (isTRUE(lowmae)) other_flags <- c(other_flags, "HM does not improve MAE by > 5%")
             if (isTRUE(aicc)) other_flags <- c(other_flags, "AICc of LM is best")
             
             if (length(other_flags) == 0) return("HM meets all criteria")
             str_c(other_flags, collapse = " | ")
           }
         )
  ) %>% 
  mutate(best.flux=case_when(best_model=="LM"~LM.flux,
                             best_model=="HM"~HM.flux,
                             best_model=="total.flux"~total.flux,
                             best_model=="None appropriate"~NA_real_))


co2_bestflux %>%
  separate_rows(best_model_flags, sep = " \\| ") %>%
  count(best_model, best_model_flags)

co2_bestflux %>% 
  count(best_model)



##3.2. Final format ------
names(co2_bestflux)
#provisionally, only best.flux to be able to work with the data. 
co2_bestflux_formated<- co2_bestflux %>% 
  #general
  select(UniqueID, duration_s,flux.term,MDF.lim,
         #totalflux results
         C0.10s,C0.10s.sd,Cf.10s,Cf.10s.sd,total.flux, total.flux.se, 
         #LM results
         LM.flux, LM.SE, LM.MAE, LM.RMSE, LM.AICc, LM.r2, LM.p.val,
         #HM.results
         HM.flux, HM.SE, HM.MAE, HM.RMSE, HM.AICc, HM.r2, HM.k, k.max,
         #Fluxchoice results
         g.fact, best_model, best_model_flags,
  ) %>% 
  rename(duration=duration_s,
         LM.flux.se=LM.SE,
         HM.flux.se=HM.SE,
         best.model=best_model, best.model.flags=best_model_flags)

names(co2_bestflux_formated)

#save co2_bestflux_formated
write.csv(co2_bestflux_formated, file = paste0(results_path, "co2_bestflux.csv"),row.names = F)


#Check CO2 flux distribution according to models
co2_bestflux %>% 
  filter(!is.na(best.flux)) %>% 
  ggplot(aes( x=abs(best.flux), fill=best_model))+
  geom_histogram(bins = 100)+
  scale_x_log10()+
  geom_vline(xintercept=25)+
  facet_wrap(~best_model, nrow=3)

co2_bestflux %>% 
  filter(abs(best.flux)<1) %>% 
  filter(!is.na(best.flux)) %>% 
  ggplot(aes( x=(best.flux), fill=best_model))+
  geom_histogram(bins = 100)+
  # scale_x_log10()+
  # geom_vline(xintercept=25)+
  facet_wrap(~best_model, nrow=2, scales="free_y")

co2_bestflux_formated %>% 
  filter(best.model!="None appropriate") %>% 
  ggplot(aes(x=LM.flux.se, col=best.model))+
  geom_histogram()


#Quality check: inspect highest incubations with highest uncertainty for estimate

error_incub<- c()
ok<- c("s1-cu-r2-16-v-t-09:36","s1-da-p1-14-v-t-11:32","s1-da-p1-14-v-d-11:36","s3-cu-r2-11-v-d-11:28","s3-cu-r2-12-v-t-11:39","s3-da-p1-10-v-d-10:30","s4-da-a1-7-v-t-09:22","s4-du-a1-11-v-d-11:47","s4-du-a2-4-v-t-08:34","s4-du-r2-9-v-t-10:09","s4-du-r2-10-v-t-10:26","s1-da-p1-13-v-t-11:16","s3-cu-p1-2-v-t-06:41","s3-da-p1-14-v-t-11:49","s4-da-r2-4-v-t-07:21","s1-da-p1-13-v-d-11:22")

co2_bestflux_formated %>%   
  filter(best.model%in%c("HM","LM")) %>% 
  filter(LM.flux.se>quantile(LM.flux.se, 0.995)) %>% 
  filter(!UniqueID%in%c(ok,error_incub)) %>% 
  select(UniqueID, LM.flux.se,g.fact, best.model,best.model.flags)

