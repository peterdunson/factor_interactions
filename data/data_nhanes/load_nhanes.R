## Reading in the NHANES datasets

# Read in the Phtalate information
library(SASxport)
library(tidyverse)

setwd("~/Desktop/rechemicalsfriends")

years = c("0910", "1112", "1314") 

get_files <- function(str, years){
  files = system(paste("ls", str, sep=" "), intern=TRUE)
  tmp = files[map_dbl(years, function(y) ifelse(length(str_which(files, y)) == 1, str_which(files, y), NA))]
  tmp[!is.na(tmp)]
}

get_data <- function(files){
  map_dfr(files, function(f) as.data.frame(as.matrix(read.xport(f))))
}

# Get the different dataset files
strings = c("PHT*.XPT", "BMX*.XPT", "DEMO*.XPT","ALB*.XPT",
            "PP*.XPT","UHM*.XPT","EPH*.XPT","TCHO*.XPT")

names = c("Phthlates", "Body_Measures", "Demographics","Creatinine",
          "Organochlorine_Pesticides", "Metals","Phenols","Cholesterol")

data = tibble(names, strings) %>% 
  mutate(files = map(strings, function(s) get_files(s, years)),
         data = map(files, get_data))

phthalates_nhanes = c("URXMBP", "URXMIB", "URXMEP", 
                      "URXMZP", "URXMCP", "URXECP", 
                      "URXMHH", "URXMOH", "URXMHP")
# URXMBP: Mono-n-butyl phthalate (ng/mL)
# URXMIB: Mono-isobutyl phthalate
# URXMEP: Mono-ethyl phthalate (ng/mL)
# URXMZP: Mono-benzyl phthalate (ng/mL)
# URXMCP: Mono-cyclohexyl phthalate (ng/mL)
# URXECP: Mono-2-ethyl-5-carboxypentyl phthalate
# URXMHH: Mono-(2-ethyl-5-hydroxyhexyl) phthalate
# URXMOH: Mono-(2-ethyl-5-oxohexyl) phthalate
# URXMHP: Mono-(2-ethyl)-hexyl phthalate (ng/mL)

chem_nhanes = data %>%
  filter(names == "Phthlates") %>%
  select(data) %>%
  unnest() %>%
  select(SEQN, phthalates_nhanes) %>%
  filter_all(all_vars(!is.na(.)))
# SEQN: Respondent sequence number

people_nhanes = data %>%
  filter(names == "Demographics") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(c(RIDAGEYR, RIAGENDR, RIDRETH1, INDFMPIR)) 
# RIDAGEYR: Age in years at screening
# RIAGENDR: Gender
# RIDRETH1: Race/Hispanic origin 
# INDFMPIR: Ratio of family income to poverty

bmi_nhanes = data %>%
  filter(names == "Body_Measures") %>%
  select(data) %>%
  unnest() %>%
  filter(SEQN %in% chem_nhanes$SEQN) %>%
  select(BMXBMI, BMXWAIST)
# BMXBMI: Body Mass Index (kg/m**2)
# BMXWAIST: Waist Circumference (cm)

creatinine_nhanes = data %>%
   filter(names == "Creatinine") %>%
   select(data) %>%
   unnest() %>%
   filter(SEQN %in% chem_nhanes$SEQN) %>%
   select(URXUCR)
# URXUCR: Creatinine, urine (mg/dL)

# Only 2010-2011 and 2013-2014
pesticides_nhanes = data %>%
   filter(names == "Organochlorine_Pesticides") %>%
   select(data) %>%
   unnest() %>%
   filter(SEQN %in% chem_nhanes$SEQN) %>%
   select(URX14D,URXDCB)
# URX14D: 2,5-dichlorophenol (ug/L) result
# URXDCB: 2,4-dichlorophenol (ug/L) result



# people that have Phalates do not have metals
metals_nhanes = data %>%
   filter(names == "Metals") %>%
   select(data) %>%
   unnest() %>%
   #filter(SEQN %in% chem_nhanes$SEQN) %>%
   select(URXUBA,URXUBE,URXUCD,
          URXUCO,URXUCS,URXUMO,
          URXUPB,URXUPT,URXUSB,
          URXUTL,URXUTU,URXUUR)
# URXUBA: Barium, urine (ug/L)
# URXUBE: Beryllium, urine (ug/L)
# URXUCD: Cadmium, urine (ug/L)
# URXUCO: Cobalt, urine (ug/L)
# URXUCS: Cesium, urine (ug/L)
# URXUMO: Molybdenum, urine (ug/L)
# URXUPB: Lead, urine (ug/L)
# URXUPT: Platinum, urine (ug/L)
# URXUSB: Antimony, urine (ug/L)
# URXUTL: Thallium, urine (ug/L)
# URXUTU: Tungsten, urine (ug/L)
# URXUUR: Uranium, urinary (ug/L)

# Missing year 2013-2014
phenols_nhanes = data %>%
   filter(names == "Phenols") %>%
   select(data) %>%
   unnest() %>%
   filter(SEQN %in% chem_nhanes$SEQN) %>%
   select(URXBPH,URX4TO,URXTRS,URXBP3,
          URXPPB,URXBUP,URXEPB,URXMPB)
# URXBPH: Urinary Bisphenol A
# URX4TO: Urinary 4-tert-octyl 
# URXTRS: Urinary 2,4,4’-Trichloro-2’-hydroxyphenyl ether (Triclosan)		
# URXBP3: Urinary 2-Hydroxy-4-metoxybenzophenone (Benzophenone-3)	
# URXPPB: Urinary Propyl paraben	
# URXBUP: Urinary Butyl paraben		
# URXEPB: Urinary Ethyl paraben	
# URXMPB: Urinary Methyl paraben	

chol_nhanes = data %>%
   filter(names == "Cholesterol") %>%
   select(data) %>%
   unnest() %>%
   filter(SEQN %in% chem_nhanes$SEQN) %>%
   select(LBXTC)
# LBXTC: Total Cholesterol (mg/dL)


data_nhanes = data.frame(chem_nhanes,people_nhanes,bmi_nhanes,
                         creatinine_nhanes,chol_nhanes,pesticides_nhanes)
save(data_nhanes, file = "nhanes.RData")

