################### Defining stages ####################
stages = list()
## prenatal 
stages[["s1"]] = c("4 pcw", "7 pcw") # 1 4-7 pcw Embryonic
stages[["s2a"]] = c("8 pcw","9 pcw") # 2A 8-9 pcw Early prenatal
stages[["s2b"]] = c("12 pcw") # 2B 10-12 pcw Early prenatal
stages[["s3a"]] = c("13 pcw") # 3A 13-15 pcw Early mid-prenatal
stages[["s3b"]] = c("16 pcw","17 pcw") # 3B 16-18 pcw Early mid-prenatal
stages[["s4"]] = c("19 pcw","21 pcw","24 pcw") # 4 19-24 pcw Late mid-prenatal
stages[["s5"]] = c("25 pcw","26 pcw","35 pcw","37 pcw") # 5 25-38 pcw Late prenatal
stages[["s6"]] = c("4 mos") # 6 Birth-5 months Early infancy


## Postnatal 
stages[["s7"]] = c("10 mos", "1 yrs") # 7 6-18 months Late infancy
stages[["s8"]] = c("2 yrs","3 yrs","4 yrs") # 8 19 months-5 yrs Early childhood
stages[["s9"]] = c("8 yrs","11 yrs") # 9 6-11 yrs Late childhood
stages[["s10"]] = c("13 yrs","15 yrs","18 yrs","19 yrs") # 10 12-19 yrs Adolescence
stages[["s11"]] = c("21 yrs", "23 yrs")
stages[["s12"]] = c("30 yrs", "36 yrs", "37 yrs")
stages[["s13"]] = c("40 yrs")


#stages[["s11"]] = c("21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs") # 11 20-60+ yrs Adulthood
order.stages <- c("s1", "s2a", "s2b", "s3a", "s3b", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11", "s12", "s13")

stages

################### Defining age intervals ####################

order.age <- c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw",
               "4 mos","10 mos",
               "1 yrs","2 yrs","3 yrs","4 yrs","8 yrs",
               "11 yrs","13 yrs","15 yrs","18 yrs","19 yrs",
               "21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")


age_intervals = list()
#age_intervals[["8-9pcw"]] = c() # 1 4-7 pcw Embryonic
age_intervals[["4-7pcw"]] = paste0(c(4:7), " pcw")
age_intervals[["8-9pcw"]] = paste0(c(8:9), " pcw") # 2A 8-9 pcw Early prenatal
age_intervals[["10-12pcw"]] = paste0(c(10:12), " pcw") # 2B 10-12 pcw Early prenatal
age_intervals[["13-15pcw"]] = paste0(c(13:15), " pcw") # 3A 13-15 pcw Early mid-prenatal
age_intervals[["16-18pcw"]] = paste0(c(16:18), " pcw") # 3B 16-18 pcw Early mid-prenatal
age_intervals[["19-24pcw"]] = paste0(c(19:24), " pcw") # 4 19-24 pcw Late mid-prenatal
age_intervals[["25-38pcw"]] = paste0(c(25:38), " pcw") # 5 25-38 pcw Late prenatal
age_intervals[["39-40pcw"]] = paste0(c(39:40), " pcw") # 5 25-38 pcw Late prenatal
age_intervals[["0-5mos"]] = paste0(c(0:5), " mos") # 6 Birth-5 months Early infancy
age_intervals[["6-18mos"]] = c(paste0(c(6:18)," mos"), "1 yrs") # 7 6-18 months Late infancy
age_intervals[["19mos-5yrs"]] = c(paste0(c(19:20)," mos"), paste0(c(2:5), " yrs"))  # 8 19 months-5 yrs Early childhood
age_intervals[["6-11yrs"]] = paste0(c(6:11), " yrs") # 9 6-11 yrs Late childhood
age_intervals[["12-19yrs"]] = paste0(c(12:19), " yrs") # 10 12-19 yrs Adolescence
age_intervals[["20-29yrs"]] = paste0(c(20:29), " yrs") 
age_intervals[["30-39yrs"]] = paste0(c(30:39), " yrs")
age_intervals[["40-49yrs"]] = paste0(c(40:49), " yrs")
age_intervals[["50-59yrs"]] = paste0(c(50:59), " yrs")
age_intervals[["60-69yrs"]] = paste0(c(60:69), " yrs")
age_intervals[["70-79yrs"]] = paste0(c(70:79), " yrs")
age_intervals[["80-89yrs"]] = paste0(c(80:89), " yrs")
age_intervals[["90-99yrs"]] = paste0(c(90:99), " yrs")


order.intervals = c("8-9pcw", "10-12pcw", "13-15pcw", "16-18pcw",
                    "19-24pcw", "25-38pcw", "0-5mos", "6-18mos", "19mos-5yrs", 
                    "6-11yrs", "12-19yrs", "20-29yrs", "30-39yrs", "40-49yrs")



################### Defining period ####################
period = list()
period[["Prenatal"]] = c("8-9pcw", "10-12pcw", "13-15pcw", "16-18pcw",
                         "19-24pcw", "25-38pcw")


period[["Postnatal"]] = c("0-5mos", "6-18mos", "19mos-5yrs", 
                          "6-11yrs", "12-19yrs", "20-29yrs", "30-39yrs", "40- 49yrs")

order.period =c("Prenatal", "Postnatal")



############## Defining structures #####################
structure_acronym = list()

structure_acronym[["AMY"]] = c("Brain - Amygdala")
structure_acronym[["CB"]] = c("Brain - Cerebellum")
structure_acronym[["HYP"]] = c("Brain - Hypothalamus")
structure_acronym[["SNA"]] = c("Brain - Substantia nigra")
structure_acronym[["ACC"]] = c("Brain - Anterior cingulate cortex (BA24)")
structure_acronym[["CTX"]] = c("Brain - Cortex")
structure_acronym[["NAC"]] = c("Brain - Nucleus accumbens (basal ganglia)")
structure_acronym[["CGE"]] = c("Brain - Caudate (basal ganglia)")
structure_acronym[["DLPFC"]] = c("Brain - Frontal Cortex (BA9)")
structure_acronym[["PUT"]] = c("Brain - Putamen (basal ganglia)")
structure_acronym[["CBC"]] = c("Brain - Cerebellar Hemisphere")
structure_acronym[["HIP"]] = c("Brain - Hippocampus")
structure_acronym[["SCI"]] = c("Brain - Spinal cord (cervical c-1)")


################### Defining regions ####################
regions = list()

regions[["Subcortex"]] = c("AMY", "CGE", "DTH", "HIP", "LGE", "MGE", 
                           "MD", "STR", "SNA", "PUT", "HYP", "NAC")

regions[["Cortex"]] = c("MFC", "DFC", "Ocx", "OFC", "PCx", 
                        "TCx", "VFC", "ITC", "STC", "IPC", "V1C", 
                        "M1C", "M1C-S1C", "S1C", "A1C", "ACC", "CTX", "DLPFC", 
                        "MTG", "CgG", "M1lm", "M1ul", "S1lm", "S1ul")

regions[["Cerebellum"]] = c("CBC", "CB", "URL")

regions[["Spinal Cord"]] = c("SCI")
