#' @title Fetch NHANES data 
#' @description Script to fetch data from NHANES website for later analysis
#' @author Dr. James Rooney, \email{jrooney@@rcsi.ie}
#' 
#' 

#### setup
library(tidyverse)
library(foreign)
#library(dplyr)

### Download up to date file and variable names from NHANES
# New files are released on an ongiong basis

#files <- nhanes_data_files()
#variables <- nhanes_variables()


### Extract data from 2015-2016 NHANES data release
# Demographics                                                      "DEMO_I.XPT"
# BMI filename                                                      "BMX_I.XPT"
# Blood Pb, Cd, Total Hg, Se, Mg                                    "PBCD_I.XPT"
# Blood mercury: inorganic, ethyl and methyl                        "IHGEM_I.XPT"
# Chromium & Cobalt                                                 "CRCO_I.XPT"
# Urine albumin & creatinine                                        "ALB_CR_I.XPT"
# Cholesterol - High-Density Lipoprotein (HDL)                      "HDL_I.XPT"
# Cholesterol - Total                                               "TCHOL_I.XPT"
# Complete Blood Count with 5-Part Differential - Whole Blood       "CBC_I.XPT"
# Fasting Questionnaire                                             "FASTQX_I.XPT"
# Glycohemoglobin                                                   "GHB_I.XPT"
# Insulin                                                           "INS_I .XPT"
# Oral Glucose Tolerance Test                                       "OGTT_I.XPT"
# Plasma Fasting Glucose                                            "GLU_I.XPT"
# Pregnancy Test - Urine                                            "UCPREG_I.XPT"
# Standard Biochemistry Profile                                     "BIOPRO_I.XPT"
# Urine Flow Rate                                                   "UCFLOW_I.XPT"
# Oral health Dentition                                             "OHXDEN_I.XPT"

# Iodine - Urine                                                    "UIO_I.XPT"
# Mercury - Urine                                                   "UHG_I.XPT"
# Metals - Urine                                                    "UM_I.XPT"
# Metals - Urine - Special Sample (Subsample)                       "UMS_I.XPT"
# Speciated Arsenics - Urine - Special Sample (Subsample)           "UASS_I.XPT"
# Speciated Arsenics - Urine (Subsample)                            "UAS_I.XPT"
# Fluoride - Plasma                                                 "FLDEP_I.XPT"
# Fluoride - Water                                                  "FLDEW_I.XPT"
# Folate - whole blood                                              "FOLATE_I.XPT"

##### Additional files Nov2020 ####
# Acrylamide & Glycidamide                                          "AMDGYD_I.XPT"
# Apolipoprotein B                                                  "APOB_I.XPT"
# Aromatic Diamines                                                 "UADM_I.XPT"
# Brominated Flame Retardants (BFRs) - Pooled Samples               "BFRPOL_I.XPT"
# Chlamydia - Urine                                                 "CHLMDA_I.XPT"
# Chlamydia Pgp3                                                    "SSCT_I.XPT"
# Cotinine and Hydroxycotinine - Serum                              "COT_I.XPT"
# Cotinine, Hydroxycotinine, & Other ......<SNIP>                   "UCOT_I.XPT"
# Cotinine, Hydroxycotinine, .... Special Sample (Subsample)        "UCOTS_I.XPT"
# DEET Metabolite - Urine                                           "DEET_I.XPT"
# DEET Metabolites - Urine - Surplus                                "SSDEET_I.XPT"
# Ethylene Oxide                                                    "ETHOX_I.XPT"
# Formaldehyde                                                      "FORMAL_I.XPT"
# Hepatitis A                                                       "HEPA_I.XPT"
# Hepatitis B: Core antibody, Surface antigen, and HepD antibody    "HEPBD_I.XPT"
# Hepatitis B: Surface Antibody                                     "HEPB_S_I.XPT"
# Hepatitis C: RNA (HCV-RNA) and Hepatitis C Genotype               "HEPC_I.XPT"
# Hepatitis E: IgG & IgM Antibodies                                 "HEPE_I.XPT"
# Herpes Simplex Virus Type-1 & Type-2                              "HSV_I.XPT"
# HIV Antibody Test                                                 "HIV_I.XPT"
# Human Papillomavirus (HPV) - Oral Rinse                           "ORHPV_I.XPT"
# Imidacloprid, 5-Hydroxy imidacloprid, Acetamiprid, N-desmethyl Acetamiprid, Clothianidin, and Thiacloprid in NHANES 2015-16 Surplus Urine                                     "SSNEON_I.XPT"
# Mono-2-ethyl-5-hydroxyhexyl terephthalate, mono-2-ethyl-5-carboxypentyl terephthalate, and monooxoisononyl phthalate - Urine (Surplus)                                         "SSMHHT_I.XPT"
# Non-dioxin-like Polychlorinated Biphenyls & Mono-ortho-substituted Polychlorinated Biphenyls - Pooled Samples                                                             "PCBPOL_I.XPT"
# Perchlorate, Nitrate & Thiocyanate - Urine                        "PERNT_I.XPT"
# Perchlorate, Nitrate & Thiocyanate - Urine - Special Sample       "PERNTS_I.XPT"
# Perfluoroalkyl and Polyfluoroalkyl                                "PFAS_I.XPT"
# Personal Care and Consumer Product Chemicals and Metabolites      "EPHPP_I.XPT"
# Pesticides - Organochlorine Pesticides - Serum - Pooled Samples   "PSTPOL_I.XPT"
# Phthalates and Plasticizers Metabolites - Urine                   "PHTHTE_I.XPT"
# Polycyclic Aromatic Hydrocarbons (PAH) - Urine                    "PAH_I.XPT"
# Polycyclic Aromatic Hydrocarbons (PAH) - Urine - Special Sample   "PAHS_I.XPT"
# Pooled-Sample Technical Support File                              "POOLTF_I.XPT"
# Sex Steroid Hormone - Serum                                       "TST_I.XPT"
# Transferrin Receptor                                              "TFR_I.XPT"
# Trichomonas - Urine                                               "TRICH_I.XPT"
# Urine Flow Rate                                                   "UCFLOW_I.XPT"
# Volatile Organic Compound (VOC) Metabolites - Urine               "UVOC_I.XPT"
# Volatile Organic Compound (VOC) Metabolites - Urine - Spec Samp   "UVOCS_I.XPT"
# Volatile Organic Compounds and Trihalomethanes/MTBE - Blood       "VOCWB_I.XPT"
# Volatile Organic Compounds and Trihalomethanes/MTBE – Blood – Special Sample "VOCWBS_I.XPT"

# blood pressure data                                               "BPX_I.XPT"


#### download NHANES files
extra_files <- c("AMDGYD_I.XPT", "APOB_I.XPT", "UADM_I.XPT", "BFRPOL_I.XPT", "CHLMDA_I.XPT",
                 "SSCT_I.XPT", "COT_I.XPT", "UCOT_I.XPT", "UCOTS_I.XPT", "DEET_I.XPT",
                 "SSDEET_I.XPT", "ETHOX_I.XPT", "FORMAL_I.XPT", "HEPA_I.XPT", "HEPBD_I.XPT",
                 "HEPB_S_I.XPT", "HEPC_I.XPT", "HEPE_I.XPT", "HSV_I.XPT", "HIV_I.XPT",
                 "ORHPV_I.XPT", "SSNEON_I.XPT", "PCBPOL_I.XPT", "PERNT_I.XPT", "PERNTS_I.XPT", 
                 "PFAS_I.XPT", "EPHPP_I.XPT", "PSTPOL_I.XPT", "PHTHTE_I.XPT", "PAH_I.XPT", 
                 "PAHS_I.XPT", "POOLTF_I.XPT", "TST_I.XPT", "TFR_I.XPT", "TRICH_I.XPT", 
                 "UCFLOW_I.XPT", "UVOC_I.XPT", "UVOCS_I.XPT", "VOCWB_I.XPT", "VOCWBS_I.XPT")

filestem <- "https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/"
dl_filenames <- c("DEMO_I.XPT", "BMX_I.XPT", "PBCD_I.XPT", "IHGEM_I.XPT", "CRCO_I.XPT",
                  "ALB_CR_I.XPT", "HDL_I.XPT", "TCHOL_I.XPT", "CBC_I.XPT", "FASTQX_I.XPT",
                  "GHB_I.XPT", "INS_I.XPT", "OGTT_I.XPT", "GLU_I.XPT", "UCPREG_I.XPT",
                  "BIOPRO_I.XPT", "UCFLOW_I.XPT", "OHXDEN_I.XPT", "UIO_I.XPT", "UHG_I.XPT",
                  "UM_I.XPT", "UMS_I.XPT", "UASS_I.XPT", "UAS_I.XPT", "FLDEP_I.XPT",
                  "FLDEW_I.XPT", "FOLATE_I.XPT", extra_files, "BPX_I.XPT")


for (i in 1: length(dl_filenames)){
    URL <- paste0(filestem, dl_filenames[i] )
    destfile <- paste0("Data/2015_2016/", dl_filenames[i])
    
    if(!file.exists(destfile)){
        res <- tryCatch(download.file(URL, destfile=destfile, method="libcurl"),
                        error=function(e) 1)
        print(destfile)
    }
}



# Load downloaded files
df_dem <- read.xport("Data/2015_2016/DEMO_I.XPT")
df_bmi <- read.xport("Data/2015_2016/BMX_I.XPT")
df_pbcd <- read.xport("Data/2015_2016/PBCD_I.XPT")
df_ihgem <- read.xport("Data/2015_2016/IHGEM_I.XPT")
df_crco <- read.xport("Data/2015_2016/CRCO_I.XPT")
df_alb_cr <- read.xport("Data/2015_2016/ALB_CR_I.XPT")
df_chol_hdl <- read.xport("Data/2015_2016/HDL_I.XPT")
df_tchol <- read.xport("Data/2015_2016/TCHOL_I.XPT")
df_fbc <- read.xport("Data/2015_2016/CBC_I.XPT")
#df_fast_q <- read.xport("Data/2015_2016/FASTQX_I.XPT")
df_glyco_hb <- read.xport("Data/2015_2016/GHB_I.XPT")
df_insulin <- read.xport("Data/2015_2016/INS_I.XPT")
df_ogtt <- read.xport("Data/2015_2016/OGTT_I.XPT")
df_glu <- read.xport("Data/2015_2016/GLU_I.XPT")
df_preg <- read.xport("Data/2015_2016/UCPREG_I.XPT")
df_biochem <- read.xport("Data/2015_2016/BIOPRO_I.XPT")
df_ur_flow <- read.xport("Data/2015_2016/UCFLOW_I.XPT")
df_oralhealth <- read.xport("Data/2015_2016/OHXDEN_I.XPT")

df_ur_iod <- read.xport("Data/2015_2016/UIO_I.XPT")
df_ur_hg <- read.xport("Data/2015_2016/UHG_I.XPT")
df_ur_metals <- read.xport("Data/2015_2016/UM_I.XPT")
#df_ur_metals_special <- read.xport("Data/2015_2016/UMS_I.XPT")
df_ur_spec_ars <- read.xport("Data/2015_2016/UASS_I.XPT")
#df_ur_spec_ars_sub <- read.xport("Data/2015_2016/UAS_I.XPT")
df_pl_fluor <- read.xport("Data/2015_2016/FLDEP_I.XPT")
df_water_fluor <- read.xport("Data/2015_2016/FLDEW_I.XPT")
df_folate <- read.xport("Data/2015_2016/FOLATE_I.XPT")
df_sexhorm <- read.xport("Data/2015_2016/TST_I.XPT")
df_cot <- read.xport("Data/2015_2016/COT_I.XPT")

# read blood pressure data
df_bp <- read.xport("Data/2015_2016/BPX_I.XPT")


# Need to recode sample comments before merging fields
df_pbcd$LBDBPBLC <- ifelse( is.na(df_pbcd$LBDBPBLC), "Missing",
                        ifelse(df_pbcd$LBDBPBLC == 0, "Valid",
                           ifelse(df_pbcd$LBDBPBLC == 1, "Below_Lod",
                                  NA)))
df_pbcd$LBDBCDLC <- ifelse( is.na(df_pbcd$LBDBCDLC), "Missing",
                            ifelse(df_pbcd$LBDBCDLC == 0, "Valid",
                                   ifelse(df_pbcd$LBDBCDLC == 1, "Below_Lod",
                                          NA)))
df_pbcd$LBDTHGLC <- ifelse( is.na(df_pbcd$LBDTHGLC), "Missing",
                            ifelse(df_pbcd$LBDTHGLC == 0, "Valid",
                                   ifelse(df_pbcd$LBDTHGLC == 1, "Below_Lod",
                                          NA)))
df_pbcd$LBDBSELC <- ifelse( is.na(df_pbcd$LBDBSELC), "Missing",
                            ifelse(df_pbcd$LBDBSELC == 0, "Valid",
                                   ifelse(df_pbcd$LBDBSELC == 1, "Below_Lod",
                                          NA)))
df_pbcd$LBDBMNLC <- ifelse( is.na(df_pbcd$LBDBMNLC), "Missing",
                            ifelse(df_pbcd$LBDBMNLC == 0, "Valid",
                                   ifelse(df_pbcd$LBDBMNLC == 1, "Below_Lod",
                                          NA)))
df_ihgem$LBDIHGLC <- ifelse( is.na(df_ihgem$LBDIHGLC), "Missing",
                            ifelse(df_ihgem$LBDIHGLC == 0, "Valid",
                                   ifelse(df_ihgem$LBDIHGLC == 1, "Below_Lod",
                                          NA)))
df_ihgem$LBDBGELC <- ifelse( is.na(df_ihgem$LBDBGELC), "Missing",
                             ifelse(df_ihgem$LBDBGELC == 0, "Valid",
                                    ifelse(df_ihgem$LBDBGELC == 1, "Below_Lod",
                                           NA)))
df_ihgem$LBDBGMLC <- ifelse( is.na(df_ihgem$LBDBGMLC), "Missing",
                             ifelse(df_ihgem$LBDBGMLC == 0, "Valid",
                                    ifelse(df_ihgem$LBDBGMLC == 1, "Below_Lod",
                                           NA)))
df_ur_iod$URDUIOLC <- ifelse( is.na(df_ur_iod$URDUIOLC), "Missing",
                              ifelse(df_ur_iod$URDUIOLC == 0, "Valid",
                                     ifelse(df_ur_iod$URDUIOLC == 1, "Below_Lod",
                                            NA)))
df_ur_hg$URDUHGLC <- ifelse( is.na(df_ur_hg$URDUHGLC), "Missing",
                             ifelse(df_ur_hg$URDUHGLC == 0, "Valid",
                                    ifelse(df_ur_hg$URDUHGLC == 1, "Below_Lod",
                                           NA)))
df_ur_metals$URDUBALC <- ifelse( is.na(df_ur_metals$URDUBALC), "Missing",
                                 ifelse(df_ur_metals$URDUBALC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUBALC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUCDLC <- ifelse( is.na(df_ur_metals$URDUCDLC), "Missing",
                                 ifelse(df_ur_metals$URDUCDLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUCDLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUCOLC <- ifelse( is.na(df_ur_metals$URDUCOLC), "Missing",
                                 ifelse(df_ur_metals$URDUCOLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUCOLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUCSLC <- ifelse( is.na(df_ur_metals$URDUCSLC), "Missing",
                                 ifelse(df_ur_metals$URDUCSLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUCSLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUMOLC <- ifelse( is.na(df_ur_metals$URDUMOLC), "Missing",
                                 ifelse(df_ur_metals$URDUMOLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUMOLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUMNLC <- ifelse( is.na(df_ur_metals$URDUMNLC), "Missing",
                                 ifelse(df_ur_metals$URDUMNLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUMNLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUPBLC <- ifelse( is.na(df_ur_metals$URDUPBLC), "Missing",
                                 ifelse(df_ur_metals$URDUPBLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUPBLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUSBLC <- ifelse( is.na(df_ur_metals$URDUSBLC), "Missing",
                                 ifelse(df_ur_metals$URDUSBLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUSBLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUSNLC <- ifelse( is.na(df_ur_metals$URDUSNLC), "Missing",
                                 ifelse(df_ur_metals$URDUSNLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUSNLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUSRLC <- ifelse( is.na(df_ur_metals$URDUSRLC), "Missing",
                                 ifelse(df_ur_metals$URDUSRLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUSRLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUTLLC <- ifelse( is.na(df_ur_metals$URDUTLLC), "Missing",
                                 ifelse(df_ur_metals$URDUTLLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUTLLC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUTULC <- ifelse( is.na(df_ur_metals$URDUTULC), "Missing",
                                 ifelse(df_ur_metals$URDUTULC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUTULC == 1, "Below_Lod",
                                               NA)))
df_ur_metals$URDUURLC <- ifelse( is.na(df_ur_metals$URDUURLC), "Missing",
                                 ifelse(df_ur_metals$URDUURLC == 0, "Valid",
                                        ifelse(df_ur_metals$URDUURLC == 1, "Below_Lod",
                                               NA)))
df_crco$LBDBCRLC <- ifelse( is.na(df_crco$LBDBCRLC), "Missing",
                                 ifelse(df_crco$LBDBCRLC == 0, "Valid",
                                        ifelse(df_crco$LBDBCRLC == 1, "Below_Lod",
                                               NA)))
df_crco$LBDBCOLC <- ifelse( is.na(df_crco$LBDBCOLC), "Missing",
                            ifelse(df_crco$LBDBCOLC == 0, "Valid",
                                   ifelse(df_crco$LBDBCOLC == 1, "Below_Lod",
                                          NA)))
df_alb_cr$URDUMALC <- ifelse( is.na(df_alb_cr$URDUMALC), "Missing",
                            ifelse(df_alb_cr$URDUMALC == 0, "Valid",
                                   ifelse(df_alb_cr$URDUMALC == 1, "Below_Lod",
                                          NA)))
df_alb_cr$URDUCRLC <- ifelse( is.na(df_alb_cr$URDUCRLC), "Missing",
                              ifelse(df_alb_cr$URDUCRLC == 0, "Valid",
                                     ifelse(df_alb_cr$URDUCRLC == 1, "Below_Lod",
                                            NA)))
df_insulin$LBDINLC <- ifelse( is.na(df_insulin$LBDINLC), "Missing",
                              ifelse(df_insulin$LBDINLC == 0, "Valid",
                                     ifelse(df_insulin$LBDINLC == 1, "Below_Lod",
                                            NA)))
df_sexhorm$LBDTSTLC <- ifelse( is.na(df_sexhorm$LBDTSTLC), "Missing",
                      ifelse(df_sexhorm$LBDTSTLC == 0, "Valid",
                             ifelse(df_sexhorm$LBDTSTLC == 1, "Below_Lod",
                                    NA)))
df_sexhorm$LBDESTLC <- ifelse( is.na(df_sexhorm$LBDESTLC), "Missing",
                               ifelse(df_sexhorm$LBDESTLC == 0, "Valid",
                                      ifelse(df_sexhorm$LBDESTLC == 1, "Below_Lod",
                                             NA)))
df_sexhorm$LBDSHGLC <- ifelse( is.na(df_sexhorm$LBDSHGLC), "Missing",
                               ifelse(df_sexhorm$LBDSHGLC == 0, "Valid",
                                      ifelse(df_sexhorm$LBDSHGLC == 1, "Below_Lod",
                                             NA)))
df_cot$LBDCOTLC <- ifelse( is.na(df_cot$LBDCOTLC), "Missing",
                           ifelse(df_cot$LBDCOTLC == 0, "Valid",
                                  ifelse(df_cot$LBDCOTLC == 1, "Below_Lod",
                                        ifelse(df_cot$LBDCOTLC == 2, "OutOfRange",
                                         NA))))
df_cot$LBDHCTLC <- ifelse( is.na(df_cot$LBDHCTLC), "Missing",
                           ifelse(df_cot$LBDHCTLC == 0, "Valid",
                                  ifelse(df_cot$LBDHCTLC == 1, "Below_Lod",
                                         ifelse(df_cot$LBDHCTLC == 2, "OutOfRange",
                                                NA))))

## Merge downloaded data
df_nhanes <- dplyr::left_join(df_dem, df_bmi, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_pbcd, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_ihgem, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_crco, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_alb_cr, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_tchol, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_fbc, by="SEQN")
#df_nhanes <- dplyr::left_join(df_nhanes, df_fast_q, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_glyco_hb, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_insulin, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_ogtt, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_glu, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_preg, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_biochem, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_ur_flow, by="SEQN")

df_nhanes <- dplyr::left_join(df_nhanes, df_ur_iod, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_ur_hg, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_ur_metals, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_pl_fluor, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_water_fluor, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_folate, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_sexhorm, by="SEQN")
df_nhanes <- dplyr::left_join(df_nhanes, df_cot, by="SEQN")

#df_nhanes <- dplyr::left_join(df_nhanes, df_oralhealth, by="SEQN")



## Two columns shared a name despite having different data - rename now
df_nhanes <- df_nhanes %>%
    dplyr::rename("Hct_pct" = "LBXHCT.x",
        "serum_hydrocot_ng_ml" = "LBXHCT.y")


# Remove duplicate columns
df_nhanes$WTSA2YR.x <- NULL
df_nhanes$WTSA2YR.y <- NULL

df_nhanes$PHAFSTHR.x <- NULL
df_nhanes$PHAFSTHR.y <- NULL

df_nhanes$PHAFSTMN.x <- NULL
df_nhanes$PHAFSTMN.y <- NULL

df_nhanes$WTSH2YR <- df_nhanes$WTSH2YR.x
df_nhanes$WTSH2YR.x <- NULL
df_nhanes$WTSH2YR.y <- NULL

df_nhanes$WTSAF2YR <- df_nhanes$WTSAF2YR.x
df_nhanes$WTSAF2YR.x <- NULL
df_nhanes$WTSAF2YR.y <- NULL


# Reformat survey weight variables
df_nhanes <- df_nhanes %>%
    dplyr::rename("WT_full2yr" = "WTINT2YR",
                  "WT_mec2yr"= "WTMEC2YR",
                  "WT_subsampleA_2yr" = "WTSA2YR",
                  "WT_fasting_subsamp_mec2yr" = "WTSAF2YR",
                  "WT_blood_metal" = "WTSH2YR",
                  "WT_oggt_subsamp_mec2yr" = "WTSOG2YR")

df_nhanes$WT_mec2yr_comment <- factor( ifelse(df_nhanes$WT_mec2yr == 0, "Not MEC examined",
                                      ifelse(df_nhanes$WT_mec2yr > 0, "Examined",
                                             NA)))
df_nhanes$WT_subsampleA_2yr_comment <- factor( ifelse( is.na(df_nhanes$WT_subsampleA_2yr), "Not selected",
                                ifelse( df_nhanes$WT_subsampleA_2yr ==0, "Selected but no sample/result",
                                ifelse( df_nhanes$WT_subsampleA_2yr > 0 , "Valid sample",
                                NA))) )
df_nhanes$WT_blood_metal_comment <- factor( ifelse( is.na(df_nhanes$WT_blood_metal), "Not selected",
                                ifelse( df_nhanes$WT_blood_metal ==0, "Selected but no sample/result",
                                ifelse( df_nhanes$WT_blood_metal > 0, "Valid sample",
                                NA))) )
df_nhanes$WT_fasting_subsamp_mec2yr_comment <- factor(
                            ifelse( is.na(df_nhanes$WT_fasting_subsamp_mec2yr), "Not selected",
                            ifelse( df_nhanes$WT_fasting_subsamp_mec2yr ==0, "Selected but no sample/result",
                            ifelse( df_nhanes$WT_fasting_subsamp_mec2yr > 0, "Valid sample",
                                    NA))) )
df_nhanes$WT_oggt_subsamp_mec2yr_comment <- factor(
                            ifelse( is.na(df_nhanes$WT_oggt_subsamp_mec2yr), "Not selected",
                            ifelse( df_nhanes$WT_oggt_subsamp_mec2yr ==0, "Selected but no sample/result",
                            ifelse( df_nhanes$WT_oggt_subsamp_mec2yr > 0, "Valid sample",
                                    NA))) )


# Format variables of interest for clarity
df_nhanes <- df_nhanes %>%
                dplyr::rename("exam_status" = "RIDSTATR",
                       "gender" = "RIAGENDR",
                       "age_screen_yrs" = "RIDAGEYR",
                       "age_screen_mnths" = "RIDAGEMN",
                       "race1" = "RIDRETH1",
                       "race3" = "RIDRETH3",
                       "exam_month" = "RIDEXMON",
                       "age_exam_mnths" = "RIDEXAGM",
                       "military"= "DMQMILIZ",
                       "served_abroad"= "DMQADFC",
                       "country_birth"= "DMDBORN4",
                       "citizenship" = "DMDCITZN",
                       "yrs_in_US" = "DMDYRSUS",
                       "educ_children"= "DMDEDUC3",
                       "educ_adults" = "DMDEDUC2",
                       "marital_status" = "DMDMARTL",
                       "preg_status_exam" = "RIDEXPRG",
                       "language"= "SIALANG",
                       "proxy" = "SIAPROXY",
                       "interpreter" = "SIAINTRP",
                       "lang_family" = "FIALANG",
                       "proxy_family" = "FIAPROXY",
                       "interp_family" = "FIAINTRP",
                       "lang_mec" = "MIALANG",
                       "proxy_mec" = "MIAPROXY",
                       "interp_mec" = "MIAINTRP",
                       "lang_acasi" = "AIALANGA",
                       "num_house" = "DMDHHSIZ",
                       "num_family" = "DMDFMSIZ",
                       "num_under5" = "DMDHHSZA",
                       "num_6to17" = "DMDHHSZB",
                       "num_over60"= "DMDHHSZE",
                       "HH_gender" = "DMDHRGND",
                       "HH_age" = "DMDHRAGE",
                       "HH_birth_country" = "DMDHRBR4",
                       "HH_educ" = "DMDHREDU",
                       "HH_marriage" = "DMDHRMAR",
                       "HH_spouse_educ" = "DMDHSEDU",
                       "house_income" = "INDHHIN2",
                       "family_income" = "INDFMIN2",
                       "income_2poverty" = "INDFMPIR"
                        )

df_nhanes$exam_status <- factor( ifelse(df_nhanes$exam_status == 1, "Interview", "Interview + MEC") )
df_nhanes$gender <- factor( ifelse(df_nhanes$gender ==1, "Male", "Female"))
df_nhanes$race1 <- factor( ifelse(df_nhanes$race1 ==1, "Mexican American",
                                ifelse(df_nhanes$race1 ==2, "Other Hispanic",
                                ifelse(df_nhanes$race1 ==3, "Non-Hispanic White",
                                ifelse(df_nhanes$race1 ==4, "Non-Hispanic Black",
                                ifelse(df_nhanes$race1 ==5, "Other Race - Including Multi-Racial",
                                       NA))))))
df_nhanes$race3 <- factor( ifelse(df_nhanes$race3 ==1, "Mexican American",
                                ifelse(df_nhanes$race3 ==2, "Other Hispanic",
                                ifelse(df_nhanes$race3 ==3, "Non-Hispanic White",
                                ifelse(df_nhanes$race3 ==4, "Non-Hispanic Black",
                                ifelse(df_nhanes$race3 ==6, "Non-Hispanic Asian",
                                ifelse(df_nhanes$race3 ==7, "Other Race - Including Multi-Racial",
                                              NA)))))))
df_nhanes$exam_month <- factor( ifelse(df_nhanes$exam_month == 1, "November 1 through April 30",
                               "May 1 through October 31") )
df_nhanes$military <- factor( ifelse(df_nhanes$military == 1, "Yes",
                                    ifelse(df_nhanes$military == 2, "No",
                                    ifelse(df_nhanes$military == 7, "Refused",
                                    ifelse(df_nhanes$military == 9, "Don't know",
                                     NA)))))
df_nhanes$served_abroad <- factor( ifelse(df_nhanes$served_abroad == 1, "Yes",
                                    ifelse(df_nhanes$served_abroad == 2, "No",
                                    ifelse(df_nhanes$served_abroad == 7, "Refused",
                                    ifelse(df_nhanes$served_abroad == 9, "Don't know",
                                          NA)))))
df_nhanes$country_birth <- factor( ifelse(df_nhanes$country_birth == 1, "US_50states",
                                    ifelse(df_nhanes$country_birth == 2, "Others",
                                    ifelse(df_nhanes$country_birth == 77, "Refused",
                                    ifelse(df_nhanes$country_birth == 99, "Don't know",
                                          NA)))))
df_nhanes$citizenship <- factor( ifelse(df_nhanes$citizenship == 1, "US_citizen",
                                    ifelse(df_nhanes$citizenship == 2, "Not_US_citizen",
                                    ifelse(df_nhanes$citizenship == 7, "Refused",
                                    ifelse(df_nhanes$citizenship == 9, "don't know",
                                        NA)))))
df_nhanes$yrs_in_US <- factor( ifelse(df_nhanes$yrs_in_US == 1, "< 1yr",
                                    ifelse(df_nhanes$yrs_in_US == 2, "01 to 04yr",   
                                    ifelse(df_nhanes$yrs_in_US == 3, "05 to 09yr",
                                    ifelse(df_nhanes$yrs_in_US == 4, "10 to 14yr",
                                    ifelse(df_nhanes$yrs_in_US == 5, "15 to 19yr",
                                    ifelse(df_nhanes$yrs_in_US == 6, "20 to 29yr",
                                    ifelse(df_nhanes$yrs_in_US == 7, "30 to 39yr",
                                    ifelse(df_nhanes$yrs_in_US == 8, "40 to 49yr",
                                    ifelse(df_nhanes$yrs_in_US == 9, "50yrs+",
                                    ifelse(df_nhanes$yrs_in_US == 77, "Refused",
                                    ifelse(df_nhanes$yrs_in_US == 99, "Don't know",
                                      NA))))))))))))
df_nhanes$educ_children <- factor( ifelse(df_nhanes$educ_children == 0, "Never attended / kindergarten only",
                                    ifelse(df_nhanes$educ_children == 1, "1st grade",     
                                    ifelse(df_nhanes$educ_children == 2, "2nd grade",
                                    ifelse(df_nhanes$educ_children == 3, "3rd grade",
                                    ifelse(df_nhanes$educ_children == 4, "4th grade",
                                    ifelse(df_nhanes$educ_children == 5, "5th grade",
                                    ifelse(df_nhanes$educ_children == 6, "6th grade",
                                    ifelse(df_nhanes$educ_children == 7, "7th grade",
                                    ifelse(df_nhanes$educ_children == 8, "8th grade",
                                    ifelse(df_nhanes$educ_children == 9, "9th grade",
                                    ifelse(df_nhanes$educ_children == 10, "10th grade",
                                    ifelse(df_nhanes$educ_children == 11, "11th grade",
                                    ifelse(df_nhanes$educ_children == 12, "12th grade, no diploma",
                                    ifelse(df_nhanes$educ_children == 13, "High school graduate",
                                    ifelse(df_nhanes$educ_children == 14, "GED or equivalent",
                                    ifelse(df_nhanes$educ_children == 15, "More than high school",
                                    ifelse(df_nhanes$educ_children == 55, "Less than 5th grade",
                                    ifelse(df_nhanes$educ_children == 66, "Less than 9th grade",
                                    ifelse(df_nhanes$educ_children == 77, "Refused",
                                    ifelse(df_nhanes$educ_children == 99, "Don't know",
                                        NA)))))))))))))))))))))
df_nhanes$educ_adults <- factor( ifelse(df_nhanes$educ_adults == 1, "Less than 9th grade",
                                    ifelse(df_nhanes$educ_adults == 2, "9-11th grade (Includes 12th grade with no diploma)",
                                    ifelse(df_nhanes$educ_adults == 3, "High school graduate/GED or equivalent",
                                    ifelse(df_nhanes$educ_adults == 4, "Some college or AA degree",
                                    ifelse(df_nhanes$educ_adults == 5, "College graduate or above",
                                    ifelse(df_nhanes$educ_adults == 7, "Refused",
                                    ifelse(df_nhanes$educ_adults == 9, "Don't know",
                                        NA))))))))
df_nhanes$marital_status <- factor( ifelse(df_nhanes$marital_status == 1, "Married",
                                    ifelse(df_nhanes$marital_status == 2, "Widowed",
                                    ifelse(df_nhanes$marital_status == 3, "Divorced",
                                    ifelse(df_nhanes$marital_status == 4, "Separated",
                                    ifelse(df_nhanes$marital_status == 5, "Never married",
                                    ifelse(df_nhanes$marital_status == 6, "Living with partner",
                                    ifelse(df_nhanes$marital_status == 77, "Refused",
                                    ifelse(df_nhanes$marital_status == 99, "Don't know",
                                        NA)))))))))
df_nhanes$preg_status_exam <- factor( ifelse(df_nhanes$preg_status_exam == 1, "Yes, lab or self-reported",
                                        ifelse(df_nhanes$preg_status_exam == 2, "Not Pregnant",
                                        ifelse(df_nhanes$preg_status_exam == 3, "Cannot ascertain",
                                        NA))))
df_nhanes$language <- factor( ifelse(df_nhanes$language == 1, "English",
                                     ifelse(df_nhanes$language == 2, "Spanish",   
                                     NA)))
df_nhanes$proxy <- factor( ifelse(df_nhanes$proxy == 1, "Yes",
                                  ifelse(df_nhanes$proxy == 2, "No",
                                    NA)))
df_nhanes$interpreter <- factor( ifelse(df_nhanes$interpreter == 1, "Yes",
                                    ifelse(df_nhanes$interpreter == 2, "No",  
                                    NA)))
df_nhanes$lang_family <- factor( ifelse(df_nhanes$lang_family == 1, "English",
                                    ifelse(df_nhanes$lang_family == 2, "Spanish",  
                                    NA)))
df_nhanes$proxy_family <- factor( ifelse(df_nhanes$proxy_family == 1, "Yes",
                                    ifelse(df_nhanes$proxy_family == 2, "No",     
                                    NA)))
df_nhanes$interp_family <- factor( ifelse(df_nhanes$interp_family == 1, "Yes",
                                    ifelse(df_nhanes$interp_family == 2, "No",      
                                    NA)))
df_nhanes$lang_mec <- factor( ifelse(df_nhanes$lang_mec == 1, "English",
                                ifelse(df_nhanes$lang_mec == 2, "Spanish",   
                                     NA)))
df_nhanes$proxy_mec <- factor( ifelse(df_nhanes$proxy_mec == 1, "Yes",
                                ifelse(df_nhanes$proxy_mec == 2, "No",  
                                      NA)))
df_nhanes$interp_mec <- factor( ifelse(df_nhanes$interp_mec == 1, "Yes",
                                ifelse(df_nhanes$interp_mec == 2, "No",  
                                       NA)))
df_nhanes$lang_acasi <- factor( ifelse(df_nhanes$lang_acasi == 1, "English",
                                ifelse(df_nhanes$lang_acasi == 2, "Spanish",
                                ifelse(df_nhanes$lang_acasi == 3, "Asian languages",
                                       NA))))
df_nhanes$num_house <- factor( ifelse(df_nhanes$num_house == 1, "1",
                                ifelse(df_nhanes$num_house == 2, "2",
                                ifelse(df_nhanes$num_house == 3, "3",
                                ifelse(df_nhanes$num_house == 4, "4",
                                ifelse(df_nhanes$num_house == 5, "5",
                                ifelse(df_nhanes$num_house == 6, "6",
                                ifelse(df_nhanes$num_house == 7, "7 or more",
                                      NA))))))))
df_nhanes$num_family <- factor( ifelse(df_nhanes$num_family == 1, "1",
                                ifelse(df_nhanes$num_family == 2, "2",
                                ifelse(df_nhanes$num_family == 3, "3",
                                ifelse(df_nhanes$num_family == 4, "4",
                                ifelse(df_nhanes$num_family == 5, "5",
                                ifelse(df_nhanes$num_family == 6, "6",
                                ifelse(df_nhanes$num_family == 7, "7 or more",
                                       NA))))))))
df_nhanes$num_under5 <- factor( ifelse(df_nhanes$num_under5 == 0, "0",
                                ifelse(df_nhanes$num_under5 == 1, "1",
                                ifelse(df_nhanes$num_under5 == 2, "2",
                                ifelse(df_nhanes$num_under5 == 3, "3 or more",
                                       NA)))))
df_nhanes$num_6to17 <- factor( ifelse(df_nhanes$num_6to17 == 0, "0",
                                ifelse(df_nhanes$num_6to17 == 1, "1",
                                ifelse(df_nhanes$num_6to17 == 2, "2",
                                ifelse(df_nhanes$num_6to17 == 3, "3",
                                ifelse(df_nhanes$num_6to17 == 4, "4 or more",
                                       NA))))))
df_nhanes$num_over60 <- factor( ifelse(df_nhanes$num_over60 == 0, "0",
                                ifelse(df_nhanes$num_over60 == 1, "1",
                                ifelse(df_nhanes$num_over60 == 2, "2",
                                ifelse(df_nhanes$num_over60 == 3, "3 or more",
                                       NA)))))
df_nhanes$HH_gender <- factor( ifelse(df_nhanes$HH_gender == 1, "Male",
                                ifelse(df_nhanes$HH_gender == 2, "Female", 
                                      NA)))
## df_nhanes$HH_age - skip
df_nhanes$HH_birth_country <- factor( ifelse(df_nhanes$HH_birth_country == 1, "US_50states",
                                        ifelse(df_nhanes$HH_birth_country == 2, "Others",
                                        ifelse(df_nhanes$HH_birth_country == 77, "Refused",
                                        ifelse(df_nhanes$HH_birth_country == 99, "Don't know",
                                               NA)))))
df_nhanes$HH_educ <- factor( ifelse(df_nhanes$HH_educ == 1, "Less than 9th grade",
                        ifelse(df_nhanes$HH_educ == 2, "9-11th grade (Includes 12th grade with no diploma)",
                        ifelse(df_nhanes$HH_educ == 3, "High school graduate/GED or equivalent",
                        ifelse(df_nhanes$HH_educ == 4, "Some college or AA degree",
                        ifelse(df_nhanes$HH_educ == 5, "College graduate or above",
                        ifelse(df_nhanes$HH_educ == 7, "Refused",
                        ifelse(df_nhanes$HH_educ == 9, "Don't know",
                               NA))))))))
#df_nhanes$HH_marriage - skip
#df_nhanes$HH_spouse_educ - skip
df_nhanes$house_income <- factor( ifelse(df_nhanes$house_income == 1, "$0 to $4,999",
                                    ifelse(df_nhanes$house_income == 2, "$5,000 to $9,999",
                                    ifelse(df_nhanes$house_income == 3, "$10,000 to $14,999",
                                    ifelse(df_nhanes$house_income == 4, "$15,000 to $19,999",
                                    ifelse(df_nhanes$house_income == 5, "$20,000 to $24,999",
                                    ifelse(df_nhanes$house_income == 6, "$25,000 to $34,999",
                                    ifelse(df_nhanes$house_income == 7, "$35,000 to $44,999",
                                    ifelse(df_nhanes$house_income == 8, "$45,000 to $54,999",
                                    ifelse(df_nhanes$house_income == 9, "$55,000 to $64,999",
                                    ifelse(df_nhanes$house_income == 10, "$65,000 to $74,999",
                                    ifelse(df_nhanes$house_income == 12, "$20,000 and Over",
                                    ifelse(df_nhanes$house_income == 13, "Under $20,000",
                                    ifelse(df_nhanes$house_income == 14, "$75,000 to $99,999",
                                    ifelse(df_nhanes$house_income == 15, "$100,000 and Over",
                                    ifelse(df_nhanes$house_income == 77, "Refused",
                                    ifelse(df_nhanes$house_income == 99, "Don't know",
                                         NA)))))))))))))))))
df_nhanes$family_income <- factor( ifelse(df_nhanes$family_income == 1, "$0 to $4,999",
                                    ifelse(df_nhanes$family_income == 2, "$5,000 to $9,999",
                                    ifelse(df_nhanes$family_income == 3, "$10,000 to $14,999",
                                    ifelse(df_nhanes$family_income == 4, "$15,000 to $19,999",
                                    ifelse(df_nhanes$family_income == 5, "$20,000 to $24,999",
                                    ifelse(df_nhanes$family_income == 6, "$25,000 to $34,999",
                                    ifelse(df_nhanes$family_income == 7, "$35,000 to $44,999",
                                    ifelse(df_nhanes$family_income == 8, "$45,000 to $54,999",
                                    ifelse(df_nhanes$family_income == 9, "$55,000 to $64,999",
                                    ifelse(df_nhanes$family_income == 10, "$65,000 to $74,999",
                                    ifelse(df_nhanes$family_income == 12, "$20,000 and Over",
                                    ifelse(df_nhanes$family_income == 13, "Under $20,000",
                                    ifelse(df_nhanes$family_income == 14, "$75,000 to $99,999",
                                    ifelse(df_nhanes$family_income == 15, "$100,000 and Over",
                                    ifelse(df_nhanes$family_income == 77, "Refused",
                                    ifelse(df_nhanes$family_income == 99, "Don't know",
                                        NA)))))))))))))))))


df_nhanes <- df_nhanes %>%
    dplyr::rename("body_meas_status" = "BMDSTATS",
           "weight_kg" = "BMXWT",
           "weight_comment" = "BMIWT",
           "recum_length_cm" = "BMXRECUM",
           "recum_length_comment" = "BMIRECUM",
           "head_circum" = "BMXHEAD",
           "head_circum_comment" = "BMIHEAD",
           "height" = "BMXHT",
           "height_comment" = "BMIHT",
           "bmi" = "BMXBMI",
           "bmi_cat" = "BMDBMIC",
           "upp_leg_length" = "BMXLEG",
           "upp_leg_length_comment" = "BMILEG",
           "upp_arm_length" = "BMXARML",
           "upp_arm_length_comment" = "BMIARML",
           "arm_circ" = "BMXARMC",
           "arm_circ_comment" = "BMIARMC",
           "waist_circ_cm" = "BMXWAIST",
           "waist_circ_comment" = "BMIWAIST",
           "sagg_diam_1" = "BMXSAD1",
           "sagg_diam_2" = "BMXSAD2",
           "sagg_diam_3" = "BMXSAD3",
           "sagg_diam_4" = "BMXSAD4",
           "sagg_diam_avg" = "BMDAVSAD",
           "sagg_diam_comment" = "BMDSADCM"
           )

##
df_nhanes$weight_comment <- factor( ifelse(df_nhanes$weight_comment == 1, "Could not obtain",
                                    ifelse(df_nhanes$weight_comment == 3, "Clothing",
                                    ifelse(df_nhanes$weight_comment == 4, "Medical appliance",
                                           NA))))
df_nhanes$recum_length_comment <- factor( ifelse(df_nhanes$recum_length_comment == 1, "Could not obtain",
                                    ifelse(df_nhanes$recum_length_comment == 2, "2",
                                    ifelse(df_nhanes$recum_length_comment == 3, "Not straight",
                                        NA))))
df_nhanes$height_comment <- factor( ifelse(df_nhanes$height_comment == 1, "Could not obtain",
                                    ifelse(df_nhanes$height_comment == 2, "2",
                                    ifelse(df_nhanes$height_comment == 3, "Not straight",      
                                           NA))))
df_nhanes$bmi_cat <- factor( ifelse(df_nhanes$bmi_cat == 1, "Underweight",
                            ifelse(df_nhanes$bmi_cat == 2, "Normal weight",
                            ifelse(df_nhanes$bmi_cat == 3, "Overweight",
                            ifelse(df_nhanes$bmi_cat == 4, "Obese",
                                   NA)))))
df_nhanes$upp_leg_length_comment <- factor( ifelse(df_nhanes$upp_leg_length_comment == 1, "Could not obtain",
                                                   NA))
df_nhanes$upp_arm_length_comment <- factor( ifelse(df_nhanes$upp_arm_length_comment == 1, "Could not obtain",
                                                   NA))
df_nhanes$arm_circ_comment <- factor( ifelse(df_nhanes$arm_circ_comment ==1, "Could not obtain",
                                             NA))
df_nhanes$waist_circ_comment <- factor( ifelse(df_nhanes$waist_circ_comment == 1, "Could not obtain",
                                               NA))                                           
df_nhanes$sagg_diam_comment <- factor( ifelse(df_nhanes$sagg_diam_comment == 1, "Could not obtain",
                        ifelse(df_nhanes$sagg_diam_comment == 2, "SP unable to comply with exam instruction",
                        ifelse(df_nhanes$sagg_diam_comment == 3, "SP discomfort",
                        ifelse(df_nhanes$sagg_diam_comment == 4, "Use of positioning cushion",
                        ifelse(df_nhanes$sagg_diam_comment == 5, "Other",
                                              NA))))))

df_nhanes <- df_nhanes %>%
    dplyr::rename("bl_pb"= "LBXBPB",
           "bl_pb_si" = "LBDBPBSI",
           "bl_pb_comment" = "LBDBPBLC",
           "bl_cd" = "LBXBCD",
           "bl_cd_si" = "LBDBCDSI",
           "bl_cd_comment" = "LBDBCDLC",
           "bl_total_hg" = "LBXTHG",
           "bl_total_hg_si" = "LBDTHGSI",
           "bl_total_hg_comment" = "LBDTHGLC",
           "bl_sel"= "LBXBSE",
           "bl_sel_si"= "LBDBSESI",
           "bl_se_comment" = "LBDBSELC",
           "bl_mang" = "LBXBMN",
           "bl_mang_si" = "LBDBMNSI",
           "bl_mang_comment" = "LBDBMNLC",
           "bl_I_hg" = "LBXIHG",
           "bl_I_hg_si" = "LBDIHGSI",
           "bl_I_hg_comment" = "LBDIHGLC",
           "bl_Et_hg" = "LBXBGE",
           "bl_Et_hg_comment" = "LBDBGELC",
           "bl_Me_hg" = "LBXBGM",
           "bl_Me_hg_comment" = "LBDBGMLC",
           "ur_iod" = "URXUIO",
           "ur_iod_comment" = "URDUIOLC",
           "ur_hg"= "URXUHG",
           "ur_hg_comment"= "URDUHGLC",
           "ur_bar" = "URXUBA",
           "ur_bar_comment" = "URDUBALC",
           "ur_cd" = "URXUCD",
           "ur_cd_comment" = "URDUCDLC",
           "ur_co" = "URXUCO",
           "ur_co_comment" = "URDUCOLC",
           "ur_ce" = "URXUCS",
           "ur_ce_comment" = "URDUCSLC",
           "ur_mo" = "URXUMO",
           "ur_mo_comment" = "URDUMOLC",
           "ur_mang" = "URXUMN",
           "ur_mang_comment" = "URDUMNLC",
           "ur_pb" = "URXUPB",
           "ur_pb_comment" = "URDUPBLC",
           "ur_ant"= "URXUSB",
           "ur_ant_comment"= "URDUSBLC",
           "ur_sn" = "URXUSN",
           "ur_sn_comment" = "URDUSNLC",
           "ur_sr" = "URXUSR",
           "ur_sr_comment" = "URDUSRLC",
           "ur_tl" = "URXUTL",
           "ur_tl_comment" = "URDUTLLC",
           "ur_tung" = "URXUTU",
           "ur_tung_comment" = "URDUTULC",
           "ur_uran" = "URXUUR",
           "ur_uran_comment" = "URDUURLC"
           )

df_nhanes <- df_nhanes %>%
    dplyr::rename("bl_cr_ug_L" = "LBXBCR",
                  "bl_cr_si" = "LBDBCRSI",
                  "bl_cr_comment" = "LBDBCRLC",
                  "bl_co_ug_L" = "LBXBCO", 
                  "bl_co_si" = "LBDBCOSI",
                  "bl_co_comment" = "LBDBCOLC")

     
df_nhanes <- df_nhanes %>%
    dplyr::rename("ur_alb"= "URXUMA",
           "ur_alb_comment" = "URDUMALC",
           "ur_alb_mgL" = "URXUMS",
           "ur_cr_mgdL" = "URXUCR",
           "ur_cr_comment" = "URDUCRLC",
           "ur_cr_umolL" = "URXCRS",
           "ur_alb_cr" = "URDACT")

df_nhanes$ur_alb_comment <- factor( ifelse(df_nhanes$ur_alb_comment == 0, "Valid",
                                               ifelse(df_nhanes$ur_alb_comment == 1, "Below_Lod",
                                                      NA)))
df_nhanes$ur_cr_comment <- factor( ifelse(df_nhanes$ur_cr_comment == 0, "Valid",
                                              ifelse(df_nhanes$ur_cr_comment == 1, "Below_Lod",
                                                     NA)))
           
df_nhanes <- df_nhanes %>%
    dplyr::rename("serum_alb_gdL" = "LBXSAL",
           "serum_alb_gL_si" = "LBDSALSI",
           "alkphos_si" = "LBXSAPSI",
           "AST_si" = "LBXSASSI",
           "ALT_si"= "LBXSATSI",
           "BUN_mgdL" = "LBXSBU",
           "BUN_mmolL_si" = "LBDSBUSI",
           "bicarb_mmolL_si" = "LBXSC3SI",
           "bl_tot_ca_mg_dL" = "LBXSCA",
           "bl_tot_ca_mmolL_si" = "LBDSCASI",
#           "serum_chol_mgdL" = "LBXSCH",
#           "serum_chol_mmolL_si" = "LBDSCHSI",
           "serum_tot_chol_mgdl" = "LBXTC",  # these chol values recommended as preferred on nhanes website
           "serum_tot_chol_si" = "LBDTCSI",
           "CPK" = "LBXSCK",
           "cl_mmolL_si" = "LBXSCLSI",
           "serum_creat_mg_dL" = "LBXSCR",
           "serum_creat_umolL_si" = "LBDSCRSI",
           "globulin_gdL" = "LBXSGB",
           "globulin_gL_si" = "LBDSGBSI",
           "bl_glucose_mgdL" = "LBXSGL",
           "bl_glucose_mmolL_si" = "LBDSGLSI",
           "GGT_si" = "LBXSGTSI",
           "serum_iron_ugdL" = "LBXSIR",
           "serum_iron_umolL" = "LBDSIRSI",
           "serum_potassion_si" = "LBXSKSI",
           "LDH_si" = "LBXSLDSI",
           "serum_sodium_si" = "LBXSNASI",
           "serum_osmol_si" = "LBXSOSSI",
           "serum_phos_mgdL" = "LBXSPH",
           "serum_phos_mmolL_si" = "LBDSPHSI",
           "serum_bili_mgdL" = "LBXSTB",
           "serum_bili_umolL" = "LBDSTBSI",
           "serum_prot_gdL" = "LBXSTP",
           "serum_prot_gL_si" = "LBDSTPSI",
           "serum_triglyc_mgdL" = "LBXSTR",
           "serum_triglyc_mmolL_si" = "LBDSTRSI",
           "bl_uric_acid_mgdL" = "LBXSUA",
           "bl_uric_acid_umolL_si" = "LBDSUASI"
    )


df_nhanes <- df_nhanes %>%
    dplyr::rename("insulin_uUmL" = "LBXIN",
           "insulin_pmolL_si" = "LBDINSI",
           "insulin_comment" = "LBDINLC",
           "foodfast_hrs" = "PHAFSTHR",
           "foodfast_mins" = "PHAFSTMN",
           "fast_glucose_mgdL" = "LBXGLU",
           "fast_glucose_mmolL_si" = "LBDGLUSI",
           "glycohemo_pct" = "LBXGH"
           )
df_nhanes <- df_nhanes %>%
    dplyr::rename("ogtt2hr_mg_dl" = "LBXGLT",
                  "ogtt2hr_mmol/L_si" = "LBDGLTSI")  

df_nhanes <- df_nhanes %>%
    dplyr::rename("pl_fluor_umolL" = "LBDPFL")


df_nhanes <- df_nhanes %>%
    dplyr::rename("wcc" = "LBXWBCSI",
           "lymph_pct" = "LBXLYPCT",
           "monocyte_pct" = "LBXMOPCT",
           "neut_pct" = "LBXNEPCT",
           "eosin_pct" = "LBXEOPCT",
           "basophil_pct" = "LBXBAPCT",
           "lymph_cells_uL" = "LBDLYMNO",
           "mono_cells_uL" = "LBDMONO",
           "neut_cells_uL" = "LBDNENO",
           "eosin_cells_uL" = "LBDEONO",
           "basophil_cells_uL" = "LBDBANO",
           "rcc" = "LBXRBCSI",
           "Hb_gdL" = "LBXHGB",
#           "Hct_pct" = "LBXHCT",
           "mcv_fL" = "LBXMCVSI",
           "mch_gdL" = "LBXMCHSI",
           "mch_pg" = "LBXMCH",
           "rdw_pct" = "LBXRDW",
           "platelets" = "LBXPLTSI",
           "mpv_fL" = "LBXMPSI",
           "folate_ng_ml" = "LBDRFO",
           "folate_ng_ml_si" = "LBDRFOSI")

df_nhanes <- df_nhanes %>%
    dplyr::rename("serum_cotin_ng_ml" = "LBXCOT",
                  "serum_cotin_ng_ml_comment" = "LBDCOTLC",
                  #"serum_hydrocot_ng_ml" = "LBXHCT",
                  "serum_hydrocot_ng_ml_comment" = "LBDHCTLC",
                  "serum_testosterone_ng_dL" = "LBXTST",
                  "serum_testosterone_ng_dL_comment" = "LBDTSTLC",
                  "serum_estradiol_pg_ml" = "LBXEST",
                  "serum_estradiol_pg_ml_comment" = "LBDESTLC",
                  "serum_SHBG_nmol_L" = "LBXSHBG",
                  "serum_SHBG_nmol_L_comment" = "LBDSHGLC",
    )


## Save data for later analysis
saveRDS(df_nhanes, "Data/Nhanes_cleaned.RDS")


# Save new names to old names as look up table
var_lookup <- data.frame( rbind(
    c("Hct_pct", "LBXHCT.x"),
       c("serum_hydrocot_ng_ml", "LBXHCT.y"),
       c("WT_full2yr", "WTINT2YR"),
       c("WT_mec2yr", "WTMEC2YR"),
       c("WT_subsampleA_2yr", "WTSA2YR"),
       c("WT_fasting_subsamp_mec2yr", "WTSAF2YR"),
       c("WT_blood_metal", "WTSH2YR"),
       c("WT_oggt_subsamp_mec2yr", "WTSOG2YR"),
       c("exam_status", "RIDSTATR"),
       c("gender", "RIAGENDR"),
       c("age_screen_yrs", "RIDAGEYR"),
       c("age_screen_mnths", "RIDAGEMN"),
       c("race1", "RIDRETH1"),
       c("race3", "RIDRETH3"),
       c("exam_month", "RIDEXMON"),
       c("age_exam_mnths", "RIDEXAGM"),
       c("military", "DMQMILIZ"),
       c("served_abroad", "DMQADFC"),
       c("country_birth", "DMDBORN4"),
       c("citizenship", "DMDCITZN"),
       c("yrs_in_US", "DMDYRSUS"),
       c("educ_children", "DMDEDUC3"),
       c("educ_adults", "DMDEDUC2"),
       c("marital_status", "DMDMARTL"),
       c("preg_status_exam", "RIDEXPRG"),
       c("language", "SIALANG"),
       c("proxy", "SIAPROXY"),
       c("interpreter", "SIAINTRP"),
       c("lang_family", "FIALANG"),
       c("proxy_family", "FIAPROXY"),
       c("interp_family", "FIAINTRP"),
       c("lang_mec", "MIALANG"),
       c("proxy_mec", "MIAPROXY"),
       c("interp_mec", "MIAINTRP"),
       c("lang_acasi", "AIALANGA"),
       c("num_house", "DMDHHSIZ"),
       c("num_family", "DMDFMSIZ"),
       c("num_under5", "DMDHHSZA"),
       c("num_6to17", "DMDHHSZB"),
       c("num_over60", "DMDHHSZE"),
       c("HH_gender", "DMDHRGND"),
       c("HH_age", "DMDHRAGE"),
       c("HH_birth_country", "DMDHRBR4"),
       c("HH_educ", "DMDHREDU"),
       c("HH_marriage", "DMDHRMAR"),
       c("HH_spouse_educ", "DMDHSEDU"),
       c("house_income", "INDHHIN2"),
       c("family_income", "INDFMIN2"),
       c("income_2poverty", "INDFMPIR"),
       c("body_meas_status", "BMDSTATS"),
       c("weight_kg", "BMXWT"),
       c("weight_comment", "BMIWT"),
       c("recum_length_cm", "BMXRECUM"),
       c("recum_length_comment", "BMIRECUM"),
       c("head_circum", "BMXHEAD"),
       c("head_circum_comment", "BMIHEAD"),
       c("height", "BMXHT"),
       c("height_comment", "BMIHT"),
       c("bmi", "BMXBMI"),
       c("bmi_cat", "BMDBMIC"),
       c("upp_leg_length", "BMXLEG"),
       c("upp_leg_length_comment", "BMILEG"),
       c("upp_arm_length", "BMXARML"),
       c("upp_arm_length_comment", "BMIARML"),
       c("arm_circ", "BMXARMC"),
       c("arm_circ_comment", "BMIARMC"),
       c("waist_circ_cm", "BMXWAIST"),
       c("waist_circ_comment", "BMIWAIST"),
       c("sagg_diam_1", "BMXSAD1"),
       c("sagg_diam_2", "BMXSAD2"),
       c("sagg_diam_3", "BMXSAD3"),
       c("sagg_diam_4", "BMXSAD4"),
       c("sagg_diam_avg", "BMDAVSAD"),
       c("sagg_diam_comment", "BMDSADCM"),
       c("bl_pb", "LBXBPB"),
       c("bl_pb_si", "LBDBPBSI"),
       c("bl_pb_comment", "LBDBPBLC"),
       c("bl_cd", "LBXBCD"),
       c("bl_cd_si", "LBDBCDSI"),
       c("bl_cd_comment", "LBDBCDLC"),
       c("bl_total_hg", "LBXTHG"),
       c("bl_total_hg_si", "LBDTHGSI"),
       c("bl_total_hg_comment", "LBDTHGLC"),
       c("bl_sel", "LBXBSE"),
       c("bl_sel_si", "LBDBSESI"),
       c("bl_se_comment", "LBDBSELC"),
       c("bl_mang", "LBXBMN"),
       c("bl_mang_si", "LBDBMNSI"),
       c("bl_mang_comment", "LBDBMNLC"),
       c("bl_I_hg", "LBXIHG"),
       c("bl_I_hg_si", "LBDIHGSI"),
       c("bl_I_hg_comment", "LBDIHGLC"),
       c("bl_Et_hg", "LBXBGE"),
       c("bl_Et_hg_comment", "LBDBGELC"),
       c("bl_Me_hg", "LBXBGM"),
       c("bl_Me_hg_comment", "LBDBGMLC"),
       c("ur_iod", "URXUIO"),
       c("ur_iod_comment", "URDUIOLC"),
       c("ur_hg", "URXUHG"),
       c("ur_hg_comment", "URDUHGLC"),
       c("ur_bar", "URXUBA"),
       c("ur_bar_comment", "URDUBALC"),
       c("ur_cd", "URXUCD"),
       c("ur_cd_comment", "URDUCDLC"),
       c("ur_co", "URXUCO"),
       c("ur_co_comment", "URDUCOLC"),
       c("ur_ce", "URXUCS"),
       c("ur_ce_comment", "URDUCSLC"),
       c("ur_mo", "URXUMO"),
       c("ur_mo_comment", "URDUMOLC"),
       c("ur_mang", "URXUMN"),
       c("ur_mang_comment", "URDUMNLC"),
       c("ur_pb", "URXUPB"),
       c("ur_pb_comment", "URDUPBLC"),
       c("ur_ant", "URXUSB"),
       c("ur_ant_comment", "URDUSBLC"),
       c("ur_sn", "URXUSN"),
       c("ur_sn_comment", "URDUSNLC"),
       c("ur_sr", "URXUSR"),
       c("ur_sr_comment", "URDUSRLC"),
       c("ur_tl", "URXUTL"),
       c("ur_tl_comment", "URDUTLLC"),
       c("ur_tung", "URXUTU"),
       c("ur_tung_comment", "URDUTULC"),
       c("ur_uran", "URXUUR"),
       c("ur_uran_comment", "URDUURLC"),
       c("bl_cr_ug_L", "LBXBCR"),
       c("bl_cr_si", "LBDBCRSI"),
       c("bl_cr_comment", "LBDBCRLC"),
       c("bl_co_ug_L", "LBXBCO"),
       c("bl_co_si", "LBDBCOSI"),
       c("bl_co_comment", "LBDBCOLC"),
       c("ur_alb", "URXUMA"),
       c("ur_alb_comment", "URDUMALC"),
       c("ur_alb_mgL", "URXUMS"),
       c("ur_cr_mgdL", "URXUCR"),
       c("ur_cr_comment", "URDUCRLC"),
       c("ur_cr_umolL", "URXCRS"),
       c("ur_alb_cr", "URDACT"),
       c("serum_alb_gdL", "LBXSAL"),
       c("serum_alb_gL_si", "LBDSALSI"),
       c("alkphos_si", "LBXSAPSI"),
       c("AST_si", "LBXSASSI"),
       c("ALT_si", "LBXSATSI"),
       c("BUN_mgdL", "LBXSBU"),
       c("BUN_mmolL_si", "LBDSBUSI"),
       c("bicarb_mmolL_si", "LBXSC3SI"),
       c("bl_tot_ca_mg_dL", "LBXSCA"),
       c("bl_tot_ca_mmolL_si", "LBDSCASI"),
       c("serum_tot_chol_mgdl", "LBXTC"),
       c("serum_tot_chol_si", "LBDTCSI"),
       c("CPK", "LBXSCK"),
       c("cl_mmolL_si", "LBXSCLSI"),
       c("serum_creat_mg_dL", "LBXSCR"),
       c("serum_creat_umolL_si", "LBDSCRSI"),
       c("globulin_gdL", "LBXSGB"),
       c("globulin_gL_si", "LBDSGBSI"),
       c("bl_glucose_mgdL", "LBXSGL"),
       c("bl_glucose_mmolL_si", "LBDSGLSI"),
       c("GGT_si", "LBXSGTSI"),
       c("serum_iron_ugdL", "LBXSIR"),
       c("serum_iron_umolL", "LBDSIRSI"),
       c("serum_potassion_si", "LBXSKSI"),
       c("LDH_si", "LBXSLDSI"),
       c("serum_sodium_si", "LBXSNASI"),
       c("serum_osmol_si", "LBXSOSSI"),
       c("serum_phos_mgdL", "LBXSPH"),
       c("serum_phos_mmolL_si", "LBDSPHSI"),
       c("serum_bili_mgdL", "LBXSTB"),
       c("serum_bili_umolL", "LBDSTBSI"),
       c("serum_prot_gdL", "LBXSTP"),
       c("serum_prot_gL_si", "LBDSTPSI"),
       c("serum_triglyc_mgdL", "LBXSTR"),
       c("serum_triglyc_mmolL_si", "LBDSTRSI"),
       c("bl_uric_acid_mgdL", "LBXSUA"),
       c("bl_uric_acid_umolL_si", "LBDSUASI"),
       c("insulin_uUmL", "LBXIN"),
       c("insulin_pmolL_si", "LBDINSI"),
       c("insulin_comment", "LBDINLC"),
       c("foodfast_hrs", "PHAFSTHR"),
       c("foodfast_mins", "PHAFSTMN"),
       c("fast_glucose_mgdL", "LBXGLU"),
       c("fast_glucose_mmolL_si", "LBDGLUSI"),
       c("glycohemo_pct", "LBXGH"),
       c("ogtt2hr_mg_dl", "LBXGLT"),
       c("ogtt2hr_mmol/L_si", "LBDGLTSI"),
       c("pl_fluor_umolL", "LBDPFL"),
       c("wcc", "LBXWBCSI"),
       c("lymph_pct", "LBXLYPCT"),
       c("monocyte_pct", "LBXMOPCT"),
       c("neut_pct", "LBXNEPCT"),
       c("eosin_pct", "LBXEOPCT"),
       c("basophil_pct", "LBXBAPCT"),
       c("lymph_cells_uL", "LBDLYMNO"),
       c("mono_cells_uL", "LBDMONO"),
       c("neut_cells_uL", "LBDNENO"),
       c("eosin_cells_uL", "LBDEONO"),
       c("basophil_cells_uL", "LBDBANO"),
       c("rcc", "LBXRBCSI"),
       c("Hb_gdL", "LBXHGB"),
       c("mcv_fL", "LBXMCVSI"),
       c("mch_gdL", "LBXMCHSI"),
       c("mch_pg", "LBXMCH"),
       c("rdw_pct", "LBXRDW"),
       c("platelets", "LBXPLTSI"),
       c("mpv_fL", "LBXMPSI"),
       c("folate_ng_ml", "LBDRFO"),
       c("folate_ng_ml_si", "LBDRFOSI"),
       c("serum_cotin_ng_ml", "LBXCOT"),
       c("serum_cotin_ng_ml_comment", "LBDCOTLC"),
       c("serum_hydrocot_ng_ml_comment", "LBDHCTLC"),
       c("serum_testosterone_ng_dL", "LBXTST"),
       c("serum_testosterone_ng_dL_comment", "LBDTSTLC"),
       c("serum_estradiol_pg_ml", "LBXEST"),
       c("serum_estradiol_pg_ml_comment", "LBDESTLC"),
       c("serum_SHBG_nmol_L", "LBXSHBG"),
       c("serum_SHBG_nmol_L_comment", "LBDSHGLC")))
names(var_lookup) <- c("newname", "nhanesname")

write.csv(var_lookup, "Data/var_lookuplist.csv", row.names = FALSE)


# Blood pressure data
systolics <- c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")
diastolic <- c("BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4")
df_bp$sysBP <- rowMeans(df_bp[, systolics], na.rm=TRUE)
df_bp$diaBP <- rowMeans(df_bp[, diastolic], na.rm=TRUE)
df_bp <- df_bp %>%
    rename("HR_60s" = "BPXPLS",
           "Pulse_regularity"= "BPXPULS")
# Pulse regularity is coded aa 1 = regular and 2 = irregular
# Will subtract 1 to make this work with bernoulli marginal description
df_bp$Pulse_regularity <- df_bp$Pulse_regularity -1


# some have zero mean diastolic - why ?
df_bp[df_bp$mean_diabp ==0 & !is.na(df_bp$mean_diabp), ]

#save
saveRDS(df_bp, "Data/bp_data.RDS")

