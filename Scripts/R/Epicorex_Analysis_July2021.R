# Title: An analysis of urine and blood toxic metals concentrations versus BMI and haematological and biochemical blood panels in the NHANES dataset
# author: "James PK Rooney"
# date: "July 2020"

#### load libraries needed for analysis
library(tidyverse)
library(patchwork)
library(rcorex)
library(igraph)
library(ggraph)
library(future.apply)
library(tidygraph)


# Load data
df1 <- readRDS("Data/Nhanes_cleaned.RDS")
var_lookup <- read.csv("Data/var_lookuplist.csv", stringsAsFactors = FALSE)

# Load the blood pressure data
dfbp <- readRDS("Data/bp_data.RDS")
# Merge bp data to df1
df1 <- data.frame(df1, dfbp[ match(df1$SEQN, dfbp$SEQN),])



# data prep - categorical variables to be included need to be dummy coded
df1$gender <- ifelse(df1$gender == "Female", 0, 1)


##### Prep metal variables
# For urinary measure metals it is normal practice to correct for dilution by dividing by creatinine
# ur_iod and ur_hg are in ng/ml
# "ur_bar", "ur_cd", "ur_co", "ur_ce", "ur_mo", "ur_mang", "ur_pb", are in ug/L
# "ur_ant", "ur_sn", "ur_sr", "ur_tl", "ur_tung", "ur_uran") are in ug/L
# Convert all to g/L
df1 <- df1 %>%
    mutate(ur_iod_gL = ur_iod *1e-9 / 1e-3,
           ur_hg_gL = ur_hg *1e-9 / 1e-3,
           ur_bar_gL = ur_bar *1e-6,
           ur_cd_gL = ur_cd *1e-6,
           ur_co_gL = ur_co *1e-6,
           ur_ce_gL = ur_ce *1e-6,
           ur_mo_gL = ur_mo *1e-6,
           ur_mang_gL = ur_mang *1e-6,
           ur_pb_gL = ur_pb *1e-6,
           ur_ant_gL = ur_ant *1e-6,
           ur_sn_gL = ur_sn *1e-6,
           ur_sr_gL = ur_sr *1e-6,
           ur_tl_gL = ur_tl *1e-6,
           ur_tung_gL = ur_tung *1e-6,
           ur_uran_gL = ur_uran *1e-6,
           ur_cr_gL = ur_cr_mgdL *1e-3/1e-1
    )
# Divide metals concs by cr conc
df1 <- df1 %>%
    mutate(ur_iod_uggcrea = ur_iod_gL *1e6 / ur_cr_gL,
           ur_hg_uggcrea = ur_hg_gL *1e6 / ur_cr_gL,
           ur_bar_uggcrea = ur_bar_gL *1e6 / ur_cr_gL,
           ur_cd_uggcrea = ur_cd_gL *1e6 / ur_cr_gL,
           ur_co_uggcrea = ur_co_gL *1e6 / ur_cr_gL,
           ur_ce_uggcrea = ur_ce_gL *1e6 / ur_cr_gL,
           ur_mo_uggcrea = ur_mo_gL *1e6 / ur_cr_gL,
           ur_mang_uggcrea = ur_mang_gL *1e6 / ur_cr_gL,
           ur_pb_uggcrea = ur_pb_gL *1e6 / ur_cr_gL,
           ur_ant_uggcrea = ur_ant_gL *1e6  / ur_cr_gL,
           ur_sn_uggcrea = ur_sn_gL *1e6 / ur_cr_gL,
           ur_sr_uggcrea = ur_sr_gL *1e6 / ur_cr_gL,
           ur_tl_uggcrea = ur_tl_gL *1e6 / ur_cr_gL,
           ur_tung_uggcrea = ur_tung_gL / ur_cr_gL,
           ur_uran_uggcrea = ur_uran_gL / ur_cr_gL)

# make a variable containing variables of interest
vars4corex <- c("gender", "age_screen_yrs", "bmi",
                "bl_pb", "bl_cd", "bl_total_hg", "bl_sel", "bl_mang", "bl_I_hg",
                "bl_Et_hg", "bl_Me_hg", "bl_cr_ug_L", "bl_co_ug_L", "ur_alb_cr",
                "serum_tot_chol_mgdl", "wcc", "lymph_pct", "monocyte_pct",
                "neut_pct", "eosin_pct","basophil_pct", "rcc", "Hb_gdL", "Hct_pct", 
                "mcv_fL", "rdw_pct", "platelets", "mpv_fL", 
                "glycohemo_pct", "insulin_uUmL", "foodfast_mins", "ogtt2hr_mg_dl",
                "serum_alb_gdL", "alkphos_si", "AST_si", "ALT_si", "BUN_mgdL",
                "bicarb_mmolL_si", "bl_tot_ca_mg_dL", "CPK", "cl_mmolL_si",
                "serum_creat_mg_dL", "globulin_gdL", "bl_glucose_mgdL", "GGT_si",
                "serum_iron_ugdL", "serum_potassion_si", "LDH_si", "serum_sodium_si",
                "serum_osmol_si", "serum_phos_mgdL", "serum_bili_mgdL", "serum_prot_gdL",
                "serum_triglyc_mgdL", "bl_uric_acid_mgdL", "ur_iod_uggcrea",
                "ur_hg_uggcrea", "ur_bar_uggcrea", "ur_cd_uggcrea", 
                "ur_co_uggcrea", "ur_ce_uggcrea", "ur_mo_uggcrea", "ur_mang_uggcrea", 
                "ur_pb_uggcrea", "ur_ant_uggcrea", "ur_sn_uggcrea", "ur_sr_uggcrea", 
                "ur_tl_uggcrea", "ur_tung_uggcrea", "ur_uran_uggcrea",
                "folate_ng_ml", "serum_testosterone_ng_dL",
                "serum_estradiol_pg_ml", "serum_SHBG_nmol_L", "serum_cotin_ng_ml",
                "serum_hydrocot_ng_ml",
                "HR_60s", 
                "Pulse_regularity", "sysBP", "diaBP")

# Make a df to categorise these vars by type
df_varcats <- data.frame(names = vars4corex,
                         cat = c("demographic", "demographic", "demographic", 
                                 "metal-blood", "metal-blood", "metal-blood",
                                 "metal-blood", "metal-blood", "metal-blood",
                                 "metal-blood", "metal-blood", "metal-blood",
                                 "metal-blood", "renal", "lipid",
                                 "haem", "haem", "haem", "haem", "haem",
                                 "haem", "haem", "haem", "haem", "haem",
                                 "haem", "haem", "haem", "lipid", "biochem", "biochem", "biochem",
                                 "biochem", "biochem", "biochem", "biochem",
                                 "biochem", "biochem", "biochem", "biochem", 
                                 "biochem", "biochem", "biochem", "lipid",
                                 "biochem", "biochem", "renal", "renal", "renal",
                                 "renal", "renal", "biochem", "biochem", "lipid",
                                 "biochem",
                                 "metal-urine", "metal-urine", "metal-urine", 
                                 "metal-urine", "metal-urine", "metal-urine", 
                                 "metal-urine", "metal-urine", "metal-urine", 
                                 "metal-urine", "metal-urine", "metal-urine", 
                                 "metal-urine", "metal-urine", "metal-urine", 
                                 "biochem", "hormone", "hormone",
                                 "hormone",
                                 "smoking", "smoking",
                                 "heart", "heart", "heart", "heart"))


# Make an indicator variable for adults vs children
df1$agecat <- ifelse(df1$age_screen_yrs < 20, "child", "adult")


### Remove those not selected for ANY sub-samples as they won't have included relevant measures
df2 <- df1[ !(df1$WT_blood_metal_comment =="Not selected" &
                  df1$WT_subsampleA_2yr_comment == "Not selected"), ]

# Make 2 separate dataframes for adults and children
df_kids <- df2[df2$agecat =="child", ]
df_adults <- df2[df2$agecat =="adult", ]

# Inspect data for missing values

mis_plot_kids <- df_kids[, c(vars4corex, "agecat")] %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
                      values = c('steelblue', 'tomato3'),
                      labels = c("Present", "Missing")) +
    #scale_x_discrete(limits = levels) +
    labs(x = "Variable",
         y = "Row Number", title = "Children") +
    coord_flip()
mis_plot_kids

mis_plot_adults <- df_adults[, c(vars4corex, "agecat")] %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
                      values = c('steelblue', 'tomato3'),
                      labels = c("Present", "Missing")) +
    #scale_x_discrete(limits = levels) +
    labs(x = "Variable",
         y = "Row Number", title = "Adults") +
    coord_flip()
mis_plot_adults

######
# For main analysis will use only adults as children have many missing values

### Check missing values by % missing per variable
vars_pct_missing <- df_adults[, c(vars4corex)] %>%
    mutate(id = row_number()) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>% 
    group_by(key) %>% 
    summarise(n_miss = sum(isna),
              pct_miss = sum(isna)*100 / n()) %>% 
    arrange(-pct_miss )
vars_pct_missing

# Exlcude variables missing more than 50%
vars_keep <- vars_pct_missing[ vars_pct_missing$pct_miss <=50, ]$key
vars4corex <- vars4corex[ vars4corex %in% vars_keep ]

# Make table describing variables
var_describe <- left_join(df_varcats, var_lookup, by = c("names" = "newname"))
var_describe <- var_describe %>% arrange(cat)

write.csv(var_describe, "Results/variable_description.csv", row.names = FALSE)


# Convert variables to log when appropriate (all vars except those excluded in follow line)
vars_log_convert <- vars4corex [ !vars4corex %in% c("gender", "age_screen_yrs", "bmi",
                                                    "HR_60s", "Pulse_regularity", "sysBP", 
                                                    "diaBP")]
df_adults[, vars_log_convert] <- log(df_adults[, vars_log_convert])

# replace infinities with NA
df_adults[df_adults == -Inf] <- NA



##### Run corex layer 1 across a selection of n_hidden and dim_hidden ####

set.seed(4)
# Parallelisation using future.apply package
plan(multisession, workers=10)
# set parameters limits for corex runs
hidmax <- round(length(vars4corex)/2) # max latent groups to half total variables
dimhidmax <- 3 

# Define marginal descriptions. Only gender and pulse regularity are bernoulli - all others are gaussian
marg_desc <- c("bernoulli", rep("gaussian", 73), "bernoulli", "gaussian", "gaussian")
# run the corex layer in parallel
res_layer1 <- future_lapply(1:hidmax, function(i){
    lapply(2:dimhidmax, function(j)
        epicorex(data = df_adults[, vars4corex], n_hidden = i, dim_hidden=j,
                 marginal_description = marg_desc, smooth_marginals = FALSE,
                 max_iter = 200, repeats = 25, eps=1e-5, logpx_method = "pycorex")
    )}, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
saveRDS(res_layer1, "Results/epi_corex_results_layer1.RDS")
# Optional code to save previously run corex layer
#res_layer1 <- readRDS("Results/epi_corex_results_layer1.RDS")


#### Extract TCS from layer 1 fits
new_names <- paste0("dim", 2:dimhidmax)

tcs_all <- res_layer1 %>%
    map_depth(2, "tcs") %>%
    map_depth(2, sum) %>%
    map(as.data.frame)

tcs_all <- lapply(1:length(tcs_all),
                  function(i)
                      tcs_all[[i]] %>% rename_at(vars(starts_with('X')), ~ new_names) )
tcs_all <- do.call("rbind", tcs_all)
tcs_all$nhid <- 1:hidmax

write.csv(tcs_all, "Results/epi_layer1_tcs_summary.csv")

# Check for convergence - recorded in 'state'
state_all <- res_layer1 %>%
    map_depth(2, "state") %>%
    map(as.data.frame)

state_all <- lapply(1:length(state_all),
                    function(i)
                        state_all[[i]] %>% rename_at(vars(starts_with('X')), ~ new_names) )

state_all <- do.call("rbind", state_all)
state_all$nhid <- 1:hidmax
state_long <- state_all %>% pivot_longer(cols = starts_with("dim"), values_to = "state")

# Make tcs long format
tcs_long <- tcs_all %>% pivot_longer(cols = starts_with("dim"), values_to = "tcs")
tcs_long$dim <- as.numeric( as.character( str_extract( tcs_long$name, "[[:digit:]]+") ) )
tcs_long$n_pars <- 2* as.numeric(tcs_long$nhid) *tcs_long$dim
tcs_long$lin <- tcs_long$tcs / log(tcs_long$dim)

# Join iterations to tcs data
tcs_long <- left_join(tcs_long, state_long)
tcs_long$nhid <- as.factor(tcs_long$nhid)
tcs_long$converged <- ifelse(tcs_long$state == "Converged", TRUE, FALSE)

# remove unconverged rows
tcs_long <- tcs_long %>% filter(converged == TRUE)

lin_tcs <- tcs_long %>% pivot_wider(id_cols = nhid, names_from = dim, values_from = lin)
write.csv(lin_tcs, "Results/epi_layer1_lintcs.csv")


# Plot TCS by parameters - exclude negative tcs from plot
g_tcs <- ggplot(tcs_long %>% filter(tcs >=0), aes(x=dim, y = tcs, group=nhid, col=nhid)) +
    geom_point() + geom_line() + expand_limits(y=0) +
    xlab("Dim_hidden") + ylab("TC") + ggtitle("Total Correlation") +
    theme_minimal()

# Plot
g_tcslogdim <- ggplot(tcs_long %>% filter(tcs >=0),
                      aes(x=dim, y = lin, group=nhid, col=nhid)) +
    geom_point() + geom_line() + expand_limits(y=0) +
    xlab("Dim_hidden") + ylab("TC / log(dim_hidden)") +
    ggtitle("Total Correlation per dimensions, k, of hidden variables") +
    theme_minimal()

tiff("Graphs/epi_TCS_by_nhid_dim.tiff", 800, 530)
    print(g_tcs + g_tcslogdim + plot_layout(guides = "collect"))
dev.off()


#### choose best fit - maximise the linearised tcs and check tcs by iteration
ord <- order( tcs_long$lin, decreasing = TRUE  )
dim(res_layer1[[ tcs_long[ ord[1] ,]$nhid ]][[ tcs_long[ ord[1] ,]$dim - 1]]$log_p_y)
plot( res_layer1[[ tcs_long[ ord[1] ,]$nhid ]][[ tcs_long[ ord[1] ,]$dim - 1]] )

# The model fit appears stable
# save plot to file
tiff("Graphs/model_fit.tiff", 800, 530)
    print(plot( res_layer1[[ tcs_long[ ord[1] ,]$nhid ]][[ tcs_long[ ord[1] ,]$dim - 1]] ))
dev.off()


# Assign preferred model to layer1
layer1 <- res_layer1[[ tcs_long[ ord[1] ,]$nhid ]][[ tcs_long[ ord[1] ,]$dim - 1]]

#### run layer2 corex ####
set.seed(4)
plan(multisession, workers=14)

res_layer2 <- future_lapply(1:dim(layer1$labels)[2], function(i){
    biocorex(data = layer1$labels, n_hidden = i, dim_hidden=2,
             marginal_description = "discrete", smooth_marginals = FALSE,
             max_iter = 200, repeats = 25, eps=1e-4, logpx_method = "pycorex")
}, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))
saveRDS(res_layer2, "Results/epi_corex_results_layer2.RDS")
#res_layer2 <- readRDS("Results/epi_corex_results_layer2.RDS")

l2_names <- paste0("dim", 2:dimhidmax)

tcs_all_l2 <- res_layer2 %>%
    map_depth(1, "tcs") %>%
    map_depth(1, sum) %>%
    unlist()
state_all_l2 <- res_layer2 %>%
    map_depth(1, "state") %>%
    #map_depth(1, sum) %>%
    unlist()
table(state_all_l2)
# all are true don't need to remove any

layer2 <- res_layer2[[ which.max(tcs_all_l2) ]]
plot(layer2) # fit looks stable



#### run layer3 corex ####
set.seed(456)
res_layer3 <- future_lapply(1:dim(layer2$labels)[2], function(i){
    biocorex(data = layer2$labels, n_hidden = i, dim_hidden=2,
             marginal_description = "discrete", smooth_marginals = FALSE,
             max_iter = 200, repeats = 25, eps=1e-4, logpx_method= "mean")
    }, future.seed = TRUE, future.scheduling = structure(TRUE, ordering = "random"))

l3_names <- paste0("dim", 2:dimhidmax)

tcs_all_l3 <- res_layer3 %>%
    map_depth(1, "tcs") %>%
    map_depth(1, sum) %>%
    unlist()
state_all_l3 <- res_layer3 %>%
    map_depth(1, "state") %>%
    #map_depth(1, sum) %>%
    unlist()
table(state_all_l3)
# all are true don't need to remove any

layer3 <- res_layer3[[ which.max(tcs_all_l3) ]]
plot(layer3)


sum(layer1$tcs)
sum(layer2$tcs)
sum(layer3$tcs)

g1 <- make_corex_tidygraph(list(layer1, layer2, layer3))

g1 <- activate(g1, nodes) %>% 
    left_join(df_varcats) %>% 
    mutate(col = factor(cat))

g1 <- g1 %>%
    activate(edges) %>%
    mutate(names = .N()$names[from]) %>% 
    left_join(df_varcats ) %>% 
    mutate(col_edge = factor(cat))


g1 <- ggraph(g1, layout = "fr") +
    geom_node_point(aes(size = node_size), show.legend = FALSE) +
    geom_edge_hive(aes(width = thickness), alpha = 0.75, show.legend = FALSE) +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_text(aes(label = names), repel = TRUE) +
    theme_graph()

tiff("Graphs/epi_Nhanes_network.tiff", 1000, 650)
    print(g1)
dev.off()


# Make a table of variable names by cluster for each layer
clust_table <- data.frame(varname = names(layer1$clusters),
                          layer1 = layer1$clusters + 1)
row.names(clust_table) <- NULL
temp <- data.frame(l1 = (1: length(unique(clust_table$layer1))),
                   TC_l1 = layer1$tcs,
                   layer2 = layer2$clusters)
clust_table$layer2 <- temp[match(clust_table$layer1, temp$l1),]$layer2 + 1
clust_table$TC_l1 <- layer1$tcs[ clust_table$layer1 ]
clust_table$TC_l2 <- layer2$tcs[ clust_table$layer2 ]

clust_table <- clust_table %>% left_join(var_describe, by = c("varname" = "names"))
clust_table$TC_l1[ duplicated(clust_table$TC_l1) ] <- NA
clust_table$TC_l2[ duplicated(clust_table$TC_l2) ] <- NA

clust_table <- clust_table %>% arrange(layer1, layer2)

clust_table <- clust_table %>%
    dplyr::select(varname, layer1, TC_l1, layer2, TC_l2, cat)

write.csv(clust_table, "Results/cluster_layer_summstats.csv")




