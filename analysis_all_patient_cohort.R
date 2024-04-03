pacman::p_load(haven, tidyverse, lubridate, gdata, survival, survminer, ggsci, gtsummary, WeightIt, cobalt, smd, RISCA, missForest)
set.cobalt.options(binary = "std")

# define ggplot2 theme...
theme_jikei <- function(base_size = 10, 
                        dark_text = "#1A242F") {
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = "#38414A", lineheight = 1.1),
          plot.title = element_text(color = "#38414A", size = rel(1.2)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = "#38414A", size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = "#38414A", size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

make_ppt <- function(data) {
  data_ft <- data %>% as_flex_table()
  print(data_ft, preview = "pptx")
}

D <- LATITUDE %>%
  dplyr::mutate(
    NODE = case_when(is.na(NODE) ~ "yes", TRUE ~ NODE), # trial handles NA as yes
    SURGERY = as.numeric(0),
    EBRT = case_when(EBRT == "yes" ~ 1, TRUE ~ 0),
    LOCAL = case_when((SURGERY == 0 & EBRT == 0) ~ "no", TRUE ~ "yes"),
    NUM_BONE = NUM_BONE %>% as.integer(),
    PSA = as.numeric(PSA),
    LDH = as.numeric(LDH),
    SCR = as.numeric(SCR),
    logPSA = log(PSA),
    logLDH = log(LDH),
    ARM = as.factor(ARM) %>% factor(levels = c("ADT", "AAP"))
  ) %>% 
  select(BMA_WITHIN, PPI_WITHIN, METFORMIN_WITHIN, STATIN_WITHIN, GC_WITHIN, 
         PFS, PROG, FU, DEATH, TTSRE, TTSRE_PRTB, TTSRE_CF, TTSRE_SCC, TTSRE_STB, SRE, SRE_PRTB, SRE_CF, SRE_SCC, SRE_STB, DM, 
         AGE_CAT, PS, logPSA, BPI, BPI_CAT, GG_CAT, BONE, NUM_BONE, NUM_BONE_CAT, VISCERAL, LUNG, LIVER, NODE, HGB, SCR, logLDH, ARM) %>%
  mutate_if(is.character, as.factor) %>% mutate_if(is.logical, as.factor) %>%
  as.data.frame()

# details of BMAs
theme_gtsummary_journal(journal = "jama")
tbl_bma <- LATITUDE %>%
  filter(BMA_WITHIN == 1 | DN_WITHIN == 1) %>%
  mutate(
    BMA = case_when(
      str_detect(BP_CMDECOD, "ALENDRONIC|ALENDRONATE") ~ "Alendronic acid",
      str_detect(BP_CMDECOD, "CLODRONIC|CLODRONATE") ~ "Clodronic acid",
      str_detect(BP_CMDECOD, "IBANDRONIC|IBANDRONATE") ~ "Ibandronic acid",
      str_detect(BP_CMDECOD, "FOSAVANCE") ~ "Alendronic acid",
      str_detect(BP_CMDECOD, "DISODIUM") ~ "Disodium incadronate",
      str_detect(BP_CMDECOD, "PAMIDRONIC|PAMIDRONATE") ~ "Pamidronic acid",
      str_detect(BP_CMDECOD, "RISEDRONIC|RISEDRONATE") ~ "Risedronic acid",
      str_detect(BP_CMDECOD, "ZOLEDRONIC|ZOLEDRONATE") ~ "Zoledronic acid",
      str_detect(DN_CMDECOD, "DENOSUMAB") ~ "Denosumab",
      TRUE ~ as.character(NA)
    )
  ) %>%
  select(BMA, ARM) %>% 
  tbl_summary(
    by = ARM,
    label = list(BMA ~ "Bone modifying agents"),
    missing = "no",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall()

# types of SRE
tbl_aa_sre <- LATITUDE %>%
  filter(ARM == "AAP") %>%
  filter(SRE == 1) %>%
  select(BMA_WITHIN, SPINAL_CORD_COMPRESSION, CLINICAL_FRACTURE, PALLATIVE_RADIATION_TO_BONE, SURGERY_TO_BONE) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      SPINAL_CORD_COMPRESSION ~ "Spinal cord compression",
      CLINICAL_FRACTURE ~ "Clinical fracture",
      PALLATIVE_RADIATION_TO_BONE ~ "Palliative radiation to bone",
      SURGERY_TO_BONE ~ "Surgery to bone"
    ),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall()

tbl_adt_sre <- LATITUDE %>%
  filter(ARM == "ADT") %>%
  filter(SRE == 1) %>%
  select(BMA_WITHIN, SPINAL_CORD_COMPRESSION, CLINICAL_FRACTURE, PALLATIVE_RADIATION_TO_BONE, SURGERY_TO_BONE) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      SPINAL_CORD_COMPRESSION ~ "Spinal cord compression",
      CLINICAL_FRACTURE ~ "Clinical fracture",
      PALLATIVE_RADIATION_TO_BONE ~ "Palliative radiation to bone",
      SURGERY_TO_BONE ~ "Surgery to bone"
    ),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall()

tbl_sre <- tbl_merge(tbls = list(tbl_aa_sre, tbl_adt_sre),
                     tab_spanner = c("AAP", "ADT"))

dist_bma <- ggplot(aes(x = DIST_BMA), data = LATITUDE %>% filter(!is.na(DIST_BMA))) + 
  geom_histogram(binwidth = 15, position = "identity", alpha = 0.8, fill = "#39568CFF") + 
  geom_vline(aes(xintercept = -90), linetype = "dashed", color = "grey") + 
  geom_vline(aes(xintercept = 90), linetype = "dashed", color = "grey") + 
  scale_x_continuous(breaks = seq(-1080, 1620, by = 90)) + 
  theme_jikei() + 
  ggsci::scale_fill_lancet() + 
  labs(
    x = "Days since randomization", 
    y = "Count"
  ) + 
  facet_wrap( ~ ARM, ncol = 1)

ggsave("dist_bma.pdf", dist_bma,
       units = "in", width = 8, height = 4)

all_bma <- LATITUDE %>% group_by(ARM) %>% count(.by = BMA_ALL) %>% flextable::flextable()

dur_bma <- ggplot(aes(x = DUR_BMA), data = LATITUDE %>% filter(!is.na(DUR_BMA))) + 
  geom_histogram(binwidth = 1.5, position = "identity", alpha = 0.8, fill = "#39568CFF") + 
  scale_x_continuous(breaks = seq(0, 72, by = 12)) + 
  theme_jikei() + 
  ggsci::scale_fill_lancet() + 
  labs(
    x = "Duration of BMA exposure (months)", 
    y = "Count"
  ) + 
  facet_wrap( ~ ARM, ncol = 1)

ggsave("dur_bma.pdf", units = "in", width = 10, height = 4)

tbl_dur_bma <- LATITUDE %>% 
  filter(!is.na(DUR_BMA)) %>% 
  filter(BMA_WITHIN == TRUE) %>% 
  select(DUR_BMA, ARM) %>%
  tbl_summary(
    by = ARM, 
    label = list(
      DUR_BMA ~ "Duration of BMA exposure, month"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall() 

# random forest imputation
set.seed(123)
IMP <- missForest(D, maxiter = 3, verbose = TRUE)
DF <- IMP$ximp

### IPTWs in AAP
iptw_formula <- as.formula(BMA_WITHIN ~ 
                             AGE_CAT + PS + PPI_WITHIN + DM + 
                             logPSA + BPI_CAT + GG_CAT + NUM_BONE_CAT + VISCERAL + LIVER + HGB + SCR + logLDH + ARM)

# check between-group balance
wgt_out <- weightit(
  iptw_formula,
  data = DF, method = "ps", stabilize = TRUE
)

balance <- bal.tab(wgt_out, abs = FALSE, un = TRUE, thresholds = c(m = 0.1)) # the summary of mean/max SMDs...
balance

DF <- tibble(DF, pscores = wgt_out$ps, weight = wgt_out$weights)

# print in table
bal <- rownames_to_column(balance$Balance) %>% 
  mutate(
    Diff.Un = Diff.Un %>% round(digits = 4),
    Diff.Adj = Diff.Adj %>% round(digits = 4)
  ) %>%
  rename(variable = rowname, diff_un = Diff.Un, diff_adj = Diff.Adj) %>%
  select(variable, diff_un, diff_adj) %>%
  flextable::flextable()

# standardized mean difference plot
label <- data.frame(
  old = c("AGE_CAT_less than 75", "PS_2", "PPI_WITHIN", "DM_yes", "logPSA",
          "BPI_CAT_4 or more", "GG_CAT_GG5", "NUM_BONE_CAT_more than 10", 
          "VISCERAL_yes", "LIVER_yes", "HGB", "SCR", "logLDH", "ARM_ADT"),
  new = c("Age at randomization (75 or more vs. less than 75)", "ECOG performance status (2 vs. 0-1)", 
          "Concomitant PPI use (Yes vs. No)", "Diabetes Mellitus (Yes vs. No)",
          "Baseline PSA (log-transformed)", "Brief pain inventory-SF3 (4 or more vs. 0-3)",
          "Grade group (GG5 vs. GG1-4)", "Number of bone metastasis (more than 10 vs. 0-10)", 
          "Visceral metastasis (Yes vs. No)", "Liver metastasis (Yes vs. No)", 
          "Baseline hemoglobin", "Baseline serum creatinine", "Baseline LDH (log-transformed)", "Arm (ADT vs. AAP)")
)

smdplot <- love.plot(wgt_out, drop.distance = TRUE, abs = FALSE, 
                     colors = c("#EE0000FF", "#39568CFF"), shapes = c("diamond filled", "circle filled"), sample.names = c("Crude", "Weighted"), 
                     r.thresholds = 0.1, size = 3.3, themes = theme_grey(), var.names = label) + 
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.title = element_blank(), 
        legend.position = "top") +
  labs(
    x = "Standarized mean difference (BMA users vs. non-users)"
  ) + 
  scale_x_continuous(n.breaks = 5) + 
  geom_vline(aes(xintercept = 0.1), linetype = "dashed", color = "grey35") + 
  geom_vline(aes(xintercept = -0.1), linetype = "dashed", color = "grey35") 

ggsave("smdplot_all.pdf", smdplot,
       units = "in", width = 8, height = 4)

# propensity score distributions
theme_set(theme_jikei() + theme(legend.position = "none"))
# propensity score distributions
psdist <- ggplot() + 
  geom_density(data = filter(DF, BMA_WITHIN == TRUE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN), alpha = 0.2) + 
  geom_density(data = filter(DF, BMA_WITHIN == FALSE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN, y = -..density..), alpha = 0.2) + 
  geom_density(data = filter(DF, BMA_WITHIN == TRUE),
               aes(x = pscores, fill = BMA_WITHIN), alpha = 0.6) + 
  geom_density(data = filter(DF, BMA_WITHIN == FALSE),
               aes(x = pscores, fill = BMA_WITHIN, y = -..density..), alpha = 0.6) +
  geom_hline(yintercept = 0, color = "white", size = 0.25) + 
  annotate(geom = "label", x = 0, y = 5, label = "BMA users\n(Pseudo-population)",
           fill = "#EE0000FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = 5, label = "BMA users\n(actual)",
           fill = "#EE0000FF", alpha = 0.6, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0, y = -5, label = "Non-users\n(Pseudo-population)",
           fill = "#3B4992FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = -5, label = "Non-users\n(actual)",
           fill = "#3B4992FF", alpha = 0.6, color = "white", hjust = 0) +
  scale_fill_aaas() + 
  scale_y_continuous(label = abs) + 
  labs(x = "Probability of BMA use", y = "Density")

ggsave("psdist_all.pdf", psdist,
       units = "in", width = 8, height = 4)

# Baseline characteristics
tbl1_crude <- D %>%
  mutate(
    NUM_BONE = NUM_BONE %>% as.integer()
  ) %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, ARM, BMA_WITHIN) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      AGE_CAT ~ "Age at randomization",
      PS ~ "ECOG performance status",
      PPI_WITHIN ~ "Concomitant PPI use",
      DM ~ "Diabetes Mellitus",
      logPSA ~ "Baseline PSA (log-transformed)",
      BPI_CAT ~ "Brief pain inventory-short form 3",
      GG_CAT ~ "Grade group",
      NUM_BONE_CAT ~ "Number of bone metastasis",
      VISCERAL ~ "Evidence of visceral metastasis",
      LIVER ~ "Evidence of liver metastasis",
      HGB ~ "Baseline hemoglobin concentration",
      SCR ~ "Baseline serum creatinine",
      logLDH ~ "Baseline LDH (log-transformed)", 
      ARM ~ "Treatment arm"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

tbl1_imp <- DF %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, ARM, BMA_WITHIN) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      AGE_CAT ~ "Age at randomization",
      PS ~ "ECOG performance status",
      PPI_WITHIN ~ "Concomitant PPI use",
      DM ~ "Diabetes Mellitus",
      logPSA ~ "Baseline PSA (log-transformed)",
      BPI_CAT ~ "Brief pain inventory-short form 3",
      GG_CAT ~ "Grade group",
      NUM_BONE_CAT ~ "Number of bone metastasis",
      VISCERAL ~ "Evidence of visceral metastasis",
      LIVER ~ "Evidence of liver metastasis",
      HGB ~ "Baseline hemoglobin concentration",
      SCR ~ "Baseline serum creatinine",
      logLDH ~ "Baseline LDH (log-transformed)", 
      ARM ~ "Treatment arm"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

tbl1_merge <- tbl_merge(list(tbl1_crude, tbl1_imp),
                        tab_spanner = c("Crude", "Imputed"))

# Kaplan-Meier plots
ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF)
ttsre_plot_crude <- ggsurvplot(ttsre_fit, 
                               risk.table = TRUE, 
                               title = "Time to skeletal-related event",
                               subtitle = "Crude population, the AAP cohort",
                               xlab = "Months since randomization", 
                               ylab = "SRE probability (95% CI)",
                               legend.labs = c("Non-users", "BMA users"),
                               censor = TRUE, censor.shape = "O", censor.size = 2.2,
                               palette = "lancet", size = 0.5,  break.time.by = 12,
                               ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                               risk.table.col = "strata",
                               tables.height = 0.12, risk.table.fontsize = 3.0,
                               tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_plot_crude$table <- ttsre_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF, weights = weight)
ttsre_plot_iptw <- ggsurvplot(ttsre_fit, 
                              risk.table = TRUE, 
                              title = "Time to skeletal-related event",
                              subtitle = "Weighted population, the AAP cohort",
                              xlab = "Months since randomization", 
                              ylab = "SRE probability (95% CI)",
                              legend.labs = c("Non-users", "BMA users"),
                              censor = TRUE, censor.shape = "O", censor.size = 2.2,
                              palette = "lancet", size = 0.5,  break.time.by = 12,
                              ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                              risk.table.col = "strata",
                              tables.height = 0.12, risk.table.fontsize = 3.0,
                              tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_plot_iptw$table <- ttsre_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            subtitle = "Crude population, the AAP cohort",
                            xlab = "Months since randomization", 
                            ylab = "OS probability (95% CI)",
                            legend.labs = c("Non-users", "BMA users"),
                            censor = TRUE, censor.shape = "O", censor.size = 2.2,
                            palette = "lancet", size = 0.5,  break.time.by = 12,
                            ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                            risk.table.col = "strata",
                            tables.height = 0.12, risk.table.fontsize = 3.0,
                            tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os_plot_crude$table <- os_plot_crude$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF, weights = weight)
os_plot_iptw <- ggsurvplot(os_fit, 
                           risk.table = TRUE, 
                           title = "Overall survival",
                           subtitle = "Weighted population, the AAP cohort",
                           xlab = "Months since randomization", 
                           ylab = "OS probability (95% CI)",
                           legend.labs = c("Non-users", "BMA users"),
                           censor = TRUE, censor.shape = "O", censor.size = 2.2,
                           palette = "lancet", size = 0.5,  break.time.by = 12,
                           ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                           risk.table.col = "strata",
                           tables.height = 0.12, risk.table.fontsize = 3.0,
                           tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

os_plot_iptw$table <- os_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

iptw_merge <- arrange_ggsurvplots(list(ttsre_plot_crude, ttsre_plot_iptw, os_plot_crude, os_plot_iptw), 
                                  nrow = 2, ncol = 2)

ggsave("iptw_merge.pdf", iptw_merge,
       units = "in", width = 8, height = 8)

# calculate IPTW-adjusted RMST
rmst_os <- akm_rmst(time = DF$FU, 
                    status = DF$DEATH, 
                    group = as.factor(DF$BMA_WITHIN), 
                    weight = DF$weight, 
                    tau = 64.2) 

rmst_ttsre <- akm_rmst(time = DF$TTSRE, 
                       status = DF$SRE, 
                       group = as.factor(DF$BMA_WITHIN), 
                       weight = DF$weight, 
                       tau = 63.9) 

# calculate last observed or censored survival time
os <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

ttsre <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

# effect modifier assessment
theme_gtsummary_journal(journal = "jama")
eod_effect <- coxph(Surv(TTSRE, SRE) ~ NUM_BONE_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
visceral_effect <- coxph(Surv(TTSRE, SRE) ~ VISCERAL * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
bpi_effect <- coxph(Surv(TTSRE, SRE) ~ BPI_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
age_effect <- coxph(Surv(TTSRE, SRE) ~ AGE_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
scr_effect <- coxph(Surv(TTSRE, SRE) ~ SCR * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
arm_effect <- coxph(Surv(TTSRE, SRE) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")

tbl_sre_int <- tbl_stack(
  tbls = list(eod_effect, visceral_effect, bpi_effect, age_effect, scr_effect, arm_effect))

theme_gtsummary_journal(journal = "jama")
eod_effect <- coxph(Surv(FU, DEATH) ~ NUM_BONE_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
visceral_effect <- coxph(Surv(FU, DEATH) ~ VISCERAL * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
bpi_effect <- coxph(Surv(FU, DEATH) ~ BPI_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
age_effect <- coxph(Surv(FU, DEATH) ~ AGE_CAT * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
scr_effect <- coxph(Surv(FU, DEATH) ~ SCR * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
arm_effect <- coxph(Surv(FU, DEATH) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")

tbl_os_int <- tbl_stack(
  tbls = list(eod_effect, visceral_effect, bpi_effect, age_effect, scr_effect, arm_effect))

# SRE components
arm_prtb_effect <- coxph(Surv(TTSRE_PRTB, SRE_PRTB) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
arm_cf_effect <- coxph(Surv(TTSRE_CF, SRE_CF) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
arm_scc_effect <- coxph(Surv(TTSRE_SCC, SRE_SCC) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
arm_stb_effect <- coxph(Surv(TTSRE_STB, SRE_STB) ~ ARM * BMA_WITHIN, data = DF, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")

tbl_sre_arm_int <- tbl_stack(
  tbls = list(arm_prtb_effect, arm_cf_effect, arm_scc_effect, arm_stb_effect))