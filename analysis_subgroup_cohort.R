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
    HGB = as.numeric(HGB) / 10, # g/l to g/dl
    SCR = as.numeric(SCR) / 88.4, # umol/l to mg/dl
    logPSA = log(PSA),
    logLDH = log(LDH),
    ARM = ARM %>% factor(levels = c("ADT", "AAP"))
  ) %>% 
  select(BMA_WITHIN, PPI_WITHIN, METFORMIN_WITHIN, STATIN_WITHIN, GC_WITHIN, 
         PFS, PROG, FU, DEATH, 
         TTSRE, TTSRE_PRTB, TTSRE_CF, TTSRE_SCC, TTSRE_STB, 
         SRE, SRE_PRTB, SRE_CF, SRE_SCC, SRE_STB, 
         LOCAL, DM, 
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

# types of 1st SRE
tbl_aa_sre <- LATITUDE %>%
  filter(ARM == "AAP") %>%
  filter(SRE == 1) %>%
  select(BMA_WITHIN, SPINAL_CORD_COMPRESSION, CLINICAL_FRACTURE, PALLIATIVE_RADIATION_TO_BONE, SURGERY_TO_BONE) %>%
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
      PALLIATIVE_RADIATION_TO_BONE ~ "Palliative radiation to bone",
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

# types of all SRE
tbl_all_aa_sre <- LATITUDE %>%
  filter(ARM == "AAP") %>%
  filter(SRE == 1) %>%
  select(BMA_WITHIN, SRE_PRTB, SRE_SCC, SRE_CF, SRE_STB) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      SRE_SCC ~ "Spinal cord compression",
      SRE_CF ~ "Clinical fracture",
      SRE_PRTB ~ "Palliative radiation to bone",
      SRE_STB ~ "Surgery to bone"
    ),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall()

tbl_all_adt_sre <- LATITUDE %>%
  filter(ARM == "ADT") %>%
  filter(SRE == 1) %>%
  select(BMA_WITHIN, SRE_PRTB, SRE_SCC, SRE_CF, SRE_STB) %>%
  tbl_summary(
    by = BMA_WITHIN,
    label = list(
      SRE_SCC ~ "Spinal cord compression",
      SRE_CF ~ "Clinical fracture",
      SRE_PRTB ~ "Palliative radiation to bone",
      SRE_STB ~ "Surgery to bone"
    ),
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall()

tbl_all_sre <- tbl_merge(tbls = list(tbl_all_aa_sre, tbl_all_adt_sre),
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

all_bma <- LATITUDE %>% filter(!is.na(DIST_BMA)) %>% group_by(ARM) %>% count(.by = BMA_ALL) %>% flextable::flextable()

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
    statistic = list(all_continuous() ~ "{median} ({p25}--{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  ) %>%
  add_overall() 

# random forest imputation
set.seed(123)
IMP <- missForest(D, maxiter = 3, verbose = TRUE)
DF <- IMP$ximp

# baseline characteristics for all patients
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
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  SCR ~ c(2, 2),
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
    statistic = list(all_continuous() ~ "{median} ({p25}--{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  SCR ~ c(2, 2),
                  all_categorical() ~ c(0, 1))
  )

tbl1_merge <- tbl_merge(list(tbl1_crude, tbl1_imp),
                        tab_spanner = c("Crude", "Imputed"))

### IPTWs in AAP
DF_AA <- DF %>% filter(ARM == "AAP")

iptw_formula <- as.formula(BMA_WITHIN ~ 
                             AGE_CAT + PS + PPI_WITHIN + DM + 
                             logPSA + BPI_CAT + GG_CAT + NUM_BONE_CAT + VISCERAL + LIVER + HGB + SCR + logLDH)

# check between-group balance
wgt_aa_out <- weightit(
  iptw_formula,
  data = DF_AA, method = "ps", stabilize = TRUE
)

balance_aa <- bal.tab(wgt_aa_out, abs = FALSE, un = TRUE, thresholds = c(m = 0.1)) # the summary of mean/max SMDs...
balance_aa

DF_AA <- tibble(DF_AA, pscores = wgt_aa_out$ps, weight = wgt_aa_out$weights)

# print in table
bal_aa <- rownames_to_column(balance_aa$Balance) %>% 
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
          "VISCERAL_yes", "LIVER_yes", "HGB", "SCR", "logLDH"),
  new = c("Age at randomization (75 or more vs. less than 75)", "ECOG performance status (2 vs. 0-1)", 
          "Concomitant PPI use (Yes vs. No)", "Diabetes Mellitus (Yes vs. No)",
          "Baseline PSA (log-transformed)", "Brief pain inventory-SF3 (4 or more vs. 0-3)",
          "Grade group (GG5 vs. GG1-4)", "Number of bone metastasis (more than 10 vs. 0-10)", 
          "Visceral metastasis (Yes vs. No)", "Liver metastasis (Yes vs. No)", 
          "Baseline hemoglobin", "Baseline serum creatinine", "Baseline LDH (log-transformed)")
)

smdplot <- love.plot(wgt_aa_out, drop.distance = TRUE, abs = FALSE, 
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

ggsave("smdplot_aap.pdf", smdplot,
       units = "in", width = 8, height = 4)

# propensity score distributions
theme_set(theme_jikei() + theme(legend.position = "none"))
# propensity score distributions
psdist_aa <- ggplot() + 
  geom_density(data = filter(DF_AA, BMA_WITHIN == TRUE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN), alpha = 0.2) + 
  geom_density(data = filter(DF_AA, BMA_WITHIN == FALSE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN, y = -..density..), alpha = 0.2) + 
  geom_density(data = filter(DF_AA, BMA_WITHIN == TRUE),
               aes(x = pscores, fill = BMA_WITHIN), alpha = 0.6) + 
  geom_density(data = filter(DF_AA, BMA_WITHIN == FALSE),
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

ggsave("psdist_aa.pdf", psdist_aa,
       units = "in", width = 8, height = 4)

# Baseline characteristics
tbl1_aap_crude <- LATITUDE %>%
  filter(ARM == "AAP") %>%
  mutate(
    NUM_BONE = NUM_BONE %>% as.integer(),
    logPSA = log(PSA),
    logLDH = log(LDH)
  ) %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, BMA_WITHIN) %>%
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
      logLDH ~ "Baseline LDH (log-transformed)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

tbl1_aap_imp <- DF_AA %>%
  filter(ARM == "AAP") %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, BMA_WITHIN) %>%
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
      logLDH ~ "Baseline LDH (log-transformed)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

# time to event distribution
library(pec)
fu <- quantile(prodlim(Hist(FU, DEATH) ~ 1, data = DF_AA, reverse = TRUE))

# Kaplan-Meier plots
ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_AA)
ttsre_plot_crude <- ggsurvplot(ttsre_fit, 
                               risk.table = TRUE, 
                               title = "Time to skeletal-related event",
                               subtitle = "Crude, the AAP cohort",
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

ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_AA, weights = weight)
ttsre_plot_iptw <- ggsurvplot(ttsre_fit, 
                              risk.table = TRUE, 
                              title = "Time to skeletal-related event",
                              subtitle = "IPTW-adjusted, the AAP cohort",
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

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_AA)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            subtitle = "Crude, the AAP cohort",
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

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_AA, weights = weight)
os_plot_iptw <- ggsurvplot(os_fit, 
                           risk.table = TRUE, 
                           title = "Overall survival",
                           subtitle = "IPTW-adjusted, the AAP cohort",
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

iptw_aap_merge <- arrange_ggsurvplots(list(ttsre_plot_crude, ttsre_plot_iptw, os_plot_crude, os_plot_iptw), 
                                      nrow = 2, ncol = 2)

ggsave("iptw_aap.pdf", iptw_aap_merge,
       units = "in", width = 8, height = 8)

# calculate IPTW-adjusted RMST
rmst_aa_os <- akm_rmst(time = DF_AA$FU, 
                       status = DF_AA$DEATH, 
                       group = as.factor(DF_AA$BMA_WITHIN), 
                       weight = DF_AA$weight, 
                       tau = 64.2) 

rmst_aa_ttsre <- akm_rmst(time = DF_AA$TTSRE, 
                          status = DF_AA$SRE, 
                          group = as.factor(DF_AA$BMA_WITHIN), 
                          weight = DF_AA$weight, 
                          tau = 63.9) 

# calculate last observed or censored survival time
os_aa <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_AA, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

ttsre_aa <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_AA, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

# effect modifier assessment
theme_gtsummary_journal(journal = "jama")
eod_effect <- coxph(Surv(TTSRE, SRE) ~ NUM_BONE_CAT * BMA_WITHIN, data = DF_AA, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
visceral_effect <- coxph(Surv(TTSRE, SRE) ~ VISCERAL * BMA_WITHIN, data = DF_AA, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
bpi_effect <- coxph(Surv(TTSRE, SRE) ~ BPI_CAT * BMA_WITHIN, data = DF_AA, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
age_effect <- coxph(Surv(TTSRE, SRE) ~ AGE_CAT * BMA_WITHIN, data = DF_AA, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")
scr_effect <- coxph(Surv(TTSRE, SRE) ~ SCR * BMA_WITHIN, data = DF_AA, weights = weight) %>% tbl_regression(exponentiate = TRUE) %>% add_n(location = "level")

tbl_int <- tbl_stack(
  tbls = list(eod_effect, visceral_effect, bpi_effect, age_effect, scr_effect)
)

# details of SREs
# Kaplan-Meier plots
ttsre_prtb_fit <- survfit(Surv(TTSRE_PRTB, SRE_PRTB) ~ BMA_WITHIN, data = DF_AA, weights = weight)
ttsre_prtb_plot_iptw <- ggsurvplot(ttsre_prtb_fit, 
                                   risk.table = TRUE, 
                                   title = "Time to palliative radiation to bone",
                                   subtitle = "IPTW-adjusted, the AAP cohort",
                                   xlab = "Months since randomization", 
                                   ylab = "Event probability (95% CI)",
                                   legend.labs = c("Non-users", "BMA users"),
                                   censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                   palette = "lancet", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.12, risk.table.fontsize = 3.0,
                                   tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_prtb_plot_iptw$table <- ttsre_prtb_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ttsre_stb_fit <- survfit(Surv(TTSRE_STB, SRE_STB) ~ BMA_WITHIN, data = DF_AA, weights = weight)
ttsre_stb_plot_iptw <- ggsurvplot(ttsre_stb_fit, 
                                  risk.table = TRUE, 
                                  title = "Time to surgery to bone",
                                  subtitle = "IPTW-adjusted, the AAP cohort",
                                  xlab = "Months since randomization", 
                                  ylab = "Event probability (95% CI)",
                                  legend.labs = c("Non-users", "BMA users"),
                                  censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                  palette = "lancet", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.12, risk.table.fontsize = 3.0,
                                  tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_stb_plot_iptw$table <- ttsre_stb_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ttsre_cf_fit <- survfit(Surv(TTSRE_CF, SRE_CF) ~ BMA_WITHIN, data = DF_AA, weights = weight)
ttsre_cf_plot_iptw <- ggsurvplot(ttsre_cf_fit, 
                                 risk.table = TRUE, 
                                 title = "Time to clinical fracture",
                                 subtitle = "IPTW-adjusted, the AAP cohort",
                                 xlab = "Months since randomization", 
                                 ylab = "Event probability (95% CI)",
                                 legend.labs = c("Non-users", "BMA users"),
                                 censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                 palette = "lancet", size = 0.5,  break.time.by = 12,
                                 ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                 risk.table.col = "strata",
                                 tables.height = 0.12, risk.table.fontsize = 3.0,
                                 tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_cf_plot_iptw$table <- ttsre_cf_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ttsre_scc_fit <- survfit(Surv(TTSRE_SCC, SRE_SCC) ~ BMA_WITHIN, data = DF_AA, weights = weight)
ttsre_scc_plot_iptw <- ggsurvplot(ttsre_scc_fit, 
                                  risk.table = TRUE, 
                                  title = "Time to spinal cord compression",
                                  subtitle = "IPTW-adjusted, the AAP cohort",
                                  xlab = "Months since randomization", 
                                  ylab = "Event probability (95% CI)",
                                  legend.labs = c("Non-users", "BMA users"),
                                  censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                  palette = "lancet", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.12, risk.table.fontsize = 3.0,
                                  tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_scc_plot_iptw$table <- ttsre_scc_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

iptw_ttsre_aap_merge <- arrange_ggsurvplots(list(ttsre_prtb_plot_iptw, ttsre_scc_plot_iptw, ttsre_cf_plot_iptw, ttsre_stb_plot_iptw), 
                                            nrow = 2, ncol = 2)

ggsave("iptw_ttsre_aap.pdf", iptw_ttsre_aap_merge,
       units = "in", width = 8, height = 8)

# ADT cohort
# Baseline characteristics
tbl1_adt_crude <- LATITUDE %>%
  filter(ARM == "ADT") %>%
  mutate(
    NUM_BONE = NUM_BONE %>% as.integer(),
    logPSA = log(PSA),
    logLDH = log(LDH)
  ) %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, BMA_WITHIN) %>%
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
      logLDH ~ "Baseline LDH (log-transformed)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

tbl1_merge <- tbl_merge(tbls = list(tbl1_aap_crude, tbl1_adt_crude),
                        tab_spanner = c("AAP", "ADT"))

tbl1_adt_imp <- DF_ADT %>%
  filter(ARM == "ADT") %>%
  select(AGE_CAT, PS, PPI_WITHIN, DM, logPSA, BPI_CAT, GG_CAT, NUM_BONE_CAT, VISCERAL, LIVER, HGB, SCR, logLDH, BMA_WITHIN) %>%
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
      logLDH ~ "Baseline LDH (log-transformed)"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>% 
  add_overall(
    statistic = list(all_continuous() ~ "{median} ({p25} to {p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))
  )

tbl1_merge_imp <- tbl_merge(tbls = list(tbl1_aap_imp, tbl1_adt_imp),
                            tab_spanner = c("AAP", "ADT"))

# check between-group balance using MatchThem
### IPTWs in ADT
DF_ADT <- DF %>% filter(ARM == "ADT")

wgt_adt_out <- weightit(
  iptw_formula,
  data = DF_ADT, method = "ps", stabilize = TRUE
)

balance_adt <- bal.tab(wgt_adt_out, abs = TRUE, un = TRUE, thresholds = c(m = 0.1)) # the summary of mean/max SMDs...
balance_adt

DF_ADT <- tibble(DF_ADT, pscores = wgt_adt_out$ps, weight = wgt_adt_out$weights)

# print in table
bal_adt <- rownames_to_column(balance_adt$Balance) %>% 
  mutate(
    Diff.Un = Diff.Un %>% round(digits = 4),
    Diff.Adj = Diff.Adj %>% round(digits = 4)
  ) %>%
  rename(variable = rowname, diff_un = Diff.Un, diff_adj = Diff.Adj) %>%
  select(variable, diff_un, diff_adj) %>%
  flextable::flextable()

# standardized mean difference plot
smdplot <- love.plot(wgt_adt_out, drop.distance = TRUE, abs = FALSE, 
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

ggsave("smdplot_adt.pdf", units = "in", width = 8, height = 4)

# propensity score distributions
theme_set(theme_jikei() + theme(legend.position = "none"))
# propensity score distributions
psdist_adt <- ggplot() + 
  geom_density(data = filter(DF_ADT, BMA_WITHIN == TRUE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN), alpha = 0.2) + 
  geom_density(data = filter(DF_ADT, BMA_WITHIN == FALSE),
               aes(x = pscores, weight = weight, fill = BMA_WITHIN, y = -..density..), alpha = 0.2) + 
  geom_density(data = filter(DF_ADT, BMA_WITHIN == TRUE),
               aes(x = pscores, fill = BMA_WITHIN), alpha = 0.6) + 
  geom_density(data = filter(DF_ADT, BMA_WITHIN == FALSE),
               aes(x = pscores, fill = BMA_WITHIN, y = -..density..), alpha = 0.6) +
  geom_hline(yintercept = 0, color = "white", size = 0.25) + 
  annotate(geom = "label", x = 0.05, y = 5, label = "BMA users\n(Pseudo-population)",
           fill = "#EE0000FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = 7.5, label = "BMA users\n(actual)",
           fill = "#EE0000FF", alpha = 0.6, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.05, y = -5, label = "Non-users\n(Pseudo-population)",
           fill = "#3B4992FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = -5, label = "Non-users\n(actual)",
           fill = "#3B4992FF", alpha = 0.6, color = "white", hjust = 0) +
  scale_fill_aaas() + 
  scale_y_continuous(label = abs) + 
  labs(x = "Probability of BMA use", y = "Density")

ggsave("psdist_adt.pdf", psdist_adt,
       units = "in", width = 8, height = 4)

# Kaplan-Meier plots
ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_ADT)
ttsre_plot_crude <- ggsurvplot(ttsre_fit, 
                               risk.table = TRUE, 
                               title = "Time to skeletal-related event",
                               subtitle = "Crude, the ADT cohort",
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

ttsre_fit <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
ttsre_plot_iptw <- ggsurvplot(ttsre_fit, 
                              risk.table = TRUE, 
                              title = "Time to skeletal-related event",
                              subtitle = "IPTW-adjusted, the ADT cohort",
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

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_ADT)
os_plot_crude <- ggsurvplot(os_fit, 
                            risk.table = TRUE, 
                            title = "Overall survival",
                            subtitle = "Crude, the ADT cohort",
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

os_fit <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
os_plot_iptw <- ggsurvplot(os_fit, 
                           risk.table = TRUE, 
                           title = "Overall survival",
                           subtitle = "IPTW-adjusted, the ADT cohort",
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

iptw_adt_merge <- arrange_ggsurvplots(list(ttsre_plot_crude, ttsre_plot_iptw, os_plot_crude, os_plot_iptw), 
                                      nrow = 2, ncol = 2)
ggsave("iptw_adt.pdf", iptw_adt_merge,
       units = "in", width = 8, height = 8)

# calculate IPTW-adjusted RMST
rmst_adt_os <- akm_rmst(time = DF_ADT$FU, 
                        status = DF_ADT$DEATH, 
                        group = as.factor(DF_ADT$BMA_WITHIN), 
                        weight = DF_ADT$weight, 
                        tau = 63.9)

rmst_adt_ttsre <- akm_rmst(time = DF_ADT$TTSRE, 
                           status = DF_ADT$SRE, 
                           group = as.factor(DF_ADT$BMA_WITHIN), 
                           weight = DF_ADT$weight, 
                           tau = 63.9) 

os_adt <- survfit(Surv(FU, DEATH) ~ BMA_WITHIN, data = DF_ADT, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

ttsre_adt <- survfit(Surv(TTSRE, SRE) ~ BMA_WITHIN, data = DF_ADT, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

# Kaplan-Meier plots
ttsre_prtb_fit <- survfit(Surv(TTSRE_PRTB, SRE_PRTB) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
ttsre_prtb_plot_iptw <- ggsurvplot(ttsre_prtb_fit, 
                                   risk.table = TRUE, 
                                   title = "Time to palliative radiation to bone",
                                   subtitle = "IPTW-adjusted, the ADT cohort",
                                   xlab = "Months since randomization", 
                                   ylab = "Event probability (95% CI)",
                                   legend.labs = c("Non-users", "BMA users"),
                                   censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                   palette = "lancet", size = 0.5,  break.time.by = 12,
                                   ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                   risk.table.col = "strata",
                                   tables.height = 0.12, risk.table.fontsize = 3.0,
                                   tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_prtb_plot_iptw$table <- ttsre_prtb_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ttsre_stb_fit <- survfit(Surv(TTSRE_STB, SRE_STB) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
ttsre_stb_plot_iptw <- ggsurvplot(ttsre_stb_fit, 
                                  risk.table = TRUE, 
                                  title = "Time to surgery to bone",
                                  subtitle = "IPTW-adjusted, the ADT cohort",
                                  xlab = "Months since randomization", 
                                  ylab = "Event probability (95% CI)",
                                  legend.labs = c("Non-users", "BMA users"),
                                  censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                  palette = "lancet", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.12, risk.table.fontsize = 3.0,
                                  tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_stb_plot_iptw$table <- ttsre_stb_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ttsre_cf_fit <- survfit(Surv(TTSRE_CF, SRE_CF) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
ttsre_cf_plot_iptw <- ggsurvplot(ttsre_cf_fit, 
                                 risk.table = TRUE, 
                                 title = "Time to clinical fracture",
                                 subtitle = "IPTW-adjusted, the ADT cohort",
                                 xlab = "Months since randomization", 
                                 ylab = "Event probability (95% CI)",
                                 legend.labs = c("Non-users", "BMA users"),
                                 censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                 palette = "lancet", size = 0.5,  break.time.by = 12,
                                 ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                 risk.table.col = "strata",
                                 tables.height = 0.12, risk.table.fontsize = 3.0,
                                 tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_cf_plot_iptw$table <- ttsre_cf_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ttsre_scc_fit <- survfit(Surv(TTSRE_SCC, SRE_SCC) ~ BMA_WITHIN, data = DF_ADT, weights = weight)
ttsre_scc_plot_iptw <- ggsurvplot(ttsre_scc_fit, 
                                  risk.table = TRUE, 
                                  title = "Time to spinal cord compression",
                                  subtitle = "IPTW-adjusted, the ADT cohort",
                                  xlab = "Months since randomization", 
                                  ylab = "Event probability (95% CI)",
                                  legend.labs = c("Non-users", "BMA users"),
                                  censor = TRUE, censor.shape = "O", censor.size = 2.2,
                                  palette = "lancet", size = 0.5,  break.time.by = 12,
                                  ggtheme = theme_jikei(), risk.table.title = "Number at risk",
                                  risk.table.col = "strata",
                                  tables.height = 0.12, risk.table.fontsize = 3.0,
                                  tables.theme = survminer::theme_cleantable(), conf.int = TRUE) 

ttsre_scc_plot_iptw$table <- ttsre_scc_plot_iptw$table + 
  theme_void(base_size = 9) + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

iptw_ttsre_adt_merge <- arrange_ggsurvplots(list(ttsre_prtb_plot_iptw, ttsre_scc_plot_iptw, ttsre_cf_plot_iptw, ttsre_stb_plot_iptw), 
                                            nrow = 2, ncol = 2)

ggsave("iptw_ttsre_adt.pdf", iptw_ttsre_adt_merge,
       units = "in", width = 8, height = 8)
