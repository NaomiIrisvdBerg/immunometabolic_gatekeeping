##incidence distributions of cancer types (WHO, https://ci5.iarc.fr/ci5plus/) 
#by patient sex and age bin
#this code reproduces the figure and analysis generated in Supporting materials IIA: Tropism in peadatric cancers explored in the context of
#immunometabolic gatekeeping

##load libraries first
library(readr) ##we ran the code with v 2.1.5
library(dplyr) ##v 1.1.4
library(stringr) ##v 1.5.2
library(purrr) ##v 1.1.0
library(tibble) ##v 3.3.0
library(tidytext) ##v 0.4.3
library(ggplot2) ##v 4.0.0

##read in the incidences files (corresponding to focus tissues as listed in Table 1 of the manuscript):
Adrenal_gland_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Adrenal gland_age_sex_incidence.csv")
Bladder_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Bladder_age_sex_incidence.csv")
Brain_central_nervous_system_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Brain, central nervous system_age_sex_incidence.csv")
Colon_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Colon_age_sex_incidence.csv")
Eye_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Eye_age_sex_incidence.csv")
Kidney_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Kidney_age_sex_incidence.csv")
Liver_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Liver_age_sex_incidence.csv")
Lung_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Lung_age_sex_incidence.csv")
Melanoma_of_skin_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Melanoma of skin_age_sex_incidence.csv")
Pancreas_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Pancreas_age_sex_incidence.csv")
Stomach_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Stomach_age_sex_incidence.csv")
Oesophagus_age_sex_incidence <- read_csv("Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Oesophagus_age_sex_incidence.csv")


##refer to all files:
cancer_files <- list(
  Adrenal_gland = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Adrenal gland_age_sex_incidence.csv",
  Bladder = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Bladder_age_sex_incidence.csv",
  Brain_CNS = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Brain, central nervous system_age_sex_incidence.csv",
  Colon = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Colon_age_sex_incidence.csv",
  Eye = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Eye_age_sex_incidence.csv",
  Kidney = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Kidney_age_sex_incidence.csv",
  Liver = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Liver_age_sex_incidence.csv",
  Lung = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Lung_age_sex_incidence.csv",
  Melanoma_skin = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Melanoma of skin_age_sex_incidence.csv",
  Pancreas = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Pancreas_age_sex_incidence.csv",
  Stomach = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Stomach_age_sex_incidence.csv",
  Oesophagus = "Downloads/CI5_PAN_BROAD_INCIDENCE_CSV/Oesophagus_age_sex_incidence.csv"
)


##now compare the peak incidences (such as presented in the manuscript's Supporting materials IIA: Tropism in peadatric cancers explored in the context of
#immunometabolic gatekeeping)

get_age_mid <- function(age_range) {
  nums <- str_extract_all(age_range, "[0-9]+")
  sapply(nums, function(x) if (length(x) == 2) mean(as.numeric(x)) else as.numeric(x))
}

calc_peak_ratio <- function(df, split_age = 20) {
  df %>%
    dplyr::mutate(age_mid = get_age_mid(age_range)) %>%
    dplyr::arrange(sex, age_mid) %>%
    dplyr::group_split(sex) %>%
    purrr::map_dfr(function(g) {
      s <- unique(g$sex)
      early <- g %>% filter(age_mid < split_age)
      adult <- g %>% filter(age_mid >= split_age)
      
      early_peak      <- if (nrow(early)) early$cases[which.max(early$cases)] else NA_real_
      early_peak_age  <- if (nrow(early)) early$age_range[which.max(early$cases)] else NA_character_
      adult_peak      <- if (nrow(adult)) adult$cases[which.max(adult$cases)] else NA_real_
      adult_peak_age  <- if (nrow(adult)) adult$age_range[which.max(adult$cases)] else NA_character_
      
      tibble(
        sex = s,
        early_peak     = early_peak,
        early_peak_age = early_peak_age,
        adult_peak     = adult_peak,
        adult_peak_age = adult_peak_age,
        peak_ratio_early_over_adult = early_peak / adult_peak,
        peak_ratio_adult_over_early = adult_peak / early_peak
      )
    })
}


peak_results <- map_dfr(names(cancer_files), function(name) {
  df <- read_csv(cancer_files[[name]], show_col_types = FALSE)
  res <- calc_peak_ratio(df, split_age = 15)   ##adjustadble; we selected for bins under 15yo
  res$cancer_site <- name
  res
}) %>%
  dplyr::select(cancer_site, sex, early_peak, early_peak_age, adult_peak, adult_peak_age,
         peak_ratio_early_over_adult, peak_ratio_adult_over_early)


#View(peak_results)
#z#write_csv(peak_results, "~/Downloads/bimodality_peak_ratios_all_cancers.csv")


str(peak_results)


peak_results_plot <- peak_results %>%
  ##keep only valid positive ratios for log10
  dplyr::filter(is.finite(peak_ratio_early_over_adult),
         peak_ratio_early_over_adult > 0,
         !is.na(cancer_site), !is.na(early_peak), !is.na(early_peak_age),
         !is.na(adult_peak_age)) %>%
  dplyr::mutate(
    log10_ratio = log10(peak_ratio_early_over_adult),
    ##order cancer_site within each sex by log10 ratio
    cancer_site_ord = reorder_within(cancer_site, log10_ratio, sex),
    sex = factor(sex, levels = c(1,2), labels = c("Male","Female"))
  )

ggplot(peak_results_plot,
       aes(x = log10_ratio,
           y = cancer_site_ord,
           size = early_peak,
           colour = early_peak_age,
           shape = adult_peak_age)) +
  geom_point(alpha = 0.8, stroke = 1.1) +
  scale_y_reordered() +
  scale_size_continuous(name = "Early peak (cases)", range = c(3, 12)) +
  scale_x_continuous(name = expression(log[10]("early/adult peak ratio"))) +
  labs(
    y = "Cancer site (ordered by log10 ratio within sex)",
    colour = "Early peak age",
    shape  = "Adult peak age",
    title  = "Bimodality in incidence: young vs adult peaks by cancer site and sex"
  ) +
  facet_wrap(~ sex, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.3, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black")
  )

