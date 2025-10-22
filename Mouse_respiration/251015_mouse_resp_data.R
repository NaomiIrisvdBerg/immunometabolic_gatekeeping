##Mouse cross-tissue respiration data as visualised in Supporting materials IB – Cross-tissue mitochondrial respiration in normal mice
#source data from: Mitochondrial respiration atlas reveals differential changes in mitochondrial function across sex and age
#https://elifesciences.org/reviewed-preprints/96926v2?utm
#Link directing to Figure 1—source data 1: 
#https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvOTY5MjYvZWxpZmUtOTY5MjYtZmlnMS1kYXRhMS12MS54bHN4/elife-96926-fig1-data1-v1.xlsx?_hash=Bhx%2Bw4PQpD3o4V0qFtqLtptKWqbzi0MRdgt06rmNt20%3D


##load necesssary libraries:
library(readxl) ##we ran code with: v 1.4.5
library(purrr) ##v 1.1.0
library(dplyr) ##v 1.1.4
library(tidyr) ##v 1.3.1
library(stringr) ##v 1.5.2
library(ggplot2) ##v 4.0.0
library(forcats) ##v 1.0.1


##example of unpacking for BAT. Change directory if not in local Downloads folder:
elife_96926_fig1_data2_v1_2_a_Young_BAT <- readxl::read_excel("Downloads/elife-96926-fig1-data2-v1.xlsx", sheet = "Young-BAT")[, 4:37]

##unpack for all tissues (one per sheet)
file_path <- "Downloads/elife-96926-fig1-data2-v1 (2).xlsx"

sheets <- excel_sheets(file_path)




##34:37: Complex Linked Resp. Averages
##78:81: Complex Linked Resp. Avg/MTDR: How hard are individual mitochondria working?




elife_data_list <- map(sheets, ~ read_excel(file_path, sheet = .x)[34:37, 4:37]) %>%
  set_names(sheets)



tidy_one_sheet <- function(df, sheet_name) {
  ##separate header & data
  header_row <- df[1, ]
  data_rows  <- df[-1, ]
  
  ##first column == mito. complex
  colnames(data_rows)[1] <- "complex"
  
  colnames(data_rows)[-1] <- make.names(
    ifelse(is.na(header_row[-1]) | header_row[-1] == "",
           paste0("col_", seq_along(header_row[-1])),
           header_row[-1]),
    unique = TRUE
  )
  
  data_rows[] <- lapply(data_rows, as.character)
  
  ##transform to long format:
  tidy_df <- data_rows %>%
    tidyr::pivot_longer(
      cols = -complex,
      names_to = "ID",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      sex = case_when(
        str_detect(ID, "^YM") ~ "Male",
        str_detect(ID, "^YF") ~ "Female",
        TRUE ~ NA_character_
      ),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::select(complex, ID, sex, value)
  
  ##extract metadata from sheet name: e.g. "Young-BAT" contains age =="Young", tissue =="BAT"
  meta <- str_split_fixed(sheet_name, "-", 2)
  tidy_df <- tidy_df %>%
    dplyr::mutate(
      age = meta[1],
      tissue = meta[2]
    )
  
  return(tidy_df)
}

##Loop over all sheets in the excel file (so all tissues):
tidy_all <- purrr::imap_dfr(elife_data_list, tidy_one_sheet)



##match tissues as listed in Table 1
tidy_all <- tidy_all %>%
  dplyr::mutate(
    focus_tissue = case_when(
      ##kidney
      tissue %in% c("Kid-Cortex", "Kid-Med") ~ "Kidney",
      
      ##brain
      tissue %in% "Cortex" ~ "Brain",
      
      ##eye
      tissue == "Eye" ~ "Eye",
      
      ##stomach
      tissue == "Stomach" ~ "Stomach",
      
      ##Pancreas
      tissue == "Pancreas" ~ "Pancreas",
      
      ##Lung
      tissue == "Lung" ~ "Lung",
      
      ##Liver
      tissue == "Liver" ~ "Liver",
      
      ##Fallopian tubes (only for females)
      tissue == "Test Fallop tube" & sex == "Female" ~ "Fallopian tube",
      
      ##Skin
      tissue == "Skin" ~ "Skin",
      
      #Colon 
      tissue %in% c("Proximal colon", "Distal Colon") ~ "Colon",
      
      #everything else
      TRUE ~ "none"
    )
  )


##if focusing on terminal complex only:
tidy_all_I <- subset(tidy_all, complex == "CIV")


focus_order <- tidy_all_I %>%
  dplyr::filter(focus_tissue != "none") %>%
  dplyr::group_by(focus_tissue) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) %>%
  dplyr::arrange(desc(median_value)) %>%
  dplyr::pull(focus_tissue)

tidy_all_I %>%
  dplyr::filter(focus_tissue != "none") %>%
  dplyr::mutate(focus_tissue = factor(focus_tissue, levels = focus_order)) %>%
  ggplot(aes(x = focus_tissue, y = log10(value + 0.01), fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
  geom_jitter(aes(color = sex), width = 0.15, alpha = 0.5, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Focus tissue (ordered by median metabolic value)",
    y = "Metabolic value",
    fill = "Sex",
    color = "Sex",
    title = "Metabolic value: log10(Complex Linked Resp. Averages for CIV) 
across focus tissues in young mice (ordered by median)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )






##next; across all complexes, normalised by mitochondrial density: Complex Linked Resp. Avg/MTDR

file_path <- "Downloads/elife-96926-fig1-data2-v1 (2).xlsx"

sheets <- excel_sheets(file_path)

elife_data_list <- map(sheets, ~ read_excel(file_path, sheet = .x)[78:81, 4:37]) %>%
  set_names(sheets)



tidy_one_sheet <- function(df, sheet_name) {
  ##separate header and data
  header_row <- df[1, ]
  data_rows  <- df[-1, ]
  
  ##first column == mito. complex
  colnames(data_rows)[1] <- "complex"
  
  colnames(data_rows)[-1] <- make.names(
    ifelse(is.na(header_row[-1]) | header_row[-1] == "",
           paste0("col_", seq_along(header_row[-1])),
           header_row[-1]),
    unique = TRUE
  )
  
  data_rows[] <- lapply(data_rows, as.character)
  
  ##transform to long format
  tidy_df <- data_rows %>%
    tidyr::pivot_longer(
      cols = -complex,
      names_to = "ID",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      sex = case_when(
        str_detect(ID, "^YM") ~ "Male",
        str_detect(ID, "^YF") ~ "Female",
        TRUE ~ NA_character_
      ),
      value = suppressWarnings(as.numeric(value))
    ) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::select(complex, ID, sex, value)
  
  ##extract metadata from sheet name: e.g. "Young-BAT" contains age =="Young", tissue =="BAT"
  meta <- str_split_fixed(sheet_name, "-", 2)
  tidy_df <- tidy_df %>%
    dplyr::mutate(
      age = meta[1],
      tissue = meta[2]
    )
  
  return(tidy_df)
}

##Loop over all sheets in the excel file (so all tissues)
tidy_all <- purrr::imap_dfr(elife_data_list, tidy_one_sheet)

print(tidy_all)


tidy_all <- tidy_all %>%
  dplyr::mutate(
    focus_tissue = case_when(
      ##kidney
      tissue %in% c("Kid-Cortex", "Kid-Med") ~ "Kidney",
      
      ##brain
      tissue %in% "Cortex" ~ "Brain",
      
      ##eye
      tissue == "Eye" ~ "Eye",
      
      ##stomach
      tissue == "Stomach" ~ "Stomach",
      
      ##pancreas
      tissue == "Pancreas" ~ "Pancreas",
      
      ##lung
      tissue == "Lung" ~ "Lung",
      
      ##liver
      tissue == "Liver" ~ "Liver",
      
      ##fallopian tubes (only for females)
      tissue == "Test Fallop tube" & sex == "Female" ~ "Fallopian tube",
      
      ##skin
      tissue == "Skin" ~ "Skin",
      
      ##colon 
      tissue %in% c("Proximal colon", "Distal Colon") ~ "Colon",
      
      #everything else
      TRUE ~ "none"
    )
  )


#if focusing on terminal complex only:
#tidy_all_I <- subset(tidy_all, complex == "CIV")


focus_order <- tidy_all %>%
  dplyr::filter(focus_tissue != "none") %>%
  dplyr::group_by(focus_tissue) %>%
  dplyr::summarise(median_value = median(value, na.rm = TRUE)) %>%
  dplyr::arrange(desc(median_value)) %>%
  dplyr::pull(focus_tissue)

tidy_all %>%
  dplyr::filter(focus_tissue != "none") %>%
  dplyr::mutate(focus_tissue = factor(focus_tissue, levels = focus_order)) %>%
  ggplot(aes(x = focus_tissue, y = log10(value + 0.01), fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.7) +
  geom_jitter(aes(color = sex), width = 0.15, alpha = 0.5, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  ylim(-2.1, -1.35) +
  labs(
    x = "Focus tissue (ordered by median metabolic value)",
    y = "Metabolic value",
    fill = "Sex",
    color = "Sex",
    title = "Metabolic value: log10(Complex Linked Resp. Avg/MTDR across all complexes) 
across focus tissues in young mice (ordered by median)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )








