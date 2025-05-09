library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(zoo)
library(tidyr)

file_path <- "data/NicMetabolomics_Resultss_4.18.23.xlsx"

# Step 1: Read and transpose the first two rows
header_rows <- read_excel(file_path, range = "A1:CX2", col_names = FALSE)
header_t <- t(as.matrix(header_rows))
header_df <- as.data.frame(header_t, stringsAsFactors = FALSE)
colnames(header_df) <- c("group", "label")

# ✅ Step 2: Forward-fill group before dropping phantom columns
header_df <- header_df %>%
  mutate(group = zoo::na.locf(group, na.rm = FALSE),
         label = replace_na(label, ""),
         combined = str_trim(paste(group, label)))

# ✅ Step 3: Drop the first two columns (phantom columns A and B)
header_df <- header_df[-c(1, 2), ]

# Step 4: Read in data
data <- read_excel(file_path, skip = 2, col_names = FALSE)

# Step 5: Align and apply headers
header_df <- header_df[seq_len(ncol(data)), ]
colnames(data) <- header_df$combined

# ✅ Final check
glimpse(data)

# === STEP 2: Keep only Sample Name, Sample Type, and 'Area' columns ===
filtered_data <- data %>%
  select(
    matches("(?i)sample.*name"),
    matches("(?i)sample.*type"),
    contains("Area")
  ) %>%
  filter(str_to_lower(`Sample Type`) %in% c("cal", "sample")) %>%
  select(-matches("(?i)sample.*type")) %>%
  column_to_rownames(var="Sample Name")

head(filtered_data)

# === STEP 6: Compute ratio for each metabolite ===
# Identify all area columns
area_cols <- colnames(filtered_data)
# Extract metabolite names and abbreviations
# Updated metabolite info: strip trailing `-d3` from abbreviations for pairing
metabolite_matches <- str_match(area_cols, "^(.+?)\\s*\\(([-A-Za-z0-9]+)\\).*Results Area")

metabolite_info <- tibble(
  full_name = area_cols,
  name      = metabolite_matches[, 2],
  abbrev    = metabolite_matches[, 3]
) %>%
  filter(!is.na(abbrev)) %>%
  mutate(
    base_abbrev = str_replace(abbrev, "-d3$", ""),   # normalize e.g. CO-d3 -> CO
    is_d3 = str_detect(abbrev, "-d3")
  )

# Compute ratios
ratio_df <- purrr::map_dfc(unique(metabolite_info$base_abbrev), function(base) {
  rows <- metabolite_info %>% filter(base_abbrev == base)

  d3_row     <- rows %>% filter(is_d3)
  non_d3_row <- rows %>% filter(!is_d3)

  if (nrow(d3_row) == 1 && nrow(non_d3_row) == 1) {
    num <- filtered_data[[non_d3_row$full_name]]
    den <- filtered_data[[d3_row$full_name]]
    tibble(!!base := num / den)
  } else {
    message("Skipping: ", base, " (non_d3 = ", nrow(non_d3_row), ", d3 = ", nrow(d3_row), ")")
    NULL
  }
})

# Add sample names back as rownames
rownames(ratio_df) <- rownames(filtered_data)

# View result
glimpse(final_df)