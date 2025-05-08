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

# Helper: Extract unique metabolite prefixes from 'Area' column names
area_cols <- colnames(filtered_data)
metabolite_prefixes <- str_match(area_cols, "^(.+?)\\s*\\(([-A-Za-z0-9]+)\\).*Results Area")[, 2:3]
metabolite_info <- tibble(full_name = area_cols, name = metabolite_prefixes[, 1], abbrev = metabolite_prefixes[, 2]) %>%
  filter(!is.na(abbrev)) %>%
  distinct(abbrev, .keep_all = TRUE)

# Match d3 (ISTD) and non-d3 pairs
ratio_df <- purrr::map_dfc(metabolite_info$abbrev, function(abbr) {
  # Find matching column names
  non_d3_col <- grep(paste0("(?<!-d3)\\(", abbr, "\\).*Area"), area_cols, value = TRUE)
  d3_col     <- grep(paste0("-d3 \\(", abbr, "-d3\\).*Area"), area_cols, value = TRUE)

  if (length(non_d3_col) == 1 && length(d3_col) == 1) {
    ratio <- filtered_data[[non_d3_col]] / filtered_data[[d3_col]]
    tibble(!!abbr := ratio)
  } else {
    NULL  # Skip if one or both columns are not found
  }
})

# Final dataframe: ratios with rownames as Sample Name
final_df <- ratio_df
rownames(final_df) <- rownames(filtered_data)

# View result
print(head(final_df))