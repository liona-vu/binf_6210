##***************************
## Course: BINF6210 
##
## Assignment 01
##
## by: Liona Vu 
## Contributors: Arwa Sheheryar (Editor) 
## 13 October 2025

## Background:Tardigrades are tiny, resilient animals found worldwide.Studying their biodiversity helps us understand global patterns of species diversity and adaptation. This analysis uses BOLD database records to compare diversity and sampling across continents.

## Research Question:How do tardigrade diversity and sampling completeness differ across North America, South America, and Antarctica using data from the BOLD database?

##===========================================================================

#Package set-up 
required_packages <- c(
  "tidyverse",   # includes dplyr, ggplot2, tidyr, readr, etc.
  "janitor",     # data cleaning (e.g., clean_names)
  "vegan",       # biodiversity analysis
  "ggpubr"       # adding p-values and annotations to ggplots
)
## Install any missing packages automatically
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
## Visual set-up
theme_set(theme_light())

## Lock packages versions for reproducibility 
sessionInfo()

## ===========================================================================
# 0. Data Cleaning and Standardization to Prepare the Tardigrada Dataset for Analysis
## ===========================================================================

## Can also use BOLD API to read in the data 
## df_tardigrada <- read_delim("https://www.boldsystems.org/index.php/API_Public/combined?taxon=Tardigrada&format=tsv")

## Import, clean names, and fix empty strings as NA in one pipe
df_raw <- read_tsv(
  "../data/tardigrada.tsv",
  col_types = cols(.default = col_character()), #forced everything as character to avoid reading issues
  progress = FALSE, show_col_types = FALSE
)

glimpse(df_raw)

#Should not have any parsing warnings however can be investigated by uncommenting below
#readr::problems(raw)

## Data cleaning to snake case, removing weird characters/spaces

df_tardigrada <- df_raw |>
  janitor::clean_names() |>
  mutate( #replacing empty strings with NA values
    bin_uri       = na_if(bin_uri, ""),
    country_ocean = na_if(country_ocean, "")
  )

## CHECK: Should be TRUE if names were cleaned
identical(
  names(df_tardigrada),
  janitor::make_clean_names(names(df_raw))
)

## CHECK: These should both be 0 now, confirming if empty strings are converted to NA
sum(df_tardigrada$bin_uri == "", na.rm = TRUE)
sum(df_tardigrada$country_ocean == "", na.rm = TRUE)

## Take a look at the data
glimpse(df_tardigrada)
class(df_tardigrada) #is a dataframe

## Exploratory data analysis to gain insight into the tardigrade BOLD dataset
sort(table(df_tardigrada$country_ocean), decreasing = TRUE) # 1547 rows of unrecoverable and 889 empty
sort(table(df_tardigrada$bin_uri), decreasing = TRUE) # 2647 rows with no BOLD ID

nrow(df_tardigrada) #5492

## Simple histogram to see top 10 country counts
barplot(head(table(df_tardigrada$country_ocean), n = 10))

## ===========================================================================
# 1. Building a Reproducible Continent Lookup Table to Clean the Dataset and enable Regional Diversity Comparisons    
## ===========================================================================

## Creating continents as a variable/new column to define which countries belong to which continent
## Reusable country–continent lookup table (easier to add more countries and expand dataset later, good for reproducibility) 

continent_lu <- tibble::tibble(
  country_ocean = c(
    "Canada", "United States", "Mexico",
    "Chile", "Argentina", "Brazil", "French Guiana", "Ecuador", "Colombia",
    "Antarctica"
  ),
  continent = c(
    rep("North America", 3),
    rep("South America", 6),
    "Antarctica"
  )
)

## Join the new country/continent table to main dataset  
df_tardigrada <- df_tardigrada %>%
  dplyr::left_join(continent_lu, by = "country_ocean") %>%
  dplyr::mutate(continent = dplyr::coalesce(continent, "Other")) #replace the NAs as "Other"

## Confirming the newly made continent column
glimpse(df_tardigrada)
head(df_tardigrada$continent)

#Check counts by continent (4573 are other, 681 are Antartica, 164 NA and 74 SA)
df_tardigrada %>% count(continent, sort = TRUE)

#See what became other incase you want to explore more countries
df_tardigrada %>%
distinct(country_ocean) %>%
  anti_join(continent_lu, by = "country_ocean") %>%
  arrange(country_ocean) %>%
  head(30)

## Filter to the target groups
df_tardigrada <- df_tardigrada %>%
  dplyr::filter(continent %in% c("North America","South America","Antarctica"))

## Selecting only columns of interest for the analysis
df_tardigrada_simple <- df_tardigrada %>% 
  select(bin_uri, country_ocean, continent)

## Checking whether the df_tardigrada_simple was properly made 
glimpse(df_tardigrada_simple)

## Filtering out NAs from the dataset in continent and bin uri
df_tardigrada_filtered <- df_tardigrada_simple %>%
  filter(!is.na(bin_uri))

##Checking whether the filtering worked
glimpse(df_tardigrada_filtered)

## Checking if Na values are removed from bin_uri column
sum(is.na(df_tardigrada_filtered$bin_uri))  # Returns value of 0, no NA values 


## Count and plot unique BINs by continent using summarise(.by=)
## Reorders bars by BIN count and flips coordinates for easier comparison
df_count_bins <- df_tardigrada_filtered |>
  dplyr::summarise(record_count = dplyr::n_distinct(bin_uri), .by = continent)

## Create and store the plot

plot_bins <- ggplot(df_count_bins,
                    aes(x = reorder(continent, record_count), y = record_count, fill = continent)) +
  geom_col(show.legend = FALSE) +
  theme(panel.grid = element_blank()) +
  labs(title = "Unique BINs by Continent",
       x = "Continent", y = "Unique BINs") +
  coord_flip()

plot_bins

## Save the plot
ggsave("../figs/tardigrada_unique_bins.png", plot = plot_bins, width = 5, height = 5, dpi = 320)

## ===========================================================================
# 2.Rarefaction Workflow: Building a Community Matrix and Estimating Comparable BIN Richness Across Continents
## ===========================================================================

## Build counts matrix in one step (rows = continent, cols = BINs)
comm <- xtabs(~ continent + bin_uri, data = df_tardigrada_filtered)

## Drop any all-zero rows (just in case)
comm <- comm[rowSums(comm) > 0, , drop = FALSE]

## This is a helper to create creates up to 200 evenly spaced integer sample sizes (1..N) for rarefaction. It will help to avoid uneven comparisons with sample sizes and get smooth curves.
sizes_fun <- function(n) unique(round(seq(1, n, length.out = min(n, 200))))

## compute rarefaction curves (class-style lapply + rbind)
rare_list <- lapply(seq_len(nrow(comm)), function(i) {
  counts <- as.numeric(comm[i, ])
  n_tot  <- sum(counts)
  if (n_tot == 0) return(NULL)
  sizes <- sizes_fun(n_tot)
  es    <- as.numeric(vegan::rarefy(counts, sample = sizes))
  data.frame(
    continent   = rownames(comm)[i],
    sample_size = sizes,
    richness    = es,
    row.names   = NULL
  )
})
rare_tbl <- do.call(rbind, rare_list)

## optional endpoints at observed richness
rare_end <- data.frame(
  continent   = rownames(comm),
  sample_size = rowSums(comm),
  richness    = apply(comm, 1, function(x) vegan::specnumber(as.numeric(x))),
  row.names   = NULL
)

## plot
plot_rarefaction <- ggplot(rare_tbl, aes(sample_size, richness, color = continent)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = rare_end, aes(sample_size, richness, color = continent),
             size = 2, show.legend = FALSE) +
  labs(
    title = "Rarefaction Curves of BIN Richness by Continent",
    x = "Individuals Barcoded (sample size)",
    y = "Expected BIN Richness (E[S])",
    color = "Continent"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank())

print(plot_rarefaction)


# Save alongside your other figures
if (!dir.exists("../figs")) dir.create("../figs", recursive = TRUE)
ggsave("../figs/tardigrada_rarefaction_ggplot.png", plot = plot_rarefaction,
       width = 6, height = 5, dpi = 320)

## ===========================================================================
## 3. Computing Shannon Diversity, Assessing Statistical Significance, and Producing Comparative Boxplots
## ===========================================================================

## Calculating the diversity using Shannon's Index of Diversity for each country in each continent using the diversity function in vegan package
df_shannon_index <- df_tardigrada_filtered %>%  
  group_by(country_ocean, continent) %>% 
  reframe(shannon = diversity(table(bin_uri)))

df_shannon_index %>% 
  group_by(continent) %>% 
  summarise(mean(shannon))

## Removing Antartica from analysis due to single observation, 1.88
df_shannon_index <- df_shannon_index[-1,]
df_shannon_index

## Checking whether the data is normally distributed by performing the Shapiro-Wilk test
shapiro.test(df_shannon_index$shannon)
## p value = 0.5005, data is not significantly different from normal distribution therefore, normality is assumed.

## Calculate significance 
t.test(shannon ~ continent, data = df_shannon_index) #not significant, p = 0.38

## Graphing with ggplot
ggplot(data = df_shannon_index, 
       mapping = aes(x = continent, y = shannon, fill = continent)) + 
  geom_boxplot(show.legend = FALSE) +
## Adding the points of each countries' diversity index to the corresponding continent
  geom_jitter(width = 0.1, size = 2, show.legend = FALSE) +
  theme(panel.grid = element_blank())+ #removes grid lines
  stat_compare_means(method = "t.test",
                   label.x = 1.4) + #shifts the p value more to the right
  labs(title = "Shannon's Index between Continents",
           x = "Continents", 
           y = "Shannon's Index (H')")

## Saves the figure to the figures file
ggsave("../figs/tardigrada_boxplot_shannon.png", width = 5, height = 5, dpi = 320)

