##***************************
## Course: BINF6210 
##
## Assignment 01
##
## by: Liona Vu 
##
## 13 October 2025
##***************************

## Packages that will be used for this analysis
## Loading in the packages required for this analysis.
library("dplyr")
library("tidyverse")
library("readr")
library("janitor")
library("tidyr")
library("vegan")
library("ggplot2")
library("ggpubr") #for adding p values to ggplot2 graphs
theme_set(theme_light())

## Can also use BOLD API to read in the data 
## df_tardigrada <- read_delim("https://www.boldsystems.org/index.php/API_Public/combined?taxon=Tardigrada&format=tsv")

## Importing the data file
df_tardigrada <- read.delim("../data/tardigrada.tsv", sep = "\t") 

## Checking of the data
glimpse(df_tardigrada)
class(df_tardigrada) #is a dataframe

## Cleaning up header names to snake cases.
df_tardigrada <- df_tardigrada %>% 
  janitor::clean_names()

## Exploratory data analysis to gain insight into the tardigrade BOLD dataset
sort(table(df_tardigrada$country_ocean), decreasing = TRUE) #1547 rows of unrecoverable
## and 889 empty
sort(table(df_tardigrada$bin_uri), decreasing = TRUE) #2647 rows with no BOLD ID
nrow(df_tardigrada)

#simple histogram to see top 10 country counts
barplot(head(table(df_tardigrada$country_ocean), n = 10))

## Creating continent vectors to define which countries belong to which continent
## based on the countries listed in the database
north_america <- c("Canada", "United States", "Mexico")
south_america <- c("Chile", "Argentina", "Brazil", "French Guiana", "Ecuador", "Colombia")
antartica <- "Antarctica" #decide to include since it has so many BIN counts

## Creating another continent column and categorizing which countries belong to 
## which regions. 
df_tardigrada <- df_tardigrada %>% 
  mutate(continent = case_when(
    country_ocean %in% north_america ~ "North America",
    country_ocean %in% south_america ~ "South America",
    country_ocean %in% antartica ~ "Antarctica",
    TRUE ~ "Other"))

## Confirming the newly made continent column
glimpse(df_tardigrada)
head(df_tardigrada$continent)

## Since some cells in the bin_uri column are empty, adding NA to the empty cells
df_tardigrada$bin_uri[df_tardigrada$bin_uri == ""] <- NA

## Selecting only columns of interest for the analysis
df_tardigrada_simple <- df_tardigrada %>% 
  select(bin_uri, country_ocean, continent)

##Checking whether the df_tardigrada_simple was properly made
head(df_tardigrada_simple)

## Filtering out NAs from the dataset in continent and bin uri
df_tardigrada_filtered <- df_tardigrada_simple %>%
  filter(continent %in% c("North America", "South America", "Antarctica")) %>% 
  filter(!is.na(bin_uri))

##Checking whether the filtering worked
head(df_tardigrada_filtered)

## Checking if Na values are removed from bin_uri column
sum(is.na(df_tardigrada_filtered$bin_uri)) 
## Returns value of 0, no NA values 

## Counting how many unique bins for each continent
df_count_bins <- df_tardigrada_filtered %>% 
  group_by(continent) %>% 
  summarize(record_count = n_distinct(bin_uri))

df_count_bins

## Plotting the unique BINS for each continent in a bar graph
ggplot(data = df_count_bins) +
  geom_col(mapping = aes(x = continent, y = record_count, fill = continent),
           show.legend = FALSE) + #removes default legend
  theme(panel.grid = element_blank()) + #removes grid lines
  labs(title = "Unique Number of BINs in Different Continents",
       x = "Continent", 
       y = " Number of Unique BINs")

## Saves the figure to the figures file
ggsave("../figs/tardigrada_unique_bins.png", width = 5, height = 5, dpi = 320)

## Next, how complete are the samples collection in all 3 regions?
## Determining the number of unique BINS per continent
df_tardigrada_count <- df_tardigrada_filtered %>% 
  group_by(continent, bin_uri) %>% 
  count(bin_uri)

## Reshaping the data to be able to graph, added values fill = 0 parameter to 
## replace Na values to 0.
df_tardigrada_wide <- df_tardigrada_count %>%
  pivot_wider(names_from = bin_uri, values_from = n, values_fill = 0)

## Checking class of df_tardigrada_wide
class(df_tardigrada_wide)

## Because df_tardigrada_wide is also a tibble, will need 
## to convert to dataframe to retain the row names. Else, warning message 
## 'Setting row names on a tibble is deprecated.' will appear on console
df_tardigrada_wide <- data.frame(df_tardigrada_wide)

## Retaining the row names to label rarecurve plot later
rownames(df_tardigrada_wide) <- df_tardigrada_wide$continent

## Dropping the first column which contains the continent names because the 
## rarecurve function only takes in numerical values
df_tardigrada_wide <- df_tardigrada_wide[,-1]

## Checking the dataframe to see if it was converted properly
glimpse(df_tardigrada_wide)

## Mapping rarecurve to show rarecuve on plot screen
rarecurve(df_tardigrada_wide, 
          col = c("steelblue", "violet", "green"), #adding colour
          main = "Rarefraction curve of Different Continents",
          xlab = "Individuals Barcoded", 
          ylab = "BIN Richness")

## Setting the directory to save a copy
png("../figs/tardigrada_rarecurve.png", width = 800, height = 800, res = 150)

## repeat the same code as above to save plot to png device
rarecurve(df_tardigrada_wide, 
                          col = c("steelblue", "violet", "green"), #adding colour
                          main = "Rarefraction curve of Different Continents",
                          xlab = "Individuals Barcoded", 
                          ylab = "BIN Richness")

## optional: may need to run the following code below twice to make figure disappear
dev.off()

## Calculating the diversity using Shannon's Index of Diversity for each country
## in each continent using the diversity function in vegan package
df_shannon_index <- df_tardigrada_filtered %>%  
  group_by(country_ocean, continent) %>% 
  reframe(shannon = diversity(table(bin_uri)))

df_shannon_index %>% 
  group_by(continent) %>% 
  summarise(mean(shannon))

## Removing Antartica from analysis due to single observation, 1.88
df_shannon_index <- df_shannon_index[-1,]
df_shannon_index

## Checking whether the data is normally distributed by performing the 
## Shapiro-Wilk test
shapiro.test(df_shannon_index$shannon)
## p value = 0.5005, data is not significantly
## different from normal distribution therefore, normality is assumed.

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

rm(list = ls())
