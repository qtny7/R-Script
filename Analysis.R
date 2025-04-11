#Required Packages####
library(geomorph)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(car)
library(FactoMineR)
library(factoextra)
library(nFactors)
library(emmeans)
library(candisc)
library(reshape2)
library(plotly)
library(stats)
library(MASS)

#GPA####
file_list <- list.files("~/Desktop/Dissertation/Data/Scaled Specimen", pattern = "\\.tps$", full.names = TRUE)

# Extract file names (without path or extension) as specimen IDs
specimen_ids <- tools::file_path_sans_ext(basename(file_list))

# Read all the TPS files without auto-ID
tps_list <- lapply(file_list, function(file) {
  readland.tps(file, specID = "None")
})

# landmarks x 2D coordinates x specimens
master_array <- abind::abind(tps_list, along = 3)

# File names as specimen IDs
dimnames(master_array)[[3]] <- specimen_ids

# Function to create a master datasheet (TPS) with all landmark files
write_tps <- function(landmark_array, file_name, specimen_ids) {
  file_conn <- file(file_name, "w")
  
  for (i in 1:dim(landmark_array)[3]) {
    writeLines(paste("LM=", dim(landmark_array)[1], sep=""), file_conn)
    
    # Write landmark coordinates
    writeLines(apply(landmark_array[,,i], 1, paste, collapse=" "), file_conn)
    
    # Use file name as Specimen ID
    writeLines(paste("ID=", specimen_ids[i], sep=""), file_conn)
    
    # Blank line between specimens
    writeLines("", file_conn)
  }
  
  close(file_conn)
}

# save into a new TPS file (the "master datasheet")
write_tps(master_array, "master_Landmarks.tps", specimen_ids)

# reading the master datasheet
land <- readland.tps("master_landmarks.tps", specID = "ID")

# GPA test
GPA_results <- gpagen(land)

# run to see the means
GPA_results

# Extract the aligned landmark coordinates
aligned_data <- GPA_results$coords

# dimensions?
num_landmarks <- dim(aligned_data)[1] 
num_dims <- dim(aligned_data)[2]      
num_specimens <- dim(aligned_data)[3]

gpa_df <- data.frame(Specimen = specimen_ids)

# landmark coordinates for each specimen?
for (i in 1:num_landmarks) {
  gpa_df[[paste0("X_LM", i)]] <- aligned_data[i, 1, ]  
  gpa_df[[paste0("Y_LM", i)]] <- aligned_data[i, 2, ]
}

# add centroid size
gpa_df$Centroid_Size <- GPA_results$Csize

write.csv(gpa_df,"~/Desktop/Dissertation/Analysis/GPA_wideframe.csv", row.names = FALSE)

# Load Data
data <- read.csv("GPA_wideframe.csv")

# Extract Replicate_Population, Replicate_Wing, Sex, and Selection Regime for each specimen and add as new columns

data <- data %>%
  mutate(
    Replicate_Population = case_when(
      str_detect(Specimen, "^LF|^LM") ~ "33",  # LHM is now Replicate 33
      TRUE ~ str_extract(Specimen, "^\\d+")  # Extract number before letter for others (rep pop)
    ),
    
    Replicate_Wing = str_extract(Specimen, "(?<=\\D)(\\d+)"),  # Extract number after letter (rep wing)
    Sex = ifelse(str_detect(Specimen, "F"), "Female", "Male"),  # Assign sex based on "F" or "M"
    
    Selection_Regime = case_when(
      grepl("^(1F|2F|3F|4F|5F|6F|7F|8F|17F|18F|19F|20F|21F|22F|23F|24F)", Specimen) |
        grepl("^(1M|2M|3M|4M|5M|6M|7M|8M|17M|18M|19M|20M|21M|22M|23M|24M)", Specimen) ~ "Female-Limited",
      
      grepl("^(9F|10F|11F|12F|13F|14F|15F|16F|25F|26F|27F|28F|29F|30F|31F|32F)", Specimen) |
        grepl("^(9M|10M|11M|12M|13M|14M|15M|16M|25M|26M|27M|28M|29M|30M|31M|32M)", Specimen) ~ "Male-Limited",
      
      grepl("(LF|LM)", Specimen) ~ "Founder"
    )
  )

write.csv(data,"~/Desktop/Dissertation/Analysis/Dataset.csv", row.names = FALSE)

#RP18 troubleshooting####
## Sex
male_subset <- subset(data, Sex == "Male") # total number of males: 942
female_subset <- subset(data, Sex == "Female") # total number of females: 888

## Selection_Regime
ML_subset <- subset(data, Selection_Regime == "Male-Limited") # male-limited total: 884
FL_subset <- subset(data, Selection_Regime == "Female-Limited") # female-limited total: 886

## Sex::Selection_Regime
subset1 <- subset(ML_subset, Sex == "Male") # males in male-limited: 455
subset2 <- subset(ML_subset, Sex == "Female") # females in male-limited: 429

subset3 <- subset(FL_subset, Sex == "Male") # males in female-limited: 457
subset4 <- subset(FL_subset, Sex == "Female") # females in female-limited: 429

## Sex::Selection_Regime::Rep_Pop::Rep_Wing (Female-Limited, females)
sub1F <- subset(subset4, Replicate_Population == "1")
sub2F <- subset(subset4, Replicate_Population == "2")
sub3F <- subset(subset4, Replicate_Population == "3")
sub4F <- subset(subset4, Replicate_Population == "4")
sub5F <- subset(subset4, Replicate_Population == "5")
sub6F <- subset(subset4, Replicate_Population == "6")
sub7F <- subset(subset4, Replicate_Population == "7")
sub8F <- subset(subset4, Replicate_Population == "8")
sub17F <- subset(subset4, Replicate_Population == "17")
sub18F <- subset(subset4, Replicate_Population == "18")
sub19F <- subset(subset4, Replicate_Population == "19")
sub20F <- subset(subset4, Replicate_Population == "20")
sub21F <- subset(subset4, Replicate_Population == "21")
sub22F <- subset(subset4, Replicate_Population == "22")
sub23F <- subset(subset4, Replicate_Population == "23")
sub24F <- subset(subset4, Replicate_Population == "24")

## Sex::Selection_Regime::Rep_Pop::Rep_Wing (Female-Limited, males)
sub1M <- subset(subset3, Replicate_Population == "1")
sub2M <- subset(subset3, Replicate_Population == "2")
sub3M <- subset(subset3, Replicate_Population == "3")
sub4M <- subset(subset3, Replicate_Population == "4")
sub5M <- subset(subset3, Replicate_Population == "5")
sub6M <- subset(subset3, Replicate_Population == "6")
sub7M <- subset(subset3, Replicate_Population == "7")
sub8M <- subset(subset3, Replicate_Population == "8")
sub17M <- subset(subset3, Replicate_Population == "17")
sub18M <- subset(subset3, Replicate_Population == "18")
sub19M <- subset(subset3, Replicate_Population == "19")
sub20M <- subset(subset3, Replicate_Population == "20")
sub21M <- subset(subset3, Replicate_Population == "21")
sub22M <- subset(subset3, Replicate_Population == "22")
sub23M <- subset(subset3, Replicate_Population == "23")
sub24M <- subset(subset3, Replicate_Population == "24")


## Sex::Selection_Regime::Rep_Pop::Rep_Wing (Male-Limited, females)
sub9F <- subset(subset2, Replicate_Population == "9")
sub10F <- subset(subset2, Replicate_Population == "10")
sub11F <- subset(subset2, Replicate_Population == "11")
sub12F <- subset(subset2, Replicate_Population == "12")
sub13F <- subset(subset2, Replicate_Population == "13")
sub14F <- subset(subset2, Replicate_Population == "14")
sub15F <- subset(subset2, Replicate_Population == "15")
sub16F <- subset(subset2, Replicate_Population == "16")
sub25F <- subset(subset2, Replicate_Population == "25")
sub26F <- subset(subset2, Replicate_Population == "26")
sub27F <- subset(subset2, Replicate_Population == "27")
sub28F <- subset(subset2, Replicate_Population == "28")
sub29F <- subset(subset2, Replicate_Population == "29")
sub30F <- subset(subset2, Replicate_Population == "30")
sub31F <- subset(subset2, Replicate_Population == "31")
sub32F <- subset(subset2, Replicate_Population == "32")

## Sex::Selection_Regime::Rep_Pop::Rep_Wing (Male-Limited, males)
sub9M <- subset(subset1, Replicate_Population == "9")
sub10M <- subset(subset1, Replicate_Population == "10")
sub11M <- subset(subset1, Replicate_Population == "11")
sub12M <- subset(subset1, Replicate_Population == "12")
sub13M <- subset(subset1, Replicate_Population == "13")
sub14M <- subset(subset1, Replicate_Population == "14")
sub15M <- subset(subset1, Replicate_Population == "15")
sub16M <- subset(subset1, Replicate_Population == "16")
sub25M <- subset(subset1, Replicate_Population == "25")
sub26M <- subset(subset1, Replicate_Population == "26")
sub27M <- subset(subset1, Replicate_Population == "27")
sub28M <- subset(subset1, Replicate_Population == "28")
sub29M <- subset(subset1, Replicate_Population == "29")
sub30M <- subset(subset1, Replicate_Population == "30")
sub31M <- subset(subset1, Replicate_Population == "31")
sub32M <- subset(subset1, Replicate_Population == "32")

# error bars (standard deviation)
summary_data <- data %>%
  group_by(Selection_Regime, Sex) %>%
  summarise(
    mean_centroid = mean(Centroid_Size, na.rm = TRUE),
    sd_centroid = sd(Centroid_Size, na.rm = TRUE),  # Standard deviation
    .groups = 'drop'
  )

# box plots for all centroid size by selection regime
ggplot(data, aes(x = Selection_Regime, y = Centroid_Size, fill = Selection_Regime)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5, aes(color = Selection_Regime)) + 
  geom_errorbar(data = summary_data, aes(x = Selection_Regime, ymin = mean_centroid - sd_centroid, 
                                         ymax = mean_centroid + sd_centroid, y = mean_centroid), 
                width = 0.2, inherit.aes = FALSE) +  # Error bars fix
  scale_fill_manual(values = c("Female-Limited" = "palegreen2", "Founder" = "gray", "Male-Limited" = "plum2")) +
  scale_color_manual(values = c("Female-Limited" = "olivedrab", "Founder" = "black", "Male-Limited" = "purple4")) +
  theme_minimal() +
  labs(title = "Centroid Size Across Selection Regime",
       x = "Selection Regime",
       y = "Centroid Size",
       fill = "Selection_Regime",
       color = "Selection_Regime") +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# Centroid Size by sex and selection regime
group_colors <- c("Female-Limited" = "palegreen2", "Founder" = "gray", "Male-Limited" = "plum2")

ggplot(data, aes(x = Selection_Regime, y = Centroid_Size, fill = Selection_Regime)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  
  
  geom_jitter(
    aes(color = Selection_Regime),  # <== ADD this to map color to a variable
    width = 0.2,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    data = summary_data, 
    aes(x = Selection_Regime, 
        ymin = mean_centroid - sd_centroid, 
        ymax = mean_centroid + sd_centroid, 
        y = mean_centroid),
    color = "red", 
    width = 0.3, 
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = group_colors) +  
  scale_color_manual(values = c(
    "Female-Limited" = "palegreen4",
    "Founder" = "gray50",
    "Male-Limited" = "plum4"
  )) +
  
  theme_minimal() +
  
  labs(
    title = "Centroid Size by Sex and Selection Regime",
    x = "Selection Regime",
    y = "Centroid Size",
    fill = "Selection Regime",
    color = "Selection Regime"
  ) +
  theme(axis.text.x = element_text(size = 12, face = "bold")) +
  facet_wrap(~ Sex) 

#Alternative plot??
group_colors <- c("Female-Limited" = "palegreen2", "Founder" = "gray", "Male-Limited" = "plum2")

ggplot(data, aes(x = Sex, y = Centroid_Size, fill = Sex)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  
  
  geom_jitter(
    aes(color = Sex),  # <== ADD this to map color to a variable
    width = 0.2,
    alpha = 0.6
  ) +
  
  geom_errorbar(
    data = summary_data, 
    aes(x = Sex, 
        ymin = mean_centroid - sd_centroid, 
        ymax = mean_centroid + sd_centroid, 
        y = mean_centroid),
    color = "red", 
    width = 0.3, 
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = group_colors) +  
  scale_color_manual(values = c(
    "Female-Limited" = "palegreen4",
    "Founder" = "gray50",
    "Male-Limited" = "plum4"
  )) +
  
  theme_minimal() +
  
  labs(
    title = "Centroid Size by Sex and Selection Regime",
    x = "Sex",
    y = "Centroid Size",
    fill = "Sex",
    color = "Sex"
  ) +
  theme(axis.text.x = element_text(size = 12, face = "bold")) +
  facet_wrap(~ Selection_Regime) 


# males of female-limited, plot
FL_male <- data %>%
  filter(Sex == "Male", Selection_Regime == "Female-Limited", 
         Replicate_Population %in% c(1:8, 17:24))

ggplot(FL_male, aes(x = factor(Replicate_Population), y = Centroid_Size)) +
  geom_boxplot(aes(fill = factor(Replicate_Population)), alpha = 0.6) +  
  geom_jitter(aes(color = factor(Replicate_Population)), width = 0.2, alpha = 0.8) + 
  labs(x = "Replicate Population", y = "Centroid Size", title = "Male Flies in Female-Limited Selection Regime") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# females of female-limited, plot
FL_female <- data %>%
  filter(Sex == "Female", Selection_Regime == "Female-Limited", 
         Replicate_Population %in% c(1:8, 17:24))

ggplot(FL_female, aes(x = factor(Replicate_Population), y = Centroid_Size)) +
  geom_boxplot(aes(fill = factor(Replicate_Population)), alpha = 0.6) +  
  geom_jitter(aes(color = factor(Replicate_Population)), width = 0.2, alpha = 0.8) + 
  labs(x = "Replicate Population", y = "Centroid Size", title = "Female Flies in Female-Limited Selection Regime") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# males and females of rep pop 18
ggplot(data %>% filter(Replicate_Population == 18), aes(x = factor(Replicate_Population), y = Centroid_Size)) +
  geom_boxplot(aes(fill = factor(Sex)), alpha = 0.6) +  
  geom_jitter(aes(color = factor(Sex)), width = 0.2, alpha = 0.8) + 
  labs(x = "Replicate Population", y = "Centroid Size", title = "Male and Females of Replicate Population 18") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# make sure there are no duplicates/overlaps/etc
ggplot() +
  geom_boxplot(
    data = data %>% filter(Sex == "Female"), 
    mapping = aes(x = Replicate_Population, y = Centroid_Size)
  ) +
  geom_boxplot(
    data = data %>% filter(Sex == "Male") %>% filter(Replicate_Population==18), 
    mapping = aes(x = Replicate_Population, y = Centroid_Size), 
    fill = "black", 
    alpha = 0.5
  ) 

# males of male-limited, plot
ML_male <- data %>%
  filter(Sex == "Male", Selection_Regime == "Male-Limited", 
         Replicate_Population %in% c(9:16, 25:32))

ggplot(ML_male, aes(x = factor(Replicate_Population), y = Centroid_Size)) +
  geom_boxplot(aes(fill = factor(Replicate_Population)), alpha = 0.6) +  
  geom_jitter(aes(color = factor(Replicate_Population)), width = 0.2, alpha = 0.8) + 
  labs(x = "Replicate Population", y = "Centroid Size", title = "Male Flies in Male-Limited Selection Regime") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# females of male-limited, plot
ML_female <- data %>%
  filter(Sex == "Female", Selection_Regime == "Male-Limited", 
         Replicate_Population %in% c(9:16, 25:32))

ggplot(ML_female, aes(x = factor(Replicate_Population), y = Centroid_Size)) +
  geom_boxplot(aes(fill = factor(Replicate_Population)), alpha = 0.6) +  
  geom_jitter(aes(color = factor(Replicate_Population)), width = 0.2, alpha = 0.8) + 
  labs(x = "Replicate Population", y = "Centroid Size", title = "Female Flies in Male-Limited Selection Regime") +
  theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# distance plot (LM5 & 10 replicate pops 17 and 18)
data <- data %>%
  mutate(Distance_LM5_LM10 = sqrt((X_LM10 - X_LM5)^2 + (Y_LM10 - Y_LM5)^2))

ggplot(data, aes(x = Sex, y = Distance_LM5_LM10, fill = Sex)) +
  geom_boxplot(alpha = 0.6) +  # Boxplot for distributions
  geom_jitter(aes(color = Sex), width = 0.2, alpha = 0.8) +  # Individual points
  facet_wrap(~Replicate_Population) +  # Separate by Replicate Population
  labs(x = "Sex", y = "Distance between LM 5 & 10", 
       title = "Distance between LM 5 and 10 in Males and Females of Replicate Populations 17 & 18") + theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

# distance plot (LM7 & 12 replicate pops 17 and 18)
data <- data %>%
  mutate(Distance_LM7_LM12 = sqrt((X_LM12 - X_LM7)^2 + (Y_LM12 - Y_LM7)^2))

ggplot(data, aes(x = Sex, y = Distance_LM7_LM12, fill = Sex)) +
  geom_boxplot(alpha = 0.6) +  # Boxplot for distributions
  geom_jitter(aes(color = Sex), width = 0.2, alpha = 0.8) +  # Individual points
  facet_wrap(~Replicate_Population) +  # Separate by Replicate Population
  labs(x = "Sex", y = "Distance between LM 7 & 12", 
       title = "Distance between LM 7 and 12 in Males and Females of Replicate Populations 17 & 18") + theme_minimal() +
  scale_fill_viridis_d() +
  scale_color_viridis_d()

#ANOVA for Centroid Size####
head(data)

model1 <- lmer(Centroid_Size ~ Selection_Regime + Sex + (1 | Replicate_Population) + (1 | Replicate_Wing), data = data)

Anova(model1)
summary(model1)

#PC Analysis and Visualisation####
df_new <- data %>% dplyr::select(-Centroid_Size) # removing Centroid size from the data being analysed
df_new <- df_new %>% mutate(Group = paste(df_new$Sex, df_new$Selection_Regime, sep = "_"))

# PCA
pca <- PCA(df_new, scale.unit = FALSE, ncp = 24, quali.sup = c(1, 26, 27, 28, 29, 30, 31, 32), graph = FALSE) # runs PCA, specimen IDs are not quantifiable data
summary(pca)

# visualise and extract cumulative variance
fviz_eig(pca, addlabels = TRUE, ncp = 24, ylim = c(0, 25))
eig_values <- pca$eig

# contributions of variables to components
pca$var$contrib
pca$ind$cos2

# plot of variables
var_explained <- cumsum(eig_values[,2]) 
fviz_pca_var(pca, axes = c(4, 6), repel = TRUE) 

# by Sex (1 & 2)
fviz_pca_ind(pca, 
             axes = c(1, 2), 
             geom = "point",  
             col.ind = df_new$Sex,  
             palette = c("palegreen4", "plum3"),  # green for females, purple for males
             addEllipses = TRUE,  
             title = "Scatterplot of Dimensions 1 and 2",  
             legend.title = "Sex")

# by Selection Regime
fviz_pca_ind(pca, 
             axes = c(1, 2), 
             geom = "point",  
             col.ind = df_new$Selection_Regime,  # Colour by Selection Regime
             palette = c("palegreen4", "yellow", "plum2"),  # blue for female-limited, yellow for male-limited, pink for founder 
             addEllipses = TRUE,  
             title = "Scatterplot of Dimensions 1 and 2",  
             legend.title = "Selection Regime")  # Set legend title
 
# by Selection Regime::Sex
## Perform PCA based on grouping
pca1 <- PCA(df_new, scale.unit = FALSE, ncp = 24, quali.sup = c(1, 26, 27, 28, 29, 30, 31, 32), graph = FALSE)  # Exclude non-numeric columns

## everything against everything
fviz_pca_ind(pca1, 
             axes = c(1, 2), 
             geom = "point",  
             col.ind = df_new$Group,  # Color by combined factor
             palette = c("palegreen4", "olivedrab3", "darkolivegreen1", "plum4", "orchid3", "violetred"),
             addEllipses = TRUE,  # Add confidence ellipses
             legend.title = "Sex & Selection Regime")
#ANOVA (with 18)####
## ANOVA for Selection Regime and Sex (with 18, excluding LHm) ##
pca_data <- as.data.frame(pca$ind$coord)  # PCA scores
pca_data$Group <- df_new$Group  # (Sex_Selection_Regime)
pca_data$Sex <- df_new$Sex  # (Sex)
pca_data$Selection_Regime <- df_new$Selection_Regime  # (Selection Regime)
pca_data$Replicate_Population <- df_new$Replicate_Population  # (rep pop)
pca_data$Replicate_Wing <- df_new$Replicate_Wing  # (rep wing)
pca_data_filtered <- pca_data %>% filter(!Selection_Regime %in% c("Founder"))

# ANOVA - mixed effects #

# Store dimensions
dims <- paste0("Dim.", 1:7)

# Empty list to store models, ANOVA summaries, and post-hoc results
models <- list()
anova_summaries <- list()
posthoc_results <- list()

# Loop through dimensions 1 through 7
for (dim in dims) {
  
  # Fit the mixed-effects model for the current dimension
  formula <- as.formula(paste(dim, "~ Selection_Regime * Sex + (1 | Selection_Regime:Replicate_Population) + (1 | Replicate_Wing)"))
  
  model <- lmer(formula, data = pca_data_filtered)
  
  # Save the model object
  models[[dim]] <- model
  
  # Get ANOVA (Type 3)
  options(contrasts = c("contr.sum", "contr.poly"))
  anova_summary <- Anova(model, type = 3)
  
  # Save ANOVA summary
  anova_summaries[[dim]] <- anova_summary
  
  # Print ANOVA summary for reference
  print(paste("ANOVA Type III for", dim))
  print(anova_summary)
  
  # Get estimated marginal means
  emm <- emmeans(model, ~ Selection_Regime * Sex)
  
  # Tukey's HSD post-hoc pairwise comparisons
  posthoc <- pairs(emm, adjust = "tukey")
  
  # Save post-hoc results
  posthoc_results[[dim]] <- posthoc
  
  # Print post-hoc summary for reference
  print(paste("Post-hoc Tukey HSD for", dim))
  print(posthoc)
}

tukey_df <- as.data.frame(posthoc_results)

#MANOVA (with 18)####

# matrix of the PCs as dependent variables
response_matrix <- as.matrix(pca_data_filtered[, c("Dim.1", "Dim.2", "Dim.3", "Dim.4", "Dim.5", "Dim.6", "Dim.7")])

# Selection_Regime and Sex as the fixed effects
manova <- manova(response_matrix ~ Selection_Regime * Sex, data = pca_data_filtered)

# Summary of MANOVA (Pillai’s Trace)
summary(manova)

#All analysis without 18####
df_new_18 <- df_new %>% filter(Replicate_Population != 18) # removing pop 18 (males and females)

pca2 <- PCA(df_new_18, scale.unit = FALSE, ncp = 24, quali.sup = c(1, 26, 27, 28, 29, 30, 31, 32, 33), graph = FALSE) # runs PCA, specimen IDs are not quantifiable data
summary(pca2)

# Visualise and extract cumulative variance
fviz_eig(pca2, addlabels = TRUE, ncp = 24, ylim = c(0, 25))
eig_values <- pca2$eig

# by Sex (1 & 2)
fviz_pca_ind(pca2, 
             axes = c(1, 2), 
             geom = "point",  
             col.ind = df_new_18$Sex,  # Colour by sex
             palette = c("slategray2", "hotpink2"),  # green for females, purple for males
             addEllipses = TRUE,  
             title = "Scatterplot of Dimensions 1 and 2 (RP18 removed)",  
             legend.title = "Sex")  # Set legend title

# by Sex (1 & 3)
fviz_pca_ind(pca2, 
             axes = c(1, 3), 
             geom = "point",  
             col.ind = df_new_18$Sex,  # Colour by sex
             palette = c("slategray2", "hotpink2"),  # green for females, purple for males
             addEllipses = TRUE,  
             title = "Scatterplot of Dimensions 1 and 23 (RP18 removed)",  
             legend.title = "Sex")  # Set legend title

# by Selection Regime
fviz_pca_ind(pca2, 
             axes = c(1, 2), 
             geom = "point",  
             col.ind = df_new$Selection_Regime,  # Colour by Selection Regime
             palette = c("palegreen4", "black", "plum3"),  # blue for female-limited, yellow for male-limited, pink for founder 
             addEllipses = TRUE,  
             title = "Scatterplot of Dimensions 1 and 2",  
             legend.title = "Selection Regime")  # Set legend title

## ANOVA for Selection Regime and Sex (without 18, excluding LHm) ##
pca_data2 <- as.data.frame(pca2$ind$coord)  # PCA scores
pca_data$Group <- df_new_18$Group  # (Sex_Selection_Regime)
pca_data2$Sex <- df_new_18$Sex  # (Sex)
pca_data2$Selection_Regime <- df_new_18$Selection_Regime  # (Selection Regime)
pca_data2$Replicate_Population <- df_new_18$Replicate_Population  # (rep pop)
pca_data2$Replicate_Wing <- df_new_18$Replicate_Wing  # (rep wing)
pca_data_filtered2 <- pca_data2 %>% filter(!Selection_Regime %in% c("Founder"))

# Store dimensions
dims2 <- paste0("Dim.", 1:7)

# Empty list to store models, ANOVA summaries, and post-hoc results
models2 <- list()
anova_summaries2 <- list()
posthoc_results2 <- list()

# Loop through dimensions 1 through 7
for (dim in dims2) {
  
  # Fit the mixed-effects model for the current dimension
  formula2 <- as.formula(paste(dim, "~ Selection_Regime * Sex + (1 | Selection_Regime:Replicate_Population) + (1 | Replicate_Wing)"))
  
  model2 <- lmer(formula2, data = pca_data_filtered2)
  
  # Save the model object
  models[[dim]] <- model2
  
  # Get ANOVA (Type 3)
  options(contrasts = c("contr.sum", "contr.poly"))
  anova_summary2 <- Anova(model2, type = 3)
  
  # Save ANOVA summary
  anova_summaries2[[dim]] <- anova_summary2
  
  # Print ANOVA summary for reference
  print(paste("ANOVA Type III for", dim))
  print(anova_summary2)
  
  # Get estimated marginal means
  emm2 <- emmeans(model2, ~ Selection_Regime * Sex)
  
  # Tukey's HSD post-hoc pairwise comparisons
  posthoc2 <- pairs(emm2, adjust = "tukey")
  
  # Save post-hoc results
  posthoc_results2[[dim]] <- posthoc2
  
  # Print post-hoc summary for reference
  print(paste("Post-hoc Tukey HSD for", dim))
  print(posthoc2)
}

#Tukey's Heatmap####
# Calculate group means across PCs
group_means <- pca_data_filtered %>%
  group_by(Selection_Regime, Sex) %>%
  summarise(across(starts_with("Dim."), mean, na.rm = TRUE))

group_means_melt <- melt(group_means, id.vars = c("Selection_Regime", "Sex"))
group_means_melt <- group_means_melt %>% filter(variable %in% paste0("Dim.", 1:7))
group_means_melt$variable <- factor(group_means_melt$variable, levels = paste0("Dim.", 1:7))

letters_list <- list()

for (dim in paste0("Dim.", 1:7)) {
  emm <- emmeans(models[[dim]], ~ Selection_Regime * Sex)
  
  letters_df <- cld(emm, Letters = letters) %>%
    dplyr::select(Selection_Regime, Sex, .group) %>%
    mutate(variable = dim)
  
  letters_list[[dim]] <- letters_df
}

letters_combined <- bind_rows(letters_list)

group_means_melt <- group_means_melt %>%
  mutate(across(c(Selection_Regime, Sex, variable), as.character))

letters_combined <- letters_combined %>%
  mutate(across(c(Selection_Regime, Sex, variable), as.character))

heatmap_data <- left_join(group_means_melt, letters_combined,
                          by = c("Selection_Regime", "Sex", "variable"))

# Heatmap with Tukey's results as letters
ggplot(heatmap_data, aes(x = variable, y = interaction(Selection_Regime, Sex), fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = .group), color = "black", size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "palegreen4", mid = "white", high = "plum4") +
  labs(
    x = "PC Dimension",
    y = "Sex and Selection Regime",
    fill = "Mean PC Score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )


#Incomplete Crossvein####

## Chi Sq test for the incomplete cross veins 
crossvein_data <- matrix(c(37, 50, 23, 25), 
               nrow = 2, 
               byrow = TRUE)

rownames(crossvein_data) <- c("Female", "Male")
colnames(crossvein_data) <- c("Female_Limited", "Male_Limited")

# Observed counts in different groups
observed_ML <- c(50, 25)
observed_FL <- c(37, 23)
observed_sex <- c(87, 48)
observed_SR <- c(75, 60)
observed_females <- c(37, 50)

# Run Chi-Square Goodness-of-Fit Test
chisq.test(observed_ML)
chisq.test(observed_FL)
chisq.test(observed_sex)
chisq.test(observed_SR)
chisq.test(observed_females)

crossvein_data <- data.frame(
  Sex = c("Female", "Female", "Male", "Male"),
  Selection_Regime = c("Female-Limited", "Male-Limited", "Female-Limited", "Male-Limited"),
  Num_Incomplete_Crossveins = c(37, 50, 23, 25)
)

ggplot(crossvein_data, aes(x = Selection_Regime, y = Num_Incomplete_Crossveins, fill = Sex)) +
  geom_bar(alpha = 0.7, width = 0.5, stat = "identity", position = "dodge") +
  labs(
    title = "Incomplete Crossveins by Sex and Selection Regime",
    x = "Selection Regime",
    y = "Number of Wings with Incomplete Crossveins"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Female" = "palegreen4", "Male" = "plum4")) +
  theme(
    plot.title = element_text(size = 12),       
    axis.title.x = element_text(size = 14),       
    axis.title.y = element_text(size = 14),       
    axis.text.x = element_text(size = 13),      
    axis.text.y = element_text(size = 13)        
  )

head(crossvein_data)
