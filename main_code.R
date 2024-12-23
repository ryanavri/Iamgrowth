## Preparation ----
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(plotly)

## Simulate dataset ----
# Set simulation parameters
set.seed(212) # For reproducibility
n_plots <- 5 # Number of plots
n_trees_per_plot <- 50 # Number of trees per plot
species <- c("Species_A", "Species_B", "Species_C", "Species_D", "Species_E")

# Simulate initial dataset
simulated_data <- data.frame(
  plot_id = rep(1:n_plots, each = n_trees_per_plot),
  tree_id = paste0("Tree", sprintf("%04d", 1:(n_plots * n_trees_per_plot))), # Unique tree IDs across all plots
  species = sample(species, n_plots * n_trees_per_plot, replace = TRUE),
  DBH = round(rnorm(n_plots * n_trees_per_plot, mean = 20, sd = 5), 1), # Simulated DBH in cm
  height = round(rnorm(n_plots * n_trees_per_plot, mean = 15, sd = 3), 1), # Simulated height in m
  status = "alive" # All trees start as alive
)

# Ensure no negative values for DBH or height
simulated_data <- simulated_data %>%
  mutate(DBH = pmax(DBH, 0),
         height = pmax(height, 0))

# Simulate changes over three visits
visits <- list()
for (visit in 1:5) {
  if (visit == 1) {
    # First visit: Baseline, all trees are alive
    simulated_data$visit <- visit
    visits[[visit]] <- simulated_data
  } else {
    # Simulate some trees dying in subsequent visits
    simulated_data <- simulated_data %>%
      mutate(status = ifelse(runif(n()) < 0.05 * (visit - 1), "dead", status)) # around 5% are dead on subsequent year
    
    # Simulate growth for living trees
    simulated_data <- simulated_data %>%
      mutate(DBH = ifelse(status == "alive", DBH + rnorm(n(), mean = 1, sd = 0.5), DBH), # grow 1 cm in average for subsequent year
             height = ifelse(status == "alive", height + rnorm(n(), mean = 0.3, sd = 0.1), height)) # grow 0.3 m in average for subsequent year
    
    # Ensure no negative values for DBH or height
    simulated_data <- simulated_data %>%
      mutate(DBH = ifelse(DBH < 0, abs(DBH), DBH),
             height = ifelse(height < 0, abs(height), height))
    
    # Add visit identifier
    simulated_data$visit <- visit
    
    # Store data for the visit
    visits[[visit]] <- simulated_data
  }
}

# Combine visits into longitudinal dataset
longitudinal_data <- bind_rows(visits)

## Survival analysis ----
survival_summary <- longitudinal_data %>%
  group_by(plot_id, tree_id, species) %>%
  summarize(survival_status = ifelse(any(status == "dead"), "dead", "alive"),
            total_visits = n(), .groups = "drop") %>%
  group_by(plot_id, species) %>%
  summarize(survival_rate = mean(survival_status == "alive") * 100,
            total_trees = n(), .groups = "drop")

# Assess survivability by species
species_survival <- longitudinal_data %>%
  group_by(tree_id, species) %>%
  summarize(survival_status = ifelse(any(status == "dead"), "dead", "alive"), .groups = "drop") %>%
  group_by(species) %>%
  summarize(survival_rate = mean(survival_status == "alive") * 100,
            total_trees = n(), .groups = "drop")

# Plot survival rates by species
ggplot(species_survival, aes(x = species, y = survival_rate, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = round(survival_rate, 1)), vjust = -0.5) +
  labs(title = "Survival Rates by Species", x = "Species", y = "Survival Rate (%)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

# Statistical test: Compare survival rates between species
survival_test_data <- longitudinal_data %>%
  group_by(tree_id, species) %>%
  summarize(survival_status = ifelse(any(status == "dead"), 0, 1), .groups = "drop")

survival_model <- glm(survival_status ~ species, data = survival_test_data, family = binomial)
summary(survival_model)

## Tree Growth Analysis----
growth_data <- longitudinal_data %>%
  filter(status == "alive") %>%
  group_by(tree_id, plot_id, species) %>%
  mutate(year = visit - min(visit)) %>%
  ungroup()

# Fit mixed-effects models for growth
dbh_growth_model <- lmer(DBH ~ year * species + (1 | plot_id), data = growth_data)
height_growth_model <- lmer(height ~ year * species + (1 | plot_id), data = growth_data)

# Summaries of growth models
summary(dbh_growth_model)
summary(height_growth_model)

# Visualize DBH growth
ggplot(growth_data, aes(x = year, y = DBH, color = species)) +
  geom_smooth(aes(group = species), method = "lm", se = TRUE) +
  labs(title = "DBH Growth Over Time by Species", x = "Year", y = "DBH (cm)", color = "Species") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Visualize DBH growth using plotly
growth_plot <- ggplot(growth_data, aes(x = year, y = DBH, color = species)) +
  geom_smooth(aes(group = species), method = "lm", se = TRUE) +
  labs(title = "DBH Growth Over Time by Species", x = "Year", y = "DBH (cm)", color = "Species") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggplotly(growth_plot)

# Visualize height growth
ggplot(growth_data, aes(x = year, y = height, color = species)) +
  geom_smooth(aes(group = species), method = "lm", se = TRUE) +
  labs(title = "Height Growth Over Time by Species", x = "Year", y = "Height (m)", color = "Species") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## Simulate carbon stock gain overtime----
# Define mean wood density (meanWD) for each species
meanWD_values <- c(
  Species_A = 0.65,
  Species_B = 0.60,
  Species_C = 0.72,
  Species_D = 0.55,
  Species_E = 0.68
)

# Add meanWD to growth_data based on species
growth_data <- growth_data %>%
  mutate(meanWD = recode(species, !!!meanWD_values)) # Assign meanWD based on species

# Calculate carbon stock using the allometric equation (kg/tree)
growth_data <- growth_data %>%
  mutate(carbon_stock_kg = 0.47 * (0.151 * DBH^2.560 * meanWD^0.889))

# Adjust carbon stock to ton/ha
growth_data <- growth_data %>%
  mutate(carbon_stock_ton_ha = (4 * carbon_stock_kg) / 1000) # Convert to ton/ha

# Summarize carbon stock trend over visits
carbon_trend <- growth_data %>%
  group_by(visit) %>%
  summarize(
    total_carbon_stock_ton_ha = sum(carbon_stock_ton_ha, na.rm = TRUE),
    mean_carbon_stock_ton_ha = mean(carbon_stock_ton_ha, na.rm = TRUE),
    sd_carbon_stock_ton_ha = sd(carbon_stock_ton_ha, na.rm = TRUE),
    se_carbon_stock_ton_ha = sd_carbon_stock_ton_ha / sqrt(n()), # Standard error
    .groups = "drop"
  )

ggplot(carbon_trend, aes(x = factor(visit), y = mean_carbon_stock_ton_ha)) +
  geom_bar(stat = "identity", fill = "#00994C", alpha = 0.8) +
  geom_errorbar(aes(
    ymin = mean_carbon_stock_ton_ha - se_carbon_stock_ton_ha,
    ymax = mean_carbon_stock_ton_ha + se_carbon_stock_ton_ha
  ), width = 0.2, color = "#00994C") +
  labs(
    title = "Mean Carbon Stock Trend Overtime",
    x = "Year",
    y = "Mean Carbon Stock (ton/ha)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
