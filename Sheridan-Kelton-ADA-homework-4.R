library(tidyverse)
library(patchwork) 
library(scales) 
library(collapse)
library(BBmisc)
library(devtools)
remotes::install_github("joachim-gassen/tidycovid19")
library(tidycovid19)

merged <- download_merged_data(cached = TRUE)
merged <- merged %>%
  group_by(country, iso3c) %>%
  arrange(country, iso3c, date) %>%
  mutate(
    daily_confirmed = confirmed - lag(confirmed, n = 1),
    daily_deaths = deaths - lag(deaths, n = 1),
    daily_recovered = recovered - lag(recovered, n = 1)
  ) %>%
  mutate(
    daily_confirmed = replace_na(daily_confirmed, 0),
    daily_deaths = replace_na(daily_deaths, 0),
    daily_recovered = replace_na(daily_recovered, 0)
  ) %>%
  ungroup() %>%
  arrange(country, iso3c, date)

add_world1 <- merged %>%
  group_by(date) %>%
  arrange(date) %>%
  summarize(
    country = "World", iso3c = NA,
    confirmed = sum(confirmed, na.rm = TRUE),
    deaths = sum(deaths, na.rm = TRUE),
    recovered = sum(recovered, na.rm = TRUE),
    timestamp = fmode(timestamp)
  ) %>%
  mutate(
    daily_confirmed = confirmed - lag(confirmed, n = 1),
    daily_deaths = deaths - lag(deaths, n = 1),
    daily_recovered = recovered - lag(recovered, n = 1)
  ) %>%
  mutate(
    daily_confirmed = replace_na(daily_confirmed, 0),
    daily_deaths = replace_na(daily_deaths, 0),
    daily_recovered = replace_na(daily_recovered, 0)
  ) %>%
  ungroup() %>%
  arrange(country, iso3c, date)

add_world2 <- merged %>%
  group_by(country, iso3c) %>%
  summarize(
    population = fmode(population),
    land_area_skm = fmode(land_area_skm),
    timestamp = fmode(timestamp)
  ) %>%
  ungroup() %>%
  summarize(
    country = "World", iso3c = NA,
    population = sum(population, na.rm = TRUE),
    land_area_skm = sum(land_area_skm, na.rm = TRUE)
  ) %>%
  mutate(pop_density = population / land_area_skm)

add_world <- left_join(add_world1, add_world2, by = c("country", "iso3c"))
merged <- bind_rows(merged, add_world)

cv_data <- pivot_longer(merged,
                        cols = c(
                          "confirmed", "deaths", "recovered",
                          "daily_confirmed", "daily_deaths", "daily_recovered"
                        ),
                        names_to = "variable", values_to = "cases"
) %>%
  arrange(country, variable, date) %>%
  rename(area = land_area_skm, density = pop_density) %>%
  mutate(rate = cases / population * 10^6)

cv_summary <- function(d, country_list = "World",
                                  plot = TRUE, facet = "country",
                                  status = c("confirmed", "deaths", "recovered")) {
  
  # based on `wes_palettes()` GrandBudapest1, IsleofDogs1, IsleofDogs2
  # from the {wesanderson} package
  my_palette <- c(
    "#5B1A18", "#FD6467", "#F1BB7B", "#D67236",
    "#0F0D0E", "#9986A5", "#79402E", "#CCBA72", "#D9D0D3", "#8D8680",
    "#EAD3BF", "#AA9486", "#B6854D", "#39312F", "#1C1718"
  )
  
  if (facet == "country") {
    fill <- "variable"
    n <- length(unique(d$variable)) / 2
    # need only half of unique # of variable (3)
  }
  
  if (facet == "variable") {
    fill <- "country"
    n <- length(country_list)
    # need number of countries
  }
  
  if ("All" %in% country_list) {
    country_list <- unique(d$country)
    country_list <- setdiff(country_list, "World")
  }
  
  if ("World" %in% country_list) {
    d <- d %>% filter(country %in% country_list)
    
    totals <- d %>%
      group_by(variable) %>%
      summarize(
        country = "World",
        cases = max(cases),
        population = max(population),
        area = max(area),
        density = max(density),
        rate = max(rate, na.rm = TRUE),
        on = max(date)
      ) %>%
      select(country, variable, cases, population, area, density, rate, on) %>%
      arrange(variable) %>%
      ungroup()
  }
  
  if ("World" %nin% country_list) {
    d <- d %>% filter(country %in% country_list)
    totals <- d %>%
      group_by(country, variable) %>%
      summarize(
        cases = max(cases),
        population = max(population),
        area = max(area),
        density = max(density),
        rate = max(rate, na.rm = TRUE),
        on = max(date),
        gdp_capita = fmode(gdp_capita),
        income = fmode(income),
        life_expectancy = fmode(life_expectancy),
        max_sd = max(soc_dist),
        max_mr = max(mov_rest)
      ) %>%
      select(
        country, variable, cases, population, area, density, rate,
        gdp_capita, income, life_expectancy, max_sd, max_mr, on
      ) %>%
      arrange(country, variable) %>%
      ungroup()
  }
  
  if (plot == TRUE) {
    cc <- filter(d, variable %in% status)
    cum_cases_plot <- ggplot(
      data = cc,
      # use the tidy evaluation pronoun .data to slice the chosen fill
      # variable from the data frame
      aes(
        x = date, y = cases + 1, color = .data[[fill]],
        fill = .data[[fill]]
      )
    ) +
      geom_point(size = 0.5) +
      geom_line() +
      # use the tidy evaluation pronoun .data to slice the chosen facet_wrap
      # variable from the data frame
      facet_wrap(~ .data[[facet]], ncol = 5) +
      xlab("Date") +
      ylab("Log Cumulative Cases") +
      scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))
      ) +
      scale_color_manual(
        aesthetics = c("color", "fill"),
        name = NULL, values = my_palette
      )
    
    dc <- filter(d, variable %in% paste0("daily_", status))
    daily_cases_plot <- ggplot(
      data = dc,
      aes(
        x = date, y = cases, color = .data[[fill]],
        fill = .data[[fill]]
      )
    ) +
      geom_point(size = 0.5) +
      geom_line() +
      facet_wrap(~ .data[[facet]], ncol = 5) +
      xlab("Date") +
      ylab("Daily Cases") +
      scale_color_manual(
        aesthetics = c("color", "fill"),
        name = NULL, values = my_palette
      )
  }
  
  if (plot == TRUE) {
    return(list(
      totals = totals,
      cum_cases_plot = cum_cases_plot,
      daily_cases_plot = daily_cases_plot
    ))
  } else {
    return(list(totals = totals))
  }
}


#Challenge 1
#Use the dataset and function generated above to plot global data 
#on confirmed coronavirus infections, deaths, and recoveries.

cv_summary(cv_data)    

#Challenge 2
#Use the dataset and function generated above to plot data on confirmed 
#coronavirus infections, deaths, and recoveries for the “Group of Seven” 
#(G7) countries, which are the largest IMF-advanced economies in the world 
#(i.e., the US, United Kingdom, Canada, France, Germany, Italy, and Japan) 
#plus China, Russia, and Iran. Facet your plots first by “country” and then 
#by “variable”.

cv_summary(cv_data, facet = "country", country_list = c("Russia", "China", "US", "United
                                                        Kingdom", "Canada", "France", 
                                                        "Germany", "Italy", "Japan", "Iran"))
cv_summary(cv_data, facet = "variable", country_list = c("Russia", "China", "US", "United
                                                        Kingdom", "Canada", "France", 
                                                         "Germany", "Italy", "Japan", "Iran"))

#Challenge 3
#Use the dataset and function generated above to return summary data for ALL 
#countries in the dataset, and then filter this returned dataset to only those 
#countries with populations of over 1 million, storing this dataset as a tibble d. 
#How many countries does this tibble include?
results <- cv_summary(cv_data, country_list = "All")
results$cum_cases_plot
results$daily_cases_plot
data <- results$total

d <- filter(data, population > 1000000)
d
nrow(d)/6

#Challenge 4
#Filter d to generate two additional tibbles, overall and daily that include only 
#data on the variables “confirmed” and “daily_confirmed” cases, respectively. 
#Depending on the dataset, the case and rate variables either reflect the overall 
#(i.e., across the pandemic) or maximum daily number of cases and number of cases 
#recorded per million people in the population. Which 10 countries have experienced 
#the highest over rate of confirmed cases? Which 10 countries have experienced the 
#highest single-day rate of confirmed cases?

overall <- filter(d, variable == "confirmed")
daily <- filter(d, variable == "daily_confirmed")


#Challenge 5
#Run a linear model to evaluate how the overall infection rate (rate) is related to 
#the variables population density (density), population size (population), gross 
#domestic product per capita (gdp_capita), and overall income level (income). In 
#doing so, you should run exploratory visualizations to see whether or not the four 
#numeric variables should be transformed.

#Based on the full model, what predictors variables have slopes significantly 
#different from zero?

lm(rate ~ density + population + gdp_capita + income, data = cv_data)

ggplot(data = cv_data, aes(x = rate, y = density, z = population)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data = cv_data, aes(x = density, y = gdp_capita, z = income)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data = cv_data, aes(x = income, y = rate, z = gdp_capita)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

#Challenge 6
#Run stepwise selection using AIC to evaluate whether the full model or a nested, 
#simpler model is preferred. What is the best model (based on AIC) of the possible ones 
#involving these 4 predictors? What are the “pseudo- R2” values associated with the full 
#and “best” models?
#HINT: Check out Module 22 on module selection using AIC with the {MASS} package and 
#Module 24 on calculating “pseudo-R2” values with {MuMIn}.
#Repeat this modeling process to evaluate what combination of explanatory variables best 
#predicts maximum daily infection rate. Are the important predictors the same? What 
#additional or different variables are included?

library(MASS)
m <- lm(data = cv_data, rate ~ density + population + gdp_capita + income)
(s <- stepAIC(m, scope = . ~ ., direction = "both"))
m
summary(s)

m2 <- lm(data = cv_data, rate ~ density)
(s2 <- stepAIC(m, scope = . ~ ., direction = "both"))
m2

#Challenge 7
#To the best model you determined in CHALLENGE 6 for predicting the maximum daily 
#infection rate, add in the maximum social distancing (max_sd) and maximum movement 
#restriction (max_mr) score per country. Do either of these additional variable improve 
#the model significantly?




