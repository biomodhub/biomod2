library(tibble)
library(dplyr)
library(tidyr)

# Prep and select basic data
data("mtcars")
base_data <- mtcars |>
  rownames_to_column() |>
  mutate(gear = paste(gear, "Gears"),
         cyl = paste(cyl, "Cylinders")) |>
  select(rowname, cyl, gear, wt, mpg, qsec)

base_data |>
  pivot_longer(names_to = "per_metric",
               cols = c(wt, mpg, qsec)) |>
  group_by(cyl, gear, per_metric) |>
  summarise(value_Mean = round(mean(value), 1),
            value_Min = round(min(value), 1),
            value_Max = round(max(value), 1),
            .groups = "drop") |>
  pivot_wider(names_from = per_metric,
              values_from = starts_with("value_")) |>
  # Round the values into a nicer format where we want the weights to have two decimals
  txtRound(ends_with("_wt"), digits = 2) |>
  txtRound(starts_with("value") & !ends_with("_wt"), digits = 1) |>
  # Convert into long format
  pivot_longer(cols = starts_with("value_"), names_prefix = "value_") |>
  separate(name, into = c("summary_stat", "per_metric")) |>
  # Without sorting the row groups wont appear right
  # If the columns end up in the wrong order you may want to change the columns
  # into factors
  arrange(per_metric) |>
  addHtmlTableStyle(align = "r") |>
  tidyHtmlTable(
    header = gear,
    cgroup = cyl,
    rnames = summary_stat,
    rgroup = per_metric,
    skip_removal_warning = TRUE
  )
