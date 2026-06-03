## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 3.75
)

## ----fig.width = 6, fig.height = 3, message = FALSE---------------------------
library(ggplot2)
library(dplyr)
library(glue)

iris_cor <- iris %>% 
  group_by(Species) %>%
  summarize(r_square = cor(Sepal.Length, Sepal.Width)^2) %>%
  mutate(
    # location of each text label in data coordinates
    Sepal.Length = 8, Sepal.Width = 4.5,
    # text label containing r^2 value 
    label = glue("r^2 = {round(r_square, 2)}")
  )

iris_cor

iris_facets <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~Species) +
  theme_bw()

iris_facets + 
  geom_text(
    data = iris_cor,
    aes(label = label),
    hjust = 1, vjust = 1
  )

## ----fig.width = 6, fig.height = 3, message = FALSE---------------------------
library(ggtext)

iris_cor_md <- iris_cor %>% 
  mutate(
    # markdown version of text label
    label = glue("*r*<sup>2</sup> = {round(r_square, 2)}")
  )

iris_cor_md

iris_facets + 
  geom_richtext(
    data = iris_cor_md,
    aes(label = label),
    hjust = 1, vjust = 1
  )

## ----fig.width = 6, fig.height = 3--------------------------------------------
iris_facets + 
  geom_richtext(
    data = iris_cor_md,
    aes(label = label),
    hjust = 1, vjust = 1,
    # remove label background and outline
    fill = NA, label.color = NA,
    # remove label padding, since we have removed the label outline
    label.padding = grid::unit(rep(0, 4), "pt") 
  )

## ----fig.width = 6, fig.height = 3--------------------------------------------
iris_facets + 
  aes(colour = Species) +
  geom_richtext(
    data = iris_cor_md,
    aes(
      label = label,
      fill = after_scale(alpha(colour, .2))
    ),
    text.colour = "black",
    hjust = 1, vjust = 1
  ) +
  theme(legend.position = "none")

## ----fig.width = 6, fig.height = 3--------------------------------------------
iris_facets + 
  aes(colour = Species) +
  geom_richtext(
    data = iris_cor_md,
    aes(
      x = 7.5,
      label = label,
      fill = after_scale(alpha(colour, .2))
    ),
    text.colour = "black",
    hjust = 1, vjust = 1,
    angle = 30
  ) +
  theme(legend.position = "none")

## ----fig.width = 6, fig.height = 3--------------------------------------------
df <- data.frame(
  x = 0.1,
  y = 0.8,
  label = "*Lorem ipsum dolor sit amet,* consectetur adipiscing
elit. Quisque tincidunt eget arcu in pulvinar. Morbi varius leo
vel consectetur luctus. **Morbi facilisis justo non fringilla.**
Vivamus sagittis sem felis, vel lobortis risus mattis eget. Nam
quis imperdiet felis, in convallis elit."
)

p <- ggplot() +
  geom_textbox(
    data = df,
    aes(x, y, label = label),
    width = grid::unit(0.73, "npc"), # 73% of plot panel width
    hjust = 0, vjust = 1
  ) +
  xlim(0, 1) + ylim(0, 1)

p

## ----fig.width = 4, fig.height = 4--------------------------------------------
p

## ----fig.width = 4, fig.height = 4--------------------------------------------
ggplot() +
  geom_textbox(
    data = df,
    aes(x, y, label = label),
    width = grid::unit(0.73, "npc"), # 73% of plot panel width
    hjust = 0, vjust = 1,
    halign = 0.5 # centered text
  ) +
  xlim(0, 1) + ylim(0, 1)

## -----------------------------------------------------------------------------
df <- data.frame(
  x = 0.5,
  y = 0.5,
  label = "The quick brown fox jumps over the lazy dog.",
  orientation = c("upright", "left-rotated", "inverted", "right-rotated")
)

ggplot() +
  geom_textbox(
    data = df,
    aes(x, y, label = label, orientation = orientation),
    width = grid::unit(1.5, "in"),
    height = grid::unit(1.5, "in"),
    box.margin = grid::unit(rep(0.25, 4), "in"),
    hjust = 0, vjust = 1
  ) +
  xlim(0, 1) + ylim(0, 1) +
  scale_discrete_identity(aesthetics = "orientation")


