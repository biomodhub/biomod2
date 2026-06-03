## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 3.75
)

## ----message = FALSE----------------------------------------------------------
library(ggplot2)
library(ggtext)

base <- ggplot(data.frame(x = c(-5, 5)), aes(x)) +
  stat_function(fun = ~ .x*.x)

base

## -----------------------------------------------------------------------------
base +
  labs(
    x = "independent variable *x*",
    y = "dependent variable *y* = *x*<sup>2</sup>"
  ) +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )

## -----------------------------------------------------------------------------
base +
  labs(
    x = "independent variable *x*",
    y = "dependent variable *y* = *x*<sup>2</sup>"
  ) +
  theme(
    axis.title.x = element_markdown(color = "blue"),
    axis.title.y = element_markdown(size = rel(0.8))
  )

## -----------------------------------------------------------------------------
base +
  labs(
    x = "independent variable *x*",
    y = "dependent variable *y* = *x*<sup>2</sup>"
  ) +
  theme(
    axis.title = element_text(color = "blue", size = rel(0.8)),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )

## -----------------------------------------------------------------------------
base +
  labs(
    x = "independent variable *x*",
    y = "dependent variable *y* = *x*<sup>2</sup>"
  ) +
  scale_y_continuous(position = "right") +
  theme(
    axis.title = element_text(color = "blue", size = rel(0.8)),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown()
  )

## -----------------------------------------------------------------------------
base +
  labs(
    x = "independent variable *x*",
    y = "dependent variable *y* = *x*<sup>2</sup>"
  ) +
  scale_y_continuous(position = "right") +
  theme(
    axis.title = element_text(color = "blue", size = rel(0.8)),
    axis.title.x = element_markdown(),
    axis.title.y.right = element_markdown()
  )

## ----fig.width = 5, message = FALSE-------------------------------------------
library(dplyr)

mtcars %>%
  mutate(
    transmission = ifelse(am == 1, "automatic", "manual")
  ) %>%
  ggplot(aes(hp, mpg, color = transmission)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(automatic = "#0072B2", manual = "#D55E00"),
    guide = "none"
  ) +
  labs(
    x = "Horse power",
    y = "Miles per gallon (MPG)",
    title = "<span style = 'font-size:14pt; font-family:Helvetica;'>Transmission type impacts fuel efficiency</span><br>
MPG is higher for <span style = 'color:#0072B2;'>automatic</span>
than for <span style = 'color:#D55E00;'>manual</span> transmissions"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times"),
    plot.title.position = "plot",
    plot.title = element_markdown(size = 11, lineheight = 1.2)
  )


## ----fig.width = 5, message = FALSE-------------------------------------------
base <- mtcars %>%
  mutate(
    transmission = ifelse(am == 1, "automatic", "manual")
  ) %>%
  ggplot(aes(hp, mpg, color = transmission)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(automatic = "#0072B2", manual = "#D55E00"),
    guide = "none"
  ) +
  labs(
    x = "Horse power",
    y = "Miles per gallon (MPG)",
    title = "Transmission type impacts fuel efficiency<br>
<span style = 'font-size:10pt;'>Miles per gallon (MPG) is on average higher for cars
with <span style = 'color:#0072B2;'>automatic transmission</span> than for cars with
<span style = 'color:#D55E00;'>manual transmission.</span> However, MPG generally
declines with increasing horse power.</span>"
  ) +
  theme_bw() + theme(plot.title.position = "plot")

base +
  theme(
    plot.title = element_textbox_simple(
      size = 14, lineheight = 1, padding = margin(0, 0, 5, 0)
    )
  )

## ----fig.width = 5, message = FALSE-------------------------------------------
base +
  theme(
    plot.title = element_textbox_simple(
      size = 14, lineheight = 1,
      linetype = 1, # turn on border
      box.color = "#748696", # border color
      fill = "#F0F7FF", # background fill color
      r = grid::unit(3, "pt"), # radius for rounded corners
      padding = margin(5, 5, 5, 5), # padding around text inside the box
      margin = margin(0, 0, 10, 0) # margin outside the box
    )
  )

## ----fig.width = 5, message = FALSE-------------------------------------------
base +
  theme(
    plot.title = element_textbox_simple(
      size = 14, lineheight = 1, 
      width = grid::unit(4, "in"), # fixed width
      hjust = 1, # alignment of box relative to plot
      linetype = 1, # turn on border
      box.color = "#748696", # border color
      fill = "#F0F7FF", # background fill color
      r = grid::unit(3, "pt"), # radius for rounded corners
      padding = margin(5, 5, 5, 5), # padding around text inside the box
      margin = margin(0, 0, 10, 0) # margin outside the box
    )
  )

## -----------------------------------------------------------------------------
mtcars %>%
  mutate(
    transmission = ifelse(am == 1, "automatic", "manual")
  ) %>%
  ggplot(aes(hp, mpg, color = transmission)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(automatic = "#0072B2", manual = "#D55E00"),
    guide = "none"
  ) +
  labs(
    x = "Horse power",
    y = "Miles per gallon (MPG) is on average higher for cars with <span style =
'color:#0072B2;'>automatic transmission</span> than for cars with <span style =
'color:#D55E00;'>manual transmission.</span>",
    title = "Transmission type impacts fuel efficiency"
  ) +
  theme_bw() + 
  theme(
    plot.title.position = "plot",
    axis.title.y = element_textbox_simple(
      orientation = "left-rotated",
      width = grid::unit(2.5, "in"),
      hjust = 0, fill = "#F0F7FF",  
      padding = margin(5, 5, 5, 5),
      margin = margin(0, 0, 10, 0)
    )
  )

