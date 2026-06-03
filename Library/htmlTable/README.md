[![Downloads](https://cranlogs.r-pkg.org/badges/htmlTable)](https://cran.r-project.org/package=htmlTable)

# Basics

The **htmlTable** package is intended for generating tables using [HTML](https://en.wikipedia.org/wiki/HTML) formatting. This format is compatible with [Markdown](https://rmarkdown.rstudio.com/) when used for HTML-output. The most basic table can easily be created by just passing a `matrix` or a `data.frame` to the `htmlTable`-function:

```r
library(magrittr)
library(htmlTable)
# A simple output
output <- matrix(1:4,
                 ncol=2,
                 dimnames = list(list("Row 1", "Row 2"),
                                 list("Column 1", "Column 2")))
htmlTable(output)
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Column 1</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Column 2</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>Row 1</td>
		<td style='text-align: center;'>1</td>
		<td style='text-align: center;'>3</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>Row 2</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>4</td>
	</tr>
	</tbody>
</table>

If you are using `dplyr` and `tidyverse` a convenient wrapper is the `tidyHtmlTable` function (check out `vignette("tidyHtmlTable")`). A simple example of the `tidyHtmlTable` would look something like this:

```r
library(tidyverse)
library(glue)
mtcars |>
  as_tibble(rownames = "rnames") |>
  filter(cyl == 6 & qsec < 18) |>
  pivot_longer(names_to = "per_metric",
               cols = c(hp, mpg, qsec)) |>
  arrange(gear, rnames) |>
  mutate(gear = glue("{gear} gears")) |>
  addHtmlTableStyle(align = "r") |>
  tidyHtmlTable(header = per_metric, rnames = rnames, rgroup = gear,
                caption = "A simple <code>tidyHtmlTable</code> example using <code>mtcars</code>")
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
	<thead>
	<tr><td colspan='4' style='text-align: left;'>
	  A simple <code>tidyHtmlTable</code> example using <code>mtcars</code>
	</td></tr>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>hp</th>
		<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>mpg</th>
		<th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: right;'>qsec</th>
	</tr>
	</thead>
	<tbody>
	<tr><td colspan='4' style='font-weight: 900;'>4 gears</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;Mazda RX4</td>
		<td style='text-align: right;'>110</td>
		<td style='text-align: right;'>21</td>
		<td style='text-align: right;'>16.46</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;Mazda RX4 Wag</td>
		<td style='text-align: right;'>110</td>
		<td style='text-align: right;'>21</td>
		<td style='text-align: right;'>17.02</td>
	</tr>
	<tr><td colspan='4' style='font-weight: 900;'>5 gears</td></tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>&nbsp;&nbsp;Ferrari Dino</td>
		<td style='border-bottom: 2px solid grey; text-align: right;'>175</td>
		<td style='border-bottom: 2px solid grey; text-align: right;'>19.7</td>
		<td style='border-bottom: 2px solid grey; text-align: right;'>15.5</td>
	</tr>
	</tbody>
</table>

# Advanced

While it may be sufficient for basic tables a more advanced layout is often needed in medical publications with elements such as:

- row groups
- column spanners
- table spanners
- caption
- table footer
- zebra coloring (also know as _banding_):
  - rows
  - columns

As many journals require that a MS Word-document is submitted it is furthermore also important that the table imports correctly to a word processor, i.e. that the table doesn't only look nice in a web browser but also in the final document. The `htmlTable`-function is written for all these purposes.

**Note:** Due to GitHub CSS-styles the rows get automatically zebra-striped (in a bad way), borders get overridden and I haven't been able to figure out how to change this. See the vignette for a correct example: `vignette("general", package = "htmlTable")`

For demonstration purposes we will setup a basic matrix:

```r
mx <-
  matrix(ncol=6, nrow=8) |>
  set_rownames(paste(c("1st", "2nd", "3rd",
                       paste0(4:8, "th")),
                     "row")) |>
  set_colnames(paste(c("1st", "2nd", "3rd",
                       paste0(4:6, "th")),
                     "hdr"))

for (nr in 1:nrow(mx)){
  for (nc in 1:ncol(mx)){
    mx[nr, nc] <-
      paste0(nr, ":", nc)
  }
}
```

## Row groups

The purpose of the row groups is to group variables that belong to the same group, e.g. a factored variable with more than two levels often benefit from grouping variables together.

```r
htmlTable(mx,
          rgroup = paste("Group", LETTERS[1:3]),
          n.rgroup = c(2,4,nrow(mx) - 6))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr><td colspan='7' style='font-weight: 900;'>Group A</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
		<td style='text-align: center;'>1:3</td>
		<td style='text-align: center;'>1:4</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='text-align: center;'>2:2</td>
		<td style='text-align: center;'>2:3</td>
		<td style='text-align: center;'>2:4</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr><td colspan='7' style='font-weight: 900;'>Group B</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='text-align: center;'>3:2</td>
		<td style='text-align: center;'>3:3</td>
		<td style='text-align: center;'>3:4</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='text-align: center;'>4:2</td>
		<td style='text-align: center;'>4:3</td>
		<td style='text-align: center;'>4:4</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='text-align: center;'>5:2</td>
		<td style='text-align: center;'>5:3</td>
		<td style='text-align: center;'>5:4</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='text-align: center;'>6:2</td>
		<td style='text-align: center;'>6:3</td>
		<td style='text-align: center;'>6:4</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr><td colspan='7' style='font-weight: 900;'>Group C</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='text-align: center;'>7:2</td>
		<td style='text-align: center;'>7:3</td>
		<td style='text-align: center;'>7:4</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>&nbsp;&nbsp;8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

We can easily mix row groups with regular variables by having an empty row group name `""`:

```r
htmlTable(mx,
          rgroup = c(paste("Group", LETTERS[1:2]), ""),
          n.rgroup = c(2,4,nrow(mx) - 6))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr><td colspan='7' style='font-weight: 900;'>Group A</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
		<td style='text-align: center;'>1:3</td>
		<td style='text-align: center;'>1:4</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='text-align: center;'>2:2</td>
		<td style='text-align: center;'>2:3</td>
		<td style='text-align: center;'>2:4</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr><td colspan='7' style='font-weight: 900;'>Group B</td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='text-align: center;'>3:2</td>
		<td style='text-align: center;'>3:3</td>
		<td style='text-align: center;'>3:4</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='text-align: center;'>4:2</td>
		<td style='text-align: center;'>4:3</td>
		<td style='text-align: center;'>4:4</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='text-align: center;'>5:2</td>
		<td style='text-align: center;'>5:3</td>
		<td style='text-align: center;'>5:4</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='text-align: center;'>6:2</td>
		<td style='text-align: center;'>6:3</td>
		<td style='text-align: center;'>6:4</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='text-align: center;'>7:2</td>
		<td style='text-align: center;'>7:3</td>
		<td style='text-align: center;'>7:4</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

When mixing row groups with variables without row groups we may want to omit the bold formatting of the row group label. As of htmlTable version 2.0
you can separate the css styling using `addHtmlTableStyle`:

```r
mx |>
  addHtmlTableStyle(css.rgroup = "") |>
  htmlTable(rgroup = c(paste("Group", LETTERS[1:2]), ""),
            n.rgroup = c(2,4,nrow(mx) - 6))
```

### Row highlighting

If you want to highlight an entire row based on a condition, use `highlightRow()` before `htmlTable()`. The condition is evaluated against the table columns, and you can also use `.rowname`:

```r
df <- data.frame(
  Name = c("Max", "Eva", "Nils"),
  Score = c(10, 20, 30),
  stringsAsFactors = FALSE
)

df |>
  highlightRow(Name == "Max", style = "warning") |>
  highlightRow(.rowname == "3", style = "background-color: #d1ecf1; color: #0c5460;") |>
  htmlTable(rnames = FALSE)
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
  <thead>
  <tr>
    <th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Name</th>
    <th style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Score</th>
  </tr>
  </thead>
  <tbody>
    <tr style='background-color: #fff3cd; color: #856404;'>
    <td style='background-color: #fff3cd; color: #856404; text-align: center;'>Max</td>
    <td style='background-color: #fff3cd; color: #856404; text-align: center;'>10</td>
  </tr>
  <tr>
    <td style='text-align: center;'>Eva</td>
    <td style='text-align: center;'>20</td>
  </tr>
  <tr style='background-color: #d1ecf1; color: #0c5460;'>
    <td style='background-color: #d1ecf1; color: #0c5460; border-bottom: 2px solid grey; text-align: center;'>Nils</td>
    <td style='background-color: #d1ecf1; color: #0c5460; border-bottom: 2px solid grey; text-align: center;'>30</td>
  </tr>
  </tbody>
</table>

## Column spanners

A column spanner spans 2 or more columns:

```r
htmlTable(mx,
          cgroup = c("Cgroup 1", "Cgroup 2"),
          n.cgroup = c(2,4))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-top: 2px solid grey;'></th>
		<th colspan='2' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Cgroup 1</th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='4' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Cgroup 2</th>
	</tr>
	<tr>
		<th style='border-bottom: 1px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:3</td>
		<td style='text-align: center;'>1:4</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='text-align: center;'>2:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:3</td>
		<td style='text-align: center;'>2:4</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='text-align: center;'>3:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:3</td>
		<td style='text-align: center;'>3:4</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='text-align: center;'>4:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:3</td>
		<td style='text-align: center;'>4:4</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='text-align: center;'>5:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:3</td>
		<td style='text-align: center;'>5:4</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='text-align: center;'>6:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:3</td>
		<td style='text-align: center;'>6:4</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='text-align: center;'>7:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:3</td>
		<td style='text-align: center;'>7:4</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

It can sometimes be convenient to have column spanners in multiple levels:

```r
htmlTable(mx,
          cgroup = rbind(c("", "Column spanners", NA),
                         c("", "Cgroup 1", "Cgroup 2")),
          n.cgroup = rbind(c(1,2,NA),
                           c(2,2,2)))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-top: 2px solid grey;'></th>
		<th colspan='2' style='font-weight: 900; border-top: 2px solid grey; text-align: center;'></th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='5' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Column spanners</th>
	</tr>
	<tr>
		<th style=''></th>
		<th colspan='2' style='font-weight: 900; text-align: center;'></th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='2' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 1</th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='2' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 2</th>
	</tr>
	<tr>
		<th style='border-bottom: 1px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:3</td>
		<td style='text-align: center;'>1:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='text-align: center;'>2:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:3</td>
		<td style='text-align: center;'>2:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='text-align: center;'>3:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:3</td>
		<td style='text-align: center;'>3:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='text-align: center;'>4:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:3</td>
		<td style='text-align: center;'>4:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='text-align: center;'>5:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:3</td>
		<td style='text-align: center;'>5:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='text-align: center;'>6:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:3</td>
		<td style='text-align: center;'>6:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='text-align: center;'>7:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:3</td>
		<td style='text-align: center;'>7:4</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

Above example allows the column spanner to be a sum of the underlying cgroups (see n.cgroup), this is not required by the function:

```r
htmlTable(mx,
          cgroup = rbind(c("", "Column spanners", NA),
                         c("", "Cgroup 1", "Cgroup 2")),
          n.cgroup = rbind(c(1,5,NA),
                           c(2,1,3)))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-top: 2px solid grey;'></th>
		<th colspan='1' style='font-weight: 900; border-top: 2px solid grey; text-align: center;'></th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='7' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Column spanners</th>
	</tr>
	<tr>
		<th style=''></th>
		<th colspan='3' style='font-weight: 900; text-align: center;'></th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='1' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 1</th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='3' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 2</th>
	</tr>
	<tr>
		<th style='border-bottom: 1px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>1:4</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>2:4</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>3:4</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>4:4</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>5:4</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>6:4</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:3</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='text-align: center;'>7:4</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

## Table spanners

A table spanner is similar to rgroup but has the primary purpose of combining 2 or more tables with the same columns into one:

```r
htmlTable(mx,
          tspanner = paste("Spanner", LETTERS[1:3]),
          n.tspanner = c(2,4,nrow(mx) - 6))
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr><td colspan='7' style='font-weight: 900; text-transform: capitalize; text-align: left;'>Spanner A</td></tr>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
		<td style='text-align: center;'>1:3</td>
		<td style='text-align: center;'>1:4</td>
		<td style='text-align: center;'>1:5</td>
		<td style='text-align: center;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>2nd row</td>
		<td style='text-align: center;'>2:1</td>
		<td style='text-align: center;'>2:2</td>
		<td style='text-align: center;'>2:3</td>
		<td style='text-align: center;'>2:4</td>
		<td style='text-align: center;'>2:5</td>
		<td style='text-align: center;'>2:6</td>
	</tr>
	<tr><td colspan='7' style='font-weight: 900; text-transform: capitalize; text-align: left; border-top: 1px solid grey;'>Spanner B</td></tr>
	<tr>
		<td style='text-align: left;'>3rd row</td>
		<td style='text-align: center;'>3:1</td>
		<td style='text-align: center;'>3:2</td>
		<td style='text-align: center;'>3:3</td>
		<td style='text-align: center;'>3:4</td>
		<td style='text-align: center;'>3:5</td>
		<td style='text-align: center;'>3:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>4th row</td>
		<td style='text-align: center;'>4:1</td>
		<td style='text-align: center;'>4:2</td>
		<td style='text-align: center;'>4:3</td>
		<td style='text-align: center;'>4:4</td>
		<td style='text-align: center;'>4:5</td>
		<td style='text-align: center;'>4:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>5th row</td>
		<td style='text-align: center;'>5:1</td>
		<td style='text-align: center;'>5:2</td>
		<td style='text-align: center;'>5:3</td>
		<td style='text-align: center;'>5:4</td>
		<td style='text-align: center;'>5:5</td>
		<td style='text-align: center;'>5:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>6th row</td>
		<td style='text-align: center;'>6:1</td>
		<td style='text-align: center;'>6:2</td>
		<td style='text-align: center;'>6:3</td>
		<td style='text-align: center;'>6:4</td>
		<td style='text-align: center;'>6:5</td>
		<td style='text-align: center;'>6:6</td>
	</tr>
	<tr><td colspan='7' style='font-weight: 900; text-transform: capitalize; text-align: left; border-top: 1px solid grey;'>Spanner C</td></tr>
	<tr>
		<td style='text-align: left;'>7th row</td>
		<td style='text-align: center;'>7:1</td>
		<td style='text-align: center;'>7:2</td>
		<td style='text-align: center;'>7:3</td>
		<td style='text-align: center;'>7:4</td>
		<td style='text-align: center;'>7:5</td>
		<td style='text-align: center;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>8th row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:2</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:3</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:4</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:5</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>8:6</td>
	</tr>
	</tbody>
</table>

## Table caption

The table caption is simply the table description and can be either located above or below the table:

```r
htmlTable(mx[1:2,1:2],
          caption="A table caption above")
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_4'>
	<thead>
	<tr><td colspan='3' style='text-align: left;'>
	Table 5:  A table caption above</td></tr>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>2nd row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:2</td>
	</tr>
	</tbody>
</table>

```r
mx[1:2,1:2] |>
  addHtmlTableStyle(pos.caption = "bottom") |>
  htmlTable(caption="A table caption below")
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_5'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>2nd row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:2</td>
	</tr>
	</tbody>
	<tr><td colspan='3' style='text-align: left;'>
	Table 6:  A table caption below</td></tr>
</table>

A more interesting detail that the function allows for is table numbering, initialized by:

```r
options(table_counter = TRUE)
```

```r
htmlTable(mx[1:2,1:2],
          caption="A table caption with a numbering")
```

<table class='gmisc_table' style='border-collapse: collapse;' >
	<thead>
	<tr><td colspan='3' style='text-align: left;'>
	Table 1:  A table caption with a numbering</td></tr>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>2nd row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:2</td>
	</tr>
	</tbody>
</table>

As we often want to reference the table number in the text there are two associated functions:

```r
tblNoLast()
```

```
## [1] 1
```

```r
tblNoNext()
```

```
## [1] 2
```

## Table footer

The footer usually contains specifics regarding variables and is always located at the foot of the table:

```r
htmlTable(mx[1:2,1:2],
          tfoot="A table footer")
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_1'>
	<thead>
	<tr>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>2nd hdr</th>
	</tr>
	</thead><tbody>
	<tr>
		<td style='text-align: left;'>1st row</td>
		<td style='text-align: center;'>1:1</td>
		<td style='text-align: center;'>1:2</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>2nd row</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:1</td>
		<td style='border-bottom: 2px solid grey; text-align: center;'>2:2</td>
	</tr>
	</tbody>
	<tfoot><tr><td colspan='3'>
	A table footer</td></tr></tfoot>
</table>

## Putting it all together

Now if we want to do everything in one table it may look like this:

```r
mx |>
  addHtmlTableStyle(col.columns = c(rep("none", 2), rep("#F5FBFF", 4)),
                    col.rgroup = c("none", "#F7F7F7"),
                    css.cell = "padding-left: .5em; padding-right: .2em;",
                    align="r") |>
  htmlTable(rgroup = paste("Group", LETTERS[1:3]),
            n.rgroup = c(2, 4),
            cgroup = rbind(c("", "Column spanners", NA),
                           c("", "Cgroup 1", "Cgroup 2&dagger;")),
            n.cgroup = rbind(c(1, 2, NA), c(2, 2, 2)),
            caption="A table with column spanners, row groups, and zebra striping",
            tfoot="&dagger; A table footer commment",
            cspan.rgroup = 2)
```

<table class='gmisc_table' style='border-collapse: collapse;'  id='table_1'>
	<thead>
	<tr><td colspan='9' style='text-align: left;'>
	Table 2:  A table with column spanners, row groups, and zebra striping</td></tr>
	<tr>
		<th style='border-top: 2px solid grey;'></th>
		<th colspan='2' style='font-weight: 900; border-top: 2px solid grey; text-align: center;'></th><th style='border-top: 2px solid grey;; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='5' style='font-weight: 900; border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Column spanners</th>
	</tr>
	<tr>
		<th style=''></th>
		<th colspan='2' style='font-weight: 900; text-align: center;'></th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='2' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 1</th><th style='; border-bottom: hidden;'>&nbsp;</th>
		<th colspan='2' style='font-weight: 900; border-bottom: 1px solid grey; text-align: center;'>Cgroup 2&dagger;</th>
	</tr>
	<tr>
		<th style='border-bottom: 1px solid grey;'> </th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>1st hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>2nd hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>3rd hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>4th hdr</th>
		<th style='border-bottom: 1px solid grey;' colspan='1'>&nbsp;</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>5th hdr</th>
		<th style='border-bottom: 1px solid grey; text-align: center;'>6th hdr</th>
	</tr>
	</thead><tbody>
	<tr><td colspan='2' style='font-weight: 900;'>Group A</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right;'></td>
		<td style='font-weight: 900;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='font-weight: 900; background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;1st row</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>1:1</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>1:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>1:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>1:4</td>
		<td style='background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>1:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>1:6</td>
	</tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;2nd row</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>2:1</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>2:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>2:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>2:4</td>
		<td style='background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>2:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>2:6</td>
	</tr>
	<tr><td colspan='2' style='font-weight: 900; background-color: #f7f7f7;'>Group B</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; background-color: #f7f7f7; text-align: right;'></td>
		<td style='font-weight: 900; background-color: #f7f7f7;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #F6F9FB;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #F6F9FB;'></td>
		<td style='font-weight: 900; background-color: #F6F9FB;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #F6F9FB;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #F6F9FB;'></td></tr>
	<tr style='background-color: #f7f7f7;'>
		<td style='background-color: #f7f7f7; text-align: left;'>&nbsp;&nbsp;3rd row</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>3:1</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>3:2</td>
		<td style='background-color: #f7f7f7;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>3:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>3:4</td>
		<td style='background-color: #F6F9FB;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>3:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>3:6</td>
	</tr>
	<tr style='background-color: #f7f7f7;'>
		<td style='background-color: #f7f7f7; text-align: left;'>&nbsp;&nbsp;4th row</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>4:1</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>4:2</td>
		<td style='background-color: #f7f7f7;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>4:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>4:4</td>
		<td style='background-color: #F6F9FB;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>4:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>4:6</td>
	</tr>
	<tr style='background-color: #f7f7f7;'>
		<td style='background-color: #f7f7f7; text-align: left;'>&nbsp;&nbsp;5th row</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>5:1</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>5:2</td>
		<td style='background-color: #f7f7f7;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>5:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>5:4</td>
		<td style='background-color: #F6F9FB;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>5:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>5:6</td>
	</tr>
	<tr style='background-color: #f7f7f7;'>
		<td style='background-color: #f7f7f7; text-align: left;'>&nbsp;&nbsp;6th row</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>6:1</td>
		<td style='padding-left: .5em; padding-right: .2em; background-color: #f7f7f7; text-align: right;'>6:2</td>
		<td style='background-color: #f7f7f7;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>6:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>6:4</td>
		<td style='background-color: #F6F9FB;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>6:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #F6F9FB;'>6:6</td>
	</tr>
	<tr><td colspan='2' style='font-weight: 900;'>Group C</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right;'></td>
		<td style='font-weight: 900;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='font-weight: 900; background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td>
		<td style='padding-left: .5em; padding-right: .2em; font-weight: 900; text-align: right; background-color: #f5fbff;'></td></tr>
	<tr>
		<td style='text-align: left;'>&nbsp;&nbsp;7th row</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>7:1</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right;'>7:2</td>
		<td style='' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>7:3</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>7:4</td>
		<td style='background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>7:5</td>
		<td style='padding-left: .5em; padding-right: .2em; text-align: right; background-color: #f5fbff;'>7:6</td>
	</tr>
	<tr>
		<td style='border-bottom: 2px solid grey; text-align: left;'>&nbsp;&nbsp;8th row</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right;'>8:1</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right;'>8:2</td>
		<td style='border-bottom: 2px solid grey;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right; background-color: #f5fbff;'>8:3</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right; background-color: #f5fbff;'>8:4</td>
		<td style='border-bottom: 2px solid grey; background-color: #f5fbff;' colspan='1'>&nbsp;</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right; background-color: #f5fbff;'>8:5</td>
		<td style='padding-left: .5em; padding-right: .2em; border-bottom: 2px solid grey; text-align: right; background-color: #f5fbff;'>8:6</td>
	</tr>
	</tbody>
	<tfoot><tr><td colspan='9'>
	&dagger; A table footer comment</td></tr></tfoot>
</table>
