## -----------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = tempdir())
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(tibble.print_min = 3)

## -----------------------------------------------------------------------------
library(vroom)

## -----------------------------------------------------------------------------
# See where the example file is stored on your machine
file <- vroom_example("mtcars.csv")
file

# Read the file, by default vroom will guess the delimiter automatically.
vroom(file)

# You can also specify it explicitly, which is (slightly) faster, and safer if
# you know how the file is delimited.
vroom(file, delim = ",")

## -----------------------------------------------------------------------------
ve <- grep("mtcars-[0-9].csv", vroom_examples(), value = TRUE)
files <- sapply(ve, vroom_example)
files

## -----------------------------------------------------------------------------
vroom(files)

## -----------------------------------------------------------------------------
vroom(files, id = "path")

## -----------------------------------------------------------------------------
file <- vroom_example("mtcars.csv.gz")

vroom(file)

## -----------------------------------------------------------------------------
zip_file <- vroom_example("mtcars-multi-cyl.zip")
filenames <- unzip(zip_file, list = TRUE)$Name
filenames

# imagine we only want to read 2 of the 3 files
vroom(purrr::map(filenames[c(1, 3)], \(x) unz(zip_file, x)))

## -----------------------------------------------------------------------------
# file <- "https://raw.githubusercontent.com/tidyverse/vroom/main/inst/extdata/mtcars.csv"
# vroom(file)

## -----------------------------------------------------------------------------
# file <- "https://raw.githubusercontent.com/tidyverse/vroom/main/inst/extdata/mtcars.csv.gz"
# vroom(file)

## -----------------------------------------------------------------------------
file <- vroom_example("mtcars.csv.gz")

vroom(file, col_select = c(model, cyl, gear))

## -----------------------------------------------------------------------------
vroom(file, col_select = c(1, 3, 11))

## -----------------------------------------------------------------------------
vroom(file, col_select = starts_with("d"))

## -----------------------------------------------------------------------------
vroom(file, col_select = c(car = model, everything()))

## -----------------------------------------------------------------------------
fwf_sample <- vroom_example("fwf-sample.txt")
cat(readLines(fwf_sample))

## -----------------------------------------------------------------------------
vroom_fwf(fwf_sample, fwf_empty(fwf_sample, col_names = c("first", "last", "state", "ssn")))

## -----------------------------------------------------------------------------
vroom_fwf(fwf_sample, fwf_widths(c(20, 10, 12), c("name", "state", "ssn")))

## -----------------------------------------------------------------------------
vroom_fwf(fwf_sample, fwf_positions(c(1, 30), c(20, 42), c("name", "ssn")))

## -----------------------------------------------------------------------------
vroom_fwf(fwf_sample, fwf_cols(name = 20, state = 10, ssn = 12))

## -----------------------------------------------------------------------------
vroom_fwf(fwf_sample, fwf_cols(name = c(1, 20), ssn = c(30, 42)))

## -----------------------------------------------------------------------------
# read the 'hp' columns as an integer
vroom(vroom_example("mtcars.csv"), col_types = c(hp = "i"))

# also skip reading the 'cyl' column
vroom(vroom_example("mtcars.csv"), col_types = c(hp = "i", cyl = "_"))

# also read the gears as a factor
vroom(vroom_example("mtcars.csv"), col_types = c(hp = "i", cyl = "_", gear = "f"))

## -----------------------------------------------------------------------------
vroom(vroom_example("mtcars.csv"), col_types = c(.default = "c"))

## -----------------------------------------------------------------------------
vroom(
  vroom_example("mtcars.csv"),
  col_types = list(hp = col_integer(), cyl = col_skip(), gear = col_factor())
)

## -----------------------------------------------------------------------------
vroom(
  vroom_example("mtcars.csv"),
  col_types = list(gear = col_factor(levels = c(gear = c("3", "4", "5"))))
)

## -----------------------------------------------------------------------------
# vroom(
#   vroom_example("mtcars.csv"),
#   .name_repair = \(x) janitor::make_clean_names(x, case = "all_caps")
# )

## -----------------------------------------------------------------------------
vroom_write(mtcars, "mtcars.tsv")

## -----------------------------------------------------------------------------
unlink("mtcars.tsv")

## -----------------------------------------------------------------------------
vroom_write(mtcars, "mtcars.csv", delim = ",")

## -----------------------------------------------------------------------------
unlink("mtcars.csv")

## -----------------------------------------------------------------------------
vroom_write(mtcars, "mtcars.tsv.gz")

vroom_write(mtcars, "mtcars.tsv.bz2")

vroom_write(mtcars, "mtcars.tsv.xz")

## -----------------------------------------------------------------------------
unlink(c("mtcars.tsv.gz", "mtcars.tsv.bz2", "mtcars.tsv.xz"))

## -----------------------------------------------------------------------------
# vroom_write(mtcars, pipe("pigz > mtcars.tsv.gz"))

## -----------------------------------------------------------------------------
unlink("mtcars.tsv.gz")

