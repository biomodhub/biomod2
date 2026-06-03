library(rangeModelMetadata);

#Location for .csv templates
#dir <- "~/Dropbox/RMMmanuscript/";
dir='/Users/ctg/Dropbox/Projects/Range_Metadata/Supplements_for_R1'
#Building a base template
baseTemplate <- rmmTemplate(family = "base");
rmmToCSV(baseTemplate, file= paste(dir, "rmm_base.csv", sep = ""));

#Building the full deluxe template
fullTemplate <- rmmTemplate();
rmmToCSV(fullTemplate, file= paste(dir, "rmm_full.csv", sep = ""));
