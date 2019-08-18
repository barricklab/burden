#!/usr/bin/env Rscript

##############################################################
#
# Using RStudio to run this code? Follow these steps:
#
# 1) Define the variable 'input.prefix' in your Rstudio shell 
#    and then all of the command-line option code will be 
#    skipped and you can run this script within RStudio!
#
#    Example command
#      input.prefix = "my_exp"
#
# 2) Set the working directory to a folder on your computer
#    that contains comma separated value (CSV) input files. 
#    It should have a "measurements file saved from the 
#    platereader and metadata file that you created that
#    describes the samples (see full docs for format)
#    beginning with your prefix:
#    
#       'input.prefix.measurements.csv'
#       'input.prefix.metadata.csv'
#
#     You can also use tab-separated value (TSV) files that end in *.tsv
#
# 3) Set any of these variables ahead of time to override the defaults.
#
#    Example commands:
#      max.method = 1
#      output.prefix = "my_output"
#
##############################################################

library(tidyverse)
library(gridExtra)
library(cowplot)

##############################################################
#### Load command line options and set global parameters
##############################################################


## Defaults:


if (!exists("input.prefix")) {
  require(optparse)
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Input file prefix. Expects to find the files <input>.metadata.tsv and <input>.measurement.tsv", metavar="input.csv"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output file prefix. Output files of the form <output>.* will be created. If this option is not provided the input file prefix will be used.", metavar="output_prefix"),
    make_option(c("-m", "--max-method"), type="character", default=2, 
                help="Maximum picking method. 1=pick max growth rate time and use other rates at this same time. 2=pick maximum value of each curve at whatever time it occurs (may be different for growth rate and for fluorescence)", metavar="output_prefix")
    
    #TODO: We need to make more options accessible at the command line
  )
  
  usage_string = paste(
    "burden.R -i input -o output\n\n",
    sep = ""
  ) 
  
  opt_parser = OptionParser(usage=usage_string, option_list=option_list);
  opt = parse_args(opt_parser);
  
  if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("You must supply the -i|--input argument for the input", call.=FALSE)
  }
  input.prefix = opt$input
  
  if (!is.null(opt$output)) {
    output.prefix = opt$output
  }
  
  if (!is.null(opt$max.method)) {
    max.method = opt$max.method 
  }
  
} else {
  # RStudio mode 
  # You have manually assigned input.prefix to reach this block
  # Nothing else is necessary...
}


# Set values of these variables to command-line values if provided.
# "option.____" variable names are the ones that should be used throughout
# the code below in case the settings are changed within RStudio
# between runs of the code

option.input.prefix = input.prefix

# Default to the same output prefix as input prefix. 
if (exists("output.prefix")) {
  option.output.prefix = output.prefix
} else {
  option.output.prefix = input.prefix
}

#Default to max.method 2 (for now)
if (exists("max.method")) {
  option.max.method = max.method
} else {
  option.max.method = 2
}

##############################################################
#### Advanced options that are currently hard-coded
##############################################################

#TODO: Add an option to have a third reading (e.g., of BFP or luminescence)
readings = c("OD", "GFP") 

# Fix all initial measurements to be equal across all non-blank samples in a category
offset.readings.to.average = T
offset.readings.to.average.time.span = c(0,60)

minimum.OD = 0.0
maximum.time = 1000000
maximum.time = 250

# time.point.span should be an odd number 
# time.point.delta is the number of points forward and back to use in the window for determining rates
time.point.span = 3
time.point.delta = floor(time.point.span/2)


##############################################################
#### Print the Settings
##############################################################

#Create a dataframe of settings so that we can write them out
# For sdf (settings data frame), the columns are ["name", "value"]

add_setting <- function(sdf, n, v) {
  sdf = suppressWarnings(bind_rows(sdf, data.frame(name=as.character(n), value=as.character(v))))
  return(sdf)
}

print_settings_df <- function(sdf) {
  for (i in 1:nrow(sdf)) {
    cat(paste0(sdf$name[i], " = ", sdf$value[i],"\n"))
  }
}

# If possible (script has not been moved), get the repository/version of the code
# and output this as a setting so that we can know which version was used!
#
# @JEB: This code finds the path of the current script so that we can offset to the 
# Git file containing the currently checked out version. This is a pretty complex
# thing to do and I have not tested this across platforms. It may need more work
# so that we can allow it to gracefully fail w/o erroring out of the script.

this_file_path <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    #cat("Rscript\n")
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    # 'source'd via R console
    #cat("Sourced\n")
    return(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    # Rstudio
    #cat("RStudio\n")
    return(rstudioapi::getActiveDocumentContext()$path)
  }
}

#Returns data frame with 'version' and 'repository' fields
git_repo_version <- function() {
  git.file.path = file.path(dirname(this_file_path()), "..", ".git", "FETCH_HEAD")
  #cat(git.file.path, "\n")
  if(file.exists(git.file.path)) {
    mystring <- read_file(git.file.path)
    mystring = str_remove_all(mystring, "[\r\n]")
    mystringlist = strsplit(mystring, "\t\t")
    grv = c()
    grv$version = mystringlist[[1]][1]
    grv$repository = mystringlist[[1]][2]
    return(grv)
  } else {
    grv = c()
    grv$version = "unknown"
    grv$version = "unknown"
    return(grv)
  }
}
git.repo.version.info = git_repo_version()

sdf = data.frame()
sdf = add_setting(sdf, "Run At", Sys.time())
sdf = add_setting(sdf, "Repository", git.repo.version.info$repository)
sdf = add_setting(sdf, "Version", git.repo.version.info$version)
sdf = add_setting(sdf, "Input Prefix", option.input.prefix)
sdf = add_setting(sdf, "Output Prefix", option.output.prefix)

# Translate the method to something meaningful
option.max.method.string = "unknown"
if (option.max.method == 1) {
  option.max.method.string = "1 (Find time point with maximum growth rate. Report all rates at this time point.) [one time = max growth rate]"
} else if (option.max.method == 2) {
  option.max.method.string = "1 (Find time point with maximum rate and report that for each curve separately.) [multiple times = max of each curve]"  
}
sdf = add_setting(sdf, "Maximum Picking Method", option.max.method.string)

sdf = add_setting(sdf, "Readings", paste(readings))
sdf = add_setting(sdf, "Offset Readings to Average", offset.readings.to.average)
sdf = add_setting(sdf, "Offset Readings to Average Time Span", 
                  paste(offset.readings.to.average.time.span, collapse="-"))
sdf = add_setting(sdf, "Minimum OD", minimum.OD)
sdf = add_setting(sdf, "Maximum Time", maximum.time)
sdf = add_setting(sdf, "Time Point Span", time.point.span)
sdf = add_setting(sdf, "Time Point Delta", time.point.delta)

cat("=== SETTINGS ==\n")
print_settings_df(sdf)
write.csv(sdf, paste0(option.output.prefix, ".settings.csv"), row.names=F)

#TODO: Get the github revision from the .git directory and add this information!

##############################################################
#### Create the plot directory
##############################################################
plot.directory = paste0(option.output.prefix, "-plots")
dir.create(file.path(plot.directory), showWarnings = FALSE)

##############################################################
#### Load measurement file and tidy
##############################################################

#SPL Note: had to run "dos2unix" to clean up the file from our platereader before it could be read properly. Not sure how to fix this on the windows end.
#TODO: Check for incorrect line endings and warn user

all_data = data.frame()

if (file.exists(paste0(option.input.prefix, ".measurements.csv"))) {
  all_data <- read_csv(paste0(option.input.prefix, ".measurements.csv"), col_names=F, comment = "#" )
} else if (file.exists(paste0(option.input.prefix, ".measurements.tsv"))) {
  all_data <- read_tsv(paste0(option.input.prefix, ".measurements.tsv"), col_names=F, comment = "#" )
} else {
  stop(paste0("Could not find a valid measurements file. Tried:\n", paste0(option.input.prefix, ".measurements.csv"), "\n", paste0(option.input.prefix, ".measurements.tsv")))
}

# If there is no header row of wells, assume standard order by row from top of plate.
# Otherwise, take the first row as the well labels.
if (all_data$X1[1] == "0s") {
  names(all_data) = c("Time",
                      "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
                      "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
                      "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
                      "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
                      "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
                      "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
                      "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
                      "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"
                      )
    
    
    
} else {
  names(all_data) = all_data[1,]
  all_data = all_data[-1,]
}

#Remove row thats are NA in time (this is the trialing trash data). 
#Remove "s" and convert to integer, remove row that dont convert properly and leave NA

all_data$Time <- str_replace(all_data$Time, "s", "")
all_data$Time <- as.integer(all_data$Time)
all_data <- na.omit (all_data)

#Separate data into data frames to make tidy

#Find the zero times? (they separate the different types of readings)
zeros <- which(all_data$Time == "0")
zeros = c(zeros, nrow(all_data)+1)

#generic names for the readings
reading.names = c("OD", "GFP", "other")

tidy_all = data.frame()

## Need to keep all times in lock step despite some rounding
## so we use the first chunk times
definitive_times = c()

for(i in 1:(length(zeros)-1)) {
  
  this.data.chunk = all_data[zeros[i]:(zeros[i+1]-1),]
  tidy.data.chunk <- gather(this.data.chunk, key = "well", value = "value", -Time)
  tidy.data.chunk$reading <- reading.names[i]

  if (i==1) {
    definitive_times = floor(tidy.data.chunk$Time/60)
  }
  
  ## If the file was truncated - readings stopped in the middle of a row by the user,
  ## Then later chunks may have fewer values. In this case, we need to shorten the earlier
  ## ones and the number of definitive times.
  if (nrow(tidy.data.chunk) < nrow(tidy_all)) {
    tidy_all = tidy_all[1:nrow(tidy.data.chunk),]
    definitive_times = definitive_times[1:nrow(tidy.data.chunk)]
  }
  
  tidy.data.chunk$time.min = definitive_times
  tidy.data.chunk = tidy.data.chunk %>% select (-Time)
  
  tidy_all = tidy_all %>% bind_rows(tidy.data.chunk)
}

# Remove overflow readings
tidy_all = tidy_all %>% filter(value != "Overflow") 


#Write the tidy data
write_csv( tidy_all, paste0(option.output.prefix, ".tidy.measurements.csv"))

##############################################################
#### Load metadata file and tidy
##############################################################

metadata <- data.frame()

if (file.exists(paste0(option.input.prefix, ".metadata.csv"))) {
  metadata <- read_csv(paste0(option.input.prefix, ".metadata.csv"), col_names=T, comment = "#" )
} else if (file.exists(paste0(option.input.prefix, ".metadata.tsv"))) {
  metadata <- read_tsv(paste0(option.input.prefix, ".metadata.tsv"), col_names=T, comment = "#" )
} else {
  stop(paste0("Could not find a valid measurements file. Tried:\n", paste0(option.input.prefix, ".measurements.csv"), "\n", paste0(option.input.prefix, ".measurements.tsv")))
}

## Check for required columns
if(! ("well" %in% colnames(metadata)) ) {
  stop("Metadata does not have 'well' column (case-sensitive)", call.=FALSE)
}
if(! ("strain" %in% colnames(metadata)) ) {
  stop("Metadata does not have 'strain' column (case-sensitive)", call.=FALSE)
}

#Add optional columns if missing
if(! ("include" %in% colnames(metadata)) ) {
  metadata$include = T
}
if(! ("isolate" %in% colnames(metadata)) ) {
  metadata$isolate = NA
}
if(! ("description" %in% colnames(metadata)) ) {
  metadata$description = ""
}

#convert boolean fields
metadata$include = as.logical(metadata$include)


## Determine what wells we are ignoring and remove from metadata
ignore.wells = metadata$well[metadata$include==F]
metadata = metadata %>% filter(include==T)

# Fix the names of strains that are synonyms for blank (case-insensitive)
metadata$strain.upper = toupper(metadata$strain)
metadata$strain[grepl("^(B|BLANK)$", metadata$strain.upper, perl=T)] = "blank"
metadata = metadata %>% select(-strain.upper)

# Blanks can't have isolate numbers
metadata$isolate[metadata$strain=="blank"] = NA

metadata$strain.isolate = paste0(metadata$strain, "__", metadata$isolate)

#Fix wells that have blank back to blank
metadata$strain.isolate = sub("blank__NA", "blank", metadata$strain.isolate)

write_csv( metadata, paste0(option.output.prefix, ".tidy.metadata.csv"))


##############################################################
#### Combine tidy measurement and metadata files
##############################################################

#Now,  load the data file and wrangle it
X = tidy_all
X = X%>% filter(time.min < maximum.time)
X = X %>% filter(!(well %in% ignore.wells))
X$value=as.numeric(X$value)
Y = X %>% spread(reading, value)

#join sample data to master tidy data
Y <- full_join(metadata, Y, by = "well")
Y$well = factor(Y$well)
Y = Y %>% filter(time.min != 0)


# Analysis starts here ------->
##############################################################
#### Background Correction
#### * Subtracts the average of ALL blank wells at a given time from readings
##############################################################

Y$GFP = as.numeric(Y$GFP)
Y$OD = as.numeric(Y$OD)

BG = Y %>% filter(strain=="blank")

BG$well = droplevels(BG$well)

plot1 = ggplot(BG, aes(time.min, OD, color=well)) + geom_point() 
plot2 = ggplot(BG, aes(time.min, GFP, color=well)) + geom_point() 

if (length(readings) == 2) {
  p = grid.arrange(plot1, plot2, nrow=2)
} else if (length(readings) == 3) {
  plot3 = ggplot(BG, aes(time.min, BFP, color=well)) + geom_point() 
  p = grid.arrange(plot1, plot2, plot3, nrow=3)
}

ggsave(file.path(plot.directory, "background.pdf"), plot=p)

Z = Y
if (length(readings) == 3) {
  
  BG.means = BG %>% group_by(time.min) %>% summarize(other.bg = mean(other), GFP.bg = mean(GFP), OD.bg = mean(OD))
  
  Z = Z %>% left_join(BG.means, by="time.min")
  
  Z$OD = Z$OD - Z$OD.bg
  Z$GFP = Z$GFP - Z$GFP.bg
  Z$other = Z$other - Z$other.bg
  
} else if (length(readings) == 2) {
  
  BG.means = BG %>% group_by(time.min) %>% summarize(GFP.bg = mean(GFP), OD.bg = mean(OD))
  
  Z = Z %>% left_join(BG.means, by="time.min")
  
  Z$OD = Z$OD - Z$OD.bg
  Z$GFP = Z$GFP - Z$GFP.bg

}
#Ditch the blanks here
Z = Z %>% filter(strain != "blank")

#For debugging
#write_tsv(Z, "temp.tsv")

##############################################################
#### Average Per-Sample Offset Correction
#### * Equalizes the beginning of each curve for a given strain
####   This is a kludge to correct for well-to-well background differences
####   TODO:  figure out a way to experimentally measure these differences!
##############################################################
if (offset.readings.to.average) {
  
  cat("Offsetting to average measured values\n")
  
  #calculate the average over all times for a sample
  if (length(readings) == 2) {
    
    ZZ = Z %>% filter((time.min >= offset.readings.to.average.time.span[1]) & (time.min <= offset.readings.to.average.time.span[2]))
    
    strain.means = ZZ %>% group_by(strain.isolate) %>% summarize(mean.strain.OD = mean(OD), mean.strain.GFP = mean(GFP))
    well.means = ZZ %>% group_by(well, strain.isolate) %>% summarize(mean.well.OD = mean(OD), mean.well.GFP = mean(GFP))
  
    well.means = well.means %>% left_join(strain.means, by="strain.isolate")
    
    well.means$well.OD.offset = well.means$mean.strain.OD - well.means$mean.well.OD
    well.means$well.GFP.offset = well.means$mean.strain.GFP - well.means$mean.well.GFP
    
    well.means = well.means %>% select(well, well.OD.offset, well.GFP.offset)
    
    Z = Z %>% left_join(well.means, by="well")
    
    Z$OD = Z$OD + Z$well.OD.offset
    Z$GFP = Z$GFP + Z$well.GFP.offset
    
    
  } else if (length(readings) == 3) {
    
    #stop("Not implemented")
    #Z %>% group_by(strain, reading) %>% summarize(mean.sample.OD = mean(OD), mean.sample.GFP = mean(GFP), mean.sample.other = mean(other))
  }
}


##############################################################
#### Calculate and graph rates per strain
##############################################################


Z = Z %>% filter (OD >= minimum.OD)
final.table = data.frame()
for (strain.of.interest in unique(Z$strain.isolate) )
{
  cat("STRAIN:", strain.of.interest, "\n") 
  strain.data = Z %>% filter(strain.isolate==strain.of.interest)
  head(strain.data)
  
  ## graph
  
  plot1 = ggplot(strain.data, aes(time.min, OD, color=well)) + geom_point() 
  plot2 = ggplot(strain.data, aes(time.min, GFP, color=well)) + geom_point() 
  
  if (length(readings) == 2) {
    p = grid.arrange(plot1, plot2, nrow=2)
  } else if (length(readings) == 3) {
    plot3 = ggplot(strain.data, aes(time.min, other, color=well)) + geom_point() 
    p = grid.arrange(plot1, plot2, plot3, nrow=3)
  }
  
  ggsave(file.path(plot.directory, paste0(strain.of.interest, ".pdf")), plot=p)
  
  
  strain.data$logOD = log(strain.data$OD)
  #ggplot(strain.data, aes(time.min, logOD, color=well)) + scale_x_continuous(limits = c(0, 400))+ geom_point() 
  
  #Calculate the growth rate over a given time period and graph
  growth.rate.data = strain.data
  
  growth.rate.data = growth.rate.data %>% 
    group_by(well) %>% 
    mutate(t1 = lag(time.min, time.point.delta, order_by=well)) %>% 
    mutate( t2 = lead(time.min, time.point.delta, order_by=well)) %>% 
    mutate( logOD.t1 = lag(logOD, time.point.delta, order_by=well)) %>% 
    mutate( logOD.t2 = lead(logOD, time.point.delta, order_by=well))
  
  growth.rate.data = growth.rate.data %>% filter(!is.na(t1) & !is.na(t2))
  
  growth.rate.data = growth.rate.data %>% mutate(time.min = (t1+t2)/2)
  growth.rate.data = growth.rate.data %>% mutate(specific.growth.rate = 60 * (logOD.t2 - logOD.t1) /  (t2 - t1))
  
  plot1 = ggplot(growth.rate.data , aes(time.min, specific.growth.rate, color=well)) +  geom_point() 

  fluorescence.data = strain.data
  
  fluorescence.data = fluorescence.data %>% 
    group_by(well) %>% 
    mutate(t1 = lag(time.min, time.point.delta, order_by=well)) %>% 
    mutate( t2 = lead(time.min, time.point.delta, order_by=well)) %>% 
    mutate( GFP.t1 = lag(GFP, time.point.delta, order_by=well)) %>% 
    mutate( GFP.t2 = lead(GFP, time.point.delta, order_by=well)) 
  
  if (length(readings) == 3) {
    fluorescence.data = fluorescence.data %>% 
      mutate( other.t1 = lag(other, time.point.delta, order_by=well)) %>% 
      mutate( other.t2 = lead(other, time.point.delta, order_by=well))
  }
  
  fluorescence.data = fluorescence.data %>% filter(!is.na(t1) & !is.na(t2))
  fluorescence.data = fluorescence.data %>% mutate(time.min = (t1+t2)/2)
  fluorescence.data = fluorescence.data %>% mutate(GFP.rate = (GFP.t2 - GFP.t1) / (t2 - t1) / OD)
  
  if (length(readings) == 3) {
    fluorescence.data = fluorescence.data %>% mutate(other.rate = (other.t2 - other.t1) / (t2 - t1) / OD)
    
  }
  
  plot2 = ggplot(fluorescence.data , aes(time.min, GFP.rate, color=well)) +  geom_point() 
  if (length(readings) == 2) {
    p = grid.arrange(plot1, plot2, nrow=2)
  } else if (length(readings) == 3) {
    fluorescence.data = fluorescence.data %>% mutate(other.rate = (other.t2 - other.t1) / (t2 - t1) / OD)
    plot3 = ggplot(fluorescence.data , aes(time.min, other.rate, color=well)) +  geom_point() 
    p = grid.arrange(plot1, plot2, plot3, nrow=3)
  }
  
  ggsave(file.path(plot.directory, paste0(strain.of.interest, ".rates.pdf")), plot=p)
  
  
  #ggplot(fluorescence.data , aes(time.min, GFP.rate, color=well)) + scale_x_continuous(limits = c(0, 400)) +  geom_point() 
  
  #ggplot(fluorescence.data , aes(time.min, BFP.rate, color=well)) + scale_x_continuous(limits = c(0, 400)) +  geom_point() 

  growth.rate.data$well = droplevels(growth.rate.data$well)
  
  strain.max.values = data.frame()
  
  
  cat("STRAIN__ISOLATE:", strain.of.interest, "\n\n")
  
  for (this.well in levels(growth.rate.data$well)) {
    
    cat("Well:", this.well, "\n")
    
    replicate.growth.rate.data = growth.rate.data %>% filter(well == this.well)
    
    #identify the best growth rate row/time
    max.growth.rate.row = replicate.growth.rate.data[which.max(replicate.growth.rate.data$specific.growth.rate),]
    
    cat("Time of maximum growth rate:", max.growth.rate.row$time.min[1], "min\n")
    
    cat("Max growth rate:", max.growth.rate.row$specific.growth.rate, "per hour\n")
    
    if (option.max.method==1) {
      max.GFP.fluorescence.data.row = fluorescence.data %>% filter(well == this.well) %>% filter(time.min == max.growth.rate.row$time.min[1])
    } else {
      
      replicate.fluorescence.data = fluorescence.data %>% filter(well == this.well)
      
      max.GFP.fluorescence.data.row = replicate.fluorescence.data[which.max(replicate.fluorescence.data$GFP.rate),]
    }
    
    cat("Time of GFP production rate:", max.GFP.fluorescence.data.row$time.min[1], "min\n")

    cat("Max GFP production rate:", max.GFP.fluorescence.data.row$GFP.rate, "per hour\n")
    

    
    if (length(readings) == 3) {
      
      if (option.max.method==1) {
        max.other.fluorescence.data.row = fluorescence.data %>% filter(well == this.well) %>% filter(time.min == max.growth.rate.row$time.min[1]) 
      } else {
        
        max.other.fluorescence.data.row = replicate.fluorescence.data[which.max(replicate.fluorescence.data$other.rate),]
        
        max.GFP.fluorescence.data.row = replicate.fluorescence.data[which.max(replicate.fluorescence.data$other.rate),]
      }
      
      cat("Max BPF production:", max.other.fluorescence.data.row$other.rate, "per hour\n")
    }
    
  if (length(readings) == 2) {
      
    strain.max.values = bind_rows(strain.max.values, 
                                  data.frame(
                                    growth.rate = max.growth.rate.row$specific.growth.rate,
                                    GFP.rate = max.GFP.fluorescence.data.row$GFP.rate
                                  )
    )
    
   } else  if (length(readings) == 3) {
      
    strain.max.values = bind_rows(strain.max.values, 
                                  data.frame(
                                    growth.rate = max.growth.rate.row$specific.growth.rate,
                                    GFP.rate = max.GFP.fluorescence.data.row$GFP.rate,
                                    other.rate = max.other.fluorescence.data.row$other.rate
                                  )
    )
    }
  }
  
  cat("\n\n")
  cat("Growth rate:", mean(strain.max.values$growth.rate), "±", sd(strain.max.values$growth.rate), "\n")
  cat("GFP rate:", mean(strain.max.values$GFP.rate), "±", sd(strain.max.values$GFP.rate), "\n")
  
  if (length(readings) == 3) {
    cat("other rate:", mean(strain.max.values$other.rate), "±", sd(strain.max.values$other.rate), "\n")
  }
  cat("\n\n")
  
  
  if (length(readings) == 2) {
    
    final.table = bind_rows(final.table, 
                            data.frame(
                              strain.isolate = strain.of.interest,
                              growth.rate = mean(strain.max.values$growth.rate),
                              growth.rate.sd = sd(strain.max.values$growth.rate),
                              GFP.rate = mean(max.GFP.fluorescence.data.row$GFP.rate),
                              GFP.rate.sd = sd(strain.max.values$GFP.rate)
                            )
    )
    
  } else  if (length(readings) == 3) {
    
    final.table = bind_rows(final.table, 
                            data.frame(
                              strain.isolate = strain.of.interest,
                              growth.rate = mean(strain.max.values$growth.rate),
                              growth.rate.sd = sd(strain.max.values$growth.rate),
                              GFP.rate = mean(max.GFP.fluorescence.data.row$GFP.rate),
                              GFP.rate.sd = sd(strain.max.values$GFP.rate),
                              other.rate = mean(max.other.fluorescence.data.row$other.rate),
                              other.rate.sd = sd(strain.max.values$other.rate)
                            )
    )
  }
  

}


##Split out strain and replicate
final.table$isolate= sub("^.+__", "", final.table$strain.isolate, perl = T)
final.table$strain = sub("__.+$", "", final.table$strain.isolate, perl = T)
final.table = final.table %>% select(-strain.isolate)

write_csv(final.table, paste0(option.output.prefix, ".rates.summary.csv"))

# Need to fix output of fit in each well
#write_csv(strain.max.values, paste0(dataset.name, ".rates.all.csv"))
          
