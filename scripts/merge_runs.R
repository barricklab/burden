##############  Read in the input files

#input.path = "_rates.all.csv_"
input.rates.file.names  = list.files(path = input.path)

readings = c("growth", "GFP") 

all.data = data.frame()
for (this.file.name in input.rates.file.names) {
  # this loop iterates over each input file name and reads them into
  # data frame this.data, then combines them into one data frame all.data
  
  this.file.path = file.path(input.path,this.file.name)
  this.data <- read.csv(this.file.path)
  
  if("graph" %in% colnames(this.data)) {
    this.data = this.data %>% filter(graph==1)
  }
  
  #we have to add other.rate if it is missing to join all data together
  if(!("other.rate" %in% colnames(this.data))) {
    this.data$other.rate=c(NA)
  }
  
  if(!("max.other.rate.time" %in% colnames(this.data))) {
    this.data$max.other.rate.time=c(NA)
  }
  
  this.data$run=this.file.name
  this.data$run = sub(".rates.all.csv", "", this.data$run)
  #this.data$run = sub("exp", "", this.data$run)
  #this.data$run = as.numeric(this.data$run)
    
  all.data = rbind(all.data, this.data)
}

############## Write out the merged file

write.csv(all.data, "rates.all.merged.csv", row.names=F)

