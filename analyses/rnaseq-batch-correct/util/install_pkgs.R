# Function: Install required packages

# check if packages are available if not download 
if (!("tidyverse" %in% installed.packages())){
  install.packages("tidyverse")
}
if (!("optparse" %in% installed.packages())){
  install.packages("optparse")
}
if (!("Rtsne" %in% installed.packages())){
  install.packages("Rtsne")
}
if (!("ggpubr" %in% installed.packages())){
  install.packages("ggpubr")
}
if (!("plyr" %in% installed.packages())){
  install.packages("plyr")
}
if (!("dplyr" %in% installed.packages())){
  install.packages("dplyr")
}
if (!("data.table" %in% installed.packages())){
  install.packages("data.table")
}
if (!("reshape2" %in% installed.packages())){
  install.packages("reshape2")
}
if (!("sva" %in% installed.packages())){
  install.packages("sva")
}
if (!("grid" %in% installed.packages())){
  install.packages("grid")
}
if (!("ggthemes" %in% installed.packages())){
  install.packages("ggthemes")
}
if (!("scales" %in% installed.packages())){
  install.packages("scales")
}

