#' This funciton takes the file location and convert raw data into CNV class for downstream analysis.
#'
#' @param fileloc
#' @return a CNV class containing raw data information
#' @export load_data_bis
#' @export load_data

load_data_bis <- function(fileloc){
  library(data.table)
  library(dplyr)
  if (!grepl("/$", fileloc)) fileloc = paste0(fileloc, "/")
  file_list <- list()
  fi <- list.files(fileloc, pattern = "*.txt|*.csv")
  for (file in fi){
    sep = ifelse(grepl(".txt$", file), "\t", ",")
    file_list[[length(file_list)+1]] <- fread(paste0(fileloc, file), skip = 1, sep = sep, header = T) %>%
      .[, plate_ID := NA] %>%
      setnames("Sample Name", "Sample.Name") %>%
      setnames("Target Name", "Target.Name")
  }

  # plate_ID = vector()
  #   for(i in 1:length(fi)){
  #     file_list[[i]]$plate_ID = strsplit(fi[i], split = "-")[[1]][3]
  #     attr(file_list[[i]], "plate_ID") = strsplit(fi[i], split = "-")[[1]][3]
  #     attr(file_list[[i]], "plate_ID") -> plate_ID[i]
  #   }


  data_raw = do.call(rbind, file_list) %>% tbl_df
  data = list()
  data$data_raw = data_raw
  data$file_list = file_list
  data
}


load_data <- function(fileloc){
  file_list <- list()
  fi <- list.files(fileloc, pattern = "*.csv")
  for (file in fi){
    file_list[[length(file_list)+1]] <- read.table(paste0(fileloc, file),skip = 1, sep = ",", header = T)
  }

  plate_ID = vector()
  for(i in 1:length(fi)){
    file_list[[i]]$plate_ID = strsplit(fi[i], split = "-")[[1]][3]
    attr(file_list[[i]], "plate_ID") = strsplit(fi[i], split = "-")[[1]][3]
    attr(file_list[[i]], "plate_ID") -> plate_ID[i]
  }

  data_raw = do.call(rbind, file_list) %>% tbl_df

  data = list()
  data$data_raw = data_raw
  data$file_list = file_list
  data
}
