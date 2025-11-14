work_dir <- "/Volumes/albertine/10_filter/"

filelist <- list.files(work_dir, pattern = "*.stats", full.names = T)

stats_df <- NULL

for(i in 1:length(filelist)){
  x <- read_lines(filelist[i])
  sp_name <- strsplit(x[1], "__")[[1]][1]
  sp_name <- gsub("_", " ", sp_name)
  sample_num <- strsplit(x[1], "__")[[1]][2]
  sample_num <- gsub("_", "", sample_num)
  sum_depth <- strsplit(x[3], "=  ")[[1]][2]
  prop_dup <- x[5]
  site_geno <- x[7]
  record <- cbind(sp_name, sample_num, sum_depth, prop_dup, site_geno)
  if(i == 1){
    stats_df <- record
  }else{
    stats_df <- rbind(stats_df, record)
  }
}

stats_df |> 
  as_data_frame() |> 
  filter(sp_name == "Batis diops" | sp_name == "Chamaetylas poliophrys" | sp_name == "Cossypha archeri" | sp_name == "Phylloscopus laetus") |> 
  write_csv("~/Desktop/summary_stats_albertine.csv")
