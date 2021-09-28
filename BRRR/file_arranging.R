## For re-arranging the spectra files
setwd("/projects/FABEEG/BRRR")
data_path <- paste0(getwd(), "/data/absolute_spectra/")

subjects <- list.files(data_path)
#ages = read.csv('data/ages.csv')
#ages <- ages[,-1]


for(i in 2:length(subjects)){
  subj_dir = paste0(data_path, subjects[i], sep="/")
  files = list.files(subj_dir)
  
  for(file in files){
    # Rename files to be unique per subject
    file.rename(paste0(subj_dir, file),
                paste0(subj_dir, paste0(subjects[i], file)))
  }
  files2 = list.files(subj_dir) #TODO: change order (now n2e = 1st segment)
  dirs = c("n1_sleep/", "n2a_sleep/", "n2b_sleep/", "n2c_sleep/", "n2d_sleep/", "n2e_sleep/")
  for(f in 1:length(files2)){
    newdir=dirs[f]
    file2=files2[f]
    file.copy(paste0(subj_dir, file2), paste0("data/",newdir))
  }
}


