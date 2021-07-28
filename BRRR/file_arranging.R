## For re-arranging the spectra files
data_path <- paste0(getwd(), "/data/subjects/")

subjects <- list.files(data_path)
#ages = read.csv('data/ages.csv')
#ages <- ages[,-1]


for(i in 1:length(subjects)){
  subj_dir = paste0(data_path, subjects[i], sep="/")
  files = list.files(subj_dir)
  # Rename files to unique per subject
  file.rename(paste0(subj_dir, files[1]), 
              paste0(subj_dir, paste0(subjects[i], "_n1.csv")))
  file.rename(paste0(subj_dir, files[2]), 
              paste0(subj_dir, paste0(subjects[i], "_n2.csv")))
  files2 = list.files(subj_dir)
  
  file.copy(paste0(subj_dir, files2[1]), "data/n1_sleep/")
  file.copy(paste0(subj_dir, files2[2]), "data/n2_sleep/")
}


