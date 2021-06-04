#t1 <- upload.file("Unnormalized peak areas/Unnormalized all protein-peptides No Dup.csv", 1, 2:88)

presence.absence <- function(input_file){
  library(tidyr)
  library(dplyr)
  t1 <- input_file %>%
    na_if(0)
  pheno_groups <- unique(t1$Label); group_num <- length(pheno_groups)

  for (i in 1:group_num){
    t2 <- t1 %>%
      filter(Label == pheno_groups[i])

    t1[paste0('missing', pheno_groups[i]),] <- as.matrix(apply(t2, 2, function(x) (sum(is.na(x)))/dim(t2)[1]))
    t1 <- rbind.data.frame(slice(t1[paste0('missing', pheno_groups[i]),]),t1[-dim(t1)[1],])
  }
  t1 <- as.data.frame(t(t1[,-1]))
  keep_file <- data.frame()
  reject_file <- data.frame()
#  write(reject_file, "absentFeatures.txt", append = TRUE)
  for (i in 1:dim(t1)[1]){
    if( mean(as.numeric(unlist(t1[i,1:group_num]))) < .50 ){
      keep_file <- rbind.data.frame(t1[i,], keep_file)
    }
    else {
      reject_file <- rbind.data.frame(t1[i,], reject_file)
    }
  }
  return(keep_file)
}
