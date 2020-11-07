library(STRINGdb)

# get mouse STRING database 
cat("Started interactions_from_ensembl...\n")

args <- commandArgs(trailingOnly = TRUE)

group_id = args[1]
output_directory = args[2]
output_file <- paste(output_directory,group_id,'/ensembl_interactions.csv',sep='')

if(file.exists(output_file)){

    cat(paste("ensembl_interactions.csv already exists. Please delete existing",
         " file to run interactions_from_ensembl.\n",sep=""))

}else{

    cat("Downloading mouse STRINGdb interaction data.\n")
    string_db <- STRINGdb$new(version="10",species=10090,
                              score_threshold=0, input_directory="")
    cat("Finished downloading mouse STRINGdb interaction data.\n")

    input_file = paste(output_directory, group_id, '/all_expressed_genes.txt', sep='')
    fa1_clean <- read.csv(input_file, header = FALSE)
    names(fa1_clean) <- c('ens_id')
    fa1_levels<- fa1_clean$ens_id
    fa1_hits <- as.character(fa1_levels)
    interactions <- string_db$get_interactions(fa1_hits)
    write.csv(interactions,file=output_file, 
              row.names = FALSE) 
    cat("interactions_from_ensembl completed!\n") 
}


