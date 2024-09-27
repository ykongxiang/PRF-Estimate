folder_path <- 'C:/Users/31598/Desktop/prfect/data/Atkins/gbk'
gbk_files <- list.files(folder_path, pattern = "\\.gbk$", full.names = TRUE)
gene_data_list <- list()

for (i in gbk_files){
  gene_data <- readLines(i)
  gene_data_list <- append(gene_data_list, list(gene_data)) 
}
