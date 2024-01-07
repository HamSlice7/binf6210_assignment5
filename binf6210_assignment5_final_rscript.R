#LFQ automated end-to-end proteomics pipeline developed by Jacob Hambly. 
#Github:https://github.com/HamSlice7/binf_6210_assignment5
#The example input files for the proteinGroup.txt file and the experimental design file can be found on the Github page. 
#example_pg.txt is the proteinGroup.txt file for the example run
#example_ed.txt is the experimental design file for the example run 
#More information on the pipeline can be found under the README.md file on Github

#-----------------------------------------------------------------------------------------------

#Installing if not previously done so and loading in the necessary packages
#install.packages('promor')
#install.packages('tidyverse')
library('promor')
library('tidyverse')

#Pipeline function
LFQ_proteomics_pipeline <- function() { 
#Prompting the user to enter the name of the proteinGroups.txt file they are using for their analysis. File must be located in the users current working directory
  MQ_PG_File <- readline(prompt = "Enter the file name for your MaxQuant proteinGroup.txt file (include .txt): ")
#Prompting the user to enter the name of the experimental design file they will be using for their analysis. File must be located in the users current working directory.
  ED_File <- readline(prompt = "Enter the file name for your experimental design file (include .txt): " )
  
  cat("proteinGroups.txt file and experimental design file sucessfully read in", '\n')
  
#reading in the raw proteinGroups.txt file
  inital_data <- read.table(MQ_PG_File, sep = '\t', header = T)
  
  
#Initial number of rows for the entered proteinGroup.txt file
  cat("The entered proteinGroup.txt file contains:", nrow(inital_data), "rows", "\n")
 
#Creating a data frame called msdata which will contain the label free quantification (LFQ) intensities of each the technical replicate. The experimental design file will be read to separate each technical replicate into a experimental group for easier downstream analysis.Empty rows and columns are removes as well as reverse peptides, proteins identified by site, and potential contaminations. Next, missing values are converted to Na for downstream imputation. Lastly, the intensity values are log transformed. 

   msdata <- create_df(
    MQ_PG_File, 
    ED_File, 
    input_type = "MaxQuant",
    data_type = "LFQ",
    filter_na = TRUE,
    filter_prot = TRUE,
    uniq_pep = 0, 
    tech_reps = FALSE,
    zero_na = TRUE,
    log_tr = TRUE,
    base = 2
  )
  
#summary statistics of each column which represents LFQ intensities in each of the technical replicates. Notice number of NA's.
  print(summary(msdata))
  cat("The number of rows filtered out after removing empty rows, removing potential contaiminations, removing reverse peptides, and removing proteins identified by site is:", (nrow(inital_data) - nrow(msdata)), "\n")
  cat("Number of rows remaining:", nrow(msdata))

#creating a for loop that appends each column name to an empty vector called sample_name. Another for loop will determine the number of protein groups identified in each column and appends these numbers to an empty vector called sample_num_p.

#creating an empty vector to store the names of each technical replicates
  sample_names <- c()

#looping through the column names of each column in msdata and appending the names to sample_names
  for (i in colnames(msdata)) {
    sample_names <- append(sample_names, i)
  
  }

#creating an empty vector to store the number of proteins in each LFQ column
  sample_num_p = c()
#for loop that loops through each LFQ column and adds the sum of proteins to the lfq_num_p vector.
  for (i in sample_names) {
    sample_num_p <- append(sample_num_p, sum(!is.na(msdata[[i]])))
  }

#Creating a data frame by combining the sample_names and sample_num_p vectors. A csv file containing this data frame is written in the current working directory along with a barplot showing the number of proteins identified in each technical replicate. The number of proteins identified should be around the same in each technical replicate.
  number_of_proteins <- data.frame(Name = sample_names,
                                 Total = sample_num_p )

  write.csv(number_of_proteins, "numberofproteins.csv", row.names = F)
  
  num_p_bar_plot <- ggplot(number_of_proteins, aes(x = as.factor(sample_names), y = sample_num_p, fill = as.factor(sample_names) )) +
    geom_col() +
    scale_fill_brewer(palette = "Dark2") +
    theme(legend.position = "none") +
    labs(title = "Number of Protein Groups in each Technical Replicate", x = "Replicates", y = "Count"  )

  print(num_p_bar_plot)
  ggsave("proteingroups.pdf", plot = num_p_bar_plot, width = 7, height = 7)
  cat("Number of protein groups per technical replicate successfully calculated and saved in current working directory", '\n')
  
  

#The filterbygroup_na function will filter out rows that contain more than 75% missing data in either of the experimental groups for quality control. 
  
  ms_data_filter_by_group <- filterbygroup_na(msdata, set_na = 0.75, filter_condition = "either")
  cat("Number of rows after filtering proteins by group level missing data is: ",(nrow(ms_data_filter_by_group)), "\n")

#The impute_na function will impute any missing data using the miniProb method which imputes data from a normal distribution column wise. Imputation is needed for the downstream statisical analysis. 
  ms_data_imputed <- impute_na(
    ms_data_filter_by_group,
    method = "minProb",
    tune_sigma = 1,
    q = 0.01,
    maxiter = 10,
    ntree = 20,
    n_pcs = 2,
    seed = NULL
  )

#Print the summary statistics for each LFQ intensity column, Notice there should be no NA's.
  print(summary(ms_data_imputed))

#The norm_data function will normalize the data to correct for non-biological biases or variations in samples
  norm_data <- normalize_data(ms_data_imputed, method = 'quantile')
  
  cat("Imputation and normalization of data set complete", '\n')

#Identify differentially expressed proteins between experimentally groups using the limma package within the promor function find_dep. P-values from a t-test are adjusted using the Benjamini-Hochberg method to mitigate FDR. Cut off for p-values and adjusted p-values are 0.05.  
  dep_data <- find_dep(
    norm_data,
    save_output = F,
    save_tophits = F,
    file_path = NULL,
    adj_method = "BH",
    cutoff = 0.05,
    lfc = 1,
    n_top = 20
  )

#Creating an annotation output file for the significant proteins by determining top hits based on adjusted p-value and merging the fasta.header column from initial_data within top_hits.
  p_value <- as.data.frame(dep_data$p.value)

  top_hits <- p_value |>
    arrange(groupR) |>
    mutate(number = 1:n()) |>
    mutate(FDR = (groupR*nrow(p_value))/number) |>
    filter(FDR <= 0.05) 

  top_hits <- rownames_to_column(top_hits, var = 'Protein.IDs')

  annotation <- top_hits |>
    left_join(inital_data, by = 'Protein.IDs') |>
    select(Protein.IDs, Fasta.headers)

#Write the annotation data frame to a CSV file in the working directory
  write.csv(annotation, "tophits_file.csv", row.names = FALSE)
  cat("Identification of differentially abundant proteins complete. See tophits_file.csv.", '\n')

#PCA plot to visualize the clustering of the data to see how related each experimental group are to each other. PCA plot is saved in the users current working directory. Using the first two principle components in the plot. This video on youtube helped me create PCA plot:https://www.youtube.com/watch?v=0Jp4gsfOLMs&ab_channel=StatQuestwithJoshStarmer
  
  pca <- prcomp(t(ms_data_imputed), scale=T)
  pca.data <- data.frame(Replicates=rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
  
  pca.var <- pca$sdev^2
  pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
  
  pca_plot <- ggplot(data = pca.data, aes(x=X, y=Y, label=Replicates)) +
    geom_text() +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "" )) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "" )) +
    theme_bw() +
    ggtitle("PCA")
  
  print(pca_plot)
  ggsave("pca.pdf", plot = pca_plot, width = 7, height = 7)
  cat("PCA plot successfully generated and saved in current working directory", '\n')
  
  
#Creating and saving a volcano plot in the users current working directory using the volcano_plot function. A volcano plot is used to visualize the difference in protein expression in each of the experimental groups. Function uses the adjusted p-value data from dep_data with a significance threshold in the adjusted p-values of 0.05. 
  
  vp <- volcano_plot(
    dep_data,
    adj_method = "BH",
    sig = "adjP",
    cutoff = 0.05,
    lfc = 1,
    line_fc = TRUE,
    line_p = TRUE,
    palette = "viridis",
    text_size = 10,
    label_top = T,
    n_top = 10,
    save = T,
    file_path = '.',
    file_name = "Volcano_plot",
    file_type = "pdf",
    plot_height = 7,
    plot_width = 7,
    dpi = 80
  )
  print(vp)
  cat("Volcano plot successfully generated and saved in current working directory", '\n')

}

LFQ_proteomics_pipeline()

