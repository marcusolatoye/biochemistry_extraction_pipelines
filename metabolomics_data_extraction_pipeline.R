# Set current working directory
rm(list=ls())
setwd("~/2023/Biochemistry/extraction_app/")

# Store current working directory as a variable
cwd = getwd()

#directory_path = cwd

#time_value = 16

#file_names = c('CY22_051_CW_metabolomic_analysis_MOO_TR1.csv', 'CY22_051_CW_metabolomic_analysis_MOO_TR2.csv', 'CY22_051_CW_metabolomic_analysis_MOO_TR3.csv')

#experiment='CY22_051'

#runner='Aaron'


data_extractor <- function(directory_path, experiment, runner, file_names, time_value){
  
                  output_dir = paste(directory_path,'/', sep='')
                  exp_name = experiment
                  runner = runner
                  list_input_filnames = file_names
                  time_threshold = as.numeric(as.character(time_value))
                  # Set Date and Time
                  current_datetime <- Sys.time()
                  date_time_stamp <- format(current_datetime, "%Y%m%d%H%M%S")
                  
                  #Initiate output file as null because of loop
                  output_file_null=NULL
                  
                  # list of unwanted column names
                  drop_column_names = c('Annotation Source: Predicted Compositions',	'Annotation Source: mzVault Search',
                                        'Annotation Source: Metabolika Search',	'Annotation Source: MassList Search',	'FISh Coverage',
                                        '# mzVault Results',	'# Metabolika Pathways',	'Metabolika Pathways',	'mzVault Best Match',	
                                        'mzVault Library Match: Bamba lab 34 lipid mediators library stepped NCE 10 30 45',	
                                        'mzVault Library Match: Bamba lab 598 polar metabolites stepped NCE 10 30 45',	
                                        'mzVault Library Match: mzVault Autoprocessed May 2019',	
                                        'mzVault Library Match: mzVault Reference May 2019',	
                                        'Mass List Match: Arita Lab 6549 Flavonoid Structure Database',	
                                        'Mass List Match: EFS HRAM Compound Database',	
                                        'Mass List Match: Endogenous Metabolites database 4400 compounds',	
                                        'Mass List Match: Extractables and Leachables HRAM Compound Database')
              
              
                  # Loop through all rep files
                  for(f in 1:length(list_input_filnames)){
                    
                      print(paste("Reading input file:", list_input_filnames[f], sep=" "))
                      input_file = read.csv(paste(output_dir, list_input_filnames[f], sep=""), header = F, stringsAsFactors = F)
                      
                      names_vec = as.vector(as.character(input_file[1,])) # Extract the first row of the dataframe as name vector
                      names_vec
                      colnames(input_file) <- names_vec
                      names(input_file)
                      
                      input_file = input_file[-1,] # Remove the first row which is now the header.
                      colnames(input_file)
                      
                      # Drop unwanted columns
                      input_file = input_file[, !(names(input_file) %in% drop_column_names)]
                      
                      # Change Name column to Chemical Name
                      names(input_file)[names(input_file) == 'Name'] <- 'Chemical_Name'
                      
                      columns_notto_convert = c('Checked', 'Chemical_Name', 'Formula', 'MS2')
                      columns_to_convert = names(input_file)[! names(input_file) %in% columns_notto_convert] 
                      input_file[columns_to_convert] <- lapply(input_file[columns_to_convert], as.numeric)
                      # Create a concatenated column with MW and RT
                      input_file$`Name (MW_RT)` <- paste(round(input_file$`Molecular Weight`, 1), 
                                                         round(input_file$`RT [min]`, 2), sep='_')
                    
                      input_file = input_file[order(input_file$`Molecular Weight`, input_file$`RT [min]`, decreasing = FALSE),]
                      input_file = input_file[which(input_file$`RT [min]` < time_threshold),]
                      input_file$TechRep <- rep(f, nrow(input_file))
                      input_file$EXP = rep(exp_name, nrow(input_file))
                      input_file$File = rep(list_input_filnames[f], nrow(input_file))
                      head(input_file)
                      
                      # Re-order the columns
                      # Bring the EXP, TechRep, and Name (MW_RT) Columns forward
                      input_file = input_file[,c(ncol(input_file),(ncol(input_file)-1),(ncol(input_file)-2),(ncol(input_file)-3),2:(ncol(input_file)-4))]
                      head(input_file)
                      Standard_Chem_Name = '7-Hydroxycoumarine'
                      Standard_MWRT = c('162_6.3', '162_6.31', '162_6.32', '162_6.33')  #162.03046  6.327   #162.03045 6.322
           
                      stdrd_value_vec = c()
                      
                      non_area_colnames = c("File", "EXP", "TechRep", "Name (MW_RT)", "Chemical_Name", "Formula", "Molecular Weight", "RT [min]", "Area (Max.)", "MS2", "Area: b6.raw (F2)",	"Area: IS-4.raw (F8)", "Group CV [%]: Control", "Group CV [%]: Sample")
                      area_colnames = colnames(input_file)[!colnames(input_file) %in% non_area_colnames]
                      counter = 1
                      
                      for(c in area_colnames){
                        
                          print(paste('Trait being processed is: ', c, ' for file: ', list_input_filnames[f], sep=''))
                          stdrd_value = input_file[((input_file$Chemical_Name == Standard_Chem_Name) & (input_file$`Name (MW_RT)` %in% Standard_MWRT)),c]
                          print(paste('Standard value for column ', c, 'is: ', stdrd_value, sep=''))
                          stdrd_value_vec[counter] = stdrd_value
                          
                          counter = counter + 1
                          
                          for(r in 1:nrow(input_file)){
                              print(paste('File: ', list_input_filnames[f],'; Column: ', c, '; Row: ', r, sep=""))
                              val = input_file[r, c]
                              est_standard = (val/stdrd_value) * 100
                              input_file[r, c] <- est_standard
                          }
                      }
                      
                      write.table(input_file, paste(paste(output_dir, gsub('.csv', '', list_input_filnames[f]), sep=''), runner,'reformatted_data', date_time_stamp, 'csv', sep='.'), sep=",", quote=F, row.names=F, col.names = T)
                      
                        # Extract string pattern before -1.raw, -2.raw, and -3.raw.
                          for (col in area_colnames) {
                              index <- which(names(input_file) == col)
                              new_name <- sub("-[1-9]\\.raw.*$", "", col)
                              names(input_file)[index] <- new_name
                          }
                      
                      output_file_null = rbind(output_file_null, input_file)
                      
                      # Reset input_file before going through the loop again
                      input_file=NULL
                  }
                  
                  #Export a CSV file out
                  output_file_null = output_file_null[order(output_file_null$`Molecular Weight`, output_file_null$`RT [min]`, output_file_null$TechRep),]
                  write.table(output_file_null, paste(paste(output_dir,exp_name, sep=''), runner,'reformatted_data', date_time_stamp, 'csv', sep='.'), sep=",", quote=F, row.names=F, col.names = T)
              
              }



data_extractor(directory_path = cwd, experiment = 'CY22_051', runner = 'Aaron', file_names = c('CY22_051_CW_metabolomic_analysis_MOO_TR1.csv', 
                                                                                               'CY22_051_CW_metabolomic_analysis_MOO_TR2.csv', 
                                                                                               'CY22_051_CW_metabolomic_analysis_MOO_TR3.csv'), time_value = 16)







