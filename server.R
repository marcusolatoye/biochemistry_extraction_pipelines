library(shiny)

shinyServer(function(input, output, session) {
  
            data_extractor <- function(directory_path, experiment, runner, file_names, time_value, replications){
  
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
                                    'Mass List Match: Extractables and Leachables HRAM Compound Database', 
                                    "Area: b6.raw (F2)",	"Area: IS-4.raw (F8)", "Group CV [%]: Control", "Group CV [%]: Sample")
                  
                  
                  # Loop through all rep files
                for(f in 1:length(list_input_filnames)){
                    
                    if (!require(openxlsx)) {
                        install.packages("openxlsx")
                    }
                    
                    library(openxlsx)
                    
                    print(paste("Reading input file:", list_input_filnames[f], sep=" "))
                    
                    #input_file = read.csv(paste(output_dir, list_input_filnames[f], sep=""), header = F, stringsAsFactors = F)
                    input_file <- read.xlsx(paste(output_dir, list_input_filnames[f], sep=""), colNames = FALSE)
                    head(input_file)
                    
                    
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
                    input_file$RT_MW <- paste(round(input_file$`RT [min]`, 1), round(input_file$`Molecular Weight`, 3), sep='_')
                    
                    input_file = input_file[order(input_file$`RT [min]`, input_file$`Molecular Weight`, decreasing = FALSE),]
                    
                    # Filter by the retention time threshold
                    input_file = input_file[which(input_file$`RT [min]` < time_threshold),]
                    
                    # Filter by only chemicals with MS column having DDA for preferred ion
                    input_file = input_file[which(input_file$MS2 == 'DDA for preferred ion'),]
                    
                    input_file$TechRep <- rep(f, nrow(input_file))
                    input_file$EXP = rep(exp_name, nrow(input_file))
                    input_file$File = rep(list_input_filnames[f], nrow(input_file))
                    head(input_file)
                    
                    # Re-order the columns
                    # Bring the EXP, TechRep, and Name (MW_RT) Columns forward
                    input_file = input_file[,c(ncol(input_file),(ncol(input_file)-1),(ncol(input_file)-2),(ncol(input_file)-3),2:(ncol(input_file)-4))]
                    head(input_file)
                    Standard_Chem_Name = '7-Hydroxycoumarine'
                    Standard_MWRT = c('6.3_162.03')  #, '6.3_162.031'
                    
                    stdrd_value_vec = c()
                    
                    non_area_colnames = c("File", "EXP", "TechRep", "RT_MW", "Chemical_Name", "Formula", "Molecular Weight", "RT [min]", "Area (Max.)", "MS2")
                    area_colnames = colnames(input_file)[!colnames(input_file) %in% non_area_colnames]
                    counter = 1
                    
                    # Loop through all environments
                    for(c in area_colnames){
                    
                        print(paste('Trait being processed is: ', c, ' for file: ', list_input_filnames[f], sep=''))
                        stdrd_value = input_file[((input_file$Chemical_Name == Standard_Chem_Name) & (input_file$RT_MW %in% Standard_MWRT)),c]
                        print(paste('Standard value for column ', c, 'is: ', stdrd_value, sep=''))
                        stdrd_value_vec[counter] = stdrd_value
                        
                        counter = counter + 1
                        # Loop through all rows
                        for(r in 1:nrow(input_file)){
                            print(paste('File: ', list_input_filnames[f],'; Column: ', c, '; Row: ', r, sep=""))
                            val = input_file[r, c]
                            est_standard = (val/stdrd_value) * 100
                            input_file[r, c] <- est_standard
                        }
                    }
                    
                    #write.table(input_file, paste(paste(output_dir, gsub('.csv', '', list_input_filnames[f]), sep=''), runner,'reformatted_data', date_time_stamp, 'csv', sep='.'), sep=",", quote=F, row.names=F, col.names = T)
                    
                    # Extract string pattern before -1.raw, -2.raw, and -3.raw.
                    for (col in area_colnames) {
                        index <- which(names(input_file) == col)
                        new_name <- sub("-[1-9]\\.raw.*$", "", col)
                        names(input_file)[index] <- new_name
                    }
                    
                    output_file_null = rbind(output_file_null, input_file)
                    
                    # Reset input_file before going through the loop again
                    input_file=NULL
                
                } # End of files loop
                
                # Write out all data
                
                write.xlsx(output_file_null, paste(paste(output_dir,exp_name, sep=''), runner,'reformatted_data_all', date_time_stamp, 'xlsx', sep='.'), rowNames = FALSE)
                #write.table(output_file_null, paste(paste(output_dir,exp_name, sep=''), runner,'reformatted_data_all', date_time_stamp, 'csv', sep='.'), sep=",", quote=F, row.names=F, col.names = T)
                
                # Identify Compunds with complete three reps of RT_MW
                # Install and load dplyr package in R
                if (!require(dplyr)) {
                    install.packages("dplyr")
                }
                library(dplyr)
                
                find_multiple_reps <- function(data, id_column) {
                
                id_column <- ensym(id_column)
                
                ids_with_multiple_reps <- data %>%
                    group_by(!!id_column) %>%
                    summarise(rep_count = n()) %>%
                    filter(rep_count == as.numeric(as.character(replications))) %>%
                    pull(!!id_column)
                
                data %>%
                    filter(!!id_column %in% ids_with_multiple_reps)
                
                }
                
                # Identify IDs with multiple repetitions in the "RT_MW" column
                output_file_three_reps <- find_multiple_reps(output_file_null, RT_MW)
                
                #Export a CSV file out
                output_file_three_reps = output_file_three_reps[order(output_file_three_reps$RT_MW),]
                write.xlsx(output_file_three_reps, paste(paste(output_dir,exp_name, sep=''), runner,'reformatted_data_three_reps', date_time_stamp, 'xlsx', sep='.'), rowNames = FALSE)
                #write.table(output_file_three_reps, paste(paste(output_dir,exp_name, sep=''), runner,'reformatted_data_three_reps', date_time_stamp, 'csv', sep='.'), sep=",", quote=F, row.names=F, col.names = T)
                  
}

          
          
            observeEvent(input$analyze, {
                  directory_path <- input$directory
                  runner <- input$runner
                  experiment <- input$experiment
                  #file_names <- strsplit(input$file_names, ",")[[1]]
                  files <- input$files
                  time_value <- input$time_value
                  replication <- input$replications
                  
                  #if (length(file_names) > 0) {
                    if (length(files$name) > 0) {
                      #print(files$datapath)
                      #print(files$name)
                      file_names <- files$name #list.files(files$datapath)
                      #file_names <- file.path(files$datapath, files$name)
                      #print(file_paths)
                      data_extractor(directory_path, experiment, runner, file_names, time_value, replications)
                      output$message <- renderText("Data Conversion and Extraction complete. Check your specified directory for the output file.")
                      
                      # Reset input fields to empty
                      updateTextInput(session, "directory", value = "")
                      updateTextInput(session, "runner", value = "")
                      updateTextInput(session, "experiment", value = "")
                      updateTextInput(session, "file_names", value = "")
                      updateTextInput(session, "time_value", value = "")
                      updateTextInput(session, "replications", value = "")
                  } else {
                      output$message <- renderText("Please provide the correct directory path where files are located.")
                  }
            })
          
          observeEvent(input$exit, {
                stopApp()
          })
        })



