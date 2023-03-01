#' Upload and combine the 16S bacterial ASV reads with the metadata of the salmon farms.
#'
#' @param sources A vector with at least one element. Contains the names of the data sources that should be kept in the final data set.
#' @param taxalevel A string. Specifies the taxonomic level used for analysis.
#' @param rarefaction_rate A numerical value or NA. Specifies the rarefaction rate applied to the ASV reads.
#'
#' @return A data frame.
data_upload=function(sources, taxalevel, rarefaction_rate, anonymous = TRUE, samples = FALSE){
  # Upload the 16S data
  S16Data = upload_16S()
  
  # Upload the metadata
  MDataF = finalise_metadata(sources)
  
  # Reformat the data set; sample names as rownames, taxa till taxalevel as columns
  MDS1 = split_to_taxa(S16Data)
  MDS2 = taxa_collate(MDS1,taxalevel)
  
  # Merge 16S training data with metadata, take mean over replicates, and clean
  MDS3 = merge_16S_metadata(MDS2, MDataF)
  MDS4 = mean_of_replicates(MDS3, rarefaction_rate)
  MDS5 = data_cleaning(MDS4, rarefaction_rate, anonymous, samples)
  
  return(MDS5)
}

#Function: upload_16S
#Takes:MetaD (see above), various data-transformation parameters, raw data file location, minReads/minSamples
#Returns: MDS (via dataframe called 'Merged4')
#Calls:NULL
#Job:Read and clean-up raw data (NOT merging raw data, now), apply data transformations (inc minReads/samples),
#    merge with metadata (== WRONG, no clean-up, merges training sets, no transformations present, no metadata is read here)
#    This function does the same thing as upload_metadata but then on the training and testing sets of 16S.
upload_16S = function(TestOrTrain){
  message("Function: 'upload_16S' start")
  
  # Upload the data and set NA to 0
  S16Data = as.data.frame(read.csv("./data16s/ASV_L7_transposed.csv"))
  S16Data[is.na(S16Data)] = 0
  
  message("Function: 'upload_16S' end")
  return(S16Data)
}

#Hello
#Function: finalise_metadata
#Takes: NA
#Returns: dataframe
#Calls:upload_metadata
#Job: Metadata uploaded and formatted to correct columns and variable names required for future RF analysis.
finalise_metadata = function(sources){
  message("Function: 'finalise_metadata' start")
  
  # Upload all the metadata and keep / rename required columns
  S16MetaData = upload_metadata()
  S16MetaData = S16MetaData[, colnames(S16MetaData) %in% c("SampleID", "Long_Site_Name",
                                                           "Dataset_Ref", "IQI", "Macro_Distance_from_cage",
                                                           "Grab_sample_month", "Grab_sample_year")] %>% mutate_all(na_if, "") %>%
    dplyr::rename(SampleID = SampleID,
                  Site = Long_Site_Name,
                  Source = Dataset_Ref,
                  IQI = IQI,
                  Distance = Macro_Distance_from_cage,
                  Month = Grab_sample_month,
                  Year = Grab_sample_year)
  
  S16MetaData$Month = str_to_title(S16MetaData$Month) # Capitalize the initial letter of the values in the "Month" column
  S16MetaData$SampleID = str_replace(S16MetaData$SampleID, "-paired-woN-withQual", "") # tidy SampleID
  
  # Standardise the Month names
  Mnth = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  LongMonth = list(c("January"),c("February"),c("March"),c("April"),c("May"),
                   c("June"), c("July"),c("August"),c("September","Sept"),
                   c("October"),c("November"), c("December"))
  for (i in 1: length(Mnth)){ #replace long form Month with short form
    S16MetaData$Month[S16MetaData$Month%in%LongMonth[[i]]]=Mnth[i]
  }
  S16MetaData$Month[!S16MetaData$Month%in%Mnth]=NA # Anything not in Mnth=NA
  
  ## Rename column content ##
  # Change all _ to spaces in Site names and generalize #
  S16MetaData$Site = gsub("_", " ", S16MetaData$Site) # Replace level indicator with space
  
  # Remove unneeded sites
  S16MetaData = S16MetaData[!grepl("Greshornish|blank|unknown|na|0", S16MetaData$Site), ]
  
  # Rename the S16MetaData
  S16MetaData$Site = gsub("^Caolas(.*)", "Caolas", S16MetaData$Site)
  S16MetaData$Site = gsub("^North Shore(.*)", "North Shore", S16MetaData$Site)
  S16MetaData$Site = gsub("^Carradale(.*)", "Carradale", S16MetaData$Site)
  S16MetaData$Site = gsub("^Tanera(.*)", "Tanera", S16MetaData$Site)
  S16MetaData$Site = gsub("^Fishnish(.*)", "Fishnish", S16MetaData$Site)
  
  S16MetaData$Source = gsub("_(.*)", "", S16MetaData$Source) # Remove _Ref from source name
  S16MetaData$Source = str_replace(S16MetaData$Source, "SV", "Shraveena") # Write SV as full
  S16MetaData$Source[which(startsWith(S16MetaData$SampleID, "MOWI"))] = "MOWI" # Set MOWI as source
  S16MetaData = S16MetaData[which(S16MetaData$Source %in% sources), ] # Only keep pre-specified sources
  
  # Remove all letters from Distance column #
  S16MetaData$Distance = as.numeric(gsub("[^0-9.-]", NA, S16MetaData$Distance))
  
  # Replace all na in IQI values to NA values
  S16MetaData$IQI = gsub("na", NA, S16MetaData$IQI)
  
  message("Function: 'finalise_metadata' end")
  return(S16MetaData)
}

#Function: upload_metadata (this is called by finalise_metadata)
#Takes: NA
#Returns: Merged metadataset
#Calls: NA
#Called by: finalise_metadata
#Job: Read and merge metadata files (current version doesn't clean-up metadata,
#     see Readme file in Metadata directory)
upload_metadata=function(){
  message("Function: 'upload_metadata' start")
  
  # Upload the data and set NA to 0
  S16Metadata = as.data.frame(read.csv("./metadata/ASV_L7_metadata.csv"))
  S16Metadata[is.na(S16Metadata)] = 0
  names(S16Metadata)[1]="SampleID"
  
  message("Function: 'upload_metadata' end")
  return(S16Metadata)
}

#Function: split_to_taxa
#Takes:
#Returns: MDS1
#Calls:NULL
#Job:Split (marker dependent) taxonomy string into separate components.
#    Clean-up tax string and remove oddities (e.g. special characters)
split_to_taxa=function(S16Data){
  message("Function: 'split_to_taxa' start")
  message("Splits data into taxa levels")
  
  Run2.1=S16Data
  Run2.2=as.data.frame(t(Run2.1[,-1]))# Transpose data frame (without SampleID)
  colnames(Run2.2) = Run2.1$index # Set SampleID as column names
  Run2.2$Taxa=rownames(Run2.2) # Add column with taxa names
  Run2.3=as.data.frame(Run2.2[,"Taxa"])# Isolate the taxa names, ahead of cleaning them up.
  colnames(Run2.3)[1]="Taxa"
  
  # Split the taxa name into seperate columns
  Run2.4 = tidyr::separate(Run2.3, col =  "Taxa", sep = ".__",
                           into = c("Empty", "Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                           fill = "right") %>% mutate_all(na_if,"")
  Run2.4$Empty = NULL
  
  Run2.4=lapply(Run2.4, gsub, pattern="._", replacement="") # Replace level indicator with space
  Run2.4=lapply(Run2.4, gsub, pattern="[.]", replacement="") # Replace set of any characters with space
  Run2.4=lapply(Run2.4, gsub, pattern="Ambiguous_taxa", replacement="") # Remove ambiguous taxa
  
  # Remove categories that are not part of the taxon
  Run2.4=lapply(Run2.4, gsub, pattern="uncult",      replacement=NA,ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern="unknown",     replacement=NA,ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern=" incertae sedis",    replacement='',ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern="benzene",     replacement=NA,ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern="unidentified",replacement=NA,ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern=";_",          replacement=NA,ignore.case=TRUE)
  Run2.4=lapply(Run2.4, gsub, pattern="[^[:alnum:]]",replacement="",ignore.case=TRUE)
  Run2.4=lapply(Run2.4, function(x) str_replace(x,";",""))
  Run2.4=lapply(Run2.4, trimws, which="both") # removes leading and tailing white_spaces
  
  Run2.5=as.data.frame(Run2.4)
  Run2.5[is.na(Run2.5)]="XUnidentified"
  
  # Add taxa columns to both data sets for merging
  Run2.5$Taxa = rownames(Run2.2)
  Run2.2$Taxa=rownames(Run2.2)
  Run2X=merge(Run2.5,Run2.2)
  
  A="_" # Unique seperator
  Run2X=within(Run2X,{
    Phylum=paste(Kingdom,Phylum,sep=A)
    Class=paste(Phylum,Class,sep=A)
    Order=paste(Class,Order,sep=A)
    Family=paste(Order,Family,sep=A)
    Genus=paste(Family,Genus,sep=A)
    Species=paste(Genus,Species,sep=A)
  })#rebuild full taxonomy
  
  message("Function: 'split_to_taxa' end")
  return(Run2X) # Returns MDS1
}

#Function: taxa_collate
#Takes:MDS1 (from MergeRuns), taxalevel)
#Returns: MDS2
#Calls:NULL
#Job:Collates reads across given taxonomic level
taxa_collate=function(MDS1,taxalevel){
  message("Function: 'taxa_collate' start")
  message("Collates data across specified taxonomic level")
  
  # Select all columns that are not == taxalevel
  df=as.data.table(MDS1)
  Excludes=c("Taxa","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  Excludes2=setdiff(Excludes,taxalevel)
  
  # Remove all taxa columns that are not == taxalevel
  df2=df[,!..Excludes2] # df2 contains only the Taxon and the sample data.
  names(df2)[1]="Taxon"
  
  # Convert data to long data format
  df3=reshape2::melt(df2,id="Taxon")
  df4=data.table(df3)
  names(df4)=c("Taxon","SampleID","Total")# Give sensible names for summarisation step via data.table
  df4$Taxon=as.factor(df4$Taxon)
  
  df5=df4[,list(#this generates the summary values by taxonomic group.
    Total=sum(Total)),
    by=c("Taxon","SampleID")]
  
  message("Function: 'taxa_collate' end")
  return(df5) # Returns MDS2
}

#Function: merge_16S_metadata
#Takes: MDS2 (Reads collated across taxon), Metadata, rarefaction_rate threshold
#Returns: MDS3
#Calls:Nil
#Job:Merge MDS2 (16S data) back with collated metadata.
merge_16S_metadata=function(MDS2, MDataF, rarefaction_rate){
  message("Function: 'merge_16S_metadata' start")
  message("This merges MDS2 with MetaData")
  
  # Convert MDS2 into a wide format
  MDS2a = dcast(MDS2,SampleID~Taxon,value.var = "Total")
  MDS2a = as.data.frame(MDS2a)
  
  # Convert sampleID to characters for merging
  colnames(MDS2a)[1] = "SampleID"
  MDS2a$SampleID=as.character(MDS2a$SampleID)
  MDataF$SampleID=as.character(MDataF$SampleID)
  
  # Merge the sequence data and metadata
  MDS3=merge(MDS2a, MDataF, by="SampleID", all = FALSE) # Only keep samples present in both data frames
  
  message("Function 'merge_16S_metadata' end")
  return(MDS3)
}

#Function: mean_of_replicates
#Takes:MDS3 (the complete data frame, including taxa and metadata)
#Returns: mean values of taxa and metadata across grab replicates
#Calls:Nil
#Job: Thorsten data includes grab replicates (so NGS generated from
#grabs that are combined to generate the IQI value). This function generates
#mean taxa and metadata and returns the revised dataframe. Remove
#all data columns where no reads are present for the taxa.
" In Thorsten data, there are cross replicates A and B, for the same Site and Station,
  the function replicates the means across both values of A and B and removes the duplicates. "
mean_of_replicates=function(MDS3, rarefaction_rate){
  message("Function: 'mean_of_replicates' start")
  # Remove all rows where the IQI is NA
  if (sum(is.na(MDS3$IQI)) != 0){
    MDS3 = MDS3[-which(is.na(MDS3$IQI)), ]
  }
  
  # Set the number of segments that should be removed per source for an unique code
  duplicate = data.frame("Source" = c("Thorsten", "MeioMetBar", "MOWI", "Shraveena", "SSF", "VAW"),
                         "Remove_Segments" = c(6, 1, 0, 2, 0, 1))
  
  # Change all - in SampleID to _
  MDS3$SampleID = as.character(lapply(MDS3$SampleID, function(x) {gsub("-", "_", x)}))
  
  # Take mean over duplicate rows when duplicate = data.frame()
  if (is.data.frame(duplicate)){
    MDS3$Unique = MDS3$SampleID # Create the unique column
    for (i in 1:nrow(duplicate)){
      # Find which rows belong to the data source
      rows = which(MDS3$Source == duplicate[i, 1])
      for (j in rows){
        # Split all SampleIDs, remove segments of REF code, add to unique column
        SampleID_split = unlist(strsplit(MDS3$SampleID[j], "_"))
        segment = ifelse(length(SampleID_split) - duplicate[i, 2] <= 0, 1, length(SampleID_split) - duplicate[i, 2])
        MDS3$Unique[j] = paste(SampleID_split[1:segment], collapse = "_")
      }
    }
    # Add IQI value to all unique site names
    MDS3$Unique = paste(MDS3$Unique, MDS3$Station, MDS3$IQI, sep = "_")
    
    # Keep necessary character columns (for grouping)
    keep_col = append(colnames(select_if(MDS3, is.numeric)), c("Site", "Source", "IQI", "Month", "Year", "Unique", "Distance"))
    MDS3 = MDS3[, keep_col]
    
    # Combine the duplicated rows in data based on character columns
    MDS3 = MDS3 %>%
      group_by_at(colnames(select_if(MDS3, is.character))) %>%
      summarise_all(mean)
    
    # Round bacteria taxa to integers
    MDS3[, grepl("Bacteria", names(MDS3))] = floor(MDS3[, grepl("Bacteria", names(MDS3))])
    
    # Ungroup columns and rename Unique to SampleID
    MDS3 = ungroup(MDS3)
    colnames(MDS3)[colnames(MDS3) == "Unique"] = "SampleID"
  }
  
  message("Function: 'mean_of_replicates' end")
  return(MDS3)
}

#Function: data_cleaning
#Takes: MDS4
#Returns: cleaned data set
#Calls: NA
#Job: Remove unneeded sites, rarefaction_rate, remove columns with zeroes,
# rename the sources to full names, and anonymise.
data_cleaning = function(MDS4, rarefaction_rate, anonymous, samples){
  message("Function: 'data_cleaning' start")
  
  # REMOVE THIS #
  # Only to ensure no, 1170, 5244, 13493 use the same samples -> test rarefy difference
  if(length(samples) > 1){MDS4 = MDS4[which(MDS4$SampleID %in% samples), ]}
  
  ## Rarefy the data set
  if (!is.na(rarefaction_rate)){
    # Select and rarefy bacteria columns
    data_taxa = MDS4[, grepl("Bacteria", names(MDS4))]
    data_taxa = rrarefy(data_taxa, rarefaction_rate)
    
    # Bind metadata columns to rarefied data
    data_other =  MDS4[, !grepl("Bacteria", names(MDS4))]
    MDS4 = cbind(data_taxa, data_other)
    MDS4[is.na(MDS4)] = 0 # Where merge results in NA for some taxa for some sites, replace with zero.
  }
  
  ## Remove rows with read count < rarefaction_rate
  # Find rows with fewer samples than rarefaction_rate
  data_taxa = MDS4[, grepl("Bacteria", names(MDS4))]
  taxa_reads = rowSums(data_taxa)
  if(is.na(rarefaction_rate)){rarefaction_rate = 100} # Remove all samples below 100 reads
  EE = which(taxa_reads < rarefaction_rate)
  Removes = MDS4$SampleID[EE]
  
  # Send a message what rows are removed
  message("The following samples have been removed (<100 or <Rarefaction threshold):")
  X=sapply(Removes[1:length(EE)],function(N) message(paste(N, sep="; "))); rm(X)
  
  # Remove rows with fewer samples than rarefaction_rate
  if (length(EE) > 0){MDS4 = MDS4[-EE, ]}
  
  ## Remove rows where IQI value is zero
  if (0 %in% MDS4$IQI == TRUE){
    MDS4 = MDS4[-which(MDS4$IQI == 0), ]
  }
  MDS4$IQI = as.numeric(MDS4$IQI)
  
  ## Remove all taxa columns with only zeroes
  data_taxa = MDS4[, grepl("Bacteria", names(MDS4))]
  data_other =  MDS4[, !grepl("Bacteria", names(MDS4))]
  new_columns = data_taxa %>% select_if(colSums(.) != 0)
  MDS4 = data.frame(new_columns, data_other)
  
  # Anonymise the data
  if (anonymous == TRUE){
    MDS4 = MDS4 %>% 
      dplyr::mutate(Site = paste("Site", ifelse(as.numeric(factor(MDS4$Site)) < 10, 
                                                paste0("0", as.numeric(factor(MDS4$Site))),
                                                as.numeric(factor(MDS4$Site)))),
                    Source = str_replace_all(MDS4$Source, c("MeioMetBar" = "Source 01", 
                                                            "MOWI" = "Source 02",
                                                            "Shraveena" = "Source 03",
                                                            "Thorsten" = "Source 04", 
                                                            "SSF" = "Source 05",
                                                            "VAW" = "Source 06"))) 
  }
  
  message("Function: 'data_cleaning' end")
  return(MDS4)
}

