getUserInput <- function() {

  # What is the Area of Interest
  Strategy_landscapes <- vect(paste0(data_directory,'/data/global_data/strategy_landscapes/AllPriorityLandscapes.gdb'))

  AOI_vec <- c()
  repeat{
    # What is the Area of Interest
    cat(paste0('\n\nWhat is the Area(s) of Interest: Please list the the Wildfire Crisis Strategy Landscape(s) (as listed below) or the two letter state code of the State(s).
               Press enter to add another area or write "done" \n \n', paste(Strategy_landscapes$NAME, collapse = ', ')))

    aoi <-  scan(what = character(), sep = '\n',n = 1, quiet = TRUE)

    # check if user wants to finish
    if(tolower(aoi)== 'done'){
      break
    }
    
    if(!aoi %in% Strategy_landscapes$NAME){
      cat(paste0('"', aoi, '" is not part of the  Wildfire Crisis Strategy Landscape(s) or two letter state code. (Please check your spelling)'))
    }else{
      AOI_vec <- c(AOI_vec, aoi)
    }

  }

  # Prompt what date range
  cat('What is the desired date range for fire occurence: put the minimum year as "1:" and the maximum year as "2:"')
  YearRange <- scan(what = numeric(), n = 2, quiet = TRUE)

  # Prompt for Previously Selected Fires or all Historic Data
  cat("Do you have previously selected fire perimeters?
  \n If YES: The follow will prompt the input MTBS fire ID's
  \n IF NO: The following will select all historic fires with available treatment data")
  AllFiresYorN <- scan(what = character(), n = 1, quiet = TRUE)

  # If YES: have users type in Selected Fire IDs,
  # IF NO: End user input
  if(AllFiresYorN == "YES"){

    # Prompt User to upload csv in the global_data folder
    cat("Please Upload a CSV with the names of the selected fires. The CSV must have three columns: 'MTBS_ID' as MTBS ID numer, 'Name' with all selected fire names and 'Ig_Date' as the select fires' ignition date in the format MM/DD/YYYY.
      \n Type the name of the csv when finished (ex. fires.list.csv)")

    FireListfp <- scan(what = character(), n = 1, quiet = TRUE)

    # Double Check that user uploaded correct file and in the correct location

    while (!file.exists(paste0(data_directory,'/data/global_data/treatments/fire_list/', FireListfp))) {

      # Prompt User to upload csv in the global_data folder
      cat(paste0(data_directory,'/data/global_data/', FireListfp," does not Exist. Please check your spelling and placement.
      \n Please try to upload a CSV with the names of the selected fires again. The CSV must have three columns: 'MTBS_ID', 'Name' with all selected fire names and 'Ig_Date' as the select fires' ignition date in the format MM/DD/YYYY.
      \n Type the name of the csv when finished (ex. fire.list.csv)"))

      FireListfp <- scan(what = character() , n = 1, quiet = TRUE)

    }

    # User Update
    cat(paste0('Reading in ', FireListfp))

    # Read in selected Fires
    FireList <- read.csv(paste0(data_directory,'/data/global_data/treatments/fire_list/', FireListfp))

    # Store finalize file path
    FireIDvect <- FireList$MTBS_ID


  }else{

    #Return a character that says "All Historic Fires"
    FireIDvect  <- "All Historic Fires"
  }


  return(list(AOI = AOI_vec, Year_range = YearRange, Fires = FireIDvect))


}
