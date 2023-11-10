getPlateInfo <- function(MBKfilePath)
{
    #' Extract file annotation from MBK file exported from MetaXpress
    #'
    #' Function to extract ImageXpress image file annotations from files exported  
    #'  MetaXpress server software (Molecular Devices)
    #' @param MBKfilePath character or path to file
    #' 
    #' @return data.frame of file annotations.

    message(cat("\nProcessing:",MBKfilePath,"\n"))
    d <- tryCatch(
        {
            parseMBK(paste0(MBKfilePath,'/PlateInfo.MBK'))
        },
        error=function(cond)
        {
            message(cat('error',cond[[1]],'in:',MBKfilePath,"\n"))
            return(NA)
        }
    )
    pid <- gsub('Plate','',basename(MBKfilePath))
    d$well <- fixWellName(paste0(LETTERS[d$WELL_Y],d$WELL_X))
    
    # if z-stack, capture z positions
    if(any(d$Z_INDEX > 0))
    {
        ZSTACK <- TRUE
        z_pos <- as.integer(gsub('ZStep_','',d$Z_INDEX))
    }

    # info common to all images in plate
    plate.info <- do.call(rbind, lapply(seq(nrow(d)), function(i) unlist(strsplit(d$DIRECTORY[i],'|',fixed=TRUE))))
    plate.info <- as.data.frame(plate.info[,!apply(plate.info,2,function(x) all(x==""))])

    # expecting 6 columns of data
    pinames <- c('expt_class','expt_id','plate_name','date','plate_id','timepoint')
    
    colnames(plate.info) <- pinames

    temp <- d$T_POSITION
    d <- d[,c('OBJ_SERVER_NAME','well','SOURCE_DESCRIPTION')]
    colnames(d) <- c('file_name','well','channel')
    if(ZSTACK) d$z_pos <- z_pos

    d$time <- strptime(temp, format='%Y-%m-%d %H:%M:%S')
    
    return(cbind(d,plate.info))
}

findPlateDir <- function(dirpath)
{
    #' Recursively search for directories containing "Plate" in their names
    #'
    #' @param dirpath character or path to top directory to search
    #' 
    #' @return character of directory paths.

    # find directories in top directory
    dl <- tryCatch({list.dirs(dirpath,recursive=FALSE)},error=NA)
    # do not include Segmentation directory if it alredy exists
    dl <- dl[!grepl('Segmentation',dl)]
    if(length(dl) != 1 && is.na(dl[1]))
    {
        message(paste('Could not find any directories in',dirpath))
        return(NA)
    }
    
    if(!any(grepl('Plate',dl))) 
    {
        # check in subdirectories
        subdir <- sapply(dl, function(x) list.dirs(x,recursive=FALSE))
        dl <- append(dl,subdir)
    }
    
    # must have 'Plate' in directory name
    return(dl[grepl('Plate',dl)])
}


makeFileInfo <- function(topdir, save=TRUE, overwrite=FALSE) {
    plate_dir <- findPlateDir(topdir)
    fi <- do.call(rbind, lapply(plate_dir, getPlateInfo))
    fi <- fi[!duplicated(fi$file_name),]    # remove any image file duplicates
    
    fi_filepath <- file.path(topdir,"ImageFileInfo.csv")
    if(save & (!file.exists(fi_filepath) | overwrite)) {
        message(cat("Writing file:",fi_filepath,"\n"))
        write.csv(fi, file=fi_filepath, row.names=FALSE)
    } else {
        message(cat(fi_filepath,"file exists and will not be overwritten\n"))
    }
    return(fi)
}