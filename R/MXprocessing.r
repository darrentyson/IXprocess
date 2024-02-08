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
            parseMBK(file.path(MBKfilePath,'PlateInfo.MBK'))
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
    # do not include Segmentation directory if it already exists
    dl <- dl[!grepl('Segmentation',dl)]
    # do not include old/older directory if it already exists
    dl <- dl[!grepl('[oO]ld',dl)]
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
    expt_name <- basename(topdir)
    plate_dir <- findPlateDir(topdir)
    fi <- do.call(rbind, lapply(plate_dir, getPlateInfo))
    fi <- fi[!duplicated(fi$file_name),]    # remove any image file duplicates

    fi_filepath <- file.path(topdir,paste0(expt_name,"ImageFileInfo.csv"))
    if(save & (!file.exists(fi_filepath) | overwrite)) {
        message(cat("Writing file:",fi_filepath,"\n"))
        write.csv(fi, file=fi_filepath, row.names=FALSE)
    } else {
        message(cat(fi_filepath,"file exists and will not be overwritten\n"))
    }
    return(fi)
}

hasEmptyPlateInfo <- function(dirvec, verbose=FALSE) {
    sapply(dirvec, function(topdir) {
        plate_dir <- findPlateDir(topdir)
        if(length(plate_dir)==0) {
            message(cat("No Plate directories found in ",topdir,"\n"))
            return(FALSE)
        }
        
        file_paths <- sapply(plate_dir, function(x) file.path(x,'PlateInfo.MBK'))
        if(verbose) message(cat(paste("Expecting PlateInfo.MBK files: ",file_paths, collapse="\n")))
        file_check <- tryCatch({file.exists(file_paths)},error=function(cond) {FALSE})
        if(any(!file.exists(file_paths))) {
            message(cat("Missing PlateInfo.MBK files in:",paste(file_paths[!file.exists(file_paths)], collapse="\n")))
            file_paths <- file_paths[file.exists(file_paths)]
        }
        
        if(length(file_paths)==0) {
            message(cat("Could not find PlateInfo.MBK in ",topdir))
            return(TRUE)
        } else {
            size <- file.info(file_paths)$size
            return(any(size==0))
        }
    })
}

rename_imdir <- function(export_dirpath, pattern=export_pattern) {
    confirmed_export_dirs <- export_dirpath[grepl(pattern,export_dirpath)]
    confirmed_export_dirs <- confirmed_export_dirs[grepl(pattern,basename(confirmed_export_dirs))]
    
    file.rename(confirmed_export_dirs,file.path(dirname(confirmed_export_dirs),"images"))
}

find_and_rename_image_dirs <- function(path, depth=3) {
    lvl_1 <- list.dirs(path, recursive=FALSE)
    lvl_2 <- lapply(lvl_1, function(d) list.dirs(d, recursive=FALSE))
    
    if(any(unlist(sapply(lvl_2, function(x) grepl(export_pattern, x))))) {
        lvl2_exports <- unlist(lvl_2)[unlist(sapply(lvl_2, function(x) grepl(export_pattern, x)))]
        rename_imdir(lvl2_exports)
        lvl_2 <- lapply(lvl_1, function(d) list.dirs(d, recursive=FALSE))
    }
    
    lvl_3 <- lapply(lvl_2, function(d) list.dirs(d, recursive=FALSE))
    
    if(any(unlist(sapply(lvl_3, function(x) grepl(export_pattern, x))))) {
        lvl3_exports <- unlist(lvl_3)[unlist(sapply(lvl_3, function(x) grepl(export_pattern, x)))]
        rename_imdir(lvl3_exports)
        lvl_3 <- lapply(lvl_2, function(d) list.dirs(d, recursive=FALSE))
    }
    
    out <- c(unlist(lvl_2)[unlist(sapply(lvl_2, function(x) grepl("images",x)))],
             unlist(lvl_3)[unlist(sapply(lvl_3, function(x) grepl("images",x)))])
    out <- out[basename(out)=="images"]
    out
}

has_segdir <- function(imdir_paths) grepl("[sS]egmentation",list.dirs(dirname(imdir_paths), recursive=FALSE))

has_fileinfo <- function(paths) {
    # FileInfo.csv files should not be in images directory but their parent directory
    paths[basename(paths)=="images"] <- dirname(paths[basename(paths)=="images"])
    out <- sapply(paths, function(p) any(grepl("FileInfo\\.csv",list.files(p, recursive=FALSE))))
    names(out) <- paths
    return(out)
}

export_pattern <- "[[:alnum:]]{8}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{4}-[[:alnum:]]{12}"

missing_imfiles <- function(finfo, topdir=NULL) {
    if(class(finfo) == "character" & tryCatch({file.exists(finfo)}, error=function(cond) {FALSE})) {
        if(is.null(topdir)) topdir <- dirname(finfo)
        finfo <- read.csv(finfo)
    }
    
    if (class(finfo) == "data.frame" & any(grepl("plate_name", colnames(finfo)))) {
        plate_dir <- findPlateDir(topdir)
        imfiles <- unlist(lapply(plate_dir, function(plate) list.files(plate, full.names=TRUE)))
        finfo_imfiles <- file.path(plate_dir,finfo$file_name)
        missing_files <- setdiff(finfo_imfiles, imfiles)
        return(missing_files)
    }  else {
        message("Could not find fileInfo to determine location of imfiles")
        return(NA)
    }
}