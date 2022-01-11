# 
# 
# bm_MAXENTprepareWorkdir <- function(Data, xy, calibLines = NULL, RunName = NULL,
#                                     VarImport = 0, evalData = NULL, evalxy =  NULL,
#                                     species.name = NULL, modeling.id = '',
#                                     background_data_dir = 'default')
# {
#   cat('\n\t\tCreating Maxent Temp Proj Data...')
#   
#   ## initialise output
#   MWD <- list()
#   class(MWD) <- "maxent_workdir_info"
#   
#   ## default parameters setting
#   if (is.null(RunName)) { RunName <- colnames(Data)[1] }
#   if (is.null(species.name)) { species.name <- colnames(Data)[1] }
#   if (is.null(calibLines)) { calibLines <- rep(TRUE, nrow(Data)) }
#   
#   ## define all paths to files needed by MAXENT.Phillips
#   nameFolder = file.path(species.name, 'models', modeling.id)
#   m_workdir <- file.path(nameFolder, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#   while (file.exists(m_workdir)) { # check wordir unicity
#     m_workdir <- file.path(nameFolder, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#   }
#   m_outdir <- file.path(nameFolder, paste0(RunName, '_MAXENT.Phillips_outputs'))
#   m_predictDir <- file.path(m_workdir, "Predictions")
#   
#   MWD$m_workdir <- m_workdir
#   MWD$m_outdir <- m_outdir
#   MWD$m_outputFile <- file.path(m_outdir, paste0(RunName, '_Pred_swd.csv'))
#   MWD$m_predictDir <- m_predictDir
#   
#   ## directories creation
#   dir.create(m_workdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
#   dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
#   dir.create(m_predictDir, showWarnings = FALSE, recursive = TRUE, mode = '777')
#   
#   
#   ## Presence Data --------------------------------------------------------------------------------
#   presLines <- which((Data[, 1] == 1) & calibLines)
#   absLines <- which((Data[, 1] == 0) & calibLines)
#   Sp_swd <- cbind(rep(RunName, length(presLines))
#                   , xy[presLines, ]
#                   , Data[presLines, 2:ncol(Data), drop = FALSE])
#   colnames(Sp_swd) <- c('species', 'X', 'Y', colnames(Data)[2:ncol(Data)])
#   
#   m_speciesFile <- file.path(m_workdir, "Sp_swd.csv")
#   write.table(Sp_swd, file = m_speciesFile, quote = FALSE, row.names = FALSE, sep = ",")
#   MWD$m_speciesFile <- m_speciesFile
#   
#   
#   ## Background Data (create background file only if needed) --------------------------------------
#   if (background_data_dir == 'default') {
#     # keep only 0 of calib lines
#     Back_swd <- cbind(rep("background", length(absLines))
#                       , xy[absLines, ]
#                       , Data[absLines, 2:ncol(Data), drop = FALSE])
#     colnames(Back_swd) <- c("background", colnames(Back_swd)[-1])
#     
#     m_backgroundFile <- file.path(m_workdir, "Back_swd.csv")
#     write.table(Back_swd, file = m_backgroundFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
#     MWD$m_backgroundFile <- m_backgroundFile
#   } else { ## use background directory given as an option
#     MWD$m_backgroundFile <- background_data_dir
#   }
#   
#   
#   ## Prediction Data ------------------------------------------------------------------------------
#   Pred_swd <- cbind(rep("predict", nrow(xy))
#                     , xy
#                     , Data[, 2:ncol(Data), drop = FALSE])
#   colnames(Pred_swd)  <- c("predict", colnames(xy), colnames(Data)[-1])
#   
#   m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
#   write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
#   MWD$m_predictFile <- m_predictFile
#   
#   
#   ## dealing with independent evaluation data -----------------------------------------------------
#   if (!is.null(evalData)) {
#     Pred_eval_swd <- cbind(rep("predictEval", nrow(evalxy))
#                            , evalxy
#                            , evalData[, 2:ncol(evalData), drop = FALSE])
#     colnames(Pred_eval_swd) <- c("predict", colnames(Back_swd)[-1])
#     
#     m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
#     write.table(Pred_eval_swd, file = m_predictEvalFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
#     MWD$m_predictEvalFile <- m_predictEvalFile
#   }
#   
#   return(MWD)
# }
# 
# 
# 
# bm_MAXENTdeleteWorkdir <- function(MWD, silent = FALSE)
# {
#   if (!silent) { cat('\n\tRemoving Maxent Temp Data..') }
#   if (inherits(MWD, "maxent_workdir_info")) {
#     unlink(unique(sub("/part([[:digit:]]+)$", "", MWD$m_workdir)), recursive = TRUE, force = TRUE)
#   } else if (!silent) {
#     cat('\n\t! Invalid maxent work dir object -> MAXENT.Phillips temp files have not been removed')
#   }
# }
# 
# 
# setGeneric(".Prepare.Maxent.Proj.WorkDir",
#            def = function(Data, ...) {
#              standardGeneric(".Prepare.Maxent.Proj.WorkDir")
#            })
# 
# 
# setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='data.frame'),
#           def = function(Data, xy, species.name =".", proj.name=".", silent=FALSE)
#           {
#             if (!silent) { cat('\n\t\tCreating Maxent Temp Proj Data...') }
#             
#             ## initialise output
#             MWD <- list()
#             class(MWD) <- "maxent_workdir_info"
#             
#             ## define all paths to files needed by MAXENT.Phillips
#             m_workdir <- file.path(species.name, proj.name, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#             while (file.exists(m_workdir)) { # check wordir unicity
#               m_workdir <- file.path(species.name, proj.name, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#             }
#             dir.create(m_workdir, recursive = TRUE, showWarnings = FALSE)
#             MWD$m_workdir <- m_workdir
#             
#             # Proj Data
#             if (is.null(xy)) { xy <- matrix(1, nrow = nrow(Data), ncol = 2, dimnames = list(NULL, c("X", "Y"))) }
#             Proj_swd <- cbind(rep("proj", nrow(xy)), xy, Data)
#             colnames(Proj_swd)  <- c("proj", "X", "Y", colnames(Data))
#             
#             m_predictFile <- file.path(m_workdir, "Pred_swd.csv")
#             write.table(Proj_swd, file = m_predictFile, quote = FALSE,  row.names = FALSE, col.names = TRUE, sep = ",")
#             MWD$m_predictFile <- m_predictFile
#             
#             return(MWD)
#           })
# 
# 
# setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data = 'RasterStack'), 
#           def = function(Data, species.name = ".", proj.name = ".", silent = FALSE, split.proj = 1)
#           {
#             if (!silent) { cat('\n\t\tCreating Maxent Temp Proj Data...') }
#             
#             ## initialise output
#             MWD <- list()
#             class(MWD) <- "maxent_workdir_info"
#             
#             ## define all paths to files needed by MAXENT.Phillips
#             m_workdir <- file.path(species.name, proj.name, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#             while (file.exists(m_workdir)) { # check wordir unicity
#               m_workdir <- file.path(species.name, proj.name, paste0('m_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE)))
#             }
#             
#             ## create the list of extent our raster will be crop at
#             pred.nrow <- nrow(Data)
#             pred.ncol <- ncol(Data)
#             seq.col <- round(seq(1, pred.ncol, length.out = split.proj + 1))
#             ext.list <- lapply(1:split.proj, function(i) { extent(Data, 1, pred.nrow, seq.col[i], seq.col[i + 1]) })
#             
#             # Proj Data
#             m_predictFile <- NULL
#             for (spl in 1:split.proj) {
#               
#               ## create tmp directory
#               m_workdirTmp <- file.path(m_workdir, paste0("part", spl))
#               dir.create(m_workdirTmp, showWarnings = FALSE, recursive = TRUE)
#               MWD$m_workdir[[paste0("part", spl)]] <- m_workdirTmp
#               
#               for (l in names(Data)) {
#                 m_predictFileTmp <- file.path(m_workdirTmp, paste0(l, '.asc'))
#                 
#                 if (!file.exists(m_predictFileTmp)) {
#                   ras = raster::subset(Data, l, drop = TRUE)
#                   
#                   if (!silent) { cat("\n\t\t\t > ", l, "\t:\t") }
#                   if (split.proj == 1) { ## no croping in this case => just write the raster as an ascii file
#                     if (grepl(".asc", filename(ras))) {
#                       if (!silent) { cat("copying ascii file") }
#                       file.copy(filename(ras), m_predictFileTmp)
#                     } else {
#                       if (!silent) { cat("creating ascii file") }
#                       writeRaster(ras, filename = m_predictFileTmp, format = 'ascii', overwrite = TRUE)
#                     }
#                   } else { ## crop the raster within parts
#                     crop(ras, ext.list[[spl]], filename = m_predictFileTmp, format='ascii', overwrite = TRUE)
#                   }
#                 } else if (!silent) { cat("\n", m_predictFileTmp, 'already created !') }
#                 m_predictFile <- c(m_predictFile, m_predictFileTmp)
#               }
#               MWD$m_predictFile[[paste0("part", spl)]] <- m_predictFile
#             }
#             return(MWD)
#           })
# 
# 
# setMethod('.Prepare.Maxent.Proj.WorkDir',
#           signature(Data = 'RasterLayer'),
#           def = function(Data,
#                          species.name = ".",
#                          proj.name = ".",
#                          silent = FALSE) {
#             .Prepare.Maxent.Proj.WorkDir(Data = stack(Data),
#                                          species.name = species.name,
#                                          proj.name = proj.name ,
#                                          silent = silent)
#           })
# 
