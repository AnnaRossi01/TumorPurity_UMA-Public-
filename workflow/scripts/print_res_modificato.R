printResults_all <- function(res,
                             ismedian=FALSE,
                             outdir="./",
                             is_pdf=F,
                             projectname="project",
                             svinput=NULL,
                             lSVinput=NULL,
                             rainbowChr=TRUE)
{
  GAMMA <- 1
  if(any(names(res)=="gamma")) GAMMA <- res$gamma
  
  # Outdir cleaning
  outdir <- gsub("/$","",outdir)  # removes any / after outdir that may be dirty
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # --- LOOPS ON SAMPLES ---
  for(i in 1:length(res$allTracks.processed))
  {
    samplename <- names(res$allTracks)[i]
    # One folder for each sample
    sample_dir <- file.path(outdir, samplename) # file.path can automatically set the correct separator for the OS
    if(!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
    
    # 1. Plotting (PNG or PDF) in the sample folder
    ext <- if(is_pdf) ".pdf" else ".png"
    plot_file <- file.path(sample_dir, paste0(samplename, "_profile", ext))
    
    if(!is_pdf) {
      png(plot_file, width = 5500, height = 2496, res=300)
    } else {
      pdf(plot_file, width=15, height=5.5)
    }
    
    try({
      suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                    purity=res$allSolutions[[i]]$purity,
                                    ploidy=res$allSolutions[[i]]$ploidy,
                                    gamma=GAMMA,
                                    ismedian=ismedian,
                                    allchr=res$chr,
                                    ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                    sol=res$allSolutions[[i]],
                                    svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                    rainbowChr=rainbowChr,
                                    is_pdf=is_pdf))
      title(samplename)
    })
    dev.off()
    
    # 2. Writing Profile (Summary) in the sample folder
    try({
      writeProfile(prof=res$allProfiles[[i]],
                   samplename=paste0(samplename, "_", projectname),
                   outdir=sample_dir)
    })
  }
  
  # --- BLOCCO REFITTED (SILENZIATO) ---
  # Per riattivarlo, cambia FALSE in TRUE
  if (FALSE && any(grepl("refitted", names(res)))) 
  {
    for(i in 1:length(res$allTracks.processed))
    {
      samplename <- names(res$allTracks)[i]
      refit_dir <- file.path(outdir, samplename, "refitted")
      if(!dir.exists(refit_dir)) dir.create(refit_dir, recursive = TRUE)
      
      if(!is_pdf)
        png(file.path(refit_dir, paste0(samplename, "_refitted.png")), width = 5500, height = 2496, res=300)
      if(is_pdf)
        pdf(file.path(refit_dir, paste0(samplename, "_refitted.pdf")), width = 15, height = 5.5)
      
      try({
        suppressWarnings(plotSolution(res$allTracks.processed[[i]],
                                      purity=res$allSolutions.refitted.auto[[i]]$purity,
                                      ploidy=res$allSolutions.refitted.auto[[i]]$ploidy,
                                      ismale=if(!is.null(res$sex)) res$sex[i]=="male" else "female",
                                      allchr=res$chr,
                                      gamma=GAMMA,
                                      ismedian=ismedian,
                                      svinput=if(!is.null(lSVinput)) lSVinput[[i]] else svinput,
                                      sol=res$allSolutions[[i]],
                                      rainbowChr=rainbowChr,
                                      is_pdf=is_pdf))
        title(paste0(samplename, "-refitted"))
      })
      dev.off()
      
      try({
        writeProfile(prof=res$allProfiles.refitted.auto[[i]],
                     samplename=paste0(samplename, "_", projectname, "_refitted"),
                     outdir=refit_dir)
      })
    }
  }
  
  # --- FUNZIONI HELPER E SUMMARY GENERALE ---
  .mytry <- function(x,retVal=NA,...) {
    out <- try(x,silent=T,...)
    if(inherits(out,"try-error")) return(retVal)
    out
  }
  getploidy <- function(tt) {
    tt <- data.frame(chromosome=as.character(tt[,"chromosome"]),
                     start=as.numeric(tt[,"start"]),
                     end=as.numeric(tt[,"end"]),
                     total_copy_number=as.numeric(tt[,"total_copy_number"]))
    sizes <- (tt$end-tt$start)/1000000
    isna <- is.na(sizes) | is.na(tt$total_copy_number)
    sum(tt$total_copy_number[!isna]*sizes[!isna],na.rm=T)/sum(sizes[!isna],na.rm=T)
  }
  
  try({
    res <- append(res,
                  list(summary=list(allSols=data.frame(samplename=names(res$allTracks),
                                                       purity=sapply(res$allSolutions,function(x) .mytry(x$purity)),
                                                       ploidy=sapply(res$allSolutions,function(x) .mytry(x$ploidy)),
                                                       ploidy.tumour=sapply(res$allProfiles,function(x) .mytry(getploidy(x)))),
                                    allSols.refitted=if(!any(grepl("refitted",names(res)))) NULL
                                    else
                                      data.frame(samplename=names(res$allTracks),
                                                 purity=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$purity)),
                                                 ploidy=sapply(res$allSolutions.refitted.auto,function(x) .mytry(x$ploidy)),
                                                 ploidy.tumour=sapply(res$allProfiles.refitted.auto,function(x) .mytry(getploidy(x)))))))
  })
  
  # Writing global tables and Rda object to root folder (outdir)
  try(write.table(res$summary$allSols,
                  file=file.path(outdir, paste0("summary_", projectname, ".txt")),
                  sep="\t", quote=F, col.names=T, row.names=T))
  
  if(!is.null(res$summary$allSols.refitted)) {
    try(write.table(res$summary$allSols.refitted,
                    file=file.path(outdir, paste0("summary_", projectname, "_refitted.txt")),
                    sep="\t", quote=F, col.names=T, row.names=T))
  }
  
  save(res, file=file.path(outdir, paste0("result_object_", projectname, ".Rda")))
  
  return(res)
}