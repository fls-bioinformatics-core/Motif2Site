
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' @title Read a bed file as Genomic Ranges
#' @description Read a bed file as Genomic Ranges.
#' @param fileName A table delimeted file in bed format
#' @return granges format of given coordinates
#' @examples
#'
#' yeastExampleFile=system.file("extdata", "YeastSampleMotif.bed",
#'      package="Motif2Site")
#' ex <- Bed2Granges(yeastExampleFile)
#' ex
#'
#' @export

Bed2Granges  <- function(fileName) {

  if (!(file.exists(fileName))){
    stop(fileName, " file does not exist")
  }

  Table <- utils::read.table(fileName, header=FALSE, stringsAsFactors=FALSE)

  if (length(Table[1,])<3)
  {
    stop("Bed files must have at least three columns")
  }

  if (!(typeof(Table[,2])=="integer")){
    stop("The second column of a bed file must be an integer")
  }

  if (!(typeof(Table[,3])=="integer")){
    stop("The third column of a bed file must be an integer")
  }

  granges <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(Table[,1]),
    ranges=IRanges::IRanges(Table[,2], Table[,3]))
  rm(Table)
  gc()
  return(granges)
}

#' @title Write Genomic Ranges in a bed file
#' @description Write Genomic Ranges in a bed file.
#' @param granges coordinates to write in granges format
#' @param fileName the name of the file that contains the coordianted
#' @return No return value

Granges2Bed  <- function(granges, fileName)
{
  df <- data.frame(chromosome=GenomeInfoDb::seqnames(granges),
                   starts=BiocGenerics::start(granges),
                   ends=BiocGenerics::end(granges)
                   )
  utils::write.table(df,
                     file=fileName,
                     quote=FALSE,
                     sep="\t",
                     row.names=FALSE,
                     col.names=FALSE
                     )
  rm(df)
  gc()
}


#' @title Delete a vector of files
#' @description Delete multiple give files as a vector of characters
#' @param files a vector of files
#' @return No return value

DeleteMultipleFiles <- function(files)
{
  fileNumber <- length(files)

  if(fileNumber>0)
  {
    for (i in seq_len(fileNumber))
    {
      if(file.exists(files[i]))
        file.remove(files[i])
    }

  }
}


#' @title Compare a set of bed files to a provided regions set
#' @description Get combined ranges of bed files and compare them to given
#' @description binding regions in terms of precision/recall.
#' @param motifName a vector of motif names
#' @param bindingRegions granges of provided binding regions
#' @return A dataframe which includes precision recall values for each motif

CompareBeds2GivenRegions <- function(motifName, bindingRegions)
{

  if(!(file.exists("combinedRanges.bed")))
  {
    stop("combinedRanges.bed does not exist")
  }

  combinedRanges <- Bed2Granges("combinedRanges.bed")
  combinedMatrix <- utils::read.table("combinedMatrix",
                                      header=FALSE,
                                      stringsAsFactors=FALSE
                                      )
  combinedOverlap <- GenomicRanges::findOverlaps(combinedRanges,bindingRegions)

  motifNameNumber <- length(motifName)
  if(motifNameNumber>0)
  {
    regionCoverage <-  seq_len(motifNameNumber) * 0
    motifCoverage <-  seq_len(motifNameNumber) * 0
    for (i in seq_len(motifNameNumber))
    {
      motifBindingInds <- which((combinedMatrix[,i]==1)==TRUE)
      rowsIndices <- which((match(S4Vectors::queryHits(combinedOverlap),
                                  motifBindingInds)>0)==TRUE)
      motifCoverage[i] <- length(unique(
        S4Vectors::queryHits(combinedOverlap[rowsIndices])))/
        length(motifBindingInds)
      regionCoverage[i] <- length(unique(
        S4Vectors::subjectHits(combinedOverlap[rowsIndices])))/
        length(bindingRegions)
    }

    graphics::plot(motifCoverage,
                   regionCoverage,
                   ylim=c(0,1),
                   xlim=c(min(motifCoverage)-(min(motifCoverage)/3),
                           (max(motifCoverage)+(max(motifCoverage)/3))),
                   ylab='Recall region based',
                   xlab='Precision motif based',
                   main='motif vs given regions'
    )
    graphics::abline(h=0.95, col="red", lty=2, lwd=3)
    graphics::text(motifCoverage,regionCoverage, motifName, pos=4)

    motifsAnalysis <- data.frame(motifName=motifName,
                                 regionCoverage=regionCoverage,
                                 motifCoverage=motifCoverage
    )

    rm(combinedRanges, combinedMatrix, combinedOverlap, regionCoverage,
       motifCoverage, motifName,motifBindingInds)
    gc() #garbage collection

  } else
  {
    motifsAnalysis <- data.frame(motifName="",
                                 regionCoverage=0,
                                 motifCoverage=0
    )
  }


  return(motifsAnalysis)

}


#' @title Combine motif bed files into a combined ranges
#' @description Get motif file names and combine them into a matrix, and keep
#' the indices of original motifs in the combined file.
#' @description If the motif type is string the bed files are deleted after
#' being combined to one matrix.
#' @param motifFileNames a vector motif file names
#' @param motifType Type of motif string or give bed
#' @return No return value

combineMotifFiles <- function(motifFileNames, motifType="BioString"){

  motifFileNamesNumber <- length(motifFileNames)

  for(i in seq_len(motifFileNamesNumber))
  {
    if(!(file.exists(motifFileNames[i])))
    {
      stop(motifFileNames[i]," does not exist")
    }
  }

  if(motifFileNamesNumber>0)
  {
    combinedRanges <- GenomicRanges::GRanges()
    for (i in seq_len(motifFileNamesNumber))
    {
      assign(paste0("granges_",i), Bed2Granges(motifFileNames[i]))
      combinedRanges <- c(combinedRanges, get(paste0("granges_",i)))
    }
    combinedRanges <- GenomicRanges::reduce(combinedRanges)
    combinedMatrix <- matrix(0,
                             nrow=length(combinedRanges),
                             ncol=length(motifFileNames)
    )
    for (i in seq_len(motifFileNamesNumber))
    {
      combinedMatrix[S4Vectors::subjectHits(
        GenomicRanges::findOverlaps(get(paste0("granges_",i)),combinedRanges)),
        i] <- 1
    }
    utils::write.table(data.frame(combinedMatrix),
                       file="combinedMatrix",
                       quote=FALSE,
                       sep="\t",
                       row.names=FALSE,
                       col.names=FALSE
    )
    Granges2Bed(combinedRanges, "combinedRanges.bed")
    if(motifType=="BioString")
    {
      DeleteMultipleFiles(motifFileNames)
    }
    rm(combinedMatrix, combinedRanges)
    gc() #garbage collection

  }

}


#' @title Compare a set of bed files to a user provided regions set
#' @description This function gets user provided bedfiles and compare them with
#'  a user provided region.
#' @description It returns this comparison to given user binding regions in
#'  terms of precision/recall.
#' @param bedfiles a vector of bed files
#' @param motifnames a vector of the names related to bed files
#' @param givenRegion granges of user provided binding regions
#' @return A dataframe which includes precision recall values for each bed file
#' @examples
#'
#'yeastExampleFile=system.file("extdata", "YeastSampleMotif.bed",
#'                               package="Motif2Site")
#' YeastRegionsChIPseq <- Bed2Granges(yeastExampleFile)
#' bed1 <- system.file("extdata", "YeastBedFile1.bed", package="Motif2Site")
#' bed2 <- system.file("extdata", "YeastBedFile2.bed", package="Motif2Site")
#' BedFilesVector <- c(bed1, bed2)
#' SequenceComparison <- compareBedFiless2UserProvidedRegions(
#'      givenRegion=YeastRegionsChIPseq,
#'      bedfiles=BedFilesVector,
#'      motifnames=c("YeastBed1", "YeastBed2")
#'      )
#' SequenceComparison
#'
#' @seealso
#' \code{\link{compareMotifs2UserProvidedRegions}}
#' @export

compareBedFiless2UserProvidedRegions <-
  function(bedfiles, motifnames, givenRegion)
{

  combineMotifFiles(bedfiles, motifType="UserProvided")
  comparison <- CompareBeds2GivenRegions(motifnames, givenRegion)
  DeleteMultipleFiles(c("combinedRanges.bed","combinedMatrix"))
  gc() #garbage collection
  return(comparison)
}


#' @title Compare a set of motifs to a user provided regions set
#' @description This function gets user provided motifs and related mismatch
#' numbers, it detects motifs and compare them with a user provided region.
#' @description It returns this comparison to given user binding regions in
#' terms of precision/recall.
#' @description The genome and build information should be provided and relevant
#'  BS genomes packages such as BSgenome.Mmusculus.UCSC.mm10 or
#'   BSgenome.Hsapiens.UCSC.hg38 must
#'  be installed for the used genome and builds.
#' @param motifs a vector of motif characters in nucleotide IUPAC format
#' @param mismatchNumbers a vector Number of mismatches allowed to match with
#' motifs
#' @param genome The genome name such as "Hsapiens", "Mmusculus",
#'  "Dmelanogaster"
#' @param genomeBuild The genome build such as "hg38", "hg19", "mm10", "dm3"
#' @param DB The database of genome build. default: "UCSC"
#' @param givenRegion granges of user provided binding regions
#' @param mainCHRs If true only the major chromosome are considered, if FALSE
#' Random, Uncharacterised, and Mithocondrial chromosomes are also considered
#' @return A dataframe which includes precision recall values for each motif
#'
#' @examples
#'
#' # Artificial example in Yeast
#' # install BSgenome.Scerevisiae.UCSC.sacCer3 prior to run this code
#'  yeastExampleFile=system.file("extdata", "YeastSampleMotif.bed",
#'                                 package="Motif2Site")
#' YeastRegionsChIPseq <- Bed2Granges(yeastExampleFile)
#' SequenceComparison <- compareMotifs2UserProvidedRegions(
#'    givenRegion=YeastRegionsChIPseq,
#'    motifs=c("TGATTSCAGGANT", "TGATTCCAGGANT", "TGATWSCAGGANT"),
#'    mismatchNumbers=c(1,0,2),
#'    genome="Scerevisiae",
#'    genomeBuild="sacCer3"
#'    )
#' SequenceComparison
#'
#'
#' @seealso
#' \code{\link{compareBedFiless2UserProvidedRegions}}
#' @export

compareMotifs2UserProvidedRegions <-
  function(motifs, mismatchNumbers, genome, genomeBuild, DB="UCSC",
           givenRegion, mainCHRs=TRUE)
{

  # Compare with Given regions
  motifsNumber <- length(motifs)
  if (length(mismatchNumbers)!=motifsNumber)
  {
    stop("Motifs and mismatchNumber vectors must have the same length")
  }

  for (i in seq_len(motifsNumber))
  {
    findMotifs(motif=motifs[i],
               mismatchNumber=mismatchNumbers[i],
               genome=genome,
               genomeBuild=genomeBuild,
               DB=DB,
               mainCHRs=mainCHRs,
               firstCHR=FALSE,
               MotifLocationName=paste0("Motif_Locations_",
                                        motifs[i], "_",mismatchNumbers[i])
               )
  }
  fileNames <- vector(mode="character", length=motifsNumber)
  for (i in seq_len(motifsNumber))
  {
    fileNames[i] <- paste0("Motif_Locations_",motifs[i], "_",mismatchNumbers[i])
  }
  combineMotifFiles(fileNames, motifType="BioString")
  comparison <- CompareMotifs2GivenRegions(motifs=motifs,
                                           mismatchNumbers=mismatchNumbers,
                                           bindingRegions=givenRegion
                                           )
  DeleteMultipleFiles(c("combinedRanges.bed","combinedMatrix"))
  gc() #garbage collection
  return(comparison)
}


#' @title Comparison motifs locations to a given regions set
#' @description  Comparison of motifs  locations to user provided binding
#' regions.
#' @description It returns this comparison to given user binding regions in
#' terms of precision/recall.
#' @param motifs a vector of motif characters in nucleotide IUPAC format
#' @param mismatchNumbers a vector Number of mismatches allowed to match with
#' motifs
#' @param bindingRegions granges of user provided binding regions
#' @return A dataframe which includes precision recall values for each motif

CompareMotifs2GivenRegions  <- function(motifs, mismatchNumbers, bindingRegions)
{

  combinedRanges <- Bed2Granges("combinedRanges.bed")
  combinedMatrix <- utils::read.table("combinedMatrix",
                                      header=FALSE,
                                      stringsAsFactors=FALSE
                                      )
  combinedOverlap <- GenomicRanges::findOverlaps(combinedRanges,bindingRegions)
  motifNumber <- length(motifs)

  if(motifNumber>0)
  {
    regionCoverage <-  seq_len(motifNumber) * 0
    motifCoverage <-  seq_len(motifNumber) * 0
    motifName <- as.character(motifs)
    for (i in seq_len(motifNumber))
    {
      motifBindingInds <- which((combinedMatrix[,i]==1)==TRUE)
      rowsIndices <- which((match(S4Vectors::queryHits(combinedOverlap),
                                  motifBindingInds)>0)==TRUE)
      motifCoverage[i] <-
        length(unique(S4Vectors::queryHits(combinedOverlap[rowsIndices])))/
        length(motifBindingInds)
      regionCoverage[i] <-
        length(unique(S4Vectors::subjectHits(combinedOverlap[rowsIndices])))/
        length(bindingRegions)
      motifName[i] <- paste0(motifName[i], "_",as.numeric(mismatchNumbers[i]),
                             "mismatch")
    }

    graphics::plot(motifCoverage,
                   regionCoverage,
                   ylim=c(0,1),
                   xlim=c(min(motifCoverage)-(min(motifCoverage)/3),
                           (max(motifCoverage)+(max(motifCoverage)/3))),
                   ylab='Recall region based',
                   xlab='Precision motif based',
                   main='motif vs given regions'
    )
    graphics::abline(h=0.95, col="red", lty=2, lwd=3)
    graphics::text(motifCoverage,regionCoverage, motifName, pos=4)

    motifsAnalysis <- data.frame(motifName=motifName,
                                 regionCoverage=regionCoverage,
                                 motifCoverage=motifCoverage
    )

    rm(combinedRanges, combinedMatrix, combinedOverlap, regionCoverage,
       motifCoverage, motifName,motifBindingInds)
    gc() #garbage collection


  } else
  {
    motifsAnalysis <- data.frame(motifName="",
                                 regionCoverage=0,
                                 motifCoverage=0
    )

  }
  return(motifsAnalysis)

}


#' @title Find motif instances with a certain mismatch number
#' @description Find motif instances in a given genome. It gets motif strings
#' and related allowed mismatchnumbers and returns genomewide motif instances.
#' @description The genome and build information should be provided and relevant
#'  BS genomes packages such as BSgenome.Mmusculus.UCSC.mm10 or
#'   BSgenome.Hsapiens.UCSC.hg38 must be installed for the used genome and
#'   builds.
#' @param motif motif characters in nucleotide IUPAC format
#' @param mismatchNumber Number of mismatch allowed to match with motif
#' @param genome The genome name such as "Hsapiens", "Mmusculus",
#'  "Dmelanogaster"
#' @param genomeBuild The genome build such as "hg38", "hg19", "mm10", "dm3"
#' @param DB The database of genome build. default: "UCSC"
#' @param mainCHRs If true only the major chromosome are considered, if FALSE
#'  Random, Uncharacterised, and Mithocondrial chromosomes are also considered
#' @param firstCHR If true only Chr1 is used to find motifs. Default is FALSE
#' @param MotifLocationName The name of the file of the motif locations
#' @param limitedRegion If specified the motifs are detected in the provided
#' granges
#' @return No return value

findMotifs  <-
  function(motif, mismatchNumber, genome, genomeBuild, DB="UCSC", mainCHRs=TRUE,
           firstCHR=FALSE, MotifLocationName="Motif_Locations",
           limitedRegion=NA)
{

  BSGstring <- paste("BSgenome.", genome,".", DB, ".", genomeBuild,sep="")
  if (!requireNamespace(BSGstring, quietly=TRUE))
  {
    stop(BSGstring, " has not been installed")
  }
  loadNamespace(BSGstring)
  genome <-  BSgenome::getBSgenome(BSGstring)

  # Get accepted crhomosome index
  # If mainCHRs==TRUE: remove chrUn, chrM, and randome choromosomes

  ChromosomeNumber <- length(GenomeInfoDb::seqnames(genome))
  chrInds <- seq_len(ChromosomeNumber)
  if(mainCHRs==TRUE)
  {
    chrInds <-
      setdiff(chrInds,
              union(
                union(grep("random", GenomeInfoDb::seqnames(genome)),
                      grep("chrUn", GenomeInfoDb::seqnames(genome))),
                grep("chrM", GenomeInfoDb::seqnames(genome))
                )
              )
  }

  # get motifs in accepted chromosomes
  motifSites <- GenomicRanges::GRanges()
  if(firstCHR==FALSE)
  {
    ChrNumber <- length(chrInds)
    if(ChrNumber>0)
    {
      for (i in seq_len(ChrNumber))
      {
        chrDNA <-
          Biostrings::maskMotif(
            eval(parse(
              text=paste("genome$",
                         as.character(
                           GenomeInfoDb::seqnames(genome)[chrInds[i]]),
                         sep='')
            )),"N")
        chrMatchedCase <- IRanges::union(
          IRanges::ranges(Biostrings::matchPattern(Biostrings::DNAString(motif),
                                                   chrDNA,
                                                   fixed=FALSE,
                                                   max.mismatch=mismatchNumber
          )
          ),
          IRanges::ranges(
            Biostrings::matchPattern(
              Biostrings::reverseComplement(Biostrings::DNAString(motif)),
              chrDNA,
              fixed=FALSE,
              max.mismatch=mismatchNumber)
          )
        )
        if(length(chrMatchedCase)>0)
        {
          chrRegions <- GenomicRanges::GRanges(
            seqnames=GenomeInfoDb::seqnames(genome)[chrInds[i]],
            ranges=IRanges::IRanges(BiocGenerics::start(chrMatchedCase),
                                      end=BiocGenerics::end(chrMatchedCase))
          )
          motifSites <- GenomicRanges::union(motifSites, chrRegions)
        }
      }
    }

  }else{
    i<-1
    chrDNA <-
      Biostrings::maskMotif(
        eval(parse(
          text=paste("genome$",
                     as.character(GenomeInfoDb::seqnames(genome)[chrInds[i]]),
                     sep='')
          )),"N")
    chrMatchedCase <-
      IRanges::union(
        IRanges::ranges(Biostrings::matchPattern(Biostrings::DNAString(motif),
                                                 chrDNA,
                                                 fixed=FALSE,
                                                 max.mismatch=mismatchNumber)
                        ),
        IRanges::ranges(Biostrings::matchPattern(
          Biostrings::reverseComplement(Biostrings::DNAString(motif)),
          chrDNA,
          fixed=FALSE,
          max.mismatch=mismatchNumber)
                        )
        )
    if(length(chrMatchedCase)>0){
      chrRegions <- GenomicRanges::GRanges(
        seqnames=GenomeInfoDb::seqnames(genome)[chrInds[i]],
        ranges=IRanges::IRanges(BiocGenerics::start(chrMatchedCase),
                                  end=BiocGenerics::end(chrMatchedCase))
        )
      motifSites <- GenomicRanges::union(motifSites, chrRegions)
    }

  }

  if(!is.na(limitedRegion)[1])
  {
    motifSites <- motifSites[unique(
      S4Vectors::queryHits(GenomicRanges::findOverlaps(motifSites,limitedRegion)
                           )
      )]
  }

  Granges2Bed(granges=motifSites, fileName=MotifLocationName)
  rm(motifSites, chrRegions, chrMatchedCase, chrDNA, genome, chrInds,
     BSGstring)
  gc() #garbage collection
}




#' @title build heurisitc distribution around the binding sites
#' @description  This function generates heuristic distribution of short reads
#' around binding sites which do not need to deconvolve, total numer of short
#' reads and window size as number of neucleotid around binding sites.
#' @description It fits a kernel to the distribution and return the distribution
#'  as output. The total sum of returned values is equal to one. It plots this
#'   kernel.
#' @description Also it calculates FRiPs (Fraction of Reads in Peaks) for each
#'
#'   ChIP-seq and returns it. FRiPs and kernel distributions are measures of
#'    goodness of ChIP-seq experiments and selected motifs.
#' @param chipSeq ChIP-seq aligned 1nt short reads
#' @param averageBindings expected short reads number aligned to a random
#'  location of genes of given size
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param acceptedRegionsOutputFile Accepted binding regions
#' @param currentDir Directory for I/O operations
#' @return FRiPs Fraction of Reads in Peaks

deriveHeuristicBindingDistribution <-
  function(chipSeq,averageBindings, windowSize,
           acceptedRegionsOutputFile="BindingRegions", currentDir)
{

    if(!(file.exists(acceptedRegionsOutputFile)))
    {
      stop(acceptedRegionsOutputFile," does not exist")
    }

    if(!(dir.exists(currentDir)))
    {
      stop(currentDir," does not exist")
    }


  acceptedRegionTable <- utils::read.table(acceptedRegionsOutputFile,
                                           header=TRUE,
                                           stringsAsFactors=FALSE)
  acceptedRegion <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(acceptedRegionTable[,1]),
    ranges=IRanges::IRanges(acceptedRegionTable[,2], acceptedRegionTable[,3])
    )
  rm(acceptedRegionTable)
  regionNumber <- length(acceptedRegion)

  oneBindingSiteIndices <-
    which(((BiocGenerics::end(acceptedRegion) -
              BiocGenerics::start(acceptedRegion))==(2*windowSize))==TRUE)
  multiBindingSiteIndices <- setdiff(seq_len(regionNumber), oneBindingSiteIndices)

  replicateNumber <- length(chipSeq$IPfiles)
  FRiPs <- seq_len(replicateNumber)*0

  for (i in seq_len(replicateNumber))
  {
    if(!(file.exists(as.character(chipSeq$IPfiles[i]))))
    {
      stop(as.character(chipSeq$IPfiles[i])," does not exist")
    }

    Table <- utils::read.table(
      file=as.character(chipSeq$IPfiles[i]),
      header=FALSE,
      stringsAsFactors=FALSE
      )
    chipReads <- GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(Table[,1]),
      ranges=IRanges::IRanges(Table[,2], Table[,3])
      )
    rm(Table)
    gc() # garbage collection
    totalChipReads <- length(chipReads) ## total reads count in ChIP

    # Get short reads of ith IP in binding regions
    chipReads <- chipReads[S4Vectors::queryHits(
      GenomicRanges::findOverlaps(chipReads,acceptedRegion))]
    gc() # Garbage collection

    # Calculate FRiPs for ith IP experiment
    FRiPs[i] <- length(chipReads)/totalChipReads

    # Find cumulative binding distribution around nondecomposing binding sites
    nonDecomposingRegions <- acceptedRegion[oneBindingSiteIndices]
    tmpOverlaps <- GenomicRanges::findOverlaps(chipReads,nonDecomposingRegions)
    bindingDistribution <- as.data.frame(
      table(
        c(seq_len(2*windowSize+1),
          (BiocGenerics::end(chipReads[S4Vectors::queryHits(tmpOverlaps)])-
             BiocGenerics::start(
               nonDecomposingRegions[S4Vectors::subjectHits(tmpOverlaps)])
           )
          )
        )
      )
    bindingDistribution$Freq <- bindingDistribution$Freq-1
    totalShortReads <-  length(tmpOverlaps)
    bckgrShortReads <- averageBindings[i]*length(oneBindingSiteIndices)
    freq <- bindingDistribution$Freq/sum(bindingDistribution$Freq)
    kernelDistribution <- fitKernelDensity(
      freq[seq_len(2*windowSize+1)], totalShortReads, windowSize)
    kernelDistribution$Freq <-
      kernelDistribution$Freq - (bckgrShortReads/(2*windowSize+1))
    kernelDistribution$Freq <-
      kernelDistribution$Freq/sum(kernelDistribution$Freq)

    KFnumber <- length(kernelDistribution$Freq)
    if(KFnumber>0)
    {
      grDevices::png(filename=file.path(currentDir,
                       paste0("Heuristic binding distribute_",i)))
      graphics::plot(y=kernelDistribution$Freq,
                     x=seq_len(KFnumber),
                     xlab="nucleotide",
                     ylab="Frequency",
                     ylim=c(0,max(kernelDistribution$Freq)+0.001),
                     main="commulative distirubtion"
      )
      grDevices::dev.off()
    }

    utils::write.table(kernelDistribution,
                       file=file.path(currentDir,
                                      paste0("HeuristicDistribution_", i)),
                       quote=FALSE,
                       sep="\t",
                       row.names=FALSE,
                       col.names=TRUE
                       )


    # make binding vecotr of multi binding sites
    multiBindingSiteIndiceNumber <- length(multiBindingSiteIndices)
    if(multiBindingSiteIndiceNumber>0)
    {
      for(j in seq_len(multiBindingSiteIndiceNumber))
      {
        # Get short reads of ith IP experiment for jth binding region
        overlappingReads <-
          chipReads[
            S4Vectors::queryHits(
              GenomicRanges::findOverlaps(
                chipReads, acceptedRegion[multiBindingSiteIndices[j]])
            )
            ]
        bindingVector <- as.data.frame(
          table(BiocGenerics::end(overlappingReads)-
                  BiocGenerics::start(
                    acceptedRegion[multiBindingSiteIndices[j]]
                    )
                )
        )
        utils::write.table(
          bindingVector,
          file=file.path(currentDir,
            paste0("bindingVector_", i, "_", multiBindingSiteIndices[j])),
          quote=FALSE,
          sep="\t",
          row.names=FALSE,
          col.names=TRUE
        )

      }
    }


  }
  rm(acceptedRegion, oneBindingSiteIndices, multiBindingSiteIndices, chipReads,
     kernelDistribution, overlappingReads,bindingVector, freq)
  gc() # garbage collection
  return(FRiPs)
}





#' @title Fit a kernel density distribution to the obersever heuristic
#'  distribution
#' @description  This function gets heuristic distribution of short reads around
#'  binding sites, total numer of short reads and window size as number of
#'   neucleotid around binding sites.
#' @description It fits a kernel to the distribution and return the distribution
#'  as output. The total sum of returned values is equal to one.
#' @param heuristicDistribution Original short distribution
#' @param totalShortReads Total number of short reads
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @return kernel returns fitted kernel distribution of short reads around
#'  binding sites

fitKernelDensity <- function(heuristicDistribution, totalShortReads, windowSize)
{

  totalNumber <- 100000
  heuristicVector <- rep(1,round(totalNumber*heuristicDistribution[1]))
  for (i in c(2:((2*windowSize)+1)))
  {
    heuristicVector <- c(heuristicVector,
                         rep(i,round(totalNumber*heuristicDistribution[i])))
  }
  kernelDistribution <- stats::density(heuristicVector,
                                       from=1,
                                       to=(2*windowSize+1),
                                       n=(2*windowSize+1)
                                       )
  kernelVector <- rep(1,round(totalNumber*kernelDistribution$y[1]))
  for (i in c(2:((2*windowSize)+1)))
  {
    kernelVector <-
      c(kernelVector, rep(i,round(totalNumber*kernelDistribution$y[i])))
  }

  kernel <- as.data.frame(table(kernelVector))
  tmpTotal <- sum(kernel$Freq)
  kernel$Freq <-  kernel$Freq * (totalShortReads/tmpTotal)

  rm(heuristicVector, kernelVector, kernelDistribution)
  gc()#garbage collection

  return(kernel)

}







#' @title Decompose binding signal among accepted motifs
#' @description Gets motif locations and related short reads and select the
#'  motifs which are non-skewed: abs(skewness) < 0.3 and more short reads binds
#'   closer to site, and show strong binding after decomposition.
#' @description Decomposition is performed by using mixtools normalmixEM command
#'  fixing mu as motif locations.
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param replicateNumber experiment replicate number
#' @param acceptedRegionsOutputFile File name contains binding regions
#'  coordinates and related motifs
#' @param acceptedMotifsOutputFile File name contains motifs coordinates and
#'  related information, Pvalue, FE, etc
#' @param currentDir Directory for I/O operations
#' @return motifStatistics Ratio of accepted motifs, rejected motifs due to
#'  skewnewss, and rejected motifs after decomposition


decomposeBindingSignal <-
  function(windowSize, replicateNumber,
           acceptedRegionsOutputFile="BindingRegions",
           acceptedMotifsOutputFile="BindingMotifsTable", currentDir)
{

    if(!(file.exists(acceptedRegionsOutputFile)))
    {
      stop(acceptedRegionsOutputFile," does not exist")
    }

    if(!(file.exists(acceptedMotifsOutputFile)))
    {
      stop(acceptedMotifsOutputFile," does not exist")
    }

    if(!(dir.exists(currentDir)))
    {
      stop(currentDir," does not exist")
    }

  acceptedRegionTable <- utils::read.table(acceptedRegionsOutputFile,
                                           header=TRUE,
                                           stringsAsFactors=FALSE)
  acceptedRegion <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(acceptedRegionTable[,1]),
    ranges=IRanges::IRanges(acceptedRegionTable[,2], acceptedRegionTable[,3])
    )
  acceptedMotifTable <- utils::read.table(acceptedMotifsOutputFile,
                                          header=TRUE,
                                          stringsAsFactors=FALSE)
  acceptedMotif <-
    GenomicRanges::GRanges(seqnames=S4Vectors::Rle(acceptedMotifTable[,1]),
                           ranges=IRanges::IRanges(acceptedMotifTable[,2],
                                                     acceptedMotifTable[,3])
                           )
  if(replicateNumber>1)
  {
    minIPCount <- min(acceptedMotifTable$NormalizedCountIP)

  }else
  {
    minIPCount <- min(acceptedMotifTable$IP1)
  }
  regionNumber <- (dim(acceptedRegionTable))[1]
  skewedMotifs <- c()
  RejectedMotifs <- c()
  RemainedMotifs <- c()

  for (i in seq_len(regionNumber))
  {
    currentRegionLength <-
      acceptedRegionTable$end[i]-acceptedRegionTable$start[i]
    if(currentRegionLength > (2*windowSize))
    {

      motifLocations <- as.numeric(
        strsplit(as.character(acceptedRegionTable$bindinSites[i]), ",")[[1]])-
        acceptedRegionTable$start[i]
      motifTableIndices <- S4Vectors::queryHits(
        GenomicRanges::findOverlaps(acceptedMotif,acceptedRegion[i]))
      readLocations <- c()
      for (j in seq_len(replicateNumber))
      {
        currentFile <- paste0("bindingVector_", j, "_", i)
        if(!(file.exists(
          file.path(currentDir, currentFile))))
        {
          stop(file.path(currentDir, currentFile), " does not exist")
        }

        tab<- utils::read.table(
          file=file.path(currentDir, currentFile),
          header=TRUE,
          stringsAsFactors=FALSE
          )

        DeleteMultipleFiles(file.path(currentDir, currentFile))
        vec <- rep(tab$Var1, tab$Freq)
        readLocations <- c(readLocations,vec)
      }

      changed <- TRUE
      # Remove motifs with U-shape distributed short read around them
      newMotifLocations <- removeNonBellShapedMotifs(motifLocations,
                                                     readLocations,
                                                     windowSize
                                                     )

      ## If all motifs are skewed: TAKE the strongest motif
      if(length(newMotifLocations)==0)
      {
        newMotifLocations <- strongestMotif(motifLocations,
                                            readLocations,
                                            windowSize
                                            )
      }

      removedSkewedMotifs <-
        motifTableIndices[!is.na(match(
          motifLocations,
          setdiff(motifLocations, newMotifLocations)
          ))]
      skewedMotifs <- union(skewedMotifs, removedSkewedMotifs)
      motifTableIndices <- setdiff(motifTableIndices, removedSkewedMotifs)
      motifLocations <- newMotifLocations

      if(length(motifLocations)<2)
      {
        changed <- FALSE
        RemainedMotifs <- union(RemainedMotifs, motifTableIndices)

      }
      ## Decompose with normal mixuter
      while(changed){

        flag <- TRUE
        while(flag)
        {
          try({
            mixture <- quiet(mixtools::normalmixEM(readLocations,
                                                   mu=motifLocations,
                                                   mean.constr=motifLocations,
                                                   maxit=10))
            flag <- FALSE
          })
        }

        smallestMotif <- motifTableIndices[which.min(mixture$lambda)]

        if(replicateNumber>1)
        {
          if(((acceptedMotifTable$NormalizedCountIP[smallestMotif]*
               min(mixture$lambda))*
              currentRegionLength/(2*windowSize))<minIPCount)
          {

            # Remove the motif
            RejectedMotifs <- union(RejectedMotifs, smallestMotif)
            motifLocations <- setdiff(motifLocations,
                                      motifLocations[which.min(mixture$lambda)]
                                      )
            motifTableIndices <- setdiff(motifTableIndices, smallestMotif)

            # If only one motif remains leave the loop
          } else {
            changed <- FALSE
            RemainedMotifs <- union(RemainedMotifs, motifTableIndices)
          }

        } else {
          if(((acceptedMotifTable$IP1[smallestMotif]*min(mixture$lambda))*
              currentRegionLength/(2*windowSize))<minIPCount)
          {

            # Remove the motif
            RejectedMotifs <- union(RejectedMotifs, smallestMotif)
            motifLocations <-
              setdiff(motifLocations, motifLocations[which.min(mixture$lambda)])
            motifTableIndices <- setdiff(motifTableIndices, smallestMotif)

            # If only one motif remains leave the loop
          } else {
            changed <- FALSE
            RemainedMotifs <- union(RemainedMotifs, motifTableIndices)
          }

        }

        if( length(motifTableIndices)==1)
        {
          changed <- FALSE
          RemainedMotifs <- union(RemainedMotifs, motifTableIndices)
        }
      }

    } else {
      motifTableIndices <- S4Vectors::queryHits(
        GenomicRanges::findOverlaps(acceptedMotif,acceptedRegion[i])
        )
      RemainedMotifs <- union(RemainedMotifs, motifTableIndices)

    }

  }

  acceptedMotifTable <- acceptedMotifTable[RemainedMotifs,]
  acceptedMotif <- acceptedMotif[RemainedMotifs]

  acceptedRegion <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(acceptedMotifTable[,1]),
    ranges=IRanges::IRanges(
      (acceptedMotifTable[,2]+acceptedMotifTable[,3])/2,
      (acceptedMotifTable[,2]+acceptedMotifTable[,3])/2
      )+windowSize
    )
  acceptedRegion <- GenomicRanges::reduce(acceptedRegion)

  bindingRegionSitesIndices <-
    GenomicRanges::findOverlaps(acceptedRegion,  acceptedMotif)

  acceptedRegionNumber <-  length(acceptedRegion)
  sitesVector <- vector(mode="character", length=acceptedRegionNumber)
  if(acceptedRegionNumber>0)
  {
    for(i in seq_len(acceptedRegionNumber))
    {
      tmpBindingIndices <- S4Vectors::subjectHits(bindingRegionSitesIndices)[
        which((S4Vectors::queryHits(bindingRegionSitesIndices)==i)==TRUE)]
      sitesVector[i] <-
        paste0(
          round((BiocGenerics::start(acceptedMotif[tmpBindingIndices])+
                   BiocGenerics::end(acceptedMotif[tmpBindingIndices]))/2
          ),
          collapse=","
        )
    }
  }

  df <- data.frame(chr=as.vector(GenomeInfoDb::seqnames(acceptedRegion)),
                   start=BiocGenerics::start(acceptedRegion),
                   end=BiocGenerics::end(acceptedRegion),
                   bindinSites=sitesVector
                   )
  utils::write.table(df,
                     file=acceptedRegionsOutputFile,
                     row.names=FALSE,
                     col.names=TRUE,
                     quote=FALSE,
                     sep="\t"
                     )

  utils::write.table(acceptedMotifTable,
                     file=acceptedMotifsOutputFile,
                     row.names=FALSE,
                     col.names=TRUE,
                     quote=FALSE,
                     sep="\t"
                     )


  rm(acceptedRegionTable, readLocations, vec, tab, acceptedRegion,
     acceptedMotifTable,acceptedMotif, motifTableIndices, motifLocations,
     bindingRegionSitesIndices, sitesVector)
  gc()

  motifStatistics <- data.frame(skewnessTestRejected=length(skewedMotifs),
                                decompositionRejected=length(RejectedMotifs),
                                accepted=length(RemainedMotifs)
                                )
  return(motifStatistics)
}



#' @title Remove non-bell shpape motifs prior to binding signal decomposition
#' @description Gets motif locations and related short reads and returns the
#' motifs which are non-skewed abs(skewness) < 0.3 and more short reads binds
#'  closer to site.
#' @description  It counts around motif with interval windowSize and
#'  windowSize/2, if the smaller window is less than half of the larger one then
#'   motif is not considered as Bell-shape
#' @param motifLocations A vector of motif locations
#' @param readLocations A vector of 1nt short reads
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @return The coordinates of accepted motifs

removeNonBellShapedMotifs <-
  function (motifLocations, readLocations, windowSize )
{

  maximumMotifReturn <- 15
  skewCutoff <- 0.3
  UshapeRatioCuttoff <- 2
  motifNumber <- length(motifLocations)
  bellShaped <- rep(TRUE, motifNumber)
  skewnessVec <- rep(0,length(motifLocations))

  for (i in seq_len(motifNumber))
  {

    motifDist <-
      readLocations[
        which(((readLocations>motifLocations[i]-(windowSize+1))&
                 (readLocations<motifLocations[i]+(windowSize+1)))==TRUE)]-
      motifLocations[i]

    m3 <- sum((motifDist-mean(motifDist))^3)/length(motifDist)
    s3 <- sqrt(stats::var(motifDist))^3
    if(s3>0){
      skew <- m3/s3
    } else
    {
      skew <- 1000000
    }

    skewnessVec[i] <- skew
    if(abs(skew)> skewCutoff){
      bellShaped[i] <- FALSE
    }

  }

  # Check more sites are in the windows/2 intervals around motifs
  for (i in seq_len(motifNumber)){
    largeInterval <- c(motifLocations[i]-windowSize,
                       motifLocations[i]+windowSize
                       )
    smallInterval <- c(motifLocations[i]-floor(windowSize/2),
                       motifLocations[i]+floor(windowSize/2)
                       )

    # Ushapes
    if((
      (length(which(((readLocations>largeInterval[1])&
                     (readLocations<largeInterval[2]))==TRUE))+1)/
        (length(which(((readLocations>smallInterval[1])&
                       (readLocations<smallInterval[2]))==TRUE))+1)
        )>UshapeRatioCuttoff)
    {
      bellShaped[i] <- FALSE
    }
  }

  motifLocations <- motifLocations[bellShaped]
  skewnessVec <- abs(skewnessVec[bellShaped])

  if(length(motifLocations)>maximumMotifReturn)
  {
    motifLocations <-
      motifLocations[order(skewnessVec)[seq_len(maximumMotifReturn)]]
  }

  return(motifLocations)

}


#' @title Returns the motif with the highest count
#' @description Gets motif locations and related short reads and returns the
#'  motif which include the highest number of short reads around it.
#' @param motifLocations A vector of motif locations
#' @param readLocations A vector of 1nt short reads
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @return The strongest motif

strongestMotif <- function (motifLocations, readLocations, windowSize){

  motifLocationNumber <- length(motifLocations)
  if(motifLocationNumber>0)
  {
    tmpCount <- rep(0, motifLocationNumber)
    for (i in seq_len(motifLocationNumber)){
      largeInterval <-
        c(motifLocations[i]-windowSize, motifLocations[i]+windowSize)
      tmpCount[i] <-
        length(which(((readLocations>largeInterval[1])&
                        (readLocations<largeInterval[2]))==TRUE)
        )
    }

    return(motifLocations[which.max(tmpCount)])
  } else {
    return(0)
  }

}




#' @title Model IP and Input count values with negative Binomal
#' @description Using edgeR TMM normalization and estimating dispersion as well
#'  as Adapting exact test function from edgeR to model IP vs Input counts.
#' @description To make this function memory effcient motifs into smaller sets
#'  and compute them seperately and combine them at the end.
#' @param countTableFile Table of counts which contains all IP and Input value
#'  raw counts
#' @param replicateNumber experiment replicate number
#' @param outputFile The name of the output file generated by this function
#' @param currentDir Directory for I/O operations
#' @return A dataframe includes fold enrichment, pvalue, and normalized count
#'  values

motifBindingNegativeBinomialCount <-
  function(countTableFile, replicateNumber, outputFile, currentDir)
{

    if(!(file.exists(countTableFile)))
    {
      stop(countTableFile," does not exist")
    }

    if(!(dir.exists(currentDir)))
    {
      stop(currentDir, " does not exist")
    }

  countTable <- utils::read.table(countTableFile,
                                  header=FALSE,
                                  stringsAsFactors=FALSE
                                  )
  totalMotifNumber <- nrow(countTable)
  breakSize <- round(480000/replicateNumber)
  if(breakSize<totalMotifNumber)
  {
    iterationPoints <- 
      c(0,seq_len((round(totalMotifNumber/breakSize)-1))*breakSize)
    iterationNumber <- length(iterationPoints)-1
    iterationPoints[iterationNumber+ 1] <- totalMotifNumber

  } else {
    iterationPoints <- c(0, totalMotifNumber)
    iterationNumber <- 1
  }

  outputDF <- data.frame(logFE=seq_len(totalMotifNumber)*0,
                         PValue=seq_len(totalMotifNumber)*0,
                         normalizedCountIP=seq_len(totalMotifNumber)*0,
                         normalizedCountInput=seq_len(totalMotifNumber)*0
                         )

  for(i in seq_len(iterationNumber))
  {
    # make DGE list table of edgeR
    tableOfCounts <- edgeR::DGEList(
      counts=
        countTable[(iterationPoints[i]+1):iterationPoints[i+1],
                   seq_len((2*replicateNumber))],
      group=c(rep("IP",replicateNumber),rep("Input",replicateNumber)),
      remove.zeros=FALSE
      )
    gc()

    TotalCountsTable <- utils::read.table(
      file.path(currentDir, "TotalCounts"),
      header=TRUE,
      stringsAsFactors=FALSE
      )

    tableOfCounts$samples$lib.size[seq_len(replicateNumber)] <-
      TotalCountsTable$ipTotalCount

    tableOfCounts$samples$lib.size[(replicateNumber+1):(2*replicateNumber)] <-
      TotalCountsTable$inputTotalCount

    tableOfCounts <- edgeR::calcNormFactors(tableOfCounts, method="TMM")
    gc()

    tableOfCounts <- edgeR::estimateCommonDisp(tableOfCounts)
    gc()

    tableOfCounts <- edgeR::estimateTagwiseDisp(tableOfCounts)
    gc()

    testWithReplicate <- NegativeBinomialTestWithReplicate(tableOfCounts)
    gc()

    outputDF[(iterationPoints[i]+1):iterationPoints[i+1],] <-
      testWithReplicate$table

  }
  rm(countTable)
  gc()

  utils::write.table(outputDF,
                     file=outputFile,
                     col.names=TRUE,
                     row.names=FALSE,
                     quote=FALSE
                     )
  df <- data.frame(commonDispersion=tableOfCounts$common.dispersion,
                   pseudoLibSize=tableOfCounts$pseudo.lib.size
                   )
  rm(testWithReplicate,tableOfCounts)
  gc()
  return(df)
}


#' @title Negative binomial test of binding using all replicates
#' @description Adapted exact test function from edgeR to compare IP vs Input
#'  with replicates. Input is a DGELIST with common and tag-wise dispression has
#'   been already caluclated by edgeR commands.
#' @description It calculates abundaces with mglmOneGroup identical to edgeR.
#'  logFE was  calculated identiacl to edgeR. For the pvalue test negative
#'   binomial test is performed on the calculated abundance.
#' @param object Table of counts which contains all IP and Input value counts,
#'  TMM normalized and contains dispersion values
#' @param prior.count edgeR prior value
#' @return log fold enrichment, pvalue, and normalized count values

NegativeBinomialTestWithReplicate <- function(object, prior.count=0.125)
{
  # Check inputo
  if(!methods::is(object,"DGEList"))
    stop("Currently only supports DGEList objects as the object argument.")

  # Get group names
  group <- as.factor(object$samples$group)
  pair <- c("Input", "IP" )
  dispersion <- edgeR::getDispersion(object)
  if(is.null(dispersion)) stop("specified dispersion not found in object")
  if(is.na(dispersion[1])) stop("dispersion is NA")
  ldisp <- length(dispersion)
  ntags <- nrow(object$counts)
  if(ldisp!=1 && ldisp!=ntags)
    stop("Dispersion provided by user must have length either 1 or the number
         of tags in the DGEList object.")
  if(ldisp==1) dispersion <- rep(dispersion,ntags)

  # Reduce to two groups
  group <- as.character(group)
  j <- group %in% pair
  y <- object$counts[,j,drop=FALSE]
  lib.size <- object$samples$lib.size[j]
  norm.factors <- object$samples$norm.factors[j]
  group <- group[j]
  if(is.null(rownames(y))) rownames(y) <- paste("tag", seq_len(ntags), sep=".")

  # Normalized library sizes
  lib.size <- lib.size * norm.factors
  offset <- log(lib.size)
  lib.size.average <- exp(mean(offset))

  # logFE
  prior.count <- prior.count*lib.size/mean(lib.size)
  offset.aug <- log(lib.size+2*prior.count)
  j1 <- group==pair[1]
  n1 <- sum(j1)
  if(n1==0) stop("No libraries for",pair[1])
  y1 <- y[,j1,drop=FALSE]
  abundance1 <- edgeR::mglmOneGroup(
    y1+matrix(prior.count[j1],ntags,n1,byrow=TRUE),
    offset=offset.aug[j1],
    dispersion=dispersion
    )
  j2 <- group==pair[2]
  n2 <- sum(j2)
  if(n1==0) stop("No libraries for",pair[2])
  y2 <- y[,j2,drop=FALSE]
  abundance2 <- edgeR::mglmOneGroup(
    y2+matrix(prior.count[j2],ntags,n2,byrow=TRUE),
    offset=offset.aug[j2],
    dispersion=dispersion
    )
  logFE <- (abundance2-abundance1)/log(2)

  normalizedCountIP <- exp(abundance2)*object$pseudo.lib.size
  normalizedCountInput <- exp(abundance1)*object$pseudo.lib.size
  rm(abundance1, abundance2)
  gc() # garbage collection

  # Find size and mean with MASS package and calculate nbinomial distribution
  normalizedCountIP <- round(normalizedCountIP)

  # Don't consider upbinding to find nb distribution
  eps <- 0.001
  log2normalizedCountIP <- log2(normalizedCountIP + eps)
  upperbound <-
    ceiling(2^(stats::quantile(log2normalizedCountIP, 0.75)+
                 1.5*stats::IQR(log2normalizedCountIP))
            )

  gc()

  normalizedCountIPNumber <- length(normalizedCountIP)
  if(normalizedCountIPNumber<1)
    stop("no data to fit by negative binomial")

  tmpDist <- normalizedCountIP[
    setdiff(seq_len(normalizedCountIPNumber),
            which((normalizedCountIP>upperbound)==TRUE))
    ]
  meanTmpDist <- mean(tmpDist)
  varTmpDist <- stats::var(tmpDist)
  size <-
    if(varTmpDist > meanTmpDist) meanTmpDist^2/(varTmpDist - meanTmpDist)
    else  meanTmpDist*2
  flag <- TRUE
  try({
#    distNB <- quiet(suppressWarnings(MASS::fitdistr(tmpDist,
#                             densfun="negative binomial",
#                             start=list(size=size, mu=meanTmpDist)
#                             )))
    distNB <- quiet(MASS::fitdistr(tmpDist,
                                   densfun="negative binomial",
                                   start=list(size=size, mu=meanTmpDist)
                                   )
                    )
    flag <- FALSE
  })

  if(flag)
  {
 #   distNB <- quiet(suppressWarnings(
#      MASS::fitdistr(tmpDist, densfun="negative binomial")))
    distNB <- quiet(MASS::fitdistr(tmpDist, densfun="negative binomial"))
  }

  rm(tmpDist,meanTmpDist)
  gc()
  pvals <- stats::pnbinom(normalizedCountIP,
                          size=(distNB$estimate)[1],
                          mu=(distNB$estimate)[2],
                          lower.tail=FALSE
                          )
  normalizedCountIP <- normalizedCountIP/100

  de.out <- data.frame(logFE=logFE,
                       PValue=pvals,
                       normalizedCountIP=normalizedCountIP,
                       normalizedCountInput=normalizedCountInput
                       )
  rm(distNB, logFE, pvals, normalizedCountIP,
     normalizedCountInput, log2normalizedCountIP)
  gc() # garbage collection
  rn <- rownames(object$counts)
  if(!is.null(rn)) rownames(de.out) <- make.unique(rn)
  methods::new("DGEExact",
               list(table=de.out, comparison=pair, genes=object$genes)
               )
}




#' @title FDR cut-off detection Benjamini Hochberg method
#' @description Return FDR cut-off for a user provided fdrvalue using Benjamini
#'  Hochberg on main motif test data
#' @param TestTableFile test table which contains pvalues
#' @param fdrValue FDR cut-off
#' @return pvalue cut-off

DetectFdrCutoffBH  <- function(TestTableFile="TestResults", fdrValue=0.05) {

  if(!(file.exists(TestTableFile)))
  {
    stop(TestTableFile, " does not exist")
  }

  TestTable <-
    utils::read.table(TestTableFile, header=TRUE, stringsAsFactors=FALSE)
  # Remove NA lines
  pvals <- TestTable$Pvalue[which(!is.na(TestTable$Pvalue)==TRUE)]
  rm(TestTable)
  gc()

  pvalsAdjusted <- stats::p.adjust(pvals, method="BH")
  accepted <- which((pvalsAdjusted<fdrValue)==TRUE)
  if(length(accepted)==0)
  {
    stop("No motif has been accepted. Use higher fdrValue")
  }
  return(max(pvals[accepted], na.rm=NA))

}




#' @title count short reads around motifs for all ChIP-seq experiments
#' @description count short reads related to each motif for all ChIPseq files
#'  both IP and Input.
#' @param motifFile File contains motifs
#' @param chipSeq dataframe of ChIP-seq 1nt alignment location
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param outputName Name of the output table
#' @param currentDir Directory for I/O operations
#' @return No return value

motifCount  <- function(motifFile, chipSeq, windowSize, outputName, currentDir)
{

  if(!(file.exists(motifFile)))
  {
    stop(motifFile, " does not exist")
  }

  if(!(dir.exists(currentDir)))
  {
    stop(currentDir, " does not exist")
  }

  if(!(file.exists(as.character(chipSeq$IPfiles[1]))))
  {
    stop("IP file does not exist")
  }


  if(!(file.exists(as.character(chipSeq$BackgroundFiles[1]))))
  {
    stop("background file does not exist")
  }


  # short reads count for motifs - each IP and Input
  replicateNumber <- length(chipSeq$IPfiles)
  ipTotalCount <- seq_len(replicateNumber)*0
  inputTotalCount <- seq_len(replicateNumber)*0
  for (i in seq_len(replicateNumber))
  {

    ipTotalCount[i] <- motifChipCount(
      motifFile=motifFile,
      chipFile=as.character(chipSeq$IPfiles[i]),
      windowSize=windowSize,
      outputName=file.path(currentDir,  paste0("count_IP", i))
      )

    inputTotalCount[i] <-
      motifChipCount(
        motifFile=motifFile,
        chipFile=as.character(chipSeq$BackgroundFiles[i]),
        windowSize=windowSize,
        outputName=file.path(currentDir,  paste0("count_Input", i))
        )

    computeFoldEnrichment(
      ipCountFile=file.path(currentDir,  paste0("count_IP", i)),
      inputCountFile=
        file.path(currentDir, paste0("count_Input", i)),
      ipTotalCount=ipTotalCount[i],
      inputTotalCount=inputTotalCount[i],
      outputName=file.path(currentDir, paste0("FE",i))
      )

  }

  df <- data.frame(ipTotalCount=ipTotalCount, inputTotalCount=inputTotalCount)

  combine2Table(outputName=outputName,
                replicateNumber=replicateNumber,
                currentDir=currentDir)

  utils::write.table(
    df,
    file=file.path(currentDir, "TotalCounts"),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=TRUE
    )
  gc()

}




#' @title  count short reads related to each motif for a given ChIPseq file
#' @description count 1nt short reads related to each motif for a given ChIPseq
#'  file.
#' @param motifFile File contains motifs
#' @param chipFile ChIP-seq 1nt alignment locations in bed format
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param outputName Name of the output table
#' @return Total number of short reads in motif reagions

motifChipCount  <- function(motifFile, chipFile, windowSize, outputName)
{

  if(!(file.exists(motifFile)))
  {
    stop(motifFile, " does not exist")
  }

  if(!(file.exists(chipFile)))
  {
    stop(chipFile, " does not exist")
  }

  # 1nt alignmetn bed to granges
  Table <-
    utils::read.table(file=chipFile, header=FALSE, stringsAsFactors=FALSE)
  chipReads <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(Table[,1]),
    ranges=IRanges::IRanges(Table[,2], Table[,3]))
  # Remove Table data
  rm(Table)
  gc() # garbage collection
  totalChipReads <- length(chipReads) ## total reads count in ChIP

  # motif bed to granges of windowsize
  Table <- utils::read.table(motifFile, header=FALSE, stringsAsFactors=FALSE)
  motifWindow <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(Table[,1]),
    IRanges::IRanges(((Table[,2]+Table[,3])/2),
                     ((Table[,2]+Table[,3])/2)
                     )+windowSize
    )
  rm(Table)
  gc() # garbage collection
  motifCounts <- GenomicRanges::countOverlaps(motifWindow, chipReads)
  df <- data.frame(motifCounts=motifCounts)
  utils::write.table(df,
                     file=outputName,
                     col.names=FALSE,
                     row.names=FALSE,
                     quote=FALSE
                     )
  rm(motifWindow, chipReads,  motifCounts)
  gc() # garbage collection
  return(totalChipReads)
}





#' @title  compute fold enrichment values for an experiment
#' @description Open raw counts IP and Inut files and with given total counts
#'  calculate Fold Enrichment values for the motifs
#' @param ipCountFile File contains motifs count values for IP experiment
#' @param inputCountFile File contains motifs count values for Input experiment
#' @param ipTotalCount Total short reads number in IP experiment
#' @param inputTotalCount Total short reads number in Input experiment
#' @param outputName Name of the output table
#' @return No return value

computeFoldEnrichment  <-
  function(ipCountFile, inputCountFile, ipTotalCount,
           inputTotalCount, outputName)
{

    if(!(file.exists(ipCountFile)))
    {
      stop(ipCountFile, " does not exist")
    }

    if(!(file.exists(inputCountFile)))
    {
      stop(inputCountFile, " does not exist")
    }

  ipCounts <- (utils::read.table(file=ipCountFile,
                                 header=FALSE,
                                 stringsAsFactors=FALSE)
               )[,1]
  inputCounts <- (utils::read.table(file=inputCountFile,
                                    header=FALSE,
                                    stringsAsFactors=FALSE)
                  )[,1]
  eps <- 0.00001 # to avoid division by zero
  FEs <- ((ipCounts+eps)/ipTotalCount)/((inputCounts+eps)/inputTotalCount)
  df <- data.frame(FEs=FEs)
  utils::write.table(df,
                     file=outputName,
                     col.names=FALSE,
                     row.names=FALSE,
                     quote=FALSE
                     )
}




#' @title  Combine all IP and Input count table files
#' @description Open raw counts IP and Inut files and with given total counts
#'  calculate Fold Enrichment values, and combine them into one file
#' @param replicateNumber Number of the replicates
#' @param outputName Name of the output table
#' @param currentDir Directory for I/O operations
#' @return No return value

combine2Table  <- function(outputName,replicateNumber,currentDir)
{

  if(!(dir.exists(currentDir)))
  {
    stop(currentDir, " does not exist")
  }

  if(!(file.exists( file.path(currentDir, paste0("count_IP", 1)))))
  {
    stop("combine2Table count_IP1 does not exist")
  }

  ipCounts <-
    matrix(0,
           nrow=length((utils::read.table(
             file=file.path(currentDir, paste0("count_IP", 1)),
             header=FALSE,
             stringsAsFactors=FALSE))[,1]),
           ncol=replicateNumber
           )
  inputCount <- ipCounts
  FE <- ipCounts
  colnames(ipCounts) <- 
    paste0(rep("IP",replicateNumber), seq_len(replicateNumber))
  colnames(inputCount)<-
    paste0(rep("Input",replicateNumber), seq_len(replicateNumber))
  colnames(FE)<- paste0(rep("FE",replicateNumber), seq_len(replicateNumber))
  for (i in seq_len(replicateNumber))
  {

    if(!(file.exists( file.path(currentDir, paste0("count_IP", i)))))
    {
      stop("combine2Table count_IP", i, "does not exist")
    }
    ipCounts[,i] <- (utils::read.table(
      file=file.path(currentDir, paste0("count_IP", i)),
      header=FALSE,
      stringsAsFactors=FALSE)
      )[,1]

    if(!(file.exists( file.path(currentDir, paste0("count_Input", i)))))
    {
      stop("combine2Table count_Input", i, "does not exist")
    }
    inputCount[,i] <- (utils::read.table(
      file=file.path(currentDir, paste0("count_Input", i)),
      header=FALSE,
      stringsAsFactors=FALSE)
      )[,1]


    if(!(file.exists( file.path(currentDir, paste0("FE", i)))))
    {
      stop("combine2Table FE", i, "does not exist")
    }
    FE[,i] <- (utils::read.table(
      file.path(currentDir, paste0("FE",i)),
      header=FALSE,
      stringsAsFactors=FALSE)
      )[,1]
    DeleteMultipleFiles(
      c(file.path(currentDir, paste0("count_IP", i)),
        file.path(currentDir, paste0("count_Input", i)),
        file.path(currentDir, paste0("FE",i)))
      )
  }

  combinedMatrix <- c(data.frame(ipCounts),
                      data.frame(inputCount),
                      data.frame(FE)
                      )

  utils::write.table(combinedMatrix,
                     file=outputName,
                     col.names=TRUE,
                     row.names=FALSE,
                     quote=FALSE
                     )

}


#' @title Process count data and perform negative binomial test
#' @description Remove unmmaped regions, low and high binding regions and
#'  regions without fold change, and call negative binomial or nb test for
#'   the remaining regions.
#' @param countTableFile Tabl of count values around motifs for all ChIP-seq
#'  experiments
#' @param outFile The name of the output file
#' @param currentDir Directory for I/O operations
#' @return sequencingStatitics A dataframe consists of the ratio of
#'  non-sequenced, low-sequenced, ang high-sequenced regions.

motifTablePreProcess <- function(countTableFile,  outFile, currentDir)
{

  if(!(dir.exists(currentDir)))
  {
    stop(currentDir, " does not exist")
  }

  if(!(file.exists(countTableFile)))
  {
    stop(countTableFile," does not exist")
  }

  outVectorName <- file.path(currentDir, "motifBindingTestVector")
  outTableName <- file.path(currentDir, "reducedTableMotifs")
  countTable <-
    utils::read.table(countTableFile, header=TRUE, stringsAsFactors=FALSE)
  replicateNumber <- ncol(countTable)/3
  motifNumber <- nrow(countTable)


  ## Mark regions with no short reads mapped as nonseqeunced
  nonSequencedRegions <- which((countTable[,1]==0)==TRUE)
  for (i in seq(from=2, to=(2*replicateNumber)))
  {
    nonSequencedRegions <-
      intersect(nonSequencedRegions, which((countTable[,i]==0)==TRUE))
  }
  acceptableIndices <- setdiff(seq_len(motifNumber), nonSequencedRegions)


  ## Compute low and high binding threshold with quantiles in log2 scale
  upperbound <- seq_len(2*replicateNumber)*0
  lowerbound <- seq_len(2*replicateNumber)*0
  extremeUpperbound <- seq_len(2*replicateNumber)*0
  extremeLowerbound <- seq_len(2*replicateNumber)*0
  eps <- 1
  log2countTable <-
    log2(countTable[acceptableIndices, seq_len(2*replicateNumber)]+eps)
  for (i in seq_len(2*replicateNumber))
  {
    upperbound[i] <-
      stats::quantile(log2countTable[,i], 0.75)+
      1.5*stats::IQR(log2countTable[,i])
  }


  ## Remove high Input in any replicate
  upBindingInds <-
    which(
      (countTable[,replicateNumber+1]>(2^upperbound[replicateNumber+1])
       )==TRUE
      )
  for (i in c((replicateNumber+1):(2*replicateNumber)))
  {
    upBindingInds <-
      union(upBindingInds, which((countTable[,i]>(2^upperbound[i]))==TRUE))
  }
  acceptableIndices <- setdiff(acceptableIndices, upBindingInds)


  ## Remove extreme low in IP ChIP-seqs
  combineRowSum <-
    log2(rowSums(countTable[acceptableIndices, seq_len(2*replicateNumber)]))
  lowerbound <-
    stats::quantile(combineRowSum, 0.25)-1.5*stats::IQR(combineRowSum)
  lowBindingInds <-
    which(
      (rowSums(countTable[,seq_len(2*replicateNumber)])<(2^lowerbound))==TRUE
      )
  acceptableIndices <- setdiff(acceptableIndices, lowBindingInds)

  lowBindingInds <- setdiff(lowBindingInds, nonSequencedRegions)

  testVector <- vector(mode="character", length=motifNumber)
  testVector[seq_len(motifNumber)] <- "Tested"
  testVector[lowBindingInds] <- "UnderBound"
  testVector[upBindingInds] <- "OverBound"
  testVector[nonSequencedRegions] <- "NonSequenced"

  df <- data.frame(testVector=testVector)
  utils::write.table(df,
                     file=outVectorName,
                     col.names=FALSE,
                     row.names=FALSE,
                     quote=FALSE
                     )
  rm(testVector, df, log2countTable,combineRowSum)
  gc()
  utils::write.table(
    data.frame(countTable[acceptableIndices, seq_len(2*replicateNumber)]),
    file=outTableName,
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
    )
  rm(countTable, acceptableIndices)
  gc()

  sequencingStatitics <- data.frame(
    nonSequenced=length(nonSequencedRegions)/motifNumber,
    underBinding=length(lowBindingInds)/motifNumber,
    overBinding=length(upBindingInds)/motifNumber
    )
  rm(upBindingInds, lowBindingInds, nonSequencedRegions)
  gc()
  outputFile <- file.path(currentDir, "BinomialCount")

  if(replicateNumber>1)
  {

    NBstat <- motifBindingNegativeBinomialCount(countTableFile=outTableName,
                                                replicateNumber=replicateNumber,
                                                outputFile=outputFile,
                                                currentDir=currentDir
                                                )

    testVector <- (utils::read.table(outVectorName,header=FALSE,
                                     stringsAsFactors=FALSE))[,1]
    pvalNumber <- length(testVector)
    pvals <- vector(mode="numeric", length=pvalNumber)

    if(pvalNumber>0)
    {
      pvals[seq_len(pvalNumber)] <- NA
    }

    FEs <- pvals
    normalizedCountIP <- pvals
    normalizedCountInput <- pvals
    acceptableIndices <- which((testVector=="Tested")==TRUE)
    testResults <- utils::read.table(outputFile,header=TRUE,
                                     stringsAsFactors=FALSE)
    FEs[acceptableIndices] <- 2^testResults$logFE
    pvals[acceptableIndices] <- testResults$PValue
    normalizedCountIP[acceptableIndices] <- testResults$normalizedCountIP
    normalizedCountInput[acceptableIndices] <- testResults$normalizedCountInput
    df <- data.frame(Test=testVector,
                     FoldEnrichment=FEs,
                     Pvalue=pvals,
                     NormalizedCountIP=normalizedCountIP,
                     NormalizedCountInput=normalizedCountInput
                     )
    utils::write.table(df, file=outFile, col.names=TRUE,
                       row.names=FALSE, quote=FALSE)
    rm(df, FEs, pvals, testVector, acceptableIndices,
       normalizedCountInput, normalizedCountIP)
    gc()
    DeleteMultipleFiles(c(outputFile, outVectorName, outTableName))
    sequencingStatitics <- data.frame(
      nonSequenced=sequencingStatitics$nonSequenced,
      underBinding=sequencingStatitics$underBinding,
      overBinding=sequencingStatitics$overBinding
      )
  } else {
    countTable <- utils::read.table(
      file.path(currentDir, "reducedTableMotifs"),
      header=FALSE,
      stringsAsFactors=FALSE
      )

    eps <- 1
    log2normalizedCountIP <- log2(countTable[,1] + eps)
    upperbound <- 2^(stats::quantile(log2normalizedCountIP, 0.75)
                     +1.5*stats::IQR(log2normalizedCountIP))

    dataPointNumber <-length(countTable[,1])
    if(dataPointNumber<2)
      stop("enough motifs does not exist to fit negative binomial curve")

#    distNB <- quiet(suppressWarnings(MASS::fitdistr(
#      countTable[setdiff(seq_len(dataPointNumber),
#                         which((countTable[,1]>upperbound)==TRUE)),1],
#      densfun="negative binomial"
#      )))
    
    
    distNB <- quiet(MASS::fitdistr(
      countTable[setdiff(seq_len(dataPointNumber),
                         which((countTable[,1]>upperbound)==TRUE)),1],
      densfun="negative binomial"
    ))
    
    testedPvals <- stats::pnbinom(
      countTable[,1],
      size=(distNB$estimate)[1],
      mu=(distNB$estimate)[2],
      lower.tail=FALSE
      )
    rm(countTable, distNB, log2normalizedCountIP)
    gc()

    testVector <- (utils::read.table(
      outVectorName,header=FALSE, stringsAsFactors=FALSE))[,1]

    pvalNumber <- length(testVector)
    pvals <- vector(mode="numeric", length=pvalNumber)
    if(pvalNumber>0)
    {
      pvals[seq_len(pvalNumber)] <- NA
    }

    acceptableIndices <- which((testVector=="Tested")==TRUE)
    pvals[acceptableIndices] <- testedPvals
    df <- data.frame(Test=testVector, Pvalue=pvals)
    utils::write.table(
      df, file=outFile, col.names=TRUE, row.names=FALSE, quote=FALSE)
    rm(df,  pvals, testVector, acceptableIndices, testedPvals)
    gc()
    DeleteMultipleFiles(c(outputFile, outVectorName, outTableName))
    sequencingStatitics <- data.frame(
      nonSequenced=sequencingStatitics$nonSequenced,
      underBinding=sequencingStatitics$underBinding,
      overBinding=sequencingStatitics$overBinding
      )
  }
  return(sequencingStatitics)
}



#' @title  Combine count Table and statistics table
#' @description Combine count table and pvalue FE statistics into one file for
#'  motifs and regions seperately.
#' @param motifFile File contains motifs
#' @param acceptedMotifsOutputFile File name of accepted motif table inforation
#' @param acceptedRegionsOutputFile File name of accepteted region information
#' @param countTableFile Table of count values file name
#' @param testTableFile negative binomial test table file name
#' @param fdrCutoff Pvalue cut-off related to the used FDR
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @return The average binding intensity for each ChIP-seq

combineTestResults  <-
  function(motifFile, acceptedMotifsOutputFile, acceptedRegionsOutputFile,
           countTableFile, testTableFile, fdrCutoff, windowSize)
{

    if(!(file.exists(testTableFile)))
    {
      stop(testTableFile," does not exist")
    }

    if(!(file.exists(countTableFile)))
    {
      stop(countTableFile," does not exist")
    }

    if(!(file.exists(motifFile)))
    {
      stop(motifFile," does not exist")
    }

  testResults <- utils::read.table(
    file=testTableFile, header=TRUE, stringsAsFactors=FALSE)
  motifNumber <- nrow(testResults)
  testedIndices <- which((testResults$Test=="Tested")==TRUE)
  #select the accepted motifs rows
  acceptedMotifsIndices <- which((testResults$Pvalue<=fdrCutoff)==TRUE)
  testedIndices <- setdiff(testedIndices, acceptedMotifsIndices)
  testResults <- testResults[acceptedMotifsIndices,]
  gc()#garbage collection
  # Read the count table
  countTable <- utils::read.table(
    file=countTableFile, header=TRUE, stringsAsFactors=FALSE)
  replicateNumber <- ncol(countTable)/3
  means <- seq_len(replicateNumber) * 0
  for (i in seq_len(replicateNumber)){
    means[i] <- mean(countTable[testedIndices,i])
  }
  #select the accepted motifs rows
  countTable <- countTable[acceptedMotifsIndices,]
  gc()#garbage collection
  # Read the motifs table
  motifs <- utils::read.table(
    file=motifFile, header=FALSE, stringsAsFactors=FALSE)
  #select the accepted motifs rows
  motifs  <- motifs[acceptedMotifsIndices,]
  gc()#garbage collection
  if(replicateNumber>1)
  {
    df <- data.frame(
      chr=motifs[,1],
      start=motifs[,2],
      end=motifs[,3],
      pval=testResults$Pvalue,
      FE=testResults$FoldEnrichment,
      NormalizedCountIP=testResults$NormalizedCountIP,
      NormalizedCountInput=testResults$NormalizedCountInput
      )

  } else {
    df <- data.frame(
      chr=motifs[,1],
      start=motifs[,2],
      end=motifs[,3],
      pval=testResults$Pvalue
      )
  }
  df <- as.data.frame(c(df, countTable))
  utils::write.table(
    df, file=acceptedMotifsOutputFile,
    row.names=FALSE,
    col.names=TRUE,
    quote=FALSE
    )
  rm(df,countTable,acceptedMotifsIndices,testResults, testedIndices)
  gc()#garbage collection
  bindingSites <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(motifs[,1]),
    ranges=IRanges::IRanges(motifs[,2], motifs[,3])
    )
  bindingRegions <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(motifs[,1]),
    ranges=IRanges::IRanges(
      (motifs[,2]+motifs[,3])/2, (motifs[,2]+motifs[,3])/2)+windowSize
    )
  rm(motifs)
  gc() #garbage collection

  # Find binding sites and binding regions and build the site-region table
  bindingRegions <- GenomicRanges::reduce(bindingRegions)
  bindingRegionSitesIndices <-
    GenomicRanges::findOverlaps(bindingRegions, bindingSites)
  sitesVector <- rep("", length(bindingRegions))

  bindingRegionNumber <- length(bindingRegions)
  if(bindingRegionNumber>0)
  {
    for(i in seq_len(bindingRegionNumber))
    {
      tmpBindingIndices <-
        S4Vectors::subjectHits(bindingRegionSitesIndices)[
          which((S4Vectors::queryHits(bindingRegionSitesIndices)==i)==TRUE)]
      sitesVector[i] <- paste(round(
        (BiocGenerics::start(bindingSites)[tmpBindingIndices]+
           BiocGenerics::end(bindingSites)[tmpBindingIndices])/2),
        collapse=",")
    }
  }


  df <- data.frame(
    chr=as.vector(GenomeInfoDb::seqnames(bindingRegions)),
    start=BiocGenerics::start(bindingRegions),
    end=BiocGenerics::end(bindingRegions),
    bindinSites=sitesVector
    )
  utils::write.table(df,
                     file=acceptedRegionsOutputFile,
                     row.names=FALSE,
                     col.names=TRUE,
                     quote=FALSE,
                     sep="\t"
                     )

  return(means)

}



#' @title  Detect binding sites from motif
#' @description DETECT Binding sites with given motif and mismatch number as
#'  well genome/build, False Discovery Rate for a given experiment name.
#' @description This function is called by both
#'  \code{\link{DetectBindingSitesBed}} and
#'   \code{\link{DetectBindingSitesMotif}} with different input.
#' @param From Type of motif dataset either "Motif" or "Bed"
#' @param BedFile Motif locations in bed format file
#' @param motif motif characters in nucleotide IUPAC format
#' @param mismatchNumber Number of mismatches allowed to match with motifs
#' @param chipSeq ChIP-seq alignment both IP and background in 1nt bed format
#'  files
#' @param genome The genome name such as "Hsapiens", "Mmusculus",
#'  "Dmelanogaster"
#' @param genomeBuild The genome build such as "hg38", "hg19", "mm10", "dm3"
#' @param DB The database of genome build. default: "UCSC"
#' @param fdrValue FDR value cut-off
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param GivenRegion granges of user provided binding regions
#' @param currentDir Directory for I/O operations
#' @return A list of FRiPs, sequence statistics, and Motif statistics

DetectBindingSites  <-
  function(From, BedFile, motif, mismatchNumber,chipSeq, genome, genomeBuild,
           DB="UCSC", fdrValue=0.05, windowSize=100, GivenRegion=NA, currentDir)
{

  replicateNumber <-  length(chipSeq$IPfiles)
  if(From=="Motif")
  {
    motifFile <- file.path(currentDir, "Metadata", "Motif_Locations")
    if(is.na(GivenRegion)[1])
    {
      findMotifs(motif=motif,
                 mismatchNumber=mismatchNumber,
                 genome=genome,
                 genomeBuild=genomeBuild,
                 DB=DB,
                 mainCHRs=TRUE,
                 firstCHR=FALSE,
                 MotifLocationName=
                   file.path(currentDir, "Metadata", "Motif_Locations")
                 ) # only main chromosomes otherwise FALSE, FALSE

    } else {
      findMotifs(motif=motif,
                 mismatchNumber=mismatchNumber,
                 genome=genome,
                 genomeBuild=genomeBuild,
                 DB=DB,
                 mainCHRs=TRUE,
                 firstCHR=FALSE,
                 MotifLocationName=
                   file.path(currentDir, "Metadata", "Motif_Locations"),
                 limitedRegion=GivenRegion
                 ) # only main chromosomes otherwise FALSE, FALSE
    }

  } else if(From=="Bed"){
    motifFile <- file.path(currentDir, "Metadata", "Motif_Locations")
    file.copy(from=BedFile, to=motifFile)
  } else
  {
    stop("Motif centeric peak calling can be done on either a string motif or
         bed file")
  }

  motifCount(motifFile=motifFile, chipSeq=chipSeq, windowSize=windowSize,
             outputName=file.path(currentDir, "Metadata", "motifCountTable"),
             currentDir=file.path(currentDir, "Metadata")
             )

  sequencingStatitic <- motifTablePreProcess(
    countTableFile=file.path(currentDir, "Metadata","motifCountTable"),
    outFile=file.path(currentDir, "Metadata","TestResults"),
    currentDir=file.path(currentDir, "Metadata")
    )

  FDRcutoff <- DetectFdrCutoffBH(
    TestTableFile=file.path(currentDir, "Metadata", "TestResults"),
    fdrValue=fdrValue
    )

  averageBindings <- combineTestResults(
    motifFile=motifFile,
    acceptedMotifsOutputFile=
      file.path(currentDir, "Metadata", "BindingMotifsTable"),
    acceptedRegionsOutputFile=
      file.path(currentDir, "Metadata", "BindingRegions"),
    countTableFile=file.path(currentDir, "Metadata", "motifCountTable"),
    testTableFile=file.path(currentDir, "Metadata", "TestResults"),
    fdrCutoff=FDRcutoff,
    windowSize=windowSize
    )

  FRiP <- deriveHeuristicBindingDistribution(
    chipSeq=chipSeq,
    averageBindings=averageBindings,
    windowSize=windowSize,
    acceptedRegionsOutputFile=
      file.path(currentDir, "Metadata", "BindingRegions"),
    currentDir=file.path(currentDir, "Metadata")
    )

  motifStatistics <- decomposeBindingSignal(
    windowSize=windowSize,
    replicateNumber=replicateNumber,
    acceptedRegionsOutputFile=
      file.path(currentDir, "Metadata", "BindingRegions"),
    acceptedMotifsOutputFile=
      file.path(currentDir, "Metadata", "BindingMotifsTable"),
    currentDir=file.path(currentDir, "Metadata")
    )

  DeleteMultipleFiles(
    c(as.character(chipSeq$IPfiles), as.character(chipSeq$BackgroundFiles)))

  file.copy(file.path(currentDir, "Metadata", "BindingRegions"), currentDir)
  if(replicateNumber==1)
  {
    file.copy(file.path(currentDir, "Metadata", "BindingMotifsTable"),
              currentDir)
  }else
  {
    tmpTable <- utils::read.table(
      file.path(currentDir, "Metadata", "BindingMotifsTable"),
      header=TRUE,
      stringsAsFactors=FALSE
      )

    utils::write.table(
      tmpTable[,c(seq_len(5), seq(from=8, to=(dim(tmpTable)[2])))],
      file=file.path(currentDir, "BindingMotifsTable"),
      row.names=FALSE,
      col.names=TRUE,
      quote=FALSE
      )
  }

  return(list(
    FRiP=FRiP,
    sequencingStatitic=sequencingStatitic,
    motifStatistics=motifStatistics))
}


#' @title  Convert bam and bed files to 1 nucleotide bed
#' @description Take alignment files in bam or bed fomat and convert them to 1
#'  nucleotide bed file
#' @param InputFile Original alignment file name
#' @param bedFile Name of output 1nt bed file
#' @param format alignment format and should be one of these: "BAMPE", "BAMSE",
#'  "BEDPE", "BEDSE"
#' @return No return value

generate1ntBedAlignment  <- function(InputFile, bedFile, format=""){

  if(!((format=="BAMPE")||(format=="BAMSE")||(format=="BEDPE")||
       (format=="BEDSE")))
  {
    stop("Input format should be to one of BAMPE, BAMSE, BEDPE, or BEDSE")
  }

  if(format=="BAMSE"){
    gal <- GenomicAlignments::readGAlignments(
      InputFile,
      use.names=FALSE,
      param=NULL,
      with.which_label=FALSE
      )
    gc()
    starts <- BiocGenerics::start(gal)
    ends <- BiocGenerics::end(gal)
    shortReads <- GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(GenomeInfoDb::seqnames(gal)),
      ranges=IRanges::IRanges(round((starts+ends)/2),
                                end=round((starts+ends)/2)
                                )
      )
    rm(starts,ends)
    gc()
    Granges2Bed(shortReads, bedFile)
    rm(shortReads)
    gc()
  } else  if(format=="BAMPE"){
    gal <-
      GenomicAlignments::readGAlignmentPairs(InputFile,
                                             use.names=FALSE,
                                             param=NULL,
                                             with.which_label=FALSE,
                                             strandMode=1
                                             )
    gc()
    starts<-
      pmin(
        pmin(
          BiocGenerics::start(GenomicAlignments::first(gal)),
          BiocGenerics::end(GenomicAlignments::first(gal))
          ),
        pmin(BiocGenerics::start(GenomicAlignments::last(gal)),
             BiocGenerics::end(GenomicAlignments::last(gal))
             )
        )
    gc()
    ends<-
      pmax(
        pmax(
          BiocGenerics::start(GenomicAlignments::first(gal)),
          BiocGenerics::end(GenomicAlignments::first(gal))
          ),
        pmax(BiocGenerics::start(GenomicAlignments::last(gal)),
             BiocGenerics::end(GenomicAlignments::last(gal))
             )
        )
    gc()
    shortReads <- GenomicRanges::GRanges(
      seqnames=
        S4Vectors::Rle(as.character(
          GenomeInfoDb::seqnames(GenomicAlignments::first(gal)))),
      ranges=IRanges::IRanges(
        round((starts+ends)/2), end=round((starts+ends)/2))
      )

    rm(gal, starts, ends)
    gc()
    Granges2Bed(shortReads, bedFile)
    rm(shortReads)
    gc()
  } else  if(format=="BEDSE"){
    Table <- utils::read.table(InputFile, header=FALSE, stringsAsFactors=FALSE)
    if(dim(Table)[2]<3)
    {
      stop("BEDSE table must have at least three columns")
    }
    chrs <- as.character(Table[,1])
    starts<- as.numeric(Table[,2])
    ends<- as.numeric(Table[,3])
    rm(Table)
    gc()
    shortReads <- GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(chrs),
      ranges=
        IRanges::IRanges(round((starts+ends)/2), end=round((starts+ends)/2))
      )

    rm(chrs, starts, ends)
    gc()
    Granges2Bed(shortReads, bedFile)
    rm(shortReads)
    gc()

  } else  if(format=="BEDPE"){
    Table <- utils::read.table(InputFile, header=FALSE, stringsAsFactors=FALSE)
    if(dim(Table)[2]<6)
    {
      stop("BEDPE table must have at least six columns")
    }
    chrs <- as.character(Table[,1])
    starts<- pmin(as.numeric(Table[,2]), as.numeric(Table[,5]))
    ends<- pmax(as.numeric(Table[,3]), as.numeric(Table[,6]))
    rm(Table)
    gc()
    shortReads <- GenomicRanges::GRanges(
      seqnames=S4Vectors::Rle(chrs),
      ranges=
        IRanges::IRanges(round((starts+ends)/2), end=round((starts+ends)/2))
      )

    rm(chrs,starts, ends)
    gc()
    Granges2Bed(shortReads, bedFile)
    rm(shortReads)
    gc()

  }

}




#' @title  Detect binding sites from bed motif input
#' @description  Takes user provied bed regions, and check for validity of them.
#'  Read bam or bed alignment files and convert to 1 nt bed and call detect
#'   binding site from 1nt bed.
#' @param BedFile Motif locations in bed format file
#' @param IPfiles IP ChIP-seq alignment files
#' @param BackgroundFiles Background ChIP-seq alignment files. Can be Input
#'  experimetn, DNA whole exctract, etc.
#' @param genome The genome name such as "Hsapiens", "Mmusculus",
#'  "Dmelanogaster"
#' @param genomeBuild The genome build such as "hg38", "hg19", "mm10", "dm3"
#' @param DB The database of genome build. default: "UCSC"
#' @param fdrValue FDR value cut-off
#' @param expName The name of the output table
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param format alignment format and should be one of these: "BAMPE", "BAMSE",
#'  "BEDPE", "BEDSE"
#' @return peakCallingStatistics A list FRiPs, sequence statistics, and Motif
#'  statistics
#'
#' @examples
#'
#' # FUR candidate motifs in NC_000913 E. coli
#' FurMotifs=system.file("extdata", "FurMotifs.bed", package="Motif2Site")
#'
#' # ChIP-seq datasets in bed single end format
#' IPFe <- c(system.file("extdata", "FUR_fe1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_fe2.bed", package="Motif2Site"))
#' Inputs <- c(system.file("extdata", "Input1.bed", package="Motif2Site"),
#'             system.file("extdata", "Input2.bed", package="Motif2Site"))
#' FURfeBedInputStats <- 
#'   DetectBindingSitesBed(BedFile=FurMotifs,
#'                         IPfiles=IPFe, 
#'                         BackgroundFiles=Inputs, 
#'                         genome="Ecoli",
#'                         genomeBuild="20080805",
#'                         DB="NCBI",
#'                         expName="FUR_Fe_BedInput",
#'                         format="BEDSE"
#'                        )
#'
#'
#' @seealso
#' \code{\link{DetectBindingSitesMotif}}
#'
#' @export

DetectBindingSitesBed <-
  function(BedFile, IPfiles, BackgroundFiles, genome, genomeBuild, DB="UCSC",
           fdrValue=0.05, expName="Motif_Centric_Peaks",windowSize=100,
           format="")
{

  if(!(file.exists(BedFile)))
  {
    stop(BedFile, " does not exist")
  }
  BedFile <- normalizePath(BedFile)


  replicateNumber <- length(IPfiles)

  if(length(IPfiles)!=length(BackgroundFiles))
  {
    stop("Length of IPfiles and BackgroundFiles must be identical")
  }

  for(i in seq_len(replicateNumber))
  {
    if(!(file.exists(IPfiles[i])))
    {
      stop(IPfiles[i], " does not exist")
    }
    if(!(file.exists(BackgroundFiles[i])))
    {
      stop(BackgroundFiles[i], " does not exist")
    }

    IPfiles[i] <- normalizePath(IPfiles[i])
    BackgroundFiles[i] <- normalizePath(BackgroundFiles[i])

  }

  chipSeq <- data.frame(IPfiles=IPfiles, BackgroundFiles=BackgroundFiles)

  if(!((format=="BAMPE")||(format=="BAMSE")||(format=="BEDPE")
       ||(format=="BEDSE")))
  {
    stop("Input format should be to one of BAMPE, BAMSE, BEDPE, or BEDSE")
  }

  BSGstring <- paste("BSgenome.", genome,".", DB, ".", genomeBuild,sep="")
  if (!requireNamespace(BSGstring, quietly=TRUE))
  {
    stop(BSGstring, " has not been installed")
  }

  currentDir <- getwd()

  if (dir.exists(file.path(currentDir, expName))){
    unlink(x=(file.path(currentDir, expName)), recursive=TRUE)
  }

  dir.create(file.path(currentDir, expName))
  
  if (!file.exists(file.path(currentDir, expName, "Metadata"))){
    dir.create(file.path(currentDir, expName, "Metadata"))
  }


  for(i in seq_len(replicateNumber))
  {
    generate1ntBedAlignment(
      InputFile=as.character(chipSeq$IPfiles[i]),
      bedFile=file.path(currentDir, expName, "Metadata",paste0("IP_",i,".bed")),
      format=format
      )

    generate1ntBedAlignment(
      InputFile=as.character(chipSeq$BackgroundFiles[i]),
      bedFile=file.path(currentDir, expName, "Metadata",
                        paste0("Background_",i,".bed")
                        ),
      format=format
      )
  }
  gc()


  chipSeq <- data.frame(
    IPfiles=file.path(currentDir, expName, "Metadata",
                      paste0("IP_",seq_len(replicateNumber),".bed")
                      ),
    BackgroundFiles=
      file.path(currentDir, expName, "Metadata",
                paste0("Background_",seq_len(replicateNumber),".bed")
    )
  )


  peakCallingStatistics <- DetectBindingSites(
    From="Bed", BedFile=BedFile,  chipSeq=chipSeq, genome=genome,
    genomeBuild=genomeBuild, DB=DB, fdrValue=fdrValue, windowSize=windowSize,
    currentDir=file.path(currentDir, expName)
    )

  # Write ChIP-seq statistics in a file
  FRiPs <- ""
  for(i in seq_len(replicateNumber))
  {
    FRiPs <- paste0(FRiPs, " ",peakCallingStatistics$FRiP[i])
  }

  dfStats <- data.frame(
    replicateNumber=replicateNumber,
    NonSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$nonSequenced,
    underSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$underBinding,
    overSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$overBinding,
    skewedMotifsRejectionNumber=
      peakCallingStatistics$motifStatistics$skewnessTestRejected,
    decomposedMotifsRejectionNumber=
      peakCallingStatistics$motifStatistics$decompositionRejected,
    acceptedMotifsNumber=peakCallingStatistics$motifStatistics$accepted,
    FRiPs=FRiPs
  )
  utils::write.table(
    dfStats,
    file=file.path(currentDir, expName, "Statistics.tsv"),
    row.names=FALSE,
    col.names=TRUE,
    quote=FALSE,
    sep="\t"
    )

  return(peakCallingStatistics)

}






#' @title  Detect binding sites from sequence motif sequence and mismatchNumber
#' @description DETECT Binding sites with given motif and mismatch number as
#' well genome/build, False Discovery Rate for a given experiment name. Read bam
#'  or bed alignment files and convert to 1 nt bed and detect binding site among
#'   motifs from 1nt bed alignment.
#' @param motif motif characters in nucleotide IUPAC format
#' @param mismatchNumber Number of mismatches allowed to match with motifs
#' @param IPfiles IP ChIP-seq alignment files
#' @param BackgroundFiles Background ChIP-seq alignment files. Can be Input
#'  experimetn, DNA whole exctract, etc.
#' @param genome The genome name such as "Hsapiens", "Mmusculus",
#'  "Dmelanogaster"
#' @param genomeBuild The genome build such as "hg38", "hg19", "mm10", "dm3"
#' @param DB The database of genome build. default: "UCSC"
#' @param fdrValue FDR value cut-off
#' @param windowSize Window size around binding site. The total region would be
#'  2*windowSize+1
#' @param expName The name of the output table
#' @param format alignment format and should be one of these: "BAMPE", "BAMSE",
#'  "BEDPE", "BEDSE"
#' @param GivenRegion granges of user provided binding regions
#' @return A list FRiPs, sequence statistics, and Motif statistics
#'
#' @examples
#'
#' # ChIP-seq datasets in bed single end format
#' IPFe <- c(system.file("extdata", "FUR_fe1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_fe2.bed", package="Motif2Site"))
#' Inputs <- c(system.file("extdata", "Input1.bed", package="Motif2Site"),
#'             system.file("extdata", "Input2.bed", package="Motif2Site"))
#'
#'  # Granages region for motif search           
#'    NC_000913_Coordiante <-
#'      GenomicRanges::GRanges(seqnames=S4Vectors::Rle("NC_000913"),
#'                             ranges=IRanges::IRanges(1, 4639675))           
#'             
#' FURfeStringInputStats <- 
#'   DetectBindingSitesMotif(motif="GWWTGAGAA",
#'    mismatchNumber=1,
#'    IPfiles=IPFe, 
#'    BackgroundFiles=Inputs, 
#'    genome="Ecoli",
#'    genomeBuild="20080805",
#'    DB="NCBI",
#'    expName="FUR_Fe_StringInput",
#'    format="BEDSE",
#'    GivenRegion=NC_000913_Coordiante 
#'    )
#'
#'                                       
#' @seealso
#' \code{\link{DetectBindingSitesBed}}
#'
#' @export

DetectBindingSitesMotif <-
  function(motif, mismatchNumber, IPfiles, BackgroundFiles, genome, genomeBuild,
           DB="UCSC", fdrValue=0.05, expName="Motif_Centric_Peaks",
           windowSize=100, format="",GivenRegion=NA)
{

  replicateNumber <- length(IPfiles)

  if(length(IPfiles)!=length(BackgroundFiles))
  {
    stop("Length of IPfiles and BackgroundFiles must be identical")
  }

  for(i in seq_len(replicateNumber))
  {
    if(!(file.exists(IPfiles[i])))
    {
      stop(IPfiles[i], " does not exist")
    }
    if(!(file.exists(BackgroundFiles[i])))
    {
      stop(BackgroundFiles[i], " does not exist")
    }

    IPfiles[i] <- normalizePath(IPfiles[i])
    BackgroundFiles[i] <- normalizePath(BackgroundFiles[i])

  }

  chipSeq <- data.frame( IPfiles=IPfiles, BackgroundFiles=BackgroundFiles)

  if(!is.na(GivenRegion)[1])
  {
    if(!methods::is(GivenRegion,"GRanges")){
      stop("The user should provided Regions as GenomicRanges")
    }
  }
  if(!((format=="BAMPE")||(format=="BAMSE")||(format=="BEDPE")||
       (format=="BEDSE")))
  {
    stop("Input format should be to one of BAMPE, BAMSE, BEDPE, or BEDSE ")
  }

  BSGstring <- paste("BSgenome.", genome,".", DB, ".", genomeBuild,sep="")
  if (!requireNamespace(BSGstring, quietly=TRUE))
  {
    stop(BSGstring, " has not been installed")
  }


  currentDir <- getwd()

  if (dir.exists(file.path(currentDir, expName))){
    unlink(x=(file.path(currentDir, expName)), recursive=TRUE)
  }
  
  dir.create(file.path(currentDir, expName))
  
  if (!file.exists(file.path(currentDir, expName, "Metadata"))){
    dir.create(file.path(currentDir, expName, "Metadata" ))
  }


  for(i in seq_len(replicateNumber))
  {
    generate1ntBedAlignment(
      InputFile=as.character(chipSeq$IPfiles[i]),
      bedFile=file.path(currentDir, expName, "Metadata",paste0("IP_",i,".bed")),
      format=format)
    generate1ntBedAlignment(
      InputFile=as.character(chipSeq$BackgroundFiles[i]),
      bedFile=file.path(currentDir, expName, "Metadata",
                                  paste0("Background_",i,".bed")
      ),
      format=format
      )
  }

  gc()
  chipSeq <- data.frame(
    IPfiles=file.path(currentDir, expName, "Metadata",
                      paste0("IP_", seq_len(replicateNumber), ".bed")
    ),
    BackgroundFiles=
      file.path(currentDir, expName, "Metadata",
                paste0("Background_", seq_len(replicateNumber), ".bed")
    )
  )
  if(!is.na(GivenRegion)[1])
  {
    peakCallingStatistics <- DetectBindingSites(
      From="Motif", motif=motif, mismatchNumber=mismatchNumber,chipSeq=chipSeq,
      genome=genome, genomeBuild=genomeBuild, DB=DB, fdrValue=fdrValue,
      windowSize=windowSize, GivenRegion=GivenRegion,
      currentDir=file.path(currentDir, expName)
      )
  } else {
    peakCallingStatistics <- DetectBindingSites(
      From="Motif", motif=motif, mismatchNumber=mismatchNumber,chipSeq=chipSeq,
      genome=genome, genomeBuild=genomeBuild, DB=DB, fdrValue=fdrValue,
      windowSize=windowSize, currentDir=file.path(currentDir, expName)

      )
  }

  # Write ChIP-seq statistics in a file
  FRiPs <- ""
  for(i in seq_len(replicateNumber))
  {
    FRiPs <- paste0(FRiPs, " ",peakCallingStatistics$FRiP[i])
  }

  dfStats <- data.frame(
    replicateNumber=replicateNumber,
    NonSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$nonSequenced,
    underSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$underBinding,
    overSequencedMotifsRatio=
      peakCallingStatistics$sequencingStatitic$overBinding,
    skewedMotifsRejectionNumber=
      peakCallingStatistics$motifStatistics$skewnessTestRejected,
    decomposedMotifsRejectionNumber=
      peakCallingStatistics$motifStatistics$decompositionRejected,
    acceptedMotifsNumber=peakCallingStatistics$motifStatistics$accepted,
    FRiPs=FRiPs
  )
  utils::write.table(
    dfStats,
    file=file.path(currentDir, expName, "Statistics.tsv"),
    row.names=FALSE,
    col.names=TRUE,
    quote=FALSE,
    sep="\t"
    )

  return(peakCallingStatistics)

}


#' @title  Detect differential motifs
#' @description Take combined matrix of motif counts generated by
#'  \code{\link{recenterBindingSitesAcrossExperiments}}, and experiment names.
#'  It detect differential motifs using edgeR TMM nomralizaiton with Generalized
#'   linear model
#' @param tableOfCountsDir Directory which conatins the combined motifs and
#'  ChIP-seq count file
#' @param exp1 Experiment name which will be compared in pairwise comparison
#' @param exp2 Experiment name which will be compared in pairwise comparison
#' @param FDRcutoff FDR cutoff applies on pvalue distribution
#' @param logFCcuttoff log fold change cutoff
#' @return A list of differential motifs, motif1 and motif2 as well as a table
#'  of total motifs and log fold changes
#'
#' @examples
#'
#' # FUR candidate motifs in NC_000913 E. coli
#' FurMotifs=system.file("extdata", "FurMotifs.bed", package="Motif2Site")
#'
#' # ChIP-seq datasets fe in bed single end format
#' IPFe <- c(system.file("extdata", "FUR_fe1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_fe2.bed", package="Motif2Site"))
#' Inputs <- c(system.file("extdata", "Input1.bed", package="Motif2Site"),
#'             system.file("extdata", "Input2.bed", package="Motif2Site"))
#' FURfeBedInputStats <- 
#'   DetectBindingSitesBed(BedFile=FurMotifs,
#'    IPfiles=IPFe, 
#'    BackgroundFiles=Inputs, 
#'    genome="Ecoli",
#'    genomeBuild="20080805",
#'    DB="NCBI",
#'    expName="FUR_Fe_BedInput",
#'    format="BEDSE"
#'    )
#'
#' # ChIP-seq datasets dpd in bed single end format
#' IPDpd <- c(system.file("extdata", "FUR_dpd1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_dpd2.bed", package="Motif2Site"))
#' FURdpdBedInputStats <- 
#'   DetectBindingSitesBed(BedFile=FurMotifs,
#'    IPfiles=IPDpd, 
#'    BackgroundFiles=Inputs, 
#'    genome="Ecoli",
#'    genomeBuild="20080805",
#'    DB="NCBI",
#'    expName="FUR_Dpd_BedInput",
#'    format="BEDSE"
#'    )
#'                        
#'
#' # Combine all FUR binding sites into one table
#' corMAT <- recenterBindingSitesAcrossExperiments(
#'   expLocations=c("FUR_Fe_BedInput","FUR_Dpd_BedInput"),
#'   experimentNames=c("FUR_Fe","FUR_Dpd"),
#'   expName="combinedFUR",
#'   )
#'
#' # Differential binding sites across FUR conditions fe vs dpd
#' diffFUR <- pairwisDifferential(tableOfCountsDir="combinedFUR",
#'    exp1="FUR_Fe",
#'    exp2="FUR_Dpd",
#'    FDRcutoff=0.05,
#'    logFCcuttoff=1
#'    )
#' 
#' FeUp <- diffFUR[[1]]
#' DpdUp <- diffFUR[[2]]
#' TotalComparison <- diffFUR[[3]]
#' head(TotalComparison)
#'
#' @seealso
#' \code{\link{recenterBindingSitesAcrossExperiments}}
#'
#' @export




pairwisDifferential  <- function(tableOfCountsDir="", exp1, exp2,
                                 FDRcutoff=0.05, logFCcuttoff=1)
{

  tableOfCountsFile <- file.path(tableOfCountsDir,"CombinedMatrix")
  if (!(file.exists(tableOfCountsFile))){
    stop("table of counts file does not exist")
  }

  tableOfCounts <- utils::read.table(tableOfCountsFile,
                                     header=TRUE,
                                     check.names=FALSE,
                                     stringsAsFactors=FALSE
                                     )

  if(length(match(exp1, colnames(tableOfCounts)))<0){
    stop("experiment 1 does not exists in the combined Matrix")
  }

  if(length(match(exp2, colnames(tableOfCounts)))<0){
    stop("experiment 2 does not exists in the combined Matrix")
  }

  Motifs <- GenomicRanges::GRanges(
    seqnames=S4Vectors::Rle(tableOfCounts[,1]),
    ranges=IRanges::IRanges(tableOfCounts[,2], tableOfCounts[,3]))

  exp1Inds <- which((colnames(tableOfCounts)==exp1))
  exp2Inds <- which((colnames(tableOfCounts)==exp2))
  combinedInds <- c(exp1Inds, exp2Inds)
  tableOfCounts <- tableOfCounts[,combinedInds]

  # Remove lines with NA
  notNA <- which(!is.na(rowSums(tableOfCounts)))
  Motifs <-  Motifs[notNA]
  tableOfCounts <- tableOfCounts[notNA,]
  rownames(tableOfCounts) <- as.character(seq_len(dim(tableOfCounts)[1]))

  tableOfCounts <- edgeR::DGEList(
    counts=tableOfCounts,
    group=c(rep(exp1,length(exp1Inds)),rep(exp2,length(exp2Inds))))

  tableOfCounts <- edgeR::calcNormFactors(tableOfCounts, method="TMM")

  # Desing matrix
  Group <- stats::relevel(factor(tableOfCounts$samples$group),
                          ref=as.character(tableOfCounts$samples$group[1]))
  design <- stats::model.matrix(~Group)
  colnames(design) <- levels(Group)

  # GLM test
  tableOfCounts <- edgeR::estimateDisp(tableOfCounts, design=design)
  gc()
  fit <- edgeR::glmFit(tableOfCounts, design=design)
  dgeLRTtest <- edgeR::glmLRT(fit, coef=2)

  resLRT <- edgeR::topTags(dgeLRTtest, n=nrow(tableOfCounts$counts))
  selectedLRT <-
    resLRT$table$FDR < FDRcutoff & abs(resLRT$table$logFC) > logFCcuttoff
  selectedLRT <- resLRT$table[selectedLRT, ]
  selectedLRT$updown <- factor(ifelse(selectedLRT$logFC > 0, "up", "down"))

  motifs2 <- Motifs[as.numeric(row.names(selectedLRT)
                               [which((selectedLRT$updown=="up")==TRUE)])]

  motifs1 <- Motifs[as.numeric(row.names(selectedLRT)
                               [which((selectedLRT$updown=="down")==TRUE)])]

  totalResults <- 
    as.data.frame(c(as.data.frame(Motifs)[, seq_len(3)], dgeLRTtest$table))

  diffs <- list(Motif1Up=motifs1, Motif2Up=motifs2, TotalResults=totalResults)
  return(diffs)
}








#' @title  Combine binding sites across experiments
#' @description Take experiment folder locations and experiment names and
#'  combine them into a combined matrix of motifs and ChIP-seq counts
#' @description Experiment folders must be generated either by
#' \code{\link{DetectBindingSitesBed}} or \code{\link{DetectBindingSitesMotif}}.
#' @param expLocations The path to the experiment folders
#' @param experimentNames Name of the experiment to be used in combined ChIP-seq
#' @param expName Name of the combined matrix
#' @param fdrValue FDR cut-off to accept binding in each ChIP-seq experiments
#' @param fdrCrossExp If no experiment fullfill this cutoff, the motif is not
#'  considered
#' @return A pariwise Pearson correlation matrix across experiments
#'
#' @examples
#'
#' # FUR candidate motifs in NC_000913 E. coli
#' FurMotifs=system.file("extdata", "FurMotifs.bed", package="Motif2Site")
#'
#' # ChIP-seq datasets fe in bed single end format
#' IPFe <- c(system.file("extdata", "FUR_fe1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_fe2.bed", package="Motif2Site"))
#' Inputs <- c(system.file("extdata", "Input1.bed", package="Motif2Site"),
#'             system.file("extdata", "Input2.bed", package="Motif2Site"))
#' FURfeBedInputStats <- 
#'   DetectBindingSitesBed(BedFile=FurMotifs,
#'    IPfiles=IPFe, 
#'    BackgroundFiles=Inputs, 
#'    genome="Ecoli",
#'    genomeBuild="20080805",
#'    DB="NCBI",
#'    expName="FUR_Fe_BedInput",
#'    format="BEDSE"
#'    )
#'
#' # ChIP-seq datasets dpd in bed single end format
#' IPDpd <- c(system.file("extdata", "FUR_dpd1.bed", package="Motif2Site"),
#'         system.file("extdata", "FUR_dpd2.bed", package="Motif2Site"))
#' FURdpdBedInputStats <- 
#'   DetectBindingSitesBed(BedFile=FurMotifs,
#'    IPfiles=IPDpd, 
#'    BackgroundFiles=Inputs, 
#'    genome="Ecoli",
#'    genomeBuild="20080805",
#'    DB="NCBI",
#'    expName="FUR_Dpd_BedInput",
#'    format="BEDSE"
#'    )
#'                        
#'
#' # Combine all FUR binding sites into one table
#' corMAT <- recenterBindingSitesAcrossExperiments(
#'     expLocations=c("FUR_Fe_BedInput","FUR_Dpd_BedInput"),
#'     experimentNames=c("FUR_Fe","FUR_Dpd"),
#'     expName="combinedFUR",
#'     )
#' corMAT
#'
#' @seealso
#' \code{\link{pairwisDifferential}}
#'
#' @export


recenterBindingSitesAcrossExperiments <-
  function(expLocations, experimentNames, expName="combinedData", fdrValue=0.05,
           fdrCrossExp=0.001)
{

  currentDir <- getwd()
  datasetNumber <- length(expLocations)

  if(datasetNumber!=length(experimentNames))
  {
    stop("Length of expLocations, experimentNames,
         and replicateNumbers vectors must be identical")
  }


  if( datasetNumber<2)
  {
    stop("At least two different experiments are needed to be combined")
  }

  ## read replicate numbers from the folers' statistic files
  replicateNumbers <- rep(0, datasetNumber)


  for(i in seq_len(datasetNumber))
  {
    if (file.exists(expLocations[i]))
    {
      expLocations[i] <- normalizePath(expLocations[i])
      if (!(file.exists(file.path(as.character(expLocations[i]),
                               "Statistics.tsv"))))
        {

        stop("In directory ", as.character(expLocations[i]),
                    " Statistics.tsv does not exists")

      }else
      {
        StatTable <- file.path(as.character(expLocations[i]), "Statistics.tsv")
        replicateNumbers[i] <- (
          utils::read.table(StatTable,
                            header=TRUE,
                            sep="\t",
                            stringsAsFactors=FALSE)
          )$replicateNumber
      }
    }else
    {
      stop("Direcotry ", as.character(expLocations[i]), " does not exists")
    }

  }

  df_datasets <- data.frame(
    locations=expLocations,
    experiments=experimentNames,
    replicates=replicateNumbers
    )

  ## read motifs and count values
  if (file.exists(as.character(df_datasets$locations[1]))){
    if (!(file.exists(file.path(as.character(df_datasets$locations[1]),
                             "Metadata"))))
    {
      stop("Metadata directory does not exist in ",
           df_datasets$locations[1], " folder")
    }

    if (!(file.exists(file.path(as.character(df_datasets$locations[1]),
                             "Metadata", "BindingMotifsTable"))))
    {
      stop("Binding motif table does not exists in the first experiment
           Metadata folder")
    } else {

      Table <- utils::read.table(
        file=file.path(as.character(df_datasets$locations[1]),
                    "Metadata", "BindingMotifsTable"), header=TRUE,
        stringsAsFactors=FALSE)

      Motifs <- GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(Table[,1]),
        ranges=IRanges::IRanges(Table[,2], Table[,3]))

    }

  } else{
    stop("Direcotry ", df_datasets$locations[1], " does not exists")
  }

  for(i in c(2:datasetNumber))
  {
    if (file.exists(as.character(df_datasets$locations[i]))){
      if (!(file.exists(file.path(as.character(df_datasets$locations[i]),
                               "Metadata"))))
      {
        stop("Metadata directory does not exist in ",
             df_datasets$locations[1], " folder")
      }

      if (!(file.exists(file.path(as.character(df_datasets$locations[i]),
                               "Metadata", "BindingMotifsTable"))))
        {
        stop("Binding motif table does not exists in the ", i,
             "th experiment Metadata folder")
      } else {
        Table <- utils::read.table(
          file=file.path(as.character(df_datasets$locations[i]),
                      "Metadata", "BindingMotifsTable"),
          header=TRUE, stringsAsFactors=FALSE)
        tmpGranges <- GenomicRanges::GRanges(
          seqnames=S4Vectors::Rle(Table[,1]),
          ranges=IRanges::IRanges(Table[,2], Table[,3]))
        Motifs <- GenomicRanges::reduce(c(Motifs,tmpGranges))

      }

    } else{
      stop("Direcotry ",df_datasets$locations[i], " does not exists")
    }
  }

  rm(Table, tmpGranges)
  gc()

  motifNumber <- length(Motifs)

  acceptedMotifs <- rep(FALSE, motifNumber)

  countMatrix <- matrix(NA, nrow=motifNumber,
                        ncol=sum(df_datasets$replicates))
  colnames(countMatrix) <- rep(df_datasets$experiments, df_datasets$replicates)

  pvalMatrix <- matrix(NA, nrow=motifNumber, ncol=datasetNumber)
  colnames(pvalMatrix) <- df_datasets$experiments

  for(i in seq_len(datasetNumber))
  {
    if (!(file.exists(file.path(as.character(df_datasets$locations[i]),
                             "Metadata", "Motif_Locations"))))
     {
      stop("Motif_Locations file does not exists in the ", i,
           "th experiment Metadata folder")
    }
    if (!(file.exists(file.path(as.character(df_datasets$locations[i]),
                             "Metadata", "TestResults"))))
    {
      stop("TestResults file does not exists in the ", i,
           "th experiment Metadata folder")
    }

    if (!(file.exists(file.path(as.character(df_datasets$locations[i]),
                             "Metadata", "motifCountTable")))){
      stop("motifCountTable file does not exists in the ", i,
           "th experiment Metadata folder")
    }


    motifLocations <- 
      Bed2Granges(file.path(as.character(df_datasets$locations[i]),
                            "Metadata",
                            "Motif_Locations"
                            )
                  )
    tmpInds <- GenomicRanges::findOverlaps(motifLocations,Motifs)
    rm(motifLocations)
    gc()

    testTable <- utils::read.table(
      file=file.path(as.character(df_datasets$locations[i]),
                     "Metadata",
                     "TestResults"
                     ),
      header=TRUE,
      stringsAsFactors=FALSE
      )
    pvalMatrix[S4Vectors::subjectHits(tmpInds),i] <-
      testTable$Pvalue[S4Vectors::queryHits(tmpInds)]
    rm(testTable)
    gc()

    base <- 0
    if(i>1){ base <- sum(df_datasets$replicates[seq_len(i-1)])  }


    motifCountTable <- utils::read.table(
      file=file.path(as.character(df_datasets$locations[i]),
                  "Metadata", "motifCountTable"), header=TRUE,
      stringsAsFactors=FALSE)
    countMatrix[S4Vectors::subjectHits(tmpInds),
                (base+1):(base+df_datasets$replicates[i])]<-
      as.matrix(motifCountTable[S4Vectors::queryHits(tmpInds),
                                seq_len(df_datasets$replicates[i])])
    rm(motifCountTable)
    gc()

  }

  bindingMatrix <- matrix("nonBinding", nrow=motifNumber,
                          ncol=datasetNumber)
  for(i in seq_len(datasetNumber))
  {
    bindingMatrix[which(is.na(pvalMatrix[,i])),i] <- "NotTested"
    pvalsAdjusted <- stats::p.adjust(pvalMatrix[,i], method="BH")
    bindingMatrix[which((pvalsAdjusted<fdrValue)==TRUE),i] <- "Binding"
    acceptedMotifs[which((pvalsAdjusted<fdrCrossExp)==TRUE)] <- TRUE
  }
  colnames(bindingMatrix) <- paste0(df_datasets$experiments,"_binding")

  currentDir <- getwd()

  if (dir.exists(file.path(currentDir, expName))){
    unlink(x=(file.path(currentDir, expName)), recursive=TRUE)
  }
  
  dir.create(file.path(currentDir, expName))
  
  dfMotif<- data.frame(
    chr=as.vector(GenomeInfoDb::seqnames(Motifs)),
    start=BiocGenerics::start(Motifs),
    end=BiocGenerics::end(Motifs)

  )

  df <- as.data.frame(c(dfMotif,as.data.frame(bindingMatrix),
                        as.data.frame(countMatrix)),check.names=FALSE)
  df <- df[acceptedMotifs,]
  utils::write.table( df,
                      file=file.path(currentDir, expName, "CombinedMatrix"),
                      quote=FALSE,
                      row.names=FALSE,
                      col.names=TRUE,
                      sep="\t"
                      )

  countMatrix <- countMatrix[acceptedMotifs,]
  corMat <- matrix(1, ncol=dim(countMatrix)[2], nrow=dim(countMatrix)[2])
  colnames(corMat) <- colnames(countMatrix)
  rownames(corMat) <- colnames(countMatrix)

  for(i in seq_len(dim(countMatrix)[2]-1)){
    for(j in seq(from=(i+1), to=(dim(countMatrix)[2]))){
      corMat[i,j] <- 
        stats::cor(log2(countMatrix[,i]+1), log2(countMatrix[,j]+1),
                                use="pairwise.complete.obs"
                   )
      corMat[j,i] <- corMat[i,j]
    }
  }

  return(corMat)
}


#' @title Suppress messages generated by in external package
#' @description mixtools and MASS::fitdistr generates warning by cat which is
#' suppressed by this funcitons
#' @param func functional input call for which cat messages should be supressed
#' @return No return value

quiet <- function(func) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(func))
}

