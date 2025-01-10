if (!"optparse" %in% installed.packages()){
  stop('There is no package called "optparse"', call.=FALSE)
}
library("optparse")
option_list = list(
  make_option(c("-f", "--fastq"), type = "character", default = NULL,
              help = "The path of fastq_pass floder, which contains one or more [barcodeXX]."),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "The path of output floder"),
  make_option("--match", type = "logical", default = TRUE,
              help = "The logical variable of the match operation, the default is TRUE"),
  make_option("--stat", type = "logical", default = TRUE,
              help = "The logical variable of the statistics operation, the default is TRUE"),
  make_option("--index", type = "character", default = NULL,
              help = "The index ID,split by comma, for example: R1,R2"),
  make_option("--min_mis", type = "integer", default = 0,
              help = "The minimum number of mismatched bases allowed, the default is 0"),
  make_option("--max_mis", type = "integer", default = 0,
              help = "The maximum number of mismatched bases allowed, the default is 0"),
  make_option(c("-t", "--cores"), type = "integer", default = 1,
              help = "Number of CPUs, the default is 1")
  
);

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pheatmap))

##### input ######
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
fq_path <- "C:/Users/Zhangyj/Desktop/RData/20240509"
output <- "C:/Users/Zhangyj/Desktop/RData/20240509"
index <- "R5"

NT <- 30
min_mis <- 0
max_mis <- 0
match <- T
statistics <- T

## test
# fq_path <- "/mnt/raid5/Personal/minion/cooperation/Geng_sMdetec/data/nanopore/"
# output <- "/mnt/raid5/Personal/minion/cooperation/Geng_sMdetec/output/nanopore"
# index<-"R1,R2,R3"
# NT <- 30
# min_mis <- 0
# max_mis <- 0
# match <- T
# statistics <- T


#### output #####
input_id_patta <-
  data.table(data.frame(
    index = c("S1", "R1", "R2", "R3", "R4", "R5", "R6", "R7","R8","R9","R10","R11"),
    seqs = c(
      "TACCAGGTCCTA",
      "TAGTCATCTCTA",
      "TATCCATCCTTA",
      "TAAGCTCGCATA",
      "ATCCTCTCCTCA",
      "CGTCTACGATGC",
      "CACGAAGTGGAA",
      "GTTCCCTGTCCC",
      "ACGTTAAGGCCA",
      "TCCTCGTGAGGT",
      "CATCTAACCTAG",
      "TCGGAATTGGCT"
    )
  ))

index_id <- strsplit(index,",")[[1]]
patta_ori <-input_id_patta[index %in% c("S1",index_id)]
print(paste0("^_^ Detecting ",paste0(patta_ori$index,collapse = " "),"......................................."))

print("^_^ Reading fastq(s) now.......................................")

fq_dirs <- list.dirs(fq_path, full.names = TRUE, recursive = FALSE)
fq_files_direct <- list.files(fq_path, pattern = "barcode[0-9]+\\.fastq$", full.names = TRUE)

all_fq_files <- list()

for (file in fq_files_direct) {
  all_fq_files[[gsub(".fastq","",basename(file))]] <- file
}

for (dir in fq_dirs) {
  fq_files_in_dir <- list.files(dir, pattern = ".*\\.fastq$", full.names = TRUE)
  

  if (length(fq_files_in_dir) > 0) {
 
    all_fq_files[[basename(dir)]] <- fq_files_in_dir
  }
}


fq_raw <- list()
for (barcode in names(all_fq_files)) {
  fq_list <- lapply(all_fq_files[[barcode]], function(x) readFastq(x))
  fq_raw[[barcode]] <- do.call(c, fq_list)
}


fq_tab <- rbindlist(lapply(names(fq_raw), function(barcode) {
  short_reads <- fq_raw[[barcode]]
  if (length(short_reads) == 1) {
    dss <- sread(short_reads[[1]])
    return(data.frame(barcode = barcode, width = width(dss), seq = as.character(dss)))
  } else {

    return(do.call(rbind, lapply(short_reads, function(sr) {
      dss <- sread(sr)
      data.frame(barcode = barcode, width = width(dss), seq = as.character(dss))
    })))
  }
}))

rm(fq_raw)
gc()

##### match ######
if (match==F) {
  print("^_^ Writing the raw sequence now................................")
  res <- split(fq_tab, fq_tab$barcode)
  lapply(names(res), function(x){
    write.csv(res[[x]],paste0(output,"/",x,"_fq_tab.csv"),row.names = F)
  })
  
}else{
  print(paste0("^_^ Matching the ",index," and S1 sequences now............................"))
  match_l <- lapply(patta_ori$seqs, function(patta){
    count_p <- parallel::mclapply(fq_tab$seq, function(seq) {
      Biostrings::vcountPattern(
        pattern = patta,
        subject = seq,
        min.mismatch = min_mis,
        max.mismatch = max_mis,
        with.indels = F
      )
    },mc.cores=1)
    c_p <- do.call(c,count_p)
    return(c_p)
  })
  match_m <- data.frame(do.call(cbind,match_l))
  colnames(match_m) <- patta_ori$index
  result <- cbind(match_m,fq_tab)
  
  print("^_^ Writing the raw sequence now...................................")
  res <- split(result,result$barcode)
  lapply(names(res), function(x){
    write.csv(res[[x]],paste0(output,"/",x,"_fq_tab.csv"),row.names = F)
  })
  
  if (statistics==T) {
    print("^_^ Making match statistics now...................................")
    match_df <- data.table(match_m)
    colnames(match_df) <- patta_ori$index
    match_df$barcode <- fq_tab$barcode
    fun <- function (x) {
      sum(x>0,na.rm = T)
    }
    summary <- match_df[, lapply(.SD[, ], fun), by = barcode]
    summary_counts <- summary
    summary_counts$Sum <- match_df[,.N,by=barcode]$N
    write.csv(summary_counts,paste0(output,"/","summary_counts.csv"),row.names = F)
    calculate_percentage <- function(df) {
      if (!"S1" %in% names(df)) {
        stop("S1 column not found")
      }
      cols_to_process <- setdiff(names(df), c("barcode","S1"))
      for (col in cols_to_process) {
        df[[col]] <- (df[[col]] / df$S1) * 100
        names(df)[names(df) == col] <- paste0(col, "/S1%")
      }
      return(df)
    }
    summary_per <- calculate_percentage(df = summary)
    write.csv(summary_per,paste0(output,"/","summary_percent.csv"),row.names = F)
    #summary_per <- read.csv(paste0(output,"/","summary_percent.csv"))
    summary_per <- data.frame(summary_per)
    rownames(summary_per) <- summary_per$barcode
    summary_per$barcode <- NULL
    summary_per$S1 <- NULL
    pdf(paste0(output,"/","summary_heatmap.pdf"),width = 8,height = 12)
    pheatmap(
      summary_per,
      cluster_cols = F,
      cluster_rows = F,
      scale = "row",
      border_color = "white",
      main = "Detected sacle row percentage",
      cellwidth = 40,
      cellheight = 30,
      display_numbers = TRUE
    )
    pheatmap(
      summary_per,
      cluster_cols = F,
      cluster_rows = F,
      border_color = "white",
      main = "Detected percentage",
      cellwidth = 40,
      cellheight = 30,
      display_numbers = TRUE
    )
    dev.off()
  }
}

