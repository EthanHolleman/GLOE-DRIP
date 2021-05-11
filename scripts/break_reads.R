
read_bed_file <- function(bed.path){

    data.frame(read.table(bed.path, header=F))

}

