library(metagene)
data(promoters_hg19)



# bam files should be specified with name in command line args



main <- function(){

    args <- commandArgs(trailingOnly=TRUE)
    output_path <- args[1]
    args.length <- length(args)-1  # should be even first half is bam files second is names 
    bam_files <- args[2:args.length/2]
    bam_names <- args[args.length/2:args.length]
    names(bam_files) <- bam_names
    mg <- metagene$new(regions=promoters_hg19, bam_files=bam_files, flip_regions=TRUE)
    plt <- mg$plot(title='2kb promotor aligned GLOE-seq reads')
    ggsave(output_path, plt, dpi=500)
}

metagene_plot <- function(bed)



# get the bam files to plot
# name the bam files based on command line inputs 

