library(metagene)
data(promoters_hg19)



# bam files should be specified with name in command line args
main <- function(){

    args <- commandArgs(trailingOnly=TRUE)
    output_path <- args[1]
    args.length <- length(args)  # should be even first half is bam files second is names 
    args.no_out <- args[2:length(args)]
    midpoint <- length(args.no_out) / 2
    bam_files <- args.no_out[1:midpoint]
    bam_names <- args.no_out[(midpoint+1):length(args.no_out)]
    print(bam_names)
    names(bam_files) <- bam_names
    mg <- metagene$new(regions=promoters_hg19, bam_files=bam_files)
    plt <- mg$plot(title='2kb promotor aligned GLOE-seq reads')
    ggsave(output_path, plt, dpi=500)
}

if (! interactive()){
    main()
}


