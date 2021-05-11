library(metagene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



read_design_table <- function(table.path, named_bam_files){
    print('Reading table')
    design.df <- as.data.frame(read.table(table.path))
    samples <- rownames(design.df)
    rownames(design.df) <- named_bam_files[samples]
    design.df$Samples <- rownames(design.df)
    design.df <- design.df[,order(ncol(design.df):1)]
    rownames(design.df) <- c()
    design.df

}

plot_metagene <- function(metagene.df){
    print('Plotting')
    metagene.df$design <- as.factor(metagene.df$design)

    ggplot(metagene.df, aes(x=bin, y=value, ymin=qinf, ymax=qsup)) + 
    geom_ribbon(aes(fill=design), alpha=0.3) +
    geom_line(aes(color=design), size=1) + facet_wrap(~group) + 
    scale_color_brewer(palette="Dark2") +
    theme_pubr()
           
}

# bam files should be specified with name in command line args
main <- function(){

    args <- commandArgs(trailingOnly=TRUE)
    output_path <- args[1]
    cores <- as.numeric(args[3])
    bed_regions.path <- args[4]
    #genes <- genes(txdb)
    args.length <- length(args)  # should be even first half is bam files second is names 
    args.no_out <- args[5:length(args)]
    midpoint <- length(args.no_out) / 2

    bam_files <- args.no_out[1:midpoint]
    bam_names <- args.no_out[(midpoint+1):length(args.no_out)]
    named_bam_files <- bam_files
    names(named_bam_files) <- bam_names
    for (bam in bam_files){
        if(! file.exists(bam)){
            print(bam)
            print('Does not exist')
        }
    }
    
    print(named_bam_files)
    design_df <- read_design_table(args[2], named_bam_files)
    print('adding design table')
    #design_df$Samples <- paste0(system.file("extdata", package="metagene"), "/",
    #                    design_df$Samples)
    print(design_df)
    print('added design')
    print(bam_files)
    mg <- metagene$new(regions=bed_regions.path, bam_files=bam_files, cores=cores)
    print('made metagene promotors')
    mg$produce_table(design=design_df)
    mg$produce_data_frame()
    mg.df <- mg$get_data_frame()
    print('got df')
    print(head(mg.df))
    write.table(mg.df, paste(output_path, 'tsv', sep='.'))
    plt <- plot_metagene(mg.df)
    print(output_path)
    ggsave(output_path, plt, dpi=500, height=10, width=12)
}





if (! interactive()){
    main()
}


