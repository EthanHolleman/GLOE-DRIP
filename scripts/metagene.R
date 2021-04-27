library(metagene)
data(promoters_hg19)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)



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
    
    args.length <- length(args)  # should be even first half is bam files second is names 
    args.no_out <- args[3:length(args)]
    midpoint <- length(args.no_out) / 2

    bam_files <- args.no_out[1:midpoint]
    bam_names <- args.no_out[(midpoint+1):length(args.no_out)]
    names(bam_files) <- bam_names
    

    design_df <- read_design_table(args[2], bam_files)
    print('adding design table')
    design_df$Samples <- paste0(system.file("extdata", package="metagene"), "/",
                        design_df$Samples)
    print('added design')
    print(design_df)
    print(bam_files)
    mg <- metagene$new(regions=promoters_hg19, bam_files=bam_files)
    print('made metagene promotors')
    mg.df <- mg$get_data_frame()
    print('got df')
    print(head(mg.df))
    write.table(mg.df, paste(args.output_path, 'tsv', sep='.'))
    plt <- plot_metagene(mg.df)
    ggsave(output_path, plt, dpi=500, height=10, width=12)
}





if (! interactive()){
    main()
}


