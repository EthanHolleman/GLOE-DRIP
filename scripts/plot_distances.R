library(ggplot2)
library(ggpubr)

read_dist <- function(dist.path){

    as.data.frame(read.table(dist.path, header=F))

}

plot_distance <- function(dist.df.fwd, dist.df.rev){

    dist.df.fwd$strand <- '-'
    dist.df.rev$strand <- '+'
    dist.df.all <- rbind(dist.df.fwd, dist.df.rev)
    dist <- ggplot(dist.df.all, aes(y=`V9`, x=as.factor(strand), fill=strand)) + geom_boxplot(outlier.shape=NA) +
    #scale_y_continuous(limits = quantile(dist.df.all$signal_diff, c(0.1, 0.9), na.rm=TRUE)) +
    theme_pubr() + stat_compare_means(method='t.test') + labs(x='', y='Distance', main='Distance to closest replicate break')

    print(dist.df.all$V4)
    dist.df.all$signal_diff <- abs(as.numeric(dist.df.all$V4) - as.numeric(dist.df.all$V8))

    signal <- ggplot(dist.df.all, aes(y=signal_diff, x=as.factor(strand), fill=strand)) + geom_boxplot(outlier.shape=NA) +
    labs(x='', y='Read count difference', main='Read count differences') +  stat_compare_means(method='t.test') + 
    #scale_y_continuous(limits = quantile(dist.df.all$signal_diff, c(0.1, 0.9), na.rm=TRUE)) +
    theme_pubr()

    ggarrange(dist, signal, nrow=1, ncol=2)
}



main <- function(){

    args = commandArgs(trailingOnly=TRUE)
    dist.df.fwd <- read_dist(args[1])
    dist.df.rev <- read_dist(args[2])
    plt <- plot_distance(dist.df.fwd, dist.df.rev)
    ggsave(args[3], plt)

}

if (! interactive()){
    main()
}