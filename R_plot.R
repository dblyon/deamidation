##### Documentation for plots
library(ggplot2)

fn <- "/Users/dblyon/Downloads/Output/Deamidation.txt"
fn_plot <- gsub("txt$", "pdf", fn)
df <- read.csv(fn, sep="\t", header = TRUE)
df <- df[with(df, order(mean)), ]
df$sort_order <- seq(nrow(df))
df$Raw.File <- factor(df$Raw.File, levels=unique(df$Raw.File[df$sort_order]), ordered=TRUE)
dodge <- position_dodge(width=0.9)
limits <- aes(ymin=(mean-std), ymax=(mean+std))
fn <- "/Users/dblyon/Downloads/Output/Number_of_Peptides_per_RawFile.txt"
df_text <- read.csv(fn, sep="\t", header = TRUE)
df_text$Raw.File <- factor(df_text$Raw.File, levels=unique(df$Raw.File[df$sort_order]), ordered=TRUE)
plt <- ggplot(aes(y=mean, x=N_Q, fill=N_Q), data=df) + 
  geom_bar(position=dodge, stat="identity") + geom_errorbar(limits, position=dodge, width=0.25) + facet_grid(. ~ Raw.File) + 
  geom_text(data=df_text, aes(x=N_Q, y=105, label=num_peptides), size=2.5) + 
  labs(title="Deamidation per RawFile STD as error bars") + xlab("Asparagine Glutamine") + ylab("percent deamidation") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#7fcdbb", "#2c7fb8"),
                    name="Amino acid",
                    breaks=c("N", "Q"),
                    labels=c("Asparagine", "Glutamine"))
plt
ggsave(fn_plot, plot=last_plot(), width=12, height=9)

