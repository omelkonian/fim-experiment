library(ggplot2)
library(scales)

allE <- list(0.1, 0.05, 0.01)
allD <- list(0.01, 0.001)
allDat <- list(
  "chess.dat",
  "mushroom.dat",
  "T10I4D100K.txt",
  "retail.dat",
  "connect.dat",
  "pumsb_star.dat",
  "pumsb.dat",
  "T40I10D100K.dat",
  "kosarak.dat",
  "accidents.dat"
)
ext <- "png"
for (db in allDat) {
  path <- sprintf("output/%s", db)
  
  # Load experiment
  count <- 1
  for (e in allE) {
    for (d in allD) {
      varName <- sprintf("(ε=%g)(δ=%g)", e, d)
      load(sprintf("%s/%s", path, varName))
      assign(paste(list("c"), count, sep=""), get(varName))
      count <- count + 1
    }}
  
  #### SAMPLE SIZES ####
  pPath <- sprintf("%s/sampleSizes.%s", path, ext)
  if (!file.exists(pPath)) {
    s <- sprintf("Method Configuration SampleSizes
    Toivonen c1 %d
    Riondato c1 %d
    Toivonen c2 %d
    Riondato c2 %d
    Toivonen c3 %d
    Riondato c3 %d
    Toivonen c4 %d
    Riondato c4 %d
    Toivonen c5 %d
    Riondato c5 %d
    Toivonen c6 %d
    Riondato c6 %d
    ", c1$T$Σ, c1$R$Σ, c2$T$Σ, c2$R$Σ, c3$T$Σ, c3$R$Σ, c4$T$Σ, c4$R$Σ, c5$T$Σ, c5$R$Σ, c6$T$Σ, c6$R$Σ)
    
    SAMPLE.SIZES <- read.table(header=TRUE, text=s)
    ggsave(pPath, ggplot(SAMPLE.SIZES, aes(Configuration, SampleSizes, fill=Method)) +
      geom_bar(stat="identity", position = "dodge") +
      scale_fill_brewer(palette = "Set1"))
  }

  #### +- ####
  # Toivonen
  for (i in 1:6) {
    cc <- get(sprintf("c%d", i))
    df <- data.frame(Itemsets = c("true positives", "false positives", "false negatives"), value = c(cc$T$`+`, cc$T$`-+`, cc$T$`--`))
    total <- sum(df$value)
    
    pPath <- sprintf("%s/toivonen-c%d.%s", path, i, ext)
    if (!file.exists(pPath)) {
      ggsave(pPath, ggplot(df, aes(x="", y=value/total, fill=Itemsets)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=c("#00bfff", "#FF0000", "#7CFF00")) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position="none",
              plot.title=element_text(size=14, face="bold")))
    }
  }
  # Riondato
  for (i in 1:6) {
    cc <- get(sprintf("c%d", i))
    df <- data.frame(Itemsets = c("true positives", "false positives", "false negatives"), value = c(cc$R$`+`, cc$R$`-+`, cc$R$`--`))
    total <- sum(df$value)
    
    pPath <- sprintf("%s/riondato-c%d.%s", path, i, ext)
    if (!file.exists(pPath)) {
      ggsave(pPath, ggplot(df, aes(x="", y=value/total, fill=Itemsets)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        scale_fill_manual(values=c("#00bfff", "#FF0000", "#7CFF00")) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position="none",
              plot.title=element_text(size=14, face="bold")))
    }
  }
  
  #### PR ####
  
  # Precision
  pPath <- sprintf("%s/precision.%s", path, ext)
  if (!file.exists(pPath)) {
    df2 <- data.frame(Metric=rep(c("Toivonen", "Riondato"), each=6),
                      Configuration=rep(c("c1", "c2", "c3", "c4", "c5", "c6"),2),
                      Precision=c(c1$T$precision,c2$T$precision,c3$T$precision,c4$T$precision,c5$T$precision,c6$T$precision,
                                  c1$R$precision,c2$R$precision,c3$R$precision,c4$R$precision,c5$R$precision,c6$R$precision)
    )
    ggsave(pPath, ggplot(df2, aes(x=Configuration, y=Precision, group=Metric, color=Metric)) +
             geom_line(aes(linetype=Metric)) +
             geom_point(aes(shape=Metric), size=3) + 
             theme(legend.position="top",
                   legend.key.size = unit(1, "cm"),
                   legend.text = element_text(size=12, face="bold"),
                   legend.title = element_text(size=12)))
  }
  # Precision
  pPath <- sprintf("%s/recall.%s", path, ext)
  if (!file.exists(pPath)) {
    df2 <- data.frame(Metric=rep(c("Toivonen", "Riondato"), each=6),
                      Configuration=rep(c("c1", "c2", "c3", "c4", "c5", "c6"),2),
                      Recall=c(c1$T$recall,c2$T$recall,c3$T$recall,c4$T$recall,c5$T$recall,c6$T$recall,
                                c1$R$recall,c2$R$recall,c3$R$recall,c4$R$recall,c5$R$recall,c6$R$recall)
                      )
    ggsave(pPath, ggplot(df2, aes(x=Configuration, y=Recall, group=Metric, color=Metric)) +
      geom_line(aes(linetype=Metric)) +
      geom_point(aes(shape=Metric), size=3) + 
      theme(legend.position="top",
            legend.key.size = unit(1, "cm"),
            legend.text = element_text(size=12, face="bold"),
            legend.title = element_text(size=12)))
  }
  
  #### FREQ ####
  pPath <- sprintf("%s/freq.%s", path, ext)
  if (!file.exists(pPath)) {
    df2 <- data.frame(Metric=rep(c("Toivonen-max", "Toivonen-avg", "Riondato-max", "Riondato-avg"), each=6),
                      Configuration=rep(c("c1", "c2", "c3", "c4", "c5", "c6"),2),
                      FrequencyError=c(c1$T$maxFreq,c2$T$maxFreq,c3$T$maxFreq,c4$T$maxFreq,c5$T$maxFreq,c6$T$maxFreq,
                               c1$T$avgFreq,c2$T$avgFreq,c3$T$avgFreq,c4$T$avgFreq,c5$T$avgFreq,c6$T$avgFreq,
                               c1$R$maxFreq,c2$R$maxFreq,c3$R$maxFreq,c4$R$maxFreq,c5$R$maxFreq,c6$R$maxFreq,
                               c1$R$avgFreq,c2$R$avgFreq,c3$R$avgFreq,c4$R$avgFreq,c5$R$avgFreq,c6$R$avgFreq)
    )
    ggsave(pPath, ggplot(df2, aes(x=Configuration, y=FrequencyError, group=Metric, color=Metric)) +
     geom_line(aes(linetype=Metric)) +
     geom_point(aes(shape=Metric), size=3) +
     scale_color_manual(values=c("#ff8080", "#cc0000", "#99e699", "#1f7a1f")) +
     theme(legend.position="top",
           legend.key.size = unit(1, "cm"),
           legend.text = element_text(size=12, face="bold"),
           legend.title = element_text(size=12)))
  }
}