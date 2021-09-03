#' Protected Connected (ProtConn)
#'
#' Estimate Protected Connected (ProtConn) indicator and fractions for one region.
#' @param DataProtconn data.frame
#' @param d distance threshold
#' @importFrom ggplot2 ggplot geom_bar aes position_dodge labs rel theme_bw theme element_blank element_text scale_fill_manual geom_hline scale_linetype_manual guide_legend margin
#' @importFrom purrr compact
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
#' @keywords internal
plotprotconn <- function(DataProtconn, d){
  plot_protconn <- list()
  Data_plot <- as.data.frame(t(round(DataProtconn, 3)))
  Data_plot$name <- row.names(Data_plot)
  names(Data_plot) <- c("Value", "Indicator")
  rownames(Data_plot)<- NULL
  Data_plot$Indicator <- as.character(Data_plot$Indicator)
  Data_plot_1 <- Data_plot[which(Data_plot$Indicator %in% c("Unprotected", "Prot","ProtConn")),]
  Data_plot_1$Indicator[which(Data_plot_1$Indicator == "ProtConn")] <-  "Protected connected"
  Data_plot_1$Indicator[which(Data_plot_1$Indicator == "Prot")] <-  "Protected"

  Data_plot_1$Indicator <- factor(Data_plot_1$Indicator, levels = c("Unprotected", "Protected", "Protected connected"))
  names(Data_plot_1) <- c("Values", "name")
  Data_plot_1$col <- c("#53A768", "#C34D51", "#4C72AF")
  Data_plot_1[which(Data_plot_1$Values == 0), ] <- NULL

  if(nrow(Data_plot_1) > 1){
    plot_protconn1 <- ggplot(Data_plot_1, aes(x = Data_plot_1$name, y = Data_plot_1$Values,
                                              fill = Data_plot_1$name)) +
      geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
      labs(title = paste0("ProtConn Indicators: ", d), x = "", y = "Percentage (%)", size = rel(1.2)) +
      theme_bw()  +
      theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.5, face = "bold"),
            axis.title = element_text(color = "#252525", size = rel(1.2)),
            legend.title = element_blank(),
            legend.text = element_text(colour = "#252525", size = rel(1.2)),
            axis.text= element_text(colour = "#525252", size = rel(1)))+
      scale_fill_manual(values = Data_plot_1$col) +
      geom_hline(aes(yintercept = 17, linetype = "Aichi Target (17%)"), colour = 'black', size = 1.2) +
      scale_linetype_manual(name = " Aichi Target", values = c(2, 2),
                            guide = guide_legend(override.aes = list(color = c("black"), size = 0.8)))

    plot_protconn[[1]] <- plot_protconn1
  }

  Data_plot_2 <- Data_plot[which(Data_plot$Indicator %in% c("ProtConn_Trans", "ProtConn_Unprot", "ProtConn_Within", "ProtConn_Contig")),]
  Data_plot_2$Indicator <-  c("ProtConn[Trans]", "ProtConn[Unprot]", "ProtConn[Within]", "ProtConn[Contig]")
  Data_plot_2$Indicator <- factor(Data_plot_2$Indicator, levels = c(Data_plot_2$Indicator[3], Data_plot_2$Indicator[4], Data_plot_2$Indicator[2], Data_plot_2$Indicator[1]))
  names(Data_plot_2) <- c("Values", "name")
  Data_plot_2$col <- c("#253494", "#2c7fb8", "#41b6c4", "#7fcdbb")
  Data_plot_2 <- Data_plot_2[which(Data_plot_2$Values > 0), ]

  if(nrow(Data_plot_2) > 1){
    plot_protconn2 <- ggplot(Data_plot_2, aes(x = Data_plot_2$name, y = Data_plot_2$Values,
                                              fill = Data_plot_2$name)) +
      geom_bar(position = position_dodge(), colour = "black", stat = "identity", show.legend = FALSE, size = 0.2) +
      labs(title = "Protected connected fraction", x = "", y = "Percentage (%)", size = rel(1.2)) +
      theme_bw()  +
      theme(plot.title = element_text(color = "#252525", size = rel(1.4), hjust = 0.45, face = "bold"),
            axis.title= element_text(color = "#252525", size = rel(1.2)),
            plot.margin = margin(0, 5.3, 0, 0.3, "cm"),
            legend.title= element_blank(),
            legend.text = element_text(colour = "#252525", size = rel(1.2)),
            axis.text= element_text(colour = "#525252", size = rel(1)))+
      scale_fill_manual(values = Data_plot_2$col)
    plot_protconn[[2]] <- plot_protconn2
  }
  plot_protconn <- compact(plot_protconn)

  if(length(plot_protconn) == 2){
    figure <- ggarrange(plot_protconn[[1]], plot_protconn[[2]],
                        ncol = 1, nrow = 2)
  } else if (length(plot_protconn) == 1) {
    figure <- plot_protconn[[1]]
  } else {
    figure <- "There are insufficient data to plot"
  }

  return(figure)
}

