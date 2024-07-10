
tree_theme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                    panel.grid.minor=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.border = element_rect(fill='transparent', color='black'),
                    plot.title = element_text(size = 16, hjust = 0),
                    legend.key = element_rect( fill = "white"),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks = element_blank(),
                    legend.text = element_text(vjust = 0.4, size = 12, colour = 'black'),
                    legend.title = element_text(vjust = 0.4, size = 12, colour = 'black'),
                    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


graph.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                      panel.border = element_rect(fill='transparent', color='black'),
                      plot.title = element_text(size = 10, hjust = 0), # title
                      plot.subtitle = element_text(size = 10, hjust = 0), # subtitle
                      legend.key = element_rect( fill = "white"),
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      axis.text = element_blank(),
                      legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
                      legend.title = element_text(vjust = 0.4, size = 13, colour = 'black'),
                      legend.key.size = unit(0.2, "cm") )


traj.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                     #panel.grid.major=element_line(colour='transparent', color='white'), # background line
                     #strip.background = element_rect(color = "black"),
                     legend.position="none",
                     strip.text = element_text(size = 10),
                     panel.grid.minor=element_blank(),
                     panel.grid.major=element_blank(),
                     panel.border = element_rect(fill='transparent', color='black'),
                     plot.title = element_text(size = 14, hjust = 0), # title
                     plot.subtitle = element_text(size = 13, hjust = 0), # subtitle
                     legend.key = element_rect( fill = "white"),
                     axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
                     axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
                     #axis.line = element_line(colour = 'black'),
                     axis.ticks = element_blank(),
                     #axis.text.y = element_blank(),
                     axis.text.x = element_text(vjust = -0.5, size = 13, colour = 'black'),
                     axis.text.y = element_text(vjust = 0.5, size = 13, colour = 'black'),
                     legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                     legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                     legend.key.size = unit(0.9, "cm") )

rj.graph.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                   #panel.grid.major=element_line(colour='transparent', color='white'), # background line
                   #strip.background = element_rect(color = "black"),
                   strip.text = element_text(size = 20),
                   # panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.border = element_rect(fill='transparent', color='black'),
                   plot.title = element_text(size = 25, hjust = 0), # title
                   plot.subtitle = element_text(size = 25, hjust = 0), # subtitle
                   legend.key = element_rect( fill = "white"),
                   axis.title = element_blank(),
                   # axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
                   # axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
                   #axis.line = element_line(colour = 'black'),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   # axis.text.x = element_text(vjust = -0.5, size = 15, colour = 'black'),
                   # axis.text.y = element_text(vjust = 0.5, size = 15, colour = 'black'),
                   legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                   legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                   legend.key.size = unit(0.9, "cm") )


rj.graph.ftheme2 <- theme(panel.background = element_rect(fill='transparent', color="black"),
                         #panel.grid.major=element_line(colour='transparent', color='white'), # background line
                         #strip.background = element_rect(color = "black"),
                         strip.text = element_text(size = 20),
                         panel.grid.minor=element_blank(),
                         panel.grid.major=element_blank(),
                         panel.border = element_rect(fill='transparent', color='black'),
                         plot.title = element_text(size = 15, hjust = 0), # title
                         plot.subtitle = element_text(size = 15, hjust = 0), # subtitle
                         legend.key = element_rect( fill = "white"),
                         axis.title = element_blank(),
                         # axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'), # face = "bold"
                         # axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'), # face = "bold"
                         #axis.line = element_line(colour = 'black'),
                         axis.ticks = element_blank(),
                         axis.text = element_blank(),
                         # axis.text.x = element_text(vjust = -0.5, size = 15, colour = 'black'),
                         # axis.text.y = element_text(vjust = 0.5, size = 15, colour = 'black'),
                         legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                         legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                         legend.key.size = unit(0.9, "cm") )

rj.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                   strip.text = element_text(size = 18),
                   panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.border = element_rect(fill='transparent', color='black'),
                   plot.title = element_text(size = 16, hjust = 0),
                   plot.subtitle = element_text(size = 16, hjust = 0),
                   legend.key = element_rect( fill = "white"),
                   axis.title.x = element_text(vjust = -1.5, size = 15, colour = 'black'),
                   axis.title.y = element_text(vjust = 1.5, size = 15, colour = 'black'),
                   axis.ticks = element_blank(),
                   axis.text.x = element_text(vjust = -0.5, size = 14, colour = 'black'),
                   axis.text.y = element_text(vjust = 0.5, size = 14, colour = 'black'),
                   legend.text = element_text(vjust = 0.4, size = 14, colour = 'black'),
                   legend.title = element_text(vjust = 0.4, size = 14, colour = 'black'),
                   # legend.key.size = unit(0.9, "cm"),
                   plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

rj.ftheme.small <- theme(panel.background = element_rect(fill='transparent', color="black"),

                   strip.text = element_text(size = 18),
                   panel.grid.minor=element_blank(),
                   panel.grid.major=element_blank(),
                   panel.border = element_rect(fill='transparent', color='black'),
                   plot.title = element_text(size = 18, hjust = 0),
                   plot.subtitle = element_text(size = 18, hjust = 0),
                   legend.key = element_rect( fill = "white"),
                   axis.title.x = element_text(vjust = -1.5, size = 14, colour = 'black'),
                   axis.title.y = element_text(vjust = 1.5, size = 14, colour = 'black'),
                   axis.ticks = element_blank(),
                   axis.text.x = element_text(vjust = -0.5, size = 13, colour = 'black'),
                   axis.text.y = element_text(vjust = 0.5, size = 13, colour = 'black'),
                   legend.text = element_text(vjust = 0.4, size = 13, colour = 'black'),
                   legend.title = element_text(vjust = 0.4, size = 13, colour = 'black'),
                   # legend.key.size = unit(0.9, "cm"),
                   plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"))

rj.ftheme.merge <- theme(panel.background = element_rect(fill='transparent', color="black"),
                         #panel.grid.major=element_line(colour='transparent', color='white'), # background line
                         #strip.background = element_rect(color = "black"),
                         strip.text = element_text(size = 20),
                         panel.grid.minor=element_blank(),
                         panel.grid.major=element_blank(),
                         panel.border = element_rect(fill='transparent', color='black'),
                         plot.title = element_text(size = 25, hjust = 0), # title
                         plot.subtitle = element_text(size = 25, hjust = 0), # subtitle
                         legend.key = element_rect( fill = "white"),
                         axis.title.x = element_text(vjust = -1.5, size = 13, colour = 'black'), # face = "bold"
                         axis.title.y = element_text(vjust = 1.5, size = 13, colour = 'black'), # face = "bold"
                         #axis.line = element_line(colour = 'black'),
                         axis.ticks = element_blank(),
                         #axis.text.y = element_blank(),
                         axis.text.x = element_text(vjust = -0.5, size = 20, colour = 'black'),
                         axis.text.y = element_text(vjust = 0.5, size = 20, colour = 'black'),
                         legend.text = element_text(vjust = 0.4, size = 20, colour = 'black'),
                         legend.title = element_text(vjust = 0.4, size = 20, colour = 'black'),
                         legend.key.size = unit(0.9, "cm") )

rj.graph.ggplot.ftheme <- theme(panel.background = element_rect(fill='transparent', color="black"),
                                #panel.grid.major=element_line(colour='transparent', color='white'), # background line
                                #strip.background = element_rect(color = "black"),
                                #legend.position="none",
                                strip.text = element_text(size = 10),
                                panel.grid.minor=element_blank(),
                                panel.grid.major=element_blank(),
                                panel.border = element_rect(fill='transparent', color='black'),
                                plot.title = element_text(size = 14, hjust = 0), # title
                                plot.subtitle = element_text(size = 13, hjust = 0), # subtitle
                                legend.key = element_rect( fill = "white"),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank(),
                                axis.text.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks = element_blank(),
                                legend.text = element_text(vjust = 0.4, size = 15, colour = 'black'),
                                legend.title = element_text(vjust = 0.4, size = 15, colour = 'black'),
                                legend.key.size = unit(0.9, "cm") )

rj.ftheme.venn <- theme(panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.title = element_text(size = 15, hjust = 0.5, vjust = -25), # title
                        plot.subtitle = element_text(size = 25, hjust = 0),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank(),
                        legend.key.size = unit(0.9, "cm") )

rj.ftheme.venn2 <- theme(panel.background = element_blank(),
                        panel.border = element_blank(),
                        plot.title = element_text(size = 15, hjust = 0.5, vjust = -5), # title
                        plot.subtitle = element_text(size = 25, hjust = 0),
                        axis.title = element_blank(),
                        axis.ticks = element_blank(),
                        axis.text = element_blank(),
                        legend.key.size = unit(0.9, "cm") )




