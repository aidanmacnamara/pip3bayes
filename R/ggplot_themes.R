themeS <- theme(legend.position="none",          
                axis.title.x = element_text(size=12),
                axis.title.y = element_text(size=12),
                axis.text.x = element_text(size=10),
                axis.text.y = element_text(size=10), 
                plot.margin = unit(c(3,3,3,3), "mm")
)

themeHT <- theme(legend.position="none",          
                 axis.title.x = element_blank(),
                 axis.title.y = element_text(size=12),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(size=10), 
                 axis.ticks.x = element_blank(),
                 plot.margin = unit(c(3,3,3,7.5), "mm")
)

themeHR <- theme(legend.position="none",          
                 axis.title.y = element_blank(),
                 axis.title.x = element_text(size=12),
                 axis.text.x = element_text(size=10),
                 axis.text.y = element_blank(), 
                 axis.ticks.y = element_blank(),
                 plot.margin = unit(c(3,3,3,3), "mm")
)

histEmp <- ggplot() + geom_point(aes(1,1), colour="white") +
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank()
  )
