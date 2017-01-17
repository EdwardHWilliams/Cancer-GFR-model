# ggtheme for plots
my_theme <- function(){
  theme <- theme_bw() +
    theme(
      panel.border = element_rect(colour = NA),
      text = element_text(size = 10, colour="black"), #family="Arial"),
      line = element_line(size = 0.5),
      #                  plot.title = element_text(size = 10),
      #                  axis.title.x=element_text(size = 10),
      #                  axis.title.y=element_text(size = 10),
      #                  plot.title = element_text(size = 10),
      #                  legend.text=element_text(size = 10),
      #                  legend.title=element_text(size = 10), 
      axis.line = element_line(size = 0.5, color = "black"),
      axis.line.x = element_line(size = .5),
      axis.line.y = element_line(size = .5),
      axis.title.y = element_text(vjust =0.1),
      panel.background = element_rect(fill = 'white', colour ='white'),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      panel.grid.major = element_line(colour="#f0f0f0", size = rel(.5)),
      panel.grid.minor = element_blank())
  return(theme)
}



scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}



ggNormQQ <- function(y){ 
  a <- mean(y)
  b <- sd(y)
  q <- 
  ggplot(data.frame(y), aes())}

  

ggQQ <- function(vec, ...) # argument: a linear model
{
  y <- quantile(vec, c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  p <- ggplot(data.frame(vec), aes(sample=vec)) +
    stat_qq(...) +
    geom_abline(slope = slope, intercept = int, color="blue") + 
    xlab("Theoretical quantiles") + ylab("Sample quantiles")
  return(p)
}

