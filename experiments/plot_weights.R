library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)




#######################################
#         Plot parameters             #
#######################################
L_thres = 0.6
intv = c(1.95,3.5)
grid_size = 0.1

fill1 <- "#0033CC"
fill2 <- "#66CCFF"
col1 <- "#0033CC"
col2 <- "#66CCFF"
col3 <- "#333333"
border <- "#666666"

line_size = 0.75
line_alpha = 0.6
fill_alpha = 0.3

plot_width = 10
plot_height = 3.5
legend_size = 23
font_size = 16
axis_size = 15

fig.dir = "./pics"

#######################################
#         Plot the CLfdr              #
#######################################
grid = seq(intv[1],intv[2], grid_size)
L1 <- function(x){
  dnorm(x,0,1)/(0.9*dnorm(x,0,1)+0.1*dnorm(x,4,1))
}
L2 <- function(x){
  dnorm(x,0,1)/(0.1*dnorm(x,0,1)+0.9*dnorm(x,4,1))
}
x1 <- uniroot(function(x) L1(x)-L_thres, intv)$root
x2 <- uniroot(function(x) L2(x)-L_thres, intv)$root

L <- data.frame(
  x = grid,
  l1 = L1(grid), 
  l2 = L2(grid)
)

L_plot <- ggplot(L, aes(x = x)) + 
          geom_line(aes(y = l1, color = "l1"), size = line_size, alpha = line_alpha) +
          geom_line(aes(y = l2, color = "l2"), size = line_size, alpha = line_alpha) +
          geom_hline(yintercept=L_thres, linetype="dashed", 
                    color = col3, size=line_size, alpha = line_alpha) +
          annotate('segment', x=x1, y=0, xend=x1, yend=L_thres, linetype="dashed", 
                     color = col3, size=line_size, alpha = line_alpha) +
          annotate('segment',x=x2, y=0, xend=x2, yend=L_thres, linetype="dashed", 
                     color = col3, size=line_size, alpha = line_alpha) +
          scale_color_manual(
                    values = c(col1, col2),
                    labels = c(expression(L[1]), expression(L[2]))
                  ) +
          theme_light()+
          theme(legend.position = c(0.9, 0.9), legend.justification = c(1, 1),
                legend.key.size = unit(legend_size, "pt"),
                legend.title = element_blank(),
                legend.box.background = element_rect(color = border),
                legend.box.margin = margin(0.4,0.4,0.4,0.4,"pt"),
                axis.title = element_text(size = font_size),
                axis.text = element_text(size = axis_size),
                legend.text = element_text(size = axis_size)) +
          labs(x = expression(t), y="CLfdr(t)") 



#######################################
#         Plot the weights            #
#######################################
w_pos <- data.frame(
  x = grid,
  y = dnorm(grid)
)

w_plot <- ggplot(w_pos, aes(x = x, y = y)) + 
          geom_line(size=line_size, alpha = line_alpha) + 
          geom_area(data = subset(w_pos, x >= x1),
                    aes(y=y, fill="x1"),color = NA, alpha = fill_alpha) +
          geom_area(data = subset(w_pos, x >= x2),
                    aes(y=y, fill="x2"),color = NA, alpha = fill_alpha) +
          scale_fill_manual(
                    values = c(fill1,fill2),
                    labels = c(expression(w[1]), expression(w[2]))
                    ) +
          theme_light() +
          theme(legend.position = c(0.9, 0.9), legend.justification = c(1, 1),
                legend.key.size = unit(legend_size, "pt"),
                legend.title = element_blank(),
                legend.box.background = element_rect(color = border),
                legend.box.margin = margin(0.4,0.4,0.4,0.4,"pt"),
                axis.title = element_text(size = font_size),
                axis.text = element_text(size = axis_size),
                legend.text = element_text(size = axis_size)) +
          labs(x = expression(t), y=expression(f[0](t)))



#######################################
#              Save plots             #
#######################################
weight_vis <- arrangeGrob(L_plot, w_plot, nrow=1, ncol = 2)
grid.draw(weight_vis)
ggsave(filename = sprintf("%s/weight_visualization.pdf", fig.dir), plot = weight_vis,
        dpi = 300, device=NULL, width=plot_width, height=plot_height)

