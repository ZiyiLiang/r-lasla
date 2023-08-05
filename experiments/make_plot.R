library(tidyverse)
library(kableExtra)

options(width=160)

#######################################
#          Plots for paper            #
#######################################
fig.dir = "./pics"
results.dir = "./results"

key.values <- c("FDR", "Power")
key.labels <- c("FDR", "Power")

Method.values <- c("PV.OR", "LASLA.OR", "LAWS.OR")
Method.labels <- c("PV.OR", "LASLA.OR", "LAWS.OR")

color.scale <- c("#3366CC", "#66CCFF", "#CC79A7", "orange", "red")
shape.scale <- c(15, 2, 4, 1, 18)


plot_width = 10
plot_height = 3
legend_size = 14
font_size = 17
axis_size = 15
pt_size = 2
pt_alpha = 0.8
l_alpha = 0.95
l_size = 0.7

#######################################
#       Asymmetric setting            #
#######################################
plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
              mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.OR") %>%
            mutate(Method = factor(Method, Method.values, Method.labels)) %>%
            mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
      pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
      mutate(Key = factor(Key, key.values, key.labels))  %>%
ggplot(aes(x=Gamma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
geom_hline(data=df.nominal, aes(yintercept=Value)) +
# sets the axis limit
geom_point(data=df.ghost, aes(x=1,y=0.03), alpha=0) +
geom_point(data=df.ghost, aes(x=1,y=0.07), alpha=0) +
scale_color_manual(values=color.scale) +
scale_shape_manual(values=shape.scale) +
scale_alpha_manual(values=alpha.scale) +
theme_bw()+
theme(
  strip.text = element_text(size = font_size, color = "black"),
  axis.title = element_text(size = font_size),
  axis.text = element_text(size = axis_size),
  legend.text = element_text(size = legend_size),
  legend.title = element_text(size = legend_size),
  plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
)+
facet_wrap(.~Key, scales="free") +
xlab(expression(gamma)) +
ylab("")+
ggtitle("(a)")
pp

ggsave(filename = sprintf("%s/asymmetric.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Alt-shape setting            #
#######################################
plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.OR") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=1,y=0.03), alpha=0) +
  geom_point(data=df.ghost, aes(x=1,y=0.07), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(sigma)) +
  ylab("")+
  ggtitle("(b)")



ggsave(filename = sprintf("%s/alt_dist.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Network1  setting            #
#######################################
Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Mean, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=2.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.5,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(mu[1])) +
  ylab("")+
  ggtitle("(a)")


ggsave(filename = sprintf("%s/network1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Network2  setting            #
#######################################
Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Distance, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=1,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(mu[2])) +
  ylab("")+
  ggtitle("(b)")


ggsave(filename = sprintf("%s/network2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#         Latent1  setting            #
#######################################
Method.values <- c("BH", "LASLA.OR","LASLA.DD", "SABHA", "WBH")
Method.labels <- c("BH", "LASLA.OR","LASLA.DD", "SABHA", "WBH")
color.scale <- c("#CC79A7","#66CCFF", "#3366CC", "orange", "red")


plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(sigma)) +
  ylab("")+
  ggtitle("(a)")


ggsave(filename = sprintf("%s/latent1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#         Latent2  setting            #
#######################################
plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Mean, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=3.0,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=4.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(mu)) +
  ylab("")+
  ggtitle("(b)")


ggsave(filename = sprintf("%s/latent2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#       Latent-Mult setting 1        #
#######################################
Method.values <- c("BH", "LASLA.OR","LASLA.DD","AVG")
Method.labels <- c("BH", "LASLA.OR","Mahalanobis","Average")
color.scale <- c("#CC79A7", "#66CCFF", "#3366CC","orange")


plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.OR") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale, labels=c("BH"="BH", "LASLA.OR"="LASLA.OR","LASLA.DD"="Mahalanobis","AVG"="Average"),,breaks=Method.values) +
  scale_shape_manual(values=shape.scale, labels=c("BH"="BH", "LASLA.OR"="LASLA.OR","LASLA.DD"="Mahalanobis","AVG"="Average"),breaks=Method.values) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(sigma)) +
  ylab("")+
  ggtitle("(a)")
  pp

ggsave(filename = sprintf("%s/latent_mult1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#       Latent-Mult setting 2        #
#######################################
Method.values <- c("BH", "LASLA.OR","LASLA.DD","AVG")
Method.labels <- c("BH", "LASLA.OR","Mahalanobis","Average")
color.scale <- c("#CC79A7", "#66CCFF", "#3366CC", "orange")


plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.OR") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale, labels=c("BH"="BH", "LASLA.OR"="LASLA.OR","LASLA.DD"="Mahalanobis","AVG"="Average"),,breaks=Method.values) +
  scale_shape_manual(values=shape.scale, labels=c("BH"="BH", "LASLA.OR"="LASLA.OR","LASLA.DD"="Mahalanobis","AVG"="Average"),breaks=Method.values) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(sigma)) +
  ylab("")+
  ggtitle("(b)")
  pp

ggsave(filename = sprintf("%s/latent_mult2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Regression setting 1         #
#######################################
Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")
color.scale <- c("#3366CC", "#66CCFF", "#CC79A7", "orange", "red")


plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.05,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=0.1,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(sigma)) +
  ylab("")+
  ggtitle("(a)")


ggsave(filename = sprintf("%s/regression1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Regression setting 2         #
#######################################
Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")
color.scale <- c("#3366CC", "#66CCFF", "#CC79A7", "orange", "red")


plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Mean, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.25,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=0.35,y=0), alpha=0) +
  scale_color_manual(values=color.scale) +
  scale_shape_manual(values=shape.scale) +
  theme_bw()+
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),
    legend.title = element_text(size = legend_size),
    plot.title = element_text(hjust = -0.11, vjust = -4, size=font_size),
    
  )+
  facet_wrap(.~Key, scales="free") +
  xlab(expression(mu)) +
  ylab("")+
  ggtitle("(b)")


ggsave(filename = sprintf("%s/regression2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)


