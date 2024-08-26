library(tidyverse)
library(kableExtra)
library(dplyr)
library(patchwork)

options(width=160)

#######################################
#          Plots for paper            #
#######################################
setwd("C:/Users/liang/Documents/GitHub/r-lasla")
fig.dir = "./pics"
results.dir = "./results"

key.values <- c("FDR", "Power")
key.labels <- c("FDR", "Power")

Method.values <- c("PV.OR", "LASLA.OR", "LAWS.OR")
Method.labels <- c("PV.OR", "LASLA.OR", "LAWS.OR")

color.scale <- c("#3366CC", "#66CCFF", "#CC79A7", "yellow", "orange", "red")
shape.scale <- c(15, 2, 4, 3, 1, 18)


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
#    Hetero-asymmetry setting         #
#######################################
load(sprintf("%s/heterogeneous_asymmetric.RData", results.dir))

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
  geom_point(data=df.ghost, aes(x=0.5,y=0.03), alpha=0) +
  geom_point(data=df.ghost, aes(x=0.5,y=0.07), alpha=0) +
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
  xlab(expression(r)) +
  ylab("")
pp

ggsave(filename = sprintf("%s/hetero-asymmetric.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#      Hetero-density setting         #
#######################################
load(sprintf("%s/heterogeneous_distribution.RData", results.dir))

# Define common elements for both plots
common_theme <- theme_bw() +
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = axis_size),
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.text.y = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),  # Adjust legend text size
    legend.title = element_text(size = legend_size),
    plot.margin = margin(5, 5, 5, 5)
  )

# Filter and plot for FDP
pp_FDP <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "FDP") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  coord_cartesian(ylim = c(0, 0.13)) +  # Set y-axis limit for FDP
  labs(x = "Method", y = "") +
  facet_wrap(. ~Key)+
  common_theme +
  theme(legend.position = "none")  # Remove legend for individual plot

# Filter and plot for Power
pp_Power <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "Power") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  labs(x = "Method", y = "") + 
  facet_wrap(. ~Key)+
  common_theme

# Combine the two plots and share a single legend
pp <- pp_FDP + pp_Power + plot_layout(guides = "collect") & theme(legend.position = "right")

# Display the combined plot
pp

ggsave(filename = sprintf("%s/hetero-alt.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width-1, height=plot_height)



#######################################
#      Hetero-sparsity setting        #
#######################################
key.values <- c("FDR", "Power")
key.labels <- c("FDR", "Power")

Method.values <- c("LASLA.OR", "AWBH", "WBH")
Method.labels <- c("LASLA.OR", "Adj-WBH.OR", "WBH.OR")

load(sprintf("%s/heterogeneous_sparsity.RData", results.dir))

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

# Define common elements for both plots
common_theme <- theme_bw() +
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = axis_size),
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.text.y = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),  # Adjust legend text size
    legend.title = element_text(size = legend_size),
    plot.margin = margin(5, 5, 5, 5)
  )

# Filter and plot for FDP
pp_FDP <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "FDP") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  coord_cartesian(ylim = c(0, 0.06)) +  # Set y-axis limit for FDP
  labs(x = NULL, y = NULL) +
  facet_wrap(. ~Key)+
  common_theme +
  theme(legend.position = "none")  # Remove legend for individual plot

# Filter and plot for Power
pp_Power <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "Power") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  coord_cartesian(ylim = c(0.2, 0.85)) +
  labs(x = NULL, y = NULL)+
  facet_wrap(. ~Key)+
  common_theme

# Combine the two plots and share a single legend
pp <- pp_FDP + pp_Power + plot_layout(guides = "collect") & theme(legend.position = "right")

# Display the combined plot
pp
ggsave(filename = sprintf("%s/hetero-sparsity.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width-1, height=plot_height)



#######################################
#        Network1  setting            #
#######################################
load(sprintf("%s/network1_sigma0.7.RData", results.dir))

Method.values <- c("BH", "LASLA.DD", "adapt_cluster_5", "adapt_cluster_10", "adapt_cluster_20")
Method.labels <- c("BH", "LASLA.DD", "5 clusters", "10 clusters", "20 clusters")
color.scale <- c("#3366CC","#66CCFF", "#99CC66", "orange", "#FFCC00")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

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
load(sprintf("%s/network2_sigma0.7.RData", results.dir))

Method.values <- c("BH", "LASLA.DD", "adapt_cluster_5", "adapt_cluster_10", "adapt_cluster_20")
Method.labels <- c("BH", "LASLA.DD", "5 clusters", "10 clusters", "20 clusters")
color.scale <- c("#3366CC","#66CCFF", "#99CC66", "orange", "#FFCC00")

results <- results %>% 
  filter(Method != "adapt_cluster_2") %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

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
load(sprintf("%s/latent1.RData", results.dir))

Method.values <- c("BH", "LASLA.OR","LASLA.DD", "ADAPTMT","SABHA", "WBH", "AWBH")
Method.labels <- c("BH", "LASLA.OR","LASLA.DD", "AdaPT", "SABHA", "WBH", "Adj-WBH")
color.scale5 <- c("#CC79A7","#66CCFF", "#3366CC", "purple", "orange")
color.scale3 <- c("#3366CC", "#CC79A7","#66CCFF")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))

# Create the first plot for the first 5 methods
results_first5 <- results %>%
  filter(Method %in% c("BH", "LASLA.OR", "LASLA.DD", "AdaPT", "SABHA"))

pp_first5 <- results_first5 %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  geom_point(data=df.ghost, aes(x=0.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale5) +  # Limit to first 5 colors
  scale_shape_manual(values=shape.scale[1:5]) +  # Limit to first 5 shapes
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

# Create the second plot for LASLA.DD, WBH, and AWBH
results_lasla_wbh_awbh <- results %>%
  filter(Method %in% c("LASLA.DD", "WBH", "Adj-WBH"))

pp_lasla_wbh_awbh <- results_lasla_wbh_awbh %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Sigma, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  geom_point(data=df.ghost, aes(x=0.5,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=2.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale3) +  # Limit to specific colors for these methods
  scale_shape_manual(values=shape.scale[1:3]) +  # Limit to specific shapes for these methods
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

# Display the two plots separately
pp_first5
pp_lasla_wbh_awbh

ggsave(filename = sprintf("%s/latent1_weight_comp.pdf", fig.dir), plot = pp_first5,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)
ggsave(filename = sprintf("%s/latent1_thres_comp.pdf", fig.dir), plot = pp_lasla_wbh_awbh,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)


#######################################
#         Latent2  setting            #
#######################################
load(sprintf("%s/latent2.RData", results.dir))

Method.values <- c("BH", "LASLA.OR","LASLA.DD", "ADAPT","SABHA", "WBH", "AWBH")
Method.labels <- c("BH", "LASLA.OR","LASLA.DD", "AdaPT", "SABHA", "WBH", "Adj-WBH")
color.scale5 <- c("#CC79A7","#66CCFF", "#3366CC", "purple", "orange")
color.scale3 <- c("#3366CC", "#CC79A7","#66CCFF")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))

# Create the first plot for the first 5 methods
results_first5 <- results %>%
  filter(Method %in% c("BH", "LASLA.OR", "LASLA.DD", "AdaPT", "SABHA"))

pp_first5 <- results_first5 %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Mean, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  geom_point(data=df.ghost, aes(x=3.0,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=4.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale5) +  # Limit to first 5 colors
  scale_shape_manual(values=shape.scale[1:5]) +  # Limit to first 5 shapes
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

# Create the second plot for LASLA.DD, WBH, and AWBH
results_lasla_wbh_awbh <- results %>%
  filter(Method %in% c("LASLA.DD", "WBH", "Adj-WBH"))

pp_lasla_wbh_awbh <- results_lasla_wbh_awbh %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Mean, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  geom_point(data=df.ghost, aes(x=3.0,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=4.0,y=0), alpha=0) +
  scale_color_manual(values=color.scale3) +  # Limit to specific colors for these methods
  scale_shape_manual(values=shape.scale[1:3]) +  # Limit to specific shapes for these methods
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

# Display the two plots separately
pp_first5
pp_lasla_wbh_awbh

ggsave(filename = sprintf("%s/latent2_weight_comp.pdf", fig.dir), plot = pp_first5,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)
ggsave(filename = sprintf("%s/latent2_thres_comp.pdf", fig.dir), plot = pp_lasla_wbh_awbh,
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
load(sprintf("%s/regression1.RData", results.dir))

Method.values <- c("BH", "LASLA.DD", "ADAPT")
Method.labels <- c("BH", "LASLA.DD", "AdaPT")
color.scale <- c("#3366CC","#66CCFF", "#FFCC00")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

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
  geom_point(data=df.ghost, aes(x=0.1,y=0.1), alpha=0) +
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

pp

ggsave(filename = sprintf("%s/regression1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Regression setting 2         #
#######################################
load(sprintf("%s/regression2.RData", results.dir))

Method.values <- c("BH", "LASLA.DD", "ADAPT")
Method.labels <- c("BH", "LASLA.DD", "AdaPT")
color.scale <- c("#3366CC","#66CCFF", "#FFCC00")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

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
  geom_point(data=df.ghost, aes(x=0.25,y=0), alpha=0) +
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

pp

ggsave(filename = sprintf("%s/regression2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Dependent setting 1          #
#######################################
load(sprintf("%s/dependent1.RData", results.dir))

Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")
color.scale <- c("#3366CC","#66CCFF")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Rho, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.8,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=0.8,y=0), alpha=0) +
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
  xlab(expression(rho)) +
  ylab("")
pp

ggsave(filename = sprintf("%s/dependent1.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#        Dependent setting 2          #
#######################################
load(sprintf("%s/dependent2.RData", results.dir))

Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")
color.scale <- c("#3366CC","#66CCFF")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

plot.alpha <- 0.05
df.nominal <- tibble(Key=c("FDR"), Value=plot.alpha) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
df.ghost <- tibble(Key=c("FDR","FDR"), Value=c(0.03,0.07), Method="LASLA.DD") %>%
  mutate(Method = factor(Method, Method.values, Method.labels)) %>%
  mutate(Key = factor(Key, key.values, key.labels))    
pp <- results  %>%
  pivot_longer(cols=c("FDR", "Power"), names_to='Key', values_to='Value') %>%
  mutate(Key = factor(Key, key.values, key.labels))  %>%
  ggplot(aes(x=Adjustment, y=Value, color=Method, shape=Method)) +
  geom_point(alpha=pt_alpha, size=pt_size) +
  geom_line(size=l_size, alpha=l_alpha) +
  geom_hline(data=df.nominal, aes(yintercept=Value)) +
  # sets the axis limit
  geom_point(data=df.ghost, aes(x=0.8,y=0.1), alpha=0) +
  geom_point(data=df.ghost, aes(x=0.8,y=0), alpha=0) +
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
  xlab(expression(a)) +
  ylab("")
pp

ggsave(filename = sprintf("%s/dependent2.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width, height=plot_height)



#######################################
#         Dependent setting 3         #
#######################################
load(sprintf("%s/dependent3.RData", results.dir))

Method.values <- c("BH", "LASLA.DD")
Method.labels <- c("BH", "LASLA.DD")
color.scale <- c("#3366CC","#66CCFF")

results <- results %>%
  mutate(Method = factor(Method, levels = Method.values, labels = Method.labels))

# Define common elements for both plots
common_theme <- theme_bw() +
  theme(
    strip.text = element_text(size = font_size, color = "black"),
    axis.title = element_text(size = axis_size),
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.ticks.x = element_blank(),  # Hide x-axis ticks
    axis.text.y = element_text(size = axis_size),
    legend.text = element_text(size = legend_size),  # Adjust legend text size
    legend.title = element_text(size = legend_size),
    plot.margin = margin(5, 5, 5, 5)
  )

# Filter and plot for FDP
pp_FDP <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "FDP") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  coord_cartesian(ylim = c(0, 0.07)) +  # Set y-axis limit for FDP
  labs(x = NULL, y = NULL) +
  facet_wrap(. ~Key)+
  common_theme +
  theme(legend.position = "none")  # Remove legend for individual plot

# Filter and plot for Power
pp_Power <- results %>%
  pivot_longer(cols = c("FDP", "Power"), names_to = 'Key', values_to = 'Value') %>%
  filter(Key == "Power") %>%
  ggplot(aes(x = Method, y = Value, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, lwd = 0.5) +
  scale_fill_manual(values = color.scale) +
  coord_cartesian(ylim = c(0.2, 0.9)) +
  labs(x = NULL, y = NULL)+
  facet_wrap(. ~Key)+
  common_theme

# Combine the two plots and share a single legend
pp <- pp_FDP + pp_Power + plot_layout(guides = "collect") & theme(legend.position = "right")

# Display the combined plot
pp
ggsave(filename = sprintf("%s/dependent3.pdf", fig.dir), plot = pp,
       dpi = 300, device=NULL, width=plot_width-1, height=plot_height)


