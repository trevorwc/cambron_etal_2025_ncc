# This reproduces the final figures from Cambron et al. (2025)

library(ggplot2)
library(meta)
library(dplyr)
library(cowplot)
library(tidyr)

###### Figure 2 ................................................................

custom_colors <- c("AGB" = "#117733",
                   "BGB" = "#999933",
                   "SOC" = "#582707",
                   "Ecosystem" = '#0077BB')
sm_text <- 20

f2pa <- read.csv('f_2_p_a.csv')

fig2_panela <- ggplot(f2pa, aes(x = factor(Experiment), y = dX, fill = `Carbon.Pool`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 1) +
  geom_errorbar(aes(ymin = dX - SE.diff, ymax = dX + SE.diff), position = position_dodge(width = 0.75), width = 0.5) +  
  labs(x = NULL, y = expression(paste("Effect of ", eCO[2], " (gC/",m^2, ")", sep = ""))) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) + 
  theme_bw() +
  theme(legend.position = c(0.95, 0.95), 
        legend.justification = c(1, 1),
        axis.text.x = element_text(angle = 45, size = sm_text,hjust=1), 
        axis.text.y = element_text(size = sm_text),  
        axis.title = element_text(size = sm_text),  
        plot.title = element_text(size = sm_text, hjust = 0.5),  
        legend.text = element_text(size = sm_text),  
        legend.title = element_text(size = sm_text),
        panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(l = 3, r =3, b = 1.5, t =1.5, unit = "lines"))

fig2_panela

f2pb <- read.csv('/Users/trevor/Desktop/Research/Nutrient Modulation of C Sink Review/figures/final/f_2_p_b.csv')

metagen(
  TE = yi, 
  seTE = vi,
  data = f2pb,
  studlab = Experiment,
  comb.fixed = FALSE, 
  comb.random = TRUE,
  hakn = TRUE,
  method.tau = "REML"
)

fig2_panelb <- ggplot(f2pb, aes(x = yi, y = Experiment)) +
  geom_vline(xintercept = 0, color = "black")+
  geom_vline(xintercept = 0.1114, color = '#0077BB')+
  geom_vline(xintercept = 0.0674, color = '#0077BB', linetype = "dashed") +
  geom_vline(xintercept = 0.1555, color = '#0077BB', linetype = "dashed") +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0)+
  scale_color_manual(values = custom_colors) +  # Use custom colors
  labs(x = "Log response ratio of ecosystem C", y = NULL) +
  geom_label(aes(x = 0.11, y = 15, label = "Average Response"), color = '#0077BB', fill = "white", size = 8) +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0,  size = sm_text),  
        axis.text.y = element_text(size = sm_text),  
        axis.title = element_text(size = sm_text),  
        plot.title = element_text(size = sm_text, hjust = 0.5),  
        legend.text = element_text(size = sm_text),  
        legend.title = element_text(size = sm_text),
        panel.grid.minor.x = element_blank(), 
        plot.margin = margin(r = 3, b = 1.5, t = 1.5,  unit = "lines"))+ 
  scale_y_discrete(limits = rev(levels(f2pa$Experiment)))

fig2_panelb

##### Figure 3 .................................................................

transparent <- rgb(0, 0, 0, alpha = 0)
sm_text <- 15
lg_text <- 25
hg_green <- "#69a889"
N_colors <- c("Fixation" = "#EE6677", 
              "Mineralization" = '#FD9A44',
              "Deposition" = "lightgray", 
              "Leaching" = "#F9D576")

n_col = "#3050F8"
c_col = "black"
chem_palette <- c("CN" = n_col, "C" = c_col)

get_prediction_error <- function(x, linmod, x_mean = NA, x_var = NA, use_method = "default"){
  
  if (use_method == "default"){
    # 95% CI (shown by geom_smooth())
    out <- predict(linmod,
                   newdata = data.frame(cLand = x),
                   se.fit = TRUE,
                   interval = "none",
                   type = "response")$se.fit
    
  } else if (use_method == "cox18"){
    # Calculate prediction error following Cox et al., 2018, but omitting the '1 + ...' in the square root
    sse <- sum(linmod$residuals^2)
    s_squared <- (1 / (length(linmod$residuals) - 2)) * sse  # Eq. 9
    out <- sqrt(s_squared) * sqrt( 1/length(linmod$residuals) + (x - x_mean)^2 / (length(linmod$residuals) * x_var) )
  }
  return(out)
}

predict_range_upper_cox18 <- function(x, linmod){
  y <- predict(linmod, data.frame(cLand = x)) + 
    get_prediction_error(x, 
                         linmod,
                         x_mean = mean(x),
                         x_var = coxvar(x),
                         use_method = "cox18")
  return(y)
}

predict_range_lower_cox18 <- function(x, linmod){
  y <- predict(linmod, data.frame(cLand = x)) -
    get_prediction_error(x, 
                         linmod,
                         x_mean = mean(x),
                         x_var = coxvar(x),
                         use_method = "cox18")
  return(y)
}

coxvar <- function(vec){
  vec <- vec[!is.na(vec)]
  sum((vec - mean(vec))^2)/length(vec)
}


landsink_gens <- read.csv('landsink_c.csv') %>% subset(generation != "IPCC TAR (2001)")
landsink_gens$generation<- factor(landsink_gens$generation, levels = c("CMIP5 (2010-2014)", "CMIP6 (current generation)"))

filtered_data = landsink_gens %>% filter(generation == "CMIP6 (current generation)")

mean_cn <- mean(subset(filtered_data, nutrient == "CN")$cLand)

linmod <- lm(nLand ~ cLand, filtered_data)

b0 <- coef(linmod)[1]
b1 <- coef(linmod)[2]

cumulative_N <- 5.77

Csink <- (cumulative_N - as.numeric(b0)) / as.numeric(b1)


fig3_panela <- ggplot(data = landsink_gens %>% subset(generation == "CMIP6 (current generation)"), aes(x = cLand, y = nLand, color = nutrient)) +
  geom_point(size = 3, show.legend = TRUE) +
  geom_smooth(data = filtered_data, aes(x = cLand, y = nLand), method = "lm", color = "darkgray", se = FALSE, 
              show.legend = FALSE, fullrange = TRUE) +
  geom_function(fun = function(x) predict_range_upper_cox18(x, linmod = linmod), color = "darkgray", linetype = "dashed") +
  geom_function(fun = function(x) predict_range_lower_cox18(x, linmod = linmod), color = "darkgray", linetype = "dashed") +
  geom_hline(yintercept = cumulative_N, linetype = "solid", color = hg_green, size = 1) +
  geom_vline(xintercept = Csink, linetype = "dashed", color = hg_green, size = 1) +
  scale_color_manual(values = chem_palette) +
  labs(title = NULL,
       x = expression(Delta~"Land Sink (Pg C)"),
       y = "Predicted N Required (Pg N)",
       color = "Nutrient") +
  scale_x_continuous(breaks = seq(100, 600, by = 100), limits = c(150, 600)) + 
  scale_y_continuous(breaks = c(0, 10, 20, 30), limits = c(0, NA)) + 
  theme_light()+
  theme(legend.position = c(0.3, 0.95), 
        legend.justification = c(1, 1), 
        legend.background = element_rect(color = "#ececec", linewidth = 0.5),
        legend.box.margin = margin(5, 5, 5, 5),
        axis.text.x = element_text( size = sm_text), 
        axis.text.y = element_text( size = sm_text), 
        axis.title.x = element_text( size = lg_text, margin = margin(t=15)), 
        axis.title.y = element_text( size = lg_text, margin = margin(r=15)),
        plot.title = element_text(size = sm_text, hjust = 0.5,), 
        legend.text = element_text(size = sm_text),
        legend.title = element_text(size = lg_text)) +
  guides(color = guide_legend(title = "Nutrient Cycle")) +
  annotate("text", 
           x = 210, 
           y = 6.2, 
           label = paste0("N supply (", round(cumulative_N, 2), " Pg)"), 
           size = 5, 
           color = hg_green#, 
           #fontface = "bold"
  )  +
  annotate("text", 
           x = 380, 
           y = 13.5, 
           size = 5,
           label = paste0("C Sink\n(",round(Csink,2)," Pg)\nsupported by\nN supply"),
           color = hg_green,
           #fontface = "bold",
           hjust = 0)

fig3_panela

f_3_p_b <- read.csv('/Users/trevor/Desktop/Research/Nutrient Modulation of C Sink Review/figures/final/f_3_p_b.csv')

fig3_panelb <- ggplot(f_3_p_b, aes(x = year, y = Tg, fill = Flux)) +
  geom_area() +
  scale_fill_manual(values = N_colors) +
  theme_light() +
  xlab('Year') +
  ylab('Tg N') +
  annotate('text', x = 43.75, y = -15, label = "-0.179% / ppm", size = sm_text * 0.35,) +
  annotate('text', x = 43.75, y = 275, label =  "0.123% / ppm", size = sm_text * 0.35,) +
  annotate('text', x = 43.75, y = 125, label =  "0.057% / ppm", size = sm_text * 0.35,) +
  annotate('text', x = 43.75, y = 15, label =  "Assume constant", size = sm_text * 0.35,) +
  annotate("text", 
           x = 43.75, y = 375, 
           label = bquote("Cumulative " * Delta * " N:" ~ .(round(cumulative_N, 2)) ~ "Pg"), 
           size = sm_text * 0.35) +
  theme(
    legend.position = "bottom",           
    legend.text = element_text(size = sm_text), 
    legend.title = element_text(size = lg_text), 
    axis.title = element_text(size = lg_text),
    axis.text = element_text(size = sm_text), 
    plot.title = element_text(size = lg_text, face = "bold"),
    strip.text = element_text(size = lg_text))+  
  guides(fill=guide_legend(ncol=2))

fig3_panelb

fig3_panelc <- ggplot(data = landsink_gens %>% 
                        subset(generation == "CMIP6 (current generation)"), 
                      aes(x = generation, y = cLand, color = nutrient)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.25, show.legend = FALSE) +
  geom_hline(yintercept = Csink, linetype = "dashed", color = hg_green, size = 1) +
  geom_boxplot(fill = transparent, position = position_dodge(width = 0.5), show.legend = FALSE, outlier.shape = NA) +
  scale_color_manual(values = chem_palette) +
  labs(title = NULL,
       x = "",
       y = expression(Delta~"Land Sink (Pg C)"),
       color = "Nutrient") +
  theme_light()+
  scale_x_discrete(labels = c("CMIP6 (current generation)" = "CMIP6")) +
  theme(legend.position = c(0.95, 0.95), 
        legend.justification = c(1, 1), 
        legend.background = element_rect(color = "#ececec", linewidth = 0.5), 
        legend.box.margin = margin(5, 5, 5, 5),
        axis.text.x = element_text( size = sm_text),  
        axis.text.y = element_text( size = sm_text), 
        axis.title.x = element_text( size = sm_text, margin = margin(t=15)), 
        axis.title.y = element_text( size = sm_text, margin = margin(r=15)),
        plot.title = element_text(size = sm_text, hjust = 0.5,), 
        legend.text = element_text(size = sm_text), 
        legend.title = element_text(size = sm_text),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
  guides(color = guide_legend(title = "Nutrient Cycle"))

fig3_panelc

fig3_paneld <- ggplot(data = landsink_gens %>% 
                       subset(generation == "CMIP6 (current generation)"), 
                     aes(x = generation, y = nLand, color = nutrient)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.25, show.legend = FALSE) +
  geom_hline(yintercept = cumulative_N, linetype = "solid", color = hg_green, size = 1) +
  geom_boxplot(fill = transparent, position = position_dodge(width = 0.5), show.legend = FALSE, outlier.shape = NA) +
  scale_color_manual(values = chem_palette) +
  labs(title = NULL,
       x = "",
       y = "Land Sink N Requirement (Pg N)",
       color = "Nutrient") +
  theme_light()+
  scale_x_discrete(labels = c("CMIP6 (current generation)" = "CMIP6")) +
  theme(legend.position = c(0.95, 0.95),  
        legend.justification = c(1, 1), 
        legend.background = element_rect(color = "#ececec", linewidth = 0.5), 
        legend.box.margin = margin(5, 5, 5, 5),
        axis.text.x = element_text( size = sm_text), 
        axis.text.y = element_text( size = sm_text),
        axis.title.x = element_text( size = sm_text, margin = margin(t=15)), 
        axis.title.y = element_text( size = sm_text, margin = margin(r=15)),
        plot.title = element_text(size = sm_text, hjust = 0.5,), 
        legend.text = element_text(size = sm_text),
        legend.title = element_text(size = sm_text),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  guides(color = guide_legend(title = "Nutrient Cycle"))

fig3_paneld

(Csink - mean_cn) / mean_cn



