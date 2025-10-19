# box plots of posterior means and overlaid posterior density plots

# we focus on oracle and DRDB posterior: Bsparse, Bridge and BART 
#CHANGED 10/06/2025: change the color darkseagreen2 to deepskyblue2 

# size of the Rplot 4.84 x 5.26 when we are saving
remove(list = ls())

# Load necessary libraries
library(ggplot2); library(reshape2) ;library(patchwork) # For combining plots
library(latex2exp); library(tidyr)

psamples = readRDS("posterior_for_mu/SimPlots/ATE_psamples_oracle_normal_n1000p10s3d2Kfold5model_type_L-Lniter20alpha1_5alpha2_3coeff_0.528270543795374_2025-09_28_03-PM.RDS")
niter = length(psamples)

# obtaining x and y values of the density plots
x1 = matrix(rep(0, niter*512), nrow = niter)
y1 = matrix(rep(0, niter*512), nrow = niter)
x2 = matrix(rep(0, niter*512), nrow = niter)
y2 = matrix(rep(0, niter*512), nrow = niter)
x3 = matrix(rep(0, niter*512), nrow = niter)
y3 = matrix(rep(0, niter*512), nrow = niter)
xorac = matrix(rep(0, niter*512), nrow = niter) 
yorac = matrix(rep(0, niter*512), nrow = niter)

for (i in 1:niter) {
  denB = density(psamples[[i]]$mu_DR_barts_psamples)
  x1[i, ] = denB$x
  y1[i, ] = denB$y
  denS = density(psamples[[i]]$mu_DR_nlp_psamples)
  x3[i, ] = denS$x
  y3[i, ] = denS$y
  denR = density(psamples[[i]]$mu_DR_brid_psamples)
  x2[i, ] = denR$x
  y2[i, ] = denR$y
  den = density(psamples[[i]]$mu_DR_oracle_psamples)
  xorac[i, ] = den$x
  yorac[i, ] = den$y
}
# generate data frame to store the values for DRDB
df_xy <- data.frame(x = c(t(x1), t(x2), t(x3)), y = c(t(y1), t(y2), t(y3)), 
                    sim = rep(rep(1:niter, each = 512), 3),
                    method = rep(c("DRDB-B","DRDB-R","DRDB-S"), each = 512 * niter))

# generate the corresponding oracle data frame 
df3 = data.frame(x = c(t(xorac)), y = c(t(yorac)), 
                 sim = rep(rep(1:niter, each = 512), 3),
                 method = rep(c("DRDB-B","DRDB-R","DRDB-S"), each = 512 * niter))


# posterior mean box plots
out = readRDS("posterior_for_mu/SimPlots/ATE_n1000p10s3d2model_type=L-L_ratio3niter500alpha1_5alpha2_3coeff_0.447213595499958_2025-05_08_12-AM.RDS")
niter = length(out)
mu_true = out[[1]]$mu_true

DR_oracle = matrix(rep(0, niter*5), nrow = niter)
DR_nlp = matrix(rep(0, niter*5), nrow = niter)
DR_bridge = matrix(rep(0, niter*5), nrow = niter)
DR_bart = matrix(rep(0, niter*5), nrow = niter)

for (i in 1:niter){
  DR_oracle[i, ] = out[[i]]$mu_oracle
  DR_nlp[i, ] = out[[i]]$mu_DR_nlp
  DR_bridge[i,] = out[[i]]$mu_DR_brid
  DR_bart[i, ] = out[[i]]$mu_DR_barts
}

# posterior density plots

d = ggplot() +
  geom_line(aes(x = x, y = y, color = "Oracle", group = sim), data = df3) +
  geom_line(aes(x = x, y = y, color = "DRDB", group = sim), data = df_xy) +
  geom_vline(xintercept = mu_true, linetype = "dashed", color = "red") + # Add vertical line at x = 5
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(
    values = c("DRDB" = "deepskyblue2", "Oracle" = "orchid2"),
    name = ""
  ) +
  facet_wrap(.~method) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(),
    legend.text = element_text(size = 10),
    legend.position = "bottom", # Position legend at the top
    legend.title = element_blank(), # Remove legend title
    legend.box.margin = margin(0, 0, 0, 0), # Adjust margin to move legend closer to plot
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)) +
  labs(title = "Posterior curves overlaid (20 replicates)", color = "", x = TeX("$\\Delta$"), y = NULL) # Remove x and y axis labels
d


# Let's add 4 boxplots of the posterior means together #####
# Combine vectors into a data frame
df_pmean <- data.frame(
  `A_1` = DR_oracle[, 1],
  `A_2` = DR_bart[, 1],
  `A_3` = DR_bridge[, 1],
  `A_4` = DR_nlp[, 1]
)

# Reshape data into long format for ggplot
df_long <- pivot_longer(df_pmean, cols = everything(), names_to = "Vectors", values_to = "Values")

# Create boxplots using ggplot with filled colors, specific names, and updated x-axis labels
pm = ggplot(df_long, aes(x = Vectors, y = Values)) +
  geom_boxplot(aes(fill = Vectors), show.legend = FALSE) + # Fill the inside of the boxplots with colors
  scale_fill_manual(values = c("A_1" = "orchid2", "A_2" = "deepskyblue2", "A_3" = "deepskyblue2", "A_4" = "deepskyblue2"), guide = "none") + # Custom colors for each boxplot
  scale_x_discrete(labels = c(
    expression(widehat(Delta)[Oracle]),
    #expression(widehat(theta)[DRDB-B][","][s]),
    expression(widehat(Delta)[DRDB-B]),
    expression(widehat(Delta)[DRDB-R]),
    expression(widehat(Delta)[DRDB-S])
  )) + # Set x-axis labels with math expressions
  geom_hline(aes(yintercept = mu_true, color = "True"), linetype = "dashed", size = 0.7, show.legend = TRUE) + 
  scale_color_manual(
    name = expression(""),
    labels = expression(paste("True ", Delta^"\u2020")),  # expression("True " * Delta[0]),
    values = c("True" = "red")) +
  labs(title = "Box plot of posterior means (500 replicates)", color = "") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 1.13), # Adjust position as needed
    legend.key.size = unit(0.15, 'cm'), # Adjust legend key size
    legend.text = element_text(size = 9), # Adjust legend text size
    legend.title = element_text(size = 1), # Adjust legend title size
    legend.background = element_rect(fill = alpha('white', 1), color = 'grey60', linewidth = 0.1), #color = NA),
    axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 8)) +
  labs(title = "Box plot of posterior means (500 replicates)", color = "", x = "", y = "")
pm

# Step 3: Stack the plots on top of each other
combined_plot <- pm / d +
  plot_layout(guides = "auto", heights = c(2, 3)) 
# + # Use patchwork to stack plots vertically
#   plot_annotation(theme = theme(legend.position = "bottom"))
# Display the combined plot
print(combined_plot) 
# library(ggpubr)
# ggarrange(pm,d, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
out[[1]]$input
