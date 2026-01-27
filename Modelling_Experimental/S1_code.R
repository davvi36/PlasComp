################################################################################
###                                                                          ###
### Supplementary Code to:                                                   ###
###                                                                          ###
### CRISPR-Cas is beneficial in plasmid competition, but limited by          ###
### competitor toxin-antitoxin activity when horizontally transferred.       ###
###                                                                          ###
### David Sünderhauf, Jahn R. Rinnger, Leighton J. Payne,                    ###
### Rafael Pinilla-Redondo, William H. Gaze, Sam P. Brown, Stineke van Houte ###
###                                                                          ###
### available at:             ; correspondence to: david@sunderhauf.net      ###
###                                                                          ###
################################################################################

rm(list=ls())

library(tidyverse)
library(ggnewscale)
library(readxl)
library(flextable)
library(scales)

#this code is to be used in conjunction with the supplementary data file S1_data.xslx
#code sections for each figure can be run independently of each other


#Figure 1####
rm(list = ls())

## Equations E1 - E4; transition rates ####
E1 <- function(xr, yi) {(xr * si * ((1 - yi) + yi * fr))} #T(CT->C)
E2 <- function(xr, yi) {(sr + xr * si * yi * fi)}         #T(CT->T)

E3 <- function(xi, yr) {(xi * sr * ((1 - yr) + yr * fi))} #T(TC->C)
E4 <- function(xi, yr) {(si + xi * sr * yr * fr)}         #T(TC->T)


#Equations E5 and E6: delta.
E5 <- function(xr, yi) {(xr * si * ((1 - yi) + yi * fr)) - (sr + xr * si * yi * fi)}
E6 <- function(xi, yr) {(xi * sr * ((1 - yr) + yr * fi)) - (si + xi * sr * yr * fr)}

#E7 and E8: derivatives of E5 and E6
E7 <- function(yi) {(si * yi * (fr - fi - 1) + si)}
E8 <- function(yr) {(sr - sr * yr * (fr - fi + 1))}


##set parameters####
#baseline segregation rates; sr < si
sr <- 1
si <- 1.49

#PSK benefit to kin; fi < fr < 1
fr <- 0.2
fi <- 0.04

## x/ymin to max - n.b. only need this for axes
xi_min <- 1
xi_max <- 900
xr_min <- 1
xr_max <- 1300

yi_min <- 0
yi_max <- 0.66
yr_min <- 0
yr_max <- 0.99


#list all possible values xi and yr can adapt in offense, note that we will be plotting on a log scale - hence in three bins at different resolutions
grid_offense_1 <- expand.grid(xi = seq(1, 10, length.out = 1000), yr = seq(0, 1, length.out = 1000))
grid_offense_2 <- expand.grid(xi = seq(10, 100, length.out = 1000), yr = seq(0, 1, length.out = 1000))
grid_offense_3 <- expand.grid(xi = seq(100, 1300, length.out = 1300), yr = seq(0, 1, length.out = 1000))
grid_offense <- bind_rows(grid_offense_1, grid_offense_2, grid_offense_3) %>% mutate(Z = (E6(xi, yr)))

grid_offense <- grid_offense[!duplicated(grid_offense), ] #remove duplicated rows, this makes later things faster.
grid_offense_neg <- filter(grid_offense, Z < 0)
grid_offense_pos <- filter(grid_offense, Z > 0)

#and the same for defense
grid_defense_1 <- expand.grid(xr = seq(1, 10, length.out = 1000), yi = seq(0, 1, length.out = 1000))
grid_defense_2 <- expand.grid(xr = seq(10, 100, length.out = 1000), yi = seq(0, 1, length.out = 1000))
grid_defense_3 <- expand.grid(xr = seq(100, 1300, length.out = 1300), yi = seq(0, 1, length.out = 1000))
grid_defense <- bind_rows(grid_defense_1, grid_defense_2, grid_defense_3) %>% mutate(Z = (E5(xr, yi)))

grid_defense <- grid_defense[!duplicated(grid_defense), ] #remove duplictaed rows, this makes later things faster.
grid_defense_neg <- filter(grid_defense, Z < 0)
grid_defense_pos <- filter(grid_defense, Z > 0)


##plot panel C: offensive CRISPR-Cas####
ggplot() +
  geom_contour_filled(data = grid_offense_neg, 
                      aes(x = xi, y = yr, z = Z, fill = after_stat(level_mid))) +
  geom_blank(aes(fill = -210)) +
  scale_fill_steps2(low = '#b35806', mid = 'white', high = '#1f78b4', midpoint = 0, n.breaks = 100) +
  new_scale_fill()+
  geom_contour_filled(data = grid_offense_pos, 
                      aes(x = xi, y = yr, z = Z, fill = after_stat(level_mid))) +
  geom_blank(aes(fill = 2000)) +
  geom_function(fun = function(xi) {(-((-xi + si)/(1.16*xi)))}) +  #0 line
  scale_fill_steps2(low = '#b35806', mid = 'white', high = '#1f78b4', midpoint = 0, n.breaks = 100) +
  geom_path(aes(x = c(xi_min, xi_max), y = (1/(1-fi+fr))), linetype = 'dashed') +
  geom_path(aes(x = c(xi_min, xi_max), y = yr_max)) +
  geom_path(aes(x = xi_min, y = c(yr_min, yr_max))) +
  geom_path(aes(x = c(xi_min, si), y = yr_min)) +
  geom_path(aes(x = xi_max, y = c((1/(1-fi+fr)), yr_max))) +
  theme_classic() +
  ylim(limits = c(yr_min, yr_max)) +
  scale_x_log10(limits = c(xi_min, xi_max)) +
  annotation_logticks(sides = 'b') +
  labs(x = 'CRISPR-Cas strength xi', y = 'TA strength yr') +
  NULL

#ggsave('Figure_1c.svg', width = 5, height = 4, units = 'in')



##plot panel B: defensive CRISPR-Cas####
ggplot() +
  geom_contour_filled(data = grid_defense_neg, 
                      aes(x = xr, y = yi, z = Z, fill = after_stat(level_mid))) +
  geom_blank(aes(fill = -210)) +
  scale_fill_steps2(low = '#b35806', mid = 'white', high = '#1f78b4', midpoint = 0, n.breaks = 100) +
  new_scale_fill()+
  geom_contour_filled(data = grid_defense_pos, 
                      aes(x = xr, y = yi, z = Z, fill = after_stat(level_mid))) +
  geom_blank(aes(fill = 2000)) +
  geom_function(fun = function(xr) {(-((-si * xr + sr)/(1.2516*xr)))}) +  #0 line
  scale_fill_steps2(low = '#b35806', mid = 'white', high = '#1f78b4', midpoint = 0, n.breaks = 100) +
  geom_path(aes(x = c(xr_min, 1/0.663944), y = yi_max)) +
  geom_path(aes(x = xr_min, y = c(0.3914989, yi_max))) +
  theme_classic() +
  ylim(limits = c(yi_min, yi_max)) +
  scale_x_log10(limits = c(xr_min, xr_max)) +
  annotation_logticks(sides = 'b') +
  labs(x = 'CRISPR-Cas strength xr', y = 'TA strength yi') +
  NULL

#ggsave('Figure_1b.svg', width = 5, height = 4, units = 'in')


#Figure 2####
rm(list = ls())

raw <- read_xlsx("S1_data.xlsx", sheet = 2) 

raw$CFU <- ((raw$Count_1 * 10^raw$Dilution_1) + (raw$Count_2 * 10^raw$Dilution_2) + (raw$Count_3 * 10^raw$Dilution_3)) * (1000/15) #calculate CFU/mL by adding up all counts, normalised to dilution

cfu_df <- select(raw, Replicate, CRISPR, TA, Plate, CFU)

#quick raw data plot
ggplot() + 
  geom_point(data = cfu_df, aes(x = Plate, y = CFU)) +
  facet_grid(CRISPR ~ TA) +
  scale_y_log10()

#widen and calculate competitive ratios
cfu_wide <- pivot_wider(cfu_df, names_from = Plate, values_from = CFU)

cfu_wide$comprat_defence <- log10(cfu_wide$GT / cfu_wide$GK)
cfu_wide$comprat_offence <- log10(cfu_wide$CT / cfu_wide$CK)

comprat_df <- pivot_longer(select(cfu_wide, Replicate, CRISPR, TA, comprat_defence, comprat_offence), cols = c(comprat_defence, comprat_offence), names_to = 'CRISPR_action', values_to = 'comprat')

comprat_df$CRISPR_action <- factor(comprat_df$CRISPR_action, levels = c('comprat_defence', 'comprat_offence'), labels = c('defence', 'offence'))

#make plot df, summarise, and plot
plot_df <- comprat_df

plot_sum <- group_by(plot_df, CRISPR, TA, CRISPR_action) %>%
  summarise(mean = mean(comprat),
            sd = sd(comprat),
            se = sd(comprat) / sqrt(length(comprat))) %>%
  ungroup()

##plot Figure 2CD####

ggplot() +
  geom_pointrange(data=plot_sum,
                  aes(x=CRISPR, y = mean, ymin=mean-se, ymax=mean + se, 
                      shape = TA, fill = TA, group = TA),
                  size = 1, position = position_dodge(width = 0.5)) +
  geom_point(data = plot_df,
             aes(x = CRISPR, y = comprat, fill = TA, group = TA),
             shape = 21, size = 1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  geom_blank(data = plot_df,
             aes(x = CRISPR, y = -comprat, fill = TA, group = TA)) +
  geom_errorbar(data=plot_sum,
                aes(x=CRISPR, y = mean, ymin=mean-se, ymax=mean + se, 
                    group = TA),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_grid( ~ CRISPR_action) +
  theme_bw() +
  scale_shape_manual(values=c(22, 23)) +
  # scale_fill_manual(values = c('#1f78b4', '#b35806')) +
  #  scale_fill_manual(values = c('#1f78b4', '#97CAED')) +
  scale_fill_manual(values = c('#f0f0f0', '#636363')) + #from 3-class greys
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x='CRISPR-Cas status', y = 'competitive ratio [ pKJK5 / RP4 ]', fill = 'TA status', shape = 'TA status') +
  geom_vline(xintercept = 1.5, colour = 'light grey', alpha = 0.5) +
  theme(panel.grid.major.x = element_blank()) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  NULL

#ggsave('plascomp_2CD.svg', width = 6, height = 3, units = 'in')


##stats - t test####
hist((plot_df$comprat))


a<- group_by(plot_df, CRISPR, TA, CRISPR_action) %>%
  summarize(unadjusted_p_value = t.test(comprat, mu = 0)$p.value) %>%
  ungroup()

b <- flextable(a) %>%
  bold(~unadjusted_p_value <= (0.05 / nrow(a)), j = 'unadjusted_p_value') %>%
  set_header_labels(unadjusted_p_value = 'p value (unadjusted)') %>%
  add_footer_row(top = FALSE, colwidths = ncol(a), values = paste('Bonferonni-adjusted significance threshold α =', round(0.05 / nrow(a), digits = 4))) %>%
  autofit()
b

#save_as_image(b, 'figure2_tmatrix.png', webshot = 'webshot2')
#save_as_docx(b, path = 'figure2_tmatrix.docx')


##stats - models####
hist(filter(plot_df, CRISPR_action == 'defence')$comprat)

model_defence <- glm(comprat ~ CRISPR*TA, #Adding replicate is NOT significant
                     family = gaussian(link = 'identity'),
                     data = filter(plot_df, CRISPR_action == 'defence'))
#plot(model_defence) 

summary(model_defence)

drop1(model_defence, test = 'Chi') #n.b. this should technically be reduced to CRISPR only. Not surprising, as TA has no impact



hist(filter(plot_df, CRISPR_action == 'offence')$comprat)

model_offence <- glm(comprat ~ CRISPR*TA, #adding replicate is NOT significant
                     family = gaussian(link = 'identity'),
                     data = filter(plot_df, CRISPR_action == 'offence'))
#plot(model_offence) #better fit - these data are more normal

summary(model_offence)


#Figure 3B and Figure S6####
#note that the base data is collated in "Bioinformatics Analyses" supplementary code.

rm(list = ls())

`%>%` <- magrittr::`%>%`

stats <- readxl::read_xlsx("S1_data.xlsx", sheet =  15)

theme_heatmap <- function() {
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27/96*0.5),
    rect = ggplot2::element_rect(linewidth = 72.27/96*0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27/96*0.5),
    axis.ticks = ggplot2::element_blank(),
    legend.position = "left"
  )
}

theme_col_base <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 9, colour = "black"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      line = ggplot2::element_line(linewidth = 72.27/96*0.5),
      rect = ggplot2::element_rect(linewidth = 72.27/96*0.5),
      panel.border = ggplot2::element_rect(linewidth = 72.27/96*0.5),
    )
}

theme_col_x <- function() {
  theme_col_base() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}

theme_col_y <- function() {
  theme_col_base() + 
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}

##Plot Figure S6####

plot_contingency_scaled <- function(stats_table) {
  
  plot_01 <- stats_table %>%
    dplyr::mutate(
      label = paste0(
        scales::comma(observed_count), 
        "\n(", scales::comma(round(expected_count)), " ±", round(standard_error), ")"
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = targeted, y = has_ta)) +
    ggplot2::geom_tile(ggplot2::aes(fill = residual), colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5) +
    ggplot2::scale_fill_gradient2(
      low = "#dc6465",
      mid = "white", 
      high = "#5496ce",
      limits = c(-scale_limits, scale_limits)
    ) +
    ggplot2::labs(x = "Plasmid is targeted", y = "Plasmid has TA", fill = "Standardised residual") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
      text = ggplot2::element_text(size = 9, colour = "black"),
      axis.text = ggplot2::element_text(size = 9, colour = "black"),
      line = ggplot2::element_line(linewidth = 72.27/96*0.5),
      rect = ggplot2::element_rect(linewidth = 72.27/96*0.5),
      panel.border = ggplot2::element_rect(linewidth = 72.27/96*0.5),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "left",
      legend.title = ggplot2::element_text(angle = 90),
      legend.title.position = "left"
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  plot_02 <- stats_table %>%
    dplyr::summarise(count = sum(observed_count), .by = has_ta) %>%
    ggplot2::ggplot(ggplot2::aes(x = count, y = has_ta)) +
    ggplot2::geom_col(fill = "#e5e5e9", colour = "black", linewidth = 72.27/96*0.5) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0)), limits = c(0, 19000), n.breaks = 4) +
    ggplot2::labs(x = "Count", y = NULL) +
    theme_col_x() +
    ggplot2::coord_cartesian(clip = "off")
  
  plot_02
  
  plot_03 <- stats_table %>%
    dplyr::summarise(count = sum(observed_count), .by = targeted) %>%
    ggplot2::ggplot(ggplot2::aes(y = count, x = targeted)) +
    ggplot2::geom_col(fill = "#e5e5e9", colour = "black", linewidth = 72.27/96*0.5) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)), limits = c(0, 19000), n.breaks = 4) +
    ggplot2::labs(y = "Count", y = NULL) +
    theme_col_y() +
    ggplot2::coord_cartesian(clip = "off")
  
  plot_03
  
  layout <- "
    CC#
    AAB
    AAB
  "
  
  plot_out <- patchwork::wrap_plots(plot_01, plot_02, plot_03, design = layout)
  
  plot_out
}

scale_limits <- max(abs(stats$residual))

stats_all <- stats %>%
  dplyr::filter(group == "All")

stats_only_iv <- stats %>%
  dplyr::filter(group == "Type IV")

stats_no_iv <- stats %>%
  dplyr::filter(group == "Non-type IV")

plot_contingency_all <- stats_all %>%
  plot_contingency_scaled()

plot_contingency_only_iv <- stats_only_iv %>%
  plot_contingency_scaled()

plot_contingency_no_iv <- stats_no_iv %>%
  plot_contingency_scaled()


layout <- "
A#
BC
"

plot_out <- patchwork::wrap_plots(
  plot_contingency_all, 
  plot_contingency_only_iv,
  plot_contingency_no_iv,
  design = layout,
  guides = "collect"
)

plot_out


#ggplot2::ggsave('contingency_tables.svg', width = 7.7, height = 7, units = 'in')



##Plot Figure 3B####


x_order <- c("All", "Type IV", "Non-type IV")

stats_long <- stats %>%
  tidyr::pivot_longer(cols = c(observed_count, expected_count), names_to = "category", values_to = "n") %>%
  dplyr::mutate(category = dplyr::if_else(category == "observed_count", "Observed", "Expected")) %>%
  dplyr::mutate(standard_error = dplyr::if_else(category == "Observed", NA, standard_error)) %>%
  dplyr::mutate(direction = dplyr::if_else(category == "Observed", direction, "Expected")) %>%
  dplyr::mutate(
    lower = n - standard_error,
    upper = n + standard_error
  ) %>%
  dplyr::mutate(
    quadrant = dplyr::case_when(
      has_ta == TRUE & targeted == TRUE ~ "TA+ & targeted",
      has_ta == FALSE & targeted == TRUE ~ "TA- & targeted",
      has_ta == TRUE & targeted == FALSE ~ "TA+ & not targeted",
      has_ta == FALSE & targeted == FALSE ~ "TA- & not targeted"
    )
  )

format_p <- function(p) {
  ifelse(
    p > 0.001,
    formatC(p, format = "f", digits = 4),
    gsub("e\\-0*(\\d+)", " × 10<sup>−\\1</sup>", formatC(p, format = "e", digits = 2))
  )
}

annotations <- stats_long %>%
  dplyr::mutate(label_y = max(n) * 1.2, .by = quadrant) %>%
  dplyr::filter(category == "Observed") %>%
  dplyr::mutate(
    formatted_p = format_p(p_value),
    label = paste0(
      "Residual = ", round(residual, 2), "<br>",
      "<i>p</i> = ", formatted_p
    )
  ) %>%
  dplyr::select(group, quadrant, label, label_y)

dodge_width <- 0.9
bar_offset <- dodge_width / 2

lines <- stats_long %>%
  dplyr::mutate(y = max(n) * 1.1, .by = quadrant) %>%
  dplyr::mutate(x_pos = match(group, x_order)) %>%
  dplyr::reframe(
    x_start = unique(x_pos) - bar_offset / 2,
    x_end   = unique(x_pos) + bar_offset / 2,
    .by = c(group, quadrant, y)
  )

ticks <- lines %>%
  tidyr::pivot_longer(cols = c(x_start, x_end), names_to = "end", values_to = "x") %>%
  dplyr::mutate(
    yend = y - 0.02 * y
  )

plot_01 <- dplyr::filter(stats_long, quadrant == 'TA+ & targeted') %>%
  ggplot2::ggplot(ggplot2::aes(x = factor(group, x_order), y = n, fill = direction)) +
  ggplot2::geom_col(colour = "black", linewidth = 0.24, position = ggplot2::position_dodge(0.9)) +
  ggplot2::scale_fill_manual(values = c("#e5e5e9", "#B35806", "#B35806")) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), position = ggplot2::position_dodge(0.9), width = 0.2, linewidth = 0.24) +
  ggtext::geom_richtext(
    data = dplyr::filter(annotations, quadrant == 'TA+ & targeted'),
    ggplot2::aes(x = group, y = label_y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    label.color = NA,
    size = 5 / ggplot2::.pt,
    hjust = 0.5
  ) +
  ggplot2::geom_segment(
    data = dplyr::filter(lines, quadrant == 'TA+ & targeted'),
    ggplot2::aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::geom_segment(
    data = dplyr::filter(ticks, quadrant == 'TA+ & targeted'),
    ggplot2::aes(x = x, xend = x, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 0.24,
    color = "black"
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 9, colour = "black"),
    legend.position = "top",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27/96*0.5),
    rect = ggplot2::element_rect(linewidth = 72.27/96*0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27/96*0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    y = "Count",
    x = "Type of targeting CRISPR-Cas",
    fill = "plasmid count"
  ) +
  ggplot2::coord_cartesian(clip = "off") #+
# ggplot2::facet_wrap(~quadrant, scales = "free")

plot_01

#ggplot2::ggsave('fig3B.svg', height = 3, width = 3.75, units = 'in')

##Plot Figure 3C####
rm(list = ls())

stats <- readxl::read_xlsx("S1_data.xlsx", sheet = 16)

format_p <- function(p) {
  ifelse(
    p > 0.001,
    formatC(p, format = "f", digits = 4),
    gsub(
      "e\\-0*(\\d+)",
      " × 10<sup>-\\1</sup>",
      formatC(p, format = "e", digits = 2)
    )
  )
}

annotations <- stats |>
  dplyr::filter(category == "Observed") %>%
  dplyr::mutate(
    label_y = max(count) * 1.25,
    formatted_p = format_p(p_value),
    label = paste0(
      "Residual = ",
      round(residual, 2),
      "<br><i>p</i> = ",
      formatted_p
    )
  ) %>%
  dplyr::select(group, label, label_y)

annotations

x_order <- c("Type IV", "Non-type IV")

plot_fisher <- stats %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = factor(group, x_order),
      y = count,
      fill = category
    )
  ) +
  ggplot2::geom_col(
    colour = "black",
    linewidth = 0.24,
    position = ggplot2::position_dodge(0.9)
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Expected" = "#e5e5e9",
      "Observed" = "#ad4c19"
    )
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    position = ggplot2::position_dodge(0.9),
    width = 0.2,
    linewidth = 0.24
  ) +
  ggtext::geom_richtext(
    data = annotations,
    ggplot2::aes(x = group, y = label_y, label = label),
    inherit.aes = FALSE,
    fill = NA,
    label.color = NA,
    size = 5 / ggplot2::.pt,
    hjust = 0.5
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.1))
  ) +
  ggplot2::labs(
    y = "Count with compatible antitoxin\n(TA+ targets)",
    x = "Type of targeting CRISPR-Cas"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      size = 9,
      colour = "black"
    ),
    legend.position = "none",
    text = ggplot2::element_text(size = 9, colour = "black"),
    axis.text = ggplot2::element_text(size = 9, colour = "black"),
    line = ggplot2::element_line(linewidth = 72.27 / 96 * 0.5),
    rect = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.border = ggplot2::element_rect(linewidth = 72.27 / 96 * 0.5),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
  ) +
  ggplot2::coord_cartesian(clip = "off")


plot_fisher

#ggplot2::ggsave(plot_fisher, filename = "fig_ta_compatibility.svg", width = 91.2, height = 60, units = "mm")

#Figure S1####
rm(list=ls())

S_plotting <- read_xlsx('S1_data.xlsx', sheet = 7, col_types = c('text', 'text', 'text', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'))
S_plotting$CRISPR_status <- factor(S_plotting$CRISPR_status)
S_plotting$par_status <- factor(S_plotting$par_status)

S_sum <- group_by(S_plotting, par_status, CRISPR_status) %>%
  summarise(mean_sprop = mean(S_prop),
            sd_sprop = sd(S_prop),
            se_sprop = sd(S_prop) / sqrt(length(S_prop))) %>%
  ungroup()


##Plot Figure S1####

ggplot() +
  geom_pointrange(data = S_sum, 
                  aes(x = CRISPR_status, y = mean_sprop, ymin = ifelse(mean_sprop-se_sprop<0, 0, mean_sprop-se_sprop), ymax = mean_sprop+se_sprop,
                      shape = par_status, fill = par_status, group = par_status),
                  size = 1, position = position_dodge(width = 0.5)) +
  geom_point(data = S_plotting,
             aes(x = CRISPR_status, y = S_prop, 
                 fill = par_status),
             shape = 21, size = 1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  geom_errorbar(data=S_sum,
                aes(x=CRISPR_status, y = mean_sprop, ymin = ifelse(mean_sprop-se_sprop<0, 0, mean_sprop-se_sprop), ymax = mean_sprop+se_sprop, 
                    group = par_status),
                width = 0.07,
                position = position_dodge(width = 0.5)) +
  theme_bw() +
  scale_shape_manual(values=c(22, 23)) +
  # scale_fill_manual(values = c('#1f78b4', '#b35806')) +
  #  scale_fill_manual(values = c('#1f78b4', '#97CAED')) +
  scale_fill_manual(values = c('#f0f0f0', '#636363')) + #from 3-class greys
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x='CRISPR-Cas status', y = 'proportion of streptomycin resistant C.F.U.', fill = 'TA status', shape = 'TA status') +
  geom_vline(xintercept = 1.5, colour = 'light grey', alpha = 0.5) +
  geom_hline(yintercept=1, colour="black", linetype='dashed') +
  theme(panel.grid.major.x = element_blank()) +
  scale_y_log10(labels = scientific) +
  coord_cartesian(ylim = c(1e-02, 2)) +
  annotation_logticks(sides = 'l') +
  #guides(fill = 'none') +
  theme(legend.position = 'none') +
  NULL
#ggsave('s_prop_droplets_bw.svg', width = 3, height = 3, units = 'in')


#Figure S2####
rm(list = ls())

##read in proportional data ####
raw <- read_xlsx('S1_data.xlsx', sheet = 4)
raw$S_host <- factor(raw$S_host)
raw$Bystander <- factor(raw$Bystander)

#calculate cfu/mL
raw$LB <- 200 * (raw$LB_count / (10^raw$LB_dilution))
raw$C <- 200 * (raw$Cm_count / (10^raw$Cm_dilution))
raw$CK <- 200 * (raw$CmKm_count / (10^raw$CmKm_dilution))
raw$CT <- 200 * (raw$CmTmp_count / (10^raw$CmTmp_dilution))
raw$CTK <- 200 * (raw$CmTmpKm_count / (10^raw$CmTmpKm_dilution))
raw$G <- 200 * (raw$Gm_count / (10^raw$Gm_dilution))
raw$GK <- 200 * (raw$GmKm_count / (10^raw$GmKm_dilution))
raw$GT <- 200 * (raw$GmTmp_count / (10^raw$GmTmp_dilution))
raw$GTK <- 200 * (raw$GmTmpKm_count / (10^raw$GmTmpKm_dilution))

limit <- 200 * (1 / 10^0)

#extract only these useful ones
cfu_df <- select(raw, Sample:Bystander, LB:GTK)

cfu_long <- pivot_longer(cfu_df, LB:GTK, names_to = 'selective_plate', values_to = 'cfu_mL')
cfu_long$selective_plate <- factor(cfu_long$selective_plate, levels = c('LB', 'C', 'CK', 'CT', 'CTK', 'G', 'GK', 'GT', 'GTK'))


ggplot(cfu_long, aes(x = selective_plate, y = cfu_mL, fill = selective_plate)) + 
  geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 0, ymax = limit), fill = 'light grey', alpha = 0.2) +
  geom_boxplot() +
  geom_vline(xintercept = c(1.5, 5.5)) +
  geom_vline(xintercept = c(2.5, 3.5, 4.5, 6.5, 7.5, 8.5), linetype = 'dotted', alpha = 0.5) +
  facet_grid(S_host ~ Bystander) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  theme_bw() +
  scale_fill_manual(values = c('#f7f7f7', 
                               '#f7f7f7', '#b35806', '#1f78b4', '#69685D',
                               '#f7f7f7', '#b35806', '#1f78b4', '#69685D')) + #according to plasmid content, brown is RP4, blue is pKJK5, grey both
  theme(legend.position = 'none')
#ggsave('replate_raw_data.svg', width = 6, height = 6, units = 'in')

#make proportional dataframe
temp <- cfu_df
temp$C_prop <- temp$C / temp$LB
temp$C_RP4_prop <- temp$CK / temp$C
temp$C_pKJK5_prop <- temp$CT / temp$C
temp$C_comaint_prop <-temp$CTK / temp$C

temp$G_prop <- temp$G / temp$LB
temp$G_RP4_prop <- temp$GK / temp$G
temp$G_pKJK5_prop <- temp$GT / temp$G
temp$G_comaint_prop <- temp$GTK / temp$G

prop_df <- select(temp, Sample:Bystander, C_prop:G_comaint_prop)
prop_df$S_prop <- 1 - (prop_df$C_prop + prop_df$G_prop)

#test plot 
prop_long <- pivot_longer(prop_df, cols = c(C_prop:S_prop), names_to = 'proportion_of', values_to = 'proportion_value')
prop_long$proportion_of <- factor(prop_long$proportion_of, levels = c('C_prop', 'C_RP4_prop', 'C_pKJK5_prop', 'C_comaint_prop',
                                                                      'G_prop', 'G_RP4_prop', 'G_pKJK5_prop', 'G_comaint_prop',
                                                                      'S_prop'))
ggplot(prop_long, aes(x = proportion_of, y = proportion_value, fill = proportion_of)) +
  geom_boxplot() +
  geom_vline(xintercept = c(4.5, 8.5)) +
  geom_vline(xintercept = c(1.5, 5.5), linetype = 'dotted') +
  facet_grid(S_host ~ Bystander) +
  theme_bw() +
  scale_fill_manual(values = c('#f7f7f7', '#b35806', '#1f78b4', '#69685D',
                               '#f7f7f7', '#b35806', '#1f78b4', '#69685D',
                               '#f7f7f7')) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = 'none') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') +
  #  coord_cartesian(ylim = c(0, 1)) +
  NULL

#ggsave('replate_raw_proportions_log.svg', width = 6, height = 7, units = 'in')

#bring into more generic format and make plasmid_content dataframe
plasmid_content <- filter(prop_long, !proportion_of %in% c('C_prop', 'G_prop', 'S_prop'))
plasmid_content$pKJK5_type <- plasmid_content$S_host
plasmid_content$S_host <- NULL
plasmid_content$C_host <- NULL
plasmid_content$G_host <- NULL

plasmid_content <- separate(plasmid_content, proportion_of, into = c('host', 'plasmid', 'rubbish'),  sep = '_')
plasmid_content$rubbish <- NULL

#summarise 
plasmid_sum <- group_by(plasmid_content, Bystander, host, plasmid, pKJK5_type) %>%
  summarise(mean = mean(proportion_value),
            sd = sd(proportion_value)) %>%
  ungroup()


##read in stamp data ####
strep <- read_xlsx('S1_data.xlsx', sheet = 5, col_types = c('text', 'numeric', 'numeric', 'numeric', 'text', 'logical', 'logical', 'logical','logical', 'logical', 'logical', 'text', 'text', 'text', 'logical', 'logical', 'logical'))
strep$Sample <- factor(strep$Sample)


#first, tidy data. eliminate all rows in which LB is false
fdf <- filter(strep, LB == TRUE)

#now, check consistency. extract 'weird_df' where both C and G are positive
weird_df <- filter(fdf, C == TRUE & G == TRUE)

#exclude these as well
fdf <- filter(fdf, !(C== TRUE & G == TRUE))

#ok! now try to get a few more descriptive columns

fdf_sum <- group_by(fdf, Sample, batch, position, type, pKJK5_type, Bystander, pKJK5_present, RP4_present, empty_host_present) %>%
  summarise(across(LB:`T`, list(true = sum, false = function(x){sum(!x)})), .groups = 'drop')

#work with long format data instead
fdf_long <- pivot_longer(fdf, LB:`T`, names_to = 'selection', values_to = 'growth')

fdf_longsum <- group_by(fdf_long, Sample, batch, position, type, pKJK5_type, Bystander, pKJK5_present, RP4_present, empty_host_present, selection) %>%
  summarise(true = sum(growth),
            false = sum(!growth)) %>%
  ungroup()

#change data formats for nicer order
fdf_longsum$Sample <- factor(fdf_longsum$Sample, levels = c('1-1', '1-2', '1-3', '1-4', '1-5', '1-6', '1-7', '1-8', '1-9', '1-10',
                                                            '2-1', '2-2', '2-3', '2-4', '2-5', '2-6', '2-7', '2-8', '2-9', '2-10',
                                                            '3-1', '3-2', '3-3', '3-4', '3-5', '3-6', '3-7', '3-8'))
fdf_longsum$selection <- factor(fdf_longsum$selection, levels = c('LB', 'C', 'G', 'S', 'K', 'T'))


###make proportional data#
fdf_sum$C_prop <- fdf_sum$C_true / fdf_sum$LB_true
fdf_sum$G_prop <- fdf_sum$G_true / fdf_sum$LB_true
fdf_sum$S_host_prop <- (fdf_sum$LB_true - (fdf_sum$C_true + fdf_sum$G_true)) / fdf_sum$LB_true
fdf_sum$S_res_prop <- fdf_sum$S_true / fdf_sum$LB_true
fdf_sum$T_res_prop <- fdf_sum$T_true / fdf_sum$LB_true
fdf_sum$K_res_prop <- fdf_sum$K_true / fdf_sum$LB_true
fdf_sum$comprat_cw <- log10(fdf_sum$K_res_prop / fdf_sum$T_res_prop)


#make long 
temp_df <- select(fdf_sum, Sample:empty_host_present, LB_true, C_prop:comprat_cw)

prop_stamp_long <- pivot_longer(temp_df, C_prop:comprat_cw, names_to = 'proportion_of', values_to = 'proportion')


###calculate plasmid proportions from these data
C_df <- filter(fdf, C == TRUE)
C_df$commaint <- ifelse(C_df$K == TRUE & C_df$`T` == TRUE, TRUE, FALSE)
C_df$none <- ifelse(C_df$K == FALSE & C_df$`T` == FALSE, TRUE, FALSE)

G_df <- filter(fdf, G == TRUE)
G_df$commaint <- ifelse(G_df$K == TRUE & G_df$`T` == TRUE, TRUE, FALSE)
G_df$none <- ifelse(G_df$K == FALSE & G_df$`T` == FALSE, TRUE, FALSE)

S_df <- filter(fdf, C == FALSE & G == FALSE)
S_df$commaint <- ifelse(S_df$K == TRUE & S_df$`T` == TRUE, TRUE, FALSE)
S_df$none <- ifelse(S_df$K == FALSE & S_df$`T` == FALSE, TRUE, FALSE)


C_df_sum <- group_by(C_df, Sample, batch, position, type, pKJK5_type, Bystander, pKJK5_present, RP4_present, empty_host_present) %>%
  summarise(host = sum(C),
            RP4 = sum(K),#toggle which row is active to plot real data / limit of detection
            # RP4 = 1,      # this row for limit of detection   
            pKJK5 = sum(`T`),
            comaintenance = sum(commaint),
            empty = sum(none)) %>%
  ungroup()

G_df_sum <- group_by(G_df, Sample, batch, position, type, pKJK5_type, Bystander, pKJK5_present, RP4_present, empty_host_present) %>%
  summarise(host = sum(G),
            RP4 = sum(K),#toggle which row is active to plot real data / limit of detection
            #RP4 = 1,      # this row for limit of detection   
            pKJK5 = sum(`T`),
            comaintenance = sum(commaint),
            empty = sum(none)) %>%
  ungroup()

S_df_sum <- group_by(S_df, Sample, batch, position, type, pKJK5_type, Bystander, pKJK5_present, RP4_present, empty_host_present) %>%
  summarise(host = sum(LB),  #needs to be LB rather than S, as I don't have a direct S count in the dataframe. I previously filtered by LB == TRUE so this counts all present in the df, just what I want
            RP4 = sum(K),#toggle which row is active to plot real data / limit of detection
            #RP4 = 1,      # this row for limit of detection   
            pKJK5 = sum(`T`),
            comaintenance = sum(commaint),
            empty = sum(none)) %>%
  ungroup()

#paste dfs back together, calculate proportions, and make long (probably in several steps)
C_df_sum$pKJK5_C <- C_df_sum$pKJK5 / C_df_sum$host
C_df_sum$RP4_C <- C_df_sum$RP4 / C_df_sum$host
C_df_sum$comaint_C <- C_df_sum$comaintenance / C_df_sum$host
C_df_sum$empty_C <- C_df_sum$empty / C_df_sum$host

G_df_sum$pKJK5_G <- G_df_sum$pKJK5 / G_df_sum$host
G_df_sum$RP4_G <- G_df_sum$RP4 / G_df_sum$host
G_df_sum$comaint_G <- G_df_sum$comaintenance / G_df_sum$host
G_df_sum$empty_G <- G_df_sum$empty / G_df_sum$host

S_df_sum$pKJK5_S <- S_df_sum$pKJK5 / S_df_sum$host
S_df_sum$RP4_S <- S_df_sum$RP4 / S_df_sum$host
S_df_sum$comaint_S <- S_df_sum$comaintenance / S_df_sum$host
S_df_sum$empty_S <- S_df_sum$empty / S_df_sum$host

temp_df <- left_join(S_df_sum, C_df_sum, by = c('Sample', 'batch', 'position', 'type', 'pKJK5_type', 'Bystander'))
temp2_df <- left_join(temp_df, G_df_sum, by = c('Sample', 'batch', 'position', 'type', 'pKJK5_type', 'Bystander'))

plasmid_stamp_df <- select(temp2_df, Sample:Bystander, pKJK5_S:empty_S, pKJK5_C:empty_C, pKJK5_G:empty_G)
plasmid_stamp_df <- filter(plasmid_stamp_df, type == 'treatment')

plasmid_stamp_long <- pivot_longer(plasmid_stamp_df, cols = pKJK5_S:empty_G, names_to = 'temp', values_to = 'proportion')
plasmid_stamp_long <- separate(plasmid_stamp_long, temp, sep = '_', into = c('plasmid_content', 'host'), remove = FALSE)

#adjust factors/levels
plasmid_stamp_long$pKJK5_type <- factor(plasmid_stamp_long$pKJK5_type, levels = c('targeting', 'non-targeting'), labels = c('aphA99', 'nt2'))
plasmid_stamp_long$Bystander <- factor(plasmid_stamp_long$Bystander, levels = c('pCDF1b', 'pCDF1b_par'))
plasmid_stamp_long$plasmid_content <- factor(plasmid_stamp_long$plasmid_content, levels = c('RP4', 'pKJK5', 'comaint', 'empty'))
plasmid_stamp_long$host <- factor(plasmid_stamp_long$host, levels = c('C','G', 'S'))

plasmid_stamp_sum <- group_by(plasmid_stamp_long, pKJK5_type, Bystander, plasmid_content, host) %>%
  summarise(mean = mean(proportion, na.rm = TRUE),
            sd = sd(proportion, na.rm = TRUE)) %>%
  ungroup()

##unite dataframes####
str(plasmid_content)
str(plasmid_stamp_long)

#adjust data types
plasmid_content$host <- factor(plasmid_content$host)
plasmid_content$plasmid_content <- factor(plasmid_content$plasmid)
plasmid_content$plasmid <- NULL
plasmid_content$readout <- 'direct_selection'
plasmid_content$proportion <- plasmid_content$proportion_value
plasmid_content$proportion_value <- NULL

plasmid_stamp_long$Sample <- chartr('-', '_', as.character(plasmid_stamp_long$Sample))
plasmid_stamp_long$batch <- NULL
plasmid_stamp_long$position <- NULL
plasmid_stamp_long$type <- NULL
plasmid_stamp_long$temp <- NULL
plasmid_stamp_long$readout <- 'stamp_plating'

str(plasmid_content)
str(plasmid_stamp_long)


#combine dataframes
combined_plasmid <- bind_rows(plasmid_content, plasmid_stamp_long)
str(combined_plasmid)

combined_plasmid$plasmid_content <- factor(combined_plasmid$plasmid_content, levels = c('RP4', 'pKJK5', 'comaint', 'empty'))

combined_sum <- group_by(combined_plasmid, Bystander, host, pKJK5_type, plasmid_content, readout) %>%
  summarise(mean = mean(proportion, na.rm = TRUE),
            sd = sd(proportion, na.rm = TRUE)) %>%
  ungroup()

combined_across_readouts <- group_by(combined_plasmid, Bystander, host, pKJK5_type, plasmid_content) %>%
  summarise(mean = mean(proportion, na.rm = TRUE),
            sd = sd(proportion, na.rm = TRUE)) %>%
  ungroup()

#test plot
ggplot() +
  geom_pointrange(data = combined_sum, 
                  aes(x = plasmid_content, y = mean, ymin = ifelse(mean-sd<0, 0, mean-sd), ymax = mean+sd,
                      shape = Bystander, fill = plasmid_content),
                  position = position_dodge(width = 1),
                  size = 1) +
  geom_point(data = combined_plasmid,
             aes(x = plasmid_content, y = proportion, 
                 group = Bystander, fill = plasmid_content),
             position = position_dodge(width = 1),
             shape = 21) +
  #  geom_label(data = combined_plasmid,
  #             aes(x = plasmid_content, y = proportion, label = Sample, group = Bystander),
  #             position = position_dodge(width = 1)) +
  facet_grid(pKJK5_type ~ host*readout) +
  theme_bw() +
  scale_shape_manual(values = c(23, 22)) +
  scale_fill_manual(values = c('#b35806', '#1f78b4', '#69685D', 'white')) + #according to plasmid content, orange is RP4, blue is pKJK5, grey both, white neither
  geom_hline(yintercept = 1, linetype = 'dashed') +
  geom_vline(xintercept = c(1.5:3.5), linetype = 'dotted') +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = '...plasmid', y = 'proportion of hosts carrying...') +
  #  scale_y_log10() +
  #  annotation_logticks(sides = 'l') +
  NULL

#ggsave('plasmids_in_hosts.svg', width = 15, height = 6, units = 'in')



###calculate competitive ratio ####
combined_plasmid_wide <- pivot_wider(combined_plasmid, names_from = plasmid_content, values_from = proportion)

combined_plasmid_wide$comprat <- ifelse(is.finite(log10(combined_plasmid_wide$RP4 / combined_plasmid_wide$pKJK5)), log10(combined_plasmid_wide$RP4 / combined_plasmid_wide$pKJK5), NA) #competitive ratio: RP4 over pKJK5

#filter to correct readout method for each host
pnc_comprat <- filter(combined_plasmid_wide, 
                      host == 'C' & readout == 'direct_selection' |
                        host == 'G' & readout == 'direct_selection' |
                        host == 'S' & readout == 'stamp_plating')

pnc_compratsum <- group_by(pnc_comprat, Bystander, host, pKJK5_type) %>%
  summarise(comprat_mean = mean(comprat, na.rm = TRUE),
            comprat_sd = sd(comprat, na.rm = TRUE),
            comprat_se = sd(comprat, na.rm = TRUE) / sqrt(length(comprat))) %>%
  ungroup()


plascomp <- pnc_comprat
plascomp$par_status <- factor(case_when(plascomp$Bystander == 'pCDF1b' ~ 'ON',
                                        plascomp$Bystander == 'pCDF1b_par' ~ 'OFF'))
plascomp$CRISPR_status <- factor(case_when(plascomp$pKJK5_type == 'aphA99' ~ 'ON',
                                           plascomp$pKJK5_type == 'nt2' ~ 'OFF'))

plascomp_sum <- pnc_compratsum
plascomp_sum$par_status <- factor(case_when(plascomp_sum$Bystander == 'pCDF1b' ~ 'ON',
                                            plascomp_sum$Bystander == 'pCDF1b_par' ~ 'OFF'))
plascomp_sum$CRISPR_status <- factor(case_when(plascomp_sum$pKJK5_type == 'aphA99' ~ 'ON',
                                               plascomp_sum$pKJK5_type == 'nt2' ~ 'OFF'))

#add opacity column
plascomp_sum$N <- ifelse(plascomp_sum$host == 'S' & plascomp_sum$par_status == 'OFF', 
                         1, 
                         5)
plascomp_sum$N <- ifelse(plascomp_sum$host == 'S' & plascomp_sum$par_status == 'ON' & plascomp_sum$CRISPR_status == 'ON',
                         4,
                         plascomp_sum$N)

##plot####

p <- ggplot() +
  geom_pointrange(data=plascomp_sum,
                  aes(x=CRISPR_status, y = -comprat_mean, ymin=-(comprat_mean-comprat_se), ymax=-(comprat_mean + comprat_se), 
                      shape = par_status, fill = par_status, group = par_status, alpha = N/5),
                  size = 1, position = position_dodge(width = 0.5)) +
  geom_point(data = plascomp,
             aes(x = CRISPR_status, y = -comprat, fill = par_status, group = par_status),
             shape = 21, size = 1.5, alpha = 0.8, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5)) +
  geom_errorbar(data=plascomp_sum,
                aes(x=CRISPR_status, y = -comprat_mean, ymin=-(comprat_mean-comprat_se), ymax=-(comprat_mean + comprat_se), 
                    group = par_status),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  facet_grid( ~ factor(host, levels = c('S', 'C', 'G'))) +
  theme_bw() +
  scale_shape_manual(values=c(22, 23)) +
  # scale_fill_manual(values = c('#1f78b4', '#b35806')) +
  #  scale_fill_manual(values = c('#1f78b4', '#97CAED')) +
  scale_fill_manual(values = c('#f0f0f0', '#636363')) + #from 3-class greys
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x='CRISPR-Cas status', y = 'competitive ratio [ pKJK5 / RP4 ]', fill = 'TA status', shape = 'TA status') +
  geom_vline(xintercept = 1.5, colour = 'light grey', alpha = 0.5) +
  theme(panel.grid.major.x = element_blank()) +
  NULL
p

#ggsave('Figure S2.svg', width = 7, height = 3, units = 'in')


#Figure S3####
rm(list = ls())

raw <- read_xlsx('S1_data.xlsx', sheet = 9, col_types = c('numeric', 'numeric', 'numeric', 'text', 'text', 'text', 'text', 'numeric','numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'))

#calculate cfu/mL
raw$LB_cfu_1 <- (raw$LB_count1 / (10^-raw$LB_dilution1)) * 200
raw$LB_cfu_2 <- (raw$LB_count2 / (10^-raw$LB_dilution2)) * 200
raw$LB_cfu_3 <- (raw$LB_count3 / (10^-raw$LB_dilution3)) * 200

raw$C_cfu_1 <- (raw$C_count1 / (10^-raw$C_dilution1)) * 200
raw$C_cfu_2 <- (raw$C_count2 / (10^-raw$C_dilution2)) * 200
raw$C_cfu_3 <- (raw$C_count3 / (10^-raw$C_dilution3)) * 200

raw$CT_cfu_1 <- (raw$CT_count1 / (10^-raw$CT_dilution1)) * 200
raw$CT_cfu_2 <- (raw$CT_count2 / (10^-raw$CT_dilution2)) * 200
raw$CT_cfu_3 <- (raw$CT_count3 / (10^-raw$CT_dilution3)) * 200

raw$CK_cfu_1 <- (raw$CK_count1 / (10^-raw$CK_dilution1)) * 200
raw$CK_cfu_2 <- (raw$CK_count2 / (10^-raw$CK_dilution2)) * 200
raw$CK_cfu_3 <- (raw$CK_count3 / (10^-raw$CK_dilution3)) * 200

##n.b. this is not accurate for data from panel A, where full bead plates were used rather than droplet plates. analyse separately afterwards.

#calculate average cfu
cfu_df <- raw %>%
  transmute(Replicate, Batch, Position, Treatment_type, pKJK5_variant, Target_plasmid, Bystander_plasmid, Donor_OD_preadjust, Donor_OD_postadjust, Recipient_OD_preadjust, Recipient_OD_postadjust,
            LB_cfu = rowMeans(across(LB_cfu_1:LB_cfu_3), na.rm = TRUE),
            C_cfu = rowMeans(across(C_cfu_1:C_cfu_3), na.rm = TRUE),
            CT_cfu = rowMeans(across(CT_cfu_1:CT_cfu_3), na.rm = TRUE),
            CK_cfu = rowMeans(across(CK_cfu_1:CK_cfu_3), na.rm = TRUE))

#calculate competitive ratio in the C host
cfu_df$comprat <- log10(cfu_df$CT_cfu/cfu_df$CK_cfu)


##rough plot of comprat

ggplot() +
  geom_boxplot(data = filter(cfu_df, Treatment_type == 'treatment'), 
               aes(x = pKJK5_variant, y = comprat, colour = pKJK5_variant, fill = Target_plasmid)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_colour_manual(values = c('black', 'grey', 'red')) +
  facet_grid(Bystander_plasmid~ Target_plasmid) +
  theme_bw()


plot_df <- filter(cfu_df, Bystander_plasmid == 'NA')
plot_df$Plasmid_backbone <- case_when(plot_df$Target_plasmid %in% c('RP4') ~ 'RP4',
                                      plot_df$Target_plasmid %in% c('pHERD99', 'pHERD99_par') ~ 'pHERD99',
                                      plot_df$Target_plasmid %in% c('pOGG99', 'pOGG99_par') ~ 'pOGG99',
                                      plot_df$Target_plasmid %in% c('pSEVA251-99', 'pSEVA251-99_par') ~ 'pSEVA251-99')
plot_df$CRISPR_status <- case_when(plot_df$pKJK5_variant == 'aphA99' ~ 'ON',
                                   plot_df$pKJK5_variant == 'nt2' ~ 'OFF')
plot_df$TA_status <- case_when(plot_df$Target_plasmid %in% c('RP4', 'pHERD99_par', 'pOGG99_par', 'pSEVA251-99_par') ~ 'ON',
                               plot_df$Target_plasmid %in% c('pHERD99', 'pOGG99', 'pSEVA251-99') ~ 'OFF')

sum_df <- group_by(plot_df, Target_plasmid, Plasmid_backbone, CRISPR_status, TA_status) %>%
  summarise(mean = mean(comprat),
            sd = sd(comprat),
            se = sd(comprat) / sqrt(length(comprat))) %>%
  ungroup()

##plot panels b-e####
ggplot() +
  geom_pointrange(data = sum_df, aes(x = CRISPR_status, y = mean, ymin = mean-se, ymax = mean+se, 
                                     shape = TA_status, fill = TA_status, group = TA_status),
                  size = 1, position = position_dodge(width = 0.5)) +
  geom_point(data = plot_df, aes(x = CRISPR_status, y = comprat, fill = TA_status, group = TA_status),
             shape = 21, size = 1.5, alpha = 0.8, position = position_dodge(width = 0.5)) +
  geom_blank(data = plot_df, aes(y = -comprat)) +
  geom_errorbar(data = sum_df, aes(x = CRISPR_status, y = mean, ymin = mean-se, ymax = mean-se,
                                   group = TA_status),
                width = 0.07, position = position_dodge(width = 0.5)) +
  facet_wrap(~ factor(Plasmid_backbone, levels = c('RP4', 'pOGG99', 'pHERD99', 'pSEVA251-99')), 
             ncol = 2, scales = 'free_y') +
  theme_bw() +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_manual(values = c('#f0f0f0', '#636363')) + #from 3-class greys
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x = 'CRISPR-Cas status', y = 'competitive ratio [ pKJK5 / competitor ]', fill = 'TA status', shape = 'TA status') +
  geom_vline(xintercept = 1.5, colour = 'light grey', alpha = 0.5) +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.position = 'none') + 
  NULL
#ggsave('Figure_S3B-E.svg', width = 5, height = 4, units = 'in')


##statistical testing of panels B-E####
a <- group_by(plot_df, Plasmid_backbone, CRISPR_status, TA_status) %>%
  summarise(unadjusted_p_value = t.test(comprat, mu = 0)$p.value) %>%
  ungroup()

b <- flextable(a) %>%
  bold(~unadjusted_p_value <= (0.05/nrow(a)), j = 'unadjusted_p_value') %>%
  set_header_labels(unadjusted_p_value = 'p value (unadjusted)') %>%
  add_footer_row(top = FALSE, colwidths = ncol(a), values = paste('Bonferonni-adjusted significance threshold ⍺ =', round(0.05 / nrow(a), digits = 4))) %>%
  autofit()
b

#save_as_image(b, 'comprat_tmatrix.png', webshot = 'webshot2')
#save_as_docx(b, path = 'comprat_tmatrix.docx')

#fit a GLM overall and test inclusion of additional variables.
model2 <- glm(comprat ~ Plasmid_backbone*CRISPR_status*TA_status +
                #  Replicate + 
                #  Position,
                # Batch + 
                # Donor_OD_preadjust,# +
                #Donor_OD_postadjust +
                #Recipient_OD_preadjust +
                #Recipient_OD_postadjust,
                data = plot_df,
              family = gaussian(link = 'identity'))

summary.lm(model2)
drop1(model2, test = 'Chi') #drop replicate, batch. position is still sig cause of correllation with treatment. drop all the ODs as well

###separate models for each plasmid
modelrp4 <- glm(comprat ~ CRISPR_status, 
                data = filter(plot_df, Plasmid_backbone == 'RP4'),
                family = gaussian(link = 'identity'))
summary.lm(modelrp4)
#plot(modelrp4)

model_estimates <- as_flextable(modelrp4) %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'center', part = 'body') %>%
  align(j = 1, align = 'right', part = 'body') %>%
  autofit
model_estimates
#save_as_docx(model_estimates, path = 'model_estimates_rp4.docx', align = 'center')

modelpogg <- glm(comprat ~ CRISPR_status*TA_status, 
                 data = filter(plot_df, Plasmid_backbone == 'pOGG99'),
                 family = gaussian(link = 'identity'))
summary.lm(modelpogg)
#plot(modelpogg)

model_estimates <- as_flextable(modelpogg) %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'center', part = 'body') %>%
  align(j = 1, align = 'right', part = 'body') %>%
  autofit
model_estimates
#save_as_docx(model_estimates, path = 'model_estimates_pogg.docx', align = 'center')

temp <- TukeyHSD(aov(modelpogg), 'CRISPR_status:TA_status', ordered = TRUE)

temp2 <- as.data.frame(temp$`CRISPR_status:TA_status`)

c <- flextable(temp2 %>% rownames_to_column('comparison (CRISPR:TA status)')) %>%
  bold(~`p adj` <= (0.05), j = 'p adj') %>%
  set_header_labels(`p adj` = 'p value (adjusted)') %>%
  autofit()
c
#save_as_docx(c, path = 'modelpogg_tukey_table.docx')

modelpseva <- glm(comprat ~ CRISPR_status*TA_status, 
                  data = filter(plot_df, Plasmid_backbone == 'pSEVA251-99'),
                  family = gaussian(link = 'identity'))
summary.lm(modelpseva)
#plot(modelpseva)

model_estimates <- as_flextable(modelpseva) %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'center', part = 'body') %>%
  align(j = 1, align = 'right', part = 'body') %>%
  autofit
model_estimates
#save_as_docx(model_estimates, path = 'model_estimates_psesva.docx', align = 'center')

temp <- TukeyHSD(aov(modelpseva), 'CRISPR_status:TA_status', ordered = TRUE)

temp2 <- as.data.frame(temp$`CRISPR_status:TA_status`)

c <- flextable(temp2 %>% rownames_to_column('comparison (CRISPR:TA status)')) %>%
  bold(~`p adj` <= (0.05), j = 'p adj') %>%
  set_header_labels(`p adj` = 'p value (adjusted)') %>%
  autofit()
c
#save_as_docx(c, path = 'modelpseva_tukey_table.docx')



modelpherd <- glm(comprat ~ CRISPR_status*TA_status, 
                  data = filter(plot_df, Plasmid_backbone == 'pHERD99'),
                  family = gaussian(link = 'identity'))
summary.lm(modelpherd)
#plot(modelpherd)

model_estimates <- as_flextable(modelpherd) %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'center', part = 'body') %>%
  align(j = 1, align = 'right', part = 'body') %>%
  autofit
model_estimates
#save_as_docx(model_estimates, path = 'model_estimates_pherd.docx', align = 'center')

temp <- TukeyHSD(aov(modelpherd), 'CRISPR_status:TA_status', ordered = TRUE)

temp2 <- as.data.frame(temp$`CRISPR_status:TA_status`)

c <- flextable(temp2 %>% rownames_to_column('comparison (CRISPR:TA status)')) %>%
  bold(~`p adj` <= (0.05), j = 'p adj') %>%
  set_header_labels(`p adj` = 'p value (adjusted)') %>%
  autofit()
c
#save_as_docx(c, path = 'modelpherd_tukey_table.docx')


##structured panel A analysis####
raw_a <- filter(raw, Bystander_plasmid != 'NA')

raw_a$CRISPR_status <- case_when(raw_a$pKJK5_variant == 'aphA99' ~ 'ON',
                                 raw_a$pKJK5_variant == 'nt2' ~ 'OFF')
raw_a$TA_status <- case_when(raw_a$Bystander_plasmid =='pCDF1b' ~ 'ON',
                             raw_a$Bystander_plasmid == 'pCDF1b_par' ~ 'OFF')

raw_a$LB_cfu <- (raw_a$LB_count1 / (10^ -raw_a$LB_dilution1)) * 20
raw_a$C_cfu <- (raw_a$C_count1 / (10^ -raw_a$C_dilution1)) * 20
raw_a$CK_cfu <- (raw_a$CK_count1 / (10^ -raw_a$CK_dilution1)) * 20
raw_a$CT_cfu <- (raw_a$CT_count1 / (10^ -raw_a$CT_dilution1)) * 20


#pivot longer
raw_red <- select(raw_a, -c(12:47))

cfu_long <- pivot_longer(raw_red, cols = c(14:17), names_to = 'selective_plate', values_to = 'cfu_mL')

#plot cfu/mL
ggplot() +
  geom_boxplot(data = cfu_long, 
               aes(x = selective_plate, y = cfu_mL)) +
  theme_bw() +
  facet_grid(pKJK5_variant ~ Bystander_plasmid) +
  scale_y_log10() +
  NULL

#calculate comprat
raw_red$comprat <- log(raw_red$CT_cfu / raw_red$CK_cfu)


comprat_sum <- group_by(raw_red, CRISPR_status, TA_status) %>%
  summarise(mean_comprat = mean(comprat, na.rm = TRUE),
            sd_comprat = sd(comprat),
            se_comprat = sd(comprat) / sqrt(length(na.omit(comprat)))) %>%
  ungroup()

##plot panel A####

p <- ggplot() +
  geom_pointrange(data=comprat_sum,
                  aes(x=CRISPR_status, y = mean_comprat, ymin= mean_comprat-se_comprat, ymax= mean_comprat + se_comprat, 
                      shape = TA_status, fill = TA_status, group = TA_status),
                  size = 1, position = position_dodge(width = 0.5)) +
  geom_point(data = raw_red,
             aes(x = CRISPR_status, y = comprat, fill = TA_status, group = TA_status),
             shape = 21, size = 1.5, alpha = 0.8, position = position_dodge(width = 0.5)) +
  geom_errorbar(data=comprat_sum,
                aes(x=CRISPR_status, y = mean_comprat, ymin= mean_comprat-se_comprat, ymax= mean_comprat + se_comprat, 
                    group = TA_status),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_blank(data = raw_red,
             aes(x = CRISPR_status, y = -comprat)) +
  theme_bw() +
  scale_shape_manual(values=c(22, 23)) +
  # scale_fill_manual(values = c('#1f78b4', '#b35806')) +
  #  scale_fill_manual(values = c('#1f78b4', '#97CAED')) +
  scale_fill_manual(values = c('#f0f0f0', '#636363')) + #from 3-class greys
  geom_hline(yintercept = 0, colour = 'black', linetype = 'dashed') +
  labs(x='CRISPR-Cas status', y = 'competitive ratio [ pKJK5 / RP4 ]', fill = 'TA status', shape = 'TA status') +
  geom_vline(xintercept = 1.5, colour = 'light grey', alpha = 0.5) +
  theme(panel.grid.major.x = element_blank(), legend.position = 'none') +
  NULL
p

#ggsave('Figure S3_A.svg', height = 2.5, width = 3, units = 'in')



##statistical analyses


modelrp4bs <- glm(comprat ~ CRISPR_status*TA_status, 
                  data = raw_red,
                  family = gaussian(link = 'identity'))
summary.lm(modelrp4bs)
#plot(modelrp4bs)

model_estimates <- as_flextable(modelrp4bs) %>%
  align(align = 'center', part = 'header') %>%
  align(align = 'center', part = 'body') %>%
  align(j = 1, align = 'right', part = 'body') %>%
  autofit
model_estimates
#save_as_docx(model_estimates, path = 'model_estimates_rp4bs.docx', align = 'center')

temp <- TukeyHSD(aov(modelrp4bs), 'CRISPR_status:TA_status', ordered = TRUE)

temp2 <- as.data.frame(temp$`CRISPR_status:TA_status`)

c <- flextable(temp2 %>% rownames_to_column('comparison (CRISPR:TA status)')) %>%
  bold(~`p adj` <= (0.05), j = 'p adj') %>%
  set_header_labels(`p adj` = 'p value (adjusted)') %>%
  autofit()
c
#save_as_docx(c, path = 'modelrp4bs_tukey_table.docx')


#Figure S4####
rm(list = ls())
raw <- read_xlsx('S1_data.xlsx', sheet = 8, col_types = c('text', 'numeric', 'text', 'text'))

raw$host <- factor(raw$host, levels = c('S', 'C', 'G'))

raw_sum <- group_by(raw, CRISPR_status, host, Selective_Agent) %>%
  summarise(mean = mean(RP4_proportion, na.rm = TRUE),
            sd = sd(RP4_proportion, na.rm = TRUE),
            se = sd(RP4_proportion, na.rm = TRUE) / sqrt(length(na.omit(RP4_proportion)))) %>%
  ungroup()



ggplot(raw, aes(x=CRISPR_status, y=RP4_proportion, fill=Selective_Agent)) +
  geom_pointrange(data=raw_sum, 
                  aes(x=CRISPR_status, y=mean, ymin=ifelse(mean-se < 0, 0, mean-se), ymax=mean + se, 
                      fill=Selective_Agent, shape = Selective_Agent), 
                  size = 1, position=position_dodge(width = 0.5)) +
  geom_point(shape = 21, size=1.5, alpha=0.8, position=position_dodge(width=0.5)) +
  facet_grid(. ~ host) +
  theme_bw() +
  #theme(legend.position = 'none') +
  scale_y_log10(labels = function(x) format(x, scientific = TRUE)) +
  geom_hline(yintercept=1, colour="black", linetype='dashed') +
  scale_fill_manual(values = c('#fdb462', '#fccde5')) +
  scale_shape_manual(values=c(25, 24)) +
  annotation_logticks(sides="l") +
  coord_cartesian(ylim = c(1e-04, 2)) +
  labs(x = 'CRISPR-Cas status', y = 'proportion of hosts carrying RP4', fill = 'Selective Agent', shape = 'Selective Agent')

#ggsave('Figure_S4', height = 3, width = 8, units = 'in')


#Figure S5 and S6####
#n.b. Figure S5 code is found in the separate Supplementary Bioinformatics Analyses code 
#n.b. Figure S6 code is found above in the same section as Figure 3 code, as the underlying data are the same.