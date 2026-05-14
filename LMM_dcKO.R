# Title: Linear mixed-effects model (LMM) for P1 Vgf:vGluT2dCKO #
# Author: Daven Rock                                            #
# Date: 1/26/26                                                 #
# Project: Thalamus-derived glutamate is required for early     #
# specification of layer 4 neurons in the sensory cortex        #
# Description: Shiny app for analyzing log-transformed          #
# cell counts using linear mixed-effects model (lmer).          #
# Includes EMMeans, Tukey-adjusted pairwise tests, and          #
# visualization of genotype effects (WT, CKO, DCKO) across      #
# brain regions and cell populations.                           #


library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(lme4)
library(emmeans)
library(ggsignif)

# Load csv file
csv_file <- list.files(pattern = "\\.csv$")[1]
df <- read_csv(csv_file)

# Data wrangling
df_clean <- df %>%
  filter(!grepl("/", z)) %>%
  select(-any_of(c("edu", "value")))

df_pair <- df_clean %>%
  group_by(pair, genotype, region, z) %>%
  summarise(mean_n = mean(n, na.rm = TRUE), .groups = "drop")

df_pair <- df_pair %>%
  mutate(
    genotype = factor(genotype, levels = c("WT", "CKO", "DCKO")),
    pair     = factor(pair),
    region   = factor(region),
    log_n    = log(mean_n)
  )

# LMM analysis
plot_lmm <- function(data, z_value, region_value) {
  
  plot_data <- data %>%
    filter(z == z_value, region == region_value)
  
  if (nrow(plot_data) == 0) return(NULL)
  
  ### EDIT LMM MODEL HERE
  model <- lmer(log_n ~ genotype + (1 | pair), data = plot_data)
  
  # Estimated marginal means
  emm <- emmeans(model, ~ genotype)
  emm_df <- as.data.frame(emm)
  
  # Pairwise comps
  pairwise <- as.data.frame(
    pairs(emm, adjust = "tukey")
  )
  
  # P-values
  p_wt_cko <- pairwise$p.value[pairwise$contrast == "WT - CKO"]
  p_wt_dcko <- pairwise$p.value[pairwise$contrast == "WT - DCKO"]
  p_cko_dcko <- pairwise$p.value[pairwise$contrast == "CKO - DCKO"]
  
  comparisons <- list(
    c("WT", "CKO"),
    c("WT", "DCKO"),
    c("CKO", "DCKO")
  )
  
  annotations <- c(
    paste0("p = ", signif(p_wt_cko, 3)),
    paste0("p = ", signif(p_wt_dcko, 3)),
    paste0("p = ", signif(p_cko_dcko, 3))
  )
  
  # P-value spacing
  y_base <- max(plot_data$log_n, na.rm = TRUE)
  step <- 0.15  ### EDIT P-VALUE BRACKET SPACING HERE
  y_positions <- y_base + seq(step, step * 3, length.out = 3)
  
  # Plotting
  ggplot(plot_data, aes(genotype, log_n, group = pair)) +
    
    geom_line(alpha = 0.4, linewidth = 0.7) +
    
    # Raw data points
    geom_jitter(
      width = 0,
      size = 2.5,
      alpha = 0.8
    ) +
    
    # Marginal means
    geom_point(
      data = emm_df,
      aes(x = genotype, y = emmean),
      inherit.aes = FALSE,
      color = "red",
      size = 4
    ) +
    
    # 95% CI
    geom_errorbar(
      data = emm_df,
      aes(x = genotype, ymin = lower.CL, ymax = upper.CL),
      inherit.aes = FALSE,
      width = 0.2,
      color = "red",
      linewidth = 1
    ) +
    
    # Significance brackets
    geom_signif(
      comparisons = comparisons,
      annotations = annotations,
      y_position = y_positions,
      tip_length = 0.015,
      textsize = 4.5,
      size = 0.8        
    ) +
    
    # Prevent clipping of p-value brackets
    coord_cartesian(
      ylim = c(
        min(plot_data$log_n, na.rm = TRUE),
        max(y_positions) + 0.2
      )
    ) +
    
    theme_classic(base_size = 14) +
    
  ### EDIT PLOT FONT SIZES HERE
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13)
  ) +
    
    labs(
      y = "log(mean cell count)",
      x = "",
      title = paste("Cell counts for", z_value, "in", region_value),
      subtitle = "Red = LMM marginal means ± 95% CI"
    )
}

# Shiny UI
ui <- fluidPage(
  titlePanel("LMM Cell Counts Plot"),
  
  sidebarPanel(
    selectInput(
      "region",
      "Select Region:",
      choices = sort(unique(df_pair$region))
    ),
    
    selectInput(
      "z_val",
      "Select Cell Population (z):",
      choices = sort(unique(df_pair$z))
    ),
    
    hr(),
    
    textInput("plot_name", "File name:", value = ""),
    
    numericInput(
      "plot_width",
      "Width (inches):",
      value = 6,
      min = 2
    ),
    
    numericInput(
      "plot_height",
      "Height (inches):",
      value = 6,
      min = 2
    ),
    
    selectInput(
      "plot_format",
      "File format:",
      choices = c("png", "svg")
    ),
    
    actionButton("save_plot", "Save plot")
  ),
  
  mainPanel(
    plotOutput("lmmPlot")
  )
)

# Shiny server
server <- function(input, output, session) {
  
  current_plot <- reactive({
    plot_lmm(df_pair, input$z_val, input$region)
  })
  
  output$lmmPlot <- renderPlot({
    current_plot()
  })
  
  observeEvent(input$save_plot, {
    
    p <- current_plot()
    
    if (is.null(p)) return()
    
    ## Build filename
    fname <- ifelse(
      nchar(input$plot_name) > 0,
      input$plot_name,
      paste(input$region, input$z_val, sep = "_")
    )
    
    filename <- paste0(fname, ".", input$plot_format)
    filepath <- file.path(getwd(), filename)
    
    ggsave(
      filename = filepath,
      plot = p,
      width = input$plot_width,
      height = input$plot_height,
      units = "in",
      dpi = if (input$plot_format == "png") 300 else NULL
    )
    
    message("✅ Plot saved to: ", normalizePath(filepath))
  })
}

shinyApp(ui, server)