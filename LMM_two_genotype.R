# Title: Linear mixed-effects model (LMM) for 2-genotype data   #
# Author: Daven Rock                                            #
# Date: 1/26/26                                                 #
# Project: Thalamus-derived glutamate is required for early     #
# specification of layer 4 neurons in the sensory cortex        #
# Description: Shiny app that performs linear mixed-effects     #
# analysis of pair-level cell counts for 2 genotypes and plots  #
# the resulting estimated means and confidence intervals.       #

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

# LMM and plotting
theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "black", linewidth = 0.8),
      panel.background = element_rect(fill = "white"),
      axis.text        = element_text(colour = "black", size = base_size * 0.85),
      axis.title.y     = element_text(margin = margin(r = 8)),
      axis.title.x     = element_blank(),
      plot.title       = element_text(size = base_size, face = "plain", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size * 0.85, colour = "grey40", hjust = 0.5)
    )
}

plot_lmm <- function(data, z_value, region_value, g1, g2) {
  
  plot_data <- data %>%
    filter(
      z == z_value,
      region == region_value,
      genotype %in% c(g1, g2)
    ) %>%
    mutate(genotype = factor(genotype, levels = c(g1, g2)))
  
  if (nrow(plot_data) < 2) return(NULL)
  
  model <- tryCatch(
    ### EDIT LMM MODEL HERE
    lmer(log(mean_n) ~ genotype + (1 | pair), data = plot_data),
    error = function(e) NULL
  )
  if (is.null(model)) return(NULL)
  
  # Estimated marginal means
  emm_df <- as.data.frame(emmeans(model, ~ genotype))
  
  # P-values
  pw <- as.data.frame(pairs(emmeans(model, ~ genotype)))
  p_val <- pw$p.value[1]
  p_label <- paste0("p = ", signif(p_val, 3))
  
  y_range <- range(log(plot_data$mean_n), na.rm = TRUE)
  y_pad   <- diff(y_range) * 0.2
  y_sig   <- y_range[2] + y_pad * 0.3
  y_top   <- y_sig + diff(y_range) * 0.15
  
  ggplot(plot_data, aes(x = genotype, y = log(mean_n))) +
    
    geom_line(
      aes(group = pair),
      colour = "black",
      linewidth = 0.7,
      alpha = 0.4
    ) +
    
    #Raw data points
    geom_point(
      colour = "black",
      size = 2,
      alpha = 0.8
    ) +
    
    # Estimated marginal means
    geom_point(
      data = emm_df,
      aes(x = genotype, y = emmean),
      inherit.aes = FALSE,
      colour = "red",
      size = 3
    ) +
    
    # 95% CI
    geom_errorbar(
      data = emm_df,
      aes(x = genotype, ymin = lower.CL, ymax = upper.CL),
      inherit.aes = FALSE,
      colour = "red",
      width = 0.12,
      linewidth = 0.9
    ) +
    
    # P-value brackets
    geom_signif(
      comparisons = list(c(g1, g2)),
      annotations = p_label,
      y_position = y_sig,
      tip_length = 0.015,
      textsize = 3.6,
      colour = "black",
      size = 0.45
    ) +
    
    scale_y_continuous(name = "log(mean cell count)") +
    
    coord_cartesian(ylim = c(NA, y_top), clip = "off") +
    
    theme_pub(base_size = 14)
}

# Shiny UI
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  titlePanel("LMM Cell Counts: Flexible Genotype Comparison"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      selectInput("region", "Region:",
                  choices = sort(unique(df_pair$region))),
      
      selectInput("z_val", "Cell population (z):",
                  choices = sort(unique(df_pair$z))),
      
      hr(),
      
      selectInput("geno1", "Genotype 1:",
                  choices = sort(unique(df_pair$genotype))),
      
      selectInput("geno2", "Genotype 2:",
                  choices = sort(unique(df_pair$genotype))),
      
      hr(),
      
      numericInput("plot_width", "Plot width (inches):",
                   value = 6, min = 2, step = 0.5),
      
      numericInput("plot_height", "Plot height (inches):",
                   value = 4.5, min = 2, step = 0.5),
      
      selectInput("file_type", "Save format:",
                  choices = c("png", "svg")),
      
      downloadButton("dl_plot", "Save plot")
    ),
    
    mainPanel(
      plotOutput("lmmPlot", height = "480px")
    )
  )
)

# Shiny server
server <- function(input, output, session) {
  
  current_plot <- reactive({
    req(input$geno1, input$geno2)
    validate(need(input$geno1 != input$geno2, "Pick two different genotypes"))
    
    plot_lmm(df_pair,
             z_value = input$z_val,
             region_value = input$region,
             g1 = input$geno1,
             g2 = input$geno2)
  })
  
  output$lmmPlot <- renderPlot({
    current_plot()
  }, res = 150)
  
  observeEvent(input$region, {
    z_choices <- df_pair %>%
      filter(region == input$region) %>%
      pull(z) %>% unique() %>% sort()
    
    updateSelectInput(session, "z_val", choices = z_choices)
  })
  
  output$dl_plot <- downloadHandler(
    filename = function() {
      paste0("lmm_", input$region, "_", input$z_val, ".", input$file_type)
    },
    
    content = function(file) {
      ggsave(
        filename = file,
        plot = current_plot(),
        width = input$plot_width,
        height = input$plot_height,
        units = "in",
        dpi = if (input$file_type == "png") 300 else NULL
      )
    }
  )
}

shinyApp(ui = ui, server = server)