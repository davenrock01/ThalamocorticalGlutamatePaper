# Author: Daven Rock
# Title: dCKO t-test Shiny app (log ratio, variability, outliers, genotype comparison)
# Date: 5/14/26
# Description: This Shiny application analyzes double cKO cell count data using
# log-ratio normalization to WT, variability metrics, outlier detection, and
# DCKO vs CKO comparisons.

library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

# Load csv file and data wrangle
files <- list.files(pattern = "\\.csv$")
all_data <- lapply(files, read.csv)
data <- bind_rows(all_data, .id = "source") %>%
  select(file, genotype, pair, region, age, z, n)

pairwise <- data %>%
  group_by(pair, genotype, z, region) %>%
  summarise(mean_n = mean(n, na.rm = TRUE), .groups = "drop")

relative_data <- pairwise %>%
  group_by(pair, z, region) %>%
  mutate(
    WT_mean = mean_n[genotype == "WT"],
    rel_to_WT = mean_n / WT_mean,
    log_rel = log(rel_to_WT)
  ) %>%
  ungroup()

# Generate fine 0.1-step tick labels for log plots
log_breaks <- seq(0.5, 2.0, by = 0.1)

# Shiny UI
ui <- fluidPage(
  titlePanel("Cell Counts: Log-Ratio, Variability, Outliers, and DCKO vs CKO"),
  
  tabsetPanel(
    
    #WT-normalized log plot
    tabPanel("Log-Ratio Plot (vs WT)",
             sidebarLayout(
               sidebarPanel(
                 checkboxGroupInput("z_select", "Select Markers (z):",
                                    choices = sort(unique(relative_data$z)),
                                    selected = unique(relative_data$z)[1]),
                 checkboxGroupInput("region_select", "Select Regions:",
                                    choices = sort(unique(relative_data$region)),
                                    selected = sort(unique(relative_data$region))),
                 
                 # PNG export controls at bottom of sidebar
                 hr(),
                 h4("Save this log plot (PNG)"),
                 numericInput("wt_png_w", "Width (px):", value = 1600, min = 200),
                 numericInput("wt_png_h", "Height (px):", value = 1000, min = 200),
                 textInput("wt_png_name", "File name (no extension):", value = "wt_log_plot"),
                 downloadButton("save_wt_png", "Save PNG")
               ),
               mainPanel(
                 plotOutput("cellPlot", height = "600px")
               )
             )),
    
    # Variability plot
    tabPanel("Variability Plot",
             sidebarLayout(
               sidebarPanel(
                 selectInput("z_var", "Select Marker (z):",
                             choices = sort(unique(data$z)),
                             selected = unique(data$z)[1]),
                 checkboxGroupInput("region_var", "Select Regions:",
                                    choices = sort(unique(data$region)),
                                    selected = sort(unique(data$region)))
               ),
               mainPanel(plotOutput("varPlot", height = "600px"))
             )),
    
    # Outliers
    tabPanel("Outliers",
             DTOutput("outlierTable")
    ),
    
    # DCKO vs CKO plot
    tabPanel("DCKO vs CKO (single-bar log plot)",
             sidebarLayout(
               sidebarPanel(
                 checkboxGroupInput("z_dvck", "Select Markers (z):",
                                    choices = sort(unique(relative_data$z)),
                                    selected = unique(relative_data$z)[1]),
                 checkboxGroupInput("region_dvck", "Select Regions:",
                                    choices = sort(unique(relative_data$region)),
                                    selected = sort(unique(relative_data$region))),
                 
                 # PNG export controls at bottom of sidebar
                 hr(),
                 h4("Save this log plot (PNG)"),
                 numericInput("cd_png_w", "Width (px):", value = 1600, min = 200),
                 numericInput("cd_png_h", "Height (px):", value = 1000, min = 200),
                 textInput("cd_png_name", "File name (no extension):", value = "dcko_vs_cko_log_plot"),
                 downloadButton("save_cd_png", "Save PNG")
               ),
               mainPanel(
                 plotOutput("dckoCkoPlot", height = "600px")
               )
             ))
  )
)

server <- function(input, output, session) {
  
  # Log ratio plot vs. WT
  filtered_data <- reactive({
    relative_data %>%
      filter(z %in% input$z_select,
             region %in% input$region_select)
  })
  
  # P-values
  pval_table <- reactive({
    df <- filtered_data()
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      filter(genotype %in% c("CKO", "DCKO")) %>%
      group_by(z, region, genotype) %>%
      summarise(
        p_raw = tryCatch({
          x <- log_rel
          if (sum(!is.na(x)) < 2) return(NA_real_)
          t.test(x, mu = 0)$p.value
        }, error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(
        pval = case_when(
          is.na(p_raw) ~ NA_character_,
          TRUE ~ {
            out <- format(signif(p_raw, 2), scientific = FALSE, trim = TRUE)
            out <- sub("\\.?0+$", "", out)
            out
          }
        )
      )
  })
  
  # Plotting
  logPlotWT <- reactive({
    df_plot <- filtered_data() %>%
      filter(genotype %in% c("CKO", "DCKO")) %>%
      group_by(genotype, z, region) %>%
      summarise(
        mean_raw = mean(rel_to_WT, na.rm = TRUE),
        se_raw   = sd(rel_to_WT, na.rm = TRUE)/sqrt(sum(!is.na(rel_to_WT))),
        .groups = "drop"
      )
    
    df_points <- filtered_data() %>%
      filter(genotype %in% c("CKO", "DCKO"))
    
    pvals <- pval_table() %>%
      mutate(label = ifelse(is.na(pval), "", paste0("p = ", pval)))
    
    plot_labels <- df_plot %>%
      left_join(pvals, by = c("z", "region", "genotype")) %>%
      mutate(y_pos = mean_raw + se_raw * 1.5)
    
    y_breaks <- seq(0.1, 4, by = 0.1)
    
    ggplot(df_plot, aes(x = z, y = mean_raw, fill = genotype, colour = genotype)) +
      geom_col(position = position_dodge(0.9), width = 0.8, alpha = 0.6, color = "black") +
      geom_errorbar(data = df_plot,
                    aes(ymin = mean_raw - se_raw, ymax = mean_raw + se_raw),
                    position = position_dodge(0.9), width = 0.1) +
      geom_jitter(data = df_points,
                  aes(x = z, y = rel_to_WT, colour = genotype),
                  size = 3, alpha = 0.6,
                  position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
      geom_text(data = plot_labels,
                aes(x = z, y = y_pos, group = genotype, label = label),
                position = position_dodge(0.9),
                vjust = 0, size = 6, color = "black", fontface = "bold") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
      facet_wrap(~region, nrow = 1, strip.position = "top") +
      scale_fill_manual(values = c("CKO" = alpha("orange1", 0.25),
                                   "DCKO" = alpha("royalblue3", 0.25))) +
      scale_colour_manual(values = c("CKO" = "orange1",
                                     "DCKO" = "royalblue3")) +
      labs(y = "Cell Count (Relative to WT)", x = NULL,
           fill = "Genotype", colour = "Genotype") +
      scale_y_continuous(
        trans = "log10",
        breaks = y_breaks,
        labels = y_breaks
      ) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            legend.position = "right") +
      coord_cartesian(clip = "off")
  })
  
  output$cellPlot <- renderPlot({
    req(nrow(filtered_data()) > 0)
    print(logPlotWT())
  })
  
  # Download plots
  output$save_wt_png <- downloadHandler(
    filename = function() {
      fname <- input$wt_png_name
      if (nzchar(fname)) paste0(fname, ".png") else "wt_log_plot.png"
    },
    content = function(file) {
      plot_obj <- logPlotWT()
      # create png at requested pixel resolution
      png(filename = file, width = input$wt_png_w, height = input$wt_png_h, res = 150)
      print(plot_obj)
      dev.off()
    }
  )
  
  # Variability plot
  filtered_var <- reactive({
    df <- data %>%
      filter(z == input$z_var,
             region %in% input$region_var) %>%
      mutate(
        genotype = factor(genotype, levels = c("WT", "CKO", "DCKO")),
        region_pair = factor(paste0(region, "_", pair),
                             levels = unique(paste0(region, "_", pair)))
      )
    df
  })
  
  output$varPlot <- renderPlot({
    df <- filtered_var()
    req(nrow(df) > 0)
    
    summary_df <- df %>%
      group_by(region_pair, genotype) %>%
      summarise(
        mean_n = mean(n, na.rm = TRUE),
        se_n   = sd(n, na.rm = TRUE)/sqrt(n()),
        .groups = "drop"
      )
    
    ggplot() +
      geom_errorbar(
        data = summary_df,
        aes(x = region_pair, ymin = mean_n - se_n, ymax = mean_n + se_n, group = genotype),
        width = 0.2,
        position = position_dodge(width = 0.7),
        color = "black"
      ) +
      geom_jitter(
        data = df,
        aes(x = region_pair, y = n, colour = genotype),
        size = 3,
        alpha = 0.8,
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)
      ) +
      scale_colour_manual(values = c("WT" = "gray40", "CKO" = "orange1", "DCKO" = "royalblue3")) +
      labs(y = "Raw Cell Count", x = "Region_Pair", colour = "Genotype") +
      theme_minimal(base_size = 15) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Outlier table
  outlier_data <- reactive({
    data %>%
      group_by(genotype, region, pair, z) %>%
      mutate(
        Q1 = quantile(n, 0.25, na.rm = TRUE),
        Q3 = quantile(n, 0.75, na.rm = TRUE),
        IQR_val = Q3 - Q1,
        outlier = (n < Q1 - 3*IQR_val) | (n > Q3 + 3*IQR_val)
      ) %>%
      filter(outlier) %>%
      select(file, genotype, region, pair, age, z, n, Q1, Q3, IQR_val) %>%
      arrange(region, pair, genotype, z)
  })
  
  output$outlierTable <- renderDT({
    datatable(
      outlier_data(),
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE)
    )
  })
  
  # DCKO vs CKO pairwise comparison
  dcko_vs_cko_data <- reactive({
    df <- pairwise %>%
      filter(z %in% input$z_dvck, region %in% input$region_dvck) %>%
      group_by(pair, z, region) %>%
      summarise(
        CKO_mean  = mean(mean_n[genotype == "CKO"], na.rm = TRUE),
        DCKO_mean = mean(mean_n[genotype == "DCKO"], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        rel = ifelse(is.na(CKO_mean) | CKO_mean == 0, NA_real_, DCKO_mean / CKO_mean),
        log_rel = ifelse(is.na(rel), NA_real_, log(rel))
      ) %>%
      filter(!is.na(rel))
    
    df
  })
  
  dvc_pvals <- reactive({
    df <- dcko_vs_cko_data()
    if (nrow(df) == 0) return(tibble(z = character(), region = character(), pval = character()))
    
    df %>%
      group_by(z, region) %>%
      summarise(
        n_pairs = n(),
        p_raw = tryCatch({
          x <- log_rel
          if (sum(!is.na(x)) < 2) return(NA_real_)
          t.test(x, mu = 0)$p.value
        }, error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      mutate(
        pval = case_when(
          is.na(p_raw) ~ NA_character_,
          TRUE ~ {
            out <- format(signif(p_raw, 2), scientific = FALSE, trim = TRUE)
            out <- sub("\\.?0+$", "", out)
            out
          }
        )
      ) %>%
      select(z, region, pval, n_pairs)
  })
  
  dckoPlot <- reactive({
    df <- dcko_vs_cko_data()
    req(nrow(df) > 0)
    
    summary_df <- df %>%
      group_by(z, region) %>%
      summarise(
        mean_rel = mean(rel, na.rm = TRUE),
        se_rel   = sd(rel, na.rm = TRUE) / sqrt(sum(!is.na(rel))),
        n_pairs  = n(),
        .groups = "drop"
      )
    
    ptab <- dvc_pvals() %>%
      mutate(label = ifelse(is.na(pval), "", paste0("p = ", pval)))
    
    plot_labels <- summary_df %>%
      left_join(ptab, by = c("z", "region")) %>%
      mutate(y_pos = mean_rel + se_rel * 1.5)
    
    y_breaks <- seq(0.1, 4, by = 0.1)
    
    ggplot(summary_df, aes(x = z, y = mean_rel)) +
      geom_col(aes(fill = "DCKO_vs_CKO"), position = position_dodge(0.9),
               width = 0.8, alpha = 0.35, color = "orchid4") +
      geom_errorbar(aes(ymin = mean_rel - se_rel, ymax = mean_rel + se_rel),
                    position = position_dodge(0.9), width = 0.1, color = "black") +
      geom_jitter(data = df, aes(x = z, y = rel),
                  size = 3, color = "orchid4",
                  position = position_jitter(width = 0.08, height = 0)) +
      geom_text(data = plot_labels,
                aes(x = z, y = y_pos, label = label),
                vjust = 0, size = 6, color = "black", fontface = "bold") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
      facet_wrap(~region, nrow = 1, strip.position = "top") +
      scale_fill_manual(name = NULL, values = c("DCKO_vs_CKO" = alpha("orchid1", 0.4))) +
      labs(y = "Rel. Cell Count (DCKO/CKO)", x = NULL) +
      scale_y_continuous(trans = "log10",
                         breaks = y_breaks,
                         labels = y_breaks) +
      theme(text = element_text(size = 15),
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            legend.position = "none") +
      coord_cartesian(clip = "off")
  })
  
  output$dckoCkoPlot <- renderPlot({
    print(dckoPlot())
  })
  
  # Download plot
  output$save_cd_png <- downloadHandler(
    filename = function() {
      fname <- input$cd_png_name
      if (nzchar(fname)) paste0(fname, ".png") else "dcko_vs_cko_log_plot.png"
    },
    content = function(file) {
      plot_obj <- dckoPlot()
      png(filename = file, width = input$cd_png_w, height = input$cd_png_h, res = 150)
      print(plot_obj)
      dev.off()
    }
  )
  
}

shinyApp(ui, server)
