library(shiny)
library(shinythemes)
library(tidymodels)
library(protr)
library(Biostrings)
library(DT)
library(tidyverse)
library(Peptides)
library(plotly)
library(htmlwidgets)
library(cluster)     
library(factoextra)  

options(shiny.maxRequestSize = 50 * 1024^2)

# ============================================
# 1. MODEL LOADING
# ============================================
models_bundle <- tryCatch({
  readRDS("models/Ensemble_Bundle.rds") 
}, error = function(e) return(NULL))

# ============================================
# 2. HELPER FUNCTIONS
# ============================================

# Clean sequences: remove spaces, newlines, tabs, asterisks
# and keep only standard amino acids (A-Z, no B, J, O, U, X, Z)
clean_sequences <- function(seqs) {
  seqs <- toupper(gsub("[ \n\r\t*]", "", seqs))
  standard_aa <- "^[ACDEFGHIKLMNPQRSTVWY]+$"
  valid_idx <- grepl(standard_aa, seqs)
  return(seqs[valid_idx])
}

# Calculate cysteine percentage in a sequence
calc_cys_pct <- function(seq) {
  s <- strsplit(seq, "")[[1]]
  round((sum(s == "C") / length(s)) * 100, 2)
}

# Extract features from sequences for model prediction
# Uses sliding windows (N-terminal, middle, C-terminal) and global descriptors
extract_features_app <- function(seqs) {
  if(length(seqs) == 0) return(NULL)
  
  # 1. Padding to avoid errors with short sequences
  x_padded <- str_pad(seqs, width = 100, side = "right", pad = "A")
  
  # 2. N-terminal window definitions
  z1 <- substr(x_padded, 1, 30)
  z2 <- substr(x_padded, 31, 60)
  z3 <- substr(x_padded, 61, 100)
  
  # 3. AAC extraction per window
  # Force as.numeric() to prevent colMeans errors in downstream steps
  f_z1 <- as.data.frame(t(sapply(z1, function(x) as.numeric(extractAAC(x)))))
  colnames(f_z1) <- paste0("Z1_", names(extractAAC("ACDEF")))
  
  f_z2 <- as.data.frame(t(sapply(z2, function(x) as.numeric(extractAAC(x)))))
  colnames(f_z2) <- paste0("Z2_", names(extractAAC("ACDEF")))
  
  f_z3 <- as.data.frame(t(sapply(z3, function(x) as.numeric(extractAAC(x)))))
  colnames(f_z3) <- paste0("Z3_", names(extractAAC("ACDEF")))
  
  # 4. Global descriptors
  f_gl_c <- as.data.frame(t(sapply(seqs, function(x) as.numeric(extractCTDC(x)))))
  colnames(f_gl_c) <- paste0("GLC_", names(extractCTDC("ACDEF")))
  
  f_gl_t <- as.data.frame(t(sapply(seqs, function(x) as.numeric(extractCTDT(x)))))
  colnames(f_gl_t) <- paste0("GLT_", names(extractCTDT("ACDEF")))
  
  # Final combination
  all_features <- cbind(f_z1, f_z2, f_z3, f_gl_c, f_gl_t)
  
  # Clean column names to match tidymodels (replace ":" or " " with ".")
  colnames(all_features) <- make.names(colnames(all_features))
  
  # Ensure everything is strictly numeric
  all_features[] <- lapply(all_features, function(x) as.numeric(as.character(x)))
  
  return(all_features)
}  

# ============================================
# 3. UI
# ============================================
ui <- fluidPage(
  div(style = "text-align: center;",
      img(src = "ExoBanner.png", width = "45%", height = "auto")
       ),
  theme = shinytheme("cosmo"),
  titlePanel(span(strong("ExoproteoML"), 
                  ": Ensemble Learning Framework for Eukaryotic Exoproteome Prediction")),
  
  sidebarLayout(
    sidebarPanel(
      h4("Input Data"),
      fileInput("fasta_in", "Upload FASTA file", accept = c(".fasta", ".fa")),
      textAreaInput("manual_in", "Or Paste Sequences", 
                    placeholder = ">ID\nSEQUENCE...", height = "100px"),
      hr(),
      actionButton("predict_btn", "Run Full Pipeline", class = "btn-primary btn-block"),
      hr(),
      h4("Global Export"),
      downloadButton("download_tsv", "Full Report (TSV - Safe Format)", class = "btn-info btn-block"),
      helpText("Note: TSV prevents ID commas from breaking columns.")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Full Results Table", br(), DTOutput("results_table")),
        
        tabPanel("Global Analytics", br(),
                 fluidRow(
                   column(6, 
                          h4("CAZy Potential"), 
                          plotlyOutput("cazy_bar"),
                          br(),
                          downloadButton("down_cazy_bar_html", "HTML", class="btn-xs"),
                          downloadButton("down_cazy_bar_svg", "SVG", class="btn-xs")),
                   column(6, 
                          h4("SSP/Effector Distribution"), 
                          plotlyOutput("ssp_bar"),
                          br(),
                          downloadButton("down_ssp_bar_html", "HTML", class="btn-xs"),
                          downloadButton("down_ssp_bar_svg", "SVG", class="btn-xs"))
                 )),
        
        tabPanel("Bio-Discovery Space", br(),
                 h4("Enzymatic Profile"),
                 plotlyOutput("cazy_scatter", height = "400px"),
                 downloadButton("down_cazy_scat_html", "Download HTML", class="btn-xs"),
                 hr(),
                 h4("Cysteine-Rich Profile"),
                 plotlyOutput("ssp_scatter", height = "400px"),
                 downloadButton("down_ssp_scat_html", "Download HTML", class="btn-xs")),
        
        tabPanel("3D Secretome Space", br(),
                 h4("Physicochemical PCA"),
                 plotlyOutput("pca_3d_plot", height = "600px"),
                 div(align="right", downloadButton("down_pca_3d_html", "Export 3D Space (HTML)", class="btn-success")),
                 hr(),
                 h4("Cluster Summary Stats"),
                 DTOutput("cluster_summary_table"))
      )
    )
  )
)

# ============================================
# 4. SERVER
# ============================================

server <- function(input, output) {
  
  # Reactive values for plot objects (used in download handlers)
  p1 <- reactiveVal(); p2 <- reactiveVal(); p3 <- reactiveVal(); p4 <- reactiveVal(); p5 <- reactiveVal()
  
  # Main reactive: runs the full prediction pipeline
  results_data <- eventReactive(input$predict_btn, {
    if(is.null(models_bundle)) return(NULL)
    
    # --- Input parsing ---
    seqs <- NULL; ids <- NULL
    if (!is.null(input$fasta_in)) {
      raw <- readAAStringSet(input$fasta_in$datapath)
      ids <- names(raw)
      seqs <- as.character(raw)
    } else if (nchar(input$manual_in) > 5) {
      txt <- trimws(input$manual_in)
      if (substr(txt, 1, 1) == ">") {
        tmp <- tempfile()
        writeLines(txt, tmp)
        raw <- readAAStringSet(tmp)
        ids <- names(raw)
        seqs <- as.character(raw)
      } else {
        raw_list <- unlist(strsplit(txt, "\n"))
        seqs <- raw_list
        ids <- paste0("Manual_", seq_along(seqs))
      }
    }
    req(seqs)
    
    # Clean sequences and filter valid ones
    clean_flags <- grepl("^[ACDEFGHIKLMNPQRSTVWY*]+$", toupper(gsub("[ \n\r\t]", "", seqs)))
    seqs <- seqs[clean_flags]
    ids <- ids[clean_flags]
    seqs_final <- clean_sequences(seqs)
    ids_final <- ids[1:length(seqs_final)]
    
    # --- Pipeline with progress bar ---
    withProgress(message = 'Deep Discovery Pipeline...', value = 0, {
      
      # 1. Feature extraction
      incProgress(0.2)
      features <- extract_features_app(seqs_final)
      
      # 2. Run ensemble predictions
      incProgress(0.4)
      preds <- map(models_bundle, ~predict(.x, features, type = "prob"))
      
      # 3. Ensemble averaging
      prob_sp  <- rowMeans(do.call(cbind, map(preds, ~.x$.pred_Classical)))
      prob_ups <- rowMeans(do.call(cbind, map(preds, ~.x$.pred_NonClassical)))
      prob_neg <- rowMeans(do.call(cbind, map(preds, ~.x$.pred_NonSecreted)))
      
      probs_matrix <- cbind(prob_sp, prob_ups, prob_neg)
      best_idx <- apply(probs_matrix, 1, which.max)
      
      # 4. Additional calculations (CAZy, MW, etc.)
      incProgress(0.2)
      cazy_idx <- sapply(seqs_final, function(x) {
        aac <- extractAAC(x)
        round(sum(aac[c("F","W","Y","G","P")], na.rm = TRUE) * 100, 1)
      })
      
      mw_vals <- round(sapply(seqs_final, function(x) {
        aa_mass <- c(
          A = 71.08, R = 156.19, N = 114.10, D = 115.09, C = 103.14,
          E = 129.12, Q = 128.13, G = 57.05, H = 137.14, I = 113.16,
          L = 113.16, K = 128.17, M = 131.19, F = 147.18, P = 97.12,
          S = 87.08, T = 101.11, W = 186.21, Y = 163.18, V = 99.13
        )
        sum(aa_mass[strsplit(x, "")[[1]]], na.rm = TRUE) / 1000 + 0.018 
      }), 2)
      
      # 5. Build final dataframe
      final_df <- data.frame(
        seqID = ids_final,
        Veredict = case_when(
          grepl("[KH]DEL$", seqs_final) & best_idx == 1 ~ "SP+ (ER-Retained)",
          best_idx == 1 ~ "Conventional (SP+)",
          best_idx == 2 ~ "Non-Conventional (UPS)",
          TRUE ~ "Non-Secreted"
        ),
        Confidence = round(apply(probs_matrix, 1, max) * 100, 2),
        Prob_SP  = round(prob_sp, 4),
        Prob_UPS = round(prob_ups, 4),
        Prob_Neg = round(prob_neg, 4),
        Is_SSP = ifelse(
          nchar(seqs_final) < 250 & sapply(seqs_final, calc_cys_pct) > 3 & best_idx == 1,
          "★ Candidate SSP/Effector",
          "Standard Profile"
        ),
        Cazy_Index = cazy_idx,
        CAZy_Potential = ifelse(cazy_idx > 25, "High Potential", "Standard"),
        Length = nchar(seqs_final),
        Cys_Percent = sapply(seqs_final, calc_cys_pct),
        N_Hydro = round(sapply(seqs_final, hydrophobicity), 3),
        N_Charge = round(sapply(seqs_final, charge, pH = 7), 3),
        MW_kDa = mw_vals,
        Sequence = seqs_final
      )
      
      # --- PCA and Clustering ---
      pca_data <- final_df %>% 
        select(MW_kDa, Cazy_Index, Confidence) %>% 
        mutate(across(everything(), ~as.numeric(as.character(.)))) %>%
        na.omit()
      
      # Only compute PCA if at least 2 rows and columns are not constant
      if(nrow(pca_data) > 1) {
        pca_res <- prcomp(pca_data, scale. = TRUE)
        final_df$PC1 <- round(pca_res$x[,1], 3)
        final_df$PC2 <- round(pca_res$x[,2], 3)
        final_df$PC3 <- round(pca_res$x[,3], 3)
        
        # K-means clustering
        set.seed(42)
        km <- kmeans(scale(pca_data), centers = min(3, nrow(pca_data)))
        final_df$Cluster <- as.factor(paste0("Group_", km$cluster))
      } else {
        # Default values if PCA cannot be computed
        final_df$PC1 <- 0
        final_df$PC2 <- 0
        final_df$PC3 <- 0
        final_df$Cluster <- as.factor("Group_1")
      }
      
      final_df
    })
  })
  
  # ============================================
  # OUTPUTS & PLOTS
  # ============================================
  
  output$results_table <- renderDT({
    req(results_data())
    df_display <- results_data() %>% 
      select(-PC1, -PC2, -PC3) %>%
      as.data.frame() 
    
    datatable(df_display, options = list(scrollX = TRUE), rownames = FALSE)
  })
  
  output$cazy_bar <- renderPlotly({
    req(results_data())
    p <- ggplot(results_data(), aes(x = Veredict, fill = CAZy_Potential)) + 
      geom_bar(position = "dodge") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    plt <- ggplotly(p)
    p1(plt)
    plt
  })
  
  output$ssp_bar <- renderPlotly({
    req(results_data())
    p <- ggplot(results_data(), aes(x = Veredict, fill = Is_SSP)) + 
      geom_bar(position = "dodge") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    plt <- ggplotly(p)
    p2(plt)
    plt
  })
  
  output$cazy_scatter <- renderPlotly({
    req(results_data())
    p <- ggplot(results_data(), aes(x = MW_kDa, y = Cazy_Index, color = CAZy_Potential, text = seqID)) + 
      geom_point() + 
      theme_minimal()
    plt <- ggplotly(p, tooltip = "text")
    p3(plt)
    plt
  })
  
  output$ssp_scatter <- renderPlotly({
    req(results_data())
    p <- ggplot(results_data(), aes(x = Length, y = Cys_Percent, color = Is_SSP, text = seqID)) + 
      geom_point() + 
      theme_minimal()
    plt <- ggplotly(p, tooltip = "text")
    p4(plt)
    plt
  })
  
  output$pca_3d_plot <- renderPlotly({
    req(results_data())
    plt <- plot_ly(
      results_data(),
      x = ~PC1, y = ~PC2, z = ~PC3,
      color = ~Veredict,
      type = 'scatter3d',
      mode = 'markers',
      text = ~paste0("ID: ", seqID, "<br>Prob SP: ", Prob_SP),
      hoverinfo = 'text'
    )
    p5(plt)
    plt
  })
  
  output$cluster_summary_table <- renderDT({
    req(results_data())
    summary <- results_data() %>%
      group_by(Cluster) %>%
      summarise(
        n   = n(),
        MW  = round(mean(MW_kDa), 1),
        Cys = round(mean(Cys_Percent), 1),
        CAZy = round(mean(Cazy_Index), 1)
      )
    datatable(summary, options = list(dom = 't'), rownames = FALSE)
  })
  
  # ============================================
  # DOWNLOAD HANDLERS
  # ============================================
  
  # TSV Export (safe for commas in IDs)
  output$download_tsv <- downloadHandler(
    filename = function() { paste0("ExoproteoML_Report_", Sys.Date(), ".tsv") },
    content = function(file) {
      write.table(results_data(), file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # HTML Exports (interactive Plotly)
  output$down_cazy_bar_html <- downloadHandler(
    filename = "CAZy_Bar.html",
    content = function(file) saveWidget(p1(), file)
  )
  output$down_ssp_bar_html <- downloadHandler(
    filename = "SSP_Bar.html",
    content = function(file) saveWidget(p2(), file)
  )
  output$down_cazy_scat_html <- downloadHandler(
    filename = "CAZy_Space.html",
    content = function(file) saveWidget(p3(), file)
  )
  output$down_ssp_scat_html <- downloadHandler(
    filename = "SSP_Space.html",
    content = function(file) saveWidget(p4(), file)
  )
  output$down_pca_3d_html <- downloadHandler(
    filename = "PCA_3D.html",
    content = function(file) saveWidget(p5(), file)
  )
  
  # SVG Exports (static from Plotly)
  output$down_cazy_bar_svg <- downloadHandler(
    filename = "CAZy_Bar.svg",
    content = function(file) export(p1(), file = file)
  )
  output$down_ssp_bar_svg <- downloadHandler(
    filename = "SSP_Bar.svg",
    content = function(file) export(p2(), file = file)
  )
}

shinyApp(ui, server)
