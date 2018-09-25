##########################################################################################
### Modes:
###      VIEW: geneA
###    SEARCH: geneA similarity search
###   COMPARE: geneA vs. geneB
### SABV: include sex-as-biological-variable analysis
### Dissim: dissimilarity search instead
##########################################################################################
### GGC = gene-gene associations
### EPS = Expression Profiles
##########################################################################################
### Jeremy Yang
##########################################################################################
library(readr)
library(wCorr)
library(shiny, quietly = T)
library(DT, quietly = T) #datatable
library(dplyr, quietly = T)
library(plotly, quietly = T)

###
# This code runs once for all sessions.
###
APPNAME <- "Ex-files"
###
t0 <- proc.time()
if (file.exists("exfiles.Rdata")) { # Read Rdata if present
  load("exfiles.Rdata")
} else { # else read files and write Rdata.
  ###
  # tissue_id, tissue_name
  ###
  tissue <- read_delim("exfiles_tissue_order.tsv", "\t", col_names=c("id","name"))
  ###
  # ENSG, NCBI, HGNCID, symbol, name
  ###
  gene <- read_delim("gtex_gene_xref.tsv", "\t") # gene attributes
  ###
  # ENSG, SEX, tissue.1, tissue.2, etc.
  ###
  eps <- read_delim("exfiles_eps.tsv", "\t", col_types=cols(SEX=col_character())) # expression profiles
  ###
  # ENSGA, ENSGB, Cluster, wRho, Ruzicka
  ###
  ggc <- read_delim("exfiles_ggc.tsv.gz", "\t", col_types="cccdd")
  ###
  save(tissue, gene, eps, ggc, file="exfiles.Rdata")
}
tissue_missing <- setdiff(tissue$name, colnames(eps))
message(sprintf("NOTE: TISSUE_MISSING: %d. %s\n", 1:length(tissue_missing), tissue_missing))
tissue <- tissue[tissue$name %in% colnames(eps),]
eps <- eps[,c("ENSG","SEX",tissue$name)]
#
message(sprintf("%d. %s\n", tissue$id, tissue$name))
#
ensgs <- intersect(eps$ENSG, c(ggc$ENSGA,ggc$ENSGB))
ensgs <- intersect(ensgs, gene$ENSG)
#
ensg_dups <- gene$ENSG[duplicated(gene$ENSG)]
message(sprintf("Duplicated/ambiguous gene ID: %s\n", ensg_dups))
#
gene <- gene[!is.na(gene$symbol),]
gene <- gene[gene$ENSG %in% ensgs,]
#
symb_dups <- gene$symbol[duplicated(gene$symbol)]
message(sprintf("Duplicated/ambiguous gene SYMBs: %s\n", symb_dups))
#
gene <- gene[!duplicated(gene$ENSG),]
gene <- gene[!duplicated(gene$symbol),]
#
message(sprintf("Gene unique ENSG count: %d", length(unique(gene$ENSG))))
message(sprintf("Gene unique SYMB count: %d", length(unique(gene$symbol))))
###
eps <- merge(eps, gene[,c("ENSG","symbol")], by="ENSG", all=T)
ggc$Combo <- round(ggc$wRho*ggc$Ruzicka, digits=3)
#
t_elapsed <- (proc.time()-t0)[3]
#
db_htm <- sprintf("<B>Dataset:</B> Genes: %d ; tissues: %d ; comparisons: %d (t_load: %.1fs)", nrow(gene), nrow(tissue), nrow(ggc), t_elapsed)
###
#
qryArand <- sample(gene$symbol, 1) # initial random query
#
#############################################################################
RMSD <- function(A,B) {
  sqrt(sum((A-B)^2)/length(A))
}
Tanimoto <- function(A,B) {
  (A %*% B) / (A %*% A + B %*% B - A %*% B)
}
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  wCorr::weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
Ruzicka <- function(A,B) {
  sum(pmin(A,B))/sum(pmax(A,B))
}
###
# Create list of choices for input autocomplete.
gene_choices <- list()
for (i in 1:nrow(gene))
{
  gene_choices[[paste(gene$symbol[i], gene$name[i])]] <- gene$symbol[i]
}
#
#############################################################################
HelpHtm <- function() {(
"<P><B>Ex-files</B> allows exploration and analysis of co-expression patterns via gene expression profiles.
With <B>GTEx</B> as the data source, gene expression profiles are computed as real valued vectors of expression levels 
across the defined tissue types.</P>
<P>
<B>Inputs</B> are query geneA, and <I>optionally</I>, geneB.
Ex-files modes of operation:
<UL>
<LI><B>View</B> - view profile for one gene
<LI><B>Search</B> - search for genes based on profile similarity
<LI><B>Compare</B> - compare input genes via profiles
</UL>
<B>SABV:</B> This checkbox invokes <B><I>Sex As a Biological Variable</I></B> analysis.</P>
<P>
<B>Metrics:</B>
<UL>
<LI><B>Ruzicka</B> - Similarity measure, size-normalizing, hence advantageous over Euclidean.
<LI><B>wRho</B> - Weighted Pearson correlation coefficient, weighted by average values for each tissue, to mitigate noise.
<LI><B>Combo</B> - Product wRho*Ruzicka, balancing correlation and similarity.
</UL>
<B>Dissimilarity:</B> This checkbox affects Search results, reversing the chosen metric.</P>
<P>
<B>Results:</B>
<UL>
<LI>Search results include the top hit, displayed against the query in the
plot, and other hits displayed in a table and downloadable as CSV.  <B>Note: current dataset includes only highly
correlated or anti-correlated profiles.</B>
<LI>Compare results consist of the plot of GeneA vs. GeneB.
<LI>SABV results include both sexes for both genes, thus four profiles.
<LI>Expression units are <B>LOG<SUB>10</SUB>(1 + TPM)</B>, where <B>TPM</B> = RNA-seq median Transcripts Per Million-kilobase.
<LI>Groups denote correlations among F-only, M-only or F+M.
</UL></P>
<P>
Notes on data preparation: This version is focused on SABV knowledge discovery, thus reproductive and 
breast tissues not considered. Also we restrict to protein-encoding genes mapped to HUGO gene symbols,
for scientific comprehensibility, so only data associated with these genes is retained.
Currently also genes are ignored with mapping ambiguity between Ensembl ENSG to HUGO symbols. 
Note that the term metric is used here in the common and not mathematically rigorous sense.
</P>
<B>Algorithms:</B> Giovanni Bocci, Oleg Ursu, Cristian Bologa, Steve Mathias, Jeremy Yang &amp; Tudor Oprea<BR/>
<B>Web app:</B> Jeremy Yang<BR/>
Data from <A HREF=\"https://www.gtexportal.org/\" TARGET=\"_blank\">GTEx, The Genotype-Tissue Expression Project</A>.<BR/>
Built with R-Shiny &amp; Plotly.<BR/>
")}
#
#############################################################################
ui <- fluidPage(
  titlePanel(h2(sprintf("%s, GTEx expression-profile exploration", APPNAME), em("(BETA)")), windowTitle=APPNAME),
  fluidRow(
    column(4, 
        wellPanel(
          actionButton("randGene", "GeneA (click for random)", style='padding:4px; background-color:#DDDDDD; font-size:100%; font-weight:bold'),
          selectizeInput("qryA", label=NULL, choices = gene_choices, selected=qryArand),
          selectizeInput("qryB", label="GeneB (optional)", choices = c(list('None'='none'), gene_choices)),
          #selectizeInput("qryB", label="GeneB (optional)", choices = NULL), #server-side not working
          radioButtons("mode", "Mode", choices=c("View", "Search", "Compare"), selected="View", inline=T),
          checkboxInput("sabv", span("SABV", icon("venus",lib="font-awesome"),icon("mars", lib="font-awesome")), value=T),
          radioButtons("metric", "Metric", choices=c("Ruzicka", "wRho", "Combo"), selected="Combo", inline=T),
          checkboxGroupInput("plot_opts", "Plot", choices=c("Annotate"), selected=NULL, inline=F),
          br(),
          actionButton("showhelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
          )),
    column(8, conditionalPanel(condition="true", plotlyOutput("plot", height = "580px")))
  ),
  fluidRow(column(12, 
      wellPanel(htmlOutput(outputId = "result_htm", height = "60px")))
  ),
  
  # How to test in JS for R NULL?
  conditionalPanel(condition="(input.mode=='Search' && typeof output.datarows !== 'undefined')",
    wellPanel(
      fluidRow(column(12, DT::dataTableOutput("datarows"))),
      fluidRow(column(12, downloadButton("hits_file", label="Download")))),
    width=12),
  fluidRow(column(12, wellPanel(
      htmlOutput(outputId = "log_htm", height = "60px")))),
  fluidRow(
    column(12, em(strong(sprintf("%s", APPNAME)), " web app built with R-Shiny and Plotly, from ", tags$a(href="http://datascience.unm.edu", "UNM Translational Informatics Division")))
  )
)

#############################################################################
server <- function(input, output, session) {
  observeEvent(input$showhelp, {
    showModal(modalDialog(
      title = HTML(sprintf("<H2>%s Help</H2>", APPNAME)),
      HTML(HelpHtm()),
      easyClose = T,
      footer = tagList(modalButton("Dismiss"))
    ))
  })
  
# observe({ 
#   updateSelectizeInput(session, 'qryB', choices = c(list('None'='none'), gene_choices), server=T) #NOT WORKING
# })
  
  Sys.sleep(1)
  randGeneA_previous <- 0 # initialize once per session

  message(sprintf("NOTE: genes: %d ; correlations = %d", nrow(gene), nrow(ggc)))
  observe({
    message(sprintf("NOTE: mode: %s ; sabv: %s ; metric: %s", input$mode, input$sabv, input$metric))
    message(sprintf("NOTE: qryA = %s \"%s\"", qryA(), gene$name[gene$symbol==qryA()]))
    message(sprintf("NOTE: qryB = %s \"%s\"", qryB(), ifelse(is.null(qryB()), "(None)", gene$name[gene$symbol==qryB()])))
  })
  
  qryA <- reactive({
    if (input$randGene>randGeneA_previous) {
      randGeneA_previous <<- input$randGene # Must assign to up-scoped variable.
      qryArand <- sample(gene$symbol, 1)
      updateTextInput(session, "qryA", value=qryArand) #Works better than updateSelectizeInput ??
    }
    if (input$qryA=="") { return(NULL) }
    toupper(input$qryA)
  })
  ensgA <- reactive({
    if (is.null(qryA())) { return(NULL) }
    gene$ENSG[gene$symbol==qryA()]
  })

  qryB <- reactive({
    if (input$qryB %in% c("","NONE","none")) { return(NULL) }
    toupper(input$qryB)
  })
  ensgB <- reactive({
    if (is.null(qryB())) { return(NULL) }
    gene$ENSG[gene$symbol==qryB()]
  })
  
  hits <- reactive({
    if (input$mode!="Search" | is.null(qryA())) { return(NULL) }
    ggc_hits <- ggc[ggc$ENSGA==ensgA()|ggc$ENSGB==ensgA(),]
    if (input$metric=="wRho") {
      ggc_hits["Similarity"] <- round(ggc_hits$wRho, digits=3)
    } else if (input$metric=="Ruzicka") {
      ggc_hits["Similarity"] <- round(ggc_hits$Ruzicka, digits=3)
    } else {
      ggc_hits["Similarity"] <- round(ggc_hits$Combo, digits=3)
    }
    ggc_hits["EnsemblID"] <- NA #Populate with non-query gene.
    ggc_hits$EnsemblID[ggc_hits$ENSGA==ensgA()] <- ggc_hits$ENSGB[ggc_hits$ENSGA==ensgA()]
    ggc_hits$EnsemblID[ggc_hits$ENSGB==ensgA()] <- ggc_hits$ENSGA[ggc_hits$ENSGB==ensgA()]
    ggc_hits <- merge(ggc_hits, gene[,c("ENSG","symbol","name")], by.x="EnsemblID", by.y="ENSG", all.x=T, all.y=F)
    ggc_hits <- ggc_hits[,c("EnsemblID","symbol","name","Cluster","Similarity")]
    ggc_hits <- ggc_hits[order(-ggc_hits$Similarity),]
    return(ggc_hits)
  })
  
  hit <- reactive({
    if (is.null(qryA()) | input$mode!="Search") { return(NULL) }
    if (is.null(hits())) { return(NULL) }
    hit_best <- hits()$EnsemblID[1]
    if (!is.na(hit_best)) {
      sim <- hits()$Similarity[1]
      message(sprintf("NOTE: best hit [sim=%.2f]: %s:%s (%s)", sim, hit_best, gene$symbol[gene$ENSG==hit_best], gene$name[gene$ENSG==hit_best]))
    } else {
      message(sprintf("ERROR: search failed."))
    }
    hit_best
  })
  hit_symbol <- reactive({
    if (is.null(hit())) { return("") }
    gene$symbol[gene$ENSG==hit()]
  })
  
  output$result_htm <- reactive({
    if (input$mode == "Compare" & !is.null(qryB())) {
      htm <- sprintf("<B>Results (%s):</B> Gene queryA: %s \"%s\" ; queryB: %s \"%s\"", input$mode, qryA(), gene$name[gene$ENSG==ensgA()], qryB(), gene$name[gene$ENSG==ensgB()])
    } else if (input$mode=="Search" & !is.null(hit())) {
      htm <- sprintf("<B>Results (%s):</B> Gene query: %s \"%s\" ; metric: %s ; profiles found: %d ; top hit: %s \"%s\"", input$mode, qryA(), gene$name[gene$ENSG==ensgA()], input$metric, nrow(hits()), hit_symbol(), gene$name[gene$ENSG==hit()])
      sim <- hits()$Similarity[1]
    } else if (input$mode=="View" & !is.null(qryA())) {
      htm <- sprintf("<B>Results (%s):</B> Gene query: %s \"%s\"", input$mode, qryA(), gene$name[gene$ENSG==ensgA()])
    } else {
      htm <- ""
    }
    return(htm)
  })

  output$log_htm <- reactive({
    htm <- db_htm
    htm
  })

  ### Creates input$datarows_rows_selected
  output$datarows <- renderDataTable({
    if (is.null(hits())) { return(NULL) }
    DT::datatable(data=hits(), rownames=F, selection="multiple", class="cell-border stripe", style="bootstrap",
	options=list(autoWidth=T), colnames=c("Ensembl", "Symbol", "Name", "Group", input$metric)) %>%
        formatRound(digits=3, columns=5:ncol(hits()))
  }, server=T)
  
  
  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    ggc_hits <- hits()
    ggc_hits["Query"] <- qryA()
    ggc_hits <- ggc_hits[,c(6,1,2,3,4,5)] #reorder cols
    names(ggc_hits) <- c("Query", "EnsemblID","GeneSymbol","GeneName","Group", input$metric)
    return(ggc_hits)
  })

  output$hits_file <- downloadHandler(
    filename = function() { "exfiles_hits.csv" },
    content = function(file) {
      if (is.null(hits_export())) { return(NULL) }
      write.csv(hits_export(), file, row.names=F) 
  })

  output$plot <- renderPlotly({
    if (is.null(ensgA())) { return(NULL) }
    
    qryA_profile_f <- log10(as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="F",][1,tissue$name])+1)
    qryA_profile_m <- log10(as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="M",][1,tissue$name])+1)
    qryA_profile <- (qryA_profile_m + qryA_profile_f)/2

    wrhoAfm <- wPearson(qryA_profile_f, qryA_profile_m )
    ruzAfm <- Ruzicka(qryA_profile_f, qryA_profile_m )
    
    if (input$mode=="Compare" & !is.null(qryB())) {
      qryB_profile_f <- log10(as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="F",][1,tissue$name])+1)
      qryB_profile_m <- log10(as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="M",][1,tissue$name])+1)
      qryB_profile <- (qryB_profile_m + qryB_profile_f)/2
      #
      wrho <- wPearson(qryA_profile, qryB_profile)
      wrhoBfm <- wPearson(qryB_profile_f, qryB_profile_m)
      wrhoFab <- wPearson(qryA_profile_f, qryB_profile_f)
      wrhoMab <- wPearson(qryA_profile_m, qryB_profile_m)
      #
      ruz <- Ruzicka(qryA_profile, qryB_profile)
      ruzBfm <- Ruzicka(qryB_profile_f, qryB_profile_m)
      ruzFab <- Ruzicka(qryA_profile_f, qryB_profile_f)
      ruzMab <- Ruzicka(qryA_profile_m, qryB_profile_m)
      #
    } else if (input$mode=="Search" & !is.null(hit())) { ## Hit only if Search
      hit_profile_f <- log10(as.numeric(eps[eps$ENSG==hit() & eps$SEX=="F",][1,tissue$name])+1)
      hit_profile_m <- log10(as.numeric(eps[eps$ENSG==hit() & eps$SEX=="M",][1,tissue$name])+1)
      hit_profile <- (hit_profile_m + hit_profile_f)/2
      #
      # "B" = hit
      wrho <- wPearson(qryA_profile, hit_profile)
      wrhoBfm <- wPearson(hit_profile_f, hit_profile_m)
      wrhoFab <- wPearson(qryA_profile_f, hit_profile_f)
      wrhoMab <- wPearson(qryA_profile_m, hit_profile_m)
      #
      ruz <- Ruzicka(qryA_profile, hit_profile)
      ruzBfm <- Ruzicka(hit_profile_f, hit_profile_m)
      ruzFab <- Ruzicka(qryA_profile_f, hit_profile_f)
      ruzMab <- Ruzicka(qryA_profile_m, hit_profile_m)
    } #Else: View 

    ### PLOT:

    if (input$mode=="Compare" & !is.null(qryB())) {
      titletxt = sprintf("GTEx Gene-Tissue Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), qryB(), gene$name[gene$ENSG==ensgA()], gene$name[gene$ENSG==ensgB()])
    } else if (input$mode=="Search" & !is.null(hit())) {
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), hit_symbol(), gene$name[gene$ENSG==ensgA()], gene$name[gene$ENSG==hit()])
    } else { #View
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s<BR>\"%s\"", qryA(), gene$name[gene$ENSG==ensgA()])
    }
    
    xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10), categoryorder = "array", categoryarray = tissue$name)
    yaxis = list(title="Expression: LOG<SUB>10</SUB>(1+TPM)")
    
    p <- plot_ly() %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
         title = titletxt,
         margin = list(t=100, r=80, b=160, l=60),
         legend = list(x=.9, y=1),
         font = list(family="Arial", size=14)
      )

    if (input$sabv) {
      p <-  add_trace(p, name = paste("(F)", qryA()), x = tissue$name, y = qryA_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$name))
      p <- add_trace(p, name = paste("(M)", qryA()), x = tissue$name, y = qryA_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$name))
      if (input$mode=="Compare" & !is.null(qryB())) {
        annotxt <- paste(sep="<BR>",
            sprintf("Fab: rho = %.2f; ruz = %.2f", wrhoFab, ruzFab),
            sprintf("Mab: rho = %.2f; ruz = %.2f", wrhoMab, ruzMab),
            sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm),
            sprintf("Bfm: rho = %.2f; ruz = %.2f", wrhoBfm, ruzBfm))
        p <- add_trace(p, name = paste("(F)", qryB()), x = tissue$name, y = qryB_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$name))
        #
        p <- add_trace(p, name = paste("(M)", qryB()), x = tissue$name, y = qryB_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$name))
      } else if (input$mode=="Search" & !is.null(hit())) {
        annotxt <- paste(sep="<BR>",
            sprintf("Fab: rho = %.2f; ruz = %.2f", wrhoFab, ruzFab),
            sprintf("Mab: rho = %.2f; ruz = %.2f", wrhoMab, ruzMab),
            sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm),
            sprintf("Bfm: rho = %.2f; ruz = %.2f", wrhoBfm, ruzBfm))
        p <- add_trace(p, name = paste("(F)", hit_symbol()), x = tissue$name, y = hit_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit_symbol(), ": ", tissue$name))
        p <- add_trace(p, name = paste("(M)", hit_symbol()), x = tissue$name, y = hit_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit_symbol(), ": ", tissue$name))
        ##
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          for (i in rows_selected) {
            if (hits()$EnsemblID[i]==hit()) { next; }
            profile_f <- log10(as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="F",][1,tissue$name])+1)
            profile_m <- log10(as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="M",][1,tissue$name])+1)
            p <- add_trace(p, name = paste("(M)", hits()$symbol[i]), x = tissue$name, y = profile_m,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$name))
            p <- add_trace(p, name = paste("(F)", hits()$symbol[i]), x = tissue$name, y = profile_f,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$name))
          }
        }
      } else { #View
        annotxt <- sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm)
      }
    } else { #NOT_SABV
      p <-  add_trace(p, name = qryA(), x = tissue$name, y = qryA_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$name))
      if (input$mode=="Compare" & !is.null(qryB())) {
        annotxt <- sprintf("rho = %.2f; ruz = %.2f", wrho, ruz)
        p <- add_trace(p, name = qryB(), x = tissue$name, y = qryB_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$name))
      } else if (input$mode=="Search" & !is.null(hit())) {
        annotxt <- sprintf("rho = %.2f; ruz = %.2f", wrho, ruz)
        p <- add_trace(p, name = hit_symbol(), x = tissue$name, y = hit_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissue$name))
        ###
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          for (i in rows_selected) {
            if (hits()$EnsemblID[i]==hit()) { next; }
            profile_f <- log10(as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="F",][1,tissue$name])+1)
            profile_m <- log10(as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="M",][1,tissue$name])+1)
            profile <- (profile_f + profile_m)/2
            p <- add_trace(p, name = hits()$symbol[i], x = tissue$name, y = profile,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$name))
          }
        }
      } else { #View
        annotxt <- ""
      }
    }
    #
    if ("Annotate" %in% input$plot_opts) {
      p <- add_annotations(p, text=annotxt, showarrow=F, x=.1, y=1, xref="paper", yref="paper")
    }
    p$elementId <- NULL #Hack to suppress spurious warnings.
    return(p)
  })
}

###
shinyApp(ui, server)
