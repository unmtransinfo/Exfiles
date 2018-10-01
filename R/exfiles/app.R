##########################################################################################
### Modes:
###      VIEW: geneA
###    SEARCH: geneA similarity search
###   COMPARE: geneA vs. geneB
### SABV: sex-as-biological-variable analysis
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
if (file.exists("exfiles.Rdata")) {
  message(sprintf("Loading dataset from Rdata..."))
  load("exfiles.Rdata")
} else {
  message(sprintf("Loading dataset from files, writing Rdata..."))
  ###
  # i, SMTS, SMTSD
  ###
  tissue <- read_delim("exfiles_tissue_order.tsv", "\t")
  ###
  # ENSG, NCBI, HGNCID, chr, uniprot, symbol, name
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
#
tissue_missing <- setdiff(tissue$SMTSD, colnames(eps))
if (length(tissue_missing)>0) {
  message(sprintf("NOTE: TISSUE_MISSING: %d. %s\n", 1:length(tissue_missing), tissue_missing))
} else {
  message(sprintf("All tissues found."))
}
tissue <- tissue[tissue$SMTSD %in% colnames(eps),]
eps <- eps[,c("ENSG","SEX",tissue$SMTSD)]
#
message(sprintf("Tissue count: %d",nrow(tissue)))
message(sprintf("%d. %s : %s\n", tissue$i, tissue$SMTS, tissue$SMTSD))
#
ensgs <- intersect(eps$ENSG, c(ggc$ENSGA,ggc$ENSGB))
ensgs <- intersect(ensgs, gene$ENSG)
message(sprintf("Gene count: %d",length(ensgs)))
#
ensg_dups <- gene$ENSG[duplicated(gene$ENSG)]
message(sprintf("Duplicated/ambiguous gene IDs: %d", length(ensg_dups)))
#message(sprintf("Duplicated/ambiguous gene ID: %s\n", ensg_dups))
#
gene <- gene[gene$ENSG %in% ensgs,]
#
message(sprintf("Unknown/unmapped gene SYMBs: %d", sum(is.na(gene$symbol))))
#
gene <- gene[!is.na(gene$symbol),]
#
symb_dups <- gene$symbol[duplicated(gene$symbol)]
message(sprintf("Duplicated/ambiguous gene SYMBs: %d", length(symb_dups)))
#message(sprintf("Duplicated/ambiguous gene SYMBs: %s\n", symb_dups))
#
gene <- gene[!duplicated(gene$ENSG),]
gene <- gene[!duplicated(gene$symbol),]
#
message(sprintf("Gene unique ENSG count: %d", length(unique(gene$ENSG))))
message(sprintf("Gene unique SYMB count: %d", length(unique(gene$symbol))))
###
ggc$Combo <- round(ggc$wRho*ggc$Ruzicka, digits=2)
#
t_elapsed <- (proc.time()-t0)[3]
#
db_htm <- sprintf("<B>Dataset:</B> Genes: %d ; tissues: %d ; comparisons: %d (t_load: %.1fs)", nrow(gene), nrow(tissue), nrow(ggc), t_elapsed)
###
#
qryArand <- sample(gene$symbol, 1) # initial random query
#
#############################################################################
#
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  ok <- !is.na(A) & !is.na(B)
  A <- A[ok]
  B <- B[ok]
  wCorr::weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
}
###
Ruzicka <- function(A,B) {
  sum(pmin(A,B), rm.na=T)/sum(pmax(A,B), rm.na=T)
}
###
EpLogIf <- function(ep, condition) {
  if (condition) { return(log10(ep+1)) }
  else { return(ep) }
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
<B>Score:</B>
<UL>
<LI><B>Ruzicka</B> - Similarity measure, size-normalizing, hence advantageous over Euclidean.
<LI><B>wRho</B> - Weighted Pearson correlation coefficient, weighted by average values for each tissue, to mitigate noise.
<LI><B>Combo</B> - Product wRho*Ruzicka, scoring function balancing correlation and similarity.
</UL>
<B>Results:</B>
<UL>
<LI>Search results include the top hit, displayed against the query in the
plot, and other hits displayed in a table and downloadable as CSV.  <B>Note: current dataset includes only highly
correlated or anti-correlated profiles.</B>
<LI>Compare results consist of the plot of GeneA vs. GeneB.
<LI>SABV results include both sexes for both genes, thus four profiles.
<LI>Expression units: <B>TPM</B> = RNA-seq median Transcripts Per Million-kilobase or <B>LOG<SUB>10</SUB>(1+TPM)</B>.
<LI>Groups denote comparisons among F-only, M-only or F vs M.
<LI>To find anti-correlated genes, search by wRho and reverse-sort results table.
<LI><B>Sex-linked genes</B> may be identified via chromosomal location (X*/Y*).
</UL></P>
<P>
Notes on data preparation: This version is focused on SABV knowledge discovery, thus reproductive and 
breast tissues not considered. Also we restrict to protein-encoding genes mapped to HUGO gene symbols,
for scientific comprehensibility, so only data associated with these genes is retained.
Currently also genes are ignored with mapping ambiguity between Ensembl ENSG to HUGO symbols. 
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
          radioButtons("score", "Score", choices=c("Ruzicka", "wRho", "Combo"), selected="Combo", inline=T),
          checkboxGroupInput("plot_opts", "Plot", choices=c("LogY", "Annotate"), selected=c("LogY"), inline=T),
          br(),
          actionButton("goRefresh", "Refresh", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
          actionButton("showhelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
          )),
    column(8, conditionalPanel(condition="true", plotlyOutput("plot", height = "580px")))
  ),
  # How to test in JS for R NULL?
  conditionalPanel(condition="(input.mode=='Search' && typeof output.datarows !== 'undefined')",
    wellPanel(
      fluidRow(column(12, DT::dataTableOutput("datarows"))),
      fluidRow(column(12, downloadButton("hits_file", label="Download")))),
    width=12),
  fluidRow(column(12, wellPanel(
    htmlOutput(outputId = "result_htm", height = "60px")))),
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
    message(sprintf("NOTE: mode: %s ; sabv: %s ; score: %s", input$mode, input$sabv, input$score))
    message(sprintf("NOTE: qryA = %s \"%s\"", qryA(), gene$name[gene$symbol==qryA()]))
    message(sprintf("NOTE: qryB = %s \"%s\"", qryB(), ifelse(is.null(qryB()), "(None)", gene$name[gene$symbol==qryB()])))
  })
  
  qryA <- reactive({
    input$goRefresh # Re-run this and downstream on action button.
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
  chrA <- reactive({
    if (is.null(ensgA())) { return(NULL) }
    gene$chr[gene$ENSG==ensgA()]
  })
  
  qryB <- reactive({
    if (input$qryB %in% c("","NONE","none")) { return(NULL) }
    toupper(input$qryB)
  })
  ensgB <- reactive({
    if (is.null(qryB())) { return(NULL) }
    gene$ENSG[gene$symbol==qryB()]
  })
  chrB <- reactive({
    if (is.null(ensgB())) { return(NULL) }
    gene$chr[gene$ENSG==ensgB()]
  })
  
  hits <- reactive({
    if (input$mode!="Search" | is.null(qryA())) { return(NULL) }
    ggc_hits <- ggc[ggc$ENSGA==ensgA()|ggc$ENSGB==ensgA(),]
    if (input$score=="wRho") {
      ggc_hits["Score"] <- round(ggc_hits$wRho, digits=2)
    } else if (input$score=="Ruzicka") {
      ggc_hits["Score"] <- round(ggc_hits$Ruzicka, digits=2)
    } else {
      ggc_hits["Score"] <- round(ggc_hits$Combo, digits=2)
    }
    ggc_hits["EnsemblID"] <- NA #Populate with non-query gene.
    ggc_hits$EnsemblID[ggc_hits$ENSGA==ensgA()] <- ggc_hits$ENSGB[ggc_hits$ENSGA==ensgA()]
    ggc_hits$EnsemblID[ggc_hits$ENSGB==ensgA()] <- ggc_hits$ENSGA[ggc_hits$ENSGB==ensgA()]
    ggc_hits <- merge(ggc_hits, gene[,c("ENSG","symbol","name","chr")], by.x="EnsemblID", by.y="ENSG", all.x=T, all.y=F)
    ggc_hits <- ggc_hits[,c("EnsemblID","symbol","name","chr","Cluster","Score")] #Must match datatable() and hits_export.
    ggc_hits <- ggc_hits[order(-ggc_hits$Score),]
    return(ggc_hits)
  })
  
  hit <- reactive({
    if (is.null(qryA()) | input$mode!="Search") { return(NULL) }
    if (is.null(hits())) { return(NULL) }
    hit_best <- hits()$EnsemblID[1]
    if (!is.na(hit_best)) {
      sim <- hits()$Score[1]
      message(sprintf("NOTE: best hit [sim=%.3f]: %s:%s (%s)", sim, hit_best, gene$symbol[gene$ENSG==hit_best], gene$name[gene$ENSG==hit_best]))
    } else {
      message(sprintf("ERROR: search failed."))
    }
    hit_best
  })
  hit_symbol <- reactive({
    if (is.null(hit())) { return("") }
    gene$symbol[gene$ENSG==hit()]
  })
  hit_chr <- reactive({
    if (is.null(hit())) { return(NULL) }
    gene$chr[gene$ENSG==hit()]
  })
  
  output$result_htm <- reactive({
    htm <- sprintf("<B>Results (%s):</B>", input$mode)
    if (!is.null(qryA())) {
      htm <- paste0(htm, sprintf(" GeneA: %s \"%s\"", qryA(), gene$name[gene$ENSG==ensgA()]))
      if (!is.null(chrA()) & grepl("^[XY]", chrA())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrA())) }
    }
    if (!is.null(qryB())) {
      htm <- paste0(htm, sprintf("; geneB: %s \"%s\"", qryB(), gene$name[gene$ENSG==ensgB()]))
      if (!is.null(chrB()) & grepl("^[XY]", chrB())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrB())) }
    }
    if (!is.null(hit())) {
      htm <- paste0(htm, sprintf("; found: %d ; top hit: %s \"%s\"; score: %s=%.2f", nrow(hits()), hit_symbol(), gene$name[gene$ENSG==hit()], input$score, hits()$Score[1]))
      if (!is.null(hit_chr()) & grepl("^[XY]", hit_chr())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", hit_chr())) }
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
    DT::datatable(data=hits(), rownames=F, selection=list(target="row", mode="multiple", selected=c(1)),
	class="cell-border stripe", style="bootstrap",
	options=list(autoWidth=T), colnames=c("Ensembl", "Symbol", "Name", "Chr", "Group", input$score)) %>%
        formatRound(digits=2, columns=6:ncol(hits()))
  }, server=T)
  
  
  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    ggc_hits <- hits()
    ggc_hits["Query"] <- qryA()
    ggc_hits <- ggc_hits[,c(ncol(ggc_hits),1:(ncol(ggc_hits)-1))] #Query col 1st.
    names(ggc_hits) <- c("Query", "EnsemblID","GeneSymbol","GeneName","Chr", "Group", input$score)
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
    
    qryA_profile_f <- as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="F",][1,tissue$SMTSD])
    qryA_profile_m <- as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="M",][1,tissue$SMTSD])
    wrhoAfm <- wPearson(qryA_profile_f, qryA_profile_m)
    ruzAfm <- Ruzicka(qryA_profile_f, qryA_profile_m)
    
    qryA_profile <- (qryA_profile_m + qryA_profile_f)/2

    if (input$mode=="Compare" & !is.null(qryB())) {
      qryB_profile_f <- as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="F",][1,tissue$SMTSD])
      qryB_profile_m <- as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="M",][1,tissue$SMTSD])
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
      hit_profile_f <- as.numeric(eps[eps$ENSG==hit() & eps$SEX=="F",][1,tissue$SMTSD])
      hit_profile_m <- as.numeric(eps[eps$ENSG==hit() & eps$SEX=="M",][1,tissue$SMTSD])
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
    
    xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10), categoryorder = "array", categoryarray = tissue$SMTSD)
    yaxis = list(title=ifelse("LogY" %in% input$plot_opts, "Expression: LOG<SUB>10</SUB>(1+TPM)", "Expression: TPM"))
    
    p <- plot_ly() %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
         title = titletxt,
         margin = list(t=100, r=80, b=160, l=60),
         legend = list(x=.9, y=1),
         font = list(family="Arial", size=14)
      )

    if (input$sabv) {
      p <-  add_trace(p, name = paste("(F)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_f, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$SMTSD))
      p <- add_trace(p, name = paste("(M)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_m, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$SMTSD))
      if (input$mode=="Compare" & !is.null(qryB())) {
        annotxt <- paste(sep="<BR>",
            sprintf("Fab: rho = %.2f; ruz = %.2f", wrhoFab, ruzFab),
            sprintf("Mab: rho = %.2f; ruz = %.2f", wrhoMab, ruzMab),
            sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm),
            sprintf("Bfm: rho = %.2f; ruz = %.2f", wrhoBfm, ruzBfm))
        p <- add_trace(p, name = paste("(F)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_f, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$SMTSD))
        #
        p <- add_trace(p, name = paste("(M)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_m, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$SMTSD))
      } else if (input$mode=="Search" & !is.null(hit())) {
        annotxt <- paste(sep="<BR>",
            sprintf("Fab: rho = %.2f; ruz = %.2f", wrhoFab, ruzFab),
            sprintf("Mab: rho = %.2f; ruz = %.2f", wrhoMab, ruzMab),
            sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm),
            sprintf("Bfm: rho = %.2f; ruz = %.2f", wrhoBfm, ruzBfm))
        ##
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          for (i in rows_selected) {
            profile_f <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="F",][1,tissue$SMTSD])
            profile_m <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="M",][1,tissue$SMTSD])
            p <- add_trace(p, name = paste("(F)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(profile_f, ("LogY" %in% input$plot_opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
            p <- add_trace(p, name = paste("(M)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(profile_m, ("LogY" %in% input$plot_opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
          }
        }
      } else { #View
        annotxt <- sprintf("Afm: rho = %.2f; ruz = %.2f", wrhoAfm, ruzAfm)
      }
    } else { #NOT_SABV
      p <-  add_trace(p, name = qryA(), x = tissue$SMTSD, y = EpLogIf(qryA_profile, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissue$SMTSD))
      if (input$mode=="Compare" & !is.null(qryB())) {
        annotxt <- sprintf("rho = %.2f; ruz = %.2f", wrho, ruz)
        p <- add_trace(p, name = qryB(), x = tissue$SMTSD, y = EpLogIf(qryB_profile, ("LogY" %in% input$plot_opts)),
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissue$SMTSD))
      } else if (input$mode=="Search" & !is.null(hit())) {
        annotxt <- sprintf("rho = %.2f; ruz = %.2f", wrho, ruz)
        ###
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          for (i in rows_selected) {
            profile_f <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="F",][1,tissue$SMTSD])
            profile_m <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="M",][1,tissue$SMTSD])
            profile <- (profile_f + profile_m)/2
            p <- add_trace(p, name = hits()$symbol[i], x = tissue$SMTSD, y = EpLogIf(profile, ("LogY" %in% input$plot_opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
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
