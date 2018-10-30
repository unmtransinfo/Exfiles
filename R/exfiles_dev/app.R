##########################################################################################
### Modes:
###      VIEW: geneA
### SIMSEARCH: via profile similarity
### TXTSEARCH: via symbols or names
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
library(shiny, quietly=T)
library(shinyBS, quietly=T)
library(DT, quietly=T) #datatable
library(dplyr, quietly=T)
library(plotly, quietly=T)

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
  idg <- read_delim("gtex_gene_idg.tsv", "\t") # IDG gene/protein attributes
  idg <- idg[,c("accession", "idgTDL", "idgFamily")]
  colnames(idg) <- c("uniprot", "idgTDL", "idgDTO")
  ###
  # ENSG, SEX, tissue.1, tissue.2, etc.
  ###
  eps <- read_delim("exfiles_eps.tsv", "\t", col_types=cols(SEX=col_character())) # expression profiles
  ###
  # ENSGA, ENSGB, Group, wRho, Ruzicka
  ###
  ggc <- read_delim("exfiles_ggc.tsv.gz", "\t", col_types="cccdd")
  ggc$Group[ggc$Group=="C"] <- "N"
  ###
  save(tissue, gene, idg, eps, ggc, file="exfiles.Rdata")
}
#
tissue <- tissue[tissue$SMTSD %in% colnames(eps),]
eps <- eps[,c("ENSG","SEX",tissue$SMTSD)]
#
ensgs <- intersect(eps$ENSG, c(ggc$ENSGA,ggc$ENSGB))
ensgs <- intersect(ensgs, gene$ENSG)
#
gene <- gene[gene$ENSG %in% ensgs,]
gene <- gene[!is.na(gene$symbol),]
gene <- gene[!duplicated(gene$ENSG),]
gene <- gene[!duplicated(gene$symbol),]
gene <- merge(gene, idg, by="uniprot", all.x=T, all.y=F)
#
message(sprintf("Tissue count: %d",nrow(tissue)))
message(sprintf("Gene count: %d", nrow(gene)))
message(sprintf("Gene unique ENSG count: %d", length(unique(gene$ENSG))))
message(sprintf("Gene unique SYMB count: %d", length(unique(gene$symbol))))
message(sprintf("Gene unique UniProt count: %d", length(unique(gene$uniprot))))
#
#
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
OrderGeneSymbols <- function(symbols) {
  genes <- data.frame(symbol=symbols, prefix=NA, i=NA)
  are_prefixed <- grepl("^[A-Z]+[0-9]+$", genes$symbol)
  if (sum(!are_prefixed)>0) {
    genes$prefix[!are_prefixed] <- as.character(genes$symbol[!are_prefixed])
  }
  if (sum(are_prefixed)>0) {
   genes$prefix[are_prefixed] <- sub("[0-9]+$", "", genes$symbol[are_prefixed])
   genes$i[are_prefixed] <- sub("^[A-Z]+", "", genes$symbol[are_prefixed])
  }
  genes$i <- as.integer(genes$i)
  return(order(genes$prefix, genes$i, na.last=T))
}
#############################################################################
HelpHtm <- function() {(
"<P><B>Ex-files SABV</B> allows exploration and analysis of co-expression patterns via gene expression profiles,
from <B>GTEx</B> RNA-seq data, with <B>Sex As a Biological Variable (SABV)</B>.
Gene expression profiles are computed as real valued vectors of expression levels 
across the defined tissue types.</P>
<P>
<B>Inputs</B> are query GeneA, and <I>optionally</I>, GeneB.
Ex-files modes of operation:
<UL>
<LI><B>View</B> - view profile for one gene
<LI><B>Compare</B> - compare input genes via profiles
<LI><B>SimSearch</B> - search for genes based on profile similarity
<LI><B>TxtSearch</B> - search for genes based on symbols or names
</UL>
<P>
<B>Score:</B>
<UL>
<LI><B>Ruzicka</B> - Similarity measure, size-normalizing, hence advantageous over Euclidean, RMSD, etc.
<LI><B>wRho</B> - Weighted Pearson correlation coefficient, weighted by average values for each tissue, to mitigate noise.
<LI><B>Combo</B> - Product wRho*Ruzicka, scoring function balancing correlation and similarity.
</UL>
<B>Results:</B>
<UL>
<LI><B>View:</B> plots query GeneA with sex-specific profiles.
<LI><B>Compare:</B> plots query GeneA vs. GeneB with sex-specific profiles.
<LI><B>SimSearch:</B> returns hits based on the GeneA query. The top hit is displayed against GeneA in the plot, and other hits displayed in a table, plotted by manual selection, and downloadable as TSV.
<LI><B>TxtSearch:</B> returns hits based on the query, which may be a substring, regular expression, or list of exact symbols. Hits are displayed in the table, plotted by manual selection, and downloadable as TSV.
<LI><B>Groups</B> denote comparisons among F-only, M-only or N (Non-sexed) representing un-stratified subjects, computed as mean of F and M profiles.
<LI>Expression units: <B>TPM</B> = RNA-seq Transcripts Per Million-kilobase or <B>LOG<SUB>10</SUB>(1+TPM)</B>.
<LI>Expression profiles are computed as medians for GTEx samples, by tissue and sex.
<LI>Option Dissim[ilarity] reverse-sorts search results. To find anti-correlated genes, search by wRho.
<LI><B>Sex-linked genes</B> may be identified via chromosomal location (X*/Y*).
<LI><B>IDG (Illuminating the Druggable Genome)</B> provides target development level (TDL), Drug Target Ontology (DTO) classification,
and links to IDG data portal Pharos.
</UL></P>
<P>
<B>Notes on data preparation:</B> This version is focused on SABV knowledge discovery, thus reproductive and 
breast tissues not included. Also limited to protein-encoding genes unambiguously mapped to HUGO gene symbols. Search results where ruz&lt;0.50 and |wRho|&lt;0.50 are excluded from dataset.
</P>
<B>Reference:</B>\"Ex-files: Sex-Specific Gene Expression Profiles Explorer\", Bocci et al., <i>[submitted]</i>.<BR/>
<B>Authors:</B> Giovanni Bocci, Oleg Ursu, Cristian Bologa, Steve Mathias, Jeremy Yang &amp; Tudor Oprea, <A HREF=\"http://datascience.unm.edu\" target=\"_blank\">Translational Informatics Division, University of New Mexico</A><BR/>
<B>Correspondence</B> from users of this app is welcome, and should be directed to <a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
Data from <A HREF=\"https://www.gtexportal.org/\" TARGET=\"_blank\">GTEx, The Genotype-Tissue Expression Project</A>.<BR/>
Built with R-Shiny &amp; Plotly.<BR/>
This work was supported by the National Institutes of Health grants OT3-OD025464 and U24-CA224370.<BR/>
")}
#
#############################################################################
ui <- fluidPage(
  titlePanel(h2(sprintf("%s, GTEx expression-profile exploration", APPNAME), em("SABV"), 
                span(icon("venus",lib="font-awesome"), icon("mars", lib="font-awesome"), icon("lightbulb", lib="font-awesome"))), 
             windowTitle=paste(APPNAME, "SABV")),
  fluidRow(
    column(4, 
        wellPanel(
    radioButtons("mode", "Mode", choices=c("View", "Compare", "SimSearch", "TxtSearch"), selected="View", inline=T),
	conditionalPanel(condition="input.mode == 'TxtSearch'",
	      textInput("qryTxt", "Query (list or regex)", width='80%'),
	      checkboxInput("iCase","IgnoreCase", value=T, width='20%')),
	conditionalPanel(condition="input.mode != 'TxtSearch'",
          selectizeInput("qryA", label="GeneA", choices = gene_choices, selected=qryArand)),
	conditionalPanel(condition="input.mode == 'Compare'",
          selectizeInput("qryB", label="GeneB", choices = c(list('None'='none'), gene_choices))),
	  conditionalPanel(condition="input.mode == 'SimSearch'",
            radioButtons("score", "Score", choices=c("Ruzicka", "wRho", "Combo"), selected="Combo", inline=T)),
	  conditionalPanel(condition="input.mode == 'TxtSearch'",
            radioButtons("field", "Field", choices=c("Symbol", "Name"), selected="Symbol", inline=T)),
          checkboxGroupInput("groups", "Groups", choices=c("F","M","N"), selected=c("N"), inline=T),
          checkboxGroupInput("opts", "Output", choices=c("IDG", "Dissim","LogY","Annplot"), selected=c("IDG","LogY"), inline=T),
          br(),
          actionButton("randGene", "Demo", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
          actionButton("goRefresh", "Refresh", style='padding:4px; background-color:#DDDDDD; font-weight:bold'),
          actionButton("showhelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
          )),
    column(8, conditionalPanel(condition="true", plotlyOutput("plot", height = "580px")))
  ),
  conditionalPanel(condition="true",
    wellPanel(fluidRow(column(12, DT::dataTableOutput("datarows"))))),
  conditionalPanel(condition="output.datarows_exist", #NOT WORKING. WHY?
    wellPanel(fluidRow(column(12, downloadButton("hits_file", label="Download"))))),
  fluidRow(column(12, wellPanel(
    htmlOutput(outputId = "result_htm", height = "60px")))),
  fluidRow(column(12, wellPanel(
      htmlOutput(outputId = "log_htm", height = "60px")))),
  fluidRow(
    column(12, em(strong(sprintf("%s", APPNAME)), " web app from ", 
	tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
	" and ",
	tags$a(href="https://nihdatacommons.us", target="_blank", span("DCPPC", tags$img(id="dcppc_logo", height="60", valign="bottom", src="dcppc_logo_only.png"))),
	" data from ",
	tags$a(href="https://gtexportal.org", target="_blank", span("GTEx", tags$img(id="gtex_logo", height="50", valign="bottom", src="GTEx_logo_only.png"))),
	" links to ",
	tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png")))
	))),
  bsTooltip("qryA", "Needed for all modes.", "top"),
  bsTooltip("qryB", "Needed for Compare mode only", "top"),
  bsTooltip("qryTxt", "Enter space-separated list of gene symbols, or regular expression.", "top"),
  bsTooltip("mode", "View 1 gene, Compare 2 genes, or Search for similar genes.", "top"),
  bsTooltip("score", "Ruzicka similarity, Pearson weighted correlation, or combination of both.", "top"),
  bsTooltip("groups", "Search F-only, M-only, or Non-sexed comparisons.", "top"),
  bsTooltip("opts", "Output options affecting plot and datatable but not query logic.", "top"),
  bsTooltip("randGene", "Random GeneA query.", "top"),
  bsTooltip("goRefresh", "Refresh plot, output.", "top"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("dcppc_logo", "NIH Data Commons Pilot Phase Consortium", "right"),
  bsTooltip("gtex_logo", "GTEx, Genotype-Tissue Expression project", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right")
)

#############################################################################
server <- function(input, output, session) {
  observeEvent(input$showhelp, {
    showModal(modalDialog(easyClose=T, footer=tagList(modalButton("Dismiss")),
      title=HTML(sprintf("<H2>%s Help</H2>", paste(APPNAME, "SABV"))),
      HTML(HelpHtm())
    ))
  })
  
  Sys.sleep(1)
  randGeneA_previous <- 0 # initialize once per session
  message(sprintf("NOTE: genes: %d ; correlations = %d", nrow(gene), nrow(ggc)))
  observe({
    message(sprintf("NOTE: mode: %s", input$mode))
    if (input$mode=="SimSearch") { message(sprintf("NOTE: score: %s", input$score)) }
    if (!is.null(qryA())) { message(sprintf("NOTE: qryA = %s \"%s\"", qryA(), gene$name[gene$symbol==qryA()])) }
    if (!is.null(qryB())) { message(sprintf("NOTE: qryB = %s \"%s\"", qryB(), gene$name[gene$symbol==qryB()])) }
  })

  qryA <- reactive({
    input$goRefresh # Re-run this and downstream on action button.
    if (input$mode=="TxtSearch") { return(NULL) }
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
    ensg <- gene$ENSG[gene$symbol==qryA()]
    if (length(ensg)>1) {
      message(sprintf("ERROR: multiple ENSGs for qryA=%s: %s", qryA(), paste(collapse=",", ensg)))
      return(ensg[1])
    } else {
      return(ensg)
    }
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
    ensg <- gene$ENSG[gene$symbol==qryB()]
    if (length(ensg)>1) {
      message(sprintf("ERROR: multiple ENSGs for qryB=%s: %s", qryB(), paste(collapse=",", ensg)))
      return(ensg[1])
    } else {
      return(ensg)
    }
  })
  chrB <- reactive({
    if (is.null(ensgB())) { return(NULL) }
    gene$chr[gene$ENSG==ensgB()]
  })

  # If space-separated list, convert input query to regular expression.
  qryRex <- reactive({
    if (is.null(input$qryTxt)) { return(NULL) }
    if (gsub(" ","",input$qryTxt)=="") { return(NULL) }
    qtxt <- sub("^ *(.*[^ ]) *$","\\1",input$qryTxt)
    if (nchar(qtxt)<3) { return(NULL) }
    if (grepl("[, ]", qtxt)) {
      vals <- strsplit(qtxt, "[, ]+")[[1]]
      qrex <- sprintf("(^%s$)", paste0(vals, collapse="$|^"))
    } else {
      qrex <- qtxt
    }
    #message(sprintf("DEBUG: qrex = \"%s\"", qrex))
    return(qrex)
  })
  
  hits <- reactive({
    if (input$mode == "SimSearch" ) {
      if (is.null(qryA())) { return(NULL) }

      if (sum(ggc$ENSGA==ensgA()|ggc$ENSGB==ensgA())==0) { return(NULL) }
      ggc_hits <- ggc[ggc$ENSGA==ensgA()|ggc$ENSGB==ensgA(),]
      if (sum(ggc_hits$Group %in% input$groups)==0) { return(NULL) }
      ggc_hits <- ggc_hits[ggc_hits$Group %in% input$groups,]
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
      gene_cols <- c("ENSG", "uniprot", "symbol","name","chr","idgTDL","idgDTO")
      ggc_hits <- merge(ggc_hits, gene[,gene_cols], by.x="EnsemblID", by.y="ENSG", all.x=T, all.y=F)
      hits_cols <- c("EnsemblID","uniprot","symbol","name","chr","idgTDL","idgDTO","Group","Score") #datatable() ref by # (0+)
      ggc_hits <- ggc_hits[,hits_cols]
      if ("Dissim" %in% input$opts) {
        ggc_hits <- ggc_hits[order(ggc_hits$Score),]
      } else {
        ggc_hits <- ggc_hits[order(-ggc_hits$Score),]
      }
      return(ggc_hits)
    } else if (input$mode == "TxtSearch") {
      if (is.null(qryRex())) { return(NULL) }
      if (input$field == "Symbol") {
        gene_hits <- gene[grepl(qryRex(), gene$symbol, ignore.case=input$iCase),]
      } else if (input$field == "Name") {
        gene_hits <- gene[grepl(qryRex(), gene$name, ignore.case=input$iCase),]
      }
      gene_hits <- rename(gene_hits, EnsemblID = ENSG)
      hits_cols <- c("EnsemblID","uniprot","symbol","name","chr","idgTDL","idgDTO") #datatable() ref by # (0+)
      gene_hits <- gene_hits[,hits_cols]
      gene_hits <- gene_hits[OrderGeneSymbols(gene_hits$symbol),]
      return(gene_hits)
    } else { return(NULL) }
  })
  
  hit <- reactive({
    if (input$mode == "SimSearch" ) {
      if (is.null(qryA())) { return(NULL) }
      if (is.null(hits())) { return(NULL) }
      hit_best <- hits()$EnsemblID[1]
      if (!is.na(hit_best)) {
        sim <- hits()$Score[1]
        message(sprintf("NOTE: best hit [sim=%.3f]: %s:%s (%s)", sim, hit_best, gene$symbol[gene$ENSG==hit_best], gene$name[gene$ENSG==hit_best]))
      } else {
        message(sprintf("ERROR: search failed."))
      }
      return(hit_best)
    } else if (input$mode == "TxtSearch") {
      hit_best <- hits()$EnsemblID[1]
      return(hit_best)
    } else { return(NULL) }
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
      htm <- paste0(htm, sprintf(" GeneA: %s \"%s\"%s", qryA(), gene$name[gene$ENSG==ensgA()],
	ifelse(("IDG" %in% input$opts), sprintf(" <A HREF=\"https://pharos.nih.gov/idg/targets/%s\" target=\"_blank\">%s</A>", gene$uniprot[gene$ENSG==ensgA()], icon("external-link",lib="font-awesome")), "")))
      if (!is.null(chrA()) & grepl("^[XY]", chrA())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrA())) }
    }
    if (!is.null(qryB())) {
      htm <- paste0(htm, sprintf("; geneB: %s \"%s\"%s", qryB(), gene$name[gene$ENSG==ensgB()],
	ifelse(("IDG" %in% input$opts), sprintf(" <A HREF=\"https://pharos.nih.gov/idg/targets/%s\" target=\"_blank\">%s</A>", gene$uniprot[gene$ENSG==ensgB()], icon("external-link",lib="font-awesome")), "")))
      if (!is.null(chrB()) & grepl("^[XY]", chrB())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrB())) }
    }
    if (!is.null(qryRex())) {
      htm <- paste0(htm, sprintf(" QueryText: \"%s\"", input$qryTxt))
    }
    if (!is.null(hits())) {
      htm <- paste0(htm, sprintf("; found: %d", nrow(hits())))
      #if (!is.null(hit_chr()) & grepl("^[XY]", hit_chr())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", hit_chr())) }
    }
    return(htm)
  })

  output$log_htm <- reactive({
    htm <- db_htm
    htm
  })

  ### Assigns input$datarows_rows_selected
  ### Hide ENSG, Uniprot columns.
  ### Note col #s start with 0.
  output$datarows <- renderDataTable({
    if (is.null(hits())) { return(NULL) }
    if ("IDG" %in% input$opts) {
      invis_cols <- c(0,1)
    } else {
      invis_cols <- c(0,1,5,6)
    }
    if (input$mode=="SimSearch") {
      colnames=c("EnsemblID", "uniprot", "Symbol", "Name", "Chr", "idgTDL","idgDTO", "Group", input$score)
      center_cols <- c(2,4,5,6,7,8)
    } else if (input$mode=="TxtSearch") {
      colnames=c("EnsemblID", "uniprot", "Symbol", "Name", "Chr", "idgTDL","idgDTO")
      center_cols <- c(2,4,5,6)
    }
    DT::datatable(data=hits(), rownames=F, 
        selection=list(target="row", mode="multiple", selected=c(1)),
	      class="cell-border stripe", style="bootstrap",
	      options=list(dom='tip', #dom=[lftipr]
		autoWidth=T,
		columnDefs = list(
			list(className='dt-center', targets=center_cols),
			list(visible=F, targets=invis_cols)
			)
		), 
	      colnames=colnames) %>%
        formatRound(digits=2, columns=ncol(hits()))
  }, server=T)
  output$datarows_exist <- reactive({ #For conditionalPanel, to test for R NULL.
    return(as.logical(!is.null(hits())))
  })

  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    hits_out <- hits()
    if (input$mode=="SimSearch") {
      hits_out["Query"] <- qryA()
      hits_out <- hits_out[,c(ncol(hits_out),1:(ncol(hits_out)-1))] #Query col 1st.
      names(hits_out) <- c("Query", "EnsemblID","UniProt","GeneSymbol","GeneName","Chr", "idgTDL","idgDTO", "Group", input$score)
      return(hits_out)
    } else if (input$mode=="TxtSearch") {
      names(hits_out) <- c("EnsemblID","UniProt","GeneSymbol","GeneName","Chr", "idgTDL","idgDTO")
      return(hits_out)
    }
  })

  output$hits_file <- downloadHandler(
    filename = function() { "exfiles_hits.tsv" },
    content = function(file) {
      if (is.null(hits_export())) { return(NULL) }
      write_delim(hits_export(), file, delim="\t") 
  })

  plotTitle <- reactive({
    if (input$mode=="Compare" & !is.null(qryB())) {
      titletxt = sprintf("GTEx Gene-Tissue Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), qryB(), gene$name[gene$ENSG==ensgA()], gene$name[gene$ENSG==ensgB()])
    } else if (input$mode=="SimSearch" & !is.null(hit())) {
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s <BR>\"%s\"", qryA(), gene$name[gene$ENSG==ensgA()])
    } else if (input$mode=="TxtSearch" & !is.null(hits())) {
      titletxt = sprintf("Ex-files, GTEx-Profiles: \"%s\"", input$qryTxt)
    } else { #View
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s<BR>\"%s\"", qryA(), gene$name[gene$ENSG==ensgA()])
    }
    return(titletxt)
  })
    
  output$plot <- renderPlotly({
    if (is.null(ensgA()) & is.null(hits())) { return(NULL) }
    if (length(input$groups)==0) { return(NULL) }

    if (!is.null(ensgA())) {
      qryA_profile_f <- as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="F",][1,tissue$SMTSD])
      qryA_profile_m <- as.numeric(eps[eps$ENSG==ensgA() & eps$SEX=="M",][1,tissue$SMTSD])
      rhoAfm <- wPearson(qryA_profile_f, qryA_profile_m)
      ruzAfm <- Ruzicka(qryA_profile_f, qryA_profile_m)
      qryA_profile_n <- (qryA_profile_f+qryA_profile_m)/2 #Non-sexed (F+M)/2
    }

    if (input$mode=="Compare" & !is.null(qryB())) {
      qryB_profile_f <- as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="F",][1,tissue$SMTSD])
      qryB_profile_m <- as.numeric(eps[eps$ENSG==ensgB() & eps$SEX=="M",][1,tissue$SMTSD])
      qryB_profile_n <- (qryB_profile_f+qryB_profile_m)/2
      #
      rhoNab <- wPearson(qryA_profile_n, qryB_profile_n)
      rhoBfm <- wPearson(qryB_profile_f, qryB_profile_m)
      rhoFab <- wPearson(qryA_profile_f, qryB_profile_f)
      rhoMab <- wPearson(qryA_profile_m, qryB_profile_m)
      #
      ruzNab <- Ruzicka(qryA_profile_n, qryB_profile_n)
      ruzBfm <- Ruzicka(qryB_profile_f, qryB_profile_m)
      ruzFab <- Ruzicka(qryA_profile_f, qryB_profile_f)
      ruzMab <- Ruzicka(qryA_profile_m, qryB_profile_m)
      #
    }

    ### PLOT:

    xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10), categoryorder = "array", categoryarray = tissue$SMTSD)
    yaxis = list(title=ifelse("LogY" %in% input$opts, "Expression: LOG<SUB>10</SUB>(1+TPM)", "Expression: TPM"))
    
    p <- plot_ly() %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
         title = plotTitle(),
         margin = list(t=100, r=80, b=160, l=60),
         legend = list(x=.9, y=1),
         showlegend=T,
         font = list(family="Arial", size=14)
      )

    annos <- c()

    if (!is.null(ensgA())) {
      if ("F" %in% input$groups) {
        p <-  add_trace(p, name = paste("(F)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_f, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryA(), ": ", tissue$SMTSD))
    }
      if ("M" %in% input$groups) {
        p <- add_trace(p, name = paste("(M)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_m, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryA(), ": ", tissue$SMTSD))
    }
      if ("N" %in% input$groups) {
        p <- add_trace(p, name = paste("(N)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_n, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryA(), ": ", tissue$SMTSD))
      }
    }
    if (input$mode=="Compare" & !is.null(qryB())) {
      if ("F" %in% input$groups) {
        p <- add_trace(p, name = paste("(F)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_f, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryB(), ": ", tissue$SMTSD))
        annos <- c(annos, sprintf("Fab: rho = %.2f; ruz = %.2f", rhoFab, ruzFab))
      }
      if ("M" %in% input$groups) {
        p <- add_trace(p, name = paste("(M)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_m, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryB(), ": ", tissue$SMTSD))
        annos <- c(annos, sprintf("Mab: rho = %.2f; ruz = %.2f", rhoMab, ruzMab))
      }
      if ("N" %in% input$groups) {
        p <- add_trace(p, name = paste("(N)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_n, ("LogY" %in% input$opts)),
          type = 'scatter', mode = 'lines+markers',
          marker = list(symbol="circle", size=10),
          text = paste0(qryB(), ": ", tissue$SMTSD))
        annos <- c(annos, sprintf("Nab: rho = %.2f; ruz = %.2f", rhoNab, ruzNab))
      }
    } else if (input$mode %in% c("SimSearch","TxtSearch") & !is.null(hit())) {
      ##
      # Include genes selected via interactive table.
      # Each row has hits()$Group[i] F|M|N so how to handle that?
      ##
      rows_selected <- input$datarows_rows_selected
      if (!is.null(rows_selected)) {
        for (i in rows_selected) {
          hit_profile_f <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="F",][1,tissue$SMTSD])
          hit_profile_m <- as.numeric(eps[eps$ENSG==hits()$EnsemblID[i] & eps$SEX=="M",][1,tissue$SMTSD])
          hit_profile_n <- (hit_profile_f+hit_profile_m)/2
          #
          if ("F" %in% input$groups) {
            p <- add_trace(p, name = paste("(F)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(hit_profile_f, ("LogY" %in% input$opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
            if (!is.null(ensgA())) { annos <- c(annos, sprintf("Fab: rho = %.2f; ruz = %.2f", wPearson(qryA_profile_f, hit_profile_f), Ruzicka(qryA_profile_f, hit_profile_f))) }
          }
          if ("M" %in% input$groups) {
            p <- add_trace(p, name = paste("(M)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(hit_profile_m, ("LogY" %in% input$opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
            if (!is.null(ensgA())) { annos <- c(annos, sprintf("Mab: rho = %.2f; ruz = %.2f", wPearson(qryA_profile_m, hit_profile_m), Ruzicka(qryA_profile_m, hit_profile_m))) }
          }
          if ("N" %in% input$groups) {
            hit_profile_n <- (hit_profile_f + hit_profile_m)/2
            p <- add_trace(p, name = paste("(N)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(hit_profile_n, ("LogY" %in% input$opts)),
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(hits()$symbol[i], ": ", tissue$SMTSD))
            if (!is.null(ensgA())) { annos <- c(annos, sprintf("Nab: rho = %.2f; ruz = %.2f", wPearson(qryA_profile_n, hit_profile_n), Ruzicka(qryA_profile_n, hit_profile_n))) }
          }
        }
      }
    } else { #View
      annos <- c(sprintf("Afm: rho = %.2f; ruz = %.2f", rhoAfm, ruzAfm))
    }
    #
    if ("Annplot" %in% input$opts) {
      p <- add_annotations(p, text=paste0(collapse="<br>", annos), showarrow=F, x=.1, y=1, xref="paper", yref="paper")
    }
    p$elementId <- NULL #Hack to suppress spurious warnings.
    return(p)
  })
}
###
shinyApp(ui, server)
#
