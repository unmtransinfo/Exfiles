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
### exfiles_tissue_order.tsv - manually curated, from GDoc
### gtex_gene_xref.tsv - from gtex_gene_xref.R
### tcrd_targets.tsv - BioClients.idg.Client 
### exfiles_eps.tsv - expression profiles; gtex_rnaseq_prep_app.py
### exfiles_ggc.tsv - gene-gene comparisons; exfiles_similarity_post.py
##########################################################################################
library(readr)
library(wCorr)
library(shiny, quietly=T)
library(shinyBS, quietly=T)
library(shinysky) #textInput.typeahead
library(DT, quietly=T)
library(data.table, quietly=T)
library(plotly, quietly=T)
#
pkgs <- names(sessionInfo()$otherPkgs)
pkgVerTxt <- paste(sprintf("%s %s", pkgs, sapply(pkgs, function(p){paste(packageVersion(p), collapse=".")})), collapse="; ")
message(pkgVerTxt)
###
# This code runs once for all sessions.
###
APPNAME <- "Exfiles"
APPNAME_FULL <- "Exfiles: Expression profile analytics SABV"
GTEX_RELEASE <- "v8 (2017)"
#
IGNORE_SEX_SPECIFIC_TISSUES <- T
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
  tissue <- read_delim("data/exfiles_tissue_order.tsv", "\t", col_types=cols(.default=col_character(), SEX_SPECIFIC=col_logical()))
  setDT(tissue)
  ###
  # ENSG, NCBI, HGNCID, chr, uniprot, symbol, name
  ###
  gene <- read_delim("data/gtex_gene_xref.tsv", "\t") # gene attributes
  setDT(gene)
  ###
  idg <- read_delim("data/tcrd_targets.tsv", "\t") # IDG gene attributes
  setDT(idg)
  idg <- idg[, .(accession, idgTDL, idgFamily)]
  idg <- idg[!duplicated(accession)]
  setnames(idg, c("uniprot", "idgTDL", "idgDTO"))
  ###
  # ENSG, SEX, tissue.1, tissue.2, etc.
  ###
  eps <- read_delim("data/exfiles_eps.tsv", "\t", col_types=cols(.default=col_double(), ENSG=col_character(), SEX=col_character()))
  setDT(eps)
  ###
  # ENSGA, ENSGB, Group, wRho, Ruzicka
  ###
  ggc <- read_delim("data/exfiles_ggc.tsv", "\t", col_types="cccdd")
  setDT(ggc)
  #
  ensgs <- intersect(eps$ENSG, c(ggc$ENSGA, ggc$ENSGB))
  ensgs <- intersect(ensgs, gene$ENSG)
  eps <- eps[ENSG %in% ensgs]
  ggc <- ggc[(ENSGA %in% ensgs) & (ENSGB %in% ensgs)]
  #
  gene <- gene[ENSG %in% ensgs]
  gene <- merge(gene, idg, by="uniprot", all.x=F, all.y=F)
  gene <- gene[!is.na(symbol)]
  gene <- gene[!duplicated(ENSG)]
  gene <- gene[!duplicated(symbol)]
  ###
  save(tissue, gene, idg, eps, ggc, file="exfiles.Rdata")
}
#
message(sprintf("Gene count (ENSG): %d", uniqueN(gene$ENSG)))
message(sprintf("Gene count (SYMB): %d", uniqueN(gene$symbol)))
message(sprintf("Gene count (UniProt): %d", uniqueN(gene$uniprot)))
message(sprintf("Tissue count (profiles): %d", ncol(eps)-2))
message(sprintf("Tissue count (shown): %d", uniqueN(tissue$SMTSD)))
#
hiddenTissues <- setdiff(names(eps)[3:ncol(eps)], tissue$SMTSD)
for (tis in hiddenTissues) {
  message(sprintf("Tissue hidden (SMTSD): %s", tis))
}
if (IGNORE_SEX_SPECIFIC_TISSUES) {
  tissue <- tissue[(!SEX_SPECIFIC)]
  tissue <- tissue[SMTSD != "Kidney - Medulla"] #No female TPMs in GTEx. Why?
}
tissue <- tissue[SMTSD %in% colnames(eps)]
TAGS_THIS <- c("ENSG", "SEX", tissue$SMTSD)
eps <- eps[, ..TAGS_THIS]
#
###
ggc[, Combo := round(wRho*Ruzicka, digits=2)]
#
db_htm <- sprintf("<B>Dataset:</B> Genes: %d ; tissues: %d ; comparisons: %d", nrow(gene), nrow(tissue), nrow(ggc))
message(sprintf("t_load: %.1fs", (proc.time()-t0)[3]))
message(sprintf("IGNORE_SEX_SPECIFIC_TISSUES: %s", IGNORE_SEX_SPECIFIC_TISSUES))
###
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
  genes <- data.table(symbol=symbols, prefix=NA, i=NA)
  genes[, are_prefixed := grepl("^[A-Z]+[0-9]+$", symbol)]
  genes[(!are_prefixed), prefix := as.character(symbol)]
  genes[(are_prefixed), prefix := sub("[0-9]+$", "", symbol)]
  genes[(are_prefixed), i := as.integer(sub("^[A-Z]+", "", symbol))]
  return(order(genes$prefix, genes$i, na.last=T))
}
#############################################################################
HelpHtm <- function() {
  htm <- sprintf(
"<P><B>%s</B> allows exploration and analysis of co-expression patterns via gene expression profiles,
from <B>GTEx</B> RNA-seq data, with <B>Sex As a Biological Variable (SABV)</B>.
Gene expression profiles are computed as real valued vectors of expression levels 
across the defined tissue types.
<B>Expression data source:</B> <a href=\"https://gtexportal.org\">GTEx</a> [%s] RNA-Seq download file, 
plus samples and subjects metadata/annotations.
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
<LI><B>TxtSearch:</B> returns hits based on the query, which may be a substring, regular expression, or space-separated list of exact symbols. Hits are displayed in the table, plotted by manual selection, and downloadable as TSV.
<LI><B>Groups</B> denote comparisons among F-only, M-only or C (Combined-sex) representing un-stratified subjects, computed as mean of F and M profiles.
<LI>Expression units: <B>TPM</B> = RNA-seq Transcripts Per Million-kilobase or <B>LOG<SUB>10</SUB>(1+TPM)</B>.
<LI>Expression profiles are computed as medians for GTEx samples, by tissue and sex.
<LI>Option Dissim[ilarity] reverse-sorts search results. To find anti-correlated genes, search by wRho.
<LI><B>Sex-linked genes</B> may be identified via chromosomal location (X*/Y*).
<LI><B>IDG (Illuminating the Druggable Genome)</B> provides target development level (TDL), Drug Target Ontology (DTO) classification,
and links to IDG data portal Pharos.
</UL></P>
<P>
<B>Notes on data preparation:</B> This version is focused on SABV knowledge discovery, thus reproductive and 
breast tissues not included in profile comparison computations. Also limited to protein-encoding genes unambiguously 
mapped to HUGO gene symbols. Search results where ruz&lt;0.50 and |wRho|&lt;0.50 are excluded from dataset.
Note that tissue \"Kidney - Medulla\" is not biologically but artifactually sex-specific since in GTEx there are no
female TPM data.
</P>
<B>Reference:</B>\"Exfiles: Sex-Specific Gene Expression Profiles Analytics\", Bocci et al., <i>[in preparation]</i>.<BR/>
<B>Authors:</B> Giovanni Bocci, Oleg Ursu, Cristian Bologa, Steve Mathias, Jeremy Yang &amp; Tudor Oprea, <A HREF=\"http://datascience.unm.edu\" target=\"_blank\">Translational Informatics Division, University of New Mexico</A><BR/>
<B>Correspondence</B> from users of this app is welcome, and should be directed to <a href=\"mailto:jjyang_REPLACE_WITH_ATSIGN_salud.unm.edu\">Jeremy Yang</a>.<br/>
Data from <A HREF=\"https://www.gtexportal.org/\" TARGET=\"_blank\">GTEx, The Genotype-Tissue Expression Project</A>.<BR/>
This work was supported by the National Institutes of Health grants OT3-OD025464 and U24-CA224370.<BR/>
", APPNAME_FULL, GTEX_RELEASE)
  htm <- paste(htm, sprintf("<hr>\nPowered by: <tt>%s; %s</tt>", R.version.string, pkgVerTxt), sep="\n")
  return(htm)
  }
#
#############################################################################
ui <- fluidPage(
  titlePanel(
    tags$table(width="100%", tags$tr(tags$td(
    h2(sprintf("%s", APPNAME_FULL), span(icon("venus", lib="font-awesome"), icon("mars", lib="font-awesome")))),
    tags$td(align="right",
      actionButton("goRefresh", "Refresh", style="padding:4px; background-color:#DDDDDD;font-weight:bold"),
      actionButton("showhelp", "Help", style="padding:4px; background-color:#DDDDDD;font-weight:bold")))),
    windowTitle=APPNAME_FULL),
  fluidRow(
    column(5, 
        wellPanel(
    radioButtons("mode", "Mode", choices=c("View", "Compare", "SimSearch", "TxtSearch"), selected="View", inline=T),
	conditionalPanel(condition="input.mode == 'TxtSearch'",
	      textInput("qryTxt", "Query (list or regex)", width='80%'),
	      checkboxInput("iCase", "IgnoreCase", value=T, width='20%')),
	conditionalPanel(condition="input.mode != 'TxtSearch'",
        shinysky::textInput.typeahead(
		id="qryA"
		,placeholder="GeneA..."
		,local=gene
		,valueKey="symbol"
		,tokens=gene$symbol
		,template=HTML("<p class='repo-name'>{{symbol}}</p> <p class='repo-description'>{{name}}</p>")
		,limit=10)
        ),
	conditionalPanel(condition="input.mode == 'Compare'",
          #selectizeInput("qryB", label="GeneB", choices = c(list('None'='none'), gene_choices))
        shinysky::textInput.typeahead(
		id="qryB"
		,placeholder="GeneB..."
		,local=gene
		,valueKey="symbol"
		,tokens=gene$symbol
		,template=HTML("<p class='repo-name'>{{symbol}}</p> <p class='repo-description'>{{name}}</p>")
		,limit=10)
	),
	  conditionalPanel(condition="input.mode == 'SimSearch'",
            radioButtons("score", "Score", choices=c("Ruzicka", "wRho", "Combo"), selected="Combo", inline=T)),
	  conditionalPanel(condition="input.mode == 'TxtSearch'",
            radioButtons("field", "Field", choices=c("Symbol", "Name"), selected="Symbol", inline=T)),
          checkboxGroupInput("groups", "Groups", choiceValues=c("F", "M", "C"), choiceNames=c("F", "M", "Combined"), selected=c("F", "M"), inline=T),
          checkboxGroupInput("opts", "Output", choices=c("IDG", "Dissim", "LogY", "Annplot"), selected=c("IDG", "LogY", "Annplot"), inline=T)),
	wellPanel(htmlOutput(outputId = "log_htm", height = "120px")),
        wellPanel(htmlOutput(outputId = "result_htm", height = "120px"))
	),
    column(7, conditionalPanel(condition="true", plotlyOutput("plot", height = "580px")))
  ),
  conditionalPanel(condition="output.hits_exist=='TRUE'",
    wellPanel(fluidRow(column(12, DT::dataTableOutput("datarows"))),
              fluidRow(column(12, downloadButton("hits_file", label="Download"))))),
  fluidRow(
    column(12, em(strong(sprintf("%s", APPNAME)), " web app from ", 
	tags$a(href="http://datascience.unm.edu", target="_blank", span("UNM", tags$img(id="unm_logo", height="60", valign="bottom", src="unm_new.png"))),
	" and ",
	tags$a(href="https://druggablegenome.net", target="_blank", span("IDG", tags$img(id="idg_logo", height="60", valign="bottom", src="IDG_logo_only.png"))),
	" built from ",
	tags$a(href="https://gtexportal.org", target="_blank", span("GTEx", tags$img(id="gtex_logo", height="50", valign="bottom", src="GTEx_logo_only.png")))
	))),
  bsTooltip("qryA", "Needed for all modes.", "top"),
  bsTooltip("qryB", "Needed for Compare mode only", "top"),
  bsTooltip("qryTxt", "Enter space-separated list of gene symbols, or regular expression.", "top"),
  bsTooltip("mode", "View 1 gene, Compare 2 genes, or Search for similar genes.", "top"),
  bsTooltip("score", "Ruzicka similarity, Pearson weighted correlation, or combination of both.", "top"),
  bsTooltip("groups", "Search F-only, M-only, or Non-sexed comparisons.", "top"),
  bsTooltip("opts", "Output options affecting plot and datatable but not query logic.", "top"),
  bsTooltip("goRefresh", "Refresh plot, output.", "top"),
  bsTooltip("unm_logo", "UNM Translational Informatics Division", "right"),
  bsTooltip("idg_logo", "IDG, Illuminating the Druggable Genome project", "right"),
  bsTooltip("gtex_logo", "GTEx, Genotype-Tissue Expression project", "right"),
  conditionalPanel(condition="true", span(style="color:white", textOutput("hits_exist"))) #fails_if_hidden
)

#############################################################################
server <- function(input, output, session) {
  
  observeEvent(input$showhelp, {
    showModal(modalDialog(easyClose=T, footer=tagList(modalButton("Dismiss")),
      title=HTML(sprintf("<H2>%s Help</H2>", APPNAME)),
      HTML(HelpHtm())
    ))
  })
  
  Sys.sleep(1)
  #outputOptions(output, "hits_exist", suspendWhenHidden=F) #Error: "hits_exist is not in list of output objects"
  
  message(sprintf("NOTE: genes: %d ; correlations = %d", nrow(gene), nrow(ggc)))
  observe({
    message(sprintf("NOTE: mode: %s", input$mode))
    if (input$mode=="SimSearch") { message(sprintf("NOTE: score: %s", input$score)) }
    if (!is.null(qryA())) { message(sprintf("NOTE: qryA = %s \"%s\"", qryA(), gene[symbol==qryA(), name])) }
    if (!is.null(qryB())) { message(sprintf("NOTE: qryB = %s \"%s\"", qryB(), gene[symbol==qryB(), name])) }
  })

  qryA <- reactive({
    input$goRefresh # Re-run this and downstream on action button.
    if (input$mode=="TxtSearch") { return(NULL) }
    if (input$qryA=="") { return(NULL) }
    toupper(input$qryA)
  })
  ensgA <- reactive({
    if (is.null(qryA())) { return(NULL) }
    ensg <- gene[symbol==qryA(), ENSG]
    if (length(ensg)>1) {
      message(sprintf("ERROR: multiple ENSGs for qryA=%s: %s", qryA(), paste(collapse=",", ensg)))
      return(ensg[1])
    } else {
      message(sprintf("DEBUG: qryA=%s; ensgA=%s", qryA(), paste(collapse=",", ensg)))
      return(ensg)
    }
  })
  chrA <- reactive({
    if (is.null(ensgA())) { return(NULL) }
    return(gene[ENSG==ensgA(), chr])
  })
  
  qryB <- reactive({
    if (input$qryB %in% c("", "NONE", "none")) { return(NULL) }
    toupper(input$qryB)
  })
  ensgB <- reactive({
    if (is.null(qryB())) { return(NULL) }
    ensg <- gene[symbol==qryB(), ENSG]
    if (length(ensg)>1) {
      message(sprintf("ERROR: multiple ENSGs for qryB=%s: %s", qryB(), paste(collapse=",", ensg)))
      return(ensg[1])
    } else {
      return(ensg)
    }
  })
  chrB <- reactive({
    if (is.null(ensgB())) { return(NULL) }
    return(gene[ENSG==ensgB(), chr])
  })

  # If space-separated list, convert input query to regular expression.
  qryRex <- reactive({
    if (is.null(input$qryTxt)) { return(NULL) }
    if (gsub(" ", "", input$qryTxt)=="") { return(NULL) }
    qtxt <- sub("^ *(.*[^ ]) *$", "\\1", input$qryTxt)
    if (nchar(qtxt)<3) { return(NULL) }
    if (grepl(" ", qtxt)) {
      vals <- strsplit(qtxt, " +")[[1]]
      qrex <- sprintf("(^%s$)", paste0(vals, collapse="$|^"))
    } else {
      qrex <- qtxt
    }
    return(qrex)
  })
  
  hits <- reactive({
    if (input$mode == "SimSearch" ) {
      if (is.null(qryA())) { return(NULL) }
      if (sum(ggc$ENSGA==ensgA() | ggc$ENSGB==ensgA())==0) { return(NULL) }
      ggc_hits <- ggc[ENSGA==ensgA() | ENSGB==ensgA()]
      if (sum(ggc_hits$Group %in% input$groups)==0) { return(NULL) }
      ggc_hits <- ggc_hits[Group %in% input$groups]
      if (input$score=="wRho") {
        ggc_hits[, Score := round(wRho, digits=2)]
      } else if (input$score=="Ruzicka") {
        ggc_hits[, Score := round(Ruzicka, digits=2)]
      } else {
        ggc_hits[, Score := round(Combo, digits=2)]
      }
      ggc_hits[, EnsemblID := as.character(NA)] #Populate with non-query gene.
      ggc_hits[ENSGA==ensgA(), EnsemblID := ENSGB]
      ggc_hits[ENSGB==ensgA(), EnsemblID := ENSGA]
      gene_cols <- c("ENSG", "uniprot", "symbol", "name", "chr", "idgTDL", "idgDTO")
      ggc_hits <- merge(ggc_hits, gene[, ..gene_cols], by.x="EnsemblID", by.y="ENSG", all.x=T, all.y=F)
      hits_cols <- c("EnsemblID", "uniprot", "symbol", "name", "chr", "idgTDL", "idgDTO", "Group", "Score") #datatable() ref by # (0+)
      ggc_hits <- ggc_hits[, ..hits_cols]
      if ("Dissim" %in% input$opts) {
        setorder(ggc_hits, Score)
      } else {
        setorder(ggc_hits, -Score)
      }
      return(ggc_hits)
    } else if (input$mode == "TxtSearch") {
      message(sprintf("DEBUG: is.null(qryRex()): %s", as.character(is.null(qryRex()))))
      if (is.null(qryRex())) { return(NULL) }
      if (input$field == "Symbol") {
        gene_hits <- gene[grepl(qryRex(), symbol, ignore.case=input$iCase)]
      } else if (input$field == "Name") {
        gene_hits <- gene[grepl(qryRex(), name, ignore.case=input$iCase)]
      }
      message(sprintf("DEBUG: nrow(gene_hits)=%d", nrow(gene_hits)))
      if (nrow(gene_hits)==0) { return(NULL) }
      setnames(gene_hits, old="ENSG", new="EnsemblID")
      hits_cols <- c("EnsemblID", "uniprot", "symbol", "name", "chr", "idgTDL", "idgDTO") #datatable() ref by # (0+)
      gene_hits <- gene_hits[, ..hits_cols]
      gene_hits <- gene_hits[OrderGeneSymbols(symbol)]
      return(gene_hits)
    } else { return(NULL) }
  })
  
  output$hits_exist <- renderText({
    message(sprintf("DEBUG: hits_exist: %s", as.character(!is.null(hits()))))
    if (!grepl("Search", input$mode)) { return("FALSE") }
    return(as.character(!is.null(hits())))
  })
  
  hit <- reactive({
    if (is.null(hits())) { return(NULL) }
    if (input$mode == "SimSearch" ) {
      if (is.null(qryA())) { return(NULL) }
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
    gene[ENSG==hit(), symbol]
  })
  hit_chr <- reactive({
    if (is.null(hit())) { return(NULL) }
    gene[ENSG==hit(), chr]
  })

  output$result_htm <- reactive({
    htm <- sprintf("<B>Results (%s):</B>", input$mode)
    if (!is.null(qryA())) {
      htm <- paste0(htm, sprintf(" GeneA: %s \"%s\"%s", qryA(), gene[ENSG==ensgA(), name],
	ifelse(("IDG" %in% input$opts), sprintf(" <A HREF=\"https://pharos.nih.gov/idg/targets/%s\" target=\"_blank\">%s</A>", gene[ENSG==ensgA(), uniprot], icon("external-link",lib="font-awesome")), "")))
      if (!is.null(chrA()) & grepl("^[XY]", chrA())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrA())) }
    }
    if (!is.null(qryB())) {
      htm <- paste0(htm, sprintf("; geneB: %s \"%s\"%s", qryB(), gene[ENSG==ensgB(), name],
	ifelse(("IDG" %in% input$opts), sprintf(" <A HREF=\"https://pharos.nih.gov/idg/targets/%s\" target=\"_blank\">%s</A>", gene[ENSG==ensgB(), uniprot], icon("external-link",lib="font-awesome")), "")))
      if (!is.null(chrB()) & grepl("^[XY]", chrB())) { htm <- paste0(htm, sprintf(", sex-linked, location %s", chrB())) }
    }
    if (!is.null(qryRex())) {
      htm <- paste0(htm, sprintf(" QueryText: \"%s\"", input$qryTxt))
    }
    if (grepl("Search", input$mode)) {
      htm <- paste0(htm, sprintf("; found: %d", ifelse(!is.null(hits()), nrow(hits()), 0)))
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

  hits_export <- reactive({
    if (is.null(hits())) { return(NULL) }
    hits_out <- hits()
    if (input$mode=="SimSearch") {
      hits_out["Query"] <- qryA()
      hits_out <- hits_out[, c(ncol(hits_out), 1:(ncol(hits_out)-1))] #Query col 1st.
      setnames(hits_out, c("Query", "EnsemblID", "UniProt", "GeneSymbol", "GeneName", "Chr",  "idgTDL", "idgDTO", "Group", input$score))
      return(hits_out)
    } else if (input$mode=="TxtSearch") {
      setnames(hits_out, c("EnsemblID", "UniProt", "GeneSymbol", "GeneName", "Chr", "idgTDL", "idgDTO"))
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
      titletxt = sprintf("GTEx Gene-Tissue Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), qryB(), gene[ENSG==ensgA(), name], gene[ENSG==ensgB(),name])
    } else if (input$mode=="SimSearch" & !is.null(hit())) {
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s <BR>\"%s\"", qryA(), gene[ENSG==ensgA(), name])
    } else if (input$mode=="TxtSearch" & !is.null(hits())) {
      titletxt = sprintf("Ex-files, GTEx-Profiles: \"%s\"", input$qryTxt)
    } else { #View
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s<BR>\"%s\"", qryA(), gene[ENSG==ensgA(), name])
    }
    return(titletxt)
  })
    
  output$plot <- renderPlotly({
    if (is.null(ensgA()) & is.null(hits())) { return(NULL) }
    if (length(input$groups)==0) { return(NULL) }

    if (!is.null(ensgA())) {
      smtsds <- tissue$SMTSD
      qryA_profile_f <- as.numeric(eps[ENSG==ensgA() & SEX=="F"][1, ..smtsds])
      qryA_profile_m <- as.numeric(eps[ENSG==ensgA() & SEX=="M"][1, ..smtsds])
      rhoAfm <- wPearson(qryA_profile_f, qryA_profile_m)
      ruzAfm <- Ruzicka(qryA_profile_f, qryA_profile_m)
      qryA_profile_n <- (qryA_profile_f+qryA_profile_m)/2 #Non-sexed (F+M)/2
    }

    if (input$mode=="Compare" & !is.null(qryB())) {
      smtsds <- tissue$SMTSD
      qryB_profile_f <- as.numeric(eps[ENSG==ensgB() & SEX=="F"][1, ..smtsds])
      qryB_profile_m <- as.numeric(eps[ENSG==ensgB() & SEX=="M"][1, ..smtsds])
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
      if ("C" %in% input$groups) {
        p <- add_trace(p, name = paste("(C)", qryA()), x = tissue$SMTSD, y = EpLogIf(qryA_profile_n, ("LogY" %in% input$opts)),
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
      if ("C" %in% input$groups) {
        p <- add_trace(p, name = paste("(C)", qryB()), x = tissue$SMTSD, y = EpLogIf(qryB_profile_n, ("LogY" %in% input$opts)),
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
      smtsds <- tissue$SMTSD
      if (!is.null(rows_selected)) {
        for (i in rows_selected) {
          hit_profile_f <- as.numeric(eps[ENSG==hits()$EnsemblID[i] & SEX=="F"][1, ..smtsds])
          hit_profile_m <- as.numeric(eps[ENSG==hits()$EnsemblID[i] & SEX=="M"][1, ..smtsds])
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
          if ("C" %in% input$groups) {
            hit_profile_n <- (hit_profile_f + hit_profile_m)/2
            p <- add_trace(p, name = paste("(C)", hits()$symbol[i]), x = tissue$SMTSD, y = EpLogIf(hit_profile_n, ("LogY" %in% input$opts)),
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
