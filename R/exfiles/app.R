##########################################################################################
### See exfiles_prep.R to generate datafiles.
##########################################################################################
### Modes:
###      VIEW: geneA
###    SEARCH: geneA similarity search
###   COMPARE: geneA vs. geneB
### SABV: include sex-as-biological-variable analysis
### Dissim: dissimilarity search instead
##########################################################################################
### PROBLEM: Memory limit 1G on shinyapps.io for Starter subscriber.  This app exceeding
### limit with additional metrics.
##########################################################################################
### Jeremy Yang
##########################################################################################
library(readr)
library(wCorr)
library(shiny, quietly = T)
library(DT, quietly = T) #datatable
library(dplyr, quietly = T)
library(plotly, quietly = T)

#############################################################################
### This code runs once for all sessions.
### GGC = gene-gene associations
### GGC cols: "Gi","Gj","Cluster","rho","RMSD","ABC","wrho","tmoto"
APPNAME <- "Ex-files"
###
# Tissues now in eps column header.
#tissues <- read_delim("exfiles_tissues.csv", ",") # tissue metadata
###
gene <- read_delim("exfiles_gene.csv", ",") # gene attributes
eps <- read_delim("exfiles_eps.csv", ",") # expression profiles
###
tissues <- data.frame(id = 1:(ncol(eps)-4), name = colnames(eps)[5:ncol(eps)])
n_tissues <- nrow(tissues)
###
# Select cols to limit memory.  Unfortunately cols_only() seems to fail at that.
# "Gi", "Gj", "Cluster", "rho", "RMSD", "ABC", "wrho", "tmoto"
###
#ggc <- read_delim("exfiles_ggc.csv.gz", ",", col_types="cccddddd")
ggc <- read_delim("exfiles_ggc.csv.gz", ",", col_types=cols_only(Gi="c", Gj="c", Cluster="c", wrho="d", ABC="d"))
###
#
qryArand <- sample(union(ggc$Gi, ggc$Gj), 1) # initial random query
#
#Area_Between_Curves_(ABC)
ABC <- function(A,B) {
  abc <- 0
  if (length(A)!=length(B)) { return(abc); }
  for (i in 1:(length(A)-1)) {
    Amid <- mean(A[i],A[i+1])
    Bmid <- mean(B[i],B[i+1])
    if (A[i]>=B[i]) {
      if (A[i+1]>=B[i+1]) {
        abc <- abc + (AULS(A[i], A[i+1], 1) - AULS(B[i], B[i+1], 1))
      } else {
        abc <- abc + (AULS(A[i], Amid, .5) - AULS(B[i], Bmid, .5))
        abc <- abc + (AULS(Bmid, B[i+1], .5) - AULS(Amid, A[i+1], .5))
      }
    } else { 
      if (A[i+1]<B[i+1]) {
        abc <- abc + (AULS(B[i], B[i+1], 1) - AULS(A[i], A[i+1], 1))
      } else {
        abc <- abc + (AULS(B[i], Bmid, .5) - AULS(A[i], Amid, .5))
        abc <- abc + (AULS(Amid, A[i+1], .5) - AULS(Bmid, B[i+1], .5))
      }
    }
  }
  return(abc)
}
ABC_sim <- function(A,B) {
  n <- length(A)
  abc <- ABC(A,B)
  sim <- (1 / (1 + abc/n))
  return(sim)
}
#Area_Under_Line_Segment_(AULS)
AULS <- function(y1, y2, w) {
  auls <- min(y1,y2) * w
  auls <- auls + 0.5 * w * abs(y1-y2)
  return(auls)
}
#
RMSD <- function(A,B) {
  sqrt(sum((A-B)^2)/length(A))
}
Tanimoto <- function(A,B) {
  (A %*% B) / (A %*% A + B %*% B - A %*% B)
}
### Unweight smaller, noise-dominated expression values.
wPearson <- function(A,B) {
  wCorr::weightedCorr(A, B, method="Pearson", weights=(A+B)/2)
  #boot::corr(matrix(c(A,B), ncol=2), w=(A+B)/2)
}
###
# Create list of choices for input autocomplete.
gene_choices <- list()
for (i in 1:nrow(gene))
{
  gene_choices[[paste(gene$sym[i], gene$name[i])]] <- gene$sym[i]
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
<LI><B>ABC</B> - Area Between the Curves, a novel metric defined by the absolute value of the difference in the profile plots.
<LI><B>ABC_sim</B> - ABC transformed to [0-1] similarity by formula (1 / (1 + ABC/N_tissues))
<LI><B>wRho</B> - Weighted Pearson correlation coefficient, weighted by average values for each tissue, to mitigate noise.
<LI><B>Combo</B> - Product wRho*Tanimoto, balancing correlation and similarity.
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
<B>Algorithms:</B> Giovanni Bocci, Oleg Ursu, Cristian Bologa, Steve Mathias, Jeremy Yang &amp; Tudor Oprea<BR/>
<B>Web app:</B> Jeremy Yang<BR/>
Data from <A HREF=\"https://www.gtexportal.org/\" TARGET=\"_blank\">GTEx, The Genotype-Tissue Expression Project</A>.<BR/>
Built with R-Shiny &amp; Plotly.<BR/>
")}
#
#############################################################################
ui <- fluidPage(
  titlePanel(h2(sprintf("%s, GTEx expression-profile exploration", APPNAME), em("(PROTOTYPE)")), windowTitle=APPNAME),
  fluidRow(
    column(4, 
        wellPanel(
          actionButton("randGene", "GeneA (click for random)", style='padding:4px; background-color:#DDDDDD; font-size:100%; font-weight:bold'),
          selectizeInput("qryA", label=NULL, choices = gene_choices, selected=qryArand),
          selectizeInput("qryB", label="GeneB (optional)", choices = c(list('None'='none'), gene_choices)),
          #selectizeInput("qryB", label="GeneB (optional)", choices = NULL), #server-side not working
          radioButtons("mode", "Mode", choices=c("View", "Search", "Compare"), selected="Search", inline=T),
          checkboxInput("sabv", span("SABV", icon("venus",lib="font-awesome"),icon("mars", lib="font-awesome")), value=T),
          radioButtons("metric", "Metric", choices=c("ABC", "wRho", "Combo"), selected="Combo", inline=T),
          checkboxInput("dissim", "Dissimilarity", value=F),
#          checkboxGroupInput("searchgroups", "Searchgroups", choiceValues=list("MF","F","M"), 
#                             choiceNames = list(
#                               span("MF", icon("mars",lib="font-awesome"),icon("venus", lib="font-awesome")),
#                               span("F", icon("venus",lib="font-awesome")),
#                               span("M", icon("mars", lib="font-awesome"))
#                             ),
#                             selected=c("MF","F","M"), inline=T),
#          div(style="display:inline-block;vertical-align:top;", sliderInput("simcutoff", "Sim_cutoff", 0, 1, .8, step=.1, width="120px")),
#          div(style="display:inline-block;vertical-align:top;", sliderInput("nmax", "Max_hits", 0, 100, 10, step=10, width="120px")),
          br(),
          actionButton("showhelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
          )),
    column(8, plotlyOutput("plot", height = "580px"))
  ),
  fluidRow(
    column(12, wellPanel(
        htmlOutput(outputId = "result_htm", height = "60px")))
  ),
  conditionalPanel(condition = "input.mode == 'Search'",
    wellPanel(
      fluidRow(column(12, DT::dataTableOutput("datarows"))),
      fluidRow(column(12, downloadButton("hits_file", label="Download")))),
    width=12),
  hr(),
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
    message(sprintf("NOTE: mode: %s ; sabv: %s ; metric: %s ; dissim: %s", input$mode, input$sabv, input$metric, input$dissim))
    message(sprintf("NOTE: qryA = %s \"%s\"", qryA(), gene$name[gene$sym==qryA()]))
    message(sprintf("NOTE: qryB = %s \"%s\"", qryB(), ifelse(qryB()=="NONE", "(None)", gene$name[gene$sym==qryB()])))
  })
  
  qryA <- reactive({
    if (input$randGene>randGeneA_previous) {
      randGeneA_previous <<- input$randGene # Must assign to up-scoped variable.
      qryArand <- sample(union(ggc$Gi, ggc$Gj), 1)
      updateTextInput(session, "qryA", value=qryArand) #Works better than updateSelectizeInput ??
    }
    toupper(input$qryA)
  })

  qryB <- reactive({
    toupper(input$qryB)
  })
  
  hits <- reactive({
    if (input$mode != "Search") { return(NA) }
    genehits <- ggc[ggc$Gi == qryA() | ggc$Gj == qryA(),]
    
    #if (length(intersect(genehits$Cluster,input$searchgroups))>0) {
    #  genehits <- genehits[genehits$Cluster %in% input$searchgroups,]
    #}

    #genehits$Combo <- round(genehits$wrho*genehits$tmoto, digits=3)
    genehits$Combo <- round(genehits$wrho/(1 + genehits$ABC/n_tissues), digits=3)

    if (input$metric=="wRho") {
      genehits$sim <- genehits$wrho
    } else if (input$metric=="ABC") {
      genehits$sim <- (1 / (1 + genehits$ABC/n_tissues))
    } else if (input$metric=="Combo") {
      genehits$sim <- genehits$Combo
    }
    
    genehits$sim <- round(genehits$sim, digits=3)
    
    genehits$Gene <- NA
    genehits$Gene[genehits$Gi == qryA()] <- genehits$Gj[genehits$Gi == qryA()]
    genehits$Gene[genehits$Gj == qryA()] <- genehits$Gi[genehits$Gj == qryA()]
    
    genehits <- merge(genehits, gene[,c("sym","name")], by.x="Gene", by.y="sym", all.x=T, all.y=F)

    genehits <- genehits[,c("Gene","name","Cluster","sim")] #remove and reorder cols
    
    if (input$dissim) {
      genehits <- genehits[order(genehits$sim),]
    } else {
      genehits <- genehits[order(-genehits$sim),]
    }
    
    return(genehits)
  })
  
  hit <- reactive({
    if (qryA()=="" | input$mode != "Search") { return(NA) } # NA vs. NULL??
    hit_best <- hits()$Gene[1]
    sim <- hits()$sim[1]
    if (!is.na(hit_best)) {
      message(sprintf("NOTE: best hit [sim=%.2f]: %s (%s)", sim, hit_best, gene$name[gene$sym==hit_best]))
    } else {
      message(sprintf("ERROR: search failed."))
    }
    hit_best
  })
  
  output$result_htm <- reactive({
    if (input$mode == "View") {
      htm <- sprintf("<B>Results (%s):</B> Gene query: %s \"%s\"", input$mode, qryA(), gene$name[gene$sym==qryA()])
    } else if (input$mode == "Compare") {
      htm <- sprintf("<B>Results (%s):</B> Gene queryA: %s \"%s\" ; queryB: %s \"%s\"", input$mode, qryA(), gene$name[gene$sym==qryA()], qryB(), gene$name[gene$sym==qryB()])
    } else if (input$mode == "Search") {
      htm <- sprintf("<B>Results (%s):</B> Gene query: %s \"%s\" ; metric: %s ; profiles found: %d ; top hit: %s \"%s\"", input$mode, qryA(), gene$name[gene$sym==qryA()], input$metric, nrow(hits()), hit(), gene$name[gene$sym==hit()])
      sim <- hits()$sim[1]
      #if ((input$dissim & sim>input$simcutoff) | (!input$dissim & sim<input$simcutoff)) {
      #  htm <- paste(htm, sprintf("<B>WARNING:</B> top hit Sim (%.2f) does not satisfy Sim_cutoff (%.2f).", sim, input$simcutoff))
      #}
      #searchgroup <- hits()$Cluster[1]
      #if (!(searchgroup %in% input$searchgroups)) {
      #  htm <- paste(htm, sprintf("<B>WARNING:</B> top hit searchgroup: %s ; none in specified searchgroups: %s.", searchgroup, paste(collapse=", ", input$searchgroups)))
      #}
    } else {
      htm <- sprintf("<B>Results:</B> ERROR: invalid mode: %s", input$mode)
    }
    if (input$sabv) {
      htm <- paste(htm, "; SABV analysis ON.")
    }
    return(htm)
  })
  
  #Without-DT method:
  #output$datarows <- renderDataTable(hits(), options = list(pageLength=min(10,nrow(hits())), lengthChange=F, pagingType="simple", searching=F, info=T))
  
#  observe({
#      rows_selected <- input$datarows_rows_selected
#      if (!is.null(rows_selected)) {
#        message(sprintf("DEBUG: input$datarows_rows_selected: %s (%s)", 
#                        paste0(rows_selected, collapse=","),
#                        paste0(hits()$Gene[rows_selected], collapse=",")))
#      }
#  })

  ### Creates input$datarows_rows_selected
  output$datarows <- renderDataTable({
    dt = hits()
    DT::datatable(data = dt,
	rownames = F,
	selection = "multiple",
	class = "cell-border stripe",
	style = "bootstrap",
	options = list(autoWidth=T),
	colnames = c("Gene", "Name", "Group", "Similarity")) %>%
        formatRound(digits=3, columns=4:ncol(dt))
  }, server=T)
  
  hits_export <- reactive({
    genehits <- hits()
    genehits$Query <- qryA()
    genehits <- genehits[,c(5,1,2,3,4)] #reorder cols
    names(genehits) <- c("Query", "GeneSym", "GeneName", "Group", paste0(input$metric, "_Similarity"))
    return(genehits)
  })
  
  output$hits_file <- downloadHandler(
    filename = function() { "exfiles_hits.csv" },
    content = function(file) { write.csv(hits_export(), file, row.names=F) }
  )

  output$plot <- renderPlotly({
    
    qryA_profile_m <- as.numeric(eps[eps$gene==qryA() & eps$sex=="M",][1,5:44])
    qryA_profile_f <- as.numeric(eps[eps$gene==qryA() & eps$sex=="F",][1,5:44])
    qryA_profile <- (qryA_profile_m + qryA_profile_f)/2

    rhoAfm <- cor(qryA_profile_f, qryA_profile_m, method = "spearman")
    wrhoAfm <- wPearson(qryA_profile_f, qryA_profile_m )
    
    #message(sprintf("DEBUG: qryA_profile_f = %s", paste0(as.character(qryA_profile_f), collapse=", ")))
    #message(sprintf("DEBUG: qryA_profile_m = %s", paste0(as.character(qryA_profile_m), collapse=", ")))
    #message(sprintf("DEBUG: qryA_profile = %s", paste0(as.character(qryA_profile), collapse=", ")))
    
    abcAfm <- ABC_sim(qryA_profile_f, qryA_profile_m)
    
    if (input$mode == "Compare") {
      qryB_profile_m <- as.numeric(eps[eps$gene==qryB() & eps$sex=="M",][1,5:44])
      qryB_profile_f <- as.numeric(eps[eps$gene==qryB() & eps$sex=="F",][1,5:44])
      qryB_profile <- (qryB_profile_m + qryB_profile_f)/2
      rho <- cor(qryA_profile, qryB_profile, method = "spearman")
      wrho <- wPearson(qryA_profile, qryB_profile)

      abc <- ABC_sim(qryA_profile, qryB_profile)

      rhoBfm <- cor(qryB_profile_f, qryB_profile_m, method = "spearman")
      rhoFab <- cor(qryA_profile_f, qryB_profile_f, method = "spearman")
      rhoMab <- cor(qryA_profile_m, qryB_profile_m, method = "spearman")

      wrhoBfm <- wPearson(qryB_profile_f, qryB_profile_m)
      wrhoFab <- wPearson(qryA_profile_f, qryB_profile_f)
      wrhoMab <- wPearson(qryA_profile_m, qryB_profile_m)

      abcBfm <- ABC_sim(qryB_profile_f, qryB_profile_m)
      abcFab <- ABC_sim(qryA_profile_f, qryB_profile_f)
      abcMab <- ABC_sim(qryA_profile_m, qryB_profile_m)

    } else if (input$mode == "Search") { ## Hit only if Search
      hit_profile_m <- as.numeric(eps[eps$gene==hit() & eps$sex=="M",][1,5:44])
      hit_profile_f <- as.numeric(eps[eps$gene==hit() & eps$sex=="F",][1,5:44])
      hit_profile <- (hit_profile_m + hit_profile_f)/2
      rho <- cor(qryA_profile, hit_profile, method = "spearman")
      wrho <- wPearson(qryA_profile, hit_profile)

      #
      abc <- ABC_sim(qryA_profile, hit_profile)

      ### "B" = hit
      rhoBfm <- cor(hit_profile_f, hit_profile_m, method = "spearman")
      rhoFab <- cor(qryA_profile_f, hit_profile_f, method = "spearman")
      rhoMab <- cor(qryA_profile_m, hit_profile_m, method = "spearman")

      wrhoBfm <- wPearson(hit_profile_f, hit_profile_m)
      wrhoFab <- wPearson(qryA_profile_f, hit_profile_f)
      wrhoMab <- wPearson(qryA_profile_m, hit_profile_m)

      abcBfm <- ABC_sim(hit_profile_f, hit_profile_m)
      abcFab <- ABC_sim(qryA_profile_f, hit_profile_f)
      abcMab <- ABC_sim(qryA_profile_m, hit_profile_m)

    } #Else: View mode
      

    ### PLOT:
    xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10), categoryorder = "array", categoryarray = tissues$name)
    yaxis = list(title="Expression: LOG<SUB>10</SUB>(1 + TPM)")

    if (input$mode == "View") {
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s<BR>\"%s\"", qryA(), gene$name[gene$sym==qryA()])
    } else if (input$mode == "Search") {
      titletxt = sprintf("Ex-files, GTEx-Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), hit(), gene$name[gene$sym==qryA()], gene$name[gene$sym==hit()])
    } else if (input$mode == "Compare") {
      titletxt = sprintf("GTEx Gene-Tissue Profiles: %s vs %s<BR>\"%s\" vs \"%s\"", qryA(), qryB(), gene$name[gene$sym==qryA()], gene$name[gene$sym==qryB()])
    }

    p <- plot_ly() %>%
      layout(xaxis = xaxis, yaxis = yaxis, 
         title = titletxt,
         margin = list(t=100, r=80, b=160, l=60),
         legend = list(x=.9, y=1),
         font = list(family="Arial", size=14)
      )
    # %>% add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=.2, xref="paper", yref="paper")
    #
    if (!input$sabv) {
      p <-  add_trace(p, name = qryA(), x = tissues$name, y = qryA_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$name))
      if (input$mode == "View") {
        annotxt <- ""
      } else if (input$mode == "Compare") {
        annotxt <- sprintf("rho = %.2f; sim = %.2f", rho, abc)
        p <- add_trace(p, name = qryB(), x = tissues$name, y = qryB_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$name))
      } else if (input$mode == "Search") {
        annotxt <- sprintf("rho = %.2f; sim = %.2f", rho, abc)
        p <- add_trace(p, name = hit(), x = tissues$name, y = hit_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$name))
        ###
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          gsyms <- hits()$Gene[rows_selected]
          for (gsym in gsyms) {
            if (gsym == hit()) { next; }
            profile_m <- as.numeric(eps[eps$gene==gsym & eps$sex=="M",][1,5:44])
            profile_f <- as.numeric(eps[eps$gene==gsym & eps$sex=="F",][1,5:44])
            profile <- (profile_f + profile_m)/2
            p <- add_trace(p, name = gsym, x = tissues$name, y = profile,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(gsym, ": ", tissues$name))
          }
        }
      }
    } else if (input$sabv) {
      p <-  add_trace(p, name = paste("(F)", qryA()), x = tissues$name, y = qryA_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$name))
    #
      p <- add_trace(p, name = paste("(M)", qryA()), x = tissues$name, y = qryA_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$name))
    #
      if (input$mode == "View") {
        annotxt <- sprintf("Afm: rho = %.2f; sim = %.2f", wrhoAfm, abcAfm)
      } else if (input$mode == "Search") {
        annotxt <- paste(sep="<BR>",
            sprintf("Uab: rho = %.2f; sim = %.2f", wrho, abc),
            sprintf("Fab: rho = %.2f; sim = %.2f", wrhoFab, abcFab),
            sprintf("Mab: rho = %.2f; sim = %.2f", wrhoMab, abcMab),
            sprintf("Afm: rho = %.2f; sim = %.2f", wrhoAfm, abcAfm),
            sprintf("Bfm: rho = %.2f; sim = %.2f", wrhoBfm, abcBfm))
        p <- add_trace(p, name = paste("(F)", hit()), x = tissues$name, y = hit_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$name))
        p <- add_trace(p, name = paste("(M)", hit()), x = tissues$name, y = hit_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$name))
        ##
        # Include genes selected via interactive table.
        rows_selected <- input$datarows_rows_selected
        if (!is.null(rows_selected)) {
          gsyms <- hits()$Gene[rows_selected]
          for (gsym in gsyms) {
            if (gsym == hit()) { next; }
            profile_m <- as.numeric(eps[eps$gene==gsym & eps$sex=="M",][1,5:44])
            profile_f <- as.numeric(eps[eps$gene==gsym & eps$sex=="F",][1,5:44])
            p <- add_trace(p, name = paste("(M)", gsym), x = tissues$name, y = profile_m,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(gsym, ": ", tissues$name))
            p <- add_trace(p, name = paste("(F)", gsym), x = tissues$name, y = profile_f,
                type = 'scatter', mode = 'lines+markers',
                marker = list(symbol="circle", size=10),
                text = paste0(gsym, ": ", tissues$name))
          }
        }
      } else if (input$mode == "Compare") {
        annotxt <- paste(sep="<BR>",
            sprintf("Uab: rho = %.2f; sim = %.2f", wrho, abc),
            sprintf("Fab: rho = %.2f; sim = %.2f", wrhoFab, abcFab),
            sprintf("Mab: rho = %.2f; sim = %.2f", wrhoMab, abcMab),
            sprintf("Afm: rho = %.2f; sim = %.2f", wrhoAfm, abcAfm),
            sprintf("Bfm: rho = %.2f; sim = %.2f", wrhoBfm, abcBfm))
        p <- add_trace(p, name = paste("(F)", qryB()), x = tissues$name, y = qryB_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$name))
        #
        p <- add_trace(p, name = paste("(M)", qryB()), x = tissues$name, y = qryB_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$name))
      }
    } else { return(NULL) } #ERROR
    #
    #p <- add_annotations(p, text=annotxt, showarrow=F, x=.1, y=1, xref="paper", yref="paper")
    p$elementId <- NULL #Hack to suppress spurious warnings.
    return(p)
  })
  
}

###
shinyApp(ui, server)
