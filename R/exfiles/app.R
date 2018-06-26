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
### Maybe: Search for SABV-dissimilar?
##########################################################################################
### PROBLEM: Memory limit 1G on shinyapps.io for Starter subscriber.  This app exceeding
### limit with additional metrics.
##########################################################################################
### Jeremy Yang
##########################################################################################
library(readr)
library(wCorr)
#library(boot)
library(shiny, quietly = T)
library(DT, quietly = T) #datatable
library(dplyr, quietly = T)
library(plotly, quietly = T)

#############################################################################
### This code runs once for all sessions.
### GGC = gene-gene associations
### GGC cols: "Gi","Gj","Cluster","rho","RMSD","ABC","wrho","tmoto"
APPNAME <- "Ex-files"
gene <- read_delim("exfiles_gene.csv", ",") # gene attributes
eps <- read_delim("exfiles_eps.csv", ",") # expression profiles
tissues <- read_delim("exfiles_tissues.csv", ",") # tissue metadata
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
###
ABC <- function(A,B) {
  abc <- 0
  for (i in 1:(length(A)-1)) {
    Amid <- mean(A[i],A[i+1])
    Bmid <- mean(B[i],B[i+1])
    if (A[i]>=B[i]) {
      if (A[i+1]>=B[i+1]) {
        abc <- abc + (AreaUnderLineSegment(A[i], A[i+1], 1) - AreaUnderLineSegment(B[i], B[i+1], 1))
      } else {
        abc <- abc + (AreaUnderLineSegment(A[i], Amid, .5) - AreaUnderLineSegment(B[i], Bmid, .5))
        abc <- abc + (AreaUnderLineSegment(Bmid, B[i+1], .5) - AreaUnderLineSegment(Amid, A[i+1], .5))
      }
    } else { 
      if (A[i+1]<B[i+1]) {
        abc <- abc + (AreaUnderLineSegment(B[i], B[i+1], 1) - AreaUnderLineSegment(A[i], A[i+1], 1))
      } else {
        abc <- abc + (AreaUnderLineSegment(B[i], Bmid, .5) - AreaUnderLineSegment(A[i], Amid, .5))
        abc <- abc + (AreaUnderLineSegment(Amid, A[i+1], .5) - AreaUnderLineSegment(Bmid, B[i+1], .5))
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
#
AreaUnderLineSegment <- function(y1, y2, w) {
  a <- min(y1,y2) * w
  a <- a + 0.5 * w * abs(y1-y2)
  return(a)
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
"<P><B>Ex-files</B> allows exploration and analysis of genes represented by their expression profiles.
With <B>GTEx</B> as the data source, expression profiles are real valued vectors of expression levels across the
defined tissue types.</P>
<P>
<B>Inputs</B> are query geneA, and <I>optionally</I>, geneB.
Ex-files has two modes of operation:
<UL>
<LI><B>View</B> - view one gene: input geneA
<LI><B>Search</B> - search for genes based on input geneA and selected metric
<LI><B>Compare</B> - compare input geneA vs. geneB based on selected metric
</UL>
<B>SABV:</B> This checkbox invokes <B><I>Sex As a Biological Variable</I></B> analysis.</P>
<P>
<B>Metrics:</B>
<UL>
<LI><B>Tanimoto</B> - Similarity metric defined by dot products (AdotB)/(AdotA + BdotB - AdotB)
<LI><B>RMSD</B> - Root Mean Square Distance, a standard metric comparing two sets of points in a space of any dimensionality.  Each of the 
expression values in each profile is considered a point in a 1D space.
<LI><B>ABC</B> - Area Between the Curves, a novel metric defined by the absolute value of the difference in the profile plots.
<LI><B>wRho</B> - Weighted Pearson correlation coefficient, weighted by average values for each tissue, to mitigate noise.
<LI><B>Rho</B> - Spearman rank correlation coefficient.
<LI><B>Combo</B> - Product wRho*Tanimoto, balancing correlation and similarity.
</UL>
<B>Dissimilarity:</B> This checkbox affects Search results, reversing the chosen metric.</P>
<P>
<B>Searchgroups:</B>
<UL>
<LI><B>M</B> - search males-only 
<LI><B>F</B> - search females-only 
<LI><B>MF</B> - search both sexes
</UL></P>
<P>
<B>Results:</B>
<UL>
<LI>Search results include the top hit, displayed against the query in the
plot, and other hits displayed in a table and downloadable as CSV.
<LI>Compare results consist of the plot of GeneA vs. GeneB.
<LI>SABV results include both sexes for both genes, thus four profiles.
<LI>Expression units are <B>LOG<SUB>10</SUB>(1 + TPM)</B>, where<BR>
<B>TPM</B> = RNA-seq median Transcripts Per Million-kilobase
</UL></P>
<B>Abbrevations:</B>
<UL>
<LI><B>Uab</B> - Unisex geneA vs. geneB
<LI><B>Fab</B> - Female geneA vs. geneB
<LI><B>Mab</B> - Male geneA vs. geneB
<LI><B>Afm</B> - GeneA Female vs. Male
<LI><B>Bfm</B> - GeneB Female vs. Male
</UL>
<B>Algorithms:</B> Giovanni Bocci, Oleg Ursu, Cristian Bologa, Steve Mathias, Jeremy Yang &amp; Tudor Oprea<BR/>
<B>Web app:</B> Jeremy Yang<BR/>
Data from <A HREF=\"https://www.gtexportal.org/\">GTEx, The Genotype-Tissue Expression Project</A>.<BR/>
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
          radioButtons("mode", "Mode", choices=c("View", "Search", "Compare"), selected="Search", inline=T),
          checkboxInput("sabv", span("SABV", icon("venus",lib="font-awesome"),icon("mars", lib="font-awesome")), value=T),
          radioButtons("metric", "Metric", choices=c("ABC", "wRho", "Combo"), selected="ABC", inline=T),
          checkboxInput("dissim", "Dissimilarity", value=F),
          checkboxGroupInput("searchgroups", "Searchgroups", choiceValues=list("MF","F","M"), 
                             choiceNames = list(
                               span("MF", icon("mars",lib="font-awesome"),icon("venus", lib="font-awesome")),
                               span("F", icon("venus",lib="font-awesome")),
                               span("M", icon("mars", lib="font-awesome"))
                             ),
                             selected=c("MF"), inline=T),
          div(style="display:inline-block;vertical-align:top;", sliderInput("simcutoff", "Sim_cutoff", 0, 1, .8, step=.1, width="120px")),
          #div(style="display:inline-block;vertical-align:top;", sliderInput("nmax", "Max_hits", 0, 100, 10, step=10, width="120px")),
          br(),
          actionButton("showhelp", "Help", style='padding:4px; background-color:#DDDDDD; font-weight:bold')
          )),
    column(8, plotlyOutput("plot", height = "640px"))
  ),
  fluidRow(
    column(12, wellPanel(
        htmlOutput(outputId = "result_htm", height = "60px")))
  ),
  fluidRow(column(12, dataTableOutput("datarows"))),
  fluidRow(column(12, downloadButton("hits_file", label="Download"))),
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
  
  #updateSelectizeInput(session, 'qryB', choices = gene_choices, server=T) #server-side, NOT WORKING
  
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
  
  #CFAP44
  
  hits <- reactive({
    if (input$mode != "Search") { return(NA) }
    ggc_this <- ggc[ggc$Gi == qryA() | ggc$Gj == qryA(),]
    
    if (length(intersect(ggc_this$Cluster,input$searchgroups))>0) {
      ggc_this <- ggc_this[ggc_this$Cluster %in% input$searchgroups,]
    }

    #ggc_this$Combo <- round(ggc_this$wrho*ggc_this$tmoto, digits=3)
    ggc_this$Combo <- round(ggc_this$wrho/(1 + ggc_this$ABC/n_tissues), digits=3)

    if (input$metric=="wRho") {
      ggc_this$sim <- ggc_this$wrho
#    } else if (input$metric=="Rho") {
#      ggc_this$sim <- ggc_this$rho
#    } else if (input$metric=="RMSD") {
#      ggc_this$sim <- -ggc_this$RMSD
    } else if (input$metric=="ABC") {
      ggc_this$sim <- (1 / (1 + ggc_this$ABC/n_tissues))
#    } else if (input$metric=="Tanimoto") {
#      ggc_this$sim <- ggc_this$tmoto
    } else if (input$metric=="Combo") {
      ggc_this$sim <- ggc_this$Combo
    }
    
    ggc_this$sim <- round(ggc_this$sim, digits=3)
    
    ggc_this <- ggc_this[order(-ggc_this$sim),]
    
    if (input$dissim) df<- ggc_this[seq(dim(ggc_this)[1],1),] #reverse order
    
    # Only filter if some hits pass threshold.
    # Maybe should do similarity cutoff instead.
#    if (input$dissim) {
#      if (input$simcutoff>=min(ggc_this$sim)) {
#        ggc_this <- ggc_this[ggc_this$sim<=input$simcutoff,]
#      }
#    } else {
#      if (input$simcutoff<=max(ggc_this$sim)) {
#        ggc_this <- ggc_this[ggc_this$sim>=input$simcutoff,]
#      }
#    }
    #ggc_this <- ggc_this[1:input$nmax,]
    return(ggc_this)
  })
  
  hit <- reactive({
    if (qryA()=="" | input$mode != "Search") { return(NA) } # NA vs. NULL??
    hit_best <- ifelse(hits()$Gi[1]==qryA(), hits()$Gj[1], hits()$Gi[1])
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
      htm <- sprintf("<B>Results (%s):</B> Gene queryA: %s \"%s\"", input$mode, qryA(), gene$name[gene$sym==qryA()])
    } else if (input$mode == "Compare") {
      htm <- sprintf("<B>Results (%s):</B> Gene queryA: %s \"%s\" ; queryB: %s \"%s\" ; expression profiles found: %d", input$mode, qryA(), gene$name[gene$sym==qryA()], qryB(), gene$name[gene$sym==qryB()], nrow(hits()))
    } else if (input$mode == "Search") {
      htm <- sprintf("<B>Results (%s):</B> Gene query: %s \"%s\" ; metric: %s ; profiles found: %d ; top hit: %s \"%s\"", input$mode, qryA(), gene$name[gene$sym==qryA()], input$metric, nrow(hits()), hit(), gene$name[gene$sym==hit()])
      sim <- hits()$sim[1]
      if ((input$dissim & sim>input$simcutoff) | (!input$dissim & sim<input$simcutoff)) {
        htm <- paste(htm, sprintf("<B>WARNING:</B> top hit Sim (%.2f) does not satisfy Sim_cutoff (%.2f).", sim, input$simcutoff))
      }
      searchgroup <- hits()$Cluster[1]
      if (!(searchgroup %in% input$searchgroups)) {
        htm <- paste(htm, sprintf("<B>WARNING:</B> top hit searchgroup: %s ; none in specified searchgroups: %s.", searchgroup, paste(collapse=", ", input$searchgroups)))
      }
    } else {
      htm <- sprintf("<B>Results:</B> ERROR: invalid mode: %s", input$mode)
    }
    if (input$sabv) {
      htm <- paste(htm, "; SABV analysis ON.")
    }
    return(htm)
  })
  
  output$datarows <- renderDataTable(hits(), options = list(pageLength=min(10,nrow(hits())), lengthChange=F, pagingType="simple", searching=F, info=T))
  
  output$hits_file <- downloadHandler(
    filename = function() { paste0("genefile_hits-", Sys.Date(), ".csv") },
    content = function(file) { write.csv(hits(), file, row.names=F) }
  )

  output$plot <- renderPlotly({
    
    qryA_profile_m <- as.numeric(eps[eps$gene==qryA() & eps$sex=="M",][1,5:44])
    qryA_profile_f <- as.numeric(eps[eps$gene==qryA() & eps$sex=="F",][1,5:44])
    qryA_profile <- (qryA_profile_m + qryA_profile_f)/2

    rhoAfm <- cor(qryA_profile_f, qryA_profile_m, method = "spearman")
    wrhoAfm <- wPearson(qryA_profile_f, qryA_profile_m )

    abcAfm <- ABC_sim(qryA_profile_f, qryA_profile_m)
    #rmsdAfm <- RMSD(qryA_profile_f, qryA_profile_m)
    #tmotoAfm <- Tanimoto(qryA_profile_f, qryA_profile_m)
    
    if (input$mode == "Compare") {
      qryB_profile_m <- as.numeric(eps[eps$gene==qryB() & eps$sex=="M",][1,5:44])
      qryB_profile_f <- as.numeric(eps[eps$gene==qryB() & eps$sex=="F",][1,5:44])
      qryB_profile <- (qryB_profile_m + qryB_profile_f)/2
      rho <- cor(qryA_profile, qryB_profile, method = "spearman")
      wrho <- wPearson(qryA_profile, qryB_profile)

      abc <- ABC_sim(qryA_profile, qryB_profile)
      #rmsd <- RMSD(qryA_profile, qryB_profile)
      #tmoto <- Tanimoto(qryA_profile, qryB_profile)

      rhoBfm <- cor(qryB_profile_f, qryB_profile_m, method = "spearman")
      rhoFab <- cor(qryA_profile_f, qryB_profile_f, method = "spearman")
      rhoMab <- cor(qryA_profile_m, qryB_profile_m, method = "spearman")

      wrhoBfm <- wPearson(qryB_profile_f, qryB_profile_m)
      wrhoFab <- wPearson(qryA_profile_f, qryB_profile_f)
      wrhoMab <- wPearson(qryA_profile_m, qryB_profile_m)

      abcBfm <- ABC_sim(qryB_profile_f, qryB_profile_m)
      abcFab <- ABC_sim(qryA_profile_f, qryB_profile_f)
      abcMab <- ABC_sim(qryA_profile_m, qryB_profile_m)

      #rmsdBfm <- RMSD(qryB_profile_f, qryB_profile_m)
      #rmsdFab <- RMSD(qryA_profile_f, qryB_profile_f)
      #rmsdMab <- RMSD(qryA_profile_m, qryB_profile_m)

      #tmotoBfm <- Tanimoto(qryB_profile_f, qryB_profile_m)
      #tmotoFab <- Tanimoto(qryA_profile_f, qryB_profile_f)
      #tmotoMab <- Tanimoto(qryA_profile_m, qryB_profile_m)

    } else { ## Hit only if Search
      hit_profile_m <- as.numeric(eps[eps$gene==hit() & eps$sex=="M",][1,5:44])
      hit_profile_f <- as.numeric(eps[eps$gene==hit() & eps$sex=="F",][1,5:44])
      hit_profile <- (hit_profile_m + hit_profile_f)/2
      rho <- cor(qryA_profile, hit_profile, method = "spearman")
      wrho <- wPearson(qryA_profile, hit_profile)

      abc <- ABC_sim(qryA_profile, hit_profile)
      #rmsd <- RMSD(qryA_profile, hit_profile)
      #tmoto <- Tanimoto(qryA_profile, hit_profile)

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

      #rmsdBfm <- RMSD(hit_profile_f, hit_profile_m)
      #rmsdFab <- RMSD(qryA_profile_f, hit_profile_f)
      #rmsdMab <- RMSD(qryA_profile_m, hit_profile_m)

      #tmotoBfm <- Tanimoto(hit_profile_f, hit_profile_m)
      #tmotoFab <- Tanimoto(qryA_profile_f, hit_profile_f)
      #tmotoMab <- Tanimoto(qryA_profile_m, hit_profile_m)
    }

    ### PLOT:
    xaxis = list(tickangle=45, tickfont=list(family="Arial", size=10))
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
      ) %>%
      add_annotations(text=format(Sys.time(), "%Y-%m-%d"), showarrow=F, x=1.0, y=.2, xref="paper", yref="paper")
    #
    if (!input$sabv) {
      p <-  add_trace(p, name = qryA(), x = tissues$tissue_name, y = qryA_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$tissue_name))
      if (input$mode == "View") {
        annotxt <- ""
      } else if (input$mode == "Compare") {
        annotxt <- sprintf("rho = %.2f; rmsd = %.2f", rho, rmsd)
        p <- add_trace(p, name = qryB(), x = tissues$tissue_name, y = qryB_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$tissue_name))
      } else if (input$mode == "Search") {
        annotxt <- sprintf("rho = %.2f; rmsd = %.2f", rho, rmsd)
        p <- add_trace(p, name = hit(), x = tissues$tissue_name, y = hit_profile,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$tissue_name))
      }
    } else if (input$sabv) {
      p <-  add_trace(p, name = paste("(F)", qryA()), x = tissues$tissue_name, y = qryA_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$tissue_name))
    #
      p <- add_trace(p, name = paste("(M)", qryA()), x = tissues$tissue_name, y = qryA_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryA(), ": ", tissues$tissue_name))
    #
      if (input$mode == "View") {
        #annotxt <- sprintf("Afm: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoAfm, rhoAfm, tmotoAfm, rmsdAfm)
        annotxt <- sprintf("Afm: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoAfm, rhoAfm, abcAfm)
      } else if (input$mode == "Search") {
        annotxt <- paste(sep="<BR>",
            #sprintf("Uab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrho, rho, tmoto, rmsd),
            #sprintf("Fab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoFab, rhoFab, tmotoFab, rmsdFab),
            #sprintf("Mab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoMab, rhoMab, tmotoMab, rmsdMab),
            #sprintf("Afm: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoAfm, rhoAfm, tmotoAfm, rmsdAfm),
            #sprintf("Bfm: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoBfm, rhoBfm, tmotoBfm, rmsdBfm))
            sprintf("Uab: wrho = %.2f; rho = %.2f; abc = %.2f", wrho, rho, abc),
            sprintf("Fab: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoFab, rhoFab, abcFab),
            sprintf("Mab: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoMab, rhoMab, abcMab),
            sprintf("Afm: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoAfm, rhoAfm, abcAfm),
            sprintf("Bfm: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoBfm, rhoBfm, abcBfm))
        p <- add_trace(p, name = paste("(F)", hit()), x = tissues$tissue_name, y = hit_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$tissue_name))
        #
        p <- add_trace(p, name = paste("(M)", hit()), x = tissues$tissue_name, y = hit_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(hit(), ": ", tissues$tissue_name))
      } else if (input$mode == "Compare") {
        annotxt <- paste(sep="<BR>",
            sprintf("Uab: wrho = %.2f; rho = %.2f; abc = %.2f", wrho, rho, abc),
            sprintf("Fab: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoFab, rhoFab, abcFab),
            sprintf("Mab: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoMab, rhoMab, abcMab),
            sprintf("Afm: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoAfm, rhoAfm, abcAfm),
            sprintf("Bfm: wrho = %.2f; rho = %.2f; abc = %.2f", wrhoBfm, rhoBfm, abcBfm))
            #sprintf("Uab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrho, rho, tmoto, rmsd),
            #sprintf("Fab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoFab, rhoFab, tmotoFab, rmsdFab),
            #sprintf("Mab: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoMab, rhoMab, tmotoMab, rmsdMab),
            #sprintf("Afm: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoAfm, rhoAfm, tmotoAfm, rmsdAfm),
            #sprintf("Bfm: wrho = %.2f; rho = %.2f; tmoto = %.2f; rmsd = %.2f", wrhoBfm, rhoBfm, tmotoBfm, rmsdBfm))
        p <- add_trace(p, name = paste("(F)", qryB()), x = tissues$tissue_name, y = qryB_profile_f,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$tissue_name))
        #
        p <- add_trace(p, name = paste("(M)", qryB()), x = tissues$tissue_name, y = qryB_profile_m,
            type = 'scatter', mode = 'lines+markers',
            marker = list(symbol="circle", size=10),
            text = paste0(qryB(), ": ", tissues$tissue_name))
      }
    } else { return(NULL) } #ERROR
    #
    p <- add_annotations(p, text=annotxt, showarrow=F, x=.1, y=1, xref="paper", yref="paper")
    p$elementId <- NULL #Hack to suppress spurious warnings.
    return(p)
  })
  
}

###
shinyApp(ui, server)
