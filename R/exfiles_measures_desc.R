library(readr)
library(plotly)
#
plots <- list()
#
MAX_ANTICOR <- -0.5
MIN_COR <- 0.5
MIN_SIM <- 0.5
#
colorScale <- data.frame(z=(0:(length(cols)-1)), col=c("#999999","#ff3333","#33ff33"))
###
ggc_cor <- read_delim("/home/data/GTEx/exfiles/gtex_rnaseq_profiles_WPearson_sample.tsv", "\t")
cor_tag <- "wRho"
#
n_cor_all <- nrow(ggc_cor)
n_cor_ok <- sum(((ggc_cor[[cor_tag]]>MIN_COR) | (ggc_cor[[cor_tag]]<MAX_ANTICOR)), na.rm=T)
#
writeLines(sprintf("Pct NA: %d / %d (%.1f%%)", sum(is.na(ggc_cor[[cor_tag]])), nrow(ggc_cor), 100*sum(is.na(ggc_cor[[cor_tag]]))/nrow(ggc_cor)))
ggc_cor <- ggc_cor[!is.na(ggc_cor[[cor_tag]]),]
qtl <- quantile(ggc_cor[[cor_tag]], c(seq(0,.09,.01), seq(.1,.9,.1), seq(.91,1,.01)))
writeLines(sprintf("%4s-ILE: %6.3f", names(qtl), qtl))
#
plots[[cor_tag]] <- plot_ly(type="histogram", x=ggc_cor[[cor_tag]], name=cor_tag,
		marker=list(color=c(rep(0,50),rep(1,100), rep(2,50)), cmin=0, cmax=2, colorscale=colorScale), 
		xbins = list(end=1.0, start = -1.0, size = 0.01)) %>%
	layout(title=sprintf("Ex-files correlation: %s",cor_tag), showlegend=F,
		margin=list(t=100), font=list(size=16, titlefont=list(size=22)))
#
###
ggc_sim <- read_delim("/home/data/GTEx/exfiles/gtex_rnaseq_profiles_Ruzicka_sample.tsv", "\t")
sim_tag <- "Ruzicka"
#
n_sim_all <- nrow(ggc_sim)
n_sim_ok <- sum((ggc_sim[[sim_tag]]>MIN_SIM), na.rm=T)
#
writeLines(sprintf("Pct NA: %d / %d (%.1f%%)", sum(is.na(ggc_sim[[sim_tag]])), nrow(ggc_sim), 100*sum(is.na(ggc_sim[[sim_tag]]))/nrow(ggc_sim)))
ggc_sim <- ggc_sim[!is.na(ggc_sim[[sim_tag]]),]
qtl <- quantile(ggc_sim[[sim_tag]], c(seq(0,.09,.01), seq(.1,.9,.1), seq(.91,1,.01)))
writeLines(sprintf("%4s-ILE: %6.3f", names(qtl), qtl))
#
plots[[sim_tag]] <- plot_ly(type="histogram", x=ggc_sim[[sim_tag]], name=sim_tag,
		marker=list(color=c(rep(1,50),rep(2,50)), cmin=0, cmax=2, colorscale=colorScale), 
		xbins = list(end=1.0, start = 0.0, size = 0.01)                            ) %>%
	layout(title=sprintf("Ex-files similarity: %s",sim_tag), showlegend=F,
		margin=list(t=100), font=list(size=16, titlefont=list(size=22)))
###
anns <- c(sprintf("N_ok / N_all = %d / %d (%.1f%%)", n_cor_ok, n_cor_all,  100*n_cor_ok/n_cor_all), 
	sprintf("N_ok / N_all = %d / %d (%.1f%%)", n_sim_ok, n_sim_all,  100*n_sim_ok/n_sim_all))
###
p <- subplot(plots, nrows=2, margin=0.05, shareX=F, shareY=F, titleX=F, titleY=F) %>%
	layout(title="Ex-files: Correlation and Similarity Measures<br><i>(Based on random sampling of computed values.)</i>", margin=list(t=80, b=60, l=30),
		font=list(family="Arial", size=14), showlegend=F) %>%
	add_annotations(text=names(plots), x=c(.1,.1), y=c(.8,.3), xref="paper", yref=
		"paper", align="center", font=list(family="Arial", size=20), showarrow=F) %>%
	add_annotations(text=anns, x=c(.5,.5), y=c(.7,.2), xref="paper", yref=
		"paper", align="center", font=list(family="Arial", size=14), showarrow=F)
#
p
###
