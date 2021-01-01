library(KEGGREST)
library(tidyverse)

N_KO_norm_cov_tab <- read.table("KO-normalized-master-tab.tsv", sep="\t", header=TRUE)

N_KO_info_tab <- read.table("All_N_KO_info.tsv", sep = "\t", header = TRUE)

get_KO_info <- function( target_KO, KO_KEGG_tab = N_KO_norm_cov_tab ) {

    # bulding link
    curr_KO_link <- paste0("https://www.genome.jp/dbget-bin/www_bget?ko:", target_KO)

    # getting info on KEGG term
    try(curr_info <- keggGet(target_KO), silent = TRUE)

    # checking if it was found at KEGG
    if ( ! exists("curr_info") ) {

        cat(paste0("\n  It seems '", target_KO, "' wasn't found at KEGG. It may have been removed, or may have never existed.\n\n  This should be the link if it were there if you wanna take a look:\n\n      ", curr_KO_link, "\n\n"))
        # couldn't figure out a better way to just report this and not give an error while exiting
        stop_no_error <- function() {
            opt <- options(show.error.messages = FALSE)
            on.exit(options(opt))
            stop()
        }

        stop_no_error()
    }

    # parsing some info
    if ( length(curr_info[[1]]$ENTRY) > 0 ) { curr_KO_ID <- curr_info[[1]]$ENTRY %>% as.vector } else { curr_KO_ID <- NA }
    if ( length(curr_info[[1]]$NAME) > 0 ) { curr_KO_name <- curr_info[[1]]$NAME %>% as.vector } else { curr_KO_name <- NA }
    if ( length(curr_info[[1]]$DEFINITION) > 0 ) { curr_KO_def <- curr_info[[1]]$DEFINITION %>% as.vector } else { curr_KO_def <- NA }

    # seeing if it's in our table
    if ( target_KO %in% (KO_KEGG_tab %>% pull(KO_ID)) ) {

      curr_in_data <- "Yes"

    } else {

      curr_in_data <- "No"

    }

    # reporting info
    cat("\n  KO ID          :  ", curr_KO_ID, "\n")
    cat("  KO name        :  ", curr_KO_name, "\n")
    cat("  KO definition  :  ", curr_KO_def, "\n")
    cat("  KO link        :  ", curr_KO_link, "\n")
    cat("  In our N data? :  ", curr_in_data, "\n\n")

}

plot_KO <- function( target_KO, tab = N_KO_norm_cov_tab ) {

    # making sure specified KO term is present/was assembled/annotated
    if ( ! target_KO %in% ( tab %>% pull(KO_ID) ) ) {

        cat("\n   ", target_KO, "was not detected in our N data.\n\n")

        # couldn't figure out a better way to just report this and not give an error while exiting
        stop_no_error <- function() {
            opt <- options(show.error.messages = FALSE)
            on.exit(options(opt))
            stop()
        }
        stop_no_error()
    }

    # subsetting to target KO
    sub_tab <- tab %>% filter(KO_ID == target_KO) %>% select(c(1:5, 10, 11, 12))

    # making long formatted version
    sub_long <- sub_tab %>% select(1:5) %>% pivot_longer(-KO_ID, names_to = "Depth", values_to = "norm_cov") %>% data.frame(check.names = FALSE)

    # making labels
    num_uniq_genes <- sub_tab$num_uniq_genes
    ko_name <- sub_tab$KO_name
    ko_def <- sub_tab$KO_def
    main_lab <- paste0(target_KO, " (", num_uniq_genes, ")")
    sub_lab <- paste0(ko_name, " - ", ko_def)

    # making plot
    plot <- ggplot() + geom_bar(data = sub_long, aes(y = norm_cov, x = Depth, fill = Depth), stat = "identity", width = 0.75) +
    theme_bw() + labs(title = main_lab, subtitle = sub_lab, y = "Coverage per Million (CPM)") +
    theme(axis.title.y = element_text(face = "bold"), axis.text.y = element_text(size=12)) +
    theme(axis.title.x = element_text(face = "bold"), axis.text.x = element_text(size=12)) +
    theme(legend.position = "none")

    return(plot)
}

plot_all_N_KOs <- function( tab = N_KO_norm_cov_tab ) {

    # adding unique number of genes to KO ID
    ids <- paste0(tab$KO_ID, " (", tab$num_uniq_genes, ")")
    tab$KO_ID <- ids

    # converting to long format
    long_tab <- tab %>% select(1:5) %>% pivot_longer(-KO_ID, names_to = "Depth", values_to = "norm_cov")

    # plotting

    plot <- ggplot() +
      geom_bar(data=long_tab, aes(y=norm_cov, x=Depth, fill=Depth), stat="identity", width=0.75) +
      facet_wrap(. ~ KO_ID, scales="free_y") + theme_bw() + labs(title = "Nitrogen-related KO coverages", y = "Coverage per Million (CPM)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank()) +
      theme(legend.position="bottom", legend.title=element_text(face="bold")) + labs(fill = "Depth") +
      theme(strip.text = element_text(size = 10, face="bold"))

    return(plot)
}
