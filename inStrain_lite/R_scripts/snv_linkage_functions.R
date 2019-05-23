## Plot and summarize output from strainRep2.py program
## strainRep2.py created by Alex Crits-Christoph: https://github.com/alexcritschristoph/strains_analysis

## R_scripts created by Keith Bouma-Gregson

#### LOG FILES #####
summarize_log_files <- function(path, filename= FALSE){
  require(tidyverse)

  message("Log filenames must have format SampleID.ref_id.pidvalue.log (e.g. Sample1.Genome1.pid98.log) \n
          SampleID= sample reads used in mapping \n
          ref_id= reference genome samples were mapped to \n
          pid= percent identity filter of the mapped reads \n
          Errors will be produced if the filename does not match this format")


  ## List files
  if(filename == FALSE){
    log.files <- list.files(file.path(path), pattern= ".log")
  } else{
    log.files <- filename
  }

  ## Initialize empty data frame
  log_summary <- data.frame(matrix(rep(NA, length(log.files)*13), ncol= 13)) %>%
    as_tibble() %>%
    rename(sample= X1,
           ref_id= X2,
           pid= X3,
           cov_expect_final= X4,
           cov_mean= X5,
           reads_count_initial= X6,
           reads_count_final= X7,
           reads_percent_final= X8,
           total_sites= X9,
           snv_sites= X10,
           snv_bases= X11,
           breadth= X12,
           clonality= X13)

  ## Loop over log files
  counter <-  1
  for(file in 1:length(log.files)){

    ## Read log file
    log.file <- suppressMessages(read_tsv(file.path(path, log.files[file]),
                                          col_names = FALSE,
                                          progress= FALSE)) %>%
      rename(rows= X1)

    ## Extract information from each row
    sample_name <- str_split(log.files[file], "\\.")[[1]][1]
    ref_ID <- str_split(log.files[file], "\\.")[[1]][2]
    PID <- str_split(log.files[file], "\\.")[[1]][3] %>%
      str_replace("pid", "")

    fasta_length <- log.file %>%
      filter(str_detect(rows, "total fasta length")) %>%
      str_replace("^.*: ", "")

    cov_expect_final <- log.file %>%
      filter(str_detect(rows, "\\(final\\) expected")) %>%
      str_replace("^.*: ", "")

    cov_mean <- log.file %>%
      filter(str_detect(rows, "Mean coverage")) %>%
      str_replace("^.*: ", "")

    reads_count_initial <- log.file %>%
      filter(str_detect(rows, "total reads")) %>%
      str_replace("^.*: ", "")

    reads_count_final <- log.file %>%
      filter(str_detect(rows, "\\(final\\) reads")) %>%
      str_replace("^.*: ", "") %>%
      str_replace(" \\(.*$", "")

    reads_percent_final <- log.file %>%
      filter(str_detect(rows, "\\(final\\) reads")) %>%
      str_replace("^.*: ", "") %>%
      str_replace("^.* \\(", "") %>%
      str_replace("%\\)", "")

    total_sites <- log.file %>%
      filter(str_detect(rows, "Total sites")) %>%
      str_replace("^.*: ", "")

    snv_sites <- log.file %>%
      filter(str_detect(rows, "SNVs-sites")) %>%
      str_replace("^.*: ", "")

    snv_bases <- log.file %>%
      filter(str_detect(rows, "SNV-bases")) %>%
      str_replace("^.*: ", "")

    breadth <- as.numeric(total_sites) / as.numeric(fasta_length)

    clonality <- log.file %>%
      filter(str_detect(rows, "clonality")) %>%
      str_replace("^.*: ", "")

    ## Add information to data frame
    log_summary[counter, ] <- c(sample_name, ref_ID, PID, cov_expect_final, cov_mean, reads_count_initial,
                                reads_count_final, reads_percent_final, total_sites,
                                snv_sites, snv_bases, breadth, clonality)

    counter <- counter + 1
  }

  ## Re-format column classes
  log_summary <- log_summary %>%
    mutate_at(vars(cov_expect_final:clonality), funs(round(as.numeric(.), 5)))

  return(log_summary)
}


#### LINKAGE AND FREQ FILES #####
analyse_linkage_data <- function(freq_threshold= 0.05, ref_id, pid, path, output, linkage_plot= TRUE, freq_plot= TRUE, min.breadth= 0.6){
  ## LIBRARIES #################################################################
  require(tidyverse)
  require(progress)
  require(ggplot2)
  require(RColorBrewer)

  ## CHECK FOR ERRORS ############################################################
  if(str_detect(as.character(pid), "\\.")){
    stop("pid must be an integer with no decimals (e.g. 0.99 = 99)")
  }

  if(!is.character(ref_id)){
    stop("ref_id must be a character string identifying reference sequence identity \n
         corresponding to the .log file name")
  }


  ## DEFINE VARIABLES ############################################################
  dir_input <- path
  dir_output <- output
  #species.num <- paste0("species_", species)
  pid_value <- as.character(pid)
  linkage.files <- list.files(file.path(dir_input), pattern= str_c(pid_value, "linkage", sep= "."))
  linkage.sample <- str_replace(linkage.files, "_0\\..*$", "")
  freq.files <- list.files(file.path(dir_input), pattern= str_c(pid_value, "freq", sep= "."))
  freq.sample <- str_replace(freq.files, "_0\\..*$", "")

  ##### PLOTTING PARAMETERS ####################################################
  x_axis_format_distance <- scale_x_continuous(expand= c(0, 5))
  x_axis_format_window <- scale_x_discrete(labels= NULL, expand= c(0.02, 0))
  y_axis_format <- scale_y_continuous(limits= c(0, 1),
                                      breaks= seq(0, 1, by= 0.25),
                                      expand= c(0.02, 0))

  ## ggplot theme for snv linkage
  theme_snv <- theme(panel.grid = element_blank(),
                     plot.margin = unit(c(1, 1, 1, 1), "cm"),
                     text = element_text(size= 14),
                     plot.background = element_rect(fill = "transparent"), # bg of the plot
                     panel.background = element_rect(fill= "transparent", color="black"),
                     axis.text = element_text(colour="black"),
                     axis.title.x = element_text(vjust = -0.75),
                     axis.title.y = element_text(vjust = 1.5),
                     legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                     legend.key = element_blank(),
                     strip.background=element_rect(fill="transparent", color="transparent"),
                     legend.position = "top")
  ################################################################################


  ##### LOOP THROUGH .LINKAGE FILES AND MAKE PLOTS #########################################

  if(linkage_plot == TRUE){

    ## progress bar
    pb.linkage <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
    pb.linkage$tick(0)

    rm(counter)

    ## Initialize data frame to store model coefficients
    mod.coef <- data.frame(sample= as.character(NULL),
                           ref_id= as.character(NULL),
                           pid= as.character(NULL),
                           df= as.character(NULL),
                           intercept= as.numeric(NULL),
                           slope= as.numeric(NULL),
                           r2_adj= as.numeric(NULL),
                           p_value= as.numeric(NULL),
                           stringsAsFactors= FALSE)

    #### THE LOOP ####
    message("Looping over .linkage files")

    for(file in 1:length(linkage.files)){
      pb.linkage$tick() ## Progress Bar

      #### SAMPLE ID INFORMATION FOR OUTPUT FILES
      sample.name <- str_replace(linkage.files[file], "_0\\..*$", "")
      plot.name <- str_c(sample.name, ref_id, str_c("PID=0.", pid_value), sep= "-")

      ## CHECK TO MAKE SURE BREADTH > 0.6
      ## if not, then skip this sample and move to the next iteration of the loop
      log.file <- str_c(sample.name, ".", ref_id, ".pid", pid_value, ".log")
      samp.breadth <- summarize_log_files(path= file.path(path, "log_files"), filename= log.file) %>%
        select(breadth)

      if(samp.breadth < min.breadth){
        message("Breadth <0.6, skipping sample: ", sample.name)
        next()
      }

      #### READ IN THE DATA

      ## Coverage data from log file
      snvs.coverage <- try(suppressMessages(read_tsv(file.path(dir_input, "log_files", str_c(linkage.sample[file], ".", ref_id, ".pid", pid_value, ".log")),
                                                     col_names = FALSE,
                                                     progress= FALSE)) %>%
                             pull() %>%
                             .[str_detect(., "Mean coverage")] %>%
                             str_replace("Mean coverage: ", "") %>%
                             as.numeric(.) %>%
                             round(., 1))

      if(class(snvs.coverage) != "try-error"){
        ## Linkage data
        snvs <- suppressMessages(read_tsv(file.path(dir_input, linkage.files[file]), progress= FALSE)) %>%
          mutate(pid= pid_value) %>%
          filter(total > snvs.coverage*freq_threshold, complete.cases(.)) %>%
          rename(row_num= X1) %>%
          mutate(r2= ifelse(r2 > 1, 1, r2), # some r2 values seem to be just above 1
                 dist.bin= cut(Distance, breaks= 10),
                 scaffold= str_replace(Window, ":.*$", ""))

        #### EXPONENTIAL DECAY MODEL
        ## use try() to make sure lm function did not error, otherwise loop crashes
        exp.model <- try(lm(log(r2) ~ Distance, data= subset(snvs, r2 > 0)))

        ## Make dataframes of predicted data for ggplots below
        if(class(exp.model) != "try-error"){
          exp.model.df <- data.frame(pid= pid_value,
                                     dist.values= snvs$Distance,
                                     pred.values= exp(predict(exp.model, list(Distance= snvs$Distance))))
        }


        ## Export model coefficients to a table
        if(exists("exp.model.df")){

          if(exists("counter") == FALSE){
            counter <- 1
          }

          mod.coef[counter, ] <- c(sample.name, ref_id, pid_value, summary(exp.model)$df[2],
                                   round(exp.model$coefficients[1], 5), round(exp.model$coefficients[2], 5), round(summary(exp.model)$adj.r.squared, 4), round(summary(exp.model)$coef["Distance", "Pr(>|t|)"], 4))
          counter <- counter + 1
        }

        #### MAKE PLOTS ################################################################


        ## DISTANCE x R^2
        ## Only plot regression line if linear model did not error
        if(exists("exp.model.df")){
          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            geom_line(data = exp.model.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            theme_snv
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)

        } else {

          ggplot(data= snvs, aes(x= Distance, y= r2)) +
            geom_hline(yintercept = 0, size= 0.25) +
            geom_point() +
            #geom_line(data = exp.model.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
            labs(x= "Pairwise distance (bp)", y= expression(r^2), title= plot.name) +
            y_axis_format +
            x_axis_format_distance +
            theme_snv
          #ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".pdf"), width= 8, height= 6, units= "in", path= dir_output, device= cairo_pdf)
          ggsave(last_plot(), filename= str_c(plot.name, "_dist", ".jpg"), width= 8, height= 6, units= "in", dpi= 320,  path= dir_output)
        }

        ## Remove model dataframes for next iteration of the loop
        rm(exp.model.df)

      }
    }


    ## Clean up mod.coef table
    message("Writing model coefficients table")
    mod.coef <- mod.coef %>%
      mutate_at(vars(df:p_value), funs(as.numeric)) %>%
      as_tibble()
    write_tsv(mod.coef, path= file.path(dir_output_table, str_c("model_coef_", ref_id, ".tsv")))
  }

  #### LOOP OVER .FREQ FILES ###################################################
  if(freq_plot == TRUE){

    ## progress bar
    pb.freq <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
    pb.freq$tick(0)

    message("Looping over .freq files")
    for(file in 1:length(freq.files)){
      pb.freq$tick()

      #### SAMPLE ID INFORMATION FOR OUTPUT FILES
      sample.name <- str_replace(freq.files[file], "_0\\..*$", "")
      plot.name <- str_c(sample.name, ref_id, str_c("PID=0.", pid_value), sep= "-")


      ## CHECK TO MAKE SURE BREADTH > 0.6
      ## if not, then skip this sample and move to the next iteration of the loop
      log.file <- str_c(sample.name, ".", ref_id, ".pid", pid_value, ".log")
      samp.breadth <- summarize_log_files(path= file.path(path, "log_files"), filename= log.file) %>%
        select(breadth)

      if(samp.breadth < min.breadth){
        message("Breadth <0.6, skipping sample: ", sample.name)
        next()
      }
      #### READ IN THE DATA
      freq.df <- suppressMessages(read_tsv(file.path(dir_input, freq.files[file]), progress= FALSE)) %>%
        rename(row_num= X1)

      ## Count histogram frequencies
      freq.counts <- freq.df %>%
        mutate(freq.rounded= round(freq, 3)) %>%
        count(freq.rounded)

      ggplot(freq.counts, aes(x= freq.rounded, y= n)) +
        geom_vline(xintercept= 0.5, size= 0.25, color= "gray") +
        geom_point(size= 2, alpha= 0.5) +
        geom_smooth(method= "loess", span= 0.1, se= FALSE, size= 1, color= "tomato") +
        labs(x= "SNV frequency", y= "Count",  title= plot.name) +
        scale_x_continuous(breaks= seq(0, 1, 0.1), limits= c(0, 1), expand= c(0.01, 0)) +
        scale_y_continuous(breaks= scales::pretty_breaks()) +
        #facet_wrap(~sample, ncol= 3, scales= "free_y") +
        theme_snv
      ggsave(last_plot(), filename= str_c(plot.name, "_freq", ".jpg"), width= 8, height= 6, units= "in", dpi= 320, path= dir_output)
    }
  }

  return(mod.coef)
}


#### MULTI PANEL FIGURES ######
make_multi_panel_fig <- function(path, file.pattern=NULL, file.list=NULL, output.name){
  require(tidyverse)
  require(progress)
  require(multipanelfigure)

  if(all(!missing(file.pattern), !missing(file.list))){
    stop("Only define either file.pattern or file.list, not both")
  }

  if(!missing(file.pattern)){
    # Get number of figures in folder
    num.figures <- length(list.files(path, pattern= file.pattern))
    last_file <- list.files(path, pattern= file.pattern)[num.figures]

    # Make blank multi-panel figure
    if(exists("multi_fig") == FALSE){
      multi_fig <- multi_panel_figure(width= 300, height= 30*ceiling(num.figures/6), unit= "mm", rows= ceiling(num.figures/6), columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")

    }

    # Loop to fill panels
    message("Filling panels")
    pb.panels <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = num.figures)
    pb.panels$tick(0)

    for(file in list.files(path, pattern= file.pattern)){
      pb.panels$tick() ## Progress Bar

      multi_fig <- suppressMessages(fill_panel(figure= multi_fig, panel= file.path(path, file), scaling= "shrink"))
    }
  }


  if(!missing(file.list)){
    # Get number of figures in folder
    num.figures <- length(file.list)
    last_file <- file.list[num.figures]

    # Make blank multi-panel figure
    if(exists("multi_fig") == FALSE){
      multi_fig <- multi_panel_figure(width= 300, height= 30*ceiling(num.figures/6), unit= "mm", rows= ceiling(num.figures/6), columns= 6, row_spacing= 0, column_spacing= 0, panel_label_type = "none")

    }

    # Loop to fill panels
    message("Filling panels")
    pb.panels <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = num.figures)
    pb.panels$tick(0)

    for(file in file.list){
      pb.panels$tick() ## Progress Bar

      multi_fig <- suppressMessages(fill_panel(figure= multi_fig, panel= file.path(path, file), scaling= "shrink"))
    }
  }


  # Save figure
  message("Saving file")
  save_multi_panel_figure(multi_fig, filename= file.path(path, str_c(output.name, ".jpg")), dpi= 320)
  rm(multi_fig, last_file)
}


#### FILTER R^2 FILE BY SNV FREQ  ####

filter_snv_freq_window <- function(path, output, ref_id= NULL, pid, max.freq= 0, min.freq= 0, min.snv.sites= 100, overwrite= TRUE){
  ## LIBRARIES #################################################################
  require(tidyverse)
  require(progress)
  require(ggplot2)
  require(hexbin)
  require(RColorBrewer)
  require(multipanelfigure)

  ## CHECK FOR ERRORS ############################################################
  if(str_detect(as.character(pid), "\\.")){
    stop("pid must be an integer with no decimals (e.g. 0.99 = 99)")
  }

  if(!any(c(str_detect(min.freq, "\\."), str_detect(max.freq, "\\.")))){
    stop("freq filter values must be < 1, and include a zero and decimal. \n(e.g. 0.2)")
  }

  if(any(min.freq > max.freq)){
    stop("Minimum frequency value must be < than maximum freq value")
  }

  if(any(length(min.freq) != length(max.freq))){
    stop("Same number of minimum and maximum freq values must be supplied")
  }



  ## DEFINE VARIABLES ############################################################

  message(paste("filtering by >", str_c(as.character(min.freq), " and < ", str_c(as.character(max.freq), "\n"))))

  pid_value <- as.character(pid)
  linkage.files <- list.files(file.path(path), pattern= str_c(pid_value, "linkage", sep= "."))
  linkage.sample <- str_replace(linkage.files, "_0\\..*$", "")
  freq.files <- list.files(file.path(path), pattern= str_c(pid_value, "freq", sep= "."))
  freq.sample <- str_replace(freq.files, "_0\\..*$", "")

  ## Initialize data frame to store model coefficients and filtering summary statistics
  summary.table <- data.frame(sample= as.character(NULL),
                              ref_id= as.character(NULL),
                              pid= as.character(NULL),
                              filter_window= as.character(NULL),
                              df= as.character(NULL),
                              intercept= as.numeric(NULL),
                              slope= as.numeric(NULL),
                              r2_adj= as.numeric(NULL),
                              p_value= as.numeric(NULL),
                              prop_r2_equals_1= as.numeric(NULL),
                              stringsAsFactors= FALSE)



  ## ggplot theme for snv linkage
  theme_snv <- theme(panel.grid = element_blank(),
                     plot.margin = unit(c(1, 1, 1, 1), "cm"),
                     text = element_text(size= 14),
                     plot.background = element_rect(fill = "transparent"), # bg of the plot
                     panel.background = element_rect(fill= "transparent", color="black"),
                     axis.text = element_text(colour="black"),
                     axis.title.x = element_text(vjust = -0.75),
                     axis.title.y = element_text(vjust = 1.5),
                     legend.background = element_rect(size=0.25, color="black", fill= "transparent"),
                     legend.key = element_blank(),
                     strip.background=element_rect(fill="transparent", color="transparent"),
                     legend.position = "top")
  hex.color.gradient <- scale_fill_gradient(name = "count", trans = "log10", low= "black", high= "deepskyblue", breaks= c(1, 10, 100, 1000), limits= c(1, 1000))


  ## Progress bar
  message("Looping over files")
  pb.freq.filt <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(linkage.files))
  pb.freq.filt$tick(0)
  #file=3
  ## LOOP OVER LINKAGE FILES ############################################################
  for(file in 1:length(linkage.files)){
    pb.freq.filt$tick() ## Progress Bar

    sample.name <- str_replace(linkage.files[file], "_0\\..*$", "")

    ## CHECK TO MAKE SURE AT LEAST 100 SNV SITES
    ## if not, then skip this sample and move to the next iteration of the loop
    log.file <- str_c(sample.name, ".", ref_id, ".pid", pid_value, ".log")
    snv.sites <- summarize_log_files(path= file.path(path, "log_files"), filename= log.file) %>%
      select(snv_sites)

    if(snv.sites < min.snv.sites){
      message("<100 SNV sites, skipping sample: ", sample.name)
      next()
    }

    freq.file <- suppressMessages(read_tsv(file.path(path, freq.files[file]), progress = FALSE)) %>%
      mutate(pid= pid_value) %>%
      filter(complete.cases(.)) %>%
      rename(row_num= X1)

    link.file <- suppressMessages(read_tsv(file.path(path, linkage.files[file]), progress= FALSE)) %>%
      mutate(pid= pid_value) %>%
      filter(complete.cases(.)) %>%
      rename(row_num= X1) %>%
      mutate(r2= ifelse(r2 > 1, 1, r2)) %>%  # some r2 values seem to be just above 1
      separate(total_A, sep="'", into= c("sep1", "snv_A", "sep2"), remove= FALSE) %>%
      separate(total_a, sep="'", into= c("sep1", "snv_a", "sep2"), remove= FALSE) %>%
      separate(total_B, sep="'", into= c("sep1", "snv_B", "sep2"), remove= FALSE) %>%
      separate(total_b, sep="'", into= c("sep1", "snv_b", "sep2"), remove= FALSE) %>%
      select(-sep1, -sep2)

    link.file.long <- link.file %>%
      gather(key= "allele", value= "SNV", c(snv_A, snv_a, snv_B, snv_b)) %>%
      left_join(subset(freq.file, select= -c(row_num, Window, pid)))

    link.file.wide <- link.file.long %>%
      select(-SNV) %>%
      spread(key= allele, value= freq) %>%
      rename(freq_A= snv_A, freq_a= snv_a, freq_B= snv_B, freq_b= snv_b)

    ## EXPONENTIAL DECAY MODEL
    ## use try() to make sure lm function did not error, otherwise loop crashes
    exp.model.unfiltered <- try(lm(log(r2) ~ Distance, data= subset(link.file, r2 > 0)))

    ## Make dataframes of predicted data for ggplots below
    if(class(exp.model.unfiltered) != "try-error"){
      exp.model.unfiltered.df <- data.frame(pid= pid_value,
                                            dist.values= link.file$Distance,
                                            pred.values= exp(predict(exp.model.unfiltered, list(Distance= link.file$Distance))))
      p.value <- try(round(summary(exp.model.unfiltered)$coefficients["Distance", "Pr(>|t|)"], 5))

      ## CALCULATE PROPORTION OF R^2 = 1
      prop_r2_1_unfiltered <- link.file %>%
        count(r2 == 1) %>%
        mutate(freq= n / sum(n)) %>%
        select(freq) %>%
        slice(2) %>%
        pull()
    }
    if(class(p.value) == "try-error"){
      p.value <- NA
      prop_r2_1_unfiltered <- NA
    }


    ## Initialize list to store plots
    plot.list <- rep(list(NULL), length(min.freq)+1)
    names(plot.list)[1] <- "unfiltered_freq"
    for(id in 2:length(plot.list)){
      names(plot.list)[id] <- str_c("filtered_freq_", str_c(min.freq[id-1], max.freq[id-1], sep= ":"))
    }

    plot.list[[1]] <- ggplot(data= link.file, aes(x= Distance, y= r2)) +
      geom_hline(yintercept = 0, size= 0.25) +
      geom_point(color= "transparent") +
      geom_hex(bins= 50) +
      geom_line(data = exp.model.unfiltered.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
      labs(x= "Pairwise distance (bp)", y= expression(r^2), title= str_c(sample.name, ref_id, str_c("PID=0.", pid_value), str_c("freq=ALL"), sep= "-")) +
    #  y.axis.format +
      hex.color.gradient +
      theme_snv


    if(exists("counter") == FALSE){
      counter <- 1
    }

    ## Add unfiltered SNV data to summary.table
    summary.table[counter, ] <- c(sample.name, ref_id, pid_value, "unfiltered",
                                  summary(exp.model.unfiltered)$df[2],
                                  round(exp.model.unfiltered$coefficients[1], 5),
                                  round(exp.model.unfiltered$coefficients[2], 5),
                                  round(summary(exp.model.unfiltered)$adj.r.squared, 4),
                                  p.value,
                                  round(prop_r2_1_unfiltered, 4))

    counter <- counter + 1

    ## LOOP OVER FREQUENCY FILTER VALUES #######################################
    for(f in 1:length(min.freq)){
      plot.name <- str_c(sample.name, ref_id, str_c("PID=0.", pid_value), str_c("freq=", min.freq[f], "-", max.freq[f]), sep= "-")

      filter.window <- str_c(min.freq[f], max.freq[f], sep= "_")

      link.file.filtered <- link.file.wide %>%
        filter_at(vars(starts_with("freq")), all_vars((. > min.freq[f] & . < max.freq[f]) | (. < 1-min.freq[f] & . > 1-max.freq[f]) ))


      ## EXPONENTIAL DECAY MODEL
      ## use try() to make sure lm function did not error, otherwise loop crashes
      exp.model.filtered <- try(lm(log(r2) ~ Distance, data= subset(link.file.filtered, r2 > 0)))
      #summary(exp.model.filtered)

      ## Make dataframes of predicted data for ggplots below
      if(class(exp.model.filtered) != "try-error"){
        exp.model.filtered.df <- data.frame(pid= pid_value,
                                            filt= filter.window,
                                            dist.values= link.file.filtered$Distance,
                                            pred.values= exp(predict(exp.model.filtered, list(Distance= link.file.filtered$Distance))))
        df.model <- summary(exp.model.filtered)$df[2]
        p.value <- try(round(summary(exp.model.filtered)$coefficients["Distance", "Pr(>|t|)"], 5))
        intercept <- try(round(exp.model.filtered$coefficients[1], 5))
        slope.distance <- try(round(exp.model.filtered$coefficients[2], 5))
        r2_adjusted <- round(summary(exp.model.filtered)$adj.r.squared, 4)

        ## CALCULATE PROPORTION OF R^2 = 1
        prop_r2_1 <- link.file.filtered %>%
          count(r2 == 1) %>%
          mutate(freq= n / sum(n)) %>%
          select(freq) %>%
          slice(2) %>%
          pull()
      } else {
        df.model <- NA
        p.value <- NA
        intercept <- NA
        slope.distance <- NA
        r2_adjusted <- NA
      }

      if(length(prop_r2_1) == 0){
        prop_r2_1 <- NA
      }

      if(class(p.value) == "try-error" | is.nan(p.value)){
        p.value <- NA
        intercept <- NA
        slope.distance <- NA
      }

      ## Export model coefficients to a table
      #if(exists("exp.model.filtered.df")){
      summary.table[counter, ] <- c(sample.name, ref_id, pid_value, filter.window,
                                    df.model,
                                    intercept,
                                    slope.distance,
                                    r2_adjusted,
                                    p.value,
                                    round(prop_r2_1, 4))
      counter <- counter + 1
      #}


      ## PLOT THE FILTERED DATA
      if(exists("exp.model.filtered.df")){
        plot.list[[f+1]] <- ggplot(data= link.file.filtered, aes(x= Distance, y= r2)) +
          geom_hline(yintercept = 0, size= 0.25) +
          geom_point(color= "transparent") +
          geom_hex(bins= 50) +
          geom_line(data = exp.model.filtered.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
          labs(x= "Pairwise distance (bp)", y= expression(r^2), title= str_c(sample.name, ref_id, str_c("PID=0.", pid_value), str_c("freq=",filter.window), sep= "-")) +
          hex.color.gradient +
          theme_snv
      } else {
        plot.list[[f+1]] <- ggplot(data= link.file.filtered, aes(x= Distance, y= r2)) +
          geom_hline(yintercept = 0, size= 0.25) +
          geom_point(color= "transparent") +
          geom_hex(bins= 50) +
          #geom_line(data = exp.model.filtered.df, aes(x= dist.values, y = pred.values), color= "tomato", size= 1) +
          geom_label(aes(x= 10, y= 0.5, label= "No SNVs meet filtering threshold"), size= 8) +
          labs(x= "Pairwise distance (bp)", y= expression(r^2), title= str_c(sample.name, ref_id, str_c("PID=0.", pid_value), str_c("freq=",filter.window), sep= "-")) +
          hex.color.gradient +
          theme_snv
      }
      ## Remove model results for the next round of the loop
      rm(exp.model.filtered.df)
    }


    # Make blank multi-panel figure
    num.figs <- length(plot.list)
    multi_fig <- multi_panel_figure(width= 200*num.figs, units= "mm", rows= 1, columns= num.figs, row_spacing= 0, column_spacing= 0, panel_label_type = "none")
    # Loop to fill panels
    for(fig in 1:num.figs){
      multi_fig <- suppressMessages(fill_panel(figure= multi_fig, panel= plot.list[[fig]], panel_clip= "on"))
    }

    # Save figure
    multi_plot_name <- str_c(sample.name,"-" , ref_id, "-pid", pid_value, "-filtered_window", ".jpg")

    if(overwrite == FALSE){
      if(length(list.files(file.path(output), pattern= multi_plot_name)) > 0){
        num.files <- length(list.files(file.path(output), pattern= multi_plot_name))
        multi_plot_name <- str_c(sample.name,"-" ,ref_id, "-pid", pid_value, "-filtered_window-", num.files+1, ".jpg")
      }
    }

    save_multi_panel_figure(multi_fig, filename= file.path(output, multi_plot_name), dpi= 320, limitsize= FALSE)
    rm(multi_fig)
    str_c(sample.name,"-" ,ref_id, "-pid", pid_value, "-filtered_window-", ".jpg")

  }

  ## Clean up summary table and write to TSV file
  message("Writing summary table to .tsv")
  summary.table <- summary.table %>%
    mutate_at(vars(intercept:prop_r2_equals_1), funs(as.numeric)) %>%
    as_tibble() %>%
    slice(-1)

  if(overwrite == TRUE){
    write_tsv(summary.table, path= file.path(output, str_c("filter_summary-pid", pid_value, "-", ref_id, ".tsv")))
  } else {
    num.tables <- length(list.files(file.path(output), pattern= str_c("filter_summary-pid", pid_value, "-", ref_id)))
    write_tsv(summary.table, path= file.path(output, str_c("filter_summary-pid", pid_value, "-", ref_id, "-", num.tables+1, ".tsv")))
  }

  return(summary.table)

}



