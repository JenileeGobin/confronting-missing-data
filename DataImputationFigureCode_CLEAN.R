
# Confonting missing data in the age of pandemic lockdown - Hossie et al. 
# Code for generating Figure 1
# Author: Jenilee Gobin
# Date: June 30, 2020

#==========================================================================

# LOAD LIBRARIES -----------

  library (ggplot2) # package for plotting
  library (tidyverse) # package to create tidy code
  library (zoo) # package for mean imputation
  library (fANCOVA) #package for loess regression and smoothing spline
  library (mice)
  library (broom)
  library (lavaan)
  library (viridis)
  library (kableExtra)
  library (cowplot)
  
  set.seed (87) # to generate reproducible results


# CREATE DATA - create time series data with a linear trend and a gap ----------

  dat <- data.frame(year = 2006:2021) 
  
  trend_dat <- dat %>% 
    mutate(randnum = sample(-30:30, length(year)), x = year - 2006, complete = 3 * x + 10 + randnum) # specify normally distributed random number to add stochasticity, an x value based on the year, and calculate a y value that is linearly related to x while incorporating stochasticity  
  
  #write.csv(trend_dat, "trend_dat.csv") # commented out because I forgot to set seed beforehand and therefore wrote data to a csv file!
  
  
  # READ DATA - read in the data I used (because I forgot to set seed when initially generating the data) 
  
  input_file<-read.csv("trend_dat.csv")
  trend_dat<-input_file
  
  
  # CREATE GAP - create a series with a gap
  
  trend_dat$gap<-trend_dat$complete
  trend_dat$gap[length(trend_dat$gap)-1]<-NA
  

# FILL GAP/FIT LINEAR MODEL (COMMON) - fill gap and fit lm using common approaches ----------

  # MEAN IMPUTATION - gaps filled using mean imputation 
  
  mean_imput<-zoo(trend_dat$gap, order.by = trend_dat$year) # create zoo object with data containing gap
  mean_imput<-as.data.frame(na.approx(mean_imput)) # fill in NA's with mean of values above and below and output as data.frame to be able to add to trend_dat later
  names(mean_imput)<-c("mean_imput") # renaming column for incorporation with trend_dat 

  # LOWESS REGRESSION AND SPLINE SMOOTHING - gaps filled using lowess regression and smoothing spline
  
  trend_dat_nas_omitted<-na.omit(trend_dat) # omit NAs because fANCOVA can't deal with missing values and does not accept na.rm argument
  
  # LOWESS (GCV) - loess regression using generalized cross validation criterion for smoothing parameter selection
  
  loess_gcv_model <- loess.as(x = trend_dat_nas_omitted$year, y = trend_dat_nas_omitted$gap, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
  
  # LOWESS (AICc) - loess regression using AICc criterion for smoothing parameter selection
  
  loess_aicc_model <- loess.as(x = trend_dat_nas_omitted$year, y = trend_dat_nas_omitted$gap, degree = 1, criterion = c("aicc", "gcv")[1], user.span = NULL, plot = F)
  
  # SMOOTHING SPLINE - cubic smoothing spline 
  
  spline_model<-smooth.spline(x = trend_dat_nas_omitted$year, y = trend_dat_nas_omitted$gap)

  # FIT LINEAR MODELS (COMMON) - fit linear models for common approaches above
  
    # add common approaches output to trend_dat dataframe 
    
    trend_dat <- trend_dat %>% 
      mutate(mean_imput = mean_imput$mean_imput, loess_gcv = gap, loess_aicc = gap, spline = gap) %>% 
      replace_na(list(loess_gcv = predict(loess_gcv_model, 2020), 
                      loess_aicc = predict(loess_aicc_model, 2020), 
                      spline = predict(spline_model, x=2020)$y))
    
    # reshape to long format
    
    trend_dat_long<-trend_dat %>% 
      gather(approach, value, c("complete", "gap", "mean_imput", "loess_gcv", "loess_aicc", "spline")) 
    
    # fit linear models for each approach
    
    common_fits<-trend_dat_long %>% 
      group_by(approach) %>%
      do(common_lm = lm(value ~ year, data=.))
    
    # create dataframe with linear model output
    
    common_stats<-common_fits %>% 
      tidy(common_lm) %>% 
      merge(glance(common_fits,common_lm), by="approach")
    

# FILL GAPS/ESTIMATE PARAMETERS (ROBUST) - fill gaps/estimate parameters using robust approaches -----

  # MULTIPLE IMPUTATION (MI) - fill gaps using MI
  
    # create MI dataframe
    
    MI_dat<-trend_dat %>% 
      dplyr::select (c(year, gap))
    
    # visualize missing data
    
    md.pattern(MI_dat)
    
    # impute with predictive mean matching (ppm method)
    
    imputed_dat <- mice(MI_dat, m=5, maxit = 50, method = 'pmm', seed = 500)
    
    # check imputed values
    
    str(imputed_dat)
    imputed_dat$imp$gap
    
    # build predictive model
    
    MI_fit <- with(data = imputed_dat, exp = lm(gap ~ year))
    summary(MI_fit)
    
    # pool fits to combine estimates and output  stats
    
    combine_fit <- pool(MI_fit) # combine fit estimates by Rubin's rules
    str(combine_fit)
    
    # extract stats output
    
    MI_output <- summary(combine_fit) 
    
    (MI_rsquared<-pool.r.squared(MI_fit, adjusted = T)[1]) 
    
    # compile MI stats output in a dataframe
    
    MI_output_complete <- MI_output %>% 
      mutate(approach = "MI", adj.r.squared = MI_rsquared)

  # ML-BASED APPROACH (FIML) - fill gaps using full information maximum likelihood (FIML)
  
    # create FIML dataframe
    
    FIML_dat <- trend_dat %>% 
      dplyr::select (c(year, gap))
    
    # fit structural equation model with fiml for missing values 
    
    FIML_fit <- sem("gap~year", FIML_dat, missing="fiml") 
    
    # visualize pattern and coverage of missing data
    
    inspect(FIML_fit, 'patterns') 
    inspect(FIML_fit, 'coverage')
    
    # extraxt stats output 
    
    FIML_output <- summary(FIML_fit, fit.measures=TRUE, rsquare=TRUE)
    
    # compile FIML stats output in a dataframe
    
    FIML_output_complete <- data.frame(approach = rep("FIML",2),
                                     term = c(FIML_output$PE$rhs[1], "(Intercept)"),
                                     estimate = c(FIML_output$PE$est[1], FIML_output$PE$est[4]), 
                                     std.error = c(FIML_output$PE$se[1], FIML_output$PE$se[1]), 
                                     adj.r.squared = rep(FIML_output$PE$est[6],2))
    

# SENSITIVITY ANALYSIS - conduct sensitivity analysis for high and low scenarios if data are MNAR -----
    
    # create sensitivity analysis dataframe
    
    sensitivity_dat<-trend_dat %>% 
      dplyr::select (c(year, gap))
    
    # estimate mean and sd to determine high and low values
    
    (gap_mean<-mean(sensitivity_dat$gap, na.rm=T))
    (gap_sd<-sd(sensitivity_dat$gap, na.rm=T))
    
    sensitivity_dat<-sensitivity_dat %>% 
      mutate(high = gap, low = gap) %>% 
      replace_na(list(high = gap_mean+3*gap_sd, 
                 low = gap_mean-3*gap_sd))
    
    # FIT LINEAR MODELS (SENSITIVITY) - fit linear models for sensitivity analyses
    
    # reshape to long format
    
    sensitivity_dat_long<-sensitivity_dat %>% 
      gather(approach, value, c("gap", "high", "low"))
    
    # fit linear models for each scenario
    
    sensitivity_fits<-sensitivity_dat_long %>% 
      filter(approach != "gap") %>% 
      group_by(approach) %>%
      do(sensitivity_lm = lm(value ~ year, data=.))
    
    # create dataframe with linear model output
    
    sensitivity_stats<-sensitivity_fits %>% 
      tidy(sensitivity_lm) %>% 
      merge(glance(sensitivity_fits,sensitivity_lm), by="approach")
    

# CREATE TABLES - generate tables for common approaches, robust approaches, sensitivity analysis and complete dataset ----------
    
    # TABLE FUNCTIONS ----------
  
      # function to prepare data frames for tables
      
      make_table <- function(df, row_names){
        
        t<-df %>%
          filter(term == "year") %>%
          select(c(approach, symbol, term, estimate, std.error, adj.r.squared)) %>% 
          mutate(approach = factor(approach, levels = row_names)) %>% 
          arrange(approach) %>% 
          rename(slope = estimate, SE = std.error, R2 = adj.r.squared) %>% 
          mutate_if(is.numeric, round, digits=2) %>% 
          `row.names<-`(row_names) %>% # need to pass and re-assign row names because not supported by dplyr
          select(-c(approach,term))
          
        return(t)
        
      }
    
      # function to generate tables with kableExtra
      
      make_kable <- function (t, pal, row_names, header){
        
        k <- t %>% 
          mutate(symbol = cell_spec(symbol, color = pal, bold=T, font_size = "large", escape = F)) %>%
          `row.names<-`(row_names) %>% # need to pass and re-assign row names because not supported by dplyr
          kable(escape = F, align = c("c", rep("l", 3)), row.names = T) %>%
          kable_styling("basic",font_size = 14, full_width = F) %>% 
          add_header_above(header, bold=T, font_size = 18, align = "l")
        
        return(k)
        
      }

  # generate colour palette for symbols
  
  v_palette<-viridis(8, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
  palette<-c("#000000FF", v_palette, "grey85")

  # COMMON APPROACHES TABLE ----------
  
    # prepare data for common approaches table
      
    common_stats_clean <- common_stats %>% 
      filter(approach != "complete") %>% 
      mutate(symbol = approach, 
             symbol = replace(symbol, symbol == "gap", "&#9679;"), 
             symbol = replace(symbol, symbol == "mean_imput", "&#9670;"), 
             symbol = replace(symbol, symbol == "loess_aicc", "&#9660;"), 
             symbol = replace(symbol, symbol == "loess_gcv", "&#9632;"), 
             symbol = replace(symbol, symbol == "spline", "&#9650;"),
             approach = replace(approach, approach == "gap", "case-wise deletion"), 
             approach = replace(approach, approach == "mean_imput", "mean imputation"), 
             approach = replace(approach, approach == "loess_aicc", "loess regression (AICc)"), 
             approach = replace(approach, approach == "loess_gcv", "loess regression (GCV)"), 
             approach = replace(approach, approach == "spline", "smoothing spline"))  
    
    common_row_names <- c("case-wise deletion", "mean imputation", "loess regression (GCV)", "loess regression (AICc)", "smoothing spline") # the names and order of approaches to be listed in the table
    common_palette <- c(palette[c(1,5,3,9,7)])
    common_header <- c("a. Common approaches" = 5)
    
    # make common approaches table
    
    common_table <- make_table(common_stats_clean, common_row_names)
    (common_kable<-make_kable(common_table, common_palette, common_row_names, common_header))
    
  # ROBUST APPROACHES TABLE ----------

    # prepare data for robust approaches table
    
    robust_stats_clean <- MI_output_complete %>% 
      relocate(approach) %>% 
      select(-c(statistic, df, p.value)) %>% 
      rbind(FIML_output_complete) %>% 
      mutate(symbol = approach, 
             symbol = replace(symbol, symbol == "MI", "&#9675;"),
             symbol = replace(symbol, symbol == "FIML", "&#9549;"))
    
    robust_row_names <- c("MI", "FIML") # the names and order of approaches to be listed in the table
    robust_palette <- c(palette[c(2,6)])
    robust_header <- c("b. Robust approaches" = 5)
    robust_table_filename <- "Table_RobustApproaches.pdf"
    
    # make robust approaches table
    
    robust_table <- make_table(robust_stats_clean, robust_row_names)
    (robust_kable<-make_kable(robust_table, robust_palette, robust_row_names, robust_header))
    
  # SENSITIVITY ANALYSIS TABLE ----------

    # prepare data for sensitivity analyses table
    
    sensitivity_stats_clean <- sensitivity_stats %>% 
      mutate(symbol = approach, 
             symbol = replace(symbol, symbol == "high", "&#9651;"),
             symbol = replace(symbol, symbol == "low", "&#9661;"))
    
    sensitivity_row_names <-c ("high", "low") # the names and order of approaches to be listed in the table
    sensitivity_palette <- c(palette[c(4,8)])
    sensitivity_header <- c("c. Sensitivity analysis" = 5)
    sensitivity_table_filename <- "Table_SensitivityAnalysis.pdf"
      
    # make sensitivity analysis table
    
    sensitivity_table <- make_table(sensitivity_stats_clean, sensitivity_row_names)
    (sensitivity_kable<-make_kable(sensitivity_table, sensitivity_palette, sensitivity_row_names, sensitivity_header))

  # COMPLETE DATASET TABLE ----------
    
    # make complete dataset table
    
    complete_stats_clean <- common_stats %>% 
      filter(approach == "complete") %>% 
      mutate(symbol = NA)
    
    complete_row_names <- "complete" # the names and order of approaches to be listed in the table
    
    complete_table <- make_table(complete_stats_clean, complete_row_names)
    
    (complete_kable <- complete_table %>% 
      select(-symbol) %>% 
      kable(escape = F, align = c(rep("l", 3)), row.names = F) %>%
      kable_styling("basic",font_size = 14, full_width = F) %>% 
      add_header_above(c("Complete dataset" = 3), bold=T, font_size = 18, align = "l"))
    

# CREATE PLOTS- generate plots for common approaches, robust approaches and sensitivity analysis ----------  

    # PLOT FUNCTIONS ----------
    
      # function to prepare dataframe for plot
      
      make_stat_df <- function (df){
        
        d <- df %>% 
          select(c(approach, term, estimate)) %>% 
          spread(term, estimate) %>% 
          rename(int='(Intercept)', slope = year)
        
        return (d)
        
      }
    
    # function to make plot
      
      make_plot<-function(df, lev, pal, fl, shp){
        
        df$approach<-factor(df$approach, levels = lev)
        
        (p <- ggplot(df, aes(x = year, y = value, col = approach))
         + geom_rect(xmin = 2019.5, xmax = 2020.5, ymin = -70, ymax = 130, fill = "grey85", col = "grey85") 
         + geom_abline(aes(slope = slope, intercept = int, col = approach), lty = 2, lwd=1.5)
         + geom_path(data=subset(df, approach == "gap"), lty=1, lwd=1)
         + geom_point(data=subset(df, year == 2020), cex = 2, stroke = 1, aes(pch = approach, fill = approach))
         + geom_point(data=subset(df, approach == "gap"), cex = 3, aes(pch = approach, fill = approach))
         + theme_classic()
         + theme(axis.title = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size=10, face="bold", colour = "black"))
         + theme(legend.position = "none")
         + scale_color_manual(values = pal)
         + scale_fill_manual(values = fl)
         + scale_shape_manual(values = shp)
         + scale_y_continuous(limits=c(-60,120))
         + scale_x_continuous(limits=c(2005,2022))
        )
        
        return(p)
      }

  # COMMON PLOT - generate plot for common approaches ----------
      
    # prepare common approach data 
    
    common_stats_for_plot <- make_stat_df(common_stats)
    
    common_plot_dat<-trend_dat_long %>% 
      filter(approach != "complete") %>% 
      select(-x) %>% 
      mutate(approach = factor(approach, levels = c("gap", "loess_gcv", "mean_imput", "spline", "loess_aicc"))) %>% # adjusted levels to points/lines plotted in a particular order
      left_join(common_stats_for_plot)
    
    # make common approach plot
    
    common_palette <- c(palette[c(1,3,5,7,9)])
    common_fill <- c(palette[c(1,3,5,7,9)])
    common_shape <- c(21:25)
    common_levels <- c("gap", "loess_gcv", "mean_imput" , "spline", "loess_aicc")
    
    (common_plot<-make_plot(common_plot_dat, common_levels, common_palette, common_fill, common_shape))

  # ROBUST APPROACHES PLOT - generate plot for robust approaches ----------
    
    # prepare robust data

    robust_stats_for_plot <- make_stat_df(robust_stats_clean)
    
    MI_plot_dat <- data.frame(year = rep(2020, 5),
                              approach = rep("MI", 5), 
                              value = gather(imputed_dat$imp$gap)$value)
    
    FIML_plot_dat <- data.frame(year = 2020, 
                                approach = "FIML", 
                                value = NA)
    
    robust_plot_dat <- trend_dat %>%
      select(c(year, gap)) %>% 
      gather(approach, value, c("gap")) %>% 
      rbind(MI_plot_dat, FIML_plot_dat) %>% 
      mutate(approach = factor(approach, levels = c("gap", "MI", "FIML"))) %>% 
      left_join(robust_stats_for_plot)
    
    # make robust approaches plot ----------

    robust_palette <- c(palette[c(1,6,2)])
    robust_fill <- c(palette[c(1,10,10)])
    robust_shape <- c(21, 0, 1)
    robust_levels <- c("gap", "FIML", "MI")
    
    (robust_plot<-make_plot(robust_plot_dat, robust_levels, robust_palette, robust_fill, robust_shape))

  # SENSITIVITY ANALYSIS PLOT - generate plot for sensitivity analysis ---------- 
    
    # prepare sensitivity analysis data 
    
    sensitivity_stats_for_plot <- make_stat_df(sensitivity_stats)
    
    sensitivity_plot_dat<-sensitivity_dat_long %>% 
      mutate(approach = factor(approach, levels = c("gap", "low", "high"))) %>% 
      left_join(sensitivity_stats_for_plot)
    
    # make sensitivity analysis plot
    
    sensitivity_palette <- c(palette[c(1,4,8)])
    sensitivity_fill <- c(palette[c(1,10,10)])
    sensitivity_shape <- c(21, 24, 25)
    sensitivity_levels <- c("gap", "high", "low")
    
    (sensitivity_plot<-make_plot(sensitivity_plot_dat, sensitivity_levels, sensitivity_palette, sensitivity_fill, sensitivity_shape))

# CREATE MULTI-PANEL FIGURE ----------

  y_axis_label<-ggdraw()+ draw_text("Response variable", fontface="bold", angle=90, size=14)
  
  (figure1_final<-plot_grid(y_axis_label,plot_grid(common_plot, robust_plot, sensitivity_plot, labels=c("a", "b", "c"), hjust = -2, vjust = 2, nrow=1, ncol=3), nrow=1, ncol=2, rel_widths = c(0.2,3)))


  



























