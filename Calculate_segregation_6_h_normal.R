  #library(plyr)
  library(tidyverse)
  library(dplyr) ##for merge
  library(seg)
  library(stringr) # for adding leading zeroes
  
  library(Rcpp) ##For the sf package
  library(sf)
  
  library(ggplot2)
  
  #  Load agents
  setwd("E:/PhD/Ablauf ABM/Calculate Segregation/Different_weights_1")
  
  
  ### map : https://www.youtube.com/watch?v=GMi1ThlGFMo
  shp <- st_read("E:/PhD/Ablauf ABM/Calculate Segregation/LOR_SHP_EPSG_25833/Lor_data.shp", stringsAsFactors = FALSE)
  
 # ring <- st_read("E:/PhD/Ablauf ABM/Calculate Segregation/Ring-Bahn/ring-bahn.shp", stringsAsFactors = FALSE)
  
  #mauer <- st_read("E:/PhD/Ablauf ABM/Calculate Segregation/Berliner_Mauer-shp/Berliner_Mauer_Hinterlandmauer.shp", stringsAsFactors = FALSE)
  #setwd("E:/PhD/Ablauf ABM/Calculate Segregation")
  
#   ggplot(new.shp) + 
#     geom_sf(aes(fill = h_i))+
#     scale_fill_gradient(limits = c(quant_range[2], quant_range[4]),low = "#edf8b1", high = "#2c7fb8", oob = scales::squish) + ## squish is used to include the values that are outsinde of limits
#     geom_sf(data = ring, color = "black", size = 1) +
#     geom_sf(data = mauer, color = "black", size = 1) +
#     labs(fill = "H_i")+
# #    scale_fill_manual(values = "black", name = "Tower Location") +    
# #    guides(fill=guide_legend(title=NULL))+    
#     coord_sf()

  
  files <- Sys.glob("E:/PhD/Ablauf ABM/Calculate Segregation/input_files/Different_weights_1/*.csv")
  print(files)
  
  
  

#  filepath <- "E:/PhD/Ablauf ABM/Calculate Segregation/input_files/Different_weights/nrAgents_5_radius_5_tick_15_weight_nw1_weight_age0.75_new_true.csv"  
#  agents <- read.csv(filepath, header = TRUE, sep = ",", skip = 12)  

  res.table <- rbind("name", "entire entropy", "min H-Index", "mean H-Index", "median H-Index", "max H-Index", "min local entropy", "mean local entropy", "median local entropy", "max local entropy")
    
  i <- 1
  
   while (i <= length(files)) {
   
     filepath <- files[i]
  #  filepath <- "E:/PhD/Ablauf ABM/Calculate Segregation/input_files/Different_weights/nrAgents_5_radius_5_tick_15_weight_nw0.75_weight_age0.5_new_true.csv"
    agents <- read.csv(filepath, header = TRUE, sep = ",", skip = 12)
    data <- agents[complete.cases(agents), ]
    
    use.agents <- subset(data, select = c(age, gender, mgbg, mar_stat, fake.breed, lor)) 
    selection <- subset(data, select = c(fake.breed, lor))
    
    count_selection <- count(selection, fake.breed, lor)
    
    data_long <- count_selection %>% group_by(lor, fake.breed) %>% summarise(n = sum(n))
    data_wide <- spread(data_long, key = "fake.breed", value = "n", fill = 0)
    
    data_wide$lor<- str_pad(data_wide$lor, 8, pad = "0")
    
    ######################following: https://rpubs.com/corey_sparks/254870################################
    
    #We need the lor-level totals for the total population and each AT
    data_wide$sum.lor <- rowSums(data_wide[,2:12])
    
    ## sum of each AT
    sumATs <- t(colSums(data_wide[,2:13]))
    sumAll <- cbind("0", sumATs)
    colnames(sumAll) <- colnames(data_wide) ###weil die Namen nicht richtig übereinstimmen
    
    
    data_wide <- rbind(as.data.frame(data_wide), sumAll)
    
    
    #####unlist everything and transform it as numeric
    data_wide[1] <- as.numeric(data_wide[[1]])
    data_wide[2] <- as.numeric(data_wide[[2]])
    data_wide[3] <- as.numeric(data_wide[[3]])
    data_wide[4] <- as.numeric(data_wide[[4]])
    data_wide[5] <- as.numeric(data_wide[[5]])
    data_wide[6] <- as.numeric(data_wide[[6]])
    data_wide[7] <- as.numeric(data_wide[[7]])
    data_wide[8] <- as.numeric(data_wide[[8]])
    data_wide[9] <- as.numeric(data_wide[[9]])
    data_wide[10] <- as.numeric(data_wide[[10]])
    data_wide[11] <- as.numeric(data_wide[[11]])
    data_wide[12] <- as.numeric(data_wide[[12]])
    data_wide[13] <- as.numeric(data_wide[[13]])
    
    #For the multi-group measure of segregation, we also need population proportions
    prop_AT1 <- data_wide[nrow(data_wide),2]/data_wide[nrow(data_wide),13]
    prop_AT10 <- data_wide[nrow(data_wide),3]/data_wide[nrow(data_wide),13]
    prop_AT11 <- data_wide[nrow(data_wide),4]/data_wide[nrow(data_wide),13]
    prop_AT2 <- data_wide[nrow(data_wide),5]/data_wide[nrow(data_wide),13]
    prop_AT3 <- data_wide[nrow(data_wide),6]/data_wide[nrow(data_wide),13]
    prop_AT4 <- data_wide[nrow(data_wide),7]/data_wide[nrow(data_wide),13]
    prop_AT5 <- data_wide[nrow(data_wide),8]/data_wide[nrow(data_wide),13]
    prop_AT6 <- data_wide[nrow(data_wide),9]/data_wide[nrow(data_wide),13]
    prop_AT7 <- data_wide[nrow(data_wide),10]/data_wide[nrow(data_wide),13]
    prop_AT8 <- data_wide[nrow(data_wide),11]/data_wide[nrow(data_wide),13]
    prop_AT9 <- data_wide[nrow(data_wide),12]/data_wide[nrow(data_wide),13]
  
    
    ############# Multi- group segregation ###############################################
    name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filepath))  

    ## E = sum[^M _m=1] pi_m log 1/pi_m
    c_entropy <- prop_AT1*(log(1/prop_AT1)) + prop_AT2*(log(1/prop_AT2)) + prop_AT3*(log(1/prop_AT3)) + prop_AT4*(log(1/prop_AT4)) + 
      prop_AT5*(log(1/prop_AT5)) + prop_AT6*(log(1/prop_AT6)) + prop_AT7*(log(1/prop_AT7)) + prop_AT8*(log(1/prop_AT8)) + 
      prop_AT9*(log(1/prop_AT9)) + prop_AT10*(log(1/prop_AT10)) + prop_AT11*(log(1/prop_AT11))
      
    #Just like for the counties, we also need population proportions at the tract level (local proportion)
    h.AT1 <- data_wide$`"ats1"`/data_wide$sum.lor
    h.AT2 <- data_wide$`"ats2"`/data_wide$sum.lor
    h.AT3 <- data_wide$`"ats3"`/data_wide$sum.lor
    h.AT4 <- data_wide$`"ats4"`/data_wide$sum.lor
    h.AT5 <- data_wide$`"ats5"`/data_wide$sum.lor
    h.AT6 <- data_wide$`"ats6"`/data_wide$sum.lor
    h.AT7 <- data_wide$`"ats7"`/data_wide$sum.lor
    h.AT8 <- data_wide$`"ats8"`/data_wide$sum.lor
    h.AT9 <- data_wide$`"ats9"`/data_wide$sum.lor
    h.AT10 <- data_wide$`"ats10"`/data_wide$sum.lor
    h.AT11 <- data_wide$`"ats11"`/data_wide$sum.lor
    
    #This is the tract-level entropy measure
    local.ent <- h.AT1*(log(1/h.AT1)) + h.AT2*(log(1/h.AT2)) + h.AT3*(log(1/h.AT3)) + h.AT4*(log(1/h.AT4)) + h.AT5*(log(1/h.AT5)) +  h.AT6*(log(1/h.AT6)) +
      h.AT7*(log(1/h.AT7)) + h.AT8*(log(1/h.AT8)) + h.AT9*(log(1/h.AT9)) + h.AT10*(log(1/h.AT10)) + h.AT11*(log(1/h.AT11)) 
    
    local.ent <- ifelse(is.na(local.ent)==T, 0, local.ent)
    local.ent <- as.data.frame(local.ent)
  
    #Now I calculate each tract's contribution to the H index
    H_i <- 1 - local.ent/c_entropy
    
    h.calc <- (data_wide$sum.lor * (c_entropy - local.ent)) / (data_wide[nrow(data_wide),ncol(data_wide)] * c_entropy)
    
    h.calc <- unlist(h.calc)
    
    h.table <- cbind(data_wide[,1], h.calc, local.ent)
    colnames(h.table) <- c("Missing_0", "h_i", "local_ent")
    
    new.shp <- inner_join(shp, h.table, copy = TRUE) 
    
    new.shp$h_i[is.na(new.shp$h_i)] <- 0 
  
    quant_range <- quantile(new.shp$h_i, na.rm = TRUE)

    b <- c(quant_range[2], quant_range[3], quant_range[4])
    
   # c_entropy
  
    # ggplot(new.shp) +
    #   geom_sf(aes(fill = h_i))+
    #   scale_fill_gradient(limits = c(quant_range[2], quant_range[4]), breaks = b, low = "#edf8b1", high = "#2c7fb8", oob = scales::squish) +
    #   labs(title="H Index", x = "", y = "")  +  
    #   geom_text(x=3, y=3, label="Scatter plot")

    setwd("E:/PhD/Ablauf ABM/Calculate Segregation/output_r/h_index")
    
    graph_title_entropy <- paste(name, "H_index", c_entropy, sep = "_")
    png(paste(graph_title_entropy, ".png", sep = ""))
    myplot <- ggplot(new.shp) +
      geom_sf(aes(fill = h_i))+
      scale_fill_gradient(limits = c(quant_range[2], quant_range[4]),low = "#edf8b1", high = "#2c7fb8", oob = scales::squish) + ## squish is used to include the values that are outsinde of limits
      theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
            axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
            )+
#      geom_sf(data = mauer, color = "black", size = 1) + ####New
      labs(title="H Index", x = "", y = "")
    print(myplot)
    dev.off()
    
    setwd("E:/PhD/Ablauf ABM/Calculate Segregation/output_r/local_ent") 
    
    graph_title_entropy <- paste(name, "local entropy", c_entropy, sep = "_")
    png(paste(graph_title_entropy, ".png", sep = ""))
    myplot <- ggplot(new.shp) +
      geom_sf(aes(fill = local_ent))+
      scale_fill_gradient(low = "#edf8b1", high = "#2c7fb8") + ## squish is used to include the values that are outsinde of limits
      theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
            axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
      )+      
#      geom_sf(data = mauer, color = "black", size = 1) + ####New
      labs(title="Local Entropy", x = "", y = "")
    print(myplot)
    dev.off()    
    
    
    col1 <- rbind(name, c_entropy, min(h.table$h_i), mean(h.table$h_i), median(h.table$h_i), max(h.table$h_i), min(h.table$local_ent), mean(h.table$local_ent), median(h.table$local_ent), max(h.table$local_ent))
    res.table <- cbind(res.table, col1)

    i <- i + 1
    
    i
    
  }
  
  write.table(res.table, file = "E:/PhD/Ablauf ABM/Calculate Segregation/output_r/hTable_entropy.csv", sep = ",", row.names = FALSE, col.names = TRUE)
    