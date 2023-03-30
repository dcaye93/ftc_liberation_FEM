#-------------------------------------------------------
# SIMULATE CARBON GAINS FROM FOCAL TREE LIBERATION FROM LIANAS
# 
#  This script models theoretical growth increases resulting from 
#  liberation from lianas in five mature trees per hectare. 
#
# Code: D.T. Cayetano, E.P. Belair   Last Updated: 3/23/2023
#--------------------------------------------------------

####---- References ----
# Chave 2014: https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12629
# Finlayson et al 2022: https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.8758
# Mills et al 2019: https://doi.org/10.1016/j.foreco.2019.02.023
# Mokany etal 2006: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2486.2005.001043.x
# Schnitzer 2006: https://onlinelibrary.wiley.com/doi/full/10.1111/j.1744-7429.2006.00187.x

rm(list = ls()) # clear environment

####---- Load packages ----####
library(tidyverse)
library(ggpattern)
library(patchwork)
library(ggforce)


####---- Define initial simulation parameters ----####
tree_ht <- 25.0 # tree height in meters
tree_dbh <- 40.0 # starting DBH in cm
tree_dens <- 0.5 # tree density g/cm3
tree_gr <- 0.4 # control growth rate cm per year
liana_gr <- 0.14 # liana growth rate cm per year
root_factor <- 0.235 # root biomass as  a factor of AGB from Mokany etal 2006
n_trees_ha <- 5 # number of trees to be liberated per ha
c_cont <- 0.47 # proportion carbon content of trees and lianas
yr <- seq(from = 0, to = 30, by = 1) # simulation period years


####---- Define function to calculate tree & liana biomass in Mg ----####

biomass_Mg <- function(form, dbh) {
  
  if(form == "tree") { # If the form argument is "tree" then biomass is calculated using Chave et al 2014
    
    # Calculate tree biomass using the Chave et al. 2014 allometric equation
    chave_2014 <- (0.0673* (tree_dens * dbh^2 * tree_ht)^0.976) / 1000
    return(chave_2014)
    
  } else if(form == "liana") { # If the form argument is "liana" then biomass is calculated using Schnitzer et al 2014
    
    # Calculate liana biomass using the Schnitzer et al. 2006 allometric equation
    schnitzer_2006 <- exp(-1.484 + 2.657 * log(dbh)) / 1000
    return(schnitzer_2006)
  }
}


####---- Calculate biomass of lianas (AGB + BGB) of different diameters over the simulation period ----#####

##-- Calculate diameter growth of lianas over simulation period
liana1_growth <- cumsum(case_when(yr == yr[which.min(yr)] ~ 2, TRUE ~ liana_gr)) # 2 cm diameter liana
liana2_growth <- cumsum(case_when(yr == yr[which.min(yr)] ~ 3, TRUE ~ liana_gr)) # 3 cm diameter liana
liana3_growth <- cumsum(case_when(yr == yr[which.min(yr)] ~ 4, TRUE ~ liana_gr)) # 4 cm diameter liana

##-- Calculate AGB of lianas over simulation period
liana1_AGB <- biomass_Mg(form = "liana", dbh = liana1_growth) # 2 cm diameter liana
liana2_AGB <- biomass_Mg(form = "liana", dbh = liana2_growth) # 3 cm diameter liana
liana3_AGB <- biomass_Mg(form = "liana", dbh = liana3_growth) # 4 cm diameter liana

##-- Calculate total biomass of lianas over simulation period
liana_AGB <- liana1_AGB + liana2_AGB + liana3_AGB # Add biomass of each liana for total AGB 
liana_BGB <- liana_AGB * root_factor # Calculate BGB from toal liana AGB using Mokany et al 2006 root factor
liana_BM <- liana_AGB + liana_BGB # Sum AGB and BGB for total liana biomass
liana_CO2e <- (liana_BM * c_cont) * (44/12) # Convert liana biomass to CO2e


####---- Calculate the biomass of trees (AGB + BGB) of different diameters over the simulation period ----#####

##-- Calculate control tree diameter growth over simulation period
control_growth <- cumsum(case_when(yr == 0 ~ tree_dbh, # Start with initial DBH
                                 yr >= 1 ~ tree_gr)) # Growth rate is at pre-treatment levels at all time steps

##-- Calculate treated tree diameter growth over simulation period
# Increased growth for treated tree, relative to untreated tree, are adapted from Finlayson et al. 2022.
treated_growth <- cumsum(case_when(yr == 0 ~ tree_dbh, # Start with initial DBH
                                   yr >= 1 & yr <= 10 ~ tree_gr * 2, # Double growth rate 
                                   yr >= 11 & yr <= 13 ~ tree_gr * 1.75, # Decrease growth rate by 25% every 2 yrs. 
                                   yr >= 14 & yr <= 16 ~ tree_gr * 1.5, 
                                   yr >= 17 & yr <= 19 ~ tree_gr * 1.25, 
                                   yr >= 20 ~ tree_gr)) # Growth rate is at pre-treated levels


##--- Calculate total biomass from treated and control tree diameter for each timestep of simulation period
control_BM <- biomass_Mg(form = "tree", dbh = control_growth) # control tree biomass
treated_BM_1 <- biomass_Mg(form = "tree", dbh = treated_growth) # treated tree biomass


##-- Modify treated tree biomass to include biomass lost from liana cutting
# Here liana biomass is simply subtracted from year 1 biomass increment of treated tree
treated_BM_mod <- c(treated_BM_1[1], diff(treated_BM_1)) # Calculate biomass increment, yr 0 is total starting biomass (treated_BM_1[1])
treated_BM_mod[2] <- treated_BM_mod[2] - liana_BM[1]  # Subtract liana biomass at yr0 from year 1 tree biomass increment and replace it with new value

treated_BM <- cumsum(treated_BM_mod) # recalculate cumulative sum of biomass considering liana biomass lost in year 1


####---- Convert biomass to carbon and CO2e over the simulation period ----####
# Create a dataframe with year, control biomass, and treated biomass columns.
carbon_df <- data.frame(yr, control_BM, treated_BM) %>%
  
  # Calculate carbon content as 47% of biomass
  mutate(control_C = control_BM * c_cont,
         treated_C = treated_BM * c_cont) %>%
  
  # Convert carbon content to CO2 equivalent
  mutate(treated_CO2e = (treated_C * (44/12)), 
         control_CO2e = (control_C * (44/12)), 
         
  # Calculate additional CO2e as the difference between treated and control at any given timestep
  # NOTE: additional CO2e values are cumulative across time, so should NOT be added together
         additional_CO2e = treated_CO2e - control_CO2e) # Calculate additional CO2e at each year


#### ---- Expand to CO2e estimate to per ha values & global estimate ---- ####
# Calculate Mg CO2e/ha as product of individual tree additional CO2e and the number of trees per hectare
(per_Mgha <- carbon_df$additional_CO2e[length(yr)] * n_trees_ha) # Multiply the additional CO2e by the number of trees to be treated per hectare.

# Global mitigation potential is per_Mgha * global hectares of selectively logged forest * conversion factor
(global <- (per_Mgha * 250e+6)/1e+9) # Calculate the global additional CO2e over 250 million hectares in petagrams.

# Annual mitigation potential is estimated as global potential / length of simulation
(per_yr <- global/ (length(yr)-1)) # Calculate the annual additional CO2e based on the number of years

# Total global cost of treatment is estimated based on per tree cost adapted from Mills et al. 2019, and the 
#  estimated number of trees that would be treated across selectively logged forests
(cost_per_MgCO2e <- (0.20 *( n_trees_ha * 250e+6)) / (global * 1e+9)) # Total global cost if per tree cost is $0.20


####---- LINEGRAPH: Plot treated and control CO2e against time based on number of trees per ha ----####

##-- Determine area between the lines - used to add a fill to graph
additionality <- carbon_df %>%
  mutate(treated_CO2e = treated_CO2e * n_trees_ha, 
         control_CO2e = control_CO2e * n_trees_ha) %>% 
  mutate(ymax = pmax(treated_CO2e, control_CO2e), #Calc upper and lower bounds of fill area
         ymin = pmin(treated_CO2e, control_CO2e),
         fill = treated_CO2e > control_CO2e) %>%
  mutate(fill = case_when(fill == "TRUE" ~ "Additional Carbon",
                          TRUE ~ "Lost Carbon")) %>%
  add_row(yr = 0.5, ymax = 9.2, ymin = 9.19, fill = "Additional Carbon") %>% 
  filter(fill == "Additional Carbon")


##-- Create linegraph 


# Transform the input dataframe from wide to long format and rename columns.
# Convert the tree_type column into a factor variable with appropriate levels.
# Multiply the value column by the number of trees per hectare.
(A <- carbon_df %>%  
   pivot_longer(cols = c(treated_CO2e, control_CO2e), names_to = "tree_type") %>%
   mutate(tree_type = case_when(tree_type == "control_CO2e" ~ paste("Control Trees"), TRUE ~ paste("Liberated Trees"))) %>%
   mutate(tree_type = factor(tree_type, levels = c(paste("Liberated Trees"), paste("Control Trees")))) %>%
   mutate(value = value *n_trees_ha) %>% # multiply by number of trees per ha
   
    
   # Start a ggplot2 plot with the long format dataframe as the data source.
   ggplot(data=.) +
   #geom_ribbon(data= additionality, aes(yr,ymax = ymax, ymin=ymin), fill = "chartreuse", show.legend = F) +
   geom_line(aes(x=yr, y=value, color = tree_type, linetype = tree_type), size = 1.2, show.legend = T) +
   geom_point(aes(x = 30, y = 15.26), size = 2, color = "black") +
   geom_point(aes(x = 30, y = 18.6), size = 2, color = "black") +
   geom_segment(aes(x = 30, y = 18.53, xend = 30, yend = 15.32), 
               arrow = arrow(length = unit(0.4, "cm"), type = "closed", ends = "both"), 
              linetype = "solid", size =1, color = "grey50") +
   geom_text(x = 35, y = 17, label = expression(italic("Additional carbon\nsequestered\nby year 30 =")), size = 6) +
   geom_text(x = 35.5, y = 16.4, label = expression("3.3 MgCO"[2]*"ha"^-1), size = 6) +
   scale_color_manual(values = c("black", "black")) +
   scale_linetype_manual(values = c("longdash", "solid")) +
   scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
   scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
   theme_classic() +
   coord_cartesian(clip = "off") +
   theme(legend.title = element_blank(),
         legend.position = c(0.3, 0.8),
         legend.text = element_text(size = 20, 
                                    face = "bold", 
                                    margin = margin(t = 8, b = 8, unit = "pt")),
         legend.key.width = unit(1.5, "cm"),
         axis.title = element_text(face = "bold", 
                                   size = 20),
         axis.text = element_text(size = 20),
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.line = element_line(size =1.5),
         axis.line.y = element_line(lineend = "round"),
         axis.line.x = element_line(arrow = arrow(length = unit(0.4, "cm"),
                                                  type = "closed", angle = 25)),
         plot.margin = margin(r = 6, unit = "cm")
   ) +
   xlab("Years After Liana Cutting") +
   ylab(expression(bold("Cumulative Carbon Sequestered (MgCO"[2]*"ha"^-1*")"))) 
)


# export graph
#ggsave(paste0("figures/sim_biomass_graph_fill_", Sys.Date(),".png"),
#       bg="white", dpi = 600, width = 10.6, height = 7.5)

