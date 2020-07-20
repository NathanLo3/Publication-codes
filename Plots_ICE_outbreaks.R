#********************************************************************************#
#02/07/2020
#Authors: Sindiso Nyathi, Nathan Lo.
#Goal: Generate Descriptive Plots
#********************************************************************************#

#********************************************************************************#
#Load the required libraries. This may require you to install some of these. 
these_packages <- c("naniar", "magrittr", "lubridate", "ggpubr", "reshape2", "ggplot2",  "lme4", 
                     "gridExtra", "zoo", "tidyverse")
#Load the required packages.
lapply(these_packages, require, character.only = TRUE)

#Workign directory. 
setwd("~/Box/Sindiso Nyathi's Files/RIF")

#Load the raw datafile
rif <- read.csv('rif.csv')

#Remove nonn-IHSC run centers. 
rif <- rif[-which(rif$Facility.at.Diagnosis == "Stewart Detention Center"),]
rif <- rif[-which(rif$Facility.at.Diagnosis == " Berks County Family Shelter"),]
#********************************************************************************#

#********************************************************************************#
#Some formatting.
#Reformat dates. 
rif$Date <- format(strptime(rif$Date.of.Diagnosis, format = "%m/%d/%y"), format = "%m/%d/%Y")
rif$MnYr <- format(strptime(rif$Date.of.Diagnosis, format = "%m/%d/%y"), format = "%Y/%m")
rif$Month <- month(mdy(rif$Date), label = T)
rif$Year <- year(mdy(rif$Date))

#Format Age.
rif$Age <- as.numeric(as.character(rif$Age.at.Diagnosis))

#Add a single case column. 
rif$Case = 1
#********************************************************************************#

#********************************************************************************#
#Initial Data examination. 
head(rif)
summary(rif)

#Range. 
range(rif$Date)

#ICD Groups. 
unique(rif$ICD.Group)

#Facilities 
unique(rif$Facility.at.Diagnosis)

#Set ICD as factor.
rif$ICD <- factor(rif$ICD, levels = c('Influenza', 'Varicella', 'Mumps'))

#********************************************************************************#

#********************************************************************************#
#Generate plots. 
#Aggregate by Month-year and ICD 
rif_month_icd <- aggregate(rif$Case, 
                          by = list(rif$MnYr, rif$ICD.Group, rif$Month), 
                          FUN = sum)
colnames(rif_month_icd) <- c('MnYr', 'ICD', 'Month', 'Cases')


#Add time pioints for missing months. Do this by creating a dataset witth all months, and years, and merge
#it with the original data. 
#Months
these_months <- c(1:12)
these_years <- c("2017", "2018", "2019")
ICD <- c("Varicella", "Influenza", "Mumps")
df <- tibble(
  these_years <- c("2017", "2018", "2019")
)
df2 <- df %>% expand(these_years, these_months, ICD)
df2$MnYr <- paste(df2$these_years, '/', str_pad(df2$these_months, 2, pad = '0'), '/', '01', sep = '')

#Reformat dates in df2
df2$MnYr <- format(strptime(df2$MnYr, format = "%Y/%m/%d"), format = "%Y/%m")

#Add cumulative cases variable.
rif_month_icd <- rif_month_icd[order(rif_month_icd$MnYr),]

#Merge with the ICD data. 
rif_month_icd <- merge(df2, rif_month_icd, by = c('MnYr', 'ICD'), all.x = T)

#Add in the missing months. 
rif_month_icd$Month <- na.locf(rif_month_icd$Month)
rif_month_icd$Cases[is.na(rif_month_icd$Cases)] <- 0

#Add cumulative cases. 
rif_month_icd <- rif_month_icd %>% group_by(ICD) %>% mutate(CumCase = cumsum(Cases))

#Aggregate by Month-year
rif_month <- aggregate(rif$Case, 
                       by = list(rif$MnYr, rif$Month), 
                       FUN = sum)
colnames(rif_month) <- c('MnYr', 'Month', 'Cases')

#Add cumulative cases variable.
rif_month <- rif_month[order(rif_month$MnYr),]
rif_month <- rif_month %>% mutate(CumCase = cumsum(Cases))

rif_month_icd$ICD <- factor(rif_month_icd$ICD, levels = c('Influenza', 'Varicella', 'Mumps'))

#Plot Incident Cases
incident_cases_graph <- ggplot(rif_month_icd, aes(x = MnYr, y = Cases, fill = ICD)) + 
  geom_bar(stat = 'identity', position = 'dodge', size = 0.1) + 
  scale_x_discrete(breaks = rif_month_icd$MnYr, labels = rif_month_icd$Month) + 
  scale_fill_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Incident\n Cases by ICD Code')  +
  xlab("") +
  ylab("Incident Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 150), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -20, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -20, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -20, label = "2019" , color = "black", size = 8) 
  
setEPS()
postscript("IncidentCases.eps",  width = 13, height = 10)
plot(incident_cases_graph)
dev.off()
#********************************************************************************#

#********************************************************************************#
#Cumulative Cases.
cumulative_cases_graph <- ggplot(rif_month_icd, aes(x = MnYr, y = CumCase, group = ICD, colour = ICD)) + 
  geom_point() + 
  geom_line(aes(colour = ICD)) + 
  scale_x_discrete(breaks = rif_month_icd$MnYr, labels = rif_month_icd$Month) +
  scale_colour_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Cumulative\n Cases by ICD Code')  +
  xlab("") +
  ylab("Cumulative Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 1500), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 36.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -200, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -200, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -200, label = "2019" , color = "black", size = 8) 

setEPS()
postscript("CumulativeCases.eps",  width = 13, height = 10)
plot(cumulative_cases_graph)
dev.off()
#********************************************************************************#

#********************************************************************************#
#Cases by Locatiton.
rif_location <- aggregate(rif$Case, 
                          by = list(rif$Facility.at.Diagnosis, rif$ICD.Group), 
                          FUN = sum)

colnames(rif_location) <- c('Location', 'ICD','Cases')

rif_location$ICD <- factor(rif_location$ICD, levels = c('Influenza', 'Varicella', 'Mumps'))
total_cases_location <- ggplot(rif_location, aes(x = reorder(Location, -Cases), y = Cases, fill = ICD)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Total Cases \n by Location')  +
  xlab("") +
  ylab("Total Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 0.5, 0.8),"cm"))

setEPS()
postscript("LocationCases.eps",  width = 13, height = 10)
plot(total_cases_location)
dev.off()
#********************************************************************************#

#********************************************************************************#
#Ages. 
rif$ICD <- factor(rif$ICD, levels = c('Influenza', 'Varicella', 'Mumps'))
total_cases_age <- ggplot(rif, aes(x = Age, fill = ICD)) + 
  geom_histogram( colour = 'grey20', size = 0.1, binwidth = 2.5) + 
  scale_fill_manual(name = "ICD Code", values = c("indianred2",  "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Total Cases \n by Age')  +
  xlab("Age") +
  ylab("Total Cases") +
  ylim(0, 350) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 0.5, 0.5),"cm"))

setEPS()
postscript("AgeCases.eps",  width = 13, height = 10)
plot(total_cases_age)
dev.off()

#********************************************************************************#

#********************************************************************************#
#Look at a few specific locations. 
rif_south_tx <- rif[rif$Facility.at.Diagnosis == ' South Texas Family RC',]

rif_south_tx <- aggregate(rif_south_tx$Case, 
                           by = list(rif_south_tx$MnYr, rif_south_tx$ICD.Group, rif_south_tx$Month), 
                           FUN = sum)
colnames(rif_south_tx) <- c('MnYr', 'ICD', 'Month', 'Cases')
rif_south_tx <- rif_south_tx[order(rif_south_tx$MnYr),]
rif_south_tx <- rif_south_tx %>% group_by(ICD) %>% mutate(CumCase = cumsum(Cases))

#Plot Incident Cases
incident_cases_graph_south_tx <- ggplot(rif_south_tx, aes(x = MnYr, y = Cases, group = ICD, colour = ICD)) + 
  geom_point() + 
  geom_line(aes(colour = ICD)) + 
  scale_x_discrete(breaks = rif_month_icd$MnYr, labels = rif_month_icd$Month) +
  scale_colour_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Incident Cases by ICD Code\n South Texas Family RC')  +
  xlab("") +
  ylab("Incident Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 150), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 36.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -20, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -20, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -20, label = "2019" , color = "black", size = 8) 

setEPS()
postscript("IncidentCases. South Texas.eps",  width = 13, height = 10)
plot(incident_cases_graph_south_tx)
dev.off()

#Cumulative Cases.
cumulative_cases_graph_south_tx <- ggplot(rif_south_tx, aes(x = MnYr, y = CumCase, group = ICD, colour = ICD)) + 
  geom_point() + 
  geom_line(aes(colour = ICD)) + 
  scale_x_discrete(breaks = rif_month_icd$MnYr, labels = rif_month_icd$Month) +
  scale_colour_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Cumulative Cases by ICD Code\n South Texas Family RC')  +
  xlab("") +
  ylab("Cumulative Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 1500), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 36.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -200, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -200, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -200, label = "2019" , color = "black", size = 8) 

setEPS()
postscript("CumulativeCases. South Texas.eps",  width = 13, height = 10)
plot(cumulative_cases_graph_south_tx)
dev.off()


#Look at a few specific locations. 
rif_port_isa <- rif[rif$Facility.at.Diagnosis == ' Port Isabel SPC',]

rif_port_isa <- aggregate(rif_port_isa$Case, 
                          by = list(rif_port_isa$MnYr, rif_port_isa$ICD.Group, rif_port_isa$Month), 
                          FUN = sum)
colnames(rif_port_isa) <- c('MnYr', 'ICD', 'Month', 'Cases')
rif_port_isa <- rif_port_isa[order(rif_port_isa$MnYr),]
rif_port_isa <- rif_port_isa %>% group_by(ICD) %>% mutate(CumCase = cumsum(Cases))

#Plot Incident Cases
incident_cases_graph_port_isa <- ggplot(rif_port_isa, aes(x = MnYr, y = Cases, group = ICD, colour = ICD)) + 
  geom_point() + 
  geom_line(aes(colour = ICD)) + 
  scale_x_discrete(breaks = rif_port_isa$MnYr, labels = rif_port_isa$Month) +
  scale_colour_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Incident Cases by ICD Code\n Port Isabel')  +
  xlab("") +
  ylab("Incident Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 150), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 36.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -20, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -20, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -20, label = "2019" , color = "black", size = 8) 

setEPS()
postscript("IncidentCases. Port Isabel.eps",  width = 13, height = 10)
plot(incident_cases_graph_port_isa)
dev.off()

#Cumulative Cases.
cumulative_cases_graph_port_isa <- ggplot(rif_port_isa, aes(x = MnYr, y = CumCase, group = ICD, colour = ICD)) + 
  geom_point() + 
  geom_line(aes(colour = ICD)) + 
  scale_x_discrete(breaks = rif_port_isa$MnYr, labels = rif_port_isa$Month) +
  scale_colour_manual(name = "ICD Code", values = c("indianred2", "royalblue1", "mediumseagreen")) +
  theme_bw() +
  ggtitle('ICE Center Cumulative Cases by ICD Code \n Port Isabel')  +
  xlab("") +
  ylab("Cumulative Cases") +
  #ylim(0, 150) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 30, hjust = 0.5), 
        axis.title.x = element_text(size = 26), 
        axis.title.y = element_text(size = 26), 
        panel.grid.minor.y = element_line(), 
        panel.grid.major.x = element_blank(), 
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 26), 
        plot.margin=unit(c(0.2, 0.5, 1, 0.5),"cm")) +
  coord_cartesian(ylim = c(0, 1500), clip = "off", expand = T) +
  geom_vline(xintercept = 12.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 24.5, linetype = 6, color = "gray50", size = 0.5) +
  geom_vline(xintercept = 36.5, linetype = 6, color = "gray50", size = 0.5) +
  annotate(geom = "text", x = 6.5, y = -200, label = "2017" , color = "black", size = 8) +
  annotate(geom = "text", x = 18.5, y = -200, label = "2018" , color = "black", size = 8) +
  annotate(geom = "text", x = 30.5, y = -200, label = "2019" , color = "black", size = 8) 

setEPS()
postscript("CumulativeCases. Port Isabel.eps",  width = 13, height = 10)
plot(cumulative_cases_graph_port_isa)
dev.off()

#Fin. 
#********************************************************************************#

#********************************************************************************#


# Additional descriptive analysis



# Upload data
setwd("~/Dropbox/Lo collaborations/ICE/Plots/")

ICE_data <- read.csv(file="ICE_data.csv")
ICE_data$Age <- as.numeric(ICE_data$Age.at.Diagnosis)
ICE_data$Facility <- as.character(ICE_data$Facility.at.Diagnosis)

# Remove “Stewart detention center” and “Berks County”, as they are not IHSC run facilitie
ICE_data <- ICE_data[ICE_data$Facility!=" Berks County Family Shelter" & ICE_data$Facility!="Stewart Detention Center",]

mean(ICE_data$Age[ICE_data$ICD.Group=="Influenza"])
quantile(ICE_data$Age[ICE_data$ICD.Group=="Influenza"])
sd(ICE_data$Age[ICE_data$ICD.Group=="Influenza"])
unique(ICE_data$Facility[ICE_data$ICD.Group=="Influenza"])

mean(ICE_data$Age[ICE_data$ICD.Group=="Varicella"])
quantile(ICE_data$Age[ICE_data$ICD.Group=="Varicella"])
sd(ICE_data$Age[ICE_data$ICD.Group=="Varicella"])
unique(ICE_data$Facility[ICE_data$ICD.Group=="Varicella"])

mean(ICE_data$Age[ICE_data$ICD.Group=="Mumps"])
quantile(ICE_data$Age[ICE_data$ICD.Group=="Mumps"])
sd(ICE_data$Age[ICE_data$ICD.Group=="Mumps"])
unique(ICE_data$Facility[ICE_data$ICD.Group=="Mumps"])


# Extra
jpeg("Fig_flu_1.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" Port Isabel SPC" & ICE_data$ICD.Group=="Influenza"], 100,main="Port Isabel SPC, Flu", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()
jpeg("Fig_flu_2.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" South Texas Family Residential Center" & ICE_data$ICD.Group=="Influenza"], 100,main="South Texas, Flu", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()

jpeg("Fig_mumps_1.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" Port Isabel SPC" & ICE_data$ICD.Group=="Mumps"], 100,main="Port Isabel SPC, Mumps", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()
jpeg("Fig_mumps_2.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" South Texas Family Residential Center" & ICE_data$ICD.Group=="Mumps"], 100,main="South Texas, Mumps", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()

jpeg("Fig_varicella_1.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" Port Isabel SPC" & ICE_data$ICD.Group=="Varicella"], 100,main="Port Isabel SPC, Varicella", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()
jpeg("Fig_varicella_2.jpg")
hist(ICE_data$Age.at.Diagnosis[ICE_data$Facility.at.Diagnosis==" South Texas Family Residential Center" & ICE_data$ICD.Group=="Varicella"], 100,main="South Texas, Varicella", xlab="Age", xlim=c(0,100), ylim=c(0,75))
dev.off()
