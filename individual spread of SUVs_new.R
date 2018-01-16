setwd("//zkh/dfs/Gebruikers17/KahleX/desktop/R-projects/MYC-PET")


library(ggplot2)
library(stringr)
library(readxl)
library(utils)
library(dplyr)
library(extrafont)
library(survival)

SUV_data <- read_excel("SUVs_multiple lesions.xlsx")

colnames(SUV_data)[c(2,3,4,5,6,7,10,12,14)] <- c("maxs", "means", "mins", "MATV", "MYC", "range_SUV", "code_OS", "code_PFS", "code_DSS")

SUV_data <- SUV_data %>% mutate(OS.month = OS / (365/12)) %>%
                          mutate(PFS.month = PFS / (365/12)) %>% 
                          mutate(DSS.month = DSS / (365/12))


SUV_data_myc_pos <-  SUV_data %>% filter(MYC %in% 'pos')
SUV_data_myc_neg <-  SUV_data %>% filter(MYC %in% 'neg')

summary(SUV_data)  

cor.test(SUV_data$MATV,SUV_data$range_SUV, data = SUV_data, method = "spearman", continuity = F, conf.level = 0.95)

wilcox.test(SUV_data_myc_neg$range_SUV, SUV_data_myc_pos$range_SUV, paired = F)

cor.test(SUV_data$num_lesions, SUV_data$range_SUV, data = SUV_data, method = "spearman", continuity = F, conf.level = 0.95)
cor.test(SUV_data$num_lesions, SUV_data$MATV, data = SUV_data, method = "spearman", continuity = T, conf.level = 0.95)
cor.test(SUV_data$MATV, SUV_data$means, method = "spearman", continuity = T, conf.level = 0.95)
cor.test(SUV_data$DSS.month,SUV_data$range_SUV, data = SUV_data, method = "spearman", continuity = F, conf.level = 0.95)

ggplot(SUV_data,aes(y = MATV, x=range_SUV))+
  geom_point()+
  coord_trans(y = "log10")+
  geom_smooth(method = "lm", se = F)


ggplot(SUV_data)+
  geom_density(aes(range_SUV))
  
summary(SUV_data)

# FINAL PLOT-------------------

brk <- 10^(-10:10)
minor_brk <- rep(1:9, 21)*(10^rep(-10:10, each=9))

MYC_labels <- c("neg" = "MYC-negative", "pos" = "MYC-positive")

dotsandbars <- ggplot(SUV_data) +
  geom_point(aes(x = MATV, y = means))+
  geom_errorbar(aes(x = MATV, ymin = mins, ymax = maxs), alpha = 0.5)+
  facet_grid(.~MYC, labeller = as_labeller(MYC_labels))+
  scale_x_log10(breaks= brk, minor_breaks= minor_brk) +
  scale_color_manual(values = c("royalblue3", "red"))+
  labs(x = expression(paste("MATV in mL ",(log [10]))),
       y = expression(paste(italic("SUV"[max])," (mean + range)")),
       title = expression(paste("Intra-individual spread of ", italic("SUV"[max]),
                                " according to ", italic("MYC-"), "status in patients with multiple lesions")))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        text=element_text(size=12, face = "plain", family="Arial"))

final <- dotsandbars + annotation_logticks(sides = "b", scaled = T)

final


pdf("Fig1.pdf", height = 5, width = 8,
    family = "Arial", paper = "special",
    onefile = F)
final
dev.off()

tiff("Fig1.tiff", height= 5, width = 8, units = 'in', family = "Arial", res = 1200)
final
dev.off()

Sys.setenv(R_GSCMD="C:/Program Files/gs/gs9.22/bin/gswin64.exe")

embed_fonts("./Fig1.pdf")

loadfonts(device = "win")

postscript("Fig1.eps", height = 10, width = 16,
           family = "Arial", paper = "special",
           onefile = F, horizontal = F)
final
dev.off


# OLD versions--------------

dotsandbars <- ggplot(SUV_data, aes(col = MYC)) +
  geom_point(aes(x = log10(MATV), y = means))+
  geom_errorbar(aes(x = log10(MATV), ymin = mins, ymax = maxs))+
  #geom_tile (aes(fill= MYC))+
  facet_grid(.~MYC)+
  scale_x_continuous(name = "MATV (log10)", labels = scales::math_format(10^.x)) +
  scale_color_manual(values = c("royalblue3", "red"))+
  scale_y_continuous("mean SUVmax (+ range)")+
  ggtitle("Intra-individual spread of SUVmax \naccording to MYC-status in patients with multiple lesions ")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.position = "none", plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))

dotsandbars + annotation_logticks(sides = "b", scaled = T)


Lines <- ggplot(SUV_data, aes(fill= MYC, color = MYC)) +
  geom_line(aes(MATV, means), size = 1) +
  geom_errorbar(aes(x = MATV, ymin = mins, ymax = maxs))+
  # geom_point(aes(MATV, means))+
  geom_ribbon(aes(x = MATV, ymin = mins, ymax = maxs),alpha= 0.4) +
  coord_trans(x = "log10")+
  facet_grid(.~MYC)+
  scale_color_manual(values = c("royalblue3", "red"))+
  scale_fill_manual(values = c("royalblue3", "red"))+
  ylab("mean SUVmax (+ range)")+
  xlab("MATV (log10)")+
  ggtitle("Intra-individual spread of SUVmax according to MYC-status", subtitle = "only cases with multiple lesions are displayed (n = 43), each vertical line represents the spread of SUVmax across lesions within 1 patient")+
  theme_bw()+
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 1))

Lines

# SURVIVAL ANALYSIS---------

SUV_data$range_SUV_cat <- ifelse(SUV_data$range_SUV >= 12, ">=12", "<12")

uni.cox.os <- coxph(Surv(OS.month, code_OS) ~ range_SUV_cat, data =  SUV_data)
summary(uni.cox.os)

uni.cox.pfs <- coxph(Surv(PFS.month, code_PFS) ~ range_SUV_cat, data =  SUV_data)
summary(uni.cox.pfs)

uni.cox.dss <- coxph(Surv(PFS.month, code_DSS) ~ range_SUV_cat, data =  SUV_data)
summary(uni.cox.dss)


surv.curve.SUV.os <- survfit(Surv(SUV_data$OS.month,SUV_data$code_OS) ~ SUV_data$range_SUV_cat)

A <- ggsurvplot(
  surv.curve.SUV.os,
  data = SUV_data,
  break.y.by = 0.1,
  break.time.by = 6,
  linetype = c(1,4),
  size = 1.2,
  xlim = c(0,64),
  ylim = c(0,1.05),
  legend = c(0.9,0.2),
  legend.title = expression(paste("range of",italic(" SUV"))),
  legend.labs = c("< 12",">= 12"),
  palette = c("black","black"),
  pval = T,
  pval.method = F,
  axes.offset = F,
  censor = T,
  censor.shape = 124,
  censor.size = 3,
  title = expression(paste("Overall survival according to range of ",italic("SUVmax"))),
  xlab = "Time (months)",
  ylab = "Survival probability",
  font.family = "Arial") 

A$plot <-  A$plot + theme(text=element_text(size=12, face = "plain", family="Arial"),
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size = 10),
                          axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 10),
                          axis.title.x = element_text(size = 12),
                          axis.title.y = element_text(size = 12))

A

surv.curve.SUV.pfs <- survfit(Surv(SUV_data$PFS.month,SUV_data$code_PFS) ~ SUV_data$range_SUV_cat)

B <- ggsurvplot(
  surv.curve.SUV.pfs,
  data = SUV_data,
  break.y.by = 0.1,
  break.time.by = 6,
  linetype = c(1,4),
  size = 1.2,
  xlim = c(0,64),
  ylim = c(0,1.05),
  legend = c(0.9,0.2),
  legend.title = expression(paste("range of",italic(" SUV"))),
  legend.labs = c("< 12",">= 12"),
  palette = c("black","black"),
  pval = T,
  pval.method = F,
  axes.offset = F,
  censor = T,
  censor.shape = 124,
  censor.size = 3,
  title = expression(paste("Progression free survival according to range of ",italic("SUVmax"))),
  xlab = "Time (months)",
  ylab = "Survival probability",
  font.family = "Arial") 

B$plot <-  B$plot + theme(text=element_text(size=12, face = "plain", family="Arial"),
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size = 10),
                          axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 10),
                          axis.title.x = element_text(size = 12),
                          axis.title.y = element_text(size = 12))

B

surv.curve.SUV.dss <- survfit(Surv(SUV_data$DSS.month,SUV_data$code_DSS) ~ SUV_data$range_SUV_cat)

C <- ggsurvplot(
  surv.curve.SUV.dss,
  data = SUV_data,
  break.y.by = 0.1,
  break.time.by = 6,
  linetype = c(1,10),
  size = 1.0,
  xlim = c(0,64),
  ylim = c(0,1.05),
  legend = c(0.9,0.2),
  legend.title = expression(paste("range of",italic(" SUVmax"))),
  legend.labs = c("< 12",">= 12"),
  palette = c("black","black"),
  pval = T,
  pval.method = F,
  axes.offset = F,
  censor = T,
  censor.shape = 124,
  censor.size = 3,
  title = expression(paste("Disease specific survival according to range of ",italic("SUVmax"))),
  xlab = "Time (months)",
  ylab = "Survival probability",
  font.family = "Arial") 

C$plot <-  C$plot + theme(text=element_text(size=12, face = "plain", family="Arial"),
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size = 10),
                          axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 10),
                          axis.title.x = element_text(size = 12),
                          axis.title.y = element_text(size = 12))

C

