#### install and load packages ----
install.packages("pacman")
pacman::p_load(tidyverse, cowplot, car, ggpubr, ggpmisc, ggtext, segmented, reshape2, nlme, plantecophys, hms, mgcv, DescTools, gridExtra)

# install.packages("remotes")
# remotes::install_github("OnofriAndreaPG/aomisc", force=T)

#### import data ----
Species_full_list <- read.csv("Data//Species_names.csv")
Species <- c("AB", "EF", "AF", "ER", "AE", "PL", "HD")
Species_names <- data.frame(Full=c("Angophora bakeri", "Acacia elongata", "Acacia falcata", "Eucalyptus fibrosa", "Eucalyptus racemosa", "Hakea dactyloides", "Persoonia linearis"), Spp=c("AB", "AE", "AF", "EF", "ER", "HD", "PL"))
diel_12hr <- read.csv("Data//diel_12hr.csv")
time_values <- read.csv("Data//time_values.csv")
dat <- read.csv("Data//CNR_data.csv")
dat$Tau <- dat$Cs*((dat$Height)^2)/(dat$Ks)/60
timecurves_values <- read.csv("Data//timecurves_values.csv")
timecurves_values$Spp <- as.factor(timecurves_values$Spp)
time_values$Spp <- as.factor(time_values$Spp)
diel_12hr$Spp <- as.factor(diel_12hr$Spp)
diel_12hr$Day <- as.factor(diel_12hr$Day)
diel_12hr$hhmmss_12hr <- as.POSIXct(diel_12hr$hhmmss, format="%H:%M:%S")
diel_12hr2 <- diel_12hr

#### ggplot theme and colour palette ----
myggplottheme <- theme_bw(base_size=20) +
  theme(axis.title=element_text(size=35),
        axis.text=element_text(size=30),
        legend.position="None",
        # legend.position=c(0.8,0.2),
        # legend.title=element_text(size=30),
        # legend.text=element_text(size=27),
        axis.ticks.length=unit(-0.5, "cm"), 
        panel.grid=element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent"),
        panel.border=element_rect(linewidth=2))

cbp <- c("AB"="#D55E00", "ER"="#E69F00", "AF"="#F0E442",
         "EF"="#009E73", "AE"="#56B4E9", "PL"="#0072B2",
         "HD"="#A020F0") # pink="#CC79A7", black="#000000", grey="#999999"; colourblind-friendly palette
cbp_number <- c("1"="#CC79A7", "2"="#999999", "3"="#F0E442",
         "4"="#009E73", "5"="#56B4E9", "6"="#A020F0")
# pink="#CC79A7", black="#000000", grey="#999999"; colourblind-friendly palette

p_legend <- ggplot(time_values, aes(x=Cs, y=minPsileaf_time_gam, color=Spp))+
  geom_point()+
  lims(x = c(0,0), y = c(0,0))+
  scale_fill_manual(values=cbp, aesthetics="colour", name="species", labels=c("*Angophora bakeri*", "*Acacia elongata*", "*Acacia falcata*", "*Eucalyptus fibrosa*", "*Eucalyptus racemosa*", "*Hakea dactyloides*", "*Persoonia linearis*")) +
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_markdown(size=27),
        legend.title = element_text(size=30))+
  guides(colour=guide_legend(override.aes=list(size=8))); p_legend
ggsave("Output//species-legend-ms.png", plot=p_legend, width=250, height=200, units="mm")

#### Range of values ----
(Cs <- c(min(time_values$Cs, na.rm=T), max(time_values$Cs, na.rm=T), max(time_values$Cs, na.rm=T)/min(time_values$Cs, na.rm=T)))
(Tau <- c(min(time_values$Tau, na.rm=T), max(time_values$Tau, na.rm=T), max(time_values$Tau, na.rm=T)/min(time_values$Tau, na.rm=T)))
(d13C <- c(min(dat$d13C, na.rm=T), max(dat$d13C, na.rm=T), min(dat$d13C, na.rm=T)/max(dat$d13C, na.rm=T)))

#### Fig. 1: data fitting example ----
p1 <- ggplot(diel_12hr2[diel_12hr2$Spp=="AB",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), linewidth=1.5, colour="black", linetype="solid") +
  # geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["AB"]]) +
  labs(title=expression(italic("Angophora bakeri")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  # coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25),
        legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p1
ggsave("Output//diurnal_fits_AB.png", width=250, height=200, units="mm", limitsize=F)

p2 <- ggplot(diel_12hr2[diel_12hr2$Spp=="PL",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour="black", linetype="dashed") +
  # geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["PL"]]) +
  labs(title=expression(italic("Persoonia linearis")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-4, hjust=0.95, size=25)); p2
ggsave("Output//diurnal_fits_PL.png", width=250, height=200, units="mm", limitsize=F)

p3 <- ggplot(diel_12hr2[diel_12hr2$Spp=="HD",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["HD"]]) +
  labs(title=expression(italic("Hakea dactyloides")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  # coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25),
        legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p3
ggsave("Output//diurnal_fits_HD.png", width=250, height=200, units="mm", limitsize=F)

p4 <- ggplot(diel_12hr2[diel_12hr2$Spp=="ER",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["ER"]]) +
  labs(title=expression(italic("Eucalyptus racemosa")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25)); p4
ggsave("Output//diurnal_fits_ER.png", width=250, height=200, units="mm", limitsize=F)

p5 <- ggplot(diel_12hr2[diel_12hr2$Spp=="EF",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["EF"]]) +
  labs(title=expression(italic("Eucalyptus fibrosa")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  # coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25),
        legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p5
ggsave("Output//diurnal_fits_EF.png", width=250, height=200, units="mm", limitsize=F)

p6 <- ggplot(diel_12hr2[diel_12hr2$Spp=="AF",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["AF"]]) +
  labs(title=expression(italic("Acacia falcata")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25)); p6
ggsave("Output//diurnal_fits_AF.png", width=250, height=200, units="mm", limitsize=F)

p7 <- ggplot(diel_12hr2[diel_12hr2$Spp=="AE",], aes(y=A, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["AE"]]) +
  labs(title=expression(italic("Acacia elongata")), y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  # coord_cartesian(ylim=c(0,18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25),
        legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p7
ggsave("Output//diurnal_fits_AE.png", width=250, height=200, units="mm", limitsize=F)

## diurnal curves for stomatal conductance
# p3 <- ggplot(diel_12hr2[diel_12hr2$Spp=="AB",], aes(y=gsw, x=hhmmss_int)) +
#   geom_point(size=3, aes(colour=Day)) +
#   geom_smooth(method="gam", formula=y~s(x), size=1.5, colour="green") +
#   labs(y=expression(g[SW]~"mol/m"^2*"/s)"), x="time of day (h)") + # μ
#   scale_x_continuous(breaks=c(6,9,12,15,18)) +
#   myggplottheme; p3
# 
# p4 <- ggplot(diel_12hr2[diel_12hr2$Spp=="PL",], aes(y=gsw, x=hhmmss_int)) +
#   geom_point(size=3, aes(colour=Day)) +
#   geom_smooth(method="gam", formula=y~s(x), size=1.5, colour="green") +
#   labs(y=expression(g[SW]~"mol/m"^2*"/s)"), x="time of day (h)") + # μ
#   scale_x_continuous(breaks=c(6,9,12,15,18)) +
#   myggplottheme; p4

p3 <- ggplot(diel_12hr2[diel_12hr2$Spp=="AB",], aes(y=Psi_xylem, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour="black", linetype="solid") +
  # geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["AB"]]) +
  labs(title=expression(italic("Angophora bakeri")), y=expression(Psi["X"]~(MPa)), x="time of day (h)") + # μ
  coord_cartesian(ylim=c(-1.6,0)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-8, hjust=0.95, size=25)); p3

p4 <- ggplot(diel_12hr2[diel_12hr2$Spp=="PL",], aes(y=Psi_xylem, x=hhmmss_int)) +
  geom_point(size=3, aes(colour=Day), alpha=0.8) +
  geom_smooth(method="gam", formula=y~s(x), size=1.5, colour="black", linetype="dashed") +
  # geom_smooth(method="gam", formula=y~s(x), size=1.5, colour=cbp[["PL"]], linetype="dashed") +
  labs(title=expression(italic("Persoonia linearis")), y=expression(Psi["X"]~(MPa)), x="time of day (h)") + # μ
  coord_cartesian(ylim=c(-1.6,0)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  # scale_shape_manual(values=c(0,1,2,5,6,10))+
  scale_colour_manual(values=cbp_number)+
  myggplottheme +
  theme(plot.title = element_text(vjust=-4, hjust=0.95, size=25)); p4

ggarrange(p1, p2, p3, p4,
          # p1+rremove("xlab"), p2+rremove("xlab"),
          ncol=2, nrow=2, labels=c("C","D","E","F"),
          label.x=0.018, label.y=0.95, common.legend=T, legend = "bottom",
          font.label=list(size=30), align="hv")

ggsave("Output//diurnal_fits_examples.png", width=600, height=400, units="mm", limitsize=F)

#### Fig. 2: time curves comparing two species: AB and PL ----
p1 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=A_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p1

p2 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=A, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("A (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p2

p3 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=gsw_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(g[S]*" (mol/m"^2*"/s)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p3

p4 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=gsw, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(g[S]*" (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p4

p5 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=1000*E_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("E (mmol/m"^2*"/s)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p5

p6 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=E, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("E (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p6

p7 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=VPDleaf_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("*D* (kPa)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  scale_y_continuous(breaks=c(0.5,1,1.5,2,2.5)) +
  myggplottheme +
  theme(axis.title.y = ggtext::element_markdown()); p7

p8 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=VPDleaf, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("*D* (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(axis.title.y = ggtext::element_markdown()); p8

p9 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=WUE_abs/1000, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("WUE (ppm)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p9

p10 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=WUE, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression("WUE (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p10

p11 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=Psi_xylem_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(Psi[X]*" (MPa)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p11

p12 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=Psi_xylem, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(Psi[X]*" (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p12

p13 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=Psi_leaf_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(Psi[L]*" (MPa)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p13

p14 <- ggplot(subset(timecurves_values, Spp==c("AB", "PL")), aes(y=Psi_leaf, x=time)) +
  geom_line(size=2, aes(colour=Spp, linetype=Spp)) + 
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(Psi[L]*" (scaled)"), x="time of day (h)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p14

ggarrange(p1+rremove("xlab"), p2+rremove("xlab"),
          p3+rremove("xlab"), p4+rremove("xlab"),
          p5+rremove("xlab"), p6+rremove("xlab"),
          p7+rremove("xlab"), p8+rremove("xlab"),
          p9+rremove("xlab"), p10+rremove("xlab"),
          p11+rremove("xlab"), p12+rremove("xlab"),
          p13, p14,
          common.legend=T, legend="top",
          ncol=2, nrow=7, labels="AUTO", label.x=0.018, label.y=0.95,
          font.label=list(size=30), align="hv")

ggsave("Output//diurnal_traits_2spp_ms.png", width=500, height=1400, units="mm", limitsize=F)

#### Fig. 3: max/min time vs. Cs ----
nlmod1<-nls(maxA_time_gam~s*Cs/(t+Cs), data=time_values, start=list(t=20, s=12))
# nlmod1<-nls(maxA_time_gam~d*(1-exp(-a*Cs/d))-b, data=time_values, start=list(d=.1, a=-.2, b=2))
summary(nlmod1)
modelr::rsquare(nlmod1, time_values)

out.lm_1 <- lm(maxA_time_gam~Cs, data=time_values)
o_1 <- segmented(out.lm_1, seg.Z = ~Cs)
slope1_1 <- slope(o_1)$Cs[1]
slope2_1 <- slope(o_1)$Cs[2]
CI_1 <- o_1$psi[2]-c(1.96,-1.96)*o_1$psi[3]
inflection_1 <- data.frame(x=c(min(time_values$Cs),o_1$psi[2]), xend=c(o_1$psi[2],max(time_values$Cs)), y=c(min(time_values$Cs)*slope1_1+intercept(o_1)$Cs[1],o_1$psi[2]*slope2_1+intercept(o_1)$Cs[2]), yend=c(o_1$psi[2]*slope1_1+intercept(o_1)$Cs[1],max(time_values$Cs)*slope2_1+intercept(o_1)$Cs[2]))

p1 <- ggplot(time_values, aes(y=(maxA_time_gam), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method="lm", colour="black", se=F) +
  # geom_smooth(method="nls",formula=y~s*x/(F+x),method.args=list(start = c(F=20,s=12)), colour="black", se=FALSE) +
  ggpmisc::stat_poly_eq(formula=y~x,aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  # ggpmisc::stat_poly_eq(aes(formula=y~x,label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
  #                       parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  # annotate(geom="text", label="pseudo R² = 0.68", x=400, y=8.8, size=9) +
  geom_segment(data = inflection_1, aes(x = x, xend = xend, y = y, yend = yend), size=1.5, linetype="dashed") +
  labs(y=expression(time["A-max"]~"(o'clock)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  scale_fill_manual(values=cbp,aesthetics="colour") +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p1

nlmod2<-nls(maxgs_time_gam~s*Cs/(F+Cs), data=time_values, start=list(F=20, s=12))
summary(nlmod2)
modelr::rsquare(nlmod2, time_values)

out.lm_2 <- lm(maxgs_time_gam~Cs, data=time_values)
o_2 <- segmented(out.lm_2, seg.Z = ~Cs)
slope1_2 <- slope(o_2)$Cs[1]
slope2_2 <- slope(o_2)$Cs[2]
CI_2 <- o_2$psi[2]-c(1.96,-1.96)*o_2$psi[3]
inflection_2 <- data.frame(x=c(min(time_values$Cs),o_2$psi[2]), xend=c(o_2$psi[2],max(time_values$Cs)), y=c(min(time_values$Cs)*slope1_2+intercept(o_2)$Cs[1],o_2$psi[2]*slope2_2+intercept(o_2)$Cs[2]), yend=c(o_2$psi[2]*slope1_2+intercept(o_2)$Cs[1],max(time_values$Cs)*slope2_2+intercept(o_2)$Cs[2]))

p2 <- ggplot(time_values, aes(y=(maxgs_time_gam), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  # geom_smooth(method="nls",formula=y~s*x/(F+x),method.args=list(start = c(F=20,s=12)), colour="black", se=FALSE) +
  ggpmisc::stat_poly_eq(formula=y~x,aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  # annotate(geom="text", label="pseudo R² = 0.59", x=400, y=9, size=9) +
  geom_segment(data = inflection_2, aes(x = x, xend = xend, y = y, yend = yend), size=1.5, linetype="dashed") +
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(time["g"[S]*"-max"]~"(o'clock)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p2

p3 <- ggplot(time_values, aes(y=(maxE_time_gam), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  # geom_smooth(method = "lm", colour="black", se=F, linetype="dotted") +
  # geom_smooth(method="nls",formula=y~s*x/(F+x),method.args=list(start = c(F=20,s=12)), colour="black", se=FALSE) +
  # ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
  #                       parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(time["E-max"]~"(o'clock)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p3

nlmod4<-nls(minPsixylem_time_gam~s*Cs/(F+Cs), data=time_values, start=list(F=20, s=12))
summary(nlmod4)
modelr::rsquare(nlmod4, time_values)

out.lm_4 <- lm(minPsixylem_time_gam~Cs, data=time_values)
o_4 <- segmented(out.lm_4, seg.Z = ~Cs)
slope1_4 <- slope(o_4)$Cs[1]
slope2_4 <- slope(o_4)$Cs[2]
CI_4 <- o_4$psi[2]-c(1.96,-1.96)*o_4$psi[3]
inflection_4 <- data.frame(x=c(min(time_values$Cs),o_4$psi[2]), xend=c(o_4$psi[2],max(time_values$Cs)), y=c(min(time_values$Cs)*slope1_4+intercept(o_4)$Cs[1],o_4$psi[2]*slope2_4+intercept(o_4)$Cs[2]), yend=c(o_4$psi[2]*slope1_4+intercept(o_4)$Cs[1],max(time_values$Cs)*slope2_4+intercept(o_4)$Cs[2]))

p4 <- ggplot(time_values, aes(y=(minPsixylem_time_gam), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  # geom_smooth(method="nls",formula=y~s*x/(F+x),method.args=list(start = c(F=20,s=12)), colour="black", se=FALSE) +
  ggpmisc::stat_poly_eq(formula=y~x,aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  # annotate(geom="text", label="pseudo R² = 0.84", x=400, y=11.8, size=9) +
  geom_segment(data = inflection_4, aes(x = x, xend = xend, y = y, yend = yend), size=1.5, linetype="dashed") +
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(time[Psi[X]*"-min"]~"(o'clock)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p4

p5 <- ggplot(time_values, aes(y=(minPsileaf_time_gam), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  # geom_smooth(method = "lm", colour="black", se=F, linetype="dotted") +
  # geom_smooth(method="nls",formula=y~s*x/(F+x),method.args=list(start = c(F=20,s=12)), colour="black", se=FALSE) +
  # ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
  #                       parse = TRUE,label.x="right", label.y=0.1, size=12) + # method="spearman",
  scale_fill_manual(values=cbp,aesthetics="colour") +
  labs(y=expression(time[Psi[L]*"-min"]~"(o'clock)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim=c(50,800), ylim=c(11,14.2)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p5

ggarrange(p1+rremove("xlab"), p2+rremove("xlab"),
          p3+rremove("xlab"), p4, p5,
          ncol=2, nrow=3, legend="none",
          labels="AUTO", label.x=0.018, label.y=0.95,
          font.label=list(size=30), align="hv")
ggsave("Output//max_min_time-Cs-ms.png", width=500, height=600, units="mm")

#### Fig. 4: hydraulic time constant vs. sapwood capacitance ----
component_Cs <- lm(log(Tau)~log(Cs), data=dat); summary(component_Cs)
p9 <- ggplot(dat, aes(x=Cs, y=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F, size=1.5) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(eq.label),sep="")),
                        formula=y~x,label.x=0.05, label.y=0.95, size=12, parse=T) + # method="spearman",
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x=0.05, label.y=0.85, size=12, parse=T) + # method="spearman",
  labs(y=expression(tau~"(minutes)"), x=expression(C["S"]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(ylim=c(0.5, 100)) +
  scale_x_log10(breaks=c(10,20,30,40,50,60)) +
  scale_y_log10(breaks=c(1,10,100)) +
  myggplottheme; p9

ggsave("Output//tau-Cs.png", width=250, height=200, units="mm")

#### Fig. 5: (relative scale) stomatal diurnal dynamics ----

p4 <- ggplot(time_values, aes(y=100*(relgsw_11-relgsw_12), x=Tau)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x="right", label.y=0.9, size=12, parse=T) + # method="spearman",
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Delta*g["S,11-12"]~"(%)"), x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p4

p5 <- ggplot(time_values, aes(y=100*(absAarea_11-absAarea_12)/absAarea_11, x=Tau)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x="right", label.y=0.9, size=12, parse=T) + # method="spearman",
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Delta*A["11-12"]~"(%)"), x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p5

# p6 <- ggplot(time_values, aes(y=100*(relWUE_12-relWUE_13), x=Tau)) +
#   geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
#   geom_smooth(method = "lm", colour="black", se=F) +
#   ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
#                         formula=y~x,label.x="left", label.y=0.9, size=12, parse=T) + # method="spearman",
#   # labs(y=expression(g[s12]-g[s11]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
#   scale_fill_manual(values=cbp, aesthetics="colour") +
#   labs(y=expression(Delta*WUE["12-13"]~"(%)"), x=expression(tau~"(minutes)")) +
#   scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +  # scale_x_continuous(trans = 'log2', breaks=c(2, 4, 8, 15, 30, 60)) +
#   # scale_y_log10() +
#   # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
#   myggplottheme; p6

p7 <- ggplot(time_values, aes(y=100*(gsw_VPD1.01-gsw_VPD1), x=Tau)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x="right", label.y=0.9, size=12, parse=T) + # method="spearman",
  # labs(y=expression(g["s,12"]-g["s,11"]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("∂"*g["S"]*" / ∂D (mol/m"^2*"/s/kPa)"), x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p7

ggarrange(p4+rremove("xlab"), p5+rremove("xlab"),
          p7, ncol=2, nrow=2, legend="none",
          labels="AUTO", label.x=0.07, label.y=0.95,
          font.label=list(size=30), align="hv")
ggsave("Output//scaled_stomatal_responsiveness-tau-ms.png", width=500, height=400, units="mm")

#### (absolute scale) stomatal diurnal dynamics ----

# p4 <- ggplot(time_values, aes(y=(absgsw_12-absgsw_11), x=Tau)) +
#   geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
#   geom_smooth(method = "lm", colour="black", se=F) +
#   ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
#                         label.x="right", label.y=0.1, size=12) + # method="spearman",
#   # labs(y=expression(g["s,12"]-g["s,11"]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
#   labs(y=expression(Delta*g["S,12-11"]~"(mol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
#   scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
#   # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
#   myggplottheme; p4
# 
# p5 <- ggplot(time_values, aes(y=absAarea_12-absAarea_11, x=Tau)) +
#   geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
#   geom_smooth(method = "lm", colour="black", se=F) +
#   ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
#                         label.x="right", label.y=0.1, size=12) + # method="spearman",
#   labs(y=expression(Delta*A["12-11"]~"("*mu*"mol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
#   scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
#   # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
#   myggplottheme; p5
# 
# p6 <- ggplot(time_values, aes(y=(absWUE_12-absWUE_13), x=Tau)) +
#   geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
#   geom_smooth(method = "lm", colour="black", se=F) +
#   ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
#                         label.x="right", label.y=0.1, size=12) + # method="spearman",
#   # labs(y=expression(g[s12]-g[s11]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
#   labs(y=expression(Delta*WUE["12-13"]~"(ppm)"), x=expression(tau~"(minutes)")) +
#   scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +  # scale_x_continuous(trans = 'log2', breaks=c(2, 4, 8, 15, 30, 60)) +
#   # scale_y_log10() +
#   # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
#   myggplottheme; p6
# 
# ggarrange(p4+rremove("xlab"), p5+rremove("xlab"),
#           p6, ncol=1,
#           labels="AUTO", label.x=0.018, label.y=0.95,
#           font.label=list(size=30), align="hv")
# ggsave("Output//absolute_stomatal_responsiveness-tau.png", width=250, height=600, units="mm")

#### Fig. 6: iWUE vs. tau ----

p7 <- ggplot(time_values, aes(y=iWUE, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  # geom_errorbar(aes(ymin=iWUE-iWUE_se, ymax=iWUE+iWUE_se),
  #               width=.05, size=1.2) + # ,position=position_dodge(0.05)
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        label.x="left", label.y="top", size=12) + # method="spearman",
  # labs(y=expression(g[s12]-g[s11]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
  labs(y="iWUE (daily average; ppm)", x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(2, 4, 8, 15, 30, 60)) +
  scale_y_log10(breaks=c(60, 80, 100)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p7

p8 <- ggplot(dat, aes(y=iWUE, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        label.x="left", label.y="top", size=12) + # method="spearman",
  # labs(y=expression(g[s12]-g[s11]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
  labs(y="iWUE (morning average; ppm)", x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(1, 2, 4, 8, 15, 30, 60)) +
  scale_y_log10(breaks=c(60, 80, 100, 120, 140)) +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p8

p9 <- ggplot(dat, aes(y=iWUE_from_isotope, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        label.x="right", label.y="bottom", size=12) + # method="spearman",
  # labs(y=expression(g[s12]-g[s11]~"(μmol/m"^2*"/s)"), x=expression(tau~"(minutes)")) +
  labs(y="iWUE (isotope; ppm)", x=expression(tau~"(minutes)")) +
  scale_x_log10(breaks=c(1, 2, 4, 8, 15, 30, 60)) +
  scale_y_log10() +
  # directlabels::geom_dl(aes(label = Spp), method = list("last.points", cex = 0.8)) +
  myggplottheme; p9

ggarrange(p7+rremove("xlab"), p8, p9,
          ncol=2, nrow=2,
          labels="AUTO", label.x=0.07, label.y=0.95,
          font.label=list(size=30), align="hv")
ggsave("Output//iWUE-tau.png", width=500, height=400, units="mm")

#### Fig. S2: time curves comparing all species ----
p1 <- ggplot(timecurves_values, aes(y=A_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("A ("*mu*"mol/m"^2*"/s)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p1

p2 <- ggplot(timecurves_values, aes(y=A, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("A (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p2

p3 <- ggplot(timecurves_values, aes(y=gsw_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(g[S]*" (mol/m"^2*"/s)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p3

p4 <- ggplot(timecurves_values, aes(y=gsw, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(g[S]*" (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p4

p5 <- ggplot(timecurves_values, aes(y=1000*E_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("E (mmol/m"^2*"/s)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p5

p6 <- ggplot(timecurves_values, aes(y=E, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("E (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p6

p7 <- ggplot(timecurves_values, aes(y=VPDleaf_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("*D* (kPa)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  scale_y_continuous(breaks=c(0.5,1,1.5,2,2.5)) +
  myggplottheme +
  theme(axis.title.y = ggtext::element_markdown()); p7

p8 <- ggplot(timecurves_values, aes(y=VPDleaf, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("*D* (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(axis.title.y = ggtext::element_markdown()); p8

p9 <- ggplot(timecurves_values, aes(y=WUE_abs/1000, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("WUE (ppm)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p9

p10 <- ggplot(timecurves_values, aes(y=WUE, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("WUE (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p10

p11 <- ggplot(timecurves_values, aes(y=Psi_xylem_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Psi[X]*" (MPa)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p11

p12 <- ggplot(timecurves_values, aes(y=Psi_xylem, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Psi[X]*" (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p12

p13 <- ggplot(timecurves_values, aes(y=Psi_leaf_abs, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Psi[L]*" (MPa)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p13

p14 <- ggplot(timecurves_values, aes(y=Psi_leaf, x=time)) +
  geom_line(size=2, aes(colour=Spp)) + 
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Psi[L]*" (scaled)"), x="time (o'clock)") + # μ
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme; p14

ggarrange(p1+rremove("xlab"), p3+rremove("xlab"),
          p5+rremove("xlab"), p7+rremove("xlab"),
          p9+rremove("xlab"), p11+rremove("xlab"),
          p13+rremove("xlab"), p2+rremove("xlab"),
          p4+rremove("xlab"), p6+rremove("xlab"),
          p8, p10+rremove("xlab"),
          p12+rremove("xlab"), p14+rremove("xlab"),
          common.legend=T,
          ncol=2, nrow=7, labels="AUTO", label.x=0.018, label.y=0.95,
          font.label=list(size=30), align="hv")
ggsave("Output//diurnal_traits-7spp-ms.png", width=500, height=1400, units="mm", limitsize=F)

#### Fig. S3: deltaPsi (X-L) vs. Cs ----
p1 <- ggplot(timecurves_values, aes(y=Psi_diff, x=time)) +
  geom_line(size=1.5, aes(colour=Spp)) +
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression(Delta*Psi["X-L"]~"(MPa)"), x="time of day (hour)") + # μ
  coord_cartesian(xlim = c(6, 18)) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  myggplottheme +
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p1

p2 <- ggplot(time_values, aes(y=Psi_diff_mean, x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(p.value.label))),
                        formula=y~x, parse=T, label.x="right", label.y=0.1, size=12) + # method="spearman",
  labs(y=expression("mean "*Delta*Psi["X-L"]~"(MPa)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  scale_fill_manual(values=cbp, aesthetics="colour") +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  myggplottheme; p2

p3 <- ggplot(time_values, aes(y=sqrt(Psi_diff_var), x=Cs)) +
  geom_point(size=7, stroke=2, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\", \"*")),
                        formula=y~x, parse=T, label.x="right", label.y=0.1, size=12) + # method="spearman",
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y=expression("sd "*Delta*Psi["X-L"]~"(MPa)"), x=expression(C[S]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim = c(50, 800)) +
  scale_x_log10(breaks=c(50, 100, 200, 400, 800)) +
  myggplottheme; p3

ggarrange(p1, p2, p3,
          common.legend=T, legend="none",
          ncol=2, nrow=2, labels="AUTO",
          label.x=0.018, label.y=0.95,
          font.label=list(size=30), align="hv")

ggsave("Output//DeltaPsi-Cs-ms.png", width=500, height=400, units="mm", limitsize=F)

#### Fig. S4: iWUE vs. time of day ----
p9 <- ggplot(subset(diel_12hr2, iWUE<200), aes(y=iWUE, x=hhmmss_int)) +
  geom_point(size=4, aes(colour=Spp)) + # change pch for solid/filled circles
  geom_smooth(method="lm", size=2, aes(colour=Spp), linetype="dashed", se=F) +
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y="iWUE", x="time of day (h)") +
  guides(color=guide_legend(override.aes=list(size=5))) +
  scale_x_continuous(breaks=c(6,9,12,15,18)) +
  coord_cartesian(xlim = c(6,18)) +
  myggplottheme +
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.title=element_text(size=30),
        legend.text=element_text(size=27)); p9

p10 <- ggplot(diel_12hr2, aes(y=iWUE, x=Spp)) +
  geom_boxplot(aes(colour=Spp)) + # change pch for solid/filled circles
  geom_jitter(aes(colour=Spp), size=4) +
  scale_fill_manual(values=cbp, aesthetics="colour") +
  labs(y="iWUE", x="Species") +
  myggplottheme; p10

ggarrange(p9, p10+rremove("ylab"),
          ncol=2, nrow=1, labels="AUTO",
          label.x=0.01, label.y=0.95, 
          common.legend=T, legend="right",
          font.label=list(size=30), align="hv")

ggsave("Output//iWUE-time-ms.png", width=550, height=250, units="mm")
#### Fig. S6: tau and components ----
component_Cs <- lm(log(Cs)~log(Tau), data=dat); summary(component_Cs)
p9 <- ggplot(dat, aes(y=Cs, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F, size=1.5) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(eq.label),sep="")),
                        formula=y~x,label.x=0.05, label.y=0.95, size=12, parse=T) + # method="spearman",
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x=0.05, label.y=0.85, size=12, parse=T) + # method="spearman",
  labs(x=expression(tau~"(minutes)"), y=expression(C["S"]~"(kg/m"^3*"/MPa)")) +
  coord_cartesian(xlim=c(0.5, 100)) +
  scale_y_log10(breaks=c(10,20,30,40,50,60)) +
  scale_x_log10(breaks=c(1,10,100)) +
  myggplottheme; p9

component_Ks <- lm(log(Ks_max)~log(Tau), data=dat); summary(component_Ks)
p10 <- ggplot(dat, aes(y=Ks_max, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F, size=1.5, linetype="dotted") +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(eq.label),sep="")),
                        formula=y~x,label.x=0.95, label.y=0.15, size=12, parse=T) + # method="spearman",
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x=0.95, label.y=0.05, size=12, parse=T) + # method="spearman",
  labs(x=expression(tau~"(minutes)"), y=expression(K["S"]~"(kg/m/s/MPa)")) +
  coord_cartesian(xlim = c(0.5, 100)) +
  scale_y_log10() +
  scale_x_log10(breaks=c(1,10,100)) +
  myggplottheme; p10

component_H <- lm(log(Height)~log(Tau), data=dat); summary(component_H)
p11 <- ggplot(dat, aes(y=Height, x=Tau)) +
  geom_point(size=7, stroke=2) + # change pch for solid/filled circles
  geom_smooth(method = "lm", colour="black", se=F, size=1.5) +
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(eq.label),sep="")),
                        formula=y~x,label.x=0.05, label.y=0.95, size=12, parse=T) + # method="spearman",
  ggpmisc::stat_poly_eq(aes(label=paste(after_stat(rr.label),after_stat(p.value.label),sep="*\",\"~")),
                        formula=y~x,label.x=0.05, label.y=0.85, size=12, parse=T) + # method="spearman",
  labs(x=expression(tau~"(minutes)"), y="H (m)") +
  coord_cartesian(xlim=c(0.5, 100)) +
  scale_y_log10(breaks=c(1.5,3,5,10)) +
  scale_x_log10(breaks=c(1,10,100)) +
  myggplottheme; p11

ggarrange(p9, p10, p11,
          ncol=3, nrow=1, labels="AUTO",
          label.x=0.018, label.y=0.95,
          font.label=list(size=30), align="hv")

ggsave("Output//tau-components-ms.png", width=750, height=200, units="mm")