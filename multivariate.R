# Multivariate analysis

# data set
setwd("F:/Synthese/R scripts")
df <- read.csv("ES.csv", header = T)

# Packages
library(dplyr)
library(ggplot2)
library(egg)
library(MuMIn)
library(cowplot)
library(car)
library(plyr)
library(usdm)
library(corrplot)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)
library(ggdist)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(egg)
library(ggpubr)
library(sjstats)
library(performance)

## BRMS design variables
# seperate functions with scale
bf_tax <- bf(scale(Tax.rich) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_spe <- bf(scale(Prop.spe) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_phy <- bf(scale(Phy.div) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_fun <- bf(scale(Fun.div) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_car <- bf(scale(ctot) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_bio <- bf(scale(agb.stem) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_nec <- bf(scale(nectar) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_use <- bf(scale(usable) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_mic <- bf(scale(summer.air) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_nut <- bf(scale(foliar.CN) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_tre <- bf(scale(tree.cover) ~  scale(lat) + type + scale(plot.num) + (1|region/transect))

fit1 <- brm(
  bf_tax + bf_spe + bf_phy + bf_fun + bf_car + bf_bio + bf_nec + bf_use + bf_mic + bf_nut + bf_tre + set_rescor(FALSE),
  data = df,  iter = 2000, warmup = 1000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.design.scale",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

# seperate functions without scale
bf_tax <- bf(Tax.rich ~  scale(plot.num) + (1|region/transect), family = poisson(link = "log"))
bf_spe <- bf(Prop.spe ~  scale(plot.num) + (1|region/transect), zi ~  scale(lat) + type + scale(plot.num) + (1|region/transect), family = zero_inflated_beta(link="logit",link_phi="log",link_zi = "logit"))
bf_phy <- bf(Phy.div ~  scale(plot.num) + (1|region/transect))
bf_fun <- bf(Fun.div ~  scale(plot.num) + (1|region/transect))
bf_car <- bf(ctot ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_bio <- bf(agb.stem ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_nec <- bf(nectar ~  scale(lat) + type + scale(plot.num) + (1|region/transect), hu ~  scale(lat) + type + scale(plot.num) + (1|region/transect), family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"))
bf_use <- bf(usable ~  scale(lat) + type + scale(plot.num) + (1|region/transect), zi ~  scale(lat) + type + scale(plot.num) + (1|region/transect), family = zero_inflated_beta(link="logit",link_phi="log",link_zi = "logit"))
bf_mic <- bf(summer.air ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_nut <- bf(foliar.CN ~  scale(lat) + type + scale(plot.num) + (1|region/transect))
bf_tre <- bf(tree.cover ~  scale(lat) + type + scale(plot.num) + (1|region/transect), zi ~  scale(lat) + type + scale(plot.num) + (1|region/transect), family = zero_inflated_beta(link="logit",link_phi="log",link_zi = "logit"))

fit2 <- brm(
  bf_tax + bf_spe + bf_phy + bf_fun + bf_car + bf_bio + bf_nec + bf_use + bf_mic + bf_nut + bf_tre + set_rescor(FALSE),
  data = df,  iter = 2000, warmup = 1000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.design2",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

fit2 <- readRDS("mult.design2.rds")
bayes_R2(fit2, re_form = NA)

summary(fit1)
summary(fit2)
pp_check(fit1, resp = "scaleTaxrich", plotfun = "scatter_avg")
r2_bayes(fit1)
r2_bayes(fit2)

## Model 1 Distance to forest edge ###########################################################################################
# Z-score transformation
bf_new <- bf(mvbind(scale(Tax.rich),scale(Prop.spe),scale(Fun.div),scale(Phy.div), scale(ctot),scale(agb.stem),
             scale(nectar),scale(usable),scale(summer.air),scale(foliar.CN),scale(tree.cover)) 
            ~ scale(plot.num) + (1|p|region/transect)) + set_rescor(TRUE)

fit3 <- brm(
  bf_new,
  data = df,  iter = 4000, warmup = 2000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.design.scale2",
  control = list(adapt_delta = .99, max_treedepth = 12)
)


fit3 <- readRDS("mult.design.scale2.rds")
table3 <- data.frame(summarise_draws(fit3))
write.csv(table3,"design_table.csv")
r2_bayes(fit3)

# Log- and Z-score transformation
bf_new <- bf(mvbind(scale(Tax.rich),scale(Prop.spe),scale(log(Fun.div+1)),scale(log(Phy.div+1)),scale(ctot),scale(log(agb.stem+1)),
                    scale(log(nectar+1)),scale(usable),scale(summer.air),scale(foliar.CN),scale(log(tree.cover+1))) 
             ~ scale(plot.num) + (1|region/transect)) + set_rescor(TRUE)

fit3 <- brm(
  bf_new,
  data = df,  iter = 4000, warmup = 2000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.design.scale4",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

fit3 <- readRDS("mult.design.scale4.rds")
table3 <- data.frame(summarise_draws(fit3))
write.csv(table3,"design_table.csv")
bayes_R2(fit3)

## Posterior predictive checks
p1 <- brms::pp_check(fit3, resp = "scaleTaxrich", ndraws = 100) + ggtitle("Taxonomic richness")
p2 <- brms::pp_check(fit3, resp = "scalePropspe", ndraws = 100) + ggtitle("% forest specialists")
p3 <- brms::pp_check(fit3, resp = "scalelogPhydiv1", ndraws = 100) + ggtitle("Phylogenetic diversity")
p4 <- brms::pp_check(fit3, resp = "scalelogFundiv1", ndraws = 100) + ggtitle("Functional diversity")
p5 <- brms::pp_check(fit3, resp = "scalectot", ndraws = 100) + ggtitle("Soil C storage")
p6 <- brms::pp_check(fit3, resp = "scalelognectar1", ndraws = 100) + ggtitle("Nectar production")
p7 <- brms::pp_check(fit3, resp = "scalesummerair", ndraws = 100) + ggtitle("Summer offset")
p8 <- brms::pp_check(fit3, resp = "scalefoliarCN", ndraws = 100) + ggtitle("Foliar C:N")
p9 <- brms::pp_check(fit3, resp = "scalelogagbstem1", ndraws = 100) + ggtitle("Stemwood biomass")
p10 <- brms::pp_check(fit3, resp = "scaleusable", ndraws = 100) + ggtitle("Usable plants")
p11 <- brms::pp_check(fit3, resp = "scalelogtreecover1", ndraws = 100) + ggtitle("Tree seedling cover")

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nrow = 4, ncol=3, align='hv')
png(file="ppchecks_m1.png",width=20,height=25,units="cm",res=300, pointsize=12)
p
dev.off()

# Z-score transformation and skewnormal
bf_tax <- bf(scale(Tax.rich) ~  plot.num + (1|region/transect))
bf_spe <- bf(scale(Prop.spe) ~  plot.num + (1|region/transect))
bf_phy <- bf(scale(Phy.div) ~  plot.num + (1|region/transect)) + lf(sigma ~ 0 + plot.num) + skew_normal()
bf_fun <- bf(scale(Fun.div) ~  plot.num + (1|region/transect)) + lf(sigma ~ 0 + plot.num) + skew_normal()
bf_car <- bf(scale(ctot) ~  plot.num + (1|region/transect))
bf_bio <- bf(scale(agb.stem) ~  plot.num + (1|region/transect)) + lf(sigma ~ 0 + plot.num) + skew_normal()
bf_nec <- bf(scale(nectar) ~  plot.num + (1|region/transect)) + lf(sigma ~ 0 + plot.num) + skew_normal()
bf_use <- bf(scale(usable) ~  plot.num + (1|region/transect))
bf_mic <- bf(scale(summer.air) ~  plot.num + (1|region/transect))
bf_nut <- bf(scale(foliar.CN) ~  plot.num + (1|region/transect))
bf_tre <- bf(scale(tree.cover) ~  plot.num + (1|region/transect)) + lf(sigma ~ 0 + plot.num) + skew_normal()

fit3 <- brm(
  bf_tax + bf_spe + bf_phy + bf_fun + bf_car + bf_bio + bf_nec + bf_use + bf_mic + bf_nut + bf_tre + set_rescor(F),
  data = df,  iter = 4000, warmup = 2000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.design.skew",
  control = list(adapt_delta = .99, max_treedepth = 12)
)


## Posterior predictive checks
p1 <- brms::pp_check(fit3, resp = "scaleTaxrich", ndraws = 100) + ggtitle("Taxonomic richness")
p2 <- brms::pp_check(fit3, resp = "scalePropspe", ndraws = 100) + ggtitle("% forest specialists")
p3 <- brms::pp_check(fit3, resp = "scalePhydiv", ndraws = 100) + ggtitle("Phylogenetic diversity")
p4 <- brms::pp_check(fit3, resp = "scaleFundiv", ndraws = 100) + ggtitle("Functional diversity")
p5 <- brms::pp_check(fit3, resp = "scalectot", ndraws = 100) + ggtitle("Soil C storage")
p6 <- brms::pp_check(fit3, resp = "scalenectar", ndraws = 100) + ggtitle("Nectar production")
p7 <- brms::pp_check(fit3, resp = "scalesummerair", ndraws = 100) + ggtitle("Summer offset")
p8 <- brms::pp_check(fit3, resp = "scalefoliarCN", ndraws = 100) + ggtitle("Foliar C:N")
p9 <- brms::pp_check(fit3, resp = "scaleagbstem", ndraws = 100) + ggtitle("Stemwood biomass")
p10 <- brms::pp_check(fit3, resp = "scaleusable", ndraws = 100) + ggtitle("Usable plants")
p11 <- brms::pp_check(fit3, resp = "scaletreecover", ndraws = 100) + ggtitle("Tree seedling cover")

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nrow = 4, ncol=3, align='hv')
png(file="ppchecks_mskew.png",width=20,height=25,units="cm",res=300, pointsize=12)
p
dev.off()

## Fig. 1B posterior means plots
c_eff <- conditional_effects(fit2, "plot.num", resp = "Taxrich")
df_eff <- as.data.frame(c_eff$Taxrich.Taxrich_plot.num)
p1<-ggplot(df,aes(x=plot.num,y=Tax.rich)) +
  geom_point(size = 1, alpha = 0.3, color ="#99CC00") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#99CC00") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#99CC00") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Taxonomic richness")

c_eff <- conditional_effects(fit2, "plot.num", resp = "Propspe")
df_eff <- as.data.frame(c_eff$Propspe.Propspe_plot.num)
p2<-ggplot(df,aes(x=plot.num,y=100*Prop.spe)) +
  geom_point(size = 1, alpha = 0.3, color ="#99CC00") +
  geom_line(data=df_eff,aes(x=plot.num,y=100*estimate__), size=1, color="#99CC00") + 
  geom_ribbon(data=df_eff,aes(ymin=100*lower__, ymax=100*upper__), alpha=0.2, fill = "#99CC00") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Forest specialists (%)")

c_eff <- conditional_effects(fit2, "plot.num", resp = "Phydiv")
df_eff <- as.data.frame(c_eff$Phydiv.Phydiv_plot.num)
p3<-ggplot(df,aes(x=plot.num,y=Phy.div)) +
  geom_point(size = 1, alpha = 0.3, color ="#99CC00") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#99CC00") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#99CC00") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Phylogenetic diversity")

c_eff <- conditional_effects(fit2, "plot.num", resp = "nectar")
df_eff <- as.data.frame(c_eff$nectar.nectar_plot.num)
p4<-ggplot(df,aes(x=plot.num,y=nectar)) +
  geom_point(size = 1, alpha = 0.3, color ="#0099FF") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#0099FF") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#0099FF") +
  scale_y_continuous(limits = c(0,25)) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y=bquote("Nectar ("*g~m^-2~year^-1*")"))

c_eff <- conditional_effects(fit2, "plot.num", resp = "summerair")
df_eff <- as.data.frame(c_eff$summerair.summerair_plot.num)
p5<-ggplot(df,aes(x=plot.num,y=summer.air)) +
  geom_point(size = 1, alpha = 0.3, color ="#0099FF") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#0099FF") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#0099FF") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Summer offset (°C)")

c_eff <- conditional_effects(fit2, "plot.num", resp = "foliarCN")
df_eff <- as.data.frame(c_eff$foliarCN.foliarCN_plot.num)
p6<-ggplot(df,aes(x=plot.num,y=foliar.CN)) +
  geom_point(size = 1, alpha = 0.3, color ="#0099FF") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#0099FF") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#0099FF") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Foliar C:N")

c_eff <- conditional_effects(fit2, "plot.num", resp = "agbstem")
df_eff <- as.data.frame(c_eff$agbstem.agbstem_plot.num)
p7<-ggplot(df,aes(x=plot.num,y=agb.stem)) +
  geom_point(size = 1, alpha = 0.3, color ="#FF9900") +
  geom_line(data=df_eff,aes(x=plot.num,y=estimate__), size=1, color="#FF9900") + 
  geom_ribbon(data=df_eff,aes(ymin=lower__, ymax=upper__), alpha=0.2, fill = "#FF9900") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_text(size=7), axis.text.x = element_text(size=6), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y=bquote("Stem biomass ("*Mg~ha^-1*")"))

c_eff <- conditional_effects(fit2, "plot.num", resp = "treecover")
df_eff <- as.data.frame(c_eff$treecover.treecover_plot.num)
p8<-ggplot(df,aes(x=plot.num,y=100*tree.cover)) +
  geom_point(size = 1, alpha = 0.3, color ="#FF9900") +
  geom_line(data=df_eff,aes(x=plot.num,y=100*estimate__), size=1, color="#FF9900") + 
  geom_ribbon(data=df_eff,aes(ymin=100*lower__, ymax=100*upper__), alpha=0.2, fill = "#FF9900") +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_text(size=7), axis.text.y = element_text(size=6), axis.title.x = element_text(size=7), axis.text.x = element_text(size=6), panel.grid.major =  element_blank(), panel.grid.minor = element_blank()) + labs(x="Distance to edge (m)",y="Tree seedling cover (%)")

png(file="fig1b.png",width=9,height=12,units="cm",res=300, pointsize=12)
p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4, ncol=2, align='v')
p
dev.off()

## Fig. 1A plot draws
png(file="design_model.png",width=9,height=12,units="cm",res=300, pointsize=12)

p1<-fit3 %>%
  gather_draws(c(b_scaleTaxrich_scaleplot.num, b_scalePropspe_scaleplot.num, b_scalelogPhydiv1_scaleplot.num, b_scalelogFundiv1_scaleplot.num, b_scalefoliarCN_scaleplot.num, b_scalectot_scaleplot.num, b_scalelognectar1_scaleplot.num, b_scalesummerair_scaleplot.num, b_scalelogagbstem1_scaleplot.num, b_scaleusable_scaleplot.num, b_scalelogtreecover1_scaleplot.num)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scaleplot.num" | .variable=="b_scalePropspe_scaleplot.num" | .variable=="b_scalelogPhydiv1_scaleplot.num" | .variable=="b_scalelogFundiv1_scaleplot.num", "G1", ifelse(.variable=="b_scalectot_scaleplot.num" | .variable=="b_scalelognectar1_scaleplot.num" | .variable=="b_scalesummerair_scaleplot.num" | .variable=="b_scalefoliarCN_scaleplot.num", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scaleplot.num", "b_scalePropspe_scaleplot.num", "b_scalelogPhydiv1_scaleplot.num","b_scalelogFundiv1_scaleplot.num", "b_scalectot_scaleplot.num", "b_scalelognectar1_scaleplot.num", "b_scalesummerair_scaleplot.num", "b_scalefoliarCN_scaleplot.num", "b_scalelogagbstem1_scaleplot.num", "b_scaleusable_scaleplot.num", "b_scalelogtreecover1_scaleplot.num")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Forest edge effect") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p1

dev.off()

## Model 2 BRMS environmental variables ###########################################################################################
bf_tax <- bf(scale(Tax.rich) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_spe <- bf(scale(Prop.spe) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_phy <- bf(scale(Phy.div) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_fun <- bf(scale(Fun.div) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep))+ (1|region/transect))
bf_car <- bf(scale(ctot) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_bio <- bf(scale(agb.stem) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_nec <- bf(scale(nectar) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_use <- bf(scale(usable) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_mic <- bf(scale(summer.air) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_nut <- bf(scale(foliar.CN) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_tre <- bf(scale(tree.cover) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))

fit4 <- brm(
  bf_tax + bf_spe + bf_phy + bf_fun + bf_car + bf_bio + bf_nec + bf_use + bf_mic + bf_nut + bf_tre + set_rescor(TRUE),
  data = df,  iter = 4000, warmup = 2000, chains = 2,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.env.scale",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

# Log- and Z-score transformation
bf_new <- bf(mvbind(scale(Tax.rich),scale(Prop.spe),scale(log(Fun.div+1)),scale(log(Phy.div+1)),scale(ctot),scale(log(agb.stem+1)),
                    scale(log(nectar+1)),scale(usable),scale(summer.air),scale(foliar.CN),scale(log(tree.cover+1))) 
             ~ scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect)) + set_rescor(TRUE)

fit4 <- brm(
  bf_new,
  data = df,  iter = 4000, warmup = 2000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.env.scale",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

fit4 <- readRDS("mult.env.scale.rds")
table4 <- data.frame(summarise_draws(fit4))
write.csv(table4,"env_table.csv")
bayes_R2(fit4)

## Posterior predictive checks
p1 <- brms::pp_check(fit4, resp = "scaleTaxrich", ndraws = 100) + ggtitle("Taxonomic richness")
p2 <- brms::pp_check(fit4, resp = "scalePropspe", ndraws = 100) + ggtitle("% forest specialists")
p3 <- brms::pp_check(fit4, resp = "scalelogPhydiv1", ndraws = 100) + ggtitle("Phylogenetic diversity")
p4 <- brms::pp_check(fit4, resp = "scalelogFundiv1", ndraws = 100) + ggtitle("Functional diversity")
p5 <- brms::pp_check(fit4, resp = "scalectot", ndraws = 100) + ggtitle("Soil C storage")
p6 <- brms::pp_check(fit4, resp = "scalelognectar1", ndraws = 100) + ggtitle("Nectar production")
p7 <- brms::pp_check(fit4, resp = "scalesummerair", ndraws = 100) + ggtitle("Summer offset")
p8 <- brms::pp_check(fit4, resp = "scalefoliarCN", ndraws = 100) + ggtitle("Foliar C:N")
p9 <- brms::pp_check(fit4, resp = "scalelogagbstem1", ndraws = 100) + ggtitle("Stemwood biomass")
p10 <- brms::pp_check(fit4, resp = "scaleusable", ndraws = 100) + ggtitle("Usable plants")
p11 <- brms::pp_check(fit4, resp = "scalelogtreecover1", ndraws = 100) + ggtitle("Tree seedling cover")

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nrow = 4, ncol=3, align='hv')
png(file="ppchecks_m2.png",width=20,height=25,units="cm",res=300, pointsize=12)
p
dev.off()

# Z-score transformation and skewnormal
bf_tax <- bf(scale(Tax.rich) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_spe <- bf(scale(Prop.spe) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_phy <- bf(scale(Phy.div) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))  + skew_normal()
bf_fun <- bf(scale(Fun.div) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))  + skew_normal()
bf_car <- bf(scale(ctot) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_bio <- bf(scale(agb.stem) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))  + skew_normal()
bf_nec <- bf(scale(nectar) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))  + skew_normal()
bf_use <- bf(scale(usable) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_mic <- bf(scale(summer.air) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_nut <- bf(scale(foliar.CN) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))
bf_tre <- bf(scale(tree.cover) ~  scale(sand) + scale(ph) + scale(log(weight_oh_of_ol)) + scale(pai) + scale(sca) + scale(BIO5) +scale(forcov_500) + scale(spei21_May2018) + scale(log(ndep)) + (1|region/transect))  + skew_normal()

fit4 <- brm(
  bf_tax + bf_spe + bf_phy + bf_fun + bf_car + bf_bio + bf_nec + bf_use + bf_mic + bf_nut + bf_tre + set_rescor(F),
  data = df,  iter = 4000, warmup = 2000, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.env.skew",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

bayes_R2(fit4)

## Posterior predictive checks
p1 <- brms::pp_check(fit4, resp = "scaleTaxrich", ndraws = 100) + ggtitle("Taxonomic richness")
p2 <- brms::pp_check(fit4, resp = "scalePropspe", ndraws = 100) + ggtitle("% forest specialists")
p3 <- brms::pp_check(fit4, resp = "scalePhydiv", ndraws = 100) + ggtitle("Phylogenetic diversity")
p4 <- brms::pp_check(fit4, resp = "scaleFundiv", ndraws = 100) + ggtitle("Functional diversity")
p5 <- brms::pp_check(fit4, resp = "scalectot", ndraws = 100) + ggtitle("Soil C storage")
p6 <- brms::pp_check(fit4, resp = "scalenectar", ndraws = 100) + ggtitle("Nectar production")
p7 <- brms::pp_check(fit4, resp = "scalesummerair", ndraws = 100) + ggtitle("Summer offset")
p8 <- brms::pp_check(fit4, resp = "scalefoliarCN", ndraws = 100) + ggtitle("Foliar C:N")
p9 <- brms::pp_check(fit4, resp = "scaleagbstem", ndraws = 100) + ggtitle("Stemwood biomass")
p10 <- brms::pp_check(fit4, resp = "scaleusable", ndraws = 100) + ggtitle("Usable plants")
p11 <- brms::pp_check(fit4, resp = "scaletreecover", ndraws = 100) + ggtitle("Tree seedling cover")

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, nrow = 4, ncol=3, align='hv')
png(file="ppchecks_m2skew.png",width=20,height=25,units="cm",res=300, pointsize=12)
p
dev.off()

#pdf(file="F:/Synthese/R scripts/multi_env.pdf")

png(file="multi_sand.png",width=10,height=10,units="cm",res=300, pointsize=12)

p1<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalesand, b_scalePropspe_scalesand, b_scalelogPhydiv1_scalesand, b_scalelogFundiv1_scalesand, b_scalefoliarCN_scalesand, b_scalectot_scalesand, b_scalelognectar1_scalesand, b_scalesummerair_scalesand, b_scalelogagbstem1_scalesand, b_scaleusable_scalesand, b_scalelogtreecover1_scalesand)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalesand" | .variable=="b_scalePropspe_scalesand" | .variable=="b_scalelogPhydiv1_scalesand" | .variable=="b_scalelogFundiv1_scalesand", "G1", ifelse(.variable=="b_scalectot_scalesand" | .variable=="b_scalelognectar1_scalesand" | .variable=="b_scalesummerair_scalesand" | .variable=="b_scalefoliarCN_scalesand", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalesand", "b_scalePropspe_scalesand", "b_scalelogPhydiv1_scalesand","b_scalelogFundiv1_scalesand", "b_scalectot_scalesand", "b_scalelognectar1_scalesand", "b_scalesummerair_scalesand", "b_scalefoliarCN_scalesand", "b_scalelogagbstem1_scalesand", "b_scaleusable_scalesand", "b_scalelogtreecover1_scalesand")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("% sand") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p1

dev.off()

png(file="multi_ph.png",width=10,height=10,units="cm",res=300, pointsize=12)

p2<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scaleph, b_scalePropspe_scaleph, b_scalelogPhydiv1_scaleph, b_scalelogFundiv1_scaleph, b_scalefoliarCN_scaleph, b_scalectot_scaleph, b_scalelognectar1_scaleph, b_scalesummerair_scaleph, b_scalelogagbstem1_scaleph, b_scaleusable_scaleph, b_scalelogtreecover1_scaleph)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scaleph" | .variable=="b_scalePropspe_scaleph" | .variable=="b_scalelogPhydiv1_scaleph" | .variable=="b_scalelogFundiv1_scaleph", "G1", ifelse(.variable=="b_scalectot_scaleph" | .variable=="b_scalelognectar1_scaleph" | .variable=="b_scalesummerair_scaleph" | .variable=="b_scalefoliarCN_scaleph", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scaleph", "b_scalePropspe_scaleph", "b_scalelogPhydiv1_scaleph","b_scalelogFundiv1_scaleph", "b_scalectot_scaleph", "b_scalelognectar1_scaleph", "b_scalesummerair_scaleph", "b_scalefoliarCN_scaleph", "b_scalelogagbstem1_scaleph", "b_scaleusable_scaleph", "b_scalelogtreecover1_scaleph")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("pH") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p2

dev.off()

png(file="multi_om.png",width=10,height=10,units="cm",res=300, pointsize=12)

p3<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalelogweight_oh_of_ol, b_scalePropspe_scalelogweight_oh_of_ol, b_scalelogPhydiv1_scalelogweight_oh_of_ol, b_scalelogFundiv1_scalelogweight_oh_of_ol, b_scalefoliarCN_scalelogweight_oh_of_ol, b_scalectot_scalelogweight_oh_of_ol, b_scalelognectar1_scalelogweight_oh_of_ol, b_scalesummerair_scalelogweight_oh_of_ol, b_scalelogagbstem1_scalelogweight_oh_of_ol, b_scaleusable_scalelogweight_oh_of_ol, b_scalelogtreecover1_scalelogweight_oh_of_ol)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalelogweight_oh_of_ol" | .variable=="b_scalePropspe_scalelogweight_oh_of_ol" | .variable=="b_scalelogPhydiv1_scalelogweight_oh_of_ol" | .variable=="b_scalelogFundiv1_scalelogweight_oh_of_ol", "G1", ifelse(.variable=="b_scalectot_scalelogweight_oh_of_ol" | .variable=="b_scalelognectar1_scalelogweight_oh_of_ol" | .variable=="b_scalesummerair_scalelogweight_oh_of_ol" | .variable=="b_scalefoliarCN_scalelogweight_oh_of_ol", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalelogweight_oh_of_ol", "b_scalePropspe_scalelogweight_oh_of_ol", "b_scalelogPhydiv1_scalelogweight_oh_of_ol","b_scalelogFundiv1_scalelogweight_oh_of_ol", "b_scalectot_scalelogweight_oh_of_ol", "b_scalelognectar1_scalelogweight_oh_of_ol", "b_scalesummerair_scalelogweight_oh_of_ol", "b_scalefoliarCN_scalelogweight_oh_of_ol", "b_scalelogagbstem1_scalelogweight_oh_of_ol", "b_scaleusable_scalelogweight_oh_of_ol", "b_scalelogtreecover1_scalelogweight_oh_of_ol")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Mass OS") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p3

dev.off()

png(file="multi_PAI.png",width=10,height=10,units="cm",res=300, pointsize=12)

p4<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalepai, b_scalePropspe_scalepai, b_scalelogPhydiv1_scalepai, b_scalelogFundiv1_scalepai, b_scalefoliarCN_scalepai, b_scalectot_scalepai, b_scalelognectar1_scalepai, b_scalesummerair_scalepai, b_scalelogagbstem1_scalepai, b_scaleusable_scalepai, b_scalelogtreecover1_scalepai)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalepai" | .variable=="b_scalePropspe_scalepai" | .variable=="b_scalelogPhydiv1_scalepai" | .variable=="b_scalelogFundiv1_scalepai", "G1", ifelse(.variable=="b_scalectot_scalepai" | .variable=="b_scalelognectar1_scalepai" | .variable=="b_scalesummerair_scalepai" | .variable=="b_scalefoliarCN_scalepai", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalepai", "b_scalePropspe_scalepai", "b_scalelogPhydiv1_scalepai","b_scalelogFundiv1_scalepai", "b_scalectot_scalepai", "b_scalelognectar1_scalepai", "b_scalesummerair_scalepai", "b_scalefoliarCN_scalepai", "b_scalelogagbstem1_scalepai", "b_scaleusable_scalepai", "b_scalelogtreecover1_scalepai")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Plant area index") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p4

dev.off()

png(file="multi_SCA.png",width=10,height=10,units="cm",res=300, pointsize=12)

p5<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalesca, b_scalePropspe_scalesca, b_scalelogPhydiv1_scalesca, b_scalelogFundiv1_scalesca, b_scalefoliarCN_scalesca, b_scalectot_scalesca, b_scalelognectar1_scalesca, b_scalesummerair_scalesca, b_scalelogagbstem1_scalesca, b_scaleusable_scalesca, b_scalelogtreecover1_scalesca)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalesca" | .variable=="b_scalePropspe_scalesca" | .variable=="b_scalelogPhydiv1_scalesca" | .variable=="b_scalelogFundiv1_scalesca", "G1", ifelse(.variable=="b_scalectot_scalesca" | .variable=="b_scalelognectar1_scalesca" | .variable=="b_scalesummerair_scalesca" | .variable=="b_scalefoliarCN_scalesca", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalesca", "b_scalePropspe_scalesca", "b_scalelogPhydiv1_scalesca","b_scalelogFundiv1_scalesca", "b_scalectot_scalesca", "b_scalelognectar1_scalesca", "b_scalesummerair_scalesca", "b_scalefoliarCN_scalesca", "b_scalelogagbstem1_scalesca", "b_scaleusable_scalesca", "b_scalelogtreecover1_scalesca")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Shade-casting ability") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p5

dev.off()

png(file="multi_BIO5.png",width=10,height=10,units="cm",res=300, pointsize=12)

p6<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scaleBIO5, b_scalePropspe_scaleBIO5, b_scalelogPhydiv1_scaleBIO5, b_scalelogFundiv1_scaleBIO5, b_scalefoliarCN_scaleBIO5, b_scalectot_scaleBIO5, b_scalelognectar1_scaleBIO5, b_scalesummerair_scaleBIO5, b_scalelogagbstem1_scaleBIO5, b_scaleusable_scaleBIO5, b_scalelogtreecover1_scaleBIO5)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scaleBIO5" | .variable=="b_scalePropspe_scaleBIO5" | .variable=="b_scalelogPhydiv1_scaleBIO5" | .variable=="b_scalelogFundiv1_scaleBIO5", "G1", ifelse(.variable=="b_scalectot_scaleBIO5" | .variable=="b_scalelognectar1_scaleBIO5" | .variable=="b_scalesummerair_scaleBIO5" | .variable=="b_scalefoliarCN_scaleBIO5", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scaleBIO5", "b_scalePropspe_scaleBIO5", "b_scalelogPhydiv1_scaleBIO5","b_scalelogFundiv1_scaleBIO5", "b_scalectot_scaleBIO5", "b_scalelognectar1_scaleBIO5", "b_scalesummerair_scaleBIO5", "b_scalefoliarCN_scaleBIO5", "b_scalelogagbstem1_scaleBIO5", "b_scaleusable_scaleBIO5", "b_scalelogtreecover1_scaleBIO5")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Tmax of warmest month") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p6

dev.off()

png(file="multi_forcov.png",width=10,height=10,units="cm",res=300, pointsize=12)

p7<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scaleforcov_500, b_scalePropspe_scaleforcov_500, b_scalelogPhydiv1_scaleforcov_500, b_scalelogFundiv1_scaleforcov_500, b_scalefoliarCN_scaleforcov_500, b_scalectot_scaleforcov_500, b_scalelognectar1_scaleforcov_500, b_scalesummerair_scaleforcov_500, b_scalelogagbstem1_scaleforcov_500, b_scaleusable_scaleforcov_500, b_scalelogtreecover1_scaleforcov_500)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scaleforcov_500" | .variable=="b_scalePropspe_scaleforcov_500" | .variable=="b_scalelogPhydiv1_scaleforcov_500" | .variable=="b_scalelogFundiv1_scaleforcov_500", "G1", ifelse(.variable=="b_scalectot_scaleforcov_500" | .variable=="b_scalelognectar1_scaleforcov_500" | .variable=="b_scalesummerair_scaleforcov_500" | .variable=="b_scalefoliarCN_scaleforcov_500", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scaleforcov_500", "b_scalePropspe_scaleforcov_500", "b_scalelogPhydiv1_scaleforcov_500","b_scalelogFundiv1_scaleforcov_500", "b_scalectot_scaleforcov_500", "b_scalelognectar1_scaleforcov_500", "b_scalesummerair_scaleforcov_500", "b_scalefoliarCN_scaleforcov_500", "b_scalelogagbstem1_scaleforcov_500", "b_scaleusable_scaleforcov_500", "b_scalelogtreecover1_scaleforcov_500")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Forest cover") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p7

dev.off()

png(file="multi_spei.png",width=10,height=10,units="cm",res=300, pointsize=12)

p8<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalespei21_May2018, b_scalePropspe_scalespei21_May2018, b_scalelogPhydiv1_scalespei21_May2018, b_scalelogFundiv1_scalespei21_May2018, b_scalefoliarCN_scalespei21_May2018, b_scalectot_scalespei21_May2018, b_scalelognectar1_scalespei21_May2018, b_scalesummerair_scalespei21_May2018, b_scalelogagbstem1_scalespei21_May2018, b_scaleusable_scalespei21_May2018, b_scalelogtreecover1_scalespei21_May2018)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalespei21_May2018" | .variable=="b_scalePropspe_scalespei21_May2018" | .variable=="b_scalelogPhydiv1_scalespei21_May2018" | .variable=="b_scalelogFundiv1_scalespei21_May2018", "G1", ifelse(.variable=="b_scalectot_scalespei21_May2018" | .variable=="b_scalelognectar1_scalespei21_May2018" | .variable=="b_scalesummerair_scalespei21_May2018" | .variable=="b_scalefoliarCN_scalespei21_May2018", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalespei21_May2018", "b_scalePropspe_scalespei21_May2018", "b_scalelogPhydiv1_scalespei21_May2018","b_scalelogFundiv1_scalespei21_May2018", "b_scalectot_scalespei21_May2018", "b_scalelognectar1_scalespei21_May2018", "b_scalesummerair_scalespei21_May2018", "b_scalefoliarCN_scalespei21_May2018", "b_scalelogagbstem1_scalespei21_May2018", "b_scaleusable_scalespei21_May2018", "b_scalelogtreecover1_scalespei21_May2018")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("SPEI") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p8

dev.off()

png(file="multi_ndep.png",width=10,height=10,units="cm",res=300, pointsize=12)

p9<-fit4 %>%
  gather_draws(c(b_scaleTaxrich_scalelogndep, b_scalePropspe_scalelogndep, b_scalelogPhydiv1_scalelogndep, b_scalelogFundiv1_scalelogndep, b_scalefoliarCN_scalelogndep, b_scalectot_scalelogndep, b_scalelognectar1_scalelogndep, b_scalesummerair_scalelogndep, b_scalelogagbstem1_scalelogndep, b_scaleusable_scalelogndep, b_scalelogtreecover1_scalelogndep)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaleTaxrich_scalelogndep" | .variable=="b_scalePropspe_scalelogndep" | .variable=="b_scalelogPhydiv1_scalelogndep" | .variable=="b_scalelogFundiv1_scalelogndep", "G1", ifelse(.variable=="b_scalectot_scalelogndep" | .variable=="b_scalelognectar1_scalelogndep" | .variable=="b_scalesummerair_scalelogndep" | .variable=="b_scalefoliarCN_scalelogndep", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaleTaxrich_scalelogndep", "b_scalePropspe_scalelogndep", "b_scalelogPhydiv1_scalelogndep","b_scalelogFundiv1_scalelogndep", "b_scalectot_scalelogndep", "b_scalelognectar1_scalelogndep", "b_scalesummerair_scalelogndep", "b_scalefoliarCN_scalelogndep", "b_scalelogagbstem1_scalelogndep", "b_scaleusable_scalelogndep", "b_scalelogtreecover1_scalelogndep")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("N deposition") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p9

dev.off()

## Edge-to-interior patterns in environmental variables
p1<-ggplot(df, aes(plot.num, sand)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="% sand") + geom_smooth(method="lm", colour = "black") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2<-ggplot(df, aes(plot.num, ph)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="pH") + geom_smooth(method="lm", colour = "black") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p3<-ggplot(df, aes(plot.num, weight_oh_of_ol)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="Mass OS (g)") + geom_smooth(method="lm", colour = "black") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p4<-ggplot(df, aes(plot.num, pai)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="PAI") + geom_smooth(method="lm", colour = "black") 
p5<-ggplot(df, aes(plot.num, sca)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="SCA") + geom_smooth(method="lm", colour = "black") 
p6<-ggplot(df, aes(plot.num, BIO5)) + geom_point() + theme_bw() + labs(x="Distance to forest edge", y="TmaxWm") + geom_smooth(method="lm", colour = "black") 

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, ncol=3, align='hv')
png(file="env_edge.png",width=16,height=8,units="cm",res=300, pointsize=12)
p
dev.off()

## Overstorey included
bf_new <- bf(mvbind(scale(total.ric),scale(Prop.spe),scale(log(Fun.div+1)),scale(log(total.phy+1)),scale(ctot),scale(log(agb.stem+1)),
                    scale(log(total.nectar+1)),scale(usable),scale(summer.air),scale(foliar.CN),scale(log(tree.cover+1))) 
             ~ scale(plot.num) + (1|region/transect)) + set_rescor(TRUE)

fit4 <- brm(
  bf_new,
  data = df,  iter = 1000, warmup = 500, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/mult.overstorey",
  control = list(adapt_delta = .99, max_treedepth = 12)
)

png(file="overstorey_model.png",width=9,height=12,units="cm",res=300, pointsize=12)

p1<-fit4 %>%
  gather_draws(c(b_scaletotalric_scaleplot.num, b_scalePropspe_scaleplot.num, b_scalelogtotalphy1_scaleplot.num, b_scalelogFundiv1_scaleplot.num, b_scalefoliarCN_scaleplot.num, b_scalectot_scaleplot.num, b_scalelogtotalnectar1_scaleplot.num, b_scalesummerair_scaleplot.num, b_scalelogagbstem1_scaleplot.num, b_scaleusable_scaleplot.num, b_scalelogtreecover1_scaleplot.num)) %>%
  #mean_qi(.width = c(.95, .8)) %>%
  mutate(group = ifelse(.variable=="b_scaletotalric_scaleplot.num" | .variable=="b_scalePropspe_scaleplot.num" | .variable=="b_scalelogtotalphy1_scaleplot.num" | .variable=="b_scalelogFundiv1_scaleplot.num", "G1", ifelse(.variable=="b_scalectot_scaleplot.num" | .variable=="b_scalelogtotalnectar1_scaleplot.num" | .variable=="b_scalesummerair_scaleplot.num" | .variable=="b_scalefoliarCN_scaleplot.num", "G2", "G3"))) %>%
  ggplot(aes(y = .variable, x = .value, color= group, fill= group)) +
  #ggplot(aes(y = .variable, x = .value, xmin = .lower, xmax = .upper, color = group, fill= group)) +
  stat_halfeye() +
  #geom_pointinterval() +
  scale_color_manual(values = c("#99CC00","#0099FF","#FF9900")) +
  scale_fill_manual(values = alpha(c("#99CC00","#0099FF","#FF9900"),0.5)) +
  geom_vline(xintercept=0, linetype ="dashed") +
  scale_y_discrete(limits=rev(c("b_scaletotalric_scaleplot.num", "b_scalePropspe_scaleplot.num", "b_scalelogtotalphy1_scaleplot.num","b_scalelogFundiv1_scaleplot.num", "b_scalectot_scaleplot.num", "b_scalelogtotalnectar1_scaleplot.num", "b_scalesummerair_scaleplot.num", "b_scalefoliarCN_scaleplot.num", "b_scalelogagbstem1_scaleplot.num", "b_scaleusable_scaleplot.num", "b_scalelogtreecover1_scaleplot.num")),labels = rev(c("Taxonomic richness", "% forest specialists","Phylogenetic diversity","Functional diversity", "Soil C storage", "Nectar production", "Summer offset", "Foliar C:N", "Stemwood biomass", "Usable plants", "Tree seedling cover"))) +
  labs(x = "Estimate", y = "") +
  ggtitle("Forest edge effect") +
  theme_bw() + theme(legend.position = "none", panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  axis.text.x = element_text(size=7), axis.title.x = element_blank())
p1

dev.off()

## biomass without plot 2
df <- df %>% filter(plot.num != 4.5)

fit <- brm(
  scale(log(agb.stem+1)) ~ scale(plot.num) + (1|region/transect),
  data = df,  iter = 1000, warmup = 500, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/biomass.new",
  control = list(adapt_delta = .99, max_treedepth = 12)
)
summarise_draws(fit)

bf_new <- bf(mvbind(scale(Tax.rich),scale(Prop.spe),scale(log(Fun.div+1)),scale(log(Phy.div+1)),scale(ctot),scale(log(agb.stem+1)),
                    scale(log(nectar+1)),scale(usable),scale(summer.air),scale(foliar.CN),scale(log(tree.cover+1))) 
             ~ scale(plot.num) + (1|region/transect)) + set_rescor(TRUE)

fit <- brm(
  bf_new,
  data = df,  iter = 1000, warmup = 500, chains = 4,
  seed = 190831,
  file = "F:/Synthese/R scripts/biomass.new2",
  control = list(adapt_delta = .99, max_treedepth = 12)
)
summarise_draws(fit)
