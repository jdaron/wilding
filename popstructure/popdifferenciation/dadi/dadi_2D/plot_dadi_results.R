library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gtable)
library(RColorBrewer)
library(strucchange)
library(plyr)
library(ggrepel)
library(geosphere)
library("gridExtra")
library(measurements)
library(ggpubr)
library(tidyr)
library(scales)
library(reshape2)

mat = read.table("~/bioInf/wilding/popstructure/popdifferenciation/dadi/dadi_2D/wilding_3_LPdom_LBVdom.dadi2D8Model.txt", h=F, dec = ".")
colnames(mat) = c("rep", "model", "optm", "stats", "value")
head(mat)

#####
## 1. Find the best model
#####

### 1.1 look at optimization
head(mat)
mat_ll_model = mat[which(mat$stats=="ll_model"),]
head(mat_ll_model)

mat_ll_model$optm = factor(mat_ll_model$optm, levels=c("anneal_hot", "anneal_cold", "BFGS"))

ggplot() +
  geom_boxplot(data=mat_ll_model, aes(x=optm, y=abs(value), color=model)) +
  facet_wrap(~model) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  scale_y_log10()

### 1.2 plot model by ordered likelihood
mat = mat[which(mat$optm=="BFGS"),]
mat.c = dcast(mat, rep+model~stats)
head(mat.c)

cols = brewer.pal(8, "Paired")
mat.c = mat.c[order(mat.c$optm_ll),]
mat.c$order = 1:dim(mat.c)[1]

p1 = ggplot() +
#  geom_bar(data=mat.c, aes(x=order, y=log10(abs(optm_ll)), color=model, fill=model), stat = "identity") +
  geom_bar(data=mat.c, aes(x=order, y=AIC, color=model, fill=model), stat = "identity") +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols)
p1 

ggplot() +
  geom_histogram(data=mat.c, aes(x=theta, color=model, fill=model)) +
  geom_vline(xintercept = 78000) +
  facet_wrap(~model)
  

pdf("~/bioInf/wilding/popstructure/popdifferenciation/dadi/dadi_2D/wilding_3_LPdom_LBVdom.dadi2D8Model.aic.barplot.pdf", width = conv_unit(150, "mm", "inch"), height = conv_unit(100, "mm", "inch"), useDingbats=FALSE) # export PDF 7x9
grid.arrange(p1, ncol = 1, nrow = 1, top = "")
dev.off()

ggplot() +
  geom_jitter(data=mat.c, aes(x=model, y=optm_ll, color=model, fill=model), size=2) +
  scale_shape_manual(values = c(21,22,24)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols)+
  xlab("Model") +
  ylab("Optimized_log likelihood") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  guides(color=guide_legend(ncol=2))

### 1.3 get best model and convert metrics
#step1: get best rep per model 
mat.c.best = ddply(mat.c, .(model), plyr::summarize, optm_ll=max(optm_ll))
head(mat.c.best)

# step2: select model with deltaAIC <= 10
mat.step2 = mat.c[which(mat.c$model %in% mat.c.best$model & mat.c$optm_ll %in% mat.c.best$optm_ll),]
mat.step2$deltaAIC = mat.step2$AIC-min(mat.step2$AIC)
mat.step2[which(mat.step2$deltaAIC<=10),]

mat.step2[,c("model", "optm_ll", "AIC", "deltaAIC")]

# step3: calc weightAIC
exp((-1*mat.step2$deltaAIC)/2)/sum(exp((-1*mat.step2$deltaAIC)/2))

# step4: choose best model and convert metrics
mat.best = mat.step2[which(mat.step2$optm_ll==max(mat.step2$optm_ll)),]
mat.best$model

# set up variable
mu = 3.5e-9 # mutation rate 
L = 62000000
#L = 30484401
g = 11 # generation per years

## Effective pop size
# nref, nu1, nu2: effective pop size of ancestral population, pop1 and pop2
# theta = 4*mu*L*Nref
# Nref = theta/(4*mu*L)
# Nref = theta/(mu*L)/4

nref = mat.best$theta/(4*mu*L)
nu1 = mat.best$nu1*nref
nu2 = mat.best$nu2*nref

## Time
# Tsc number of generation since the two pop have entered into a secondaru contact
Ts = 2*nref*mat.best$Ts*g
Tsc = 2*nref*mat.best$Tsc*g

# m12, m21 effective migration rates expressed in 2.Nrefm units per generation.
# m12 proportion of pop made migrants form pop2 into pop1
# m12, m21 = 2.Nrefm units per generation

r_m12 = mat.best$m12/(2*nref) #which means m12 individuals each generation are migrating from population 2 to population 1.
nb_m12 = (mat.best$nu1*r_m12)/2 # number of migrant: nu1*M12/2
r_m21 = mat.best$m21/(2*nref) 
nb_m21 = (mat.best$nu2*r_m21)/2

# b1 b2 growth

mat.best$b1
mat.best$b2

### bootstrap sd
boot_sd = c(9.31919869e-01, 2.79958976e-01, 1.14452344e-02, 5.94743686e-02, 1.29838950e+00, 8.23823258e-01, 4.37477503e-02, 7.57683549e-02, 7.40859765e-03, 6.75724102e+03)
names(boot_sd) = c("nu1", "nu2", "b1", "b2", "m12", "m21", "Ts", "Tsc", "O", "theta")

N = 1000
boot_ci = 1.960*(boot_sd/sqrt(N))

boot_ci["theta"]/(4*mu*L)

boot_ci["nu1"]*nref
boot_ci["nu2"]*nref

2*nref*boot_ci["Ts"]*g
2*nref*boot_ci["Tsc"]*g

boot_ci["m12"]/(2*nref)
boot_ci["m21"]/(2*nref)


boot_ci["b1"]
boot_ci["b2"]


