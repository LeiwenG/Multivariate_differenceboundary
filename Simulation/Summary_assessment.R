
# Import assessment results
WAIC_ARDP_dagar_indc = readRDS("WAIC_ARDP_dagar_indc25_v1.rds")
D_ARDP_dagar_indc = readRDS("D_ARDP_dagar_indc25_v1.rds")
KL_ARDP_dagar_indc = readRDS("KL_ARDP_dagar_indc25_v1.rds")

WAIC_ARDP_car_indc = readRDS("WAIC_ARDP_car_indc25.rds")
D_ARDP_car_indc = readRDS("D_ARDP_car_indc25.rds")
KL_ARDP_car_indc = readRDS("KL_ARDP_car_indc25.rds")

WAIC_ARDP_Jdagarc = readRDS("WAIC_ARDP_Jdagarc25.rds")
D_ARDP_Jdagarc = readRDS("D_ARDP_Jdagarc25.rds")
KL_ARDP_Jdagarc = readRDS("KL_ARDP_Jdagarc25.rds")

WAIC_ARDP_Jcarc = readRDS("WAIC_ARDP_Jcarc25_1.rds")
D_ARDP_Jcarc = readRDS("D_ARDP_Jcarc25_1.rds")
KL_ARDP_Jcarc = readRDS("KL_ARDP_Jcarc25_1.rds")

# Plots for WAIC
WAIC_value = c(WAIC_ARDP_dagar_indc, WAIC_ARDP_Jdagarc, WAIC_ARDP_car_indc, WAIC_ARDP_Jcarc)
Type = c(rep("DAGAR",60),rep("CAR",60))
Model = c(rep("Independent-disease",30), rep("Joint",30), rep("Independent-disease",30), rep("Joint",30))
Full_model = paste(Type, "-", Model, sep="")
df = data.frame(Model)
df$WAIC = WAIC_value
df$Type = Type
df$Full_model= Full_model

library(plyr)
WAIC_mu <- ddply(df, "Full_model", summarise, WAIC.mean=mean(WAIC))
df1 = merge(df, WAIC_mu, by = "Full_model")

pdf("ARDP_WAIC25.pdf", height = 5, width = 12)
ggplot(df1, aes(x=WAIC, color=Model,fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.5) + 
  geom_vline(data=df1, aes(xintercept=WAIC.mean, color=Model),
             linetype="dashed") +
  facet_wrap(~Type) +
  xlab("WAIC") + ylab("Density")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.position = "none",
        strip.text.x = element_text(size = 14))
dev.off()

# Plots for D score
D_value = c(D_ARDP_dagar_indc, D_ARDP_Jdagarc,D_ARDP_car_indc, D_ARDP_Jcarc)
df = data.frame(Model)
df$D = D_value
df$Type = Type
df$Full_model= Full_model
D_mu <- ddply(df, "Full_model", summarise, D.mean=mean(D))
df1 = merge(df, D_mu, by = "Full_model")


pdf("ARDP_D25.pdf", height = 5, width = 12)
ggplot(df1, aes(x=D, color=Model,fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.5) + 
  geom_vline(data=df1, aes(xintercept=D.mean, color=Model),
             linetype="dashed") +
  facet_wrap(~Type) +
  xlab("D score") + ylab("Density")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        legend.position = "none",
        strip.text.x = element_text(size = 14))
dev.off()

# Plots for Kullback-Leibler divergence
KL_value = c(rowMeans(KL_ARDP_dagar_indc), rowMeans(KL_ARDP_Jdagarc),rowMeans(KL_ARDP_car_indc), rowMeans(KL_ARDP_Jcarc))
df = data.frame(Model)
df$KL = KL_value
df$Type = Type
df$Full_model= Full_model
KL_mu <- ddply(df, "Full_model", summarise, KL.mean=mean(KL))
df1 = merge(df, KL_mu, by = "Full_model")


pdf("ARDP_KL25.pdf", height = 5, width = 14)
ggplot(df1, aes(x=KL, color=Model,fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.5) + 
  geom_vline(data=df1, aes(xintercept=KL.mean, color=Model),
             linetype="dashed") +
  facet_wrap(~Type) +
  xlab("KL") + ylab("Density")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.title=element_text(size=14),
        legend.text=element_text(size=14),
        strip.text.x = element_text(size = 14))
dev.off()

