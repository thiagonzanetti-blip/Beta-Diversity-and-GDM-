
# Script Bray Beta diversity and GDMs PARNA-CD


setwd("/Users/thiagozanetti/Documents/Mestrado/Análises/Dados")


# Betadiversity Bray-curtis and partitioning
##############


library(betapart)
library(ggplot2)
library(vegan)
library(scales)
library(magick)
library(geosphere)
library(reshape2)
library(broom)
library(ggpubr)
library(sjPlot)

#Beta diversity Bray-Curtis


matrix<-read.csv("AbundCota.csv", stringsAsFactors=FALSE, encoding="UTF-8",sep=",")
attach(matrix) 

beta_core_abund<-betapart.core.abund(matrix[,-1])

beta.multi.abund<-beta.multi.abund(matrix[,-1], index.family="bray")

matrix_samp_abund<-beta.multi.abund(beta_core_abund, index.family="bray")

beta.sample.abund <- beta.sample.abund(matrix[,-1], index.family = "bray", sites=5, samples = 171)

dist<-beta.sample.abund$sampled.values

beta.pair.abund<-beta.pair.abund(matrix[,-1],"bray")

#Elevational distance x Dissimilarity

dissim <- vegdist(matrix[,-1], method = "bray")
matriz_dissim <- as.matrix(dissim)  

coordenadas_altitude <- data.frame(
  Area = c("1", "2", "3", "4", "5"),
  Longitude = c(-41.34039, -41.39905, -41.32481, -41.32744, -41.48656),
  Latitude = c(-12.75886, -12.54414, -12.88654, -12.93656, -12.59905),
  Altitude = c(400, 600, 800, 1100, 1300)  #meters
)

pares <- melt(matriz_dissim)
names(pares) <- c("Area1", "Area2", "Dissimilaridade")

pares <- merge(pares, coordenadas_altitude[, c("Area", "Altitude")], 
               by.x = "Area1", by.y = "Area")
pares <- merge(pares, coordenadas_altitude[, c("Area", "Altitude")], 
               by.x = "Area2", by.y = "Area", suffixes = c("1", "2"))


pares$Diff_Altitude <- abs(pares$Altitude1 - pares$Altitude2)


cor.test(pares$Dissimilaridade, pares$Diff_Altitude)

modelo <- lm(Dissimilaridade ~ Diff_Altitude, data = pares)
summary(modelo)


tidy(modelo)


sjPlot::tab_model(modelo, 
                  show.se = TRUE, 
                  show.stat = TRUE,
                  digits = 3,
                  pred.labels = c("Intercept", "Elevation difference"),
                  dv.labels = "Dissimilarity")


ggplot(pares, aes(x = Diff_Altitude, y = Dissimilaridade)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Relação entre Dissimilaridade e Diferença de Altitude",
    x = "Altitude difference",
    y = "Dissimilarity index"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


plot <- ggplot(pares, aes(x = Diff_Altitude, y = Dissimilaridade)) +
  geom_point(size = 3, alpha = 0.7, shape = 19) + 
  geom_smooth(method = "lm", color = "red", se = TRUE, fill = "lightpink") +
  labs(
    x = "Altitude difference (m)", 
    y = "Dissimilarity index",
    caption = "Error band represents 95% confidence interval"
  ) +
  theme_minimal(base_size = 12) + 
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"), 
    aspect.ratio = 0.8
  ) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()   
  )


plot_with_stats <- plot + 
  stat_regline_equation(
    aes(label = ..rr.label..),  
    label.x = 0.7,       
    label.y = 0.9,          
    size = 4,               
    color = "black",      
  )

plot_with_stats

ggsave("dissimilarity_vs_altitude.png", plot_with_stats, 
       width = 8, height = 6, dpi = 300)




#dendroclust + partitioning

#####
novos_rotulos <- c("1"= "400", "2"="600", "3"="800","4"= "1100","5"= "1300") 
hcB <- hclust(beta.pair.abund$beta.bray.bal, method = "average")
hcG <- hclust(beta.pair.abund$beta.bray.gra, method = "average")
hc <- hclust(beta.pair.abund$beta.bray, method = "average")

hc$labels <- novos_rotulos

dend <- as.dendrogram(hcB)

cores <- c("#FDE725FF", "#21908CFF", "#440154FF")  

dend <- color_branches(dend, k = 3, col = cores)

plot(dend, 
     xlab = "Elevation (m)", ylab = "Dissimilarity (%)")

png("dendrograma.png", width = 3000, height = 2000, res = 300)  
plot(dend, 
     xlab = "Elevation (m)", ylab = "Dissimilarity (%)")

rect.dendrogram(dend, k = 3, border = "gray30", lty = 2, lwd = 1.5, upper_rect = 0.1 )
grupos <- cutree(dend, k = 3)
posicoes_x <- tapply(order.dendrogram(dend), grupos, mean)
text(
  x = posicoes_x,
  y = 0.71 * max(get_nodes_attr(dend, "height")),
  labels = paste("Group", 1:3),
  pos = 4, col = "black", font = 2, cex = 0.8
)



dados <- matrix[,-1]
beta_result <- beta.pair.abund(dados, index.family = "bray")

balanced <- as.matrix(beta_result$beta.bray.bal)
gradient <- as.matrix(beta_result$beta.bray.gra) 
total_beta <- as.matrix(beta_result$beta.bray)
                        

contrib <- data.frame(
  Componente = factor(c("balanced", "gradient"), 
                      levels = c("balanced", "gradient")),
  Proporção = c(0.86, 0.14))

Particao <- ggplot(contrib, aes(x = Componente, y = Proporção, fill = Componente)) +  
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = scales::percent(Proporção)), 
            position = position_stack(vjust = 0.5), 
            size = 3.5, color = "black") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.9)) +
  scale_fill_manual(
    name = NULL,
    values = c("#7F7F7F", "#D62728"),
    labels = c(bquote(beta["Bray"]^"Bal"), bquote(beta["Bray"]^"Grad"))
  ) +
  labs(
    x = NULL,
    y = "Beta diversity partitioning",
  ) +
  theme_minimal(base_size = 14) +
  theme(
    
    legend.position = "right", 
    legend.text = element_text(size = 14, margin = margin(r = 8)),  
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(family = "Times")
  ) +
  guides(fill = guide_legend(override.aes = list(color = NA)))  

ggsave("PartiçãoBDiversidade.png", Particao, 
       width = 8, height = 6, dpi = 300)



img1 <- image_read("/Users/thiagozanetti/Documents/Mestrado/Análises/Dados/PartiçãoBDiversidade.png")
img2 <- image_read("/Users/thiagozanetti/Documents/Mestrado/Análises/Dados/dendrograma.png")

img1 <- image_resize(img1, "x1500")  
img2 <- image_resize(img2, "x1500")                        

GraphDissimili <- image_append(c(img2, img1))
print(GraphDissimili)
image_write(GraphicDiss, "GraphDissimili.png")

#####

rm(list=ls())

# GDM 
############

library(gdm)
library(ggplot2)
library (farver)
library(magick)


matrix<-read.csv("AbundCota.csv", stringsAsFactors=FALSE, 
                 encoding="UTF-8",sep=",")
attach(matrix)

EnvData <-read.csv("Envmatrix.csv", stringsAsFactors=FALSE, encoding="UTF-8", 
                   header=T, sep=",")
attach(EnvData)

gdmData_Bray <-gdm::formatsitepair(bioData = matrix, 
                                   bioFormat = 1, 
                                   dist = "bray",
                                   abundance = TRUE, 
                                   siteColumn= "Plot", 
                                   XColumn = "Lon", 
                                   YColumn = "Lat", 
                                   predData = EnvData)


gdm_Model<-gdm::gdm(data = gdmData_Bray, geo=TRUE, splines=NULL)

gdm_Model$explained/100 
gdm_Model$gdmdeviance/gdm_Model$nulldeviance 
gdm_Model$nulldeviance


gdm_Model_splineDat <- gdm::isplineExtract(gdm_Model)

png("GDMImportance.png", width = 6, height = 10, units = "in", res = 300)

par(mfcol = c(3, 1), oma = c(4, 4, 2, 2)) 

plot(gdm_Model_splineDat$x[,"Treecover"], 
     gdm_Model_splineDat$y[,"Treecover"], 
     lwd=2, type="l", col="red", ylim=c(0, 1),
     xlab="Treecover Mean", 
     ylab="",
     main="",
     cex.lab=1.7)

plot(gdm_Model_splineDat$x[,"Annual_Temperature"], 
     gdm_Model_splineDat$y[,"Annual_Temperature"],  
     lwd=3, type="l", col="red", ylim=c(0, 1),
     xlab="Annual Temperature (°C)", 
     ylab="",
     main="",
     cex.lab=1.7)

plot(gdm_Model_splineDat$x[,"Annual_precipitation"], 
     gdm_Model_splineDat$y[,"Annual_precipitation"],  
     lwd=3, type="l", col="red", ylim=c(0, 1),
     xlab="Annual Precipitation (mm)", 
     ylab="",
     main="",
     cex.lab=1.7)

mtext("Partial ecological distance", side = 2, line = 1, outer = TRUE, 
      cex = 1.2)

dev.off()


names(gdmData_Bray)
names(EnvData)
EnvData_filtered<-data.frame(EnvData[,-1])


EnvData_Transformed <- gdm::gdm.transform(model = gdm_Model, 
                                          data = EnvData_filtered)

summary(EnvData_Transformed)
PredMinMax<-Rfast::colMinsMaxs(EnvData_Transformed) 
PredRange<-PredMinMax[2,] - PredMinMax[1,] 
names(PredRange)<-names(EnvData_filtered) 

EnvDataList<-list()
EnvDataList[[1]]<-c("Lat","Lon",  "Annual_precipitation", "Solar_Radiation", 
                    "Annual_Temperature"
                    ,"Treecover")

names(EnvDataList)<-"Environment"

VarPartition<-gdm.partition.deviance(sitePairTable = gdmData_Bray,
                                     varSets=EnvDataList, partSpace=TRUE)
VarPartition


Barplot<-data.frame(Category = c("Pure space", "Shared", "Environment predictors"),
                    Value = c(VarPartition[5,]$DEVIANCE,
                              (VarPartition[3,]$DEVIANCE - VarPartition[5,]$DEVIANCE - VarPartition[4,]$DEVIANCE), 
                              VarPartition[4,]$DEVIANCE),
                    Group = c("","",""))


Barplot$Category<-factor(Barplot$Category,
                         levels=c("Pure space", "Shared", "Environment predictors"))


MyPlot1<-ggplot(data=Barplot, aes(x=Group, y=Value, fill=Category)) +
  geom_bar(position = "stack", stat="identity") +
  xlab("") + ylab("Variation explained") +
  scale_fill_manual(values=c("black", "gray50", "gray80")) +
  coord_flip() +  
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=10),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=10),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=12, colour="black", face="bold",
                                  margin=margin(t=5, r=0, b=5, l=0)), 
        legend.position="top",
        legend.title=element_blank(),
        plot.background=element_rect(fill="white"),
  ) +
  guides(fill=guide_legend(title.position=NULL, reverse = T, hjust=0, nrow=1,
                           byrow=F)); MyPlot1

ggsave("BarplotsVarPartition.png", plot=MyPlot1, width=6, height=2, 
       units="in", bg = "transparent")

dev.off()


gdm_VariableImportance<-gdm::gdm.varImp(spTable = gdmData_Bray,
                                        geo = TRUE,
                                        splines = NULL,
                                        knots = NULL,
                                        predSelect = FALSE, 
                                        nPerm = 1000,
                                        pValue = 0.05,
                                        parallel = TRUE,
                                        cores = 5,
                                        sampleSites = 1, 
                                        sampleSitePairs = 1, 
                                        outFile = NULL)

gdm_VariableImportance

PredContribImport<-merge(data.frame(Predictor=row.names(gdm_VariableImportance[[2]]),
                                    Importance=gdm_VariableImportance[[2]][,1],
                                    Importance_p=gdm_VariableImportance[[3]][,1]),  
                         
                         data.frame(Predictor=row.names(as.data.frame(PredRange)), 
                                    Contrib=PredRange),
                         
                         by="Predictor", all.x=TRUE)

ggplot(data=PredContribImport, aes(x=Importance, y=Contrib)) + 
  geom_point(shape=21, aes(size=(1/(Importance_p+0.001))), color="black", fill="gray") + 
  geom_point(data=PredContribImport[PredContribImport$Importance_p<=0.1,], 
             shape=21, fill="red", color="black", aes(size=(1/(Importance_p+0.001)))) +
  ggrepel::geom_text_repel(label=PredContribImport$Predictor, 
                           nudge_x=0.1, nudge_y=0.01) + 
  xlab("Variable Importance") + ylab("Variable Contribution") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"), 
        axis.text.y = element_text(hjust=1, vjust=0.5, angle=0, size=12),
        axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=0, size=12),
        axis.title = element_text(hjust=0.5, vjust=0.5, angle=0, size=14),
        legend.position = "top"
  )

cor(EnvData[ c("Annual_precipitation", "Treecover", "Annual_Temperature",
               "Solar_Radiation")])


img1 <- image_read("/Users/thiagozanetti/Documents/Mestrado/Análises/Dados/GDMImportance.png")
img2 <- image_read("/Users/thiagozanetti/Documents/Mestrado/Análises/Dados/BarplotsVarPartition.png")

img1 <- image_resize(img1, "x1500")  
img2 <- image_resize(img2, "x1500")                        

GraphGDM <- image_append(c(img2, img1))
print(GraphGDM)
image_write(GraphicDiss, "GraphGDM.png")


### The end ###
