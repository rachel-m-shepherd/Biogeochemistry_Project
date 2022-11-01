### Assessing variation in Variables
# read in mapping file
mapping <-read.csv("META_Panama_B.csv", comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

jpeg(file="Fe.jpeg")
hist(x=mapping$Fe,
     main = "",
     ylab = "",
     xlab="Fe",
     ylim=c(0,308),
     col="#CC6677",
     freq=TRUE)
dev.off()
jpeg(file="Zn.jpeg")
hist(x=mapping$Zn,
     main = "",
     ylab = "",
     xlab="Zn",
     ylim=c(0,308),
     col="#6699CC",
     freq=TRUE)
dev.off()
## ALL OTHER
jpeg(file="pH.jpeg")
hist(x=mapping$pH,
     main = "",
     ylab = "",
     xlab="pH",
     ylim=c(0,308),
     col="#888888",
     freq=TRUE)
dev.off()
jpeg(file="Al.jpeg")
hist(x=mapping$Al,
     main = "",
     ylab = "",
     xlab="Al",
     ylim=c(0,308),
     col="#332288",
     freq=TRUE)
dev.off()
### Assessing correlation between pH and variables
library("ggpubr")
jpeg(file="Al_by_pH.jpeg")
ggscatter(mapping, x = "pH", y = "Al", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "pH", ylab = "Al",
          color =  "#332288")
dev.off()

jpeg(file="Fe_by_pH.jpeg")
ggscatter(mapping, x = "pH", y = "Fe", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "pH", ylab = "Fe",
          color =  "#CC6677")
dev.off()

jpeg(file="Zn_by_pH.jpeg")
ggscatter(mapping, x = "pH", y = "Zn", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "pH", ylab = "Zn",
          color =  "#6699CC")
dev.off()

#### Community composition by variables
#### MANTEL
otu_taxa <- read.csv("Bact_otu_B.csv", comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
otu_only <- otu_taxa[,1:308]##select for only the count data
table<- otu_only[,order(colnames(otu_only))]

##import map file
mapping <-read.csv("META_Panama_B.csv", comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
map <- mapping[order(rownames(mapping)),]

## MAKE SURE ROW NAMES AND COLUMN NAMES MATCH
rownames(map) == colnames(table)

otu <-table
#### Make distances --ASV Table
bray = vegdist(t(otu), method = "bray")
#### Make Distances --Variables
pH = dist(map$pH, method = "euclidean")
Al = dist(map$Al, method = "euclidean")
Fe = dist(map$Fe, method = "euclidean")
Zn = dist(map$Zn, method = "euclidean")
### RUN MANTELS
vegan::mantel(bray, pH, method="spearman", permutations=999)
vegan::mantel(bray, Al, method="spearman", permutations=999)
vegan::mantel(bray, Fe, method="spearman", permutations=999)
vegan::mantel(bray, Zn, method="spearman", permutations=999)

### Record observation # (Rho) for each variable in Excel File
## Use this to then plot the results
#### PLOT MANTEL RESULTS
library(ggplot2)
mant <- read.csv("MantelNumbers.csv", comment.char="", header=T, stringsAsFactors=T, check.names=FALSE)
mant$Variable <- factor(mant$Variable, levels = c("pH","Al", "Fe","Zn"))

p <- ggplot(data=mant, aes(x=Variable, y=Rho, fill=Variable)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal()
p + theme(axis.title.x = element_blank())+ylab("Community Composition Association with Variable (Rho)")+
  guides(fill=guide_legend(element_blank()))+
  scale_fill_manual("legend", values = c("pH" = "#888888", "Al" = "#332288","Fe" = "#CC6677","Zn" = "#6699CC"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

