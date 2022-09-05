install.packages("vegan")
library(vegan)
install_github('bioBakery/Maaslin2', force = TRUE)
install_github('omicsEye/omicsArt', force = TRUE)
setwd("/Users/Rano/Desktop/Single_Cell_Wound/")


##This code is to generate the distance matrix 
data_dir <- '/Users/Rano/Desktop/Single_Cell_Wound/'
list.files(data_dir)

wound1 <- read.table(file = 'analysis/data/data_Wound1.tsv', sep = '\t', header = TRUE)
meta_wound1 <- read.table(file = 'analysis/data/meta_data_Wound1.tsv', sep = '\t', header = TRUE)
#wound1 <- wound1[apply(wound1[,-1], 1, function(x) !all(x==0)),]
#wound1 <-wound1[, colSums(wound1 == 0)/nrow(wound1) < .9, drop = FALSE]
#wound1 <- wound1[rowSums(wound1[])>0,]
wound1 <- wound1[1:100,]


rownames(meta_nonwound1) <- meta_nonwound1[,1]
meta_nonwound1[,1] <- NULL

rownames(nonwound1) <- nonwound1[,1]
nonwound1[,1] <- NULL

final <- rownames(nonwound1)
meta_data_nonwound <- subset(meta_nonwound1, rownames(meta_nonwound1)%in% final)

final <- rownames(meta_data)
data_wound1 <- subset(wound1, rownames(wound1)%in% final)


wound1_samp<-nrow(wound1_samp[wound1_samp<0,])
wound1_samp[wound1_samp<0] <- 0
dist_nonwound <- as.matrix(vegdist(nonwound1, method="bray"))
write.table(dist_nonwound, file = "dist_nonwound1.txt", sep = "\t",
            row.names = TRUE, col.names = T)
write.table(meta_data_nonwound, file = "dist_meta_nonwound1.txt", sep = "\t",
            row.names = TRUE, col.names = T)



nonwound1 <- read.table(file = 'analysis/data/data_Nonwound1.tsv', sep = '\t', header = TRUE)
meta_nonwound1 <- read.table(file = 'analysis/data/meta_data_Nonwound1.tsv', sep = '\t', header = TRUE)
#wound1 <- wound1[apply(wound1[,-1], 1, function(x) !all(x==0)),]
#wound1 <-wound1[, colSums(wound1 == 0)/nrow(wound1) < .9, drop = FALSE]
#wound1 <- wound1[rowSums(wound1[])>0,]
nonwound1 <- nonwound1[1:100,]


rownames(meta_nonwound1) <- meta_wound1[,1]
meta_wound1[,1] <- NULL

rownames(wound1) <- wound1[,1]
wound1[,1] <- NULL

final <- rownames(wound1)
meta_data <- subset(meta_wound1, rownames(meta_wound1)%in% final)

final <- rownames(meta_data)
data_wound1 <- subset(wound1, rownames(wound1)%in% final)


wound1_samp<-nrow(wound1_samp[wound1_samp<0,])
wound1_samp[wound1_samp<0] <- 0
dist <- as.matrix(vegdist(data_wound1, method="bray"))
write.table(dist, file = "dist_wound1.txt", sep = "\t",
            row.names = TRUE, col.names = T)
write.table(dist, file = "dist_meta_wound1.txt", sep = "\t",
            row.names = TRUE, col.names = T)