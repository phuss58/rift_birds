

library(vegan)
library(ecodist)
library(fossil)
options(scipen=999)

popmap <- read.table("popmap_MRM.txt", header=T)
popmap[,1] <- substr(popmap[,1], 1, nchar(popmap[,1]) - 1)
popmap[,1] <- paste0(substr(sapply(strsplit(popmap[,1], "_"), "[[", 1), 1, 1), "_",
	substr(sapply(strsplit(popmap[,1], "_"), "[[", 2), 1, 1), "_",
	substr(sapply(strsplit(popmap[,1], "__"), "[[", 2), 1, nchar(sapply(strsplit(popmap[,1], "__"), "[[", 2))))

#########################################################
#########################################################
#########################################################
#cossypha
#########################################################
#########################################################
#########################################################
dist_table <- read.table("cossypha_nothin_distmat.phy", header=F, skip=1, row.names=1)
colnames(dist_table) <- rownames(dist_table)

ind1 <- c()
ind2 <- c()
gen_dist <- c()
geo_dist <- c()
biogeo <- c()
time <- c()
# loop for each row in the distance matrix
for(a in 1:(nrow(dist_table) - 1)) {
	a_rep <- dist_table[(a+1):nrow(dist_table),a]
	a_ind1 <- rep(rownames(dist_table)[a], length(a_rep))
	a_ind2 <- colnames(dist_table)[(a+1):nrow(dist_table)]
	
	# loop for each comparison
	for(b in 1:length(a_ind1)) {
		b_ind1 <- popmap[popmap[,2] %in% a_ind1[b],]
		b_ind2 <- popmap[popmap[,2] %in% a_ind2[b],]
		ind1 <- c(ind1, b_ind1$id)
		ind2 <- c(ind2, b_ind2$id)
		# genetic distance
		gen_dist <- c(gen_dist, a_rep[b])
		
		# geographic distance
		# long1, lat1, long2, lat2
		geo_dist <- c(geo_dist, deg.dist(b_ind1$long, b_ind1$lat, b_ind2$long, b_ind2$lat))
    # temporal
		if(popmap[popmap[,1] %in% b_ind1, 3] == popmap[popmap[,1] %in% b_ind2, 3]) {
			time <- c(time, 0)
		} else {
			time <- c(time, 1)
		}
		
		# biogeographic barriers
		if(b_ind1$geography == "Kabogo" & b_ind1$geography == "Kabogo") {
			biogeo <- c(biogeo, 0)
		} else if (b_ind1$geography == "Kabogo" & b_ind2$geography == "West") {
			biogeo <- c(biogeo, 1)
		} else if (b_ind1$geography == "Kabogo" & b_ind2$geography == "East") {
			biogeo <- c(biogeo, 2)
		} else if (b_ind1$geography == "West" & b_ind2$geography == "Kabogo") {
			biogeo <- c(biogeo, 1)
		} else if (b_ind1$geography == "West" & b_ind2$geography == "East") {
			biogeo <- c(biogeo, 1)
		} else if (b_ind1$geography == "West" & b_ind2$geography == "West") {
			biogeo <- c(biogeo, 0)
		} else if (b_ind1$geography == "East" & b_ind2$geography == "Kabogo") {
			biogeo <- c(biogeo, 2)
		} else if (b_ind1$geography == "East" & b_ind2$geography == "West") {
			biogeo <- c(biogeo, 1)
		} else if (b_ind1$geography == "East" & b_ind2$geography == "East") {
			biogeo <- c(biogeo, 0)
		}
	}
}

output <- data.frame(ind1=as.character(ind1), ind2=as.character(ind2), gen_dist=as.numeric(gen_dist), 
geo_dist=as.numeric(geo_dist), biogeo=as.numeric(biogeo), time=as.numeric(time))

# BELOW HERE ON PERSONAL MACHINE

# output_mod <- MRM(output$gen_dist ~ output$geo_dist + output$biogeo, nperm=100000)

output_mod <- MRM(gen_dist ~ geo_dist + biogeo, data = output, nperm=100000)

output_mod

varpart_mod <- varpart(output$gen_dist, output$geo_dist, output$biogeo, output$time)
varpart_mod
plot(varpart_mod)