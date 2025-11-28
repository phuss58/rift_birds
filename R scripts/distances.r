library(vcfR)
library(adegenet)
library(StAMPP)

# loop for each of the 8 species to read in the vcf, create a distance matrix, and output in a
# format compatible with splitstree



x_files <- list.files(pattern="*nothin.vcf")
for(a in 1:length(x_files)) {
	# read vcf
	a_rep <- read.vcfR(x_files[a])
	
	# convert to genlight
	a_rep <- vcfR2genlight(a_rep)
	
	# give fake population names (not used anyway)
	pop(a_rep) <- a_rep@ind.names
	
	# calculate Nei's distance 
	a_rep <- stamppNeisD(a_rep, pop = FALSE)
	
	# shorten names for phylip output
    rownames_file <- rownames(a_rep)
    
    for(j in 1:length(rownames_file)){
        nr1 <- rownames_file[j]
        nr2 <- strsplit(nr1,"__")[[1]][2]
        nr3 <- gsub("FMNH","F",nr2)
        nr4 <- gsub("_","",nr3)
        rownames_file[j] <- nr4
    }
    
	# rownames(a_rep) <- paste0(substr(sapply(strsplit(rownames(a_rep), "_"), "[[", 1), 1, 1), "_",mkdir
	# substr(sapply(strsplit(rownames(a_rep), "_"), "[[", 2), 1, 1), "_",
	# substr(sapply(strsplit(rownames(a_rep), "__"), "[[", 2), 1, nchar(sapply(strsplit(rownames(a_rep), "__"), "[[", 2)) - 1))
    
    rownames(a_rep) <- rownames_file
	
	# output name
	outname <- paste0(strsplit(x_files[a], "\\.")[[1]][1], "_distmat.phy")
	
	# write output distance matrix in phylip format
	stamppPhylip(distance.mat=a_rep, file= outname)
	
}