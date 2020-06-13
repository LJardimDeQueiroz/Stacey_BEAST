## Supplementary data 6 – Script to process STACEY output using speciesDA
# Adapted from Jones et al. (2015)
# If you use this script, cite Jones et al. (2015) and us (Jardim de Queiroz et al. 2020)

setwd("path/")

workdir <- "path/"

# Reading speciesDA output:
x <- read.table(paste("clusterings.txt", sep=""), header=TRUE, check.names=FALSE)

#Some adjustes:
x1 <- x[,1:4] # x1 contains only the first 4 columns

x <- x[-(1:4)] # x contains now only the columns related to the samples

# Before uploading your phylogenetic tree (the annotated tree obtained after running STACEY), order the clades  according to your wish; you'll use such to re-order the speciesDA output. Than save it in nexus format.

# Reading nexus tree.
library(ape)

tree <- read.nexus("tree.nexus") 
tree$tip.label #this is the order we want our speciesDA output

tips <- data.frame(cols = c(tree$tip))

# Minimal cluster names (samples)
mincl.names <- colnames(x)

# Check if the names in the tree and in the table are equal:
matrix.names <- as.data.frame(sort(colnames(x)))
colnames(matrix.names) <- c("names")
tree.names <- as.data.frame(sort(tree$tip))
colnames(tree.names) <- c("names")

all.equal(tree.names,matrix.names) #TRUE=everything is correct

# If typos are detect, to identify them, use:
tree.names==matrix.names
#AND/OR
tree.names[matrix.names != tree.names]


# To order the matrix according to the order in the tree

# Requrires data.table
library(data.table)

# Convert data.frame to data.table
x <- setDT(x)

# Order the columns acoording to the order of the tips in the tree
x  <- setcolorder(x, as.character(tips$cols))

# Put it back in the data.frame format
x <- as.data.frame(x)
is.data.frame(x)

# Check now if your matrix has the correct order (according to your nexus trees)
# Requires 'compare'
require(compare)
mincl.names <- colnames(x)
	compare(mincl.names,tips$cols, allowAll=TRUE)
		#OR/AND
mincl.names==tips$cols
			
# Merging x1 and x in order to have the data.frame as the original
x <- cbind(x1,x, all=TRUE)
x$all <- NULL #to eliminate a not useful column that was added to the data.frame

#Prepare 'renames' (I do it in excel):
#mincl.names <- colnames(x)[-(1:4)]
#cat(mincl.names,sep="\n") # This lists the samples I will copy to prepare the list in excel.

#Actually, I will not rename the samples, i.e., I'm keeping the original name. But I kept this step because it was in the original scritp.
renames <- matrix(c(
"sample1",	"sample1",
"sample2",	"sample2",
"sample3",	"sample3",
nrow=ncol(x)-4, ncol=2, byrow=TRUE)

# Check for typos, etc
for (i in 1:length(mincl.names)) {
	stopifnot(mincl.names[i] == renames[i,1])
}

# Make the similarity matrix
displaynames <- renames[,2]
nmincls <- length(displaynames)

sim <-matrix(0,ncol=nmincls,nrow=nmincls,dimnames=list(displaynames,displaynames))

for (i in 1:nmincls) {
	for (j in 1:nmincls) {
		coli <- x[,mincl.names[i]]
		colj <- x[,mincl.names[j]]
		w <- coli == colj
		sim[i,j] <- sum(x[w,"fraction"])
	}
}

# Ensure rounding errors don't make probabilities sum to more than 1.
sim <- pmin(sim,1)

# Change the order of minimal clusters
	# We don't need to re-order because we have already ordered the matrix before; so you want just to assume the currently order, i.e., from 1 to n. 		
neworder <- c(1:(ncol(x)-4)) #actually, no change is needed here.

# Currently recognised groups
# dividers <- c(0, 6, 8, 12, 13)

plot.rectangle <- function(v1,v2,...)
{
	polygon(c(v1[1],v2[1],v2[1],v1[1]), c(v1[2],v1[2],v2[2],v2[2]), ...)
}

# Main plotting routine;
# To change size of the font, alter cex option.
# To better place the sample names, alter "line" in 'axis'
plot.simmatrix <- function() {
  par(mar= c(0,5,5,0)+.1, cex=.75)
  plot(NULL, xlim=c(0,nmincls), ylim=c(nmincls,0), axes=FALSE, ylab="", xlab="")
  axis(3, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-2)
  axis(2, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-2)
  for (i in 1:nmincls) {
    for (j in 1:nmincls) {
      d <- 1 - sim[neworder[i],neworder[j]]
      plot.rectangle(c(i-1,j-1), c(i,j), col=rgb(d,d,d), border="white")
  } }
}
#  for (b in dividers) {
#    lines(x=c(-.5,nmincls), y=c(b,b))
#    lines(x=c(b,b), y=c(-.5,nmincls))
# } }

# Display as text, on screen, and to PDF file
matrix <- print(sim[neworder,neworder], digits=2)
similarity_matrix <- as.data.frame(as.matrix(matrix))
write.csv(similarity_matrix, file = "similarity_matrix.csv")
plot.simmatrix()
#pdf(file=paste(workdir, "simmatrix.pdf", sep=""))
plot.simmatrix()
#dev.off()





