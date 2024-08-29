# 10 March 2021
# Nancy Chen
# get Beadchip SNP positions
# run on file for Dovetail assembly 27 Sept 2021

args = commandArgs(trailingOnly=TRUE)

# get actual flanking seq
#beadchip <- read.csv('finalFSJbeadchip.csv', header = TRUE, 
	#stringsAsFactors = FALSE)
beadchip <- read.csv(args[1], header = TRUE,
        stringsAsFactors = FALSE)

beadchip$LeftFlank <- nchar(matrix(unlist(strsplit(beadchip$Sequence, "\\[")),
	ncol = 2, byrow = TRUE)[, 1])
beadchip$RightFlank <- nchar(matrix(unlist(strsplit(beadchip$Sequence, "\\]")),
	ncol = 2, byrow = TRUE)[, 2])

# read in locations file
# loc <- read.table('FSJbeadchipSeqLocDovetail02Jun2018SNP.txt', header = TRUE,
#	stringsAsFactors = FALSE, sep = "\t")
#loc <- read.table('FSJbeadchipSeqLoc.txt', header = TRUE,
	#stringsAsFactors = FALSE, sep = "\t") # note first edit SNP names

loc <- read.table(args[2], header = TRUE,
        stringsAsFactors = FALSE, sep = "\t") # note first edit SNP names

data <- merge(loc[, 1:4], beadchip[, c('Locus_Name', 'LeftFlank', 'RightFlank')], 
	by.x = 'SNPname', by.y = 'Locus_Name')

# fix SNP positions
data$SNPpos <- -1
data[data$ProbeOrientation == 0, 'SNPpos'] <- data[data$ProbeOrientation == 0,
	'ProbePos'] + data[data$ProbeOrientation == 0, 'LeftFlank']

# rev comp
data[data$ProbeOrientation >= 16, 'SNPpos'] <- data[data$ProbeOrientation >= 16, 
	'ProbePos'] + data[data$ProbeOrientation >= 16, 'RightFlank']

# write.table(data, file = 'FSJbeadchipSeqLocDovetail02Jun2018SNPfixed.txt', 
#	quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
#write.table(data, file = 'FSJbeadchipSeqLocFSJgenomeV2_06March2023Fixed.txt',
	#quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(data, file = args[3],
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
