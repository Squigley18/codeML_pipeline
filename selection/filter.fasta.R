library(seqinr)
args = commandArgs(trailingOnly=TRUE)

IDs <- read.delim(args[1],header = FALSE,sep=" ")
print(IDs)
print(str(IDs))

human.cds <- read.fasta("CDS.filtered/human.CDS.fa",as.string = TRUE)
sheep.cds <- read.fasta("CDS.filtered/sheep.CDS.fa",as.string = TRUE)
baboon.cds <- read.fasta("CDS.filtered/baboon.CDS.fa",as.string = TRUE)
mouse.cds <- read.fasta("CDS.filtered/mouse.CDS.fa",as.string = TRUE)
cow.cds <- read.fasta("CDS.filtered/cow.CDS.fa",as.string = TRUE)
opossum.cds <- read.fasta("CDS.filtered/opossum.CDS.fa",as.string = TRUE)
pig.cds <- read.fasta("CDS.filtered/pig.CDS.fa",as.string = TRUE)
rat.cds <- read.fasta("CDS.filtered/rat.CDS.fa",as.string = TRUE)
chicken.cds <- read.fasta("CDS.filtered/chicken.CDS.fa",as.string = TRUE)
chimp.cds <- read.fasta("CDS.filtered/chimp.CDS.fa",as.string = TRUE)
gibbon.cds <- read.fasta("CDS.filtered/gibbon.CDS.fa",as.string = TRUE)

human.transcript <- IDs[,6]
sheep.transcript <- IDs[,11]
baboon.transcript <- IDs[,1]
mouse.transcript <- IDs[,7]
cow.transcript <- IDs[,4]
opossum.transcript <- IDs[,8]
pig.transcript <- IDs[,9]
rat.transcript <- IDs[,10]
chicken.transcript <- IDs[,2]
chimp.transcript <- IDs[,3]
gibbon.transcript <- IDs[,5]


human <- human.cds[[human.transcript, exact = FALSE]]
human.seq <- human[1]
write.fasta(sequences = human.seq, name = "human", file.out = paste0(human.transcript,".human.fasta.txt"))

sheep <- sheep.cds[[sheep.transcript, exact = FALSE]]
sheep.seq <- sheep[1]
write.fasta(sequences = sheep.seq, name = "sheep", file.out = paste0(human.transcript,".sheep.fasta.txt"))

baboon <- baboon.cds[[baboon.transcript, exact = FALSE]]
baboon.seq <- baboon[1]
write.fasta(sequences = baboon.seq, name = "baboon", file.out = paste0(human.transcript,".baboon.fasta.txt"))


chimp <- chimp.cds[[chimp.transcript, exact = FALSE]]
chimp.seq <- chimp[1]
write.fasta(sequences = chimp.seq, name ="chimp", file.out = paste0(human.transcript,".chimp.fasta.txt"))


gibbon <- gibbon.cds[[gibbon.transcript, exact = FALSE]]
gibbon.seq <- gibbon[1]
write.fasta(sequences = gibbon.seq, name = "gibbon", file.out = paste0(human.transcript,".gibbon.fasta.txt"))


mouse <- mouse.cds[[mouse.transcript, exact = FALSE]]
mouse.seq <- mouse[1]
write.fasta(sequences = mouse.seq, name = "mouse", file.out = paste0(human.transcript,".mouse.fasta.txt"))


cow <- cow.cds[[cow.transcript, exact = FALSE]]
cow.seq <- cow[1]
write.fasta(sequences = cow.seq, name = "cow", file.out = paste0(human.transcript,".cow.fasta.txt"))


chicken <- chicken.cds[[chicken.transcript, exact = FALSE]]
chicken.seq <- chicken[1]
write.fasta(sequences = chicken.seq, name = "chicken", file.out = paste0(human.transcript,".chicken.fasta.txt"))


pig <- pig.cds[[pig.transcript, exact = FALSE]]
pig.seq <- pig[1]
write.fasta(sequences = pig.seq, name = "pig", file.out = paste0(human.transcript,".pig.fasta.txt"))

opossum <- opossum.cds[[opossum.transcript, exact = FALSE]]
opossum.seq <- opossum[1]
write.fasta(sequences = opossum.seq, name = "opossum", file.out = paste0(human.transcript,".opossum.fasta.txt"))

rat <- rat.cds[[rat.transcript, exact = FALSE]]
rat.seq <- rat[1]
write.fasta(sequences = rat.seq, name = "rat", file.out = paste0(human.transcript,".rat.fasta.txt"))


