args = commandArgs(trailingOnly=TRUE)

lnL_vals <- read.table( file = args[2], sep= " ", stringsAsFactors = FALSE, 
                        header = FALSE )
row.names( lnL_vals ) <- c( "M0-branch", "Human-branch", "Baboon-branch", "Mouse-branch", "Rat-branch", "Sheep-branch", "Cow-branch", "Pig-branch", "Chicken-branch", "Chimp-branch", "Gibbon-branch", "Opossum-branch", "HumanBaboon-branch","MouseRat-branch","CowSheep-branch", "CowSheepPig-branch", "Primate-branch","Mammal-branch", "Placental-branch")


calculate.branch.stats <- function(branch){
  df <- data.frame()
  diff <- 2*(lnL_vals[1,] - lnL_vals[branch,])
  pval <- pchisq(diff, df =1, lower.tail=F)
  df <- cbind(diff,pval)
}

human.branch <- calculate.branch.stats(2)
baboon.branch <- calculate.branch.stats(3)
mouse.branch <- calculate.branch.stats(4)
rat.branch <- calculate.branch.stats(5)
sheep.branch <- calculate.branch.stats(6)
cow.branch <- calculate.branch.stats(7)
pig.branch <- calculate.branch.stats(8)
chicken.branch <- calculate.branch.stats(9)
chimp.branch <- calculate.branch.stats(10)
gibbon.branch <- calculate.branch.stats(11)
opossum.branch <- calculate.branch.stats(12)
humanbaboon.branch <- calculate.branch.stats(13)
mouserat.branch <- calculate.branch.stats(14)
cowsheep.branch <- calculate.branch.stats(15)
cowsheeppig.branch <- calculate.branch.stats(16)
primate.branch <- calculate.branch.stats(17)
mammal.branch <- calculate.branch.stats(18)
placental.branch <- calculate.branch.stats(19)

tab <- as.data.frame(cbind(human.branch,baboon.branch,mouse.branch,rat.branch,sheep.branch,cow.branch,pig.branch,chicken.branch,chimp.branch,gibbon.branch,opossum.branch,humanbaboon.branch,mouserat.branch,cowsheep.branch,cowsheeppig.branch,primate.branch,mammal.branch,placental.branch))
names(tab) <- c("human.difference","human.p-value","baboon.difference","baboon.p-value","mouse.difference","mouse.p-value","rat.difference","rat.p-value","sheep.difference","sheep.p-value","cow.difference","cow.p-value","pig.difference","pig.p-value","chicken.difference","chicken.p-value","chimp.difference","chimp.p-value","gibbon.difference","gibbon.p-value","opossum.difference","opossum.p-value","humanbaboon.difference","humanbaboon.p-value","mouserat.difference","mouserat.p-value","cowsheep.difference","cowsheep.p-value","cowsheeppig.difference","cowsheeppig.p-value","primate.difference","primate.p-value","mammal.difference","mammal.p-value","placental.difference","placental.p-value")

w_vals <- read.delim(file = args[3], stringsAsFactors = FALSE, 
                        header = TRUE )
majority.data <- cbind(tab,w_vals)

majority.data$gene <- args[1]


write.table(majority.data,file = paste0(args[1],".stats.table"), row.names = FALSE, quote = FALSE, sep = "\t")
