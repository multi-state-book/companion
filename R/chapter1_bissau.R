# Read data
bissau <- data.frame(read.csv("data/bissau.csv"))

#------------------------------------------------------------------#
#---------------- Table 1.1 ---------------------------------------#
#------------------------------------------------------------------#

with(bissau, table(bcg, dead))
rowSums(with(bissau, table(bcg, dead)))
colSums(with(bissau, table(bcg, dead)))
with(bissau, table(bcg, dead)) / rowSums(with(bissau, table(bcg, dead)))
colSums(with(bissau, table(bcg, dead))) / (sum(colSums(with(bissau, table(bcg, dead)))))
