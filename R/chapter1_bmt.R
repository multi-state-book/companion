# Read data
bmt <- data.frame(read.csv("data/bmt.csv"))

#------------------------------------------------------------------#
#---------------- Table 1.4 ---------------------------------------#
#------------------------------------------------------------------#

with(bmt, table(rel))
with(bmt, table(rel)) / sum(with(bmt, table(rel)))

with(bmt, table(death))
with(bmt, table(death)) / sum(with(bmt, table(death)))

with(bmt, table(rel,death))
with(bmt, table(rel,death)) / rowSums(with(bmt, table(rel,death)))

with(bmt, table(gvhd))
with(bmt, table(gvhd)) / sum(with(bmt, table(gvhd)))

with(bmt, table(gvhd, rel))
with(bmt, table(gvhd, rel)) / rowSums(with(bmt, table(gvhd, rel)))

with(bmt, table(gvhd, death))
with(bmt, table(gvhd, death)) / rowSums(with(bmt, table(gvhd, death)))


