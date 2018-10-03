

# Give the number +1 of protein with cleavage site (lower case letter in their sequence) == 608 -1 = 607
cat data/viral_proteins/Viruses_protein_db.faa | grep [a-z] | cut -d'>' -f1 | uniq | grep [A-Z] -v -c

# Give the number of protein annotated with a polyprotein outline = 607
cat results/stat_viral_protein/stat_proteins_Viruses.csv | cut -f3 | grep True -c

#Sanity check the number number need to be equal !! 
