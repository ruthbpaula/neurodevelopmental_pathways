################################
#### Code by: Ruth De Paula ####
#### 2022 ######################
#### ruthbpaula@gmail.com ######
################################


##### QUERYING REACTOME #####

# Export Reactome reference file as .tab into the computer

# To filter Homo sapiens results:
grep sapiens reactome_with_gene_names-chr_human-mouse-etc.tab > reactome_with_gene_names-chr_human.tab

# To get genes pathway information:
sed 's/\t/"\t"/g' reactome_with_gene_names-chr_human.tab | sed 's/^/"/g' | sed 's/$/"/g' > reactome_new

for i in *_genes ; do sed -i $'s/\r$//' $i ; done
for i in *_genes ; do grep -w -f $i reactome_new > ${i}_reactome ; done

# To get list of unique pathway ids:
for i in *_genes_reactome ; do cut -f 4 $i | sort -u > ${i}_ids ; done    # per mutation
cat *_genes_reactome_ids | sort -u > reactome_concat_all.txt              # for all mutations

# To create header in files (necessary to create matrix in R):
for i in *_genes_reactome_ids ; do sed -i -e '1i\"baseMean"' $i ; done
sed -i -e '1i\"baseMean"' reactome_concat_all.txt

# Intersect pathway IDs. To do this, see Count matrix R script.

# Create an Excel spreadsheet to organize data. Rescue original pathway names using vlookup function on Excel ( example: =VLOOKUP(A2;reactome_full_list!D:F;3;FALSE) ).

# Count number of times a certain pathway shows up in each causal mutation:
## Create files with all pathway ids per causal gene:
for i in *_genes_reactome ; do cut -f 2,4 $i | sort -u | cut -f 2 > ${i}_ids_all ; done    # per mutation

# Count occurrences:
## Run code once per gene:
for i in $(cat FMR1_genes_reactome_ids) ; do grep -x $i FMR1_genes_reactome_ids_all | wc -l >> FMR1_counts ; done
paste FMR1_genes_reactome_ids FMR1_counts > FMR1_reactome_counts

# Again, Rescue original pathway names using vlookup function on Excel ( example: =VLOOKUP(A2;reactome_full_list!D:F;3;FALSE) ).
