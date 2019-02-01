# ISF_add_ons

Listed below all the features that are alrdeay or will
be impleted in this repository:

## ISF Batch run
	
-> Run ISF in batch to aggregates sequences to as many
input family of genes (seeds)

## ISF to MultiTwin

-> Given a list of family of genes aggregated by ISF,
parse the fasta file that comprise all the sequences 
of each family with the aim to retrieve all species-sequence
relationships of families.
This allows to associate to each family name a list of
representative species.

-> Create a tabular file whose primary key is (family_name, species_name)
and that can be taken as input by MultiTwin to construct bipartite graphs.

## ISF to phylogenetic tree

-> Given a sequence similarity network (SSN) built by ISF, boostrap the SSN,
draw as many Minimal Spanning Trees (MST) as there are iterations, learn
a consensus MST and turn it into a phylogenetic tree.

## ISF to functional annotations

-> Given a list of family of genes, run rpsblast with each family as query
against COG in order to assign COG categories to each sequence of each family.
Ultimately, COG category enrichment statistics are eventually computed for each family of genes.

-> Take COG enrichment statistics to cluster families of genes according to
their functional annotations. Clustering is agglomerative and appropriate
number of cluster is guided through consensus clustering. A shiny heatmap is 
eventually plotted to summarize all these information.

## ISF to diversity visualisation

-> Given a list of family of genes and a reference phylogenetic tree 
(the tree of life), anchor each representative species of a family (step 2) to
a leaf of the tree. 

-> Collapse the comprehensive tree to a lower-level tree, optimal for visualization.
(leaves that were species would be collapsed into leaves that are, for instance, phylums).
For each collapsed group of leaves, compute the number of species that have been anchored
to it for each family of genes. Turn systematically this information to a histogram that
summarize the distribution of each family of genes along the lower-level tree. Scale the 
histrogram independantly for each family (maximum value becomes 1, minimum 0). 
Ultimately, use itol, to represent on the left, the annotated tree of life, and on the
right, the corresponding distribution of each family of genes along the leaves of the tree.

## ISF to putative age

-> Given the reference tree, and the species anchored to its leaves (step 5), infer the evolutionary
history using the different methodologies implemented by the software "count". This allows
to compute several putative age for each family of genes.  

## ISF with HMM instead of blast

-> For each sown family of genes, do a multiple sequence alignment, learn its HMM profile.
Retrieve sequences based on this HMM profile. Include retrieved sequences to multiple alignment,
update the HMM profile and repeat the process until no sequences are retrieved anymore.
