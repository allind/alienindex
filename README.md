# alienindex

Calculate alien index as described in https://doi.org/10.1093/molbev/msw073

Input files are blastp with taxonomy information against a large database (recommended: nr_cluster_seq), blastp against self (to know what the maximum possible bitscore a protein can have is), a 'recipient' taxa (the group into which HGT is being tested to have occured. Most accurate if this is a higher taxonomic level like 'eukaryotes' or 'fungi'.) and a 'self' taxa (NCBI taxid for the organism being tested, so that blast hits against this group are not considerd possible in-group matches. This is usually species-level).

Optional inputs are a location for an ete3 taxonomy database location (if you're using different versions of NCBI databases and interacting with that through ete3, this is necsesary), and a file with a list of NCBI taxids to exclude as hits (e.g. if there is a genome assembly in a potential donor group that has been contaminated and has hits to your taxon of interest, this artificially inflates the alien index and you want to remove these).

I recommend using diamond for blast search as it's faster than ncbi blastp. The command that should be run for diamond blastp is `diamond blastp --query [qfile] --db [db] --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen slen qlen nident positive staxids qcovhsp --threads [threads] -b12 -c1 --out [output]`.
