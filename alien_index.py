#! /usr/bin/env python
#usage: script.py [bls vs nr] [bls vs self] [recipient taxid]
#uses diamond
#diamond cmd is diamond blastp --query [qfile] --db [db] --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen slen qlen nident positive staxids qcovhsp --threads [threads] -b12 -c1 --out [output]
#for bls vs self, same command except db is the same as query
#skip and recip file is structured: species_name skip_taxid recipient_name recipient_taxid 
#diamond NR

from ete3 import NCBITaxa
import sys, os, re


def parse_blast(bls, skip, recipient, ncbi):
	max_bitscore = 0
	
	bitscores = []
	
	in_group_hits = []
	out_group_hits = []
	for line in bls:
		line = line.strip("\n")
		qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,\
		sstart, send, evalue, bitscore, saccver, slen, qlen, nident, positive,\
		staxids, qcovs = line.split('\t')
		lineage_taxids = []
		
		for taxid in staxids.split(';'):
			if taxid == "N/A" or len(taxid) == 0:
				#skip this one, get out of the bls iterator
				lineage_taxids.append("N/A")
				break
			for l in ncbi.get_lineage(int(taxid)):
				if l not in lineage_taxids:
					lineage_taxids.append(l)
					
		if float(bitscore) > max_bitscore:
			max_bitscore = float(bitscore)

		if skip not in lineage_taxids and "N/A" not in lineage_taxids:
		    #make stuff to save
		    bitscores.append([float(bitscore), staxids.split(';')[0]])
		    info = [sseqid, float(bitscore), float(evalue), staxids.split(';')[0]]
		    if recipient in lineage_taxids:
		        in_group_hits.append(info)
		    else:
		        out_group_hits.append(info)

	sorted_bitscores = sorted(bitscores, key = lambda x: float(x[0]), reverse=True)
	if len(sorted_bitscores) < 5:
		top_bitscores = sorted_bitscores
	else:
		top_bitscores = sorted_bitscores[0:4]
	
	classes = []
	for hit in top_bitscores:
		lineage = ncbi.get_lineage(hit[1])
		ranks = ncbi.get_rank(lineage)
		#rank_ids = dict((v,k) for k,v in ranks.iteritems())

		rank_ids = {ranks[v]:v for v in ranks}

		#if "class" in rank_ids:
		#	classes.append(list(ncbi.get_taxid_translator([rank_ids["class"]]).values())[0])

			#classes.append(ncbi.get_taxid_translator([rank_ids['class']]))
		#elif "order" in rank_ids:
		#	classes.append(list(ncbi.get_taxid_translator([rank_ids["order"]]).values())[0])
			#classes.append(ncbi.get_taxid_translator(rank_ids['order']))
		if "phylum" in rank_ids:
			classes.append(list(ncbi.get_taxid_translator([rank_ids["phylum"]]).values())[0])
			#classes.append(ncbi.get_taxid_translator(rank_ids['phylum']))
		elif "order" in rank_ids:
			classes.append(list(ncbi.get_taxid_translator([rank_ids["order"]]).values())[0])
		else:
			classes.append("No phlyum or order in NCBI taxid")
				
	classes = set(classes)
	return float(max_bitscore), classes, in_group_hits, out_group_hits
	
	#returns: float(max_bitscore), classes of 4 best bitscore hits, all in_group_hits, all out_group_hits

def calculate_alien_info(max_bitscore, in_lineage_hits, out_lineage_hits, ncbi):
	best_inlineage = ["", 0, float("inf")]
	best_outlineage = ["", 0, float("inf")]
	for hit in in_lineage_hits:
		if hit[1] > best_inlineage[1]:
			best_inlineage = hit
	
	for hit in out_lineage_hits:
		if hit[1] > best_outlineage[1]:
			best_outlineage = hit
	
	alien = float(float(best_outlineage[1]/max_bitscore) - float(best_inlineage[1]/max_bitscore))
	if len(in_lineage_hits) == 0:
		best_inlineage = "No in group hits"
		best_inlineage_species = "-"
		in_class_name = '-'
	else:
		best_inlineage_species = list(ncbi.get_taxid_translator([best_inlineage[3]]).values())[0]
		in_lineage = ncbi.get_lineage(best_inlineage[3])
		in_ranks = ncbi.get_rank(in_lineage)
		#in_rank_ids = dict((v,k) for k,v in in_ranks.iteritems())
		in_rank_ids = {in_ranks[v]: v for v in in_ranks}
		#if "class" in in_rank_ids:
		#	in_class_name = list(ncbi.get_taxid_translator([in_rank_ids["class"]]).values())[0]
		#elif "order" in out_rank_ids:
		#	in_class_name = list(ncbi.get_taxid_translator([in_rank_ids["order"]]).values())[0]
		if "phylum" in in_rank_ids:
			in_class_name= list(ncbi.get_taxid_translator([in_rank_ids["phylum"]]).values())[0]
		elif "order" in in_rank_ids:
			in_class_name= list(ncbi.get_taxid_translator([in_rank_ids["order"]]).values())[0]
		else:
			in_class_name = "No phlyum in NCBI taxid"

	if len(out_lineage_hits) == 0:
		best_outlineage = "No non-group hits"
		best_outlineage_species = "-"
		out_class_name = "-"
	else:

		best_outlineage_species = list(ncbi.get_taxid_translator([best_outlineage[3]]).values())[0]

		out_lineage = ncbi.get_lineage(best_outlineage[3])
		out_ranks = ncbi.get_rank(out_lineage)
		#out_rank_ids = dict((v,k) for k,v in out_ranks.iteritems())
		out_rank_ids = {out_ranks[v]: v for v in out_ranks}
		#if "class" in out_rank_ids:
		#	out_class_name = list(ncbi.get_taxid_translator([out_rank_ids["class"]]).values())[0]
		#elif "order" in out_rank_ids:
		#	out_class_name = list(ncbi.get_taxid_translator([out_rank_ids["order"]]).values())[0]
		if "phylum" in out_rank_ids:
			out_class_name=list(ncbi.get_taxid_translator([out_rank_ids["phylum"]]).values())[0]
		elif "order" in out_rank_ids:
			out_class_name=list(ncbi.get_taxid_translator([out_rank_ids["order"]]).values())[0]
		else:
			#just grab the species, whatever
			out_class_name = "No phlyum in NCBI taxid"

	return alien, best_inlineage, best_inlineage_species, in_class_name, best_outlineage, best_outlineage_species, out_class_name
def main(argv):
	
	ncbi = NCBITaxa("/pollard/data/projects/alind/eukdetect/ncbi_met_and_arch/busco/alien/nr_dl/taxdump-current/taxa.sqlite")

	#
	#test file just has one sequence
	
	#structure:
	#qseqid(0) sseqid(1) pident(2) length(3) mismatch(4) gapopen(5) qstart(6) qend(7) sstart(8) send(9) evalue(10) bitscore(11) saccver(12) slen(13) qlen(14) nident(15) positive(16) staxids(17) qcovs(18)
	#
	
	skip = 12967 #Blastocystis
	recipient = 2759 #Euks
	recipient_name = "Eukaryotes"
	#skip_and_recip = {line.split('\t')[0]: line.strip('\n').split('\t')[1:] for line in open(sys.argv[2])} #open("alien_above_0_species_skip_recipient.txt")}

	in_lineage_hits = []
	out_hits = []
	
	self_bitscores = {}

	#parse self blast:
	for line in open(sys.argv[2]):
		line = line.strip('\n')
		q = line.split('\t')[0]
		s = line.split('\t')[1]
		if q == s:
			self_bitscores[q] = float(line.split('\t')[-1])


	#parse blast
	parse = []
	currgene = ""
	dest = open(sys.argv[3], 'w')
	bls_linecount = len(open(sys.argv[1]).readlines())
	#print("here")
	linecounter = 0
	dest.write("query\talien_score\trecipient_name\trecipient_taxid\tbest_ingroup_hit\tbest_ingroup_species\tbest_ingroup_class\tbest_ingroup_bitscore\tbest_outgroup_hit\tbest_outgroup_species\tbest_outgroup_class\tbest_outgroup_bitscore\tclasses_of_top_4_hits\tself_bitscore\n")
	for line in open(sys.argv[1]):
		linecounter += 1
		#print(line)
		line = line.strip('\n')
		gene = line.split('\t')[0]
		
		if currgene == "":
			currgene = gene
		elif currgene == gene and linecounter != bls_linecount:
			parse.append(line)

		else:
			print("parsing:", currgene)

			#species = '-'.join(re.split('-\d*at\d*-', currgene)[0].split('-')[1:])
			species = "BT1" #hardcode for now
			#skip = int(skip_and_recip[species][0]) #donor
			#recipient = int(skip_and_recip[species][2]) #ingroup
			#recipient_name = skip_and_recip[species][1]

			max_bitscore = self_bitscores[currgene]

			max_blast_bitscore, best_lineages, in_lineage_hits, out_hits = parse_blast(parse, skip, recipient, ncbi)
			#iterate over in hits and out hits
			alien, best_inlineage, best_inlineage_species, in_class_name, \
			best_outlineage, best_outlineage_species, out_class_name = calculate_alien_info(max_bitscore, in_lineage_hits, out_hits, ncbi)
			
			dest.write(currgene + '\t' + str(alien) + '\t' + recipient_name + '\t' + str(recipient) + '\t' + best_inlineage[0] + \
			"\t" + best_inlineage_species + '\t' + in_class_name + '\t' + \
			str(best_inlineage[1])+ '\t' + best_outlineage[0] + '\t' + \
			best_outlineage_species+ '\t' + out_class_name + '\t' + \
			str(best_outlineage[1]) + '\t' + ",".join(best_lineages) + '\t' + str(max_bitscore) + '\n')
			
			#reset
			parse = []
			currgene = gene
			parse.append(line)

	print("parsing:", currgene)
	
	#once more for the last guy
	alien, best_inlineage, best_inlineage_species, in_class_name, \
	best_outlineage, best_outlineage_species, out_class_name = calculate_alien_info(max_bitscore, in_lineage_hits, out_hits, ncbi)
			
	dest.write(currgene + '\t' + str(alien) + '\t' + recipient_name + '\t' + str(recipient)  + '\t' + best_inlineage[0] + \
	"\t" + best_inlineage_species + '\t' + in_class_name + '\t' + \
	str(best_inlineage[1])+ '\t' + best_outlineage[0] + '\t' + \
	best_outlineage_species+ '\t' + out_class_name + '\t' + \
	str(best_outlineage[1]) + '\t' + ",".join(best_lineages) +'\t' + str(max_bitscore) +  '\n')

if __name__ == "__main__":
  main(sys.argv)
