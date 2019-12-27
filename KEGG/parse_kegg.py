
import sys
import os
import argparse
import itertools

## other utility functions
import file_utils 
from convert_utils import *  ## map_namespace(), c(), convert()

## Biopython modules to interact with KEGG
# https://www.kegg.jp/kegg/xml/docs/
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG.KGML import KGML_pathway


def main(args):
	"""
	Main function.

	Parameters
	---------
	args: ArgumentParser object

	"""

	## get all the pathways listed for the species if --graph or --list is specified.
	if args.graph or args.list:
		pathways = REST.kegg_list('pathway',org=args.species)
	else: # pathways is simply the single --graph_single value.
		pathways = [args.graph_single]
		
	## If --list is specified, print all pathways to console.
	if args.list: 
		for p in pathways:
			print(p.strip())
		print()

	## if --graph or --graph_single is specified, 
	## get interactions and make graph for each pathway
	if args.graph or args.graph_single:

		## get namespace mapper. We will always map to SOME namespace.
		## map_namespace also takes care of filtering IDs if --filter is specified.
		kegg2id,id2kegg = map_namespace(args)

		num = 0
		for p in pathways:
			num+=1

			if args.graph:
				name = p.split()[0]
				short_name = name.split(':')[1]
			else: # --graph_single was specified
				name = args.graph_single
				short_name = args.graph_single

			## get KGML file if it's not already in the output directory.
			print('processing pathway #%d: %s' % (num,short_name))
			kgml_file = '%s/%s.kgml' % (args.outdir,short_name)
			if not os.path.isfile(kgml_file):
				print('KGML file does not exist. Pull it down from KEGG...')
				kgml = REST.kegg_get(name,option='kgml')
				file_utils.write_kgml(kgml_file,kgml)
			
			# parse the pathway.
			pathway = KGML_parser.read(open(kgml_file))

			print(' %s "%s": %d entries (incl. genes & groups) & %d relations' % (pathway.name,pathway.title,len(pathway.entries),len(pathway.relations)))

			# retain gene entries & map keggIDs to namespace.
			to_delete = set()
			pathway.gene_entries = {g:pathway.entries[g] for g in pathway.entries if pathway.entries[g].type == 'gene'}
			for node,entry in pathway.gene_entries.items():
				pathway.gene_entries[node].mapped_name = convert(pathway.gene_entries[node].name,kegg2id)
				if pathway.gene_entries[node].mapped_name == None:
					to_delete.add(node)
			print(' deleting %d gene entries with no mapping' % (len(to_delete)))
			for n in to_delete:
				del pathway.gene_entries[n]

			# retain gene groups & (a) add component IDs, (b) add component keggIDs, and (c) add map keggIDs to namespace.
			to_delete = set()
			pathway.gene_groups = {g:pathway.entries[g]for g in pathway.entries if pathway.entries[g].type == 'group'}
			for node,entry in pathway.gene_groups.items():
				pathway.gene_groups[node].ids = [component.id for component in entry.components]
				pathway.gene_groups[node].kegg_name = [pathway.gene_entries[i].name for i in pathway.gene_groups[node].ids if i in pathway.gene_entries]
				pathway.gene_groups[node].mapped_name = convert(pathway.gene_groups[node].kegg_name,kegg2id)
				if pathway.gene_groups[node].mapped_name == None:
					to_delete.add(node)
			print(' deleting %d gene groups with no mapping' % (len(to_delete)))
			for n in to_delete:
				del pathway.gene_groups[n]
			
			# pathway_ids are all the pathway IDs in genes & groups.
			pathway_ids = set(pathway.gene_entries.keys()).union(set(pathway.gene_groups.keys()))

			# retain relations that are among gene or group entries only and 
			pathway.gene_relations = [r for r in pathway.relations if r.entry1.id in pathway_ids and r.entry2.id in pathway_ids and not ignore(r)]

			print(' %d entries, %d groups, & %d relations after retaining genes & groups and removing ignored edges.' % (len(pathway.gene_entries),len(pathway.gene_groups),len(pathway.gene_relations)))

			# write entries, groups, and relations files (just for 'gene' and 'group' entities and relations)
			entries_file = '%s/%s-gene-entries.txt' % (args.outdir,short_name)
			file_utils.write_kgml_entries(entries_file,pathway)

			groups_file = '%s/%s-gene-groups.txt' % (args.outdir,short_name)
			file_utils.write_kgml_groups(groups_file,pathway)

			relations_file = '%s/%s-gene-relations.txt' % (args.outdir,short_name)
			file_utils.write_kgml_relations(relations_file,pathway)
			
			## generate graphs
			relation_counts = {'dir':0,'undir':0}

			# instead of writing edges directly, keep dictionaries that are keyed
			# by the edge identifiers.  Sometimes there are duplicate edges for various
			# reasons - this guarantees that we will only write unique edges to the file at the end.
			collapse_edges = {} # dictionary of collapsed edges
			expand_edges = {} # dictionary of expanded edges
			expanded_groups = set() # this will keep track of the groups that we have already expanded.
			for entry in pathway.gene_relations:

				## get node names, types, and whether the interaction is directed.
				n1,n2,t1,t2,is_directed = get_relation_entry_info(entry,pathway)

				if is_directed:
					relation_counts['dir']+=1
				else:
					relation_counts['undir']+=1

				## store collapsed edges
				add_to_dictionary(collapse_edges,(c(n1),c(n2),t1,t2),[e[0] for e in entry.subtypes])
				if not is_directed:
					add_to_dictionary(collapse_edges,(c(n2),c(n1),t2,t1),[e[0] for e in entry.subtypes])

				# expand edges
				expanded, expanded_groups = expand_entry_edges(n1,n2,t1,t2,c([e[0] for e in entry.subtypes]),expanded_groups)

				## store expanded edges
				for n1,n2,t in expanded:
					add_to_dictionary(expand_edges,(n1,n2),t)
					if not is_directed:
						add_to_dictionary(expand_edges,(n2,n1),t)

			print('Processed %d directed and %d undirected KEGG relations' % (relation_counts['dir'],relation_counts['undir']))

			## write edge files 
			collapse_file = '%s/%s-collapsed-edges.txt' % (args.outdir,short_name)
			expand_file = '%s/%s-expanded-edges.txt' % (args.outdir,short_name)
			file_utils.write_edge_files(collapse_file,collapse_edges,expand_file,expand_edges)

		print('Done making graph for each pathway.')

	return 

def add_to_dictionary(d,key,value):
	"""
	Utility function to add an entry to a dictionary that stores edges.

	Parameters
	--------------
	d: dict
	   dictionary to add key/value pair
	key: string
	   key
	value: string, list, or set
	   add all elements of value to the set keyed by the key in the dictionary.

	"""

	if key not in d: # create a <key,set()> pair
		d[key] = set()

	# if the value is a list or a set, add all elements to the set.
	# if the value is a string, just add the single element to the set.
	if type(value) == list or type(value) == set:
		d[key].update(set(value))
	else:
		d[key].add(value)

	return

def get_relation_entry_info(entry,pathway):
	"""
	Extracts the nodes and node types for a relation entry.

	Parameters
	-------------
	entry: Bio.KEGG.KGML.KGML_pathway.Entry object
	pathway: Bio.KEGG.KGML.KGML_pathway.Pathway object

	Returns
	-------------
	list
	   node1 mapped name(s)
	list
	   node2 mapped names(s)
	string
	   node1 entity type ('gene' or 'group')
	string
	   node2 entity type ('gene' or 'group')
	bool
	   True if the relation is directed; False otherwise

	"""
	is_directed = directed(entry)
	e1 = entry.entry1.id
	if e1 in pathway.gene_entries:
		n1  = pathway.gene_entries[e1].mapped_name
		t1 = pathway.gene_entries[e1].type
	else:
		n1 = pathway.gene_groups[e1].mapped_name
		t1 = pathway.gene_groups[e1].type
	e2 = entry.entry2.id
	if e2 in pathway.gene_entries:
		n2  = pathway.gene_entries[e2].mapped_name
		t2 = pathway.gene_entries[e2].type
	else:
		n2 = pathway.gene_groups[e2].mapped_name
		t2 = pathway.gene_groups[e2].type
	return n1,n2,t1,t2,is_directed

def expand_entry_edges(n1,n2,t1,t2,rel_type,expanded_groups):
	"""
	Take a collapsed edge and "expand" it by adding edges for certain
	pairs of elements.  Entities labeled as "groups" are protein complexes.

	Parameters
	------------
	n1: string, string with spaces, or list
	  IDs for the first node in the collapsed edge
	n2: string, string with spaces, or list
	  IDs for the second node in the collapsed edge
	t1: string ('gene' or 'group')
	  entity type of the first node
	t2: string ('gene' or 'group')
	  entity type of the second node
	rel_type: string
	  Relation type of the collapsed edge (e.g. 'activation')
	expanded_groups: set
	  Set of group entities that have already been expanded (no need to re-process them)

	Returns
	--------------
	set
	  expanded edges, represented as 3-tuples of (n1,n2,edge_type)
	  edge_type may be one of 'group_expansion','one_to_one_mapping:rel_type', or 'mult_mapping_expansion:rel_type'
	set
	  expanded_groups set

	"""

	expanded = set() # set of 3-tuples that will conain (n1,n2,edge_type)

	## if n1 is a group and we haven't expanded it yet, introduce all vs. all edges.
	if t1 == 'group' and c(n1) not in expanded_groups:
		for u1,u2 in itertools.permutations(n1,2): 
			expanded.add((u1,u2,'group_expansion'))
		expanded_groups.add(c(n1))

	## if n2 is a group and we haven't expanded it yet, introduce all vs. all edges.
	if t2 == 'group' and c(n2) not in expanded_groups:
		for u1,u2 in itertools.permutations(n2,2):
			expanded.add((u1,u2,'group_expansion'))
		expanded_groups.add(c(n2))

	## if n1 and n2 are single nodes, we have a one-to-one mapping of the collapsed edge.
	if len(n1) == 1 and len(n2) == 1:
		expanded.add((n1.pop(),n2.pop(),'one_to_one_mapping:%s' % (rel_type)))
	else: 
		# otherwise, the edges are expanded from a multiple mappings (many-to-one, one-to-many, or many-to-many).
		# for now, these are all considered as "multiple mappings".
		for u1,u2 in itertools.product(n1,n2):
			expanded.add((u1,u2,'mult_mapping_expansion:%s' % (rel_type)))

	return expanded,expanded_groups


def ignore(entry):
	"""
	Checks if the relation entry should be ignored.

	Relation entries have an interaction type and interaction subtypes.
	Ignore GErel (gene expression ineraction) types and some interaction
	subtypes that denote vague or expression-based relationships.

	Parameters
	-------------
	entry: Bio.KEGG.KGML.KGML_pathway.Entry object

	Returns
	-------------
	bool
	   True if the relation entry should be ignored, False otherwise.

	"""

	# ignore GErel (gene expression interaction, TF -> Target Gene).
	if entry.type == 'GErel':
		return True

	subtypes = set(s[0] for s in entry.subtypes)

	# ignore relations with no subtypes.
	if len(subtypes) == 0: 
		return True

	# ignore 'state-change', 'missing-interaction', or 'expression' subtypes
	if 'state change' in subtypes or 'missing interaction' in subtypes or 'expression' in subtypes:
		return True

	# if we get this far, edge is not ignored. 
	return False

def directed(entry):
	"""
	Checks if the relation entry is directed according to user-specified rules.  

	Relations have multiple subtypes; this function processes the subtypes in a 
	specific order. 

	Parameters
	-----------
	entry: Bio.KEGG.KGML.KGML_pathway.Entry object

	Returns
	-----------
	bool
	  True if the relation should be treated as a directed edge, False otherwise.
	  Note that the function exits with an error if direction cannot be determined.

	"""

	## get the subtype names of this entry.
	subtypes = set([s[0] for s in entry.subtypes])
	
	# If edge has 'activation' or 'inhibition' edge type, let this set the direction.
	# activation/inhibition: 'positive and negative effects which may be associated 
	# with molecular information below'
	if 'activation' in subtypes or 'inhibition' in subtypes:
		return True
		
	# If the edge is a molecular event (phos,dephos,glyco,ubiq,or methyl),
	# it is a directed edge.  Doesn't matter if it also has an undirected subtype.
	if 'phosphorylation' in subtypes or 'dephosphorylation' in subtypes or \
			'glycosylation' in subtypes or 'ubiquitination' in subtypes or \
			'methylation' in subtypes:
		return True
		
	# If edge is 'indirect effect', then it is directed.
	elif 'indirect effect' in subtypes:
		return True

	# Compounds are tricky.  For now, add directed edges for compound.  
	# compound: "shared with two successive reactions (ECrel) or intermediate of two 
	# interacting proteins (PPrel)"
	elif 'compound' in subtypes:
		return True

	# Otherwise, the edge is an undirected edge.
	# binding/association, dissociation
	# group-entry: KEGG entry that is labeled as 'group', which is interpreted as a complex.
	# group-entry is the ONLY label that is NOT from the KEGG manual (we add it when making edges
	# to represent the complex)
	elif 'binding/association' in subtypes or 'dissociation' in subtypes \
			or 'group' in subtypes:
		return False
		
	else: 
		print('ERROR! Edge direction cannot be established with subtypes %s' % (c(subtypes)))
		print(entry)
		sys.exit()
	
	return

def parse_arguments():
	""" 
	Argument Parser for parse_kegg.py.

	Returns
	-----------
	ArgumentParser object

	"""
	parser = argparse.ArgumentParser('KEGG Pathway Processor. At least one of --list, --graph, or --graph_single must be specified.')
	parser.add_argument('--list',action='store_true',help='list pathways to stdout.')
	parser.add_argument('--graph',action='store_true',help='make graph for all pathways from the specified species.')
	parser.add_argument('--graph_single',help='make graph of a single pathway. Pass in the pathway identifier (e.g. hsa04310).')
	parser.add_argument('-s','--species',default='hsa',help='species/taxon identifier. Default is hsa.')
	parser.add_argument('-c','--convert',default='uniprot',help='convert kegg id to this case insensitive id/namespace (ncbi-geneid | uniprot). Default is uniprot')
	parser.add_argument('-f','--filter',help='filter converted IDs by single-column file of ids. Only IDs that appear in this file will be used.')
	parser.add_argument('-o','--outdir',help='outfile directory.')
	args = parser.parse_args()

	## one of --list, -graph, or --graph_single must be specified.
	if not(args.list or args.graph or args.graph_single):
		sys.exit('ERROR: --list, --graph, or --graph_single must be specified. Exiting.')
	
	## only one of --graph and --graph_single can be specified.
	if args.graph and args.graph_single:
		sys.exit('ERROR: --graph and --graph_single cannot both be specified. Exiting.')
	
	## outdir must be specified if we are generating a graph
	if (args.graph or args.graph_single) and not args.outdir:
		sys.exit('ERROR: --species and --outdir must be specified to make graphs. Exiting.')

	## if a filter file is specified, it must exist.
	if args.filter and not os.path.isfile(args.filter):
		sys.exit('ERROR: namespace file filter "%s" does not exist. Exiting.' % (args.filter))

	## make output directory if it does not exist.
	if (args.graph or args.graph_single) and not os.path.isdir(args.outdir):
		print('making output directory %s...' % (args.outdir))
		os.makedirs(args.outdir)
	
	return args

if __name__ == '__main__':
	args = parse_arguments()
	main(args)