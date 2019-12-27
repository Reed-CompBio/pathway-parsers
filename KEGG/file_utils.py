## IO utilities for reading/writing KGML and related files.
import sys
from convert_utils import *

def write_kgml(kgml_file,kgml):
	"""
	write KGML file to output file.

	Parameters
	--------------
	kgml_file: str
	   KGML output file name
	kgml: IOTextWrapper 
	   File handle
	"""

	out = open(kgml_file,'w')
	for line in kgml:
		out.write(line)
	out.close()

	print(' wrote to %s' % (kgml_file))
	return

def write_kgml_entries(entries_file,pathway):
	"""
	Writes KGML gene entries (id, mapped_names, kegg_names). Pathway object must 
	have been modified to add gene_entries and gene_entries.mapped_name.  

	Parameters
	-------------------
	entries_file: string
	   Entries output file name 
	pathway: Bio.KEGG.KGML.KGML_pathway object

	"""
	
	out = open(entries_file,'w')
	out.write('#id\tmapped_name\tkegg_name\n')
	for node,entry in pathway.gene_entries.items():
		out.write('%s\t%s\t%s\n' % (c(node),c(entry.mapped_name),c(entry.name)))
	out.close()

	print(' wrote to %s' % (entries_file))
	return

def write_kgml_groups(groups_file,pathway):
	"""
	Writes KGML gene groups (id, component_ids, mapped_names, kegg_names). Pathway 
	object must have been modified to add gene_groups, gene_groups.ids, 
	gene_groups.kegg_name, and gene_groups.mapped_name.  

	Parameters
	-------------------
	groups_file: string
	   Groups output file name 
	pathway: Bio.KEGG.KGML.KGML_pathway object

	"""
	out = open(groups_file,'w')
	out.write('#id\tcomponent_ids\tmapped_names\tkegg_names\n')
	for node,entry in pathway.gene_groups.items():
		out.write('%s\t%s\t%s\t%s\n' % (c(node),c(entry.ids),c(entry.mapped_name),c(entry.kegg_name)))
	out.close()

	print(' wrote to %s' % (groups_file))
	return

def write_kgml_relations(relations_file,pathway):
	"""
	Writes KGML relations (id1, id2, relation_type, relation_subtype). Pathway 
	object must have been modified to add gene_relations.  

	Parameters
	-------------------
	relations_file: string
	   Relations output file name 
	pathway: Bio.KEGG.KGML.KGML_pathway object

	"""

	out = open(relations_file,'w')
	out.write('#id1\tid2\ttype\tsubtype\n')
	for entry in pathway.gene_relations: ## NOTE: this is the gene_relations, different than ALL relations.
		out.write('%s\t%s\t%s\t%s\n' % (entry.entry1.id,entry.entry2.id,c(entry.type),c(['%s' % s[0] for s in entry.subtypes])))
	out.close()

	print(' wrote to %s' % (relations_file))
	return

def write_edge_files(collapse_file,collapse_edges,expand_file,expand_edges):
	"""
	Write two sets of edge files: one of "collapsed" edges, and one of "expanded" edges.

	Parameters
	---------------
	collapse_file: string
	   Output file of collapsed graph
	collapse_edges: dict
	   Dictionary of (edge,relation_type) key/value pairs
	expand_file: string
	   Output file of expanded graph
	expand_edges: dict
	   Dictionary of (edge,edge_types) key/value pairs. Edge_types include relation_types from 
	   collapsed versions, as well as an indication of why the edge was expanded.

	"""
	## write collapsed file 
	out_collapse = open(collapse_file,'w')
	out_collapse.write('#node1\tnode2\tnode1type\tnode2type\trelation_type\n')
	for (n1,n2,t1,t2),relation_types in collapse_edges.items():
		# nodes are already in text form at this point.
		out_collapse.write('%s\t%s\t%s\t%s\t%s\n' % (n1,n2,t1,t2,c(relation_types))) 
	out_collapse.close()

	print(' wrote %d (collapsed) edges to %s' % (len(collapse_edges),collapse_file))

	## write expanded file
	out_expand = open(expand_file,'w')
	out_expand.write('#node1\tnode2\tedge_expansion:relation_type\n')
	for (n1,n2) in sorted(expand_edges.keys()):
		edge_types = expand_edges[(n1,n2)]
		# nodes shouldn't need to be collapsed with the c() function - they are singletons!
		out_expand.write('%s\t%s\t%s\n' % (n1,n2,c(edge_types)))
	out_expand.close()
	
	print(' wrote %d (expanded) edges to %s' % (len(expand_edges),expand_file))
	return