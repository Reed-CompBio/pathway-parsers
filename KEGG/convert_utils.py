## IO utilities for converting, mapping, and compressing identifiers.
from Bio.KEGG import REST
import sys

def map_namespace(args):
	"""
	Uses BioPython's KEGG REST api to convert kegg IDs to a different namespace.

	Parameters
	-------------------
	args: ArgumentParser object
	  Contains args.convert (namespace to conver to) and args.filter (file to filter identifiers).

	Returns
	---------------------------
	dict
	   dictionary of kegg IDs to namespace IDs (treated as a set, in case one keggID maps to multiple namespace IDs)
	dict
	   dictionary of namespace IDs to kegg IDs (treated as a set, in case one namespace ID maps to multiple keggIDs)

	"""

	kegg2id = {} # dictionary of kegg IDs to namespace IDs
	id2kegg = {} # dictionary of namespace IDs to kegg IDs

	## convert all the keggIDs to the args.convert namespace for the args.species species
	response = REST.kegg_conv(args.convert,args.species)
	
	## if there is a filter file (args.filter), read the file as a single column and store
	## identifiers as a set. Any identifier NOT in this file will subseqently be ignored.
	to_filter = None
	if args.filter:
		to_filter = set()
		with open(args.filter) as fin:
			for line in fin:
				to_filter.add(line.strip())
		print('retaining only ids that are in the filter file (%d total)' % (len(to_filter)))

	## for every keggID - to - namespaceID, add it to the dictionaries.
	for e in response:

		## get the keggID and the namespaceID from the response string
		row = e.strip().split()
		kegg = row[0]
		new_id = row[1].split(':')[1]

		# if filtered file was present, only continue
		# if the mapped ID is in the filter file.
		if to_filter and new_id not in to_filter:
			continue
		
		## update the dictionaries
		if kegg not in kegg2id:
			kegg2id[kegg] = set()
		kegg2id[kegg].add(new_id)
		
		if new_id not in id2kegg:
			id2kegg[new_id] = set()
		id2kegg[new_id].add(kegg)
		
	return kegg2id, id2kegg

def convert(kegg,kegg2id):
	"""
	Convert the keggID(s) to the namespace IDs stored in the kegg2id dictionary.

	Parameters
	------------
	kegg: string, string with spaces, or list/set
	   If a single keggID, this is a string. This can also be a string with multiple
	   keggIDs separated by spaces, or a list or set of keggIDs.
	kegg2id: dict
	   A kegg-to-namespace mapping dictionary (provided by a call to map_namespace() above)

	Returns
	-------------
	set
	   A set of mapped identifiers (which may be a different size than the original keggID input).  
	   If there are no mapped IDs then the function returns None. This is handled in other functions.
	"""

	if type(kegg) == str and ' ' in kegg: #if there are spaces in the string, assume something like 'hsa:1234 hsa:5432 hsa:9432'
		kegg_ids = kegg.split()
	elif type(kegg) == str: # if there are no spaces in the string, assume something like 'hsa:1234'
		kegg_ids = [kegg]
	elif type(kegg) == list or type(kegg) == set: # if this is a list/set, assume something like ['hsa:1234','hsa:5423','hsa:9432']
		kegg_ids = kegg

	## take each id separately and map to the namespace
	ids = set()
	for kegg in kegg_ids:
		if kegg in kegg2id: # if there was a filter file, this might not be true.
			ids.update(kegg2id[kegg])

	## return None if there are no ids.
	if len(ids) == 0:
		return None

	return ids

def c(s,delim='|'):
	"""
	Small function to "compess" an entry with the specified delimiter.

	Parameters:
	------------------
	s: string, int, float, list, or set
	   This input can be anything. Spaces are interpreted as multiples. (e.g. 'hsa:1234 hsa:5432 hsa:9432', 'hsa:1234',
	   ['hsa:1234','hsa:5423','hsa:9432'], and 1234 are all appropriate)
	delim: str
	   Delimiter to use (default = '|')

	Returns
	------------------
	string
	  String representation of compressed fields, using the delimiter specified. Function exits with an Error
	  if the value cannot be compressed.
	"""

	if type(s) == str: 
		return s.replace(' ',delim)
	elif type(s) == int or type(s) == float: 
		return str(s)
	elif type(s) == list or type(s) == set: 
		return delim.join(sorted([str(a) for a in s]))
	else:
		print("Error: cannot compress ",s)
		sys.exit()
		