# KEGG Database Pathway Parser

This parser converts [KEGG KGML files](https://www.kegg.jp/kegg/xml/docs/) into two versions of pathways: one that retains family and complex information, and one that "expands" the protein families and complexes into individual proteins.  It was heavily adapted from the parser developed in TM Murali's group at Virginia Tech, which was authored by myself, Allison Tegge, and Richard Rodrigues at Virginia Tech (around 2013 or so).

This parser is written in Python3 using the [BioPython module](https://biopython.org/) for accessing and traversing [KEGG KGML objects](https://www.kegg.jp/kegg/xml/docs/) through [KEGG's REST API](https://www.kegg.jp/kegg/rest/keggapi.html).  We currently use it to convert KEGG pathways into UniProtKB identifiers; entities that are not mapped to the appropriate namespace are ignored.  Note that the edges here are not necessarily represented in the large interactomes that our group uses, since the most recent parsing might have introduced new interactions that are not part of the original [PathLinker](https://github.com/Murali-group/PathLinker) interactomes.  An example of a recently developed interactome is included in the [Localized PathLinker (LocPL)](https://github.com/annaritz/localized-pathlinker) repository.

The main function is located within `parse_kegg.py`. Running `python3 parse_kegg.py -h` will print usage information:

```
usage: KEGG Pathway Processor. At least one of --list, --graph, or --graph_single must be specified.
       [-h] [--list] [--graph] [--graph_single GRAPH_SINGLE] [-s SPECIES]
       [-c CONVERT] [-f FILTER] [-o OUTDIR]

optional arguments:
  -h, --help            show this help message and exit
  --list                list pathways to stdout.
  --graph               make graph for all pathways from the specified
                        species.
  --graph_single GRAPH_SINGLE
                        make graph of a single pathway. Pass in the pathway
                        identifier (e.g. hsa04310).
  -s SPECIES, --species SPECIES
                        species/taxon identifier. Default is hsa.
  -c CONVERT, --convert CONVERT
                        convert kegg id to this case insensitive id/namespace
                        (ncbi-geneid | ncbi-proteinid | uniprot). Default is
                        uniprot
  -f FILTER, --filter FILTER
                        filter converted IDs by single-column file of ids.
                        Only IDs that appear in this file will be used.
  -o OUTDIR, --outdir OUTDIR
                        outfile directory.
```

## Requirements
* Python3 (most recently tested with `python 3.8.2`)
* BioPython (most recently tested with `biopython==1.77`)

## Options and Output Files

There are three main functions of this script: (a) list all files (`--list`), (b) parse a single pathway (`--graph_single`), and (c) parse all pathways for a species (`--graph`).  If one or more pathways are parsed, files are placed in an output directory specified with `-o` or `--outdir`. If no directory exists, then one is automatically created.

* `pathway.kgml`: KGML file from KEGG.  If this file exists, the program reads from the file instead of queries KEGG through the REST API.
* `pathway-gene-entries.txt`: tab-delimited file of gene entries in the pathway.
* `pathway-gene-groups.txt`: tab-delimited file of gene groups (complexes) in the pathway.
* `pathway-gene-relations.txt`: tab-delimited file of entity relations (interactions) in the pathway.
* `pathway-collapsed-edges.txt`: graph with "collapsed" edges.
* `pathway-expanded-edges.txt`: graph with "expanded" edges.

See Parsing Details for more information about the intermediate and final output files.

## Quick Start

### Listing KEGG Pathways

List all human pathways:
```
python3 parse_kegg.py --list
```

List all yeast pathways:
```
python3 parse_kegg.py --list -s sce
```
### Parse a single KEGG Pathway

Parse the human Wnt signaling pathway, converted to UniProtKB IDs (default). Store files in `output/` directory:
```
python3 parse_kegg.py --graph_single hsa04310 -o output/
```

Parse the human Wnt signaling pathway with reviewed UniProtKB IDs (passed as a filter file).  See the Filter File section below for more details.
```
python3 parse_kegg.py --graph_single hsa04310 -o output/ -f uniprot-swissprot-ids.txt
```

Parse the human Wnt signaling pathway with NCBI gene IDs. *TODO:* This is a little buggy.
```
python3 parse_kegg.py --graph_single hsa04310 -o output -c ncbi-geneid
```

### Parse all KEGG Pathways
Parse all human signaling pathways, converted to UniProtKB IDs (default). Store files in `output/` directory:
```
python3 parse_kegg.py --graph -o output/
```

Parse all yeast signaling pathways, converted to UniProtKB IDs (default). Store files in `output/` directory:
```
python3 parse_kegg.py --graph -s cse -o output/
```
## Filter File

I downloaded the filter file of UniProtKB reviewed proteins (SwissProt) from the [UniProt Database website](https://www.uniprot.org/).  
![download all human reviewed proteins](uniprot-reviewed.png)

The filter file passed with the `-f` or `--filter` option is a single-column file that contains all *allowed* namespace identifiers.  Passing this filter file in will require that all UniProtKB identifiers are reviewed.

## Parsing Details

Here, I use the parsed, reviewed-filtered [Wnt](https://www.genome.jp/kegg/pathway/hsa/hsa04310.html) to illustrate parsing details.  The files can be gnerated with the automatically

```
python3 parse_kegg.py --graph_single hsa04310 -o output/ -f uniprot-swissprot-ids.txt
```

### Gene Entries

In some cases, a `gene-entry` is a collection of proteins (e.g. a protein family) that is considered as a single node in a KEGG pathway.  Using , gene entry 49 is `P36402|Q9HCS4|Q9NQB0|Q9UJU2`,
which represents members of [TCF/LEF transcription factors](https://www.genome.jp/dbget-bin/www_bget?hsa:51176+hsa:6932+hsa:6934+hsa:83439).   Similarly, gene entry 50 is `Q09472|Q92793` represents [CREB-binding proteins](https://www.genome.jp/dbget-bin/www_bget?hsa:1387+hsa:2033).

## Gene Groups

Gene groups represent protein complexes from the `gene-entry` list.  A `gene-group` element consists of multiple `gene-entry` identifiers.  An example from `hsa05418` is

```
#id	component_ids	mapped_names	kegg_names
228	68|70|82	P33151|P35222|P35968	hsa:1003|hsa:1499|hsa:3791
```

### From Relations to Edges

The `gene-relations` file relies on gene entry and gene group IDs.  In our example above, there is an edge from gene entry 49 to gene entry 50:

```
#id1	id2	type	subtype
50	49	PPrel	binding/association
```

In the **collapsed** edges file, the gene entry and gene group IDs are considered the nodes and the concatenated UniProt IDs are passed as node names.   Thus, we would have these bidirected edges in `output/hsa04310-collapsed-edges.txt`.

```
#node1	node2	node1type	node2type	relation_type
Q09472|Q92793	P36402|Q9HCS4|Q9NQB0|Q9UJU2	gene	gene	binding/association
P36402|Q9HCS4|Q9NQB0|Q9UJU2	Q09472|Q92793	gene	gene	binding/association
```

However, we often want to have the UniProt IDs as the nodes. In this case, in the **expanded** edges file we create edges for all (node1,node2) UniProt IDs.  In `output/hsa04310-expanded-edges.txt`, the single gene relation becomes

```
#node1	node2	edge_expansion:relation_type
P36402	Q09472	mult_mapping_expansion:binding/association
P36402	Q92793	mult_mapping_expansion:binding/association
Q09472	P36402	mult_mapping_expansion:binding/association
Q09472	Q9HCS4	mult_mapping_expansion:binding/association
Q09472	Q9NQB0	mult_mapping_expansion:binding/association
Q09472	Q9UJU2	mult_mapping_expansion:binding/association
...
```

When there are `gene-group` entries for a pathway, these are expanded to include edges among all pairs of UniProtIDs in the complex.  While the Wnt pathway doesn't have any group entries, the `hsa05418` group entry above appears as the following in `output/hsa05418-expanded-edges.txt`:

```
#node1	node2	edge_expansion:relation_type
P33151	P35222	group_expansion
P33151	P35968	group_expansion
...
```
