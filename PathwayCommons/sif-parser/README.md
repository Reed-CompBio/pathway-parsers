This parser parses Pathway Commons SIF files.   All files downloaded from [PathwayCommons v11](https://www.pathwaycommons.org/archives/PC2/v11/).

## NetPath

All pathways are used for the graphlet project.

```
python3 parse_pc.py -i PathwayCommons11.netpath.hgnc.txt -o netpath
```

## KEGG

NO pathways are used for the graphlet project. _NOTE:_ PathwayCommons includes mainly metabolic pathways. Use the KEGG parser instead.
```
python3 parse_pc.py -i PathwayCommons11.kegg.hgnc.txt -o kegg -d KEGG
```

## NCI-PID

A subset of pathways are used for the graphlet project.

```
python3 parse_pc.py -i PathwayCommons11.pid.hgnc.txt -o pid
```

After generating the pathway files for all pathways, make a hand-curated list of signaling pathways to use as `pid_pathways.txt`.  Then,

```
mkdir pid_subset
cat pid_pathways.txt | awk '{print "cp pid/"$1" pid_subset/"$1}' | bash
ls pid_subset | wc -l
```

## Reactome

A subset of pathways are used for the graphlet project.

```
python3 parse_pc.py -i PathwayCommons11.reactome.hgnc.txt -o reactome
```

After generating the pathway files for all pathways, make a hand-curated list of signaling pathways to use as `reactome_pathways.txt`.  Then,

```
mkdir reactome_subset
cat reactome_pathways.txt | awk '{print "cp reactome/"$1" reactome_subset/"$1}' | bash
ls reactome_subset | wc -l
```

**TODO:** I bet that some pathways are actually combinations of other pathways (e.g. Signaling By Wnt is definitely a Reactome pathway).  Look into this more closely.
