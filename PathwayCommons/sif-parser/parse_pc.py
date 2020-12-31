import sys
import os
import argparse
import itertools
import glob

def main(args):

    proteins = read_proteins(args.infile)
    print('%d proteins processed' % (len(proteins)))

    interactions_by_pathways = read_interactions(args.infile,proteins)
    print('%d pathways processed' % (len(interactions_by_pathways)))

    for pathway in interactions_by_pathways.keys():
        print('Pathway "%s" has %d edges' % (pathway,len(interactions_by_pathways[pathway])))
        if len(interactions_by_pathways[pathway]) > args.thres:
            outfile = '%s/%s-edges.txt' % (args.outdir,pathway.replace(' ','-').replace('/','-or-').replace('(','').replace(')',''))
            write_file(interactions_by_pathways[pathway],proteins,outfile)
        else:
            print('  not writing %s -- not enough edges.' % (pathway))

    print('done!')
    return

'''
Reads all the nodes in the file.
Ignores any entries that are not ProteinReferences (e.g. small molecules, RNA, DNA, etc.).
Ignores any entries that do not have a uniprot ID.
'''
def read_proteins(infile):
    proteins = {} # common name to uniprot.
    num_missed = 0
    found_header = False
    with open(infile) as fin:
        for line in fin:
            if not found_header and 'PARTICIPANT_TYPE' not in line:
                continue
            elif not found_header and 'PARTICIPANT_TYPE' in line:
                found_header = True
                continue
            else: # check entry.
                if 'ProteinReference' not in line:
                    num_missed +=1
                    continue
                row = line.strip().split('\t')
                if len(row) < 4 or 'uniprot' not in row[3]:
                    num_missed+=1
                    continue
                proteins[row[0]] = row[3].split(':')[-1]
    print('%d of %d (%.2f) missed' % (num_missed,num_missed+len(proteins),num_missed/(num_missed+len(proteins))))
    return proteins

'''
Parsed undirected edges from SIF.
Sorts each edge to ensure that only one of (u,v) or (v,u) are added.
Only considers edges where both nodes are in the proteins dictionary.
Does not handle the MEDIATOR_IDS.
'''
def read_interactions(infile,proteins):
    pathways = {} # pathway to edge tuple
    DELIM = ';'
    num_missed = 0
    with open(infile) as fin:
        for line in fin:
            if 'PARTICIPANT_TYPE' in line:
                # done reading interactions; return.
                return pathways
            if 'INTERACTION_TYPE' in line:
                # skip header
                continue
            row = line.strip().split('\t')
            if row[0] not in proteins or row[2] not in proteins:
                num_missed +=1
                continue
            edge = tuple(sorted([row[0],row[2]]))
            for p in row[5].split(';'):
                p = p.strip()
                if p == '':
                    continue
                if p not in pathways:
                    pathways[p] = set()
                pathways[p].add(edge)
    ## this should never happen (we should return when we hit PARTICIPANT_TYPE).
    ## return None in this case, it's an error.
    return None

def write_file(edges,proteins,outfile):
    out = open(outfile,'w')
    for e in edges:
        out.write('\t'.join([e[0],e[1],proteins[e[0]],proteins[e[1]]])+'\n')
    out.close()
    print('  wrote to %s' % (outfile))
    return

def parse_arguments():
    """
    Argument Parser for parse_pc.py.

    Returns
    -----------
    ArgumentParser object

    """
    parser = argparse.ArgumentParser('PathwayCommons Parser.  Converts pathways into undirected graphs with UniProtKB identifiers.')
    parser.add_argument('-i','--infile',
        help='BioPAX file (.owl format). Required.',required=True)
    parser.add_argument('-o','--outdir',
        help='outfile directory. Default = out.',default='out/')
    #parser.add_argument('-f','--filter',
    #    help='Filter converted IDs by single-column file of ids. Only IDs that appear in this file will be used.')
    parser.add_argument('-t','--thres',type=int,default=10,
        help='Do not write pathways with fewer than THRES edges. Default 10.')
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        print('making output directory %s...' % (args.outdir))
        os.makedirs(args.outdir)
    return args

if __name__ == '__main__':
    main(parse_arguments())
