# seek_uniq_seq.py, author: vandijk@evolbio.mpg.de

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from datetime import datetime
import sys, getopt, os
import glob, shutil

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size. (https://biopython.org/wiki/Split_large_file)

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                #entry = iterator.next() // Depricated due to python 3
                entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
            
            
def main(argv):

    argsgiven = 0
    query = ''
    subject = ''
    build_DB = True;    
    usage = 'seq_uniq_seek.py -q <queryfile>.fasta -s <subjectfile>.fasta -B [build database true/false]'
    verbal=True
    opts, args = getopt.getopt(argv,"xmhq:s:o:",["subject=","query="])
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt == '-x':
            build_DB = False;
        elif opt in ("-q", "--query"):
            query = arg
            argsgiven+=1
        elif opt in ("-m", "--mute"):
            verbal = False
            argsgiven+=1
        elif opt in ("-s", "--subject"):
            subject = arg
            argsgiven+=1
        elif opt in ("-o", "--output"):
            output = arg
            argsgiven+=1
    if(argsgiven < 3):
        print(usage)
        sys.exit(2)

            
    if(verbal): print("\n ---- ==== SEEK UNIQ SEQ ==== ---- \nFinding sequences occuring in "+query+" that are not occuring in "+subject+" and saving in "+output+".fasta\n")
    
    if(build_DB):
        if(verbal): print("Building blast database for subject file ("+subject+")")
        makedb = NcbimakeblastdbCommandline(cmd='makeblastdb', input_file=subject, dbtype='nucl', parse_seqids=True)
        makedb()
        if(verbal): print("Done.\n")
    else:
        if(verbal): print("Not building database. Hoping for the best") 
    
    if(verbal): print("Blasting query ("+query+") against subject database ("+subject+")")
    
    if(verbal): print("Splitting query into multiple files to save memory.")
    
    shutil.rmtree("chunks", ignore_errors=True)
    os.mkdir("chunks")    
    
    record_iter = SeqIO.parse(open(query), "fasta")
    for i, batch in enumerate(batch_iterator(record_iter, 10000)):
        filename = "chunks/chunk_%i.fasta" % (i + 1)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
    
    
    
    if(verbal): print("Building query index dictionary")
    q_dict = SeqIO.index(query,"fasta") 
    hits = []
    
    
    chunks = glob.glob('chunks/chunk*')
    for i,file in enumerate(chunks):
        now = datetime.now()
        dt_string = now.strftime("%d-%m_%H:%M:%S")
        print("[xenoseq_blast "+dt_string+"] So anyway... I'm busy blasting... "+str(round(i/len(chunks)*100,2))+"%")
        
        blastn_cline = NcbiblastnCommandline(cmd='blastn',query=file, db=subject, num_threads=8, evalue=1e-5, perc_identity=90, outfmt=5, out="reads_all_vs_all.xml")
       
        blastn_cline()
        
        
        
        # Bit below is from: https://biopython.org/wiki/Retrieve_nonmatching_blast_queries
        
        for record in NCBIXML.parse(open("reads_all_vs_all.xml")):
            for alignment in record.alignments:
                if(alignment.length > 100):                    
                    hits.append(record.query.split()[0])
        os.remove("reads_all_vs_all.xml") 
        
    shutil.rmtree("chunks")  
    
    if(verbal): print("Subtracting hits from query dict keys")
    misses = set(q_dict.keys()) - set(hits)
    orphans = [q_dict[name] for name in misses]    
    if(verbal): print("%i out of %i records in query are unique" % (len(misses), len(q_dict)))
    if(verbal): print("Writing to file %s" % (output))
    SeqIO.write(orphans, output, 'fasta')
    if(verbal): print("Done. Hoping for the best.\n")
    
    
if __name__ == "__main__":
    main(sys.argv[1:])


