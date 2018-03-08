"""
parses raw TCGA SNP txt file
returns SNP call if confidence is 0.05 or lower (birdseed call)
otherwise sets SNP value to NA
txt file output should be ready for easy loading in R

(to run)
cd into directory with SNP txt file and run in command line: 
python parse_TCGA_SNP.py FILENAME.txt  

"""

from sys import argv
import re
import time

def format_header(line):
    """ parses first line in file """
    l = line.strip("\n").split("\t")
    l = [l[0]] + l[1::2]
    return "\t".join(l) + "\n"

def get_snp_vals(line, thresh):
    """ parses line, returns SNP call if confidence <= thresh else NA
    only valid for 3rd to last lines in txt file """
    l = line.strip("\n").split("\t")
    zipped = zip(l[1::2], l[2::2])
    snpvals = [x[0] if float(x[1]) <= float(thresh) else 'NA' for x in zipped]
    snpvals.insert(0, l[0])
    return "\t".join(snpvals) + "\n"
    
start = time.time()
script, filename = argv
match = re.match(r'^[A-Z]*', filename)
output = match.group() + "_SNP.txt"
with open(filename, "r") as infile:
    print "opening file %r:" % filename
    with open(output, "w") as outfile:
        for i, line in enumerate(infile):
            if i == 0:
                outfile.write(format_header(line = line))
            elif i == 1:
                pass
            else: 
                outfile.write(get_snp_vals(line = line, thresh = 0.05))
print "parsing complete on %r:" %output
print "script took %r seconds" %(time.time() - start) 

    
