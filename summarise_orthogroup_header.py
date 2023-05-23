###########################################################
# Parse OrthoFinder's output file and 

from Bio import SeqIO, Entrez
import pandas as pd
Entrez.email = "elian.strozyk@etu.umontpellier.fr"


def parseInput(file) :
    orthogroup_ids = {} # stores for each orthogroup the list of IDs affiliated
    with open(file, "r") as infile :
        for line in infile :
            element = line.split()
            key = element[0].replace(":","")
            orthogroup_ids[key] = element[1:]
    return orthogroup_ids

def queryOnline(dictionary) : 
headers = []
for key in dictionary :
    handle = Entrez.efetch(db = "protein", id = ",".join(dictionary[key]), rettype = "fasta")
    for record in SeqIO.parse(handle, "fasta") :
        description = record.description.split()
        species = " ".join(description[-2:]).replace("]","").replace("[","")
        protein_name = " ".join(description[1:-2])
        headers.append([key, species, record.id, protein_name])
    handle.close()
    return headers

def makeOutput(list, outfile) :
#csv table
    df = pd.DataFrame(list, columns = ["orthogroup", "species", "NCBI_ID", "protein_name"])
    df.to_csv(outfile, index = False)

def main() :
    if len(sys.argv) == 1 :
        print("This script parses OrthoFinder's output file named Orthogroups.txt")
        print("It returns a summary file as a csv for biologists")
        print("The input file is located /Results_<date>/Orthogroups/Orthogroups.txt after running OrthoFinder")
        print("Input: a txt file named Orthogroups.txt")
        print("Output: a csv file")
