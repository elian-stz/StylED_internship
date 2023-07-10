#!/usr/bin/env python3

###########################################################
# Parse OrthoFinder's output file named Orthogroups.txt   #
# obtained after running OrthoFinder. It summarises data  #
# from fasta's header as a csv file for biologists.       #
#                                                         #
#  Input: a txt file (named Orthogroups.txt)              #
# Output: two csv files                                   #
#               one containing all orthogroups            #
#               one containing selected orthogroups       #
###########################################################
"""
Additional information:
    Orthogroups.txt location: /Results_<date>/Orthogroups/
    Snippet of Orthogroups.txt:
        OG0000000: CAH1716057.1 KAE9536424.1 KAF0773754.1 NP_001165859.1
        OG0000001: AWF70758.1 KAF0766326.1 NP_001155786.1 NP_001156143.1
        OG0000002: AAZ20451.1 AGT98613.1
    The script assumes the fasta header is written like this:
        <ID> <protein_name> <species>
        Example:
        KAE9524596.1 hypothetical protein AGLY_014646 [Aphis glycines]
"""

from Bio import SeqIO, Entrez
import pandas as pd
import sys 

def parseInput(txtfile) :
#create a dictionary: orthogroup(key); list of corresponding IDs(value)
    orthogroup_ids = {} 
    with open(txtfile, "r") as file :
        for line in file :
            element = line.split() # 1st element is orthogroup, others are IDs
            key = element[0].replace(":","")
            orthogroup_ids[key] = element[1:]
    return orthogroup_ids

def queryNCBI(dictionary) :
#query NCBI to retrieve the header of every sequence ID
    Entrez.email = "elian.strozyk@etu.umontpellier.fr"
    headers = [] #store a list of headers, each one split by info 

    for key in dictionary :
        handle = Entrez.efetch(db = "protein",
                               id = ",".join(dictionary[key]),
                               rettype = "fasta"
                              )
        for record in SeqIO.parse(handle, "fasta") :
            #assuming the header is written like this:
                # <ID> <protein_name> <species_name>
                # After split():
                # First field is ID
                # Last two fields are species e.g. [Myzus persicae]
                # Other fields are protein name
            description = record.description.split()
            species = " ".join(description[-2:]).replace("]","").replace("[","")
            protein_name = " ".join(description[1:-2])
            headers.append([key, species, record.id, protein_name])
        handle.close()
    return headers

def makeOutput(list, csvfile) :
    df = pd.DataFrame(list,
                      columns = ["orthogroup", "species", "NCBI_ID", "protein_name"]
                     )

    # all stylins identified for Acyrthosiphon pisum
    stylin_protein_ID = ["NP_001155786.1",
                        "NP_001156143.1",
                        "NP_001155431.1",
                        "NP_001165731.1",
                        "NP_001165739.1",
                        "NP_001156724.1"
                        ]
    # add an attribute: if NCBI_ID row contains the stylin ID, add a 1(True)
    df["acypi_identified_stylin"] = df["NCBI_ID"].isin(stylin_protein_ID).astype(int)
    df.to_csv(csvfile, index = False) # write the full table

    #create a subfile containing only orthogroups with A. pisum's identified stylins
    groups = df.loc[df["acypi_identified_stylin"] == 1, "orthogroup"]
    df2 = df[df["orthogroup"].isin(groups.tolist())]
    df2.to_csv("selected_" + csvfile, index = False) # write the subfile 

def main() :
    if len(sys.argv) != 5 :
        print("This script parses OrthoFinder's output file named Orthogroups.txt")
        print("It returns two csv files summarising header's fasta")
        print("The input file is located in /Results_<date>/Orthogroups/ after running OrthoFinder")
        print("Input: a txt file (named Orthogroups.txt)")
        print("Output: two csv files")
        print("        one containing all orthogroups")
        print("        one containing only selected orthogroups with proteins of interest (identified stylins)")
        print("To run this script, type:")
        print("python3", sys.argv[0], "--i <Orthogroups.txt> --o <csv file>")
    else :
        print("Parsing the file")
        data = parseInput(sys.argv[2])
        print("Querying NCBI. This might take a few time...")
        headers = queryNCBI(data)
        print("Writing the output")
        makeOutput(headers, sys.argv[4])
        print("Summary completed")

if __name__ == "__main__" :
    main()
