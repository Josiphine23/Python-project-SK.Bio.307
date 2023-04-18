
"""
Script that designs a forward and a reverse primer for a Fasta file as input
Writer: Josephine Franke
Last updated: 18.04.2023

"""

import os

# User should give a fasta file as input. If no fasta-file is given an error message is displayed
error = True
while error == True:
    file_path = input("Please type the file path of your fasta file: ")
    name, file_type = os.path.splitext(file_path)
    if file_type != ".fa" and file_type != ".fasta":
        print("Please give a fasta-file as input.")
    else:
        error = False

# content of given fasta file is saved in variable gene_content
file = open(file_path, "r")
gene_content_list = file.readlines()
# every line is saved as an element of the list gene_content
file.close()

# the first line which is always the header line in a fasta file is deleted
gene_content_list = gene_content_list[1:]
# the join-method takes every element of gene_content_list and concatenates it to the string gene_content_str
gene_content_str = "".join(gene_content_list)

# deletes all \n characters
gene_content_strip = gene_content_str.rstrip("\n")

# saving substring starting with ATG until the first appearance of a stop codon in the new variable coding_region
start = gene_content_strip.find("ATG")
stop_list = ["TGA", "TAG", "TAA"]
index_list = []
for elm in stop_list:
    index = gene_content_strip.find(elm)
    # only saves index of stop codon in index_list when stop codon exists and it appears after the first ATG
    if index != -1 and index >start:
        index_list.append(index)

# prints an error message when no stop codon could be found
codon = False
for elm in index_list:
    if elm != -1:
        codon = True
if codon == False:
    print("No Stop codon was found. Please enter a file containing a complete reading frame.")
    exit()

# finding the first appearance of a stop codon
first_stop = min(index_list)
# saving the substring from start codon ATG until the first stop codon in new variable
coding_region = gene_content_strip[start:first_stop + 3]

# printing an error message when the coding region is too short for the primer design
if len(coding_region) < 20:
    print("The coding region is shorter than 20 nucleotides and no primer design is possible.")
    exit()

# Now we have a variable with only the sequence of the coding region and we can start design the primers

# defining a function design_primer_fwd that returns the complementary strand for a given substring from the coding region
def design_primer_fwd(substring):
    primer_fwd = ""
    for elm in substring:
        if elm == "A":
            primer_fwd = primer_fwd + "A"
        elif elm == "T":
            primer_fwd = primer_fwd + "T"
        elif elm == "G":
            primer_fwd = primer_fwd + "G"
        elif elm == "C":
            primer_fwd = primer_fwd + "C"
    return primer_fwd
#defining a function design_primer_rev that returns the complementary strand for a reversed substring from the coding region
def design_primer_rev(substring):
    primer_rev = ""
    #reverses the given substring
    for elm in substring[::-1]:
        if elm == "A":
            primer_rev = primer_rev + "T"
        elif elm == "T":
            primer_rev = primer_rev + "A"
        elif elm == "G":
            primer_rev = primer_rev + "C"
        elif elm == "C":
            primer_rev = primer_rev + "G"
    return primer_rev
#defining a function that calculates the annealing temperature of a given substring(the primer)
def calc_AT(primer):
    numb_A = primer.count("A")
    numb_T = primer.count("T")
    numb_G = primer.count("G")
    numb_C = primer.count("C")

    AT = 2 * (numb_A + numb_T) + 4 * (numb_G + numb_C)

    return AT

#generating the substrings which will be the input for the design_primer function
sub_fwd_list =[coding_region[0:20], coding_region[0:21], coding_region[0:22], coding_region[0:23]]
sub_rev_list =[coding_region[-20:], coding_region[-21:], coding_region[-22:], coding_region[-23:]]

#saves fwd primer in new list primer_fwd_list
primer_fwd_list = []
for elm in sub_fwd_list:
    primer_fwd_list.append(design_primer_fwd(elm))
#saves rev primers in new list primer_rev_list
primer_rev_list = []
for elm in sub_rev_list:
    primer_rev_list.append(design_primer_rev(elm))
#saves annealing temperatures of all fwd primer in new list AT_fwd_list
AT_fwd_list = []
for elm in primer_fwd_list:
    AT_fwd_list.append(calc_AT(elm))
#saves annealing temperatures of all rev primer in new list AT_rev_list
AT_rev_list = []
for elm in primer_rev_list:
    AT_rev_list.append(calc_AT(elm))

# combining each primer sequence with its corresponding annealing temperature in a dictionary
# the function lambda combines two elements with the same index from two lists into a tuple and then saves it as key and value in a dictionary
# the function map iterates through the lists primer_fwd_list and AT_fwd_list and applies the function lambda
dict_fwd = dict(map(lambda i,j : (i,j), primer_fwd_list, AT_fwd_list))
dict_rev = dict(map(lambda i,j : (i,j), primer_rev_list, AT_rev_list))

# deleting primers from the dictionaries which have an AT which is lower than 55°C or higher than 62°C
# the method items() allows for iteration through keys and values in a dictionary
for key, value in dict(dict_fwd).items():
    if value < 55 or value > 62:
        del dict_fwd[key]
for key, value in dict(dict_rev).items():
    if value < 55 or value > 62:
        del dict_rev[key]

# when dict is empty because all the AT are <55°C or >62°C, an error message is printed
if bool(dict_fwd) == False:
    print("No forward primer sequences with annealing temperatures between 55°C and 62 °C exist for this file")
    exit()
if bool(dict_rev) == False:
    print("No reverse primer sequences with annealing temperatures between 55°C and 62 °C exist for this file")
    exit()

# extracting the two sequences from dict_fwd and dict_rev whose AT does not diverge more than 4 degrees
solution = False
for (k_f,v_f) in dict_fwd.items():
    for (k_r,v_r) in dict_rev.items():
        # when solution is true the first combination that has been found is printed. If this would be deleted, all combinations would be printed.
        if solution == True:
            continue
        if abs(dict_fwd[k_f]-dict_rev[k_r]) <=4:
            print("Forward primer: ", k_f, " | ", "Annealing temperature: ", dict_fwd[k_f], "°C")
            print("Reverse primer: ", k_r, " | ", "Annealing temperature: ", dict_rev[k_r], "°C")
            solution = True

# an error message is printed when the primer sequences diverge more than 4 degrees from each other
if solution == False:
    print("The determined primer sequences diverge more than 4°C from each other")