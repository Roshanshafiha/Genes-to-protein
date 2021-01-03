#open the file
infile=open("lactasegene.txt","r")
dnasequence=""
for ln in infile.readlines():
    dnasequence = dnasequence+ln.strip()#this strip() function is used to remove the spaces
finaldna=dnasequence
#load in exon interval (start and end values)
def assignexonintervals(file_name):
    open_file=open(file_name)
    exon_intervals=[]
    for ln in open_file.readlines():
        exon_values=ln.split()
        start=int(exon_values[0])
        end=int(exon_values[1])
        dictionary={'start':start,'end':end}
        exon_intervals.append(dictionary)
    return exon_intervals
#create mrna
def mrnasequence(sequence,exon_range):
    mrna=''
    for dic in exon_range:
        start=dic['start']
        end=dic['end']
        mrna=mrna+sequence[start:end].replace('T','U')
    return(mrna)
#protein dictionary code
def load_protein_dic(file_name):
    new_file=open (file_name)
    protein_dict={}
    for ln in new_file.readlines():
        toks=ln.split()
        protein_dict[toks[0]] = {'1-letter':toks[1], '3-letter':toks[2], 'amino acid':toks[3]}
    return (protein_dict)
#to get the protein sequence :
def proteinsequence (mrna_sequence,protein_dictionary):
    proteinstring=''
    start_codon='AUG'
    start_location=mrna_sequence.find(start_codon)
    while start_location <len(mrna_sequence):
        codon=mrna_sequence[start_location:start_location+3]
        protein=protein_dictionary[codon]['1-letter']
        if protein != 'X':
            proteinstring += protein
            start_location += 3
        else:    
            break
    return (proteinstring)

#load the exon range function first
file_name=('lactase_exon.tsv')
exon_range=assignexonintervals(file_name)
#then assign the dna sequence 
sequence=finaldna
#load the mrnasequence and assign it to the mrnaseq 
mrna_sequence=mrnasequence(sequence,exon_range)
#load the protein dictionary 
file_name=('genetic_code.tsv')
protein_dictionary=load_protein_dic(file_name)
#obtain the protein sequence
protein_sequence_result=proteinsequence(mrna_sequence,protein_dictionary)
#find the length of the protein sequence 
length_proteinseq=len(protein_sequence_result)
#check the last ten amino acid 
last_ten_Aminoacid=protein_sequence_result[-10:]
#print the statement and compare the result with the mutated sequence 
print("The length of non mutated protein sequence: ",length_proteinseq,"The last 10 aminoacid of non mutated sequence:",last_ten_Aminoacid)

