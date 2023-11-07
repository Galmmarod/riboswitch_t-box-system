import pandas as pd
from collections import Counter
import os.path
codon_table={
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
    'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W',
}
lista_codones = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT'],
    'X': ['TAA', 'TGA', 'TAG']
}

name_aa= {'A' : 'Alanina', 'C' : 'Cisteína', 'D' : 'Aspartato', 'E': 'Glutamato',
     'F' : 'Fenilalanina','G' : 'Glicina','H' : 'Histidina','I' : 'Isoleucina',
     'K' : 'Lisina','L': 'Leucina','M' : 'Metionina','N' : 'Asparagina',
     'P' : 'Prolina','Q' : 'Glutamina','R' : 'Arginina','S' : 'Serina',
     'T' : 'Treonina','V' : 'Valina','W' : 'Triptófano','Y': 'Tirosina','X': 'STOP'}

carpeta= "/home/fastas_2206/"
ids= "/home/galmanza/Riboswitches/Programs/orgs_S21_con_Tbox"
conteo=0
with open(ids, "r") as s:
    orgs=s.read().split()
    f1 = pd.DataFrame()
    f1['IDs']=(orgs)
df = {}
ids_faa = os.popen('ls /home/fastas_2206/*.faa | cut -f4 -d "/" | cut -f1 -d "."').read()
ids_faa= ids_faa.split('\n')
error = ids_faa.pop()
ids_fna = os.popen('ls /home/fastas_2206/*.fna | cut -f4 -d "/" | cut -f1 -d "."').read()
ids_fna= ids_fna.split('\n')
error2 = ids_fna.pop()
for amino_acid, codons in lista_codones.items():
    name = name_aa[amino_acid]
    letras=list(name_aa)
    first_column={'ÓRG':[orgs]}
    df[name]= pd.DataFrame(columns=codons)
    df[name]=f1.join(df[name])
    #Iterar entre todos los archivos y ejecutar el código
#print(ids_faa)
for ides in orgs:
    #print(orgs)
    if ides in ids_faa:
        #Extraer el código del organismo de los archivos en la carpeta
            conteo+=1
            #print('>>',codigo_fna)
            fna= carpeta + ides + ".fna"
            faa= carpeta + ides + ".faa"
            with open(fna, "r") as f:
                sequence=f.read().strip()
            with open(faa, "r") as prot:
                idprot=prot.read().strip()
                #Separación de genes por id y secuencia
            codons_by_aa= {}
            gene_data= sequence.split(">")[1:]
            prot_data= idprot.split(">")[1:]
            prot_id_list=[]
            for prot in prot_data:
                prot_id= prot.split("\n")[0]
                prot_id_list.append(prot_id)
            for gene in gene_data:
                gene_id, seq= gene.split("\n")[0], "".join(gene.split("\n")[1:]).upper()
                gene_length= len(seq)
                # Separa la secuencia en codones y los asigna a los aminoácidos correspondientes
                if gene_id in prot_id_list:
                    for i in range(0, len(seq)-2, 3):
                        codon= seq[i:i+3]
                        if codon in codon_table:
                            aa= codon_table[codon]
                            if aa not in codons_by_aa:
                                codons_by_aa[aa]= []
                            codons_by_aa[aa].append(codon)
            # Cuenta el número de codones por aminoácido
            counts_by_aa= {}
            for aa, codons in codons_by_aa.items():
                counts_by_aa[aa]= Counter(codons)
            # Calcula las frecuencias de los codones por aminoácido
            for aa, counts in sorted(counts_by_aa.items()):
                    total_codons= sum(counts.values())
                    frequencies= {codon: count/total_codons*100 for codon, count in sorted(counts.items())}
                    name=name_aa[aa]
                    for codones in frequencies:
                        df[name].loc[conteo-1:conteo,codones]=frequencies[codones]
                        
                        

# Nombre del archivo HDF5 donde se guardarán los DataFrames
hdf5_file = 'FCO.h5'
# Iterar a través del diccionario de DataFrames
for name, dataframe in df.items():
    if name != 'Triptófano' and name != 'Metionina' and name != 'STOP':
        # Crear un archivo HDF5
        with pd.HDFStore(hdf5_file, mode='a') as hdf5:
        # Guardar cada DataFrame en una tabla dentro del archivo HDF5
            hdf5[name]=dataframe
            print(hdf5[name])