#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 15:27:57 2023

@author: ringo
"""
import re
import pandas as pd
#from collections import Counter
codon_table={
    'UUU':'F', 'UCU':'S', 'UAU':'Y', 'UGU':'C',
    'UUC':'F', 'UCC':'S', 'UAC':'Y', 'UGC':'C',
    'UUA':'L', 'UCA':'S', 'UAA':'X', 'UGA':'X',
    'UUG':'L', 'UCG':'S', 'UAG':'X', 'UGG':'W',                
    'CUU':'L', 'CCU':'P', 'CAU':'H', 'CGU':'R',
    'CUC':'L', 'CCC':'P', 'CAC':'H', 'CGC':'R',
    'CUA':'L', 'CCA':'P', 'CAA':'Q', 'CGA':'R',
    'CUG':'L', 'CCG':'P', 'CAG':'Q', 'CGG':'R',
    'AUU':'I', 'ACU':'T', 'AAU':'N', 'AGU':'S',
    'AUC':'I', 'ACC':'T', 'AAC':'N', 'AGC':'S',
    'AUA':'I', 'ACA':'T', 'AAA':'K', 'AGA':'R',
    'AUG':'M', 'ACG':'T', 'AAG':'K', 'AGG':'R',
    'GUU':'V', 'GCU':'A', 'GAU':'D', 'GGU':'G',
    'GUC':'V', 'GCC':'A', 'GAC':'D', 'GGC':'G',
    'GUA':'V', 'GCA':'A', 'GAA':'E', 'GGA':'G',
    'GUG':'V', 'GCG':'A', 'GAG':'E', 'GGG':'G',}

def frecuencia_codones(codones_encontrados, total_modelos):
    frecuencias = {}
    for codon in codones_encontrados:
        if codon not in frecuencias:
            frecuencias[codon] = 0
        frecuencias[codon] += 1
    for codon in frecuencias:
        frecuencias[codon] = frecuencias[codon] 
    frecuencias_ordenadas = sorted(frecuencias.items())
    return frecuencias_ordenadas

aa_list = {
    'A': 'Alanine',
    'C': 'Cysteine',
    'D': 'Aspartate',
    'E': 'Glutamate',
    'F': 'Phenylalanine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'K': 'Lysine',
    'L': 'Leucine',
    'M': 'Methionine',
    'N': 'Asparagine',
    'P': 'Proline',
    'Q': 'Glutamine',
    'R': 'Arginine',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
    'X': 'STOP'
}

print('Tabla de codones:\nA = Alanina     ','\tC = Cisteína ','\tD = Aspartato','\tE= Glutamato\n'
     'F = Fenilalanina','\tG = Glicina  ','\tH = Histidina    ','\tI = Isoleucina\n'
     'K = Lisina      ','\tL = Leucina  ','\tM = Metionina    ','\tN = Asparagina\n'
     'P = Prolina     ','\tQ = Glutamina','\tR = Arginina     ','\tS = Serina\n'
     'T = Treonina    ','\tV = Valina   ','\tW = Triptófano   ','\tY = Tirosina\nX = STOP')

archivo='./COGS/COG2269'
aminoacido=input('Escribe la letra del aminoácido correspondiente al T-Box: ')
count=0
df={}

#df[T_box]=pd.DataFrame()
with open(archivo) as f:
    codones = []
    print('\n')
    print('*************',archivo,'*************:')
    print('*************',aa_list[aminoacido],'\n')
    ids=[]
    codoones=[]
    for i, line in enumerate(f):
        if line.startswith(">>"):
            line_id= line.split(">> ")
            idd="  ".join(line_id[0:])
            id_gen=line.strip('\n')
            for n in range(5):
                f.readline()
            line2=f.readline()
            line3=f.readline()
            none= line3.split()
            non= none[2].startswith('<') 
            if re.search(r"\s*::", line2) and (non==False):
                line4 = f.readline()
                start = line2.find('->')
                end = line2.find('>,', start)+1
                #Caracteres
                matching_1 = line3[start:end].upper()
                start2= matching_1.find('GAA')
                end2= matching_1.find('CCCC', start2)
                match1= matching_1[start2:end2]
                non1= line3[:start].find('[')
                if (non1== -1):
                    matching_0=line2[start:end]
                    match0= matching_0[start2:end2]
                    #Modelo
                    line5= f.readline()
                    matching_3= line5[start:end]
                    match3= matching_3[start2:end2].upper()
                    seq="".join(match3.split("-")[0:]).upper()
                    codones_aminoacido = [codon for codon, amino in codon_table.items() if amino == aminoacido]
                    # Buscamos el primer codón que aparece en la secuencia y lo guardamos en una variable
                    encontrado = False
                    i = 0
                    while not encontrado and i < len(seq)-2:
                        codon = seq[i+3:i+6]
                        if codon in codones_aminoacido:
                            encontrado = True
                            codones.append(codon)
                        else:
                            i += 1     
                    # Si se ha encontrado un codón, imprimimos la secuencia y el codón encontrado
                    if encontrado:
                        codigo=idd.split('-')
                        codigo=codigo[0].split('  ')
                        cod=codigo[1]
                        print(id_gen)
                        print(match0)
                        print(match1)
                        print(seq)
                        print(' '*(i+3) + codon)
                        ids.append(cod)
                        codoones.append(codon)
                        count+=1

archivo_name=archivo.split('./COGS/')[1]
T_box= archivo_name + '-' + aa_list[aminoacido] + '-tRNA_sintetasa' 
df[T_box]= pd.DataFrame({'IDs':ids, 'codon':codoones})
print('******************************************************')
print(T_box)
print('Cantidad de modelos analizados:', count)
frecuencias_codones = frecuencia_codones(codones, count)
unique_codones=[]
for codones in codones:
  if codones not in unique_codones:
    unique_codones.append(codones)
print('Codones encontrados:',sorted(unique_codones))
print('Frecuencias de codones:')
for codon, frecuencia in frecuencias_codones:
    frecuencias= codon , ':' , frecuencia
    print(codon,':' , frecuencia)
    

# Nombre del archivo HDF5 donde se guardarán los DataFrames
hdf5_file = 'SLS.h5'
# Iterar a través del diccionario de DataFrames
for name, dataframe in df.items():
    # Crear un archivo HDF5
    with pd.HDFStore(hdf5_file, mode='a') as hdf5:
    # Guardar cada DataFrame en una tabla dentro del archivo HDF5
        hdf5[name]=dataframe
        print(hdf5[name])
