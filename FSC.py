#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 12:51:43 2023

@author: ringo
"""
import pandas as pd

# Abre los archivos HDF5
SLS = pd.HDFStore('SLS.h5')
FCO = pd.HDFStore('FCO.h5')
SLS.open()
FCO.open()
# Crear un nuevo archivo HDF5 para almacenar el resultado
FSC = pd.HDFStore('FSC.h5')
count=0
#----> PRIMER PASO: Unir los dataframes de FCO y SLS
# Recorre los dataframes en SLS y busca sus contrapartes en FCO
for name_A in SLS:
    # Separar el nombre en FCO y tomar el valor [1]
    name_B = name_A.split('-')[1] if '-' in name_A else name_A
    if name_B in FCO:
        # Lee los dataframes correspondientes
        df_A = SLS[name_A]
        df_B = FCO[name_B]
        # Realiza la fusión basada en la columna 'IDs'
        df1= pd.merge(df_A, df_B, on='IDs', how='inner')
        #Podemos evaluar el progreso imprimiendo df1
        #print(df1)
        
        #----> SEGUNDO PASO: Encontrar la frequencia correspondiente al codon del loop de especificidad
        #                    y ordenarlos en una nueva columna 'codon_freq'
        
        # Inicializa una lista para almacenar los valores de 'codon_freq'
        codon_freq = []
        # Itera a través de las filas
        for index, row in df1.iterrows():
            codon_value = row['codon']
            col_name = codon_value  # Suponiendo que el valor en 'codon' coincide con el nombre de la columna
            freq_value = row[col_name]
            codon_freq.append(freq_value)
        # Agrega la lista 'codon_freq' como una nueva columna 'codon_freq'
        df1['codon_freq'] = codon_freq
        # Reorganiza el orden de las columnas para poner 'codon_freq' como la tercera columna
        columns = df1.columns.tolist()
        columns = columns[:2] + ['codon_freq'] + columns[2:-1]
        df2 = df1[columns]
        #Podemos evaluar el progreso imprimiendo df2
        #print(df2)
        
        #----> TERCER PASO: Ordenar las frecuencias de menor a mayor en nuevas columnas P#:
        
        # Obtén todas las columnas después de 'codon_freq' en un DataFrame temporal
        columnas_data = df2.columns[3:]
        df_temp = df1[columnas_data]
        # Ordena las columnas en el DataFrame temporal
        df_temp = df_temp.apply(lambda x: x.sort_values().values, axis=1, result_type='broadcast')
        # Renombra las columnas ordenadas en el DataFrame temporal
        df_temp.columns = [f'P{i}' for i in range(1, len(columnas_data) + 1)]
        # Combina el DataFrame original y el DataFrame temporal
        df3 = pd.concat([df2, df_temp], axis=1)
        #Podemos evaluar el progreso imprimiendo df3
        #print(df3)
        
        #----> CUARTO PASO: Calcular el promedio yentre el porcentaje mejor y mayor
        # Obtén las columnas P1 y la última columna P# (la más alta)
        columna_p1 = df3['P1']
        columna_pn = df3.iloc[:,-1]
        # Calcula el promedio entre P1 y la columna P# más alta
        df3['promedio'] = (columna_p1 + columna_pn) / 2
        #Podemos evaluar el progreso imprimiendo df3
        #print(df3)
        
        #----> QUINTO PASO: Calcular el ratio entre la frecuencia correspondiente al loop y el promedio
        columna_cf=df3['codon_freq']
        columna_prom= df3['promedio']
        df3['coeficiente']= (columna_cf / columna_prom)
        #Podemos evaluar el progreso imprimiendo df3
        #print(df3)
        
        #----> SEXTO PASO: Ordenar de menor a mayor las frequencias en codon_freq
        df4 = df3.sort_values('codon_freq')
        print(df4)
        
        # Almacena el resultado en el nuevo archivo HDF5
        name=name_A.split('-')[0]
        namedf= name + '-' + name_B
        namedf=namedf.split('/')[1]
        print(namedf)
        FSC.put(namedf, df4)
        count+=1
        
# Cierra los archivos HDF5
SLS.close()
FCO.close()
FSC.close()
print('Se guardaron:',count,'correlaciones')
