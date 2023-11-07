#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:21:02 2023

@author: ringo
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

FSC=pd.HDFStore('FSC.h5')
FSC.open()

IDs=[]
amino=[]
for dataframe in FSC:
    print(dataframe)
    name1=dataframe.split('-')[1]
    name2=dataframe.split('-')[0]
    name2=name2.split('/COG')[1]
    name=name1+ '-' +name2
    print(name)
    df=FSC[dataframe]
    amino.append(name)
    amino=sorted(amino)
    if 'IDs' in df.columns:
        ids=df['IDs'].tolist()
        IDs.extend(ids)
    list_ids=sorted(set(IDs))    
data={'IDs':list_ids}
for aa in amino:
    data[aa]=None
df1=pd.DataFrame(data)
    
for dataframe in FSC:
    df=FSC[dataframe]
    name1=dataframe.split('-')[1]
    name2=dataframe.split('-')[0]
    name2=name2.split('/COG')[1]
    name=name1+ '-' +name2
    colum=[]
    orgs=[]
    for org in list_ids:
        if org in df['IDs'].tolist():
            valor= df.loc[df['IDs']==org, 'ratio']
            colum.extend(valor.tolist())
            orgs.append(org)
    for codigo,dato in zip(orgs, colum):
        if codigo in df1['IDs'].tolist():
            df1.loc[df1['IDs'] == codigo, name] = dato
df_ids = df1[['IDs']]
# Crear un DataFrame con las demás columnas
df_calculo = df1.iloc[:,1:]
#print(df_calculo)
# Calcula la cantidad de números en cada columna en el rango de 0.001 a 2
counts = df_calculo.apply(lambda col: col.between(0.001, 2).sum())
# Ordena las columnas en función de los recuentos
valores= counts.sort_values(ascending=False)
print("****"*30)
print(valores)
sorted_columns = counts.sort_values(ascending=False).index
# Aplica la nueva secuencia de columnas al DataFrame
print(sorted_columns)
df_sorted = df1[sorted_columns]
df_new = df_sorted.fillna(0)
df_final = df_ids.join(df_new)
df_prueba = df_ids.join(df_sorted)
#print(df_final)
#print(df_new)
#df_final.to_csv("DataFrame_Temporal.csv")

df=df_final
# Asignar un número entero único a cada ID
df['ID_encoded'] = df['IDs'].astype('category').cat.codes
# Establecer la columna 'ID_encoded' como índice y eliminar la columna 'ID'
df.set_index('ID_encoded', inplace=True)
df.drop('IDs', axis=1, inplace=True)
# Crear un mapa de colores personalizado
colors = ['#FFFFFF' , '#0000FF', '#0000FF','#0000FF', '#0000FF','#0000FF','#0000FF','#0000FF','#FF0000','#FF0000','#FF0000' ,'#FF0000','#FF0000','#FF0000','#FF0000','#FF0000' ]
cmap = mcolors.LinearSegmentedColormap.from_list('custom_cmap', colors)
# Crear el heatmap con el mapa de colores personalizado
plt.figure(figsize=(8, 6))
heatmap = plt.pcolor(df, cmap=cmap, vmin=0, vmax=df.values.max())
# Agregar etiquetas en el eje x
plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation='vertical')
# Agregar etiquetas en el eje y
yticks = np.arange(0, len(df.index), 50)
yticks = np.append(yticks, len(df.index))  # Agregar el número total (580) al final
plt.yticks(yticks)
# Agregar una barra de color a la derecha del heatmap
plt.colorbar(heatmap)
#Guardar heatmap en una imagen svg
#plt.savefig("heatmap.svg", dpi=300, bbox_inches='tight')
# Mostrar el heatmap
#plt.show()
FSC.close()