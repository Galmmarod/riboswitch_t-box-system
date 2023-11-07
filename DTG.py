#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 00:47:38 2023

@author: ringo
"""

import pandas as pd
import matplotlib.pyplot as plt
import math

plt.style.use('seaborn-v0_8-darkgrid')

FSC = pd.HDFStore('FSC.h5')
FSC.open()
for dataframe in FSC:
    print(dataframe)
    nom=dataframe.split('/')[1]
    f1=FSC[dataframe]
    #Elimina todos los organismos duplicados, sin importar si tieen diferente (anti-anti)-codon
    f1= f1.drop_duplicates(['IDs'], keep='last')
    
    # Obtener las columnas que comienzan con 'P'
    columnas_P = [col for col in f1.columns if col.startswith('P')]
    
    # Crear una lista de colores para los puntos de dispersión
    colores = ['cyan', 'red', 'green', 'magenta', 'darkorange', 'blue']
    
    # Datos de ejemplo
    x = list(range(0, len(f1)))  # Datos para el eje x
    
    y0 = f1['codon_freq']  # Datos para la línea de tendencia
    y1 = f1["promedio"]
    
    # Invertir el orden de las columnas en la lista
    columnas_P.reverse()
    
    # Definir colores para P1 y el resto de las columnas "P#"
    color_p1 = 'blue'
    color_p_max = 'red'
    
    # Colores para las otras columnas "P#"
    colores_P = {
        'P2': 'green',
        'P3': 'cyan',
        'P4': 'orange',
        'P5': 'purple',
        'P6': 'pink'
    }
    
   # Calcular el espacio entre las marcas en el eje X
    num_marcas = 10  # Número deseado de marcas en el eje X
    espacio_entre_marcas = max(1, len(x) // num_marcas)
    
    # Crear la figura y los ejes
    fig, ax = plt.subplots()
    
    
    # Graficar los puntos de dispersión y la línea de tendencia para cada columna 'P#'
    for i, columna in enumerate(columnas_P):
        y = f1[columna]
        marker = 'o' if i == len(columnas_P) - 1 else 'o'  # Marcar la última columna con '^' para mantener la propiedad
        size=8 #Tamaño de los puntos
        label = columna
        if columna == 'P1':
            color = color_p1  # Asignar color azul a P1
        elif columna == columnas_P[0]:
            color = color_p_max  # Asignar color rojo al P# más alto
        else:
            color = colores_P.get(columna)  # Asignar colores definidos
        ax.scatter(x, y, marker=marker, color=color, label=label, s=size)
    
    
    # Graficar la línea de promedio
    ax.plot(x, y1, linewidth=1, color='black', linestyle='-', label='Promedio')
    # Establecer el grosor de línea de la línea de tendencia
    ax.plot(x, y0, linewidth=1, color='yellow', linestyle='-', label='Tendencia')
    
    # Configurar los límites de los ejes
    # Ajustar el límite derecho del eje X para que se muestren todas las marcas
    ax.set_xlim(0, len(x))
    ax.set_ylim(0, math.ceil(y0.max()) + 5)
    
    # Establecer las marcas de los ejes x y y con saltos de 10 en 10
    # Establecer las marcas en el eje x con espaciado igual
    ax.set_xticks(range(0, len(x) + 1, espacio_entre_marcas))
    ax.set_yticks(range(0, math.ceil(y0.max()) + 9, 10))
    
    # Agregar un encabezado
    ax.set_title(nom)
    
    # Agregar etiquetas a los ejes
    ax.set_xlabel('Organismos Firmicutes')
    ax.set_ylabel('Porcentaje de Frecuencia')
    

    # Agregar una leyenda personalizada con el número de elementos de cada columna "P#" en 'codon_freq'
    legend_labels = []  # Lista para almacenar las etiquetas de la leyenda
    freq_menores= []
    freq_mayores=[]
    for columna in columnas_P:
        elementos_P_en_codon_freq = f1[f1[columna].isin(f1['codon_freq'])]
        num_elementos_P_en_codon_freq = len(elementos_P_en_codon_freq)
        label = f"{columna} ({num_elementos_P_en_codon_freq})"
        legend_labels.append(label)
        
    #Agregar número de frecuencias bajas y altas
    frecuencias=f1['codon_freq']
    promedio=f1['promedio']
    for porcentaje, prom in zip(frecuencias,promedio):
        if porcentaje < prom:
            freq_menores.append(porcentaje)
        elif porcentaje > prom:
            freq_mayores.append(porcentaje)
        num_freq_menores= len(freq_menores)
        num_freq_mayores= len(freq_mayores)
    menores= f"FCO Bajas ({num_freq_menores})"
    mayores= f"FCO Altas ({num_freq_mayores})"
    
    # Calcular el número total de elementos en 'codon_freq'
    num_elementos_total_codon_freq = len(f1['codon_freq'])
    
    # Agregar la etiqueta de "Tendencia" a la leyenda
    legend_labels.append('Promedio')
    legend_labels.append('Tendencia')
    legend_labels.append(menores)
    legend_labels.append(mayores)
    
    # Configurar la leyenda con las etiquetas personalizadas
    ax.legend(legend_labels, loc='center left', bbox_to_anchor=(1, 0.6))
    
    # Agregar la leyenda "Total: " con el número total de elementos en 'codon_freq'
    ax.text(1.03, 0.01, f'Total: {num_elementos_total_codon_freq}', transform=ax.transAxes, fontsize=10, verticalalignment='bottom')
    ax.text(1.03, 0.25, mayores, transform=ax.transAxes, fontsize=10, verticalalignment='bottom')
    ax.text(1.03, 0.20, menores, transform=ax.transAxes, fontsize=10, verticalalignment='bottom')
    
    # Guardar la gráfica en formato SVG
    nombre_graf = nom + '.svg'
    plt.savefig(nombre_graf, dpi=300, bbox_inches='tight')
    # Mostrar la gráfica
    plt.show()
    print('Hecho')
  
FSC.close()