
Laboratorio 6 - Informática II
Cálculo numérico con interfaz gráfica (Tkinter).

Repositorio GitHub: https://github.com/seba5680181-eng/Laboratorio6

import tkinter as tk
from tkinter import ttk, messagebox
import math


# Funciones

def evaluar_polinomio(coeficientes, x):
    # Uso el método de Horner para evaluar el polinomio sin calcular cada potencia.
    # Los coeficientes se ingresan desde el mayor grado hasta el término independiente.
    "Evalúa un polinomio (coeficientes están ordenados de mayor a menor grado.)"
    resultado = 0.0
    for coef in coeficientes:
        resultado = resultado * x + coef
    return resultado

    
def derivada_analitica(coeficientes, x):
    # Calculo los coeficientes de la derivada del polinomio.
    grado = len(coeficientes) - 1
    derivada = []

    # Aplico la regla de la potencia a cada término.
    for i in range(len(coeficientes) - 1):
        coef = coeficientes[i] * (grado - i)
        derivada.append(coef)

    # Evaluo la derivada en el valor de x que ingreso.
    return evaluar_polinomio(derivada, x)
