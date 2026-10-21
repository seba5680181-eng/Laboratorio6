
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


def metodo_biseccion(coeficientes, a, b, tolerancia, max_iteraciones):
    # Calculo el valor del polinomio en los extremos.
    fa = evaluar_polinomio(coeficientes, a)
    fb = evaluar_polinomio(coeficientes, b)

    # Verifico que haya un cambio de signo en el intervalo.
    if fa * fb > 0:
        raise ValueError("El intervalo no contiene un cambio de signo.")

    # Repito el cálculo hasta alcanzar la tolerancia o el máximo.
    for i in range(max_iteraciones):
        medio = (a + b) / 2
        fm = evaluar_polinomio(coeficientes, medio)

        # Compruebo si encontré una raíz aproximada.
        if abs(fm) < tolerancia or (b - a) / 2 < tolerancia:
            return medio

        # Elijo la mitad del intervalo que contiene el cambio de signo.
        if fa * fm < 0:
            b = medio
            fb = fm
        else:
            a = medio
            fa = fm

    return (a + b) / 2
