
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

def metodo_newton_raphson(coeficientes, x0, tolerancia, max_iteraciones):
    # Empiezo con el valor inicial que ingreso.
    x = x0

    # Repito el cálculo hasta alcanzar la tolerancia.
    for i in range(max_iteraciones):
        fx = evaluar_polinomio(coeficientes, x)
        dfx = derivada_analitica(coeficientes, x)

        # Verifico que la derivada no sea cero.
        if abs(dfx) < 1e-12:
            raise ValueError("La derivada es cero o demasiado pequeña.")

        # Calculo la siguiente aproximación de la raíz.
        x_nuevo = x - fx / dfx

        # Compruebo si la aproximación ya es suficientemente buena.
        if abs(x_nuevo - x) < tolerancia:
            return x_nuevo

        x = x_nuevo

    return x

def metodo_secante(coeficientes, x0, x1, tolerancia, max_iteraciones):
    # Calculo el valor inicial del polinomio en los dos puntos.
    f0 = evaluar_polinomio(coeficientes, x0)
    f1 = evaluar_polinomio(coeficientes, x1)

    # Repito el cálculo hasta encontrar una aproximación a la raíz.
    for i in range(max_iteraciones):
        # Verifico que no se produzca una división por cero.
        if abs(f1 - f0) < 1e-12:
            raise ValueError("No se puede continuar con el método de la secante.")

        # Calculo una nueva aproximación usando la fórmula de la secante.
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

        # Compruebo si la aproximación cumple la tolerancia.
        if abs(x2 - x1) < tolerancia:
            return x2

        # Actualizo los puntos para la siguiente iteración.
        x0 = x1
        f0 = f1
        x1 = x2
        f1 = evaluar_polinomio(coeficientes, x1)

    return x1