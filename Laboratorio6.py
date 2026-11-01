# Laboratorio 6 - Informática II
# Alumnos: Sebastian Vallejos - Larrea Toftum Cristian
# Cálculo numérico
# Repositorio GitHub:
# https://github.com/seba5680181-eng/Laboratorio6
# Se inicia con "python Laboratorio6.py"



# FUNCIONES MATEMATICAS

def evaluar_polinomio(coeficientes, x):
    # Uso el método de Horner para evaluar el polinomio sin calcular cada potencia.
    # Los coeficientes se ingresan desde el mayor grado hasta el término independiente.

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

    # Evalúo la derivada en el valor de x.
    if not derivada:
        return 0.0

    return evaluar_polinomio(derivada, x)


def derivada_numerica(coeficientes, x, h):
    # Calculo la derivada mediante diferencia finita centrada.

    if h <= 0:
        raise ValueError("El paso h debe ser mayor que cero.")

    f_adelante = evaluar_polinomio(coeficientes, x + h)
    f_atras = evaluar_polinomio(coeficientes, x - h)

    derivada = (f_adelante - f_atras) / (2 * h)

    return derivada


# METODO DE BISECCION


def metodo_biseccion(coeficientes, a, b, tolerancia, max_iteraciones):

    # Calculo el valor del polinomio en los extremos.
    fa = evaluar_polinomio(coeficientes, a)
    fb = evaluar_polinomio(coeficientes, b)

    if a >= b:
        raise ValueError("El límite inferior debe ser menor que el superior.")

    # Verifico que exista un cambio de signo.
    if fa == 0:
        return a, 0, 0.0, fa

    if fb == 0:
        return b, 0, 0.0, fb

    if fa * fb > 0:
        raise ValueError("El intervalo no contiene un cambio de signo.")

    anterior = None

    # Repito el cálculo hasta alcanzar la tolerancia.
    for i in range(1, max_iteraciones + 1):

        medio = (a + b) / 2
        fm = evaluar_polinomio(coeficientes, medio)

        # Calculo el error aproximado relativo porcentual.
        if anterior is None:
            error = None
        elif medio != 0:
            error = abs((medio - anterior) / medio) * 100
        else:
            error = abs(medio - anterior)

        # Compruebo si encontré la raíz exacta.
        if fm == 0:
            return medio, i, 0.0, fm

        # Compruebo si se alcanzó la tolerancia.
        if error is not None and error <= tolerancia:
            return medio, i, error, fm

        # Elijo el intervalo que contiene el cambio de signo.
        if fa * fm < 0:
            b = medio
            fb = fm
        else:
            a = medio
            fa = fm

        anterior = medio

    raise ValueError(
        "Se alcanzó el máximo de iteraciones sin cumplir la tolerancia."
    )

  # METODO NEWTON-RAPHSON
def metodo_newton_raphson(coeficientes, x0, tolerancia, max_iteraciones):

    # Empiezo con el valor inicial ingresado.
    x = x0

    for i in range(1, max_iteraciones + 1):

        fx = evaluar_polinomio(coeficientes, x)

        # Compruebo si el valor inicial ya es una raíz.
        if fx == 0:
            return x, i - 1, 0.0, fx

        # Calculo la derivada analítica.
        dfx = derivada_analitica(coeficientes, x)

        if abs(dfx) < 1e-12:
            raise ValueError(
                "La derivada es cero o demasiado pequeña."
            )

        # Calculo la siguiente aproximación.
        x_nuevo = x - fx / dfx

        # Calculo el error aproximado relativo porcentual.
        if x_nuevo != 0:
            error = abs((x_nuevo - x) / x_nuevo) * 100
        else:
            error = abs(x_nuevo - x)

        fx_nuevo = evaluar_polinomio(coeficientes, x_nuevo)

        # Compruebo si se encontró la raíz o se alcanzó la tolerancia.
        if fx_nuevo == 0 or error <= tolerancia:
            return x_nuevo, i, error, fx_nuevo

        x = x_nuevo

    raise ValueError(
        "Se alcanzó el máximo de iteraciones sin cumplir la tolerancia."
    )



# METODO DE LA SECANTE


def metodo_secante(coeficientes, x0, x1, tolerancia, max_iteraciones):

    # Calculo el valor inicial del polinomio en los dos puntos.
    f0 = evaluar_polinomio(coeficientes, x0)
    f1 = evaluar_polinomio(coeficientes, x1)

    for i in range(1, max_iteraciones + 1):

        # Evito una división por cero.
        if abs(f1 - f0) < 1e-12:
            raise ValueError(
                "No se puede continuar: los valores de la función son iguales."
            )

        # Calculo una nueva aproximación.
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

        # Calculo el error aproximado relativo porcentual.
        if x2 != 0:
            error = abs((x2 - x1) / x2) * 100
        else:
            error = abs(x2 - x1)

        f2 = evaluar_polinomio(coeficientes, x2)

        # Compruebo si se encontró la raíz o se alcanzó la tolerancia.
        if f2 == 0 or error <= tolerancia:
            return x2, i, error, f2

        # Actualizo los puntos para la siguiente iteración.
        x0 = x1
        f0 = f1

        x1 = x2
        f1 = f2

    raise ValueError(
        "Se alcanzó el máximo de iteraciones sin cumplir la tolerancia."
    )



# CALCULO DE TOLERANCIA


def obtener_tolerancia(modo, valor):

    # Si se ingresan cifras significativas, calculo el error permitido.
    if modo == "Cifras significativas":

        cifras = int(valor)

        if cifras < 1:
            raise ValueError("Las cifras significativas deben ser positivas.")

        tolerancia = 0.5 * (10 ** (2 - cifras))

        return tolerancia, cifras

    # Si se ingresa directamente el error porcentual.
    else:

        epsilon = float(valor)

        if epsilon <= 0:
            raise ValueError("El error porcentual debe ser mayor que cero.")

        return epsilon, None



# METODO DE GAUSS


def metodo_gauss(A, b, pivoteo=True):

    # Trabajo con copias para no modificar los datos originales.
    A = [fila[:] for fila in A]
    b = b[:]

    n = len(b)

    if n == 0 or len(A) != n or any(len(fila) != n for fila in A):
        raise ValueError("La matriz debe ser cuadrada y compatible.")

    # Eliminación hacia adelante.
    for k in range(n - 1):

        # Pivoteo parcial: busco el elemento de mayor valor absoluto.
        if pivoteo:

            fila_max = max(
                range(k, n),
                key=lambda i: abs(A[i][k])
            )

            A[k], A[fila_max] = A[fila_max], A[k]
            b[k], b[fila_max] = b[fila_max], b[k]

        # Verifico que el pivote no sea cero.
        if abs(A[k][k]) < 1e-12:
            raise ValueError("Se encontró un pivote nulo.")

        # Hago ceros debajo del pivote.
        for i in range(k + 1, n):

            factor = A[i][k] / A[k][k]

            A[i][k] = 0.0

            for j in range(k + 1, n):
                A[i][j] -= factor * A[k][j]

            b[i] -= factor * b[k]

    # Verifico el último pivote.
    if abs(A[n - 1][n - 1]) < 1e-12:
        raise ValueError(
            "Se encontró un pivote nulo. El sistema no tiene solución única."
        )

    # Sustitución hacia atrás.
    x = [0.0] * n

    for i in range(n - 1, -1, -1):

        suma = sum(
            A[i][j] * x[j]
            for j in range(i + 1, n)
        )

        x[i] = (b[i] - suma) / A[i][i]

    return x



# METODO DE SIMPSON 


def metodo_simpson(funcion, a, b, n):

    # La cantidad de subintervalos debe ser positiva y par.
    if n <= 0 or n % 2 != 0:
        raise ValueError(
            "La cantidad de subintervalos debe ser positiva y par."
        )

    # Calculo el ancho de cada subintervalo.
    h = (b - a) / n

    # Sumo los valores de los extremos.
    suma = funcion(a) + funcion(b)

    # Sumo los términos interiores.
    for i in range(1, n):

        x = a + i * h

        if i % 2 == 0:
            suma += 2 * funcion(x)
        else:
            suma += 4 * funcion(x)

    # Aplico la fórmula de Simpson 1/3.
    integral = (h / 3) * suma

    return integral



# VALORES PREDETERMINADOS DEL LABORATORIO


# Polinomio:
# f(x) = -25 + 82x - 8x^3 + x^5

coeficientes_prueba = [1, 0, -8, 0, 82, -25]

# Matriz del sistema de cuatro ecuaciones.
matriz_prueba = [
    [2,  2,  1, -1],
    [2, -2,  6, -4],
    [2,  2,  1,  1],
    [1,  5, -4, -2]
]

# Vector de términos independientes.
vector_prueba = [5, 0, 13, -9]


# Función para utilizar el polinomio en Simpson.
def funcion_prueba(x):
    return evaluar_polinomio(coeficientes_prueba, x)



# INTERFAZ GRAFICA - TKINTER


import tkinter as tk
from tkinter import ttk, messagebox


def iniciar_interfaz():

    ventana = tk.Tk()
    ventana.title("Laboratorio 6 - Informática II")
    ventana.geometry("850x720")
    ventana.resizable(False, False)

    # FUNCIONES AUXILIARES
   

    def crear_entrada(parent, texto, valor, fila, columna=0):

        tk.Label(parent, text=texto).grid(
            row=fila, column=columna, padx=5, pady=6, sticky="e"
        )

        entrada = tk.Entry(parent, width=18)
        entrada.grid(row=fila, column=columna + 1, padx=5, pady=6)
        entrada.insert(0, str(valor))

        return entrada

    def leer_float(entrada, nombre):

        try:
            return float(entrada.get())
        except ValueError:
            raise ValueError(nombre + " debe ser un número válido.")

    def leer_entero(entrada, nombre):

        try:
            return int(entrada.get())
        except ValueError:
            raise ValueError(nombre + " debe ser un número entero.")

    def leer_coeficientes():

        texto = entrada_coeficientes.get().strip()

        if not texto:
            raise ValueError("Debe ingresar los coeficientes.")

        try:
            coeficientes = [
                float(valor.strip())
                for valor in texto.split(",")
            ]
        except ValueError:
            raise ValueError(
                "Ingrese los coeficientes separados por comas."
            )

        if len(coeficientes) > 7:
            raise ValueError(
                "El polinomio no puede superar el grado 6."
            )

        return coeficientes

    def mostrar_resultado(area, texto):

        area.config(state="normal")
        area.delete("1.0", tk.END)
        area.insert(tk.END, texto)
        area.config(state="disabled")

    def limpiar_resultado(area):

        mostrar_resultado(area, "")

    def manejar_error(error):

        messagebox.showerror("Error", str(error))

  
    # VENTANA PRINCIPAL
   

    tk.Label(
        ventana,
        text="LABORATORIO 6 - CALCULO NUMERICO",
        font=("Arial", 17, "bold")
    ).pack(pady=10)

    tk.Label(
        ventana,
        text="Informática II",
        font=("Arial", 11)
    ).pack()

   
    # INGRESO DEL POLINOMIO
  

    marco_polinomio = tk.LabelFrame(
        ventana,
        text="Polinomio",
        padx=10,
        pady=10
    )
    marco_polinomio.pack(fill="x", padx=15, pady=10)

    tk.Label(
        marco_polinomio,
        text="Coeficientes de mayor a menor grado:"
    ).pack(side="left", padx=5)

    entrada_coeficientes = tk.Entry(
        marco_polinomio,
        width=60
    )
    entrada_coeficientes.pack(side="left", padx=5)

    entrada_coeficientes.insert(
        0, "1, 0, -8, 0, 82, -25"
    )

     # PESTAÑAS DE OPERACIONES
    

    pestañas = ttk.Notebook(ventana)
    pestañas.pack(fill="both", expand=True, padx=15, pady=10)

    pestaña_raices = ttk.Frame(pestañas)
    pestaña_derivada = ttk.Frame(pestañas)
    pestaña_gauss = ttk.Frame(pestañas)
    pestaña_simpson = ttk.Frame(pestañas)

    pestañas.add(pestaña_raices, text="Raíces")
    pestañas.add(pestaña_derivada, text="Derivada")
    pestañas.add(pestaña_gauss, text="Gauss")
    pestañas.add(pestaña_simpson, text="Simpson")

  
    # PESTAÑA 1 - RAICES
    

    marco_raices = tk.LabelFrame(
        pestaña_raices,
        text="Métodos para calcular raíces",
        padx=10,
        pady=10
    )
    marco_raices.pack(fill="x", padx=10, pady=10)

    tk.Label(
        marco_raices,
        text="Método:"
    ).grid(row=0, column=0, padx=5, pady=5)

    combo_metodo = ttk.Combobox(
        marco_raices,
        values=["Bisección", "Newton-Raphson", "Secante"],
        state="readonly",
        width=20
    )
    combo_metodo.grid(row=0, column=1, padx=5, pady=5)
    combo_metodo.current(0)

    entrada_a = crear_entrada(
        marco_raices, "Extremo a:", "0", 1, 0
    )

    entrada_b = crear_entrada(
        marco_raices, "Extremo b:", "1", 1, 2
    )

    entrada_x0 = crear_entrada(
        marco_raices, "Valor x0:", "0", 2, 0
    )

    entrada_x1 = crear_entrada(
        marco_raices, "Valor x1:", "-1", 2, 2
    )

    tk.Label(
        marco_raices,
        text="Condición de parada:"
    ).grid(row=3, column=0, padx=5, pady=5)

    combo_tolerancia = ttk.Combobox(
        marco_raices,
        values=["Error porcentual", "Cifras significativas"],
        state="readonly",
        width=22
    )
    combo_tolerancia.grid(row=3, column=1, padx=5, pady=5)
    combo_tolerancia.current(0)

    entrada_valor_tolerancia = crear_entrada(
        marco_raices, "Valor:", "0.5", 3, 2
    )

    entrada_max_iteraciones = crear_entrada(
        marco_raices, "Máximo iteraciones:", "100", 4, 0
    )

    def calcular_raiz():

        try:

            coeficientes = leer_coeficientes()

            metodo = combo_metodo.get()

            valor = leer_float(
                entrada_valor_tolerancia,
                "La tolerancia"
            )

            if combo_tolerancia.get() == "Cifras significativas":
                modo = "Cifras significativas"
            else:
                modo = "Error porcentual"

            tolerancia, cifras = obtener_tolerancia(
                modo, valor
            )

            max_iteraciones = leer_entero(
                entrada_max_iteraciones,
                "El máximo de iteraciones"
            )

            if max_iteraciones <= 0:
                raise ValueError(
                    "El máximo de iteraciones debe ser positivo."
                )

            if metodo == "Bisección":

                a = leer_float(entrada_a, "El extremo a")
                b = leer_float(entrada_b, "El extremo b")

                resultado = metodo_biseccion(
                    coeficientes,
                    a,
                    b,
                    tolerancia,
                    max_iteraciones
                )

            elif metodo == "Newton-Raphson":

                x0 = leer_float(entrada_x0, "El valor x0")

                resultado = metodo_newton_raphson(
                    coeficientes,
                    x0,
                    tolerancia,
                    max_iteraciones
                )

            else:

                x0 = leer_float(entrada_x0, "El valor x0")
                x1 = leer_float(entrada_x1, "El valor x1")

                resultado = metodo_secante(
                    coeficientes,
                    x0,
                    x1,
                    tolerancia,
                    max_iteraciones
                )

            raiz, iteraciones, error, fx = resultado

            if cifras is not None:
                raiz_formateada = f"{raiz:.{cifras}g}"
            else:
                raiz_formateada = f"{raiz:.10g}"

            texto = (
                "RESULTADOS\n"
                "------------------------------\n"
                f"Método: {metodo}\n"
                f"Raíz encontrada: {raiz_formateada}\n"
                f"Iteraciones: {iteraciones}\n"
                f"Error aproximado relativo: {error:.10g} %\n"
                f"f(xr): {fx:.10g}\n"
                f"Tolerancia: {tolerancia:.10g} %"
            )

            mostrar_resultado(resultado_raices, texto)

        except (ValueError, ZeroDivisionError, OverflowError) as error:
            manejar_error(error)

    tk.Button(
        pestaña_raices,
        text="Calcular raíz",
        command=calcular_raiz,
        width=20
    ).pack(pady=10)

    resultado_raices = tk.Text(
        pestaña_raices,
        height=10,
        width=90,
        state="disabled"
    )
    resultado_raices.pack(padx=10, pady=10)

    
    # PESTAÑA 2 - DERIVADA NUMERICA
  

    marco_derivada = tk.LabelFrame(
        pestaña_derivada,
        text="Diferencia finita centrada",
        padx=10,
        pady=10
    )
    marco_derivada.pack(fill="x", padx=15, pady=15)

    entrada_x = crear_entrada(
        marco_derivada, "Punto x:", "1", 0
    )

    entrada_h = crear_entrada(
        marco_derivada, "Paso h:", "0.1", 1
    )

    def calcular_derivada():

        try:

            coeficientes = leer_coeficientes()

            x = leer_float(entrada_x, "El punto x")
            h = leer_float(entrada_h, "El paso h")

            resultado = derivada_numerica(
                coeficientes, x, h
            )

            texto = (
                "RESULTADO DE LA DERIVADA\n"
                "------------------------------\n"
                f"Punto x: {x}\n"
                f"Paso h: {h}\n"
                f"Derivada aproximada: {resultado:.10g}"
            )

            mostrar_resultado(resultado_derivada, texto)

        except (ValueError, ZeroDivisionError, OverflowError) as error:
            manejar_error(error)

    tk.Button(
        pestaña_derivada,
        text="Calcular derivada",
        command=calcular_derivada,
        width=20
    ).pack(pady=10)

    resultado_derivada = tk.Text(
        pestaña_derivada,
        height=10,
        width=90,
        state="disabled"
    )
    resultado_derivada.pack(padx=10, pady=10)

 
    # PESTAÑA 3 - SISTEMA DE ECUACIONES GAUSS
  

    marco_gauss = tk.LabelFrame(
        pestaña_gauss,
        text="Matriz de coeficientes 4x4",
        padx=10,
        pady=10
    )
    marco_gauss.pack(pady=10)

    entradas_matriz = []

    for i in range(4):

        fila_entradas = []

        for j in range(4):

            entrada = tk.Entry(
                marco_gauss,
                width=9,
                justify="center"
            )

            entrada.grid(
                row=i,
                column=j,
                padx=5,
                pady=5
            )

            entrada.insert(
                0, str(matriz_prueba[i][j])
            )

            fila_entradas.append(entrada)

        entradas_matriz.append(fila_entradas)

    marco_vector = tk.LabelFrame(
        pestaña_gauss,
        text="Vector de términos independientes",
        padx=10,
        pady=10
    )
    marco_vector.pack(pady=5)

    entradas_vector = []

    for i, valor in enumerate(vector_prueba):

        entrada = tk.Entry(
            marco_vector,
            width=9,
            justify="center"
        )

        entrada.grid(
            row=0,
            column=i,
            padx=5,
            pady=5
        )

        entrada.insert(0, str(valor))
        entradas_vector.append(entrada)

    tk.Label(
        pestaña_gauss,
        text="Método de resolución:"
    ).pack(pady=5)

    combo_pivoteo = ttk.Combobox(
        pestaña_gauss,
        values=["Con pivoteo parcial", "Sin pivoteo"],
        state="readonly",
        width=25
    )
    combo_pivoteo.pack()
    combo_pivoteo.current(0)

    def calcular_gauss():

        try:

            A = []
            b = []

            # Leo los valores de la matriz.
            for i in range(4):

                fila = []

                for j in range(4):

                    valor = leer_float(
                        entradas_matriz[i][j],
                        "La matriz"
                    )

                    fila.append(valor)

                A.append(fila)

            # Leo el vector independiente.
            for entrada in entradas_vector:

                valor = leer_float(
                    entrada,
                    "El vector independiente"
                )

                b.append(valor)

            pivoteo = (
                combo_pivoteo.get() == "Con pivoteo parcial"
            )

            solucion = metodo_gauss(A, b, pivoteo)

            texto = (
                "RESULTADO DEL SISTEMA\n"
                "------------------------------\n"
                f"x1 = {solucion[0]:.10g}\n"
                f"x2 = {solucion[1]:.10g}\n"
                f"x3 = {solucion[2]:.10g}\n"
                f"x4 = {solucion[3]:.10g}\n"
                f"\nMétodo: {combo_pivoteo.get()}"
            )

            mostrar_resultado(resultado_gauss, texto)

        except (ValueError, ZeroDivisionError, OverflowError) as error:
            manejar_error(error)

    tk.Button(
        pestaña_gauss,
        text="Resolver sistema",
        command=calcular_gauss,
        width=20
    ).pack(pady=10)

    resultado_gauss = tk.Text(
        pestaña_gauss,
        height=9,
        width=90,
        state="disabled"
    )
    resultado_gauss.pack(padx=10, pady=10)

    
    # PESTAÑA 4 - INTEGRAL SIMPSON 1/3
    

    marco_simpson = tk.LabelFrame(
        pestaña_simpson,
        text="Integral definida",
        padx=10,
        pady=10
    )
    marco_simpson.pack(fill="x", padx=15, pady=15)

    entrada_integral_a = crear_entrada(
        marco_simpson, "Límite inferior a:", "0", 0
    )

    entrada_integral_b = crear_entrada(
        marco_simpson, "Límite superior b:", "2", 1
    )

    entrada_n = crear_entrada(
        marco_simpson, "Cantidad de segmentos n:", "2", 2
    )

    def calcular_integral():

        try:

            coeficientes = leer_coeficientes()

            a = leer_float(
                entrada_integral_a,
                "El límite inferior"
            )

            b = leer_float(
                entrada_integral_b,
                "El límite superior"
            )

            n = leer_entero(
                entrada_n,
                "La cantidad de segmentos"
            )

            if a >= b:
                raise ValueError(
                    "El límite inferior debe ser menor que el superior."
                )

            # Defino el polinomio que se va a integrar.
            def funcion(x):
                return evaluar_polinomio(coeficientes, x)

            resultado = metodo_simpson(
                funcion, a, b, n
            )

            texto = (
                "RESULTADO DE LA INTEGRAL\n"
                "------------------------------\n"
                f"Límite inferior: {a}\n"
                f"Límite superior: {b}\n"
                f"Cantidad de segmentos: {n}\n"
                f"Integral aproximada: {resultado:.10g}"
            )

            mostrar_resultado(resultado_simpson, texto)

        except (ValueError, ZeroDivisionError, OverflowError) as error:
            manejar_error(error)

    tk.Button(
        pestaña_simpson,
        text="Calcular integral",
        command=calcular_integral,
        width=20
    ).pack(pady=10)

    resultado_simpson = tk.Text(
        pestaña_simpson,
        height=10,
        width=90,
        state="disabled"
    )
    resultado_simpson.pack(padx=10, pady=10)

   
    # PIE DE VENTANA
    

    tk.Label(
        ventana,
        text="Laboratorio 6 - Informática II",
        font=("Arial", 9)
    ).pack(pady=5)

    ventana.mainloop()



# EJECUCION DEL PROGRAMA


if __name__ == "__main__":
    iniciar_interfaz()

