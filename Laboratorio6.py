

# Laboratorio 6 - Informática II
# Cálculo numérico
# Repositorio GitHub:
# https://github.com/seba5680181-eng/Laboratorio6


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
