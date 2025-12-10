# client.cpp – Guía rápida

## Propósito
Cliente interactivo: genera una matriz `N x N`, la guarda en `matrix.bin`, y la envía al jefe para obtener la SVD distribuida. Recibe U, S, V^T y los guarda en disco.

## Flujo principal
- Menú:
  - `1) Generar matriz aleatoria`: crea `matrix.bin` con `N x N` floats uniformes en [-1,1], semilla variable, actualiza `N_global`.
  - `2) Cargar matriz desde binario`: lee un archivo cuadrado y lo copia a `matrix.bin`, actualiza `N_global`.
  - `3) Enviar matriz`: pide `k`, `p`, calcula `k_total = k+p`, y ejecuta `enviar_matriz_y_recibir_svd`.
  - `4) View matrix head`: muestra cabeceras de `matrix.bin`, `U.bin`, `VT.bin`, `S.bin` (si existen).
  - `5) Exit`.

## `send_matrix_and_receive_svd(k_total, k_target)`
- Pre‑check de workers (mínimo 2):
  - Envía `ID_H`; si hay 1 worker, pregunta si desea continuar; reintenta hasta 3 veces con espera de 5s si hay menos de 2.
  - Si no hay suficientes, aborta sin enviar la matriz.
- Envío:
  - Header `ID_A` con `a = N_global`, `b = k_total`.
  - Envía `matrix.bin` en bloques de tamaño `CHUNK`.
  - Envía `ID_K` con `k_target` (k deseado sin oversampling).
- Recepción:
  - Espera secuencia: `ID_UT` (U), `ID_S` (Sigma), `ID_VT` (V^T), `ID_DONE`.
  - Cada bloque recibido se guarda en `U.bin`, `S.bin`, `VT.bin` usando `mmap`/`ftruncate`.

## Utilidades
- `generar_matriz(N)`: crea `matrix.bin` con floats aleatorios uniformes en [-1,1], semilla variable.
- `cargar_matriz_desde_bin(ruta)`: valida que el binario sea cuadrado y lo copia a `matrix.bin`.
- `mostrar_matriz/ mostrar_sigma`: leen con `mmap`, si el archivo no existe sólo informan.

## Protocolos usados
- `ID_H`: consulta disponibilidad de workers al jefe/server.
- `ID_A`: anuncia envío de A (dimensiones n, k_total).
- `ID_K`: k objetivo sin oversampling.
- `ID_UT`, `ID_S`, `ID_VT`, `ID_DONE`: resultados finales desde el server.
