# jefe.cpp – Guía rápida

## Propósito
Actúa como intermediario: recibe la matriz del cliente, consulta disponibilidad de workers al server, reenvía la matriz y parámetros, y devuelve U, S, V^T al cliente. También genera la semilla si el cliente no la fija.

## Flujo de `handle_client`
1. Puede recibir primero un `ID_H` del cliente para consultar disponibilidad:
   - Se conecta al server, envía `ID_H`, recibe `available` y lo reenvía al cliente.
   - Luego espera el verdadero `ID_A` del cliente.
2. Recibe `ID_A` (n, k_full) y la matriz completa; la almacena en `jefe_matrix_<cid>.bin` con `mmap`.
3. Recibe opcional `ID_K` (k objetivo sin oversampling) del cliente.
4. Conecta al server (`BOSS_SERVER_PORT`).
5. Handshake de disponibilidad con el server (`ID_H`, `a = min_workers`):
   - Mínimo interno = 1 (el cliente decide si procede con 1 o espera).
   - Si el server responde menos de `min_workers`, reintenta 3 veces con 5s; si falla, envía `DONE` al cliente y termina.
6. Reenvía al server:
   - `ID_A` con n, k_full y la matriz (stream en `CHUNK`).
   - `ID_K` con k_target.
   - Semilla `ID_S` (generada con `choose_seed`).
7. Reenvía resultados del server al cliente:
   - Pasa headers y payloads para `ID_UT`, `ID_S`, `ID_VT`, `ID_DONE`.

## Funciones clave
- `conectar_servidor_trabajo()`: conecta al server de workers.
- `reenviar_matriz(serverSock, mm, n, k)`: manda header `ID_A` y la matriz.
- `retransmitir_resultados(serverSock, clientSock)`: túnel U, S, V^T y `DONE`.
- `elegir_semilla()`: genera semilla aleatoria de 64 bits.

## Protocolos
- `ID_H`: consulta workers disponibles.
- `ID_A`: dimensiones y matriz A.
- `ID_K`: k objetivo sin oversampling.
- `ID_S`: semilla.
- `ID_UT`, `ID_VT`, `ID_S`, `ID_DONE`: resultados.
