# server/server.cpp – Guía rápida

## Propósito
Orquestador distribuido: recibe A desde el jefe, reparte trabajo a los workers, ensambla resultados y devuelve U, S, V^T al jefe. Usa `mmap` para manejar matrices en disco.

## Flujo de `manejar_cliente`
1. Handshake opcional: si llega `ID_H`, responde con número de workers conectados.
2. Recibe `ID_A` con `n`, `k_full` (k+p). Mappea `server_matrix.bin`, recibe A completa en disco.
3. Recibe opcional `ID_K` (k objetivo) y `ID_S` (semilla; si no, usa valor por defecto).
4. Toma la lista de workers (`worker_sockets`). Si `W==0`, aborta y libera recursos.
5. Particiona filas/columnas entre `W` workers: `start/nrows` y `start_cols/ncols`.
6. Envíos y cómputo:
   - Envía Ai (por filas) a cada worker (`ID_A`), luego semilla y k (`ID_S`).
   - Recibe R_i (`ID_R`), hace TSQR → `Qr.bin` y `Rglobal`.
   - Envía Qr_i (`ID_Q`) a cada worker.
   - Recibe B_i (`ID_B`) y arma `B_final` (k x n).
   - Envía bloques Bj por columnas (`ID_D`).
   - Recibe C_j (`ID_C`) y arma `C_final` (k x k).
   - Eigendecompone C → Utilde, Lambda; calcula Sigma y SigmaInv.
   - Envía Utilde y SigmaInv a workers; recibe V_j (`ID_VT` payload) y arma `Vt_final`.
   - Recibe U_i (`ID_UI`) y arma `U_final`.
7. Trunca a `k_req` si `k_target < k_full` (genera U_k, Vt_k, Sigma_k).
8. Envía al cliente U (`ID_UT`), S (`ID_S`), V^T (`ID_VT`), y `ID_DONE`.

## Funciones principales
- `hilo_aceptador_workers()`: acepta conexiones de workers y guarda sockets.
- `enviar_Ai_a_worker`, `enviar_semilla_a_worker`.
- `recibir_R_de_worker_mmap`, `construir_Rstack_mmap`, `descomponer_qr_mmap` (TSQR).
- `enviar_Qr_a_worker`.
- `recibir_Bi_de_worker_mmap`, `ensamblar_B_mmap`.
- `enviar_bloque_B_a_worker_mmap` (Bj por columnas).
- `recibir_Cj_de_worker_mmap`, ensamblar C con `ensamblar_B_mmap`.
- `descomponer_C_mmap`, `calcular_sigma_y_inv_mmap`, `enviar_Util_Sinv_a_worker_mmap`.
- `recibir_Vj_de_worker_mmap`, `ensamblar_Vt_mmap`.
- `recibir_Ui_de_worker_mmap`, `ensamblar_U_mmap`.
- `enviar_Vt_a_cliente`, `enviar_Sigma_mmap`.

## Protocolos y tablas
- Usa `print_protocol_table` para mostrar campos de mensajes (A, S, R_i, B_i, V_j, U_i, finales).

## Notas
- A se materializa en disco (`server_matrix.bin`) vía `mmap`.
- Abortará si no hay workers (`W==0`) tras recibir headers.
- Usa `CHUNK` para envíos/recepciones grandes.
