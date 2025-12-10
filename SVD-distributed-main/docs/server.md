# server/server.cpp – Guía rápida

## Propósito
Orquestador distribuido: recibe A desde el jefe, reparte trabajo a los workers, ensambla resultados y devuelve U, S, V^T al jefe. Usa `mmap` para manejar matrices en disco.

## Flujo de `clientHandler`
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
- `worker_acceptor_thread()`: acepta conexiones de workers y guarda sockets.
- `send_Ai_to_worker`, `send_seed_to_worker`.
- `recv_R_from_worker_mmap`, `build_Rstack_mmap`, `qr_mmap` (TSQR).
- `send_Qr_to_worker`.
- `recv_Bi_from_worker_mmap`, `assemble_B_mmap`.
- `send_B_block_to_worker_mmap` (Bj por columnas).
- `recv_Cj_from_worker_mmap`, ensamblar C con `assemble_B_mmap`.
- `eigendecompose_C_mmap`, `sigma_and_inv_mmap`, `send_Util_Sinv_to_worker_mmap`.
- `recv_Vj_from_worker_mmap`, `assemble_Vt_mmap`.
- `recv_Ui_from_worker_mmap`, `assemble_U_mmap`.
- `send_Vt_to_client`, `send_Sigma_mmap`.

## Protocolos y tablas
- Usa `print_protocol_table` para mostrar campos de mensajes (A, S, R_i, B_i, V_j, U_i, finales).

## Notas
- A se materializa en disco (`server_matrix.bin`) vía `mmap`.
- Abortará si no hay workers (`W==0`) tras recibir headers.
- Usa `CHUNK` para envíos/recepciones grandes.
