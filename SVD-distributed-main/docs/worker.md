# worker/worker.cpp – Guía rápida

## Propósito
Cada worker recibe su bloque A_i, ejecuta pasos locales de RSVD y devuelve resultados parciales (R_i, B_i, C_j, V_j, U_i) al server.

## Flujo en `handle_connection`
- Espera headers del server en bucle.
- `ID_A`: recibe Ai (rows_i x n) en `Ai.bin`; luego recibe `ID_S` (seed, k).
  - Genera `Omega`, calcula `Y = A_i * Omega`.
  - QR local → `Qi.bin`, `Ri.bin`; envía `R_i` (`ID_R`).
  - Recibe `Qr` (`ID_Q`) desde server; forma `Q_final = Qi * Qr`.
  - Calcula `B_i = Q_final * A_i`; envía `B_i` (`ID_B`).
  - Recibe bloque Bj (`ID_D`) con columnas asignadas.
  - Calcula `C_j = Bj * Bj^T`; envía `C_j` (`ID_C`).
  - Recibe `Utilde` y `SigmaInv`; calcula `V_j` y lo envía (payload, ID_VT).
  - Calcula `U_i = Q_final * Utilde`; envía `U_i` (`ID_UI`).
- Logs informan cada etapa y archivos intermedios.

## Archivos intermedios (prefijo `worker_<pid>_`)
- `Ai.bin`, `Omega.bin`, `Y.bin`, `Qi.bin`, `Ri.bin`, `Qr.bin`, `Qfinal.bin`, `Bi.bin`, `Bj.bin`, `Cj.bin`, `Utilde.bin`, `SigmaInv.bin`, `Vj.bin`, `Ui.bin`.

## Protocolos
- `ID_A` (Ai), `ID_S` (seed/k), `ID_R`, `ID_Q`, `ID_B`, `ID_D`, `ID_C`, `ID_VT` (payload V_j), `ID_UI`.
