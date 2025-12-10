# algebra_matriz.hpp – Guía rápida

## Propósito
Rutinas numéricas usadas por workers/server para RSVD.

## Funciones
- `generar_omega_mmap(seed, n, k, fileOmega)`: crea matriz aleatoria `n x k` en disco.
- `multiplicar_A_por_Omega_mmap(fileA, rows, cols, fileOmega, k, fileY)`: `Y = A * Omega` en mmap.
- `descomponer_qr_mmap(fileY, rows, k, fileQi, fileRi)`: QR (Gram-Schmidt) sobre `Y` → `Qi`, `Ri`.
- `jacobi_autovalores_inplace(A, n, Util, Lambda, max_iter, tol)`: eigen-decomp simétrica (Jacobi) in-place sobre copia de `A`.
- `calcular_Vj_mmap(UtilFile, SigmaInvFile, BiFile, VjOutFile, k, ncols)`: `V_j = SigmaInv * Util^T * B_i` (almacena k x ncols).
- `calcular_Ui_mmap(Qi_file, Utilde_file, Ui_file, mi, k)`: `U_i = Q_i * Utilde`.

## Notas
- Todas operan con archivos mapeados; asegúrate de cerrar con `mmap_close`.
- Supone `float` y matrices densas.
