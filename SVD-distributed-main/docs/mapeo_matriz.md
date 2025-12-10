# mapeo_matriz.hpp – Guía rápida

## Propósito
Helpers para mapear matrices en archivos usando `mmap` (floats).

## Struct
- `MMapMatrix { float* data; uint64_t filas, columnas; size_t bytes; int fd; string path; }`

## Funciones
- `mmap_create(path, filas, columnas)`: crea/trunca archivo, `ftruncate` al tamaño, `mmap` RW compartido. Lanza `runtime_error` si falla.
- `mmap_open_read(path, filas, columnas)`: abre archivo solo lectura y lo mapea RO compartido.
- `mmap_close(MMapMatrix&)`: `munmap` + `close(fd)`.

## Uso típico
```cpp
auto mm = mmap_create("A.bin", n, n);
// escribir mm.data ...
mmap_close(mm);
```
