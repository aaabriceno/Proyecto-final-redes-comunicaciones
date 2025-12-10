# protocolo.hpp – Guía rápida

## Propósito
Constantes de puertos, tamaño de chunk, IDs de mensajes y utilidades de envío/recepción.

## Constantes
- Puertos: `CLIENT_PORT=45000`, `WORKER_PORT=45001`, `BOSS_SERVER_PORT=45002`.
- `CHUNK = 64 KB`.
- IDs: `ID_K`, `ID_H`, `ID_A`, `ID_S`, `ID_R`, `ID_Q`, `ID_B`, `ID_D`, `ID_C`, `ID_UT`, `ID_V`, `ID_VT`, `ID_UI`, `ID_DONE`.

## Struct
- `MsgHeader { char id; uint64_t a; uint64_t b; }` (helper ctor incluido).

## Funciones inline
- `send_all(sock, buf, len)`: envía hasta completar o falla.
- `recv_all(sock, buf, len)`: recibe hasta completar o falla.
- `htonll_u64/ntohll_u64`: conversión endian 64 bits.
