#ifndef PROTOCOLO_HPP
#define PROTOCOLO_HPP

#include <arpa/inet.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <cstdint>
#include <cstring>
#include <iostream>

static const int CLIENT_PORT = 45000;
static const int WORKER_PORT = 45001;
static const int BOSS_SERVER_PORT = 45002;
static const size_t CHUNK = 64*1024; // 64 KB

// message IDs (single char)
#define ID_K 'K'  // client -> server: target k (without oversampling)
#define ID_H 'H'  // handshake: availability / worker count
#define ID_A 'A'  // server -> worker: block Ai
#define ID_S 'S'  // server -> worker: seed + k, and Sigma for client
#define ID_R 'R'  // worker -> server: R_i (k x k)
#define ID_Q 'Q'  // server -> worker: Qr_i (k x k)
#define ID_B 'b'  // worker -> server: B_i (k x n) OR server->worker block
#define ID_D 'D'  // server -> worker: send block of B (k x colsblock)
#define ID_C 'c'  // worker -> server: C_j (k x k)
#define ID_UT 'U' // server -> worker: Utilde (k x k) and Sigma^{-1} (k)
#define ID_V 'v'  // worker -> server: V_j^T (k x colsblock)
#define ID_VT 'V'  // server -> client: V^T (k x n)
#define ID_UI 'u' // worker -> server: U_i (rows_i x k)
#define ID_DONE 'X' // generic done

struct MsgHeader {
    char id;
    uint64_t a;
    uint64_t b;

    MsgHeader(char id_ = 0, uint64_t a_ = 0, uint64_t b_=0): id(id_), a(a_), b(b_) {}
};

inline bool send_all(int sock, const void* buf, size_t len) {
    const char* p = (const char*)buf;
    while (len > 0) {
        ssize_t n = send(sock, p, len, 0);
        if (n <= 0) return false;
        p += n; len -= n;
    }
    return true;
}

inline bool recv_all(int sock, void* buf, size_t len) {
    char* p = (char*)buf;
    while (len > 0) {
        ssize_t n = recv(sock, p, len, 0);
        if (n <= 0) return false;
        p += n; len -= n;
    }
    return true;
}

inline uint64_t htonll_u64(uint64_t v) {
    static const int num = 42;
    if (*(const char*)&num == 42) {
        return ((uint64_t)htonl((uint32_t)(v & 0xffffffffULL)) << 32) | htonl((uint32_t)(v >> 32));
    } else return v;
}

inline uint64_t ntohll_u64(uint64_t v) { return htonll_u64(v); }

#endif
