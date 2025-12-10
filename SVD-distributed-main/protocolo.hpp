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

#define ID_K 'K'  
#define ID_H 'H'  
#define ID_A 'A'  
#define ID_S 'S'  
#define ID_R 'R' 
#define ID_Q 'Q'  
#define ID_B 'b'  
#define ID_D 'D'  
#define ID_C 'c'  
#define ID_UT 'U' 
#define ID_V 'v'  
#define ID_VT 'V' 
#define ID_UI 'u' 
#define ID_DONE 'X' 

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
