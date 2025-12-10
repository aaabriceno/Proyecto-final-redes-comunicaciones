
#include "../protocolo.hpp"
#include "../mapeo_matriz.hpp"
#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <vector>

using namespace std;

namespace {
atomic<uint64_t> client_counter{0};

int conectar_servidor_trabajo() {
    int s = socket(AF_INET, SOCK_STREAM, 0);
    if (s < 0) { perror("Conexion jefe --> servidor"); return -1; }

    sockaddr_in serv{};
    serv.sin_family = AF_INET;
    serv.sin_port = htons(BOSS_SERVER_PORT);
    inet_pton(AF_INET, "127.0.0.1", &serv.sin_addr);
    if (connect(s, (sockaddr*)&serv, sizeof(serv)) < 0) {
        perror("Conectar jefe --> servidor");
        close(s);
        return -1;
    }
    return s;
}

bool retransmitir_payload(int src, int dst, size_t bytes) {
    vector<char> buf(min(CHUNK, bytes > 0 ? bytes : CHUNK));
    size_t done = 0;
    while (done < bytes) {
        size_t chunk = min(buf.size(), bytes - done);
        if (!recv_all(src, buf.data(), chunk)) return false;
        if (!send_all(dst, buf.data(), chunk)) return false;
        done += chunk;
    }
    return true;
}

bool reenviar_matriz(int serverSock, MMapMatrix& mm, uint64_t n, uint64_t k) {
    MsgHeader h(ID_A, n, k);
    if (!send_all(serverSock, &h, sizeof(h))) return false;

    size_t bytes = mm.bytes;
    size_t sent = 0;
    const char* ptr = reinterpret_cast<const char*>(mm.data);
    while (sent < bytes) {
        size_t chunk = min(CHUNK, bytes - sent);
        if (!send_all(serverSock, ptr + sent, chunk)) return false;
        sent += chunk;
    }
    return true;
}

bool retransmitir_resultados(int serverSock, int clientSock) {
    while (true) {
        MsgHeader h;
        if (!recv_all(serverSock, &h, sizeof(h))) {
            cerr << "[jefe] fallo recibiendo encabezado desde servidor\n";
            return false;
        }
        if (!send_all(clientSock, &h, sizeof(h))) return false;

        size_t payload = 0;
        if (h.id == ID_UT || h.id == ID_VT) {
            payload = (size_t)h.a * h.b * sizeof(float);
        } else if (h.id == ID_S) {
            payload = (size_t)h.a * sizeof(float);
        } else if (h.id == ID_DONE) {
            break;
        } else {
            cerr << "[jefe] mensaje desconocido del servidor id=" << h.id << "\n";
            return false;
        }

        if (payload && !retransmitir_payload(serverSock, clientSock, payload))
            return false;
    }
    return true;
}

uint64_t elegir_semilla() {
    uint64_t seed = (uint64_t)
        chrono::steady_clock::now().time_since_epoch().count();
    seed ^= ((uint64_t)random_device{}()) + 0x9e3779b97f4a7c15ULL
            + (seed << 6) + (seed >> 2);
    return seed;
}
}

void manejar_cliente(int cs) {
    const uint64_t min_workers = 1;
    auto send_done = [&](int sock){
        MsgHeader d(ID_DONE,0,0);
        send_all(sock, &d, sizeof(d));
    };

    uint64_t cid = client_counter.fetch_add(1, memory_order_relaxed);
    MsgHeader h;
    if (!recv_all(cs, &h, sizeof(h))) {
        cerr << "[jefe] encabezado invalido del cliente\n";
        close(cs);
        return;
    }

    if (h.id == ID_H) {
        int ss_chk = conectar_servidor_trabajo();
        if (ss_chk < 0) { send_done(cs); close(cs); return; }
        MsgHeader req(ID_H, 0, 0), resp;
        if (!send_all(ss_chk, &req, sizeof(req)) ||
            !recv_all(ss_chk, &resp, sizeof(resp)) ||
            resp.id != ID_H) {
            cerr << "[jefe] fallo en handshake de disponibilidad con el server\n";
            send_done(cs);
            close(ss_chk);
            close(cs);
            return;
        }
        send_all(cs, &resp, sizeof(resp));
        close(ss_chk);

        if (!recv_all(cs, &h, sizeof(h)) || h.id != ID_A) {
            cerr << "[jefe] encabezado invalido del cliente tras handshake\n";
            send_done(cs);
            close(cs);
            return;
        }
    } else if (h.id != ID_A) {
        cerr << "[jefe] encabezado invalido del cliente\n";
        close(cs);
        return;
    }
    uint64_t n = h.a;
    uint64_t k_full = h.b; // k + p
    cerr << "[jefe] cliente " << cid << " n=" << n << " k_full=" << k_full << "\n";

    string matrixFile = "jefe_matrix_" + to_string(cid) + ".bin";
    size_t bytes = (size_t)n * n * sizeof(float);
    MMapMatrix mm = mmap_crear(matrixFile, n, n);

    size_t rec = 0;
    char* base = reinterpret_cast<char*>(mm.data);
    while (rec < bytes) {
        size_t chunk = min(CHUNK, bytes - rec);
        if (!recv_all(cs, base + rec, chunk)) {
            cerr << "[jefe] fallo recibiendo matriz del cliente\n";
            mmap_cerrar(mm);
            close(cs);
            unlink(matrixFile.c_str());
            return;
        }
        rec += chunk;
    }
    cerr << "[cliente->jefe] recibida matriz (" << bytes << " bytes)\n";

    uint64_t k_target = k_full;
    MsgHeader hk_from_client;
    if (recv_all(cs, &hk_from_client, sizeof(hk_from_client)) && hk_from_client.id == ID_K) {
        k_target = hk_from_client.a;
    } else {
        cerr << "[jefe] advertencia: cliente no envió ID_K, usando k_full\n";
    }

    int ss = conectar_servidor_trabajo();
    if (ss < 0) {
        mmap_cerrar(mm);
        close(cs);
        unlink(matrixFile.c_str());
        return;
    }

    MsgHeader availReq(ID_H, min_workers, 0);
    MsgHeader availResp{};
    const int max_attempts = 3;
    int attempt = 0;
    bool enough_workers = false;
    while (attempt < max_attempts) {
        if (!send_all(ss, &availReq, sizeof(availReq)) ||
            !recv_all(ss, &availResp, sizeof(availResp)) ||
            availResp.id != ID_H) {
            cerr << "[jefe] fallo en handshake de disponibilidad con el server\n";
            send_done(cs);
            mmap_cerrar(mm);
            close(ss);
            close(cs);
            unlink(matrixFile.c_str());
            return;
        }
        if (availResp.a >= min_workers) {
            enough_workers = true;
            break;
        }
        ++attempt;
        if (attempt < max_attempts) {
            cerr << "[jefe] workers disponibles ("<<availResp.a<<") < min ("<<min_workers<<"), reintentando en 5s...\n";
            this_thread::sleep_for(std::chrono::seconds(5));
        }
    }
    if (!enough_workers) {
        cerr << "[jefe] no se alcanzó el mínimo de workers tras "<<max_attempts<<" intentos (disp="<<availResp.a<<")\n";
        send_done(cs);
        mmap_cerrar(mm);
        close(ss);
        close(cs);
        unlink(matrixFile.c_str());
        return;
    }

    if (!reenviar_matriz(ss, mm, n, k_full)) {
        cerr << "[jefe] fallo enviando matriz al servidor\n";
        send_done(cs);
        mmap_cerrar(mm);
        close(ss);
        close(cs);
        unlink(matrixFile.c_str());
        return;
    }

    // Enviar k solicitado al server
    MsgHeader hk(ID_K, k_target, 0);
    send_all(ss, &hk, sizeof(hk));

    uint64_t seed = elegir_semilla();
    MsgHeader seedMsg(ID_S, seed, k_full);
    send_all(ss, &seedMsg, sizeof(seedMsg));
    cerr << "[jefe->server] enviado seed=" << seed << "\n";

    mmap_cerrar(mm);
    unlink(matrixFile.c_str());

    if (!retransmitir_resultados(ss, cs)) {
        cerr << "[jefe] relay fallido para cliente " << cid << "\n";
    } else {
        cerr << "[jefe] completado cliente " << cid << "\n";
    }

    close(ss);
    close(cs);
}

int main() {
    int ls = socket(AF_INET, SOCK_STREAM, 0);
    if (ls < 0) { perror("socket jefe"); return 1; }
    int opt = 1;
    setsockopt(ls, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(CLIENT_PORT);
    addr.sin_addr.s_addr = INADDR_ANY;
    if (bind(ls, (sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("bind jefe");
        return 1;
    }
    if (listen(ls, 16) < 0) { perror("listen jefe"); return 1; }
    cerr << "[jefe] escuchando clientes en puerto " << CLIENT_PORT << "\n";

    while (true) {
        int cs = accept(ls, nullptr, nullptr);
        if (cs < 0) { perror("accept jefe"); continue; }
        thread(manejar_cliente, cs).detach();
    }

    return 0;
}
