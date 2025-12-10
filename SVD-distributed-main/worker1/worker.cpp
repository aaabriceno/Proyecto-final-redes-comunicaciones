#include "../protocolo.hpp"
#include "../mapeo_matriz.hpp"
#include "../algebra_matriz.hpp"
#include <arpa/inet.h>
#include <sys/socket.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <thread>
#include <fcntl.h>
#include <sys/mman.h>

using namespace std;

void manejar_conexion(int sock) {
    const std::string prefix = "worker_" + std::to_string(getpid()) + "_";
    auto fn = [&](const std::string &name){ return prefix + name; };
    while (true) {
        MsgHeader h;
        if (!recv_all(sock, &h, sizeof(h))) { 
            cout << "[worker] server closed\n"; 
            break; 
        }
        if (h.id == ID_A) {
            uint64_t rows_i = h.a;
            uint64_t n = h.b;
            cout << "[server->worker] Received Ai dims: " << rows_i << " x " << n << "\n";

            string fname = fn("Ai.bin");
            size_t bytes = rows_i * n * sizeof(float);
            int fd = open(fname.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
            ftruncate(fd, bytes);
            void* mapv = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            float* M = (float*)mapv;

            size_t rec = 0;
            char* base = (char*)mapv;
            while (rec < bytes) {
                size_t chunk = (bytes - rec) > CHUNK ? CHUNK : (bytes - rec);
                if (!recv_all(sock, base + rec, chunk)) {
                    cerr << "recv matrix failed\n";
                    munmap(mapv, bytes); close(fd); close(sock); return; 
                }
                rec += chunk;
            }

            MsgHeader hs;
            recv_all(sock, &hs, sizeof(hs));
            uint64_t seed = hs.a;
            uint64_t k = hs.b;
            cout << "[server->worker] seed="<<seed<<" k="<<k<<"\n";

            generar_omega_mmap(seed, n, k, fn("Omega.bin"));
            multiplicar_A_por_Omega_mmap(fn("Ai.bin"), rows_i, n, fn("Omega.bin"), k, fn("Y.bin"));

            descomponer_qr_mmap(fn("Y.bin"), rows_i, k, fn("Qi.bin"), fn("Ri.bin"));

            {
                auto Rm = mmap_abrir_lectura(fn("Ri.bin"), k, k);
                MsgHeader hr(ID_R, k, k);
                send_all(sock, &hr, sizeof(hr));
                send_all(sock, Rm.data, (size_t)k * k * sizeof(float));
                mmap_cerrar(Rm);
                cout << "[worker->server] sent R_i to server\n";
            }

            {
                MsgHeader hq;
                recv_all(sock, &hq, sizeof(hq));
                auto Qr = mmap_crear(fn("Qr.bin"), hq.a, hq.b);
                recv_all(sock, Qr.data, sizeof(float)*hq.a*hq.b);
                mmap_cerrar(Qr);
                cout << "[server->worker] received Qr\n";
            }

            {
                auto Qi  = mmap_abrir_lectura(fn("Qi.bin"), rows_i, k);
                auto Qr  = mmap_abrir_lectura(fn("Qr.bin"), k, k);
                auto Qf  = mmap_crear(fn("Qfinal.bin"), rows_i, k);

                for (uint64_t r = 0; r < rows_i; ++r)
                    for (uint64_t c = 0; c < k; ++c) {
                        float sum = 0;
                        for (uint64_t kk = 0; kk < (uint64_t)k; ++kk)
                            sum += Qi.data[r*(uint64_t)k + kk] * Qr.data[kk*(uint64_t)k + c];

                        Qf.data[r*(uint64_t)k + c] = sum;
                    }
                mmap_cerrar(Qi);
                mmap_cerrar(Qr);
                mmap_cerrar(Qf);
                
                cout << "[worker] Calculate Q_final\n";
            }

            {
                auto Qf = mmap_abrir_lectura(fn("Qfinal.bin"), rows_i, k);
                auto Ai = mmap_abrir_lectura(fn("Ai.bin"), rows_i, n);
                auto Bi = mmap_crear(fn("Bi.bin"), k, n);

                for (uint64_t i = 0; i < (uint64_t)k; ++i)
                    for (uint64_t j = 0; j < (uint64_t)n; ++j) {
                        float sum = 0.0f;
                        for (uint64_t t = 0; t < (uint64_t)rows_i; ++t)
                            sum += Qf.data[t*(uint64_t)k + i] * Ai.data[t*(uint64_t)n + j];
                        Bi.data[i*(uint64_t)n + j] = sum;
                    }
                mmap_cerrar(Qf);
                mmap_cerrar(Ai);
                mmap_cerrar(Bi);

                cout << "[worker] Calculate B_i = Q_final x A_i\n";
            }

            {
                auto Bi = mmap_abrir_lectura(fn("Bi.bin"), k, n);
                MsgHeader hb;
                hb.id = ID_B;
                hb.a = k;
                hb.b = n;

                send_all(sock, &hb, sizeof(hb));
                send_all(sock, Bi.data, sizeof(float)*k*n);
                mmap_cerrar(Bi);

                cout << "[worker->server] sent B_i to server\n";
            }
            
            int cols_j = 0;
            {
                MsgHeader hd;
                recv_all(sock, &hd, sizeof(hd));
                cols_j = static_cast<int>(hd.b);
                auto Bj = mmap_crear(fn("Bj.bin"), k, cols_j);
                recv_all(sock, Bj.data, sizeof(float) * k * cols_j);
                mmap_cerrar(Bj);
                cout << "[server->worker] received B_j from server (" << cols_j << " cols)\n";                
            }

            {
                auto Bj = mmap_abrir_lectura(fn("Bj.bin"), k, cols_j);
                auto Cj = mmap_crear(fn("Cj.bin"), k, k);
                for (int i = 0; i < k; ++i) {
                    for (int j = 0; j < k; ++j) {
                        float sum = 0.0f;
                        for (int c = 0; c < cols_j; ++c) {
                            float bi = Bj.data[(uint64_t)i * cols_j + c];
                            float bj = Bj.data[(uint64_t)j * cols_j + c];
                            sum += bi * bj;
                        }

                        Cj.data[(uint64_t)i * k + j] = sum;
                    }
                }

                mmap_cerrar(Bj);
                mmap_cerrar(Cj);

                cout << "[worker] Calculate C_j (k x k)\n";  
            }

            {
                MsgHeader h(ID_C, k, k);
                send_all(sock, &h, sizeof(h));

                auto Cj = mmap_abrir_lectura(fn("Cj.bin"), k, k);
                send_all(sock, Cj.data, sizeof(float) * k * k);
                mmap_cerrar(Cj);

                cout << "[worker->server] sent C_j to server\n";
            }
            
            {
                auto Utilmm = mmap_crear(fn("Utilde.bin"), k, k);
                recv_all(sock, Utilmm.data, sizeof(float) * k * k);

                auto SigmaInvmm = mmap_crear(fn("SigmaInv.bin"), k, 1);
                recv_all(sock, SigmaInvmm.data, sizeof(float) * k);

                mmap_cerrar(Utilmm);
                mmap_cerrar(SigmaInvmm);
                cout << "[server->worker] received Utilde and SigmaInv from server\n";           
            }

            calcular_Vj_mmap(fn("Utilde.bin"),fn("SigmaInv.bin"),fn("Bj.bin"),fn("Vj.bin"),k,cols_j);

            {
                auto Vj = mmap_abrir_lectura(fn("Vj.bin"), k, cols_j);
                send_all(sock, Vj.data, Vj.bytes);
                mmap_cerrar(Vj);
                cout << "[worker->server] send V_j to server\n";
            }

            calcular_Ui_mmap(fn("Qfinal.bin"),fn("Utilde.bin"),fn("Ui.bin"),rows_i,k);

            send_all(sock, &h, sizeof(h));
            MMapMatrix Ui = mmap_abrir_lectura(fn("Ui.bin"), rows_i, k);
            send_all(sock, Ui.data, Ui.bytes);
            mmap_cerrar(Ui);

            std::cout << "[worker -> server] sent U_i\n";                

        } else {
            cerr << "[worker] unknown header id: " << h.id << "\n";
            break;
        }
    }
    close(sock);
}

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "Uso: " << argv[0] << " <IP_del_servidor> [puerto]\n";
        return 1;
    }

    const char* ip = argv[1];
    int port = (argc >= 3 ? atoi(argv[2]) : WORKER_PORT);

    int s = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in serv{};
    serv.sin_family = AF_INET;
    serv.sin_port = htons(port);
    inet_pton(AF_INET, ip, &serv.sin_addr);

    if (connect(s, (sockaddr*)&serv, sizeof(serv)) < 0) {
        perror("connect worker");
        return 1;
    }

    cout << "[worker] connected to server on port " << port << "\n";
    manejar_conexion(s);
    return 0;
}
