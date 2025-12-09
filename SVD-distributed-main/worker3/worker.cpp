#include "../common.hpp"
#include "../matrix_io.hpp"
#include "../matrix_math.hpp"


#include <arpa/inet.h>
#include <sys/socket.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include <thread>
#include <fcntl.h>
#include <sys/mman.h>

using namespace std;

void handle_connection(int sock) {
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

            // receive matrix A_i
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

            // Seed proccess
            MsgHeader hs;
            recv_all(sock, &hs, sizeof(hs));
            uint64_t seed = hs.a;
            uint64_t k = hs.b;
            cout << "[server->worker] seed="<<seed<<" k="<<k<<"\n";

            generate_omega_mmap(seed, n, k, fn("Omega.bin"));
            matmul_A_Omega_mmap(fn("Ai.bin"), rows_i, n, fn("Omega.bin"), k, fn("Y.bin"));

            qr_mmap(fn("Y.bin"), rows_i, k, fn("Qi.bin"), fn("Ri.bin"));

            {
                auto Rm = mmap_open_read(fn("Ri.bin"), k, k);
                MsgHeader hr(ID_R, k, k);
                send_all(sock, &hr, sizeof(hr));
                send_all(sock, Rm.data, (size_t)k * k * sizeof(float));
                mmap_close(Rm);
                cout << "[worker->server] sent R_i to server\n";
            }

            {
                MsgHeader hq;
                recv_all(sock, &hq, sizeof(hq));
                auto Qr = mmap_create(fn("Qr.bin"), hq.a, hq.b);
                recv_all(sock, Qr.data, sizeof(float)*hq.a*hq.b);
                mmap_close(Qr);
                cout << "[server->worker] received Qr\n";
            }

            {
                auto Qi  = mmap_open_read(fn("Qi.bin"), rows_i, k);
                auto Qr  = mmap_open_read(fn("Qr.bin"), k, k);
                auto Qf  = mmap_create(fn("Qfinal.bin"), rows_i, k);

                for (uint64_t r = 0; r < rows_i; ++r)
                    for (uint64_t c = 0; c < k; ++c) {
                        float sum = 0;
                        for (uint64_t kk = 0; kk < (uint64_t)k; ++kk)
                            sum += Qi.data[r*(uint64_t)k + kk] * Qr.data[kk*(uint64_t)k + c];

                        Qf.data[r*(uint64_t)k + c] = sum;
                    }
                mmap_close(Qi);
                mmap_close(Qr);
                mmap_close(Qf);
                
                cout << "[worker] Calculate Q_final\n";
            }

            {
                auto Qf = mmap_open_read(fn("Qfinal.bin"), rows_i, k);
                auto Ai = mmap_open_read(fn("Ai.bin"), rows_i, n);
                auto Bi = mmap_create(fn("Bi.bin"), k, n);

                for (uint64_t i = 0; i < (uint64_t)k; ++i)
                    for (uint64_t j = 0; j < (uint64_t)n; ++j) {
                        float sum = 0.0f;
                        for (uint64_t t = 0; t < (uint64_t)rows_i; ++t)
                            sum += Qf.data[t*(uint64_t)k + i] * Ai.data[t*(uint64_t)n + j];
                        Bi.data[i*(uint64_t)n + j] = sum;
                    }
                mmap_close(Qf);
                mmap_close(Ai);
                mmap_close(Bi);

                cout << "[worker] Calculate B_i = Q_final x A_i\n";
            }

            {
                auto Bi = mmap_open_read(fn("Bi.bin"), k, n);
                MsgHeader hb;
                hb.id = ID_B;
                hb.a = k;
                hb.b = n;

                send_all(sock, &hb, sizeof(hb));
                send_all(sock, Bi.data, sizeof(float)*k*n);
                mmap_close(Bi);

                cout << "[worker->server] sent B_i to server\n";
            }
            
            int cols_j = 0;
            {
                MsgHeader hd;
                recv_all(sock, &hd, sizeof(hd));
                cols_j = static_cast<int>(hd.b);
                auto Bj = mmap_create(fn("Bj.bin"), k, cols_j);
                recv_all(sock, Bj.data, sizeof(float) * k * cols_j);
                mmap_close(Bj);
                cout << "[server->worker] received B_j from server (" << cols_j << " cols)\n";                
            }

            {
                auto Bj = mmap_open_read(fn("Bj.bin"), k, cols_j);
                auto Cj = mmap_create(fn("Cj.bin"), k, k);
                // BB^T
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

                mmap_close(Bj);
                mmap_close(Cj);

                cout << "[worker] Calculate C_j (k x k)\n";  
            }

            {
                MsgHeader h(ID_C, k, k);
                send_all(sock, &h, sizeof(h));

                auto Cj = mmap_open_read(fn("Cj.bin"), k, k);
                send_all(sock, Cj.data, sizeof(float) * k * k);
                mmap_close(Cj);

                cout << "[worker->server] sent C_j to server\n";
            }
            
            {
                auto Utilmm = mmap_create(fn("Utilde.bin"), k, k);
                recv_all(sock, Utilmm.data, sizeof(float) * k * k);

                auto SigmaInvmm = mmap_create(fn("SigmaInv.bin"), k, 1);
                recv_all(sock, SigmaInvmm.data, sizeof(float) * k);

                mmap_close(Utilmm);
                mmap_close(SigmaInvmm);
                cout << "[server->worker] received Utilde and SigmaInv from server\n";           
            }

            compute_Vj_mmap(fn("Utilde.bin"),fn("SigmaInv.bin"),fn("Bj.bin"),fn("Vj.bin"),k,cols_j);

            {
                auto Vj = mmap_open_read(fn("Vj.bin"), k, cols_j);
                send_all(sock, Vj.data, Vj.bytes);
                mmap_close(Vj);
                cout << "[worker->server] send V_j to server\n";
            }

            compute_Ui_mmap(fn("Qfinal.bin"),fn("Utilde.bin"),fn("Ui.bin"),rows_i,k);

            send_all(sock, &h, sizeof(h));
            MMapMatrix Ui = mmap_open_read(fn("Ui.bin"), rows_i, k);
            send_all(sock, Ui.data, Ui.bytes);
            mmap_close(Ui);

            std::cout << "[worker -> server] sent U_i\n";                

            /*
            while (true) {

                MsgHeader mh;
                if (!recv_all(sock, &mh, sizeof(mh))) { cerr << "[worker] server closed mid-flow\n"; return; }
                if (mh.id == ID_D) {
                    // server sends (k x colsblock) as header a=k b=colsblock
                    int krecv = (int)mh.a;
                    int colsblock = (int)mh.b;
                    vector<float> Bj((size_t)krecv * colsblock);
                    recv_all(sock, Bj.data(), Bj.size()*sizeof(float));
                    // compute Cj = Bj * Bj^T (k x k)
                    vector<float> Cj = compute_Cj(Bj, krecv, colsblock);
                    // send Cj back with ID_C
                    MsgHeader mc; mc.id = ID_C; mc.a = krecv; mc.b = krecv;
                    send_all(sock, &mc, sizeof(mc));
                    send_all(sock, Cj.data(), (size_t)krecv*krecv*sizeof(float));
                    cout << "[worker] processed D -> sent Cj\n";
                } else if (mh.id == ID_UT) {
                    // server sends Utilde (k x k) and Sigma^{-1 (k)}
                    uint64_t krecv = mh.a;
                    uint64_t dummy = mh.b; // maybe zero
                    vector<float> Util((size_t)krecv * krecv);
                    recv_all(sock, Util.data(), Util.size()*sizeof(float));
                    vector<float> SigmaInv((size_t)krecv);
                    recv_all(sock, SigmaInv.data(), SigmaInv.size()*sizeof(float));
                    // Compute Vj^T = Sigma^{-1} * Util^T * B_i
                    // First compute M = Util^T * B_i  => (k x n)
                    vector<float> M((size_t)krecv * Bi.size() / (size_t)krecv); // k x n
                    int ncols = Bi.size() / (size_t)krecv;
                    for (int i = 0; i < krecv; ++i)
                        for (int j = 0; j < ncols; ++j) {
                            float s = 0;
                            for (int t = 0; t < krecv; ++t)
                                s += Util[(size_t)t * krecv + i] * Bi[(size_t)t * ncols + j];
                            M[(size_t)i * ncols + j] = s;
                        }
                    // multiply by Sigma^{-1}: each row i scaled
                    for (int i = 0; i < krecv; ++i)
                        for (int j = 0; j < ncols; ++j)
                            M[(size_t)i * ncols + j] *= SigmaInv[i];

                    // send V_j^T with header ID_V: a=k b=ncols
                    MsgHeader mv; mv.id = ID_V; mv.a = krecv; mv.b = ncols;
                    send_all(sock, &mv, sizeof(mv));
                    send_all(sock, M.data(), (size_t)krecv * ncols * sizeof(float));
                    cout << "[worker] sent V_j^T to server\n";
                } else if (mh.id == ID_UI) {
                    // server requests U_i: sends Utilde (k x k)
                    uint64_t krecv = mh.a;
                    vector<float> Util((size_t)krecv * krecv);
                    recv_all(sock, Util.data(), Util.size()*sizeof(float));
                    // compute U_i = Qf * Util   => rows_i x k
                    vector<float> Ui((size_t)rows_i * krecv, 0.0f);
                    for (size_t r = 0; r < (size_t)rows_i; ++r)
                        for (size_t c = 0; c < (size_t)krecv; ++c)
                            for (size_t t = 0; t < (size_t)krecv; ++t)
                                Ui[r*krecv + c] += Qf[r*krecv + t] * Util[t*krecv + c];
                    // send Ui back
                    MsgHeader mu; mu.id = ID_UI; mu.a = rows_i; mu.b = krecv;
                    send_all(sock, &mu, sizeof(mu));
                    send_all(sock, Ui.data(), (size_t)rows_i * krecv * sizeof(float));
                    cout << "[worker] sent U_i to server\n";
                } else if (mh.id == ID_DONE) {
                    cout << "[worker] received DONE from server stage\n";
                    break; // exit inner loop and wait for next client block or close
                } else {
                    cerr << "[worker] unknown message id: " << mh.id << "\n";
                    break;
                }
            
            } */
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
    handle_connection(s);
    return 0;
}
