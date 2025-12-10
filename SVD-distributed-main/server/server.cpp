// server.cpp
// Compile: g++ -std=c++17 server.cpp -o server -pthread

#include "../protocolo.hpp"
#include "../matrix_io.hpp"
#include "../matrix_math.hpp"

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include <map>
#include <algorithm>
#include <cstring>
#include <iomanip>


using namespace std;

vector<int> worker_sockets;
mutex worker_mtx;

// Imprime en forma de tabla los campos de un protocolo para facilitar el debug.
void print_protocol_table(const string& protocol_name, const vector<vector<string>>& rows) {
    cout << "\n=================================================================================\n";
    cout << " PROTOCOLO: " << protocol_name << "\n";
    cout << "=================================================================================\n";
    cout << "| " << left << setw(12) << "Campo"
         << "| " << setw(8)  << "Tamano"
         << "| " << setw(20) << "Descripcion"
         << "| " << "Valor Real" << "\n";
    cout << "|-------------|--------|--------------------|------------------------------------\n";
    for (const auto& row : rows) {
        cout << "| " << left << setw(12) << row[0]
             << "| " << setw(8)  << row[1]
             << "| " << setw(20) << row[2]
             << "| " << row[3] << "\n";
    }
    cout << "=================================================================================\n\n";
}

void worker_acceptor_thread() {
    int ls = socket(AF_INET, SOCK_STREAM, 0);
    if (ls < 0) { perror("socket"); return; }
    int opt = 1; 
    setsockopt(ls, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    sockaddr_in addr{}; 
    addr.sin_family = AF_INET; 
    addr.sin_port = htons(WORKER_PORT); 
    addr.sin_addr.s_addr = INADDR_ANY;
    if (bind(ls, (sockaddr*)&addr, sizeof(addr)) < 0) { perror("bind"); close(ls); return; }
    if (listen(ls, 64) < 0) { perror("listen"); close(ls); return; }
    cerr << "[server] el worker_acceptor está escuchando en el puerto " << WORKER_PORT << "\n";
    while (true) {
        int ws = accept(ls, NULL, NULL);
        if (ws < 0) { perror("accept"); continue; }
        {
            lock_guard<mutex> lk(worker_mtx);
            worker_sockets.push_back(ws);
        }
        cerr << "[server] trabajador conectado fd=" << ws << " total=" << worker_sockets.size() << "\n";
    }
}

void send_Ai_to_worker(int ws, uint64_t rows_i, uint64_t n, const float* Ai) {
    vector<vector<string>> table = {
        {"B",         "1 B",  "Tipo de mensaje", "ID_A (1)"},
        {"matrix_id", "4 B",  "Id global",       "1"},
        {"block_id",  "4 B",  "Id bloque",       to_string(ws)},
        {"rows",      "2 B",  "Filas",           to_string(rows_i)},
        {"cols",      "2 B",  "Columnas",        to_string(n)},
        {"data",      "Var",  "Valores bloque",  "[" + to_string(rows_i * n) + " floats]"}
    };
    print_protocol_table("ENVIO_BLOQUE(B)", table);

    MsgHeader h; h.id = ID_A; h.a = rows_i; h.b = n;
    send_all(ws, &h, sizeof(h));
    size_t bytes = (size_t)rows_i * n * sizeof(float);
    send_all(ws, Ai, bytes);
}

void send_seed_to_worker(int ws, uint64_t seed, uint64_t k) {
    vector<vector<string>> table = {
        {"S",    "1 B", "Tipo de mensaje", "ID_S (2)"},
        {"seed", "8 B", "Semilla",         to_string(seed)}
    };
    print_protocol_table("ENVIOS(S)", table);

    MsgHeader h(ID_S,seed,k); //h.id = ID_S; h.a = k; h.b = 0;
    send_all(ws, &h, sizeof(h));
}

void recv_R_from_worker_mmap(int ws, const string &fileRi, uint64_t k){
    MsgHeader h;
    recv_all(ws, &h, sizeof(h));

    vector<vector<string>> table = {
        {"C",         "1 B", "Tipo",             "ID_UI (3)"},
        {"matrix_id", "4 B", "Id matriz",        "1"},
        {"rows",      "2 B", "Filas (k)",        to_string(h.a)},
        {"cols",      "2 B", "Columnas (k)",     to_string(h.b)},
        {"data",      "Var", "Datos de Y (R)",   "[Matriz reducida]"}
    };
    print_protocol_table("ENVIA_RESULTADOSR_i(C)", table);

    auto Ri = mmap_create(fileRi, k, k);

    size_t bytes = (size_t)k * k * sizeof(float);
    recv_all(ws, Ri.data, bytes);

    mmap_close(Ri);
}

void build_Rstack_mmap(const vector<string>& fileRis, int W, int k, const string& fileRstack){

    auto Rstack = mmap_create(fileRstack, W*k, k);

    for (int w = 0; w < W; ++w) {
        auto Ri = mmap_open_read(fileRis[w], k, k);
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                Rstack.data[(uint64_t)(w*k + i)*k + j] = Ri.data[(uint64_t)i*k + j];
        mmap_close(Ri);
    }

    mmap_close(Rstack);
}

void extract_Qr_blocks(const string fileQ, int W, int k, vector<string>& fileQri_list){
    auto Qbig = mmap_open_read(fileQ, W*k, k);

    for (int w = 0; w < W; ++w) {
        fileQri_list[w] = "Qr_" + to_string(w) + ".bin";

        auto Qr = mmap_create(fileQri_list[w], k, k);

        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j)
                Qr.data[(uint64_t)i*k + j] = Qbig.data[(uint64_t)(w*k + i)*k + j];

        mmap_close(Qr);
    }

    mmap_close(Qbig);
}

void send_Qr_to_worker(int sock, const string& fileQr, int k)
{
    MsgHeader h(ID_Q, k, k);
    send_all(sock, &h, sizeof(h));

    auto Qr = mmap_open_read(fileQr, k, k);
    send_all(sock, Qr.data, sizeof(float)*k*k);
    mmap_close(Qr);
}

void recv_Bi_from_worker_mmap(int ws, const string& fileBi, int k, int n) {
    MsgHeader h;
    recv_all(ws, &h, sizeof(h));

    vector<vector<string>> table = {
        {"D",         "1 B", "Tipo",             "ID_UI (3)"},
        {"matrix_id", "4 B", "Id matriz",        "1"},
        {"rows",      "2 B", "Filas (k)",        to_string(h.a)},
        {"cols",      "2 B", "Columnas (n)",     to_string(h.b)},
        {"data",      "Var", "Datos Z (Bi)",     "[Bloque Bi]"}
    };
    print_protocol_table("ENVIA_RESULTADOSB_i(D)", table);

    auto Bi = mmap_create(fileBi, k, n);
    recv_all(ws, Bi.data, sizeof(float)*k*n);
    mmap_close(Bi);
}

void assemble_B_mmap(const vector<string>& fileBi_list, const vector<int>& nrows,
                     int W, int k, int n, const string& outFile){
    auto Bfinal = mmap_create(outFile, k, n);
    memset(Bfinal.data, 0, sizeof(float) * (size_t)k * n);

    for (int i = 0; i < W; ++i) {
        if (nrows[i] == 0) continue;
        auto Bi = mmap_open_read(fileBi_list[i], k, n);

        size_t elems = (size_t)k * n;
        for (size_t p = 0; p < elems; ++p)
            Bfinal.data[p] += Bi.data[p];

        mmap_close(Bi);
    }

    mmap_close(Bfinal);
}

void send_B_block_to_worker_mmap(int sock, MMapMatrix& B, int k, int n,
                                 int col_start, int col_count){
    // Enviamos cabecera para que el worker conozca cuántas columnas recibe
    MsgHeader h(ID_D, k, col_count);
    send_all(sock, &h, sizeof(h));
    for (int r = 0; r < k; ++r) {
        float* row_ptr = &B.data[(size_t)r * n + col_start];
        size_t bytes = sizeof(float) * col_count;
        send_all(sock, row_ptr, bytes);
    }
}
// Función NUEVA para recibir C_j correctamente
void recv_Cj_from_worker_mmap(int ws, const string& fileCj, int k) {
    MsgHeader h;
    recv_all(ws, &h, sizeof(h));


    vector<vector<string>> table = {
        {"C_part",    "1 B", "Tipo",             "ID_UI (Reusado)"}, // O el ID que estés usando
        {"matrix_id", "4 B", "Id matriz",        "1"},
        {"rows",      "2 B", "Filas (k)",        to_string(h.a)},
        {"cols",      "2 B", "Columnas (k)",     to_string(h.b)},
        {"data",      "Var", "Datos C_j",        "[Matriz k x k]"}
    };
    print_protocol_table("RECIBIR_C_PARCIAL(E)", table);
    // ------------------------------

    auto Cj = mmap_create(fileCj, k, k);
    // Nota: Cj es k x k
    size_t bytes = (size_t)k * k * sizeof(float);
    recv_all(ws, Cj.data, bytes);
    mmap_close(Cj);
}

void eigendecompose_C_mmap(const string &Cfile, int k, const string &UtilFile,
    const string &LambdaFile){

    MMapMatrix Cmm = mmap_open_read(Cfile, k, k);
    MMapMatrix Utilmm = mmap_create(UtilFile,k, k);
    MMapMatrix Lambdamm = mmap_create(LambdaFile, k, 1);

    vector<float> A((size_t)k * k);
    for (int i = 0; i < k * k; ++i) A[i] = Cmm.data[i]; //necesario pal jacobo

    jacobi_eigen_inplace(A.data(), k, Utilmm.data, Lambdamm.data);
    vector<pair<float,int>>order(k);
    for (int i = 0; i < k; ++i) order[i] = {Lambdamm.data[i], i};
    sort(order.begin(), order.end(),
              [](const auto &a, const auto &b){ return a.first > b.first; });

    vector<float> Util_sorted((size_t)k * k);
    vector<float> Lambda_sorted(k);
    for (int newc = 0; newc < k; ++newc) {
        int oldc = order[newc].second;
        Lambda_sorted[newc] = order[newc].first;
        for (int r = 0; r < k; ++r)
            Util_sorted[(size_t)r * k + newc] = Utilmm.data[(size_t)r * k + oldc];
    }
    memcpy(Lambdamm.data, Lambda_sorted.data(), (size_t)k * sizeof(float));
    memcpy(Utilmm.data, Util_sorted.data(), (size_t)k * k * sizeof(float));

    msync(Utilmm.data, Utilmm.bytes, MS_SYNC);
    msync(Lambdamm.data, Lambdamm.bytes, MS_SYNC);

    mmap_close(Cmm);
    mmap_close(Utilmm);
    mmap_close(Lambdamm);
}

void sigma_and_inv_mmap(const string &LambdaFile,
    const string &SigmaFile,const string &SigmaInvFile,int k){

    auto Lambdamm = mmap_open_read(LambdaFile, k, 1);
    auto Sigmamm = mmap_create(SigmaFile, k, 1);
    auto SigmaInvmm = mmap_create(SigmaInvFile, k, 1);

    for (int i = 0; i < k; ++i) {
        float lambda = Lambdamm.data[i];
        float sigma = (lambda > 0.0f) ? sqrtf(max(lambda, 0.0f)) : 0.0f;
        Sigmamm.data[i] = sigma;
        SigmaInvmm.data[i] = (sigma > 1e-12f) ? (1.0f / sigma) : 0.0f;
    }

    msync(Sigmamm.data, Sigmamm.bytes, MS_SYNC);
    msync(SigmaInvmm.data, SigmaInvmm.bytes, MS_SYNC);

    mmap_close(Lambdamm);
    mmap_close(Sigmamm);
    mmap_close(SigmaInvmm);
}

void send_Sigma_mmap(int sock, const string& sigmaFile, int k){
    MMapMatrix Sm = mmap_open_read(sigmaFile, k, 1);

    MsgHeader hs(ID_S,k,0);
    send_all(sock, &hs, sizeof(hs));
    send_all(sock, Sm.data, (size_t)k * sizeof(float));
    mmap_close(Sm);

    cout << "[server->client] se envió Sigma mediante mmap (" << sigmaFile << ")\n";
}

void send_Util_Sinv_to_worker_mmap(int sock,const string &UtilFile,
    const string &SigmaInvFile, int k){

    MMapMatrix Utilmm = mmap_open_read(UtilFile, k, k);
    MMapMatrix SigmaInvmm = mmap_open_read(SigmaInvFile, k, 1);

    send_all(sock, Utilmm.data, Utilmm.bytes);
    send_all(sock, SigmaInvmm.data, SigmaInvmm.bytes);

    mmap_close(Utilmm);
    mmap_close(SigmaInvmm);
}

void recv_Vj_from_worker_mmap(int sock, const string& fileVj, int k, int cols_j){
    vector<vector<string>> table = {
        {"V_part",    "1 B", "Tipo",             "ID_VT (Gen)"},
        {"rows",      "2 B", "Filas (k)",        to_string(k)},
        {"cols",      "2 B", "Columnas (sub)",   to_string(cols_j)},
        {"data",      "Var", "Datos V_j",        "[Bloque V_j]"}
    };
    print_protocol_table("RECIBIR_V_PARCIAL", table);

    MMapMatrix Vj = mmap_create(fileVj, k, cols_j);
    size_t bytes = (size_t)k * cols_j * sizeof(float);
    recv_all(sock, Vj.data, bytes);
    mmap_close(Vj);
}

void assemble_Vt_mmap(const vector<string>& fileVj_list, const vector<int>& start_cols,
    const vector<int>& cols_list, int W, int k, int n, const string& outFile){

    MMapMatrix Vt = mmap_create(outFile, k, n);
    memset(Vt.data, 0, (size_t)k * n * sizeof(float));

    for (int i = 0; i < W; ++i) {
        int cols_j = cols_list[i];
        if (cols_j == 0) continue;

        int start = start_cols[i];
        MMapMatrix Vj = mmap_open_read(fileVj_list[i], k, cols_j);

        for (int r = 0; r < k; ++r) {
            float* dst = &Vt.data[(size_t)r * n + start];
            float* src = &Vj.data[(size_t)r * cols_j];

            memcpy(dst, src, sizeof(float) * (size_t)cols_j);
        }

        mmap_close(Vj);
    }

    mmap_close(Vt);
}

void send_Vt_to_client(int cs, const string &VtFile, int k, int n){
    MsgHeader h(ID_VT,k,n);
    send_all(cs, &h, sizeof(h));

    MMapMatrix Vt = mmap_open_read(VtFile, k, n);
    send_all(cs, Vt.data, Vt.bytes);
    mmap_close(Vt);

    cout << "[server->client] archivo enviado (" <<VtFile<< ")\n";
}

void recv_Ui_from_worker_mmap(int sock, const string &outfile, int rows_i, int k){
    vector<vector<string>> table = {
        {"U_part",    "1 B", "Tipo",             "ID_UI (Gen)"},
        {"rows",      "2 B", "Filas (sub)",      to_string(rows_i)},
        {"cols",      "2 B", "Columnas (k)",     to_string(k)},
        {"data",      "Var", "Datos U_i",        "[Bloque U_i]"}
    };
    print_protocol_table("RECIBIR_U_PARCIAL", table);

    MMapMatrix Uimm = mmap_create(outfile, rows_i, k);
    recv_all(sock, Uimm.data, Uimm.bytes);
    mmap_close(Uimm);

    cout << "[server] se recibió Ui -> " << outfile 
              << " (" << rows_i << " x " << k << ")\n";
}

void assemble_U_mmap(const vector<string> &Ui_files, 
    const vector<int> &nrows,int W, int k, int m, const string &Ufile){

    MMapMatrix Ufinal = mmap_create(Ufile, m, k);
    int row_offset = 0;

    for (int i = 0; i < W; ++i){
        if (nrows[i] == 0) continue;

        int ri = nrows[i];
        MMapMatrix Ui = mmap_open_read(Ui_files[i], ri, k);

        for (int r = 0; r < ri; ++r){
            memcpy(&Ufinal.data[(size_t)(row_offset + r) * k],
                   &Ui.data[(size_t)r * k],
                   sizeof(float) * k);
        }

        row_offset += ri;
        mmap_close(Ui);
    }

    mmap_close(Ufinal);
    cout << "[server] ensambló U con mmap\n";
}

void clientHandler(int cs){
    cout << "[server] cliente conectado fd=" << cs << "\n";
    MsgHeader h;
    if (!recv_all(cs, &h, sizeof(h))) { close(cs); return; }

    // Optional handshake: jefe asks for available workers before sending matrix
    if (h.id == ID_H) {
        uint64_t available = 0;
        {
            lock_guard<mutex> lk(worker_mtx);
            available = worker_sockets.size();
        }
        MsgHeader resp(ID_H, available, 0);
        send_all(cs, &resp, sizeof(resp));
        if (!recv_all(cs, &h, sizeof(h))) { close(cs); return; }
    }

    // Expect header ID_A with a = n, b = k
    if (h.id != ID_A) { cerr << "se esperaba un encabezado A\n"; close(cs); return; }
    uint64_t n = h.a;
    uint64_t k_full = h.b; // k + p (oversampling)
    cerr << "[server] randomized SVD con n = "<<n<<" k_full="<<k_full<<"\n";

    // prepare file server_matrix.bin
    string fname = "server_matrix.bin";
    uint64_t elems = n * n;
    size_t bytes = (size_t)elems * sizeof(float);
    int fd = open(fname.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    ftruncate(fd, bytes);
    void* mapv = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    float* M = (float*)mapv;

    size_t rec = 0;
    char* base = (char*)mapv;
    while (rec < bytes) {
        size_t chunk = (bytes - rec) > CHUNK ? CHUNK : (bytes - rec);
        if (!recv_all(cs, base + rec, chunk)) { 
            cerr << "error al recibir la matriz\n"; 
            munmap(mapv, bytes); close(fd); close(cs); return; 
        }
        rec += chunk;
    }
    cerr << "[client->server] matriz A recibida ("<<n<<"x"<<n<<")\n";

    // receive target k (without oversampling) to truncate later
    uint64_t k_target = k_full;
    MsgHeader hk;
    if (recv_all(cs, &hk, sizeof(hk)) && hk.id == ID_K) {
        k_target = hk.a;
    } else {
        cerr << "[server] advertencia: no se encontró el encabezado ID_K, se usará k_full\n";
    }
    int k_req = (int)k_target;
    int k = (int)k_full;

    // Expect the boss to provide the random seed (ID_S). If missing, keep default.
    uint64_t seed = 1234567;
    MsgHeader seedMsg;
    if (recv_all(cs, &seedMsg, sizeof(seedMsg)) && seedMsg.id == ID_S) {
        seed = seedMsg.a;
        if (seedMsg.b && seedMsg.b != k)
            cerr << "[server] advertencia: el jefe indicó k="<<seedMsg.b<<" , pero el encabezado tenía k="<<k<<"\n";
    } else {
        cerr << "[server] advertencia: el jefe no envió un encabezado de semilla; se usará el valor por defecto "<<seed<<"\n";
    }

    vector<int> ws;
    {
        lock_guard<mutex> lk(worker_mtx);
        ws = worker_sockets;
    }
    int W = ws.size();
    vector<int> start(W), nrows(W);
    vector<int> start_cols(W), ncols(W);
    uint64_t base_rows = n / W;
    uint64_t extra = n % W;
    uint64_t cur = 0;
    for (int i = 0; i < W; ++i) {
        start[i] = cur;
        nrows[i] = base_rows + (i < (int)extra ? 1 : 0);
        cur += nrows[i];
    }
    // Partición por columnas (coherente con B y C)
    uint64_t base_cols = n / W;
    uint64_t extra_cols = n % W;
    uint64_t curc = 0;
    for (int i = 0; i < W; ++i) {
        start_cols[i] = curc;
        ncols[i] = base_cols + (i < (int)extra_cols ? 1 : 0);
        curc += ncols[i];
    }

    // 1) send Ai to each worker
    for (int i = 0; i < W; ++i) {
        uint64_t rowsi = nrows[i];
        if (rowsi == 0) continue;
        const float* Ai = M + (size_t)start[i] * n;
        send_Ai_to_worker(ws[i], rowsi, n, Ai);
    }
    cerr << "[server->workers] matriz Ai enviada a los trabajadores\n";

    // 2) send seed and k
    for (int i = 0; i < W; ++i) {
        send_seed_to_worker(ws[i], seed, k);
    }
    cerr << "[server->workers] se enviaron la semilla y el valor k\n";

    // 3) receive R_i from each worker
    vector<string> Ri_files(W);
    for (int i = 0; i < W; ++i) {
        if (nrows[i] == 0) continue;
        Ri_files[i] = "Ri_" + to_string(i) + ".bin";
        recv_R_from_worker_mmap(ws[i], Ri_files[i], k);
        cout<<"[workers->server] se recibió R_"<<i<<" del trabajador "<<i<<"\n";
    }

    
    // 4) TSQR mmap
    build_Rstack_mmap(Ri_files,W,k,"Rstack.bin");
    qr_mmap("Rstack.bin", W*k, k, "Qr.bin", "Rglobal.bin");
    cout << "[server] R_stack y TSQR completados\n";

    // 5) send Qr_i to each worker
    vector<string> fileQri(W);
    extract_Qr_blocks("Qr.bin", W, k, fileQri);
    for (int i = 0; i < W; ++i) {
        if (nrows[i] == 0) continue;
        send_Qr_to_worker(ws[i], fileQri[i], k);
        cout << "[server->workers] se envió Qglbal_"<<i<<" al trabajador " << i << "\n";
    }

    // 6) receive B_i from workers
    vector<string> fileBi_list(W);
    for (int i = 0; i < W; ++i) {
        if (nrows[i] == 0) continue;
        string fileBi = "B_i_" + to_string(i) + ".bin";
        fileBi_list[i] = fileBi;

        recv_Bi_from_worker_mmap(ws[i], fileBi, k, n);
        cout << "[workers->server] se recibió B_"<<i<<" del trabajador " << i << " (k x n)\n";
    }

    
    // 7) Construct B = [B1 B2 ...] horizontally (k x n)
    assemble_B_mmap(fileBi_list, nrows, W, k, n, "B_final.bin");
    cout << "[server] ensambló B (k x n) usando mmap\n";

    auto Bmap = mmap_open_read("B_final.bin", k, n);

    // 8) Enviar Bj a cada worker usando mmap (por columnas)
    for (int j = 0; j < W; ++j) {
        if (ncols[j] == 0) continue;
        send_B_block_to_worker_mmap(ws[j],Bmap,k, n, start_cols[j], ncols[j]);
        cout << "[server->workers] se envió Bj ("<< start_cols[j]<< " - " << start_cols[j]+ncols[j] << ")\n";
    }
    mmap_close(Bmap);


    // 9) receive C_j from workers and sum to C
    vector<string> fileCj_list(W);

    for (int j = 0; j < W; ++j) {
        if (ncols[j] == 0) continue;
         fileCj_list[j] = "Cj_" + to_string(j) + ".bin";

        recv_Cj_from_worker_mmap(ws[j], fileCj_list[j], k);
        cout << "[server] se recibió C_" << j << " del trabajador " << j << "\n";
    }
    //se usa la misma funcion que bmmap porque solo es para ensamblar
    assemble_B_mmap(fileCj_list, ncols, W, k, k,"C_final.bin");
    cout << "[server] matriz C (k x k) ensamblada usando mmap\n";
    
    eigendecompose_C_mmap("C_final.bin",k,"Utilde.bin","Lambda.bin");
    sigma_and_inv_mmap("Lambda.bin", "Sigma.bin","SigmaInv.bin", k);

    cout << "[server] se calcularon Sigma, Utilde y SigmaInv\n";

    for (int i = 0; i < W; ++i) {
        send_Util_Sinv_to_worker_mmap(ws[i],"Utilde.bin", "SigmaInv.bin",k);
        cout << "[server->Workers] se enviaron Utilde y SigmaInv al trabajador " << i << "\n";
    }

    vector<string> fileVj_list(W);

    for (int i = 0; i < W; ++i) {
        if (ncols[i] == 0) continue;
        fileVj_list[i] = "V_j_" + to_string(i) + ".bin";

        recv_Vj_from_worker_mmap(ws[i], fileVj_list[i], k, ncols[i]);

        cout << "[workers->server] se recibió V_" << i 
            << " (" << k << " x " << ncols[i] << ")\n";
    }

    assemble_Vt_mmap(fileVj_list, start_cols, ncols, W, k, n, "Vt_final.bin");
    cout << "[server] ensambló Vt (k x n) usando mmap\n";

    vector<string> Ui_files(W);
    for (int i = 0; i < W; i++){
        if (nrows[i] == 0) continue;
        Ui_files[i] = "U_i_" + to_string(i) + ".bin";

        recv_Ui_from_worker_mmap(ws[i], Ui_files[i], nrows[i], k);
    }


    assemble_U_mmap(Ui_files, nrows, W, k, n, "U_final.bin");

    // Truncate to requested k (top-k) if oversampling was used
    string U_to_send = "U_final.bin";
    string Vt_to_send = "Vt_final.bin";
    string Sigma_to_send = "Sigma.bin";
    if (k_req < k) {
        U_to_send = "U_final_k.bin";
        Vt_to_send = "Vt_final_k.bin";
        Sigma_to_send = "Sigma_k.bin";

        // U: n x k_req
        {
            MMapMatrix Ufull = mmap_open_read("U_final.bin", n, k);
            MMapMatrix Utrim = mmap_create(U_to_send, n, k_req);
            for (uint64_t r = 0; r < n; ++r) {
                memcpy(&Utrim.data[(size_t)r * k_req],
                       &Ufull.data[(size_t)r * k],
                       sizeof(float) * (size_t)k_req);
            }
            mmap_close(Ufull);
            mmap_close(Utrim);
        }

        // Vt: k_req x n (take first k_req rows)
        {
            MMapMatrix Vfull = mmap_open_read("Vt_final.bin", k, n);
            MMapMatrix Vtrim = mmap_create(Vt_to_send, k_req, n);
            for (int r = 0; r < k_req; ++r) {
                memcpy(&Vtrim.data[(size_t)r * n],
                       &Vfull.data[(size_t)r * n],
                       sizeof(float) * (size_t)n);
            }
            mmap_close(Vfull);
            mmap_close(Vtrim);
        }

        // Sigma: first k_req
        {
            MMapMatrix Sfull = mmap_open_read("Sigma.bin", k, 1);
            MMapMatrix Strim = mmap_create(Sigma_to_send, k_req, 1);
            memcpy(Strim.data, Sfull.data, sizeof(float) * (size_t)k_req);
            mmap_close(Sfull);
            mmap_close(Strim);
        }

        cout << "[server] se truncaron los resultados a k_req=" << k_req << "\n";
    }

    vector<vector<string>> tableF = {
        {"F",    "1 B", "Tipo",            "ID_DONE (99)"},
        {"U",    "Var", "Matriz U",        "Enviada (" + to_string(n) + "x" + to_string(k_req) + ")"},
        {"S",    "Var", "Vector S",        "Enviado (" + to_string(k_req) + ")"},
        {"V",    "Var", "Matriz V",        "Enviada (" + to_string(n) + "x" + to_string(k_req) + ")"}
    };
    print_protocol_table("Enviar_resultadoFinal(F)", tableF);

    // Enviar U (header + stream fichero U.bin)
    {
        MsgHeader h(ID_UT,n,k_req);
        send_all(cs,&h,sizeof(h));
        MMapMatrix Ufile = mmap_open_read(U_to_send, n, k_req);
        send_all(cs,Ufile.data,Ufile.bytes);
        mmap_close(Ufile);
    }

    // send S
    send_Sigma_mmap(cs, Sigma_to_send, k_req);

    // send V^T
    send_Vt_to_client(cs,Vt_to_send,k_req,n);

    vector<vector<string>> tableH = {
        {"H",    "1 B", "Tipo de mensaje", "ID_DONE (99)"},
        {"Info", "-",   "Estado",          "Operacion Completada"}
    };
    print_protocol_table("Enviar_sol_cliente(H)", tableH);

    MsgHeader md; md.id = ID_DONE; md.a = 0; md.b = 0;
    send_all(cs, &md, sizeof(md));
    cerr << "[server] se enviaron U, S, V^T y DONE al cliente\n";

    munmap(mapv, bytes);
    close(fd);
    close(cs);   
}

int main() {
    thread(worker_acceptor_thread).detach();

    int ls = socket(AF_INET, SOCK_STREAM, 0);
    if (ls < 0) { perror("Conexion"); return 1; }
    int opt = 1; 
    setsockopt(ls, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    sockaddr_in addr{}; 
    addr.sin_family = AF_INET; 
    addr.sin_port = htons(BOSS_SERVER_PORT); 
    addr.sin_addr.s_addr = INADDR_ANY;
    if (bind(ls, (sockaddr*)&addr, sizeof(addr)) < 0) { 
        perror("bind"); return 1; 
    }
    if (listen(ls, 10) < 0) { 
        perror("Escuchar"); return 1; 
    }
    cerr << "[server] esperando conexión del jefe en el puerto " << BOSS_SERVER_PORT << "\n";

    while (true) {
        int clientSocket = accept(ls, NULL, NULL);
        if (clientSocket < 0) { perror("Aceptar"); continue; }

        thread(clientHandler, clientSocket).detach();
    }

    return 0;
}