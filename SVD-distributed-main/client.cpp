
#include "protocolo.hpp"
#include "mapeo_matriz.hpp"
#include <arpa/inet.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <chrono>
#include <filesystem>

using namespace std;


string MATRIX_FILE = "matrix.bin";
uint64_t N_global = 0;

bool cargar_matriz_desde_bin(const string& ruta) {
    namespace fs = std::filesystem;
    if (!fs::exists(ruta)) {
        cerr << "No existe el archivo: " << ruta << "\n";
        return false;
    }
    uint64_t bytes = fs::file_size(ruta);
    if (bytes % sizeof(float) != 0) {
        cerr << "El archivo no es múltiplo de sizeof(float)\n";
        return false;
    }
    uint64_t elems = bytes / sizeof(float);
    uint64_t n = static_cast<uint64_t>(sqrt(static_cast<long double>(elems)));
    if (n * n != elems) {
        cerr << "El archivo no representa una matriz cuadrada (elementos="<<elems<<")\n";
        return false;
    }

    if (ruta == MATRIX_FILE) {
        N_global = n;
        cout << "Matriz ya presente en " << MATRIX_FILE << " ("<< n <<" x "<< n <<")\n";
        return true;
    }

    auto src = mmap_abrir_lectura(ruta, n, n);
    auto dst = mmap_crear(MATRIX_FILE, n, n);
    memcpy(dst.data, src.data, bytes);
    msync(dst.data, dst.bytes, MS_SYNC);
    mmap_cerrar(src);
    mmap_cerrar(dst);

    N_global = n;
    cout << "Matriz cargada desde " << ruta << " ("<< n <<" x "<< n <<") a " << MATRIX_FILE << "\n";
    return true;
}

void generar_matriz(uint64_t N) {
    auto mm = mmap_crear(MATRIX_FILE, N, N);
    uint32_t seed = static_cast<uint32_t>(
        std::chrono::steady_clock::now().time_since_epoch().count() ^
        std::random_device{}());
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    uint64_t elems = N * (uint64_t)N;
    for (uint64_t i = 0; i < elems; ++i) mm.data[i] = dist(rng);
    msync(mm.data, mm.bytes, MS_SYNC);
    mmap_cerrar(mm);
    N_global = N;
    cout << "Matriz Generada " << N << "x" << N << " en archivo " << MATRIX_FILE << "\n";
}

void mostrar_matriz(const string &filename,uint64_t max_rows = 5,uint64_t max_cols = 5){
    namespace fs = std::filesystem;
    if (!fs::exists(filename)) {
        cout << filename << ": no existe\n\n";
        return;
    }
    if (N_global == 0) {
        cout << filename << ": No data (N_global==0)\n\n";
        return;
    }

    cout << "========= " << filename << " (max "
         << max_rows << "x" << max_cols << ") =========\n";

    try {
        auto mm = mmap_abrir_lectura(filename, N_global, N_global);
        for (uint64_t r = 0; r < min(max_rows, N_global); ++r) {
            for (uint64_t c = 0; c < min(max_cols, N_global); ++c) {
                cout << mm.data[r * N_global + c] << "\t";
            }
            cout << "\n";
        }
        cout << "\n";
        mmap_cerrar(mm);
    } catch (const std::exception& e) {
        cout << filename << ": error al leer (" << e.what() << ")\n\n";
    }
}

void mostrar_sigma(const string &filename, uint64_t max_vals = 5){
    namespace fs = std::filesystem;
    if (!fs::exists(filename)) {
        cout << filename << ": no existe\n\n";
        return;
    }
    if (N_global == 0) {
        cout << "Sigma: No data (N_global==0)\n\n";
        return;
    }

    cout << "========= SIGMA (S) (max " << max_vals << ") =========\n";
    try {
        auto mm = mmap_abrir_lectura(filename, N_global, 1);  // vector de tamaño N
        for (uint64_t i = 0; i < min(max_vals, N_global); ++i) {
            cout << mm.data[i] << "\t";
        }
        cout << "\n\n";
        mmap_cerrar(mm);
    } catch (const std::exception& e) {
        cout << filename << ": error al leer (" << e.what() << ")\n\n";
    }
}

void enviar_matriz_y_recibir_svd(int k_total, int k_target) {
    if (N_global == 0) { cout << "Generar primero la matriz.\n"; return; }
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    sockaddr_in serv{};
    serv.sin_family = AF_INET;
    serv.sin_port = htons(CLIENT_PORT);
    inet_pton(AF_INET, "127.0.0.1", &serv.sin_addr);
    if (connect(sock, (sockaddr*)&serv, sizeof(serv)) < 0) { perror("connect"); close(sock); return; }

    const uint64_t min_workers = 2;
    const int max_attempts = 3;
    int attempt = 0;
    bool proceed = false;
    while (attempt < max_attempts) {
        MsgHeader chk(ID_H, 0, 0);
        if (!send_all(sock, &chk, sizeof(chk))) { cerr << "No se pudo enviar ID_H\n"; close(sock); return; }
        MsgHeader resp;
        if (!recv_all(sock, &resp, sizeof(resp)) || resp.id != ID_H) {
            cerr << "No se pudo recibir disponibilidad de workers\n";
            close(sock); return;
        }
        uint64_t available = resp.a;
        if (available >= min_workers) { proceed = true; break; }
        if (available == 1) {
            cout << "Solo hay 1 worker disponible. ¿Deseas continuar de todos modos? (s/n): ";
            char ans='n'; cin >> ans;
            if (ans == 's' || ans == 'S') { proceed = true; break; }
        }
        ++attempt;
        if (attempt < max_attempts) {
            cout << "Workers disponibles ("<<available<<") < "<<min_workers<<". Reintentando en 5s...\n";
            this_thread::sleep_for(std::chrono::seconds(5));
        }
    }
    if (!proceed) {
        cout << "No se enviará la matriz porque no hay suficientes workers.\n";
        close(sock);
        return;
    }

    MsgHeader h;
    h.id = ID_A; h.a = N_global; h.b = k_total;
    if (!send_all(sock, &h, sizeof(h))) { cerr << "Error al enviar el encabezado\n"; close(sock); return; }

    auto mm = mmap_abrir_lectura(MATRIX_FILE, N_global, N_global);
    size_t bytes = mm.bytes;
    size_t sent = 0;
    const char* p = (const char*)mm.data;
    while (sent < bytes) {
        size_t chunk = (bytes - sent) > CHUNK ? CHUNK : (bytes - sent);
        if (!send_all(sock, p + sent, chunk)) { cerr << "Envio fallido\n"; munmap(mm.data, mm.bytes); close(mm.fd); close(sock); return; }
        sent += chunk;
    }
    cout << "Enviar matriz (" << sent << " bytes)\n";
    munmap(mm.data, mm.bytes); close(mm.fd);

    MsgHeader hk(ID_K, (uint64_t)k_target, 0);
    send_all(sock, &hk, sizeof(hk));

    while (true) {
        MsgHeader rh;
        if (!recv_all(sock, &rh, sizeof(rh))) { cerr << "recv header failed\n"; break; }
        if (rh.id == ID_UT) { 
            uint64_t rows = rh.a, cols = rh.b;
            cout << "Receiving U ("<<rows<<"x"<<cols<<")\n";
            string fname = "U.bin";
            int fd = open(fname.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0666);
            ftruncate(fd, rows*cols*sizeof(float));
            void* map = mmap(nullptr, rows*cols*sizeof(float), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
            recv_all(sock, map, rows*cols*sizeof(float));
            msync(map, rows*cols*sizeof(float), MS_SYNC);
            munmap(map, rows*cols*sizeof(float));
            close(fd);
            cout << "Saved U to " << fname << "\n";
        } else if (rh.id == ID_S) { 
            uint64_t len = rh.a;
            cout << "Receiving singular values (" << len << ")\n";
            string fname = "S.bin";
            int fd = open(fname.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0666);
            ftruncate(fd, len * sizeof(float));
            void* map = mmap(nullptr, len*sizeof(float), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
            recv_all(sock, map, len*sizeof(float));
            msync(map, len*sizeof(float), MS_SYNC);
            munmap(map, len*sizeof(float));
            close(fd);
            cout << "Saved S to " << fname << "\n";
        } else if (rh.id == ID_VT) {
            uint64_t rows = rh.a, cols = rh.b;
            cout << "Receiving V^T ("<<rows<<"x"<<cols<<")\n";
            string fname = "VT.bin";
            int fd = open(fname.c_str(), O_CREAT | O_RDWR | O_TRUNC, 0666);
            ftruncate(fd, rows*cols*sizeof(float));
            void* map = mmap(nullptr, rows*cols*sizeof(float), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
            recv_all(sock, map, rows*cols*sizeof(float));
            msync(map, rows*cols*sizeof(float), MS_SYNC);
            munmap(map, rows*cols*sizeof(float));
            close(fd);
            cout << "Saved V^T to " << fname << "\n";
        } else if (rh.id == ID_DONE) {
            cout << "Server done. Closing.\n";
            break;
        } else {
            cout << "Unknown msg id from server: " << rh.id << "\n";
            break;
        }
    }

    close(sock);
}

void imprimirMenu(){
    cout << "\nMenu:";
    cout << "\n1) Generar matriz NxN aleatoria";
    cout << "\n2) Cargar matriz desde archivo binario";
    cout << "\n3) Enviar matriz al servidor para distribuir SVD";
    cout << "\n4) View matrix head";
    cout << "\n5) Exit";
    cout << "\nOption: ";
}

int main() {
    cout << "Cliente SVD Distribuido\n";
    while (true) {
        imprimirMenu();
        int op; if (!(cin >> op)) break;
        if (op == 1) {
            cout << "n = "; uint64_t n; cin >> n; generar_matriz(n);
        } else if (op == 2) {
            cout << "Ruta del archivo binario: ";
            string ruta; cin >> ruta;
            cargar_matriz_desde_bin(ruta);
        } else if (op == 3) {
            cout << "k (valores singulares, k<= n) = "; int k; cin >> k;
            cout << "p (oversampling) = "; int p ; cin >> p;
            int k_total = k + p;

            if (k <= 0 || (N_global && k > (int)N_global)) {
                cout << "k debe ser >0 y <= n (" << N_global << ")\n";
                continue;
            }
            if (k_total > (int)N_global) {
                cout << "k + p debe ser <= n (" << N_global << ")\n";
                continue;
            }
            enviar_matriz_y_recibir_svd(k_total, k);
        } else if (op == 4) {
            mostrar_matriz("matrix.bin");
            mostrar_matriz("U.bin");
            mostrar_matriz("VT.bin");
            mostrar_sigma("S.bin");
        } else break;
    }
    return 0;
}
