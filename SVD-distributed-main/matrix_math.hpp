#ifndef MATRIX_MATH_HPP
#define MATRIX_MATH_HPP

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cassert>
using namespace std;

void generate_omega_mmap(uint64_t seed, int n, int k, const string &fileOmega) {
    auto O = mmap_create(fileOmega, n, k);

    mt19937 rng(seed);
    uniform_real_distribution<float> dist(-1.0f, 1.0f);

    uint64_t total = (uint64_t)n * k;
    for (uint64_t i = 0; i < total; i++)
        O.data[i] = dist(rng);
    mmap_close(O);
    cout << "[worker] Omega builded: " << fileOmega << "\n";
}

void matmul_A_Omega_mmap(const string &fileA,int rows, int cols,
                         const string &fileOmega,int k,
                         const string &fileY){

    auto A = mmap_open_read(fileA, rows, cols);
    auto O = mmap_open_read(fileOmega, cols, k);
    auto Y = mmap_create(fileY, rows, k);

    for (int r = 0; r < rows; ++r) {
        const float* arow = A.data + (uint64_t)r * cols;
        for (int c = 0; c < k; ++c) {
            float sum = 0;
            for (int j = 0; j < cols; ++j)
                sum += arow[j] * O.data[(uint64_t)j * k + c];
            Y.data[(uint64_t)r * k + c] = sum;
        }
    }

    mmap_close(A);
    mmap_close(O);
    mmap_close(Y);

    cout << "[worker] Y builded: " << fileY << "\n";
}


void qr_mmap(const string &fileY,int rows, int k,
                 const string &fileQi,const string &fileRi){
    auto Y = mmap_open_read(fileY, rows, k);
    auto Q = mmap_create(fileQi, rows, k);
    auto R = mmap_create(fileRi, k, k);

    vector<float> v(rows);

    for (int j = 0; j < k; ++j) {

        // v = Y[:, j]
        for (int i = 0; i < rows; ++i)
            v[i] = Y.data[(uint64_t)i * k + j];

        // --- Proyecciones ---
        for (int i = 0; i < j; ++i) {
            float rij = 0;

            // rij = Q[:,i]^T * v
            for (int t = 0; t < rows; ++t)
                rij += Q.data[(uint64_t)t * k + i] * v[t];

            R.data[(uint64_t)i * k + j] = rij;

            // v -= rij * Q[:,i]
            for (int t = 0; t < rows; ++t)
                v[t] -= rij * Q.data[(uint64_t)t * k + i];
        }

        float norm = 0;
        for (int t = 0; t < rows; ++t)
            norm += v[t] * v[t];

        norm = sqrt(norm);

        if (norm < 1e-12f) {
            R.data[(uint64_t)j * k + j] = 0;
            for (int t = 0; t < rows; ++t)
                Q.data[(uint64_t)t * k + j] = 0;
        } else {
            R.data[(uint64_t)j * k + j] = norm;
            for (int t = 0; t < rows; ++t)
                Q.data[(uint64_t)t * k + j] = v[t] / norm;
        }
    }

    mmap_close(Y);
    mmap_close(Q);
    mmap_close(R);
}

inline void jacobi_eigen_inplace(float* A, int n, float* Util, float* Lambda,
                                 int max_iter = 100, float tol = 1e-7f){
    // Initialize Util = I
    for (int i = 0; i < n * n; ++i) Util[i] = 0.0f;
    for (int i = 0; i < n; ++i) Util[i * n + i] = 1.0f;

    // Work on A_copy in-place: copy A to a local buffer since caller may pass mmap read-only
    vector<float> M((size_t)n * n);
    for (int i = 0; i < n * n; ++i) M[i] = A[i];

    // Jacobi iterations
    for (int iter = 0; iter < max_iter; ++iter) {
        // find largest off-diagonal element
        int p = 0, q = 1;
        float maxv = 0.0f;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                float v = fabs(M[i * n + j]);
                if (v > maxv) { maxv = v; p = i; q = j; }
            }
        }
        if (maxv < tol) break;

        float app = M[p * n + p];
        float aqq = M[q * n + q];
        float apq = M[p * n + q];

        float phi = 0.5f * atan2f(2.0f * apq, (aqq - app));
        float c = cosf(phi), s = sinf(phi);

        // update matrix M: perform rotation in (p,q) plane
        for (int i = 0; i < n; ++i) {
            if (i == p || i == q) continue;
            float mip = M[i * n + p];
            float miq = M[i * n + q];
            M[i * n + p] = c * mip - s * miq;
            M[p * n + i] = M[i * n + p];
            M[i * n + q] = s * mip + c * miq;
            M[q * n + i] = M[i * n + q];
        }

        float new_pp = c*c*app - 2.0f*s*c*apq + s*s*aqq;
        float new_qq = s*s*app + 2.0f*s*c*apq + c*c*aqq;
        M[p * n + p] = new_pp;
        M[q * n + q] = new_qq;
        M[p * n + q] = 0.0f;
        M[q * n + p] = 0.0f;

        // update Util columns p and q
        for (int i = 0; i < n; ++i) {
            float up = Util[i * n + p];
            float uq = Util[i * n + q];
            Util[i * n + p] = c * up - s * uq;
            Util[i * n + q] = s * up + c * uq;
        }
    }

    // write eigenvalues (diagonal) into Lambda
    for (int i = 0; i < n; ++i) Lambda[i] = M[i * n + i];
}


void compute_Vj_mmap(const string &UtilFile,const string &SigmaInvFile,
    const string &BiFile,const string &VjOutFile,int k, int ncols){

    MMapMatrix Util = mmap_open_read(UtilFile, k, k);
    MMapMatrix SigmaInv = mmap_open_read(SigmaInvFile, k, 1);
    MMapMatrix Bi = mmap_open_read(BiFile, k, ncols);
    MMapMatrix Vj = mmap_create(VjOutFile, k, ncols);

    // Compute M = Util^T * Bi  -> result k x ncols (store into Vj temporarily)
    // Util^T: (k x k), Bi: (k x ncols)
    for (int i = 0; i < k; ++i) {          // row of Util^T (i) == column of Util
        for (int j = 0; j < ncols; ++j) {
            float s = 0.0f;
            for (int t = 0; t < k; ++t) {
                // Util^T[i,t] = Util[t,i]
                s += Util.data[(size_t)t * k + i] * Bi.data[(size_t)t * ncols + j];
            }
            Vj.data[(size_t)i * ncols + j] = s;
        }
    }

    // Multiplicar cada fila por su correspondiente valor en SigmaInv: Vj[i, :] *= SigmaInv[i]
    for (int i = 0; i < k; ++i) {
        float scale = SigmaInv.data[i];
        if (scale == 1.0f) continue;
        for (int j = 0; j < ncols; ++j)
            Vj.data[(size_t)i * ncols + j] *= scale;
    }

    msync(Vj.data, Vj.bytes, MS_SYNC);

    mmap_close(Util);
    mmap_close(SigmaInv);
    mmap_close(Bi);
    mmap_close(Vj);
}

void compute_Ui_mmap(const string &Qi_file,const string &Utilde_file,
    const string &Ui_file,int mi, int k){

    MMapMatrix Qmm = mmap_open_read(Qi_file, mi, k);
    MMapMatrix Utilde = mmap_open_read(Utilde_file, k, k);
    MMapMatrix Uimm = mmap_create(Ui_file, mi, k);

    // --- Multiplicación U_i = Q_i * Utilde ---
    // U_i (m_i x k) = Q_i (m_i x k) * Utilde (k x k)

    for (int r = 0; r < mi; ++r){
        float* qrow = &Qmm.data[(size_t)r * k];
        float* urow = &Uimm.data[(size_t)r * k];

        for (int j = 0; j < k; ++j)
        {
            float sum = 0.0f;
            for (int t = 0; t < k; ++t)
                sum += qrow[t] * Utilde.data[(size_t)t * k + j];

            urow[j] = sum;
        }
    }

    mmap_close(Qmm);
    mmap_close(Utilde);
    mmap_close(Uimm);

    cout << "[worker] computed Ui = Q_i * Utilde and stored in "
    << Ui_file << " (" << mi << " x " << k << ")\n";
}

#endif
