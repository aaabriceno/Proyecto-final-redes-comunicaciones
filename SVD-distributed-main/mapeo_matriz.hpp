#ifndef MAPEO_MATRIZ_HPP
#define MAPEO_MATRIZ_HPP

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <stdexcept>
#include <iostream>
using namespace std;

struct MMapMatrix {
    float* data;
    uint64_t filas;
    uint64_t columnas;
    size_t bytes;
    int fd;
    string path;
};

inline MMapMatrix mmap_create(const string &path, uint64_t filas, uint64_t columnas) {
    uint64_t elems = filas * columnas;
    size_t bytes = elems * sizeof(float);
    int fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) throw runtime_error("open failed");
    if (ftruncate(fd, bytes) != 0) throw runtime_error("ftruncate failed");
    void* m = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) throw runtime_error("mmap failed");
    return {(float*)m, filas, columnas, bytes, fd, path};
}

inline MMapMatrix mmap_open_read(const string &path, uint64_t filas, uint64_t columnas) {
    uint64_t elems = filas * columnas;
    size_t bytes = elems * sizeof(float);
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw runtime_error("open failed");
    void* m = mmap(nullptr, bytes, PROT_READ, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) throw runtime_error("mmap failed");
    return {(float*)m, filas, columnas, bytes, fd, path};
}

inline void mmap_close(MMapMatrix &mm) {
    munmap(mm.data, mm.bytes);
    close(mm.fd);
}

#endif
