#ifndef MATRIX_IO_HPP
#define MATRIX_IO_HPP

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <stdexcept>
#include <iostream>

struct MMapMatrix {
    float* data;
    uint64_t rows;
    uint64_t cols;
    size_t bytes;
    int fd;
    std::string path;
};

inline MMapMatrix mmap_create(const std::string &path, uint64_t rows, uint64_t cols) {
    uint64_t elems = rows * cols;
    size_t bytes = elems * sizeof(float);
    int fd = open(path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0666);
    if (fd < 0) throw std::runtime_error("open failed");
    if (ftruncate(fd, bytes) != 0) throw std::runtime_error("ftruncate failed");
    void* m = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) throw std::runtime_error("mmap failed");
    return {(float*)m, rows, cols, bytes, fd, path};
}

inline MMapMatrix mmap_open_read(const std::string &path, uint64_t rows, uint64_t cols) {
    uint64_t elems = rows * cols;
    size_t bytes = elems * sizeof(float);
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw std::runtime_error("open failed");
    void* m = mmap(nullptr, bytes, PROT_READ, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) throw std::runtime_error("mmap failed");
    return {(float*)m, rows, cols, bytes, fd, path};
}

inline void mmap_close(MMapMatrix &mm) {
    munmap(mm.data, mm.bytes);
    close(mm.fd);
}

#endif
