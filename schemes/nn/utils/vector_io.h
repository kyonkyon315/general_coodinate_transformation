#ifndef VECTOR_IO_H
#define VECTOR_IO_H
#include <vector>
#include <fstream>

template <class T>
void save_vector(const std::vector<T>& vec, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    size_t size = vec.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));      // サイズ保存
    ofs.write(reinterpret_cast<const char*>(vec.data()), sizeof(T) * size); // 本体保存
}

template <class T>
std::vector<T> load_vector(const std::string& filename) {
    std::ifstream ifs(filename, std::ios::binary);
    size_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));  // サイズ読み込み

    std::vector<T> vec(size);
    ifs.read(reinterpret_cast<char*>(vec.data()), sizeof(T) * size); // 本体読み込み
    return vec;
}

#endif //VECTOR_IO_H