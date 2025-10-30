#ifndef JACOBIAN_H
#define JACOBIAN_H
#include <tuple>
template<typename... Xi_diff_x>
class Jacobian
{
private:
    std::tuple<const std::add_lvalue_reference_t<Xi_diff_x>...> elements;
    static constexpr int num_args = sizeof...(elements);
    static constexpr int num_label = []{
    constexpr int n = sizeof...(Xi_diff_x);
    for (int i = 1; i <= 20; ++i) {
        if (i * i == n) return i;
    }
    static_assert(false, "テンプレート引数の数は平方数である必要があります。");
    return 0; // never reached
}();
public:
    Jacobian(const Xi_diff_x&... args):
        elements(args...)
    {}
    template<int I,int J>
    const auto& operator()(){
        return std::get<num_label*I+J>(elements);
    }

};
#endif //JACOBIAN_H