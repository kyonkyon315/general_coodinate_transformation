#ifndef CALC_CURRENT_IN_X_AND_Y
#define CALC_CURRENT_IN_X_AND_Y

#include <tuple>
#include <utility>

// 一次元専門
template<typename Current, typename DistFunction, typename Operators>
class CalcCurrent_1d {
public:
    static constexpr int dim      = DistFunction::get_dimension();
    static constexpr int real_dim = Current::get_dimension();
    static constexpr int velo_dim = dim - real_dim;

    static_assert(real_dim == 1,
        "this current solver is only for 1 dimensional.\n");
    static_assert(velo_dim >= 2,
        "the velocity dimension must be more than 2.\n");
    static_assert(dim == Operators::get_num_objects(),
        "total dimension and size of Operators mismatch.\n");
    static_assert(Current::shape[0] == DistFunction::shape[0],
        "the shape in real space mismatch.\n");

private:
    Current&       current;
    DistFunction& distfunction;
    Operators&    operators;

public:
    CalcCurrent_1d(Current& current,
                   DistFunction& distfunction,
                   Operators& operators)
        : current(current),
          distfunction(distfunction),
          operators(operators) {}

private:
    // Depth: 全次元走査
    // R...  : real-space indices
    // V...  : velocity-space indices
     // ---- velocity recursion ----
    template<int Depth, typename... V>
    inline void vel_loop(int x, V... v) {
        if constexpr (Depth < velo_dim) {
            constexpr int axis = real_dim + Depth;
            for (int i = 0; i < DistFunction::shape[axis]; ++i) {
                vel_loop<Depth + 1>(x, v..., i);
            }
        } else {
            // ---- leaf ----
            auto f = distfunction.at(x, v...);


            //電流はyee格子により半グリッドずれているので、以下のような実装になる。
            //j(i)=j(Δx (i+1/2)) = (v f(i) +v f(i+1))/2
            //
            //言い換えると、
            //for i in range
            //  j(i) += v f(i)/2
            //  j(i-1) += v f(i)/2
            
            current.at(x-1).x +=
                operators.template get_object<1>().translate(x, v...) * f/2.;
            current.at(x).x +=
                operators.template get_object<1>().translate(x, v...) * f/2.;

            current.at(x-1).y +=
                operators.template get_object<2>().translate(x, v...) * f/2.;
            current.at(x).y +=
                operators.template get_object<2>().translate(x, v...) * f/2.;
        }
    }
public:
    void calc() {
        for (int x = 0; x < DistFunction::shape[0]; ++x) {
            current.at(x).x = 0.0;
            current.at(x).y = 0.0;
            vel_loop<0>(x);
        }
    }
};

#endif // CALC_CURRENT_IN_X_AND_Y
