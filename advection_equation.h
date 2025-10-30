#ifndef ADVECTION_EQUATION_H
#define ADVECTION_EQUATION_H
#include "independent.h"
using Value = double;
template<typename TargetFunction,typename Operators,typename Advections,typename Jacobian>
class AdvectionEquation
{
private:
    TargetFunction& target_func;
    const Operators& operators;
    const Advections& advections;
    const Jacobian& jacobian;

    static constexpr int dimension = TargetFunction::get_dimension();

    //連鎖率を用いて、計算空間でのフラックスを計算します。
    template<int I,typename CalcAxis>
    Value flux_in_calc_space_helper(int... indices){
        if constexpr(I==dimension-1){
            if constexpr(is_same<
                                jacobian.get_calc_diff_phys<CalcAxis::label,I>(),
                                Independent
                                >){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return 0.;
            }
            return jacobian.get_calc_diff_phys<CalcAxis::label,I>().at(indices...)
                    *advections.get_object<I>().at(indices);
        }
        else{
            if constexpr(is_same<jacobian.get_calc_diff_phys<CalcAxis::label,I>(),Independent>){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return flux_in_calc_space_helper<I+1,CalcAxis>(indices);
            }
            return jacobian.get_calc_diff_phys<CalcAxis::label,I>().at(indices...)
                    *advections.get_object<I>().at(indices)
                    +flux_in_calc_space_helper<I+1,CalcAxis>(indices);
        }
    }

    template<typename CalcAxis>
    Value flux_in_calc_space(int... indices){
        return flux_in_calc_space_helper<0,CalcAxis>(indices);
    }



public:
    AdvectionEquation(
        TargetFunction& target_func,
        const Operators& operators,
        const Advections& advections,
        const Jacobian& jacobian):
        target_func(target_func),
        operators(operators),
        advections(advections),
        jacobian(jacobian)
    {
        static_assert(target_func.get_dimension()==operators.get_num_objects());
        static_assert(target_func.get_dimension()==advections.get_num_objects());
    }

    template<typename CalcAxis>
    void solve(Value dt){
        
    }
};
#endif //ADVECTION_EQUATION_H