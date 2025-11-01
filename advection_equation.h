#ifndef ADVECTION_EQUATION_H
#define ADVECTION_EQUATION_H
#include "independent.h"
#include "utils/arg_changer.h"
#include <type_traits>
#include <array>
#include <utility>

using Value = double;
template<typename TargetFunction,typename Operators,typename Advections,typename Jacobian,typename Scheme,typename BoundaryCondition>
class AdvectionEquation
{
private:
    TargetFunction& target_func;
    const Operators& operators;
    const Advections& advections;
    const Jacobian& jacobian;
    const Scheme& scheme;
    const BoundaryCondition& boundary_condition;

    static constexpr int dimension = TargetFunction::get_dimension();
    static constexpr int L = Scheme::used_id_left;
    static constexpr int R = Scheme::used_id_right;

    //連鎖率を用いて、計算空間でのフラックスを計算します。
    template<int I,int Target_Dim,typename... Ints>
    Value advection_in_calc_space_helper(Ints... indices){
        using E = typename Jacobian::element_t<Target_Dim,I>;
        if constexpr(I==dimension-1){
            if constexpr(std::is_same_v<E, Independent>){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return 0.;
            }
            return jacobian.get_element<Target_Dim,I>().at(indices...)
                    *advections.get_object<I>().at(indices...);
        }
        else{
            if constexpr(std::is_same_v<E, Independent>){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return advection_in_calc_space_helper<I+1,Target_Dim>(indices...);
            }
            return jacobian.get_element<Target_Dim,I>().at(indices...)
                    *advections.get_object<I>().at(indices...)
                    +advection_in_calc_space_helper<I+1,Target_Dim>(indices...);
        }
    }

    template<int Target_Dim,typename... Ints>
    Value advection_in_calc_space(Ints... indices){
        return advection_in_calc_space_helper<0,Target_Dim>(indices...);
    }

    template<int Depth,int Dim,int Target_Dim,typename... Ints, std::size_t... Is>
    void solve_helper(Value dt,Ints... indices,std::index_sequence<Is...>){
        if constexpr(Depth == Dim){
            constexpr int stencil_offsets[] = { (int(Is) + L) ... };
            Value advection = advection_in_calc_space<Target_Dim>(indices...);
            Value nyu = - dt * advection;
            Value df = limiter.calc_df(Utility::arg_changer<Target_Dim,stencil_offsets[Is]>(target_func.at,indices...)...,
                            nyu);
            //dfをどうするかは後で考える。
        }
        else{
            constexpr int axis_len = TargetFunction::sizes[Depth];
            for(int i=0;i<axis_len;++i){
                solve_helper<Depth+1,Dim,Target_Dim>(dt,indices...,i,std::make_index_sequence<R-L+1>{});
            }
        }
    }


public:
    AdvectionEquation(
        TargetFunction& target_func,
        const Operators& operators,
        const Advections& advections,
        const Jacobian& jacobian,
        const Scheme& scheme,
        const BoundaryCondition& boundary_condition
    ):
        target_func(target_func),
        operators(operators),
        advections(advections),
        jacobian(jacobian),
        scheme(scheme),
        boundary_condition(boundary_condition)
    {
        static_assert(target_func.get_dimension()==operators.get_num_objects());
        static_assert(target_func.get_dimension()==advections.get_num_objects());
    }

    template<typename CalcAxis>
    void solve(Value dt){
        constexpr int target_dim = CalcAxis::label;
        solve_helper<0,dimension,target_dim>(dt,std::make_index_sequence<R-L+1>{});
    }
};
#endif //ADVECTION_EQUATION_H