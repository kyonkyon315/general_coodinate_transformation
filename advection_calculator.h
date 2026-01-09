#ifndef ADVECTION_CALCULATOR_H
#define ADVECTION_CALCULATOR_H
#if defined(__GNUC__) || defined(__clang__)
#define ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#define ALWAYS_INLINE inline
#endif
#include "independent.h"
#include "utils/arg_changer.h"
#include <type_traits>
using Value = double;

template<typename TargetFunction,typename Advections,typename Jacobian> 
class AdvectionCalculator{
private:
    const Advections& advections;
    const Jacobian& jacobian;

    static constexpr int dimension = TargetFunction::get_dimension();

    //連鎖率を用いて、計算空間でのフラックスを計算します。
    template<int I,int Target_Dim,typename... Ints>
    ALWAYS_INLINE
    Value advection_in_calc_space_helper(Ints... indices)const{
        using E = typename Jacobian::element_t<Target_Dim,I>;
        if constexpr(I==dimension-1){
            if constexpr(std::is_same_v<E, Independent>){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return 0.;
            }
            return jacobian.template get_element<Target_Dim,I>().at(indices...)
                    *advections.template get_object<I>().at(indices...);
        }
        else{
            if constexpr(std::is_same_v<E, Independent>){
                //こんなことしなくても最適化で０の項はなくなるかも。
                return advection_in_calc_space_helper<I+1,Target_Dim>(indices...);
            }
            return jacobian.template get_element<Target_Dim,I>().at(indices...)
                    *advections.template get_object<I>().at(indices...)
                    +advection_in_calc_space_helper<I+1,Target_Dim>(indices...);
        }
    }

public:
    template<int Target_Dim,typename... Ints>
    ALWAYS_INLINE
    Value flux_p_half(Value dt,Ints... indices)const{
        const Value v = advection_in_calc_space_helper<0,Target_Dim>(indices...);
        const Value v_p_1 = Utility::arg_changer<Target_Dim,1>(
                            [this](auto... idxs) -> Value{
                                
                                return this->advection_in_calc_space_helper<0,Target_Dim>(idxs...);
                            },
                            indices...
                        );
        const Value v_p_half = (v + v_p_1)/2.;
        return v_p_half*dt;//一次精度
    }

    template<int Target_Dim,typename... Ints>
    ALWAYS_INLINE
    Value flux_m_half(Value dt,Ints... indices)const{
        const Value v = advection_in_calc_space_helper<0,Target_Dim>(indices...);
        const Value v_m_1 = Utility::arg_changer<Target_Dim,-1>(
                            [this](auto... idxs) -> Value{
                                return this->advection_in_calc_space_helper<0,Target_Dim>(idxs...);
                            },
                            indices...
                        );
        const Value v_m_half = (v + v_m_1)/2.;
        return v_m_half*dt;//一次精度
    }
    AdvectionCalculator(
        TargetFunction& target_func,
        const Advections& advections,
        const Jacobian& jacobian
    ):
        advections(advections),
        jacobian(jacobian)
    {}
};

#endif //ADVECTION_CALCULATOR_H
