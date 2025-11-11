#ifndef ADVECTION_EQUATION_H
#define ADVECTION_EQUATION_H
#include "independent.h"
#include "utils/arg_changer.h"
#include <type_traits>
#include <array>
#include <utility>
#include <omp.h>
#include <cmath>

using Value = double;
template<typename TargetFunction,typename Operators,typename Advections,typename Jacobian,typename Scheme,typename BoundaryCondition>
class AdvectionEquation
{
private:
    TargetFunction& target_func;
    TargetFunction func_buffer1;
    TargetFunction func_buffer2;
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

    template<int Target_Dim,typename... Ints>
    Value advection_in_calc_space(Ints... indices){
        return advection_in_calc_space_helper<0,Target_Dim>(indices...);
    }

    // 再帰ヘルパ：indices を集める
    template<int Depth,int Dim,int Target_Dim,typename... Ints>
    void solve_helper(Value dt, Ints... indices){
        if constexpr(Depth == Dim){
            const Value df = solve_leaf<Target_Dim>(dt, indices...);
            //増加分を保存
            func_buffer1.at(indices...) = df;
        }
        else if constexpr(Depth==0){
            //[@TODO]Dim==1のときはomp発動しないようにした方がいいかも。
            constexpr int axis_len = TargetFunction::shape[Depth];
            #pragma omp parallel for
            for(int i=0;i<axis_len;++i){
                // 再帰：indices に i を追加
                solve_helper<Depth+1,Dim,Target_Dim>(dt, indices..., i);
            }
        }
        else{
            constexpr int axis_len = TargetFunction::shape[Depth];
            for(int i=0;i<axis_len;++i){
                // 再帰：indices に i を追加
                solve_helper<Depth+1,Dim,Target_Dim>(dt, indices..., i);
            }
        }
    }

    // leaf: index_sequence を生成して「index-first」ヘルパを呼ぶ
    template<int Target_Dim, typename... Ints>
    Value solve_leaf(Value dt, Ints... indices){
        constexpr std::size_t stencil_size = (R - L + 1);
        return solve_leaf_impl_indices<Target_Dim>(dt, std::make_index_sequence<stencil_size>{}, indices...);
    }

    // helper: index_sequence を最初の引数に置く（※これで推論が安定）
    template<int Target_Dim, std::size_t... Is, typename... Ints>
    Value solve_leaf_impl_indices(Value dt, std::index_sequence<Is...>, Ints... indices){
        // コンパイル時に確定するオフセット配列
        static constexpr int stencil_offsets[] = { (int(Is) + L)... };

        // チェーンルールによる advection 計算
        Value advection_p_1 = Utility::arg_changer<Target_Dim,1>(
                            [this](auto... idxs) -> Value{
                                return this->advection_in_calc_space<Target_Dim>(idxs...);
                            },
                            indices...
                        );
        Value advection = Utility::arg_changer<Target_Dim,0>(
                            [this](auto... idxs) -> Value{
                                return this->advection_in_calc_space<Target_Dim>(idxs...);
                            },
                            indices...
                        );
        Value advection_m_1 = Utility::arg_changer<Target_Dim,-1>(
                            [this](auto... idxs) -> Value{
                                return this->advection_in_calc_space<Target_Dim>(idxs...);
                            },
                            indices...
                        );
        Value nyu_p_half = - dt * (advection+advection_p_1)/2.;
        Value nyu_m_half = - dt * (advection+advection_m_1)/2.;
        
        Value df = this->scheme.calc_df(
            Utility::arg_changer<Target_Dim, stencil_offsets[Is]>(
                [this](auto... idxs) -> Value {
                    return this->target_func.at(idxs...);
                },
                indices...
            )...,
            nyu_p_half,nyu_m_half
        );

        if(std::isnan(df)){
            ((std::cout<<indices<<" "),...);
            std::cout<<"nyu_p_half: "<<nyu_p_half<<"\n";
            std::cout<<"nyu_m_half: "<<nyu_m_half<<"\n";
            std::cout<<"\n";
            throw 1;
        }

        return df;
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
        //以下ルンゲクッタ2の実装
        constexpr int target_dim = CalcAxis::label;
        // dfを計算してfunc_buffer1に格納
        solve_helper<0, dimension, target_dim>(dt/2.);

        func_buffer2.copy(target_func);

        //dfをfに足し込む（関数の更新）
        target_func.add(func_buffer1);
        // dfを計算してfunc_buffer1に格納
        solve_helper<0, dimension, target_dim>(dt);

        target_func.copy(func_buffer2);
        //dfをfに足し込む（関数の更新）
        target_func.add(func_buffer1);

    }
};
#endif //ADVECTION_EQUATION_H