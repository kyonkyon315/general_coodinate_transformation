#ifndef ADVECTION_EQUATION_H
#define ADVECTION_EQUATION_H
#include "independent.h"
#include "utils/arg_changer.h"
#include "slice.h"
#include <type_traits>
#include <array>
#include <utility>
#include <omp.h>

using Value = double;
template<typename TargetFunction,typename Operators,typename Advections,typename Jacobian,typename Scheme,typename BoundaryCondition>
class AdvectionEquation
{
private:
    TargetFunction& target_func;
    TargetFunction func_buffer;
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
    template<int Depth,int Dim,typename...Slices,int Target_Dim,typename... Ints>
    void solve_helper(const Value dt, Ints... indices){
        if constexpr(Depth == Dim){
            const Value df = solve_leaf<Target_Dim>(dt, indices...);
            //増加分を保存
            func_buffer.at(indices...) = df;
        }
        else if constexpr(Depth==0){
            using CurrentSlice = std::tuple_element_t<Depth, std::tuple<Slices...>>;
            //[@TODO]Dim==1のときはomp発動しないようにした方がいいかも。
            constexpr int axis_len = TargetFunction::shape[Depth];
            #pragma omp parallel for
            for(int i=CurrentSlice::START_val;i<CurrentSlice::End_val;++i){
                // 再帰：indices に i を追加
                solve_helper<Depth+1,Dim,Target_Dim>(dt, indices..., i);
            }
        }
        else{
            using CurrentSlice = std::tuple_element_t<Depth, std::tuple<Slices...>>;
            constexpr int axis_len = TargetFunction::shape[Depth];
            for(int i=CurrentSlice::START_val;i<CurrentSlice::End_val;++i){
                // 再帰：indices に i を追加
                solve_helper<Depth+1,Dim,Slices,Target_Dim>(dt, indices..., i);
            }
        }
    }

    // leaf: index_sequence を生成して「index-first」ヘルパを呼ぶ
    template<int Target_Dim, typename... Ints>
    Value solve_leaf(const Value dt, Ints... indices){
        constexpr std::size_t stencil_size = (R - L + 1);
        return solve_leaf_impl_indices<Target_Dim>(dt, std::make_index_sequence<stencil_size>{}, indices...);
    }

    // helper: index_sequence を最初の引数に置く（※これで推論が安定）
    template<int Target_Dim, std::size_t... Is, typename... Ints>
    Value solve_leaf_impl_indices(const Value dt, std::index_sequence<Is...>, Ints... indices){
        // コンパイル時に確定するオフセット配列
        static constexpr int stencil_offsets[] = { (int(Is) + L)... };

        // チェーンルールによる advection 計算
        Value advection = this->advection_in_calc_space<Target_Dim>(indices...);
        Value nyu = - dt * advection;

        Value df = this->scheme.calc_df(
            Utility::arg_changer<Target_Dim, stencil_offsets[Is]>(
                [this](auto... idxs) -> Value {
                    return this->target_func.at(idxs...);
                },
                indices...
            )...,
            nyu,nyu//本当はnyu_minus_half, nyu_plus_halfを入れないといけない。
        );
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
        constexpr int target_dim = CalcAxis::label;
        // dfを計算してfunc_bufferに格納
        solve_helper<0, dimension, target_dim>(dt);

        //dfをfに足し込む（関数の更新）
        target_func.add(func_buffer);
    }

    template<typename CalcAxis,typename... Slices>
    void solve(const Value dt){
        constexpr int target_dim = CalcAxis::label;
        static_assert(sizeof...(Slices)==dimension,"Slices の数が次元数と一致しません。");
        
        solve_helper<0,dimension,Slices...,target_dim>(dt);
        //dfをfに足し込む（関数の更新）
        target_func.add<Slices...>(func_buffer);
    }
};
#endif //ADVECTION_EQUATION_H