#ifndef BOUNDARY_MANAGER_H
#define BOUNDARY_MANAGER_H
#include <tuple>
#include <type_traits>
#include <utility>
#include "block_id2rank.h"
#include "../utils/arg_changer.h"

#include "axis.h"
#include "BC2block_id.h"
#include "n_d_tensor_with_ghost_cell.h"

using Value = double;
struct CommInfo {
    int dst_rank;     // 送り先
    int src_rank;     // 受信元
    bool is_my_left;  // 自分の left / right
    bool is_src_left; // 相手の left / right
};


//スパコン用のboundary manager では、境界条件の交換と共に、スレッド間での通信も担当する。

template<typename TargetFunc,typename BoundaryCondition, Axis_T... Axes>
class BoundaryManager{
private:
    const BC2BlockId<BoundaryCondition, Axes...> bc2block_id;
    std::tuple<Axes...> axes_tuple;

    const BoundaryCondition& boundary_condition;
    TargetFunc& target_func;
    static constexpr int Dim = TargetFunc::shape.size();
    BlockId2Rank<Axes...> blockid2rank;
    std::array<std::array<CommInfo, 2>, Dim> comm_info;//送信先、受信元のランクなどを格納

    //Axesの入力方法が正しいか確認(開始)
    template<int I = 0>
    static constexpr bool axis_order_checker(){
        if constexpr (I == Dim) {
            return true;
        } else {
            return std::tuple_element_t<I, std::tuple<Axes...>>::label == I
                && axis_order_checker<I+1>();
        }
    }
    static_assert(
        axis_order_checker(),
        "Axes... はAxes::labelが小さい順になるように並べてください。\n"
    );
    //Axesの入力方法が正しいか確認(終了)

    template<Axis_T TargetAxis>
    void init_left_comm() {
        constexpr int target_dim = TargetAxis::label;
        int block_id = std::get<target_dim>(axes_tuple).block_id;

        auto& info = comm_info[target_dim][0]; // left

        if (block_id == 0) {
            // 境界条件あり
            auto [src_rank, is_src_left] = std::apply(
                [&](const Axes&... axes){
                    return bc2block_id.template
                        get_rcv_rank<true, TargetAxis>(axes...);
                },
                axes_tuple
            );

            auto [dst_rank, is_my_left] = std::apply(
                [&](const Axes&... axes){
                    return bc2block_id.template
                        get_snd_rank<true, TargetAxis>();
                },
                axes_tuple
            );

            info = {dst_rank, src_rank, is_my_left, is_src_left};
        }
        else {
            // 通常通信
             //右隣ブロックに送信する
            int dst_rank = std::apply(
                 [&](auto const&... axes){
                    return Utility::arg_changer<int, target_dim, 1>(
                        [this](auto... idxs) -> int {
                            return blockid2rank.template get_rank(idxs...);
                        },
                        axes.block_id...   // ← ここで初めて pack 展開できる
                    );
                },
                axes_tuple
            );

            //左隣ブロックから受信する
            int src_rank = std::apply(
                 [&](auto const&... axes){
                    return Utility::arg_changer<int, target_dim, -1>(
                        [this](auto... idxs) -> int {
                            return blockid2rank.template get_rank(idxs...);
                        },
                        axes.block_id...   // ← ここで初めて pack 展開できる
                    );
                },
                axes_tuple
            );

            info = {
                dst_rank,
                src_rank,
                false, // 自分は right 側を送る
                false  // 相手も right 側
            };
        }
    }

    
    template<Axis_T TargetAxis>
    void init_right_comm() {
        constexpr int target_dim = TargetAxis::label;
        int block_id = std::get<target_dim>(axes_tuple).block_id;

        auto& info = comm_info[target_dim][1]; // right

        if (block_id == TargetAxis::num_blocks -1) {
            // 境界条件あり
            auto [src_rank, is_src_left] = std::apply(
                [&](const Axes&... axes){
                    return bc2block_id.template
                        get_rcv_rank<false, TargetAxis>(axes...);
                },
                axes_tuple
            );

            auto [dst_rank, is_my_left] = std::apply(
                [&](const Axes&... axes){
                    return bc2block_id.template
                        get_snd_rank<false, TargetAxis>();
                },
                axes_tuple
            );

            info = {dst_rank, src_rank, is_my_left, is_src_left};
        }
        else {
            // 通常通信
             //左隣ブロックに送信する
            int dst_rank = std::apply(
                 [&](auto const&... axes){
                    return Utility::arg_changer<int, target_dim, -1>(
                        [this](auto... idxs) -> int {
                            return blockid2rank.template get_rank(idxs...);
                        },
                        axes.block_id...   // ← ここで初めて pack 展開できる
                    );
                },
                axes_tuple
            );

            //右隣ブロックから受信する
            int src_rank = std::apply(
                 [&](auto const&... axes){
                    return Utility::arg_changer<int, target_dim, 1>(
                        [this](auto... idxs) -> int {
                            return blockid2rank.template get_rank(idxs...);
                        },
                        axes.block_id...
                    );
                },
                axes_tuple
            );

            info = {
                dst_rank,
                src_rank,
                true, // 自分は left 側を送る
                true  // 相手も left 側
            };
        }
    }

    template<int I=0>
    void init_comm(){
        if constexpr(I==Dim){
            return;
        }
        else{
            init_left_comm<std::tuple_element_t<I, std::tuple<Axes...>>>();
            init_right_comm<std::tuple_element_t<I, std::tuple<Axes...>>>();

            init_comm<I+1>();
        }
    }


    template<Axis_T TargetAxis>
    void apply_helper_l(){
        constexpr int target_dim = TargetAxis::label;
        constexpr int num_ghost_grid = TargetAxis::L_ghost_length;//=R_ghost_length
        const int block_id  = std::get<target_dim>(axes_tuple).block_id;

        //左端ブロックの場合は別行動
        if(block_id == 0){
            auto [dst_rank, src_rank, is_my_left, is_src_left] = comm_info[target_dim][0];

            if(is_src_left){
                if(is_my_left){
                    target_func.template send_ghosts<TargetAxis,true,true>(dst_rank,src_rank);
                } else {
                    target_func.template send_ghosts<TargetAxis,false,true>(dst_rank,src_rank);
                }
                target_func.template set_value_sliced<std::conditional_t<Axes::label==target_dim, Slice< - num_ghost_grid, 0>, FullSlice>...>(
                    [&]<typename... Ins1>(Ins1 const &... indices_){
                        return [&]<std::size_t... Is,typename... Ins2>(std::index_sequence<Is...>, Ins2 const &... indices){
                            return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                                [&]<typename... Ins3>(Ins3... idx) -> Value {
                                    using BC = typename BoundaryCondition::template element<target_dim>;
                                    return this->target_func.template buf_at<TargetAxis,true>(
                                        BC::template left<Is>(idx...)...
                                    );
                                },
                                indices...
                            );
                        }(std::make_index_sequence<Dim>{},indices_...);
                    }
                );
            }
            else{
                if(is_my_left){
                    target_func.template send_ghosts<TargetAxis,true,false>(dst_rank,src_rank);
                } else {
                    target_func.template send_ghosts<TargetAxis,false,false>(dst_rank,src_rank);   
                }
                std::tuple<Axes...> axes_tuple_src = axis_instantiator<Axes...>(src_rank);
                target_func.template set_value_sliced<std::conditional_t<Axes::label==target_dim, Slice< - num_ghost_grid, 0>, FullSlice>...>
                (
                    [&]<typename... Ins1>(Ins1 const &... indices_){
                        return [&]<std::size_t... Is,typename... Ins2>(std::index_sequence<Is...>, Ins2 const &... indices){
                            return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                                [&]<typename... Ins3>(Ins3... idx) -> Value {
                                    using BC = typename BoundaryCondition::template element<target_dim>;
                                    return this->target_func.template buf_at<TargetAxis,false>(
                                        (BC::template left<Is>(
                                            (std::get<Is>(axes_tuple).L_id //左端idを足してグローバルインデックスに変換
                                                                    + idx
                                                                )...)
                                            -std::get<Is>(axes_tuple_src).L_id//ローカルインデックスに戻す
                                        )...
                                    );
                                },
                                indices...
                            );
                        }(std::make_index_sequence<Dim>{},indices_...);
                    }
                );
            }
        }
        else{//左端のブロックでないならば、隣の境界セルを平行移動
            auto [dst_rank, src_rank, is_my_left, is_src_left] = comm_info[target_dim][0];

            //template<typename TargetAxis,bool Is_left_send,bool Is_left_source>
            target_func.template send_ghosts<TargetAxis,false,false>(dst_rank,src_rank);

            target_func.template set_value_sliced<std::conditional_t<Axes::label==target_dim, Slice< - num_ghost_grid, 0>, FullSlice>...>(
                [&]<typename... Ins1>(Ins1 const &... indices){
                    return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                        [&]<typename... Ins2>(Ins2... idx){
                            return this->target_func.template buf_at<TargetAxis,false>(idx...);
                        },
                        indices...
                    );
                }
            );
        }
    }

    
    template<Axis_T TargetAxis>
    void apply_helper_r(){
        constexpr int target_dim = TargetAxis::label;
        constexpr int num_ghost_grid = TargetAxis::R_ghost_length;//=L_ghost_length
        const int block_id  = std::get<target_dim>(axes_tuple).block_id;

        //左端ブロックの場合は別行動
        if(block_id == TargetAxis::num_blocks - 1){
            auto [dst_rank, src_rank, is_my_left, is_src_left] = comm_info[target_dim][1];

            if(is_src_left){
                if(is_my_left){
                    target_func.template send_ghosts<TargetAxis,true,true>(dst_rank,src_rank);
                } else {
                    target_func.template send_ghosts<TargetAxis,false,true>(dst_rank,src_rank);
                }
                std::tuple<Axes...> axes_tuple_src = axis_instantiator<Axes...>(src_rank);
                target_func.template set_value_sliced<
                                        std::conditional_t<
                                            Axes::label==target_dim, 
                                            Slice<TargetAxis::num_grid, TargetAxis::num_grid + num_ghost_grid>, 
                                            FullSlice
                                        >...
                                    >
                (
                    [&]<typename... Ins1>(Ins1 const &... indices_){
                        return [&]<std::size_t... Is,typename... Ins2>(std::index_sequence<Is...>, Ins2 const &... indices){
                            return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                                [&]<typename... Ins3>(Ins3... idx) -> Value {
                                    using BC = typename BoundaryCondition::template element<target_dim>;
                                    return this->target_func.template buf_at<TargetAxis,true>(
                                        (BC::template right<Is>((std::get<Is>(axes_tuple).L_id 
                                                                    + idx
                                                                )...)
                                            -std::get<Is>(axes_tuple_src).L_id
                                        )...
                                    );
                                },
                                indices...
                            );
                        }(std::make_index_sequence<Dim>{},indices_...);
                    }
                );
            }
            else{
                if(is_my_left){
                    target_func.template send_ghosts<TargetAxis,true,false>(dst_rank,src_rank);
                } else {
                    target_func.template send_ghosts<TargetAxis,false,false>(dst_rank,src_rank);   
                }
                target_func.template set_value_sliced<
                                        std::conditional_t<
                                            Axes::label==target_dim, 
                                            Slice<TargetAxis::num_grid, TargetAxis::num_grid + num_ghost_grid>, 
                                            FullSlice
                                        >...
                                    >
                (
                    [&]<typename... Ins1>(Ins1 const &... indices_){
                        return [&]<std::size_t... Is, typename... Ins2>(std::index_sequence<Is...>, Ins2 const &... indices){
                            return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                                [&]<typename... Ins3>(Ins3 ... idx) -> Value {
                                    using BC = typename BoundaryCondition::template element<target_dim>;
                                    return this->target_func.template buf_at<TargetAxis,false>(
                                        (BC::template right<Is>((std::get<Is>(axes_tuple).L_id 
                                                                    + idx
                                                                )...)
                                            -std::get<Is>(axes_tuple).L_id
                                        )...
                                    );
                                },
                                indices...
                            );
                        }(std::make_index_sequence<Dim>{},indices_...);
                    }
                );
            }
        }
        else{//左端のブロックでないならば、隣の境界セルを平行移動
            auto [dst_rank, src_rank, is_my_left, is_src_left] = comm_info[target_dim][1];

            //template<typename TargetAxis,bool Is_left_send,bool Is_left_source>
            target_func.template send_ghosts<TargetAxis,true,true>(dst_rank,src_rank);

            target_func.template set_value_sliced<
                                        std::conditional_t<
                                            Axes::label==target_dim, 
                                            Slice<TargetAxis::num_grid, TargetAxis::num_grid + num_ghost_grid>, 
                                            FullSlice
                                        >...
                                    >
            (
                [&]<typename... Ins1>(Ins1 const &... indices){
                    return Utility::arg_changer<Value, target_dim, TargetAxis::num_grid>(
                        [&]<typename... Ins2>(Ins2... idx){
                            return this->target_func.template buf_at<TargetAxis,true>(idx...);
                        },
                        indices...
                    );
                }
            );
        }
    }




public:
    BoundaryManager(int my_world_rank,int thread_num,TargetFunc& target_func,const BoundaryCondition& boundary_condition,const Axes&... axes):
        boundary_condition(boundary_condition),
        target_func(target_func),
        axes_tuple(std::tuple<Axes...>{axes...}),
        bc2block_id(my_world_rank,thread_num)
    {
        init_comm();
    }
    //Pack<BoundaryCondition_x_,BoundaryCondition_vr,BoundaryCondition_vt,BoundaryCondition_vp>::element<1>::left(hoge,1,2,3,4);
    
    template<typename TargetAxis>
    void apply(){
        apply_helper_l<TargetAxis>();
        apply_helper_r<TargetAxis>();
    }
};
#endif //BOUNDARY_MANAGER_H