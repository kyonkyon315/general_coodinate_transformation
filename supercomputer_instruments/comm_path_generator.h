#ifndef COMM_PATH_GENERATOR
#define COMM_PATH_GENERATOR

#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>
#include "node.h"
#include "block_id2rank.h"
#include "axis_instantiator.h"
#include "tuple"

template<typename BoundaryCondition, typename... Axes>
class CommPathGenerator{
private:
    BlockId2Rank<Axes...> block_id2rank;
    static constexpr int DIM = sizeof...(Axes);
    const int thread_num;
    template<int I>
    using BoundaryElement = typename BoundaryCondition::template element<I>;

    template<bool Is_tgt_left,int TargetDim>
    std::array<int, DIM> make_ghost_indices(const Axes&... axes)const{
        if constexpr(Is_tgt_left){
            return std::array<int, DIM>{
                (Axes::label == TargetDim
                    ? axes.L_id - 1   // ghost
                    //? - 1   // ghost
                    : axes.L_id      // interior
                    //: 0      // interior
                )...
            };
        }
        else{
            return std::array<int, DIM>{
                (Axes::label == TargetDim
                    ? axes.R_id      // ghost
                    //? axes.num_grid      // ghost
                    : axes.L_id      // interior
                    //: 0      // interior
                )...
            };
        }
    }

    
    template<bool Is_tgt_left,int TargetDim>
    std::tuple<int,bool> calc_candidate(const Axes&... axes)const{//src_rank, from whether left or right data come from. 
        std::tuple<Axes...> axis_tuple(axes...);
        const std::array<int, DIM> idx = make_ghost_indices<Is_tgt_left,TargetDim>(axes...);
        using TargetAxis = std::tuple_element_t<TargetDim, std::tuple<Axes...>>;

        if constexpr(Is_tgt_left){
            std::array<int,DIM> src_indices={BoundaryElement<TargetDim>::template left<Axes::label>(idx[Axes::label]...)...};
            /*
            std::cout<<"src_indices:[";
            for(int i : src_indices){
                std::cout<<i<<" ";
            }
            std::cout<<"]\n\n";*/
            return std::make_tuple(
                    /*int*/block_id2rank.get_rank((src_indices[Axes::label]/Axes::num_grid)...),
                    /*bool*/src_indices[TargetDim] % TargetAxis::num_grid <  TargetAxis::L_ghost_length//=R_ghost_length
                );
        }
        else{
            std::array<int,DIM> src_indices={BoundaryElement<TargetDim>::template right<Axes::label>(idx[Axes::label]...)...};
            return std::make_tuple(
                    /*int*/block_id2rank.get_rank((src_indices[Axes::label]/Axes::num_grid)...),
                    /*bool*/src_indices[TargetDim] % TargetAxis::num_grid <  TargetAxis::L_ghost_length//=R_ghost_length
                );
        }
    }


public:
    CommPathGenerator(int thread_num):
        thread_num(thread_num)
    {
    }
    template<typename TargetAxis>
    std::vector<Node> get_comm_path()const{
        std::vector<Node> retval(thread_num);
        constexpr int TargetDim = TargetAxis::label;

        //まずは内部の通信パスを登録する。
        //注目している軸について、block_idを±1するだけでよい。
        for(int rank=0;rank<thread_num;rank++){
            //std::cout<<rank<<" "<<std::flush;
            std::tuple<Axes...> axes = axis_instantiator<Axes...>(rank);
            auto& target_axis = std::get<TargetDim>(axes);
            if (target_axis.block_id != target_axis.num_blocks-1){
                int dst_rank = [&]<std::size_t... Is>(std::index_sequence<Is...>) -> int 
                    {
                        return block_id2rank.get_rank(
                                    (Is == TargetDim ? 
                                        std::get<Is>(axes).block_id+1 : 
                                        std::get<Is>(axes).block_id
                                    )...);
                    }(std::make_index_sequence<DIM>{});
                    //std::cout<<rank<<"->"<<dst_rank<<" ";
                retval[rank].right = Endpoint{dst_rank, Hand::LEFT};
            }
            if (target_axis.block_id != 0){
                int dst_rank = [&]<std::size_t... Is>(std::index_sequence<Is...>) -> int 
                    {
                        return block_id2rank.get_rank(
                                    (Is == TargetDim ? 
                                        std::get<Is>(axes).block_id-1 : 
                                        std::get<Is>(axes).block_id
                                    )...);
                    }(std::make_index_sequence<DIM>{});
                //std::cout<<rank<<"->"<<dst_rank<<" ";
                retval[rank].left = Endpoint{dst_rank, Hand::RIGHT};
            }
        }
        //内部の通信パスの登録　終了

        //境界条件のための通信パスの登録
        for(int rank=0;rank<thread_num;rank++){
            //std::cout<<rank<<" "<<std::flush;
            std::tuple<Axes...> axes = axis_instantiator<Axes...>(rank);
            auto& target_axis = std::get<TargetDim>(axes);
            if (target_axis.block_id == 0){
                std::tuple<int,bool> dst_rank_and_is_left = [&]<std::size_t... Is>(std::index_sequence<Is...>) -> std::tuple<int,bool> {
                        return this->template calc_candidate<true, TargetDim>((std::get<Is>(axes))...);
                    }(std::make_index_sequence<DIM>{});
                auto& [dst_rank, is_left] = dst_rank_and_is_left;
                retval[rank].left = Endpoint{dst_rank, (is_left ? Hand::LEFT : Hand::RIGHT)};
            }
            if (target_axis.block_id == TargetAxis::num_blocks -1){
                std::tuple<int,bool> dst_rank_and_is_left = [&]<std::size_t... Is>(std::index_sequence<Is...>) -> std::tuple<int,bool>
                    {
                        return calc_candidate<false, TargetDim>((std::get<Is>(axes))...);
                    }(std::make_index_sequence<DIM>{});
                auto& [dst_rank, is_left] = dst_rank_and_is_left;
                retval[rank].right = Endpoint{dst_rank, (is_left ? Hand::LEFT : Hand::RIGHT)};
            }
        }
        return retval;
    }
};

#endif //COMM_PATH_GENERATOR