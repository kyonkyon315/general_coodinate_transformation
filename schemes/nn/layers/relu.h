#ifndef RELU_H
#define RELU_H
#include <vector>

using Index = int;

template<typename PreviousLayer>
class Relu{
public:
    using previous_layer_t = PreviousLayer;
    using Value = typename PreviousLayer::Value;

    static constexpr Index input_dimension 
        = PreviousLayer::output_dimension;
    static constexpr Index output_dimension
        = input_dimension;
    static constexpr Index output_id_left
        = PreviousLayer::output_id_left + PreviousLayer::output_dimension;
    
    static constexpr Index param_dimension = 0;

    static constexpr Index param_id_left
        = PreviousLayer::param_id_left + PreviousLayer::param_dimension;

    static void forward(
        std::vector<Value> &out_buf,
        const std::vector<Value> & w_buf
    ){
        for(int i=0;i<input_dimension;++i){
            out_buf[output_id_left + i]
                = (out_buf[PreviousLayer::output_id_left + i] > 0. ? 
                    out_buf[PreviousLayer::output_id_left + i]:
                    0.
                );
        }
    }

    static void backward(
        const std::vector<Value> & out_buf,
        std::vector<Value> & d_out_buf,
        const std::vector<Value> & w_buf,
        const std::vector<Value> & d_w_buf
    ){
        for(int i=0;i<input_dimension;++i){
            d_out_buf[PreviousLayer::output_id_left + i]
                = (out_buf[output_id_left + i] > 0. ? 
                    d_out_buf[output_id_left + i]:
                    0.
                );
        }
    }
};
#endif //RELU_H