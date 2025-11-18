#ifndef AFFINE_H
#define AFFINE_H
#include <vector>
using Index = int;
template<typename PreviousLayer, Index OUTPUT_DIM>
class Affine {
public:
    using previous_layer_t = PreviousLayer;
    using Value = typename PreviousLayer::Value;

    static constexpr Index input_dimension  = PreviousLayer::output_dimension;
    static constexpr Index output_dimension = OUTPUT_DIM;

    // 出力の開始位置
    static constexpr Index output_id_left =
        PreviousLayer::output_id_left + PreviousLayer::output_dimension;

    // パラメータ W, b のサイズ
    static constexpr Index weight_dimension = output_dimension * input_dimension;
    static constexpr Index bias_dimension   = output_dimension;
    static constexpr Index param_dimension  = weight_dimension + bias_dimension;

    // パラメータの開始位置
    static constexpr Index param_id_left =
        PreviousLayer::param_id_left + PreviousLayer::param_dimension;

    // -------------------
    // forward
    // -------------------
    static void forward(
        std::vector<Value> &out_buf,
        const std::vector<Value> &w_buf
    ){
        // 入力 x は PreviousLayer の出力
        const Index x0 = PreviousLayer::output_id_left;
        const Index y0 = output_id_left;
        const Index w0 = param_id_left;
        const Index b0 = param_id_left + weight_dimension;

        for(Index o = 0; o < output_dimension; ++o){
            // y[o] = b[o]
            Value sum = w_buf[b0 + o];

            // W[o][i] * x[i] の和
            const Index wrow = w0 + o * input_dimension;

            for(Index i = 0; i < input_dimension; ++i){
                sum += w_buf[wrow + i] * out_buf[x0 + i];
            }
            out_buf[y0 + o] = sum;
        }
    }

    // -------------------
    // backward
    // -------------------
    static void backward(
        const std::vector<Value> &out_buf,
        std::vector<Value> &d_out_buf,
        const std::vector<Value> &w_buf,
        std::vector<Value> &d_w_buf
    ){
        const Index x0  = PreviousLayer::output_id_left;
        const Index y0  = output_id_left;
        const Index w0  = param_id_left;
        const Index b0  = param_id_left + weight_dimension;

        // dx = W^T * dy
        for(Index i = 0; i < input_dimension; ++i){
            Value s = 0.0;
            for(Index o = 0; o < output_dimension; ++o){
                s += w_buf[w0 + o * input_dimension + i] * d_out_buf[y0 + o];
            }
            d_out_buf[x0 + i] = s;
        }

        // dW += dy[o] * x[i]
        for(Index o = 0; o < output_dimension; ++o){
            const Value g = d_out_buf[y0 + o];
            const Index wrow = w0 + o * input_dimension;

            for(Index i = 0; i < input_dimension; ++i){
                d_w_buf[wrow + i] += g * out_buf[x0 + i];
            }
        }

        // db += dy
        for(Index o = 0; o < output_dimension; ++o){
            d_w_buf[b0 + o] += d_out_buf[y0 + o];
        }
    }

    static std::string get_network_structure(){
        return "Affine<" + PreviousLayer::get_network_structure() + "," + std::to_string(output_dimension) + ">";
    }
};

#endif // AFFINE_H
