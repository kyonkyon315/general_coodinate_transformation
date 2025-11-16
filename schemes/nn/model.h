#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <stdexcept> // for std::runtime_error
#include <random>    // <--- [追加] 乱数生成のため
#include <cmath>     // <--- [追加] std::sqrt のため

#include "layers/input_layer.h"
#include "layers/affine.h"
#include "layers/relu.h"

// 各ファイルで共通の型エイリアス
using Index = int;

// -----------------------------------------------------------------
// 内部ヘルパー (テンプレート再帰用)
// -----------------------------------------------------------------
namespace internal {

// --- 順伝播ヘルパー ---
// 実行順: Input -> Layer1 -> ... -> FinalLayer
// (再帰の「戻り」で実行する)
template<typename Layer>
struct ForwardPass {
    static void run(
        std::vector<typename Layer::Value>& out_buf,
        const std::vector<typename Layer::Value>& w_buf
    ) {
        // 1. 先に前のレイヤーを実行
        ForwardPass<typename Layer::previous_layer_t>::run(out_buf, w_buf);
        // 2. 次に自分のレイヤーを実行
        Layer::forward(out_buf, w_buf);
    }
};

// InputLayer のための基底ケース (何もしない)
template<Index DIM>
struct ForwardPass<InputLayer<DIM>> {
    static void run(
        std::vector<double>& out_buf,
        const std::vector<double>& w_buf
    ) {
        // InputLayer には forward() がないので何もしない
        // データは Model::forward() でコピー済み
    }
};

// --- 逆伝播ヘルパー ---
// 実行順: FinalLayer -> ... -> Layer1 -> Input
// (再帰の「行き」で実行する)
template<typename Layer>
struct BackwardPass {
    static void run(
        const std::vector<typename Layer::Value>& out_buf,
        std::vector<typename Layer::Value>& d_out_buf,
        const std::vector<typename Layer::Value>& w_buf,
        std::vector<typename Layer::Value>& d_w_buf
    ) {
        // 1. 最初に自分のレイヤーを実行
        Layer::backward(out_buf, d_out_buf, w_buf, d_w_buf);
        // 2. 次に前のレイヤーを実行
        BackwardPass<typename Layer::previous_layer_t>::run(out_buf, d_out_buf, w_buf, d_w_buf);
    }
};

// InputLayer のための基底ケース (何もしない)
template<Index DIM>
struct BackwardPass<InputLayer<DIM>> {
    static void run(
        const std::vector<double>& out_buf,
        std::vector<double>& d_out_buf,
        const std::vector<double>& w_buf,
        std::vector<double>& d_w_buf
    ) {
        // InputLayer には backward() がないので何もしない
    }
};


// --- 入力層の次元を取得するヘルパー ---
template<typename Layer>
struct GetInputDim {
    static constexpr Index value = GetInputDim<typename Layer::previous_layer_t>::value;
};

// InputLayer のための基底ケース
template<Index DIM>
struct GetInputDim<InputLayer<DIM>> {
    static constexpr Index value = InputLayer<DIM>::input_dimension;
};

// <--- [ここから追加] ---

// --- 重み初期化ヘルパー ---
// デフォルト (Relu層など、パラメータを持たない層用)
template<typename Layer>
struct InitializeWeights {
    static void run(std::vector<typename Layer::Value>& w_buf, std::mt19937& gen) {
        // 前の層の初期化を呼び出す
        InitializeWeights<typename Layer::previous_layer_t>::run(w_buf, gen);
        // この層 (Reluなど) はパラメータがないので何もしない
    }
};

// Affine レイヤーのためのテンプレート特殊化
template<typename PreviousLayer, Index OUTPUT_DIM>
struct InitializeWeights<Affine<PreviousLayer, OUTPUT_DIM>> {
    
    using Layer = Affine<PreviousLayer, OUTPUT_DIM>;
    using Value = typename Layer::Value;

    static void run(std::vector<Value>& w_buf, std::mt19937& gen) {
        // 1. まず再帰的に前の層の初期化を実行
        InitializeWeights<typename Layer::previous_layer_t>::run(w_buf, gen);

        // 2. この Affine 層の W (重み) と b (バイアス) を初期化
        const Index w0 = Layer::param_id_left;
        const Index b0 = Layer::param_id_left + Layer::weight_dimension;
        const Index n_in = Layer::input_dimension;
        
        // He 初期化 (Relu を想定)
        // 平均 0, 標準偏差 sqrt(2.0 / n_in) の正規分布
        double stddev = std::sqrt(2.0 / static_cast<double>(n_in));
        std::normal_distribution<Value> d_w(0.0, stddev);

        // W (重み) を He 初期化
        for (Index i = 0; i < Layer::weight_dimension; ++i) {
            w_buf[w0 + i] = d_w(gen);
        }

        // b (バイアス) を 0.0 で初期化
        for (Index i = 0; i < Layer::bias_dimension; ++i) {
            w_buf[b0 + i] = 0.0;
        }
    }
};

// InputLayer のための基底ケース (再帰の終点)
template<Index DIM>
struct InitializeWeights<InputLayer<DIM>> {
    static void run(std::vector<double>& w_buf, std::mt19937& gen) {
        // 何もしない
    }
};

// <--- [追加ここまで] --

} // namespace internal


// -----------------------------------------------------------------
// Model クラス
// -----------------------------------------------------------------
template<typename FinalLayer>
class Model {
public:
    using Value = typename FinalLayer::Value;

    // --- ネットワーク全体の次元 ---
    
    // 入力次元 (再帰ヘルパーで取得)
    static constexpr Index input_dimension = internal::GetInputDim<FinalLayer>::value;
    
    // 出力次元 (最終レイヤーから取得)
    static constexpr Index output_dimension = FinalLayer::output_dimension;

    // 全パラメータ数 (最終レイヤーから取得)
    static constexpr Index total_param_dimension = 
        FinalLayer::param_id_left + FinalLayer::param_dimension;
    
    // 全中間出力バッファサイズ (最終レイヤーから取得)
    static constexpr Index total_output_buffer_dimension =
        FinalLayer::output_id_left + FinalLayer::output_dimension;

    // --- バッファ ---
    std::vector<Value> w_buf;     // 重み (Weights)
    std::vector<Value> d_w_buf;   // 重みの勾配 (Delta Weights)
    std::vector<Value> out_buf;   // 中間出力 (Activations)
    std::vector<Value> d_out_buf; // 中間出力の勾配 (Delta Activations)

public:
    // --- コンストラクタ ---
    Model() :
        w_buf(total_param_dimension),
        d_w_buf(total_param_dimension), // 勾配は 0 で初期化
        out_buf(total_output_buffer_dimension),
        d_out_buf(total_output_buffer_dimension)
    {
        // 重みのランダム初期化を実行
        std::mt19937 gen(std::random_device{}());
        internal::InitializeWeights<FinalLayer>::run(w_buf, gen);
    }

    // --- 順伝播 ---
    void forward(const std::vector<Value>& input) {
        if (input.size() != input_dimension) {
            throw std::runtime_error("Invalid input dimension");
        }

        // 1. 入力を out_buf の先頭にコピー
        // (InputLayer::output_id_left は 0 のはず)
        const Index x0 = 0; // InputLayer::output_id_left
        for (Index i = 0; i < input_dimension; ++i) {
            out_buf[x0 + i] = input[i];
        }

        // 2. 再帰的に順伝播を実行
        internal::ForwardPass<FinalLayer>::run(out_buf, w_buf);
    }

    // --- 逆伝播 ---
    void backward(const std::vector<Value>& output_gradient) {
        if (output_gradient.size() != output_dimension) {
            throw std::runtime_error("Invalid output gradient dimension");
        }

        // 1. 出力層の勾配を d_out_buf の該当箇所にコピー
        const Index y0 = FinalLayer::output_id_left;
        for (Index i = 0; i < output_dimension; ++i) {
            d_out_buf[y0 + i] = output_gradient[i];
        }

        // 2. 再帰的に逆伝播を実行
        internal::BackwardPass<FinalLayer>::run(out_buf, d_out_buf, w_buf, d_w_buf);
    }

    // --- ユーティリティ (結果の取得) ---

    // 最終出力を取得
    void get_output(std::vector<Value>& output) const {
        if (output.size() != output_dimension) {
            std::cout<<"output_dimension = "<<output_dimension<<"\n";
            std::cout<<"output.size()= "<<output.size()<<"\n";
            throw std::runtime_error("Invalid output dimension");
        }
        const Index y0 = FinalLayer::output_id_left;
        std::copy(out_buf.begin() + y0, 
                  out_buf.begin() + y0 + output_dimension, 
                  output.begin());
    }

    // 入力層の勾配を取得 (d_x)
    void get_input_gradient(std::vector<Value> &input_grad) const {
        if (input_grad.size() != input_dimension) {
            throw std::runtime_error("Invalid input dimension");
        }
        const Index x0 = 0; // InputLayer::input_id_left
        std::copy(d_out_buf.begin() + x0, 
                  d_out_buf.begin() + x0 + input_dimension, 
                  input_grad.begin());
    }
};

#endif // MODEL_H