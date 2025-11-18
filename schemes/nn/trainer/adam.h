#ifndef ADAM_H
#define ADAM_H

#include "../model.h"      // Model クラス
#include "../dataloader.h" // DataLoader クラス
#include "../../../Timer.h"
#include <vector>
#include <cmath>      // std::sqrt, std::pow
#include <iostream>   // std::cout (ログ出力用)
#include <numeric>    // std::fill (勾配のゼロ初期化用)

template<typename FinalLayer>
class Adam {
public:
    // 参照
    Model<FinalLayer>& model;
    DataLoader& data_loader;
    
    // 訓練パラメータ
    int batch_size;
    int num_epochs;

    // Adam ハイパーパラメータ
    double learning_rate;
    double beta1;
    double beta2;
    double epsilon;

    // Adam オプティマイザの状態
    int t; // タイムステップ (バイアス補正用)
    std::vector<double> m; // 1次モーメント (移動平均)
    std::vector<double> v; // 2次モーメント (二乗勾配の移動平均)

public:
    /**
     * @brief コンストラクタ
     */
    Adam(
        Model<FinalLayer>& model_ref,
        DataLoader& loader_ref,
        int b_size,
        int n_epochs,
        double lr = 0.001,  // 学習率 (alpha)
        double b1 = 0.9,
        double b2 = 0.999,
        double eps = 1e-8
    ) : 
        model(model_ref), 
        data_loader(loader_ref), 
        batch_size(b_size), 
        num_epochs(n_epochs),
        learning_rate(lr),
        beta1(b1),
        beta2(b2),
        epsilon(eps),
        t(0), // タイムステップを 0 で初期化
        // モーメントベクトルをモデルの総パラメータ数で、0.0 で初期化
        m(model.total_param_dimension, 0.0), 
        v(model.total_param_dimension, 0.0)
    {
        // コンストラクタ本体
    }

    /**
     * @brief 学習ループを実行
     */
    void learn() {
        
        std::vector<double> output(model.output_dimension);
        std::vector<double> output_gradient(model.output_dimension);

        // 1エポックあたりのイテレーション数 (ダミー)
        // 本来は (全データ数 / バッチサイズ)
        int iterations_per_epoch = 100; 

        std::cout << "Training started..." << std::endl;
        std::cout << "Total Parameters: " << model.total_param_dimension << std::endl;
        Timer timer;
        timer.start();
        for (int epoch = 0; epoch < num_epochs; ++epoch) {
            
            double total_epoch_loss = 0.0;

            for (int iter = 0; iter < iterations_per_epoch; ++iter) {
                
                // 1. データを取得
                auto [batch_inputs, batch_labels] = data_loader.get_batch(batch_size);

                // 2. バッチの勾配をゼロ初期化
                std::fill(model.d_w_buf.begin(), model.d_w_buf.end(), 0.0);
                
                double current_batch_loss = 0.0;

                // 3. バッチ内の各データで勾配を計算・蓄積
                for (int i = 0; i < batch_size; ++i) {
                    const auto& input = batch_inputs[i];
                    const auto& true_label = batch_labels[i]; // true_label[0] に正解値

                    // 3a. 順伝播
                    model.forward(input);
                    model.get_output(output); // output[0] に予測値

                    // 3b. 損失 (MSE) と出力層の勾配を計算
                    // 損失 L = 0.5 * (output[0] - true_label[0])^2
                    // 損失の微分 dL/d_output = output[0] - true_label[0]
                    double loss = 0.;
                    for(int j=0;j<model.output_dimension;++j){
                        loss += 0.5 * std::pow(output[j] - true_label[j], 2);
                    }
                    current_batch_loss += loss;

                    std::vector<double> output_gradient(model.output_dimension);

                    for(int j=0;j<model.output_dimension;++j){
                        output_gradient[j] = output[j] - true_label[j];
                    }

                    // 3c. 逆伝播 (勾配は model.d_w_buf に蓄積される)
                    model.backward(output_gradient);
                }

                // 4. バッチの平均勾配を計算
                const double inv_batch_size = 1.0 / batch_size;
                for (double& grad : model.d_w_buf) {
                    grad *= inv_batch_size;
                }
                total_epoch_loss += (current_batch_loss * inv_batch_size);

                // 5. Adam によるパラメータ更新
                t++; // タイムステップをインクリメント

                // バイアス補正済み学習率
                double lr_t = learning_rate * std::sqrt(1.0 - std::pow(beta2, t)) 
                                           / (1.0 - std::pow(beta1, t));

                // 全パラメータを更新
                for (Index i = 0; i < model.total_param_dimension; ++i) {
                    // 勾配
                    double g = model.d_w_buf[i];

                    // モーメントの更新
                    m[i] = beta1 * m[i] + (1.0 - beta1) * g;
                    v[i] = beta2 * v[i] + (1.0 - beta2) * (g * g);

                    // パラメータ (重み) の更新
                    model.w_buf[i] -= lr_t * m[i] / (std::sqrt(v[i]) + epsilon);
                }

            } // end iteration loop

            // エポックの平均損失を出力
            if ((epoch + 1) % 100 == 0 || epoch == 0) {
                std::cout << "Epoch " << (epoch + 1) << "/" << num_epochs 
                          << ", Average Loss: " << (total_epoch_loss / iterations_per_epoch) 
                          << std::endl;
            }

        } // end epoch loop
        timer.stop();
        std::cout << "Training finished." << std::endl;
        std::cout<<timer<<"\n";
    }
};

#endif // ADAM_H