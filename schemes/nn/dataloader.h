#ifndef DATALOADER_H
#define DATALOADER_H

#include <vector>
#include <random>
#include <utility> // for std::pair

class DataLoader {
private:
    std::mt19937 gen; // 乱数生成器
    std::uniform_real_distribution<> dis; // 0.0 から 1.0 までの一様分布

public:
    DataLoader() : gen(std::random_device{}()), dis(0.0, 1.0) {
        // コンストラクタ
    }

    /**
     * @brief バッチデータを取得する (ダミー)
     * @param batch_size バッチサイズ
     * @return { (入力バッチ), (ラベルバッチ) } のペア
     * 入力: [batch_size][2]
     * ラベル: [batch_size][1]
     */
    std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
    get_batch(int batch_size) {
        
        std::vector<std::vector<double>> inputs(batch_size, std::vector<double>(2));
        std::vector<std::vector<double>> labels(batch_size, std::vector<double>(1));

        for (int i = 0; i < batch_size; ++i) {
            double x1 = dis(gen);
            double x2 = dis(gen);
            inputs[i][0] = x1;
            inputs[i][1] = x2;

            // ダミーの学習タスク: AND ゲート
            // (0.5以上を 1, 0.5未満を 0 とする)
            bool x1_bool = (x1 > 0.5);
            bool x2_bool = (x2 > 0.5);
            labels[i][0] = (x1_bool && x2_bool) ? 1.0 : 0.0;
        }
        
        // C++17 の構造化束縛で受け取れるようペアで返す
        return {inputs, labels};
    }
};

#endif // DATALOADER_H
