#ifndef SLICE_H
#define SLICE_H

using Index = int;

// --- スライス用ヘルパー型 ---

/**
 * @brief スライス指定: この次元の全物理領域 (:) を対象とすることを示すタグ
 */
struct FullSlice {};

/**
 * @brief スライス指定: [START, END) の半開区間を対象とすることを示すタグ
 * (Pythonの [START:END] と同じ)
 */
template<Index START, Index END>
struct Slice {
    static_assert(START <= END, "Slice: START must be less than or equal to END");
    static constexpr Index START_val = START;
    static constexpr Index END_val = END;
};

#endif //SLICE_H