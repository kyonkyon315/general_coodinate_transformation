#ifndef AXIS_H
#define AXIS_H

template<int LABEL,int LENGTH,int NUM_BLOCKS,int LEN_GHOST=3>
class Axis {
    public:
    static constexpr int label = LABEL;
    static constexpr int num_grid = LENGTH;
    static constexpr int L_ghost_length = LEN_GHOST;
    static constexpr int R_ghost_length = REN_GHOST;
    static constexpr int num_blocks = NUM_BLOCKS;
    const int L_id;
    const int R_id;
    const int block_id;
    explicit Axis(int block_id):
        L_id(block_id*num_grid),
        R_id((block_id+1)*num_grid),
        block_id(block_id)
    {}
};

#endif //AXIS_H
