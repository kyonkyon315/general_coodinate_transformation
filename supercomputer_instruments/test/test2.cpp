#include <mpi.h>
#include <iostream>
#include "../new_n_d_tensor_with_ghost_cell.h"
#include "../axis.h"

#include "../axis_instantiator.h"
#include "../block_id2rank.h"
// 1D
using X = Axis<0,4,2,2>;


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0)
            std::cerr << "This test requires exactly 2 MPI ranks\n";
        MPI_Finalize();
        return 0;
    }

    const auto [axis_x] = axis_instantiator<X>(rank);
    BlockId2Rank blockid_2_rank(axis_x);

    NdTensorWithGhostCell<double,X> tensor;

    tensor.set_value([&](int i){return rank*10.+i;});

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){
        std::cout << "[Rank 0] physical before exchange:\n";
        for(int i=0;i<X::num_grid;++i){
            std::cout<<tensor.at(i)<<" ";
        }
        std::cout<<"\n";
    }

    if(rank == 1){
        std::cout<< "[Rank 1] physical before exchange:\n";
        for(int i=0;i<X::num_grid;++i){
            std::cout<<tensor.at(i)<<" ";
        }
        std::cout << "\n";
    }


    if(rank == 0){
        tensor.send_ghosts<X,false,true>(1,1);
    }
    if(rank == 1){
        tensor.send_ghosts<X,true,false>(0,0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ----------------------------------
    // 3. buf_at を使って ghost を反転コピー
    // ----------------------------------
    if (rank == 0) {
        // right ghost <- reversed recv
        for (int k = 0; k < X::L_ghost_length; ++k) {
            tensor.at(X::num_grid + k) =
                tensor.buf_at<X, false>(X::num_grid - 1 - k);
        }
    }
    if (rank == 1) {
        // left ghost <- reversed recv
        for (int k = 0; k < X::L_ghost_length; ++k) {
            tensor.at(-k - 1) =
                tensor.buf_at<X, true>(k);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ----------------------------------
    // 4. 結果表示
    // ----------------------------------
    if (rank == 0) {
        std::cout << "[Rank 0] full including ghost:\n";
        for (int i = -X::L_ghost_length; i < X::num_grid + X::R_ghost_length; ++i)
            std::cout << tensor.at(i) << " ";
        std::cout << "\n";
    }
    if (rank == 1) {
        std::cout << "[Rank 1] full including ghost:\n";
        for (int i = -X::L_ghost_length; i < X::num_grid + X::R_ghost_length; ++i)
            std::cout << tensor.at(i) << " ";
        std::cout << "\n";
    }

    MPI_Finalize();

    return 0;
}


