#ifndef MPI_SENDRECV_BYTES
#define MPI_SENDRECV_BYTES
//#include <ostream>
#include <vector>
#include <mpi.h>
#include <assert.h>
template<typename T>
void mpi_sendrecv_bytes(const std::vector<T>& send_buf,
                        std::vector<T>& recv_buf,
                        int count,
                        int dest, int source,
                        int send_tag, 
                        int recv_tag,
                        MPI_Comm comm)
{
    static_assert(std::is_trivially_copyable_v<T>);
    //destination_world_rank==-1 のときは送信しない
    //source_world_rank==-1 のときは受信しない
    assert(send_buf.size() >= count);
    assert(recv_buf.size() >= count);
    /*std::cout << "rank " << "?"
          << " send_count=" << count
          << " send_size=" << send_buf.size()
          << " recv_size=" << recv_buf.size()
          << std::flush;*/
    if(dest!=-1 && source!=-1){
        MPI_Sendrecv(
            send_buf.data(), count * sizeof(T), MPI_BYTE, dest, send_tag,
            recv_buf.data(), count * sizeof(T), MPI_BYTE, source, recv_tag,
            comm, MPI_STATUS_IGNORE
        );
    }
    else if(dest != -1){
        // 送信だけ行う
        MPI_Send(
            send_buf.data(),
            count * sizeof(T),
            MPI_BYTE,
            dest,
            send_tag,
            comm
        );
    }
    else if(source != -1){
        // 受信だけ行う
        MPI_Recv(
            recv_buf.data(),
            count * sizeof(T),
            MPI_BYTE,
            source,
            recv_tag,
            comm,
            MPI_STATUS_IGNORE
        );
    }
    else{
        //なにもしない。
    }
}
#endif //MPI_SENDRECV_BYTES