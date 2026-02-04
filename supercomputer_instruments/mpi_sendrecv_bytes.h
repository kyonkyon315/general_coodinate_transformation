#ifndef MPI_SENDRECV_BYTES
#define MPI_SENDRECV_BYTES
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
    assert(send_buf.size() >= count);
    assert(recv_buf.size() >= count);
    if(dest == -1){
        std::cout<<"mpi_sendrecv_bytes rank=-1\n";
        throw 1;
    }
    
    MPI_Sendrecv(
        send_buf.data(), count * sizeof(T), MPI_BYTE, dest, send_tag,
        recv_buf.data(), count * sizeof(T), MPI_BYTE, source, recv_tag,
        comm, MPI_STATUS_IGNORE
    );
}
#endif //MPI_SENDRECV_BYTES