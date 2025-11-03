#ifndef PROJECTED_SAVER_2D_H
#define PROJECTED_SAVER_2D_H

#include <fstream>
#include <cstdint>
using Index = int;
template<
    typename Tensor,  // NdTensorWithGhostCell<T, AxisA, AxisB>
    typename PhysX,   // Physic_x class: translate(int,int)
    typename PhysY,   // Physic_y class: translate(int,int)
    typename AxisA,   // Axis_0
    typename AxisB    // Axis_1
>
class ProjectedSaver2D {
private:
    void validated_index(Index& idA,Index& idB){

    }
public:
    ProjectedSaver2D(Tensor& tensor,
                     const PhysX& phys_x,
                     const PhysY& phys_y,
                     AxisA, AxisB
                    )
        : tensor_(tensor), phys_x(phys_x), phys_y(phys_y)
    {}

    // save to binary file
    void save(const std::string& filename) const {
        std::ofstream fout(filename, std::ios::binary);
        if (!fout) throw std::runtime_error("Cannot open file: " + filename);

        const int64_t Nx = AxisA::num_grid;
        const int64_t Ny = AxisB::num_grid;

        // write dimensions
        fout.write((char*)&Nx, sizeof(int64_t));
        fout.write((char*)&Ny, sizeof(int64_t));

        // loop cells
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {

                int vr0 = i;
                int vr1 = i+1;
                int vt0 = j;
                int vt1 = j+1;
                double vx[4];
                double vy[4];
                // CCW 順序
                if (i+1<Nx && j+1<Ny){
                    vx[0] = phys_x.translate(vr0, vt0);
                    vy[0] = phys_y.translate(vr0, vt0);

                    vx[1] = phys_x.translate(vr1, vt0);
                    vy[1] = phys_y.translate(vr1, vt0);

                    vx[2] = phys_x.translate(vr1, vt1);
                    vy[2] = phys_y.translate(vr1, vt1);

                    vx[3] = phys_x.translate(vr0, vt1);
                    vy[3] = phys_y.translate(vr0, vt1);
                }
                else{
                    vx[0] = phys_x.honestly_translate(vr0, vt0);
                    vy[0] = phys_y.honestly_translate(vr0, vt0);

                    vx[1] = phys_x.honestly_translate(vr1, vt0);
                    vy[1] = phys_y.honestly_translate(vr1, vt0);

                    vx[2] = phys_x.honestly_translate(vr1, vt1);
                    vy[2] = phys_y.honestly_translate(vr1, vt1);

                    vx[3] = phys_x.honestly_translate(vr0, vt1);
                    vy[3] = phys_y.honestly_translate(vr0, vt1);
                }
                double fval = tensor_.at(i,j);
                //std::cout<<fval<<"\n";

                fout.write((char*)vx, sizeof(double)*4);
                fout.write((char*)vy, sizeof(double)*4);
                fout.write((char*)&fval, sizeof(double));
            }
        }

        fout.close();
    }

private:
    Tensor& tensor_;
    const PhysX& phys_x;
    const PhysY& phys_y;
};

#endif // PROJECTED_SAVER_2D_H
