#ifndef FDTD_SOLVER_1D
#define FDTD_SOLVER_1D
#include "../parameters.h"

using Value = double;
using Index = int;
//軸の方向はz
template<typename ElectricField,typename MagneticField>
class FDTD_solver_1d{
    using MagneticField_t = typename MagneticField::Element_t;
    static_assert(ElectricField::shape.size()==1,"electric field must be 1d.\n");
    static_assert(MagneticField_t::shape.size()==1,"magnetic field must be 1d.\n");
    static_assert(ElectricField::shape[0]==MagneticField_t::shape[0],"sizes of electric field and magnetic field mismatch.\n");

    static constexpr Index num_grid = ElectricField::shape[0];

    private:
    ElectricField& e_field;
    MagneticField& m_field;
    void develop_m(Value dt_per_dz){
        swap(m_field.p_half,m_field.m_half);
        for(Index i=0;i<num_grid;++i){
            m_field.p_half.at(i).x = m_field.m_half.at(i).x 
                + dt_per_dz*(e_field.at(i).y - e_field.at(i-1).y);
            m_field.p_half.at(i).y = m_field.m_half.at(i).y 
                - dt_per_dz*(e_field.at(i).x - e_field.at(i-1).x);
        }
    }

    void develop_e(Value dt_per_dz){
        dt_per_dz*=(Parameters::c2);
        for(Index i=0;i<num_grid;++i){
            e_field.at(i).x -= dt_per_dz*(m_field.p_half.at(i+1).y - m_field.p_half.at(i).y);
            e_field.at(i).y += dt_per_dz*(m_field.p_half.at(i+1).x - m_field.p_half.at(i).x);
        }
    }

    public:
    FDTD_solver_1d(ElectricField& e_field,MagneticField& m_field):
        e_field(e_field),
        m_field(m_field)
    {}

    void develop(Value dt_per_dz){
        develop_e(dt_per_dz);
        develop_m(dt_per_dz);
    }
};

#endif //FDTD_SOLVER_1D