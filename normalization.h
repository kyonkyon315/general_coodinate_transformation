#ifndef NORMALIZATION_H
#define NORMALIZATION_H
#include "utils/constexpr_sqrt.h"
#include <cmath>
namespace Norm{
namespace Param{
    //素電化
    static constexpr double e = 1.602e-19; //[C]

    //電子質量
    static constexpr double m = 9.109e-31; //[kg]

    //平均電子密度
    static constexpr double Ne = 1e8;

    //光速度
    static constexpr double c = 3e8;

    //誘電率
    static constexpr double epsilon0 = 8.854e-12; //[F/m]

    //透磁率
    static constexpr double mu0 = 1./(epsilon0*c*c);

    //プラズマ周波数
    static constexpr double omega_pe = Utils::ConstExpr::sqrt(Ne*e*e/(m*epsilon0));

    //サイクロトロン周波数 
    /*
    Katoh, Y.,and Y. Omura (2007), Computer simulation of chorus wave
    generation in the Earth’s inner magnetosphere, Geophys. Res.
L   ett., 34, L03102, doi:10.1029/2006GL028594.
    にならってomega_ce = omega_pe/4.
    */
    static constexpr double omega_ce = omega_pe/4.;

    //背景磁場
    static constexpr double B = omega_ce*m/e;

    //熱速度
    static constexpr double v_thermal = 0.003 * c;

    //デバイ長
    static constexpr double debye_length = v_thermal / omega_pe;

    //熱速度電子のジャイロ半径
    static constexpr double R_ce = v_thermal / omega_ce;


}

namespace Base{
    //static constexpr double t0 = 1./Param::omega_pe;
    static constexpr double t0 = 1./Param::omega_ce;
    //static constexpr double x0 = Param::debye_length;
    static constexpr double x0 = Param::R_ce;
    static constexpr double v0 = x0/t0;
    static constexpr double q0 = Param::e;
    static constexpr double m0 = Param::m;
    static constexpr double B0 = m0/(q0*t0);
    static constexpr double E0 = B0*v0;
    static constexpr double j0 = q0/(x0*x0*t0);
    static constexpr double n0 = 1./(x0*x0*x0);
    static constexpr double f0 = 1./(x0*x0*x0*v0*v0*v0);
}

namespace Coef{//無次元化された各種物理定数および物理環境定数

    static constexpr double c_tilde = Param::c/Base::v0;
    static constexpr double Ne_tilde = Param::Ne/Base::n0;

    static constexpr double mu_tilde 
        = Param::mu0 
          / ( (Base::m0 * Base::x0) / (Base::q0 * Base::q0) );

    static constexpr double epsilon_tilde
        =  Param::epsilon0
            /(Base::q0*Base::q0*Base::t0*Base::t0
                            /(Base::m0*Base::x0*Base::x0*Base::x0)
            );
    static constexpr double B_tilde = Param::B/Base::B0;

    // Maxwell方程式の curl B 項の係数（無次元化）
    static constexpr double maxwell_curlB_coef = c_tilde*c_tilde;

    // Maxwell方程式の 電流項の係数（無次元化）
    static constexpr double maxwell_current_coef = 1./epsilon_tilde;

    // ポアソン方程式の無次元化係数
    // ∇·E = rho / epsilon0 -> ∇̃·Ẽ = poisson_coef * ne_tilde
    static constexpr double poisson_coef = 1./epsilon_tilde;
}//Coef

}//Norm
#endif //NORMALIZATION_H