#ifndef NORMALIZATION_H
#define NORMALIZATION_H
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
    //static constexpr double omega_pe = std::sqrt(Ne*e*e/(m*epsilon0));
    static constexpr double omega_pe = 564102.5;

    //熱速度
    static constexpr double v_thermal = 0.003 * c;

    //デバイ長
    static constexpr double debye_length = v_thermal / omega_pe;

}

namespace Base{
    static constexpr double t0 = 1./Param::omega_pe;
    static constexpr double x0 = Param::debye_length;
    static constexpr double v0 = x0/t0;
    static constexpr double q0 = Param::e;
    static constexpr double m0 = Param::m;
    static constexpr double B0 = m0/(q0*t0);
    static constexpr double E0 = B0*v0;
    static constexpr double j0 = q0/(x0*x0*t0);
    static constexpr double n0 = 1./(x0*x0*x0);
    static constexpr double f0 = 1./(x0*x0*x0*v0*v0*v0);
}

namespace Coef{

    static constexpr double c_tilde = Param::c/Base::v0;
    static constexpr double Ne_tilde = Param::Ne/Base::n0;
    static constexpr double mu_tilde =  1./(c_tilde*c_tilde*Ne_tilde);
    static constexpr double epsilon_tilde =  Ne_tilde;

    // Maxwell方程式の curl B 項の係数（無次元化）
    static constexpr double maxwell_curlB_coef = c_tilde*c_tilde;

    // Maxwell方程式の 電流項の係数（無次元化）
    static constexpr double maxwell_current_coef = 1./Ne_tilde;

    // ポアソン方程式の無次元化係数
    // ∇·E = rho / epsilon0 -> ∇̃·Ẽ = poisson_coef * ne_tilde
    static constexpr double poisson_coef = 1./epsilon_tilde;
}//Coef

}//Norm
#endif //NORMALIZATION_H