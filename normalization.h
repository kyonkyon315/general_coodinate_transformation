#ifndef NORMALIZATION_H
#define NORMALIZATION_H
namespace Norm{

struct Base{
    static constexpr double q0 = 1.602e-19; //[C]
    static constexpr double m0 = 9.109e-31; //[kg]
    static constexpr double t0 = 1.;        //[s]
    static constexpr double x0 = 1.;        //[m]
};

struct SI:Base{
    static constexpr double t0 = 1.;
    static constexpr double x0 = 1.;
};

struct Plasma : Base{
    static constexpr double omega_pe = 1.;
    static constexpr double c = 3e8;
    static constexpr double t0 = 1./ omega_pe;
    static constexpr double x0 = c/ omega_pe;
};

template<typename Units>
struct Derived{
    static constexpr double v0 = Units::x0 / Units::t0;
    static constexpr double E0 = Units::m0 * v0 /(Units::q0 * Units::t0);
    static constexpr double B0 = Units::m0 /(Units::q0 * Units::t0);

    static constexpr double n0 = 1.;
    static constexpr double J0 = n0 * Units::q0 * v0;

    static constexpr double curlB_coef = 
        (Units::x0 * Units::x0) / (Units::t0 * Units::t0);

    static constexpr double current_coef = 
        J0 * Units::t0 / E0;
};

}
#endif //NORMALIZATION_H