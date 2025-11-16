#ifndef UMEDA_2008_5_H
#define UMEDA_2008_5_H
#include <algorithm> // std::min / std::max

using Value = double;

class Umeda2008_5 {
private:
    static constexpr Value inv120 = 1.0 / 120.0;

    // calc_flux_rightward5 は中心 f0 を用いる (fm2,fm1,f0,fp1,fp2)
    // nyu >= 0 (左→右 情報) を想定した計算
    static inline Value calc_flux_rightward5(
        Value fm2, Value fm1, Value f0, Value fp1, Value fp2,
        Value nyu
    ){
        // Fortran の pic5 の流れを踏襲
        // hmax/hmin 計算（ここでは f を用いる）
        Value fmax1 = std::max(std::max(fm1, f0),
                               std::min(2.0*fm1 - fm2, 2.0*f0 - fp1));
        Value fmin1 = std::min(std::min(fm1, f0),
                               std::max(2.0*fm1 - fm2, 2.0*f0 - fp1));
        Value fmax2 = std::max(std::max(fp1, f0),
                               std::min(2.0*fp1 - fp2, 2.0*f0 - fm1));
        Value fmin2 = std::min(std::min(fp1, f0),
                               std::max(2.0*fp1 - fp2, 2.0*f0 - fm1));

        Value fmax = std::max(fmax1, fmax2);
        Value fmin = std::max(0.0, std::min(fmin1, fmin2));

        // Fortran 内の変数名に合わせる（nu 相当が nyu）
        Value nu = nyu;
        Value nup1 = nu + 1.0;
        Value nup2 = nu + 2.0;
        Value num1 = nu - 1.0;
        Value num2 = nu - 2.0;
        Value num3 = nu - 3.0;

        // ep2, ep3, ep1, ep4 を Fortran と同じロジックで決定
        Value ep3, ep2, ep1, ep4;

        if (f0 >= fm1) {
            ep3 = std::min(f0 - fm1, 2.0 * (fmax - f0));
        } else {
            ep3 = -std::min(fm1 - f0, 2.0 * (f0 - fmin));
        }

        if (fp1 >= f0) {
            ep2 = std::min(fp1 - f0, 2.0 * (f0 - fmin));
        } else {
            ep2 = -std::min(f0 - fp1, 2.0 * (fmax - f0));
        }

        if (ep2 >= ep3) {
            if (fp2 >= fp1) {
                ep1 = std::min(fp2 - fp1, 2.0 * (f0 - fmin));
            } else {
                ep1 = -std::min(fp1 - fp2, 2.0 * (f0 - fmin));
            }
            if (fm1 >= fm2) {
                ep4 = std::min(fm1 - fm2, 2.0 * (f0 - fmin));
            } else {
                ep4 = -std::min(fm2 - fm1, 2.0 * (f0 - fmin));
            }
        } else {
            if (fp2 >= fp1) {
                ep1 = std::min(fp2 - fp1, 2.0 * (fmax - f0));
            } else {
                ep1 = -std::min(fp1 - fp2, 2.0 * (fmax - f0));
            }
            if (fm1 >= fm2) {
                ep4 = std::min(fm1 - fm2, 2.0 * (fmax - f0));
            } else {
                ep4 = -std::min(fm2 - fm1, 2.0 * (fmax - f0));
            }
        }

        // Fortran の b2,b3,b1 計算を再現
        Value b2 = num1 * num2 * num3;
        Value b3 = nup2 * nup1 * num1;
        Value b1 = b2 * ( nup1 * ep1 - (3.0 * nu + 8.0) * ep2 )
                 - b3 * ( num2 * ep4 - (3.0 * nu - 11.0) * ep3 );

        // 最終的な流束（Fortran pic5 の式）
        Value U = nu * (f0 + b1 * inv120);
        return U;
    }

public:
    // Umeda のインターフェースに合わせる
    static const int used_id_left  = -3;
    static const int used_id_right =  3;

    // f_im3 .. f_ip3 の 7 点を受け取り、i に対する増分を返す
    Value calc_df(
        Value f_im3, Value f_im2, Value f_im1,
        Value f_i,
        Value f_ip1, Value f_ip2, Value f_ip3,
        Value nyu_m_half, Value nyu_p_half
    ) const {
        // 元コードに合わせて符号反転
        nyu_m_half = -nyu_m_half;
        nyu_p_half = -nyu_p_half;

        Value U_im_half, U_ip_half;

        // 左境界側フラックス U_{i-1/2}
        if (nyu_m_half >= 0.0) {
            // center を f_im1 (i-1) に合わせる配置：
            // fm2 = f_im3, fm1 = f_im2, f0 = f_im1, fp1 = f_i, fp2 = f_ip1
            U_im_half = calc_flux_rightward5(
                /*fm2*/ f_im3,
                /*fm1*/ f_im2,
                /*f0*/  f_im1,
                /*fp1*/ f_i,
                /*fp2*/ f_ip1,
                nyu_m_half
            );
        } else {
            // 負の場合は向きを反転して計算（符号も反転）
            U_im_half = - calc_flux_rightward5(
                /*fm2*/ f_ip2,
                /*fm1*/ f_ip1,
                /*f0*/  f_i,
                /*fp1*/ f_im1,
                /*fp2*/ f_im2,
                -nyu_m_half
            );
        }

        // 右境界側フラックス U_{i+1/2}
        if (nyu_p_half >= 0.0) {
            // center = f_i
            U_ip_half = calc_flux_rightward5(
                /*fm2*/ f_im2,
                /*fm1*/ f_im1,
                /*f0*/  f_i,
                /*fp1*/ f_ip1,
                /*fp2*/ f_ip2,
                nyu_p_half
            );
        } else {
            U_ip_half = - calc_flux_rightward5(
                /*fm2*/ f_ip3,
                /*fm1*/ f_ip2,
                /*f0*/  f_ip1,
                /*fp1*/ f_i,
                /*fp2*/ f_im1,
                -nyu_p_half
            );
        }

        Value delta_f_i = U_im_half - U_ip_half;
        return delta_f_i;
    }
};

#endif // UMEDA_2008_5_H
