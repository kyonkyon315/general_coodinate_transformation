#ifndef CONFIG_H
#define CONFIG_H
#include <cmath>
#include "axis.h"
#include "n_d_tensor.h"
#include "vec3.h"

const double Q = 1.602e-19;
const double m = 9.101e-31;

using Value = double;
using namespace std;

//計算空間の座標を設定します。
//Axis<ここには軸の通し番号をintで入力します。,ここには座標のグリッドの数をintで入力します>
//全体をなめる計算においては、通し番号が小さいものほど、より外側のループを担当することになります。
//通し番号は重複することなく、互いに隣り合った0以上の整数である必要があります。また、0を含む必要があります。
//計算空間の軸なので、一律Δx=1であり、軸同士は直交しています。
using Axis_x_ = Axis<0,1000>;
using Axis_vr = Axis<1,100>;
using Axis_vt = Axis<2,100>;
using Axis_vp = Axis<3,100>;

//電子分布関数の型を定義
//先頭に入力する型はテンソルの値の型です。その後に続く軸は、通し番号が小さいものほど左に入力してください。
using DistributionFunction = NdTensor<Value,Axis_x_,Axis_vr,Axis_vt,Axis_vp>;

//電場の型を定義
using ElectricField = NdTensor<Vec3<Value>,Axis_x_>;

//磁場の型を定義
using MagneticField = NdTensor<Vec3<Value>,Axis_x_>;

//グローバル変数としてインスタンス化しておく。
namespace Global{
    DistributionFunction dist_function;
    ElectricField e_field;
    MagneticField m_field;
}

/***********************************************
 * 物理空間と計算空間の関係を表す関数を書きます(始)*
 ***********************************************/

// --- グローバル定数とヘルパー関数の定義 ---

const Value grid_size_x_ = 3.;
const Value grid_size_vr = 4.;
const Value grid_size_vt = M_PI / (double)(Axis_vt::get_length()+1);
const Value grid_size_vp = 2.*M_PI / (double)(Axis_vp::get_length()+1);

Value vr(int calc_vr){ return grid_size_vr * (0.5 + (double)calc_vr);}
Value vt(int calc_vt){ return grid_size_vt * (0.5 + (double)calc_vt);}
Value vp(int calc_vp){ return grid_size_vp * (0.5 + (double)calc_vp);}


// --- 物理量クラス ---
//honestly_translateで計算座標↔物理座標の変換の式を定義します。
//それを用いてコンストラクタで各場所での値を事前計算してテーブルに格納します。（table.set_value(honestly_translate))
//シミュレーション中はテーブルを参照します。
//こちらも計算軸クラスと同様に通し番号を設定します。

class Physic_x_
{
public:
    Physic_x_(){}
    Value translate(int calc_x,int calc_vr,int calc_vt,int calc_vp)const{
        return grid_size_x_*calc_x;
    }
    static const int label = 0;
};

class Physic_vx
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_translate(int calc_vr,int calc_vt,int calc_vp){
        // v_x = vr * sin(vt) * cos(vp)
        return vr(calc_vr) * sin(vt(calc_vt)) * cos(vp(calc_vp));
    }
public:
    Physic_vx(){table.set_value(honestly_translate);}
    Value translate(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);    
    }
    static const int label = 1;
};

class Physic_vy
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_translate(int calc_vr,int calc_vt,int calc_vp){
        // v_y = vr * sin(vt) * sin(vp)
        return vr(calc_vr) * sin(vt(calc_vt)) * sin(vp(calc_vp));
    }
public:
    Physic_vy(){table.set_value(honestly_translate);}
    Value translate(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);    
    }
    static const int label = 2;
};

class Physic_vz
{
private:
    NdTensor<Value,Axis_vr,Axis_vt> table;
    static Value honestly_translate(int calc_vr,int calc_vt){
        // v_z = vr * cos(vt)
        return vr(calc_vr) * cos(vt(calc_vt));
    }
public:
    Physic_vz(){table.set_value(honestly_translate);}
    Value translate(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt);    
    }
    static const int label = 3;
};

//グローバル変数としてインスタンス化しておく。
namespace Global{
    Physic_x_ physic_x_;
    Physic_vx physic_vx;
    Physic_vy physic_vy;
    Physic_vz physic_vz;
}
/***********************************************
 * 計算軸を物理軸で微分した値の関数を書きます　(始)*
 ***********************************************/



class X__diff_x_
{
public:
    X__diff_x_(){}
    static constexpr Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp){
        return 1./grid_size_x_;
    }
};

class Vr_diff_vx
{
private:
    NdTensor<Value,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vt,int calc_vp){
        // d(vr) / d(vx) = sin(theta) * cos(phi)
        return sin(vt(calc_vt)) * cos(vp(calc_vp)) / grid_size_vr;
    }
public:
    Vr_diff_vx(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vt,calc_vp);
    }
};

class Vr_diff_vy
{
private:
    NdTensor<Value,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vt,int calc_vp){
        // d(vr) / d(vy) = sin(theta) * sin(phi)
        return sin(vt(calc_vt)) * sin(vp(calc_vp)) / grid_size_vr;
    }
public:
    Vr_diff_vy(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vt,calc_vp);
    }
};

class Vr_diff_vz
{
private:
    NdTensor<Value,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vt,int calc_vp){
        // d(vr) / d(vz) = cos(theta)
        return cos(vt(calc_vt)) / grid_size_vr;
    }
public:
    Vr_diff_vz(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vt,calc_vp);
    }
};

class Vt_diff_vx
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vr,int calc_vt,int calc_vp){
        // d(v_theta) / d(vx) = cos(theta) * cos(phi) / vr
        return cos(vt(calc_vt)) * cos(vp(calc_vp)) / (vr(calc_vr) * grid_size_vt);
    }
public:
    Vt_diff_vx(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);
    }
};

class Vt_diff_vy
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vr,int calc_vt,int calc_vp){
        // d(v_theta) / d(vy) = cos(theta) * sin(phi) / vr
        return cos(vt(calc_vt)) * sin(vp(calc_vp)) / (vr(calc_vr) * grid_size_vt);
    }
public:
    Vt_diff_vy(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);
    }
};

class Vt_diff_vz
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vr,int calc_vt,int calc_vp){
        // d(v_theta) / d(vz) = -sin(theta) / vr
        return -sin(vt(calc_vt)) / (vr(calc_vr) * grid_size_vt);
    }
public:
    Vt_diff_vz(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);
    }
};

class Vp_diff_vx
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vr,int calc_vt,int calc_vp){
        // d(v_phi) / d(vx) = -sin(phi) / (vr * sin(theta))
        return -sin(vp(calc_vp)) / (vr(calc_vr) * sin(vt(calc_vt)) * grid_size_vp);
    }
public:
    Vp_diff_vx(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);
    }
};

class Vp_diff_vy
{
private:
    NdTensor<Value,Axis_vr,Axis_vt,Axis_vp> table;
    static Value honestly_calc(int calc_vr,int calc_vt,int calc_vp){
        // d(v_phi) / d(vy) = cos(phi) / (vr * sin(theta))
        return cos(vp(calc_vp)) / (vr(calc_vr) * sin(vt(calc_vt)) * grid_size_vp);
    }
public:
    Vp_diff_vy(){table.set_value(honestly_calc);}
    Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp)const{
        return table(calc_vr,calc_vt,calc_vp);
    }
};

//グローバル変数としてインスタンス化
namespace Global{
    const Independent independent;
    const X__diff_x_ x__diff_x_;
    const Vr_diff_vx vr_diff_vx;
    const Vr_diff_vy vr_diff_vy;
    const Vr_diff_vz vr_diff_vz;
    const Vt_diff_vx vt_diff_vx;
    const Vt_diff_vy vt_diff_vy;
    const Vt_diff_vz vt_diff_vz;
    const Vp_diff_vx vp_diff_vx;
    const Vp_diff_vy vp_diff_vy;
}

/*******************************************************************
 * Jacobian行列を定義します。上で作成したクラスを行列風に並べてください。*
 * 互いに独立な軸の箇所（微分が０）はIndependent classを入れてください。*
 *                                                                 *
 * 具体的には、Jacobian[I,J]には「通し番号Iの計算軸」を「通し番号Jの物理*
 * 軸」で微分したものを入れてください。                               *
 *******************************************************************/
#include "jacobian.h"
Jacobian jacobian(
    Global::x__diff_x_,  Global::independent, Global::independent, Global::independent,
    Global::independent, Global::vr_diff_vx, Global::vr_diff_vy, Global::vr_diff_vz,
    Global::independent, Global::vt_diff_vx, Global::vt_diff_vy, Global::vt_diff_vz,
    Global::independent, Global::vp_diff_vx, Global::vp_diff_vy, Global::independent
);

/*******************************************
 * 解くべき移流方程式を定義します。           *
 * df/dt + v_x df/dx + q/m(E+v*B)・∇_v f = 0 *
 * を例に定義の仕方を解説                    *
 *******************************************/

//移流項の定義
//------------------------------------------
// 1. v_x * df/dx
//------------------------------------------
class Fx_ {
public:
    Fx_(){}
    Value at(int calc_x_, int calc_vr, int calc_vt, int calc_vp) const {
        return Global::physic_vx.translate(calc_x_, calc_vr, calc_vt, calc_vp);
    }
};

//------------------------------------------
// 2. q/m (E + v×B)_x
//------------------------------------------
class Fvx {
private:
public:
    Fvx(){}

    Value at(int calc_x_, int calc_vr, int calc_vt, int calc_vp) const {
        const Value vx = Global::physic_vx.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vy = Global::physic_vy.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vz = Global::physic_vz.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        return Q/m * (      Global::e_field(calc_x_).x 
                     + vy * Global::m_field(calc_x_).z
                     - vz * Global::m_field(calc_x_).y );
    }
};

//------------------------------------------
// 3. q/m (E + v×B)_y
//------------------------------------------
class Fvy {

public:
    Fvy(){}

    Value at(int calc_x_, int calc_vr, int calc_vt, int calc_vp) const {
        const Value vx = Global::physic_vx.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vy = Global::physic_vy.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vz = Global::physic_vz.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        return Q/m * (      Global::e_field(calc_x_).y 
                     + vz * Global::m_field(calc_x_).x
                     - vx * Global::m_field(calc_x_).z );
    }
};

//------------------------------------------
// 4. q/m (E + v×B)_z
//------------------------------------------
class Fvz {
public:
    Fvz(){}

    Value at(int calc_x_, int calc_vr, int calc_vt, int calc_vp) const {
        const Value vx = Global::physic_vx.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vy = Global::physic_vy.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        const Value vz = Global::physic_vz.translate(calc_x_, calc_vr, calc_vt, calc_vp);
        return Q/m * (      Global::e_field(calc_x_).z 
                     + vx * Global::m_field(calc_x_).y
                     - vy * Global::m_field(calc_x_).x );
    }
};
namespace Global{
    Fx_ flux_x_;
    Fvx flux_vx;
    Fvy flux_vy;
    Fvz flux_vz;
}
/****************************************************************************
 *物理演算子はOperatorsに、物理移流項はAdvectionsに、それぞれclass Packを用いてまとめる。
 *ただし、Operatorsの順番とAdvectionsの順番は式の順番と同じにしてください。
 *
 *さらにそれらOperatorsとAdvections、および発展させたい関数（ここではdist_func）を
 *用いてAdvectionEquationをインスタンス化します。これが、本シミュレーションのメインと
 *なるソルバーとして働きます。
 ****************************************************************************/
#include "pack.h"
#include "advection_equation.h"
namespace Global{
    Pack operators(Global::physic_x_,Global::physic_vx,Global::physic_vy,Global::physic_vz);
    Pack advections(Global::flux_x_,Global::flux_vx,Global::flux_vy,Global::flux_vz);
    AdvectionEquation equation(Global::dist_function,operators,advections,jacobian);
}

/****************************************************************************
 * 次に、数値振動を抑えるためのlimiterを選択します。今回は、Umeda2008を用います。
 ****************************************************************************/
using Limiter = Umeda2008;

#endif //CONFIG_H
