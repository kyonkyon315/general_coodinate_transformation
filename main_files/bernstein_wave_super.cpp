#include "mpi.h"
#include <cmath>
#include <random>

//#include "../include_super.h"
//#include "../supercomputer_instruments/axis_instantiator.h"
int rank;

using Value = double;
using namespace std;

//計算空間の座標を設定します。
//Axis<ここには軸の通し番号をintで入力します。,ここには座標のグリッドの数をintで入力します,3,3>
//全体をなめる計算においては、通し番号が小さいものほど、より外側のループを担当することになります。
//また、x,v空間で、∂x/∂v = 0である必要があります。（電流の計算を簡単に行うための措置です。）
//∂v/∂x = 0 は要求されていません。（例えば背景磁場に沿って速度空間の向きを変えたい時など）
//通し番号は重複することなく、互いに隣り合った0以上の整数である必要があります。また、0を含む必要があります。
//計算空間の軸なので、一律Δx=1であり、軸同士は直交しています。
//最後の3,3 >はゴーストセルのグリッド数です。
//物理空間↔計算空間の写像は、全単射である必要があります。
#include "../supercomputer_instruments/axis.h"
using Axis_z_ = Axis<0,256/16,16,3>;
using Axis_vr = Axis<1,512/32,32,3>;
using Axis_vt = Axis<2,128,1,3>;

//担当するブロックの各軸の左端インデックス
//main関数内で設定される。グローバル変数で使うので、ここで定義している。
int z__start_id;
int vr_start_id;
int vt_start_id;

//電子分布関数の型を定義
//先頭に入力する型はテンソルの値の型です。その後に続く軸は、通し番号が小さいものほど左に入力してください。
#include "../supercomputer_instruments/n_d_tensor_with_ghost_cell.h"
using DistributionFunction = NdTensorWithGhostCell<Value,Axis_z_,Axis_vr,Axis_vt>;

//磁場の型を定義
#include "../vec3.h"
using MagneticField = NdTensorWithGhostCell<Vec3<Value>,Axis_z_>;
//B(i,j).z=Bz(x=Δx i     ,t=Δt(j+1/2))
//B(i,j).x=Bx(x=Δx(i+1/2),t=Δt(j+1/2))
//B(i,j).y=By(x=Δx(i+1/2),t=Δt(j+1/2))

//電場の型を定義
using ElectricField = NdTensorWithGhostCell<Vec3<Value>,Axis_z_>;
//E(i,j).z=Ez(x=Δx(i+1/2),t=Δt j)
//E(i,j).x=Ex(x=Δx i     ,t=Δt j)
//E(i,j).y=Ey(x=Δx i     ,t=Δt j)

#include "../pack.h"
using VeloPack = Pack<Axis_vr,Axis_vt>;
//電流の型を定義
//電流は実空間のみのグリッドを持つので、Axis_vxは与えない。
//ただし、電流計算用の足し合わせで速度空間の情報が必要なので、Pack<Axis_vx>を与える。
#include "../supercomputer_instruments/current.h"
using Current_type = Current<Vec3<Value>,VeloPack,Axis_z_>;
//current.at(i).x = j_x(x=Δx i)
//current.at(i).y = j_y(x=Δx i)
//current.at(i).z = j_z(x=Δx(i+1/2))

//電流計算が不要の時（磁場固定のときなど）はCurrentをNone_currentにしておく
//using Current = None_current;

//グローバル変数としてインスタンス化しておく。
namespace Global{
    DistributionFunction dist_function;
    ElectricField e_field;
    MagneticField m_field;
}


/***********************************************
 * 物理空間と計算空間の関係を表す関数を書きます(始)*
 ***********************************************/
#include "../normalization.h"
// --- グローバル定数とヘルパー関数の定義 ---
constexpr Value grid_size_z_ = 0.5*3.3;
//0.3 * lambda_D

constexpr Value v_max = 5.*3.3* Norm::Param::v_thermal/Norm::Base::v0;
constexpr Value grid_size_vr = v_max / Axis_vr::num_global_grid;

constexpr Value grid_size_vt = 2.*M_PI / (double)(Axis_vt::num_global_grid);

Value vr(const int calc_vr){ return grid_size_vr * (0.5 + (double)(vr_start_id+calc_vr));}
Value vt(const int calc_vt){ return grid_size_vt * (0.5 + (double)(vt_start_id+calc_vt));}


// --- 物理量クラス ---
//honestly_translateで計算座標↔物理座標の変換の式を定義します。
//それを用いてコンストラクタで各場所での値を事前計算してテーブルに格納します。（table.set_value(honestly_translate))
//シミュレーション中はテーブルを参照します。
//こちらも計算軸クラスと同様に通し番号を設定します。

class Physic_z_
{
public:
    Physic_z_(){}
    Value honestly_translate(const int calc_z,const int calc_vr,const int calc_vt)const{
        return grid_size_z_ * (z__start_id + calc_z);
    }
    Value translate(const int calc_z,const int calc_vr,const int calc_vt)const{
        return grid_size_z_ * (z__start_id + calc_z);
    }
    static const int label = 0;
};

class Physic_vx
{
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_x = vr * cos(vt)
        return vr(calc_vr) * cos(vt(calc_vt));
    }

    Physic_vx(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    Value translate(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);    
    }
    static const int label = 1;
};

class Physic_vz
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_y = vr * sin(vt) 
        return vr(calc_vr) * sin(vt(calc_vt));
    }
    Physic_vz(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    Value translate(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);    
    }
    static const int label = 2;
};

using Operators = Pack<Physic_z_,Physic_vx,Physic_vz>;


//グローバル変数としてインスタンス化しておく。
namespace Global{
    Physic_z_ physic_z_;
    Physic_vx physic_vx;
    Physic_vz physic_vz;
    Operators operators(physic_z_,physic_vx,physic_vz);
}
/***********************************************
 * 計算軸を物理軸で微分した値の関数を書きます　(始)*
 ***********************************************/



class Z__diff_z_
{
public:
    Z__diff_z_(){}
    Value at(const int calc_z,const int calc_vr,const int calc_vt)const{
        return 1./grid_size_z_;
    }
};

class Vr_diff_vx
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_y = vr * sin(vt) 
        return std::cos(vt(calc_vt))/(double)grid_size_vr;
    }
    Vr_diff_vx(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    
    Value at(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);
    }
};

class Vr_diff_vz
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_y = vr * sin(vt) 
        return std::sin(vt(calc_vt))/(double)grid_size_vr;
    }
    Vr_diff_vz(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    
    Value at(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);
    }
};

class Vt_diff_vx
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_y = vr * sin(vt) 
        return  - std::sin(vt(calc_vt))/(vr(calc_vr)*(double)grid_size_vt);
    }
    Vt_diff_vx(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
   
    Value at(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);
    }
};

class Vt_diff_vz
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(const int calc_vr,const int calc_vt){
        // v_y = vr * sin(vt) 
        return  std::cos(vt(calc_vt))/(vr(calc_vr)*(double)grid_size_vt);
    }
    Vt_diff_vz(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
  
    Value at(const int calc_z,const int calc_vr,const int calc_vt)const{
        return table.at(calc_vr,calc_vt);
    }
};

#include "../independent.h"
//グローバル変数としてインスタンス化
namespace Global{
    const Independent independent;
    const Z__diff_z_ z__diff_z_;
    const Vr_diff_vx vr_diff_vx;
    const Vr_diff_vz vr_diff_vz;
    const Vt_diff_vx vt_diff_vx;
    const Vt_diff_vz vt_diff_vz;
}

/*******************************************************************
 * Jacobian行列を定義します。上で作成したクラスを行列風に並べてください。*
 * 互いに独立な軸の箇所（微分が０）はIndependent classを入れてください。*
 *                                                                 *
 * 具体的には、Jacobian[I,J]には「通し番号Iの計算軸」を「通し番号Jの物理*
 * 軸」で微分したものを入れてください。                               *
 *******************************************************************/
#include "../jacobian.h"
namespace Global{
Jacobian jacobian(
    z__diff_z_ , independent, independent, 
    independent, vr_diff_vx , vr_diff_vz, 
    independent, vt_diff_vx , vt_diff_vz 
);
}
Value jacobi_det(const int calc_z,const int calc_vr,const int calc_vt){
    return grid_size_z_*grid_size_vr*grid_size_vt*vr(calc_vr);
}
/**********************************************
 * 解くべき移流方程式を定義します。              *
 * df/dt + v_x df/dx + q/m(E+v*B)・∇_v f = 0 *
 * を例に定義の仕方を解説                       *
 **********************************************/

//移流項の定義
//------------------------------------------
// 1. v_x * df/dx
//------------------------------------------
class Fz_ {
public:
    Fz_(){}
    Value at(const int calc_z,const int calc_vr,const int calc_vt) const {
        return Global::physic_vz.translate(calc_z, calc_vr, calc_vt);
    }
};

//B(i,j).z=Bz(x=Δx i     ,t=Δt(j+1/2))
//B(i,j).x=Bx(x=Δx(i+1/2),t=Δt(j+1/2))
//B(i,j).y=By(x=Δx(i+1/2),t=Δt(j+1/2))

//E(i,j).z=Ez(x=Δx(i+1/2),t=Δt j)
//E(i,j).x=Ex(x=Δx i     ,t=Δt j)
//E(i,j).y=Ey(x=Δx i     ,t=Δt j)


//------------------------------------------
// 2. q/m (E + v×B)_x
//------------------------------------------
namespace Global{
    bool is_velo_left_edge;
    bool is_velo_right_edge;
}
class Fvx {
private:
public:
    Fvx(){}

    Value at(const int calc_z,const int calc_vr,const int calc_vt) const {
        //速度空間の境界でフラックスが0になるように、移流を反対称にする。
        if(Global::is_velo_right_edge && calc_vr == Axis_vr::num_grid){
            return - at(calc_z,  Axis_vr::num_grid-1,calc_vt);
        }
        else{
            //Yee格子を採用しているため、電場はｘ方向に、磁場はｔ方向に補間しなければならない。
            const Value Ex = Global::e_field.at(calc_z).x;
            const Value By = (  Global::m_field.at(calc_z-1).y +
                                Global::m_field.at(calc_z).y)/2.;
            const Value vz = Global::physic_vz.translate(calc_z, calc_vr, calc_vt);

            return - (Ex - vz*By);//電子の電荷が負なので - がつく
            //return Parameters::Q/Parameters::m * Ex;
            //規格化したので移流項はExのみ
        }
    }
};

class Fvz {
private:
public:
    Fvz(){}

    Value at(const int calc_z,const int calc_vr,const int calc_vt) const {
        //速度空間の境界でフラックスが0になるように、移流を反対称にする。
        if(Global::is_velo_right_edge && calc_vr == Axis_vr::num_grid){
            return - at(calc_z,  Axis_vr::num_grid-1,calc_vt);
        }
        else{
            //Yee格子を採用しているため、電場はｘ方向に、磁場はｔ方向に補間しなければならない。
            const Value Ez = (  Global::e_field.at(calc_z-1).y +
                                Global::e_field.at(calc_z).y)/2.;
            const Value By = (  Global::m_field.at(calc_z-1).y +
                                Global::m_field.at(calc_z).y)/2.;
            const Value vx = Global::physic_vx.translate(calc_z, calc_vr, calc_vt);
            
            return - (Ez + vx*By);//電子の電荷が負なので - がつく
            //return Parameters::Q/Parameters::m * Ex;
            //規格化したので移流項はExのみ
        }
    }
};

namespace Global{
    Fz_ flux_z_;
    Fvx flux_vx;
    Fvz flux_vz;
}

/****************************************************************************
 * 次に、Flux計算機を選択します。今回は、Umeda2008を用います。
 ****************************************************************************/
#include "../schemes/umeda_2008_fifth_order.h"
#include "../schemes/umeda_2008.h"
//using Scheme = Umeda2008;
using Scheme = Umeda2008_5;
namespace Global{
    Scheme scheme;
}
/*****************************************************************************
 * 境界条件を設定します。
 * ここで、境界条件を設定することと、ゴーストセルの更新方法を設定することは同値です。
 * ユーザーの目的を満たすようなゴーストセルの更新方法を設定してください。
 * 
 * 例えば、x方向に周期境界条件を用いたいならば、ゴーストセルは次のように更新すべきだ
 * ということは自明でしょう。
 * f(-x) ← f(x_length-x)
 * f(x_length+x) ← f(x)
 * 
 * また、反射境界を設定する場合は次のようになります。
 * f(-x,v) ← f(x,-v)
 * f(x_length+x,v) ← f(x_length-x,-v)
 * 
 * Axes と同様、labelを付けます。
 * ここではグローバルインデックスで記述してください。
 * 
 ****************************************************************************/
//xは周期境界条件
class BoundaryCondition_z_
{
public:
    static const int label = 0;

    //吸収項などを実装するとき、ゴーストセルにほかのセルの値を代入するだけではなくなる。このときtrueにする。->その場合の動作は未定義
    static constexpr bool not_only_comm = false;

    //ghost_cell[calc_x,calc_vx] <- cell[calc_x_ + Axis_x_::num_grid,calc_vx] 
    template<int Index>
    static int left(const int calc_z,const int calc_vr,const int calc_vt){
        if constexpr(Index == 0){
            return calc_z + Axis_z_::num_global_grid;
        }
        else if constexpr(Index == 1){
            return calc_vr;
        }
        else if constexpr(Index == 2){
            return calc_vt;
        }
        else return 0;
    }
    //ghost_cell[calc_x,calc_vx] <- cell[calc_x_ - Axis_x_::num_grid,calc_vx] 
    template<int Index>
    static int right(const int calc_z,const int calc_vr,const int calc_vt){
        if constexpr(Index == 0){
            return calc_z - Axis_z_::num_global_grid;
        }
        else if constexpr(Index == 1){
            return calc_vr;
        }
        else if constexpr(Index == 2){
            return calc_vt;
        }
        else return 0;
    }
};
class BoundaryCondition_vr
{
public:
    static const int label = 1;
    static constexpr bool not_only_comm = false;

    template<int Index>
    static int left(const int calc_z,const int calc_vr,const int calc_vt){
        static_assert(Axis_vt::num_grid%2 == 0,"v_theta空間のグリッド数は偶数である必要がある");
        constexpr int vt_half_num_grid = Axis_vt::num_global_grid/2;
        if constexpr(Index == 0){
            return calc_z;
        }
        else if constexpr(Index == 1){
            return -calc_vr-1;
        }
        else if constexpr(Index == 2){
            const int index_vt=(
                calc_vt < vt_half_num_grid ? 
                calc_vt+vt_half_num_grid:
                calc_vt-vt_half_num_grid
            );
            return index_vt;
        }
        else return 0;
    }

    template<int Index>
    static int right(const int calc_z,const int calc_vr,const int calc_vt){
        if constexpr(Index == 0){
            return calc_z;
        }
        else if constexpr(Index == 1){
            return calc_vr;
        }
        else if constexpr(Index == 2){
            return Axis_vr::num_global_grid - 1 - (calc_vr - Axis_vr::num_global_grid);
        }
        else return 0;
    }
};

//theta は周期境界条件
class BoundaryCondition_vt
{
public:
    static const int label = 2;
    static constexpr bool not_only_comm = false;

    template<int Index>
    static int left(const int calc_z,const int calc_vr,const int calc_vt){
        if constexpr(Index == 0){
            return calc_z;
        }
        else if constexpr(Index == 1){
            return calc_vr;
        }
        else if constexpr(Index == 2){
            return calc_vt + Axis_vt::num_global_grid;
        }
        else return 0;
    }

    template<int Index>
    static int right(const int calc_z,const int calc_vr,const int calc_vt){
        if constexpr(Index == 0){
            return calc_z;
        }
        else if constexpr(Index == 1){
            return calc_vr;
        }
        else if constexpr(Index == 2){
            return calc_vt - Axis_vt::num_global_grid;
        }
        else return 0;
    }
};

/*--------------------------------------
 * Pack を用いて境界条件をまとめます。
 *----------------------------------------------*/
using BoundaryCondition = Pack<BoundaryCondition_z_, BoundaryCondition_vr, BoundaryCondition_vt>;
namespace Global{
    BoundaryCondition_z_ boundary_condition_z_;
    BoundaryCondition_vr boundary_condition_vr;
    BoundaryCondition_vt boundary_condition_vt;

    Pack boundary_condition(
        boundary_condition_z_,
        boundary_condition_vr,
        boundary_condition_vt
    );
}

//電場、磁場についても境界条件をまとめる
class BoundaryCondition_M_z_
{
public:
    static const int label = 0;

    //吸収項などを実装するとき、ゴーストセルにほかのセルの値を代入するだけではなくなる。このときtrueにする。->その場合の動作は未定義
    static constexpr bool not_only_comm = false;

    //ghost_cell[calc_x,calc_vx] <- cell[calc_x_ + Axis_x_::num_grid,calc_vx] 
    template<int Index>
    static int left(const int calc_z){
        if constexpr(Index == 0){
            return calc_z + Axis_z_::num_global_grid;
        }
        else return 0;
    }
    //ghost_cell[calc_x,calc_vx] <- cell[calc_x_ - Axis_x_::num_grid,calc_vx] 
    template<int Index>
    static int right(const int calc_z){
        if constexpr(Index == 0){
            return calc_z - Axis_z_::num_global_grid;
        }
        else return 0;
    }
};

namespace Global{
    BoundaryCondition_M_z_ boundary_condition_M_z_;

    Pack boundary_condition_em(boundary_condition_M_z_);
}

/*----------------------------------------------------------------------------
 * ターゲットとなる関数とboundary_conditionを用いてboundary_managerを作成します。
 *---------------------------------------------------------------------------*/



/****************************************************************************
 * 最後に、解くべき移流方程式を定義します。
 *
 *Advections、および発展させたい関数（ここではdist_func）を
 *用いてAdvectionEquationをインスタンス化します。これが、本シミュレーションにおける
 *ブラソフソルバーとして働きます。
 ****************************************************************************/

/*
fdtd 関連
*/


template<typename Field>
void apply_periodic_1d(Field& f)
{
    constexpr int N  = Axis_z_::num_grid;
    constexpr int GL = Axis_z_::L_ghost_length;
    constexpr int GR = Axis_z_::R_ghost_length;

    for(int g = 1; g <= GL; ++g){
        f.at(-g) = f.at(N - g);
    }
    for(int g = 0; g < GR; ++g){
        f.at(N + g) = f.at(g);
    }
}
/*
fdtd関連の設定終わり。
*/

/*分布関数の初期化関数の設定*/
Value fM(const Value vr_tilde/*無次元量が入る*/){
    return Norm::Coef::Ne_tilde * std::exp(-vr_tilde * vr_tilde /2.)
           /(M_PI * std::sqrt(2*M_PI));
    //Ne_tilde = int f_tilde dv_tilde^3
}


void initialize_distribution(int seed)
{
    constexpr Value eps = 1e-3;

    std::mt19937 rng(12345 + seed);
    std::uniform_real_distribution<Value> uni(-1.0,1.0);

    for(int iz=0; iz<Axis_z_::num_grid; iz++){
        Value eta = uni(rng);  // x 依存ノイズ
        //Value eta = std::sin(30.*2.*M_PI*(Value)ix/(Value)Axis_x_::num_grid);
        Value base = 1.;
        //if(ix>Axis_x_::num_grid/4 && ix<3*Axis_x_::num_grid/4)base = 0.01;

        for(int ivr=0; ivr<Axis_vr::num_grid; ivr++){
            for(int ivt=0; ivt<Axis_vt::num_grid; ivt++){
                Value vr_tilde = vr(ivr);
            Global::dist_function.at(iz,ivr,ivt)
                = jacobi_det(iz,ivr,ivt)
                    *fM(vr_tilde) * (base + eps * eta);
                //やこびあんで計算空間にスケールする
            }
        }
    }
}

void solve_poisson_1d_periodic() {
    int N = Axis_z_::num_grid;
    
    // イオン密度を計算（電子密度の平均値
    Value n_e_tilde_avg = 0.0;
    for(int i=0;i<Axis_z_::num_grid;++i){
        for(int j=0;j<Axis_vr::num_grid;++j){
            for(int k=0;k<Axis_vt::num_grid;++k){
                n_e_tilde_avg += Global::dist_function.at(i,j,k);
            }
        }
    }
    n_e_tilde_avg /= (Value)N;

    // 電場の積分
    Global::e_field.at(0).z = 0.0; // 基準値（ポテンシャル自由度）
    for(int i=0;i<Axis_z_::num_grid;i++){
        Value n_e_tilde=0.;
        for(int j=0;j<Axis_vr::num_grid;++j){
            for(int k=0;k<Axis_vt::num_grid;++k){
                n_e_tilde += Global::dist_function.at(i,j,k);
            }
        }
        Global::e_field.at(i+1).z 
            = Global::e_field.at(i).z 
            + Norm::Coef::poisson_coef * (n_e_tilde_avg - n_e_tilde)/* *gridsize(=1)*/;

    }
    
    // 平均を引いて、周期境界条件を調整
    double E_mean = 0.0;
    for(int i=0;i<N;i++) E_mean += Global::e_field.at(i).z;
    E_mean /= (N);

    for(int i=0;i<N;i++) Global::e_field.at(i).z -= E_mean;

    apply_periodic_1d(Global::e_field);
}

#include "../supercomputer_instruments/advection_equation.h"
#include "../supercomputer_instruments/FDTD/fdtd_solver_1d.h"
#include "../supercomputer_instruments/axis_instantiator.h"
#include "../supercomputer_instruments/boundary_manager.h"
#include "../utils/Timer.h"
#include "../calc_current_in_x.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    Global::dist_function.comm_init(world_rank);
    Global::e_field.comm_init(world_rank);
    Global::m_field.comm_init(world_rank);
    
    Current_type current(world_rank);
    
    Pack advections(Global::flux_z_,Global::flux_vx,Global::flux_vz);
    AdvectionEquation equation(Global::dist_function,advections,Global::jacobian,Global::scheme, current);
    
    FDTD_solver_1d fdtd_solver(Global::e_field,Global::m_field,current);
    CalcCurrent_1d
    current_calculator(
        current,
        Global::dist_function,
        Global::operators,
        grid_size_z_);

    //担当するブロックの各軸の左端インデックスをrankから計算
    auto [axis_z_, axis_vr, axis_vt] = axis_instantiator<Axis_z_,Axis_vr,Axis_vt>(world_rank);
    z__start_id = axis_z_.L_id;
    vr_start_id = axis_vr.L_id;
    vt_start_id = axis_vt.L_id;
    Global::is_velo_left_edge = (axis_vr.block_id == 0);
    Global::is_velo_right_edge = (axis_vr.block_id == Axis_vr::num_blocks-1);

    BoundaryManager boundary_manager(world_rank,world_size,Global::dist_function,Global::boundary_condition,axis_z_,axis_vr,axis_vt);
    BoundaryManager boundary_manager_e(world_rank,world_size,Global::e_field,Global::boundary_condition_em,axis_z_,axis_vr, axis_vt);
    BoundaryManager boundary_manager_m(world_rank,world_size,Global::m_field,Global::boundary_condition_em,axis_z_,axis_vr, axis_vt);
    
    //初期化
    initialize_distribution(axis_z_.block_id);
    //背景磁場の設定
    for(int i=0;i<Axis_z_::num_grid;i++){
        Global::m_field.at(i).y = Norm::Coef::B_tilde;
    }
    // 初期化後のゴーストセル更新（重要）
    boundary_manager.apply<Axis_z_>();
    boundary_manager.apply<Axis_vr>();
    boundary_manager.apply<Axis_vt>();
    boundary_manager_e.apply<Axis_z_>();
    boundary_manager_m.apply<Axis_z_>();

    solve_poisson_1d_periodic();

    Value dt = 0.01 ;

    int num_steps = 100000;
    std::ofstream ex_log("Ez_t_blockid_z_" + std::to_string(axis_z_.block_id) + ".dat");
    Timer timer;
    timer.start();
    for(int i=0;i<num_steps;i++){
        //if(world_rank==0 && i%100==0)std::cout<<i<<std::endl;
        //v(0), x(0), E(0), B(1/2), J(-1/2)
        
        equation.solve<Axis_vr>(dt/2.);
        boundary_manager.apply<Axis_vr>();

        equation.solve<Axis_vt>(dt/2.);
        boundary_manager.apply<Axis_vt>();

        //v(1/2), x(0), E(0), B(1/2), J(-1/2)

        equation.solve<Axis_z_>(dt);
        boundary_manager.apply<Axis_z_>();
        
        current.clear();
        current_calculator.calc();
        current.compute_global_current();

        //v(1/2), x(1), E(0), B(1/2), J(1/2)

        fdtd_solver.develop_e(dt , grid_size_z_);
        boundary_manager_e.apply<Axis_z_>();

        //v(1/2), x(1), E(1), B(1/2), J(1/2)
        

        equation.solve<Axis_vt>(dt/2.);
        boundary_manager.apply<Axis_vt>();

        equation.solve<Axis_vr>(dt/2.);
        boundary_manager.apply<Axis_vr>();
        //v(1), x(1), E(1), B(1/2), J(1/2)
        
        fdtd_solver.develop_m(dt , grid_size_z_);
        boundary_manager_m.apply<Axis_z_>();
        //v(1), x(1), E(1), B(3/2), J(1/2)


        //if(i%20 == 0){
        if(false){
            Global::dist_function.save_physical_fast("../output/two_stream/rank_" 
                                    + std::to_string(world_rank) 
                                    + "__step_" 
                                    + std::to_string(i/20) 
                                    + ".bin");
        }

        if(axis_vr.block_id==0 && axis_vt.block_id==0){
        //if(false){

            for(int ix=0; ix<Axis_z_::num_grid; ix++){
                ex_log << Global::e_field.at(ix).z << " ";
            }
            ex_log << "\n";
        }
    }
    timer.stop();
    if(world_rank==0)std::cout<<timer<<"\n";
    MPI_Finalize();
    return 0;
}