#include <cmath>
#include "../include.h"
const double Q = 1.602e-19;
const double m = 9.101e-31;

using Value = double;
using namespace std;

//計算空間の座標を設定します。
//Axis<ここには軸の通し番号をintで入力します。,ここには座標のグリッドの数をintで入力します,3,3>
//全体をなめる計算においては、通し番号が小さいものほど、より外側のループを担当することになります。
//通し番号は重複することなく、互いに隣り合った0以上の整数である必要があります。また、0を含む必要があります。
//計算空間の軸なので、一律Δx=1であり、軸同士は直交しています。
//最後の3,3 >はゴーストセルのグリッド数です。
using Axis_vr = Axis<0,100,3,3>;
using Axis_vt = Axis<1,50,3,3>;

//電子分布関数の型を定義
//先頭に入力する型はテンソルの値の型です。その後に続く軸は、通し番号が小さいものほど左に入力してください。
using DistributionFunction = NdTensorWithGhostCell<Value,Axis_vr,Axis_vt>;

//磁場の型を定義
using MagneticField = Vec3<Value>;
//電流計算が不要の時（磁場固定のときなど）はCurrentをNone_currentにしておく
using Current = None_current;

//グローバル変数としてインスタンス化しておく。
namespace Global{
    Current current;
    DistributionFunction dist_function;
    MagneticField m_field;
}

/***********************************************
 * 物理空間と計算空間の関係を表す関数を書きます(始)*
 ***********************************************/

// --- グローバル定数とヘルパー関数の定義 ---
const Value grid_size_vr = 0.6;
const Value grid_size_vt = 2.*M_PI / (double)(Axis_vt::num_grid);

Value vr(int calc_vr){ return grid_size_vr * (0.5 + (double)calc_vr);}
Value vt(int calc_vt){ return grid_size_vt * (0.5 + (double)calc_vt);}



// --- 物理座標クラス ---
//honestly_translateで計算座標↔物理座標の変換の式を定義します。
//それを用いてコンストラクタで各場所での値を事前計算してテーブルに格納します。（table.set_value(honestly_translate))
//シミュレーション中はテーブルを参照します。
//こちらも計算軸クラスと同様に通し番号を設定します。

class Physic_vx
{
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(int calc_vr,int calc_vt){
        // v_x = vr * cos(vt)
        return vr(calc_vr) * cos(vt(calc_vt));
    }

    Physic_vx(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    Value translate(int calc_vr,int calc_vt)const{
        return table.at(calc_vr,calc_vt);    
    }
    static const int label = 0;
};

class Physic_vy
{
private:
    NdTensorWithGhostCell<Value,Axis_vr,Axis_vt> table;
public:
    static Value honestly_translate(int calc_vr,int calc_vt){
        // v_y = vr * sin(vt) 
        return vr(calc_vr) * sin(vt(calc_vt));
    }
    Physic_vy(){table.set_value_sliced<FullSlice,FullSlice>(honestly_translate);}
    Value translate(int calc_vr,int calc_vt)const{
        return table.at(calc_vr,calc_vt);    
    }
    static const int label = 1;
};


//グローバル変数としてインスタンス化しておく。
namespace Global{
    Physic_vx physic_vx;
    Physic_vy physic_vy;
}
/***********************************************
 * 計算軸を物理軸で微分した値の関数を書きます　(始)*
 ***********************************************/



class Vr_diff_vx
{
public:
    Vr_diff_vx(){}
    Value at(int calc_vr,int calc_vt)const{
        return std::cos(vt(calc_vt))/(double)grid_size_vr;
    }
};

class Vr_diff_vy
{
public:
    Vr_diff_vy(){}
    Value at(int calc_vr,int calc_vt)const{
        return std::sin(vt(calc_vt))/(double)grid_size_vr;
    }
};

class Vt_diff_vx
{
public:
    Vt_diff_vx(){}
    Value at(int calc_vr,int calc_vt)const{
        return - std::sin(vt(calc_vt))/(vr(calc_vr)*(double)grid_size_vt);
    }
};

class Vt_diff_vy
{
public:
    Vt_diff_vy(){}
    Value at(int calc_vr,int calc_vt)const{
        return  std::cos(vt(calc_vt))/(vr(calc_vr)*(double)grid_size_vt);
    }
};


//グローバル変数としてインスタンス化
namespace Global{
    const Independent independent;
    const Vr_diff_vx vr_diff_vx;
    const Vr_diff_vy vr_diff_vy;
    const Vt_diff_vx vt_diff_vx;
    const Vt_diff_vy vt_diff_vy;
}

/*******************************************************************
 * Jacobian行列を定義します。上で作成したクラスを行列風に並べてください。*
 * 互いに独立な軸の箇所（微分が０）はIndependent classを入れてください。*
 *                                                                 *
 * 具体的には、Jacobian[I,J]には「通し番号Iの計算軸」を「通し番号Jの物理*
 * 軸」で微分したものを入れてください。                               *
 *******************************************************************/

namespace Global{
Jacobian jacobian(
    vr_diff_vx , vr_diff_vy,  
    vt_diff_vx , vt_diff_vy 
);
}
/*******************************************
 * 解くべき移流方程式を定義します。           *
 * df/dt + q/m(v*B)・∇_v f = 0 *
 * を例に定義の仕方を解説                    *
 *******************************************/

//移流項の定義
//------------------------------------------
// 1. q/m (E + v×B)_x
//------------------------------------------
class Fvx {
public:
    Fvx(){}
    Value at(int calc_vr, int calc_vt) const {
        const Value vy = Global::physic_vy.translate(calc_vr, calc_vt);
        return Q/m*vy*Global::m_field.z;
    }
};

//------------------------------------------
// 2. q/m (E + v×B)_x
//------------------------------------------
class Fvy {
public:
    Fvy(){}
    Value at(int calc_vr, int calc_vt) const {
        const Value vx = Global::physic_vx.translate(calc_vr, calc_vt);
        return -Q/m*vx*Global::m_field.z;
    }
};
namespace Global{
    Fvx flux_vx;
    Fvy flux_vy;
}

/****************************************************************************
 * 次に、Flux計算機を選択します。今回は、Umeda2008を用います。
 ****************************************************************************/

using Scheme = Umeda2008;
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
 * 
 ****************************************************************************/

class BoundaryCondition_vr
{
public:
    static const int label = 0;
    template<typename Func>
    static Value left(Func func,int calc_vr, int calc_vt){
        static_assert(Axis_vt::num_grid%2 == 0,"v_theta空間のグリッド数は偶数である必要がある");
        constexpr int vt_half_num_grid=Axis_vt::num_grid/2;
        const int index_vt=(
            calc_vt < vt_half_num_grid ? 
            calc_vt+vt_half_num_grid
            :calc_vt-vt_half_num_grid
        );
        //std::cout<<func(-calc_vr-1, index_vt);
        //std::cout<<"("<<-calc_vr-1<<","<< index_vt<<") "<<std::flush;
        return func(-calc_vr-1, index_vt);
    }

    //物理的におかしいけど、とりあえずこうしておく
    template<typename Func>
    static Value right(Func func,int calc_vr, int calc_vt){
        return 0.;
        //return func(calc_vr, calc_v_-Axis_x_::num_grid);
    }
};

//thetaは周期境界条件
class BoundaryCondition_vt
{
public:
    static const int label = 1;
    template<typename Func>
    static Value left(Func func,int calc_vr, int calc_vt){
        return func(calc_vr,Axis_vt::num_grid+calc_vt);
    }
    template<typename Func>
    static Value right(Func func,int calc_vr, int calc_vt){
        return func(calc_vr,calc_vt - Axis_vt::num_grid);
    }
};

/*--------------------------------------
 * Pack を用いて境界条件をまとめます。
 *----------------------------------------------*/
namespace Global{
    BoundaryCondition_vr boundary_condition_vr;
    BoundaryCondition_vt boundary_condition_vt;
    
    Pack boundary_condition(
        boundary_condition_vr,
        boundary_condition_vt
    );
}
/*----------------------------------------------------------------------------
 * ターゲットとなる関数とboundary_conditionを用いてboundary_managerを作成します。
 *---------------------------------------------------------------------------*/


namespace Global{
    BoundaryManager boundary_manager(dist_function,boundary_condition);
}


/****************************************************************************
 * 最後に、解くべき移流方程式を定義します。
 * 
 *物理演算子はOperatorsに、物理移流項はAdvectionsに、それぞれclass Packを用いてまとめる。
 *ただし、Operatorsの順番とAdvectionsの順番は式の順番と同じにしてください。
 *
 *さらにそれらOperatorsとAdvections、および発展させたい関数（ここではdist_func）を
 *用いてAdvectionEquationをインスタンス化します。これが、本シミュレーションのメインと
 *なるソルバーとして働きます。
 ****************************************************************************/

namespace Global{
    Pack operators(physic_vx,physic_vy);
    Pack advections(flux_vx,flux_vy);
    AdvectionEquation equation(dist_function,operators,advections,jacobian,scheme,boundary_condition,current);
}
Value gauss(Value x,Value y){
    Value sigma = 25.;
    Value m_x =  0.;
    Value m_y =  0.;
    return std::exp(-((x-m_x)*(x-m_x)+(y-m_y)*(y-m_y))/sigma);
}

ProjectedSaver2D saver(
    Global::dist_function,
    Global::physic_vx, Global::physic_vy,
    Axis_vr{}, Axis_vt{}
);

int main(){
    
    Global::dist_function.set_value(
        [](int calc_vr,int calc_vt){
            Value vx = Global::physic_vx.translate(calc_vr,calc_vt);
            Value vy = Global::physic_vy.translate(calc_vr,calc_vt);
            return gauss(vx,vy);
        }
    );

    Value dt = 0.1;

    int num_steps = 100000;

    for(int i=0;i<num_steps;i++){
        if(i%(30)==0)saver.save("../data/0D2V/"+std::to_string(i/(30))+".bin");
        if(i%(30)==0)std::cout<<i<<"../data/0D2V/"+std::to_string(i/(30))+".bin"<<std::endl;
        Global::equation.solve<Axis_vr>(dt/2.);
        Global::boundary_manager.apply<Axis_vr>();
        Global::equation.solve<Axis_vt>(dt);
        Global::boundary_manager.apply<Axis_vt>();
        Global::equation.solve<Axis_vr>(dt/2.);
        Global::boundary_manager.apply<Axis_vr>();
        //saver.save("data/0D2V/"+std::to_string(i)+".bin");
    }
    return 0;
}