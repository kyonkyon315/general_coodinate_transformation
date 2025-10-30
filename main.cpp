#include <iostream>
#include <vector>

double PI=3.14159265;
#include <array>
#include <iostream>
#include <functional>
#include <omp.h>
#include <cmath>
#include "Timer.h"
#include "axis.h"
#include "n_d_tensor.h"

template<typename T>
class Vec3{
    T x,y,z;
};



int main(){
    using Value = double;

    //計算空間の座標を設定します。
    //Axis<ここには座標のグリッドの数をintで入力します>
    //計算空間の軸なので、一律Δx=1であり、軸同士は直交しています。
    using Axis_vr = Axis<100>;
    using Axis_vt = Axis<100>;
    using Axis_vp = Axis<100>;
    using Axis_x_ = Axis<1000>;

    //電子分布関数作成
    NdTensor<Value,Axis_x_,Axis_vr,Axis_vt,Axis_vp> f;
    
    //電場作成
    NdTensor<Vec3<Value>,Axis_x_> E;
    
    //磁場作成
    NdTensor<Vec3<Value>,Axis_x_> B;
    
}