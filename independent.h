// --- (微分ゼロ) を表すクラス ---
using Value = double;
class Independent
{
public:
    Independent(){}
    // 常に 0.0 を返す
    static constexpr Value operator()(int calc_x_,int calc_vr,int calc_vt,int calc_vp){
        return 0.0;
    }
};
