// --- (微分ゼロ) を表すクラス ---
class Independent
{
public:
    Independent(){}
    // 常に 0.0 を返す
    static constexpr Value at(int calc_x_,int calc_vr,int calc_vt,int calc_vp){
        return 0.0;
    }
};
