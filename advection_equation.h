#ifndef ADVECTION_EQUATION_H
#define ADVECTION_EQUATION_H
using Value = double;
template<typename TargetFunction,typename Operators,typename Advections>
class AdvectionEquation
{
public:
    AdvectionEquation(TargetFunction& target_func,const Operators& operators,const Advections&){
        
    }
    template<typename CalcAxis>
    void solve(Value dt){

    }
};
#endif //ADVECTION_EQUATION_H