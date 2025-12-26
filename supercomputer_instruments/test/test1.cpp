#include "../axis.h"
#include "../axis_instantiator.h"
#include <iostream>
using Axis_x_ = Axis<0,10,3,3,3>;
using Axis_y_ = Axis<0,10,4,3,3>;
using Axis_z_ = Axis<0,10,5,3,3>;
int main(){
    int rank =100;
    for(;;){
        std::cin>>rank;
        std::cout<<"\n";
        try{
            auto [axis_x,axis_y,axis_z] = axis_instantiator<Axis_x_,Axis_y_,Axis_z_>(rank);
            std::cout<<axis_x.block_id<<"\n";
            std::cout<<axis_y.block_id<<"\n";
            std::cout<<axis_z.block_id<<"\n";
        }
        catch(const std::out_of_range& e){
            std::cout<<e.what()<<"\n";
        }
    }
}