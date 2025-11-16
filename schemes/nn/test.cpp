#include "model.h"
#include "trainer/adam.h"
#include "dataloader.h"
#include <vector>
using Value = double;
using Input = InputLayer<2>;
using L1 = Relu<Affine<Input,10>>;
using L2 = Relu<Affine<L1,10>>;
using Output = Affine<L2,1>;

int main(){
    Model<Output> model;
    DataLoader data_loader;
    int butch_size = 64;
    int num_epock = 1000;
    Adam trainer(model,data_loader,butch_size,num_epock);

    std::vector<Value> input(2);
    trainer.learn();

#include "../../Timer.h"
    Timer timer;
    std::vector<Value> ans(1);
    std::vector<Value> test1 = {0.2,0.8};
    timer.start();
    for(int i=0;i<100000000;i++){
        //std::cout<<"input: ["<<test1[0]<<", "<<test1[1]<<"]\n";
        model.forward(test1);
        model.get_output(ans);
        //std::cout<<"pred: "<<ans[0]<<"\n";
    }
    timer.stop();
    std::cout<<timer<<"\n";
}