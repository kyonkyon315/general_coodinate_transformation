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
    std::cout<<model.get_network_structure()<<"\n";
    DataLoader data_loader;
    int butch_size = 64;
    int num_epock = 1000;
    Adam trainer(model,data_loader,butch_size,num_epock);

    std::vector<Value> input(2);
    std::vector<Value> ans(1);
    model.load_param("parameters/" + model.get_network_structure());
    for(;;){
        std::cin>>input[0]>>input[1];
        model.forward(input);
        model.get_output(ans);
        std::cout<<ans[0]<<"\n";
    }
    /*
    trainer.learn();

    model.save_param("parameters/" + model.get_network_structure());

#include "../../Timer.h"
    Timer timer;
    std::vector<Value> ans(1);
    std::vector<Value> test1 = {0.2,0.8};
    timer.start();
    for(int i=0;i<100000000;i++){
        model.forward(test1);
        model.get_output(ans);
    }
    timer.stop();
    std::cout<<timer<<"\n";
    */
}