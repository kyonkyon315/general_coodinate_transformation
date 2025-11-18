#ifndef INPUT_LAYER_H
#define INPUT_LAYER_H
#include <vector>
#include <string>

using Index = int;

template<Index OutputDimension>
class InputLayer{
public:
    using Value = double;

    static constexpr Index input_dimension 
        = OutputDimension;
    static constexpr Index output_dimension
        = input_dimension;
    static constexpr Index output_id_left
        = 0;
    
    static constexpr Index param_dimension = 0;

    static constexpr Index param_id_left
        = 0;

    static std::string get_network_structure(){
        return "InputLayer<" + std::to_string(input_dimension) + ">";
    }
};
#endif //INPUT_LAYER_H