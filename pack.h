#ifndef PACK_H
#define PACK_H
#include <tuple>
template<typename... Objects>
class Pack
{
private:
    std::tuple<std::add_lvalue_reference_t<Objects>...> objects;
public:
    Pack(Objects& ...objects):
        objects(objects...)
    {}
    template<int I>
    auto& get_object(){
        return std::get<I>(objects);
    }
};

#endif //PACK_H