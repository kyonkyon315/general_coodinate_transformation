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
    static constexpr int get_num_objects(){
        return sizeof...(Objects);
    }
    template<int I>
    using element = typename std::tuple_element<I, std::tuple<Objects...>>::type;
};

#endif //PACK_H