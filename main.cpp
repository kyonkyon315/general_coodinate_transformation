#include <iostream>
#include <vector>
using Value = double;
double PI=3.14159265;
#include <array>
#include <iostream>
#include "Timer.h"

template<int LENGTH>
struct Coordinate {
    private:
    int val;
    public:
    static constexpr int length = LENGTH;
    int operator=(int r){
        return val=r;
    }
    int operator()(){
        return val;
    }
    int get_length(){
        return length;
    }
};


template<typename... Axes>
class Variable {
private:
    // 軸ごとの長さをコンパイル時に収集
    static constexpr std::array<int, sizeof...(Axes)> shape = {Axes::length...};

    // 総要素数をコンパイル時計算
    static constexpr int total_size = []() constexpr {
        int prod = 1;
        for (auto s : shape) prod *= s;
        return prod;
    }();

    std::vector<Value> data;

    template<size_t I = 0, typename... Idx>
    constexpr int flatten_index(int i, Idx... rest) const noexcept {
        if constexpr (I == sizeof...(Axes) - 1)
            return i;
        else
            return i + shape[I] * flatten_index<I + 1>(rest...);
    }

public:
    Variable(Axes ...args){
        data.resize(total_size);
    }
    
    template<typename... Idx>
    inline Value& operator()(Idx... indices) noexcept {
        static_assert(sizeof...(Idx) == sizeof...(Axes));
        constexpr bool valid = (sizeof...(Idx) == sizeof...(Axes));
        return data[flatten_index(0, indices...)];
    }

    template<typename... Idx>
    inline const Value& operator()(Idx... indices) const noexcept {
        static_assert(sizeof...(Idx) == sizeof...(Axes));
        return data[flatten_index(0, indices...)];
    }
};

int main(){
    Coordinate<100> xi,vr,vt;
    Coordinate<1000> vp;
    Variable p(xi,vr,vt,vp);
    p(1,2,3,4)=4;
    std::cout<<p(1,2,3,4)<<" "<<p(3,3,3,4)<<"\n";
    Timer timer;
    timer.start();
    for(int i=0;i<100;++i){
    for(int j=0;i<100;++i){
    for(int k=0;i<100;++i){
    for(int l=0;i<1000;++i){
        p(i,j,k,l) =0.;
        p(i,j,k,l)*=(double)i;
        p(i,j,k,l)*=(double)j;
        p(i,j,k,l)*=(double)k;
        p(i,j,k,l)*=(double)l;
    }
    }
    }
    }
    timer.stop();
    std::cout<<timer<<"\n";
}