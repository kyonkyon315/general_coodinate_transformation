#ifndef AXIS_H
#define AXIS_H

template<int LABEL,int LENGTH>
class Axis {
    private:
    int val;
    public:
    static constexpr int label = LABEL;
    static constexpr int length = LENGTH;
    int operator=(int r){
        return val=r;
    }
    int operator()(){
        return val;
    }
    static int get_length(){
        return length;
    }
};

#endif //AXIS_H
