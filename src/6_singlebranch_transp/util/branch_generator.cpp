#include<iostream>
#include <vector>

class branch{
    std::vector<int> p1={0,0,0};
    std::vector<int> p2={0,0,0};
    int ndiv=0;
    generate(void);
    public:   
    branch (std::vector<int>, std::vector<int>,int );
    };

branch::branch (std::vector<int> up1, std::vector<int> up2,int div): 
            p1(up1),p2(up2),ndiv(div)
    {
     generate();
    }
    
branch::generate(){
    std::cout<< "Branch generation"<<std::endl;
    
    }
int main (){
    std::cout<<"Generator of branches"<<std::endl;
    
    std::vector<int> p1={0,0,0};
    std::vector<int> p2={0,1,0};
    branch b1(p1,p2,10);
    
    std::cout<<"Endl of Generator of branches"<<std::endl;
    
    return 1;}
