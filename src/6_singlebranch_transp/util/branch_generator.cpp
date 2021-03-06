#include<iostream>
#include <vector>
#include <iterator> // for ostream_iterator
#include <fstream>

class branch{
    std::vector<double> p1={0,0,0};
    std::vector<double> p2={0,0,0};
    int ndiv=0;
    void generate(void);
    void print(void);    
public:   
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    double val1=0;
    double val2=0;
    branch (std::vector<double>, std::vector<double>,int, double, double );
    };

branch::branch (std::vector<double> up1, std::vector<double> up2,int div,
                double valpt1,double valpt2): 
            p1(up1),p2(up2),ndiv(div), val1(valpt1), val2(valpt2)
    {
     generate();
    }
    
void branch::generate(){
    std::cout<< "Branch generation"<<std::endl;
    double dx= (p2[0]-p1[0])/ndiv;
    double dy= (p2[1]-p1[1])/ndiv;
    double dz= (p2[2]-p1[2])/ndiv;
    x.resize(ndiv+1);
    y.resize(ndiv+1);
    z.resize(ndiv+1);
    for(int i=0;i<ndiv+1;i++){
      x[i]=dx*i + p1[0];y[i]=dy*i+ p1[1];z[i]=dz*i+ p1[2];
      }
    }

void print(std::vector<branch> branches){
std::ofstream out;
out.open("./seg.pts");
out<<"BEGIN_LIST" <<std::endl;
for(int b=0;b<branches.size();b++){
out<<"BEGIN_ARC" <<std::endl;
out<<"BC DIR  " << branches[b].val1 <<std::endl;
out<<"BC DIR  " << branches[b].val2 <<std::endl;
out<<"\t" << 1<<"\t"<<branches[b].x[0]<< "\t" << branches[b].y[0]<< "\t" << branches[b].z[0]<<"\t"<<"begin"<<std::endl;
out<<"\t"  << 1<<"\t"
          << branches[b].x[branches[b].x.size()-1] << "\t"  
          << branches[b].y[branches[b].x.size()-1] << "\t" 
          << branches[b].z[branches[b].x.size()-1] << "\t"
          << "end" <<std::endl;
for (int i=1; i< branches[b].x.size()-1;i++)
            out<<"\t" << 1<<"\t"
                      << branches[b].x[i]<< "\t" 
                      << branches[b].y[i]<< "\t" 
                      << branches[b].z[i]<< "\t" 
                     << "point" <<std::endl;
out<<"END_ARC"<<std::endl;
}
out<<"END_LIST"<<std::endl;
out.close();
}

int main (){
    std::cout<<"Generator of branches"<<std::endl;
    
    std::vector<double> p1={0.7,0.33,0}; double pt1d=0.77;
    std::vector<double> p2={0.7,0.33,0.1}; double pt2d=1.;
    branch b1(p1,p2,10, pt1d,pt2d);

    std::vector<double> p3={0.65,0.33,0};double pt3d=0.77;
    std::vector<double> p4={0.1,0.33,0.1};double pt4d=0.;
    branch b2(p3,p4,10, pt3d,pt4d);

std::vector<branch> branches;
branches.push_back(b1);
branches.push_back(b2);
print(branches); 
    
    std::cout<<"Endl of Generator of branches"<<std::endl;
    
    return 1;}
