#include "geometry/geometry.hpp"

#include <list>
#include <iostream>

using namespace carpio; 

typedef PointChain_<double, 2> PointChain;
typedef PolygonPartition_<double, 2> PP;

typedef IOFile_Geometry_<double, 2> IOF;


int main(int argc, char** argv) {
    std::string filname(argv[1]);
    PointChain pc;
    IOF::ReadPointChain("./" + filname, pc);
    pc.set_close();
    std::cout << "File = " << filname   << std::endl;
    std::cout << "size = " << pc.size() << std::endl;
    std::list<PointChain> lres;
    std::cout << "EerClipping ...\n";
    int r = PP::EerClipping(pc, lres);
    std::cout << "return code = " << r << std::endl;
    int count = 0;
    for(auto& pc: lres){
        IOF::WritePointChain("./res_"+ filname +"/" + fmt::format("{0:03d}", count)
                             , pc);
        count++;
    }
    
}
