#include "bmp.h"
using namespace std;

int main() {
    // load the file. The constructor now does most of the work
    char* c = "test.bmp";
    char* c1 = "test2.bmp";
    char* c2 = "test3.bmp";
    BitMap example_bmp(c);
    BitMap example_bmp1(c1); 
    BitMap example_bmp2(c2);
    /*for (int i = 0; i < 3; i++){
        example_bmp.getPixel(i,0);
    }*/
    int opc;
    cin>>opc;
    // get the vector <R,G,B> for the pixel at (1,1)
    //example_bmp.dispPixelData();
    if(opc==1){
        example_bmp.encrypt("ad",3);
        example_bmp1.encrypt("ad",3);
        example_bmp2.encrypt("ad",3);
    }
    else{
        example_bmp.decrypt("ad",3);
        example_bmp1.decrypt("ad",3);
        example_bmp2.decrypt("ad",3);
    }
    //example_bmp.writePixel(10,10,255,0,0);
}