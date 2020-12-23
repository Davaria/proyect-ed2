/*
 *****************************************************************************
 *                                                                           *
 *                          Platform Independent                             *
 *                     Bitmap Image Reader Writer Library                    *
 *                                                                           *
 * Author: Arash Partow - 2002                                               *
 * URL: http://partow.net/programming/bitmap/index.html                      *
 *                                                                           *
 * Note: This library only supports 24-bits per pixel bitmap format files.   *
 *                                                                           *
 * Copyright notice:                                                         *
 * Free use of the Platform Independent Bitmap Image Reader Writer Library   *
 * is permitted under the guidelines and in accordance with the most current *
 * version of the MIT License.                                               *
 * http://www.opensource.org/licenses/MIT                                    *
 *                                                                           *
 *****************************************************************************
*/


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "bitmap_image.h"


void test01()
{
   std::string file_name1("imagenes/img1.bmp"); //1
   std::string file_name2("imagenes/img2.bmp");  //2
   std::string file_name3("imagenes/img3.bmp");  //3
   std::string file_name4("imagenes/img4.bmp");  //4
   std::string file_name5("imagenes/img5.bmp");  //5
   std::string file_name6("imagenes/img6.bmp");  //6

   std::cout<<"Imagen"<<std::endl;
   std::string n;
   std::cin>>n;

   std::string file_name = "imagenes/img" + n + ".bmp";
   bitmap_image image(file_name);
   if (!image)
   {
      printf("test01() - Error - Failed to open '%s'\n",file_name.c_str());
      return;
   }
   int k;
   std::cout<<"(1) Encriptacion por permutacion"<<std::endl;
   std::cout<<"(2) Encriptacion por funciones caoticas"<<std::endl;
   std::cin>>k;
   switch (k){
   case 1:
      std::cout<<"Semilla para encriptar"<<std::endl;
      std::cin>>k;
      image.encriptar_permutacion(k);
      std::cout<<"Semilla para desencriptar"<<std::endl;
      std::cin>>k;
      image.desencriptar_permutacion(k);
      break;
   case 2:   
      std::cout<<"Clave para encriptar"<<std::endl;
      std::cin>>n;
      image.encriptar_caotico(n);
      std::cout<<"Clave para desencriptar"<<std::endl;
      std::cin>>n;
      image.desencriptar_caotico(n);
      break;
   default:
      std::cout<<"Error de opcion"<<std::endl;
      break;
   }
}


int main()
{
   test01();
   return 0;
}
