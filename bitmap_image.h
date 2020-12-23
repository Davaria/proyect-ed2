/*
 *****************************************************************************
 *                                                                           *
 *                          Platform Independent                             *
 *                    Bitmap Image Reader Writer Library                     *
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


#ifndef INCLUDE_BITMAP_IMAGE_HPP
#define INCLUDE_BITMAP_IMAGE_HPP

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>
#include <random>


class bitmap_image
{
public:

   void encriptar_permutacion(int k){
      unsigned char matriz[height_][width_*3];

      int cont=0,n=0;
      unsigned char R,G,B;
      //iniciar valores
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            get_pixel(j/3,i,R,G,B);
            matriz[i][j+2] = R;//R
            matriz[i][j+1] = G;//G
            matriz[i][j] = B;//B
         }
      }

      //Encriptacion horizontal
      cont=0;n=0;
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            n=(cont*k)%width_;
            set_pixel(n,i,matriz[i][j+2],matriz[i][j+1],matriz[i][j]);
            cont++;
         }
      }

      //iniciar valores
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            get_pixel(j/3,i,R,G,B);
            matriz[i][j+2] = R;//R
            matriz[i][j+1] = G;//G
            matriz[i][j] = B;//B
         }
      }

      //Encriptacion vertical
      int _i=0;cont=0;n=0;
      for (int i = 0; i < width_*3; i=i+3){
         for (int j = 0; j < height_; j++){
            n=(cont*k)%height_;
            set_pixel(_i,n,matriz[j][i+2],matriz[j][i+1],matriz[j][i]);
            cont++;
         }
         _i++;
      }

      std::string encriptado = "encriptado/perm_" + name_s_e();
      encrypt_name=encriptado;
      this->save_image(encriptado);
      
   }

   void desencriptar_permutacion(int k){
      unsigned char matriz[height_][width_*3];

      int cont=0,n=0;
      //Iniciar valores
      unsigned char R,G,B;
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            get_pixel(j/3,i,R,G,B);
            matriz[i][j+2] = R;//R
            matriz[i][j+1] = G;//G
            matriz[i][j] = B;//B
         }
      }

      //Desencriptar horizontal
      cont=0;n=0;
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            n=(j*k)%(width_*3);
            set_pixel(cont,i,matriz[i][n+2],matriz[i][n+1],matriz[i][n]);
            cont=(cont+1)%width_;
         }
      }

      //Copiar valores
      for (int i = 0; i < height_; i++){
         for (int j = 0; j < width_*3; j=j+3){
            get_pixel(j/3,i,R,G,B);
            matriz[i][j+2] = R;//R
            matriz[i][j+1] = G;//G
            matriz[i][j] = B;//B
         }
      }

      //Desencriptar vertical
      cont=0;n=0;
      for (int i = 0; i < width_*3; i=i+3){
         for (int j = 0; j < height_; j++){
            n=(j*k)%height_;
            set_pixel(cont,j,matriz[n][i+2],matriz[n][i+1],matriz[n][i]);
         }
         cont++;
      }

      std::string desencriptado = "desencriptado/perm_" + name_s_e();
      decrypt_name=desencriptado;
      this->save_image(desencriptado);
   }

   //Funcion para generar un numero aleatorio decimal
   double fRand(double fMin, double fMax, double semilla){
      srand(fmod(semilla,10e5));
      double f = (double)rand() / RAND_MAX;
      return fMin + f * (fMax - fMin);
   }

   //Para generar la clave
   void generar_clave(double& lambda, double& X0,std::string clave){
      std::string caracteres="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789~`!@#$%^&*()-_=+[{]}\\|;:\'\",<.>/?";
      double aux=1;
      for (int i = 0; i < clave.size(); i++){
         for (int j = 0; j < caracteres.size(); j++){
            if (clave[i]==caracteres[j]){
               aux=aux*(j+1);
               break;
            }
         }
      }
      lambda = fRand(3.57,4,aux);
      X0 = fRand(0.01, 0.99, aux);
   }
   
   //Para encriptar por medio de una funcion caotica
   void encriptar_caotico(std::string clave){
      if(clave.size()<5){
         std::cout<<"clave de minimo 5 caracteres"<<std::endl;
         return;
      }
      //Inicializar los parametros
      int L = data_.size();
      int d = floor(log10(L)) + 3;
      int _R = floor(128/(log2(9*pow(10,d-1)))) + 1;
      int r = 1;
      double _lambda;
      double X0;
      generar_clave(_lambda,X0,clave);
      std::cout<<_lambda<<" "<<X0<<std::endl;
      double _alfa = pow(10,d);
      int i,U,F=L-2,j,Q1,Q2;
      double X;
      
      while(r<=_R){
         i=0;U=L-1;X=X0;
         while (i<F){
            X = _lambda*X*(1-X);
            j = i+1+(int)(_alfa*X)%U;
            Q1 = data_[i];
            Q2 = data_[j];
            data_[i] = Q2;
            data_[j] = Q1;
            i = i + 1;
            U = U - 1;
         }
         r = r + 1;
      }
      
      std::string encriptado = "encriptado/caot_" + name_s_e();
      encrypt_name = encriptado;
      this->save_image(encriptado);
   }

   //Desencriptar
   void desencriptar_caotico(std::string clave){
      if(clave.size()<5){
         std::cout<<"clave de minimo 5 caracteres"<<std::endl;
         return;
      }
      //Inicializar los parametros
      int L = data_.size();
      int d = floor(log10(L)) + 3;
      int _R = floor(128/(log2(9*pow(10,d-1)))) + 1;
      int r = _R;
      double _lambda;
      double X0;
      generar_clave(_lambda,X0,clave);
      std::cout<<_lambda<<" "<<X0<<std::endl;
      double _alfa = pow(10,d);
      int i,U,F=L-2,j,Q1,Q2;
      double X;
      int W[L];

      while(r>0){
         i=0;U=L-1;X=X0;
         while (i<F){
            X = _lambda*X*(1-X);
            j = i+1+(int)(_alfa*X)%U;
            W[i] = j;
            i = i + 1;
            U = U - 1;
         }
         j = F - 1;
         while (j>=0){
            i = W[j];
            Q1 = data_[j];
            Q2 = data_[i];
            data_[j]=Q2;
            data_[i]=Q1;
            j = j - 1;
         }
         r = r - 1;
      }

      std::string desencriptado = "desencriptado/caot_" + name_s_e();
      decrypt_name = desencriptado;
      this->save_image(desencriptado);
   }

   std::string get_encrypt_name(){
      return encrypt_name;
   } 

   std::string get_decrypt_name(){
      return decrypt_name;
   }

   std::string name_s_e(){
      std::string n=file_name_;

      for (int i = 0; i < n.size(); i++){
         if (file_name_[i-1]=='/'){
            n.erase(0,i);
         }
      }

      return n;
   }

   enum channel_mode {
                        rgb_mode = 0,
                        bgr_mode = 1
                     };

   enum color_plane {
                       blue_plane  = 0,
                       green_plane = 1,
                       red_plane   = 2
                    };

   struct rgb_t
   {
      unsigned char   red;
      unsigned char green;
      unsigned char  blue;
   };

   bitmap_image()
   : file_name_(""),
     width_          (0),
     height_         (0),
     row_increment_  (0),
     bytes_per_pixel_(3),
     channel_mode_(bgr_mode)
   {}

   bitmap_image(const std::string& filename)
   : file_name_(filename),
     width_          (0),
     height_         (0),
     row_increment_  (0),
     bytes_per_pixel_(0),
     channel_mode_(bgr_mode)
   {
      load_bitmap();
   }

   bitmap_image(const unsigned int width, const unsigned int height)
   : file_name_(""),
     width_ (width ),
     height_(height),
     row_increment_  (0),
     bytes_per_pixel_(3),
     channel_mode_(bgr_mode)
   {
      create_bitmap();
   }

   bitmap_image(const bitmap_image& image)
   : file_name_(image.file_name_),
     width_    (image.width_    ),
     height_   (image.height_   ),
     row_increment_  (0),
     bytes_per_pixel_(3),
     channel_mode_(bgr_mode)
   {
      create_bitmap();
      data_ = image.data_;
   }

   void set_file(const std::string& filename){
      file_name_ = filename;
      load_bitmap();
   }

   bitmap_image& operator=(const bitmap_image& image)
   {
      if (this != &image)
      {
         file_name_       = image.file_name_;
         bytes_per_pixel_ = image.bytes_per_pixel_;
         width_           = image.width_;
         height_          = image.height_;
         row_increment_   = 0;
         channel_mode_    = image.channel_mode_;
         create_bitmap();
         data_ = image.data_;
      }

      return *this;
   }

   inline bool operator!()
   {
      return (data_.size()   == 0) ||
             (width_         == 0) ||
             (height_        == 0) ||
             (row_increment_ == 0);
   }

   inline void clear(const unsigned char v = 0x00)
   {
      std::fill(data_.begin(), data_.end(), v);
   }

   inline unsigned char red_channel(const unsigned int x, const unsigned int y) const
   {
      return data_[(y * row_increment_) + (x * bytes_per_pixel_ + 2)];
   }

   inline unsigned char green_channel(const unsigned int x, const unsigned int y) const
   {
      return data_[(y * row_increment_) + (x * bytes_per_pixel_ + 1)];
   }

   inline unsigned char blue_channel (const unsigned int x, const unsigned int y) const
   {
      return data_[(y * row_increment_) + (x * bytes_per_pixel_ + 0)];
   }

   inline void red_channel(const unsigned int x, const unsigned int y, const unsigned char value)
   {
      data_[(y * row_increment_) + (x * bytes_per_pixel_ + 2)] = value;
   }

   inline void green_channel(const unsigned int x, const unsigned int y, const unsigned char value)
   {
      data_[(y * row_increment_) + (x * bytes_per_pixel_ + 1)] = value;
   }

   inline void blue_channel (const unsigned int x, const unsigned int y, const unsigned char value)
   {
      data_[(y * row_increment_) + (x * bytes_per_pixel_ + 0)] = value;
   }

   inline unsigned char* row(unsigned int row_index) const
   {
      return const_cast<unsigned char*>(&data_[(row_index * row_increment_)]);
   }

   inline void get_pixel(const unsigned int x, const unsigned int y,
                         unsigned char& red,
                         unsigned char& green,
                         unsigned char& blue) const
   {
      const unsigned int y_offset = y * row_increment_;
      const unsigned int x_offset = x * bytes_per_pixel_;
      const unsigned int offset   = y_offset + x_offset;

      blue  = data_[offset + 0];
      green = data_[offset + 1];
      red   = data_[offset + 2];
   }

   template <typename RGB>
   inline void get_pixel(const unsigned int x, const unsigned int y, RGB& colour) const
   {
      get_pixel(x, y, colour.red, colour.green, colour.blue);
   }

   inline rgb_t get_pixel(const unsigned int x, const unsigned int y) const
   {
      rgb_t colour;
      get_pixel(x, y, colour.red, colour.green, colour.blue);
      return colour;
   }

   inline void set_pixel(const unsigned int x, const unsigned int y,
                         const unsigned char red,
                         const unsigned char green,
                         const unsigned char blue)
   {
      const unsigned int y_offset = y * row_increment_;
      const unsigned int x_offset = x * bytes_per_pixel_;
      const unsigned int offset   = y_offset + x_offset;

      data_[offset + 0] = blue;
      data_[offset + 1] = green;
      data_[offset + 2] = red;
   }

   template <typename RGB>
   inline void set_pixel(const unsigned int x, const unsigned int y, const RGB& colour)
   {
      set_pixel(x, y, colour.red, colour.green, colour.blue);
   }

   inline bool copy_from(const bitmap_image& image)
   {
      if (
           (image.height_ != height_) ||
           (image.width_  != width_ )
         )
      {
         return false;
      }

      data_ = image.data_;

      return true;
   }

   inline bool copy_from(const bitmap_image& source_image,
                         const unsigned int& x_offset,
                         const unsigned int& y_offset)
   {
      if ((x_offset + source_image.width_ ) > width_ ) { return false; }
      if ((y_offset + source_image.height_) > height_) { return false; }

      for (unsigned int y = 0; y < source_image.height_; ++y)
      {
         unsigned char* itr1           = row(y + y_offset) + x_offset * bytes_per_pixel_;
         const unsigned char* itr2     = source_image.row(y);
         const unsigned char* itr2_end = itr2 + source_image.width_ * bytes_per_pixel_;

         std::copy(itr2, itr2_end, itr1);
      }

      return true;
   }

   inline unsigned int width() const
   {
      return width_;
   }

   inline unsigned int height() const
   {
      return height_;
   }

   inline unsigned int bytes_per_pixel() const
   {
      return bytes_per_pixel_;
   }

   inline unsigned int pixel_count() const
   {
      return width_ *  height_;
   }

   inline void setwidth_height(const unsigned int width,
                               const unsigned int height,
                               const bool clear = false)
   {
      data_.clear();
      width_  = width;
      height_ = height;

      create_bitmap();

      if (clear)
      {
         std::fill(data_.begin(), data_.end(), static_cast<unsigned char>(0x00));
      }
   }

   void save_image(const std::string& file_name) const
   {
      std::ofstream stream(file_name.c_str(),std::ios::binary);

      if (!stream)
      {
         std::cerr << "bitmap_image::save_image(): Error - Could not open file "
                   << file_name << " for writing!" << std::endl;
         return;
      }

      bitmap_information_header bih;

      bih.width            = width_;
      bih.height           = height_;
      bih.bit_count        = static_cast<unsigned short>(bytes_per_pixel_ << 3);
      bih.clr_important    = 0;
      bih.clr_used         = 0;
      bih.compression      = 0;
      bih.planes           = 1;
      bih.size             = bih.struct_size();
      bih.x_pels_per_meter = 0;
      bih.y_pels_per_meter = 0;
      bih.size_image       = (((bih.width * bytes_per_pixel_) + 3) & 0x0000FFFC) * bih.height;

      bitmap_file_header bfh;

      bfh.type             = 19778;
      bfh.size             = bfh.struct_size() + bih.struct_size() + bih.size_image;
      bfh.reserved1        = 0;
      bfh.reserved2        = 0;
      bfh.off_bits         = bih.struct_size() + bfh.struct_size();

      write_bfh(stream,bfh);
      write_bih(stream,bih);

      unsigned int padding = (4 - ((3 * width_) % 4)) % 4;
      char padding_data[4] = { 0x00, 0x00, 0x00, 0x00 };

      for (unsigned int i = 0; i < height_; ++i)
      {
         const unsigned char* data_ptr = &data_[(row_increment_ * (height_ - i - 1))];

         stream.write(reinterpret_cast<const char*>(data_ptr), sizeof(unsigned char) * bytes_per_pixel_ * width_);
         stream.write(padding_data,padding);
      }

      stream.close();
   }

   inline void convert_to_grayscale()
   {
      double r_scaler = 0.299;
      double g_scaler = 0.587;
      double b_scaler = 0.114;

      if (rgb_mode == channel_mode_)
      {
         std::swap(r_scaler, b_scaler);
      }

      for (unsigned char* itr = data(); itr < end(); )
      {
         unsigned char gray_value = static_cast<unsigned char>
                       (
                         (r_scaler * (*(itr + 2))) +
                         (g_scaler * (*(itr + 1))) +
                         (b_scaler * (*(itr + 0)))
                       );

         *(itr++) = gray_value;
         *(itr++) = gray_value;
         *(itr++) = gray_value;
      }
   }

   inline const unsigned char* data() const
   {
      return data_.data();
   }

   inline unsigned char* data()
   {
      return const_cast<unsigned char*>(data_.data());
   }

   inline void reverse()
   {
      unsigned char* itr1 = data();
      unsigned char* itr2 = end() - bytes_per_pixel_;

      while (itr1 < itr2)
      {
         for (std::size_t i = 0; i < bytes_per_pixel_; ++i)
         {
            unsigned char* citr1 = itr1 + i;
            unsigned char* citr2 = itr2 + i;

            std::swap(*citr1,*citr2);
         }

         itr1 += bytes_per_pixel_;
         itr2 -= bytes_per_pixel_;
      }
   }

private:

   inline const unsigned char* end() const
   {
      return data_.data() + data_.size();
   }

   inline unsigned char* end()
   {
      return const_cast<unsigned char*>(data() + data_.size());
   }

   struct bitmap_file_header
   {
      unsigned short type;
      unsigned int   size;
      unsigned short reserved1;
      unsigned short reserved2;
      unsigned int   off_bits;

      unsigned int struct_size() const
      {
         return sizeof(type     ) +
                sizeof(size     ) +
                sizeof(reserved1) +
                sizeof(reserved2) +
                sizeof(off_bits ) ;
      }

      void clear()
      {
         std::memset(this, 0x00, sizeof(bitmap_file_header));
      }
   };

   struct bitmap_information_header
   {
      unsigned int   size;
      unsigned int   width;
      unsigned int   height;
      unsigned short planes;
      unsigned short bit_count;
      unsigned int   compression;
      unsigned int   size_image;
      unsigned int   x_pels_per_meter;
      unsigned int   y_pels_per_meter;
      unsigned int   clr_used;
      unsigned int   clr_important;

      unsigned int struct_size() const
      {
         return sizeof(size            ) +
                sizeof(width           ) +
                sizeof(height          ) +
                sizeof(planes          ) +
                sizeof(bit_count       ) +
                sizeof(compression     ) +
                sizeof(size_image      ) +
                sizeof(x_pels_per_meter) +
                sizeof(y_pels_per_meter) +
                sizeof(clr_used        ) +
                sizeof(clr_important   ) ;
      }

      void clear()
      {
         std::memset(this, 0x00, sizeof(bitmap_information_header));
      }
   };

   inline bool big_endian() const
   {
      unsigned int v = 0x01;

      return (1 != reinterpret_cast<char*>(&v)[0]);
   }

   inline unsigned short flip(const unsigned short& v) const
   {
      return ((v >> 8) | (v << 8));
   }

   inline unsigned int flip(const unsigned int& v) const
   {
      return (
               ((v & 0xFF000000) >> 0x18) |
               ((v & 0x000000FF) << 0x18) |
               ((v & 0x00FF0000) >> 0x08) |
               ((v & 0x0000FF00) << 0x08)
             );
   }

   template <typename T>
   inline void read_from_stream(std::ifstream& stream,T& t)
   {
      stream.read(reinterpret_cast<char*>(&t),sizeof(T));
   }

   template <typename T>
   inline void write_to_stream(std::ofstream& stream,const T& t) const
   {
      stream.write(reinterpret_cast<const char*>(&t),sizeof(T));
   }

   inline void read_bfh(std::ifstream& stream, bitmap_file_header& bfh)
   {
      read_from_stream(stream,bfh.type     );
      read_from_stream(stream,bfh.size     );
      read_from_stream(stream,bfh.reserved1);
      read_from_stream(stream,bfh.reserved2);
      read_from_stream(stream,bfh.off_bits );

      if (big_endian())
      {
         bfh.type      = flip(bfh.type     );
         bfh.size      = flip(bfh.size     );
         bfh.reserved1 = flip(bfh.reserved1);
         bfh.reserved2 = flip(bfh.reserved2);
         bfh.off_bits  = flip(bfh.off_bits );
      }
   }

   inline void write_bfh(std::ofstream& stream, const bitmap_file_header& bfh) const
   {
      if (big_endian())
      {
         write_to_stream(stream,flip(bfh.type     ));
         write_to_stream(stream,flip(bfh.size     ));
         write_to_stream(stream,flip(bfh.reserved1));
         write_to_stream(stream,flip(bfh.reserved2));
         write_to_stream(stream,flip(bfh.off_bits ));
      }
      else
      {
         write_to_stream(stream,bfh.type     );
         write_to_stream(stream,bfh.size     );
         write_to_stream(stream,bfh.reserved1);
         write_to_stream(stream,bfh.reserved2);
         write_to_stream(stream,bfh.off_bits );
      }
   }

   inline void read_bih(std::ifstream& stream,bitmap_information_header& bih)
   {
      read_from_stream(stream,bih.size            );
      read_from_stream(stream,bih.width           );
      read_from_stream(stream,bih.height          );
      read_from_stream(stream,bih.planes          );
      read_from_stream(stream,bih.bit_count       );
      read_from_stream(stream,bih.compression     );
      read_from_stream(stream,bih.size_image      );
      read_from_stream(stream,bih.x_pels_per_meter);
      read_from_stream(stream,bih.y_pels_per_meter);
      read_from_stream(stream,bih.clr_used        );
      read_from_stream(stream,bih.clr_important   );

      if (big_endian())
      {
         bih.size          = flip(bih.size               );
         bih.width         = flip(bih.width              );
         bih.height        = flip(bih.height             );
         bih.planes        = flip(bih.planes             );
         bih.bit_count     = flip(bih.bit_count          );
         bih.compression   = flip(bih.compression        );
         bih.size_image    = flip(bih.size_image         );
         bih.x_pels_per_meter = flip(bih.x_pels_per_meter);
         bih.y_pels_per_meter = flip(bih.y_pels_per_meter);
         bih.clr_used      = flip(bih.clr_used           );
         bih.clr_important = flip(bih.clr_important      );
      }
   }

   inline void write_bih(std::ofstream& stream, const bitmap_information_header& bih) const
   {
      if (big_endian())
      {
         write_to_stream(stream,flip(bih.size            ));
         write_to_stream(stream,flip(bih.width           ));
         write_to_stream(stream,flip(bih.height          ));
         write_to_stream(stream,flip(bih.planes          ));
         write_to_stream(stream,flip(bih.bit_count       ));
         write_to_stream(stream,flip(bih.compression     ));
         write_to_stream(stream,flip(bih.size_image      ));
         write_to_stream(stream,flip(bih.x_pels_per_meter));
         write_to_stream(stream,flip(bih.y_pels_per_meter));
         write_to_stream(stream,flip(bih.clr_used        ));
         write_to_stream(stream,flip(bih.clr_important   ));
      }
      else
      {
         write_to_stream(stream,bih.size            );
         write_to_stream(stream,bih.width           );
         write_to_stream(stream,bih.height          );
         write_to_stream(stream,bih.planes          );
         write_to_stream(stream,bih.bit_count       );
         write_to_stream(stream,bih.compression     );
         write_to_stream(stream,bih.size_image      );
         write_to_stream(stream,bih.x_pels_per_meter);
         write_to_stream(stream,bih.y_pels_per_meter);
         write_to_stream(stream,bih.clr_used        );
         write_to_stream(stream,bih.clr_important   );
      }
   }

   inline std::size_t file_size(const std::string& file_name) const
   {
      std::ifstream file(file_name.c_str(),std::ios::in | std::ios::binary);
      if (!file) return 0;
      file.seekg (0, std::ios::end);
      return static_cast<std::size_t>(file.tellg());
   }

   void create_bitmap()
   {
      row_increment_ = width_ * bytes_per_pixel_;
      data_.resize(height_ * row_increment_);
   }

   void load_bitmap()
   {
      std::ifstream stream(file_name_.c_str(),std::ios::binary);

      if (!stream)
      {
         std::cerr << "bitmap_image::load_bitmap() ERROR: bitmap_image - "
                   << "file " << file_name_ << " not found!" << std::endl;
         return;
      }

      width_  = 0;
      height_ = 0;

      bitmap_file_header bfh;
      bitmap_information_header bih;

      bfh.clear();
      bih.clear();

      read_bfh(stream,bfh);
      read_bih(stream,bih);

      if (bfh.type != 19778)
      {
         std::cerr << "bitmap_image::load_bitmap() ERROR: bitmap_image - "
                   << "Invalid type value " << bfh.type << " expected 19778." << std::endl;

         bfh.clear();
         bih.clear();

         stream.close();

         return;
      }

      if (bih.bit_count != 24)
      {
         std::cerr << "bitmap_image::load_bitmap() ERROR: bitmap_image - "
                   << "Invalid bit depth " << bih.bit_count << " expected 24." << std::endl;

         bfh.clear();
         bih.clear();

         stream.close();

         return;
      }

      if (bih.size != bih.struct_size())
      {
         std::cerr << "bitmap_image::load_bitmap() ERROR: bitmap_image - "
                   << "Invalid BIH size " << bih.size
                   << " expected " << bih.struct_size() << std::endl;

         bfh.clear();
         bih.clear();

         stream.close();

         return;
      }

      width_  = bih.width;
      height_ = bih.height;

      bytes_per_pixel_ = bih.bit_count >> 3;

      unsigned int padding = (4 - ((3 * width_) % 4)) % 4;
      char padding_data[4] = { 0x00, 0x00, 0x00, 0x00 };

      std::size_t bitmap_file_size = file_size(file_name_);

      std::size_t bitmap_logical_size = (height_ * width_ * bytes_per_pixel_) +
                                        (height_ * padding)                   +
                                         bih.struct_size()                    +
                                         bfh.struct_size()                    ;

      if (bitmap_file_size != bitmap_logical_size)
      {
         std::cerr << "bitmap_image::load_bitmap() ERROR: bitmap_image - "
                   << "Mismatch between logical and physical sizes of bitmap. "
                   << "Logical: "  << bitmap_logical_size << " "
                   << "Physical: " << bitmap_file_size    << std::endl;

         bfh.clear();
         bih.clear();

         stream.close();

         return;
      }

      create_bitmap();

      for (unsigned int i = 0; i < height_; ++i)
      {
         unsigned char* data_ptr = row(height_ - i - 1); // read in inverted row order

         stream.read(reinterpret_cast<char*>(data_ptr), sizeof(char) * bytes_per_pixel_ * width_);
         stream.read(padding_data,padding);
      }
   }

   std::string  file_name_;
   std::string  encrypt_name;
   std::string  decrypt_name;
   unsigned int width_;
   unsigned int height_;
   unsigned int row_increment_;
   unsigned int bytes_per_pixel_;
   channel_mode channel_mode_;
   std::vector<unsigned char> data_;
};

typedef bitmap_image::rgb_t rgb_t;

inline bool operator==(const rgb_t& c0, const rgb_t& c1)
{
   return (c0.red   == c1  .red) &&
          (c0.green == c1.green) &&
          (c0.blue  == c1 .blue) ;
}

inline bool operator!=(const rgb_t& c0, const rgb_t& c1)
{
   return (c0.red   != c1  .red) ||
          (c0.green != c1.green) ||
          (c0.blue  != c1 .blue) ;
}

inline std::size_t hamming_distance(const rgb_t& c0, const rgb_t& c1)
{
   std::size_t result = 0;

   if (c0.red   != c1  .red) ++result;
   if (c0.green != c1.green) ++result;
   if (c0.blue  != c1 .blue) ++result;

   return result;
}

inline rgb_t make_colour(const unsigned int& red, const unsigned int& green, const unsigned int& blue)
{
   rgb_t result;

   result.red   = static_cast<unsigned char>(red  );
   result.green = static_cast<unsigned char>(green);
   result.blue  = static_cast<unsigned char>(blue );

   return result;
}

#endif
