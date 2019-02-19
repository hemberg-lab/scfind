#pragma once

#include "EliasFano.h"

#include <Rcpp.h>
#include <iostream>
#include <cassert>
#include <vector>


class SerializationDB
{
public:
  friend class EliasFanoDB;
  unsigned long byte_pointer;
  std::vector<unsigned char> serialized_bytestream;
  
  template<typename T>
    int read(T& val)
  {
    memcpy(&val, &this->serialized_bytestream[byte_pointer], sizeof(T));
    byte_pointer += sizeof(T);
    return 0;
  }

  template<typename T>
  int expandStream(const T& new_length)
  {
    if (this->byte_pointer + new_length >= this->serialized_bytestream.size())
    {
      // TODO(Nikos) Is there a chance the buffer to be bigger than 16MB ?
      // Grow the buffer by 16MB
      // std::cout << "expanding stream" << std::endl;
      this->serialized_bytestream.resize(this->serialized_bytestream.size() + (1 << 24));
    }
    return 0;
  }

  template<typename T>
    int write(const T& val)
  {
    // const unsigned char* ptr = reinterpret_cast<const unsigned char*>(&val);
    expandStream(sizeof(T));
    memcpy(&this->serialized_bytestream[this->byte_pointer], &val, sizeof(T));
    byte_pointer += sizeof(T);
    return 0;
  }

  SerializationDB();

  int loadByteStream(const Rcpp::RawVector& stream);
  
  Rcpp::RawVector getByteStream(const EliasFanoDB&);

  int writeBuffer(const char* buf, int buf_len);

  int readBuffer(void* buf, int buf_len);
  
  int readFile(const std::string& filename);
  
  void deserializeEliasFano(EliasFano& ef, int quantization_bits);

  void binarizeEliasFano(const EliasFano& ef);

  void deserializeDB(EliasFanoDB&);

  void serialize(const EliasFanoDB&);


};
