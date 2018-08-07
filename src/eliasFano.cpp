#include <Rcpp.h>
// #include <unistd.h>
#include <iostream>
#include <bitset>
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <set>
#include <cmath>
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "scfind_types.h"
#include "fp_growth.hpp"


// the bits used for the encoding
#define BITS 32
#define SERIALIZATION_VERSION 4
// #define DEBUG

typedef std::pair<unsigned short, std::bitset<BITS> > BitSet32;
typedef std::vector<bool> BoolVec;


class CellType 
{
 public:
  std::string name;
  int total_cells;
  size_t operator()() const
  {
    return std::hash<std::string>{}(name);
  }

};



struct CellTypeCompare
{
  using is_transparent = std::true_type;
  bool operator()(const CellType& lhs, const CellType& rhs) const
  {
    return lhs.name < rhs.name;
  }
  
  bool operator()(const CellType& lhs, const std::string& name) const
  {
    return lhs.name < name;
  }
  
  bool operator()(const std::string& name, const CellType& rhs) const
  {
    return name < rhs.name;
  }
};


typedef float IDFtype;



// struct that holds the quantization vector
typedef struct{
  double mu;
  double sigma;
  std::vector<bool> quantile;
} Quantile;

typedef struct
{
  BoolVec H;
  BoolVec L;
  int l;
  IDFtype idf; // tfidf
  Quantile expr;
} EliasFano;

typedef struct
{
  int gene;
  int cell_type;
  int index;
} IndexRecord;


struct Cell_ID
{
  unsigned int num;
  const int cell_type;
  // Hashing
  size_t operator()() const
  {
    return std::hash<int>{}(cell_type) ^ std::hash<unsigned int>{}(num);
  }
  
  size_t operator==(const struct Cell_ID& obj) const 
  {
    return (num == obj.num) && (cell_type == obj.cell_type);
  }
};


typedef struct Cell_ID CellID;

namespace std
{
  template<>
  struct hash<CellID>
  {
    size_t operator()(const CellID& obj) const
    {
      // Return overloaded operator
      return obj();
    }
  };

  template<>
  struct hash<CellType>
  {
    size_t operator()(const CellType& obj) const
    {
      return obj();
    }
  };
  
}

// Highly Recommended
// TODO(Nikos) refactor code with these to avoid nasty bug that will misalign the whole bytestream
// ;)


std::string str_join( const std::vector<std::string>& elements, const char* const separator)
{
  switch (elements.size())
  {
    case 0:
      return "";
    case 1:
      return *(elements.begin());
    default:
      std::ostringstream os; 
      std::copy(elements.cbegin(), elements.cend() - 1, std::ostream_iterator<std::string>(os, separator));
      os << *elements.crbegin();
      return os.str();
  }
}

inline std::bitset<BITS> int2bin_core(const unsigned int id)
{
  // we have a fixed BITS bitset so we have to offset 
  // the leading zeros whenever we need to access the elements
  return std::bitset<BITS>(id);
}


inline BitSet32 int2bin_bounded(const unsigned int id, unsigned int min_bit_length)
{
  // Check if number is larger than the desired size in bits
  unsigned short bit_set_length = id == 0 ? 0 : __builtin_clz(id);
  bit_set_length = bit_set_length > min_bit_length ? bit_set_length : min_bit_length;
  return make_pair(bit_set_length, int2bin_core(id));
}


// These functions work with 1-based indexes on arrays, 
// __builtin_clz has undefined behavior with 0
inline BitSet32 int2bin(unsigned int id)
{
  return make_pair(__builtin_clz(id), int2bin_core(id));
}



inline double normalCDF(double x, double mu, double sigma)
{
  // this is an inline function for the cdm of the normal distribution
  // it depends on the cmath library where it contains the erfc function
  // it return a value ranging from zero to one 
  
  return 1 - (0.5 * erfc( (x - mu)/ (sigma * M_SQRT1_2) ));
   
}


// Accepts a vector, transforms and returns a quantization logical vector
// This function aims for space efficiency of the expression vector
Quantile lognormalcdf(std::vector<int> ids, const Rcpp::NumericVector& v, unsigned int bits)
{
  Quantile expr;
  expr.mu = std::accumulate(ids.begin(),ids.end(), 0, [&v](const double& mean, const int& index){
      return  mean + v[index - 1];
    }) / ids.size();
  
  expr.sigma = sqrt(std::accumulate(ids.begin(), ids.end(), 0, [&v, &expr](const double& variance, const int& index){
        return pow(expr.mu - v[index - 1], 2);
      }) / ids.size());
  // initialize vector with zeros
  expr.quantile.resize(ids.size() * bits, 0);
  //std::cerr << "Mean,std" << expr.mu << "," << expr.sigma << std::endl;
  //std::cerr << "ids size " << ids.size() << " v size " << v.size() << std::endl;
  int expr_quantile_i = 0;
  for (auto const& s : ids)
  {
    unsigned int t = round(normalCDF(v[s], expr.mu, expr.sigma) * (1 << bits));  
    std::bitset<BITS> q = int2bin_core(t);
    for (int i = 0; i < bits; ++i)
    {
       expr.quantile[expr_quantile_i++] = q[i];
      
    }
  }
  return expr;
 }


int byteToBoolVector(const std::vector<char> buf, std::vector<bool>& bool_vec)
{
  bool_vec.resize((buf.size() << 3), 0);
  int c = 0;
  for (auto const& b : buf)
  {
    for ( int i = 0; i < 8; i++)
    {
      bool_vec[c++] = ((b >> i) & 1);
    }
  }
    
}


class EliasFanoDB;
RCPP_EXPOSED_CLASS(EliasFanoDB)

typedef int EliasFanoID;
typedef int CellTypeID;

class EliasFanoDB
{

 public:
  typedef std::unordered_map<CellTypeID, EliasFanoID> EliasFanoIndex;
  // gene -> cell type -> eliasFano
  typedef std::unordered_map<std::string, EliasFanoIndex > CellTypeIndex;
  typedef std::deque<EliasFano> ExpressionMatrix;
  
  // Store the gene metadata, gene support in cells at the index
  typedef std::map<std::string, unsigned int> GeneIndex;

 // private:
  CellTypeIndex metadata;
  ExpressionMatrix ef_data;
  std::map<CellType,int,CellTypeCompare> cell_types_id;
  std::deque<CellType> inverse_cell_type;
  GeneIndex gene_counts;
  unsigned int total_cells;
  long byte_pointer;
  unsigned char quantization_bits;
  std::vector<unsigned char> serialized_bytestream;
  
  
  bool global_indices;
  int warnings;

  int getSizeBoolVector(const std::vector<bool>& v)
  {
    int size = v.size() / 8;
    if (v.size() % 8 != 0)
    {
      ++size;
    }
    return size;
  }
  
  template<typename T>
  int read(T& val)
  {
    memcpy(&val, &this->serialized_bytestream[byte_pointer], sizeof(T));
    // SLOOOOOOOOOOOOOW
    // this->serialized_bytestream.erase(this->serialized_bytestream.begin(), this->serialized_bytestream.begin() + sizeof(T));
    byte_pointer += sizeof(T);
    // if(read(fd, &val, sizeof(T)) < 0)
    //   {
        // std::cerr << "Something went wrong" << std::endl;
    //     return 1;
    //   }
    return 0;
  }

  template<typename T>
  int expandStream(const T& new_length)
  {
    if(this->byte_pointer + new_length >= this->serialized_bytestream.size())
    {
      // Is there a chance the buffer to be bigger than 16MB
      // Grow the buffer by 16MB
      std::cout << "expanding stream" << std::endl;
      this->serialized_bytestream.resize(this->serialized_bytestream.size() + (1 << 24));
    }
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


  int writeBuffer(const char* buf, int buf_len)
  {
    
    expandStream(buf_len);
    memcpy(&this->serialized_bytestream[this->byte_pointer], buf, buf_len);
    byte_pointer += buf_len;
    
    // this->serialized_bytestream.insert(this->serialized_bytestream.end(), buf, buf + buf_len);
    return 0;
  }
 
  int readBuffer(void* buf, int buf_len)
  {
    
    memcpy(buf, &this->serialized_bytestream[this->byte_pointer], buf_len);
    this->byte_pointer += buf_len;
    // erase is cheap constant in a vector so i do not think this needs optimization
    // this->serialized_bytestream.erase(this->serialized_bytestream.begin(), this->serialized_bytestream.begin() + buf_len);
    return 0;
  }
  
  
  int readFile(const std::string& filename)
  {
    FILE *fp = fopen(filename.c_str(), "rb");
    fseek(fp, 0, SEEK_END);
    long bytes = ftell(fp);
    // Rewind the fp to the start of the file so we can read the contents
    fseek(fp, 0, SEEK_SET); 
    this->serialized_bytestream.resize(bytes);
    fread(&this->serialized_bytestream[0], bytes, 1, fp);
    fclose(fp);
    return 0;
  }
  
  void deserializeEliasFano(EliasFano& ef)
  {
    int cells;
    
    read(cells);
    read(ef.l);
    read(ef.idf);
    
    int H_size, L_size, quant_size;
    unsigned char test;
    read(test);
    assert(test == 0xFF);
    read(H_size);
    read(L_size);
    read(quant_size);
    
    // Multiply by eight, shift three positions
    // std::cout << "l = " << ef.l << std::endl;
    // std::cout << "Read " << " H " << H_size << std::endl;
    // std::cout << "Read " << " L " << L_size << std::endl;
    // std::cout << "Read " << " quant " << quant_size << std::endl;
    
    
    std::vector<char> buffer;

    // Read H
    buffer.resize(H_size, 0);
    readBuffer(&buffer[0], buffer.size());
    byteToBoolVector(buffer, ef.H);
    // Trim padding to last 1
    for(int i = ef.H.size() - 1; i > 0; i--)
    {
      if(ef.H[i])
      {
        // Trim vector size and exit 
        ef.H.resize(i + 1, 0);
        break;
      }
    }

    
    // Read L
    buffer.resize(L_size, 0);
    readBuffer(&buffer[0], buffer.size());
    
    byteToBoolVector(buffer, ef.L);
    // Trim buffer to the right size
    ef.L.resize(cells * ef.l);
    
    
    // Read expression values
    buffer.resize(quant_size, 0);
    readBuffer(&buffer[0], buffer.size());
    byteToBoolVector(buffer, ef.expr.quantile);
    
    
    // Take care of the byte quantization
    // Resize container to the right size
    ef.expr.quantile.resize(this->quantization_bits * cells);
   
    read(ef.expr.mu);
    read(ef.expr.sigma);

  }

  void binarizeEliasFano(const EliasFano& ef)
  {
    // int i = 0;
    int cells = ef.L.size() / ef.l;
    write(cells);
    write(ef.l);
    write(ef.idf);
    int H_size = getSizeBoolVector(ef.H);
    int L_size = getSizeBoolVector(ef.L);
    int quant_size = getSizeBoolVector(ef.expr.quantile);
    unsigned char t = 0xFF;
    write(t);
    write(H_size);
    write(L_size);
    write(quant_size);
    
    std::vector<char> H_buf, L_buf, expr;
    // Initialize memory to zeros
    H_buf.resize(H_size, 0);
    L_buf.resize(L_size, 0);
    expr.resize(quant_size, 0);

    for (size_t i = 0; i < ef.H.size(); i++)
    {
      H_buf[i / 8] |= ef.H[i] << (i % 8);
    }

    for ( size_t i = 0; i < ef.L.size(); i++)
    {
      L_buf[i / 8] |= ef.L[i] << (i % 8);
    }

    
    for (size_t i = 0; i < ef.expr.quantile.size(); i++)
    {
      expr[i / 8] |= ef.expr.quantile[i] << (i % 8);
    }

    writeBuffer(&H_buf[0], H_buf.size());
    writeBuffer(&L_buf[0], L_buf.size());
    writeBuffer(&expr[0], expr.size());
    write(ef.expr.mu);
    write(ef.expr.sigma);
    
    return;
  }
  
  void dumpGenes()
  {
    // for (auto const g : gene_counts)
    // {
    //   // std::cout << g.first << " " ;
    // }
    
    std::cout << "Total Genes:" << gene_counts.size() << std::endl;
    for (auto const c : this->cell_types_id)
    {
      std::cout << c.first.name << " ";
    }
    std::cout << "Total Cell types" << std::endl;
    
    for (auto const& m : metadata)
    {
      std::cout << "Gene:" << m.first << std::endl;
      long total = 0;
      for (auto const& t : m.second)
      {
          auto ct = inverse_cell_type[t.first];
          // str.size();
          total += ct.name.size() + 4;
      }
      std::cout << total <<std:: endl;
    }
  }


  void deserializeDB()
  {

    // Read the database version
    int version;
    read(version);
    std::cout << "Version " << version << std::endl;
    if (version != SERIALIZATION_VERSION)
    {
      std::cerr << "The model that the database was stored is not matching the of the existing version that the database will be built. Exiting now" << std::endl;
      return;
    }

    
    // Read the total cells stored in the database
    read(this->total_cells);
    
    // Read the number of bits used for the quantization of the expression level quantiles
    
    read(this->quantization_bits);

    
    int genes_present;
    read(genes_present);
    

    if (not(genes_present > 0 and genes_present < 100000))
    {
      std::cerr << "something went wrong" << std::endl;
      return;
    }

    // Read gene names
    char buffer[256];
    std::vector<std::string> gene_ids;
    gene_ids.reserve(genes_present);
    for (int i = 0; i < genes_present; ++i)
    {
      int cell_support;
      unsigned char gene_name_length;
      read(gene_name_length);
      readBuffer(buffer, gene_name_length);
      // Insert end of character string
      buffer[gene_name_length] = '\0';
      read(cell_support);
      this->gene_counts[buffer] = cell_support;
      metadata[buffer] = EliasFanoIndex();
      // std::cout << "gene" << buffer << std::endl;
      gene_ids.push_back(buffer);
    }

    // Read cell type names
    int cell_types_present;
    read(cell_types_present);

    if(not(cell_types_present > 0 and cell_types_present < 50000))
    {
      std::cerr << "something went wrong with the cell types" << std::endl;
      return;
    }
    

    unsigned int cell_support;
    std::vector<CellTypeID> cell_type_ids;
    cell_type_ids.reserve(cell_types_present);
    for (int i = 0; i < cell_types_present; ++i)
    {
      unsigned char cell_type_length;
      int total_cells;
      read(cell_type_length);
      readBuffer(buffer, cell_type_length);
      read(total_cells);
      // Insert end of string character
      buffer[cell_type_length] = '\0';
      CellType ct;
      ct.name = buffer;
      ct.total_cells = total_cells;
      int cell_type_id = this->cell_types_id.size();
      this->cell_types_id[ct] = cell_type_id;
      this->inverse_cell_type.push_back(ct);
      
    }

    int index_size;
    read(index_size);
    
    std::vector<IndexRecord> records;
    for (int i = 0; i < index_size; ++i)
    {
      IndexRecord record;
      read(record);
      records.push_back(record);
    }
    
    
    for (int i = 0; i < index_size; i++)
    {
      EliasFano ef;
      deserializeEliasFano(ef);
      ef_data.push_back(ef);
       
    }
    
    
    // Build database
    for (auto const& r : records)
    {
      //std::cerr << r.gene << " " <<r.cell_type << " " << r.index << std::endl;
      metadata[gene_ids[r.gene]][r.cell_type] = r.index;
    }
    
#ifdef DEBUG
    std::cout << "Total Cells " << this->total_cells << std::endl;
    std::cout << "Quantization bits " << (unsigned int)this->quantization_bits << std::endl;
    std::cout << "Present genes " << genes_present << std::endl;
    std::cout << "Present cell_types " << cell_types_present << std::endl;
    std::cout << "Genes support " << gene_counts.size() << std::endl;
    std::cout << "Index size " << index_size << std::endl;
    std::cout << "Bytes left to read:" << this->serialized_bytestream.size() - byte_pointer << std::endl;
#endif
    this->serialized_bytestream.clear();
    
  }

  void clearDB()
  {
    // Clear the database
    metadata.clear();
    cell_types_id.clear();
    inverse_cell_type.clear();
    gene_counts.clear();
    serialized_bytestream.clear();
  }

  

  void loadFromFile(const std::string& filename)
  {
    clearDB();
    readFile(filename);
    deserializeDB();
  }


  int loadByteStream(const Rcpp::RawVector& stream)
  {
    this->serialized_bytestream = std::vector<unsigned char>(stream.begin(), stream.end());
    
    deserializeDB();
    std::cout << "Database deserialized!" << std:: endl;
    return 1;
    
  }


  void serialize()
  {
    this->serialized_bytestream.clear();
    this->byte_pointer = 0;
    std::map<const EliasFano*, int> ef_ids;
    int ef_id = 0;
    for (auto const& ef : ef_data)
    {
      ef_ids[&ef] = ef_id++;
    }
    
    int version = SERIALIZATION_VERSION;
    write(version);

    // Dump the total amount of cells
    write(this->total_cells);
    // Dump the quantization bits for storing the expression level of the genes
    write(this->quantization_bits);
    
    // Number of genes present in the database
    int genes_present = gene_counts.size();
    write(genes_present);

    
    int gene_id = 0;
    // Read the gene names
    std::map<std::string, int> gene_ids;
    
    for (auto const& g : gene_counts)
    {
      // Size of cell type can be from 0 to 255
      unsigned char gene_name_size = g.first.size();
      
      write(gene_name_size);
      writeBuffer(&g.first[0], gene_name_size);
      // Dump gene idf
      write(g.second);
      gene_ids[g.first] = gene_id++;
    }
    // Dump cell types

    int cell_type_id = 0;
    int cell_types_present = this->cell_types_id.size();
    write(cell_types_present);
    
    for (auto const& ct : this->inverse_cell_type)
    {
      //assign unique id
      unsigned char cell_type_name_size = ct.name.size();
      
      write(cell_type_name_size);
      writeBuffer(&ct.name[0], cell_type_name_size);
      write(ct.total_cells);
      //std::cout << ct << std::endl;
    }
    
    // Dump records
    // Write size of index, if it is 1:1 relation it should be consistent
    int index_size = ef_data.size();
    write(index_size);
    std::cout << "Index size " << index_size << std::endl;
    
    for (auto const& g : metadata)
    {

      for (auto const& ct : g.second)
      {
        IndexRecord record;
        record.gene = gene_ids[g.first];
        record.cell_type = ct.first;
        record.index = ct.second;
        write(record);
      }
    }

    // Dump raw data
    for (auto const& ef : ef_data)
    {
      binarizeEliasFano(ef);
    }

    return;
  }
  
  Rcpp::RawVector getByteStream()
  {
    serialize();
    Rcpp::RawVector r_obj(this->byte_pointer);
    memcpy(&r_obj[0], &this->serialized_bytestream[0], this->byte_pointer);
    
    //Rcpp::wrap(this->serialized_bytestream);
    this->serialized_bytestream.clear();
    std::cout << r_obj.size() <<" is the size of the stream " << std::endl;
    for( int i = 0 ; i < 13; i++)
    {
      std::cerr << (unsigned int)r_obj[i] << std::endl;
    }
    return r_obj;
  }
  void serializeToFile(const std::string& filename)
  {
    serialize();
    FILE* fp = fopen(filename.c_str(), "wb");
    fwrite(&this->serialized_bytestream[0], this->byte_pointer, 1, fp);
    std::cout << "Database was dumped on " << filename << std::endl;
    fclose(fp);
    this->serialized_bytestream.clear();
    return;
  }

  void insertToDB(int ef_index, const std::string& gene_name, const CellType& cell_type)
  {
    if (ef_index == -1)
    {
      // Something went wrong so do not do anything
      this->warnings++;
      return;
    }

    EliasFano* ef = &(this->ef_data[ef_index]);
    
    if (metadata.find(gene_name) == metadata.end())
    {
      metadata[gene_name] = EliasFanoIndex();
    }
    
    auto celltype_record = this->cell_types_id.find(cell_type);
    if (celltype_record == this->cell_types_id.end())
    {
      int id = this->cell_types_id.size();
      this->cell_types_id[cell_type] = id;
      this->inverse_cell_type.push_back(cell_type);
      metadata[gene_name][id] = ef_index;
    }else
    {
      metadata[gene_name][this->cell_types_id[cell_type]] = ef_index;
    }
    
    
  } 

  long eliasFanoCoding(const std::vector<int>& ids, const Rcpp::NumericVector& values) 
  {
    
    if(ids.empty())
    {
      return -1;
    }
    int items = values.size();

    
    EliasFano ef;
    ef.l = int(log2(items / (float)ids.size()) + 0.5) + 1;
    ef.idf = log2(items / (float)ids.size());
    int l = ef.l;

    int prev_indexH = 0;
    ef.L.resize(l * ids.size(), false);

    BoolVec::iterator l_iter = ef.L.begin();
    ef.expr = lognormalcdf(ids, values, this->quantization_bits);
    
    
    for (auto expr = ids.begin(); expr != ids.end(); ++expr)
    {
      BitSet32 c = int2bin_bounded(*expr, l);
    
      for( int i = 0; i < l; i++, ++l_iter)
      {
        // TODO(Nikos) optimize this ? check if c.second[i] is false and THEN assign
        *l_iter = c.second[i];
      }
      //Use a unary code for the high bits
      unsigned int upper_bits = (*expr >> l);
      unsigned int m =  ef.H.size() + upper_bits - prev_indexH + 1;
      prev_indexH = upper_bits;
      ef.H.resize(m, false);
      ef.H[m - 1] = true;
    }
    ef_data.push_back(ef);
    // std::cerr <<"New index" << ef_data.size() - 1 << std::endl;
    // return the index of the ef_data in the deque
    return ef_data.size() - 1;
  }

  std::vector<int> eliasFanoDecoding(const EliasFano& ef)
  {
    
    std::vector<int> ids(ef.L.size() / ef.l);

    // This step inflates the vector by a factor of 8
    std::vector<char> H;
    H.reserve(ef.H.size());
    H.insert(H.end(), ef.H.begin(), ef.H.end());
    

    unsigned int H_i = 0;
    // Warning: Very very dodgy I might want to replace this with a check in the loop
    auto prev_it = H.begin() - 1;
    int i = 0;
    for (auto true_it = std::find(H.begin(), H.end(), true); 
         true_it != H.end() && i < ids.size(); 
         true_it = std::find(true_it + 1, H.end(), true), ++i)
    {
      size_t offset  = std::distance(prev_it, true_it);
      prev_it = true_it;
      H_i += offset - 1;
      int id = H_i << ef.l;
      for (unsigned short k = 0; k < ef.l; ++k)
      {
        id |= (ef.L[(i * ef.l) + k] << k); 
      }
      ids[i] = id;
    }
    return ids;
  }



 public:

  // constructor
  EliasFanoDB(): 
    global_indices(false), 
    warnings(0), 
    total_cells(0), 
    quantization_bits(2), 
    byte_pointer(0)
  {
    
  }

  // This is invoked on slices of the expression matrix of the dataset 
  long encodeMatrix(const std::string& cell_type_name, const Rcpp::NumericMatrix& gene_matrix)
  {
    int items = gene_matrix.ncol();
    Rcpp::CharacterVector genes = Rcpp::rownames(gene_matrix);

    CellType cell_type;
    cell_type.name = cell_type_name;
    cell_type.total_cells = gene_matrix.ncol();

    std::vector<std::string> gene_names;
    gene_names.reserve(genes.size());


    // Increase the cell number present in the index
    this->total_cells += gene_matrix.ncol();
    
    for(Rcpp::CharacterVector::iterator gene_it = genes.begin(); gene_it != genes.end(); ++gene_it)
    {
      std::string gene_name = Rcpp::as<std::string>(*gene_it);
      gene_names.push_back(gene_name);
      // If this is the first time this gene occurs initialize a counter
      if (this->gene_counts.find(gene_name) == this->gene_counts.end())
      {
        this->gene_counts[gene_name] = 0;
      }
    }
    
    for(unsigned int gene_row = 0; gene_row < gene_matrix.nrow(); ++gene_row)
    {
      
      const Rcpp::NumericVector& expression_vector = gene_matrix(gene_row, Rcpp::_);
      std::deque<int> sparse_index;
      int i = 0;
      for (Rcpp::NumericVector::const_iterator expr = expression_vector.begin(); expr != expression_vector.end(); ++expr)
      {
        i++; 
        if (*expr > 0)
        {
          sparse_index.push_back(i);
        }
      }
      if ( i == 0)
      {
        warnings++;
        continue;
      }
      std::vector<int> ids(sparse_index.begin(), sparse_index.end());
      this->gene_counts[gene_names[gene_row]] += ids.size();
      
      insertToDB(eliasFanoCoding(ids, expression_vector), gene_names[gene_row], cell_type);
    }
    return 0;
    //std::cerr << "Total Warnings: "<<warnings << std::endl;
  }

  Rcpp::List total_genes()
  {
    Rcpp::List t;
    for(auto & d : metadata)
    {
      t.push_back(Rcpp::wrap(d.first));
    }
    return t;
  }

  Rcpp::NumericVector getTotalCells()
  {
    return Rcpp::wrap(this->total_cells);
  }

  Rcpp::NumericVector getCellTypeSupport(Rcpp::CharacterVector& cell_types)
  {
    std::vector<std::string> cts = Rcpp::as<std::vector<std::string>>(cell_types);
    std::vector<int> ct_support;
    ct_support.reserve(cts.size());
    for (auto const& ct : cts)
    {
      // TODO(Nikos) fix this, otherwise we will get a nice segfault no error control
      auto cit = this->cell_types_id.find(ct);
      if(cit != this->cell_types_id.end())
      {
        ct_support.push_back(this->inverse_cell_type[cit->second].total_cells);
      }
      else
      {
        ct_support.push_back(0);
      }
    }
    
    
    return Rcpp::wrap(ct_support);
  }
  

  Rcpp::List queryGenes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active)
  {
    Rcpp::List t;
    for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
    {
      
      std::string gene_name = Rcpp::as<std::string>(*it);
      //t.add(Rcpp::wrap(gene_names[i]), Rcpp::List::create());
      Rcpp::List cell_types;
      
      if (metadata.find(gene_name) == metadata.end())
      {
        
        std::cout << "Gene " << gene_name << " not found in the index " << std::endl;
        continue;
      }
      std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);
      const auto& gene_meta = metadata[gene_name];
      for (auto const& dat : gene_meta)
      {
        CellType current_cell_type = this->inverse_cell_type[dat.first];
        
        std::string dataset = current_cell_type.name.substr(0, current_cell_type.name.find("."));
        auto ct_find = std::find(datasets.begin(), datasets.end(), dataset);
        
        if (ct_find == datasets.end())
        {
          continue;
        }
        std::vector<int> ids = eliasFanoDecoding(ef_data[dat.second]);
        cell_types[current_cell_type.name] = Rcpp::wrap(ids);
      }
      t[gene_name] = cell_types;
      
    }

    return t;
  }
  
  size_t dataMemoryFootprint()
  {
    size_t bytes = 0;
    for(auto & d : ef_data)
    {
      bytes += int((d.H.size() / 8) + 1);
      bytes += int((d.L.size() / 8) + 1);
      bytes += int((d.expr.quantile.size() / 8) + 1);
      
    }
    bytes += ef_data.size() * 32; // overhead of l idf and deque struct
    return bytes;
  }

  size_t dbMemoryFootprint()
  {
    size_t bytes = dataMemoryFootprint();

    std::cout << "Raw elias Fano Index size " << bytes / (1024 * 1024) << "MB" << std::endl;
    
    for(auto& d : metadata)
    {
      bytes += d.first.size();
      bytes += d.second.size() * 8;
    }
    return bytes;
  }


  // And query
  Rcpp::List findCellTypes(const Rcpp::CharacterVector& gene_names, const Rcpp::CharacterVector& datasets_active)
  {
    
    std::unordered_map<CellTypeID, std::set<std::string> > cell_types;
    std::vector<std::string> genes;

    std::vector<std::string> datasets = Rcpp::as<std::vector<std::string>>(datasets_active);

    // Fast pruning if there is not an entry we do not need to consider
    for (Rcpp::CharacterVector::const_iterator it = gene_names.begin(); it != gene_names.end(); ++it)
    {
      std::string gene_name = Rcpp::as<std::string>(*it);
      
      // check if gene exists in the database
      auto db_it = metadata.find(gene_name);
      if (db_it == metadata.end())
      {
        std::cerr << gene_name << " is ignored, not found in the index"<< std::endl;
        continue;
      }

      // iterate cell type
      for (auto const& ct_it : db_it->second)
      {
        CellType ct = this->inverse_cell_type[ct_it.first];
        
        // Remove cells if not in the selected datasets
        std::string ct_dataset = ct.name.substr(0, ct.name.find("."));
        auto find_dataset = std::find(datasets.begin(), datasets.end(), ct_dataset);
        // check if the cells are in active datasets
        if (find_dataset == datasets.end())
        {
          continue;
        }

        if (cell_types.find(ct_it.first) == cell_types.end())
        {
          cell_types[ct_it.first] = std::set<std::string>();
        }
        cell_types[ct_it.first].insert(gene_name);
      }
      genes.push_back(gene_name);
    }
    
    // Store the results here
    Rcpp::List t;

    for (auto const& ct : cell_types)
    {
      bool empty_set = false;
      bool initial_set = true;
      if (ct.second.size() != genes.size())
      {
        continue;
      }
      
      auto g_it = ct.second.begin();
      std::vector<int> ef = eliasFanoDecoding(this->ef_data[metadata[*g_it][ct.first]]);

      for (++g_it; g_it != ct.second.end();++g_it)
      {
        auto cells = eliasFanoDecoding(this->ef_data[metadata[*g_it][ct.first]]);
        
        std::vector<int> intersected_cells;
        std::set_intersection(ef.begin(), 
                              ef.end(), 
                              cells.begin(), 
                              cells.end(), 
                              std::back_inserter(intersected_cells));
        ef = intersected_cells;
        if (ef.empty())
        {
          break;
        }
      }
      
      if(!ef.empty())
      {
        t[this->inverse_cell_type[ct.first].name] = Rcpp::wrap(ef);
      }
    }
    return t;
  }


  // TODO(Nikos) this function can be optimized.. It uses the native quering mechanism
  // that casts the results into native R data structures
  Rcpp::DataFrame findMarkerGenes(const Rcpp::CharacterVector& gene_list, const Rcpp::CharacterVector datasets_active,unsigned int min_support_cutoff = 5)
  {
    // Store the results in this list
    Rcpp::List results;
    
    std::vector<std::string> query;
    std::vector<double> query_scores;
    std::vector<int> query_cell_type_cardinality;
    std::vector<int> query_cell_cardinality;
    std::vector<int> query_gene_cardinality;

    std::map<std::string, std::map<int, Transaction> > cells;


    int cells_present = 0;
    
    // Perform an OR query on the database as a first step
    const Rcpp::List genes_results = queryGenes(gene_list, datasets_active);

    // Start inversing the list to a cell level
    const Rcpp::CharacterVector& gene_names = genes_results.names();
    for (auto const& gene : gene_names)
    {
      const auto gene_name = Rcpp::as<std::string>(gene);
      // Gene hits contains the cell type hierarchy
      const Rcpp::List& gene_hits = genes_results[gene_name];
      const Rcpp::CharacterVector& cell_type_names = gene_hits.names();
      for (auto const& _ct : cell_type_names)
      {
        const auto ct = Rcpp::as<std::string>(_ct);
        
        if (cells.find(ct) == cells.end())
        {
          cells[ct] = std::map<int, Transaction>();
        }
        
        //

        std::vector<unsigned int> ids  = Rcpp::as<std::vector<unsigned int> >(gene_hits[ct]);
        auto& cell_index = cells[ct];
        // For all the hits
        for (auto const& id : ids)
        {
          // search for the id in the cell type space!
          if (cell_index.find(id) == cell_index.end())
          {
            // insert new cell
            cell_index[id] = Transaction();
            cells_present++;
          }
          // Add gene hit in the cell
          cell_index[id].push_back(gene_name);
        }
      }
    }

    std::cerr << "Query Done: found " << cells_present << " rules" << std::endl;
    
    // Collect all transactions for fp-growth
    std::vector<Transaction> transactions;
    transactions.reserve(cells_present);
    for (auto & ct : cells)
    {
      for (auto & cl : ct.second)
      {
        // Maybe sort?        
        // std::sort(cl.second.begin(), cl.second.end());
        if (cl.second.size() != 1)
        {
          transactions.push_back(std::vector<std::string>(cl.second.begin(), cl.second.end()));
        }
      }
    }
    
    std::cerr << transactions.size() << " transactions" << std::endl;
    // Run fp-growth algorithm
    const FPTree fptree{transactions, min_support_cutoff};
    const std::set<Pattern> patterns = fptree_growth(fptree);

    // Patterns contain genesets
    std::vector<std::pair<std::string, double> > tfidf;
    
    // Iterate through the calculated frequent patterns
    for (auto const& item : patterns)
    {
      Rcpp::List gene_query;
      const auto& gene_set = item.first;
      int fp_support = item.second;
      
      if (gene_set.size() == 1)
      {
        continue;
      }
      // Support scaled by genes
      // double query_base_score = 1 + (log(item.second) * gene_set.size());
      
      // Get the cell type intersection of the involved cells
      // TODO(Nikos) optimize wrapping 
      std::map<std::string, std::vector<int> > ct_map;
      
      // First gene init set for intersections
      auto git = gene_set.begin();
     
      const Rcpp::List& first_gene = genes_results[*git];
      const std::vector<std::string> initial_cell_types = Rcpp::as<std::vector<std::string> >(first_gene.names());
      for (auto const& ct : initial_cell_types)
      {
        ct_map[ct] = Rcpp::as< std::vector<int> >(first_gene[ct]);
      }
      
      // std::cout << *git << ", cell types " << ct_map.size() << std::endl;

      for (++git; git != gene_set.end(); ++git)
      {
        const Rcpp::List& g_res = genes_results[*git];
        std::vector<std::string> cell_types = Rcpp::as<std::vector<std::string> >(g_res.names());
        std::vector<std::string> cell_types_in_ct_map;
        cell_types_in_ct_map.reserve(ct_map.size());
        for (auto const& ct : ct_map)
        {
          cell_types_in_ct_map.push_back(ct.first);
        }
        
        std::sort(cell_types_in_ct_map.begin(), cell_types_in_ct_map.end());
        std::sort(cell_types.begin(), cell_types.end());
        
        
        std::vector<std::string> ct_intersection;
        std::set_intersection(cell_types_in_ct_map.begin(),
                              cell_types_in_ct_map.end(),
                              cell_types.begin(),
                              cell_types.end(),
                              std::back_inserter(ct_intersection));

        // Update ct_map with the values that have been thrown out cause of the intersection
        for(auto m_it = ct_map.begin(); m_it != ct_map.end();)
        {
          auto f_it = std::find(ct_intersection.begin(),ct_intersection.end(), m_it->first);
          if (f_it == ct_intersection.end())
          {
            m_it = ct_map.erase(m_it);
          }
          else
          {
            ++m_it;
          }
        }
        
      }
      
      if (ct_map.empty())
      {
        continue;
      }
      // std::cout << " Common cell types accross genes" << ct_map.size() << std::endl;
      for (auto const& _ct : ct_map)
      {
        auto ct = _ct.first;
        for(auto& gene : gene_set){
          auto& current_cells = ct_map[ct];
          const Rcpp::List& g_res = genes_results[gene];
          std::vector<int> cells;
          std::vector<int> cells_in_gct = Rcpp::as<std::vector<int>>(g_res[ct]);
          // std::sort(current_cells.begin(), current_cells.end());
          // std::sort(cells_in_gct.begin(), cell_in_gct.end
          std::set_intersection(
            current_cells.begin(),
            current_cells.end(),
            cells_in_gct.begin(),
            cells_in_gct.end(),
            std::back_inserter(cells));
          
          // This invalidates the current cells reference, careful there
          ct_map[ct] = cells;
          
          if (cells.empty())
          {
            ct_map.erase(ct_map.find(ct));
            break;
          }
        }
      }
      
      
      

      // Remove all empty records
      int total_support = 0;
      for (auto it = ct_map.begin(); it != ct_map.end();)
      {
        if(it->second.empty())
        {
          // std::cerr << "deleted" << it->first << std::endl;
          it = ct_map.erase(it);
          
        }
        else
        {
          total_support+= it->second.size();
          ++it;
        }
      }
      
      // Continue to the next pattern
      if (ct_map.empty())
      {
        continue;
      }

      // std::cout << "FP-Growth support " << fp_support << ", Discovered support: " << total_support << std::endl;
     
      // Give query a score
      std::string view_string = str_join(std::vector<Item>(gene_set.begin(), gene_set.end()), ",");
      
      
      double query_score = 0;
      int cells_in_query = 0;
      std::map<std::string, double> ct_tfidf;
      for ( auto const& _ct : ct_map)
      {
        const std::string& ct = _ct.first;
        auto cct =  this->cell_types_id.find(_ct.first);
        double tfidf = 0;
        for (auto const&  g : gene_set)
        {

          tfidf += this->ef_data[this->metadata[g][cct->second]].idf;
        }
        query_score += tfidf * log(_ct.second.size());
        cells_in_query += _ct.second.size();
      }
      // Get the mean ct size
      query_score /= ct_map.size();
      
      
      query_cell_cardinality.push_back(cells_in_query);
      query_gene_cardinality.push_back(gene_set.size());
      query.push_back(view_string);
      query_scores.push_back(query_score);
      query_cell_type_cardinality.push_back(ct_map.size());
      
    }

    std::vector<int> query_rank(query_scores.size());    
    std::iota(query_rank.begin(), 
              query_rank.end(), 
              1);
    
    std::sort(query_rank.begin(), 
              query_rank.end(), 
              [&query_scores](const int& i1, const int& i2){
                return query_scores[i1] < query_scores[i2];
              });

    // Dump the list
    return Rcpp::DataFrame::create(Rcpp::Named("Genes") = Rcpp::wrap(query_gene_cardinality),
                                   Rcpp::Named("Query") = Rcpp::wrap(query),
                                   Rcpp::Named("Rank") = Rcpp::wrap(query_rank),
                                   Rcpp::Named("Cells") = Rcpp::wrap(query_cell_cardinality),
                                   Rcpp::Named("Cell Types") = Rcpp::wrap(query_cell_type_cardinality));
  }


  int dbSize()
  {
    std::cout << metadata.size() << "genes in the DB" << std::endl;
    return ef_data.size();
    
  }

  int sample(int index)
  {
    auto iter = metadata.begin();
    for(int i = 0; i < index; i++, ++iter);

    std::cout << "Gene: " << iter->first << std::endl;
    for(auto const& ct : iter->second)
    {
      std::cout << "Cell Type:" << ct.first << std::endl;
      auto v = eliasFanoDecoding(ef_data[ct.second]);
      for( auto const& cell : v)
      {
        std::cout << cell << ", ";
      }
      std::cout << std::endl;
      
    }
    return 0;
    
  }

  std::vector<int> decode(int index)
  {
    if(index >= dbSize())
    {
      std::cerr << "Invalid index for database with size "<< dbSize() << std::endl;
      return std::vector<int>();
    }
    return eliasFanoDecoding(ef_data[index]);
  }
  
  int insertNewCellType(const CellType& cell_type)
  {
    // Fetch the last id
    int id = this->inverse_cell_type.size();

    if ( this->cell_types_id.find(cell_type) != this->cell_types_id.end())
    {

      std::cerr << "This should not happen!! Duplicate Cell Type: " << cell_type.name << std::endl;
      id = this->cell_types_id[cell_type];
    }
    else
    {
      id = this->inverse_cell_type.size();
      this->inverse_cell_type.push_back(cell_type);
      this->cell_types_id[cell_type] = id;
    }
    return id;
    
  }

  int mergeDB(const EliasFanoDB& db)
  {    
    EliasFanoDB extdb(db);
    
    // the DB will grow by this amount of cells
    this->total_cells += extdb.total_cells;

    // Insert new cell types in the database
    for (auto const& ct : extdb.inverse_cell_type)
    {
      insertNewCellType(ct);
    }
    
    // Iterate through the data model
    for ( auto& gene: extdb.metadata)
    {
      // Update cell counts for the individual gene
      if (this->gene_counts.find(gene.first) == this->gene_counts.end())
      {
        this->gene_counts[gene.first] = 0;
      }
      
      this->gene_counts[gene.first] += extdb.gene_counts[gene.first];
      
      
      // if gene does not exist yet initialize entry in metadata
      if(metadata.find(gene.first) == metadata.end())
      {
        metadata[gene.first] = EliasFanoIndex();
      }
      
      for( auto& ct : gene.second)
      {
        int new_id = ef_data.size();
        // Push the new elias fano index in the database
        ef_data.push_back(extdb.ef_data[ct.second]);
        // Update with the new entry
        int cell_type_id = this->cell_types_id[extdb.inverse_cell_type[ct.first]];
        metadata[gene.first][cell_type_id] = new_id;
      }
    }
    return 0;
  }

};

RCPP_MODULE(EliasFanoDB)
{
  Rcpp::class_<EliasFanoDB>("EliasFanoDB")
    .constructor()
    .method("indexMatrix", &EliasFanoDB::encodeMatrix)
    .method("queryGenes", &EliasFanoDB::queryGenes)
    .method("dbSize", &EliasFanoDB::dbSize)
    .method("decode", &EliasFanoDB::decode)
    .method("mergeDB", &EliasFanoDB::mergeDB)
    .method("sample", &EliasFanoDB::sample)
    .method("findCellTypes", &EliasFanoDB::findCellTypes)
    .method("efMemoryFootprint", &EliasFanoDB::dataMemoryFootprint)
    .method("dbMemoryFootprint", &EliasFanoDB::dbMemoryFootprint)
    .method("findMarkerGenes", &EliasFanoDB::findMarkerGenes)
    .method("serializeToFile", &EliasFanoDB::serializeToFile)
    .method("loadFromFile",&EliasFanoDB::loadFromFile)
    .method("dumpGenes", &EliasFanoDB::dumpGenes)
    .method("getByteStream", &EliasFanoDB::getByteStream)
    .method("loadByteStream", &EliasFanoDB::loadByteStream)
    .method("getCellsInDB", &EliasFanoDB::getTotalCells)
    .method("getCellTypeSupport", &EliasFanoDB::getCellTypeSupport);
  
}




