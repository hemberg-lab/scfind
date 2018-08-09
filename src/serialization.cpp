#include "Serialization.hpp"


SerializationDB::SerializationDB(): byte_pointer(0)
{}


int SerializationDB::loadByteStream(const Rcpp::RawVector& stream)
{
  this->serialized_bytestream = std::vector<unsigned char>(stream.begin(), stream.end());
  std::cerr << this->serialized_bytestream.size();
  return 1;
    
}

  
Rcpp::RawVector SerializationDB::getByteStream(const EliasFanoDB& efdb)
{
  
  serialize(efdb);
  Rcpp::RawVector r_obj(this->byte_pointer);
  memcpy(&r_obj[0], &this->serialized_bytestream[0], this->byte_pointer);
    
  //Rcpp::wrap(this->serialized_bytestream);
  this->serialized_bytestream.clear();
  std::cout << r_obj.size() <<" is the size of the stream " << std::endl;
  return r_obj;
}

int SerializationDB::writeBuffer(const char* buf, int buf_len)
{
    
  expandStream(buf_len);
  memcpy(&this->serialized_bytestream[this->byte_pointer], buf, buf_len);
  byte_pointer += buf_len;
    
  // this->serialized_bytestream.insert(this->serialized_bytestream.end(), buf, buf + buf_len);
  return 0;
}


int SerializationDB::readBuffer(void* buf, int buf_len)
{
    
  memcpy(buf, &this->serialized_bytestream[this->byte_pointer], buf_len);
  this->byte_pointer += buf_len;
  // erase is cheap constant in a vector so i do not think this needs optimization
  // this->serialized_bytestream.erase(this->serialized_bytestream.begin(), this->serialized_bytestream.begin() + buf_len);
  return 0;
}
  
  
int SerializationDB::readFile(const std::string& filename)
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
  
void SerializationDB::deserializeEliasFano(EliasFano& ef, int quantization_bits)
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
  ef.expr.quantile.resize(quantization_bits * cells);
   
  read(ef.expr.mu);
  read(ef.expr.sigma);

}

void SerializationDB::binarizeEliasFano(const EliasFano& ef)
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


void SerializationDB::deserializeDB(EliasFanoDB& efdb)
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
  read(efdb.total_cells);
    
  // Read the number of bits used for the quantization of the expression level quantiles
    
  read(efdb.quantization_bits);

    
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
    efdb.gene_counts[buffer] = cell_support;
    efdb.metadata[buffer] = EliasFanoDB::EliasFanoIndex();
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
    int cell_type_id = efdb.cell_types_id.size();
    efdb.cell_types_id[ct] = cell_type_id;
    efdb.inverse_cell_type.push_back(ct);
      
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
    deserializeEliasFano(ef, efdb.quantization_bits);
    efdb.ef_data.push_back(ef);
       
  }
    
    
  // Build database
  for (auto const& r : records)
  {
    //std::cerr << r.gene << " " <<r.cell_type << " " << r.index << std::endl;
    efdb.metadata[gene_ids[r.gene]][r.cell_type] = r.index;
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
  // return efdb;
}



void SerializationDB::serialize(const EliasFanoDB& efdb)
{
  this->serialized_bytestream.clear();
  this->byte_pointer = 0;
  std::map<const EliasFano*, int> ef_ids;
  int ef_id = 0;
  for (auto const& ef : efdb.ef_data)
  {
    ef_ids[&ef] = ef_id++;
  }
  
  int version = SERIALIZATION_VERSION;
  write(version);

  // Dump the total amount of cells
  write(efdb.total_cells);
  // Dump the quantization bits for storing the expression level of the genes
  write(efdb.quantization_bits);
    
  // Number of genes present in the database
  int genes_present = efdb.gene_counts.size();
  write(genes_present);

    
  int gene_id = 0;
  // Read the gene names
  std::map<std::string, int> gene_ids;
    
  for (auto const& g : efdb.gene_counts)
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
  int cell_types_present = efdb.cell_types_id.size();
  write(cell_types_present);
    
  for (auto const& ct : efdb.inverse_cell_type)
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
  int index_size = efdb.ef_data.size();
  write(index_size);
  std::cout << "Index size " << index_size << std::endl;
    
  for (auto const& g : efdb.metadata)
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
  for (auto const& ef : efdb.ef_data)
  {
    binarizeEliasFano(ef);
  }

  std::cout << "Database deserialized!" << std:: endl;
}

