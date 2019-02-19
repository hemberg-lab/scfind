#include "Serialization.h"


SerializationDB::SerializationDB(): byte_pointer(0)
{}


int SerializationDB::loadByteStream(const Rcpp::RawVector& stream)
{
  this->serialized_bytestream = std::vector<unsigned char>(stream.begin(), stream.end());
  Rcpp::Rcout << this->serialized_bytestream.size() / (1 << 20)  << " MB "<< std::endl;
  return 1;
    
}

  
Rcpp::RawVector SerializationDB::getByteStream(const EliasFanoDB& efdb)
{
  
  serialize(efdb);
  Rcpp::RawVector r_obj(this->byte_pointer);
  memcpy(&r_obj[0], &this->serialized_bytestream[0], this->byte_pointer);
    
  //Rcpp::wrap(this->serialized_bytestream);
  this->serialized_bytestream.clear();
  Rcpp::Rcout << r_obj.size() <<" is the size of the stream " << std::endl;
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
  // Rcpp::Rcout << "l = " << ef.l << std::endl;
  // Rcpp::Rcout << "Read " << " H " << H_size << std::endl;
  // Rcpp::Rcout << "Read " << " L " << L_size << std::endl;
  // Rcpp::Rcout << "Read " << " quant " << quant_size << std::endl;
    
    
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
  Rcpp::Rcout << "Version " << version << std::endl;
  if (version != SERIALIZATION_VERSION)
  {
    Rcpp::Rcerr << "The model that the database was stored is not matching the of the existing version that the database will be built. Exiting now" << std::endl;
    return;
  }

    
  // Read the total cells stored in the database
  read(efdb.total_cells);
    
  // Read the number of bits used for the quantization of the expression level quantiles
    
  read(efdb.quantization_bits);

    
  int genes_present;
  read(genes_present);
    

  if (not (genes_present > 0) )
  {
    Rcpp::Rcerr << "something went wrong" << std::endl;
    return;
  }

  // Read gene names
  char buffer[256];
  std::vector<std::string> gene_ids;
  gene_ids.reserve(genes_present);
  for (int i = 0; i < genes_present; ++i)
  {
    GeneMeta gene_meta;
    unsigned char gene_name_length;
    read(gene_name_length);
    readBuffer(buffer, gene_name_length);
    // Insert end of character string
    buffer[gene_name_length] = '\0';
    read(gene_meta);
    efdb.genes[buffer] = gene_meta;
    efdb.index[buffer] = EliasFanoDB::GeneContainer();
    // Rcpp::Rcout << "gene" << buffer << std::endl;
    gene_ids.push_back(buffer);
  }


  int number_of_cells;
  read(number_of_cells);
  
  for (int i = 0; i < number_of_cells; i++)
  {
    CellID cell_id(0,0);
    read(cell_id);
    CellMeta cell;
    read(cell);
    efdb.cells.insert({cell_id, cell});
  }
  


  // Read cell type names
  int cell_types_present;
  read(cell_types_present);

  if(not(cell_types_present > 0))
  {
    Rcpp::Rcerr << "something went wrong with the cell types" << std::endl;
    return;
  }
    
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
    int cell_type_id = efdb.cell_types.size();
    efdb.cell_types[buffer] = cell_type_id;
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
    //Rcpp::Rcerr << r.gene << " " <<r.cell_type << " " << r.index << std::endl;
    efdb.index[gene_ids[r.gene]][r.cell_type] = r.index;
  }
    
#ifdef DEBUG
  Rcpp::Rcout << "Total Cells " << this->total_cells << std::endl;
  Rcpp::Rcout << "Quantization bits " << (unsigned int)this->quantization_bits << std::endl;
  Rcpp::Rcout << "Present genes " << genes_present << std::endl;
  Rcpp::Rcout << "Present cell_types " << cell_types_present << std::endl;
  Rcpp::Rcout << "Index size " << index_size << std::endl;
  Rcpp::Rcout << "Bytes left to read:" << this->serialized_bytestream.size() - byte_pointer << std::endl;
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
  int genes_present = efdb.genes.size();
  write(genes_present);

    
  int gene_id = 0;
  // Read the gene names
  std::map<std::string, int> gene_ids;
    
  for (auto const& g : efdb.genes)
  {
    // Size of cell type can be from 0 to 255
    unsigned char gene_name_size = g.first.size();
      
    write(gene_name_size);
    writeBuffer(&g.first[0], gene_name_size);
    // Dump gene idf
    write(g.second);
    gene_ids[g.first] = gene_id++;
  }

  int number_of_cells = efdb.cells.size();
  write(number_of_cells);
  
  for ( auto const& cell : efdb.cells)
  {
    write(cell.first);
    write(cell.second);
  }

  // Dump cell types
  int cell_types_present = efdb.cell_types.size();
  write(cell_types_present);
    
  for (auto const& ct : efdb.inverse_cell_type)
  {
    //assign unique id
    unsigned char cell_type_name_size = ct.name.size();
      
    write(cell_type_name_size);
    writeBuffer(&ct.name[0], cell_type_name_size);
    write(ct.total_cells);
    //Rcpp::Rcout << ct << std::endl;
  }
    
  // Dump records
  // Write size of index, if it is 1:1 relation it should be consistent
  int index_size = efdb.ef_data.size();
  write(index_size);
  Rcpp::Rcout << "Index size " << index_size << std::endl;
    
  for (auto const& g : efdb.index)
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

  Rcpp::Rcout << "Database deserialized!" << std:: endl;
}

