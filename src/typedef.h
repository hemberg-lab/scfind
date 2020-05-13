#pragma once
#include <vector>
#include <bitset>

#include "const.h"


typedef std::pair<unsigned short, std::bitset<BITS> > BitSet32;
typedef std::vector<bool> BoolVec;
typedef int EliasFanoID;
typedef int CellTypeID;


using Item = std::string;
using Transaction = std::vector<Item>;
using Pattern = std::pair<std::set<Item>, uint64_t>;

typedef struct
{
  int gene;
  int cell_type;
  int index;
} IndexRecord;

class GeneMeta
{
public:
  int total_reads;
  GeneMeta();
  void merge(const GeneMeta& other);
};

class CellMeta
{
public:                                               
  int reads;
  int features;
  int getReads() const 
  {
    return reads;
  }
  
  int getFeatures() const
  {
    return features;
  }
  CellMeta();
};


class CellType
{
public:
  std::string name;
  int total_cells;
  int getTotalCells()const 
  {
    return total_cells;
  }
};



typedef struct
{
  double mu;
  double sigma;
  std::vector<bool> quantile;
} Quantile;


class EliasFano
{
public:
  BoolVec H;
  BoolVec L;
  int l;
  float idf; // tfidf
  Quantile expr;
  int getSize() const
  {
    return L.size() / l;
  }
};


class CellID
{
public:
  CellTypeID cell_type;
  int cell_id;
  CellID(CellTypeID, int);
  bool operator==(const CellID& obj) const
  {
    return (obj.cell_type == cell_type) && (obj.cell_id == cell_id);

  }
};

namespace std
{
  template<>
  struct hash<CellID>
  {
    inline size_t operator()(const CellID& cid) const
    {
      return hash<CellTypeID>()(cid.cell_type) ^ hash<int>()(cid.cell_id);
    }
  };

}

typedef struct
{
  int tp;
  int fp;
  int tn;
  int fn;
  float inv_precision() const
  {
    return  (tp + fp) / float(tp);
  }
  float inv_recall() const
  {
    return (fn + tp) /  float(tp);
  }

  float recall() const
  {
    return 1 / inv_recall();
  }

  float precision() const
  {
    return 1 / inv_precision();
  }

  float f1() const
  {
    return 2/(inv_precision() + inv_recall());
  }
} CellTypeMarker;




