#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<bool> int2bin(int n) {
    std::vector<bool> b;
    int rem;
    while (n>=0) {
        if (n==0) {
            b.insert(b.begin(), false);
            break;
        }
        else if (n==1) {
            b.insert(b.begin(), true);
            break;
        }
        else {
            b.insert(b.begin(), (n % 2)==1);
            n /= 2;
        }
    }
    return b;
}

// [[Rcpp::export]]
Rcpp::List eliasFanoCodingCpp(NumericVector ids, int l) {
    std::vector<bool> L, H;
    //Append 0 to the start of the list
    for (int i=0; i<ids.size(); i++) {
        std::vector<bool> c = int2bin(ids[i]);
        if (c.size()<l) { //pre-pend zeros to make length l
            std::vector<bool> tmp (l-c.size(), false);
            c.insert(c.begin(), tmp.begin(), tmp.end());
        }
        L.insert(L.end(), c.begin()+c.size()-l, c.end()); //append the lower bits
        //Use a unary code for the high bits
        int m = (int)(ids[i]/pow(2, l)) - (int)(ids[i-1]/pow(2, l));
        std::vector<bool> h (m+1, false);
        h[m] = true;
        H.insert(H.end(), h.begin(), h.end());
    }
    return List::create(_["H"] = H, _["L"] = L);
}

// [[Rcpp::export]]
NumericVector eliasFanoDecodingCpp(LogicalVector H, LogicalVector L, int m, int l, int ih, int il) {
    NumericVector ids(m);
    int nZeros = 0, nCellsFound = 0, i = ih, j = 0, prevH = 0;
    while (nCellsFound<m) {
        if (!H[i])
            nZeros++;
        else {
            for (int k=0; k<l; k++) {//Calculate the value stored in the low bits
                if (L[il+j*l+k])
                    ids[nCellsFound] += pow(2, l-k-1);
            }
            //Calculate the value stored in the high bits
            std::vector<bool> h = int2bin(nZeros + prevH);
            for (int k=0; k<h.size(); k++) {
                if (h[k])
                    ids[nCellsFound] += pow(2, h.size()-k+l-1);
            }
            prevH = (int)(ids[nCellsFound]/pow(2, l));
            j++;
            nZeros = 0;
            nCellsFound++;
        }
        i++;
    }
    return ids;
}

//////////////////////////////////////////////////////////////////


//Rcpp::List createIndex(NumericMatrix d, int u) {
//Create Elias-Fano index for the genes in the matrix d. It is assumed that they all come from the same cell-type
//NumericVector m(d.nrow());
//for (int i=0; i<d.nrow(); i++) {
//Find the indexes of the cells expressing the gene

//if (m[i]>0) {
//int l = floor(log2(u/m[i]));
