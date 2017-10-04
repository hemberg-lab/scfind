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
List eliasFanoCoding(List inds_list, NumericVector l) {
    StringVector feature_names;
    NumericVector inds;
    List H_all;
    List L_all;
    
    int m;
    
    std::vector<bool> c;
    std::vector<bool> tmp;
    std::vector<bool> L, H;
    
    for(int j = 0; j < inds_list.length(); j++) {
        inds = inds_list(j);
        for (int i = 0; i < inds.size(); i++) {
            c = int2bin(inds[i]);
            if (c.size() < l[j]) { //pre-pend zeros to make length l
                std::vector<bool> tmp (l[j] - c.size(), false);
                c.insert(c.begin(), tmp.begin(), tmp.end());
            }
            L.insert(L.end(), c.begin() + c.size() - l[j], c.end()); //append the lower bits
            //Use a unary code for the high bits
            if(i > 0) {
                m = (int)(inds[i]/pow(2, l[j])) - (int)(inds[i - 1]/pow(2, l[j]));
            } else {
                m = (int)(inds[i]/pow(2, l[j]));
            }
            std::vector<bool> h (m + 1, false);
            h[m] = true;
            H.insert(H.end(), h.begin(), h.end());
        }
        H_all.push_back(H);
        L_all.push_back(L);
        H.clear();
        L.clear();
    }
    return List::create(_["H"] = H_all, _["L"] = L_all);
}

// [[Rcpp::export]]
NumericVector eliasFanoDecoding(LogicalVector H, LogicalVector L, int l) {
    NumericVector ids;
    int nZeros = 0, nCellsFound = 0, i = 0,  j = 0, prevH = 0;
    while (i < H.length()) {
        if (!H[i]) {
            nZeros++;
        }
        else {
            ids.push_back(0);
            //Calculate the value stored in the low bits
            for (int k = 0; k < l; k++) {
                if (L[j*l+k])
                    ids[nCellsFound] += pow(2, l - k - 1);
            }
            //Calculate the value stored in the high bits
            std::vector<bool> h = int2bin(nZeros + prevH);
            for (int k = 0; k < h.size(); k++) {
                if (h[k])
                    ids[nCellsFound] += pow(2, h.size() + l - k - 1);
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
