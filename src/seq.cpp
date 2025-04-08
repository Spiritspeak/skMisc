#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector seq_composite(IntegerVector from, IntegerVector to){
  int fl=from.length();
  if(fl != to.length()){
    stop("from and to must be same length");
  }
  if(is_true(any(is_na(from))) || is_true(any(is_na(to))) || 
     is_true(any(is_infinite(from))) || is_true(any(is_infinite(to)))){
    stop("No NA or infinite values are allowed");
  }
  int ol=sum(abs(from-to)+1);
  IntegerVector out(ol);
  int offset=0;
  for(int i=0; i<fl; i++){
    IntegerVector currseq = seq(from[i],to[i]);
    IntegerVector currrange = seq(offset,offset+currseq.length()-1);
    out[currrange] = currseq;
    offset+=currseq.length();
  }
  return out;
}
