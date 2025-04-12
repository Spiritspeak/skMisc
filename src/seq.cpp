#include <Rcpp.h>
using namespace Rcpp;

//' Generate and concatenate multiple integer sequences
//' 
//' This generates multiple integer sequences and concatenates them.
//'
//' @param from,to the starting and (maximal) end values of the sequences. 
//' Multiple can be given.
//'
//' @returns An integer vector of multiple concatenated integer sequences.
//' @export
//' @author Sercan Kahveci
//'
//' @examples
//' seq_composite(from=c(1,7),to=c(3,8))
//' 
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

// [[Rcpp::export]]
NumericVector carryforward_numeric(NumericVector x){
  int tl=x.length();
  if(tl<2){
    return x;
  }
  NumericVector out = clone(x);
  LogicalVector nas = is_na(x);
  for(int i=1; i<tl; i++){
    if(nas[i]){
      out[i]=out[i-1];
    }
  }
  return out;
}

// [[Rcpp::export]]
List EnumerateFrom(List sequences, List times, String target){
  int ll = sequences.length();
  List out (ll);
  for(int i = 0; i < ll; ++i){
    CharacterVector currauth = sequences[i];
    int cl = currauth.length();
    IntegerVector currtimes = times[i];
    
    bool found = false;
    int flippoint = cl+1;
    for(int j = 0; j < cl; ++j){
      if(currauth[j] == target){
        flippoint = currtimes[j];
        found = true;
        break;
      }
    }
    
    NumericVector newtimes (cl);
    if(found){
      newtimes = currtimes - flippoint;
    }else{
      newtimes = rep(0, cl);
    }
    out[i] = newtimes;
  }
  return out;
}
