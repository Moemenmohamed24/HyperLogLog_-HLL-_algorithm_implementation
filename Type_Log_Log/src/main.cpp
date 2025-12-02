#include<iostream>
#include<cmath>
#include<vector>
#include<cstdint>
#include<functional>

using namespace std;

template<typename T>
class HyperLogLog
{
private:
  int _inetial_zeros;
  int Number_ofbuckets;
  uint64_t HashNumber;
  int bucket;    
  vector<int>registers;
  
  int Make_Numofbuckets(int numofzeros)
  {
    return 1 << numofzeros;
  }
  
public:
    HyperLogLog(const int inetial_zeros)
    {
      _inetial_zeros = inetial_zeros;
      Number_ofbuckets = Make_Numofbuckets(_inetial_zeros);
      registers.resize(Number_ofbuckets, 0);
    }
    
    void Add_Element(const T& value)
    {
      string _value = to_string(value);
      HashNumber = static_cast<uint64_t>(hash<string>{}(_value));      
      
      
      bucket = HashNumber >> (64 - _inetial_zeros);
      
      
      uint64_t suffix = HashNumber << _inetial_zeros;
      
      
      int count = 1;
      
          uint64_t mask = 1ULL << 63;
          while ((suffix & mask) == 0) {
              count++;
              suffix <<= 1;  
              if (suffix == 0) break;
          }
      
      
      registers[bucket] = max(registers[bucket], count);
    }
    
    double ComputeCardinality()
    {
      int m = registers.size();
      double summingnum = 0.0;
      int zero_count = 0;
      
      for(int i = 0; i < m; i++)
      {
        summingnum += pow(2.0, -registers[i]);
        if (registers[i] == 0) zero_count++;
      }
      
      
      double alpha = 0.7213 / (1.0 + 1.079 / m);
      double estimate = (alpha * m * m) / summingnum;
      
      
      if (estimate <= 2.5 * m && zero_count > 0) {
          estimate = m * log((double)m / zero_count);
      }
      
      return estimate;
    }  
    
    double GetCardinality()
    {
      return ComputeCardinality();
    } 
};


template<typename T>
class HyperLogLogPresto
{
private:
  int _inetial_zeros;
  int Number_ofbuckets;
  uint64_t HashNumber;
  int bucket;    
  vector<int>registers;
  vector<int>overflow_bucket;

  int Make_Numofbuckets(int numofzeros)
  {
    return 1 << numofzeros;
  }
  
public:
  HyperLogLogPresto(const int inetial_zeros)
  {
    _inetial_zeros = inetial_zeros;
    Number_ofbuckets = Make_Numofbuckets(_inetial_zeros);
    registers.resize(Number_ofbuckets, 0);
    overflow_bucket.resize(Number_ofbuckets, 0);
  } 
  
  void Add_Element(const T& value)
  {
    string _value = to_string(value);
    HashNumber = static_cast<uint64_t>(hash<string>{}(_value));
    
    
    bucket = HashNumber >> (64 - _inetial_zeros);
    
    
    
    int rho;    
    rho = __builtin_ctzll(HashNumber) + 1;    
    
    
    int dense_bits = _inetial_zeros;
    int dense_value = rho > ((1 << _inetial_zeros) - 1 ) ? (1 << _inetial_zeros) - 1 : rho ;
    int overflow_value = rho > ((1 << _inetial_zeros) - 1 ) ? (rho - (1 << _inetial_zeros))  : 0 ;
    
    registers[bucket] = max(registers[bucket], dense_value);
    overflow_bucket[bucket] = max(overflow_bucket[bucket], overflow_value);
  }
  
  double ComputeCardinality()
{
    int m = registers.size();     
    double summingnum = 0.0;      
    int zero_count = 0;
    
    for(int i = 0; i < m; i++)
    {
        
        int rho = registers[i] + (overflow_bucket[i] << _inetial_zeros);
        
        if (rho == 0) {
            zero_count++;
        }
        
        
        summingnum += pow(2.0, -rho);
    }
    
    
    double alpha = 0.7213 / (1.0 + 1.079 / m);
    double estimate = (alpha * m * m) / summingnum;
    
    
    if (estimate <= 2.5 * m) {
        if (zero_count > 0) {
            estimate = m * log((double)m / zero_count);
        }
    }
    
    
    return floor(estimate);
}
  
  double GetCardinality()
  {
    return ComputeCardinality();
  }
};

int main()
{
    HyperLogLog<int> H1(8);
    for (int i = 1; i <= 20; i++) {
        H1.Add_Element(i);
    }
    cout << "20 elements: " << H1.GetCardinality() << endl;
    
    HyperLogLog<int> H2(8);
    for (int i = 1; i <= 100; i++) {
        H2.Add_Element(i);
    }
    cout << "100 elements: " << H2.GetCardinality() << endl;
    
    
    HyperLogLog<int> H3(8);
    for (int i = 1; i <= 10000; i++) {
        H3.Add_Element(i);
    }
    cout << "1000 elements: " << H3.GetCardinality() << endl;
    
    return 0;
}