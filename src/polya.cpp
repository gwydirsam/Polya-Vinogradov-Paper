#include <iostream>
#include <cstdio>
#include <gmpxx.h>
// #include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

#define PI 3.14159265358

// From NTL
int Legendre(const NTL::ZZ& aa, const NTL::ZZ& nn)
{
  NTL::ZZ a, n;
  long t, k;
  long d;

  a = aa;
  n = nn;
  t = 1;

  while (a != 0) {
    k = MakeOdd(a);
    d = trunc_long(n, 3);
    if ((k & 1) && (d == 3 || d == 5)) t = -t;

    if (trunc_long(a, 2) == 3 && (d & 3) == 3) t = -t;
    swap(a, n);
    rem(a, a, n);
  }

  if (n == 1)
    return t;
  else
    return 0;
}

int main(int argc, char *argv[])
{
  // NTL::SetNumThreads(4);

  // set the number of bits in our modulus using first argument of program
  int numbit = atoi(argv[1]);

  // generate prime of numbits
  NTL::ZZ p = NTL::GenPrime_ZZ(numbit);
  NTL::ZZ_p::init(p);

  // set our range to 512 numbers
  // NTL::ZZ n = NTL::ZZ(32768);
  // NTL::ZZ n = p-1;
  // Generate a random starting point 
  // std::srand(std::time(0)); // use current time as seed for random generator
  // NTL::ZZ m = (NTL::RandomBits_ZZ(std::rand())) % p; // use standard random function as seed for "big int" random number
  // NTL::ZZ m = NTL::ZZ(0);

  // optionally on command line
  NTL::ZZ m = NTL::ZZ(atoi(argv[2]));
  NTL::ZZ n = NTL::ZZ(atoi(argv[3]));

  std::cerr << "Modulus: "<< p << std::endl;
  std::cerr << "Length of Interval (N): "<< n << std::endl;
  std::cerr << "Start of Interval (M): "<< m << std::endl;
  std::cerr << "Length of char: "<< sizeof(char) << std::endl;

  unsigned char c = 0;
  int j;
  unsigned int k = 0;
  NTL::ZZ sum = NTL::ZZ(0);

  for (NTL::ZZ i = NTL::ZZ(m + 1); i <= (n+m); ++i) {
    // // std::cout << i << "\t" << NTL::Jacobi(i, p) << std::endl;
    // // std::cout << NTL::Jacobi(i, p) << std::endl;
    // int j = Legendre(i, p);
    // // if (j == -1) std::cout << '0';
    // // else std::cout << '1';

    // if (j == -1) j = 0;
    // std::cout << j;

    // calculate legendre
    j = Legendre(i, p);
    sum += j;
    // if you get 1, add 1
    // else j is -1. Add nothing, shifting in a 0 as our interpretation of -1
    if (j == 1) {
      c |= 1;
    }

    // if our char if full, print it, else shift it for the next bit
    if (k % 8 == 0) {
      std::cout << c;
      c = 0;
    } else {
      c <<= 1;
    }
    k++;

  }
  std::cerr << "Sum: " << NTL::abs(sum) << std::endl;
  // print Polya-Vinogradov Bound
  std::cerr << "Bounds" << std::endl;

  std::cerr << "Polya-Vinogradov " << "\t" << n/2+NTL::SqrRoot(p)*NTL::log(p) << std::endl;
  std::cerr << "Polya-Vinogradov (lower bound)" << "\t" << NTL::SqrRoot(p) << std::endl;
  std::cerr << "Polya-Vinogradov (upper bound)" << "\t" << NTL::SqrRoot(p)*NTL::log(p) << std::endl;
  std::cerr << "Montgomery-Vaughan" << "\t" << NTL::SqrRoot(p)*log(NTL::log(p)) << std::endl;
  std::cerr << "Tao w/ Burgess Trick" << "\t" << NTL::SqrRoot(p)*NTL::log(p)*NTL::log(p) << std::endl;
  std::cerr << "Schur (1918)" << "\t" << NTL::SqrRoot(p)/(2*PI) << std::endl;
  return 0;
}
