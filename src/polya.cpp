#include <iostream>
#include <gmpxx.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

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

  // set the number of bits in our modulus using first argument of program
  int numbit = atoi(argv[1]);

  // generate prime of numbits
  NTL::ZZ p = NTL::GenPrime_ZZ(numbit);
  NTL::ZZ_p::init(p);

  // set our range using second and third argument of program
  NTL::ZZ m = NTL::ZZ(0);
  NTL::ZZ n = p -1;

  // optionally on command line
  // NTL::ZZ m = NTL::ZZ(atoi(argv[2]));
  // NTL::ZZ n = NTL::ZZ(atoi(argv[3]));

  for (NTL::ZZ i = NTL::ZZ(m + 1); i <= (n+m); ++i) {
    // std::cout << i << "\t" << NTL::Jacobi(i, p) << std::endl;
    // std::cout << NTL::Jacobi(i, p) << std::endl;
    unsigned int j = Legendre(i, p);
    if (j == (unsigned int)-1) std::cout << "0";
    else std::cout << "1";
  }
  return 0;
}
