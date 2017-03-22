#include <iostream>
#include <gmpxx.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

int main(int argc, char *argv[])
{
  mpz_t rand1,rand2;
  mpz_init (rand1);
  mpz_init (rand2);

  mpz_random (rand1,8);
  //mpz_random (rand2,512);

  mpz_nextprime ( rand2, rand1 );
  gmp_printf("random %Zd\n", rand2);
  //free the big ints
  mpz_clear(rand1);
  mpz_clear(rand2);

  return 0;
}
