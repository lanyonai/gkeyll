#include <acutest.h>

#include <gkyl_fvec.h>

void
test_null(void)
{
  int *a = 0;
  TEST_CHECK( 0 == gkyl_fvec_size(a) );
  TEST_CHECK( 0 == gkyl_fvec_capacity(a) );  
  gkyl_fvec_free(a);
}

void
test_int(void)
{
  int *a = 0;
  TEST_CHECK( 0 == gkyl_fvec_size(a) );
  
  for (int i=0; i<100; ++i)
    gkyl_fvec_push(a, i+10);

  TEST_ASSERT( gkyl_fvec_size(a) == 100 );

  for (int i=0; i<100; ++i)
    TEST_ASSERT( a[i] == (i+10) );

  gkyl_fvec_free(a);
}

TEST_LIST = {
  { "test_null", test_null },
  { "test_int", test_int },
  { NULL, NULL },
};
