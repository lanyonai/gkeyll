#pragma once

#include <gkyl_alloc.h>

struct gkyl_fvec_header { size_t len, capacity; };
#define _gkyl_fvec_header_(arr) ((struct gkyl_fvec_header*) (arr) -  1)
#define GKYL_FVEC_INIT_SIZE 256
#define _gkyl_fvec_resize_(hdr, arr)                                          \
    do {                                                                \
      if (hdr->len >= hdr->capacity) {                                  \
        size_t cap = hdr->capacity = 2*hdr->capacity;                   \
        hdr = gkyl_realloc(hdr, sizeof(struct gkyl_fvec_header) + cap*sizeof(*arr)); \
        arr = (void*) (hdr + 1);                                        \
      }                                                                 \
    } while(0)

/**
 * Push value @a val at the end of array @a arr.
 */
#define gkyl_fvec_push(arr, val)                                        \
    do {                                                                \
      if (0 == arr) {                                                   \
        size_t cap = GKYL_FVEC_INIT_SIZE;                               \
        struct gkyl_fvec_header *hdr = gkyl_malloc(sizeof(struct gkyl_fvec_header) + cap*sizeof(*arr)); \
        hdr->len = 0; hdr->capacity = cap;                              \
        arr = (void*) (hdr + 1);                                        \
      }                                                                 \
      struct gkyl_fvec_header *hdr = _gkyl_fvec_header_(arr);           \
      _gkyl_fvec_resize_(hdr, arr);                                     \
    (arr)[hdr->len++] = val;                                            \
    } while (0)

/**
 * Returns capacity of vector (amount of space allocated) in @a
 * arr. This is a multiple of GKYL_FVEC_INIT_SIZE
 */
#define gkyl_fvec_capacity(arr) (0 == arr ? 0 : _gkyl_fvec_header_(arr)->capacity)

/**
 * Returns number of elements in vector @a arr
 */
#define gkyl_fvec_size(arr) (0 == arr ? 0 :_gkyl_fvec_header_(arr)->len)

/**
 * Frees the array.
 */
#define gkyl_fvec_free(arr)                                             \
    do {                                                                \
      if (0 != arr) gkyl_free( _gkyl_fvec_header_(arr) );               \
    } while (0)
