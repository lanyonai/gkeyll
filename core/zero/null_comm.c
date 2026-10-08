#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_comm_priv.h>
#include <gkyl_elem_type_priv.h>
#include <gkyl_null_comm.h>
#include <gkyl_null_comm_priv.h>

#include <math.h>
#include <string.h>

#include <gkyl_range.h>

// simple hash-map from longs to skin_ghost_ranges
typedef struct sgr_map_entry {
  long key;
  struct skin_ghost_ranges value;
} sgr_map_entry;

typedef struct sgr_map {
  sgr_map_entry *entries;
  unsigned char *used; /* 1 if slot occupied */
  size_t capacity;     /* always a power of two (or 0) */
  size_t count;
} sgr_map;

/* Returns 0 on success, -1 on allocation failure. */
static int sgr_map_init(sgr_map *m, size_t expected_count);
static void sgr_map_free(sgr_map *m);
static void sgr_map_clear(sgr_map *m);
/* Insert or overwrite. 0 on success, -1 on alloc failure. */
static int sgr_map_put(sgr_map *m, long key,
                       const struct skin_ghost_ranges *val);
/* Pointer to the stored value, or NULL if absent. The pointer is
 * invalidated by any subsequent insert/put/remove. */
static struct skin_ghost_ranges *sgr_map_get(const sgr_map *m, long key);
/* Returns 1 if removed, 0 if the key was not present. */
static int sgr_map_remove(sgr_map *m, long key);
// number of elements in map
static size_t sgr_map_size(const sgr_map *m);

#define SGR_MAP_MIN_CAP 16

/* splitmix64 finalizer: good avalanche for integer keys */
static size_t
sgr__hash(long key)
{
  unsigned long long x = (unsigned long long)key;
  x ^= x >> 30;
  x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27;
  x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return (size_t)x;
}

static size_t
sgr__round_pow2(size_t n)
{
  size_t c = SGR_MAP_MIN_CAP;
  while (c < n)
    c <<= 1;
  return c;
}

static int
sgr__alloc(sgr_map *m, size_t cap)
{
  m->entries = (sgr_map_entry *)gkyl_malloc(cap * sizeof *m->entries);
  m->used = (unsigned char *)calloc(cap, 1);
  if (!m->entries || !m->used) {
    gkyl_free(m->entries);
    gkyl_free(m->used);
    m->entries = NULL;
    m->used = NULL;
    return -1;
  }
  m->capacity = cap;
  m->count = 0;
  return 0;
}

/* Place into a table known to have room and not contain key. */
static void
sgr__place(sgr_map *m, long key, const struct skin_ghost_ranges *val)
{
  size_t mask = m->capacity - 1;
  size_t i = sgr__hash(key) & mask;
  while (m->used[i])
    i = (i + 1) & mask;
  m->used[i] = 1;
  m->entries[i].key = key;
  m->entries[i].value = *val;
  m->count++;
}

static int
sgr__grow(sgr_map *m, size_t new_cap)
{
  sgr_map old = *m;
  size_t i;
  if (sgr__alloc(m, new_cap) != 0) {
    *m = old;
    return -1;
  }
  for (i = 0; i < old.capacity; i++)
    if (old.used[i])
      sgr__place(m, old.entries[i].key, &old.entries[i].value);
  gkyl_free(old.entries);
  gkyl_free(old.used);
  return 0;
}

/* Keep load factor <= 0.75 */
static int
sgr__reserve_one(sgr_map *m)
{
  if (m->capacity == 0)
    return sgr__alloc(m, SGR_MAP_MIN_CAP);
  if ((m->count + 1) * 4 > m->capacity * 3)
    return sgr__grow(m, m->capacity * 2);
  return 0;
}

static size_t
sgr__find(const sgr_map *m, long key, int *found)
{
  size_t mask, i;
  *found = 0;
  if (m->capacity == 0)
    return 0;
  mask = m->capacity - 1;
  i = sgr__hash(key) & mask;
  while (m->used[i]) {
    if (m->entries[i].key == key) {
      *found = 1;
      return i;
    }
    i = (i + 1) & mask;
  }
  return i;
}

static int
sgr_map_init(sgr_map *m, size_t expected_count)
{
  m->entries = NULL;
  m->used = NULL;
  m->capacity = 0;
  m->count = 0;
  if (expected_count == 0)
    return 0;
  return sgr__alloc(m, sgr__round_pow2(expected_count + expected_count / 3 + 1));
}

static void
sgr_map_free(sgr_map *m)
{
  gkyl_free(m->entries);
  gkyl_free(m->used);
  m->entries = NULL;
  m->used = NULL;
  m->capacity = 0;
  m->count = 0;
}

void sgr_map_clear(sgr_map *m)
{
  if (m->used)
    memset(m->used, 0, m->capacity);
  m->count = 0;
}

static int
sgr_map_insert(sgr_map *m, long key, const struct skin_ghost_ranges *val)
{
  if (sgr__reserve_one(m) != 0)
    return -1;
  sgr__place(m, key, val);
  return 0;
}

int
sgr_map_put(sgr_map *m, long key, const struct skin_ghost_ranges *val)
{
  int found;
  size_t i = sgr__find(m, key, &found);
  if (found) {
    m->entries[i].value = *val;
    return 0;
  }
  return sgr_map_insert(m, key, val);
}

struct skin_ghost_ranges *sgr_map_get(const sgr_map *m, long key) {
  int found;
  size_t i = sgr__find(m, key, &found);
  return found ? &m->entries[i].value : NULL;
}

size_t sgr_map_size(const sgr_map *m) { return m->count; }

// Create ghost and skin sub-ranges given a parent range
static void skin_ghost_ranges_init(struct skin_ghost_ranges *sgr, const struct gkyl_range *parent,  const int *ghost) {
#define G_MAX(a, b) (a) > (b) ? (a) : (b)

  int ndim = parent->ndim;
  long max_vol = 0;

  for (int d = 0; d < ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d], d,
      GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d], d,
      GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

// Create ghost and skin sub-ranges given a parent range: includes
// corners
static void
skin_ghost_ranges_with_corners_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent,
  const int *ghost)
{
#define G_MAX(a, b) (a) > (b) ? (a) : (b)

  int ndim = parent->ndim;
  long max_vol = 0;

  for (int d = 0; d < ndim; ++d) {
    gkyl_skin_ghost_with_corners_ranges(&sgr->lower_skin[d],
      &sgr->lower_ghost[d], d,
      GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);

    gkyl_skin_ghost_with_corners_ranges(&sgr->upper_skin[d],
      &sgr->upper_ghost[d], d,
      GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);

  sgr_map_free(null_comm->l2sgr);
  gkyl_free(null_comm->l2sgr);
  
  sgr_map_free(null_comm->l2sgr_wc);
  gkyl_free(null_comm->l2sgr_wc);
  
  gkyl_rect_decomp_release(null_comm->decomp);
  gkyl_mem_buff_release(null_comm->pbuff);
  gkyl_free(null_comm);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = 0;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  *sz = 1;
  return 0;
}

static int
allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp,
  void *out)
{
  struct null_comm *null_comm =
      container_of(comm, struct null_comm, priv_comm.pub_comm);
  if (null_comm->use_gpu)
    gkyl_cu_memcpy(out, inp, gkyl_elem_type_size[type] * nelem,
      GKYL_CU_MEMCPY_D2D);
  else
    memcpy(out, inp, gkyl_elem_type_size[type] * nelem);
  return 0;
}

static int
allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);
  memcpy(out, inp, gkyl_elem_type_size[type] * nelem);
  return 0;
}

static int
array_allgather(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *global,
  const struct gkyl_array *array_local,
  struct gkyl_array *array_global)
{
  gkyl_array_copy(array_global, array_local);
  return 0;
}

static int
array_bcast(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  gkyl_array_copy(arecv, asend);
  return 0;
}


static int
array_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  return 0;
}

// apply periodic BCs
static void
apply_periodic_bc(const struct skin_ghost_ranges *sgr, char *data,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(data, f, &(sgr->lower_skin[dir]));
  gkyl_array_copy_from_buffer(f, data, &(sgr->upper_ghost[dir]));

  gkyl_array_copy_to_buffer(data, f, &(sgr->upper_skin[dir]));
  gkyl_array_copy_from_buffer(f, data, &(sgr->lower_ghost[dir]));
}

static int
array_per_no_corners_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);

  int nghost[GKYL_MAX_DIM];
  for (int d = 0; d < null_comm->decomp->ndim; ++d)
    nghost[d] = local_ext->upper[d] - local->upper[d];

  long lkey = gkyl_range_idx(&null_comm->grange, nghost);

  if (!sgr_map_get(null_comm->l2sgr, lkey)) {
    struct skin_ghost_ranges sgr;
    skin_ghost_ranges_init(&sgr, local_ext, nghost);
    sgr_map_put(null_comm->l2sgr, lkey, &sgr);
  }

  const struct skin_ghost_ranges *val = sgr_map_get(null_comm->l2sgr, lkey);
  long max_vol_esnz = val->max_vol * array->esznc;

  if (max_vol_esnz > gkyl_mem_buff_size(null_comm->pbuff))
    gkyl_mem_buff_resize(null_comm->pbuff, max_vol_esnz);

  char *data = gkyl_mem_buff_data(null_comm->pbuff);

  for (int d = 0; d < nper_dirs; ++d)
    apply_periodic_bc(val, data, per_dirs[d], array);

  return 0;
}

static int
array_per_with_corners_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);

  int nghost[GKYL_MAX_DIM];
  for (int d = 0; d < null_comm->decomp->ndim; ++d)
    nghost[d] = local_ext->upper[d] - local->upper[d];

  long lkey = gkyl_range_idx(&null_comm->grange, nghost);

  if (!sgr_map_get(null_comm->l2sgr_wc, lkey)) {
    struct skin_ghost_ranges sgr;
    skin_ghost_ranges_with_corners_init(&sgr, local_ext, nghost);
    sgr_map_put(null_comm->l2sgr_wc, lkey, &sgr);
  }

  const struct skin_ghost_ranges *val = sgr_map_get(null_comm->l2sgr_wc, lkey);
  long max_vol_esnz = val->max_vol * array->esznc;

  if (max_vol_esnz > gkyl_mem_buff_size(null_comm->pbuff))
    gkyl_mem_buff_resize(null_comm->pbuff, max_vol_esnz);

  char *data = gkyl_mem_buff_data(null_comm->pbuff);

  for (int d = 0; d < nper_dirs; ++d)
    apply_periodic_bc(val, data, per_dirs[d], array);

  return 0;
}

static int
array_per_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext, int nper_dirs,
  const int *per_dirs, struct gkyl_array *array)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);
  array_per_no_corners_sync(comm, local, local_ext, nper_dirs, per_dirs, array);
  if (null_comm->sync_corners)
    array_per_with_corners_sync(comm, local, local_ext, nper_dirs, per_dirs, array);

  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{ return 0; }

static int
array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  return gkyl_grid_sub_array_write(grid, range, meta, arr, fname);
}

static int
array_read(struct gkyl_comm *comm, const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, struct gkyl_array *arr,
  const char *fname)
{
  struct gkyl_rect_grid fgrid;
  int status = gkyl_grid_sub_array_read(&fgrid, range, arr, fname);
  if (status == 0) {
    if (!gkyl_rect_grid_cmp(grid, &fgrid))
      status = 1;
  }
  return status;
}

static struct gkyl_comm *
extend_comm(const struct gkyl_comm *comm,
  const struct gkyl_range *erange)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);
  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp =
    gkyl_rect_decomp_extended_new(erange, null_comm->decomp);
  struct gkyl_comm *ext_comm = gkyl_null_comm_inew(
    &(struct gkyl_null_comm_inp){
      .decomp = ext_decomp,
      .use_gpu = null_comm->use_gpu,
      .sync_corners = null_comm->sync_corners
    }
  );
  gkyl_rect_decomp_release(ext_decomp);

  return ext_comm;
}

static struct gkyl_comm *
split_comm(const struct gkyl_comm *comm, int color,
  struct gkyl_rect_decomp *new_decomp)
{
  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);

  return gkyl_null_comm_inew(
      &(struct gkyl_null_comm_inp){
        .use_gpu = null_comm->use_gpu,
        .sync_corners = null_comm->sync_corners,
        .decomp = new_decomp
      }
  );
}

static struct gkyl_comm *
create_comm_from_ranks(const struct gkyl_comm *comm, int nranks,
  const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  if (nranks > 1) {
    *is_valid = false;
    return 0;
  }

  *is_valid = true;

  struct null_comm *null_comm =
    container_of(comm, struct null_comm, priv_comm.pub_comm);
  return gkyl_null_comm_inew(
      &(struct gkyl_null_comm_inp){
        .use_gpu = null_comm->use_gpu,
        .sync_corners = null_comm->sync_corners,
        .decomp = new_decomp
      }
  );
}

struct gkyl_comm *
gkyl_null_comm_inew(const struct gkyl_null_comm_inp *inp)
{
  struct null_comm *comm = gkyl_malloc(sizeof *comm);
  strcpy(comm->priv_comm.pub_comm.id, "null_comm");

  comm->priv_comm.pub_comm.has_decomp = true;
  if (0 == inp->decomp) {
    comm->priv_comm.pub_comm.has_decomp = false;

    // construct a dummy decomposition
    comm->decomp =
        gkyl_rect_decomp_new_from_cuts_and_cells(1, (int[]){1}, (int[]){1});
  } else {
    comm->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  }

  // construct range to hash ghost layout
  int lower[GKYL_MAX_DIM] = {0};
  int upper[GKYL_MAX_DIM];
  for (int d = 0; d < comm->decomp->ndim; ++d)
    upper[d] = GKYL_MAX_NGHOST;
  gkyl_range_init(&comm->grange, comm->decomp->ndim, lower, upper);

  comm->use_gpu = inp->use_gpu;
  comm->sync_corners = inp->sync_corners;

  comm->l2sgr = gkyl_malloc(sizeof(sgr_map));
  sgr_map_init(comm->l2sgr, 0);
  comm->l2sgr_wc = gkyl_malloc(sizeof(sgr_map));
  sgr_map_init(comm->l2sgr_wc, 0);

  if (comm->use_gpu)
    comm->pbuff = gkyl_mem_buff_cu_new(1024); // will be reallocated
  else
    comm->pbuff = gkyl_mem_buff_new(1024); // will be reallocated

  comm->priv_comm.get_rank = get_rank;
  comm->priv_comm.get_size = get_size;
  comm->priv_comm.allreduce = allreduce;
  comm->priv_comm.allreduce_host = allreduce_host;
  comm->priv_comm.gkyl_array_allgather = array_allgather;
  comm->priv_comm.gkyl_array_allgather_host = array_allgather;
  comm->priv_comm.gkyl_array_bcast = array_bcast;
  comm->priv_comm.gkyl_array_bcast_host = array_bcast;
  comm->priv_comm.gkyl_array_sync = array_sync;
  comm->priv_comm.gkyl_array_per_sync = array_per_sync;
  comm->priv_comm.barrier = barrier;
  comm->priv_comm.gkyl_array_write = array_write;
  comm->priv_comm.gkyl_array_read = array_read;
  comm->priv_comm.extend_comm = extend_comm;
  comm->priv_comm.split_comm = split_comm;
  comm->priv_comm.create_comm_from_ranks = create_comm_from_ranks;

  comm->priv_comm.pub_comm.ref_count = gkyl_ref_count_init(comm_free);

  return &comm->priv_comm.pub_comm;
}
