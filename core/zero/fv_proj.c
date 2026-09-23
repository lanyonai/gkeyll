#include <math.h>
#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_fv_proj.h>
#include <gkyl_gauss_quad_data.h>
#include <gkyl_quad_type.h>


GKYL_CU_DH
void static inline
fv_eval_1d_ser_p0(const double *z, double *b )
{
  b[0] = 7.0710678118654757e-01;
}

GKYL_CU_DH
void static inline
fv_eval_2d_ser_p0(const double *z, double *b )
{
  b[0] = 5.0000000000000000e-01;
}

GKYL_CU_DH
void static inline
fv_eval_3d_ser_p0(const double *z, double *b )
{
  b[0] = 3.5355339059327379e-01;
}

GKYL_CU_DH
void static inline
fv_eval_4d_ser_p0(const double *z, double *b )
{
  b[0] = 2.5000000000000000e-01;
}

GKYL_CU_DH
void static inline
fv_eval_5d_ser_p0(const double *z, double *b )
{
  b[0] = 1.7677669529663689e-01;
}

GKYL_CU_DH
void static inline
fv_eval_6d_ser_p0(const double *z, double *b )
{
  b[0] = 1.2500000000000000e-01;
}

// Types for various kernels.
typedef void (*fv_proj_on_basis_c2p_t)(const double *xcomp, double *xphys, void *ctx);

// input packaged as a struct
struct gkyl_fv_proj_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project

  enum gkyl_quad_type qtype; // quadrature to use
  
  int num_quad; // number of quadrature points
  int num_ret_vals; // number of return values in eval function
  evalf_t eval; // function to project
  void *ctx; // function context

  fv_proj_on_basis_c2p_t c2p_func; // Function that transforms a set of ndim
                                // computational coordinates to physical ones.
  void *c2p_func_ctx; // Context for c2p_func.
};

struct gkyl_fv_proj {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context

  int num_basis; // number of basis functions
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates

  fv_proj_on_basis_c2p_t c2p; // Function transformin comp to phys coords.
  void *c2p_ctx; // Context for the c2p mapping.
};

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_rect_grid *grid = ctx;
  int ndim = grid->ndim;
  for (int d=0; d<ndim; d++) xphys[d] = xcomp[d];
}


static struct gkyl_fv_proj*
gkyl_fv_proj_priv_inew(const struct gkyl_fv_proj_inp *inp)
{
  struct gkyl_fv_proj *up = gkyl_malloc(sizeof(struct gkyl_fv_proj));

  up->grid = *inp->grid;
  int num_quad = up->num_quad = inp->num_quad == 0 ? 1 : inp->num_quad;
  int num_ret_vals = up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;
  up->num_basis = 1;

  if (inp->c2p_func == 0) {
    up->c2p = c2p_identity;
    up->c2p_ctx = &up->grid; // Use grid as the context since all we need is ndim.
  }
  else {
    up->c2p = inp->c2p_func;
    up->c2p_ctx = inp->c2p_func_ctx;
  }

  double ordinates1[num_quad], weights1[num_quad];

  if (inp->qtype == GKYL_GAUSS_QUAD) {
    if (num_quad <= gkyl_gauss_max) {
      // use pre-computed values if possible (these are more accurate
      // than computing them on the fly)
      memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
      memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
    }
    else {
      gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
    }
  }
  else if (inp->qtype == GKYL_GAUSS_LOBATTO_QUAD) {
    assert( (num_quad > 1) && (num_quad <= gkyl_gauss_max) );
    
    // Gauss-Lobatto quadrature
    memcpy(ordinates1, gkyl_gauss_lobatto_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_lobatto_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    fprintf(stderr, "Quadrature type not available. Exiting... \n");
    assert(false);
  }

  // create range to loop over quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<inp->grid->ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, inp->grid->ndim, qshape);

  int tot_quad = up->tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature  
  up->ordinates = gkyl_array_new(GKYL_DOUBLE, inp->grid->ndim, tot_quad);
  up->weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    long node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(up->ordinates, node);
    for (int i=0; i<inp->grid->ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(up->weights, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  void (*eval_func[])(const double *, double *) = {
    0,
    fv_eval_1d_ser_p0,
    fv_eval_2d_ser_p0,
    fv_eval_3d_ser_p0,
    fv_eval_4d_ser_p0,
    fv_eval_5d_ser_p0,
    fv_eval_6d_ser_p0
  };

  // pre-compute basis functions at ordinates
  up->basis_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    eval_func[inp->grid->ndim](gkyl_array_fetch(up->ordinates, n),
      gkyl_array_fetch(up->basis_at_ords, n));

  return up;
}

static struct gkyl_fv_proj*
gkyl_fv_proj_priv_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  return gkyl_fv_proj_priv_inew( &(struct gkyl_fv_proj_inp) {
      .grid = grid,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = num_quad,
      .num_ret_vals = num_ret_vals,
      .eval = eval,
      .ctx = ctx,
      .c2p_func = 0,
      .c2p_func_ctx = NULL,
    }
  );
}


static inline void
log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  // Convert logical to computational coordinates.
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static void
gkyl_fv_proj_priv_quad(const struct gkyl_fv_proj *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;
  int num_ret_vals = up->num_ret_vals;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  // arrangement of f is as:
  // c0[0], c0[1], ... c1[0], c1[1], ....
  // where c0, c1, ... are components of f (num_ret_vals)
  int offset = 0;
  for (int n=0; n<num_ret_vals; ++n) {
    for (int k=0; k<num_basis; ++k) f[offset+k] = 0.0;

    for (int imu=0; imu<tot_quad; ++imu) {
      double tmp = weights[imu]*func_at_ords[n+num_ret_vals*imu];
      for (int k=0; k<num_basis; ++k)
        f[offset+k] += tmp*basis_at_ords[k+num_basis*imu];
    }
    offset += num_basis;
  }
}

static void
gkyl_fv_proj_priv_advance(const struct gkyl_fv_proj *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int tot_quad = up->tot_quad;
  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, tot_quad);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i=0; i<tot_quad; ++i) {
      log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->ordinates, i),
        up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    gkyl_fv_proj_priv_quad(up, fun_at_ords, gkyl_array_fetch(arr, lidx));
  }

  gkyl_array_release(fun_at_ords);
}

static void
gkyl_fv_proj_priv_release(struct gkyl_fv_proj* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_free(up);
}

// functions below are not static
gkyl_fv_proj*
gkyl_fv_proj_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  return gkyl_fv_proj_priv_new(grid, num_quad, num_ret_vals, eval, ctx);
}

void
gkyl_fv_proj_advance(const gkyl_fv_proj *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  gkyl_fv_proj_priv_advance(pob, tm, update_rng, out);

  // from projections, compute cell average
  double denorm = 1.0/sqrt(pow(2, update_rng->ndim));
  gkyl_array_scale_range(out, denorm, update_rng);
}

void
gkyl_fv_proj_release(gkyl_fv_proj* pob)
{
  gkyl_fv_proj_priv_release(pob);
}
