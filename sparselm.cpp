#include "splm.h"
#include "sparselm.h"

int sparselm::sparselm_dercrs(void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
	void(*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
	double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
	int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
	return ::sparselm_dercrs(func, fjac, p, x, nvars, nconvars, nobs, Jnnz, JtJnnz, itmax, opts, info, adata); // CRS Jacobian
}

int sparselm::sparselm_derccs(
	void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
	void(*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
	double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
	int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
	return ::sparselm_derccs(func, fjac, p, x, nvars, nconvars, nobs, Jnnz, JtJnnz, itmax, opts, info, adata); // CCS Jacobian
}

int sparselm::sparselm_difcrs(
	void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
	void(*fjac)(double *p, struct splm_crsm *jac, int nvars, int nobs, void *adata),
	double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
	int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
	return ::sparselm_difcrs(func, fjac, p, x, nvars, nconvars, nobs, Jnnz, JtJnnz, itmax, opts, info, adata); // CRS Jacobian
}

int sparselm::sparselm_difccs(
	void(*func)(double *p, double *hx, int nvars, int nobs, void *adata),
	void(*fjac)(double *p, struct splm_ccsm *jac, int nvars, int nobs, void *adata),
	double *p, double *x, int nvars, const int nconvars, int nobs, int Jnnz, int JtJnnz,
	int itmax, double opts[SPLM_OPTS_SZ], double info[SPLM_INFO_SZ], void *adata)
{
	return ::sparselm_difccs(func, fjac, p, x, nvars, nconvars, nobs, Jnnz, JtJnnz, itmax, opts, info, adata); // CCS Jacobian
}

void sparselm::sparselm_chkjaccrs(
	void(*func)(double *p, double *hx, int m, int n, void *adata),
	void(*jacf)(double *p, struct splm_crsm *jac, int m, int n, void *adata),
	double *p, int m, int n, int jnnz, void *adata, double *err)
{
	::sparselm_chkjaccrs(func, jacf, p, m, n, jnnz, adata, err);
}

void sparselm::sparselm_chkjacccs(
	void(*func)(double *p, double *hx, int m, int n, void *adata),
	void(*jacf)(double *p, struct splm_ccsm *jac, int m, int n, void *adata),
	double *p, int m, int n, int jnnz, void *adata, double *err)
{
	::sparselm_chkjacccs(func, jacf, p, m, n, jnnz, adata, err);
}

/* CRS sparse matrices manipulation routines */
void sparselm::splm_crsm_alloc(struct splm_crsm *sm, int nr, int nc, int nnz)
{
	::splm_crsm_alloc(sm, nr, nc, nnz);
}

void sparselm::splm_crsm_alloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz)
{
	::splm_crsm_alloc_novalues(sm, nr, nc, nnz);
}

void sparselm::splm_crsm_alloc_values(struct splm_crsm *sm)
{
	::splm_crsm_alloc_values(sm);
}

void sparselm::splm_crsm_realloc_novalues(struct splm_crsm *sm, int nr, int nc, int nnz)
{
	::splm_crsm_realloc_novalues(sm, nr, nc, nnz);
}

void sparselm::splm_crsm_free(struct splm_crsm *sm)
{
	::splm_crsm_free(sm);
}

int sparselm::splm_crsm_elmidx(struct splm_crsm *sm, int i, int j)
{
	return ::splm_crsm_elmidx(sm, i, j);
}

int sparselm::splm_crsm_elmrow(struct splm_crsm *sm, int idx)
{
	return ::splm_crsm_elmrow(sm, idx);
}

int sparselm::splm_crsm_row_elmidxs(struct splm_crsm *sm, int i, int *vidxs, int *jidxs)
{
	return ::splm_crsm_row_elmidxs(sm, i, vidxs, jidxs);
}

int sparselm::splm_crsm_row_maxnelms(struct splm_crsm *sm)
{
	return ::splm_crsm_row_maxnelms(sm);
}

int sparselm::splm_crsm_col_elmidxs(struct splm_crsm *sm, int j, int *vidxs, int *iidxs)
{
	return ::splm_crsm_col_elmidxs(sm, j, vidxs, iidxs);
}

void sparselm::splm_crsm2ccsm(struct splm_crsm *crs, struct splm_ccsm *ccs)
{
	::splm_crsm2ccsm(crs, ccs);
}

void sparselm::splm_crsm_row_sort(struct splm_crsm *sm)
{
	::splm_crsm_row_sort(sm);
}

/* CCS sparse matrices manipulation routines */
void sparselm::splm_ccsm_alloc(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
	::splm_ccsm_alloc(sm, nr, nc, nnz);
}

void sparselm::splm_ccsm_alloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
	::splm_ccsm_alloc_novalues(sm, nr, nc, nnz);
}

void sparselm::splm_ccsm_alloc_values(struct splm_ccsm *sm)
{
	::splm_ccsm_alloc_values(sm);
}

void sparselm::splm_ccsm_realloc_novalues(struct splm_ccsm *sm, int nr, int nc, int nnz)
{
	::splm_ccsm_realloc_novalues(sm, nr, nc, nnz);
}

void sparselm::splm_ccsm_free(struct splm_ccsm *sm)
{
	::splm_ccsm_free(sm);
}

int sparselm::splm_ccsm_elmidx(struct splm_ccsm *sm, int i, int j)
{
	return ::splm_ccsm_elmidx(sm, i, j);
}

int sparselm::splm_crsm_elmcol(struct splm_ccsm *sm, int idx)
{
	return ::splm_crsm_elmcol(sm, idx);
}

int sparselm::splm_ccsm_row_elmidxs(struct splm_ccsm *sm, int i, int *vidxs, int *jidxs)
{
	return ::splm_ccsm_row_elmidxs(sm, i, vidxs, jidxs);
}

int sparselm::splm_ccsm_col_elmidxs(struct splm_ccsm *sm, int j, int *vidxs, int *iidxs)
{
	return ::splm_ccsm_col_elmidxs(sm, j, vidxs, iidxs);
}

int sparselm::splm_ccsm_col_maxnelms(struct splm_ccsm *sm)
{
	return ::splm_ccsm_col_maxnelms(sm);
}

void sparselm::splm_ccsm2crsm(struct splm_ccsm *ccs, struct splm_crsm *crs)
{
	::splm_ccsm2crsm(ccs, crs);
}

int sparselm::splm_ccsm_drop_cols(struct splm_ccsm *A, int ncols)
{
	return ::splm_ccsm_drop_cols(A, ncols);
}

void sparselm::splm_ccsm_restore_cols(struct splm_ccsm *A, int ncols, int ncnnz)
{
	::splm_ccsm_restore_cols(A, ncols, ncnnz);
}

void sparselm::splm_ccsm_col_sort(struct splm_ccsm *sm)
{
	::splm_ccsm_col_sort(sm);
}

double sparselm::splm_gettime(void)
{
	return ::splm_gettime();
}

void sparselm::splm_stm_alloc(struct splm_stm *sm, int nr, int nc, int maxnnz)
{
	::splm_stm_alloc(sm, nr, nc, maxnnz);
}

void sparselm::splm_stm_allocval(struct splm_stm *sm, int nr, int nc, int maxnnz)
{
	::splm_stm_allocval(sm, nr, nc, maxnnz);
}

void sparselm::splm_stm_free(struct splm_stm *sm)
{
	::splm_stm_free(sm);
}

int sparselm::splm_stm_nonzero(struct splm_stm *sm, int i, int j)
{
	return ::splm_stm_nonzero(sm, i, j);
}

int sparselm::splm_stm_nonzeroval(struct splm_stm *sm, int i, int j, double val)
{
	return ::splm_stm_nonzeroval(sm, i, j, val);
}

void sparselm::splm_stm2ccsm(struct splm_stm *st, struct splm_ccsm *ccs)
{
	::splm_stm2ccsm(st, ccs);
}

void sparselm::splm_tri2ccsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_ccsm *ccs)
{
	::splm_tri2ccsm(i, j, s, m, n, nzmax, ccs);
}

void sparselm::splm_stm2crsm(struct splm_stm *st, struct splm_crsm *crs)
{
	::splm_stm2crsm(st, crs);
}

void sparselm::splm_tri2crsm(int *i, int *j, double *s, int m, int n, int nzmax, struct splm_crsm *crs)
{
	::splm_tri2crsm(i, j, s, m, n, nzmax, crs);
}