#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]






// [[Rcpp::export]]
double eucld(const rowvec &x, const rowvec &y) {
  int n = x.size();
  double dist = 0;
  for (int i = 0; i < n; i++) {
     dist += pow(x(i)-y(i), 2.0);
  }
  dist = sqrt(dist);

  return dist;
}


// [[Rcpp::export]]
mat distMat(const mat &x) {
  int Mn = x.n_rows;
  int Mp = x.n_rows;
  int p = x.n_cols;
  rowvec a(p);
  rowvec b(p);
  mat M(Mn, Mn);
  for (int i = 0; i < Mn; i++) {
    for (int j = i+1; j < Mp+1; j++) {
      a = x.row(i);
      b = x.row(j-1);
      double d = eucld(a,b);
      M(j-1, i) = d;
      M(i, j-1) = d;
    }
  }

  return M;
}


// [[Rcpp::export]]
uvec sortcppvec(const vec &x) {
  uvec c = sort_index(x);

  return c;
}


// [[Rcpp::export]]
umat arma_rank(vec x) {
  umat d(x.n_rows,x.n_rows,fill::eye);
  umat res=index_max(d.rows(sort_index(x)));
  return(res);
}


// [[Rcpp::export]]
umat sortcpp(const mat &x) {
  umat M(x.n_rows, x.n_cols);
  uvec s(x.n_rows);
    for (uword i = 0; i < x.n_cols; i++) {
      M.col(i) = sort_index(x.col(i));
      // M.col(i) = s;
    }

  return M;
}


// [[Rcpp::export]]
bool contains(uvec X, uword z) {
  return std::find(X.begin(), X.end(), z) != X.end();
}





// [[Rcpp::export]]
uvec diff_nn(const umat &sX, const uvec &cl) {
  uvec s(sX.n_rows);
  uword i = 0;
  for (uword col = 0; col < sX.n_cols; col++) {

    i = 0;
    uvec v = sX.col(col);
    uvec l = find(cl == cl(col));


    while (contains(l, v(i))) {
      i++;
      l = find(cl == cl(col));
    }
    s(col) = i;
  }

  return s;
}





// [[Rcpp::export]]
arma::uvec num_same_class_neighbors(const arma::mat& data,
                                    const arma::uvec& labels,
                                    const arma::uword num_neighbors) {
  // Get the number of data points
  const arma::uword n = data.n_rows;

  // Initialize a vector to store the number of same-class neighbors for each point
  arma::uvec num_same_class(n);

  // Compute the number of same-class neighbors for each point
  for (arma::uword i = 0; i < n; i++) {
    // Get the label for the current point
    const arma::uword label = labels(i);

    // Compute the distance from the current point to all other points
    arma::rowvec point = data.row(i);
    arma::mat rep_point = arma::repmat(point, n, 1);
    arma::rowvec diff = data - rep_point;
    arma::rowvec sq_diff = diff % diff;
    arma::vec dist = arma::sum(sq_diff, 1);

    // Sort the distances in ascending order
    arma::uvec indices = arma::sort_index(dist);

    // Count the number of same-class neighbors
    arma::uword count = 0;
    for (arma::uword j = 1; j <= num_neighbors; j++) {
      if (labels(indices(j)) == label) {
        count++;
      }
    }

    // Store the number of same-class neighbors for the current point
    num_same_class(i) = count;
  }

  return num_same_class;
}








// [[Rcpp::export]]
umat is_in_radius(const mat &data, const mat &centroids, double r) {

  // initialization
  int n = centroids.n_rows, k = data.n_rows;
  umat outmat(n, k);
  mat distmat(n, k);

  // for each centroid: compute distance
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < n; j++) {
      rowvec x =  data.row(i) - centroids.row(j);
      double d = arma::sum(x % x);

      distmat(j,i) = d;
    }
  }

  // for (int i = 0; i < n; i++) {
  //   rowvec x = data.row(i);
  //   mat m = data.each_row() - x;
  //   vec d = arma::sqrt(arma::sum(m % m, 1));
  //   distmat.col(i) = d;
  // }



  return distmat <= r;
}


// [[Rcpp::export]]
rowvec most_frequent_row(const mat &x) {
  int n = x.n_rows;
  uvec counts(n);
  unsigned int count;
  bool is_equal = FALSE;

  for (int i = 0; i < n; i++) {
    count = 0;
    rowvec r = x.row(i);
    for (int j = 0; j < n; j++) {
      rowvec rX = x.row(j);
      is_equal = arma::approx_equal(r, rX, "absdiff", 1e-5);
      if (is_equal) count++;
    }
    counts(i) = count;
  }
  return x.row(index_max(counts));
}


// [[Rcpp::export]]
rowvec most_frequent_row2(const mat &x) {

  int n = x.n_rows, k = x.n_cols;
  uvec counts(n);
  unsigned int count;
  bool is_equal = FALSE, was_searched_before = FALSE;
  mat lookup = zeros(n, k);

  for (int i = 0; i < n; i++) {
    count = 0;
    rowvec r = x.row(i);

    // check if row was counted before; if yes, begin counting new row.
    // this is (much) faster if there are many duplicates in x.
    for (int l = 0; l < i; l++) {
      rowvec lX = lookup.row(l);
      was_searched_before = arma::approx_equal(r, lX, "absdiff", 1e-5);
      if (was_searched_before) break;
    }

    if (was_searched_before) continue;

    for (int j = 0; j < n; j++) {
      rowvec rX = x.row(j);
      is_equal = arma::approx_equal(r, rX, "absdiff", 1e-5);
      if (is_equal) count++;
    }
    counts(i) = count;
    lookup.row(i) = r;
  }
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("row") = x.row(index_max(counts)),
                                      Rcpp::Named("idx") = index_max(counts));

  return x.row(index_max(counts));
}


// [[Rcpp::export]]
mat mean_shift(const mat& data,
                      unsigned int s,
                      int maxiter = 50,
                      double etaX = 1e-5,
                      double r = 0.1) {

  // initialization
  unsigned int n = data.n_rows, p = data.n_cols, iter;
  umat in_radius_mat(s, n);
  bool converged = FALSE;
  double eta = 0;

  // choose s random data points
  uvec rand_ind = randperm(n, s);
  mat rand_subset = data.rows(rand_ind);

  // convergence check
  mat rand_subset_dupl = zeros(s, p);

  for (iter = 0; iter < maxiter; iter++) {
    rand_subset_dupl = rand_subset;
    in_radius_mat = is_in_radius(data, rand_subset, r);
    for (int l = 0; l < s; l++) {
      uvec in_radius_vec = arma::find(in_radius_mat.row(l)==1);
      rand_subset.row(l) = arma::mean(data.rows(in_radius_vec), 0);
    }
    eta = accu(pow(rand_subset - rand_subset_dupl, 2)) / (s*p);
    converged = eta < etaX;
    if (converged) break;
  }


  Rcpp::List out = Rcpp::List::create(Rcpp::Named("used_points") = in_radius_mat,
                                      Rcpp::Named("centroids0") = rand_subset_dupl,
                                      Rcpp::Named("centroids") = rand_subset,
                                      Rcpp::Named("iterations") = iter,
                                      Rcpp::Named("eta") = eta);
  return rand_subset;
}






// [[Rcpp::export]]
Rcpp::List nearest_diff_class(const arma::mat& data, const arma::vec& labels) {

  // init
  mat DATA = data;
  vec LABELS = labels;
  vec unique_labels = unique(labels);
  unsigned int c = unique_labels.size(), n = data.n_rows, p = data.n_cols;


  // start iteration loop here


  // Initialize result vector
  arma::uvec res(DATA.n_rows);
  arma::uvec nn(DATA.n_rows);


  // Loop over each data point
  for (int i = 0; i < DATA.n_rows; i++) {
    // Get current data point
    arma::rowvec x = DATA.row(i);

    // Calculate distance from current data point to all other data points
    mat m = data.each_row() - x;
    arma::vec d = arma::sqrt(arma::sum(m % m, 1));

    // Sort distances in ascending order and get the corresponding labels
    arma::uvec sorted_idx = arma::sort_index(d);
    nn(i) = sorted_idx(1);
    arma::vec sorted_labels = labels.elem(sorted_idx);

    // Loop over sorted labels and find the first label that is different from the current data point's label
    for (size_t j = 0; j < sorted_labels.n_elem; j++) {
      if (sorted_labels(j) != labels(i)) {
        res(i) = j;
        break;
      }
    }
  }

  arma::uvec nn_met_at_same = res == res(nn) && labels == labels(nn);
  res %= nn_met_at_same;
  arma::uvec uq = unique(res);
  umat histmat(uq.size(), c);


  unsigned int idx = 0;
  rowvec cent(p);
  umat rm_idx(n, c);
  mat new_cents(c, p);
  urowvec no_drop(c);


  // start for each class label
  for (int i = 0; i < c; i++) {
    uvec current_class = labels == unique_labels(i);
    uvec nnc = res(find(current_class));
    uvec uqc = unique(nnc);
    uvec histoc(uq.size());
    histoc = hist(nnc, uq);
    histmat.col(i) = histoc;

    unsigned int max_count = max(histoc.subvec(1, uq.size()-1));
    bool too_few = max_count < 2;
    no_drop(i) = too_few;
    if (too_few) continue;

    uvec sums = uq % histoc.replace(1, 0);
    idx = index_max(sums);

    unsigned int maxval = uq(idx);
    uvec rr = res == maxval && current_class;
    uvec rr_s = find(rr==1);
    mat ms = mean_shift(data.rows(rr_s), rr_s.size());
    cent = most_frequent_row2(ms);

    new_cents.row(i) = cent;
    rm_idx.col(i) = rr;
  }

  // Remove used rows & labels
  DATA.shed_rows(find( arma::sum(rm_idx, 1) == 1 ));
  LABELS.shed_rows(find( arma::sum(rm_idx, 1) == 1 ));
  unsigned int number_new_rows = c - sum(no_drop);
  n = DATA.n_rows + number_new_rows;
  DATA.resize(n, p);
  LABELS.resize(n);

  // Add new cents & new labels
  unsigned int inc = 0;
  for (int i = 0; i < c; i++) {
    if (no_drop(i) == 1) continue;
    else {
      DATA.row(n-number_new_rows+inc) = new_cents.row(i);
      LABELS.row(n-number_new_rows+inc) = i+1;
      inc++;
    }
  }

  // incorporate DATA and add new labels.


  Rcpp::List out = Rcpp::List::create(Rcpp::Named("first_diff_class") = res,
                                      Rcpp::Named("nn_idx") = nn,
                                      Rcpp::Named("nn_met_at_same") = nn_met_at_same,
                                      Rcpp::Named("uniques") = uq,
                                      Rcpp::Named("hist") = histmat,
                                      Rcpp::Named("drop") = no_drop,
                                      Rcpp::Named("cents") = new_cents,
                                      Rcpp::Named("rm_idx") = rm_idx,
                                      Rcpp::Named("data") = DATA,
                                      Rcpp::Named("labels") = LABELS);
  return out;
}




