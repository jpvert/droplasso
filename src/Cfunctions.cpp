#include <Rcpp.h>

using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int randWrapper(const int n) { return floor(unif_rand()*n);}



// SGD  drop_lasso mini batch
// [[Rcpp::export]]
NumericVector droplassoC(const NumericMatrix& x, const NumericVector& y, 
                         const std::string& family, const double& keep_prob, 
                         const double& lambda, const NumericVector& binit, 
                         const double& gamma0, const double& decay,
                         const unsigned& n_passes, 
                         const unsigned& minibatch_size ) {
  unsigned n_samples = x.nrow();
  unsigned n_features = x.ncol();
  unsigned t = 1, i, j, i_batch, i_pass;
  NumericVector b = clone(binit); // the variable to optimize
  NumericVector xi = no_init(n_features); // a dropout-ed sample
  NumericVector grad_drop = no_init(n_features); // gradient
  double gamma = 0; // learning rate
  double xb;

  // Create list of samples, to be shuffled at each pass over the training set
  std::vector<unsigned> idx = std::vector<unsigned>(n_samples);
  for (i = 0; i < n_samples; i++) {
    idx[i] = i;
  }
  
  // Main loop: passes over the training set
  for (i_pass = 0; i_pass < n_passes; ++i_pass) {

    // Shuffle training set
    std::random_shuffle(idx.begin(), idx.end(), randWrapper);

    // Pass over the batches 1
    for (i_batch=0; i_batch < (n_samples/minibatch_size) ; ++i_batch)
    {
      
      // Initialize gradient for the batch
      std::fill(grad_drop.begin(), grad_drop.end(), 0);
      
      // Loop over the samples in the batch
      for (i=(i_batch*minibatch_size); i<(i_batch*minibatch_size+minibatch_size) ; ++i){
        
        xb = 0; // initialize inner product between current b and dropout-ed sample
        for (j=0; j<n_features; j++) {
          // create dropout sample from x(idx[j],_)
          if (R::runif(0,1) < keep_prob) {
            xi[j] = x(idx[i],j) / keep_prob; // non dropout-ed feature
            xb += xi[j]*b[j]; // update inner product is the feature is not dropout-ed
          }
          else
            xi[j] = 0; // dropout-ed feature
        }
        
        // Increment gradient of the batch with the gradient of the current sample
        if (family=="gaussian") {
          grad_drop += ( xb - y[idx[i]] ) * xi;
          }
        else if (family=="binomial") {
          grad_drop += ( - y[idx[i]] + 1/(1+exp(-xb)) ) * xi ;
          }
        else {
          std::cout << "family is invalid" << std::endl;
          std::exit( EXIT_FAILURE );
        }
      }
      // Don't forget to divide by minibatch_size to account for the batch size
      grad_drop = grad_drop / minibatch_size;

      // SGD update
      gamma = gamma0 / (1+ decay * t); // update learning rate
      for (j=0; j< n_features ; ++j){
        if ( b[j] > gamma * ( grad_drop[j] + lambda ) ) {
          b[j] -= gamma * ( grad_drop[j] + lambda );
        }
        else if ( b[j] < gamma * ( grad_drop[j] - lambda ) ) {
          b[j] -= gamma * ( grad_drop[j] - lambda );
        }
        else {
          b[j] = 0;
        }
      }
    }
    t++;
  }
  return(b);
}



