

# @import methods
# @importFrom stats rnorm



#' X decomposition
#'
#' @param X
#' @param coef.d
#' @param coef.h
#'
#' @return
#' @export
#'
#' @examples
X.decomp <- function(X = NULL, coef.d = NULL, coef.h = setdiff(seq_len(ncol(X)), coef.d))
{
  Xq = qr.Q(qr(X), complete = TRUE)

  list(
  Xd = Xq[,coef.d, drop = FALSE],
  Xhe = Xq[,setdiff(seq_len(ncol(Xq)), coef.d), drop = FALSE]
  )
}




### pre-allocate memory before executing loops


# ### set following methods (see also limma's "classes.R" file):
# dim
# dimnames
# as.matrix
# "dimnames<-.RGList" <- .setdimnames

### ev. auch die [] operatoren definieren. Aber bringt das was ? Das gesamte Objekt
### müsste neu berechnet werden wenn die Anzahl an Spalten verändert wird.





#### implement contrasts and moderated statistics !





### Initialization can either be done with limma compatible objects containing $design, $E and $weights
### or by hand.
### We do not use annotation data here in order to avoid excessive copying of annotation data when the
### resampled data is executed (e.g. with lmFit or ComBat)

init.randrot <- function(Y = NULL, X = NULL, coef.h = NULL, weights = NULL)
{

  if(is.null(X)) X <- matrix(1, ncol(Y))
  if(any(dim(X) < 1)) stop("Dimensions of X (design matrix) invalid.")
  if(is.null(Y) || ncol(Y) != nrow(X)) stop("Number of rows in X (design matrix) and number of columns in Y do not match.")

  if(is.null(coef.h)) coef.h <- seq_len(ncol(X))
  if(!all(coef.h %in% seq_len(ncol(X)))) stop("coef.h not contained in X (design matrix).")
  if(length(coef.h) < 1) stop("length(coef.h) must be at least 1.")

  if(is.null(colnames(X))) colnames(X) <- seq_len(ncol(X))

  coef.h = sort(coef.h)
  coef.d = setdiff(seq_len(ncol(X)), coef.h)
  if(length(coef.d > 0) && (min(coef.h) < max(coef.d))){
    message("coef.d do not correspond to the last columns in the design matrix.\nColumns of design matrix are rearranged.")
    X = X[, c(coef.d, coef.h), drop = FALSE]
    coef.d = seq_len(length(coef.d))
    coef.h = seq_len(length(coef.h)) + length(coef.d)
  }

  if(!is.null(weights)) return(init.randrot.w(Y, X, coef.h, coef.d, weights))

  decomp = X.decomp(X, coef.d, coef.h)

  Yd = Y %*% decomp$Xd %*% t(decomp$Xd)
  dimnames(Yd) = dimnames(Y)

  Xhe.Y = t(decomp$Xhe) %*% t(Y)

  new("init.randrot", list(X = X, Xhe.Y = Xhe.Y, Yd = Yd, Xhe = decomp$Xhe,
                           coef.d = coef.d, coef.h = coef.h))
}

init.randrot.w = function(Y, X, coef.h, coef.d, weights)
{
  if(any(weights <= 0) | any(!is.finite(weights))) stop("Weights must be finite > 0.")
  w = sqrt(weights)


  #### Transform Y
  Y.w = w * Y


  decomp.list = apply(w, 1, function(w.i){
    X.decomp(w.i * X, coef.d, coef.h)
  })

  Yd = t(sapply(seq_len(nrow(Y.w)), function(i){
    Y.w[i,,drop = FALSE] %*% decomp.list[[i]]$Xd %*% t(decomp.list[[i]]$Xd)
  }))

  dimnames(Yd) = dimnames(Y)

  Xhe.Y.w = sapply(1:nrow(Y.w), function(i){
    Y.w[i,,drop = FALSE] %*% decomp.list[[i]]$Xhe
  })

  new("init.randrot.w", list(X = X, Xhe.Y.w = Xhe.Y.w, Yd = Yd, decomp.list = decomp.list,
                             coef.d = coef.d, coef.h = coef.h, w = w))
}







#' Initialised randRot class
#'
#' @return
#' @export
#' @rdname class.init.randrot
#'
#' @examples
setClass("init.randrot", contains = "list")



#' Class init.randrot.w
#'
#' @return
#' @export
#' @rdname class.init.randrot
#'
#' @examples
setClass("init.randrot.w", contains = c("init.randrot", "list"))



#' Initialisation of a random rotation Object
#'
#' Initialization of a linear model for subsequent generation of randomly rotated data (\code{\link[randRotation:randrot]{randrot}}) associated with the null hypothesis
#' \eqn{H_{0}: \beta_{coef.h} = 0}{H0: \beta_coef.h = 0}. Basics of rotation tests are found in \insertCite{Langsrud2005}{randRotation}.
#'
#'
#' @param Y a data matrix with \code{features x samples} dimensions. Missing values (\code{\link[base:NA]{NA}}) are allowed but lead to NAs for all samples of the respective features in the rotated dataset and should henced be avoided.
#' @param X the design matrix of the experiment with \code{samples x coefficients} dimensions. If no design matrix is specified, intercept only model is used (design matrix with one column where all elements are \code{1}).
#' @param coef.h single integer or vector of integers specifying the \code{H_0} coefficients. \code{coef.h} should correspond to the last columns in \code{X} (see Details).
#' @param weights numerical matrix of finite positive weights > 0. Dimensions must be equal to dimensions of \code{Y}.
#'
#' @rdname init.randrot
#' @return an initialised object containing all necessary information for generating randomly rotated data (\code{\link[randRotation:randrot]{randrot}}).
#' @author Peter Hettegger
#' @references \insertAllCited{}
#'
#' @details
#'
#' This function performs basic initial checks and preparatory calculations for random rotation data generation, see \insertCite{Langsrud2005}{randRotation}.
#' Nomenclature of variables is mainly as in \insertCite{Langsrud2005}{randRotation}. See also package vignette for application examples.
#'
#' \code{coef.h} specifies the model coefficients associated with the null hypothesis. All other model coefficients #### continue here
#' The design matrix is rearranged so that \code{coef.h} correspond to the last columns of the design matrix. This is necessary for adequate transformation of the combined null-hypothesis \eqn{H_{0}: \beta_{coef.h} = 0}{H0: \beta_coef.h = 0} by QR decomposition.
#'
#'
#' Weights must be finite positive values greater zero. This ist necessary for model (QR) decomposition and for backtransformation of the rotated data into the original variance structure, see also \code{\link[randRotation:randrot]{randrot}}.
#' Weights as estimated e.g. by voom \insertCite{Law2014}{randRotation} are suitable and can be used without further processing.
#' @export
#' @seealso \code{\link[randRotation:randrot]{randrot}}
#'
#' @examples
#'
setGeneric("init.randrot", function(Y = NULL, X = NULL, coef.h = NULL, weights = NULL) {standardGeneric("init.randrot")})



#' Title
#'
#' @param list This is a test
#' @param Y2 Another test
#'
#' @return
#' @rdname init.randrot
#' @export
#'
#' @examples
setMethod("init.randrot", "list",
          function(Y = NULL, X = Y$design, coef.h = NULL, weights = Y$weights){
            init.randrot(Y = Y$E, X = X, coef.h = coef.h, weights = weights)
          })





#' @param init.randrot
#' @export
#' @rdname init.randrot

setMethod("show", "init.randrot",
          function(object)
          {
            cat("Initialised random rotation object", if(!is.null(object$w))"(with weights)","\n\n")
            cat(dim(object$Yd)[1], "features   ", dim(object$Yd)[2], "samples\n")
            cat("Coefficients to test (coef.h):", colnames(object$X)[object$coef.h], "\n\n")

            cat("\nmodel matrix (X):\n")
            print(head(init.rand$X, n = 6))
            if(nrow(object$X) > 6)cat(".\n.\n",nrow(object$X)-6, "more rows\n")

          }
)




setGeneric("randrot", function(object, ...) standardGeneric("randrot"))





#' Random rotation of initialised object
#'
#' Perform random data rotation of a previously initialised object (\code{\link[randRotation:init.randrot]{init.randrot}}) associated with the null hypothesis \eqn{H_{0}: \beta_{coef.h} = 0}{H0: \beta_coef.h = 0}.
#'
#'
#' @param init.randrot
#'
#' @return a numerical matrix of rotated data under the specified combined null hypothesis.
#' @export
#' @rdname randrot
#'
#' @details
#'
#' This function generates a randomly rotated dataset from an initialised randrot object (\code{\link[randRotation:init.randrot]{init.randrot}}). See also package vignette for application examples.
#' Only the numerical matrix of rotated data is returned, no design matrix, weights or other info is return for efficiency purposes.
#' restricted random rotation matrix \eqn{R_t^*}{Rt*} \deqn{R_t^* = X_dX_d' + \left[X_h \quad X_e \right] R^* \left[X_h \quad X_e \right]^\prime}{Rt* = Xd Xd' + [Xh  Xe] R* [Xh  Xe]'}
#' with ' being the transposed and \code{R*} being a (reduced) random rotation matrix with reduced dimensions ### continue here
#'
#' In the case of weighted data, for each feature a separate QR decomposition and
#' restricted random orthogonal matrix is calculated with the same reduced random rotation matrix \eqn{R*} for all features.
#' Weighted is rotated by feature wise pre-multiplying \code{Y} and \code{X} with the respective weights.
#' @author Peter Hettegger
#' @examples


setMethod("randrot", "init.randrot",
          function(object){
            ### No excessive initial checks for efficiency purposes.
            R = randorth(nrow(object$Xhe.Y))

            with(object, Yd + t(Xhe %*% R %*% Xhe.Y))
          })


#' Title
#'
#' @param init.randrot.w
#'
#' @rdname randrot
setMethod("randrot", "init.randrot.w",
          function(object){
            ### No excessive initial checks for efficiency purposes.
            R = randorth(nrow(object$Xhe.Y.w))

            R.Xhe.Y.w = R %*% object$Xhe.Y.w

            n = length(object$decomp.list)
            Yhe = t(sapply(1:n, function(i){
              object$decomp.list[[i]]$Xhe %*% R.Xhe.Y.w[,i]
            }))

            # backtransform Y.w to Y with inverse weigths
            (1/object$w) * (object$Yd + Yhe)
          })






#' Random rotation matrix
#'
#' Generation of a random rotation matrix (random orthogonal matrix) distributed with haar measure.
#' @param n
#' @param type
#'
#' @return
#' @export
#' @details Adapted from pracma package (randortho function) (\code{\link[pracma:randortho]{randortho}})
#' \insertCite{Stewart1980a}{randRotation}
#'
#' @author Peter Hettegger
#' @examples
randorth <- function (n, type = c("orthonormal", "unitary"))
{
  ### this function was adapted from the pracma package (randortho function)

  stopifnot(is.numeric(n), length(n) == 1, floor(n) == ceiling(n),
            n >= 1)
  type <- match.arg(type)
  if (type == "orthonormal") {
    z <- matrix(rnorm(n^2), n)
  }
  else {
    z <- (matrix(rnorm(n^2), n) + (0+1i) * matrix(rnorm(n^2), n))
  }
  Z <- qr(z)
  q <- qr.Q(Z)
  r <- qr.R(Z)
  d <- diag(r)
  ph <- d/abs(d)
  q %*% diag(ph, length(ph))
}



#'
#' #' Generate random permutation matrix for n samples
#' #'
#' #' Generate a random permutation matrix for \code{n} samples.
#' #'
#' #' @param n Number of samples
#' #'
#' #' @return A random permutation matrix
#' #' @export
#' #' @details This methods generates an orthogonal matrix with entries with only one entry in each row and column being \code{1}, all other entries being \code{0}.
#' #' @examples
#' #' design = model.matrix(~Species, iris)
#' #' randpermut(5)
#' #' @author Peter Hettegger
#' randpermut <- function(n){
#'   m1 = matrix(0,n,n)
#'   m1[cbind(1:n,sample(n))] = 1
#'   m1
#' }



