\name{filters}
\alias{linearWindow}
\alias{squareWindow}
\alias{gaussianWindow}
\alias{parzenWindow}
\alias{welchWindow}
\alias{sineWindow}
\alias{cosineSumWindow}
\alias{tukeyWindow}
\alias{plankWindow}
\alias{exponentialWindow}
\alias{HannConstants}
\alias{HammingConstants}
\alias{BlackmanConsts}
\alias{BlackmanExact}
\alias{NuttallConstants}
\alias{BlackmanNuttallConstants}
\alias{BlackmanHarrisConstants}
\alias{FlattopConstants}
\title{Creates a window function.}
\description{

  These functions create a filter with the given window and parameters.
  The filter is a function with arguments \code{(t,tstar)}
  which provides a weight between 0 and 1.  It can be used as the value of
  \code{\link{IRTF}$wfun}.


}
\usage{
squareWindow(window = 1, ..., symmetric = TRUE)
gaussianWindow(window = 1, sigma=.5, power=1, ..., symmetric = TRUE)
linearWindow(window=1, slope=1, ..., symmetric = TRUE)
parzenWindow(window=1,...,symmetric=TRUE)
welchWindow(window=1,slope=1,...,symmetric=TRUE)
sineWindow(window=1,power=1,...,symmetric=TRUE)
BlackmanConsts(alpha=.16)
cosineSumWindow(window=1,
                a=list(HannConstants,HammingConstants,
                       BlackmanConsts(0.16),BlackmanExact,
                       NuttallConstants,BlackmanNuttallConstants,
                       BlackmanHarrisConstants,FlattopConstants),
                ...,symmetric=TRUE)
tukeyWindow(window=1,alpha=1/2,...,symmetric=TRUE)
plankWindow(window=1,epsilon=1/4,...,symmetric=TRUE)
exponentialWindow(window=1,lambda=.8,...,symmetric=TRUE)
}
\arguments{
  \item{window}{Times outside of the window will have weight 0.}
  \item{\dots}{Placeholder for tuning parameters for other filters.}
  \item{sigma}{Tuning parameter for Gaussian window.}
  \item{power}{Used in the Gaussian and sine window to raise the value
    to a power (which makes it closer to a square window).  For the
    Gaussian window, this should be an integer.}
  \item{slope}{Used in the linear and Welch windows to add a slope (note
  this may cause the window to not hit zero at the edges.}
  \item{a}{A vector of coefficients for the cosine sum rule.  See details.}
  \item{alpha}{Used in the \code{BlackmanConst()} function to determine
    the value of the constant and in the Tukey window to indicate the
    point at which to switch from square to cosine.}
  \item{epsilon}{Used in the plank window to indicate the place at which
    the window starts to move from 1.}
  \item{lambda}{Used in the exponential filter to indicate the amount by
    which observations should be discounted per time unit.  Should be
    less than one.}
  \item{symmetric}{A logical value.  If false, the the weight of observations
    with \code{t > tstar} will be zero.}
}
\details{

  These functions create a window function which has signature \code{(t,tstar)}
  where \code{tstar} is the center or right hand boundary of the window.
  The functions all have a signature \code{(window=1,...,symmetric=TRUE)},
  where \code{window} is the size of the window and \code{...} is any other
  parameters which vary by the window type.  If \code{symmetric} is true,
  then the window will be symmetrical about \code{tstar}, and only the left half
  window if false.

  Let \eqn{w(t,t_*)} be a function of the observation time \eqn{t} and the
  target time \eqn{t_*}.  Let \eqn{\Delta t = t_{*} - t}, and let \eqn{L} be the
  window width.  Expect where noted, if \eqn{|\Delta t| > L}, then \eqn{w(t,t_*)=0}.
  If \code{symmetric=FALSE}, then if \eqn{\Delta t < 0} then \eqn{w(t,t_*)=0}.

  \emph{Square Window}.  This filter is one inside the window and zero outside.

  \emph{Linear Window}. This filter goes downward from $t_*$ at slope \code{slope}.

  \emph{Gaussian Window}.  When $p=2$ this is a Gaussian distribution, but can be sharper if
   if $p$ is greater.  \eqn{w(\Delta t) = \exp \left ( -\frac{1}{2} \left ( \frac{\Delta t}{\sigma L}
  \right )^p \right), } where \eqn{\sigma \leq 0.5}.   When \eqn{p=2}, then
  the standard deviation is \eqn{\sigma L}.

  \emph{Exponential Window}.  This filter ignores the window and decays
  exponentially, \eqn{\lambda^{|\Delta t|}}.  Note that \code{lambda} should be
  less than one.

  \emph{Parzan Window}.

  \deqn{w(\Delta t) = \begin{cases} 1 - 6 \left( \frac{\Delta t}{L} \right )^2
  \left ( 1 - \frac{|\Delta t|}{L} \right )\ , \tab 0 \le |\Delta t| \le L/2\cr
  2 \left (1 - \frac{|\Delta t|}{L} \right )^3 \tab L/2 \le |\Delta t| \le L\cr
  \end{cases}}

  \emph{Welch Window}.

  \deqn{w(\Delta t) = 1 - \left ( \frac{\Delta t}{L} \right )^2}

  \emph{Sine Window}. This is the sine function raised to the \code{power}, \eqn{k}.
  (Putting 0 at the center of the window shifts the phase so that sine becomes
  cosine).  When \code{power=2} this is the Hann window.

  \deqn{w(\Delta t) = \cos^k \left ( \frac{\pi \Delta t}{L} \right )}


  \emph{Cosine Sum Window}.  This is a sum of cosines (sines after the phase shift)
  at different frequencies.  The parameter, \code{a} is a vector of coefficients.

  \deqn{w(\Delta t) = \sum_k (-1)^k a_k \sin\left(\frac{k\pi \Delta t}{2L}\right).}

  There are a number of filters that are cosine window filters with
  specific parameters.
  \describe{
    \item{\emph{Hann Window}}{\code{HannConstants=c(0.5,0.5)}.}
    \item{\emph{Hamming Window}}{\code{HammingConstants=c(0.5434783,0.4565217)}.}
    \item{\emph{Blackman Window}}{The function \code{BlackmanConsts(alpha)}
      has coefficiencts \eqn{a_0=(1-\alpha)/2},
      \eqn{a_1 = 1/2}, and \eqn{a_2 = \alpha/2}
      (default \eqn{\alpha=0.16}).  Or use \code{BlackmanExact}.}
    \item{\emph{Nuttall}}{\code{NuttallConstants=c(0.355768, 0.487396, 0.144232, 0.012604)}.}
    \item{\emph{Blackman--Nuttall}}{\code{BlackmanNuttallConstants=c(
      0.3635819, 0.4891775, 0.1365995, 0.0106411)}.}
    \item{\emph{Blackman--Harris}}{\code{BlackmanHarrisConstants=c(
      0.35875, 0.48829, 0.14128, 0.01168)}.}
    \item{\emph{Hann}}{\code{HannConstants=c()}.}
    \item{\emph{Hamming}}{\code{HammingConstants=c()}.}
    \item{\emph{Flattop Window}}{\code{FlattopConstants=c(
      0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368
    )}.}
  }


  \emph{Tukey Window}. This is a rectangular window with width
  \code{alpha*window} and then
  tapers off using a cosine window. At \eqn{\alpha=0} it is rectangular at
  \eqn{alpha=1} it is a cosine window.

  \deqn{w(\Delta t) = \begin{cases}
    \frac{1}{2} \left (1 + \cos \left ( \frac{2\pi | (\Delta t)/L - \alpha/2|}{\alpha}
    \right ) \right ) , \tab L \ge |\Delta t| > (1-\alpha/2) L \cr
    1, \tab (1-\alpha/2) L \ge |\Delta t| \ge 0 \cr
    \end{cases}.}

  \emph{Planck Window}.  Similar to the Tukey window,
  this has a flat region of size epsilon.

\deqn{w(\Delta t) = \begin{cases}
0, & |\Delta t| > L \cr
\left ( 1 + \exp \left ( \frac{2\epsilon L}{L-\Delta t} -
\frac{2\epsilon L}{(2\epsilon -1) L + \Delta t} \right ) \right)^{-1},
& L \ge |\Delta t| > (1-2 \epsilon)L \cr
1, & (1-2\epsilon) L \ge |\Delta t|\cr
\end{cases}
}


}
\value{

  A function with arguments \code{t,tstar} which gives a value between
  0 and 1 based the distance between the observation time (\code{t}) and
  the target time (\code{tstar}).  This is suitable for the
  \code{\link{IRTF}$wfun} field.
}
\references{

  A large number of window functions are listed on Wikipedia
  \url{https://en.wikipedia.org/wiki/Window_function}

}
\author{
Russell Almond
}
\seealso{
  The \code{\link{IRTF}} uses a window function as one of its fields, and the
  \code{\link{irtf}} function uses a weight function as one of its examples.
}
\examples{


curve(squareWindow()(x,0),xlim=c(-1.5,1.5))
curve(squareWindow(symmetric=FALSE)(x,0),xlim=c(-1.5,1.5))

curve(linearWindow()(x,0),xlim=c(-1.5,1.5))
curve(linearWindow(slope=.5,symmetric=FALSE)(x,0),xlim=c(-1.25,.25))

curve(exponentialWindow(lambda=.5)(x,0),xlim=c(-1.5,1.5))
curve(exponentialWindow(lambda=.5,symmetric=FALSE)(x,0),xlim=c(-3.25,.25))

curve(parzenWindow()(x,0),xlim=c(-1.5,1.5))
curve(parzenWindow(symmetric=FALSE)(x,0),xlim=c(-1.5,1.5))

curve(welchWindow()(x,0),xlim=c(-1.5,1.5))
curve(welchWindow(slope=.5,symmetric=FALSE)(x,0),xlim=c(-1.25,.25))

curve(sineWindow()(x,0),xlim=c(-1.5,1.5))
curve(sineWindow(power=2,symmetric=FALSE)(x,0),xlim=c(-1.25,.25))

curve(cosineSumWindow(a=HannConstants)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=HammingConstants)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=BlackmanConsts())(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=BlackmanExact)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=NuttallConstants)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=BlackmanNuttallConstants)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=BlackmanHarrisConstants)(x,0),xlim=c(-1.25,1.25))
curve(cosineSumWindow(a=FlattopConstants)(x,0),xlim=c(-1.25,1.25))

curve(tukeyWindow()(x,0),xlim=c(-1.5,1.5))
curve(tukeyWindow(alpha=.5,symmetric=FALSE)(x,0),xlim=c(-1.25,.25))

curve(plankWindow()(x,0),xlim=c(-1.5,1.5))
curve(plankWindow(alpha=.5,symmetric=FALSE)(x,0),xlim=c(-1.25,.25))

par(opar)
}
\keyword{ts}
\concept{filter}
\concept{window}

