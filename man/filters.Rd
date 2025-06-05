\name{filters}
\alias{linearWindow}
\alias{squareWindow}
\alias{gaussianWindow}
\alias{parzenWindow}
\alias{welchWindow}
\alias{sineWindow}
\alias{cosineSumWindow}
\alias{tukeyWindwo}
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
gaussianWindow(window = 1, sigma=1, power=1, ..., symmetric = TRUE)
linearWindow(window=1, slope=1, ..., symmetric = TRUE)
parzenWindow(window=1,...,symmetric=TRUE)
welchWindow(window=1,slope=1,...,symmetric=TRUE)
sineWindow(window=1,power=1,...,symmetric=TRUE) 
BlackmanConsts(alpha=1.6)
cosineSumWindow(window=1,
                a=list(HannConstants,HammingConstants,
                       BlackmanConsts(1.6),BlackmanExact,
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

  Let \eqn{w(t,t_*)} be a function of the observation time \eqn{t} and the
  target time \eqn{t_*}.  Let \eqn{\Delta t = t_{*} - t}, and let \eqn{L} be the
  window width.  Expect where noted, if \eqn{|\Delta t| > L}, then \eqn{w(t,t_*)=0}.
  If \code{symmetric=FALSE}, then if \eqn{\Delta t < 0} then \eqn{w(t,t_*)=0}.

  \emph{Square Window}.

  \emph{Linear Window}.

  \emph{Gaussian Window}.

  \eqn{w(\Delta t) = \exp \left ( -\frac{1}{2} \left ( \frac{\Delta t}{\sigma L}
  \right )^p \right), } where \eqn{\sigma \leq 0.5}.   When \eqn{p=2}, then
  the standard deviation is \eqn{\sigma L}.

  \emph{Parzan Window}.

  \deqn{w(\Delta t) = \begin{cases} 1 - 6 \left( \frac{\Delta t}{L} \right )^2
  \left ( 1 - \frac{|\Delta t|}{L} \right )\ , \tab 0 \le |\Delta t| \le L/2\cr
  2 \left (1 - \frac{|\Delta t|}{L} \right )^3 \tab L/2 \le |\Delta t| \le L\cr
  \end{cases}}

  \emph{Welch Window}.

  \deqn{w(\Delta t) = 1 - \left ( \frac{\Delta t}{L} \right )^2}

  \emph{Sine Window}  (Hann window when squared)

  \deqn{w(\Delta t) = \cos^k \left ( \frac{\pi \Delta t}{L} \right )}


  \emph{Blackman Window} default \eqn{\alpha=0.16}.

  \dqn{w(\Delta t) = a_0 - a_1 \sin \left(\frac{2 \pi \Delta t}{2L} \right ) +
  a_2 \sin \left(\frac{4 \pi \Delta t}{2L} \right ), } where \eqn{a_0=(1-\alpha)/2},
  \eqn{a_1 = 1/2}, and \eqn{a_2 = \alpha/2}.

  \emph{Nuttall, Blackman--Nuttall, Blackman--Harris}.

  \dqn{w(\Delta t) = a_0 - a_1 \sin \left(\frac{2 \pi \Delta t}{2L} \right ) +
  a_2 \sin \left(\frac{4 \pi \Delta t}{2L} \right ) -
  a_3 \sin \left(\frac{6 \pi \Delta t}{2L} \right ),  }

  Nuttall \tab 0.355768 \tab 0.487396 \tab 0.144232 \tab 0.012604 \cr
  Blackman--Nuttall \tab 0.3635819 \tab 0.4891775 \tab 0.1365995 \tab 0.0106411 \cr
  Blackman--Harris \tab 0.35875 \tab 0.48829 \tab 0.14128 \tab 0.01168 \cr

  \emph{Flat top Window}.

   \dqn{w(\Delta t) = a_0 - a_1 \sin \left(\frac{2 \pi \Delta t}{2L} \right ) +
  a_2 \sin \left(\frac{4 \pi \Delta t}{2L} \right ) -
  a_3 \sin \left(\frac{6 \pi \Delta t}{2L} \right ) +
  a_4 \sin \left(\frac{8 \pi \Delta t}{2L} \right )
  }

  Matlab parameters.
  a0=0.21557895, a1=0.41663158, a2=0.277263158, a3=0.083578947, a4=0.006947368

  \emph{Tukey Window}.  At \eqn{\alpha=0} it is rectangular at \eqn{alpha=1}

  \deqn{w(\Delta t) = \begin{cases}
    \frac{1}{2} \left (1 - \sin \left ( \frac{2\pi | \Delta t - \alpha L|}{2 \alpha L}
    \right ) \right ) , \tab 0 \le |\Delta t| < \alpha L \cr
    1, \tab \alpha L \le |\Delta t| \le L \cr
    \end{cases}.}


  \emph{Planck Window}.

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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (window = 1, ..., symmetric = TRUE)
{
    if (symmetric)
        function(t, tstar) {
            as.numeric(abs(tstar - t) <= window)
        }
    else function(t, tstar) {
        as.numeric((tstar - t) <= window & tstar >= t)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
