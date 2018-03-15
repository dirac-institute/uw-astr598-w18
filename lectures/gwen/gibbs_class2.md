The terms *joint distribution*, *conditional distribution*, and
*marginal distribution* get thrown around a lot. Let's make sure we know
what they mean before jumping into hierarchical Bayesian models.

### Joint Distribution

The *joint distribution* is the multidimensional probability
distribution of all possible parameter values. For example, if you have
two parameters *θ*<sub>1</sub> and *θ*<sub>2</sub> in your model, their
joint posterior distribution would be
*p*(*θ*<sub>1</sub>, *θ*<sub>2</sub>).
 You can plot this distribution in a two-dimensional plot with
*θ*<sub>1</sub> on one axis and *θ*<sub>2</sub> on the other. The
distribution is usually shown with contours or a scatter plot if you
have samples from the distribution.

For example, let's pretend we have a *joint distribution* that is a
multivariate normal:

    library(MASS)
    library(ggplot2)
    # draw 5,000 samples from a multivariance normal distribution
    samples = mvrnorm(n = 5000, mu = c(1,1), Sigma = matrix(data = c(0.3, 0.1, 0.3, 0.5), nrow = 2, ncol = 2, byrow = TRUE))
    # name the columns
    colnames(samples) = c("theta1", "theta2")

    # make samples a data frame
    samples = as.data.frame(samples)

    # see what "samples" looks like using the structure function
    str(samples)

    ## 'data.frame':    5000 obs. of  2 variables:
    ##  $ theta1: num  1.042 0.341 0.976 0.782 1.849 ...
    ##  $ theta2: num  1.3517 -0.0448 0.5137 0.1029 0.9574 ...

    # plot the samples with some transparency
    # pch = point type to use
    plot(x=samples$theta1, y = samples$theta2, pch=19, col=rgb(0,0,1, alpha = 0.1), ylab=expression(theta[2]), xlab=expression(theta[1]), main="Samples from the joint distribution of a multivariate normal", xlim=c(-1,3), ylim=c(-2,4))
    grid()

<img src="gibbs_class2_files/figure-markdown_strict/examplejoint-1.png" width="400px" />

### Conditional Distribution

The *conditional distribution* is the distribution of a parameter given
a specific value of the other parameter(s). For example, in the above
plot, we could look at the conditional distribution of *θ*<sub>2</sub>
given that *θ*<sub>1</sub> = 1.1. We denote this as
*p*(*θ*<sub>2</sub>|*θ*<sub>1</sub> = 1.1).
 If we were to write out the actual mathematical equation for the above
distribution, it would be a function of *θ*<sub>2</sub> alone because we
set *θ*<sub>1</sub> = 1.1. A way to visualize the conditional
distribution of *θ*<sub>2</sub> given *θ*<sub>1</sub> = 1.1 is to take a
slice through the joint distribution *at the value
*θ*<sub>1</sub> = 1.1* and plot the distribution/histogram of
*θ*<sub>2</sub> at that value. Let's do that next!

<img src="gibbs_class2_files/figure-markdown_strict/conditional-1.png" width="400px" />

Because we have a finite sample from the joint distribution, very few
(if any) of the values of *θ*<sub>1</sub> equal 1.1. So let's make a
small window around *θ*<sub>1</sub> = 1.1 to approximate the conditional
distribution

    # get the indices of the rows that are TRUE in the logical statement 1.08 < theta1 < 1.12
    theseones = which(samples$theta1 > 1.08 & samples$theta1 < 1.12)

    # how many are in this range? 
    length(theseones)

    ## [1] 146

Plot the samples in a histogram with the same axis as the plot above, or
as a smoothed density plot:

<img src="gibbs_class2_files/figure-markdown_strict/conditionalthree-1.png" width="400px" /><img src="gibbs_class2_files/figure-markdown_strict/conditionalthree-2.png" width="400px" />

We could have chosen any value of *θ*<sub>1</sub> and found the
conditional distribution for *θ*<sub>2</sub> given that value. For
example, choosing *θ*<sub>1</sub> = 0.2 gives a different conditional
distribution than *θ*<sub>1</sub> = 1.1:

<img src="gibbs_class2_files/figure-markdown_strict/anotherexample-1.png" width="400px" /><img src="gibbs_class2_files/figure-markdown_strict/anotherexample-2.png" width="400px" />

Because the conditional distribution depends on the value of
*θ*<sub>1</sub>, we can generalize and write out the conditional
distribution as a *function* of *θ*<sub>2</sub> and *θ*<sub>1</sub>,

*p*(*θ*<sub>2</sub>|*θ*<sub>1</sub>).
 Once we have this *function*, then we can calculate the conditional
distribution of *θ*<sub>2</sub> for any value of *θ*<sub>1</sub>.

Of course, this can all be done the other way around--- we could find
the conditional distribution of *θ*<sub>1</sub> given *θ*<sub>2</sub>,
*p*(*θ*<sub>1</sub>|*θ*<sub>2</sub>).

When you know the conditional distributions for a model, then you can
use a Gibbs Sampler. In the case of known conditional distributions, a
Gibbs Sampler is often more efficient than a regular Metropolis
algorithm. (Recall the exercise from last class!)

### Marginal Distribution

The *marginal distribution* is the distribution of a parameter
regardless of the values of the other parameter(s). Mathematically, you
can find the marginal distribution by integrating out the unwanted,
uninteresting, or "nuisance" parameters from the joint distribution:

*p*(*θ*<sub>1</sub>)=∫*p*(*θ*<sub>1</sub>, *θ*<sub>2</sub>)*d**θ*<sub>2</sub>.
 Using our example of a multivariate normal distribution, we can look at
the marginal distribution of *θ*<sub>1</sub> and *θ*<sub>2</sub> by
making histograms of those parameters. Visually, this is like collapsing
all of the points in our scatterplot onto one axis and plotting a
histogram. Below we do this for both *θ*<sub>1</sub> and
*θ*<sub>2</sub>:

    library(psych)

    ## 
    ## Attaching package: 'psych'

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha

    # this package has a nice plotting function for a scatter plot plus side histograms

    scatter.hist(x = samples$theta1, y = samples$theta2, ellipse=FALSE, correl=FALSE, density = TRUE, freq = FALSE, smooth = FALSE, pch=19, col=rgb(0,0,1,0.1), grid=TRUE, title="Marginal distributions (top and side)", xlab = expression(theta[1]), ylab=expression(theta[2]))

<img src="gibbs_class2_files/figure-markdown_strict/marginaldist-1.png" width="500px" />

Notice that the marginal distributions are centered on the expected
(i.e. the mean) values for *θ*<sub>1</sub> and *θ*<sub>2</sub>.

When you have a Markov Chain and you plot the histogram of the samples
of one of the parameters, you are looking at (an estimate of) the
marginal distribution.
