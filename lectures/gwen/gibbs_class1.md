Sampling from the Normal model using Gibbs
------------------------------------------

Here is a set of data

    y = c(0.57, 0.71, -0.45, 0.92, -0.67, 3.04, 0.32, 1.38, 1.76, -0.14, -0.37, 0.69)

In class, we saw the model
$$ \\boldsymbol{y} \\overset{\\mathrm{iid}}{\\sim} N(\\theta, \\sigma^2) $$
 with conjugate priors
*θ* ∼ *N*(*μ*<sub>*n*</sub>, *τ*<sub>*n*</sub><sup>2</sup>)

and
$$ \\sigma^2 \\sim \\mathrm{Inv\\text{-}\\chi^2}(\\nu\_0+n, \\frac{\\nu\_0\\sigma^2\_0 + nv}{\\nu\_0+n}). $$
 The mean *μ*<sub>*n*</sub> and variance *τ*<sub>*n*</sub><sup>2</sup>
depend on *n* (the number of data points), the mean of the data
$\\bar{\\boldsymbol{y}}$, the known value of *σ*<sup>2</sup>, and the
hyperparameters *μ*<sub>0</sub> and *τ*<sub>0</sub>. The two parameters
of the scaled-inverse-*χ*<sup>2</sup> distribution (the degrees of
freedom and the scale factor) depend on *n*,
$v = \\frac{1}{n}\\Sigma\_i^{n}(y\_i - \\theta)^2$, and hyperparameters
*ν*<sub>0</sub> and *σ*<sub>0</sub><sup>2</sup>. See your notes from
class.

We're interested in the joint posterior distribution of the *θ* and
*σ*<sup>2</sup> parameters given the data **y**:
$$ p(\\theta,\\sigma^2|\\boldsymbol{y}) = \\frac{p(\\boldsymbol{y}|\\theta, \\sigma^2)p(\\theta,\\sigma^2)}{p(\\boldsymbol{y})}. $$
 We want to know what *p*(*θ*, *σ*<sup>2</sup>|**y**) looks like. We
don't have *p*(**y**), but it's just a constant anyway, so we can
rewrite this as
*p*(*θ*, *σ*<sup>2</sup>|**y**)∝*p*(**y**|*θ*, *σ*<sup>2</sup>)*p*(*θ*, *σ*<sup>2</sup>)

We assume independent prior probability distributions for *θ* and
*σ*<sup>2</sup>:
*p*(*θ*, *σ*<sup>2</sup>|**y**)∝*p*(**y**|*θ*, *σ*<sup>2</sup>)*p*(*θ*)*p*(*σ*<sup>2</sup>).

To set up the conjugate prior distributions, we must decide on fixed
values for the hyperparameters
{*μ*<sub>0</sub>, *τ*<sub>0</sub><sup>2</sup>, *σ*<sub>0</sub><sup>2</sup>, *ν*<sub>0</sub>}
(until Thursday, anyway...). Let's do that right now:

    mu0 = 2
    tausq0 = 4.3
    sigsq0 = 1.2
    nu0 = 1.2

We've set up our prior distributions by defining the values above, and
we are ready to draw samples from the posterior distribution!

To draw samples from the posterior, we could use a Metropolis sampler.
However, we haven't written out the distribution
*p*(**y**|*θ*, *σ*<sup>2</sup>) (i.e. the likelihood) and for whatever
reason don't want to (maybe it seems like a lot of work, or it is
totally intractable). Instead, we can be clever and use the
*conditional* distributions for *θ* and *σ*<sup>2</sup> to get samples
from the posterior distribution via a Gibbs Sampler!

A Gibbs sampler doesn't have to be tuned like a Metropolis algorithm
does. We can use a Gibbs sampler in this example because we have closed
form equations for the conditional posteriors, i.e.
$$ p(\\theta | \\boldsymbol{y}, \\sigma^2) \\\\ p(\\sigma^2|\\boldsymbol{y}, \\theta) $$
 (see your notes). It is easy to draw values from these distributions in
a Markov chain algorithm if you have some nice statistical software
(like R! or Python...?).

### Gibbs Sampler Algorithm

The following is mostly from Carlin & Lewis (2008), *Bayesian Methods
for Data Analysis*, CRC Press.

Just like in Metropolis, we have to choose the initial parameter values
to start the Markov chain. I'm going to call the set of parameters for
the model **ϑ** (so in our example, a normal model,
**ϑ** = {*θ*, *σ*<sup>2</sup>}).

The initial parameter values for the chain are denoted
**ϑ**<sup>(1)</sup>, and the next set in the chain will be
**ϑ**<sup>(2)</sup>, then **ϑ**<sup>(3)</sup>, etc. In our normal model
example, we need to set starting values for *θ* and *σ*<sup>2</sup>.
Let's set the initial values for the chain *θ*<sup>(1)</sup> and
*σ*<sup>2(1)</sup>

    theta1 = 0
    sigmasq1 = 1

With a Gibbs Sampler, you sample only one parameter at a time, given the
current values of all the other parameters. To get the next value of *θ*
in the chain, *θ*<sup>(2)</sup>, we don't have to *suggest* a new value,
we simply need to *draw* a new value from the conditional distribution
given the current value of *σ*<sup>2</sup> and the data:
*p*(*θ*|**y**, *σ*<sup>2(1)</sup>)

Drawing a random *θ* given the data and *σ*<sup>2(1)</sup>, gives us
*θ*<sup>(2)</sup>. Now we are ready to draw *σ*<sup>2(2)</sup>, given
the data and the updated value for *θ*. We draw *σ*<sup>2(2)</sup> from
*p*(*σ*<sup>2</sup>|**y**, *θ*<sup>(2)</sup>)
.

Now we go back and sample *θ* again! Get *θ*<sup>(3)</sup> by drawing
from
*p*(*θ*|**y**, *σ*<sup>2(2)</sup>)
 ... and then get *σ*<sup>2(3)</sup> by drawing from
*p*(*σ*<sup>2</sup>|**y**, *θ*<sup>(3)</sup>)
 Repeat!

Let's put this into action with our data **y**. Remember: we already set
values for the hyperparameters
{*μ*<sub>0</sub>, *τ*<sub>0</sub><sup>2</sup>, *σ*<sub>0</sub><sup>2</sup>, *ν*<sub>0</sub>}
above.

    # make an empty Markov chain of length T=1100 and with two columns
    chain = data.frame(theta = rep(NA_real_, 1100), sigsq = rep(NA_real_, 1100))
    # have a look at the first few rows of chain
    head(chain)

    ##   theta sigsq
    ## 1    NA    NA
    ## 2    NA    NA
    ## 3    NA    NA
    ## 4    NA    NA
    ## 5    NA    NA
    ## 6    NA    NA

    # define functions to calculate the mean and variance for the conditional distribution on theta (this will be needed for each sample draw)
    mun = function(mu0, y, tausq0, sigsq, n){
      (mu0/tausq0 + n*mean(y)/sigsq) / (1/tausq0 + n/sigsq)
    }

    taun = function(tausq0, sigsq, n){
      1/(1/tausq0 + n/sigsq)
    }

    # define a function that draws n samples from the scaled inverse chi-squared distribution with df degrees of freedom, and a scale parameter
    rinvchisq <- function(n,df,scale) (df * scale)/rchisq(n, df = df)

    # n is the number of data points
    n = length(y)

    # put the initial values in the chain, in the first row
    chain[1, ] = c(theta1, sigmasq1)
    head(chain)

    ##   theta sigsq
    ## 1     0     1
    ## 2    NA    NA
    ## 3    NA    NA
    ## 4    NA    NA
    ## 5    NA    NA
    ## 6    NA    NA

    # set the current sigma^2 as "newsigsq", which will get changed at every iteration
    newsigsq = sigmasq1

    for(t in 2:nrow(chain)){

      # sample a new theta from the conditional distribution for theta, given the data, hyperparameters, and current value for sigma^2 ("newsigsq")
      newtheta = rnorm(n = 1, mean = mun(mu0 = mu0, y = y, tausq0 = tausq0, sigsq = newsigsq, n=n), sd = sqrt(taun(tausq0 = tausq0, sigsq = newsigsq, n=n)))
            
      # sample a new sigma^2 from its conditional distribution, given the data, hyperparameters, and current value for theta ("newtheta")
      newsigsq = rinvchisq(n = 1, df = nu0+n, scale = (nu0*sigsq0 + sum((y-newtheta)^2))/(nu0 + n) )

      chain[t, ] = c(newtheta, newsigsq)
      
    }

    # plot the chain as a Markov Chain using the as.mcmc() function from the coda package. Don't plot the first 100 steps because we are treating that as a burn-in
    library(coda)
    plot(as.mcmc(chain[101:1100,]))

![](gibbs_class1_files/figure-markdown_strict/makeaGibbssampler-1.png)

    # use R's summary() function to get the summary statistics about the Markov chain (i.e. the samples from the posterior distribution)
    summary(as.mcmc(chain[101:1100, ]))

    ## 
    ## Iterations = 1:1000
    ## Thinning interval = 1 
    ## Number of chains = 1 
    ## Sample size per chain = 1000 
    ## 
    ## 1. Empirical mean and standard deviation for each variable,
    ##    plus standard error of the mean:
    ## 
    ##         Mean     SD Naive SE Time-series SE
    ## theta 0.6791 0.3314  0.01048        0.01048
    ## sigsq 1.3475 0.7004  0.02215        0.02445
    ## 
    ## 2. Quantiles for each variable:
    ## 
    ##          2.5%    25%    50%    75% 97.5%
    ## theta -0.0152 0.4715 0.6884 0.8789 1.357
    ## sigsq  0.5899 0.9177 1.1747 1.5620 2.976

To summarize and generalize, here's some pseudocode for a Gibbs sampler
algorithm when you have k parameters and want to make a chain of length
T,

For (*t* = 2, ..., *T*), repeat:

**Step 1:** Draw *ϑ*<sub>1</sub><sup>(*t*)</sup> from
*p*(*ϑ*<sub>1</sub>|**y**, *ϑ*<sub>2</sub><sup>(*t* − 1)</sup>, *ϑ*<sub>3</sub><sup>(*t* − 1)</sup>, …, *ϑ*<sub>*k*</sub><sup>(*t* − 1)</sup>)

**Step 2:** Draw *ϑ*<sub>2</sub><sup>(*t*)</sup> from
*p*(*ϑ*<sub>1</sub>|**y**, *ϑ*<sub>1</sub><sup>(*t* − 1)</sup>, *ϑ*<sub>3</sub><sup>(*t* − 1)</sup>, …, *ϑ*<sub>*k*</sub><sup>(*t* − 1)</sup>)

**Step k:** Draw *ϑ*<sub>*k*</sub><sup>(*t*)</sup> from
*p*(*ϑ*<sub>1</sub>|**y**, *ϑ*<sub>1</sub><sup>(*t* − 1)</sup>, *ϑ*<sub>2</sub><sup>(*t* − 1)</sup>, …, *ϑ*<sub>*k* − 1</sub><sup>(*t* − 1)</sup>)
