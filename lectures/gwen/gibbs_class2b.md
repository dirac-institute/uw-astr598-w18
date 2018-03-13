Review of the Normal model
--------------------------

Last class, we introduced the normal model for data **y**, with mean *θ*
and variance *σ*<sup>2</sup> (following the notation and Chapter 2 of
Gelmen et al 2014, *Bayesian Data Analysis*, CRC Press). We started by
looking at two different scenarios--- one in which the mean was known
and the variance was unknown, and then where the mean was unknown and
the variance was known. In both cases, we derived the *conditional*
distribution of the unknown parameter, given the data and the known
parameter.

For example, the conditional distribution of the mean given the data and
the variance was
*p*(*θ*|**y**, *σ*<sup>2</sup>)=*p*(**y**, *σ*<sup>2</sup>|*θ*)*p*(*θ*)
 The prior for *θ* was a normal distribution with hyperparameters
*μ*<sub>0</sub> and *τ*<sub>0</sub><sup>2</sup>, and the conditional
distribution turned out to be a normal with mean *μ*<sub>*n*</sub> and
and variance *τ*<sub>*n*</sub><sup>2</sup>.

Similarly, we also found the conditional distribution for the variance,
given the data and a known mean:
*p*(*σ*<sup>2</sup>|**y**, *θ*)=*p*(**y**, *θ*|*σ*<sup>2</sup>)*p*(*σ*<sup>2</sup>),
 which turned out to be a scaled inverse-*χ*<sup>2</sup> distribution
with *ν*<sub>0</sub> + *n* degrees of freedom and scale factor
$\\frac{\\nu\_0\\sigma^2\_0 + nv}{\\nu\_0+n}$, where
$v = \\frac{1}{n}\\Sigma\_i^{n}(y\_i - \\theta)^2$, and *ν*<sub>0</sub>
and *σ*<sub>0</sub><sup>2</sup> are hyperparameters.

In Tuesday's class, we set the hyperparameters to specific values. In a
hierarchical Bayesian model, we don't do this--- instead we set
distributions for the hyperparameters, thus defining *hyperpriors*
*p*(*ν*<sub>0</sub>),*p*(*μ*<sub>0</sub>),*p*(*τ*<sub>0</sub><sup>2</sup>),and *p*(*σ*<sub>0</sub><sup>2</sup>).

Now we can use Bayes' rule recursively, so where we used to have the
joint posterior distribution of the *θ* and *σ*<sup>2</sup> parameters
given the data **y**:
*p*(*θ*, *σ*<sup>2</sup>|**y**)∝*p*(**y**|*θ*, *σ*<sup>2</sup>)*p*(*θ*)*p*(*σ*<sup>2</sup>).
 we now have
*p*(*θ*, *σ*<sup>2</sup>|**y**)∝*p*(**y**|*θ*, *σ*<sup>2</sup>)*p*(*θ*, *σ*<sub>0</sub><sup>2</sup>|*μ*<sub>0</sub>, *τ*<sub>0</sub><sup>2</sup>, *ν*<sub>0</sub>, *σ*<sub>0</sub><sup>2</sup>)*p*(*ν*<sub>0</sub>)*p*(*μ*<sub>0</sub>)*p*(*τ*<sub>0</sub><sup>2</sup>)*p*(*σ*<sub>0</sub><sup>2</sup>).

The *hyperprior distributions* will have their own set of parameters.
For example, if we set the hyperprior distribution *p*(*ν*<sub>0</sub>)
to be a uniform, truncated prior,
$$p(\\nu\_0) = \\text{unif}(a,b) \\\\
= \\left\\{ \\begin{array}{11}{ \\frac{1}{b-a}, \\nu\_0 \\in \[ a,b \] \\\\ 0, \\text{otherwise} }\\end{array} \\right.$$
 then we would have to choose values for the upper and lower bounds
(i.e. the parameters *a* and *b*) to the best of our ability.

Example of a Hierarchical Bayesian Model for Kuiper Belt Objects
----------------------------------------------------------------

*Disclaimer: I am not a solar system scientist, so please forgive my
possibly bad choices of distributions for the following example!*

Imagine we are interested in two quantities: 1) the masses of objects in
the Kuiper Belt (KB), and 2) the distribution of masses in the KB (in
particular the population mean and variance). Now suppose we have mass
estimates for *n* KB objects, and we label this data *y*<sub>*i*</sub>,
where *i* is the index of the object. With these estimates, we could
plot an empirical distribution of the masses, and estimate for the
population mean and variance. However, I think we can all agree that our
individual mass estimates are not the *true* masses of these objects,
and the uncertainty in the population mean is going to be rather large.

There is an arguably more natural way of estimating the true individual
masses and the mass distribution *at the same time*, and that is through
a hierarchical Bayesian model. In this set up, we assume that the masses
of all Kuiper Belt objects are following some underlying population
distribution, and fit for both the individual masses and the population
parameters simultaneously.

We keep our label for each KB object's mass estimate as data
*y*<sub>*i*</sub>, where *i* is the index of the object. But now we
assume that each mass measurement is drawn from a normal distribution
*y*<sub>*i*</sub> ∼ *N*(*μ*<sub>*i*</sub>, *σ*<sub>*i*</sub><sup>2</sup>)
 with mean *μ*<sub>*i*</sub> equal to the *true but unknown mass* for
object *i*, and variance *σ*<sub>*i*</sub><sup>2</sup> which is related
to our measurement uncertainty. For this example, we will assume that
our measurement uncertainty is consistent between objects, so that they
all share the same variance, i.e.
*y*<sub>*i*</sub> ∼ *N*(*μ*<sub>*i*</sub>, *σ*<sup>2</sup>).

Next, we assume that the natual log of the true masses, i.e.
log*μ*<sub>*i*</sub>, come from a normal distribution that defines the
*population* of KB objects. That is,
log*μ*<sub>*i*</sub> ∼ *N*(*μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>)
 where *μ*<sub>*p*</sub> is the true population mean and
*σ*<sub>*p*</sub><sup>2</sup> is the true population variance. We then
assign prior distributions to these hyperparameters
$$ \\mu\_p \\sim N(\\mu\_0, \\sigma^2\_0) \\\\
\\text{and} \\\\
\\sigma^2\_p \\sim \\text{Inv-Gamma}(\\alpha, \\beta).$$
 and set values for *μ*<sub>0</sub>, *σ*<sub>0</sub><sup>2</sup>, *α*,
and *β* to define these distributions.

This model has a three-structure hierarchy: the data level
(*y*<sub>*i*</sub>), the prior level (the population), and the
hyperprior level (distributions for the hyperparameters). For a visual
representation of the hierarchical model, see my sketch below.

<img src="KBhierarchy.png" width="400px" style="display: block; margin: auto;" />

Mathematically we can write this out using Bayes' Theorem. We are
interested in the true mass of each object, its variance, the population
mean, and the population variance, given the data set **y**), so the
posterior distribution is
*p*(**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>|**y**)∝*p*(**y**|**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>)*p*(**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>).
 Note that **μ** is bold because it is the vector of the individual true
masses of the objects. The variance *σ*<sup>2</sup> is assumed
independent of the population mean and variance, and is shared between
objects, so
*p*(**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>|**y**)∝*p*(**y**|**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>)*p*(*σ*<sup>2</sup>)*p*(**μ**, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>).
 Next, we use Bayes' theorem again to rewrite the joint distribution
*p*(**μ**, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>) in terms of
a conditional and prior distributions,
*p*(**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>|**y**)∝*p*(**y**|**μ**, *σ*<sup>2</sup>, *μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>)*p*(*σ*<sup>2</sup>)*p*(**μ**, *σ*<sup>2</sup>|*μ*<sub>*p*</sub>, *σ*<sub>*p*</sub><sup>2</sup>)*p*(*μ*<sub>*p*</sub>)*p*(*σ*<sub>*p*</sub><sup>2</sup>).
 Voila! We have a hierarchical Bayesian model.

### Sampling from the Posterior Distribution of a Hierarchical Model

The posterior distribution can be sampled in various ways, and this
subject could probably take up an entire week of classes. A lot of
statistical software already exists for doing hierarchical Bayesian
analysis, so before you run off and re-invent the wheel, take a look
below! :)

-   Stan <http://mc-stan.org/users/>
    -   see the "Stan Best Practices" github page for some solid advice
        that is not necessarily limited to Stan
        <https://github.com/stan-dev/stan/wiki/Stan-Best-Practices>
    -   Stan can be used via R (RStan), Python (PyStan), Julia
        (Stan.jl), Matlab (MatlabStan), the command line (CmdStan), etc.
-   BUGS - *Bayesian inference Using Gibbs Sampling*
    <https://www.mrc-bsu.cam.ac.uk/software/bugs/>

-   JAGS - *Just Another Gibbs Sampler*
    <http://mcmc-jags.sourceforge.net/>
