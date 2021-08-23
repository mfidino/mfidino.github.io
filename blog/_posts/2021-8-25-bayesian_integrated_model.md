---
layout: post
title: A gentle introduction to an integrated occupancy model that combines presence-only and detection/non-detection data, and how to fit it in `JAGS`
category: blog
---

Species distribution models (SDM) are useful tools in ecology and conservation biology. As the name implies, SDMs are used to estimate species presence or abundance across geographic space and/or through time. At their core, SDMs are a way to correlate environmental covariates--be they spatial, temporal, or both--to known locations of a speices. Generally this is done through some form of regression analysis.

SDMs are being regularly improved, and one class of model I am very excited about is integrated SDMs. These models were developed, in part, to take advantage of the ever growing data sources that ecologists have access to (e.g., data from public participitory science, museum collections, and wildlife surveys, among others). What is great about integrated SDMs is that they allow you to combine these different types of species location data within the same statistical model. And because of this, integrated SDMs usually (not always) improve the predictive accuracy and precision of your model estimates. And while models that combine data sources are being used now more than ever, I would argue that they are still not seeing widespread use. So why is that?

In my opinion, the biggest hurdle that prevents the widescale adoption of integrated SDMs is that the math associated to them can be complex and use non-standard likelihood functions that may need to be coded up by hand. And to be honest, trying to understand that math is hard, especially as papers may need to skip some algebraic steps in the model formulation to save space! And while I am no expert in integrated SDMs, I want to take the time to dispel some of the confusion for one model I have been using in my own research, which is an integrated SDM that combine opportunistic presence-only data (e.g., data from public participatory science) with detection / non-detection data (e.g., data from camera traps, repeated bird counts, etc.). 

In this post I am going to 1) walk through the model in [Koshkina *et al.* (2017)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12738), 2) show how the model can be written out in `JAGS`, and 3) simulate some data and fit the model. To keep some of the coding parts short, I've compartmentalized a lot of the code into functions, all of which is stored on this GitHub repository.

Here are some links to the other parts of the post in case you want to skip around.

1. [The model](#markdown-header-the-model)
2. [How to code up the model in `JAGS`](#markdown-header-how-to-code-up-the-model-in-jags)
3. 



## The model

---

The Koshkina *et al.* (2017) integrated SDM, which was based on the model in Dorazio (2014), uses an inhomogeneous Poisson process model for the latent state. A Poisson process is just a mathematical object (e.g., a square with an x and y axis) in which points are randomly located. An inhomogeneous Poisson process, however, means the density of points on this object depends on your location on the object. Ecologically, this means that there is some region in geographic space (e.g., the city of Chicago, Illinois) and the modeled species abundance in that region varies across space.

So, for our region, *B*, we have a set of locations that individuals of our species occupies <span>$$s_1,...,s_n$$</span> and those locations are assumed to be a inhomogenous Poisson process. Thus, the limiting expected density of our species (number of individuals per unit area) can be described by the non-negative intensity function <span>$$\lambda(s)$$</span>. Taking this one step further, the abundance of our species in region *B* is a Poisson random variable with mean:


$$\mu(B) = \int_B \lambda(s)ds$$


which just means we take the integral across our landscape (i.e., slice up the landscape into tiny little bits, count all the individuals in each slice, and add that all together). It's not difficult to make <span>$$\lambda(s)$$</span> a function of environmental variables, you use the log link and add in an intercept and slope term for each covariate <span>$$\boldsymbol{x}(s)$$</span>. The fact that <span>$$\boldsymbol{x}(s)$$</span> is bold means that there are a vector of covariates for each location, and the first element of <span>$$\boldsymbol{x}(s)$$</span> should be a 1 so that the intercept is always included. For example, if you have 3 covariates <span>$$\boldsymbol{x}(s)$$</span> would actually have a length of 4.  Thus, for *i* in 1,...,*I* covariates

$$\text{log}(\lambda(s)) = \boldsymbol{\beta}\boldsymbol{x}(s)^\intercal = \beta_1 \times x(s)_1 +\cdots+ \beta_I x(s)_I$$


therefore the abundance model across the whole study area is

$$N(B) \sim \text{Poisson}(\mu(B))$$

Some other important things to note here is that we can also model the number of individuals within any subregion of *B*, which also has a Poisson distribution. So, if you have a set of non-overlapping subregions where <span>$$r_1 \cup \cdots \cup r_k = B$$</span> (this means the union of all the subregions equals *B*) then <span>$$N(r_k) \sim \text{Poisson}(\mu(r_k))$$</span>. 

I bring up this subregion modeling for two reasons. First, when we write the model in `JAGS` we are going to have to discretize the landscape to approximate the integral in the latent state model. So, we may as well start thinking about what breaking up the region into little pieces looks like now. Second, creating subregions makes it more intuitive to think about modeling occupancy instead of abundance. Why? Because if we set out to study a species in an area, chances are we know a priori that the probability at least one individual occupies *B* is 1, but there may be variation in the different subregions. To model occupancy we need to convert<span>$$\mu(r_k)$$</span> into a probability, which can be done with the inverse complimentary log-log link (inverse cloglog)

$$\psi_k = Pr(N(r_k) > 0)  = (1 - \text{exp}({-\mu(r_k)}))$$

where <span>$$\psi_k$$</span> is the probability of occupancy for subregion <span>$$r_k$$</span>. This could then me modeled as a Bernoulli trial where <span>$$z_k = 1$$</span> if the species occupies subregion <span>$$r_k$$</span>.

$$z_k \sim \text{Bernoulli}(\psi_k)$$

And that is the latent state model! So far this should mostly remind you of Poisson regression (except for that little jaunt into occupancy modeling). The big exception, however, is that integral across space that we need to contend with in `JAGS`, and we'll get there soon enough.

Because we are fitting a Bayesian model, we get to take a small departure from the way the data model is described in Koshkina *et al.* (2017). In their model explanation, they initially combined the presence-only data detectability and probability of detection from the detection/non-detection data into one term. We, however, will keep them seperate the whole time. The detection / non-detection data model is the easier of the two, because it's equivalent to the detection process of a standard occupancy model, so let's start there.

Remember those <span>$$r_k$$</span> subregions of *B* that are non-overlapping parts of the study area? Let's assume we do not sample every point in *B* and instead sample some subset of the subregions. For *j* in 1,...,*J* subregions sampled (hereafter sites) let  <span>$$r_k[j]$$</span> represent the *jth* site sampled (<span>$$k[j]$$</span> means there is [nested indexing](https://masonfidino.com/nested_indexing/)). We then have for *w* in 1,..,*W* sampling occasions at these sites, which results in a *J* by *W* binary observation matrix <span>$$y_{j,w}$$</span>. If we detected the species at site *j* on sampling occasion *w*, <span>$$y_{j,w} = 1$$</span>, otherwise it is zero (assuming equal sampling at all sites). The probability of detecting the species given their presence <span>$$\rho_{j,w}$$</span> can then be estimated as:

$$y_{j,w}|z_{k[j]} \sim \text{Bernoulli}(\rho_{j,w} \times z_{k[j]})$$

and <span>$$\rho_{j,w}$$</span> can be made a function of covariates using the logit link, which could vary across sites or sampling occasions. For *q* in 1,...,*Q* covariates this is the linear predictor of the detection / non-detection data model is

$$\text{logit}(\rho_{j,w}) = \boldsymbol{a}\boldsymbol{v}_{j,w}^\intercal = a_1 \times v_{1,j,w} + \cdots + a_Q \times v_{Q,j,w} $$

where *a* is a vector of regression coefficients (intercept and slope terms) and *V* is a *Q* by *J* by *W*  array where the first element of <span>$$\boldsymbol{v}_{j,w} = 1$$</span> to account for the model intercept.


Modeling the opportunistic presence-only data is a little more complex because we assume the presence-only data are a thinned Poisson process. What that means is that we multiply <span>$$\lambda(s)</span> by some estimated probability (i.e., the presence-only detectabilty), which "thins" it. This thinning is helpful because we expect that the opportunistic presence-only data has some bias in where it was collected becuase there may not be a rigorous sampling design. As such, some areas may be oversampled while other areas are undersampled and we do not want that to bias our species abundance estimate.  

The probability of being detected in the presence-only data, *b(s)*, can be made a function of *g* in 1,..,*G* covariates such that

$$\text{logit}(b(s)) = \boldsymbol{c}\boldsymbol{h}(s)^\intercal = c_1 \times h(s)_1 + \cdots + c_G \times h(s)_G$$


Where <span>$$\boldsymbol(c)$$</span> is a vector that contains the interecpt and slope terms  while <span>$$\boldsymbol(h)(s)$$</span> is a vector of G covariates, where the first element is 1 to account for the intercept  Before we move on there is one really important thing to bring up. For this model to be identifiable,  <span>$$\boldsymbol{x}$$(s)</span> must have one unique covariate that <span>$$\boldsymbol{h}(s)$$</span> does not have, and <span>$$\boldsymbol{h}(s)$$</span> must have one unique covariate <span>$$\boldsymbol{x}(s)$$</span> does not have. If you put the exact same covariates on the latent state and the presence-only data model, your model will not converge. Keep that in mind while you consider hypotheses for your own data!

As this is a thinned Poisson process, the the expected number of presence only locations throughout *B* is the product of <span>$$\lambda(s)$$</span> and<span>$$b(s)$$</span>, assuming that the species was detected at *m* locations <span>$$s_1,\cdots,s_m$$</span> where m < n. In other words, we make the assumption that the number of individuals detected in the presence-only data is less than the total number of individuals of that species on the landscape. Thus,

$$\pi(B) = \int_B \lambda(s)b(s)ds$$

Following Dorazio (2014) the likelihood of the presence only data is

$$\begin{eqnarray} 
L(\boldsymbol{\beta},\boldsymbol{c}) &=& \text{exp} \big(- \int_B\lambda(s)b(s)ds \big) \prod_{i=1}^m\lambda(s_i)b(s_i)      \nonumber \\
&=& \text{exp}\Big(-\int_B\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal + \boldsymbol{c}\boldsymbol{h}(s)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s)^\intercal )}ds\Big) \prod_{i=1}^m\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal + \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}  \nonumber
\end{eqnarray}$$

And if I had to guess, this is where people can get lost with this model! This likelihood is not necisarilly that intuitive so it makes it difficult to know how to code it in your Bayesian programming language of choice. But if we go through it piece by piece, the math here is not all that scary.

To let you know we are headed, this equation is dividing the thinned Poisson process at all of the presence only data points by the thinned Poisson process over the entire region *B* (i.e., the integral). Yet, there is no fraction here, so how did I arrive at this conclusion? Whenever you see negative signs thrown around with exponentials or natural logarithms, there could be some "secret division". For example, if you remember that `exp(log(x)) = x`, let me show you two different ways to write the equation <span>$$10 / 2 = 5$$</span> using exponentials and natural logs.

$$\begin{eqnarray}
10 / 2 &=& 5 \nonumber \\
\text{exp}\big(\text{log}(10) - \text{log}(2)\big) &=& 5 \nonumber \\
\text{exp}(-\text{log}(2)) \times 10 &=& 5 \nonumber
\end{eqnarray}$$

This is also really easy to check in `R` as well

```R
10 / 2
[1] 5

exp(log(10) - log(2))
[1] 5

exp(-log(2)) * 10
[1] 5

```

If you look closely at that third example, its not too hard to recognize it as a simplified version of <span>$$\text{exp} \big(- \int_B\lambda(s)b(s)ds \big) \prod_{i=1}^m\lambda(s_i)b(s_i)$$</span>. So, because `exp(-x) * y == y/x`, the first term term in our likelihood that represents the thinned Poisson process over the region *B* is denominator of a fraction, while the second term in our likelihood that represents the thinned Poisson process over the *m* presence-only locations is the numerator. With that out of the way, we just need to figure out what the likelihood is doing with the regression coefficients from the latent state model (<span>$$\beta$$</span>) and the regression coefficients associated to the presence-only thinning probability (<span>$$\boldsymbol{c}$$</span>).

For those familiar with logistic regression, you can hopefully notice that the likelihood

$$\text{exp}\Big(-\int_B\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal + \boldsymbol{c}\boldsymbol{h}(s)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s)^\intercal )}ds\Big) \prod_{i=1}^m\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal + \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}$$

made a slight modification to the inverse logit link, <span>$$\text{logit}^{-1}(x) = \frac{\text{exp}(x)}{1 + \text{exp}(x)}$$</span>. The inverse logit link is important to know, especially if you do any amount of occupancy modeling, because it converts logit-scale coefficients back to the probability scale. It looks like this likelihood function added the log-scale coefficients from the latent state model (<span>$$\boldsymbol{\beta}$$</span>) into the numerator of the inverse logit link, while the logit-scale coefficients (<span>$$\boldsymbol{c}$$</span>) from the thinning process are in the numerator and the denominator. To see what this means, let's explore how this modified inverse logit link function works in `R`, assuming that <span>$$\lambda(s) = 20$$</span> and <span>$$b(s) = 0.75$$</span>.

```R
# The parameters
lambda <- 20
b <- 0.75

# Transform these variables to their
#  respective scales. 

# Log link for lambda_s
log_lambda <- log(lambda)

# logit link for b_s
logit_b <- log(b/(1 - b)) 

# Modified inverse logit
answer <- exp(log_lambda + logit_b) / (1 + exp(logit_b))

# compare that to lambda * b
#  both == 15
c(answer, lambda * b)
```

Therefore, this modified inverse logit is the inverse link function for the regression coefficients associated to the thinned Poisson process <span>$$\lambda(s)b(s)$$</span>! Knowing that, we could take the presesence-only data likelihood and abstract it out a little further:

$$L(\boldsymbol{\beta},\boldsymbol{c}) = \frac{\prod_{i=1}^m PO^{-1}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal, \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal)}{\text{exp}(\int_B PO^{-1}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal, \boldsymbol{c}\boldsymbol{h}(s)^\intercal ds))}$$

where

$$PO^{-1}(x,y) = \frac{\text{exp}(x + y)}{1 + \text{exp}(y)}$$


And that is the breakdown of the likelihood of the presence-only data. Hopefully this additional explanation here makes the model a bit easier to understand. In summary, this model has three linear predictors. One for the latent state that is partially informed by the presence-only data and detection/non-detection data, one to account for the bias in the opportunistic presence-only data, and one to account for false absencses in the detection / non-detection data.

<p><a href="#top" style>Back to top ⤒</a></p>

## How to code up the model in `JAGS`

For this example I'm going to focus on modeling species occupancy instead of abundance, though the model can easily be modified to estimate abundance by changing the likelihood function of the latent state model.

To code up this model in `JAGS` there are two things we need to overcome.

1. Approximating the integrals in the model with a Riemann sum.
2. Incorporate the non-standard likelihood of the presence-only data model with the Bernoulli "ones trick".

The first step is not difficult. For the region *B*, all we need to do is break it up into a "sufficiently fine grid" (Dorazio 2014) and evaluate the model in each grid cell.
![A landscape with a grid of cells]({{site.url}}/blog/images/iocm01.jpeg#center)

 What is sufficiently fine? I've seen some suggestions that you [need at least 10000](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12352), but in reality the number of cells is going to depend on the landscape heterogeneity you are trying to capture and your data. For example, a small study region roughly the size of a city may not need 10000 grid cells. In my own research, I've had success splitting Chicago, IL, which is about 606 km<sup>2</sup>, into roughly 2500 500m<sup>2</sup> grid cells. So if you are trying to use this model on your own you may need to do some trial and error to determine what is an appropriate cell count is (i.e., fit the model with different grid cell sizes to see if it has an influence on the parameters in the latent state model). 

Perhaps the trickiest part related to discretizing the region, however, means all of the data input into the model needs to be aggregated to the scale of the grid (i.e., the covariates and species data). I'll show you how to do that when I simulate some data.

Now back to the model. Let `npixel` be the total number of grid cells on our landscape and `cell_area` be the log area of the gridded cell (we have to include a log offset term into the latent state linear predictor to account for the the gridding we did). The latent state model is then:

```R
for(pixel in 1:npixel){
	# latent state linear predictor
	#
	# x_s  = covariates for latent state
	# beta = latent state regression coefficients
	# cell_area = log area of grid cell 
	#
	log(lambda[pixel]) <-inprod(x_s[pixel,1:I], beta[1:I]) + cell_area
	# Species presence in a gridcell as a Bernoulli trial
	z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
}

```

Remember, however, that the presence-only data model needs to be evaluated across the entire landscape as well (i.e., the integral that is the denominator of the presence-only data model likelihood). Therefore, we should also calcluate the thinning probability for each cell here too.

```R
for(pixel in 1:npixel){
	# latent state linear predictor
	#
	# x_s  = covariates for latent state
	# beta = latent state model regression coefficients
	# cell_area = log area of grid cell 
	#
	log(lambda[pixel]) <-inprod(x_s[pixel,], beta) + cell_area
	# Species presence in a gridcell as a Bernoulli trial
	z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
	# presence only thinning prob linear predictor
	#
	# h_s = covariates for thinning probability
	# cc  = presence-only data model regression coefficients
	#
	logit(b[pixel]) <-  inprod(h_s[pixel,] , cc)
}

```
And now we can start coding up the likelihood of the presence-only data model. While you can easily create your own likelihood functions in `Nimble`, it's not quite as simple in `JAGS`. However, there are some workarounds, and my personal favorite is the [Bernoulli "ones trick"](https://stackoverflow.com/questions/45888970/specify-a-discrete-weibull-distribution-in-jags-or-bugs-for-r/45920885#45920885). 

To do this for the presence-only data model there are a two steps. First, code up the likelihood for each of the `m` presence-only datapoints. Second, input that likelihood for each data point into a Bernoulli trial and divide it by some large constant value to ensure the likelihood is between 0 and 1 such that `ones[i] ~ dbern(likelihood[i]/CONSTANT)`, where `ones` is a vector of ones the same length as `m`. This little trick evaluates the likelihood of datapoint `i` given the current model parameters in the MCMC chain, which is exactly what we need. 

Note, however, that the presence-only data model likelihood function presented in Koshkina *et al.* (2017)is for ALL the datapoints, and we need to evaluate the likelihood of each data point on its own. To do this, all we need to do is the product from the numerator bceause we are only looking at a single datapoint from the numerator and divide the denominator by `m` (i.e., the number of presence-only datapoints). If you recall that `exp(log(x) - log(y)) = x/y` we can write the presence-only likelihood in `JAGS` using [nested indexing](https://masonfidino.com/nested_indexing/) to map the `m` datapoints back to their respective grid cell. If I wanted to I could have coded up the presence-only data model inverse link function and applied it to the regression coefficients & covariates for each cell. However, we have already calculated `lambda` and `b` in the model so we can use them instead.

```R
# The presence only data model.
#
# This part of the model uses
#  what we have calculated above (lambda
#  and b). The denominator of this likelihood
#  is a scalar so we can calculate it
#  outside of a for loop. Let's do that first,
#  which we are approximating as a Riemann sum.
#
po_denominator <- inprod(lambda[1:npixel], b[1:npixel] ) / m
#
# Loop through each presence-only datapoint
#  using Bernoulli one's trick. The numerator
#  is just the thinned poisson process for
#  the ith data point, and we use some
#  nested indexing to map the data to 
#  it's associated grid cell.
for(po in 1:m){
  ones[po] ~ dbern(
  	exp(
  	  log(lambda[po_pixel[po]]*b[po_pixel[po]]) -
  	  log(po_denominator)
	) / CONSTANT
  )
} 
```
Here, `po_pixel` denotes the grid cell of the *ith* presence only datapoint.

Finally, the detection/non-detection data model is basically what you'd see in a standard occupancy model, but again uses nested indexing to map the sampled site to it's respective grid cell. I'm making the assumption here that the `W` sampling occasions are i.i.d. random variables and so will model them as a binomial process such that

```R
for(site in 1:nsite){
  # detection/non-detection data model linear predictor
  #
  #  v = detection covariates for the entire region B
  #  a = det/non-det data model regression coefficients
  #
  logit(rho[site]) <-inprod(v[pa_pixel[site], ],a)
  # The number of detections for site is a binomial
  #  process with Pr(rho[site]) | z = 1 with
  #  W sampling occasions per site.
  y[site] ~ dbin(
    z[pa_pixel[site]] * rho[site],
    W
  )
}
```
Here is all of the model put together plus some non-informative priors for all the parameters and a derived parameter to track how many cells the species occupies. You can find this model on [this repository here].
```R
model{
  # Bayesian version of the Koshkina (2017) model.
  #
  # The latent-state model
  for(pixel in 1:npixel){
    # latent state linear predictor
    #
    # x_s  = covariates for latent state
    # beta = latent state model regression coefficients
    # cell_area = log area of grid cell 
    #
    log(lambda[pixel]) <-inprod(x_s[pixel,], beta) + cell_area
    # Species presence in a gridcell as a Bernoulli trial
    z[pixel] ~ dbern(1 - exp(-lambda[pixel]))
    # presence only thinning prob linear predictor
    #
    # h_s = covariates for thinning probability
    # cc  = presence-only data model regression coefficients
    #
    logit(b[pixel]) <-  inprod(h_s[pixel,] , cc)
  }
  # The presence only data model.
  #
  # This part of the model just uses the
  #  what we have calculated above (lambda
  #  and b). The denominator of this likelihood
  #  is actually a scalar so we can calculate it
  #  outside of a for loop. Let's do that first.
  #
  # The presence_only data model denominator, which
  #  is the thinned poisson process across the
  #  whole region (divided by the total number of 
  #  data points because it has to be 
  #  evaluated for each data point).
  po_denominator <- inprod(lambda[1:npixel], b[1:npixel] ) / m
  #
  # Loop through each presence-only datapoint
  #  using Bernoulli one's trick. The numerator
  #  is just the thinned poisson process for
  #  the ith data point.
  for(po in 1:m){
    ones[po] ~ dbern(
      exp(
        log(lambda[po_pixel[po]]*b[po_pixel[po]]) -
          log(po_denominator)
      ) / CONSTANT
    )
  } 
#
# Detection / non-detection data model
for(site in 1:nsite){
  # detection/non-detection data model linear predictor
  #
  #  v = detection covariates for the entire region B
  #  a = det/non-det data model regression coefficients
  #
  logit(rho[site]) <-inprod(v[pa_pixel[site], ],a)
  # The number of detections for site is a binomial
  #  process with Pr(rho[site]) | z = 1 with
  #  W sampling occasions per site.
  y[site] ~ dbin(
    z[pa_pixel[site]] * rho[site],
    W
  )
}
#
# Priors for latent state model
for(latent in 1:nlatent){
  beta[latent] ~ dnorm(0, 0.01)
}
# Priors for presence-only data model
for(po in 1:npar_po){
  cc[po] ~ dlogis(0, 1)
}
# Priors for det/non-det data model
for(pa in 1:npar_pa){
  a[pa] ~ dlogis(0, 1)
}
# Derived parameter, the number of cells occupied
zsum <- sum(z)
}
```

<p><a href="#top" style>Back to top ⤒</a></p>

## Simulating some data

To look at all the simulation code and try it for yourself, visit this [repository here](https://github.com/mfidino/integrated-occupancy-model). You just need to run through `simulate_poisson_process.R` to simulate the data and fit the models.

If you are not interested in that code, here are the steps I took to simulate the data:

1. Create the study region and give it a spatially autocorrelated environmental covariate.
2. Generate latent species abundance on the landscape, which is correlated to the environmental covariate from step 1.
3. Create a secondary environmental covariate that is correlated to the bias in the presence-only data.
4. Generate presence-only data, which samples from the latent species abundance in step 2 based on the environmental covariate in step 3 (i.e., some individuals are observed in the presence-only dataset).
5. Aggregate covariates and species presence on the landscape down to a coarser resolution for modeling purposes.
6. Place some sites to sample your detection / non-detection data across the landscape at the same scale as this new coarser resolution.
7. If the species is present at these sites, sample their presence with imperfect detection.


Following this we fit our simulated data to two models. The first model had the same latent-state and detection/non-detection data model, but we have removed the presence-only data model so that it is a standard occupancy model using a cloglog link. The second model is the integrated model from above fit to all of the data.

For this simulation, I assumed there was a single covariate that influenced the latent state and one covariate that influenced the bias of the presence-only data. For simplicity, I also assumed that imperfect detection at sampling sites did not vary (i.e., the data model for detection/non-detection data was intercept only). Some other general bookkeeping for this simulation:

1. One hundred sites were sampled for the detection / non-detection data.
2. The intercept and slope coefficients for the latent state were respectively `6` and `1`.
3. The intercept and slope coefficients for the presence-only thinning probability were respectively `-0.75` and `1.5`.
4. The by-survey probability fo detecting an individual given their presence at a sample site was 0.3 (about -0.85 on the logit scale). 




