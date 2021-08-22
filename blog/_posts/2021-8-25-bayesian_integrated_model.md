---
layout: post
title: So, you don't have enough data to fit a dynamic occupancy model? An introduction to auto-logistic occupancy models.
category: blog
---

Species distribution models (SDM) are useful tools in ecology and conservation biology. As the name implies, SDMs are used to estimate species presence or abundance across geographic space and/or through time. At their core, SDMs are a way to correlate environmental covariates--be they spatial, temporal, or both--to known locations of a speices. Generally this is done through some form of regression analysis.

SDMs are being regularly improved, and one class of model I am very excited about is integrated SDMs. These models were developed, in part, to take advantage of the ever growing data sources that ecologists have access to (e.g., data from public participitory science, museum collections, and wildlife surveys, among others). What is great about integrated SDMs is that they allow you to combine these different types of species location data within the same statistical model. And because of this, integrated SDMs usually (not always) improve your predictive accuracy and the precision of your model estimates. And while models that combine data sources are being used now more than ever, I would argue that they are still not seeing widespread use. So why is that?

In my opinion, the biggest hurdle that prevents the widescale adoption of integrated SDMs is that the math associated to them can be complex and use non-standard likelihood functions that may need to be coded up by hand. And to be honest, trying to understand that math is hard, especially as papers may need to skip some algebraic steps in the model formulation to save space! And while I am no expert in integrated SDMs, I do have some experience with them, and so wanted to dispell some of the confusion for at least one model I have been using in my own research, which is an integrated SDM that can combine opportunistic presence-only data (e.g., data from public participatory science) with detection / non-detection data (e.g., data from camera traps, repeated bird counts, etc.). 

So, in this post I am going to 1) walk through the model in Koshkina *et al.* (2017), 2) show how they can be written out in `JAGS`, and 3) simulate some data and show you that this model works. If you have not caught on by now from my other posts, simulating data is a sure fire way to ensure you've coded your model up correctly and therefore I do it whenever I'm trying to learn new statistical techniques. To keep some of the coding parts short, I've compartmentalized a lot of the code into functions, all of which is stored on this GitHub repository.


## The model

---

The Koshkina *et al.* (2017) integrated SDM, which was based on the model in Dorazio (2014), uses an inhomogeneous Poisson process model for the latent state. A Poisson process is just a mathematical object (e.g., a square with an x and y axis) in which points are randomly located. In our case, however, this is an inhomogeneous Poisson process, which means the density of points on this object depends on your location on the object. What this means is that we have some landscape and abundance of a species varies across that landscape.

So, for our study region, *B*, we have a set of locations that individuals of our species occupies <span>$$s_1,...,s_n$$</span> and those locations are assumed to be a inhomogenous Poisson process. Thus, the limiting expected density of our species (number of individuals per unit area) can be described by a non-negative intensity function <span>$$\lambda(s)$$</span>. Taking this one step further, the abundance of our species in region *B* is a Poisson random variable with mean:


$$\mu(B) = \int_B \lambda(s)ds$$


which just means we take the integral across our landscape (i.e., slice up the landscape into tiny little bits, count all the individuals in each slice, and add that all together). It's not difficult to make <span>$$\lambda(s)$$</span> a function of environmental variables, you use the log link and add in an intercept and slope term for each covariate <span>$$\boldsymbol{x}(s)$$</span>. The fact that <span>$$\boldsymbol{x}(s)$$</span> is bold means that there are a vector of covariates for each location, and the first element of <span>$$\boldsymbol{x}(s)$$</span> should be a 1 so that the intercept is always included. For example, if you have 3 covariates <span>$$\boldsymbol{x}(s)$$</span> would actually have a length of 4.  Thus, for *i* in 1,...,*I* covariates

$$\text{log}(\lambda(s)) = \boldsymbol{\beta}\boldsymbol{x}(s)^\intercal = \beta_1 \times x(s)_1 +\cdots+ \beta_I x(s)_I$$


and so our model of abundance across the whole study area is

$$N(B) \sim \text{Poisson}(\mu(B))$$

Some other important things to note here is that we can also model the number of individuals within any subregion of *B*, which also has a Poisson distribution. So, if you have a set of non-overlapping subregions where <span>$$r_1 \cup \cdots \cup r_k = B$$</span> (this means the union of all the subregions equals *B*) then <span>$$N(r_k) \sim \text{Poisson}(\mu(r_k))</span>. 

The reason I bring up the subregion modeling is for two reasons. First, when we write the model in `JAGS` we are going to have to discretize the landscape to approximate the integral in the latent state model. So, we may as well start thinking about what that looks like now. Second, creating subregions makes it more intuitive to think about modeling occupancy instead of abundance. Why? Because if we set out to study a species in an area, chances are we know a priori that the probability at least one individual occupies *B* is 1, but there may be variation in the different subregions. To model occupancy we need to convert<span>$$\mu(r_k)$$</span> into a probability, which can be done with the inverse complimentary log-log link (inverse cloglog)

$$\psi_k = Pr(N(r_k) > 0)  = (1 - \text{exp}({-\mu(r_k)}))$$

where <span>$$\psi_k$$</span> is the probability of occupancy for subregion <span>$$r_k$$</span>. This could then me modeled as a Bernoulli trial where <span>$$z_k = 1$$</span> if the species occupies subregion <span>$$r_k$$</span>.

$$z_k \sim \text{Bernoulli}(\psi_k)$$

And that is the latent state model! So far this should mostly remind you of Poisson regression (except for that little jaunt into occupancy modeling). The big exception, however, is that integral across space that we need to contend with in `JAGS`, and we'll get there soon enough.

Because we are fitting a Bayesian model, we get to take a small departure from the way the data model is described in Koshkina *et al.* (2017). In their model explanation, they initially combined the presence-only data detectability and probability of detection from the detection/non-detection data into one term. We, however, will keep them seperate the whole time. The detection / non-detection data model is the easier of the two, because it's equivalent to the detection process of a standard occupancy model, so let's start there.

Remember those <span>$$r_k$$</span> subregions of *B* that are non-overlapping parts of the study area? Assuming you do not sample every point in *B* and instead you sample some subset of the subregions. For *j* in 1,...,*J* subregions sampled (hereafter sites) let  <span>$$r_k[j]$$</span> represent the *jth* site sampled (<span>$$k[j]$$</span> means there is [nested indexing](https://masonfidino.com/nested_indexing/)). We then have for *w* in 1,..,*W* sampling occasions as these sites, which results in a *J* by *W* binary observation matrix <span>$$y_{j,w}$$</span>. If we detected the species at site *j* on sampling occasion *w*, <span>$$y_{j,w} = 1$$</span>, otherwise it is zero (assuming equal sampling at all sites). The probability of detecting the species given their presence <span>$$\rho_{j,w}$$</span> can then be estimated as:

$$y_{j,w}|z_{k[j]} \sim \text{Bernoulli}(\rho_{j,w} \times z_{k[j]})$$

and <span>$$\rho_{j,w}$$</span> can be made a function of covariates using the logit link, which could vary across sites or sampling occasions. For *q* in 1,...,*Q* covariates this is then

$$\text{logit}(\rho_{j,w}) = \boldsymbol{a}\boldsymbol{v}_{j,w}^\intercal = a_1 \times v_{1,j,w} + \cdots + a_Q \times v_{Q,j,w} $$

where *a* is a vector of regression coefficients (intercept and slope terms) and *V* is a *Q* by *J* by *W*  array where the first element of <span>$$\boldsymbol{v}_{j,w} = 1$$</span> to account for the model intercept.


Moving onto modeling the opportunistic presence-only data. We model this as a thinned Poisson process. What that means is that we multiply <span>$$\lambda(s)</span> by a probability (i.e., the presence-only detectabilty), which "thins" it. This thinning is helpful because we expect that the opportunistic presence-only data has some bias in where it was collected becuase there may not be a rigorous sampling design. As such, some areas may be oversampled while other areas are undersampled.  

The probability of being detected in the presence-only data, *b(s)*, can be made a function of *g* in 1,..,*G* covariates such that

$$\text{logit}(b(s)) = \boldsymbol{c}\boldsymbol{h}(s)^\intercal = c_1 \times h(s)_1 + \cdots + c_G \times h(s)_G$$


Where <span>$$\boldsymbol(c)$$</span> is a vector that contains the interecpt and slope terms  while <span>$$\boldsymbol(h)(s)$$</span> is a vector of G covariates, where the first element is 1 to account for the intercept  Before we move on there is one really important thing to bring up. For this model to be identifiable,  <span>$$\boldsymbol{x}(s)</span> must have one unique covariate that <span>$$\boldsymbol{h}(s)</span> does not have, and <span>$$\boldsymbol{h}(s)</span> must have one unique covariate <span>$$\boldsymbol{x}(s)</span> does not have. If you put the exact same covariates on the latent state and the presence-only data model, your model will not converge.

As this is a thinned Poisson process, the the expected number of presence only locations throughout *B* is just the product of <span>$$\lambda(s)$$</span> and<span>$$b(s)$$</span>, assuming that the species was detected at *m* locations <span>$$s_1,\cdots,s_m$$</span> where m < n (i.e., the number of individuals detected in the presence-only is less than the total number of individuals of that species on the landscape).

$$\pi(B) = \int_B \lambda(s)b(s)ds$$

Following Dorazio (2014) the likelihood of the presence only data is

$$\begin{eqnarray} 
L(\boldsymbol{\beta},\boldsymbol{c}) &=& \text{exp} \big(- \int_B\lambda(s)b(s)ds \big) \prod_{i=1}^m\lambda(s_i)b(s_i)      \nonumber \\
&=& \text{exp}\Big(-\int_B\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal + \boldsymbol{c}\boldsymbol{h}(s)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s)^\intercal )}ds\Big) \prod_{i=1}^m\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal + \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}  \nonumber
\end{eqnarray}$$

And if I had to guess, this is where people can get lost with this model. It's not necisarilly intuitive so it makes it difficult to know how to code it in your Bayesian programming language of choice!

But let's go through this piece by piece. To let you know we are headed, this equation is dividing the thinned Poisson process at all of the presence only points by the thinned Poisson process over the entire region *B* (i.e., the integral). Yet, there is no fraction here, so how did I arrive at this conclusion? Whenever you see negative signs thrown around with exponentials or natural logarithms, this can indicate that some division is going on. For example, if you remember that `exp(log(x)) = x`, let me show you two different ways to write the equation <span>$$10 / 2 = 5$$</span>.

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

If you look closely at that third example, its not too hard to recognize it as a simplified version of <span>$$\text{exp} \big(- \int_B\lambda(s)b(s)ds \big) \prod_{i=1}^m\lambda(s_i)b(s_i)$$</span>. So, because `exp(-x) * y == y/x`, then that first term that represents the thinned Poisson process over the region *B* is denominator of a fraction, while the second term that represents the thinned Poisson process over the *m* presence-only locations is the numerator. Therefore, the rest of the equation must take the log and logit scale values from the model, transform them a bit, and calcuate the thinned Poisson process (and that is exactly what it does).

For those familiar with logistic regression, you can hopefully notice that the likelihood

$$\text{exp}\Big(-\int_B\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal + \boldsymbol{c}\boldsymbol{h}(s)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s)^\intercal )}ds\Big) \prod_{i=1}^m\frac{\text{exp}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal + \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}{1 + \text{exp}(\boldsymbol{c}\boldsymbol{h}(s_i)^\intercal )}$$

is doing something a little different with the inverse logit link, <span>$$\text{logit}^{-1}(x) = \frac{\text{exp}(x)}{1 + \text{exp}(x)}$$</span>, which used to convert logit-scale coefficients back to the probability scale. In fact, it looks like they added the log-scale coefficients from the latent state model (<span>$$\boldsymbol{\beta}$$</span>) into the numerator of the inverse logit link, while the logit-scale coefficients (<span>$$\boldsymbol{c}$$</span>) from the thinning process are in the numerator and the denominator. To see what this means, let's explore how this modified inverse logit link function works in `R`, assuming that <span>$$\lambda(s) = 20$$</span> and <span>$$\b(s) = 0.75$$</span>.

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

Therefore, this modified inverse logit is just the inverse link function of the thinned Poisson process, <span>$$\lambda(s)b(s)</span>, which makes sense because this likelihood function is for the presence-only data! So if we wanted to abstract this likelihood out a little bit, we could write it out as:

$$L(\boldsymbol{\beta},\boldsymbol{c}) = \frac{\prod_{i=1}^m PO^{-1}(\boldsymbol{\beta}\boldsymbol{x}(s_i)^\intercal, \boldsymbol{c}\boldsymbol{h}(s_i)^\intercal)}{\text{exp}(\int_B PO^{-1}(\boldsymbol{\beta}\boldsymbol{x}(s)^\intercal, \boldsymbol{c}\boldsymbol{h}(s)^\intercal))ds}$$

where

$$PO^{-1}(x,y) = \frac{\text{exp}(x + y)}{1 + \text{exp}(y)}$$


Knowing this makes it much easier to code into `JAGS` because now we know what we are after. Now, let's code this up in `JAGS`.


## How to code up the model in `JAGS`







We are almost to the likelihood of the presence-only data model, and what I think was the trickiest part to unpack from this whole model.




One issue, however, is that integral cannot easily be evaluated in closed form as so we most approximate it as a Riemann sum. To do so, we split region *B* into a grid of cells, like you made a raster of *B* where each cell in the raster represents a portion of the total study area. I have seen suggestions that you need to break up the landscape into at least 10,000 cells or so, but in reality the number of cells is going to depend on the landscape heterogeneity you are trying to capture and your data. By discretizing space here, we need to slightly modify our linear predictor by adding a log offset term that denotes that area of a cell. Thus, for *c* in 1,...,*C* cells of region *B*, our linear predictor is now:


 in our specific case, the likelihood of our model can most easily be evaluated by discretizing our reg