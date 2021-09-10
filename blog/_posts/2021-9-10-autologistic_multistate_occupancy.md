---
layout: post
title: Two ways to write autologistic multistate occupancy models in `JAGS`
category: blog
---

This is the second post in a series I am writing on multistate occupancy models. I will be building off my previous example of a static multistate model, so if you've not read the first post be sure to check it out [here](https://masonfidino.com/intro_2_multistate_occupancy/). Likewise, if want more of an introduction to autologistic occupancy models, please see this post [here](https://masonfidino.com/autologistic_occupancy_model/).

As a reminder, multistate occupancy models partition the probability of success (i.e., species presence) into multiple discrete states. In my previous example, we estimated the probability of three states:

1. That owls were not present
2. That owls were present but not breeding
3. That owls were present and breeding.

Autologistic occupancy models, on the other hand, are a simplified version of a dynamic occupancy model that estimates species occupancy instead of local colonization / extinction dynamics. They include an additional term in their logit-linear predictor to account for temporal dependence between primary sampling periods. Therefore, an autologistic multistate model can be used to estimate multiple discrete states across space and through time, and includes additional parameters to account for temporal dependence (e.g., a site where we found owls breeding at time *t* may be more likely to have owls that breed there at time *t+1*). This class of model has never received a formal write-up, and so it has rarely been used in the literature. Nevertheless, it is a very useful model to know when you are interested in applying a multistate model to data over time.

I'm going to walk through two different ways to write a three state model for a single species. The first model is the simplest parameterization, and uses the logit-link to estimate the occupancy probability of multiple states. I like this model for two reasons. First, it uses a standard link function that people should be familiar with, which makes it more approachable. Second, it is the absolute simplest parameterization of an autologistic multistate model, and so may be your best bet if you have a small sample size. In the second model, I use the softmax function to parameterize the model and include extra autologistic terms to estimate turnover among states over time. For example, if a site was in state 2 (owls present, not breeding) at time *t-1*, the probability it transitions to state 3 (owls present & breeding) at time *t* may be different than a site in state 3 that remains in state 3 from one time period to the next. This second model is useful for two reasons. First, the way this model is parameterized is very similar to multistate occupancy models for potentially interacting species. Therefore, understanding how this model works here will hopefully help you better understand how the [Rota *et al.* (2016)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12587) or [Fidino *et al.* (2018)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13117) models work. Second, it demonstrates that you have multiple choices to make when you attempt to use such a model, and so it's important to critically think about what it is you are interested in estimating (and whether or not you may have the sample size to do so).

Links to different parts of this post

1. [The transition probability matrix. What is how and how to use it.](#the-transition-probability-matrix-what-it-is-and-how-to-use-it)
2. [An autologistic multistate occupancy model that uses the logit link](#an-autologistic-multistate-occupancy-model-that-uses-the-logit-link)
  1. [Fitting the model in JAGS]
  2. [Plotting out the results]
3. [An autologistic multistate model occupancy that uses the softmax function]
  1. [Fitting the model in JAGS]
  2. [Plotting our the results]


Again, we are going to continue with the example of modeling the distribution of great horned owls and their breeding territory, except we now have data over time as well. Therefore, assume we have conducted *j* in 1,...,*J* surveys at *s* in 1,...,*S* sites throughout Chicago, Illinois during *t* in *1,...,T* great horned owl breeding seasons. During each of these surveys we can record 1 of 3 possible values:

1. No detection of an owl
2. Detection of an adult owl but no detection of a young owl
3. Detection of a young owl (i.e., reproductive success)

We store these observations in the *S* by *J* by *T* array <span>$$\boldsymbol{Y}$$</span>. Again, this lets us estimate three separate states: owls are not present (state = 1), owls are present but not breeding (state = 2), and owls are present and breeding (state = 3). However, the probability of occupancy of the different states at site *s* and time *t* is going to depend on the the state at time *t-1*. Therefore, instead of having a vector of occupancy probabilties, we actually place them in a transition probability matrix (TPM), which stores the probability of transitioning from one state to another in a single time unit:

$$
\boldsymbol{\psi_{s,t}} = 
\begin{bmatrix}
\psi_{s,t}^1 & \psi_{s,t}^{2} & \psi_{s,t}^3 \\
\psi_{s,t}^{1|2} & \psi_{s,t}^{2|2} & \psi_{s,t}^{3|2} \\
\psi_{s,t}^{1|3} & \psi_{s,t}^{2|3} & \psi_{s,t}^{3|3}
\end{bmatrix}
 $$

One important aspect to remember about TPMs is that one dimension of the matrix represents "from this state" and the other dimension represents "to that state". In our case, the rows are the former and the columns are the latter. Because of this, each row of our TPM sums to 1. For example, the probability of transitioning from state 3 to state 2 means we need to look for the element in row 3 ("from" state 3) and column 2 ("to" state 2), which is <span>$$\psi_{s,t}^{2|3}$$</span>. I've added the conditional probabilities (e.g,. 2|3) on the 2nd and 3rd rows to remind you that we are including autologistic terms to help account for temporal dependence between time units.

Before I break down how TPMs apply to autologistic occupancy models, however, I want to acknowledge that interpretting TPMs can be difficult! As such, I'm going to work through a simple example of a TPM, explain some important qualities that they may have, and introduce you to the `markovchain` package in `R`, which has some useful functions for dealing with TPMs.

## The transition probability matrix. What it is and how to use it.

Imagine we have this TPM:

$$TPM = \begin{bmatrix}
0.3 & 0.4 & 0.3 \\
0.6 & 0.2 & 0.2 \\
0.1 & 0.5 & 0.4
\end{bmatrix}$$

It has three states, and the probability of transitioning from one state to another in a single time unit varies based on what state you are currently in. Each row of this TPM sums to one. Another feature of this TPM is that it does not have any absorbing states (i.e., a state that, once entered, cannot be left). Because of this, this TPM is an *irreducible stochastic matrix*, which has some pretty unique properties. Perhaps the most useful in our case is that you can estimate a stationary probability vector, <span>$$\boldsymbol{\pi}$$</span>, which describes the long-term distribution of the TPM (i.e., the expected proportion of time each state should be visited). In terms of autologistic multistate occupancy models, this stationary probability vector describes the expected occupancy of each of the three states!

How do you calculate <span>$$\boldsymbol{\pi}$$</span>? You need to solve this equation:

$$\boldsymbol{\pi}TPM = \boldsymbol{\pi}$$

This equation, as it turns out, looks almost identical to the equation for left eigenvectors:

$$x A = \lambda x$$

where *x*, our left eigenvector, is a non-zero row vector that satisfies the equation above, *A* is a *n* by *n* matrix (the TPM), and <span>$$\lambda$$</span> is an eigenvalue, which is a special scalar associated to a system of linear equations that ensures the equality above holds. 

For an *n* by *n* matrix we have *n* eigenvectors and *n* eigenvalues. However, in our case we are pretty much only interested in the first of each. This brings me to another unique property of TPMs: the leading eigenvalue of a TPM is always 1. This means that to calculate the left eigenvectors for a TPM we don't need to worry about the eigenvalue <span>$$\lambda$$</span>, and can just solve the equation <span>$$$$\boldsymbol{\pi}TPM = \boldsymbol{\pi}$$</span>. Here is how we can do this in R:


```R
# Our TPM
tpm <- matrix(
  c(0.3,0.4,0.3,
    0.6,0.2,0.2,
    0.1,0.5,0.4),
  ncol = 3,
  nrow = 3,
  byrow = TRUE
)

# Note: base::eigen() gives right eigenvectors, so we need
#         to calculate their inverse to get the left
#         eigenvectors.

# base::eigen() gives right eigenvectors, which satisfy the equation
#   A %*% x = lambda %*%x (which is not the equation we are trying 
#   to solve). However, we can invert the right eigenvector matrix to
#   get our left eigenvectors. The output of this function is a named list
#  with (eigen)'values' and (eigen)'vectors'.
right <- eigen(tpm)

# See that the first eigenvalue of our TPM is 1
round(right$values,2)

[1]  1.00 -0.16  0.06

# Calculate the generalized inverse of the right eigenvectors
#  to get the left eigenvectors
left_vectors <- MASS::ginv(right$vectors)

# Check really quickly that the decomposition of our
#  matrix is correct, this should return the tpm (it does).
right$vectors%*%diag(right$values)%*%left_vectors

     [,1] [,2] [,3]
[1,]  0.3  0.4  0.3
[2,]  0.6  0.2  0.2
[3,]  0.1  0.5  0.4

# All we need to do is normalize the first row of
#   the left eigenvectors (i.e., divide each value
#   by the sum of the first row). To get the
#   stationary distribution.
stationary_dist <- left_vectors[1,]/sum(left_vectors[1,])

# Check out the stationary distribution,
#   in this case the amount of time spent in
#   each state is about 0.33.
stationary_dist
[1] 0.3486239 0.3577982 0.2935780

# These values fulfill the equation pi %*% tpm = pi
as.numeric(stationary_dist %*% tpm) 
[1] 0.3486239 0.3577982 0.2935780
```

Now, you can write a function to do this for you if you'd like, but you can also get the stationary distribution using `markovchain::SteadyStates()` just as easily.

```R
library(markovchain)

# Our TPM
tpm <- matrix(
  c(0.3,0.4,0.3,
    0.6,0.2,0.2,
    0.1,0.5,0.4),
  ncol = 3,
  nrow = 3,
  byrow = TRUE
)
# Make tpm a markovchain 
markov_tpm <- new(
  "markovchain",
  states = c("a", "b", "c"),
  transitionMatrix = tpm,
  byrow = TRUE
)
# The stationary distribution from the markovchain package
stationary_distmc <- markovchain::steadyStates(markov_tpm)

# our stationary distribution from "scratch"
stationary_dist
[1]  0.3486239 0.3577982 0.293578

# our stationary distribution from markovchain package
#   They are equal.
stationary_distmc
             a         b        c
[1,] 0.3486239 0.3577982 0.293578
```

There are two reasons I am took the time to introduce TPMs. First, we are going to use TPMs in our autologistic occupancy model, and because of this we will want to derive the expected occupancy of each state based on the estimated TPM. Likewise, as we will include spatial covariates in this model, our TPM will also spatially vary (i.e., we will calculate <span>$$\pi_s$$</span>, not <span>$$\pi$$</span>). Therefore we need to derive the stationary distribution along whatever environmental gradient we include if we are interested to plotting out the expected occupancy of different states (which we probably want to do).

Second, you may not have known it. But you may have already used a similar approach while fitting dynamic occupancy models. As a reminder, dynamic occupancy models estimate local colonization (<span>$$\gamma$$</span>) and extinction (<span>$$\epsilon$$</span>) rates, and you can use this equation to calculate the expected occupancy of your species:

$$E(\psi) = \frac{\gamma}{\gamma + \epsilon}$$

This equation provides the same solution as our approach above for a 2x2 TPM. Instead of using this equation, we could still use the left eigenvector approach from above.

```R
# colonization
gamma <- 0.8
# extinction
epsilon <- 0.6

# calculate expected occupancy
ex_occ <- gamma / (gamma + epsilon)

# and tack on the complimentary probability,
#  which is the proportion of landscape the 
#  species is absent
ex_occ <- c(ex_occ, 1 - ex_occ)

# give names to vector
names(ex_occ) <- c("present", "absent")


# make a tpm for a dynamic 
#  occupancy model, each row
#  sums to 1
dynamic_tpm <- matrix(
  c(1 - epsilon, epsilon,
    gamma, 1 - gamma),
  ncol = 2,
  nrow = 2,
  byrow = TRUE
)

# Calculate steady states using the
#  markovchain package
dynamic_tpm <- new(
  "markovchain",
  states = c("present", "absent"),
  transitionMatrix = dynamic_tpm,
  byrow = TRUE
)

ex_occmc <- markovchain::steadyStates(dynamic_tpm)

# compare outputs. From markovchain
ex_occmc
  present    absent
0.5714286 0.4285714

# from gamma / (gamma + epsilon)    
ex_occ
  present    absent 
0.5714286 0.4285714 

# They are the same!
```

While it's relatively simple to derive the expected occupancy from a 2x2 TPM (i.e., a dynamic occupancy model), you will need to use the left eigenvector approach when dealing with larger TPMs (which autologistic multistate occupancy models have). As such, this is really helpful to know! Let's move on to our models now.

<p><a href="#top" style>Back to top ⤒</a></p>

## An autologistic multistate occupancy model that uses the logit link

All autologistic multistate occupancy models require two components:

1. A TPM for the latent state model
2. A detection probability matrix for the data model

However, like I demonstrated in the last post on multistate models, there are different ways you can specify the probabilities for each row vector of the latent state TPM. What I did not tell you however, is that the way you write the latent state probabilities has an influence on what link function you will use to convert your linear predictors to probabilities. For this model, let's assume there is no temporal variation in occupancy (i.e., no temporal covariates), and let the latent state TPM, <span>$$\boldsymbol{\psi_{s}}$$</span>, be:

$$
\boldsymbol{\psi_{s}} = 
\begin{bmatrix}
1 - \omega_{s} & \omega_{s} (1 - \delta_{s}) & \omega_{s} \delta_{s})  \\
1 - \omega_s^{2,3|2,3} & \omega_s^{2,3|2,3} (1 - \delta_{s}) & \omega_s^{2,3|2,3} \delta_{s} \\
1 - \omega_s^{2,3|2,3} & \omega_s^{2,3|2,3} (1 - \delta_{s}^{3|3}) & \omega_s^{2,3|2,3} \delta_{s}^{3|3}
\end{bmatrix}
 $$


The probabilities in <span>$$\boldsymbol{\psi_{s}}$$</span> can be defined as:

1. <span>$$\omega_{s}</span>: The probability of occupancy regardless of breeding status.
2. <span>$$\delta_s$$</span>: The conditional probability of breeding given presence.
3. <span>$$\omega_s^{2,3|2,3}$$</span>: The conditional probability of occupancy regardless of breeding status given that the species was present in the last time step.
4. <span>$$\delta_{s}^{3|3}$$</span>: The conditional probability of breeding given the species was present and breeding in the last time step.

For the first time step, we have no knowledge of species presence or breeding status. Because of this, all we need to do is use the first row of <span>$$\boldsymbol{\psi_{s}}$$</span> to estimate the distribution of each state at *t=1* such that:

$$z_{s,t=1} \sim \text{Categorical}(\boldsymbol{\psi_{s,1}})$$

Where <span>$$z_{s,t}$$</span> is a Categorical random variable that can take a value of either 1, 2, or 3 depending on what state site *s* is in at time *t*. Going through the rest of this equation, it's important to remember that <span>$$\boldsymbol{\psi_{s}}$$</span> is a *S* x 3 x 3 array. Therefore, that additional indexing of 1 in <span>$$\boldsymbol{\psi_{s,1}}$$</span> indicates that we selected the first row of the that site-specific TPM, which would be:

$$\boldsymbol{\psi_{s,1}} =  \begin{bmatrix}(1-\omega_s) & \omega_s(1-\delta_s) & \omega_s\delta_s\end{bmatrix}$$


After we estimate occupancy at *t=1*, we can use the latent state at <span>$$z_{s,t-1}$$</span> to index the appropriate row of <span>$$\boldsymbol{\psi_{s}}$$</span> for time *t*. Thus for *t>1* the model is:

$$z_{s,t} \sim \text{Categorical}(\boldsymbol{\psi_{s,z_{s,t-1}}}), t>1$$

So, for example, if <span>$$z_{s,t-1} = 2$$</span>, we grab the 2nd row of our TPM for site *s*. 

Now, the way I wrote out the TPM is a little unique. We have written the probability of the species presence states (i.e., 2 or 3) as a product of two probabilities (e.g., <span>$$\omega_{s} (1 - \delta_{s})$$</span>). In my last post on multistate models, I showed that all of the probabilities in this TPM can take any value between 0 and 1, and each row will still sum to 1. Because of this quality, we can use the logit link instead of the softmax function to specify the linear predictors for each probability. Why can we do that? Because the logit link maps a linear predictor back to a probability between 0 and 1. We just need to specify our logit-linear predictors, convert them back to probabilities, and fill in our TPM. So, for *n* in 1,...,*N* covariates in the *S* by *N* matrix *X* for <span>$$\omega$$</span> and *r* in 1,...,*R*  covariates in the *S* by *R* matrix *U* for <span>$$\delta$$</span>, the logit-linear predictors of the four unique probabilities in this TPM are

$$
\begin{eqnarray}
\text{logit}(\omega_s) &=& \boldsymbol{\beta}^{2,3} \boldsymbol{x}_s  &=& \beta_1^{2,3} \times x_{s,1} + \cdots + \beta_N^{2,3} \times x_{s,N}  \nonumber\\
\text{logit}(\omega_s^{2,3|2,3}) &=& \boldsymbol{\beta}^{2,3} \boldsymbol{x}_s + \theta_1  &=& \beta_1^{2,3} \times x_{s,1} + \cdots + \beta_N^{2,3} \times x_{s,N} + \theta_1 \nonumber \\
\text{logit}(\delta_s) &=& \boldsymbol{\beta}^3 \boldsymbol{u}_s &=& \beta_1^{3} \times u_{s,1} + \cdots + \beta_R^3 \times u_{s,R}  \nonumber \\
\text{logit}(\delta_s^{3|3}) &=& \boldsymbol{\beta}^3 \boldsymbol{u}_s + \theta_2 &=& \beta_1^{3} \times u_{s,1} + \cdots + \beta_R^3 \times u_{s,R} + \theta_2 \nonumber 
\end{eqnarray}
$$

Where the first column of both *X* and *U* is a vector of 1's to accomodate the model intercepts. I've indexed each of the regression coefficients by the states they apply to. For example, <span>$$\boldsymbol{\beta}^{2,3}$$</span> are regression coefficients associated to the probability of presences ignoring breeding status (i.e., states 2 or 3). Conversely, <span>$$\boldsymbol{\beta}^3$$</span> are regression coefficients associated to the conditional probability of breeding given species presence (i.e,. state 3). Looking at the equations above, the only thing that differs between the linear predictors in <span>$$\omega_s$$</span> and <span>$$\omega_s^{2,3|2,3}$$</span> is the latter linear predictor has the autologistic term <span>$$\theta_1$$</span>. Likewise, the only thing that differs between the linear predictors in <span>$$\delta_s$$</span> and <span>$$\delta_s^{3|3}$$</span> is the autologistic term <span>$$\theta_2</span>. Thus, this specification of an autologistic multistate occupancy model only adds two parameters to the model in order to account temporal dependence in species presence or breeding status from one time period to the next.

The data model for an autologistic multistate occupancy model is almost identical to a static multistate occupancy model, except it is indexed by sites (*s* in 1,..,*S*), surveys (*j* in 1,...,*J*), and time periods (*t* in 1,...,*T*). Again, we are going to assume there is no temporal variability in this level of the model (i.e, across *t*), but adding that in is not hard to do. Thus, let the detection matrix of the data model be:

$$
\boldsymbol{\eta}_{s,j} = \begin{bmatrix}
 1 & 0 & 0 \\
 (1 - \rho_{s,j}^{\omega}) & \rho_{s,j}^{\omega} & 0 \\
 (1 - \rho_{s,j}^{\omega}) & \rho_{s,j}^{\omega} (1 - \rho_{s,j}^\delta) & \rho_{s,j}^{\omega} \rho_{s,j}^\delta
 \end{bmatrix} 
 $$

and the two detection probabilities in this matrix can be defined as:

1. <span>$$\rho_{s,j}^{\omega}$$</span>: The probability of detecting the species regardless of breeding status.
2. <span>$$\rho_{s,j}^{\delta}$$</span>: The conditional probability of detecting breeding given the species is present.

Again, the rows sum to 1 in this matrix, the rows correspond to the state site *s* is in at time *t* and the columns correspond to the probability of detecting each state given the row. For example, the only state that can be detected if a site is in state 1 (species absence) is state 1. 



