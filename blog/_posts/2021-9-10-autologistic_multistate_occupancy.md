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
2. [An autologistic multistate occupancy model that uses the logit link]
  1. [Fitting the model in JAGS]
  2. [Plotting out your results]
3. [An autologistic multistate model occupancy that uses the softmax function]
  1. [Fitting the model in JAGS]
  2. [Plotting our your results]


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

It has three states, and the probability of transitioning from one state to another in a single time unit varies based on what state you are currently in. This matrix does not have any absorbing states (i.e., a state that, once entered, cannot be left). Because of this, this TPM is an irreducible stochastic matrix, which has some pretty unique properties. Perhaps the most useful in our case is that you can estimate a stationary probability vector, <span>$$\boldsymbol{\pi}$$</span>, which describes the distribution of the TPM after a "sufficently long time". In terms of autologistic multistate occupancy models, this stationary probability vector describes the expected occupancy of each of the three states!

How do you calculate <span>$$\boldsymbol{\pi}$$</span>? You need to solve this equation:

$$\boldsymbol{\pi}TPM = \boldsymbol{\pi}$$

This equation happens to have a lot of similarity to the equation for left eigenvectors:

$$x A = \lambda x$$

where *x*, our left eigenvector, is a row vector that satisfies that above equation, *A* is a matrix (the TPM), and <span>$$\lambda$$</span> is an eigenvalue, which is a special scalar associated to a system of linear equations and tells you how variable the matrix *A* is around the eigenvector. Conceptually, eigenvectors have somewhat similar properties to the variance of your data around the best fit line of a linear regression, the an eigenvector is the more spread out the matrix is from the eigenvector.

Now, this brings me to another unique property of TPMs, namely that the leading eigenvalue of a TPM is always 1. This means that to calculate the left eigenvectors for a TPM we don't need to worry about the eigenvalue <span>$$\lambda$$</span>, and can just solve the equation <span>$$$$\boldsymbol{\pi}TPM = \boldsymbol{\pi}$$</span>. Here is how we can do this in R:


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

# Get the eigenvectors of the tpm. This is named list
#  with (eigen)'values' and (eigen)'vectors'.
right <- eigen(tpm)

# See that the first eigenvalue of our TPM is 1
round(right$values,2)

[1]  1.00 -0.16  0.06

# Calculate the generalized inverse of the right eigenvectors
#  to get the left eigenvectors
left_vectors <- MASS::ginv(right$vectors)


# Check really quickly that the decomposition of our
#  matrix is correct, this should return the tpm
right$vectors%*%diag(right$values)%*%left_vectors


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

There are two reasons I am walking you through TPMs. First, we are going to use TPMs in our autologistic occupancy model, and because of this we will want to derive the expected occupancy of each state based on the estimated TPM. Likewise, as we will include spatial covariates in this model, our TPM will also spatially vary. Therefore we need to derive the stationary distribution along whatever environmental gradient we include if we are interested to plotting out the expected occupancy of different states (which we probably want to do).

Second, you may not have known it. But you may have already used a similar approach while fitting dynamic occupancy models! Dynamic occupancy models estimates local colonization (<span>$$\gamma$$</span>) and extinction (<span>$$\epsilon$$</span>) rates, and you can use this equation to calculate the expected occupancy of your species

$$E(\psi) = \frac{\gamma}{\gamma + \epsilon}$$

This equation provides the same solution as our approach above for a 2x2 TPM. Instead of using this equation, we could use the left eigenvector approach from above.

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










So, after we estimate occupancy at *t=1* (we'll get to that in a second), we can use the latent state at <span>$$z_{s,t-1}$$</span> to index the appropriate row of <span>$$\boldsymbol{\psi_{s,t}}$$</span> while modeling <span>$$z_{s,t}$$</span>. Hopefully, if you read my previous post on multistate models, you can recall that we used the latent state *z* to index the appropriate row of the data model detection probability matrix. We are just using similar logic here for the latent state model as well!