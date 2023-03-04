---
layout: post
title: A really helpful function I wrote to make the MCMC output from NIMBLE or JAGS so much easier to work with for model predictions and plotting
category: blog
---

When fitting complex models in `NIMBLE` or `JAGS` you can end up with hundreds or even thousands of parameters. Take, for example, a multi-species occupancy model. In `NIMBLE` that would look something like this:

```R
multi_species_occupancy <- nimble::nimbleCode(
  {
    # Latent state part of model
    for(species in 1:nspecies){
      for(site in 1:nsite){
        # logit-linear predictor latent state
        logit(mu_psi[site,species]) <- inprod(
          psi_beta[species,1:nbeta_psi],
          psi_design_matrix[site,1:nbeta_psi]
        )
        z[site,species] ~ dbern(
          mu_psi[site,species]
        )
      }
    }
    # Detection part of model
    for(species in 1:nspecies){
      for(site in 1:nsite){
        for(visit in 1:nvisit){
          # assuming covariates for sites and visits
          logit(mu_rho[site,species,visit]) <- inprod(
            rho_beta[species, 1:nbeta_rho],
            rho_design_matrix[site, visit, 1:nbeta_rho]
          )
          y[site,visit,species] ~ dbern(
            mu_rho[site,species,visit] * z[site,species]
          )
        }
      }
    }
    # Priors for among-species average
    for(ipsi in 1:nbeta_psi){
      psi_beta_mu[ipsi] ~ dnorm(0, sd = 1)
      psi_beta_sd[ipsi] ~ dgamma(1, 1)
    }
    for(irho in 1:nbeta_rho){
      rho_beta_mu[irho] ~ dnorm(0, sd = 1)
      rho_beta_sd[irho] ~ dgamma(1, 1)
    }
    # species-specific parameters
    for(ipsi in 1:nbeta_psi){
      for(species in 1:nspecies){
        psi_beta[species,ipsi] ~ dnorm(psi_beta_mu[ipsi], sd = psi_beta_sd[ipsi])
      }
    }
    for(irho in 1:nbeta_rho){
      for(species in 1:nspecies){
        rho_beta[species,irho] ~ dnorm(rho_beta_mu[irho], sd = rho_beta_sd[irho])
      }
    }
  }
)
```

From this model we would want to track the species-level parameters (`psi_beta` and `rho_beta`), as well as the among-species parameters (`psi_beta_mu`, `rho_beta_mu`, `psi_beta_sd`, and `rho_beta_sd`). When we track these
parameters across multiple chains, most MCMC software will return a list-style object with each element representing a MCMC matrix that has a number of rows equal to the number of samples in the chain and a number of columns equal to the total number of parameters monitored. For example, I fitted the model above to some simulated data for 7 species where `nbeta_psi = 2` and `nbeta_rho = 2`. This resulted in 36 parameters that were tracked

```R
# Assuming we have an object called my_model_ouput
colnames(
	my_model_output$chain1
)

[1] "psi_beta[1, 1]" "psi_beta[2, 1]"
[3] "psi_beta[3, 1]" "psi_beta[4, 1]"
[5] "psi_beta[5, 1]" "psi_beta[6, 1]"
[7] "psi_beta[7, 1]" "psi_beta[1, 2]"
[9] "psi_beta[2, 2]" "psi_beta[3, 2]"
[11] "psi_beta[4, 2]" "psi_beta[5, 2]"
[13] "psi_beta[6, 2]" "psi_beta[7, 2]"
[15] "psi_beta_mu[1]" "psi_beta_mu[2]"
[17] "psi_beta_sd[1]" "psi_beta_sd[2]"
[19] "rho_beta[1, 1]" "rho_beta[2, 1]"
[21] "rho_beta[3, 1]" "rho_beta[4, 1]"
[23] "rho_beta[5, 1]" "rho_beta[6, 1]"
[25] "rho_beta[7, 1]" "rho_beta[1, 2]"
[27] "rho_beta[2, 2]" "rho_beta[3, 2]"
[29] "rho_beta[4, 2]" "rho_beta[5, 2]"
[31] "rho_beta[6, 2]" "rho_beta[7, 2]"
[33] "rho_beta_mu[1]" "rho_beta_mu[2]"
[35] "rho_beta_sd[1]" "rho_beta_sd[2]"

```

You can see above that the indexing gets added to the column names. This is great as we can use it to determine what species the parameter is associated to as well as if the parameter is an intercept or slope term. However, it makes it a bit
difficult to work with the MCMC output, especially when we get around to making predictions from the model. Instead,
it would be awesome if we could make it so that each parameter type (e.g., `psi_beta`, `rho_beta`, etc.) was appropriately indexed in an array instead of storing the array information in the column name itself. A couple months back I wrote a function to do exactly this, and have been using it for the last few analyses I've been working on, and I absolutely love it. I'm regularly sending said function to collaborators and the like, and as a result I figured it would be a good idea to showcase it here so you can use it as well!

The function, `split_mcmc()`, takes in the full MCMC matrix (i.e., you have to join together all of the chains into one big matrix) and then returns a named list object, one name for every parameter type, that is indexed exactly as it would be in your Bayesian model, save for the fact that the first dimension is always your MCMC samples. The code for this 
function is:

```R
split_mcmc <- function(x){
  # get parameter names
  pars <- colnames(x)
  # unique parameters
  unq_pars <- unique(
    gsub(
      "\\[.*\\]",
      "",
      pars
    )
  )
  # make list object to store arrays in
  result_list <- vector(
    "list",
    length = length(unq_pars)
  )
  # name the list
  names(result_list) <- unq_pars
  # fill in the arrays
  for(i in 1:length(result_list)){
    # get just the parameters
    tmp <- pars[grep(
      paste0(
        "^",unq_pars[i], "\\["
      ),
      pars
    )]
    if(length(tmp) == 0){
      tmp <- pars[grep(
        paste0("^",unq_pars[i],"$"),
        pars
      )]
    }
    # and then the array dimensions
    arr_dim <- gsub(
      paste0(
        unq_pars[i],"\\[|\\]"
      ),
      "",
      tmp
    )
    arr_dim <- strsplit(
      arr_dim,
      ","
    )
    ndim <- length(arr_dim[[1]])
    npar <- length(arr_dim)
    # make a matrix
    arr_ind <- suppressWarnings(
      matrix(
        as.numeric(
          unlist(arr_dim)
        ),
        ncol = ndim,
        nrow = npar,
        byrow = TRUE
      )
    )
    if(nrow(arr_ind) == 1 & ncol(arr_ind) == 1){
      arr_ind[1,1] <- 1
    }
    # get max index for each
    max_ind <- apply(arr_ind, 2, max)
    # and then fill in the array
    result_list[[i]] <- array(
      x[,tmp],
      dim = c(nrow(x), max_ind)
    )
    
  }
  return(result_list)
}
```

Let's use that function on the `NIMBLE` output of this model and see what it looks like


```R
# combine all MCMC chains
my_mcmc <- do.call(
  "rbind",
  my_model_output
)

# use split_mcmc()
my_mcmc <- split_mcmc(
  my_mcmc
)

# check out names of list
names(my_mcmc)

[1] "psi_beta"    "psi_beta_mu"
[3] "psi_beta_sd" "rho_beta"   
[5] "rho_beta_mu" "rho_beta_sd"

```

And now let's look at the dimensions of `psi_beta`. Again, we have seven species and two parameters for each species.

```R

dim(my_mcmc$psi_beta)
[1] 10000     7     2

```
The first dimension is for the 10,000 MCMC samples, the second dimension is for species, and the third is for the two parameters associated to that species. So, if I wanted to take a peek at the first few MCMC samples of the intercept and slope for the third species, I can just subset the second dimension of `my_mcmc$psi_beta`!

```R
head(my_mcmc$psi_beta[,3,])
         [,1]      [,2]
[1,] 4.282852 0.6882747
[2,] 4.282852 0.4528713
[3,] 4.282852 0.7064803
[4,] 5.730312 0.2103931
[5,] 4.517122 0.6484852
[6,] 5.473493 0.4883727
```

Other parameter types were not matrices in `NIMBLE`, like the among-species parameters. These elements have fewer dimensions in `my_mcmc`


```R
dim(my_mcmc$psi_beta_mu)
[1] 10000     2

```
Here we just have a matrix with 10,000 MCMC samples and two parameters. I cannot stress how much easier this makes it to work with the model output / make model predictions with. See the plotting part of [this blog post here](https://masonfidino.com/interpret_rota_model/) as an example.

Anways, I just wanted to highlight this function as it has definitely made my life easier when working with more complex outputs from different analyses. And just to demonstrate that it works correctly, here is a little
matrix you can use to test it with. Enjoy!


```R
# testing split_mcmc
test_matrix <- matrix(
  rep(1:10, each = 3),
  ncol = 10,
  nrow = 3
)

# give column names like they would appear in NIMBLE
#  or JAGS.
colnames(test_matrix) <- c(
  "a", "b[1,1]", "b[2,1]", "b[1,2]", "b[2,2]",
  "c[1]", "c[2]", "c[3]", "c[4]", "c[5]"
)

# take a look at this object
test_matrix

     a b[1,1] b[2,1] b[1,2] b[2,2] c[1] c[2] c[3] c[4] c[5]
[1,] 1      2      3      4      5    6    7    8    9   10
[2,] 1      2      3      4      5    6    7    8    9   10
[3,] 1      2      3      4      5    6    7    8    9   10

# use split_mcmc!
test_mcmc <- split_mcmc(
  test_matrix
)

# should be a matrix, three rows, one column, all 1's
test_mcmc$a

      [,1]
[1,]    1
[2,]    1
[3,]    1

# should be a matrix, three rows, two columns, 2 and 3's
test_mcmc$b[,,1]

      [,1] [,2]
[1,]    2    3
[2,]    2    3
[3,]    2    3

# should be a matrix, three rows, two columns, 4 and 5's
test_mcmc$b[,,2]
      [,1] [,2]
[1,]    4    5
[2,]    4    5
[3,]    4    5

# should be amatrix, 3 rows, 5 columns, 6 through 10
test_mcmc$c
      [,1] [,2] [,3] [,4] [,5]
[1,]    6    7    8    9   10
[2,]    6    7    8    9   10
[3,]    6    7    8    9   10

```