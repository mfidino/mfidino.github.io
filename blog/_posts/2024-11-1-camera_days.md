---
layout: post
title: Please reconsider including the number of sampling days as a covariate when estimating detection probabilities in occupancy models
category: blog
---

When you flip a fair coin, the probability that you get heads at least once increases with the number of times you flip the coin. This makes intuitive sense. More chances means you are more likely to succeed at least once. If you wanted to describe this mathematically the easiest way to do so is by calculating what is called the 'complement' of the probability you are interested in. If you have probability <span>$$p$$</span>, which in this example is the probability of getting heads in one coin flip, it's complement is <span>$$1-p$$</span>. If we can calculate the probability we don't get heads in EVERY coin flip we do, it's complement represents the probability we get at least one heads. Thus, the probability you get at least one success in  <span>$$J$$</span> coin flips ends up using two complementary probabilities:

$$
P(\text{at least one success}) = 1 - (1 - p)^J
$$

where <span>$$p$$</span> is the probability of success (i.e., getting heads on a coin flip). As we are flipping a fair coin (i.e., <span>$$p = 0.5$$</span>) we can plug in various values for <span>$$J$$</span> to see how the probability of getting at least one success increases with the number of trials.

```R
p <- 0.5
J <- 1:10

pr_one_success <- 1 - (1 - p)^J

plot(
  pr_one_success ~ J,
  xlab = "Number of coin flips",
  ylab = "Pr(at least one flip is heads)",
  type = "l",
  bty ="l",
  las = 1,
  lwd = 3
)

```

Plotted out, you can see that as the number of coin flips increases, so too does the probability you get at least one heads.

![The probability of getting at least one heads with j coin flips]({{site.url}}/blog/images/cd_coin_flip.jpeg#center)

Now, it may seem like the example above has nothing to do with occupancy modeling. We're flipping coins, not deploying camera traps, setting up acoustic recording devices, or conducting some other type of ecological surveys. However, as I review a fair number of papers that use occupancy models for their analysis, one common trend I've encountered is that people will often think they need to control for the fact that the number of surveys between sites differs and they do so by including the number of survey days as a covariate in the detection part of the model. As occupancy models estimate survey-specific detection probabilities, including the number of survey days is essentially stating that as the number of survey days increases, the survey-specific detection probability changes. To bring it back to coin flips, this is like saying as the number of times you flip a coin increases, the probability you get a heads on your next flip, or any flip after that, also changes. Can you see what I am getting at here? Whether it is coin flips or survey-specific detection probabilities, the probability you get at least one success increases with the number of trails but **the individual probability of a single trial does not change**.

So, the question is then what should you do if you want to account for a variable number of survey days in your occupancy model? 

ABSOLUTELY NOTHING. Occupancy models already account for missing data, and they always have!

I'm going to go ahead and demonstrate this in two ways. First, I'll go over some of the math. Following that, I'll simulate a dataset and analyze it using the unmarked package in R. After that, I am going to circle back to why I wrote this blog post to begin with and consider if there could be a reason you may want to include the number of survey days in the detection part of your model.


### The math ####

If you take a second and look at the MacKenzie et al. (2002) [occupancy modeling paper](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/0012-9658(2002)083[2248:ESORWD]2.0.CO;2) that introduces this class of model, they demonstrate how to construct the likelihood of a given detection history (i.e., whether or not you detected your species of interest on each survey at a single location). Take, for example, the detection history [1, 0, 1]. We detected the species on survey 1 and 3, but failed to detect it on survey 2. If <span>$$\psi$$</span> is the occupancy probability and <span>$$\rho$$</span> is the conditional probability of detecting the species on a given survey given they are present, then we can write out the overall likelihood of this detection history using those probabilities and their complements.

$$
\psi\rho_1(1 - \rho_2)\rho_3
$$

On the likelihood above the subscripts on <span>$$\rho$$</span> just denote the survey. Because we detected the species at least once we know it is there, so we must include <span>$$\psi$$</span> in the likelihood. From there, we use either <span>$$\rho$$</span> or <span>$$1 - \rho$$</span> for each survey, depending on if the species was detected or not. Let's consider another example. If we had the detection history [0, 0, 0] then either the species is not there or it is and we failed to detect it. As such, the likelihood of this detection history is a bit more complex because we need to account for both.

$$
(1 - \psi) + \psi(1 - \rho_1)(1 - \rho_2)(1 - \rho_3)
$$

Again, either the species is not there (<span>$$\psi$$</span>) or it was and we did not detect it (<span>$$\psi(1 - \rho_1)(1 - \rho_2)(1 - \rho_3)$$</span>). If you are wondering why we add these probabilities together to generate this likelihood (I know I did at one time), it is because they are two mutually exclusive events (the species is either there or not there). Because of that we can use the addition rule in probability. For mutually exclusive events (<span>$$A$$</span>) and (<span>$$B$$</span>), the probability of either occurring is <span>$$Pr(A \textrm{ or } B) = Pr(A) + Pr(B)$$</span>. Looks similar to the likelihood for the [0, 0 ,0] detection history right?

Anyways, the MacKenzie et al. (2002) paper also demonstrates how to account for missing data and it's incredibly simple. If you missed a survey you just drop that survey-specific detection probability from the likelihood. Because of this, detection histories with missing data contribute less to the overall likelihood of the model because they have less data. As such, the model already accounts for the fact that the number of surveys differs! For example, the likelihood of the detection history [0, 1, NA], where NA means we did not sample on the third survey is:

$$
\psi(1 - \rho_1)\rho_2
$$

You can see above that <span>$$\rho_3$$</span> has been dropped from the likelihood. We didn't conduct that survey so we don't add it to the likelihood. That is all you need to do! But how do actually do that when it comes to modeling in R?

### In R ####

The good news is that every occupancy modeling R package I've encountered already provides a way for you to account for missing data, which is often through adding NA values to a detection history. Let's see what that looks like with the unmarked package for a single-species, single-season model. In this cherry picked example I have a very large sample size, which makes missing an observation here or there less of an issue. If you took this code and changed the number of sites to say 100, you would find that the model gives some pretty bogus estimates. Why? Mostly because occupancy is relatively low and having a survey-specific detection probability of 0.2 means that, at best we have a roughly 60% chance of detecting the species if it was present (i.e., (1 - (1 - 0.2)^4) = 0.594). As such, don't take this simulation as a general rule that you need a whole bunch of sites when you have some missing data, it is going to depend a whole lot on your species and sampling process!


```R
# Load unmarked
library(unmarked)

# set seed for reproducibility
set.seed(8675309)

# Probability of occupancy and survey-specific detection probability
psi <- 0.35
rho <- 0.2

# Number of sites sampled
n <- 150

# maximum number of surveys per site
max_j <- 4

# Number of surveys per site, sampling because it is variable. We
#  are assuming here that there are just a few sites we sampled only
#  one time and then the rest were either 2, 3, or 4 surveys.
j <- sample(
  x = 1:max_j,
  size = n,
  replace = TRUE,
  prob = c(0.1, 0.3, 0.3, 0.3)
)

# Simulate species presence.
z <- rbinom(
  n = n,
  size = 1,
  prob = psi
)

# For each site, generate a detection history, which will
#  vary based on the number of surveys per site. This can
#  be done without a for loop, but I think it is a bit
#  easier to understand with one. Filling y with NA
#  as that is what unmarked uses to denote a missed
#  observation.

y <- matrix(
  NA,
  ncol = max_j,
  nrow = n
)

for(i in 1:n){
  y[i,1:j[i]] <- rbinom(
    n = j[i],
    size = 1,
    prob = rho * z[i]
  )
}

# Look at the first few rows to see the variable detection histories
head(y)

#       [,1] [,2] [,3] [,4]
# [1,]    0    0   NA   NA
# [2,]    0    0    0   NA
# [3,]    0    0    0    0
# [4,]    0    0    0    0
# [5,]    0    0   NA   NA
# [6,]    0    0    0    0

# Prepare detection history for unmarked. We have
#  no covariates so this is very simple.
my_occu_frame <- unmarked::unmarkedFrameOccu(
  y = y
)

# Fit the single-season occuapncy model
my_model <- unmarked::occu(
  formula = ~1~1,
  data = my_occu_frame
)



# Convert the logit-scale intercept values back to 
#  occupancy and survey specific detection probability.
psi_ci <- unmarked::confint(
  my_model,
  type = "state"
)

rho_ci <- unmarked::confint(
  my_model,
  type = "det"
)

# Back transform to probability scale, round them a little bit
#  and add on 'truth' as we simulated the data.
my_estimates <- data.frame(
  parameter = c("psi", "rho"),
  estimate = plogis(
    coef(my_model)
  ),
  lower_ci = plogis(
    c(
      psi_ci[1],
      rho_ci[1]
    )
  ),
  upper_ci= plogis(
    c(
      psi_ci[2],
      rho_ci[2]
    )
  ),
  truth = c(
    psi, rho
  )
)

my_estimates <- as.data.frame(
  lapply(
    my_estimates,
    function(x){
      if(is.numeric(x)){
        round(x,2)
      } else {
        x
      }
    }
   )
)

# And look, in this cherry picked simulated dataset with a large sample
#  size, the model estimates have the truth within their 95% CI!

my_estimates
#   parameter estimate lower_ci upper_ci truth
# 1       psi     0.40     0.18     0.67  0.35
# 2       rho     0.19     0.09     0.36  0.20

```

### Is there any reason we should add number of survey days to the detection part of the model? ###

While writing this I started thinking about whether there would be any valid reason to include the number of survey days in your model. Realistically, I think the answer is no, more often then not you should have some other covariate that better captures whatever pattern you are trying to capture. For example, with camera trapping, having fewer surveys could be because the memory card on your camera filled up due to grass blowing in the wind. As such, if a camera is not sampling correctly over your surveys you probably want to account for that, and including survey days may be one way to do it. Conversely, I would argue that the better option would be to include something like the total number of photos collected over a survey in linear and quadratic terms (i.e., photos + photos^2) because we'd likely expect detection probability to go up and then start decreasing once you get a site with too many photos due to grass blowing in the wind. Or, if you tracked that while tagging your data, you may be able to include an 'over trigged' covariate into the model that equals 1 during surveys you are getting a lot of grass blowing in the wind. This is, of course, just one example that comes to mind. As always, it is important to carefully consider why you are doing something and not just listen to some person (e.g., me) tell you that 'it must be done in this way.' More often then not, including the number of survey days in the detection part of the model does not make a lot of sense even though you will often find significant differences with it. But if you have a valid reason for it, just make sure you justify it and provide some logic behind why you are doing it! It's as 'simple' as that.

