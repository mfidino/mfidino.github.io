---
layout: post
title: Linear regression&#58 The normal distribution
category: blog
---

Regression analysis is a widely used set of statistical tools that can estimate the relationship between different variables. To better understand regression analysis, I think it is important to learn about some common probability distributions that are used to fit different types of models to data. Without getting too technical, probability distributions are just mathematical functions that provide the probability of occurrence for different outcomes of an experiment, and we use these distributions because they have proven to be helpful. Probability distributions have parameters, which are numeric inputs that can be used to change how likely different outcomes are for a specific experiment. In the real world, we generally do not know what these parameter values are, which is why we collect data to estimate them.

Arguably, one of the best-known probability distributions is the normal distribution, which is one of a few different distributions that resemble a bell curve:

INCLUDE PLOT OF NORMAL DISTRIBUTION HERE.

And often follows this notation:

$$Normal(\mu, \sigma)$$

The normal distribution consists of two parameters, <span>$$\mu$$</span> and <span>$$\sigma$$</span>, which respectively represent the mean and standard deviation of this probability distribution. In more general terms, <span>$$\mu$$</span> is a location parameter, which often indicates the mean, median, or mode of a statistical distribution while <span>$$\sigma$$</span> is a scale parameter. Scale parameters indicate the spread of a distribution and the larger the scale parameter the more spread out the distribution. 

While it is possible for <span>$$\mu$$</span> to be any number between <span>$$-\infty$$</span> and <span>$$\infty$$</span>,  <span>$$\sigma \geq 0$$.

Depending on the type of data or what the experiment is, you may be encouraged to use different types of probability distributions.

These probability distributions have been found to be useful in predicting different types of   

At it's most, a probability distribution is a mathematical function that can provide the probability of occurrence 