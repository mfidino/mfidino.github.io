---
layout: post
title: Linear regression&#58 Using maximum likelihood to estimate the mean and standard deviation of a normally distributed random variable
category: blog
---

When a regression is fit to data in `R`, there is a fair bit of math that occurs *behind the scenes* between when the regression is executed and what `R` formats and returns. Hiding some of this code and math is generally useful as it does not require a person to code common statistical algorithms on their own. However, the math that is hidden away is not only interesting but also important to understand regression.

This set of posts breaks down a bit of the math that may occur when a regression model is fit to data. Starting here, I will discuss how to use maximum likelihood to estimate the mean and standard deviation of a normal distribution from data. Maximum likelihood is an incredibly useful method that determines the most likely values for the parameters in a model given the data. In very simple cases, like the example here, these values can actually be calculated by hand.


### Motivating example

We are interested in the average weight (in grams) of eastern gray squirrels (*Sciurus carolinensis*) at [Humboldt park](https://en.wikipedia.org/wiki/Humboldt_Park_(Chicago_park)) in Chicago, Illinois and how much the weight of squirrels in this specific population deviates from this central tendency. 

After procuring the necessary permits and a conducting a few days of field work, we collected the weight of 10 squirrels. Further, we know that these 10 squirrels do not represent the entire Humboldt park population because many other squirrels eluded our capture. Therefore, these squirrels represent a sample of the overall population.

Mathematically, we will represent the squirrel weight data as the vector <span>$$\mathbf{y}$$</span>. This vector has a length of 10 and each element represents the weight of a squirrel. In `R`, this data vector can be represented as:
```R
# The squirrel weight data
y <- c(557.1, 416.3, 393.6, 459.0, 588.6, 
       503.5, 507.7, 649.3, 529.2, 524.2)
```

For the sake of this example we will assume that these weights are a random sample drawn from a normally distributed population such that: 

$$y_i \sim Normal(\mu, \sigma).$$

Here, <span>$$\mu$$</span> represents the weight of an average squirrel in Humboldt park while <span>$$\sigma$$</span> is the standard deviation or amount of dispersion in the population. If <span>$$\sigma$$</span> is small, than squirrels in this population tend to be much closer to the average. Likewise, <span>$$\sigma$$</span> is large then there can be a fair bit of spread in weights. Currently, we do not know these population level parameters, but we would like to estimate them with the data we collected.

### Estimating this in `R`

In `R`, trying to estimate <span>$$\mu$$</span> and <span>$$\sigma$$</span> would be akin to fitting an intercept only linear regression. This can be done in `R` with the `lm()` function by only including a `1` on the right hand side of the linear predictor:

```R
# Fit the model
squirrel_model <- lm(y ~ 1)

```
By looking at the summary of the model we can determine our maximum likelihood estimates for <span>$$\mu$$</span> and <span>$$\sigma$$</span>, which would respectively be the `(Intercept)` estimate and the `Residual standard error`:

```R
summary(squirrel_model)

Call:
lm(formula = y ~ 1)

Residuals:
    Min      1Q  Median      3Q     Max 
-119.25  -42.73    3.10   37.27  136.45 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)   512.85      24.33   21.08 5.71e-09 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 76.93 on 9 degrees of freedom
```

This model indicates that <span>$$\hat{\mu} = 512.85$$</span>  and <span>$$\hat{\sigma} = 76.93$$</span> (we add the little hats on the parameters because they are estimates). These estimates indicate that your perfectly average squirrel in Humboldt park weighs about 512 grams. Likewise, we would expect 95% of Humboldt park squirrels to roughly weigh between 362 and 663 grams (i.e., <span>$$\hat{\mu} \pm 1.96 \times \hat{\sigma}$$</span>). Knowing this, lets move on to seeing how to estimate these values by hand.

### Doing it by hand

With a little bit of differentiation we can generate point estimates for <span>$$\mu$$</span> and <span>$$\sigma$$</span>. Remember that derivatives can be used to locate where a function is at its maximum value, which sounds like exactly what we need given we are using **maximum** likelihood estimation. All we need now is a function to differentiate that takes parameter values and data and returns a probability. Such functions are called probability density functions, and for the normal distribution, the probability of observing data point <span>$$x$$</span> given a mean <span>$$\mu$$</span> and variance <span>$$\sigma^2$$</span> is:


$$f\left(x | \mu, \sigma^2\right)=\frac{1}{\sqrt{2 \pi\sigma^2}} \textrm{exp}\left({-\frac{(x-\mu)^{2}}{2 \sigma^{2}}}\right)$$

However, this equation only calculates the probability of observing a single data point given <span>$$\mu$$</span> and <span>$$\sigma^2$$</span>. Maximum likelihood estimation requires us to calculate the joint probability of all our data. This joint probability can be written as:


$$f\left(y_{1}, y_{2}, \ldots, y_{n} | \sigma, \mu\right)=\prod_{i=1}^{n} \frac{1}{\sqrt{2 \pi\sigma^2}} \textrm{exp}\left({-\frac{\left(y_{i}-\mu\right)^{2}}{2 \sigma^{2}}}\right)$$

Knowing that we want to differentiate this, there are a couple algebraic manipulations we can use to make this a little easier to tackle.  First, we will want to remove that reciprocal in the first term and apply the product throughout the whole equation. Second, we want to avoid some confusion with that square on the variance term. However, we could represent these parameters with different symbols to get around this.  Therefore, let <span>$$\theta_1 = \mu$$</span> and <span>$$\theta_2 = \sigma^2$$</span>.

$$f\left(y_{1}, y_{2}, \ldots, y_{n} | \theta_1, \theta_2 \right)= \theta_2^{-\frac{n}{2}} \times 2\pi^{-\frac{n}{2}} \times \textrm{exp}\left(\frac{\sum_{i=1}^{n}\left(y_{i}- \theta_1\right)^2}{2 \theta_2}\right)$$

We are getting closer, but with all of these products this function is still tricky to differentiate. However, we can make this function additive by taking the natural logarithm of both sides. For simplicity, let <span> $$\mathcal{L} = f\left(y_{1}, y_{2}, \ldots, y_{n} | \theta_1, \theta_2 \right)$$</span>. Additionally, recall that <span> $$\textrm{log}\left(a^b\right) = b\:\textrm{log}(a)$$</span>.

$$\textrm{log}(\mathcal{L})= -\frac{n}{2}\textrm{log}(\theta_2) - \frac{n}{2}\textrm{log}(2\pi) - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{2 \theta_2}$$


Now this looks much easier to handle. First, we will take the partial derivative of this function with respect to <span>$$\theta_1$$</span> to get the point estimate of the mean. The first two terms do not contain <span>$$\theta_1$$</span>, so they drop out of the derivative. Additionally, we can use the chain rule, <span>$$(f \circ g)^{\prime}=\left(f^{\prime} \circ g\right) \cdot g^{\prime}$$</span>, on the numerator of the third term.

$$\frac{\partial\textrm{log}(\mathcal{L})}{\partial \theta_1}= \ccancel{red}{-\frac{n}{2}\textrm{log}(\theta_2)} \ccancel{red}{- \frac{n}{2}\textrm{log}(2\pi)}- 2\frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)(-1)}{2 \theta_2}$$

Following this we can further simplify:

$$\frac{\partial\textrm{log}(\mathcal{L})}{\partial \theta_1}= \ccancel{Gray}{-}\ccancel{red}{2}\frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)\ccancel{Gray}{(-1)}}{\ccancel{red}{2} \theta_2} = \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)}{ \theta_2}$$

Now all we need to do is set this equation equal to zero and solve for <span>$$\theta_1$$</span>:

$$\begin{aligned}
0 &= \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)}{ \theta_2}\\
 &= \sum_{i=1}^{n}\left(y_{i}-\theta_1\right) \\
 &= \sum_{i=1}^{n}y_{i}- n\theta_1 \\
 n\theta_1 &= \sum_{i=1}^{n}y_{i} \\
 \theta_1 &= \hat{\mu}= \frac{\sum_{i=1}^{n}y_{i}}{n}
\end{aligned}$$

And there we have it. The maximum likelihood estimate for <span>$$\mu$$</span> is the mean of the data. Remember, `R` estimated the intercept of the linear regression to be 512.85. The mean of the squirrel weights should be identical to our model estimate.

```R
# The squirrel weights
y <- c(557.1, 416.3, 393.6, 459.0, 588.6, 
       503.5, 507.7, 649.3, 529.2, 524.2)

# Our maximum likelihood estimate of mu for a normal distribution
mean(y)
[1] 512.85

# The estimated intercept from the linear model
coef(squirrel_model)
(Intercept) 
    512.85 

```

Moving onto the second point estimate, we need to take the partial derivative with respect to <span>$$\theta_2$$</span>. For the first term, remember that <span>$$\textrm{log}(x)^\prime = \frac{1}{x}$$</span>. The second term is a constant and therefore cancels out. And finally we use the reciprocal rule, <span>$$\frac{d}{d x} \frac{1}{f(x)}=-\frac{f^{\prime}(x)}{f(x)^{2}}$$</span>, on the last term to calculate our partial derivative. 

$$
\begin{aligned}
\frac{\partial\textrm{log}(\mathcal{L})}{\partial \theta_2}&= -\frac{n}{2}\textrm{log}(\theta_2) \ccancel{red}{- \frac{n}{2}\textrm{log}(2\pi)} - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{2 \theta_2}\\

&= -\frac{n}{2\theta_2} - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{2 \theta_2^2}
\end{aligned}$$

From here we set this equation equal to zero and solve for <span>$$\theta_2$$</span>.

$$
\begin{aligned}
0 &= \left[-\frac{n}{2\theta_2} - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{2 \theta_2^2}\right] \times 2\theta_2^2\\
&= -\frac{n2\theta_2^2}{2\theta_2} - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}2 \theta_2^2}{2 \theta_2^2}\\
&= -\frac{n\ccancel{red}{2}\theta_2\ccancel{Gray}{^2}}{\ccancel{red}{2}\ccancel{Gray}{\theta_2}} - \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}\ccancel{violet}{2 \theta_2^2}}{\ccancel{violet}{2 \theta_2^2}}\\
&=n\theta_2 - \sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2} \\
n\theta_2 &= \sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2} \\
\theta_2 &= \hat{\sigma}^2 = \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{n}
\end{aligned}
$$

Given that we already know that <span>$$\theta_1$$</span> is the mean of our observed data, the maximum likelihood estimate for <span>$$\sigma^2$$</span> is the residual sum of squares divided by the sample size. We can then take the square of this to calculate the standard deviation, which `R` estimated to be 76.93.

```R
y <- c(557.1, 416.3, 393.6, 459.0, 588.6, 
       503.5, 507.7, 649.3, 529.2, 524.2)


squirrel_variance <- sum((y-mean(y))^2)/ length(y)

squirrel_sd <- sqrt(squirrel_variance)

[1] 72.98432

```
It looks like the residual standard error is not the same as our own maximum likelihood estimate! The reason for this is because the one we derived here is a biased estimator ([link to wikipedia article](https://en.wikipedia.org/wiki/Bias_of_an_estimator#Sample_variance)), and `R` provides the unbiased estimate. To estimate the unbiased population variance or standard deviation we must actually use Bessel's correction ([link to wikipedia article](https://en.wikipedia.org/wiki/Bessel%27s_correction)). This is a pretty simple correction as it just means we modify the denominator of our sample variance estimate from <span>$$n$$</span> to <span>$$n-1$$</span>.

$$
\theta_2 = \hat{\sigma}^2 = \frac{\sum_{i=1}^{n}\left(y_{i}-\theta_1\right)^{2}}{n-1}

$$

```R
# population variance using Bessel's correction
squirrel_variance <- sum((y-mean(y))^2)/ (length(y)-1)

squirrel_sd <- sqrt(squirrel_variance)

# Estimate calculated by hand
squirell_sd
[1] 76.93222

# compare to the estimate from the linear model
summary(squirrel_model)$sigma

[1] 76.93222

```

And there we have it. By taking some partial derivatives of the normal distribution's likelihood function we were able to calculate by hand a little bit of what `R` does behind the scenes. 
