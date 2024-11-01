---
layout: post
title: An autologistic version of the Rota et al. (2016) multi-species occupancy model for interacting species so you can use data collected across multiple sampling periods!
category: blog
---

I've gotten quite a few emails about people interested in fitting the Rota et al. (2016) model to their data. Nine times out of ten, however, there is one issue that prevents them from being able to easily do this, which is that they have they collected data over time. With a large sample size, you could fit the dynamic version of that model, but often times a researcher either may not have a huge sample or they may not be all that interested in estimating the dynamics (i.e., colonization and extinction). 

However, there are actually two relatively easy ways to get around this. Either you can incorporate a random site effect into your model or you can reparameterize it so that you add first-order autologistic parameters. Given this post is about the second option, we are going to ignore the random effect route (even though random effects are not currently available for the static model in `unmarked` ). Autologistic parameters trigger when something happened in the past, and a first-order autologistic parameter is one that depends on the time period immediately before it. So, in our case, if a species is present at a site at time t-1 then we should include our autologistic term into the model, otherwise it becomes zero. To see an example of this for a standard occupancy model, check my [blog post here](https://masonfidino.com/autologistic_occupancy_model/).


 For this example we are going to use the Categorical distribution in our model, which means we will represent community states as different numbers. For a two-species system we have four possible states: (1) no species present, (2) species A present, (3) species B present, and (4) species A and B present. However, if we need to turn those autologistic parameters of and on depending on the community state in the previous time period, we have to convert these numeric states to binary variables. The easiest way to do that is to add a matrix to the model that can be used as an indicator variable. This matrix will have a number of rows equal to the number of community states and number of columns equal to the number of species. So, for two species this would be

 $$
\text{indicator matrix} = \begin{bmatrix}
0 & 0 \\
1 & 0 \\
0 & 1 \\
1 & 1
 \end{bmatrix} 
 $$

With this added to the model, you can easily 'flip the switch' on your autologistic terms if your species is present in the previous timestep.

## The model

This model is essentially identical to the Rota et al. (2016) model except we add one additional parameter per species (the autologistic terms). If you really wanted to you could create different autologistic terms for every community state, but for this example we will start as simple as possible. For two species at <span>$$s$$</span> in <span>$$1,\dots,S$$</span> sites and <span>$$t$$</span> in <span>$$1,\dots,T$$</span> sampling seasons we will have 


$$
\begin{eqnarray}
\psi_1 &=&\frac{\beta^{\psi_1}}{\beta^{\psi_1} + \beta^{\psi_2} + \beta^{\psi_3} + \beta^{\psi_4}} = \frac{1}{1 + \text{exp}(a_0) + \text{exp}(b_0) + \text{exp}(a_0 + b_0 + c_0)} \\
\psi_2 &=&\frac{\beta^{\psi_2}}{\beta^{\psi_1} + \beta^{\psi_2} + \beta^{\psi_3} + \beta^{\psi_4}} = \frac{\text{exp}(a_0)}{(1 + \text{exp}(a_0) + \text{exp}(b_0) + \text{exp}(a_0 + b_0 + c_0)} \\
\psi_3 &=&\frac{\beta^{\psi_3}}{\beta^{\psi_1} + \beta^{\psi_2} + \beta^{\psi_3} + \beta^{\psi_4}} = \frac{\text{exp}(b_0)}{(1 + \text{exp}(a_0) + \text{exp}(b_0) + \text{exp}(a_0 + b_0 + c_0)} \\
\psi_4 &=&\frac{\beta^{\psi_5}}{\beta^{\psi_1} + \beta^{\psi_2} + \beta^{\psi_3} + \beta^{\psi_4}} = \frac{\text{exp}(a_0 + b_0 + c_0)}{(1 + \text{exp}(a_0) + \text{exp}(b_0) + \text{exp}(a_0 + b_0 + c_0)}
\end{eqnarray}
$$




