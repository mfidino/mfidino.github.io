---
layout: post
title: A dynamic multistate occupancy model for one species in `JAGS`
category: blog
---

This is now the third post I've written on multistate models, and it will be the last one I write for a single species. To recap, I have covered [static multistate models](https://masonfidino.com/intro_2_multistate_occupancy/) and [autologistic multistate models](https://masonfidino.com/autologistic_multistate_occupancy/) that estimate occupancy across multiple sampling periods. This model here borrows a little bit from those previous posts, so if you've not read those you may want to start from the beginning of this multistate modeling series.

Dynamic occupancy models seperate changes in occupancy status over time as a function of local colonization and extinction rates. Essentially, if a site is not occupied in the previous time step, what is the probability it becomes occupied in the next time step? Likewise, if a site is occupied in the previous time step, what is the probability it is not occupied in the next time step? Extending such a model out to multiple states is a little bit trickier than the static occupancy model, but it's not insurmountable. In fact, I already did most of the heavy lifting by specifying the autologistic multistate model, which uses a transition probability matrix to house the different rates to estimate changes between states over time.

#### Links to different parts of the post

1. [The model](#the-model)
2. [Simulating some data and fitting it in `JAGS`](#)
3. [Plotting out the results](#)