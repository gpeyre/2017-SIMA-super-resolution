This toolbox solves the BLASSO problem for 1-D and 2-D super-resolution problem using an enhanced version of Frank-Wolfe.

The algorithm alternates between:

- a Frank-Wolfe step to add a new dirac location (which is the matching pursuit step of correlation maximization)
- a non-convex update step on the (amplitude,positions) (with an initialization using a classical LASSO).

It is guaranteed to converge even without the non-convex step, but this non-convex step typically leads to termination in a finite number of steps (which in practice matches the number of spikes in the measure to recover when noise and regularization are low).

References
-------

This algorithm was first described in:

> Kristian Bredies and Hanna Pikkarainen. Inverse problems in spaces of measures. ESAIM: Control, Optimisation and Calculus of Variations, 19(1):190-218, 2013.  

and used extensively in:

>  N. Boyd, G. Schiebinger, B. Recht [The Alternating Descent Conditional Gradient Method for Sparse Inverse Problems](http://arxiv.org/abs/1507.01562), 2015.

with a [Julia implementation](https://www.stat.berkeley.edu/~nickboyd/adcg/).

For an modern view on Frank-Wolfe (including 1/k rate on the objective), see:

> Revisiting Frank-Wolfe: Projection-Free Sparse Convex Optimization
ICML 2013: International Conference on Machine Learning (2013)

Copyright
-------

Copyright (c) 2017 Clarice Poon and Gabriel Peyré
