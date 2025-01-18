## Compensated summation

A method known as compensated summation [1], [4] is implemented for all integrators EXCEPT WHFast:

When we advance our system by $\text{d}t$, we have 

$x_{n+1} = x_n + \delta x$

Since $\delta x$ is very small compared to $x_n$, many digits of precision will be lost.
By compensated summation, we keep track of the losing digits using another variable, which
allows us to effectively eliminates round off error with very little cost.

It is not implemented for WHFast yet due to implementation
difficulty.
