# Ascertainment

Real data, especially low coverage data, often discard singletons.
This is easy to deal with in principle: given a branch length weighting method,
we need only modify the method so that any branch subtending only a single sample gets weight zero.
For instance pairwise divergence between $i$ and $j$ out of a sample of size $N$ --
with $x_1=1$ if a branch is ancestral to $i$;
with $x_2=1$ if a branch is ancestral to $j$;
and $x_3$ the number of the $N$ samples below the branch --
a weighting function is
$$\begin{aligned}
    f(x_1, x_2, x_3) &=
        \left( x_1 (1-x_2) + (1-x_1) x_2 \right) \mathbb{1}_{x_3 > 1} \mathbb{1}_{x_3 < N} .
\end{aligned}$$
Mean divergence for randomly chosen pairs from the group is given by,
with now $x$ the number of individuals subtending the branch,
$$\begin{aligned}
    f(x) &=
        \frac{2}{N(N-1)} \left( x (1-x) \right) \mathbb{1}_{x > 1} \mathbb{1}_{x < N} .
\end{aligned}$$

Now suppose you want to compute the *expected* value of ascertained divergence
for a group of size $N=5$ out of a larger sample of size $M$.
This is the chance that we pick (a) two samples, that differ, and then
(b) three more, that are also not all the same,
so
$$\begin{aligned}
    f(x) &=
        \frac{2}{M(M-1)} \left( x (1-x) \right) \frac{3}{(M-2)(M-3)(M-4)} \left( (x-1) (x-2) (M-x-1) + (x-1) (M-x-1) (M-x-2) \right) .
\end{aligned}$$

Working out these special cases could get tiring, though.


# A general method

Given a particular set of samples and a pattern,
it's easy to write down a formula for what we want to compute.
For instance, "one and three both have allele 1 but not two" is $x_1 (1-x_2) x_3$.
But then, we want to *average* such things across choices of samples and alleles:
like maybe samples one and two are drawn from group $A$, while three is from group $B$.
Here's a general method to do this.
First start with a polynomial function of $p$ variables that defines the pattern we want:
$$\begin{aligned}
    g(y_1, \ldots, y_p) = \sum_{I \subset \{1, \ldots, p\}} a_I y^I .
\end{aligned}$$
Then pick a partition of $[p]$ that says which variable corresponds to a draw from the same group of samples.
which we encode by saying that $j(i)$ is the index of the partition that variable $i$ is in, for $1 \le i \le p$.
Define some random variables $Z_i$ for $1 \le i \le p$
so that $Z^{(j)} = \{ Z_i \; : \; j(i) = j\}$ are uniform draws without replacement from the $j$th group of samples.
Suppose there are $x_j$ type 1 alleles in the $n_j$ samples of the $j$th group.
We want to compute
$$\begin{aligned}
   f(x_1, \ldots, x_p) &= \E[g(Z_1, \ldots, Z_p)] .
\end{aligned}$$

By linearity, each term can be evaluated separately.
Also, $Z_i$ and $Z_k$ are independent if $j(i) \neq j(k)$.
On the other hand, if $j(1) = \cdots = j(k) = j$ then
$$\begin{aligned}
    \E[Z_1 \cdots Z_k]
        &=
        \frac{x_j (x_j - 1) \cdots (x_j - k + 1)}{n_j (n_j - 1) \cdots (n_j - k + 1)} \\
        &=
        \frac{(x_j)_k}{(n_j)_k},
\end{aligned}$$
where $(a)_k = a (a-1) \cdots (a-k+1)$ is the falling factorial.

This says that to evaluate the function $f$, we can do the following.
First, find the coefficients of the polynomial
$$\begin{aligned}
    h(u_1, \ldots, u_j) 
        &= 
        g(u_{j(1)}, \ldots, u_{j(p)}) \\
        &=
        \sum_m a_m u^m,
\end{aligned}$$
where the sum is over multiindices $m = (m_1, \ldots, m_p)$.
Then, 
$$\begin{aligned}
    f(x_1, \ldots, x_p)
    &=
    \sum_m a_m \prod_{j=1}^p \frac{ (x_j)_{a_j} }{ (n_j)_{a_j} } .
\end{aligned}$$


