# Monte Carlo and Binomial Simulations for Option Pricing
Program that compares 3 implementations of the Binomial option pricing model, and then uses a Monte carlo simulation to price European, Asian, Barrier, and Minimum Strike call options.

Program that compares and constrasts Stock pice path models

Program to price European, Asian, Barrier, and minimum strike call options, Using both Binomial pricing and monte Carlo pricing options

## Program Overview 
This program that takes as input µ, r, σ, S0, T, ∆t and simulates the stock price from
time 0 to T in time steps of ∆t for the risk neutral dynamics given by:

∆S = µS∆t + σS∆W    and    ∆S˜ = rS˜∆t + σS˜∆W .

using each of the following Models:

(a) Binomial mode I: compute λ± from µ, σ assuming that p = 1/2, and then computing ~p (probaility in the risk nuteral world).

(b) Binomial mode II: compute λ± from µ, σ assuming that p = 2/3, and then computing ~p (probaility in the risk nuteral world).

(c) Continuous mode: using the continuous risk neutral dynamics r, σ generate at time step ∆t as if the discrete model were taken to the limit dt → 0.

Then for each of the Stock pice path models above, and each of the three time discretizations, ∆t = 0.1, 0.01, 0.0001 the program computes the price of each of the  following options using Monte carlo with 5000 samples (in all cases S0 is the initial stock price):

1. C(S0, K, T): the European Call option with strike K and exercise time T.

2. B(S0, K, B, T): the Barrier option with payoff K, barrier B and horizon T. If and
when the price hits the barrier B the holder may buy at B and sell at K.

3. A(S0, T): the average strike Asian call option with expiry T – a European call option
with strike at the average value of the stock price over [0, T].

4. M(S0, T): the minimum strike call option – a European call option with strike at the
minimum price over [0, T].

## Results

In an effort to compare and contrast the each of the three models, here are the plots of representative stock price paths for S0 = 1, µ = 0.07, r = 0.03, σ = 0.2, T = 2 using ∆t = 0.1, 0.01, 0.0001.

Binomial mode I:

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

Binomial mode II:

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

Continuous mode:
![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

For each of the Stock pice path models above, and each of the three time discretizations, ∆t = 0.1, 0.01, 0.0001 the program computes the price for each of the following options using Monte carlo simulation with 5000 samples (S0 = 1, µ = 0.07, r = 0.03, σ = 0.2, T = 2):

(a) C(1, 1, 2) and compared with the analytic formula.\
![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

(b) B(1, 1, 0.95, 2)
![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

(c) A(1, 2), for three possible definitions of “average”: the harmonic, arithmetic,
and geometric means. Explain the relative ordering of these prices.
![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

(d) M(1, 2).
![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

## Usage

```
Python3 monte_carlo_with_binomial_pricing_model.py
```


