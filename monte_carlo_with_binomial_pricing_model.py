import random
import numpy as np
#import Tkinter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
import random
import sys

import scipy.stats as si
#from scipy.stats import norm
import statistics
 

def black_scholes(S, K, T, r, sigma):

    d1 = (np.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    d2 = (np.log(S / K) + (r - 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    
    call = (S * si.norm.cdf(d1, 0.0, 1.0) - K * np.exp(-r * T) * si.norm.cdf(d2, 0.0, 1.0))
    
    return call

def simulate_paths_continuous(S_zero,mu,r,sigma,T,d_t):
    time_steps = int(T/d_t)
    S = S_zero
    stock_path = [S_zero] #tracks S at each step
    for i in range(time_steps):
        rv = np.random.normal(loc=(r-(0.5*math.pow(sigma,2)))*d_t, scale=sigma*math.sqrt(d_t))
        S = S*math.exp(rv)
        stock_path.append(S)

    return stock_path

def simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p):

    time_steps = int(T/d_t)
    
    lambda_up = math.exp(mu*d_t+sigma*math.sqrt((1-p)/p)*math.sqrt(d_t))
    lambda_down = math.exp(mu*d_t-sigma*math.sqrt(p/(1-p))*math.sqrt(d_t))
    p_tilde =  (math.exp(r*d_t)-lambda_down) / (lambda_up - lambda_down)
    S = S_zero
    stock_path = [S_zero] #tracks S at each step

    for i in range(time_steps):
        random = np.random.uniform(0, 1)
        if(random < p_tilde):
            S = S * lambda_up
            stock_path.append(S)
        else:
            S = S * lambda_down
            stock_path.append(S)
    return stock_path

'''
change in the log of the stock pice - has a normal
run with 100 run with another 100 see the deviatiions if a big number run it with more paths
'''
def monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths, B = 0.95):
    
    #calc risk nuteral world dynamics:
    if path_type == 1:
        p = 0.5
    elif path_type == 2:
        p = (2.0/3.0)
    elif path_type == 3:
        p = 0.5
    else:
        #print('invalid path type', file=sys.stderr)
        print('invalid path type: ',path_type)
        return


    cash_flows = []
    avg_discounted_cash_flow = 0.0
    A_cashflow = 0.0
    G_cashflow = 0.0
    H_cashflow = 0.0

    for i in range(num_paths):
        if path_type == 1 or path_type == 2:
            path = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
        elif path_type == 3:
            path = simulate_paths_continuous(S_zero,mu,r,sigma,T,d_t)

        cashflow = 0.0

        if(option_type == "C"):
            #calulate cash flows for european call option
            if(path[-1] > K):
                #discount it by the number of time steps (T/d_t) * d_t
                cashflow = math.exp(-1*r*T)*(path[-1] - K)
        elif(option_type == "B"):
            #calulate cash flows for barrier 
            previous = path[0]
            for j, price in enumerate(path):
                if j == 0:
                    continue
                if (previous > B  and price <= B) or (previous < B and price >= B):
                    if(K > B):
                        cashflow += math.exp(-1*r*(j*d_t)) * (K-B)
                    else:
                        cashflow = 0.0
                    break
                previous = price
                
                
        elif(option_type == "A"):
            #calulate cash flows for average strike Asian call option

            #arithmetic mean
            K = sum(path)/len(path)
            if(path[-1] > K):
                cashflow = math.exp(-1*r*T)*(path[-1] - K)
            A_cashflow += (1/num_paths) * cashflow

            #geometric mean
            K = si.gmean(path)
            if(path[-1] > K):
                cashflow = math.exp(-1*r*d_t)*(path[-1] - K)
            G_cashflow += (1/num_paths) * cashflow

            #harmonic mean
            K = statistics.harmonic_mean(path)
            if(path[-1] > K):
                cashflow = math.exp(-1*r*d_t)*(path[-1] - K)
            H_cashflow += (1/num_paths) * cashflow

        elif(option_type == "M"):
            #calulate cash flows for Minimum strike call option
            K = min(path)
            if(path[-1] > K):
                cashflow = math.exp(-1*r*d_t)*(path[-1] - K)
        else:
            #print('invalid option type', file=sys.stderr)
            print('invalid option type: ',option_type)
            return
        avg_discounted_cash_flow += (1/num_paths) * cashflow
        
    #avg_discounted_cash_flow = sum(cash_flows)/len(cash_flows)
    print_string = ''
    if(path_type != 3):
        print_string += ("binomial:  "+str(path_type)+"  ")
    else:
        print_string += ("continuous:   ")
    
    if(d_t == 0.1):
        print_string += "d_t: 0.1    "
    elif(d_t == 0.01):
        print_string += "d_t: 0.01   "
    elif(d_t == 0.001):
        print_string += "d_t: 0.001  "


    if(option_type == "C"):
        print_string += "C("+str(S_zero)+", "+str(K)+", "+str(T)+") = " + str(avg_discounted_cash_flow)
   
    elif(option_type == 'B'):
        print_string += "B("+str(S_zero)+", "+str(K)+", "+str(B)+", "+str(T)+") = " +str(avg_discounted_cash_flow)
    
    elif(option_type == 'A'):
        print(print_string + "mean_type: A " + "A("+str(S_zero)+", "+str(T)+") = " +str(A_cashflow))
        print(print_string + "mean_type: G " + "A("+str(S_zero)+", "+str(T)+") = " +str(G_cashflow))
        print(print_string + "mean_type: H " + "A("+str(S_zero)+", "+str(T)+") = " +str(H_cashflow))
        return (A_cashflow,G_cashflow,H_cashflow)

    elif(option_type == 'M'):
        print_string += "M("+str(S_zero)+", "+str(T)+") = " +str(avg_discounted_cash_flow)
    
    print(print_string)
    return avg_discounted_cash_flow

'''
Questions:
1. is our funciton for part 1 just supposed to simulate 1 path?
2. is supposed to calculate stock price path.
3. does D_S forumla calulate s at time d_t or does it 
4. for the asian call: should we set cash flow = 0 if S[i] < K? same for euro call and min strike
5. how do i calulate the cashflow for barrier option?

'''
if __name__ == "__main__":
    #def simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p):
    #simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p):
    '''
    =================binomial_1=================
    '''
    S_zero = 1
    mu = 0.07
    r = 0.03
    sigma = 0.2
    T = 2
    d_t = 0.1
    p = 0.5
    binomial_1_1 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_1_1)
    #print()
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_1_1)
    plt.legend()
    plt.title("binomial_1 dt = 0.1")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    plt.savefig("binomial_1_1.pdf")
    plt.show()

    #simulate_path(1,0.07,0.03,0.2,2,0.01,0.5)
    d_t = 0.01
    binomial_1_2 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_1_2)
    #print()
    plt.cla()
    plt.title("binomial_1 dt = 0.01")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_1_2)
    plt.savefig("binomial_1_2.pdf")
    plt.show()

    #simulate_path(1,0.07,0.03,0.2,2,0.001,0.5)
    d_t = 0.001
    binomial_1_3 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_1_3)
    #print()
    plt.cla()
    plt.title("binomial_1 dt = 0.001")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_1_3)
    plt.savefig("binomial_1_3.pdf")
    plt.show()
    '''
    =================binomial 2=================
    '''
    #simulate_path(1,0.07,0.03,0.2,2,0.1,(2/3))
    d_t = 0.1
    p = (2.0/3.0)
    binomial_2_1 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_2_1)
    #print()
    plt.cla()
    plt.title("binomial_2 dt = 0.1")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_2_1)
    plt.savefig("binomial_2_1.pdf")
    plt.show()

    #simulate_path(1,0.07,0.03,0.2,2,0.01,(2/3))
    d_t = 0.01
    binomial_2_2 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_2_2)
    #print()
    plt.cla()
    plt.title("binomial_2 dt = 0.01")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_2_2)
    plt.savefig("binomial_2_2.pdf")
    plt.show()

    #simulate_path(1,0.07,0.03,0.2,2,0.001,(2/3))
    d_t = 0.001
    binomial_2_3 = simulate_stock_path(S_zero, mu, r, sigma, T , d_t, p)
    #print(binomial_2_3)
    #print()
    plt.cla()
    plt.title("binomial_2 dt = 0.001")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), binomial_2_3)
    plt.savefig("binomial_2_3.pdf")
    plt.show()
    '''
    =================continuous=================
    '''
    #simulate_paths_continuous(S_zero,r,sigma,T,d_t):
    #simulate_paths_continuous(1,0.03,0.2,,2,0.1) 
    S_zero = 1
    mu = 0.07
    r = 0.03
    sigma = 0.2
    T = 2 
    d_t = 0.1
    p = 0.5
    continuous_1 = simulate_paths_continuous(S_zero,mu,r,sigma,T,d_t)
    #print(continuous_1)
    #print()
    plt.cla()
    plt.title("continuous dt = 0.1")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), continuous_1)
    plt.savefig("continuous_1.pdf")
    plt.show()
    
    #simulate_paths_continuous(1,0.03,0.2,,2,0.01)    
    d_t = 0.01
    continuous_2 = simulate_paths_continuous(S_zero,mu,r,sigma,T,d_t)
    #print(continuous_2)
    #print()
    plt.cla()
    plt.title("continuous dt = 0.01")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), continuous_2)
    plt.savefig("continuous_2.pdf")
    plt.show()
    
    #simulate_paths_continuous(1,0.03,0.2,,2,0.001) 
    d_t = 0.001
    continuous_3 = simulate_paths_continuous(S_zero,mu,r,sigma,T,d_t)
    #print(continuous_3)
    #print()
    plt.cla()
    plt.title("continuous dt = 0.001")
    plt.ylabel('stock Price')
    plt.xlabel('Time ')
    g, = plt.plot( np.arange(0,int(T/d_t)+1), continuous_3)
    plt.savefig("continuous_3.pdf")
    plt.show()

    '''
    ===================================================
                     Monte Carlo
    ===================================================
    '''

    '''
    =================European Call=================
    '''

    S_zero = 1
    mu = 0.07
    r = 0.03
    sigma = 0.2
    T = 2
    d_t = 0.1
    p = 0.5
    #def monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths, mean_option = 'A', B = 0.95)
    print()
    print("European Call Option:")
    
    option_type = 'C'
    K = 1
    path_type = 1
    num_paths = 5000

    bs = black_scholes(S_zero,K,T,r,sigma)
    print("black_scholes:      ",bs)

    #european_call_binomial_1_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)

    #european_call_binomial_1_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)

    #european_call_binomial_1_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 2
    #european_call_binomial_2_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #european_call_binomial_2_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #european_call_binomial_2_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 3
    #european_call_continuous_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #european_call_continuous_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #european_call_continuous_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    print()
    print("Barrier Call Option:")
    option_type = 'B'
    path_type = 1
    #barrier_binomial_1_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_binomial_1_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_binomial_1_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 2
    #barrier_binomial_2_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_binomial_2_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_binomial_2_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 3
    #barrier_continuous_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_continuous_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_continuous_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()
    
    print()
    print("Average Strike Asian Call Option:")
    option_type = 'A'
    path_type = 1
    #asian_binomial_1_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #asian_binomial_1_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #asian_binomial_1_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 2
    #asian_binomial_2_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #asian_binomial_2_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #asian_binomial_2_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 3
    #barrier_continuous_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_continuous_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #barrier_continuous_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    print()
    print("Minimum Strike Call Option:")
    option_type = 'M'
    path_type = 1
    #min_strike_binomial_1_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_binomial_1_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_binomial_1_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 2
    #min_strike_binomial_2_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_binomial_2_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_binomial_2_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()

    path_type = 3
    #min_strike_continuous_1
    d_t = 0.1
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_continuous_2
    d_t = 0.01
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    #min_strike_continuous_3
    d_t = 0.001
    monte_carlo(S_zero, mu, r, sigma, T, d_t, option_type, K, path_type, num_paths)
    print()
    