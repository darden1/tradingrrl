# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tradingrrl import TradingRRL

def main():

    fname = "../data/USDJPY30.csv"
    init_t = 6000

    T = 1000
    M = 200
    mu = 10000
    sigma = 0.04
    rho = 1.0
    n_epoch = 10000

    # RRL agent with initial weight.
    ini_rrl = TradingRRL(T, M, init_t, mu, sigma, rho, n_epoch)
    ini_rrl.load_csv(fname)
    ini_rrl.set_t_p_r()
    ini_rrl.calc_dSdw()

    rrl = TradingRRL(T, M, init_t,  mu, sigma, rho, n_epoch)
    rrl.all_t = ini_rrl.all_t
    rrl.all_p = ini_rrl.all_p
    rrl.set_t_p_r()
    
    #----GA fit
    import random

    from deap import algorithms
    from deap import base
    from deap import creator
    from deap import tools

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", np.ndarray, fitness=creator.FitnessMax)

    toolbox = base.Toolbox()

    # toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("attr_bool", random.uniform, -3, 3)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=M+2)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    def evalOneMax(individual):
        return rrl.calc_S(individual),
    #     return sum(individual),

    def cxTwoPointCopy(ind1, ind2):
        size = len(ind1)
        cxpoint1 = random.randint(1, size)
        cxpoint2 = random.randint(1, size - 1)
        if cxpoint2 >= cxpoint1:
            cxpoint2 += 1
        else: # Swap the two cx points
            cxpoint1, cxpoint2 = cxpoint2, cxpoint1

        ind1[cxpoint1:cxpoint2], ind2[cxpoint1:cxpoint2] \
            = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

        return ind1, ind2


    toolbox.register("evaluate", evalOneMax)
    toolbox.register("mate", cxTwoPointCopy)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)


    random.seed(64)

    pop = toolbox.population(n=100)
    hof = tools.HallOfFame(1, similar=np.array_equal)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    population, logbook = algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=100, stats=stats, halloffame=hof)

    # Plot results.
    # Training for initial term T.
    df =pd.DataFrame(logbook)
    plt.plot((df["gen"])*100,df["max"])
    plt.title("Sharp's ratio optimization")
    plt.xlabel("Epoch times")
    plt.ylabel("Sharp's ratio")
    plt.grid(True)
    plt.savefig("sharp's ratio optimization.png", dpi=300)
    plt.close

    fig, ax = plt.subplots(nrows=3, figsize=(15, 10))
    t = np.linspace(1, rrl.T, rrl.T)[::-1]
    ax[0].plot(t, rrl.p[:rrl.T])
    ax[0].set_xlabel("time")
    ax[0].set_ylabel("USDJPY")
    ax[0].grid(True)

    ax[1].plot(t, ini_rrl.F[:rrl.T], color="blue", label="With initial weights")
    ax[1].plot(t, rrl.F[:rrl.T], color="red", label="With optimized weights")
    ax[1].set_xlabel("time")
    ax[1].set_ylabel("F")
    ax[1].legend(loc="upper left")
    ax[1].grid(True)

    ax[2].plot(t, ini_rrl.sumR, color="blue", label="With initial weights")
    ax[2].plot(t, rrl.sumR, color="red", label="With optimized weights")
    ax[2].set_xlabel("time")
    ax[2].set_ylabel("Sum of reward[yen]")
    ax[2].legend(loc="upper left")
    ax[2].grid(True)
    plt.savefig("rrl_train.png", dpi=300)
    fig.clear()


    # Prediction for next term T with optimized weight.
    # RRL agent with initial weight.
    ini_rrl_f = TradingRRL(T, M, init_t-T, mu, sigma, rho, n_epoch)
    ini_rrl_f.all_t = ini_rrl.all_t
    ini_rrl_f.all_p = ini_rrl.all_p
    ini_rrl_f.set_t_p_r()
    ini_rrl_f.calc_dSdw()

    # RRL agent with optimized weight.
    rrl_f = TradingRRL(T, M, init_t-T, mu, sigma, rho, n_epoch)
    rrl_f.all_t = ini_rrl.all_t
    rrl_f.all_p = ini_rrl.all_p
    rrl_f.set_t_p_r()
    rrl_f.w = rrl.w
    rrl_f.calc_dSdw()

    fig, ax = plt.subplots(nrows=3, figsize=(15, 10))
    t_f = np.linspace(rrl.T+1, rrl.T+rrl.T, rrl.T)[::-1]
    ax[0].plot(t_f, rrl_f.p[:rrl_f.T])
    ax[0].set_xlabel("time")
    ax[0].set_ylabel("USDJPY")
    ax[0].grid(True)

    ax[1].plot(t_f, ini_rrl_f.F[:rrl_f.T], color="blue", label="With initial weights")
    ax[1].plot(t_f, rrl_f.F[:rrl_f.T], color="red", label="With optimized weights")
    ax[1].set_xlabel("time")
    ax[1].set_ylabel("F")
    ax[1].legend(loc="lower right")
    ax[1].grid(True)

    ax[2].plot(t_f, ini_rrl_f.sumR, color="blue", label="With initial weights")
    ax[2].plot(t_f, rrl_f.sumR, color="red", label="With optimized weights")
    ax[2].set_xlabel("time")
    ax[2].set_ylabel("Sum of reward[yen]")
    ax[2].legend(loc="lower right")
    ax[2].grid(True)
    plt.savefig("rrl_prediction.png", dpi=300)
    fig.clear()

if __name__ == "__main__":
    main()