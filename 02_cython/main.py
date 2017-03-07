# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from tradingrrl import ChartData, TradingRRL, plot_hist, save_weights, load_weights

def main():

    fname = "../data/USDJPY30.csv"
    init_t = 6000

    T = 1000
    M = 200
    q_threshold = 0.7
    mu = 10000
    sigma = 0.04
    alpha = 1.0
    n_epoch = 10000

    # Prepare chart data
    cd = ChartData()
    cd.load_csv(fname)
    all_t, all_p = cd.all_t, cd.all_p

    # Training RRL agent.
    ini_rrl = TradingRRL(T, M, init_t, q_threshold, mu, sigma, alpha, n_epoch)
    ini_rrl.set_t_p_r(all_t, all_p)
    ini_rrl.calc_dSdw()
    rrl = TradingRRL(T, M, init_t, q_threshold, mu, sigma, alpha, n_epoch)
    rrl.set_t_p_r(all_t, all_p)
    rrl.fit()

    # Plot results.
    # Training for initial term T.
    plt.plot(range(len(rrl.epoch_S)),rrl.epoch_S)
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
    ini_rrl_f = TradingRRL(T, M, init_t-T, q_threshold, mu, sigma, alpha, n_epoch)
    ini_rrl_f.set_t_p_r(all_t, all_p)
    ini_rrl_f.calc_dSdw()
    rrl_f = TradingRRL(T, M, init_t-T, q_threshold, mu, sigma, alpha, n_epoch)
    rrl_f.set_t_p_r(all_t, all_p)
    rrl_f.w = rrl.w.copy()
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