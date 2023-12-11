import numpy as np
import matplotlib.pyplot as plt

def plot_alpha(alpha, gene=None):
    alpha_sumk = []
    for g in range(len(alpha)):
        alpha_partk = 0.
        for i in range(len(alpha[g])):
            alpha_partk += np.array(alpha[g][i])
        alpha_sumk.append(alpha_partk)
    #alpha_sumk = np.array(alpha_sumk)
    #alpha_sumk = np.sum(alpha, axis=1).tolist()
    if gene is not None:
        _ = plt.plot(alpha_sumk[gene], '.')
        plt.xlabel('index')
        plt.ylabel('alpha')
    
def cal_pip(alpha, ng):
    pip = []
    for g in range(ng):
        ans = 1
        for k in range(len(alpha[g])):
            ans *= (1-np.array(alpha[g][k]))
        pip.append(1-ans)
    #pip = np.array(pip)
    return pip

def plot_pip(alpha, ng, gene=None):
    pip = cal_pip(alpha, ng)
    _ = plt.plot(pip[gene], '.')
    plt.xlabel('index')
    plt.ylabel('PIP')

def cal_cs(alpha, ng, K, rho):
    cs = []
    for g in range(ng):
        csg = []
        for k in range(K[g]):
            alphagk = alpha[g][k]
            alphagk_argsort = np.argsort(alphagk)[::-1]
            alphagk = np.sort(alphagk)[::-1]
            ans = 0
            i = 0
            while ans < rho:
                ans += alphagk[i]
                i += 1
            csg.append(alphagk_argsort[:i])
        cs.append(csg)
    return cs

    
def plot_cs(alpha, ng, K, rho, gene):
    cs = cal_cs(alpha, ng, K, rho)
    ncol = min(K[gene], 3)
    nrow = K[gene] // ncol + 1
    plt.figure(figsize=(15,5))
    for k in range(K[gene]):
        plt.subplot(nrow, ncol, k+1)
        plt.plot(alpha[gene][k], '.')
        plt.plot(cs[gene][k], np.array(alpha[gene][k])[cs[gene][k]], 'o', mfc='none')
        plt.xlabel('index')
        plt.ylabel('alpha')
    return cs

def model_plot(alpha, beta=None, beta_hat=None, rho=None, gene=None, sigma=None, sigma_hat=None):
    ng = len(alpha)
    K = [len(alpha[g]) for g in range(ng)]
    if beta is None:
        plt.figure(figsize=(15,5))
        plt.subplot(1, 2, 1)
        plot_alpha(alpha, gene)
        plt.subplot(1, 2, 2)
        plot_pip(alpha, ng, gene)
        plt.show()
    else:
        plt.figure(figsize=(15,5))
        plt.subplot(1, 3, 1)
        plot_alpha(alpha, gene)
        plt.subplot(1, 3, 2)
        plot_pip(alpha, ng, gene)
        plt.subplot(1, 3, 3)
        plt.plot(beta, beta_hat, '.')
        x = np.min([beta, beta_hat])
        y = np.max([beta, beta_hat])
        plt.plot(np.arange(x, y, 0.1), np.arange(x, y, 0.1), alpha=0.5)
        plt.ylabel('beta_hat')
        plt.xlabel('beta')
        #plt.subplot(1, 4, 4)
        #plt.plot(sigma, sigma_hat, '.')
        #x = np.min([sigma, sigma_hat])
        #y = np.max([sigma, sigma_hat])
        #plt.plot(np.arange(x, y, 0.1), np.arange(x, y, 0.1), alpha=0.5)
        #plt.ylabel('sigma_hat')
        #plt.xlabel('sigma')
        
    if gene is not None:
        cs = plot_cs(alpha, ng, K, rho, gene)
        
def simu_plot(alpha, beta=None, beta_hat=None, rho=None, gene=None, sigma=None, sigma_hat=None):
    # lg are equal
    ng = len(alpha)
    K = [len(alpha[g]) for g in range(ng)]
    alpha_sumk = np.sum(alpha, axis=1)
    pip = cal_pip(alpha, ng)
    plt.figure(figsize=(15,5))
    plt.subplot(1, 3, 1)
    _=plt.boxplot(alpha_sumk)
    plt.xlabel('index')
    plt.ylabel('alpha')
    
    plt.subplot(1, 3, 2)
    _=plt.boxplot(pip)
    plt.xlabel('index')
    plt.ylabel('pip')
    
    plt.subplot(1, 3, 3)
    plt.plot(beta, beta_hat, '.')
    x = np.min([beta, beta_hat])
    y = np.max([beta, beta_hat])
    plt.plot(np.arange(x, y, 0.1), np.arange(x, y, 0.1), alpha=0.5)
    plt.ylabel('beta_hat')
    plt.xlabel('beta')
    
    #plt.subplot(1, 4, 4)
    #plt.plot(sigma, sigma_hat, '.')
    #x = np.min([sigma, sigma_hat])
    #y = np.max([sigma, sigma_hat])
    #plt.plot(np.arange(x, y, 0.1), np.arange(x, y, 0.1), alpha=0.5)
    #plt.ylabel('sigma_hat')
    #plt.xlabel('sigma')
    
    if gene is not None:
        cs = plot_cs(alpha, ng, K, rho, gene)