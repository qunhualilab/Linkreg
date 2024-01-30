import multiprocessing
import numpy as np
from scipy.stats import norm
import copy
from qpsolvers import solve_qp
from datetime import datetime
import logging
#logger = logging.getLogger(__name__)
#logging.basicConfig(level=logging.INFO, datefmt='%Y/%m/%d %H:%M:%S', format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler()])
#ctx = multiprocessing.get_context('spawn')

#import logging
import argparse
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from data_process import process
from gene_class import gene_class

def estimate(D, ymean, alpha_sumk):
    return [np.dot(Dg, alpha_sumk[g]) + ymean[g] for g, Dg in enumerate(D)]

def cal_prior(genes, trsd=100000, maxK=30, medianK=10):
    def helper(cur, ref, maxK):
        adjusted = np.array([1-(1-1/cur)**K for K in range(1, maxK+1)])
        return np.argmin(np.abs(adjusted-ref))+1
    priors, Kg = [], np.zeros(len(genes))
    firsts = []
    for gene in genes:
        distances = np.abs(np.mean(gene.cCRE_loc, axis=1)-gene.tss)
        prior = np.exp((-distances/trsd).astype(float))
        firsts.append(np.sum(prior))
        prior /= firsts[-1]
        priors.append(prior.tolist())
    Kg[np.argsort(firsts)[len(firsts)//2]] = medianK
    Kg = [helper(firsts[i], 1-(1-1/firsts[np.argsort(firsts)[len(firsts)//2]])**medianK, maxK) for i in range(len(genes))]
    pi = [[priors[i] for _ in range(Kg[i])] for i in range(len(Kg))]
    return pi, Kg

def cal_variance(X, y, K, alpha, sigma, converge, beta_hat):
    ng = len(X)
    nc = [len(X[g]) for g in range(ng)]
    l = [len(X[g][0]) for g in range(ng)]
    converge_idx = np.arange(ng)[converge==1]
    Sigma = []
    for g in converge_idx:
        Sigmag = []
        for k in range(K[g]):
            Sigmagk = -np.outer(alpha[g][k], alpha[g][k])
            for i in range(l[g]):
                Sigmagk[i][i] = alpha[g][k][i]*(1-alpha[g][k][i])
            Sigmag.append(Sigmagk)
        Sigma.append(Sigmag)
    alpha_sumk = [np.sum(alpha[g], axis=0).tolist() for g in range(ng)]
    G = np.array([[np.dot(alpha_sumk[g], X[g][i]) for i in range(nc[g])] for g in converge_idx])
    first = 0.
    for j, g in enumerate(converge_idx):
        for i in range(nc[g]):
            soft = 0. 
            for k in range(K[g]):
                soft += np.dot(np.dot(X[g][i].transpose(), Sigma[j][k]), X[g][i])
            first += (np.outer(G[j][i], G[j][i]) + soft) / sigma[g]**2
    Xg = [np.transpose([np.dot(alpha_sumk[g], X[g][i])/sigma[g]**2 for i in range(nc[g])]) for g in converge_idx]
    var_y = []
    for j, g in enumerate(converge_idx):
        var_yg = np.zeros((nc[g], nc[g]))
        for m in range(nc[g]):
            for n in range(nc[g]):
                var_yg[m, n] = np.dot(np.dot(np.transpose(np.dot(X[g][m], beta_hat)), np.sum(np.array(Sigma[j]), axis=0)), np.dot(X[g][n], beta_hat))
                if m == n: var_yg[m, n] += sigma[g]**2
        var_y.append(var_yg)
    var_y = np.array(var_y)
    middle = 0
    for g in range(len(converge_idx)):
        middle += np.dot(np.dot(Xg[g], var_y[g]), np.transpose(Xg[g]))
    return np.dot(np.dot(np.linalg.inv(first.astype(float)), middle), np.transpose(np.linalg.inv(first.astype(float))))

def mtnm(yg, Dg, sigmag, alphag, pig, nc, Kg, lg, scaleg, scale_change, g, converge):
    alpha_sumk = np.sum(alphag, axis=0).tolist()
    n_iter = 0
    alphag_old = .001
    while True:
        ans = [np.abs(np.array(alphag_old)-np.array(alphag))/np.maximum(1e-8, np.array(alphag_old))]
        if np.max(ans) < 1e-3:
            break
        alphag_old = np.copy(alphag)
        n_iter += 1
        rg = yg - np.dot(scaleg*Dg, alpha_sumk)
        for k in range(Kg):
            rgk = (rg + np.dot(scaleg*Dg, alphag[k])).astype(float)
            lbf = np.array([np.sum(norm.logpdf(rgk, loc=scaleg*Dg[:, j], scale=sigmag) - 
                                   norm.logpdf(rgk, loc=0, scale=sigmag)) for j in range(lg)])
            lbf = lbf - np.max(lbf)
            alphag[k] = [np.exp(lbf[j])*pig[k][j]/(np.sum(np.exp(lbf)*np.array(pig[k]))) for j in range(lg)]
            rg = rgk - np.dot(scaleg*Dg, alphag[k])
        alpha_sumk = np.sum(alphag, axis=0).tolist()
        
        x = np.dot(Dg, alpha_sumk)
        eres_var = np.sum((yg-x*scaleg)**2)
        var = [(np.dot(Dg**2, alphag[k])-np.dot(Dg, alphag[k])**2).tolist() for k in range(len(alphag))]
        var = [[np.max([_var, 0.]) for _var in _vars] for _vars in var]
        eres_var = eres_var + np.sum(var)
        sigmag = np.sqrt(eres_var/nc)
        if scale_change == True:
            numerator, denominator = np.dot(yg, x), np.dot(x, x)
            if denominator == 0:
                print('%d gene, np.dot(x, x)=0' % g)
                return alpha_sumk, alphag, sigmag, scaleg, -1
            elif numerator < 0:
                print('%d gene, np.dot(yg, x)<0' % g)
                return alpha_sumk, alphag, sigmag, scaleg, -1
            scaleg = 1/denominator*numerator
        if n_iter >= 300:
            print('%d gene, SUSIE not convergent: %f' % (g, np.max(ans)))
            return alpha_sumk, alphag, sigmag, scaleg, -1
    return alpha_sumk, alphag, sigmag, scaleg, 1

def model(X, y, K=None, n_iter=200, parallel=True, genes=None, allow_scale_change=False, scale_change=False, norm_beta=False, pi=None, ridge=0):
    if parallel == True:
        ctx = multiprocessing.get_context('spawn')
    if scale_change == True and allow_scale_change == False:
        allow_scale_change = True
        print('scale_change is set to True, allow_scale_change has been changed to True.')
    if scale_change == True and allow_scale_change == True and norm_beta == False:
        norm_beta = True
        print('If both scale_change and allow_scale_change are set to be True, norm_beta should be set to be True')
    test = 0.
    for g in range(len(X)):
        for i in range(len(X[g])):
            test += np.outer(np.dot(np.ones(len(X[g][i])), X[g][i]), np.dot(np.ones(len(X[g][i])), X[g][i]))
    if np.linalg.det(test.astype(float)) == 0:
        print('Singular matrix')

    ng = len(X)
    nc = [len(X[g]) for g in range(ng)]
    l = [len(X[g][0]) for g in range(ng)]
    if (genes is not None) and (K is None):
        pi, K = cal_prior(genes)
    elif genes is not None:
        pi, _ = cal_prior(genes)
    elif pi is None:
        pi = [[[1/l[g] for _ in range(l[g])] for _ in range(K[g])] for g in range(ng)]
    if isinstance(K, list):
        pass
    else:
        K = [K for _ in range(ng)]
    for g in range(ng):
        pi[g] = [pi[g][0] for _ in range(K[g])]

    alpha = copy.deepcopy(pi)
    #alpha = [[[1/l[g] for _ in range(l[g])] for _ in range(K[g])] for g in range(ng)]
    alpha_sumk = []
    for g in range(ng):
        alpha_partk = 0.
        for i in range(len(alpha[g])):
            alpha_partk += np.array(alpha[g][i])
        alpha_sumk.append(alpha_partk)

    ymean = [np.mean(y[i]) for i in range(len(y))]
    for i in range(len(y)):
        y[i] = y[i] - np.mean(y[i])
    sigma = np.array([1. for g in range(ng)])
    Xmean = [np.mean(X[g], axis=0) for g in range(ng)]
    for g in range(ng):
        for i in range(nc[g]):
            X[g][i] = np.array(X[g][i]) - Xmean[g]

    scales = np.ones(ng)
    scales_old = np.ones(ng)*0.01
    beta_track = []
    beta_old = np.ones(len(X[0][0][0]))*0.01
    beta_hat = np.ones(len(X[0][0][0]))
    converge = np.ones(ng)
    iter_num = 0
    
    def max_diff(beta_hat, beta_old):
        return np.max([np.abs(beta_hat[i]-beta_old[i])/np.maximum(np.abs(beta_old[i]), 1e-8) for i in range(len(beta_hat))])
    while True:
        ans = [max_diff(beta_hat*scales[i], beta_old*scales_old[i]) for i in range(len(scales)) if converge[i]==1]
        now = datetime.now()
        dt_string = now.strftime("%m/%d/%Y %H:%M:%S")
        print(dt_string, 'Relative change: %.4f' % np.max(ans))
        beta_old = beta_hat
        
        # converge_idx
        converge_idx = np.arange(ng)[converge==1]
        # sigma_g,k
        Sigma = []
        for g in converge_idx:
            Sigmag = []
            for k in range(K[g]):
                Sigmagk = -np.outer(alpha[g][k], alpha[g][k])
                for i in range(l[g]):
                    Sigmagk[i][i] = alpha[g][k][i]*(1-alpha[g][k][i])
                Sigmag.append(Sigmagk)
            Sigma.append(Sigmag)

        G = np.array([[np.dot(alpha_sumk[g], scales[g]*X[g][i]) for i in range(nc[g])] for g in converge_idx])
        #pip = []
        #for g in range(ng):
        #    temp = 1
        #    for k in range(len(alpha[g])):
        #        temp *= (1-np.array(alpha[g][k]))
        #    pip.append(1-temp)
        #G = np.array([[np.sum(X[g][i][np.argsort(pip[g])[-K[g]:]]*np.array(alpha_sumk[g])[np.argsort(pip[g])[-K[g]:]].reshape(-1, 1), axis=0) for i in range(nc[g])] for g in converge_idx])
        #G = np.array([[np.sum(X[g][i][np.argsort(pip[g])[-K[g]:]], axis=0) for i in range(nc[g])] for g in converge_idx])
        #G = np.array([[np.sum([X[g][i][np.argmax(alpha[g][k])] for k in range(K[g])], axis=0) for i in range(nc[g])] for g in converge_idx])
        first, second = 0., 0.
        for j, g in enumerate(converge_idx):
            for i in range(nc[g]):
                soft = 0. 
                for k in range(K[g]):
                    soft += np.dot(np.dot(scales[g]*X[g][i].transpose(), Sigma[j][k]), scales[g]*X[g][i])
                first += (np.outer(G[j][i], G[j][i]) + soft) / sigma[g]**2
                second += y[g][i] * G[j][i] / sigma[g]**2
        
        if not (np.linalg.det(first.astype(float))):
            print('Unable to estimate beta due to non-invertible issue.')
            break
        #logger.info('Ridge Regression')
        beta_hat = np.dot(np.linalg.inv(first.astype(float)+ridge*np.diag(np.ones(len(first)))), second)
        #beta_hat = np.dot(np.linalg.inv(first.astype(float)), second)
        if norm_beta == True:
            beta_hat /= np.sqrt(np.sum(beta_hat**2))
        beta_track.append(beta_hat)

        D = []
        for g in range(ng):
            Dg = []
            for i in range(nc[g]):
                Dg.append(np.dot(X[g][i], beta_hat))
            Dg = np.array(Dg)
            D.append(Dg)
            
        if np.max(ans) < 0.001:
            if allow_scale_change==True and scale_change==False:
                scale_change=True
                print('Scale change')
            else:
                est = estimate(D, ymean, alpha_sumk)
                break
        if iter_num > n_iter:
            print('Not convergent')
            est = estimate(D, ymean, alpha_sumk)
            if ng == 1:
                converge = np.array([-1])
            break
        iter_num += 1
        
        if parallel==True:
            cores = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=cores)
            #num_arrays = len(D)
            #array_size1 = len(D[0])
            #array_size2 = len(D[0][0])
            #shared_D = create_shared_memory(D)
            
            with multiprocessing.pool.Pool(context=ctx) as pool:
                tasks = [(y[g], D[g], sigma[g], alpha[g], pi[g], nc[g], K[g], l[g], scales[g], scale_change, g, converge[g]) for g in range(ng)]
                #tasks = [(y[g], shared_matrix(shared_D, g, array_size1, array_size2), sigma[g], alpha[g], pi[g], nc[g], K[g], l[g], scales[g], scale_change, g, converge[g]) for g in range(ng)]
                res = pool.starmap(mtnm, tasks)
        else:
            res = [mtnm(y[g], D[g], sigma[g], alpha[g], pi[g], nc[g], K[g], l[g], scales[g], scale_change, g, converge[g]) for g in range(ng)]

        scales_old = scales
        alpha_sumk = [ans[0] for ans in res]
        alpha = [ans[1] for ans in res]
        sigma = np.array([ans[2] for ans in res])
        if scale_change == True:
            scales = np.array([ans[3] for ans in res])
        else:
            scales = scales
            scales_old = scales
        converge = np.array([min(converge[a], ans[4]) for a, ans in enumerate(res)])
        if ng == 1:
            converge = np.array([1])
    return alpha, beta_hat, sigma, beta_track, est, K, scales, converge, cal_variance(X, y, K, alpha, sigma, converge, beta_hat)

def create_shared_memory(D):
    flattened_D = np.array(D).ravel()
    shared_D = multiprocessing.Array('d', flattened_D)
    return shared_D

def shared_matrix(shared_D, g, array_size1, array_size2):
    shared_array = np.ctypeslib.as_array(shared_D.get_obj())
    shared_array_g = shared_array[(g*array_size1*array_size2):((g+1)*array_size1*array_size2)].reshape(array_size1, array_size2)
    return shared_array_g

def cal_pip(alpha, ng):
    pip = []
    for g in range(ng):
        ans = 1
        for k in range(len(alpha[g])):
            ans *= (1-np.array(alpha[g][k]))
        pip.append(1-ans)
    #pip = np.array(pip)
    return pip

def expand_feature(feature, distances, trsd):
    idxs = []
    trsd = trsd + [float('inf')]
    for distance in distances:
        idx = 0
        while idx < len(trsd) and distance >= trsd[idx]:
            idx += 1
        idxs.append(idx)
    zeros = [0] * len(feature[0][0])
    new_Xg = []
    for i in range(len(feature)):
        Xgi = [zeros*(idx) + _cCRE + zeros*(len(trsd)-idx-1) for idx, _cCRE in zip(idxs, feature[i])]
        new_Xg.append(Xgi)
    return new_Xg
    
def Linkreg(expression_file, track_files, output, Kg=15, distance=500000):
    genes = process(expression_file, track_files, distance)
    X, y = [], []
    pi, K = [], []
    for g in range(len(genes)):
        gene = genes[g]
        distances = np.abs(np.mean(gene.cCRE_loc, axis=1)-gene.tss)
        prior = np.exp((-distances/100000).astype(float))
        prior /= np.sum(prior)
        #idx = np.where(gene.expression!=0)[0]
        #if len(idx) <= 30:
        #    continue
        if not isinstance(gene.cCRE, list):
            cCRE = gene.cCRE.tolist()
        Xg = [[cCRE[i*gene.cell_num+j] for i in range(int(gene.cCRE_num))] for j in range(gene.cell_num)]
        trsd = [100000]
        Xg = expand_feature(Xg, distances, trsd=trsd)
        X.append(Xg)
        y.append(np.array(gene.expression))
        #y.append(np.array(gene.expression))
        K.append(Kg)
        pi.append([prior.tolist() for _ in range(K[-1])])
    alpha, beta_hat, sigma, beta_track, est, K, scales, converge, variance = model(X, y, K=K, genes=genes, parallel=True, allow_scale_change=False, scale_change=False, norm_beta=False, pi=pi, ridge=0)
    pip = cal_pip(alpha, len(genes))
    for i in range(len(genes)):
        genes[i].pip = pip[i]
        genes[i].alpha = alpha[i]
        genes[i].sigma = sigma[i]
        genes[i].beta = beta_hat
        genes[i].runk = K[i]
        genes[i].scale = scales[i]
        genes[i].converge = converge[i]
        genes[i].std = np.sqrt(np.diag(variance))
        genes[i].trsd = trsd

    results = []
    for gene in genes:
        for i in range(gene.cCRE_num):
            results.append([gene.chromosome, gene.gene_body[0], gene.gene_body[1], gene.ID, gene.strand, gene.cCRE_loc[0], gene.cCRE_loc[1]])
            results[-1] += (gene.cCRE[i*gene.cell_num+np.arange(gene.cell_num), 0]*gene.pip[i]).tolist()
    np.savetxt(output, results, delimiter='\t', fmt='%s')

def main():
    logger = logging.getLogger(__name__)
    parser = argparse.ArgumentParser(description='Linkreg')
    parser.add_argument('--expression_input', type=str, default='')
    parser.add_argument('--tracks_input', type=str, nargs='+', default=[''])
    parser.add_argument('--output', type=str, default='')
    parser.add_argument('--Kg', type=int, default=15)
    parser.add_argument('--distance', type=int, default=500000)
    args = parser.parse_args()
    Linkreg(args.expression_input, args.tracks_input, args.output, args.Kg, args.distance)

if __name__ == '__main__':
    main()
