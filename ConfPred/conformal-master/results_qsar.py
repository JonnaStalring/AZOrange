import pickle

from cp.nonconformity import ClassNC, RegrNC

from evaluate_nc import *

def report(args):
    lines = []
    header_data = ''
    header_prop = ''
    for id, arg in enumerate(args):
        nc_str, dataset, eps = arg
        nc = eval(nc_str)
        line = dataset+'\t'+str(eps)+'\t'+nc_str
        header_prop = ''
        with open('results/qsar/%d.p' % id, 'rb') as f:
            res = pickle.load(f)
            r, rep = res
            if isinstance(nc, ClassNC):
                sc = (1-r.accuracy(), r.singleton_criterion(), r.confidence(), rep, r.time())
                header_prop = '\t\t\terr\tsingleton\tconfidence\trep\ttime'
            else:
                sc = (1-r.accuracy(), r.median_range(), r.interdecile_mean(), rep, r.time())
                header_prop = '\t\t\terr\tmedian\tinterdecile\trep\ttime'
            line += ''.join(['\t%.3f' % s for s in sc])
        lines.append(header_prop+'\n')
        lines.append(line+'\n')

    with open('results/qsar/results.txt', 'w') as f:
        f.writelines(lines)

if __name__ == '__main__':
    with open('results/qsar/args.p', 'rb') as f:
        args = pickle.load(f)
    report(args)
