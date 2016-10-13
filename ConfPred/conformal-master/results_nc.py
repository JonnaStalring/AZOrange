import pickle

from cp.nonconformity import ClassNC, RegrNC

from evaluate_nc import *

def report(args, nc_type, rng=None):
    lines = []
    header_data = ''
    header_prop = ''
    for id, nc_str in enumerate(args):
        nc = eval(nc_str)
        line = nc_str
        if not isinstance(nc, nc_type): continue
        with open('results/nc/%d.p' % id, 'rb') as f:
            res = pickle.load(f)
            header_data = ''
            header_prop = ''
            for dataset in sorted(res):
                r = res[dataset]
                if isinstance(nc, ClassNC):
                    header_data += '\t%s\t\t\t' % dataset
                    sc = (1-r.accuracy(), r.singleton_criterion(), r.confidence(), r.time())
                    header_prop += '\terr\tsingleton\tconfidence\ttime'
                else:
                    header_data += '\t%s\t\t\t\t\t' % dataset
                    sc = (1-r.accuracy(),
                          r.median_range(), r.median_range()/rng[dataset],
                          r.interdecile_mean(), r.interdecile_mean()/rng[dataset],
                          r.time())
                    header_prop += '\terr\tmedian\tmed_norm\tinterdecile\tinter_norm\ttime'
                line += ''.join(['\t%.3f' % s for s in sc])
        lines.append(line+'\n')

    if nc_type == ClassNC: out = 'results/nc/class.txt'
    else: out = 'results/nc/regr.txt'
    with open(out, 'w') as f:
        f.write(header_data+'\n')
        f.write(header_prop+'\n')
        f.writelines(lines)

if __name__ == '__main__':
    with open('results/nc/args.p', 'rb') as f:
        args = pickle.load(f)
    rng = {}
    for dataset in dataset_regr:
        tab = Orange.data.Table(dataset)
        dataset_id = dataset.split('/')[-1].split('.')[0]
        rng[dataset_id] = max(tab.Y)-min(tab.Y)
    report(args, ClassNC)
    report(args, RegrNC, rng)
