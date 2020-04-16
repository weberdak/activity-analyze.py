import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import argparse
import matplotlib.gridspec as gridspec

def readData(fname):
    '''Read plate reader file'''
    return pd.read_csv(fname,skiprows=2,sep='\t',header=0)


def readData(fname):
    return pd.read_csv(fname,skiprows=2,sep='\t',header=0)


def linear(x,m,c):
    return m*x + c


def fit_linear(xData,yData):
    fit = optimization.curve_fit(linear, xData, yData)
    return fit[0][0],fit[0][1]


def hill(x,vmax,n,pkca,c):
    #return vmax / ( 1 + 10**(-n*(pkca-x)) )
    return c + (vmax / ( 1 + 10**(-n*(pkca-x)) ))


def fit_hill(xData,yData):
    fit = optimization.curve_fit(hill, xData, yData, np.array([1.0, 2.0, 6.5, 0.001]))
    #fit = optimization.curve_fit(hill, xData, yData, np.array([0.0010, 0.5, 6.0]))
    #return fit[0][0],fit[0][1],fit[0][2]
    return fit[0][0],fit[0][1],fit[0][2],fit[0][3]


# Convert A340 sec-1 to umole/min/mg
def convert(rate,path,epsilon,well_serca,well_vol):
    #path = 0.55 # Path length of absorbance reading (cm)
    #epsilon = 6.22e3 # Exinction coefficient of NADH at 340 nm (M-1 cm-1)
    #well_serca = 0.0044 # mg/mL
    #well_vol = 200 # ul
    return ((rate/(epsilon*path))*60*1000000*(200/1000000))/(well_serca*(well_vol/1000))


def analyze_row(name,plate,row,time_range,exclude,series,interval,path,epsilon,well_serca,well_vol):

    # Read plate file
    print('# Analyzing {}...'.format(name))
    print('# Reading row {} from {}'.format(row,plate))
    data = readData(plate)

    # Set columns to exclude
    badCols = [ row+str(i) for i in exclude ]
    print('# Excluding wells from pKCa fitting: {}'.format(' '.join(badCols)))

    # Set columns to analyze
    colDict = { row+str(n): c for n,c in zip(range(1,13,1), series) }
    colList = []
    for k in colDict.keys():
        if k not in badCols:
            colList.append(k)
    print('# Columns used for fitting: {}'.format(colList))
    pCaList = [ colDict[k] for k in colList ]
    print('# pCa values for wells measured: {}'.format(pCaList))

    # Set up fit linear decays for well time series
    nTimePoints = len(data[colList[0]])
    timeList = np.arange(0,nTimePoints*interval,interval)
    lo = time_range[0]
    hi = time_range[1]
    print('# Fitting decay rates for each well...')
    print('# Measurements detected for each well: {}'.format(nTimePoints))
    print('# Time interval between measurements: {} seconds'.format(interval))
    print('# Only fitting rates from measurement {} to {}'.format(lo, hi))
    
    # Linear fits
    rates = []
    intercepts = []
    for col in colList:
        fit = fit_linear(timeList[lo:hi],data[col][lo:hi])
        rates.append(fit[0]*-1)
        intercepts.append(fit[1])
        
    # Hill fit
    print('# Converting rates from AU/sec to international units (umole/min/mg) for hill fitting')
    print('# Well SERCA: {} mg/mL'.format(well_serca))
    print('# Well volume: {} uL'.format(well_vol))
    print('# Well path length: {} cm'.format(path))
    print('# NADH epsilon: {} M-1 cm-1'.format(epsilon))
    rates_c = [ convert(r,path,epsilon,well_serca,well_vol) for r in rates ]
    fit = fit_hill(pCaList,rates_c)
    vmax = fit[0]
    n = fit[1]
    pkca = fit[2]
    c = fit[3]

    # Show results
    print('# Results for {}:'.format(name))
    print('# Vmax: {0:.4f} umole/min/mg (maximum ATPase rate)'.format(vmax))
    print('# n: {0:.4f} (cooperativity coefficient)'.format(n))
    print('# pKCa: {0:.4f} (pCa at 50% Vmax)'.format(pkca))
    print('# Vmin: {0:.4f} umole/min/mg (minimum ATPase rate)'.format(c))

    # Return all results in list
    results = dict()
    results['name'] = name
    results['vmax'] = np.absolute(vmax)
    results['coop'] = np.absolute(n)
    results['pkca'] = np.absolute(pkca)
    results['vmin'] = np.absolute(c)
    results['pca_list'] = pCaList
    results['rates_au'] = rates
    results['intercepts'] = intercepts
    results['rates_iu'] = rates_c
    results['time'] = timeList
    results['wells'] = colList
    results['absorbance'] = [ data[well] for well in colList ]
    results['lo'] = lo
    results['hi'] = hi
    return results


def blank_results():
    results = dict()
    results['vmax'] = 0
    results['coop'] = 0
    results['pkca'] = 0
    results['vmin'] = 0
    results['pca_list'] = []
    results['rates_au'] = []
    results['intercepts'] = []
    results['rates_iu'] = []
    results['time'] = []
    results['wells'] = []
    results['absorbance'] = []
    return results


def parse_args():

    parser = argparse.ArgumentParser(description='Analyze activity assay.',
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        '--prefix', type=str,
        help='Prefix for output files.',
        default=''
    )
    parser.add_argument(
        '--series', type=float, nargs='+',
        help='pCa value from X1 to X12.',
        default=[ 8.0, 7.5, 7.0, 6.8, 6.6, 6.4, 6.2, 6.0, 5.8, 5.6, 5.4, 5.0 ]
    )
    parser.add_argument(
        '--interval', type=int,
        help='Time interval between measurements in seconds.',
        default=10.0
    )
    parser.add_argument(
        '--rep1_row', type=str,
        help='Row for replicate 1.',
        choices=[ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'N' ], default='N'
    )
    parser.add_argument(
        '--rep1_plate', type=str,
        help='Plate reader file for replicate 1.',
        default=''
    )
    parser.add_argument(
        '--rep1_exclude', type=int, nargs='+',
        help='Exclude these wells from replicate 1.',
        choices=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ]
    )
    parser.add_argument(
        '--rep1_range', type=int, nargs='+',
        help='First and last time point to calculate rate for replicate 1.',
        default=[ 0, 60 ]
    )
    parser.add_argument(
        '--rep2_row', type=str,
        help='Row for replicate 2.',
        choices=[ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'N' ], default='N'
    )
    parser.add_argument(
        '--rep2_plate', type=str,
        help='Plate reader file for replicate 2.',
        default=''
    )
    parser.add_argument(
        '--rep2_exclude', type=int, nargs='+',
        help='Exclude these wells from replicate 2.',
        choices=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ]
    )
    parser.add_argument(
        '--rep2_range', type=int, nargs='+',
        help='First and last time point to calculate rate for replicate 2.',
        default=[ 0, 60]
    )
    parser.add_argument(
        '--rep3_row', type=str,
        help='Row for replicate 3.',
        choices=[ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'N' ], default='N'
    )
    parser.add_argument(
        '--rep3_plate', type=str,
        help='Plate reader file for replicate 3.',
        default=''
    )
    parser.add_argument(
        '--rep3_exclude', type=int, nargs='+',
        help='Exclude these wells from replicate 3.',
        choices=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ]
    )
    parser.add_argument(
        '--rep3_range', type=int, nargs='+',
        help='First and last time point to calculate rate for replicate 3.',
        default=[ 0, 60 ]
    )
    parser.add_argument(
        '--well_serca', type=float,
        help='Concentration of SERCA in well (mg/mL).',
        default=0.0044
    )
    parser.add_argument(
        '--well_vol', type=float,
        help='Volume of well (uL).',
        default=200
    )
    parser.add_argument(
        '--well_eps', type=float,
        help='Exinction coefficient of NADH at 340 nm (M-1 cm-1).',
        default=6.22e3
    )
    parser.add_argument(
        '--well_path', type=float,
        help='Path length through well (cm).',
        default=0.55
    )
    
    args = parser.parse_args()
    return args


def main():

    # Read command line arguments
    args = parse_args()

    results_all = []
    # Analyse replicate 1
    print('# REPLICATE 1')
    print('# -----------')
    if args.rep1_row != 'N':
        result_1 = analyze_row(args.prefix+'.rep1',
                               args.rep1_plate,
                               args.rep1_row,
                               args.rep1_range,
                               args.rep1_exclude,
                               args.series,
                               args.interval,
                               args.well_path,
                               args.well_eps,
                               args.well_serca,
                               args.well_vol
        )
        results_all.append(result_1)
        print('#')
    else:
        result_1 = blank_results()
        print('# Replicate 1 not specified.')
        print('#')
        
    # Analyse replicate 2
    print('# REPLICATE 2')
    print('# -----------')
    if args.rep2_row != 'N':
        result_2 = analyze_row(args.prefix+'.rep2',
                               args.rep2_plate,
                               args.rep2_row,
                               args.rep2_range,
                               args.rep2_exclude,
                               args.series,
                               args.interval,
                               args.well_path,
                               args.well_eps,
                               args.well_serca,
                               args.well_vol
        )
        results_all.append(result_2)
        print('#')
    else:
        result_2 = blank_results()
        print('# Replicate 2 not specified.')
        print('#')
        
    # Analyse replicate 3
    print('# REPLICATE 3')
    print('# -----------')
    if args.rep3_row != 'N':
        result_3 = analyze_row(args.prefix+'.rep3',
                               args.rep3_plate,
                               args.rep3_row,
                               args.rep3_range,
                               args.rep3_exclude,
                               args.series,
                               args.interval,
                               args.well_path,
                               args.well_eps,
                               args.well_serca,
                               args.well_vol
        )
        results_all.append(result_3)
        print('#')
    else:
        result_3 = blank_results()
        print('# Replicate 3 not specified.')
        print('#')

    # Summary
    results = ( result_1, result_2, result_3 )
    print('# SUMMARY')
    print('# -------')
    print('# Rep.\t\tVmax\t\tCoop\t\tpKCa\t\tVmin')
    print('# 1\t\t{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}'.format(result_1['vmax'],
                                                                   result_1['coop'],
                                                                   result_1['pkca'],
                                                                   result_1['vmin']))
    
    print('# 2\t\t{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}'.format(result_2['vmax'],
                                                                   result_2['coop'],
                                                                   result_2['pkca'],
                                                                   result_2['vmin']))

    print('# 3\t\t{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}'.format(result_3['vmax'],
                                                                   result_3['coop'],
                                                                   result_3['pkca'],
                                                                   result_3['vmin']))

    # Get averages
    vmaxs = [ result['vmax'] for result in results if result['vmax'] != 0 ]
    avg_vmax = np.mean(vmaxs)
    std_vmax = np.std(vmaxs)
    coops = [ result['coop'] for result in results if result['coop'] != 0 ]
    avg_coop = np.mean(coops)
    std_coop = np.std(coops)
    pkcas = [ result['pkca'] for result in results if result['pkca'] != 0 ]
    avg_pkca = np.mean(pkcas)
    std_pkca = np.std(pkcas)
    vmins = [ result['vmin'] for result in results if result['vmin'] != 0 ]
    avg_vmin = np.mean(vmins)
    std_vmin = np.std(vmins)
    print('# Avg.\t\t{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}'.format(avg_vmax,
                                                                      avg_coop,
                                                                      avg_pkca,
                                                                      avg_vmin))
    print('# Std.\t\t{0:.4f}\t\t{1:.4f}\t\t{2:.4f}\t\t{3:.4f}'.format(std_vmax,
                                                                      std_coop,
                                                                      std_pkca,
                                                                      std_vmin))
    print('#')

    # Sort through ATPase rates. Fill missing values with 0.
    rates_lists = []
    for result in results:
        rates = []
        if not result['rates_iu']:
            result['rates_iu'] = [ 0 for pca in args.series ]
            result['pca_list'] = [ pca for pca in args.series ]
            
        temp_dict = { pca: rate for pca,rate in zip(result['pca_list'],result['rates_iu']) }
        for pca in args.series:
            if pca not in temp_dict.keys():
                rates.append(0)
            else:
                rates.append(temp_dict[pca])        
        rates_lists.append(rates)

    print('# PKCA CURVES')
    print('# -----------')
    print('# pCa\tRep1\tRep2\tRep3\tAvg\tStd\tAvg_norm\tStd_norm')
    i=0
    avgs = []
    stds = []
    avgs_norm = []
    stds_norm = []
    for pca in args.series:
        # Get average. Remove null values.
        l = range(3)
        temp_r = [ rates_lists[n][i] for n in l if rates_lists[n][i] != 0]
        temp_rn = [ (r-avg_vmin)/avg_vmax for r in temp_r ]
        
        avgs.append(np.mean(temp_r))
        stds.append(np.std(temp_r))
        avgs_norm.append(np.mean(temp_rn))
        stds_norm.append(np.std(temp_rn))
        
        print('{0}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}'.format(pca,
                                                                                          rates_lists[0][i],
                                                                                          rates_lists[1][i],
                                                                                          rates_lists[2][i],
                                                                                          avgs[i],
                                                                                          stds[i],
                                                                                          avgs_norm[i],
                                                                                          stds_norm[i]))
        i+=1

        
    # Simulate fits
    sim_xlist = [ x for x in np.arange(9,4,-0.1)]
    sim_ylists = []
    for result in results:
        if result['vmax'] != 0:
            s = [ hill(x,result['vmax'],result['coop'],result['pkca'],result['vmin']) for x in sim_xlist ]
        else:
            s = [ 0 for x in sim_xlist ]
        sim_ylists.append(s)

    # Simulate average fit
    s = [ hill(x,avg_vmax,avg_coop,avg_pkca,avg_vmin) for x in sim_xlist ]
    s_n = [ (y-avg_vmin)/avg_vmax for y in s ]
    sim_ylists.append(s)
    sim_ylists.append(s_n)

    
    # Plot all data
    fig = plt.figure(figsize=(6,6))
    gs = gridspec.GridSpec(3,2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[2,0])
    ax6 = fig.add_subplot(gs[2,1])
    #ax7 = fig.add_subplot(gs[3,0])
    #ax8 = fig.add_subplot(gs[3,1])

    j = -1
    for ax_r,ax_h,result in ((ax1,ax2,result_1),(ax3,ax4,result_2),(ax5,ax6,result_3)):
        j+=1
        
        if result['absorbance']:
            for i,c in enumerate(result['wells']):

                # Plot Absorbance vs. time
                ax_r.plot(result['time'],result['absorbance'][i],c='black',linewidth=1.0)
                ax_r.plot(result['time'],[linear(x,result['rates_au'][i]*-1,result['intercepts'][i]) for x in result['time']],
                         linestyle='--',color='red',linewidth=0.5)
                ax_r.axvspan(result['time'][result['lo']],
                             result['time'][result['hi']], alpha=0.1, color='grey')
                ax_r.xaxis.set_tick_params(labelsize=6)
                ax_r.yaxis.set_tick_params(labelsize=6)
                ax_r.set_xlabel('Time (seconds)',size=6)
                ax_r.set_ylabel('A340',size=6)
                ax_r.set_title('{}'.format(result['name']),size=6,loc='left')

            # Plot Hill
            # Average
            ax_h.scatter(args.series,avgs, c='black',
                         marker='o',s=10, linewidths=0.5,
                         edgecolor='black')
            ax_h.plot(sim_xlist,sim_ylists[3],linestyle='--',color='black',linewidth=0.5)
            ax_h.errorbar(args.series,avgs,yerr=stds,linewidth=0.0,ecolor='black', elinewidth=0.5, capsize=1)
            
            # Replicate
            ax_h.scatter(result['pca_list'],result['rates_iu'], c='red',
                         marker='o',s=10, linewidths=0.5,
                         edgecolor='black')
            ax_h.plot(sim_xlist,sim_ylists[j],linestyle='--',color='red',linewidth=0.5)

            
            ax_h.xaxis.set_tick_params(labelsize=6)
            ax_h.yaxis.set_tick_params(labelsize=6)
            ax_h.set_xlim(8.5,4.5)
            ax_h.set_xlabel('pCa',size=6)
            ax_h.set_ylabel('Activity ($\mu$mole/mg/min)',size=6)
            ax_h.set_title('{}'.format(result['name']),size=6,loc='left')

    plt.tight_layout()
    print('#\n# Outputing figure to {}.jpg'.format(args.prefix))
    plt.savefig('{}.jpg'.format(args.prefix),dpi=300)
    plt.show()
        
if __name__ == '__main__':
    main()
