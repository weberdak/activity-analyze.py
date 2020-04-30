import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import argparse
import matplotlib.gridspec as gridspec


def log(text,f):
    '''Very simple logger. Writes to file object and terminal.
    Parameters
    ----------
    text: str
       Text to write to file.
    f: file object
       File to write lof text to.
    '''
    f.write(text+'\n')
    print(text)

    
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


def analyze_row(name,plate,row,time_range,exclude,series,interval,path,epsilon,well_serca,well_vol,f):

    # Read plate file
    log('# Analyzing {}...'.format(name),f)
    log('# Reading row {} from {}'.format(row,plate),f)
    data = readData(plate)

    # Set columns to exclude
    badCols = [ row+str(i) for i in exclude ]
    log('# Excluding wells from pKCa fitting: {}'.format(' '.join(badCols)),f)

    # Set columns to analyze
    colDict = { row+str(n): c for n,c in zip(range(1,13,1), series) }
    colList = []
    for k in colDict.keys():
        if k not in badCols:
            colList.append(k)
    log('# Columns used for fitting: {}'.format(colList),f)
    pCaList = [ colDict[k] for k in colList ]
    log('# pCa values for wells measured: {}'.format(pCaList),f)

    # Set up fit linear decays for well time series
    nTimePoints = len(data[colList[0]])
    timeList = np.arange(0,nTimePoints*interval,interval)
    lo = time_range[0]
    hi = time_range[1]
    log('# Fitting decay rates for each well...',f)
    log('# Measurements detected for each well: {}'.format(nTimePoints),f)
    log('# Time interval between measurements: {} seconds'.format(interval),f)
    log('# Only fitting rates from measurement {} to {}'.format(lo, hi),f)
    
    # Linear fits
    rates = []
    intercepts = []
    for col in colList:
        fit = fit_linear(timeList[lo:hi],data[col][lo:hi])
        rates.append(fit[0]*-1)
        intercepts.append(fit[1])
        
    # Hill fit
    log('# Converting rates from AU/sec to international units (umole/min/mg) for hill fitting',f)
    log('# Well SERCA: {} mg/mL'.format(well_serca),f)
    log('# Well volume: {} uL'.format(well_vol),f)
    log('# Well path length: {} cm'.format(path),f)
    log('# NADH epsilon: {} M-1 cm-1'.format(epsilon),f)
    rates_c = [ convert(r,path,epsilon,well_serca,well_vol) for r in rates ]
    fit = fit_hill(pCaList,rates_c)
    vmax = fit[0]
    n = fit[1]
    pkca = fit[2]
    c = fit[3]

    # Show results
    log('# Results for {}:'.format(name),f)
    log('# Vmax: {0:.4f} umole/min/mg (maximum ATPase rate)'.format(vmax),f)
    log('# n: {0:.4f} (cooperativity coefficient)'.format(n),f)
    log('# pKCa: {0:.4f} (pCa at 50% Vmax)'.format(pkca),f)
    log('# Vmin: {0:.4f} umole/min/mg (minimum ATPase rate)'.format(c),f)

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
        '--rep4_row', type=str,
        help='Row for replicate 4.',
        choices=[ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'N' ], default='N'
    )
    parser.add_argument(
        '--rep4_plate', type=str,
        help='Plate reader file for replicate 4.',
        default=''
    )
    parser.add_argument(
        '--rep4_exclude', type=int, nargs='+',
        help='Exclude these wells from replicate 4.',
        choices=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ]
    )
    parser.add_argument(
        '--rep4_range', type=int, nargs='+',
        help='First and last time point to calculate rate for replicate 4.',
        default=[ 0, 60 ]
    )
    parser.add_argument(
        '--rep5_row', type=str,
        help='Row for replicate 5.',
        choices=[ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'N' ], default='N'
    )
    parser.add_argument(
        '--rep5_plate', type=str,
        help='Plate reader file for replicate 5.',
        default=''
    )
    parser.add_argument(
        '--rep5_exclude', type=int, nargs='+',
        help='Exclude these wells from replicate 5.',
        choices=[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ]
    )
    parser.add_argument(
        '--rep5_range', type=int, nargs='+',
        help='First and last time point to calculate rate for replicate 5.',
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

    # Start log file
    f = open(args.prefix+'.log', 'w')
    
    results_all = []
    # Analyse replicate 1
    log('# REPLICATE 1',f)
    log('# -----------',f)
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
                               args.well_vol,
                               f
        )
        results_all.append(result_1)
        log('#',f)
    else:
        result_1 = blank_results()
        log('# Replicate 1 not specified.',f)
        log('#',f)
        
    # Analyse replicate 2
    log('# REPLICATE 2',f)
    log('# -----------',f)
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
                               args.well_vol,
                               f
        )
        results_all.append(result_2)
        log('#',f)
    else:
        result_2 = blank_results()
        log('# Replicate 2 not specified.',f)
        log('#',f)
        
    # Analyse replicate 3
    log('# REPLICATE 3',f)
    log('# -----------',f)
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
                               args.well_vol,
                               f
        )
        results_all.append(result_3)
        log('#',f)
    else:
        result_3 = blank_results()
        log('# Replicate 3 not specified.',f)
        log('#',f)

    # Analyse replicate 4
    log('# REPLICATE 4',f)
    log('# -----------',f)
    if args.rep4_row != 'N':
        result_4 = analyze_row(args.prefix+'.rep4',
                               args.rep4_plate,
                               args.rep4_row,
                               args.rep4_range,
                               args.rep4_exclude,
                               args.series,
                               args.interval,
                               args.well_path,
                               args.well_eps,
                               args.well_serca,
                               args.well_vol,
                               f
        )
        results_all.append(result_4)
        log('#',f)
    else:
        result_4 = blank_results()
        log('# Replicate 4 not specified.',f)
        log('#',f)

    # Analyse replicate 3
    log('# REPLICATE 5',f)
    log('# -----------',f)
    if args.rep5_row != 'N':
        result_5 = analyze_row(args.prefix+'.rep5',
                               args.rep5_plate,
                               args.rep5_row,
                               args.rep5_range,
                               args.rep5_exclude,
                               args.series,
                               args.interval,
                               args.well_path,
                               args.well_eps,
                               args.well_serca,
                               args.well_vol,
                               f
        )
        results_all.append(result_5)
        log('#',f)
    else:
        result_5 = blank_results()
        log('# Replicate 5 not specified.',f)
        log('#',f)
    
    # Summary
    results = ( result_1, result_2, result_3, result_4, result_5 )
    log('# SUMMARY',f)
    log('# -------',f)
    log('# Rep.\tVmax\tCoop\tpKCa\tVmin',f)
    log('# 1\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(result_1['vmax'],
                                                                   result_1['coop'],
                                                                   result_1['pkca'],
                                                                   result_1['vmin']),f)
    
    log('# 2\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(result_2['vmax'],
                                                                   result_2['coop'],
                                                                   result_2['pkca'],
                                                                   result_2['vmin']),f)

    log('# 3\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(result_3['vmax'],
                                                                   result_3['coop'],
                                                                   result_3['pkca'],
                                                                   result_3['vmin']),f)

    log('# 4\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(result_4['vmax'],
                                                                   result_4['coop'],
                                                                   result_4['pkca'],
                                                                   result_4['vmin']),f)

    log('# 5\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(result_5['vmax'],
                                                                   result_5['coop'],
                                                                   result_5['pkca'],
                                                                   result_5['vmin']),f)
    
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
    log('# Avg.\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(avg_vmax,
                                                                      avg_coop,
                                                                      avg_pkca,
                                                                      avg_vmin),f)
    log('# Std.\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(std_vmax,
                                                                      std_coop,
                                                                      std_pkca,
                                                                      std_vmin),f)
    log('#',f)

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

    log('# PKCA CURVES',f)
    log('# -----------',f)
    log('# pCa\tRep1\tRep2\tRep3\tRep4\tRep5\tAvg\tStd\tAvg_norm\tStd_norm',f)
    i=0
    avgs = []
    stds = []
    avgs_norm = []
    stds_norm = []
    for pca in args.series:
        # Get average. Remove null values.
        l = range(5)
        temp_r = [ rates_lists[n][i] for n in l if rates_lists[n][i] != 0]
        temp_rn = [ (r-avg_vmin)/avg_vmax for r in temp_r ]
        
        avgs.append(np.mean(temp_r))
        stds.append(np.std(temp_r))
        avgs_norm.append(np.mean(temp_rn))
        stds_norm.append(np.std(temp_rn))
        
        log('{0}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}'.format(pca,
                                                                                                          rates_lists[0][i],
                                                                                                          rates_lists[1][i],
                                                                                                          rates_lists[2][i],
                                                                                                          rates_lists[3][i],
                                                                                                          rates_lists[4][i],
                                                                                                          avgs[i],
                                                                                                          stds[i],
                                                                                                          avgs_norm[i],
                                                                                                          stds_norm[i]),f)
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
    
    # Write fits to file
    log('#\n# Writing fits to {}.fit.txt.'.format(args.prefix),f)
    ff = open(args.prefix+'.fit.txt', 'w')
    ff.write('# pCa\tRep1\tRep2\tRep3\tRep4\tRep5\tAvg\tNorm\n')
    i=0
    for x in sim_xlist:
        ff.write('{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\n'.format(x,
                                                                                                   sim_ylists[0][i],
                                                                                                   sim_ylists[1][i],
                                                                                                   sim_ylists[2][i],
                                                                                                   sim_ylists[3][i],
                                                                                                   sim_ylists[4][i],
                                                                                                   sim_ylists[5][i],
                                                                                                   sim_ylists[6][i]))
        i+=1
    ff.close()
    
    # Plot all data
    j = -1
    for result in (result_1,result_2,result_3,result_4,result_5):
        j+=1

        # Only output plot if there is data
        if result['absorbance']:
            fig = plt.figure(figsize=(6,2))
            gs = gridspec.GridSpec(1,2)
            ax_r = fig.add_subplot(gs[0,0])
            ax_h = fig.add_subplot(gs[0,1])
            
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
            ax_h.plot(sim_xlist,sim_ylists[5],linestyle='--',color='black',linewidth=0.5)
            ax_h.errorbar(args.series,avgs,yerr=stds,linewidth=0.0,ecolor='black', elinewidth=0.5, capsize=1)
            
            # Replicate
            ax_h.scatter(result['pca_list'],result['rates_iu'], c='red',
                         marker='o',s=10, linewidths=0.5,
                         edgecolor='black')
            ax_h.plot(sim_xlist,sim_ylists[j],linestyle='--',color='red',linewidth=0.5)

            
            ax_h.xaxis.set_tick_params(labelsize=8)
            ax_h.yaxis.set_tick_params(labelsize=8)
            ax_h.set_xlim(8.5,4.5)
            ax_h.set_xlabel('pCa',size=8)
            ax_h.set_ylabel('Activity ($\mu$mole/mg/min)',size=8)
            ax_h.set_title('{}'.format(result['name']),size=8,loc='left')

            # Output replicate plot to file
            plt.tight_layout()
            log('# Outputting figure to {}.rep{}.jpg'.format(args.prefix,j+1),f)
            plt.savefig('{}.rep{}.jpg'.format(args.prefix,j+1),dpi=300)
            
    f.close()
    
        
if __name__ == '__main__':
    main()
