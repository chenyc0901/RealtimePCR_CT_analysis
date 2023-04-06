
## column index must be 'Treatment' 'Gene' 'CT' 
## Replicate default = 3

import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import statannot
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()
    parser.add_argument("-f", help="csv file")
    parser.add_argument("-c", help="control gene")
    parser.add_argument("-t",help="control treatment")
    parser.add_argument("-rep", type=int, default=3, help="replicate number")
    args = parser.parse_args()
    control_gene = args.c
    control_treatment = args.t
    csvfile = args.f
    rep = args.rep

    df = pd.read_csv(csvfile, sep=',')
    check_list = ['Treatment','Gene','CT']
    #to check column names are correct
    df_cols = df.columns.to_list()
    assert all(i in df_cols for i in check_list) , "Columns are wrong: must contain {0}".format(check_list)
    #to index the number of replicates
    df = df.sort_values(['Treatment','Gene','CT'])
    df['index'] = [i % rep for i in range(len(df))]
    final_data = FCcalculator_per_gene(df, control_gene, control_treatment)
    draw_plot(final_data)

def FCcalculator(control_gene, condition_gene, control_treatment,df):
    df1 = df.loc[df.Gene.isin([control_gene,condition_gene])].pivot_table(index=['Treatment','index'],columns=['Gene'], values='CT')
    df1['Target_gene'] = condition_gene
    df1['deltaCT'] = df1[condition_gene] - df1[control_gene]
    df1['deltadeltaCT'] = df1['deltaCT'] - df1.groupby('Treatment').mean('deltaCT').loc[control_treatment,].deltaCT
    df1['FC'] = 2**-(df1['deltadeltaCT'])
    df1 = df1.reset_index()
    result = df1[['Target_gene','Treatment','FC']]
    return result

def FCcalculator_per_gene(df, control_gene, control_treatment):
    final_data = pd.DataFrame()
    for i in set(df.Gene):
        if i == control_gene:
            continue
        else:
            condition_gene = i
            result = FCcalculator(control_gene, condition_gene, control_treatment,df) 
            final_data = pd.concat([final_data, result],ignore_index=True)
    final_data.to_csv('FC_result.csv',index=False)
    print(final_data)
    return final_data
    
def draw_plot(final_data):
    # use seaborn to draw plot
    plt.figure(figsize = (15,8) )
    sns.set(style='whitegrid', font_scale=1.2)
    ax = sns.barplot(x='Target_gene', y='FC', hue='Treatment', data=final_data, errorbar='sd')
    ax.set(xlabel='Gene', ylabel='Fold Change', title='')
    ax.legend(title='Treatment', loc='upper right')
    ax.tick_params(axis = 'x', rotation = 45)

    ## Add statistical result to plot
    tuple_list = []
    for name in set(final_data.Target_gene):
        list = tuple((name, i) for i in set(final_data[final_data.Target_gene==name].Treatment))
        tuple_list.append(list)

    statannot.add_stat_annotation(
        ax,
        data = final_data,
        x = 'Target_gene',
        y = 'FC',
        hue = 'Treatment',
        test = 't-test_ind',
        box_pairs = tuple_list,
        text_format = 'star',
        loc = 'outside',
    )
    plt.savefig('plot_result.pdf')

if __name__ == "__main__":
    main()
