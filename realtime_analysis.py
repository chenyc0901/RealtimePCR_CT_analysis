
## Usage: python realtime_analysis.py control_gene control_condition csv_file
## column must be 'Treatment' 'Gene' 'CT' 
## Replicate = 3

import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import statannot

control_gene = sys.argv[1]
control_treatment = sys.argv[2]
csvfile = sys.argv[3]

df = pd.read_csv(csvfile, sep=',')
final_data = pd.DataFrame()
df['index'] = [i % 3 for i in range(len(df))]

def FCcalculator(control_gene, condition_gene, control_treatment):
    df1 = df.loc[df.Gene.isin([control_gene,condition_gene])].pivot_table(index=['Treatment','index'],columns=['Gene'], values='CT')
    df1['Target_gene'] = condition_gene
    df1['deltaCT'] = df1[condition_gene] - df1[control_gene]
    df1['deltadeltaCT'] = df1['deltaCT'] - df1.groupby('Treatment').mean('deltaCT').loc[control_treatment,].deltaCT
    df1['FC'] = 2**-(df1['deltadeltaCT'])
    df1 = df1.reset_index()
    result = df1[['Target_gene','Treatment','FC']]
    return result

for i in set(df.Gene):
    if i == control_gene:
        continue
    else:
        condition_gene = i
        result = FCcalculator(control_gene, condition_gene, control_treatment) 
        final_data = pd.concat([final_data, result],ignore_index=True)
print(final_data)

final_data.to_csv('FC_result.csv',index=False)

# use seaborn to draw plot
sns.set(style='whitegrid', font_scale=1.2)
ax = sns.barplot(x='Target_gene', y='FC', hue='Treatment', data=final_data, errorbar='sd')
ax.set(xlabel='Gene', ylabel='Fold Change', title='')
ax.legend(title='Treatment', loc='upper right')



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
