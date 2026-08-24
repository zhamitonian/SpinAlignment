import pandas as pd

df = pd.read_csv("./fit_results/data_fit/nsig_results.csv")
total_nsig = df['nsig'].sum()
print("total signal yield", total_nsig)
1,224,727.7655 