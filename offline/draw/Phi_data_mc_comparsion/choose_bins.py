import pandas as pd
import numpy as np

txt_file = "./fit_result/nsig_results.txt"
                 
df = pd.read_csv(txt_file, sep='\s+', comment='#',
                names=['z_center', 'z_width', 'pt_center', 'pt_width', 
                        'cos_center', 'cos_width', 'nsig', 'nsig_err', 'nsig_err2'])

print(f"Loaded {len(df)} bins (before filtering)")

# 添加原始 bin index
df['bin_index'] = df.index

# 过滤掉全0的行（无效数据）
df = df[(df['z_center'] != 0) | (df['pt_center'] != 0) | (df['nsig'] != 0)]
print(f"After removing all-zero rows: {len(df)} bins")

# 按 (z_center, pt_center) 分组并计算总事例数
df_grouped = df.groupby(['z_center', 'pt_center']).agg(
    total_nsig=('nsig', 'sum'),
    n_bins=('nsig', 'count')
).reset_index()

print(f"\nFound {len(df_grouped)} unique (z, pt) subgroups")
print(f"Before filtering (nsig >= 1000):")
print(df_grouped)

# 过滤出总事例数 >= 1000 的子组
df_valid_groups = df_grouped[df_grouped['total_nsig'] >= 1000]
print(f"\nAfter filtering (total nsig >= 1000): {len(df_valid_groups)} subgroups")
print(df_valid_groups)

# 只保留在有效子组中的bins
df_filtered = df.merge(
    df_valid_groups[['z_center', 'pt_center']], 
    on=['z_center', 'pt_center'], 
    how='inner'
)

print(f"\nFinal filtered bins: {len(df_filtered)} bins")
print(f"Retained bin indices: {sorted(df_filtered['bin_index'].tolist())}")

# 保存结果
output_file = "./fit_result/nsig_results_filtered.csv"
df_filtered.to_csv(output_file, index=False)
print(f"\nSaved filtered results to {output_file}")

# 显示统计信息
print("\n=== Statistics ===")
print(f"Original bins: {len(df)}")
print(f"Filtered bins: {len(df_filtered)}")
print(f"Removed bins: {len(df) - len(df_filtered)}")
print(f"\nFiltered by (z, pt) subgroup:")
for _, row in df_valid_groups.iterrows():
    print(f"  z={row['z_center']:.4f}, pt={row['pt_center']:.4f}: {int(row['total_nsig'])} events in {int(row['n_bins'])} bins")