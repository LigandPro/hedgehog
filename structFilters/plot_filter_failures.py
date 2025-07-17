import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams.update({'font.size': 16,
                     'axes.titlesize': 22,
                     'axes.labelsize': 18,
                     'xtick.labelsize': 14,
                     'ytick.labelsize': 14,
                     'legend.fontsize': 14,
                     'figure.titlesize': 24
                   })

def analyze_filter_failures(file_path):
    print(f"Reading data from: {file_path}")
    df = pd.read_csv(file_path, low_memory=False)
    
    filter_columns = [col for col in df.columns if col.startswith('pass_filter_')]
    
    filter_failures = {}
    filter_reasons = {}
    all_detailed_reasons = {}
    
    for col in filter_columns:
        filter_name = col.replace('pass_filter_', '')
        
        failures = (df[col] == False).sum()
        total = len(df)
        failure_percentage = (failures / total) * 100
        
        filter_failures[filter_name] = {'failures': failures,
                                        'total': total,
                                        'percentage': failure_percentage
                                       }
        
        reasons_col = f'reasons_{filter_name}'
        if reasons_col in df.columns:
            failed_molecules = df[df[col] == False]
            reasons_data = failed_molecules[reasons_col].dropna()
            
            reason_counts = {}
            for reasons_str in reasons_data:
                if pd.notna(reasons_str) and str(reasons_str).strip():
                    individual_reasons = [r.strip() for r in str(reasons_str).split(';') if r.strip()]
                    for reason in individual_reasons:
                        reason_counts[reason] = reason_counts.get(reason, 0) + 1
            
            sorted_reasons = sorted(reason_counts.items(), key=lambda x: x[1], reverse=True)
            filter_reasons[filter_name] = sorted_reasons
            all_detailed_reasons[filter_name] = reason_counts
    
    create_main_filter_plot(filter_failures, file_path)
    create_individual_filter_plots(filter_failures, filter_reasons, file_path)
    create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path)
    
    print("\nFilter Failure Summary (sorted by failures):")
    print("-" * 80)
    for filter_name, stats in sorted(filter_failures.items(), key=lambda x: x[1]['failures'], reverse=True):
        print(f"{filter_name:25}: {stats['failures']:6d} failures ({stats['percentage']:5.1f}%)")
    
    print(f"\nTotal molecules analyzed: {len(df)}")
    max_failures = max(filter_failures.values(), key=lambda x: x['failures'])
    min_failures = min(filter_failures.values(), key=lambda x: x['failures'])
    print(f"Most restrictive filter: {[k for k, v in filter_failures.items() if v == max_failures][0]} ({max_failures['failures']} failures)")
    print(f"Least restrictive filter: {[k for k, v in filter_failures.items() if v == min_failures][0]} ({min_failures['failures']} failures)")
    
    print("\n" + "="*80)
    print("DETAILED REASONS FOR ALL FILTERS WITH FAILURES")
    print("="*80)
    
    for filter_name, stats in sorted(filter_failures.items(), key=lambda x: x[1]['failures'], reverse=True):
        if stats['failures'] > 0 and filter_name in filter_reasons and filter_reasons[filter_name]:
            print(f"\n{filter_name.upper()} ({stats['failures']} failures):")
            print("-" * 60)
            
            for i, (reason, count) in enumerate(filter_reasons[filter_name], 1):
                percentage = (count / stats['failures']) * 100
                print(f"{i:3d}. {reason:<50} : {count:5d} ({percentage:4.1f}%)")
            
            if i < len(filter_reasons[filter_name]):
                print()
    
    create_complete_reasons_breakdown(all_detailed_reasons, filter_failures, file_path)
    create_comprehensive_overview(filter_reasons, filter_failures, file_path)
    create_summary_table(filter_failures, filter_reasons, file_path)
    
    return filter_failures, filter_reasons, all_detailed_reasons

def create_main_filter_plot(filter_failures, file_path):
    plot_data = []
    for filter_name, stats in filter_failures.items():
        plot_data.append({'filter': filter_name,
                          'failures': stats['failures'],
                          'percentage': stats['percentage']
                        })
    
    plot_df = pd.DataFrame(plot_data)
    plot_df = plot_df.sort_values('failures', ascending=False)  

    plt.figure(figsize=(max(16, len(plot_df) * 0.6), 16))
    bars = plt.bar(range(len(plot_df)), plot_df['failures'], color='steelblue', alpha=0.8, width=0.3)
    
    plt.xlabel('Filters', fontsize=20)
    plt.ylabel('Number of Molecules Failed', fontsize=20)
    plt.title('Number of Molecules Failed by Each Filter', fontsize=26, fontweight='bold')
    plt.xticks(range(len(plot_df)), plot_df['filter'], rotation=45, ha='right', fontsize=16)
    
    for i, (_, row) in enumerate(plot_df.iterrows()):
        plt.text(i, row['failures'] + max(plot_df['failures']) * 0.01, f"{row['failures']}\n({row['percentage']:.1f}%)", ha='center', va='bottom', fontsize=14)
    
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'filter_failures_plot.png')
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Main filter plot saved to: {output_path}")
    plt.close()

def create_individual_filter_plots(filter_failures, filter_reasons, file_path):
    print(f"Creating individual plots for all filters with failures...")
    
    for filter_name, stats in filter_failures.items():
        if stats['failures'] > 0:
            reasons_data = filter_reasons.get(filter_name, [])
            
            if not reasons_data:
                continue
                
            plot_data = []
            for reason, count in reasons_data:
                plot_data.append({'Reason': reason,
                                  'Count': count,
                                  'Percentage_of_Filter_Failures': (count / stats['failures']) * 100 if stats['failures'] > 0 else 0
                                })
            
            if not plot_data:
                continue
            
            plot_df = pd.DataFrame(plot_data)
            plot_df = plot_df.sort_values('Count', ascending=False) 
            
            plt.figure(figsize=(max(16, len(plot_df) * 0.6), 20))
            bars = plt.bar(range(len(plot_df)), plot_df['Count'], color='steelblue', alpha=0.8, width=0.3)
            
            plt.xlabel('Failure Reasons', fontsize=20)
            plt.ylabel('Number of Molecules Failed', fontsize=20)
            plt.title(f'{filter_name.upper()} - Failure Reasons ({len(plot_df)} reasons, {stats["failures"]} total failures)', fontsize=26, fontweight='bold')
            plt.xticks(range(len(plot_df)), plot_df['Reason'], rotation=45, ha='right', fontsize=max(10, min(16, 300 // len(plot_df))))
            
            for i, (bar, count) in enumerate(zip(bars, plot_df['Count'])):
                plt.text(i, count + max(plot_df['Count']) * 0.01, f"{count}\n({plot_df.iloc[i]['Percentage_of_Filter_Failures']:.1f}%)", ha='center', va='bottom', fontsize=12)
            
            plt.grid(axis='y', alpha=0.3)
            plt.tight_layout()

            path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
            os.makedirs(path_to_save, exist_ok=True)
            output_path = os.path.join(path_to_save,  f'{filter_name}_reasons_plot.png')
            plt.savefig(output_path, dpi=600, bbox_inches='tight')
            print(f"Individual filter plot for {filter_name} saved to: {output_path}")
            plt.close()

def create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path):
    sorted_filters = sorted(filter_failures.items(), key=lambda x: x[1]['failures'], reverse=True)
    sorted_filters = [(name, stats) for name, stats in sorted_filters if stats['failures'] > 0]
    
    num_filters = len(sorted_filters)
    if num_filters == 0:
        print("No filters with failures found for multi-panel plot")
        return
    
    if num_filters <= 3:
        rows = 1
        cols = num_filters
    elif num_filters <= 6:
        rows = 2
        cols = 3
    elif num_filters <= 9:
        rows = 3
        cols = 3
    elif num_filters <= 12:
        rows = 3
        cols = 4
    else:
        cols = 4
        rows = (num_filters + cols - 1) // cols  
    
    plt.figure(figsize=(cols * 6, rows * 6))
    
    for i, (filter_name, stats) in enumerate(sorted_filters):
        all_reasons_data = filter_reasons.get(filter_name, [])
        
        plt.subplot(rows, cols, i + 1)
        reason_names = [r[0] for r in all_reasons_data]
        reason_counts = [r[1] for r in all_reasons_data]
        
        if len(reason_names) > 10:
            reason_names = reason_names[:10]
            reason_counts = reason_counts[:10]
            title_suffix = f"(Top 10 of {len(all_reasons_data)} reasons)"
        else:
            title_suffix = f"({len(all_reasons_data)} reasons)"
        
        plt.bar(range(len(reason_names)), reason_counts, color='steelblue', alpha=0.8, width=0.3)
        
        plt.xlabel('Reasons', fontsize=14)
        plt.ylabel('Molecules Failed', fontsize=14)
        plt.title(f'{filter_name.upper()}\n{title_suffix}', fontsize=16, fontweight='bold')
        
        truncated_names = []
        for name in reason_names:
            if len(name) > 15:
                truncated_names.append(name[:12] + "...")
            else:
                truncated_names.append(name)
        
        plt.xticks(range(len(truncated_names)), truncated_names, rotation=45, ha='right', fontsize=12)
        
        for j, count in enumerate(reason_counts):
            percentage = (count / stats['failures']) * 100 if stats['failures'] > 0 else 0
            plt.text(j, count + max(reason_counts) * 0.01 if reason_counts else 0, f"{count}", ha='center', va='bottom', fontsize=11)
        
        plt.grid(axis='y', alpha=0.3)
    
    plt.tight_layout(h_pad=1.5, w_pad=1.0)
    
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'all_filters_reasons_plot.png')
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Multi-panel plot for all filters saved to: {output_path}")
    plt.close()

def create_complete_reasons_breakdown(all_detailed_reasons, filter_failures, file_path):
    breakdown_data = []
    
    for filter_name, reasons_dict in all_detailed_reasons.items():
        total_failures = filter_failures[filter_name]['failures']
        
        for reason, count in sorted(reasons_dict.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_failures) * 100 if total_failures > 0 else 0
            breakdown_data.append({'Ruleset': filter_name,
                                   'Reason': reason,
                                   'Count': count,
                                   'Percentage_of_Filter_Failures': percentage,
                                   'Total_Filter_Failures': total_failures
                                  })
    
    breakdown_df = pd.DataFrame(breakdown_data)
    
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'complete_reasons_breakdown.csv')
    breakdown_df.to_csv(output_path, index=False)
    print(f"Complete reasons breakdown saved to: {output_path}")
    
    return breakdown_df

def create_comprehensive_overview(filter_reasons, filter_failures, file_path):
    all_reasons = {}
    
    for filter_name, reasons in filter_reasons.items():
        for reason, count in reasons:
            if reason in all_reasons:
                all_reasons[reason] += count
            else:
                all_reasons[reason] = count
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)
    
    if not top_reasons:
        return
    
    display_count = min(30, len(top_reasons))
    plt.figure(figsize=(max(16, display_count * 0.6), 16))
    
    reason_names = []
    reason_counts = []
    
    for reason, count in top_reasons[:display_count]:
        if len(reason) > 30:
            reason_short = reason[:27] + "..."
        else:
            reason_short = reason
        reason_names.append(reason_short)
        reason_counts.append(count)
    
    bars = plt.bar(range(len(reason_names)), reason_counts, color='darkgreen', alpha=0.7, width=0.3)
    
    plt.xlabel('Failure Reasons', fontsize=20)
    plt.ylabel('Total Number of Molecules Failed', fontsize=20)
    plt.title(f'Most Common Molecular Filter Failure Reasons (Top {display_count} of {len(top_reasons)})', fontsize=26, fontweight='bold')
    plt.xticks(range(len(reason_names)), reason_names, rotation=45, ha='right', fontsize=16)
    
    for i, (bar, count) in enumerate(zip(bars, reason_counts)):
        plt.text(i, count + max(reason_counts) * 0.01, f"{count}", ha='center', va='bottom', fontsize=14)
    
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'comprehensive_reasons_overview.png')
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    print(f"Comprehensive overview plot saved to: {output_path}")
    
    plt.close()
    
    print("\n" + "="*80)
    print(f"MOST COMMON FAILURE REASONS ACROSS ALL FILTERS (Top 20 of {len(top_reasons)})")
    print("="*80)
    
    for i, (reason, count) in enumerate(top_reasons[:20], 1):
        print(f"{i:3d}. {reason:<60} : {count:5d}")
    
    all_reasons_df = pd.DataFrame(top_reasons, columns=['Reason', 'Total_Count'])
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'all_reasons_summary.csv')
    all_reasons_df.to_csv(output_path, index=False)
    print(f"All reasons summary saved to: {output_path}")

def create_summary_table(filter_failures, filter_reasons, file_path):
    summary_data = []
    
    for filter_name, stats in filter_failures.items():
        row = {'Ruleset': filter_name,
               'Total_Failures': stats['failures'],
               'Failure_Percentage': stats['percentage'],
               'Total_Molecules': stats['total'],
               'Unique_Reasons_Count': len(filter_reasons.get(filter_name, []))
              }
        
        if filter_name in filter_reasons and filter_reasons[filter_name]:
            for i, (reason, count) in enumerate(filter_reasons[filter_name][:5], 1):
                row[f'Top_Reason_{i}'] = reason
                row[f'Top_Reason_{i}_Count'] = count
                row[f'Top_Reason_{i}_Percentage'] = (count / stats['failures']) * 100 if stats['failures'] > 0 else 0
        
        summary_data.append(row)
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('Total_Failures', ascending=False)
    
    path_to_save = os.path.join(os.path.dirname(file_path), 'CommonAlertsBreakdown')
    os.makedirs(path_to_save, exist_ok=True)
    output_path = os.path.join(path_to_save, 'filter_summary_table.csv')
    summary_df.to_csv(output_path, index=False)
    print(f"Summary table saved to: {output_path}")
    
    return summary_df


def create_comparison_plot(directory_path):
    extended_files = []
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file == 'CommonAlerts_extended.csv':
                extended_files.append(os.path.join(root, file))
    
    if len(extended_files) > 1:
        print(f"Found {len(extended_files)} CommonAlerts_extended.csv files:")
        for i, file in enumerate(extended_files):
            print(f"{i+1}. {file}")
        
        choice = input(f"\nEnter the number of the file to analyze (1-{len(extended_files)}): ")
        try:
            selected_file = extended_files[int(choice) - 1]
            return analyze_filter_failures(selected_file)
        except (ValueError, IndexError):
            print("Invalid choice. Using the first file.")
            return analyze_filter_failures(extended_files[0])
    elif len(extended_files) == 1:
        return analyze_filter_failures(extended_files[0])
    else:
        print("No CommonAlerts_extended.csv files found!")
        return None


if __name__ == "__main__":
    base_dir = "results"
    
    if os.path.exists(base_dir):
        result = create_comparison_plot(base_dir)
        if result:
            print("\n" + "="*80)
            print("ANALYSIS COMPLETE!")
            print("="*80)
            print("Generated files:")
            print("- filter_failures_plot.png: Main bar chart")
            print("- [filter_name]_reasons_plot.png: Individual plots for each filter showing ALL reasons")
            print("- all_filters_reasons_plot.png: Multi-panel plot showing all filters with reasons")
            print("- comprehensive_reasons_overview.png: Most common reasons")
            print("- filter_summary_table.csv: Summary statistics")
            print("- complete_reasons_breakdown.csv: Detailed breakdown")
            print("- all_reasons_summary.csv: All reasons ranked by frequency")
    else:
        print(f"Directory {base_dir} not found!")
        print("Please make sure you're running this script from the MolGenBenchmark directory.") 