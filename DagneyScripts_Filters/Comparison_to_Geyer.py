
# import all the things
import sys
import pandas as pd
pd.options.mode.chained_assignment = None
import pandas.util.testing as tm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

fnPalette = ["#D8F000", "#00E7FF","#A5A3FF","#FF4BBE","#8E9BAC","#4F2EEF",
          "#9EEFE1","#FECDD1","#AD2D96","#FF9D42","#BBD532","#808080"]

fnPalette3 = ["#D8F000", "#00E7FF","#FF9D42"]

fnPalette4 = ["#D8F000", "#00E7FF","#FF9D42","#FF4BBE"]

file1="Input/Geyer_deep proteome_curve.csv"
# file2="Blair_EVs/evidence.txt"
# file3="Blair_EVs/SampleSheet-EV.csv"

def read_in_data(file2, file3):
    # read in geyer data, study data, sample sheet data
    geyer_data = pd.read_csv(file1)
    study2_data = pd.read_csv(file2,delimiter = "\t")
    sample_sheet = pd.read_csv(file3)
    
    return geyer_data, study2_data, sample_sheet
    
def add_abundance_rank(df):
    
    df = df.sort_values(by="LFQ Intensity [log10]",ascending=False).reset_index()
    df = df.drop(columns="index").reset_index()
    df = df.rename(columns={"index":"Abundance Rank"})
    
    df["Abundance Level"] = df.apply (lambda row: label_abundance(row), axis=1)
    
    return df

def label_abundance(row):
    
    if row["Abundance Rank"] <= 240:
        return "Low"
    elif 240 < row["Abundance Rank"] <= 480:
        return "Mid"
    elif 480 < row["Abundance Rank"] <= 720:
        return "Deep"
    elif row["Abundance Rank"] > 720:
        return "Ultra Deep"
    return "Outside Range"


def explode_column(df, lst_col):
    
    x = df.assign(**{lst_col:df[lst_col].str.split(';')})
    
    test = pd.DataFrame({
              col:np.repeat(x[col].values, x[lst_col].str.len())
              for col in x.columns.difference([lst_col])
          }).assign(**{lst_col:np.concatenate(x[lst_col].values)})[x.columns.tolist()]

    return test

def add_sample_column(sample_sheet, study_data):
 
    df_list = []
    sample_names = sample_sheet["sample_id"].unique()
    # loop through all samples listed in sample sheet
    for x in sample_names:
        df = study_data.loc[study_data["Raw file"].str.contains(x)]
        new_df = df.assign(Sample=x)
        df_list.append(new_df)

    return pd.concat(df_list)

# add indicator value and change to outer join
def compare_to_geyer_LFQ(sample_sheet, study_data, geyer_data):
    df_list = []
    for x in range(len(sample_sheet)):

        test_sample = sample_sheet.iloc[x,0]

        study_df = add_sample_column(sample_sheet, study_data)
        new_df = study_df.loc[study_df["Sample"] == test_sample]

        merged_df = pd.merge(new_df, geyer_data, how="outer", left_on="Leading razor protein",right_on='Majority Protein IDs',indicator=True)
        merged_df = merged_df.rename(columns={"_merge": "Overlap"})
        merged_df["Overlap"] = merged_df["Overlap"].replace(["left_only","right_only","both"],["Sample Only","Geyer Only","Both"])
        df_list.append(merged_df)
        plt.figure(figsize=(15,5))
        ax = sns.scatterplot(data=merged_df, x="Abundance Rank",y="LFQ Intensity [log10]",hue="Overlap")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        plt.title(f'{test_sample} - LFQ Intensity [log10]')
        plt.xlabel("Geyer Abundance Rank")
        plt.savefig(f'Output/{test_sample} - LFQ Intensity [log10] - all points')
        #plt.show()
        
        only_both = merged_df[merged_df["Overlap"] == "Both"]
        plt.figure(figsize=(10,5))
        ax = sns.scatterplot(data=only_both, x="Abundance Rank",y="LFQ Intensity [log10]")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        plt.title(f'{test_sample} - LFQ Intensity [log10]')
        plt.xlabel("Geyer Abundance Rank")
        plt.savefig(f'Output/{test_sample} - LFQ Intensity [log10]')
        #plt.show()
        
    all_df = pd.concat(df_list)
    all_df_2 = all_df[["Majority Protein IDs","Abundance Rank","Sample","Protein Names","LFQ Intensity [log10]","Overlap","Abundance Level","Leading razor protein"]]
    all_df_2 = all_df_2.drop_duplicates()
    all_df_2.to_csv("Study_vs_Geyer.csv",index=False)
    
    return all_df_2

def shift_values(df):

    sample_list = df["Sample"].unique()
    count = 0
    for i in range(len(sample_list)):
        df["LFQ Intensity [log10]"] = np.where(df["Sample"] == sample_list[i],df["LFQ Intensity [log10]"]+(3*count),df["LFQ Intensity [log10]"])
        count = count + 1
        
    return df

def comparison_charts(df, metadata):
    
    final_meta = pd.merge(df, metadata, left_on="Sample",right_on="sample_id")
    
    df = final_meta[final_meta["Overlap"] == "Both"]
    
    sns.boxplot(data=df, x="sample_name", y="Abundance Rank",color="#DCDCDC")
    plt.xlabel("Samples")
    plt.savefig(f'Output/All Samples - Abundance Box Plot')
    #plt.show()
    
    shifted_df = shift_values(df)
    
    plt.figure(figsize=(10,5))
    sns.scatterplot(data=shifted_df, x="Abundance Rank",y="LFQ Intensity [log10]",hue="sample_name")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,facecolor="white") 
    plt.ylabel("Intensity")   
    plt.xlabel("Geyer Abundance Rank")
    plt.yticks(ticks=[])
    plt.savefig(f'Output/All Samples - Intensity Comparison')
    #plt.show()
    
    plt.figure(figsize=(10,5))
    sns.countplot(data=df, x="Abundance Level",hue="sample_name",order=["Low","Mid","Deep","Ultra Deep"])
    plt.ylabel("Number of IDs") 
    plt.title("Proteome Depth Analysis")
    plt.savefig(f'Output/All Samples - Abundance Levels')
    #plt.show()
    
    
def plot_unique_proteins(data_df,sample_sheet):
    
    df_sample = add_sample_column(sample_sheet, data_df)
    df = pd.merge(df_sample, sample_sheet, left_on="Sample",right_on="sample_id")

    values = []
    sample_names = df["sample_name"].unique()
    for sample in sample_names:
        df1 = df[(df["sample_name"] == sample)]
        values.append(len(df1["Leading razor protein"].unique()))

    new_df = pd.DataFrame({'Sample Name': sample_names, 'Unique IDs': values})

    s = sns.barplot(data=new_df, x="Sample Name",y="Unique IDs")
    
    for p in s.patches:
        s.annotate(format(p.get_height(), '.1f'), 
                       (p.get_x() + p.get_width() / 2., p.get_height()), 
                       ha = 'center', va = 'center', 
                       xytext = (0, 6), 
                       textcoords = 'offset points')

    plt.ylim(0, new_df["Unique IDs"].max()*1.2)
    plt.title("Identified Protein Groups")
    plt.savefig(f'Output/Unique IDs per Sample')

def main():

    if len(sys.argv) == 3:
        file2 = sys.argv[1]
        file3 = sys.argv[2]
    else:
        print("Please provide two files: evidence.txt and sample_sheet.csv")
        exit(0)

    # run functions
    geyer_data, study2_data, sample_sheet = read_in_data(file2, file3)
    geyer_df = add_abundance_rank(geyer_data)
    geyer_explode = explode_column(geyer_df, "Majority Protein IDs")
    plot_unique_proteins(study2_data, sample_sheet)
    final_dataset = compare_to_geyer_LFQ(sample_sheet, study2_data, geyer_explode)
    comparison_charts(final_dataset, sample_sheet)


if __name__ == '__main__':
    main()