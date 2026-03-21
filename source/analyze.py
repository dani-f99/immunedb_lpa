from source.helpers import name_cmap
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import source.lpa_master as LPA
from natsort import natsorted
from source import helpers
import seaborn as sns
import pandas as pd
import numpy as np
import regex as re
import json
import os


class lpa_analysis():
    def __init__(self,
                 metric : str = 'substitution_survival',
                 min_treshold : int = 0,
                 keep_zero_freq : bool = True,
                 lpa_input = None,
                 label_dict : dict = None):
    
        """
        metric : str -> Data analysis metric. "substitution_survival" or "trimers_usage".
        lpa_input : pd.DataFrame / str -> string path or dataframe variable of lpa input csv.
        min_treshold: int -> treshold for number of trimers in documents, lower values from this treshold will be dropped.
        keep_zero_freq: -> Keeping trimers with zero occurnces, the LPA algorithm drops them automaticly.
        label_dict : dict ->
        """
        
        # Setting lpa paths
        config =  helpers.read_json()
        self.metric = metric
        self.folder_name = f"{config["sql"]["database"]}-{self.metric}"
        self.lpa_input_path = os.path.join("lpa_input", self.folder_name, "lpa_documents.csv")
        self.lpa_output_path = os.path.join("lpa_results", self.folder_name)

        # Loading lpa input file
        if isinstance(lpa_input, pd.DataFrame):
            self.lpa_input = lpa_input
        
        elif isinstance(lpa_input, str):
            try:
                self.lpa_input = pd.read_csv(lpa_input, index_col=0)
            except:
                print(f"> Invalid lpa input path: {lpa_input}")
        
        else:
           self.lpa_input = pd.read_csv(self.lpa_input_path, index_col=0)


        # Removing documents with 
        if isinstance(min_treshold, int):
            doc_max = self.lpa_input.groupby("document").max().reset_index()
            self.low_docs = doc_max[doc_max.frequency_in_document < min_treshold].document.values
            self.high_docs = doc_max[doc_max.frequency_in_document >= min_treshold].document.values

            self.lpa_input = self.lpa_input[self.lpa_input.document.isin(self.high_docs)]
            print(f"> The following document havn't met the treshold ({str(min_treshold)}):")
            print(doc_max.loc[doc_max.document.isin(self.low_docs), ["document", "frequency_in_document"]])

        if keep_zero_freq:
            unique_trimers =  np.sort(self.lpa_input.element.unique())
            unique_labels = self.lpa_input.document.unique()
            template_df = pd.DataFrame(index=unique_trimers)

            temp_dfs = []
            for label in unique_labels:
                input_df = self.lpa_input[self.lpa_input.document == label].set_index("element")["frequency_in_document"]
                conc_df = pd.concat([template_df, input_df], axis=1).fillna(0).astype("int")
                conc_df["document"] = label
                temp_dfs.append(conc_df)

            self.lpa_input = pd.concat(temp_dfs, axis=0).reset_index(drop=False, names="element")[["document", "element", "frequency_in_document"]]
        

        # Replacing zeros with null values
        self.lpa_input = self.lpa_input.replace(np.inf, np.nan).dropna(axis=0, how="all")

        # Starting to perform the LPA analysis
        # creating corpus and domain 
        corpus = LPA.Corpus(freq=self.lpa_input)
        dvr = corpus.create_dvr()

        # sorting the domain by element value
        dvr_sort = dvr.sort_values(by="element").reset_index(drop=True)

        # defining epsilon for KLDe distance calculation
        epsilon_frac = 2
        epsilon = 1 / (len(dvr) * epsilon_frac)

        # creating distances signatures for each document
        signatures = corpus.create_signatures(epsilon=epsilon, sig_length=None, distance="KLDe")
        sig_list = [i.sort_index() for i in signatures]

        self.sig_df = dvr_sort.copy()["element"].to_frame()

        # Creating dataframe of signatures
        for i in sig_list:
           self.sig_df = self.sig_df.merge(right=i.to_frame(), how="left", left_on="element", right_index=True) 
            
        self.sig_df.set_index(keys="element", drop=True, inplace=True)
        self.sig_df = self.sig_df[np.sort(self.sig_df.columns)]

        # Replacing zeros with null values
        self.x_range = self.sig_df.replace(np.inf, np.nan).dropna(axis=0, how="all").index

    # Heatmap plot: KLDe distances of each dataset from the shared domain
    def heatmap(self,
                mask_input : pd.DataFrame = None,
                save_fig : bool = True) -> plt.Figure:
        
        domain_df = self.sig_df[self.sig_df.index.isin(self.x_range)]
        fhight = domain_df.shape[1]/5
        fig, ax = plt.subplots(1,1, figsize=(20,fhight))


        if mask_input is None:
            mask_df = np.full(domain_df.shape,False)

        try:
            cond_shape = (mask_input.shape == domain_df.shape)
            if cond_shape:
                mask_df = mask_input
            else:
                mask_df = np.full(domain_df.shape,False)
        except:
            mask_df = np.full(domain_df.shape,False)
        
        sns.heatmap(data=domain_df.T,
                    cbar_kws={'label': 'KLDe Distance'},
                    cmap="turbo",
                    xticklabels=True,
                    yticklabels=True,
                    mask=mask_df.T,
                    ax=ax)
        
        cbar_ax = ax.figure.axes[-1]
        cbar_ax.tick_params(labelsize=11)
        cbar_ax.yaxis.get_label().set_fontsize(15)

        ax.tick_params(axis="both", labelsize=11)
        ax.set_xlabel("Position", fontsize=15)
        ax.set_ylabel("Repertoire", fontsize=15)

        range_all = range(1,105) #FR1, CDR1, FR2, CDR2, FR3
        range_cdr = list(range(27,39)) + list(range(56,66)) #CDR1, CDR2
        range_fr = list(range(1,27))+list(range(39,56))+list(range(66,105)) #FR1,FR2,FR3
        for xtick, xcolor in zip(ax.get_xticklabels(), ["red" if i in range_cdr else "blue" for i in domain_df.T.columns]):
            xtick.set_color(xcolor)

        if save_fig == True:
            fig.savefig(os.path.join(self.lpa_output_path,"heatmap.png"), bbox_inches='tight')

    def pca(self,
            plot : bool = True):

            df_input = self.sig_df[self.sig_df.index.isin(self.x_range)].T
            
            pca_object = PCA(n_components=0.95, random_state=0)
            pca_object.fit(df_input)
            pca_output = pca_object.transform(df_input)
            pca_evar = pca_object.explained_variance_ratio_
            self.pca_df = pd.DataFrame(data = pca_output,
                                index = df_input.index,
                                columns = [f"PC{i+1}" for i in range(0,pca_output.shape[1])])
            
            self.pca_var = pd.DataFrame({"PC":self.pca_df.columns, "Variance":pca_evar})
            
            if plot:
                plt.bar(x=self.pca_df.columns, height=pca_evar)
                plt.xlabel("PC")
                plt.ylabel("Explaned Variance")
                plt.show()
                
            return self.pca_df, self.pca_var


    def pca_scatter(self,
                    n_pcs:int = 3,
                    label_index:int = None,
                    save_fig:bool = True):
        """
        n_pcs : int -> number of principle components to present.
        label_index : int -> present the data with different colors based on unique values of the label when splitted by ".", unique index is the
                             location from which the unique values will be selected.
        save_fig : bool -> save the figure.
        """
        
        #
        x = np.array([1, 2, 3, 4])
        y = np.array([1, 0.94, 0.92, 0.90])
        pf = np.polyfit(x, y, 2)
        eq_pf = np.poly1d(pf)

        self.pca_df["label_index"] = "no_index"
        if isinstance(label_index, int):
            self.pca_df["label_index"] = self.pca_df.index.str.split(".").map(lambda X : X[label_index])
        
        pc_range = range(1,n_pcs+1)
        pc_ranges = [(i-1,j-1) for i in pc_range for j in pc_range]
        n_subplots = int(len(pc_ranges)**0.5)
        #print(f"n_pcs: {n_pcs} \npc_range: {pc_range} \npc_ranges: {pc_ranges} \nn_subplots: {n_subplots}")

        fig, axs = plt.subplots(n_subplots,n_subplots, figsize=(n_subplots*5,n_subplots*5))
        unique_labels = self.pca_df["label_index"].unique()

        for rc in pc_ranges:
            pc_x = f"PC{rc[0]+1}"
            pc_y = f"PC{rc[1]+1}"

            for unique_label in unique_labels:
                index_label_cond = (self.pca_df.label_index == unique_label)

                if n_pcs > 1:
                    temp_axs = axs[rc[0], rc[1]]
                else:
                    temp_axs = axs

                temp_axs.scatter(self.pca_df.loc[index_label_cond ,pc_x],
                                        self.pca_df.loc[index_label_cond ,pc_y],
                                        s = n_subplots*20,
                                        alpha = 0.6,
                                        edgecolor="black",
                                        lw = 1,
                                        label=unique_label)
                
                temp_axs.set_xlabel(pc_x, size=10)
                temp_axs.set_ylabel(pc_y, size=10)
                            
        ax_handles, ax_labels = temp_axs.get_legend_handles_labels() 
        fig.legend(handles=ax_handles,
                   labels=ax_labels,
                   loc='outside upper center', 
                   bbox_to_anchor=(0.5, eq_pf(n_pcs)),
                   fontsize=12,
                   ncol=len(unique_labels))
            

        plt.show()