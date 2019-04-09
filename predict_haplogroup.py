#!/usr/bin/env python

#!/usr/bin/env python

#-- Diego Montiel 2019
#-- Genetics Identificacion Group @ Erasmus MC  --
#-- Haplogroup prediction


import pandas as pd
import numpy as np
import collections
import operator
import os
from argparse   import ArgumentParser

def get_arguments():

    parser = ArgumentParser(description="ERASMUS MC \n Preddicts Haplogroup from Yleaf/CleanTree output")    
    
    parser.add_argument("-input", "--Input",
        dest="Input", required=True, type=file_exists,
        help="Output file or path produced from Clean tree or Yleaf", metavar="FILE")    
    
    parser.add_argument("-int", "--Intermediate",
        dest="Intermediate", required=True, type=file_exists,
        help="Path with intermediate haplogroups root", metavar="FILE")    
        
    parser.add_argument("-out", "--Outfile",
            dest="Outputfile", required=True,                        
            help="Output file name", metavar="FILE")        
                   
    args = parser.parse_args()    
    return args

def file_exists(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def check_if_folder(path,ext):
    
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if f.endswith(ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]

if __name__ == "__main__":
    
    print("\tErasmus MC Department of Genetic Identification \n\n Haplogroup Prediction")

    args = get_arguments()                
    path_samples = args.Input
    samples = check_if_folder(path_samples,'.out')                                             
    #out_file = "yleaf_all_samples_updated.txt"    
    out_file = args.Outputfile
    #intermediate_tree_table = "Hg_Prediction_tables/Intermediates.txt"
    intermediate_tree_table = args.Intermediate+"Intermediates.txt"
    #path_hg_prediction_tables = "Hg_Prediction_tables/"    
    path_hg_prediction_tables = args.Intermediate
    h_flag = True    
        
    for sample_name in samples:
        
        df_intermediate = pd.read_csv(intermediate_tree_table, header=None, engine='python')
        intermediates = df_intermediate[0].values            
        df_haplogroup = pd.read_csv(sample_name, sep="\t", engine='python')    
        df_haplogroup = df_haplogroup.sort_values(by=['haplogroup'])        
        ## Removes the intermediate haplogroups temporally    
        df_derived = df_haplogroup[df_haplogroup["state"] == "D"]
        ## instance with only D state
        df_tmp = df_derived
        for hg in intermediates:    
            df_tmp = df_tmp.drop(df_tmp[df_tmp.haplogroup == hg].index)
        hg = df_tmp["haplogroup"].values

        """
        Choose the haplogroup based on the highest count in the first letter (E.x. J2a21 = J)
        """
        list_hg = []
        for i in hg:
            list_hg.append(i[0])
        collections_hg = collections.Counter(list_hg)
        init_hg = max(collections_hg.items(), key=operator.itemgetter(1))[0]
        init_hg = init_hg[0]
        tmp_init_hg = init_hg+"_int.txt"
        hg_intermediate_file = path_hg_prediction_tables+tmp_init_hg
        df_intermediate = pd.read_csv(hg_intermediate_file, header=None, sep="\t", engine='python')

        """
        QC.1
        Now you have the putative haplogroup. Open the file with the intermediate hg 
        and show the intermediate haplogroups again and compare if there are some where 
        does not match with the correct state and store it a error. Give the estimate of the  correct_count/total  
        """    
        
        total = 0.0
        correct_state = 0.0
        for i in df_intermediate.values:
            tmp = df_derived.loc[df_derived["haplogroup"] == i[0]]
            if not tmp.empty:        
                if "/" in i[1]:
                    #i = i[1].split("/")[-1] 
                    correct_state += len(tmp) 
                    total += len(tmp)                           
                else:        
                    correct_state += np.sum(i[1] == tmp["state"])                                                
                    total += len(tmp)                  
        #qc_one = round((correct_state / total),3)
        try:
            qc_one = float(round((correct_state / total),3)) if total != 0 else 0
        except ZeroDivisionError as error:                        
            qc_one = 0.0
        """ 
        Removes all haplogroup but the main one
        Check if the preffix of all main haplogroup with D state by allowing one haplogroup that does not match 
        to the main haplogroup until you get back to the origin. 
        If mismatch > 1 skip choose as main haplogroup the following in the largest length and continue. 
        In case haplogroup contains ~ symbol at the end ignore it temporally for preffix comparison but 
        if this is the main haplogroup store the special character at the end. 
        Ancestral and Derive state
        QC.2
        Check if the same name of the main haplogroup appears as an Ancestral State and 
        save the number of count and calculate QC2
        """
        qc_two = 0.0
        qc_tmp = 0.0
        haplogroup_threshold = 0.70
        list_main_hg = []
        for i in hg:    
            if i.startswith(init_hg):
                list_main_hg.append(i.replace("~",""))
        list_main_hg = set(list_main_hg)
        list_main_hg = sorted(list(set(list_main_hg)), reverse=False)        
        list_putative_hg = []        
        for putative_hg in list_main_hg:
            total_qctwo = float(len(df_haplogroup.loc[df_haplogroup["haplogroup"] == putative_hg]))
            Ahg = np.sum("A" == df_haplogroup.loc[df_haplogroup["haplogroup"] == putative_hg]["state"])
            try:
                qc_tmp = float(round((total_qctwo-Ahg)/total_qctwo,3))
            except ZeroDivisionError as error:                            
                qc_tmp = 0.0
            #print(qc_tmp)
            if qc_tmp > haplogroup_threshold:
                list_putative_hg.append(putative_hg)
                #print(putative_hg)
                qc_two = qc_tmp
                #print("{}: {}".format(putative_hg, qc_two))
            #else:
            #    print("{}: {}".format(putative_hg, qc_tmp))
        if len(list_putative_hg) > 0:
            putative_hg = max(list_putative_hg)
        #print("--------------")
        #print("{}: {}".format(putative_hg, qc_two))        

        """
        QC.3
        Show both Ancestral and Derived states from the main haplogroup and check the preffix 
        of each haplogroup that is before the main one until you reach the root. If there are 
        some Ancestral states which follows the pattern/preffix from the main haplogroup keep 
        the count of how many of these appear and this will give the total count. Substract 
        the ones found from the corrected and this will give the QC.3 score. 
        """
        list_main_hg_all = []
        list_hg_all = df_haplogroup["haplogroup"].values
        for i in list_hg_all:    
            if i.startswith(init_hg):
                list_main_hg_all.append(i)
        list_main_hg_all = sorted(list(set(list_main_hg_all)), reverse=True)
        df_main_hg_all = df_haplogroup.loc[df_haplogroup["haplogroup"].isin(list_main_hg_all)]
        df_main_hg_all = df_main_hg_all[["haplogroup","state","marker_name"]]
        df_main_hg_all = df_main_hg_all.sort_values(by="haplogroup", ascending=False).values
        a_match = 0 #ancestral state
        total_match = 0.0
        for i in df_main_hg_all:
            tmp_hg = i[0].replace("~","")
            if tmp_hg in putative_hg:
                total_match +=1        
                if i[1] == "A":
                    a_match += 1
        #qc_three = round((total_match - a_match) / total_match,3)
        try:
            qc_three = float(round((total_match - a_match)/total_match,3)) if total_match != 0 else 0
        except ZeroDivisionError as error:                        
            qc_three = 0.0 

        """
        Haplogroup and marker name
        if main haplogroup is a preffix from a higher resolution Ancestral state haplogroup 
        and report this haplogroup and marker name (Ex. HG: J21a Marker: L123)        

        Could be that that there are more than one contain as a preffix from an ancestral state haplogroup. 
        Should report all of them only if there are different haplogroup name with the resolution
        """
        
        df_putative_ancestral_hg = pd.DataFrame()
        for i in df_main_hg_all:    
            if i[0] == putative_hg:
                break
            if putative_hg in i[0]:     
                tmp = pd.DataFrame([i], columns=["haplogroup","state","marker"])        
                df_putative_ancestral_hg = df_putative_ancestral_hg.append(tmp)        

        putative_ancestral_hg = []        
        if not df_putative_ancestral_hg.empty:
            df_putative_ancestral_hg = df_putative_ancestral_hg.sort_values(by="haplogroup")
            df_putative_ancestral_hg = df_putative_ancestral_hg
            df_putative_ancestral_hg.index = range(len(df_putative_ancestral_hg))    
            for i in range(len(df_putative_ancestral_hg)):        
                row = df_putative_ancestral_hg.iloc[i]
                if putative_ancestral_hg == []:        
                    putative_ancestral_hg.append(row.values)                    
                else:        
                    if putative_ancestral_hg[-1][0] not in row[0]:
                        putative_ancestral_hg.append(row.values)                                  
            putative_ancestral_hg = np.array(putative_ancestral_hg)
        """
        Print outputfile
        """
        out_name = sample_name.split("/")[-1]
        out_name = out_name.split(".")[0]
        for i in df_derived.loc[df_derived["haplogroup"] == putative_hg].values:
            out_hg = ("{}-{}".format(i[3],i[2]))
            break
        header = "Sample_name\tHg\tHg_marker\tQC-score\tQC-1\tQC-2\tQC-3"
        out_marker = out_hg
        if len(putative_ancestral_hg) > 0:
            out_marker += "*("
            for pahg in putative_ancestral_hg:        
                out_marker += pahg[2]+","
            out_marker += ")"
        output = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(out_name,out_hg,out_marker,round((qc_one*qc_two*qc_three),3),qc_one,qc_two,qc_three)
        print(output)
        w_file = open(out_file, "a")    
        if h_flag:                
            h_flag = False
            w_file.write(header)
        w_file.write("\n")
        w_file.write(output)
    w_file.close()
    
    