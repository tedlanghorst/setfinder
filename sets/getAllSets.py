""" Script to create sets for all algorithms, using the sets class
"""

# Standard imports
import os
import sys
from pathlib import Path
import json

# Third-party imports
from netCDF4 import Dataset
import numpy as np

# Local imports
try:
    # from sets import Sets
    from nx_sets import Sets
except ImportError:
    # from sets.sets import Sets
    from sets.nx_sets import Sets

def generate_sets(reaches:list, continent:str, output_dir:str, algorithms:list, sword_dataset, sword_filepath:str, expanded:bool):
    """Main function for finding sets"""


    #context

    # try:
    #     shell = get_ipython().__class__.__name__
    #     if shell == 'ZMQInteractiveShell':
    #         print('Running from interatctive shell (e.g. Jupyter notebook) detected. Modifying command line args')
    #         sys.argv=sys.argv[1:]
    #         #print(len(sys.argv))
    #         #print(sys.argv)
    # except NameError:
    #     print("Not running in Jupyter notebook.")

    # if len(sys.argv) <= 2:
    #     try:
    #         index_to_run=int(sys.argv[0]) #integer
    #     except IndexError:
    #         index_to_run=-235
            
    #     continent=sys.argv[1]
        
    #     # specify local directories to run in, for debug purposes
    #     # input_dir=Path('input')
    #     # output_dir=Path('output')        
    #     input_dir=Path('postinput')
    #     output_dir=Path('postoutput')        
        
    # else:
    #     index_to_run=int(args.index)
    #     continent=continent
        
    ## temporary change - for mike's debug march 6
    #index_to_run=sys.argv[1]
    #continent=sys.argv[2]
    

    #data directories
    # INPUT_DIR = input_dir
    # OUTPUT_DIR = output_dir   
    # swordfilepath=INPUT_DIR.joinpath("sword")

    # read in file with all reaches to run
    # reach_json=INPUT_DIR.joinpath(f"reaches_{continent.lower()}.json")
    # with open(reach_json) as json_file:
    #     reaches = json.load(json_file)


        
    # # figure out which sword file to read
    # swordfile=swordfilepath.joinpath(reaches[0]['sword'])

    # # read in sword file
    # sword_dataset=Dataset(swordfile)



    #get set
    # Algorithms=['MetroMan','NeoBAM']
    # Algorithms=['MetroMan','HiVDI','SIC','NeoBAM']
    #Algorithms=['HiVDI']
    # Algorithms=['MetroMan']
    #Algorithms=['SIC']
    # Algorithms=['NeoBAM']

    all_reaches = []
    for Algorithm in algorithms:
        print('Getting set for',Algorithm)
        params = SetParameters(Algorithm, continent, expanded)
        # print(params)
        for key, value in params.items():
            print("    {} : {}.".format(key.capitalize(), value))

        algoset = Sets(params,reaches,sword_dataset)
        InversionSets=algoset.getsets()

        # output to json file        
        all_reaches.extend(algoset.write_inversion_set_data(InversionSets,output_dir, expanded))

    #close sword dataset
    sword_dataset.close()
    
    all_reaches = list(set(all_reaches))
    all_reaches.sort()    
    return all_reaches

def SetParameters(algo, cont, expanded):
    """Seting parameters for setfinder

    Parameters
    ----------
    algo: string
        Algorithm name
    cont: string
        Continent abrevation
    """
    
    """
     Note: 
     
     In some cases, the sword topology causes reaches to be up and downstream from another.
     This causes the the set finder to run infinitly when looking up and downstream for valid set reaches
     We put a limit of 1000 just in case this occures. When it does occur, it only adds the two problematic reaches to the set
     (previously set to np.inf)
    """    
    LargeNumber=1000

    RequireAllReachesInFile=not expanded #usual operation, set this to true, False makes it expand
    # RequireAllReachesInFile=False #in dev set "step 1" set this to false, so we can scrape a list of all reaches in sets

    
    params={}
    params['algo']=algo
    if algo == 'MetroMan':
        params['RequireIdenticalOrbits']=True
        params['DrainageAreaPctCutoff']=10.
        params['AllowRiverJunction']=False
        params['Filename']=f'metrosets_{cont.lower()}.json'
        params['MaximumReachesEachDirection']=2
        #params['MinimumReaches']=3
        params['MinimumReaches']=1
        params['AllowedReachOverlap']=-1 # specify -1 to just remove duplicates
        params['RequireSetReachesInput']=RequireAllReachesInFile # typically set to true: requires all reaches in set to be in input reach list
        # params['']
    elif algo == 'HiVDI':
        params['RequireIdenticalOrbits']=False
        params['DrainageAreaPctCutoff']=30.
        params['AllowRiverJunction']=False
        params['Filename']=f'hivdisets_{cont.lower()}.json'
        """
         In some cases, the sword topology causes reaches to be up and downstream from another.
         This causes the the set finder to run infinitly when looking up and downstream for valid set reaches.
         We put a limit of 1000 just in case this occures. When it does occur, it only adds the two problematic reaches to the set.
         (previously set to np.inf)
        """
        params['MaximumReachesEachDirection']=LargeNumber
        params['MinimumReaches']=1
        params['AllowedReachOverlap']=.5
        params['RequireSetReachesInput']=RequireAllReachesInFile # typically set to true: requires all reaches in set to be in input reach list
    elif algo == 'SIC':
        params['RequireIdenticalOrbits']=False
        params['DrainageAreaPctCutoff']=30.
        params['AllowRiverJunction']=False
        params['Filename']=f'sicsets_{cont.lower()}.json'
        params['MaximumReachesEachDirection']=LargeNumber
        params['MinimumReaches']=1
        params['AllowedReachOverlap']=.67
        params['RequireSetReachesInput']=RequireAllReachesInFile # typically set to true: requires all reaches in set to be in input reach list
    elif algo == 'NeoBAM':
        params['RequireIdenticalOrbits']=False
        params['DrainageAreaPctCutoff']=10.
        params['AllowRiverJunction']=False
        params['Filename']=f'neosets_{cont.lower()}.json'
        params['MaximumReachesEachDirection']=1000
        params['MinimumReaches']=3
        params['AllowedReachOverlap']=0 # specify -1 to just remove duplicates
        params['RequireSetReachesInput']=RequireAllReachesInFile # typically set to true: requires all reaches in set to be in input reach list
        
 
    return params

# if __name__ == "__main__":
#     from datetime import datetime
#     start = datetime.now()
#     main()
#     end = datetime.now()
#     print(f"Execution time: {end - start}")

