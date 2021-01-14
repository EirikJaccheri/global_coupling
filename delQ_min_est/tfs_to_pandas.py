import numpy as np
import tfs
import pandas as pd 

error_file = "output_files/errors.out"
err_dict = tfs.read(error_file,index = "NAME")
err_dict["K1SL"].to_csv("output_files/K1SL.out")
