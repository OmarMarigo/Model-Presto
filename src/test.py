

import pandas as pd

gedi = pd.read_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/src/dta.csv')

gedi = gedi.rename(columns={"x": "lat", "y": "lon"})

sentinel_1 = pd.read_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_1.csv')

sentinel_2 = pd.read_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/sentinel_2.csv')

data = pd.concat([gedi,sentinel_2, sentinel_1], ignore_index=True)
print(data.shape)
#data.to_csv('/Users/clementkm/Documents/School/TOLBI STAGE /PROJECT/Data_Collection_and_Processing/data/CSV/new_dataframe.csv')