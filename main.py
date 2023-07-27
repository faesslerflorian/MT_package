import pandas as pd
import numpy as np
from dynamotable.convention import COLUMN_NAMES
from dynamotable_autofill import read
from dynamotable.dynamotable import write


a = pd.DataFrame({"A": [1, 2], "B": [3, 4]})

b = a.to_numpy(copy=True)
print(b)
b[:,[0,1]]=0
print(b)




