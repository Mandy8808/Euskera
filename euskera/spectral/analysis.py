"""Spectral analysis implementation."""

import numpy as np

def Organize(data, Rtol=1e-02, Atol=1e-03):
   """
   Organizes the data based on a reference row.
   """
   R_row, i = Reference_row(data)  # Pivot row
   data_T = data.copy()

   Nrow = len(data_T.index)
   for ind in range(Nrow):
      if ind == i:  # Skip the pivot row
         continue
      data_T = Organize_row(ind, R_row, data_T, Rtol, Atol)

   return data_T

def Reference_row(data):
   """
   Gets the reference row (first one without NaN values).
   """
   Nrow = len(data.index)
   for i in range(Nrow):
      R_row = data.iloc[[i]]
      if not R_row.isnull().values.any():
         print(f'The reference row is {i}')
         return R_row, i
   return None, None  # In case all rows contain NaN values

def Organize_row(ind, R_row, data_T, Rtol, Atol):
   """
   Organizes row `ind` based on the reference row `R_row`.
   """
   for col in R_row.columns[1:]:  # Skip the first column if it is an index
      # Get real and imaginary parts from the reference row
      valSupI = np.imag(R_row[col].values[0])
      valSupR = np.real(R_row[col].values[0])

      # Get values from the row being sorted
      tempFrame = data_T.iloc[[ind]]
      # tempArray = tempFrame[col].values  # Extract values as an array
      tempArrayI = np.imag(tempFrame)[0] # np.imag(tempArray)
      tempArrayR = np.real(tempFrame)[0] # np.real(tempArray)

      # Comparison with tolerance
      compI = np.isclose(tempArrayI, valSupI, rtol=Rtol, atol=Atol)
      compR = np.isclose(tempArrayR, valSupR, rtol=Rtol, atol=Atol)
      filt = compI & compR  # Element-wise AND

      if filt.any():
         # Get the matching index
         j = np.where(filt)[0]  # np.where returns a tuple, extract the first element
         # val1 = tempArray[filt]  # Matching values
         val1 = (np.array(tempFrame)[0])[filt]  # Matching values
         if len(val1) > 1:
            val1 = val1[0]  # Take only the first if there are multiple
         val2 = data_T.at[ind, col]  # Original value

         # Swap values
         data_T.at[ind, col] = val1
         #data_T.at[ind, data_T.columns[j[0]]] = val2  # Ensure correct access
         data_T.loc[ind, j] = val2

   return data_T

def sep(Auto_Valores, Auto_Funciones, Imag=False):
   """
   Separando los autovalores y autovectores
   """
   ndatos = len(Auto_Valores)
   dataR, dataI = [], []
   dataRV, dataIV = [], []
   for i in range(ndatos):
        jj = np.real(Auto_Valores[i])!=0  # reales
        dataR.append(Auto_Valores[i][jj])
        dataRV.append(Auto_Funciones[i][:, jj])

        if Imag:
            gg = np.array([not(k) for k in jj])  # imag
            dataI.append(Auto_Valores[i][gg])
            dataIV.append(Auto_Funciones[i][:, gg])

   if Imag:
      return dataR, dataRV, dataI, dataIV
   else:
      return dataR, dataRV
